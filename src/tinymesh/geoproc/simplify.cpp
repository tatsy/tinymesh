#define TINYMESH_API_EXPORT
#include "simplify.h"

#include <vector>
#include <queue>
#include <set>
#include <array>

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/QR>
using Matrix4 = Eigen::Matrix4d;
using Vector4 = Eigen::Vector4d;

#include "core/vec.h"
#include "core/progress.h"
#include "core/openmp.h"
#include "polymesh/mesh.h"
#include "polymesh/face.h"
#include "polymesh/halfedge.h"
#include "polymesh/vertex.h"
#include "geoproc/smooth.h"

namespace {

struct UnionFindTree {
    UnionFindTree() = default;

    explicit UnionFindTree(size_t n)
        : n{ n }
        , values{ } {
        values.assign(n, -1);
    }

    int root(int x) {
        if (values[x] < 0) {
            return x;
        }

        values[x] = root(values[x]);
        return values[x];
    }

    void merge(int x, int y) {
        x = root(x);
        y = root(y);
        if (x == y) {
            return;
        }

        values[x] = values[y];
        values[y] = x;
    }

    bool same(int x, int y) {
        return values[x] == values[y];
    }

    size_t n;
    std::vector<int> values;
};

}  // anonymous namespace

namespace tinymesh {

struct QEMNode {
    QEMNode(double value, Halfedge *he, const Vec3 &v)
        : value{ value }
        , he{ he }
        , v{ v } {
    }

    bool operator<(const QEMNode &n) const {
        return value < n.value;
    }

    bool operator>(const QEMNode &n) const {
        return value > n.value;
    }

    double value;
    Halfedge *he;
    Vec3 v;
};

double computeQEM(const Matrix4 &m1, const Matrix4 &m2, const Vertex &v1, const Vertex &v2, Vec3 *v) {
    Matrix4 Q = Matrix4::Identity();
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 4; j++) {
            Q(i, j) = m1(i, j) + m2(i, j);
        }
    }

    Vector4 v_bar;
    const double D = Q.determinant();
    if (D < 1.0e-8) {
        Vec3 m = 0.5 * (v1.pos() + v2.pos());
        v_bar << m[0], m[1], m[2], 1.0;
    } else {
        Vector4 bb;
        bb << 0.0, 0.0, 0.0, 1.0;
        v_bar = Q.colPivHouseholderQr().solve(bb);
    }

    const double qem = (v_bar.transpose() * ((m1 + m2) * v_bar))(0, 0);
    if (v) {
        *v = Vec3(v_bar(0, 0), v_bar(1, 0), v_bar(2, 0));
    }
    return qem;
}

void simplifyIncremental(Mesh &mesh, int numTarget) {
    static const double Eps = 1.0e-12;
    const int numTargetRemove = mesh.num_vertices() - numTarget;
    if (numTarget <= 0) {
        Warn("#vertices is already less than #target: %d < %d", (int)mesh.num_vertices(), numTarget);
        return;
    }

    ProgressBar pbar(numTargetRemove);

    // Pre-smoothing
    laplace_smooth(mesh);

    // Simplification
    int numRemoved = 0;
    while (true) {
        // Current mesh status
        const int numVertices = (int)mesh.num_vertices();
        const int numHalfedges = (int)mesh.num_halfedges();
        const int numFaces = (int)mesh.num_faces();

        // Compute quadric metric tensor for current vertices.
        std::vector<Matrix4> Qs(numVertices, Matrix4::Zero());
        for (int i = 0; i < numFaces; i++) {
            Face *f = mesh.face(i);
            if (f->isBoundary()) {
                continue;
            }

            std::vector<Vertex *> vs;
            for (auto vit = f->v_begin(); vit != f->v_end(); ++vit) {
                vs.push_back(vit.ptr());
            }

            if (vs.size() != 3) {
                Warn("Non trianglar mesh is detected: #vertex = %d", (int)vs.size());
                return;
            }

            Vec3 norm = cross(vs[1]->pos() - vs[0]->pos(), vs[2]->pos() - vs[0]->pos());
            const double w = length(norm);
            norm /= (w + Eps);

            const double nx = norm.x();
            const double ny = norm.y();
            const double nz = norm.z();
            const double d = -dot(norm, vs[0]->pos());

            Matrix4 Kp;
            Kp << nx * nx, nx * ny, nx * nz, nx * d,
                  ny * nx, ny * ny, ny * nz, ny * d,
                  nz * nx, nz * ny, nz * nz, nz * d,
                   d * nx,  d * ny,  d * nz,  d * d;

            Qs[vs[0]->index()] += Kp;
            Qs[vs[1]->index()] += Kp;
            Qs[vs[2]->index()] += Kp;
        }

        // Push QEMs
        std::priority_queue<QEMNode, std::vector<QEMNode>, std::greater<QEMNode>> que;
        for (int i = 0; i < numHalfedges; i++) {
            Halfedge *he = mesh.halfedge(i);

            Vertex *v1 = he->src();
            Vertex *v2 = he->dst();

            int i1 = v1->index();
            int i2 = v2->index();
            Matrix4 &q1 = Qs[i1];
            Matrix4 &q2 = Qs[i2];
            Vec3 v;
            const double qem = computeQEM(q1, q2, *v1, *v2, &v);
            que.push(QEMNode(qem, he, v));
        }

        // Remove halfedges
        std::set<Halfedge *> removedHalfedges;
        while (!que.empty() && numRemoved < numTargetRemove) {
            QEMNode qn = que.top();
            que.pop();

            if (removedHalfedges.find(qn.he) != removedHalfedges.end()) {
                continue;
            }

            const int ii = qn.he->src()->index();
            const int jj = qn.he->dst()->index();
            const Vec3 v_bar = qn.v;

            Vertex *v_i = mesh.vertex(ii);
            Vertex *v_j = mesh.vertex(jj);

            // Skip boundary vertice3s
            if (v_i->isBoundary() || v_j->isBoundary()) {
                continue;
            }

            // Vertices with degree < 4 cannot be contracted
            if (v_i->degree() <= 3 || v_j->degree() <= 3) {
                continue;
            }

            // Check face flip
            bool isFlip = false;
            std::vector<Face*> faces;
            for (auto it = v_i->f_begin(); it != v_i->f_end(); ++it) {
                faces.push_back(it.ptr());
            }

            for (auto it = v_j->f_begin(); it != v_j->f_end(); ++it) {
                faces.push_back(it.ptr());
            }

            for (Face *f : faces) {
                bool has_i = false;
                bool has_j = false;
                std::vector<Vertex*> vs;
                std::vector<Vec3> ps;
                for (auto it = f->v_begin(); it != f->v_end(); ++it) {
                    vs.push_back(it.ptr());
                    ps.push_back(it->pos());
                    if (it.ptr() == v_i) {
                        has_i = true;
                    }
                    if (it.ptr() == v_j) {
                        has_j = true;
                    }
                }

                if (!has_i || !has_j) {
                    // This is contracted face
                    continue;
                }

                const Vec3 n0 = cross(ps[2] - ps[0], ps[1] - ps[0]);
                bool isFound = false;
                for (int i = 0; i < vs.size(); i++) {
                    if (vs[i] == v_i || vs[i] == v_j) {
                        ps[i] = v_bar;
                        isFound = true;
                        break;
                    }
                }

                if (!isFound) {
                    FatalError("Contractible vertex not found!");
                }

                const Vec3 n1 = cross(ps[2] - ps[0], ps[1] - ps[0]);

                const double cos = dot(n0, n1) / (length(n0) * length(n1));
                if (cos <= 1.0e-12) {
                    isFlip = true;
                }
            }

            // Face flips if target vertex is contracted
            if (isFlip) {
                continue;
            }

            // Check face degeneration
            std::set<Vertex*> s_i;
            for (auto it = v_i->v_begin(); it != v_i->v_end(); ++it) {
                s_i.insert(it.ptr());
            }
            std::vector<Vertex*> neighbors;
            for (auto it = v_j->v_begin(); it != v_j->v_end(); ++it) {
                if (s_i.count(it.ptr()) != 0) {
                    neighbors.push_back(it.ptr());
                }
            }

            bool isDegenerate = false;
            for (Vertex *v : neighbors) {
                if (v->degree() < 4) {
                    isDegenerate = true;
                    break;
                }
            }

            if (isDegenerate) {
                continue;
            }

            // Update remove halfedges
            removedHalfedges.insert(qn.he);
            removedHalfedges.insert(qn.he->next());
            removedHalfedges.insert(qn.he->next()->next());
            removedHalfedges.insert(qn.he->rev());
            removedHalfedges.insert(qn.he->rev()->next());
            removedHalfedges.insert(qn.he->rev()->next()->next());

            // Collapse halfedge
            qn.he->src()->setPos(qn.v);
            if (mesh.collapseHE(qn.he)) {
                numRemoved += 1;
                pbar.step();
            }
        }

        // Flip edges
        for (int i = 0; i < mesh.num_halfedges(); i++) {
            Halfedge *he = mesh.halfedge(i);
            if (he->face()->isBoundary() || he->rev()->face()->isBoundary()) {
                continue;
            }

            Vertex *v0 = he->src();
            Vertex *v1 = he->dst();
            Vertex *v2 = he->next()->dst();
            Vertex *v3 = he->rev()->next()->dst();
            const int d0 = v0->degree();
            const int d1 = v1->degree();
            const int d2 = v2->degree();
            const int d3 = v3->degree();

            const int score = std::abs(d0 - 6) + std::abs(d1 - 6) + std::abs(d2 - 6) + std::abs(d3 - 6);
            const int after = std::abs(d0 - 1 - 6) + std::abs(d1 - 1 - 6) + std::abs(d2 + 1 - 6) + std::abs(d3 + 1 - 6);
            if (score > after) {
                mesh.flipHE(he);
            }
        }

        // Smoothing
        laplace_smooth(mesh);

        if (numRemoved >= numTargetRemove) {
            break;
        }
    }
    printf("\n");

    Info("%d vertices removed (%6.2f%% achievement)", numRemoved, 100.0 * numRemoved / numTargetRemove);
}

}  // namespace tinymesh

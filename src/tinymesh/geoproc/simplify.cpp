#define TINYMESH_API_EXPORT
#include "simplify.h"

#include <vector>
#include <queue>
#include <set>
#include <array>

#include "core/vec.h"
#include "math/matrix.h"
#include "trimesh/mesh.h"
#include "trimesh/face.h"
#include "trimesh/halfedge.h"
#include "trimesh/vertex.h"

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
    QEMNode(double value, Halfedge *he, const Vec &v)
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
    Vec v;
};

double computeQEM(const Matrix &m1, const Matrix &m2, const Vertex &v1, const Vertex &v2, Vec *v) {
    Matrix Q = Matrix::identity(4, 4);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 4; j++) {
            Q.set(i, j, m1.get(i, j) + m2.get(i, j));
        }
    }

    Matrix v_bar;
    const double D = Q.det();
    if (D < 1.0e-8) {
        Vec m = 0.5 * (v1.pos() + v2.pos());
        double elems[4] = {m.x, m.y, m.z, 1.0};
        v_bar = Matrix(4, 1, elems);
    } else {
        double elems[4] = {0.0, 0.0, 0.0, 1.0};
        v_bar = Q.solve(Matrix(4, 1, elems));
    }

    const double qem = (v_bar.T() * ((m1 + m2) * v_bar)).get(0, 0);
    if (v) {
        *v = Vec(v_bar.get(0, 0), v_bar.get(1, 0), v_bar.get(2, 0));
    }
    return qem;
}

void simplify(Mesh &mesh, double ratio, int nRemain) {
    static const double Eps = 1.0e-12;
    const int nv = (int)mesh.num_vertices();

    // Determine #vertices to be removed
    int nRemove = (int)(nv * (1.0 - ratio));
    if (nRemain > 0) {
        if (nRemain <= 4) {
            FatalError("Remaining vertices must be more than 3!");
        }
        nRemove = nv - nRemain;
    }
    printf("#remove: %d\n", nRemove);

    // Compute matrix Q
    std::vector<Matrix> Qs(nv, Matrix::zeros(4, 4));
    for (int i = 0; i < mesh.num_faces(); i++) {
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

        Vec norm = cross(vs[1]->pos() - vs[0]->pos(), vs[2]->pos() - vs[0]->pos());
        const double w = length(norm);
        norm /= (w + Eps);

        const double nx = norm.x;
        const double ny = norm.y;
        const double nz = norm.z;
        const double d = -dot(norm, vs[0]->pos());

        double elems[] = {
            nx * nx, nx * ny, nx * nz, nx * d,
            ny * nx, ny * ny, ny * nz, ny * d,
            nz * nx, nz * ny, nz * nz, nz * d,
             d * nx,  d * ny,  d * nz,  d * d,
        };
        Matrix Kp = Matrix(4, 4, elems);

        Qs[vs[0]->index()] += Kp;
        Qs[vs[1]->index()] += Kp;
        Qs[vs[2]->index()] += Kp;
    }

    // Push QEMs
    std::priority_queue<QEMNode, std::vector<QEMNode>, std::greater<QEMNode>> que;
    for (int i = 0; i < mesh.num_halfedges(); i++) {
        Halfedge *he = mesh.halfedge(i);

        Vertex *v1 = he->src();
        Vertex *v2 = he->dst();
        int i1 = v1->index();
        int i2 = v2->index();
        Matrix &q1 = Qs[i1];
        Matrix &q2 = Qs[i2];
        Vec v;
        const double qem = computeQEM(q1, q2, *v1, *v2, &v);
        que.push(QEMNode(qem, he, v));
    }

    // Remove halfedges
    int removed = 0;
    auto uftree = UnionFindTree(nv);
    std::set<Halfedge *> removedHalfedges;
    while (removed < nRemove && !que.empty()) {
        QEMNode qn = que.top();
        que.pop();

        if (removedHalfedges.find(qn.he) != removedHalfedges.end()) {
            continue;
        }

        const int ii = qn.he->src()->index();
        const int jj = qn.he->src()->index();
        const Vec v_bar = qn.v;

        Vertex *v_i = mesh.vertex(ii);
        Vertex *v_j = mesh.vertex(jj);

        // Either end vertex is already contracted
        if (!v_i || !v_j) {
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
            std::vector<Vec> ps;
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

            const Vec n0 = cross(ps[2] - ps[0], ps[1] - ps[0]);
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

            const Vec n1 = cross(ps[2] - ps[0], ps[1] - ps[0]);

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
            removed += 1;
            uftree.merge(ii, jj);    
        }
    }

    printf("%d vertices removed!\n", removed);

    for (int it = 0; it < 3; it++) {
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

        // Volonoi tessellation
        for (int l = 0; l < 5; l++) {
            int index;
            const int nv = mesh.num_vertices();
            std::vector<Vec> centroids(nv);
            std::vector<Vec> normals(nv);

            // Compute centroids and tangent planes
            index = 0;
            for (int i = 0; i < mesh.num_vertices(); i++, index++) {
                Vertex *v = mesh.vertex(i);

                // Collect surrounding vertices
                Vec org = v->pos();
                std::vector<Vec> pts;
                for (auto vit = v->v_begin(); vit != v->v_end(); ++vit) {
                    pts.push_back(vit->pos());
                }

                // Compute centroids, tangents, and binormals
                Vec cent(0.0);
                Vec norm(0.0);
                for (int i = 0; i < pts.size(); i++) {
                    const int j = (i + 1) % pts.size();
                    Vec e1 = pts[i] - org;
                    Vec e2 = pts[j] - org;
                    Vec g = (org + pts[i] + pts[j]) / 3.0;

                    cent += g;
                    norm += cross(e1, e2);
                }

                cent /= pts.size();
                const double l = length(norm);

                if (l != 0.0) {
                    centroids[index] = cent;
                    normals[index] = norm / l;
                }
            }

            // Update vertex positions
            index = 0;
            for (int i = 0; i < mesh.num_vertices(); i++, index++) {
                Vertex *v = mesh.vertex(i);
                if (v->isBoundary()) {
                    continue;
                }

                if (length(normals[index]) != 0.0) {
                    const Vec pt = v->pos();
                    Vec e = centroids[index] - pt;
                    e -= normals[index] * dot(e, normals[index]);
                    v->setPos(pt + e);
                }
            }
        }
    }
}

}  // namespace tinymesh

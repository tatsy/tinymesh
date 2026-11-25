#include "remesh.h"

#include <cstdint>
#include <vector>
#include <queue>
#include <set>
#include <array>
#include <unordered_set>

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/QR>
using Matrix4 = Eigen::Matrix4d;
using Vector4 = Eigen::Vector4d;

#include "core/debug.h"
#include "core/vec.h"
#include "core/progress.h"
#include "core/openmp.h"
#include "core/mesh.h"
#include "core/face.h"
#include "core/halfedge.h"
#include "core/vertex.h"
#include "filters/filters.h"

namespace {

struct UnionFindTree {
    UnionFindTree() = default;

    explicit UnionFindTree(size_t n)
        : n{ n }
        , values{} {
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
    QEMNode(double value, std::uintptr_t key, const Vec3 &v)
        : value{ value }
        , key{ key }
        , v{ v } {
    }

    bool operator<(const QEMNode &n) const {
        return value < n.value;
    }

    bool operator>(const QEMNode &n) const {
        return value > n.value;
    }

    double value;
    std::uintptr_t key;
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
        //v_bar = Q.colPivHouseholderQr().solve(bb);
        v_bar = Q.partialPivLu().solve(bb);
    }

    const double qem = (v_bar.transpose() * ((m1 + m2) * v_bar))(0, 0);
    if (v) {
        *v = Vec3(v_bar(0, 0), v_bar(1, 0), v_bar(2, 0));
    }
    return qem;
}

void simplifyQEM(Mesh &mesh, int numTarget, int maxTrials, bool verbose) {
    static const double Eps = 1.0e-12;
    const int numTargetRemove = (int)mesh.numFaces() - numTarget;
    if (numTarget == 0) {
        Assertion(numTarget > 0, "Target face number is equal or less than 0!");
    }

    if (numTargetRemove <= 0) {
        Warn("#vertices is already less than #target: {} < {}", (int)mesh.numFaces(), numTarget);
        return;
    }

    ProgressBar pbar(numTargetRemove);

    // Pre-smoothing
    smoothLaplacian(mesh);

    // Simplification
    int numRemoved = 0;
    for (int trial = 0; trial < maxTrials; trial++) {
        // Current mesh status
        const int numVertices = (int)mesh.numVertices();
        const int numHalfedges = (int)mesh.numHalfedges();
        const int numFaces = (int)mesh.numFaces();
        std::unordered_set<Vertex *> checked;

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
                Warn("Non trianglar mesh is detected: #vertex = {}", (int)vs.size());
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
            Kp << nx * nx, nx * ny, nx * nz, nx * d, ny * nx, ny * ny, ny * nz, ny * d, nz * nx, nz * ny, nz * nz,
                nz * d, d * nx, d * ny, d * nz, d * d;

            Qs[vs[0]->index()] += Kp;
            Qs[vs[1]->index()] += Kp;
            Qs[vs[2]->index()] += Kp;
        }

        // Push QEMs
        std::priority_queue<QEMNode, std::vector<QEMNode>, std::greater<QEMNode>> que;
        std::unordered_map<std::uintptr_t, Halfedge *> activeHEs;
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

            const auto key = reinterpret_cast<std::uintptr_t>(he);
            que.push(QEMNode(qem, key, v));
            activeHEs.emplace(key, he);
        }

        // Remove halfedges
        std::set<Halfedge *> removedHalfedges;
        while (!que.empty() && numRemoved < numTargetRemove) {
            QEMNode qn = que.top();
            que.pop();

            const auto it = activeHEs.find(qn.key);
            if (it == activeHEs.end()) {
                continue;
            }

            Halfedge *he = it->second;
            Vertex *v_i = he->src();
            Vertex *v_j = he->dst();

            // Skip checked vertices
            if (checked.count(v_i) != 0 || checked.count(v_j) != 0) {
                continue;
            }

            // Skip boundary vertices
            if (v_i->isBoundary() || v_j->isBoundary() || v_i->isLocked() || v_j->isLocked()) {
                continue;
            }

            // Vertices with degree < 4 cannot be contracted
            if (v_i->degree() <= 3 || v_j->degree() <= 3) {
                continue;
            }

            // Check face degeneration
            std::set<Vertex *> s_i;
            for (auto it = v_i->v_begin(); it != v_i->v_end(); ++it) {
                s_i.insert(it.ptr());
            }
            std::vector<Vertex *> neighbors;
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

            // Collect removed halfedges
            const std::array<std::uintptr_t, 6> keyToRemove = {
                reinterpret_cast<std::uintptr_t>(he),
                reinterpret_cast<std::uintptr_t>(he->next()),
                reinterpret_cast<std::uintptr_t>(he->next()->next()),
                reinterpret_cast<std::uintptr_t>(he->rev()),
                reinterpret_cast<std::uintptr_t>(he->rev()->next()),
                reinterpret_cast<std::uintptr_t>(he->rev()->next()->next()),
            };

            if (const auto kept = mesh.collapseHE(he, qn.v)) {
                for (const auto key : keyToRemove) {
                    activeHEs.erase(key);
                }

                checked.insert(*kept);
                numRemoved += 2;
                if (verbose) pbar.step(2);
            }

            if (numRemoved >= numTargetRemove) {
                break;
            }
        }

        // Flip edges
        for (int i = 0; i < (int)mesh.numHalfedges(); i++) {
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
        smoothTaubin(mesh);

        if (numRemoved >= numTargetRemove) {
            break;
        }

        if (numFaces > mesh.numFaces()) {
            trial = 0;
        }
    }

    if (verbose) pbar.finish();

    if (verbose) {
        Info("{} faces removed ({:6.2f}% achievement)", numRemoved, 100.0 * numRemoved / numTargetRemove);
    }
}

}  // namespace tinymesh

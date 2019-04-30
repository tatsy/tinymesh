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

    explicit UnionFindTree(int n)
        : n{ n }
        , values{ n, -1 } {
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

    int n;
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

    // How many vertices are removed?
    int nRemove = (int)(nv * (1.0 - ratio));
    if (nRemain > 0) {
        if (nRemain <= 3) {
            FatalError("Remaining vertices must be more than 3!");
        }
        nRemove = nv - nRemain;
    }
    nRemove = nv - nRemain;

    // Compute matrix Q
    std::vector<Matrix> Qs(nv, Matrix::zeros(4, 4));
    for (auto fit = mesh.f_begin(); fit != mesh.f_end(); ++fit) {
        std::vector<Vertex *> vs;
        for (auto vit = fit->v_begin(); vit != fit->v_end(); ++fit) {
            vs.push_back(vit.ptr());
        }

        if (vs.size() != 3) {
            Warn("Non trianglar mesh is detected!");
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

        Matrix Q = Matrix(4, 4, elems);

        Qs[vs[0]->index()] += w * Q;
        Qs[vs[1]->index()] += w * Q;
        Qs[vs[2]->index()] += w * Q;
    }

    // Push QEMs
    std::priority_queue<QEMNode> que;
    for (auto it = mesh.he_begin(); it != mesh.he_end(); ++it) {
        Vertex *v1 = it->src();
        Vertex *v2 = it->dst();
        int i1 = v1->index();
        int i2 = v2->index();
        Matrix &q1 = Qs[i1];
        Matrix &q2 = Qs[i2];
        Vec v;
        const double qem = computeQEM(q1, q2, *v1, *v2, &v);
        que.push(QEMNode(qem, it.ptr(), v));
    }

    int removed = 0;
    auto uftree = UnionFindTree(nv);
    while (removed < nRemove) {
        QEMNode qn = que.top();
        que.pop();

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

        // Collapse halfedge
        mesh.collapseHE(qn.he);
        uftree.merge(ii, jj);

        Assertion(ii == uftree.root(jj), "Merge operation is not succeeded!");
        Assertion(mesh.vertex(ii) != nullptr, "Invalid null vertex is found!");
        Assertion(mesh.vertex(jj) == nullptr, "Invalid non-null vartex is found!");

        // Check triangle shapes
        //bool isUpdate = true;
        //update_vertices = [ v_i ] #list(chain([ v_i ], v_i.vertices()))
        //while is_update:
        //    is_update = False
        //for he in v_i.halfedges():
        //if he.face is None or he.opposite.face is None:
    }

//# Boundary halfedge
//    continue
//
//    v0 = he.next.vertex_to.position
//    v1 = he.vertex_to.position
//    v2 = he.vertex_from.position
//    v3 = he.opposite.next.vertex_to.position
//
//    e0 = v1 - v0
//    e1 = v2 - v0
//    c1 = e0.dot(e1) / (e0.norm() * e1.norm() + EPS)
//    a1 = math.acos(max(-1.0, min(c1, 1.0)))
//
//    e2 = v1 - v3
//    e3 = v2 - v3
//    c2 = e2.dot(e3) / (e2.norm() * e3.norm() + EPS)
//    a2 = math.acos(max(-1.0, min(c2, 1.0)))
//
//    if a1 + a2 > math.pi:
//    mesh.flip_halfedge(he)
//    is_update = True
//    break
//
//# Update matrix Q
//    for v in update_vertices:
//    Qs[v.index] = np.zeros((4, 4))
//    for f in v.faces():
//    vs = list(f.vertices())
//    assert len(vs) == 3
//
//    ps = [ v.position for v in vs ]
//    norm = (ps[1] - ps[0]).cross(ps[2] - ps[0])
//    w = norm.norm()
//    norm /= (w + EPS)
//
//    d = -norm.dot(ps[0])
//    pp = np.array([ norm.x, norm.y, norm.z, d ])
//    Q = pp.reshape((pp.size, 1)) * pp
//    Qs[v.index] += w * Q
//
//# Update QEMs
//    for v1 in update_vertices:
//    for v2 in v1.vertices():
//    assert v1.index != v2.index
//    assert mesh.vertices[v1.index] is not None
//    assert mesh.vertices[v2.index] is not None
//
//    if v1.degree() <= 3 or v2.degree() <= 3: continue
//    if v1.degree() >= 7 or v2.degree() >= 7: continue
//
//    Q1 = Qs[v1.index]
//    Q2 = Qs[v2.index]
//    qem, v_bar = compute_QEM(Q1, Q2, v1, v2)
//
//    pque.push(QEMNode(qem, v1.index, v2.index, v_bar))
//

}

}  // namespace tinymesh

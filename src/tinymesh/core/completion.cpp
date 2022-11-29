#define TINYMESH_API_EXPORT
#include "mesh.h"

#include <queue>
#include <unordered_map>

#include "types.h"
#include "eigen.h"
#include "vertex.h"
#include "face.h"
#include "edge.h"
#include "halfedge.h"

namespace tinymesh {

void Mesh::holeFillMinDihedral(Face *face, double dihedralBound) {
    /*
     * The triangulation algorithm is based on the following papers.
     * Barequet and Sharir, "Filling Gaps in the Boundary of a Polyhedron", 1995.
     * Liepa, "Filling Holes in Meshes", 2003.
     * 
     * Different from them, the triangularion is performed when the dihedral
     * angle in resulting triangularion is less than "dihedralBound".
     */
    using E = IndexPair;
    using T = std::tuple<uint32_t, uint32_t, uint32_t>;

    std::vector<Vertex *> vs;
    std::unordered_map<E, Halfedge *> pair2he;
    std::unordered_map<Halfedge *, E> he2pair;

    // List all boundary vertices
    for (auto it = face->v_begin(); it != face->v_end(); ++it) {
        vs.push_back(it.ptr());
    }
    const int n_verts = (int)vs.size();

    int count = 0;
    for (auto it = face->he_begin(); it != face->he_end(); ++it) {
        const int next = (count + 1) % n_verts;
        pair2he.insert(std::make_pair(E(count, next), it.ptr()));
        he2pair.insert(std::make_pair(it.ptr(), E(count, next)));
        count += 1;
    }

    Eigen::MatrixXd Warea(n_verts, n_verts);
    Eigen::MatrixXd Wangle(n_verts, n_verts);
    Eigen::MatrixXi O(n_verts, n_verts);
    Warea.setConstant(1.0e12);
    Wangle.setConstant(1.0e12);
    O.setConstant(-1);
    for (int i = 0; i < n_verts - 1; i++) {
        Warea(i, i + 1) = 0.0;
        Wangle(i, i + 1) = 0.0;
        O(i, i + 1) = -1;
        if (i < n_verts - 2) {
            const Vec3 v0 = vs[i]->pos();
            const Vec3 v1 = vs[i + 1]->pos();
            const Vec3 v2 = vs[i + 2]->pos();
            Warea(i, i + 2) = 0.5 * length(cross(v1 - v0, v2 - v0));

            const auto he01 = pair2he[E(i, i + 1)];
            const auto he12 = pair2he[E(i + 1, i + 2)];

            const Vec3 v01rev = he01->rev()->next()->dst()->pos();
            const Vec3 v12rev = he12->rev()->next()->dst()->pos();
            const double a01 = Pi - dihedral(v2, v0, v1, v01rev);
            const double a12 = Pi - dihedral(v0, v1, v2, v12rev);
            Wangle(i, i + 2) = std::max(a01, a12);

            O(i, i + 2) = i + 1;
        }
    }

    for (int j = 3; j < n_verts; j++) {
        for (int i = 0; i < n_verts - j; i++) {
            const int k = i + j;
            for (int m = i + 1; m < k; m++) {
                const Vec3 v0 = vs[i]->pos();
                const Vec3 v1 = vs[m]->pos();
                const Vec3 v2 = vs[k]->pos();
                const double F = 0.5 * length(cross(v1 - v0, v2 - v0));
                const double WareaNew = Warea(i, m) + Warea(m, k) + F;

                Vec3 v01rev;
                if (abs(m - i) == 1) {
                    const auto he01 = pair2he[E(i, m)];
                    v01rev = he01->rev()->next()->dst()->pos();
                } else {
                    v01rev = vs[O(i, m)]->pos();
                }

                Vec3 v12rev;
                if (abs(k - m) == 1) {
                    const auto he12 = pair2he[E(m, k)];
                    v12rev = he12->rev()->next()->dst()->pos();
                } else {
                    v12rev = vs[O(m, k)]->pos();
                }

                const double a01 = Pi - dihedral(v2, v0, v1, v01rev);
                const double a12 = Pi - dihedral(v0, v1, v2, v12rev);
                const double D = std::max(a01, a12);
                const double WangleNew = std::max(D, std::max(Wangle(i, m), Wangle(m, k)));

                if (WangleNew < Wangle(i, k) ||
                    (abs(WangleNew - Wangle(i, k)) < dihedralBound && WareaNew < Warea(i, k))) {
                    Warea(i, k) = WareaNew;
                    Wangle(i, k) = WangleNew;
                    O(i, k) = m;
                }
            }
        }
    }

    std::queue<E> que;
    std::vector<T> tris;

    que.push(E(0, n_verts - 1));
    while (!que.empty()) {
        const E ik = que.front();
        que.pop();

        const int i = std::get<0>(ik);
        const int k = std::get<1>(ik);
        const int m = O(i, k);

        if (i + 2 == k) {
            tris.emplace_back(i, m, k);
        } else {
            if (i + 1 != m) {
                que.push(E(i, m));
            }

            tris.emplace_back(i, m, k);

            if (m != k - 1) {
                que.push(E(m, k));
            }
        }
    }

    std::vector<Face *> new_faces;
    for (const auto &t : tris) {
        const int i0 = std::get<0>(t);
        const int i1 = std::get<1>(t);
        const int i2 = std::get<2>(t);

        Halfedge *new_he0 = nullptr;
        if (pair2he.count(E(i0, i1)) != 0) {
            new_he0 = pair2he[E(i0, i1)];
        } else {
            new_he0 = new Halfedge();
            pair2he.insert(std::make_pair(E(i0, i1), new_he0));
            he2pair.insert(std::make_pair(new_he0, E(i0, i1)));
            addHalfedge(new_he0);
        }

        Halfedge *new_he1 = nullptr;
        if (pair2he.count(E(i1, i2)) != 0) {
            new_he1 = pair2he[E(i1, i2)];
        } else {
            new_he1 = new Halfedge();
            pair2he.insert(std::make_pair(E(i1, i2), new_he1));
            he2pair.insert(std::make_pair(new_he1, E(i1, i2)));
            addHalfedge(new_he1);
        }

        Halfedge *new_he2 = nullptr;
        if (pair2he.count(E(i2, i0)) != 0) {
            new_he2 = pair2he[E(i2, i0)];
        } else {
            new_he2 = new Halfedge();
            pair2he.insert(std::make_pair(E(i2, i0), new_he2));
            he2pair.insert(std::make_pair(new_he2, E(i2, i0)));
            addHalfedge(new_he2);
        }

        new_he0->src_ = vs[i0];
        new_he1->src_ = vs[i1];
        new_he2->src_ = vs[i2];

        new_he0->next_ = new_he1;
        new_he1->next_ = new_he2;
        new_he2->next_ = new_he0;

        auto new_face = new Face();
        addFace(new_face);
        new_faces.push_back(new_face);

        new_he0->face_ = new_face;
        new_he1->face_ = new_face;
        new_he2->face_ = new_face;
        new_face->halfedge_ = new_he0;
    }

    for (auto f : new_faces) {
        Halfedge *he = f->halfedge_;
        do {
            if (he->rev_ == NULL) {
                auto e = he2pair[he];
                int i = std::get<0>(e);
                int j = std::get<1>(e);

                auto r = E(j, i);
                if (pair2he.count(r) != 0) {
                    Halfedge *rev = pair2he[r];
                    he->rev_ = rev;
                    rev->rev_ = he;

                    Edge *new_edge = new Edge();
                    he->edge_ = new_edge;
                    rev->edge_ = new_edge;
                    new_edge->halfedge_ = he;
                    addEdge(new_edge);
                } else {
                    printf("%d %d %p\n", i, j, he);
                }
            }

            he->src_->isBoundary_ = false;
            he = he->next_;
        } while (he != f->halfedge_);
    }

    removeFace(face);
}

void Mesh::holeFillAdvancingFront(Face *f) {
}

}  // namespace tinymesh

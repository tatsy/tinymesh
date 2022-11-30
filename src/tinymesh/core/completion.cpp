#define TINYMESH_API_EXPORT
#include "mesh.h"

#include <utility>
#include <queue>
#include <unordered_set>

#include "types.h"
#include "eigen.h"
#include "vertex.h"
#include "face.h"
#include "edge.h"
#include "halfedge.h"

using E = IndexPair;
using T = std::tuple<uint32_t, uint32_t, uint32_t>;

namespace tinymesh {

Face *Mesh::addNewTriangle(const std::vector<Vertex *> &boundary, const std::tuple<uint32_t, uint32_t, uint32_t> &tri,
                           std::unordered_map<IndexPair, Halfedge *> &pair2he,
                           std::unordered_map<Halfedge *, IndexPair> &he2pair) {
    const int i0 = std::get<0>(tri);
    const int i1 = std::get<1>(tri);
    const int i2 = std::get<2>(tri);

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

    new_he0->src_ = boundary[i0];
    new_he1->src_ = boundary[i1];
    new_he2->src_ = boundary[i2];

    new_he0->next_ = new_he1;
    new_he1->next_ = new_he2;
    new_he2->next_ = new_he0;

    new_he0->isBoundary_ = false;
    new_he1->isBoundary_ = false;
    new_he2->isBoundary_ = false;

    auto new_face = new Face();
    addFace(new_face);

    new_he0->face_ = new_face;
    new_he1->face_ = new_face;
    new_he2->face_ = new_face;
    new_face->halfedge_ = new_he0;

    return new_face;
}

void Mesh::holeFillMinDihedral_(Face *face, double dihedralBound) {
    // Collect boundary vertices
    std::vector<Vertex *> boundary;
    for (auto it = face->v_begin(); it != face->v_end(); ++it) {
        boundary.push_back(it.ptr());
    }
    const int n_verts = (int)boundary.size();

    // Collect boundary halfedges
    std::unordered_map<E, Halfedge *> pair2he;
    std::unordered_map<Halfedge *, E> he2pair;
    int count = 0;
    for (auto it = face->he_begin(); it != face->he_end(); ++it) {
        const int next = (count + 1) % n_verts;
        pair2he.insert(std::make_pair(E(count, next), it.ptr()));
        he2pair.insert(std::make_pair(it.ptr(), E(count, next)));
        count++;
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
            const Vec3 v0 = boundary[i]->pos();
            const Vec3 v1 = boundary[i + 1]->pos();
            const Vec3 v2 = boundary[i + 2]->pos();
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
                const Vec3 v0 = boundary[i]->pos();
                const Vec3 v1 = boundary[m]->pos();
                const Vec3 v2 = boundary[k]->pos();
                const double F = 0.5 * length(cross(v1 - v0, v2 - v0));
                const double WareaNew = Warea(i, m) + Warea(m, k) + F;

                Vec3 v01rev;
                if (abs(m - i) == 1) {
                    const auto he01 = pair2he[E(i, m)];
                    v01rev = he01->rev()->next()->dst()->pos();
                } else {
                    v01rev = boundary[O(i, m)]->pos();
                }

                Vec3 v12rev;
                if (abs(k - m) == 1) {
                    const auto he12 = pair2he[E(m, k)];
                    v12rev = he12->rev()->next()->dst()->pos();
                } else {
                    v12rev = boundary[O(m, k)]->pos();
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
        Face *new_face = addNewTriangle(boundary, t, pair2he, he2pair);
        new_faces.push_back(new_face);
    }

    for (auto f : new_faces) {
        for (auto it = f->he_begin(); it != f->he_end(); ++it) {
            if (it->rev_ == NULL) {
                auto e = he2pair[it.ptr()];
                int i = std::get<0>(e);
                int j = std::get<1>(e);

                auto r = E(j, i);
                if (pair2he.count(r) != 0) {
                    Halfedge *rev = pair2he[r];
                    it->rev_ = rev;
                    rev->rev_ = it.ptr();

                    Edge *new_edge = new Edge();
                    it->edge_ = new_edge;
                    rev->edge_ = new_edge;
                    new_edge->halfedge_ = it.ptr();
                    addEdge(new_edge);

                    it->isBoundary_ = false;
                    rev->isBoundary_ = false;
                } else {
                    printf("%d %d %p\n", i, j, it.ptr());
                }
            }
        }
    }

    removeFace(face);
}

void Mesh::holeFillAdvancingFront_(Face *face) {
    for (int hoge = 0; hoge < 50; hoge++) {
        // Collect boundary vertices
        std::vector<Vertex *> boundary;
        for (auto it = face->v_begin(); it != face->v_end(); ++it) {
            boundary.push_back(it.ptr());
        }
        const int N = static_cast<int>(boundary.size());
        printf("N = %d\n", N);

        // Collect boundary halfedges
        std::unordered_map<E, Halfedge *> pair2he;
        std::unordered_map<Halfedge *, E> he2pair;
        int count = 0;
        for (auto it = face->he_begin(); it != face->he_end(); ++it) {
            const int next = (count + 1) % N;
            pair2he.insert(std::make_pair(E(count, next), it.ptr()));
            he2pair.insert(std::make_pair(it.ptr(), E(count, next)));
            count++;
        }

        // Search vertex with minimum open angle
        double minAngle = 1.0e20;
        int minIndex = -1;
        for (int i = 0; i < N; i++) {
            const int prev = (i - 1 + N) % N;
            const int next = (i + 1) % N;
            const Vec3 norm1 = pair2he[E(prev, i)]->rev()->face()->normal();
            const Vec3 norm2 = pair2he[E(i, next)]->rev()->face()->normal();
            const Vec3 norm = normalize(norm1 + norm2);

            const Vec3 e1 = boundary[prev]->pos() - boundary[i]->pos();
            const Vec3 e2 = boundary[next]->pos() - boundary[i]->pos();
            const Vec3 outer = cross(e2, e1);
            const double sign = dot(outer, norm) >= 0.0 ? 1.0 : -1.0;
            const double inner = dot(e1, e2);
            double theta = std::atan2(sign * length(outer), inner);
            if (theta < 0.0) {
                theta = theta + 2.0 * Pi;
            }

            if (theta < minAngle) {
                minAngle = theta;
                minIndex = i;
            }
        }

        // Put new triangles
        std::vector<Face *> newFaces;
        if (minAngle <= Pi * (75.0 / 180.0)) {
            printf("minAngle = %f\n", minAngle / Pi * 180.0);
            const int prev = (minIndex - 1 + N) % N;
            const int next = (minIndex + 1) % N;

            Halfedge *new_he = new Halfedge;
            addHalfedge(new_he);

            new_he->src_ = boundary[prev];
            new_he->isBoundary_ = true;
            pair2he[E((prev - 1 + N) % N, prev)]->next_ = new_he;
            new_he->next_ = pair2he[E(next, (next + 1) % N)];

            pair2he[E(prev, next)] = new_he;
            he2pair[new_he] = E(prev, next);

            const T tri = std::make_tuple(prev, minIndex, next);
            Face *newFace = addNewTriangle(boundary, tri, pair2he, he2pair);
            newFaces.push_back(newFace);

            face->halfedge_ = new_he;

        } else if (minAngle <= Pi * (135.0 / 180.0)) {
            printf("minAngle = %f\n", minAngle / Pi * 180.0);
            const int prev = (minIndex - 1 + N) % N;
            const int next = (minIndex + 1) % N;

            const Vec3 p0 = boundary[prev]->pos();
            const Vec3 p1 = boundary[minIndex]->pos();
            const Vec3 p2 = boundary[next]->pos();
            const Vec3 e1 = p0 - p1;
            const Vec3 e2 = p2 - p1;
            const double len = (length(e1) + length(e2)) * 0.5;
            const Vec3 pnew = p1 + normalize(e1 + e2) * len;

            Vertex *new_vert = new Vertex(pnew);
            addVertex(new_vert);

            Halfedge *new_he1 = new Halfedge;
            Halfedge *new_he2 = new Halfedge;
            addHalfedge(new_he1);
            addHalfedge(new_he2);

            new_vert->halfedge_ = new_he2;

            new_he1->src_ = boundary[prev];
            new_he2->src_ = new_vert;
            new_he1->isBoundary_ = true;
            new_he2->isBoundary_ = true;

            pair2he[E((prev - 1 + N) % N, prev)]->next_ = new_he1;
            new_he1->next_ = new_he2;
            new_he2->next_ = pair2he[E(next, (next + 1) % N)];

            const int newIndex = (int)boundary.size();
            boundary.push_back(new_vert);

            pair2he[E(prev, newIndex)] = new_he1;
            pair2he[E(newIndex, next)] = new_he2;
            he2pair[new_he1] = E(prev, newIndex);
            he2pair[new_he2] = E(newIndex, next);

            const T t1 = std::make_tuple(prev, minIndex, newIndex);
            Face *new_face1 = addNewTriangle(boundary, t1, pair2he, he2pair);
            newFaces.push_back(new_face1);

            const T t2 = std::make_tuple(newIndex, minIndex, next);
            Face *new_face2 = addNewTriangle(boundary, t2, pair2he, he2pair);
            newFaces.push_back(new_face2);

            face->halfedge_ = new_he1;
        } else {
            printf("minAngle = %f\n", minAngle / Pi * 180.0);
            const int prev = (minIndex - 1 + N) % N;
            const int next = (minIndex + 1) % N;

            const Vec3 p0 = boundary[prev]->pos();
            const Vec3 p1 = boundary[minIndex]->pos();
            const Vec3 p2 = boundary[next]->pos();
            const Vec3 e1 = p0 - p1;
            const Vec3 e2 = p2 - p1;
            const double len = (length(e1) + length(e2)) * 0.5;

            const Vec3 pnew1 = p1 + normalize(e1 * 2.0 + e1) * len;
            const Vec3 pnew2 = p1 + normalize(e1 + e2 * 2.0) * len;

            Vertex *new_vert1 = new Vertex(pnew1);
            Vertex *new_vert2 = new Vertex(pnew2);
            addVertex(new_vert1);
            addVertex(new_vert2);
        }

        for (auto f : newFaces) {
            for (auto it = f->he_begin(); it != f->he_end(); ++it) {
                if (it->rev_ == nullptr) {
                    auto e = he2pair[it.ptr()];
                    int i = std::get<0>(e);
                    int j = std::get<1>(e);

                    auto r = E(j, i);
                    if (pair2he.count(r) != 0) {
                        Halfedge *rev = pair2he[r];
                        it->rev_ = rev;
                        rev->rev_ = it.ptr();
                        it->isBoundary_ = rev->isBoundary_;

                        Edge *new_edge = new Edge();
                        it->edge_ = new_edge;
                        rev->edge_ = new_edge;
                        new_edge->halfedge_ = it.ptr();
                        addEdge(new_edge);
                    } else {
                        printf("%d %d %p\n", i, j, it.ptr());
                    }
                }
            }
        }
    }
}

}  // namespace tinymesh

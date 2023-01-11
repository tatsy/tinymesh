#define TINYMESH_API_EXPORT
#include "mesh.h"

#include <iostream>
#include <fstream>
#include <utility>
#include <queue>
#include <unordered_set>

#include <Eigen/SparseCholesky>

#include "utils.h"
#include "eigen.h"
#include "vertex.h"
#include "face.h"
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
    // Compute average edge length
    const double avgEdgeLen = getMeanEdgeLength();

    // Original boundary vertices
    std::vector<Vertex *> boundary;
    for (auto it = face->v_begin(); it != face->v_end(); ++it) {
        boundary.push_back(it.ptr());
        it->lock();
    }

    // Initial Hole filling by advancing front
    while (true) {
        // Collect boundary vertices
        std::vector<Vertex *> face_vs;
        std::vector<Halfedge *> face_hes;
        for (auto it = face->he_begin(); it != face->he_end(); ++it) {
            face_vs.push_back(it->src());
            face_hes.push_back(it.ptr());
        }
        const int N = static_cast<int>(face_vs.size());
        if (N == 3) {
            for (auto it = face->he_begin(); it != face->he_end(); ++it) {
                it->isBoundary_ = false;
                it->rev_->isBoundary_ = false;
            }
            break;
        }

        // Merge adjacent vertices if they are too close
        bool isMerged = false;
        for (int i = 0; i < N; i++) {
            const int prev = (i - 1 + N) % N;
            const int next = (i + 1) % N;
            if (!face_vs[i]->isLocked() && !face_vs[i]->isLocked() &&
                length(face_vs[i]->pos() - face_vs[next]->pos()) < 0.5 * avgEdgeLen) {
                Halfedge *he = face_hes[i];
                std::vector<Halfedge *> v_hes;
                for (auto it = face_vs[i]->he_begin(); it != face_vs[i]->he_end(); ++it) {
                    v_hes.push_back(it.ptr());
                }

                Halfedge *he_prev = face_hes[prev];
                Halfedge *he_next = face_hes[next];
                Assertion(he_prev->next_ == he, "Error!!");
                Assertion(he->next_ == he_next, "Error!!");
                he_prev->next_ = he_next;

                Halfedge *he0 = he->rev_;
                Halfedge *he1 = he0->next_;
                Halfedge *he2 = he1->next_;
                Assertion(he2->next_ == he0, "Target face is not a triangle!");

                Vertex *v1 = he0->src_;  // == face_vs[next]
                Vertex *v0 = he1->src_;  // == face_vs[i]
                Vertex *v2 = he2->src_;  // opposite vertex of the edge to be removed.

                Halfedge *he1_rev = he1->rev_;
                Halfedge *he2_rev = he2->rev_;
                he1_rev->rev_ = he2_rev;
                he2_rev->rev_ = he1_rev;

                for (auto *v_he : v_hes) {
                    v_he->src_ = face_vs[next];
                }
                v1->halfedge_ = he2_rev;
                v2->halfedge_ = he1_rev;
                face->halfedge_ = he_next;

                v1->setPos(0.5 * (v0->pos() + v1->pos()));

                removeFace(he0->face_);
                removeHalfedge(he);
                removeHalfedge(he0);
                removeHalfedge(he1);
                removeHalfedge(he2);
                removeVertex(v0);

                isMerged = true;
                break;
            }
        }

        if (isMerged) {
            verify();
            continue;
        }

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

            const Vec3 e1 = face_vs[prev]->pos() - face_vs[i]->pos();
            const Vec3 e2 = face_vs[next]->pos() - face_vs[i]->pos();
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
            const int prev = (minIndex - 1 + N) % N;
            const int next = (minIndex + 1) % N;

            Halfedge *new_he = new Halfedge;
            addHalfedge(new_he);

            new_he->src_ = face_vs[prev];
            new_he->isBoundary_ = true;
            pair2he[E((prev - 1 + N) % N, prev)]->next_ = new_he;
            new_he->next_ = pair2he[E(next, (next + 1) % N)];

            pair2he[E(prev, next)] = new_he;
            he2pair[new_he] = E(prev, next);

            const T tri = std::make_tuple(prev, minIndex, next);
            Face *newFace = addNewTriangle(face_vs, tri, pair2he, he2pair);
            newFaces.push_back(newFace);

            face->halfedge_ = new_he;
            new_he->face_ = face;

        } else if (minAngle <= Pi * (135.0 / 180.0)) {
            const int prev = (minIndex - 1 + N) % N;
            const int next = (minIndex + 1) % N;

            const Vec3 p0 = face_vs[prev]->pos();
            const Vec3 p1 = face_vs[minIndex]->pos();
            const Vec3 p2 = face_vs[next]->pos();
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

            new_he1->src_ = face_vs[prev];
            new_he2->src_ = new_vert;
            new_he1->isBoundary_ = true;
            new_he2->isBoundary_ = true;

            pair2he[E((prev - 1 + N) % N, prev)]->next_ = new_he1;
            new_he1->next_ = new_he2;
            new_he2->next_ = pair2he[E(next, (next + 1) % N)];

            const int newIndex = (int)face_vs.size();
            face_vs.push_back(new_vert);

            pair2he[E(prev, newIndex)] = new_he1;
            pair2he[E(newIndex, next)] = new_he2;
            he2pair[new_he1] = E(prev, newIndex);
            he2pair[new_he2] = E(newIndex, next);

            const T t1 = std::make_tuple(prev, minIndex, newIndex);
            Face *new_face1 = addNewTriangle(face_vs, t1, pair2he, he2pair);
            newFaces.push_back(new_face1);

            const T t2 = std::make_tuple(newIndex, minIndex, next);
            Face *new_face2 = addNewTriangle(face_vs, t2, pair2he, he2pair);
            newFaces.push_back(new_face2);

            face->halfedge_ = new_he1;
            new_he1->face_ = face;
            new_he2->face_ = face;

        } else {
            const int prev = (minIndex - 1 + N) % N;
            const int next = (minIndex + 1) % N;

            const Vec3 p0 = face_vs[prev]->pos();
            const Vec3 p1 = face_vs[minIndex]->pos();
            const Vec3 p2 = face_vs[next]->pos();
            const Vec3 e1 = p0 - p1;
            const Vec3 e2 = p2 - p1;
            const double len1 = (length(e1) * 2.0 + length(e2)) / 3.0;
            const double len2 = (length(e1) + length(e2) * 2.0) / 3.0;

            const Vec3 pnew1 = p1 + normalize(e1 * 2.0 + e2) * len1;
            const Vec3 pnew2 = p1 + normalize(e1 + e2 * 2.0) * len2;

            Vertex *new_vert1 = new Vertex(pnew1);
            Vertex *new_vert2 = new Vertex(pnew2);
            addVertex(new_vert1);
            addVertex(new_vert2);

            Halfedge *new_he1 = new Halfedge;
            Halfedge *new_he2 = new Halfedge;
            Halfedge *new_he3 = new Halfedge;
            addHalfedge(new_he1);
            addHalfedge(new_he2);
            addHalfedge(new_he3);

            new_vert1->halfedge_ = new_he2;
            new_vert2->halfedge_ = new_he3;

            new_he1->src_ = face_vs[prev];
            new_he2->src_ = new_vert1;
            new_he3->src_ = new_vert2;
            new_he1->isBoundary_ = true;
            new_he2->isBoundary_ = true;
            new_he3->isBoundary_ = true;

            pair2he[E((prev - 1 + N) % N, prev)]->next_ = new_he1;
            new_he1->next_ = new_he2;
            new_he2->next_ = new_he3;
            new_he3->next_ = pair2he[E(next, (next + 1) % N)];

            const int newIndex1 = (int)face_vs.size();
            face_vs.push_back(new_vert1);
            const int newIndex2 = (int)face_vs.size();
            face_vs.push_back(new_vert2);

            pair2he[E(prev, newIndex1)] = new_he1;
            pair2he[E(newIndex1, newIndex2)] = new_he2;
            pair2he[E(newIndex2, next)] = new_he3;
            he2pair[new_he1] = E(prev, newIndex1);
            he2pair[new_he2] = E(newIndex1, newIndex2);
            he2pair[new_he3] = E(newIndex2, next);

            const T t1 = std::make_tuple(prev, minIndex, newIndex1);
            Face *new_face1 = addNewTriangle(face_vs, t1, pair2he, he2pair);
            newFaces.push_back(new_face1);

            const T t2 = std::make_tuple(newIndex1, minIndex, newIndex2);
            Face *new_face2 = addNewTriangle(face_vs, t2, pair2he, he2pair);
            newFaces.push_back(new_face2);

            const T t3 = std::make_tuple(newIndex2, minIndex, next);
            Face *new_face3 = addNewTriangle(face_vs, t3, pair2he, he2pair);
            newFaces.push_back(new_face3);

            face->halfedge_ = new_he1;
            new_he1->face_ = face;
            new_he2->face_ = face;
            new_he3->face_ = face;
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
                    } else {
                        printf("%d %d %p\n", i, j, it.ptr());
                    }
                }
            }
        }
    }

    // Correct internal vertices
    std::unordered_set<Vertex *> insideVerts;
    std::queue<Vertex *> que;
    for (auto it = face->v_begin(); it != face->v_end(); ++it) {
        if (!it->isLocked()) {
            que.push(it.ptr());
            break;
        }
    }
    while (!que.empty()) {
        Vertex *v = que.front();
        que.pop();

        if (insideVerts.count(v) != 0) {
            continue;
        }
        insideVerts.insert(v);

        for (auto it = v->v_begin(); it != v->v_end(); ++it) {
            if (!it->isLocked() && insideVerts.count(it.ptr()) == 0) {
                que.push(it.ptr());
            }
        }
    }

    // Smooth vertex normals
    std::unordered_map<Vertex *, Vec3> newNormals;
    for (Vertex *v : insideVerts) {
        newNormals[v] = v->normal();
    }

    for (int kIter = 0; kIter < 1000; kIter++) {
        for (Vertex *v : insideVerts) {
            Vec3 nn = Vec3(0.0);
            double sumWgt = 0.0;
            for (auto it = v->he_begin(); it != v->he_end(); ++it) {
                const Vec3 &n = newNormals[it->dst()];
                const double weight = it->cotWeight();
                nn += weight * n;
                sumWgt += weight;
            }
            nn = 0.5 * newNormals[v] + 0.5 * nn / sumWgt;
            newNormals[v] = normalize(nn);
        }
    }

    // Compute face normals
    std::unordered_set<Face *> insideFaces;
    for (Vertex *v : insideVerts) {
        for (auto it = v->f_begin(); it != v->f_end(); ++it) {
            Face *f = it.ptr();
            bool inside = true;
            for (auto vit = f->v_begin(); vit != f->v_end(); ++vit) {
                if (insideVerts.count(vit.ptr()) == 0 && !it->isLocked()) {
                    inside = false;
                    break;
                }
            }

            if (inside) {
                insideFaces.insert(f);
            }
        }
    }

    std::unordered_map<Face *, EigenMatrix3> rotations;
    for (Face *f : insideFaces) {
        const Vec3 n0 = f->normal();
        Vec3 n1(0.0);
        for (auto it = f->v_begin(); it != f->v_end(); ++it) {
            n1 += newNormals[it.ptr()];
        }
        n1 = normalize(n1);

        const Vec3 axis = cross(n0, n1);
        const double theta = std::atan2(length(axis), dot(n0, n1));
        rotations[f] = rotationAxisAngle(theta, axis);
    }

    // Compute rotated gradients
    std::unordered_map<Vertex *, Vec3> rotGrads;
    for (Vertex *v : insideVerts) {
        Vec3 sumGrad(0.0);
        int count = 0;
        for (auto it = v->he_begin(); it != v->he_end(); ++it) {
            Assertion(it->next()->next()->next() == it.ptr(), "Non-triangle face is detected!");
            Vertex *v0 = it->src();
            Vertex *v1 = it->next()->src();
            Vertex *v2 = it->next()->next()->src();
            const Vec3 p0 = v0->pos();
            const Vec3 p1 = v1->pos();
            const Vec3 p2 = v2->pos();

            const Vec3 l12 = p2 - p1;
            const Vec3 e10 = normalize(p0 - p1);
            const Vec3 e12 = normalize(p2 - p1);
            const Vec3 perp = cross(e12, e10);
            const Vec3 e12_rot90 = normalize(cross(e12, perp)) * length(l12);

            // NOTE: gradient should be divided by "2A" rather than "2,"
            // but here, it's just divided by the latter to make its value be
            // proportional to the edge length.
            // const double A = 0.5 * length(cross(p1 - p0, p2 - p0));
            const Vec3 grad = (Vec3)(rotations[it->face()] * (EigenVector3)e12_rot90) / 2.0;
            sumGrad += grad;
            count++;
        }
        rotGrads[v] = sumGrad / (double)count;
    }

    // Solve Poisson equation
    std::vector<EigenTriplet> tripInner;
    std::vector<EigenTriplet> tripOuter;
    std::unordered_map<Vertex *, uint32_t> uniqueInner;
    std::unordered_map<Vertex *, uint32_t> uniqueOuter;
    int countInner = 0;
    int countOuter = 0;
    for (Vertex *v : insideVerts) {
        if (uniqueInner.count(v) == 0) {
            uniqueInner[v] = countInner;
            countInner++;
        }
        const int i = uniqueInner[v];
        double sumWgt = 0.0;
        for (auto it = v->he_begin(); it != v->he_end(); ++it) {
            Vertex *u = it->dst();
            const double weight = it->cotWeight();
            if (insideVerts.count(u) != 0) {
                // inner vertex
                if (uniqueInner.count(u) == 0) {
                    uniqueInner[u] = countInner;
                    countInner++;
                }
                const int j = uniqueInner[u];
                tripInner.emplace_back(i, j, -weight);
            } else {
                // outer vertex
                if (uniqueOuter.count(u) == 0) {
                    uniqueOuter[u] = countOuter;
                    countOuter++;
                }
                const int j = uniqueOuter[u];
                tripOuter.emplace_back(i, j, -weight);
            }
            sumWgt += weight;
        }
        tripInner.emplace_back(i, i, sumWgt);
    }

    EigenMatrix bbInner(countInner, 3);
    bbInner.setZero();
    for (auto it : uniqueInner) {
        const Vec3 &rg = rotGrads[it.first];
        bbInner.row(it.second) << rg.x(), rg.y(), rg.z();
    }

    EigenMatrix xxOuter(countOuter, 3);
    xxOuter.setZero();
    for (auto it : uniqueOuter) {
        const Vec3 &p = it.first->pos();
        xxOuter.row(it.second) << p.x(), p.y(), p.z();
    }

    EigenSparseMatrix LL(countInner, countInner);
    LL.setFromTriplets(tripInner.begin(), tripInner.end());

    EigenSparseMatrix SS(countInner, countOuter);
    SS.setFromTriplets(tripOuter.begin(), tripOuter.end());

    EigenMatrix bb = bbInner - SS * xxOuter;
    Eigen::SimplicialLDLT<EigenSparseMatrix> solver(LL);
    EigenMatrix xx = solver.solve(bb);

    for (auto it : uniqueInner) {
        Vertex *v = it.first;
        const uint32_t i = it.second;
        const Vec3 newPos(xx(i, 0), xx(i, 1), xx(i, 2));
        v->setPos(newPos);
    }
}

}  // namespace tinymesh

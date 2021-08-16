#define TINYMESH_API_EXPORT
#include "abfxx.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <tuple>
#include <set>
#include <queue>
#include <map>
#include <unordered_map>

#include "core/eigen.h"
#include "polymesh/vertex.h"
#include "polymesh/face.h"
#include "polymesh/mesh.h"
#include "polymesh/halfedge.h"

#include <Eigen/IterativeLinearSolvers>
using LinearSolver = Eigen::LeastSquaresConjugateGradient<EigenSparse>;

namespace {

double calcAngle(const Vec3 &p0, const Vec3 &p1, const Vec3 &p2) {
    const Vec3 e1 = p0 - p1;
    const Vec3 e2 = p2 - p1;
    return std::atan2(length(cross(e1, e2)), dot(e1, e2));
}

}  // anonymous namespace

struct AngleData {
    AngleData() {}
    AngleData(int t, int k) : t(t), k(k) {}
    AngleData(int t, int k, int i, double a) : t(t), k(k), i(i), a(a) {}

    bool operator==(const AngleData &other) const {
        return t == other.t && k == other.k;
    }

    int t;     // Triangle ID
    int k;     // Vertex ID
    int i;     // AngleID
    double a;  // Initial angle
};

namespace std {

template <>
struct hash<AngleData> {
    size_t operator()(const AngleData &a) const {
        size_t ret = 0;
        ret = (ret < 1) & std::hash<int>()(a.t);
        ret = (ret < 1) & std::hash<int>()(a.k);
        return ret;
    }
};

}  // namespace std

namespace tinymesh {

void abfxx(Mesh &mesh, int maxiter, double epsilon) {
    // Check input mesh indices
    for (int vi = 0; vi < (int)mesh.num_vertices(); vi++) {
        if (vi != mesh.vertex(vi)->index()) {
            throw std::runtime_error("Vertex order and its index do not match!");
        }
    }

    for (int fi = 0; fi < (int)mesh.num_faces(); fi++) {
        if (fi != mesh.face(fi)->index()) {
            throw std::runtime_error("Face order and its index do not match!");
        }
    }

    // Compute average edge length
    const int nFaces = mesh.num_faces();
    const int nVertices = mesh.num_vertices();
    double avgLength = 0.0;
    for (int i = 0; i < (int)mesh.num_halfedges(); i++) {
        avgLength += mesh.halfedge(i)->length();
    }
    avgLength /= mesh.num_halfedges();
    printf("Avg: %f\n", avgLength);

    // Traverse all the angles
    std::vector<AngleData> angles;
    int count = 0;
    for (int faceID = 0; faceID < (int)mesh.num_faces(); faceID++) {
        Face *t = mesh.face(faceID);
        std::vector<Vertex*> vs;
        for (auto vit = t->v_begin(); vit != t->v_end(); ++vit) {
            vs.push_back(vit.ptr());
        }

        const int nv = (int)vs.size();
        for (int i = 0; i < (int)vs.size(); i++) {
            const int pre = (i - 1 + nv) % nv;
            const int post = (i + 1) % nv;
            const double a = calcAngle(vs[pre]->pos(), vs[i]->pos(), vs[post]->pos());
            angles.emplace_back(t->index(), vs[i]->index(), count++, std::abs(a));
        }
    }

    // Store in unordered map
    const int nAngles = (int)angles.size();
    std::unordered_map<AngleData, int> memo;
    for (int i = 0; i < nAngles; i++) {
        memo[angles[i]] = i;
    }

    // Traverse forward vertices group
    std::vector<std::vector<int>> forwardGroup(nAngles, std::vector<int>());
    std::vector<std::vector<int>> backwardGroup(nAngles, std::vector<int>());
    for (int i = 0; i < nAngles; i++) {
        Vertex *v = mesh.vertex(angles[i].k);
        for (auto it = v->ohe_begin(); it != v->ohe_end(); ++it) {
            Vertex *k = it->dst();
            Face *t = it->face();
            forwardGroup[i].push_back(memo[AngleData(t->index(), k->index())]);
        }

        for (auto it = v->ihe_begin(); it != v->ihe_end(); ++it) {
            Vertex *k = it->src();
            Face *t = it->face();
            backwardGroup[i].push_back(memo[AngleData(t->index(), k->index())]);
        }
    }

    // Initialize variables
    const int nConstraints = nFaces + nVertices * 2;

    EigenVector alpha(nAngles);
    EigenVector lambda(nConstraints);
    for (int i = 0; i < nAngles; i++) {
        alpha(i) = angles[i].a;
    }
    lambda.setOnes();

    // Optimization loop
    int constID = 0;
    EigenVector invL(nAngles);
    EigenSparse J(nConstraints, nAngles);
    EigenVector b1(nAngles);
    EigenVector b2(nConstraints);
    EigenVector deltaAlpha = EigenVector::Zero(nAngles);
    EigenVector deltaLambda = EigenVector::Zero(nConstraints);
    std::vector<EigenTriplet> triplets;

    const double stepSizeDecay = 0.5;
    double stepSize = 0.995;
    double prevNorm = 1.0e20;
    for (int loop = 0; loop < maxiter; loop++) {
        // Initialization for next loop
        constID = 0;
        triplets.clear();
        b1.setZero();
        b2.setZero();

        // Update for E
        {
            for (int i = 0; i < nAngles; i++) {
                const double b = angles[i].a;
                const double a = alpha(i);
                const double w = 1.0 / (b * b);
                b1(i) -= 2.0 * (a - b) / w;
                invL(i) = w / 2.0;
            }
        }

        // Update for C_tri
        {
            for (int fi = 0; fi < nFaces; fi++) {
                Face *t = mesh.face(fi);
                double sumAngle = 0.0;
                for (auto it = t->v_begin(); it != t->v_end(); ++it) {
                    Vertex *k = it.ptr();
                    const int angleID = memo[AngleData(t->index(), k->index())];
                    sumAngle += alpha(angleID);
                }
                sumAngle -= Pi;

                for (auto it = t->v_begin(); it != t->v_end(); ++it) {
                    Vertex *k = it.ptr();
                    const int angleID = memo[AngleData(t->index(), k->index())];

                    b1(angleID) -= lambda(constID);
                    triplets.emplace_back(constID, angleID, 1.0);
                }

                b2(constID) -= sumAngle;
                constID += 1;
            }
        }

        // Update for C_plan
        {
            for (int vi = 0; vi < nVertices; vi++) {
                Vertex *k = mesh.vertex(vi);
                double sumAngle = 0.0;
                for (auto it = k->f_begin(); it != k->f_end(); ++it) {
                    Face *t = it.ptr();
                    const int angleID = memo[AngleData(t->index(), k->index())];
                    sumAngle += alpha(angleID);
                }
                sumAngle -= 2.0 * Pi;

                for (auto it = k->f_begin(); it != k->f_end(); ++it) {
                    Face *t = it.ptr();
                    const int angleID = memo[AngleData(t->index(), k->index())];

                    b1(angleID) -= lambda(constID);
                    triplets.emplace_back(constID, angleID, 1.0);
                }

                b2(constID) -= sumAngle;
                constID += 1;
            }
        }

        // Update for C_len
        {
            for (int vi = 0; vi < nVertices; vi++) {
                Vertex *k = mesh.vertex(vi);
                double prodForward = 1.0;
                double prodBackward = 1.0;
                for (int i = 0; i < (int)forwardGroup[vi].size(); i++) {
                    const double a = alpha(forwardGroup[vi][i]);
                    prodForward *= std::sin(a);
                }

                for (int i = 0; i < (int)backwardGroup[vi].size(); i++) {
                    const double a = alpha(backwardGroup[vi][i]);
                    prodBackward *= std::sin(a);
                }

                for (int i = 0; i < (int)forwardGroup[vi].size(); i++) {
                    const int angleID = forwardGroup[vi][i];
                    const double a = alpha(angleID);
                    const double grad = prodForward / std::tan(a);

                    b1(angleID) -= lambda(constID) * grad;
                    triplets.emplace_back(constID, angleID, grad);
                }

                for (int i = 0; i < (int)backwardGroup[vi].size(); i++) {
                    const int angleID = backwardGroup[vi][i];
                    const double a = alpha(angleID);
                    const double grad = -prodBackward / std::tan(a);

                    b1(angleID) -= lambda(constID) * grad;
                    triplets.emplace_back(constID, angleID, grad);
                }

                const double C_len = prodForward - prodBackward;
                b2(constID) -= C_len;
                constID += 1;
            }
        }

        // Termination criteria
        const double gradNorm = std::sqrt(b1.squaredNorm() + b2.squaredNorm());
        printf("[#%d] gnorm = %.6f, step = %.6f\n", loop + 1, gradNorm, stepSize);
        if (gradNorm < epsilon || stepSize < epsilon) {
            break;
        }

        if (prevNorm <= gradNorm) {
            if (loop == 0) {
                throw std::runtime_error("Optimization does not work!");
            }

            alpha -= (stepSize * (1.0 - stepSizeDecay)) * deltaAlpha;
            lambda -= (stepSize * (1.0 - stepSizeDecay)) * deltaLambda;
            stepSize *= stepSizeDecay;
            loop -= 1;
            continue;
        }
        prevNorm = gradNorm;

        // Set J
        J.setFromTriplets(triplets.begin(), triplets.end());

        // Update
        EigenVector b_star = J * (invL.asDiagonal() * b1) - b2;
        EigenSparse AA = J * (invL.asDiagonal() * J.transpose());
        LinearSolver solver;
        solver.compute(AA);
        deltaLambda = solver.solve(b_star);
        if (solver.info() != Eigen::Success) {
            fprintf(stderr, "Failed to factorize matrix!\n");
            break;
        }

        deltaAlpha = invL.asDiagonal() * (b1 - J.transpose() * deltaLambda);

        alpha += stepSize * deltaAlpha;
        lambda += stepSize * deltaLambda;

        // Make angles positive
        for (int i = 0; i < nAngles; i++) {
            while (alpha(i) < 0.0) alpha(i) += 2.0 * Pi;
            while (alpha(i) > 2.0 * Pi) alpha(i) -= 2.0 * Pi;
            if (alpha(i) > Pi) {
                alpha(i) = 2.0 * Pi - alpha(i);
            }
        }
    }

    // Flattening
    // See the paper:
    // A. Sheffer and E. de Sturler, 2001
    // "Parameterization  of  Faceted  Surfaces  for  Meshing  using  Angle-BasedFlattening"
    std::unordered_map<Vec3, uint32_t> uniqueVertices;
    std::vector<Vec3> vertices;
    std::vector<uint32_t> indices;
    {
        using Pair = std::pair<int, int>;
        std::queue<Halfedge *> S;
        std::vector<int> visited(nFaces, 0);
        std::unordered_map<int, Vec3> project;

        // Set initial half-edge
        Halfedge *he = mesh.halfedge(0);
        {
            const Vec3 va = Vec3(0.0, 0.0, 0.0);
            const Vec3 vb = Vec3(he->length(), 0.0, 0.0);
            project[he->src()->index()] = va;
            project[he->dst()->index()] = vb;
            S.push(he);
            S.push(he->rev());

            uniqueVertices[va] = 0;
            uniqueVertices[vb] = 1;
            vertices.push_back(va);
            vertices.push_back(vb);
        }

        // Traverse
        while (!S.empty()) {
            Halfedge *he = S.front();
            S.pop();

            const int faceID = he->face()->index();
            if (visited[faceID]) {
                continue;
            }
            visited[faceID] = 1;

            const int a = he->src()->index();
            const int b = he->dst()->index();
            const int c = he->next()->dst()->index();
            const Vec3 va = project[a];
            const Vec3 vb = project[b];
            Vec3 vc;
            if (project.count(c) != 0) {
                vc = project[c];
            } else {
                const Vec3 e = normalize(va - vb);
                const int angleID = memo[AngleData(faceID, b)];
                const double angle = alpha(angleID);
                const double rotx = e.x() * std::cos(angle) - e.y() * std::sin(angle);
                const double roty = e.x() * std::sin(angle) + e.y() * std::cos(angle);
                vc = vb + he->next()->length() * Vec3(rotx, roty, 0.0);
                project[c] = vc;
            }

            if (length(va - vb) > avgLength * 2.0 || length(vb - vc) > avgLength * 2.0 || length(vc - va) > avgLength * 2.0) {
                continue;
            }

            if (uniqueVertices.count(vc) == 0) {
                uniqueVertices[vc] = (uint32_t)vertices.size();
                vertices.push_back(vc);
            }
            indices.push_back(uniqueVertices[va]);
            indices.push_back(uniqueVertices[vb]);
            indices.push_back(uniqueVertices[vc]);

            S.push(he->next()->rev());
            S.push(he->next()->next()->rev());
       }

        // Save
        std::ofstream writer("project.obj", std::ios::out);
        {
            for (const auto &v : vertices) {
                writer << "v " << v.x() << " " << v.y() << " " << v.z() << std::endl;
            }

            for (int i = 0; i < indices.size(); i += 3) {
                writer << "f " << (indices[i + 0] + 1)
                       << " "  << (indices[i + 1] + 1)
                       << " "  << (indices[i + 2] + 1) << std::endl;
            }
        }
        writer.close();
    }
}

}  // namespace tinymesh

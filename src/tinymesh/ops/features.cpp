#define TINYMESH_API_EXPORT
#include "ops.h"

#include <Eigen/Eigenvalues>
#include <Eigen/Cholesky>
#include <Spectra/SymEigsSolver.h>
#include <Spectra/MatOp/DenseSymMatProd.h>
#include <Spectra/MatOp/SparseSymMatProd.h>

#include "core/vertex.h"
#include "core/face.h"

namespace tinymesh {

EigenMatrix getHeatKernelSignatures(const EigenSparseMatrix &L, int K, int nTimes) {
    Assertion(K < L.rows(), "Eigenvalues less than matrix size can only be requested!");
    const int ncv = std::min((int)L.rows(), K * 2);
    Spectra::SparseSymMatProd<FloatType> op(L);
    Spectra::SymEigsSolver<Spectra::SparseSymMatProd<FloatType>> eigs(op, K, ncv);
    eigs.init();

    eigs.compute(Spectra::SortRule::SmallestMagn, 1000, 1.0e-10, Spectra::SortRule::SmallestMagn);
    if (eigs.info() != Spectra::CompInfo::Successful) {
        Error("Eigen decomposition failed!");
    }

    EigenMatrix U = eigs.eigenvectors();
    EigenVector lambda = eigs.eigenvalues();

    const double t_min = 4.0 * std::log(10) / std::max(1.0e-12, lambda(K - 1));
    const double t_max = 4.0 * std::log(10) / std::max(1.0e-12, lambda(1));

    const double log_t_min = std::log(t_min);
    const double log_t_max = std::log(t_max);
    const double inc = (log_t_max - log_t_min) / (nTimes - 1);
    EigenVector times(nTimes);
    for (int i = 0; i < nTimes; i++) {
        times(i) = std::exp(log_t_min + i * inc);
    }

    const int N = L.rows();
    EigenMatrix HKS(N, nTimes);
    HKS.setZero();
    for (int i = 0; i < K; i++) {
        const EigenVector exp = (-times * lambda(i)).array().exp().matrix();
        const EigenVector U2 = U.col(i).array().square().matrix();
        HKS += U2 * exp.transpose();
    }

    return HKS;
}

void getCurvatureTensors(const Mesh &mesh, std::vector<EigenMatrix2> &tensors, std::vector<LocalFrame> &frames) {
    // Compute per-triangle tensors
    std::vector<EigenMatrix2> faceTensors;
    std::vector<LocalFrame> faceFrames;
    for (int i = 0; i < mesh.numFaces(); i++) {
        Face *f = mesh.face(i);
        std::vector<Vertex *> vertices;
        for (auto it = f->v_begin(); it != f->v_end(); ++it) {
            vertices.push_back(it.ptr());
        }
        Assertion(vertices.size() == 3, "Non-triangle face detected!");

        Vertex *v0 = vertices[0];
        Vertex *v1 = vertices[1];
        Vertex *v2 = vertices[2];
        const Vec3 p0 = v0->pos();
        const Vec3 p1 = v1->pos();
        const Vec3 p2 = v2->pos();
        const Vec3 e0 = p2 - p1;
        const Vec3 e1 = p0 - p2;
        const Vec3 e2 = p1 - p0;
        const Vec3 n0 = v0->normal();
        const Vec3 n1 = v1->normal();
        const Vec3 n2 = v2->normal();

        const Vec3 u = normalize(p1 - p0);
        const Vec3 w = normalize(cross(u, p2 - p0));
        const Vec3 v = normalize(cross(w, u));

        EigenMatrix A(6, 3);
        EigenVector b(6);
        A.setZero();
        // e0
        A(0, 0) = dot(e0, u);
        A(0, 1) = dot(e0, v);
        b(0) = dot(n2 - n1, u);
        A(1, 1) = dot(e0, u);
        A(1, 2) = dot(e0, v);
        b(1) = dot(n2 - n1, v);
        // e1
        A(2, 0) = dot(e1, u);
        A(2, 1) = dot(e1, v);
        b(2) = dot(n0 - n2, u);
        A(3, 1) = dot(e1, u);
        A(3, 2) = dot(e1, v);
        b(3) = dot(n0 - n2, v);
        // e2
        A(4, 0) = dot(e2, u);
        A(4, 1) = dot(e2, v);
        b(4) = dot(n1 - n0, u);
        A(5, 1) = dot(e2, u);
        A(5, 2) = dot(e2, v);
        b(5) = dot(n1 - n0, v);

        EigenMatrix AA = A.transpose() * A;
        EigenVector Ab = A.transpose() * b;

        Eigen::LLT<EigenMatrix> solver(AA);
        EigenVector3 efg = solver.solve(Ab);

        EigenMatrix2 M;
        M << efg(0), efg(1),  //
            efg(1), efg(2);   //
        faceTensors.push_back(M);
        faceFrames.emplace_back(u, v, w);
    }

    tensors.clear();
    frames.clear();
    tensors.reserve(mesh.numVertices());
    frames.reserve(mesh.numVertices());
    for (int i = 0; i < mesh.numVertices(); i++) {
        Vertex *vt = mesh.vertex(i);
        const Vec3 nv = vt->normal();
        Vec3 uv, vv;
        if (std::abs(nv.x()) > 0.1) {
            uv = normalize(cross(nv, Vec3(0.0, 1.0, 0.0)));
        } else {
            uv = normalize(cross(nv, Vec3(1.0, 0.0, 0.0)));
        }
        vv = normalize(cross(nv, uv));

        EigenMatrix2 M;
        M.setZero();
        double sumWgt = 0.0;
        for (auto it = vt->f_begin(); it != vt->f_end(); ++it) {
            const LocalFrame F = faceFrames[it->index()];
            const Vec3 uf = std::get<0>(F);
            const Vec3 vf = std::get<1>(F);
            const Vec3 nf = std::get<2>(F);

            const double theta = std::atan2(length(cross(nv, nf)), dot(nv, nf));
            EigenMatrix3 R;
            if (theta < 1.0e-3) {
                R.setIdentity();
            } else {
                const Vec3 axis = normalize(cross(nv, nf));
                R = rotationAxisAngle(theta, axis);
            }

            const Vec3 rot_u = (Vec3)(R * (EigenVector3)uv);
            const Vec3 rot_v = (Vec3)(R * (EigenVector3)vv);
            const EigenMatrix2 Mf = faceTensors[it->index()];

            const double dot_u0 = dot(rot_u, uf);
            const double dot_u1 = dot(rot_u, vf);
            const double dot_v0 = dot(rot_v, uf);
            const double dot_v1 = dot(rot_v, vf);
            const double ev =
                Mf(0, 0) * dot_u0 * dot_u0 + 2.0 * Mf(0, 1) * dot_u0 * dot_u1 + Mf(1, 1) * dot_u1 * dot_u1;
            const double gv =
                Mf(0, 0) * dot_v0 * dot_v0 + 2.0 * Mf(0, 1) * dot_v0 * dot_v1 + Mf(1, 1) * dot_v1 * dot_v1;
            const double fv = Mf(0, 0) * dot_u0 * dot_v0 + Mf(0, 1) * dot_u0 * dot_v1 + Mf(1, 0) * dot_u1 * dot_v0 +
                              Mf(1, 1) * dot_u1 * dot_v1;

            EigenMatrix2 Mv;
            Mv << ev, fv,  //
                fv, gv;    //
            const double weight = it->area() / 3.0;
            M += weight * Mv;
            sumWgt += weight;
        }
        tensors.push_back(M / sumWgt);
        frames.emplace_back(uv, vv, nv);
    }
}

std::tuple<std::vector<double>, std::vector<double>, std::vector<Vec3>, std::vector<Vec3>> getPrincipalCurvatures(
    const Mesh &mesh) {
    std::vector<EigenMatrix2> tensors;
    std::vector<LocalFrame> frames;
    getCurvatureTensors(mesh, tensors, frames);

    std::vector<double> ks_max;
    std::vector<double> ks_min;
    std::vector<Vec3> ts_max;
    std::vector<Vec3> ts_min;
    for (size_t i = 0; i < mesh.numVertices(); i++) {
        const EigenMatrix2 M = tensors[i];
        const LocalFrame F = frames[i];
        const Vec3 u = std::get<0>(F);
        const Vec3 v = std::get<1>(F);

        Eigen::SelfAdjointEigenSolver<EigenMatrix2> solver(M);
        const EigenMatrix2 eigvec = solver.eigenvectors();
        const EigenVector2 eigval = solver.eigenvalues();

        double k_max = eigval(0);
        double k_min = eigval(1);
        Vec3 t_max = eigvec(0, 0) * u + eigvec(0, 1) * v;
        Vec3 t_min = eigvec(1, 0) * u + eigvec(1, 1) * v;
        if (k_max < k_min) {
            std::swap(k_max, k_min);
            std::swap(t_max, t_min);
        }

        ks_max.push_back(k_max);
        ks_min.push_back(k_min);
        ts_max.push_back(t_max);
        ts_min.push_back(t_min);
    }

    return { ks_max, ks_min, ts_max, ts_min };
}

}  // namespace tinymesh

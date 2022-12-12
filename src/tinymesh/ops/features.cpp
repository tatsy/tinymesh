#define TINYMESH_API_EXPORT
#include "ops.h"

#include <Eigen/Eigenvalues>
#include <Eigen/Cholesky>
#include <Spectra/SymEigsSolver.h>
#include <Spectra/MatOp/DenseSymMatProd.h>
#include <Spectra/MatOp/SparseSymMatProd.h>

#include "core/vertex.h"
#include "core/face.h"

namespace {

void getPrincipalCurvaturesFromTensors(std::vector<double> &ks_max, std::vector<double> &ks_min,
                                       std::vector<Vec3> &ts_max, std::vector<Vec3> &ts_min,
                                       const std::vector<EigenMatrix2> &tensors, std::vector<LocalFrame> &frames) {
    const int Nv = (int)tensors.size();
    ks_max.resize(Nv);
    ks_min.resize(Nv);
    ts_max.resize(Nv);
    ts_min.resize(Nv);

    for (size_t i = 0; i < Nv; i++) {
        const EigenMatrix2 M = tensors[i];
        const LocalFrame F = frames[i];
        const Vec3 u = std::get<0>(F);
        const Vec3 v = std::get<1>(F);

        Eigen::SelfAdjointEigenSolver<EigenMatrix2> solver(M);
        const EigenMatrix2 eigvec = solver.eigenvectors();
        const EigenVector2 eigval = solver.eigenvalues();

        double k_max = eigval(0);
        double k_min = eigval(1);
        Vec3 t_max = normalize(eigvec(0, 0) * u + eigvec(0, 1) * v);
        Vec3 t_min = normalize(eigvec(1, 0) * u + eigvec(1, 1) * v);
        if (k_max < k_min) {
            std::swap(k_max, k_min);
            std::swap(t_max, t_min);
        }

        ks_max[i] = k_max;
        ks_min[i] = k_min;
        ts_max[i] = t_max;
        ts_min[i] = t_min;
    }
}

}  // namespace

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

            EigenMatrix3 R = rotationFromTwoVectors(nv, nf);
            const Vec3 rot_u = normalize(R * uv);
            const Vec3 rot_v = normalize(R * vv);
            const EigenMatrix2 Mf = faceTensors[it->index()];

            const double dot_u0 = dot(rot_u, uf);
            const double dot_u1 = dot(rot_u, vf);
            const double dot_v0 = dot(rot_v, uf);
            const double dot_v1 = dot(rot_v, vf);

            const double ev = Mf(0, 0) * dot_u0 * dot_u0 +        //
                              2.0 * Mf(0, 1) * dot_u0 * dot_u1 +  //
                              Mf(1, 1) * dot_u1 * dot_u1;         //
            const double gv = Mf(0, 0) * dot_v0 * dot_v0 +        //
                              2.0 * Mf(0, 1) * dot_v0 * dot_v1 +  //
                              Mf(1, 1) * dot_v1 * dot_v1;         //
            const double fv = Mf(0, 0) * dot_u0 * dot_v0 +        //
                              Mf(0, 1) * dot_u0 * dot_v1 +        //
                              Mf(1, 0) * dot_u1 * dot_v0 +        //
                              Mf(1, 1) * dot_u1 * dot_v1;         //

            EigenMatrix2 Mv;
            Mv << ev, fv,  //
                fv, gv;    //

            const double weight = it->area();
            M += weight * Mv;
            sumWgt += weight;
        }
        tensors.push_back(M / sumWgt);
        frames.emplace_back(uv, vv, nv);
    }
}

std::tuple<EigenVector, EigenVector, EigenMatrix, EigenMatrix> getPrincipalCurvatures(const Mesh &mesh) {
    std::vector<EigenMatrix2> tensors;
    std::vector<LocalFrame> frames;
    getCurvatureTensors(mesh, tensors, frames);

    std::vector<double> ks_max;
    std::vector<double> ks_min;
    std::vector<Vec3> ts_max;
    std::vector<Vec3> ts_min;
    getPrincipalCurvaturesFromTensors(ks_max, ks_min, ts_max, ts_min, tensors, frames);

    const int Nv = (int)ks_max.size();
    EigenVector kv_max(Nv);
    EigenVector kv_min(Nv);
    EigenMatrix tv_max(Nv, 3);
    EigenMatrix tv_min(Nv, 3);
    for (int i = 0; i < Nv; i++) {
        kv_max(i) = ks_max[i];
        kv_min(i) = ks_min[i];
        tv_max.row(i) << ts_max[i].x(), ts_max[i].y(), ts_max[i].z();
        tv_min.row(i) << ts_min[i].x(), ts_min[i].y(), ts_min[i].z();
    }

    return { kv_max, kv_min, tv_max, tv_min };
}

std::tuple<EigenVector, EigenVector, EigenVector, EigenVector, EigenMatrix, EigenMatrix>
getPrincipalCurvaturesWithDerivatives(const Mesh &mesh) {
    std::vector<EigenMatrix2> tensors;
    std::vector<LocalFrame> frames;
    getCurvatureTensors(mesh, tensors, frames);

    std::vector<double> ks_max;
    std::vector<double> ks_min;
    std::vector<Vec3> ts_max;
    std::vector<Vec3> ts_min;
    getPrincipalCurvaturesFromTensors(ks_max, ks_min, ts_max, ts_min, tensors, frames);

    // Computer per-face curvature derivative tensors
    std::vector<EigenVector4> faceTensors;
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

        const Vec3 uf = normalize(p1 - p0);
        const Vec3 wf = normalize(cross(uf, p2 - p0));
        const Vec3 vf = normalize(cross(wf, uf));

        const EigenMatrix3 R0 = rotationFromTwoVectors(wf, n0);
        const EigenMatrix3 R1 = rotationFromTwoVectors(wf, n1);
        const EigenMatrix3 R2 = rotationFromTwoVectors(wf, n2);

        const double dot_e0_u = dot(e0, uf);
        const double dot_e0_v = dot(e0, vf);
        const double dot_e1_u = dot(e1, uf);
        const double dot_e1_v = dot(e1, vf);
        const double dot_e2_u = dot(e2, uf);
        const double dot_e2_v = dot(e2, vf);

        LocalFrame F0 = frames[v0->index()];
        LocalFrame F1 = frames[v1->index()];
        LocalFrame F2 = frames[v2->index()];

        const Vec3 uv0 = std::get<0>(F0);
        const Vec3 vv0 = std::get<1>(F0);
        const Vec3 uv1 = std::get<0>(F1);
        const Vec3 vv1 = std::get<1>(F1);
        const Vec3 uv2 = std::get<0>(F2);
        const Vec3 vv2 = std::get<1>(F2);

        EigenVector2 uf_v0, vf_v0;
        const Vec3 rot0_uf = R0 * EigenVector3(uf);
        const Vec3 rot0_vf = R0 * vf;
        uf_v0 << dot(rot0_uf, uv0), dot(rot0_uf, vv0);
        vf_v0 << dot(rot0_vf, uv0), dot(rot0_vf, vv0);
        EigenVector2 uf_v1, vf_v1;
        const Vec3 rot1_uf = R1 * uf;
        const Vec3 rot1_vf = R1 * vf;
        uf_v1 << dot(rot1_uf, uv1), dot(rot1_uf, vv1);
        vf_v1 << dot(rot1_vf, uv1), dot(rot1_vf, vv1);
        EigenVector2 uf_v2, vf_v2;
        const Vec3 rot2_uf = R2 * uf;
        const Vec3 rot2_vf = R2 * vf;
        uf_v2 << dot(rot2_uf, uv2), dot(rot2_uf, vv2);
        vf_v2 << dot(rot2_vf, uv2), dot(rot2_vf, vv2);

        const EigenMatrix2 M0 = tensors[v0->index()];
        const EigenMatrix2 M1 = tensors[v1->index()];
        const EigenMatrix2 M2 = tensors[v2->index()];

        EigenMatrix A0(4, 4);
        EigenVector b0(4);
        A0.setZero();
        b0.setZero();
        A0.row(0) << dot_e0_u, dot_e0_v, 0.0, 0.0;
        A0.row(1) << 0.0, dot_e0_u, dot_e0_v, 0.0;
        A0.row(2) << 0.0, dot_e0_u, dot_e0_v, 0.0;
        A0.row(3) << 0.0, 0.0, dot_e0_u, dot_e0_v;
        b0 << (M2 - M1) * uf_v0, (M2 - M1) * vf_v0;

        EigenMatrix A1(4, 4);
        EigenVector b1(4);
        A1.setZero();
        b1.setZero();
        A1.row(0) << dot_e1_u, dot_e1_v, 0.0, 0.0;
        A1.row(1) << 0.0, dot_e1_u, dot_e1_v, 0.0;
        A1.row(2) << 0.0, dot_e1_u, dot_e1_v, 0.0;
        A1.row(3) << 0.0, 0.0, dot_e1_u, dot_e1_v;
        b1 << (M0 - M2) * uf_v1, (M0 - M2) * vf_v1;

        EigenMatrix A2(4, 4);
        EigenVector b2(4);
        A2.setZero();
        b2.setZero();
        A2.row(0) << dot_e2_u, dot_e2_v, 0.0, 0.0;
        A2.row(1) << 0.0, dot_e2_u, dot_e2_v, 0.0;
        A2.row(2) << 0.0, dot_e2_u, dot_e2_v, 0.0;
        A2.row(3) << 0.0, 0.0, dot_e2_u, dot_e2_v;
        b2 << (M1 - M0) * uf_v2, (M1 - M0) * vf_v2;

        EigenMatrix A(12, 4);
        EigenVector b(12);
        A << A0, A1, A2;
        b << b0, b1, b2;

        EigenMatrix AA = A.transpose() * A;
        EigenVector Ab = A.transpose() * b;
        Eigen::LLT<EigenMatrix> solver(AA);
        EigenVector4 abcd = solver.solve(Ab);

        faceTensors.push_back(abcd);
        faceFrames.emplace_back(uf, vf, wf);
    }

    std::vector<double> es_max;
    std::vector<double> es_min;
    for (int i = 0; i < mesh.numVertices(); i++) {
        Vertex *vt = mesh.vertex(i);
        const LocalFrame F = frames[vt->index()];
        const Vec3 uv = std::get<0>(F);
        const Vec3 vv = std::get<1>(F);
        const Vec3 nv = std::get<2>(F);

        EigenVector4 M;
        M.setZero();
        double sumWgt = 0.0;
        for (auto it = vt->f_begin(); it != vt->f_end(); ++it) {
            const LocalFrame Ff = faceFrames[it->index()];
            const Vec3 uf = std::get<0>(Ff);
            const Vec3 vf = std::get<1>(Ff);
            const Vec3 nf = std::get<2>(Ff);

            EigenMatrix3 R = rotationFromTwoVectors(nv, nf);
            const Vec3 rot_u = normalize(R * uv);
            const Vec3 rot_v = normalize(R * vv);
            const EigenVector4 Mf = faceTensors[it->index()];

            const double dot_u0 = dot(rot_u, uf);
            const double dot_u1 = dot(rot_u, vf);
            const double dot_v0 = dot(rot_v, uf);
            const double dot_v1 = dot(rot_v, vf);

            const double av = Mf(0) * dot_u0 * dot_u0 * dot_u0 +                                     //
                              3.0 * Mf(1) * dot_u0 * dot_u0 * dot_u1 +                               //
                              3.0 * Mf(2) * dot_u0 * dot_u1 * dot_u1 +                               //
                              Mf(3) * dot_u1 * dot_u1 * dot_u1;                                      //
            const double bv = Mf(0) * dot_u0 * dot_u0 * dot_v0 +                                     //
                              Mf(1) * (dot_u0 * dot_u0 * dot_v1 + 2.0 * dot_u0 * dot_u1 * dot_v0) +  //
                              Mf(2) * (2.0 * dot_u0 * dot_u1 * dot_v1 + dot_u1 * dot_u1 * dot_v0) +  //
                              Mf(3) * dot_u1 * dot_u1 * dot_v1;                                      //
            const double cv = Mf(0) * dot_u0 * dot_v0 * dot_v0 +                                     //
                              Mf(1) * (2.0 * dot_u0 * dot_v0 * dot_v1 + dot_u1 * dot_v0 * dot_v0) +  //
                              Mf(2) * (dot_u0 * dot_v1 * dot_v1 + 2.0 * dot_u1 * dot_v0 * dot_v1) +  //
                              Mf(3) * dot_u1 * dot_v1 * dot_v1;                                      //
            const double dv = Mf(0) * dot_v0 * dot_v0 * dot_v0 +                                     //
                              3.0 * Mf(1) * dot_v0 * dot_v0 * dot_v1 +                               //
                              3.0 * Mf(2) * dot_v0 * dot_v1 * dot_v1 +                               //
                              Mf(3) * dot_v1 * dot_v1 * dot_v1;                                      //
            EigenVector4 Mv;
            Mv << av, bv, cv, dv;

            const double weight = it->area();
            M += weight * Mv;
            sumWgt += weight;
        }
        M /= sumWgt;

        const Vec3 t_max = ts_max[i];
        const Vec3 t_min = ts_min[i];

        const double u0_max = dot(t_max, uv);
        const double u1_max = dot(t_max, vv);
        const double u0_min = dot(t_min, uv);
        const double u1_min = dot(t_min, vv);
        const double e_max = M(0) * u0_max * u0_max * u0_max +        //
                             3.0 * M(1) * u0_max * u0_max * u1_max +  //
                             3.0 * M(2) * u0_max * u1_max * u1_max +  //
                             M(3) * u1_max * u1_max * u1_max;         //
        const double e_min = M(0) * u0_min * u0_min * u0_min +        //
                             3.0 * M(1) * u0_min * u0_min * u1_min +  //
                             3.0 * M(2) * u0_min * u1_min * u1_min +  //
                             M(3) * u1_min * u1_min * u1_min;         //
        es_max.push_back(e_max);
        es_min.push_back(e_min);
    }

    const int Nv = (int)mesh.numVertices();
    EigenVector kv_max(Nv), kv_min(Nv), ev_max(Nv), ev_min(Nv);
    EigenMatrix tv_max(Nv, 3), tv_min(Nv, 3);
    for (int i = 0; i < Nv; i++) {
        kv_max(i) = ks_max[i];
        kv_min(i) = ks_min[i];
        ev_max(i) = es_max[i];
        ev_min(i) = es_min[i];
        tv_max.row(i) << ts_max[i].x(), ts_max[i].y(), ts_max[i].z();
        tv_min.row(i) << ts_min[i].x(), ts_min[i].y(), ts_min[i].z();
    }

    return { kv_max, kv_min, ev_max, ev_min, tv_max, tv_min };
}

}  // namespace tinymesh

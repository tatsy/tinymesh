#define TINYMESH_API_EXPORT
#include "ops.h"

#include <Eigen/Eigenvalues>
#include <Eigen/Cholesky>
#include <Eigen/SparseCholesky>
#include <Spectra/SymEigsSolver.h>
#include <Spectra/SymEigsShiftSolver.h>
#include <Spectra/MatOp/DenseSymMatProd.h>
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/MatOp/SparseSymShiftSolve.h>

#include "core/vertex.h"
#include "core/halfedge.h"
#include "core/face.h"

namespace tinymesh {

namespace {

double secondForm(const EigenMatrix2 &M, const EigenVector2 &u, const EigenVector2 &v) {
    return ((u.transpose() * M) * v).value();
}

double thirdForm(const EigenVector4 &C, const EigenVector2 &u, const EigenVector2 &v, const EigenVector2 &w) {
    double ret = 0.0;
    ret += C(0) * u(0) * v(0) * w(0);
    ret += C(1) * (u(1) * v(0) * w(0) + u(0) * v(1) * w(0) + u(0) * v(0) * w(1));
    ret += C(2) * (u(1) * v(1) * w(0) + u(0) * v(1) * w(1) + u(1) * v(0) * w(1));
    ret += C(3) * u(1) * v(1) * w(1);
    return ret;
}

void getPrincipalCurvaturesFromTensors(std::vector<double> &ks_max, std::vector<double> &ks_min,
                                       std::vector<Vec3> &ts_max, std::vector<Vec3> &ts_min,
                                       const std::vector<EigenMatrix2> &tensors,
                                       const std::vector<LocalFrame> &frames) {
    const int Nv = (int)tensors.size();
    ks_max.resize(Nv);
    ks_min.resize(Nv);
    ts_max.resize(Nv);
    ts_min.resize(Nv);

    for (size_t i = 0; i < Nv; i++) {
        const EigenMatrix2 M = tensors[i];
        const LocalFrame F = frames[i];
        const Vec3 up = std::get<0>(F);
        const Vec3 vp = std::get<1>(F);
        const Vec3 np = std::get<2>(F);

        double c = 1.0, s = 0.0, tt = 0.0;
        if (M(0, 1) != 0.0) {
            const double h = 0.5 * (M(1, 1) - M(0, 0)) / M(0, 1);
            tt = (h < 0.0) ? 1.0 / (h - std::sqrt(1.0 + h * h)) : 1.0 / (h + std::sqrt(1.0 + h * h));
            c = 1.0 / std::sqrt(1.0 + tt * tt);
            s = tt * c;
        }

        double k_min = M(0, 0) - tt * M(0, 1);
        double k_max = M(1, 1) + tt * M(0, 1);
        Vec3 t_min, t_max;
        if (k_max > k_min) {
            t_max = c * up - s * vp;
        } else {
            std::swap(k_min, k_max);
            t_max = s * up + c * vp;
        }
        t_min = normalize(cross(np, t_max));

        // Ignore umbilical points
        if (std::abs(k_min - k_max) < 1.0e-6) {
            t_min = Vec3(0.0);
            t_max = Vec3(0.0);
        }

        ks_max[i] = k_max;
        ks_min[i] = k_min;
        ts_max[i] = t_max;
        ts_min[i] = t_min;
    }
}

void smoothingTensors(std::vector<EigenMatrix2> &tensors, std::vector<LocalFrame> &frames, const Mesh &mesh) {
    const int Nv = (int)tensors.size();
    const double avgLen = mesh.getMeanEdgeLength();
    std::vector<EigenMatrix2> smoothTensors(Nv);
    for (size_t i = 0; i < Nv; i++) {
        const Vertex *vertex = mesh.vertex(i);
        const LocalFrame &F = frames[i];
        const Vec3 &u = std::get<0>(F);
        const Vec3 &v = std::get<1>(F);
        const Vec3 &n = std::get<2>(F);

        EigenMatrix2 M = tensors[i];
        double sumWgt = 1.0;
        for (auto it = vertex->v_begin(); it != vertex->v_end(); ++it) {
            const double dist = length(vertex->pos() - it->pos()) / avgLen;
            const double weight = std::exp(-0.5 * dist * dist);

            const EigenMatrix2 &Mp = tensors[it->index()];
            const LocalFrame &Fp = frames[it->index()];
            const Vec3 &up = std::get<0>(Fp);
            const Vec3 &vp = std::get<1>(Fp);
            const Vec3 &np = std::get<2>(Fp);

            const EigenMatrix3 R = rotationFromTwoVectors(n, np);
            const Vec3 rot_u = normalize(R * u);
            const Vec3 rot_v = normalize(R * v);

            EigenVector2 new_u, new_v;
            new_u << dot(rot_u, up), dot(rot_u, vp);
            new_v << dot(rot_v, up), dot(rot_v, vp);

            const double e = (new_u.transpose() * (Mp * new_u)).value();
            const double f = (new_u.transpose() * (Mp * new_v)).value();
            const double g = (new_v.transpose() * (Mp * new_v)).value();

            EigenMatrix2 new_Mp;
            new_Mp << e, f, f, g;

            M += weight * new_Mp;
            sumWgt += weight;
        }
        smoothTensors[i] = M / sumWgt;
    }
    tensors = smoothTensors;
}

}  // anonymous namespace

EigenMatrix getHeatKernelSignatures(const EigenSparseMatrix &L, int K, int nTimes) {
    Assertion(K < L.rows(), "Eigenvalues less than matrix size can only be requested!");
    const int ncv = std::min(K * 8, (int)L.rows());

    // Spectra::SparseSymMatProd<FloatType> op(L);
    // Spectra::SymEigsSolver<Spectra::SparseSymMatProd<FloatType>> eigs(op, K, ncv);
    // eigs.init();

    // eigs.compute(Spectra::SortRule::SmallestMagn, 200, 1.0e-4, Spectra::SortRule::SmallestMagn);
    // if (eigs.info() != Spectra::CompInfo::Successful) {
    //     Error("Eigen decomposition failed!");
    // }

    Spectra::SparseSymShiftSolve<FloatType> op(L);
    Spectra::SymEigsShiftSolver<Spectra::SparseSymShiftSolve<FloatType>> eigs(op, K, ncv, -1.0e-8);
    eigs.init();
    eigs.compute(Spectra::SortRule::LargestMagn, 200, 1.0e-4, Spectra::SortRule::SmallestMagn);
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
    // Smoothing vertex normals
    std::vector<Vec3> normals(mesh.numVertices());
    for (size_t i = 0; i < mesh.numVertices(); i++) {
        normals[i] = mesh.vertex(i)->normal();
    }

    const double avgLen = mesh.getMeanEdgeLength();
    std::vector<Vec3> smoothNorms(mesh.numVertices());
    for (size_t i = 0; i < mesh.numVertices(); i++) {
        Vertex *vertex = mesh.vertex(i);
        Vec3 norm = normals[i];
        double sumWgt = 1.0;
        for (auto it = vertex->v_begin(); it != vertex->v_end(); ++it) {
            const double dist = length(vertex->pos() - it->pos()) / avgLen;
            const double weight = std::exp(-0.5 * dist * dist);
            norm += normals[it->index()] * weight;
            sumWgt += weight;
        }
        smoothNorms[i] = normalize(norm / sumWgt);
    }

    // Compute per-triangle tensors
    const int Nf = (int)mesh.numFaces();
    std::vector<EigenMatrix2> faceTensors(Nf);
    std::vector<LocalFrame> faceFrames(Nf);
    for (int i = 0; i < Nf; i++) {
        Face *f = mesh.face(i);
        std::vector<Vertex *> vertices;
        for (auto it = f->v_begin(); it != f->v_end(); ++it) {
            vertices.push_back(it.ptr());
        }
        Assertion(vertices.size() == 3, "Non-triangle face detected!");

        Vertex *v0 = vertices[0];
        Vertex *v1 = vertices[1];
        Vertex *v2 = vertices[2];
        std::array<Vec3, 3> ps = {
            v0->pos(),
            v1->pos(),
            v2->pos(),
        };
        std::array<Vec3, 3> ns = {
            v0->normal(),
            v1->normal(),
            v2->normal(),
        };
        std::array<Vec3, 3> es = {
            ps[2] - ps[1],
            ps[0] - ps[2],
            ps[1] - ps[0],
        };
        const Vec3 uf = normalize(es[0]);
        const Vec3 nf = normalize(cross(es[0], es[1]));
        const Vec3 vf = normalize(cross(nf, uf));

        EigenMatrix A = EigenMatrix::Zero(3, 3);
        EigenVector b = EigenVector::Zero(3);
        for (int k = 0; k < 3; k++) {
            const double u = dot(es[k], uf);
            const double v = dot(es[k], vf);
            A(0, 0) += u * u;
            A(0, 1) += u * v;
            A(2, 2) += v * v;

            const int kp = (k - 1 + 3) % 3;
            const int kn = (k + 1) % 3;
            const Vec3 dn = ns[kp] - ns[kn];
            const double dnu = dot(dn, uf);
            const double dnv = dot(dn, vf);
            b(0) += dnu * u;
            b(1) += dnu * v + dnv * u;
            b(2) += dnv * v;
        }
        A(1, 1) = A(0, 0) + A(2, 2);
        A(1, 2) = A(0, 1);

        EigenVector efg = A.selfadjointView<Eigen::Upper>().ldlt().solve(b);
        EigenMatrix2 M;
        M << efg(0), efg(1), efg(1), efg(2);
        faceTensors[i] = M;
        faceFrames[i] = LocalFrame(uf, vf, nf);
    }

    tensors.clear();
    frames.clear();
    tensors.resize(mesh.numVertices());
    frames.resize(mesh.numVertices());
    for (int i = 0; i < mesh.numVertices(); i++) {
        Vertex *vertex = mesh.vertex(i);
        const Vec3 np = vertex->normal();
        Vec3 up, vp;
        for (auto it = vertex->v_begin(); it != vertex->v_end(); ++it) {
            up = normalize(it->pos() - vertex->pos());
            break;
        }
        vp = normalize(cross(np, up));

        EigenMatrix2 M = EigenMatrix2::Zero();
        double sumWgt = 0.0;
        for (auto it = vertex->f_begin(); it != vertex->f_end(); ++it) {
            const LocalFrame &F = faceFrames[it->index()];
            const Vec3 &uf = std::get<0>(F);
            const Vec3 &vf = std::get<1>(F);
            const Vec3 &nf = std::get<2>(F);

            EigenMatrix3 R = rotationFromTwoVectors(np, nf);
            const Vec3 rot_up = R * up;
            const Vec3 rot_vp = R * vp;
            const EigenMatrix2 &Mf = faceTensors[it->index()];

            EigenVector2 new_up, new_vp;
            new_up << dot(rot_up, uf), dot(rot_up, vf);
            new_vp << dot(rot_vp, uf), dot(rot_vp, vf);

            const double ep = secondForm(Mf, new_up, new_up);
            const double fp = secondForm(Mf, new_up, new_vp);
            const double gp = secondForm(Mf, new_vp, new_vp);

            EigenMatrix2 Mp;
            Mp << ep, fp, fp, gp;

            const double weight = vertex->volonoiArea(it.ptr());
            M += weight * Mp;
            sumWgt += weight;
        }
        M /= sumWgt;

        tensors[i] = M;
        frames[i] = LocalFrame(up, vp, np);
    }
}

std::tuple<EigenVector, EigenVector, EigenMatrix, EigenMatrix> getPrincipalCurvatures(const Mesh &mesh,
                                                                                      bool smoothTensors) {
    std::vector<EigenMatrix2> tensors;
    std::vector<LocalFrame> frames;
    getCurvatureTensors(mesh, tensors, frames);

    if (smoothTensors) {
        smoothingTensors(tensors, frames, mesh);
    }

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
getPrincipalCurvaturesWithDerivatives(const Mesh &mesh, bool smoothTensors) {
    std::vector<EigenMatrix2> tensors;
    std::vector<LocalFrame> frames;
    getCurvatureTensors(mesh, tensors, frames);

    if (smoothTensors) {
        smoothingTensors(tensors, frames, mesh);
    }

    std::vector<double> ks_max;
    std::vector<double> ks_min;
    std::vector<Vec3> ts_max;
    std::vector<Vec3> ts_min;
    getPrincipalCurvaturesFromTensors(ks_max, ks_min, ts_max, ts_min, tensors, frames);

    // Computer per-face curvature derivative tensors
    const int Nf = (int)mesh.numFaces();
    std::vector<EigenVector4> faceTensors(Nf);
    std::vector<LocalFrame> faceFrames(Nf);
    for (int i = 0; i < Nf; i++) {
        Face *f = mesh.face(i);
        std::vector<Vertex *> vertices;
        for (auto it = f->v_begin(); it != f->v_end(); ++it) {
            vertices.push_back(it.ptr());
        }
        Assertion(vertices.size() == 3, "Non-triangle face detected!");

        Vertex *v0 = vertices[0];
        Vertex *v1 = vertices[1];
        Vertex *v2 = vertices[2];
        std::array<Vec3, 3> ps = {
            v0->pos(),
            v1->pos(),
            v2->pos(),
        };
        std::array<Vec3, 3> es = {
            ps[2] - ps[1],
            ps[0] - ps[2],
            ps[1] - ps[0],
        };
        std::array<int, 3> ids = { v0->index(), v1->index(), v2->index() };

        const Vec3 uf = normalize(es[0]);
        const Vec3 nf = normalize(cross(es[0], es[1]));
        const Vec3 vf = normalize(cross(nf, uf));

        std::array<EigenVector3, 3> Mfs;
        for (int k = 0; k < 3; k++) {
            const Vec3 &up = ts_max[ids[k]];
            const Vec3 &vp = ts_min[ids[k]];
            const Vec3 &np = std::get<2>(frames[ids[k]]);
            const EigenMatrix3 rot = rotationFromTwoVectors(nf, np);

            EigenVector2 new_uf, new_vf;
            const Vec3 rot_uf = rot * uf;
            const Vec3 rot_vf = rot * vf;
            new_uf << dot(rot_uf, up), dot(rot_uf, vp);
            new_vf << dot(rot_vf, up), dot(rot_vf, vp);

            EigenMatrix2 M;
            M << ks_max[ids[k]], 0.0, 0.0, ks_min[ids[k]];
            Mfs[k] << secondForm(M, new_uf, new_uf),  //
                secondForm(M, new_uf, new_vf),        //
                secondForm(M, new_vf, new_vf);
        }

        EigenMatrix A = EigenMatrix::Zero(4, 4);
        EigenVector b = EigenVector::Zero(4);
        for (int k = 0; k < 3; k++) {
            const double u = dot(es[k], uf);
            const double v = dot(es[k], vf);
            A(0, 0) += u * u;
            A(0, 1) += u * v;
            A(3, 3) += v * v;

            const int kp = (k - 1 + 3) % 3;
            const int kn = (k + 1) % 3;
            const EigenVector3 dfcurv = Mfs[kp] - Mfs[kn];
            b(0) += u * dfcurv(0);
            b(1) += v * dfcurv(0) + 2.0 * u * dfcurv(1);
            b(2) += 2.0 * v * dfcurv(1) + u * dfcurv(2);
            b(3) += v * dfcurv(2);
        }
        A(1, 1) = 2.0 * A(0, 0) + A(3, 3);
        A(1, 2) = 2.0 * A(0, 1);
        A(2, 2) = A(0, 0) + 2.0 * A(3, 3);
        A(2, 3) = A(0, 1);

        EigenVector4 abcd = A.selfadjointView<Eigen::Upper>().ldlt().solve(b);

        faceTensors[i] = abcd;
        faceFrames[i] = LocalFrame(uf, vf, nf);
    }

    // Compute per-vertex curvature derivatives
    std::vector<double> es_max(mesh.numVertices(), 0.0);
    std::vector<double> es_min(mesh.numVertices(), 0.0);
    for (size_t i = 0; i < mesh.numVertices(); i++) {
        Vertex *vertex = mesh.vertex(i);
        const Vec3 &up = std::get<0>(frames[i]);
        const Vec3 &vp = std::get<1>(frames[i]);
        const Vec3 &np = std::get<2>(frames[i]);

        EigenVector4 M = EigenVector4::Zero();
        double sumWgt = 0.0;
        for (auto it = vertex->f_begin(); it != vertex->f_end(); ++it) {
            const LocalFrame &F = faceFrames[it->index()];
            const Vec3 &uf = std::get<0>(F);
            const Vec3 &vf = std::get<1>(F);
            const Vec3 &nf = std::get<2>(F);

            const EigenMatrix3 R = rotationFromTwoVectors(np, nf);
            const Vec3 rot_up = R * up;
            const Vec3 rot_vp = R * vp;
            const EigenVector4 &Mf = faceTensors[it->index()];

            EigenVector2 new_up, new_vp;
            new_up << dot(rot_up, uf), dot(rot_up, vf);
            new_vp << dot(rot_vp, uf), dot(rot_vp, vf);

            EigenVector4 Mp;
            Mp(0) = thirdForm(Mf, new_up, new_up, new_up);
            Mp(1) = thirdForm(Mf, new_up, new_up, new_vp);
            Mp(2) = thirdForm(Mf, new_up, new_vp, new_vp);
            Mp(3) = thirdForm(Mf, new_vp, new_vp, new_vp);

            const double weight = vertex->volonoiArea(it.ptr());
            M += weight * Mp;
            sumWgt += weight;
        }
        M /= sumWgt;

        const Vec3 t_max = ts_max[i];
        const Vec3 t_min = ts_min[i];

        EigenVector2 t_max_2d, t_min_2d;
        t_max_2d << dot(t_max, up), dot(t_max, vp);
        t_min_2d << dot(t_min, up), dot(t_min, vp);

        es_max[i] = thirdForm(M, t_max_2d, t_max_2d, t_max_2d);
        es_min[i] = thirdForm(M, t_min_2d, t_min_2d, t_min_2d);
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

EigenMatrix getFeatureLineField(const Mesh &mesh, bool smoothTensors) {
    const auto ret = getFeatureLineFieldWithFlags(mesh, smoothTensors);
    return std::get<0>(ret);
}

std::tuple<EigenMatrix, EigenVector> getFeatureLineFieldWithFlags(const Mesh &mesh, bool smoothTensors) {
    const auto cuv = getPrincipalCurvaturesWithDerivatives(mesh, smoothTensors);
    const EigenVector kv_max = std::get<0>(cuv);
    const EigenVector kv_min = std::get<1>(cuv);
    const EigenVector ev_max = std::get<2>(cuv);
    const EigenVector ev_min = std::get<3>(cuv);
    const EigenMatrix tv_max = std::get<4>(cuv);
    const EigenMatrix tv_min = std::get<5>(cuv);

    const int NONE = 0;
    const int VALLEY = 1;
    const int RIDGE = -1;

    auto checkRV = [&](const Vertex *v1, const Vertex *v2) -> int {
        const int i1 = v1->index();
        const int i2 = v2->index();
        const Vec3 p1 = v1->pos();
        const Vec3 p2 = v2->pos();

        const auto t1 = tv_max.row(i1);
        Vec3 t1_max(t1(0), t1(1), t1(2));
        const auto t2 = tv_max.row(i2);
        Vec3 t2_max(t2(0), t2(1), t2(2));

        const double k1_max = kv_max(i1);
        const double k1_min = kv_min(i1);
        const double k2_max = kv_max(i2);
        const double k2_min = kv_min(i2);

        double e1_max = ev_max(i1);
        double e2_max = ev_max(i2);
        if (dot(t1_max, t2_max) < 0.0) {
            t2_max = -1.0 * t2_max;
            e2_max = -1.0 * e2_max;
        }

        // Ignore umbilical points
        if (std::abs(k1_max - k1_min) < 1.0e-2 || std::abs(k2_max - k2_min) < 1.0e-2) {
            return NONE;
        }

        if (e1_max * e2_max >= 0.0) {
            return NONE;
        }

        const double c1 = e1_max * dot(p2 - p1, t1_max);
        const double c2 = e2_max * dot(p1 - p2, t2_max);

        if (k1_max > std::abs(k1_min) && k2_max > std::abs(k2_min)) {
            if (c1 > 0.0 && c2 > 0.0) {
                return RIDGE;
            }
        }

        if (k1_max < std::abs(k1_min) && k2_max < std::abs(k2_min)) {
            if (c1 < 0.0 && c2 < 0.0) {
                return VALLEY;
            }
        }

        return NONE;
    };

    // Mark ridge/valley faces
    std::vector<int> faceRVFlags(mesh.numFaces());
    for (size_t i = 0; i < mesh.numFaces(); i++) {
        const Face *f = mesh.face(i);
        std::vector<const Vertex *> vs;
        for (auto it = f->v_begin(); it != f->v_end(); ++it) {
            vs.push_back(it.ptr());
        }
        Assertion(vs.size() == 3, "Non-triangle face detected!");

        const int rv0 = checkRV(vs[0], vs[1]);
        const int rv1 = checkRV(vs[1], vs[2]);
        const int rv2 = checkRV(vs[2], vs[0]);
        if (rv0 != NONE || rv1 != NONE || rv2 != NONE) {
            if (rv0 >= 0 && rv1 >= 0 && rv2 >= 0) {
                faceRVFlags[i] = VALLEY;
            }

            if (rv0 <= 0 && rv1 <= 0 && rv2 <= 0) {
                faceRVFlags[i] = RIDGE;
            }
        }
    }

    // Detect ridge/valley vertices
    std::vector<Vec3> rvDirections(mesh.numVertices(), Vec3(0.0));
    EigenVector ridgeValleyFlags(mesh.numVertices());
    ridgeValleyFlags.setZero();
    for (size_t i = 0; i < mesh.numFaces(); i++) {
        const Face *f = mesh.face(i);
        if (faceRVFlags[f->index()] != NONE) {
            continue;
        }

        int count = 0;
        int opposite = -1;
        for (auto it = f->f_begin(); it != f->f_end(); ++it) {
            if (faceRVFlags[it->index()] != NONE) {
                opposite = it->index();
                count++;
            }
        }

        if (count != 1) {
            continue;
        }

        const Face *g = mesh.face(opposite);
        Assertion(faceRVFlags[g->index()] != NONE, "Something is wrong!");

        bool valid = false;
        for (auto it = f->he_begin(); it != f->he_end(); ++it) {
            if (it->rev()->face() == g) {
                const Vertex *v0 = it->src();
                const Vertex *v1 = it->dst();
                const Vertex *v2 = it->next()->dst();

                const Vec3 p0 = v0->pos();
                const Vec3 p1 = v1->pos();
                const Vec3 p2 = v2->pos();
                const Vec3 e1 = p1 - p0;
                const Vec3 e2 = p2 - p0;
                const Vec3 n = normalize(cross(e1, e2));
                const Vec3 d = normalize(cross(n, e1));

                const int orient = faceRVFlags[g->index()];
                rvDirections[v0->index()] += (double)orient * d;
                rvDirections[v1->index()] += (double)orient * d;
                ridgeValleyFlags(v0->index()) += 1.0;
                ridgeValleyFlags(v1->index()) += 1.0;

                valid = true;
                break;
            }
        }

        Assertion(valid, "Something is wrong!");
    }

    for (int i = 0; i < mesh.numVertices(); i++) {
        if (ridgeValleyFlags[i] != 0) {
            rvDirections[i] /= ridgeValleyFlags(i);
            ridgeValleyFlags(i) = 1.0;
            if (length(rvDirections[i]) > 1.0e-8) {
                rvDirections[i] = normalize(rvDirections[i]);
            }
        }
    }

    // Solve Laplace equation
    std::vector<EigenTriplet> tripNonRV;
    std::vector<EigenTriplet> tripRV;
    std::unordered_map<uint32_t, uint32_t> uniqueNonRV;
    std::unordered_map<uint32_t, uint32_t> uniqueRV;
    int countNonRV = 0;
    int countRV = 0;
    for (size_t i = 0; i < mesh.numVertices(); i++) {
        const Vertex *v = mesh.vertex(i);
        if (ridgeValleyFlags(i) == 0.0) {
            // non-ridge-valley vertex
            if (uniqueNonRV.count(i) == 0) {
                uniqueNonRV[i] = countNonRV;
                countNonRV++;
            }
            const int row = uniqueNonRV[i];
            double sumWgt = 0.0;
            for (auto it = v->he_begin(); it != v->he_end(); ++it) {
                Vertex *u = it->dst();
                const double weight = it->cotWeight();
                if (ridgeValleyFlags(u->index()) != 0.0) {
                    // ridge-valley vertex
                    if (uniqueRV.count(u->index()) == 0) {
                        uniqueRV[u->index()] = countRV;
                        countRV++;
                    }
                    const int col = uniqueRV[u->index()];
                    tripRV.emplace_back(row, col, -weight);
                } else {
                    // non-ridge-valley vertex
                    if (uniqueNonRV.count(u->index()) == 0) {
                        uniqueNonRV[u->index()] = countNonRV;
                        countNonRV++;
                    }
                    const int col = uniqueNonRV[u->index()];
                    tripNonRV.emplace_back(row, col, -weight);
                }
                sumWgt += weight;
            }
            tripNonRV.emplace_back(row, row, sumWgt);
        }
    }

    EigenSparseMatrix LL(countNonRV, countNonRV);
    LL.setFromTriplets(tripNonRV.begin(), tripNonRV.end());
    EigenSparseMatrix SS(countNonRV, countRV);
    SS.setFromTriplets(tripRV.begin(), tripRV.end());

    EigenMatrix bbNonRV(countNonRV, 3);
    bbNonRV.setZero();
    for (auto it : uniqueNonRV) {
        const Vec3 &v = rvDirections[it.first];
        bbNonRV.row(it.second) << v.x(), v.y(), v.z();
    }

    EigenMatrix xxRV(countRV, 3);
    xxRV.setZero();
    for (auto it : uniqueRV) {
        const Vec3 &v = rvDirections[it.first];
        xxRV.row(it.second) << v.x(), v.y(), v.z();
    }

    EigenMatrix bb = -SS * xxRV;
    Eigen::SimplicialLDLT<EigenSparseMatrix> solver(LL);
    EigenMatrix xx = solver.solve(bb);

    for (auto it : uniqueNonRV) {
        const int i = it.second;
        rvDirections[it.first] = normalize(Vec3(xx(i, 0), xx(i, 1), xx(i, 2)));
    }

    // Convert std::vector to EigenMatrix
    EigenMatrix ret(mesh.numVertices(), 3);
    for (size_t i = 0; i < mesh.numVertices(); i++) {
        const Vec3 &d = rvDirections[i];
        ret.row(i) << d.x(), d.y(), d.z();
    }
    return std::make_tuple(ret, ridgeValleyFlags);
}

}  // namespace tinymesh

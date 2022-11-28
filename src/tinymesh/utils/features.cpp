#define TINYMESH_API_EXPORT
#include "utils.h"
#include <Spectra/SymEigsSolver.h>
#include <Spectra/MatOp/SparseGenMatProd.h>

namespace tinymesh {

EigenMatrix getHeatKernelSignatures(EigenSparseMatrix &L, int K, int nTimes) {
    Spectra::SparseGenMatProd<FloatType> op(L);
    Spectra::SymEigsSolver<Spectra::SparseGenMatProd<FloatType>> eigs(op, K, K * 2);
    eigs.init();

    int nconv = eigs.compute(Spectra::SortRule::SmallestMagn, 1000, 1.0e-10, Spectra::SortRule::SmallestMagn);
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

}  // namespace tinymesh

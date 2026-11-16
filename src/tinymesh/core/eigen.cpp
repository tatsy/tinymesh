#include "eigen.h"

#include <Eigen/SVD>

namespace tinymesh {

void eigenSVD(const EigenMatrix &A, EigenMatrix &U, EigenVector &sigma, EigenMatrix &V) {
    Eigen::JacobiSVD<EigenMatrix> svd(A, Eigen::ComputeFullU | Eigen::ComputeFullV);
    U = svd.matrixU();
    sigma = svd.singularValues();
    V = svd.matrixV();
}

}  // namespace tinymesh

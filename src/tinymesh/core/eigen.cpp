#include "eigen.h"

#include <Eigen/SVD>

void eigenSVD(const Eigen::MatrixXd &A, Eigen::MatrixXd &U, Eigen::VectorXd &sigma, Eigen::MatrixXd &V) {
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeFullU | Eigen::ComputeFullV);
    U = svd.matrixU();
    sigma = svd.singularValues();
    V = svd.matrixV();
}

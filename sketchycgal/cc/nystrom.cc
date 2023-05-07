
#include "third_party/eigen3/Eigen/SparseCore"
#include "third_party/eigen3/Eigen/SparseCholesky"
#include "third_party/eigen3/Eigen/Core"
#include "third_party/eigen3/Eigen/Cholesky"
#include "third_party/eigen3/Eigen/SVD"
#include "third_party/eigen3/Eigen/QR"
#include "third_party/eigen3/Eigen/LU"
#include <iostream>
#include <vector>
#include <random>
#include <string>
#include <cmath>
#include "nystrom.h"

#define EPS 2.2e-16

NystromSketch::NystromSketch(int n, int R){
  _n = n;
  _R = R;
  if (R > n){
    std::cerr << "sketch-size R=" << R << " is larger than n=" << n << std::endl;
    exit(1);
  }
  std::default_random_engine generator;
  std::normal_distribution<double> distribution(0.0, 1.0);
  std::cout << "Filling the Omega the matrix with random normal values" << std::endl;
  _Omega = new Eigen::MatrixXd(n, R);
  _Omega->setRandom();
  // for (int i = 0; i < n; ++i) {
  //     for (int j = 0; j < R; ++j) {
  //         (*_Omega)(i, j) = distribution(generator);
  //     }
  // }
  std::cout << "Filling the S matrix with zero values" << std::endl;
  _S = new Eigen::MatrixXd(n ,R);
  _S->setZero();
}

void NystromSketch::update( Eigen::ArrayXd& v, double eta){
  Eigen::MatrixXd new_S = (1-eta)*(*_S) + eta*( v.matrix() * ( v.matrix().transpose() * (*_Omega) ) );
  // std::cout << "norm of new_S - S = " << (new_S - (*_S)).norm() << std::endl;
  (*_S) = new_S;
}

std::pair<Eigen::MatrixXd, Eigen::MatrixXd > NystromSketch::reconstruct(){
  Eigen::VectorXd col_norms = _S->colwise().norm();
  double sigma = sqrt(_n) * EPS * col_norms.maxCoeff();
  // std::cout << "sigma val " << sigma << std::endl;
  (*_S) = (*_S) + sigma * (*_Omega);
  Eigen::MatrixXd B = (*_Omega).transpose() * (*_S);
  B = 0.5*( B.array() + B.array().transpose() );
  // skip zero norm check

  /* https://eigen.tuxfamily.org/dox/classEigen_1_1LLT.html
  Performance: for best performance, it is recommended to use a
  column-major storage format with the Lower triangular part (the default),
  or, equivalently, a row-major storage format with the Upper triangular part.
  Otherwise, you might get a 20% slowdown for the full factorization step, and
  rank-updates can be up to 3 times slower.*/
  Eigen::LLT<Eigen::MatrixXd> lltOfA(B); // LLT for lower triangular
  Eigen::MatrixXd L = lltOfA.matrixL();

  Eigen::MatrixXd U;
  Eigen::ArrayXd sv_array;
  // assume S and C have been initialized with appropriate values
  // std::cout << "Dimensions of L : " << L.rows() << " " << L.cols() << std::endl;
  // Eigen::BDCSVD<Eigen::MatrixXd> svd_S((*_S), Eigen::ComputeThinU | Eigen::ComputeThinV);
  // Eigen::MatrixXd temp = svd_S.transpose().solve(L.transpose()).transpose();

  Eigen::MatrixXd temp = (*_S) * L.inverse();
  // std::cout << "Dimensions of temp : " << temp.rows() << " " << temp.cols() << std::endl;

  Eigen::BDCSVD<Eigen::MatrixXd> svd(temp, Eigen::ComputeThinU | Eigen::ComputeThinV);
  U = svd.matrixU();
  sv_array = svd.singularValues();
  // std::cout << "sv array: \n" << sv_array << std::endl;
  Eigen::MatrixXd Delta = ( sv_array.square() - sigma ).cwiseMax(0).matrix().asDiagonal();
  // std::cout << "Delta matrix: \n" << Delta << std::endl;
  return std::pair<Eigen::MatrixXd, Eigen::MatrixXd>(U, Delta);

}


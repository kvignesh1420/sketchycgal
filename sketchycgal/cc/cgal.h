#ifndef SKETCHYCGAL_CC_CGAL_H
#define SKETCHYCGAL_CC_CGAL_H

#include "third_party/eigen3/Eigen/SparseCore"
#include "third_party/eigen3/Eigen/SparseCholesky"
#include "third_party/eigen3/Eigen/Core"
#include <vector>
#include <string>
#include <math.h>
#include <random>
#include "nystrom.h"


class SketchyCGAL{
  public:
    SketchyCGAL();
    ~SketchyCGAL(){
      delete _A;
      delete _C;
    };
    /* Allow the user to pass a N x 3 matrix of (row, col, wt) entries of the
     adjacency matrix and configure the state on demand via the `setup` method. */
    void setup(char* filepath);

    // run the sketchyCGAL iterations
    void run();

    Eigen::SparseMatrix<double> getAdjacencyMatrix();
    Eigen::SparseMatrix<double, Eigen::RowMajor> getLaplacian();
  private:
    Eigen::SparseMatrix<double, Eigen::RowMajor>* computeLaplacian(Eigen::SparseMatrix<double>* A);

    /* primitives for matrix - vector manipulations */
    Eigen::VectorXd primitive1(const Eigen::VectorXd& x);
    Eigen::VectorXd primitive2(const Eigen::VectorXd& y, const Eigen::VectorXd& x);
    Eigen::VectorXd primitive3(const Eigen::VectorXd& x);

    /* Lanczos approximation of eigenvector*/
    std::pair<Eigen::VectorXd, double> ApproxMinEvecLanczosSE(Eigen::VectorXd& vt, int n, int q);
    std::pair<Eigen::MatrixXd, Eigen::VectorXd> cgal_eig(const Eigen::SparseMatrix<double>& X);

    /* Get the value of max-cut */
    double getCutValue(Eigen::MatrixXd& U);

    /* private variables */
    Eigen::SparseMatrix<double>* _A = nullptr;
    // using _C to represent the laplacian for notational consistency
    Eigen::SparseMatrix<double, Eigen::RowMajor>* _C = nullptr;
    // scaling factors;
    double _SCALE_A = 1.0;
    double _SCALE_X = 1.0;
    double _SCALE_C = 1.0;
    double _RESCALE_OBJ = 1.0;
    double _RESCALE_FEAS = 1.0;
    // runner config
    int _MAX_ITERS = 1000;
    int _R = 10;
    double _beta0 = 1.0;
    double _K = INFINITY;
    double _TOLERANCE = 0.1;
    // nystrom sketch object;
    NystromSketch* _nysketch = nullptr;
};


#endif
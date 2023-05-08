#ifndef SKETCHYCGAL_CC_NYSTROM_H
#define SKETCHYCGAL_CC_NYSTROM_H

#include "third_party/eigen3/Eigen/SparseCore"
#include "third_party/eigen3/Eigen/SparseCholesky"
#include "third_party/eigen3/Eigen/Core"
#include <vector>
#include <random>
#include <string>
#include <cmath>


class NystromSketch{
  public:
    NystromSketch( int n, int R);
    ~NystromSketch(){
        delete _Omega;
        delete _S;
    };

    // rank-one updates to the sketched solution
    void update( Eigen::VectorXd& v, double eta);

    // reconstruct solution from sketches
    std::pair<Eigen::MatrixXd, Eigen::MatrixXd> reconstruct();

  private:
    Eigen::MatrixXd* _Omega = nullptr;
    Eigen::MatrixXd* _S = nullptr;

    int _n;
    int _R;
};


#endif
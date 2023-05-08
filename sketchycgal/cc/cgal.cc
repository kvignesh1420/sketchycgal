#include "third_party/eigen3/Eigen/SparseCore"
#include "third_party/eigen3/Eigen/SparseCholesky"
#include "third_party/eigen3/Eigen/Core"
#include "third_party/eigen3/Eigen/Cholesky"
#include "third_party/eigen3/Eigen/Eigenvalues"
#include "third_party/eigen3/Eigen/SVD"
#include <iostream>
#include <cmath>
#include <functional>
#include <random>
#include <utility>
#include <algorithm>
#include <chrono>
#include "cgal.h"
#include "reader.h"
#include "tracer.h"
#include "nystrom.h"

#define EPS 2.2e-16

// handle constructor
SketchyCGAL::SketchyCGAL(){}

/*================================== Private methods ======================================*/

Eigen::SparseMatrix<double, Eigen::RowMajor>* SketchyCGAL::computeLaplacian(Eigen::SparseMatrix<double>* A){
    Eigen::SparseMatrix<double>* D = new Eigen::SparseMatrix<double>(A->rows(), A->cols());
    Eigen::SparseMatrix<double, Eigen::RowMajor>* L = new Eigen::SparseMatrix<double, Eigen::RowMajor>(A->rows(), A->cols());
    Eigen::VectorXd v = Eigen::VectorXd::Ones(A->cols());
    Eigen::VectorXd res = *A * v;
    D->reserve(Eigen::VectorXi::Constant(A->cols(), 1));
    for(int i = 0; i< res.size(); i++){
        // std::cout << res[i] << std::endl;
        if (abs(res[i]) > 1e-10) D->insert(i,i) = res[i];
    }
    D->makeCompressed();
    *L = *D - *A;

    // Validate symmetry of the laplacian matrix
    bool is_symmetric = L->isApprox(L->transpose());

    if (!is_symmetric) {
        std::cerr << "The laplacian matrix is not symmetric.\n" << std::endl;
        exit(1);
    }
    delete D;
    return L;
}

Eigen::VectorXd SketchyCGAL::primitive1(const Eigen::VectorXd& x){
    // std::cout<<"primitive1" <<std::endl;
    return _SCALE_A * (*_C)*(x);
}

Eigen::VectorXd SketchyCGAL::primitive2(const Eigen::VectorXd& y,const Eigen::VectorXd& x){
    // std::cout<<"primitive2" <<std::endl;
    return _SCALE_A * (y.array().cwiseProduct(x.array()));
}

Eigen::VectorXd SketchyCGAL::primitive3(const Eigen::VectorXd& x){
    // std::cout<<"primitive3" <<std::endl;
    return _SCALE_A * x.array().cwiseProduct(x.array());
}

std::pair<Eigen::VectorXd, double> SketchyCGAL::ApproxMinEvecLanczosSE(Eigen::VectorXd& vt, int n, int q){

    // choose the minimum of iterations vs dimension
    // std::cout << "q " << q << std::endl;
    q = std::min(q, n-1);
    // std::cout << "q " << q << std::endl;
    double* aleph = new double[q]();
    double* beth = new double[q]();

    // std::cout << "initialized aleph and beth in ApproxMinEvecLanczosSE" << std::endl;

    std::default_random_engine _generator;
    std::normal_distribution<double> _distribution(0.0, 1.0);
    Eigen::VectorXd v(n);
    for (int i = 0; i < n; ++i) {
        v(i) = _distribution(_generator);
    }
    // v.setRandom();
    // std::cout <<" v vector " << v << std::endl;
    // std::cout << " v size = "<< v.size() << " v norm = " << v.matrix().norm() << std::endl;
    v /= v.norm();
    Eigen::VectorXd vi = v;
    Eigen::VectorXd vip, vim;

    // define a lambda function M
    auto M = [this](Eigen::VectorXd& a, Eigen::VectorXd& b) -> Eigen::VectorXd {
        return this->primitive1(a) + this->primitive2(b, a);
    };

    // i stores the number of completed iterations
    int i = 0;
    // std::cout << "starting "<< i << " iterations in ApproxMinEvecLanczosSE" << std::endl;
    for(i = 0; i < q; i++){
        vip = M( vi, vt );
        // std::cout << " vip size = "<< vip.size() << " vip norm = " << vip.matrix().norm() << std::endl;
        aleph[i] = vi.dot(vip);

        if(i==0){
            vip.noalias() -= aleph[i]*vi;
        }else{
            vip.noalias() -= (aleph[i]*vi + beth[i-1] * vim);
        }
        beth[i] = vip.norm();

        // if ( abs(beth[i]) < sqrt(n)*EPS ) break;

        vip /= beth[i];
        vim = vi;
        vi = vip;
    }

    // std::cout << "completed "<< i << " iterations in ApproxMinEvecLanczosSE" << std::endl;
    // prepare the tri-diagonal matrix

    std::vector<Eigen::Triplet<double> > triplets;
    for(int j = 0; j < q; j++){
        triplets.push_back(Eigen::Triplet<double>(j, j, aleph[j]));
        if(j+1 < q){
            triplets.push_back(Eigen::Triplet<double>(j, j+1, beth[j]));
            triplets.push_back(Eigen::Triplet<double>(j+1, j, beth[j]));
        }
    }
    Eigen::SparseMatrix<double> B = Eigen::SparseMatrix<double>(q, q);
    B.setFromTriplets(triplets.begin(), triplets.end());

    Eigen::MatrixXd U;
    Eigen::VectorXd D;
    std::pair<Eigen::MatrixXd, Eigen::VectorXd> eig_info = this->cgal_eig(B);
    U = eig_info.first;
    D = eig_info.second;
    // std::cout << "U : " << U << std::endl;
    // std::cout << "D : " << D << std::endl;

    int min_D_index;
    double min_D_value = D.minCoeff(&min_D_index);
    Eigen::VectorXd Uind = U.col(min_D_index);
    // std::cout << "min_D_index: "<< min_D_index << " min_D_value: " << min_D_value << std::endl;

    /* start second loop */

    // reset values
    for (int idx = 0; idx < q; idx++) {
        aleph[idx] = 0.0;
        beth[idx] = 0.0;
    }

    vi = v;
    v = Eigen::VectorXd::Zero(n);

    for(i = 1; i < Uind.size(); i++){
        v.noalias() += vi * Uind(i);
        vip = M( vi, vt);
        aleph[i] = vi.dot(vip);
        if(i==0){
            vip.noalias() -= aleph[i]*vi;
        }else{
            vip.noalias() -= (aleph[i]*vi + beth[i-1] * vim);
        }
        beth[i] = vip.norm();

        vip /= beth[i];
        vim = vi;
        vi = vip;
    }

    double norm_v = v.norm();
    min_D_value *= norm_v;
    v /= norm_v;

    return std::pair<Eigen::VectorXd, double>(v, min_D_value);

}

std::pair<Eigen::MatrixXd, Eigen::VectorXd> SketchyCGAL::cgal_eig(const Eigen::SparseMatrix<double>& X)
{
    try {
        Eigen::MatrixXd U;
        Eigen::VectorXd D;
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es;
        es.compute(X.toDense());
        // Eigen::EigenSolver<Eigen::MatrixXd> es(X.toDense());
        U = es.eigenvectors().real();
        D = es.eigenvalues().real();
        return std::pair<Eigen::MatrixXd, Eigen::VectorXd>(U, D);
    } catch (const std::exception& e) {
        std::cerr << "eig did not work. Using the svd based replacement instead." << std::endl;
        exit(1);
    }
}

double SketchyCGAL::getCutValue(Eigen::MatrixXd& U){
    int ncols = U.cols();
    double cut_val = 0.0;
    for(int i=0; i<ncols; i++){
        Eigen::VectorXd sign_evec = U.col(i).unaryExpr([](double val) {
            return (val >= 0.0) ? 1.0 : -1.0;
        });

        // std::cout << "sign_evec dims: " << sign_evec.rows() << " " << sign_evec.cols() << std::endl;
        double rank_val = -( sign_evec.dot( (*_C) * sign_evec  ) );

        // std::cout << "sign_evec norm: " << sign_evec.norm()<< std::endl;
        // std::cout << "rank _val " << rank_val << std::endl;
        cut_val = std::max(cut_val, rank_val);
    }
    return cut_val;
}

/*======================================================================================*/

Eigen::SparseMatrix<double> SketchyCGAL::getAdjacencyMatrix(){
    return *_A;
}

Eigen::SparseMatrix<double, Eigen::RowMajor> SketchyCGAL::getLaplacian(){
    return *_C;
}

// setup method to prepare the necessary matrices
void SketchyCGAL::setup(char* filepath){

    // Create Reader object and prepare the adjacency matrix
    Reader* r = new MMReader(filepath);
    _A = r->Read();
    std::cout << "Parsing " << filepath << std::endl;
    std::cout << "nnz(A) = " << _A->nonZeros() << std::endl;

    // prepare laplacian
    _C = this->computeLaplacian(_A);
    *_C = -(1.0/4.0)*(*_C);
    std::cout << "nnz(C) = " << _C->nonZeros() << std::endl;

    // compute scaling factors
    int n = _C->rows();
    _SCALE_X = 1.0/n;
    _SCALE_C = 1.0/_C->norm();
    std::cout << "SCALE_X : " << _SCALE_X << " _SCALE_C : " << _SCALE_C << std::endl;

    // setup the nystrom sketch object
    _nysketch = new NystromSketch(/*n=*/n, /*R=*/_R);
}


void SketchyCGAL::run(){

    int n = _C->rows();
    int a = n;

    Eigen::VectorXd b_orig = Eigen::VectorXd::Ones(n);
    Eigen::VectorXd b = Eigen::VectorXd::Ones(n);

    // rescale w.r.t A
    b *= _SCALE_A;
    _RESCALE_FEAS /= _SCALE_A;
    // rescale w.r.t X
    b *= _SCALE_X;
    a *= _SCALE_X;
    _RESCALE_OBJ /= _SCALE_X;
    _RESCALE_FEAS /= _SCALE_X;
    // rescale w.r.t C
    _RESCALE_OBJ /= _SCALE_C;

    // initialize dual
    Eigen::VectorXd z = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd y0 = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd y = y0;
    Eigen::VectorXd vt;
    double pobj = 0.0;

    double TRACE = 0.0;
    double beta, eta, FeasOrg, FeasCond, ObjCond;

    bool stop_tolerance = false;
    int t;
    for (t = 1; t < _MAX_ITERS + 1; t++){

        beta = _beta0 * sqrt( t + 1);
        eta = 2.0/(t + 1);

        vt = y + beta * (z - b);

        int q = ceil( pow(t, 0.25)*log(n) );
        std::pair<Eigen::VectorXd, double> ev_info = this->ApproxMinEvecLanczosSE(vt, n, q);

        Eigen::VectorXd smallest_ev = ev_info.first;
        smallest_ev = sqrt(a)*smallest_ev;

        if (!stop_tolerance){

            FeasOrg = ((z - b)*_RESCALE_FEAS).norm();
            double clipped_norm = std::max(b_orig.norm(), 1.0);
            FeasCond = FeasOrg / clipped_norm;
            Eigen::VectorXd Ahk = this->primitive3(smallest_ev);
            double ObjCond1 = pobj - this->primitive1(smallest_ev).dot( smallest_ev );
            double ObjCond2 = y.dot( b - Ahk ) + beta * ( z - b ).dot( z - Ahk );
            double ObjCond3 = -0.5*beta* pow((z-b).norm(), 2);
            double ObjCond4 = std::max( double(abs(pobj*_RESCALE_OBJ)) , 1.0);
            ObjCond = (ObjCond1 + ObjCond2 + ObjCond3)*_RESCALE_OBJ / ObjCond4;
            // if(t % 100 == 0){
                // std::cout << " Iter " << t << " FeasOrg = " << FeasOrg << " FeasCond = " << FeasCond << " ObjCond = " << ObjCond << std::endl;
                // std::cout << " Iter " << t << std::endl;
            // }

            if (FeasCond <= _TOLERANCE && ObjCond <= _TOLERANCE){
                break;
            }
        }

        Eigen::VectorXd zEvec = this->primitive3(smallest_ev);
        // std::cout << "norm of zEvec: " << zEvec.matrix().norm() << std::endl;
        Eigen::VectorXd z_new = (1-eta)*z + eta*zEvec;
        // std::cout << "norm of z_new: " << z_new.matrix().norm() << std::endl;
        // std::cout << "norm of z_new-z: " << (z_new-z).matrix().norm() << std::endl;
        z = z_new;
        TRACE = (1-eta)*TRACE + eta*a;
        double ObjEvec = smallest_ev.matrix().dot(  this->primitive1(smallest_ev).matrix() );
        pobj = (1-eta)*pobj + eta*ObjEvec;
        _nysketch->update(/*v=*/smallest_ev, /*eta=*/eta);

        double beta_p = _beta0 * sqrt(t+2);
        Eigen::VectorXd dual_update = z - b;

        double NORM_A = 1;
        double sigma = std::min(_beta0, 4 * beta_p * eta*eta * a*a * NORM_A*NORM_A / pow(dual_update.matrix().norm(),2) );

        // update dual
        Eigen::VectorXd yt1 = y + sigma * dual_update;
        if ( (yt1 - y0).matrix().norm() <= _K ){
            y = yt1;
        }

    }

    std::cout << " Iter " << t << " FeasOrg = " << FeasOrg << " FeasCond = " << FeasCond << " ObjCond = " << ObjCond << std::endl;

    std::pair<Eigen::MatrixXd, Eigen::MatrixXd> sketch_res = _nysketch->reconstruct();
    Eigen::MatrixXd U = sketch_res.first;
    Eigen::MatrixXd Delt = sketch_res.second;

    // std::cout << "Dimensions of U : " << U.rows() << " " << U.cols() << std::endl;
    // std::cout << "Dimensions of Delt : " << Delt.rows() << " " << Delt.cols() << std::endl;
    // std::cout << "Trace value is " << TRACE << std::endl;
    Delt.noalias() += (( TRACE - Delt.trace() )/_R)*Eigen::MatrixXd::Identity(Delt.rows(), Delt.cols());
    U = U*Delt.array().sqrt().matrix();
    U = (U.array()/sqrt(_SCALE_X)).matrix();

    std::cout << "Dimensions of U : " << U.rows() << " " << U.cols() << std::endl;
    std::cout << "Dimensions of Delt : " << Delt.rows() << " " << Delt.cols() << std::endl;
    double max_cut = this->getCutValue(U);
    std::cout << "max-cut value : " << max_cut << std::endl;


}

#ifndef SKETCHYCGAL_CC_READER_H
#define SKETCHYCGAL_CC_READER_H

#include "third_party/eigen3/Eigen/Core"
#include "third_party/eigen3/Eigen/SparseCore"

class Reader{
  public:
    // read data and return the Eigen matrix representation of the data
    virtual Eigen::SparseMatrix<double>* Read() = 0;
};

class MMReader : public Reader{
  public:
    // A reader to parse MM (MatrixMarket) files and prepare the adjacency matrix
    // with eigen sparse matrix data structure.
    MMReader(char* filename);
    Eigen::SparseMatrix<double>* Read() override;
  private:
    char* _filename;
    int _nrows, _ncols, _nnz;
};


#endif
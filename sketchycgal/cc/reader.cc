//
// reader code adapted from https://math.nist.gov/MatrixMarket/mmio/c/example_read.c
//
#include "third_party/eigen3/Eigen/Core"
#include "third_party/eigen3/Eigen/SparseCore"
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include "reader.h"

MMReader::MMReader(char* filename){
    _filename = filename;
}

Eigen::SparseMatrix<double>* MMReader::Read(){

    std::ifstream _file;

    _file.open(_filename);
    if (!_file.is_open()){
        std::cerr << "Error: The input file couldn't be opened" << std::endl;
        exit(1);
    }

    // Read and ignore all comment lines
    std::string line;
    while (std::getline(_file, line)) {
        if (line[0] != '%') {
            break;
        }
    }

    // Parse the header line
    std::istringstream iss(line);
    /* These fields are not present for the G67 dataset */
    // std::string objectClass, matrixType, storageScheme;
    iss >> _nrows >> _ncols >> _nnz;

    std::cout << "rows " << _nrows << " cols " << _ncols << " nnz " << _nnz << std::endl;
 
    // Read the matrix entries
    std::vector<Eigen::Triplet<double> > triplets;
    for (int i = 0; i < _nnz; i++) {
        std::getline(_file, line);
        if (line[0] == '%') {
            i--; // Ignore comment lines
            continue;
        }

        int row, col;
        double value;
        std::istringstream iss(line);
        iss >> row >> col >> value;

        // Store the matrix entries. Note that we add entries for the lower/upper half
        // as well to avoid transpose operations for sparse matrix.
        // handle so called "zero" in double
        if(abs(value) > 1e-10){
            // std::cout << "Pushing: "<< row-1 << " " << col-1 << " " << value << std::endl;
            triplets.push_back(Eigen::Triplet<double>(row-1, col-1, value));
            triplets.push_back(Eigen::Triplet<double>(col-1, row-1, value));
        }
        else{
            // std::cout << "Pushing: "<< row-1 << " " << col-1 << " 1.0" << std::endl;
            triplets.push_back(Eigen::Triplet<double>(row-1, col-1, 1.0));
            triplets.push_back(Eigen::Triplet<double>(col-1, row-1, 1.0));
        }
    }

    Eigen::SparseMatrix<double>* Adj = new Eigen::SparseMatrix<double>(_nrows, _ncols);
    Adj->setFromTriplets(triplets.begin(), triplets.end());

    // Validate symmetry of the adjacency matrix
    bool is_symmetric = Adj->isApprox(Adj->transpose());
    if (!is_symmetric) {
        std::cerr << "The adjacency matrix is not symmetric. Please check the input file.\n" << std::endl;
        exit(1);
    }

    return Adj;
}

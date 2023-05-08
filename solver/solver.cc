#include <iostream>
#include <chrono>
#include "third_party/eigen3/Eigen/SparseCore"
#include "third_party/eigen3/Eigen/SparseCholesky"
#include "sketchycgal/cc/reader.h"
#include "sketchycgal/cc/cgal.h"

#if defined(_OPENMP)
#include <omp.h>
#define EIGEN_DONT_PARALLELIZE
#endif

int main(int argc, char *argv[]){

    char* filepath = argv[1];
    #ifdef EIGEN_DONT_PARALLELIZE
    std::cout << "Default Eigen parallelization disabled, OpenMP available " << std::endl;
    #else
    std::cout << "Default Eigen parallelization enabled, OpenMP unavailable " << std::endl;
    #endif
    SketchyCGAL Runner = SketchyCGAL();
    Runner.setup(filepath, 10, 10, 0.1);
    std::cout << "Starting the SketchyCGAL loop " << std::endl;
    auto start_time = std::chrono::high_resolution_clock::now();
    Runner.run();
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    double seconds = static_cast<double>(duration.count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den;

    std::cout << "Elapsed time: " << seconds << " seconds." << std::endl;

    return 0;
}

#include <iostream>
#include <chrono>
#include "third_party/eigen3/Eigen/SparseCore"
#include "third_party/eigen3/Eigen/SparseCholesky"
#include "sketchycgal/cc/reader.h"
#include "sketchycgal/cc/cgal.h"


int main(int argc, char *argv[]){

    char* filepath = argv[1];

    SketchyCGAL Runner = SketchyCGAL();
    Runner.setup(filepath);
    auto start_time = std::chrono::high_resolution_clock::now();
    Runner.run();
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    double seconds = static_cast<double>(duration.count()) * std::chrono::microseconds::period::num / std::chrono::microseconds::period::den;

    std::cout << "Elapsed time: " << seconds << " seconds." << std::endl;

    return 0;
}

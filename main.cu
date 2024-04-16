//
// Created by Yanliang Li on 4/9/24.
//
#include "FlyingEdges.cuh"
//#include "utils/rawImageReader.h"
#include "utils/testDataReader.h"
#include "utils/Image3D.h"
#include "config/config.h"

void checkCudaErrors(cudaError_t result) {
    if (result != cudaSuccess) {
        std::cerr << "CUDA Runtime Error: " << cudaGetErrorString(result) << std::endl;
        assert(result == cudaSuccess);
    }
}

int main() {
    scalar_t isoval = 0.5;
    std::array<scalar_t, 3> spacing = {1, 1, 1};
    std::array<scalar_t, 3> zeroPos = {0, 0, 0};
    std::array<size_t, 3> dimensions = {100, 100, 100};
//    std::string filePath = "data/FullHead.raw";
    std::string filePath = "data/sphere_volume.bin";
    std::vector<scalar_t> pixels = readTestData(filePath, dimensions[0], dimensions[1], dimensions[2]);
    Image3D image(pixels, spacing, zeroPos, dimensions);
    FlyingEdges algo(image, isoval);

    cudaEvent_t start, stop;
    float milliseconds = 0;
    checkCudaErrors(cudaEventCreate(&start));
    checkCudaErrors(cudaEventCreate(&stop));
    checkCudaErrors(cudaEventRecord(start));
    algo.pass1();
    algo.pass2();
    algo.pass3();
    algo.pass4();
    checkCudaErrors(cudaEventRecord(stop));
    checkCudaErrors(cudaEventSynchronize(stop));
    checkCudaErrors(cudaEventElapsedTime(&milliseconds, start, stop));
    std::cout << "GPU FlyingEdges Algorithm running time: " << milliseconds << " ms" << std::endl;

    algo.moveOutput();
    algo.writeObj();

    checkCudaErrors(cudaEventDestroy(start));
    checkCudaErrors(cudaEventDestroy(stop));

    return 0;
}
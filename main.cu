//
// Created by Yanliang Li on 4/9/24.
//
#include "FlyingEdges.cuh"
#include "utils/rawImageReader.h"
#include "utils/Image3D.h"
#include "config/config.h"

int main() {
    scalar_t isoval = 20.0;
    std::array<scalar_t, 3> spacing = {1, 1, 1};
    std::array<scalar_t, 3> zeroPos = {0, 0, 0};
    std::array<size_t, 3> dimensions = {100, 100, 100};
//    std::string filePath = "data/FullHead.raw";
    std::string filePath = "data/sphere_volume.bin";
    std::vector<scalar_t> pixels = readRawFile(filePath, dimensions[0], dimensions[1], dimensions[2]);
    Image3D image(pixels, spacing, zeroPos, dimensions);
    FlyingEdges algo(image, isoval);

    algo.pass1();
    algo.pass2();
    algo.pass3();
    algo.pass4();
    algo.moveOutput();
    algo.writeObj();

    return 0;
}
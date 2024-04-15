//
// Created by Yanliang Li on 4/9/24.
//
#include "FlyingEdges.cuh"
#include "utils/rawImageReader.h"
#include "utils/Image3D.h"
#include "config/config.h"

int main() {
    scalar_t isoval = 20.0;
    std::array<scalar_t, 3> spacing = {0.9375, 0.9375, 1.5};
    std::array<scalar_t, 3> zeroPos = {0, 0, 0};
    std::array<size_t, 3> dimensions = {256, 256, 94};
    std::string filePath = "data/FullHead.raw";
    std::vector<unsigned short> pixels = readRawFile(filePath, dimensions[0], dimensions[1], dimensions[2]);
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
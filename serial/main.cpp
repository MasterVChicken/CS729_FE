//
// Created by Yanliang Li on 4/16/24.
//

/*
 * main.cpp
 *
 *  Created on: Feb 17, 2017
 *      Author: dbourge
 *
 * miniIsosurface is distributed under the OSI-approved BSD 3-clause License.
 * See LICENSE.txt for details.
 *
 * Copyright (c) 2017
 * National Technology & Engineering Solutions of Sandia, LLC (NTESS). Under
 * the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains
 * certain rights in this software.
 */
#include <iostream>
#include <string.h>
#include <cstdlib>
#include "util/Image3D.h"
#include "FlyingEdgesSerial.h"
#include "../utils/testDataReader.h"

int main(int argc, char *argv[]) {
    scalar_t isoval = 0.5;
    std::array<scalar_t, 3> spacing = {1, 1, 1};
    std::array<scalar_t, 3> zeroPos = {0, 0, 0};
    std::array<size_t, 3> dimensions = {100, 100, 100};
    std::string filePath = "data/sphere_volume.bin";
    std::vector<scalar_t> pixels = readTestData(filePath, dimensions[0], dimensions[1], dimensions[2]);
    Image3D image(pixels, spacing, zeroPos, dimensions);
    FlyingEdgesSerial algo(image, isoval);

    algo.pass1();
    algo.pass2();
    algo.pass3();
    algo.pass4();
    algo.writeObj();
}


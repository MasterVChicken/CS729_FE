//
// Created by Yanliang Li on 4/9/24.
//

#ifndef CS729_FE_CONFIG_H
#define CS729_FE_CONFIG_H

#pragma once

#include <array>
#include <vector>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <cstdint>
#include <cstdlib>
#include <cuda_runtime.h>

using std::size_t;

using scalar_t = float;

using uchar = unsigned char;

using ushort = unsigned short;

using cube_t = std::array<std::array<scalar_t, 3>, 8>;
using scalarCube_t = std::array<scalar_t, 8>;

// TODO: we might need to adjust the BLOCK WIDTH due to FullHead.raw is only about 256*256*94, so max dimension should be 256

//#define FE_BLOCK_WIDTH 512
//#define FE_BLOCK_WIDTH_PLUS_ONE 513

#define FE_BLOCK_WIDTH 128
#define FE_BLOCK_WIDTH_PLUS_ONE 129

// FE_BLOCK_WIDTH = FE_BLOCK_WIDTH_Y * FE_BLOCK_WIDTH_Z
#define FE_BLOCK_WIDTH_Y 8
#define FE_BLOCK_WIDTH_Z 16

#endif //CS729_FE_CONFIG_H

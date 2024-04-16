//
// Created by Yanliang Li on 4/16/24.
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
#include <cassert>

using std::size_t;

using scalar_t = float;

using uchar = unsigned char;

using ushort = unsigned short;

using cube_t = std::array<std::array<scalar_t, 3>, 8>;
using scalarCube_t = std::array<scalar_t, 8>;

#endif //CS729_FE_CONFIG_H

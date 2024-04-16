//
// Created by Yanliang Li on 4/15/24.
//

#ifndef CS729_FE_TESTDATAREADER_H
#define CS729_FE_TESTDATAREADER_H

#include "../config/config.h"

// read data from generated data for test usage
// returns vector of float
std::vector<scalar_t> readTestData(const std::string &filePath, int width, int height, int depth) {
    std::ifstream file(filePath, std::ios::binary);
    if (!file) {
        std::cerr << "Cannot open file.\n";
        return pixels;
    }

    int total_elements = width * height * depth;
    std::vector <scalar_t> pixels(total_elements);

    file.read(reinterpret_cast<char *>(pixels.data()), total_elements * sizeof(float));
    file.close();

    return pixels;
}

#endif //CS729_FE_TESTDATAREADER_H

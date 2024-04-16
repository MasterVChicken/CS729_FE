//
// Created by Yanliang Li on 4/9/24.
//

#ifndef CS729_FE_RAWIMAGEREADER_H
#define CS729_FE_RAWIMAGEREADER_H

#include  <vector>
#include "../config/config.h"

// read meta data from raw file
// return vector of unsigned short
std::vector<unsigned short> readRawFile(const std::string &filePath, int width, int height, int depth) {
    std::ifstream file(filePath, std::ios::binary);
    std::vector<unsigned short> pixels;

    if (!file) {
        std::cerr << "Cannot open file.\n";
        return pixels;
    }

    size_t totalPixels = width * height * depth;
    pixels.resize(totalPixels);

    file.read(reinterpret_cast<char *>(pixels.data()), totalPixels * sizeof(unsigned short));

    if (!file) {
        std::cerr << "Error reading file.\n";
        return std::vector<unsigned short>();
    }

    return pixels;
}

#endif //CS729_FE_RAWIMAGEREADER_H

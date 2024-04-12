//
// Created by Yanliang Li on 4/9/24.
//

#ifndef CS729_FE_RAWIMAGEREADER_H
#define CS729_FE_RAWIMAGEREADER_H

#include<vector>
#include "config/config.h"

// read meta data from raw file
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

// To make data suitable for iso-contouring application, we need to transform it to some scalar format
// transform data from RGB 565 format to grayscale
void RGB565toGreyscale(scalar_t *out, std::vector<unsigned short> pixels) {
    for (int i = 0; i < pixels.size(); i++) {
        // Get channel values for R, G and B
        int r = (pixels[i] >> 11) & 0x1F;
        int g = (pixels[i] >> 5) & 0x3F;
        int b = pixels[i] & 0x1F;
        out[i] = 0.299 * (scalar_t) r + 0.587 * g + 0.114 * b;
    }
}

#endif //CS729_FE_RAWIMAGEREADER_H

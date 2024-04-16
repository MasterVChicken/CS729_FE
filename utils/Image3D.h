//
// Created by Yanliang Li on 4/12/24.
//

#ifndef CS729_FE_IMAGE3D_H
#define CS729_FE_IMAGE3D_H

#include "../config/config.h"

class Image3D {
public:
    Image3D(const std::vector <scalar_t> &metadata,
            const std::array<scalar_t, 3> &spacing,
            const std::array<scalar_t, 3> &zeroPos,
            const std::array<size_t, 3> &dimensions)
            : metadata(metadata), spacing(spacing), zeroPos(zeroPos),
              nx(dimensions[0]), ny(dimensions[1]), nz(dimensions[2]) {
        if (nx == 0 || ny == 0 || nz == 0) {
            throw std::runtime_error("Dimensions cannot be zero.");
        }

        if (metadata.size() != nx * ny * nz) {
            throw std::runtime_error("Metadata size does not match dimensions.");
        }

        data = new scalar_t[nx * ny * nz];  // Allocating memory
        for (int i = 0; i < metadata.size(); i++) {
            data[i] = metadata[i];
        }
//        data = metadata.data();
        // If it is raw data with RGB 565 format
        // Transform data from RGB 565 format to grayscale
//        for (size_t i = 0; i < metadata.size(); i++) {
//            int r = (metadata[i] >> 11) & 0x1F;
//            int g = (metadata[i] >> 5) & 0x3F;
//            int b = metadata[i] & 0x1F;
//            data[i] = 0.299 * r + 0.587 * g + 0.114 * b;
//        }
    }

    int getX() const {
        return nx;
    }

    int getY() const {
        return ny;
    }

    int getZ() const {
        return nz;
    }

    scalar_t *getDataPointer() {
        return data;
    }

    std::array<scalar_t, 3> getZeroPos() const {
        return zeroPos;
    }

    std::array<scalar_t, 3> getSpacing() const {
        return spacing;
    }

    ~Image3D() {
        if (!data) {
            delete[] data;  // Properly delete allocated memory
        }
    }


private:
//    std::vector<ushort> metadata;       // A vector containing scalar values
    // along three-dimensional space.

    std::vector <scalar_t> metadata;

    // Grayscale data
    scalar_t *data;

    std::array<scalar_t, 3> spacing;    // The distance between two points in
    // the mesh.

    std::array<scalar_t, 3> zeroPos;    // The position at index (0, 0, 0).

    size_t nx;         //
    size_t ny;         // The dimensions
    size_t nz;         //
};

#endif //CS729_FE_IMAGE3D_H

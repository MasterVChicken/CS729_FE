//
// Created by Yanliang Li on 4/16/24.
//

#ifndef CS729_FE_FLYINGEDGESSERIAL_H
#define CS729_FE_FLYINGEDGESSERIAL_H

#include "util/Image3D.h"
#include "config/config.h"

struct FlyingEdgesSerial {
    FlyingEdgesSerial(Image3D const &image, scalar_t const &isoval)
            : image(image),
              isoval(isoval),
              nx(image.xdimension()),
              ny(image.ydimension()),
              nz(image.zdimension()),
              gridEdges(ny * nz),
              triCounter((ny - 1) * (nz - 1)),
              edgeCases((nx - 1) * ny * nz),
              cubeCases((nx - 1) * (ny - 1) * (nz - 1)) {}

    void pass1();

    void pass2();

    void pass3();

    void pass4();

    void writeObj();

private:
    struct gridEdge {
        gridEdge()
                : xl(0),
                  xr(0),
                  xstart(0),
                  ystart(0),
                  zstart(0) {}

        // trim values
        // set on pass 1
        size_t xl;
        size_t xr;

        // modified on pass 2
        // set on pass 3
        size_t xstart;
        size_t ystart;
        size_t zstart;
    };

private:
    Image3D const &image;
    scalar_t const isoval;

    size_t const nx; //
    size_t const ny; // for indexing
    size_t const nz; //

    std::vector <gridEdge> gridEdges; // size of ny*nz
    std::vector <size_t> triCounter;  // size of (ny-1)*(nz-1)

    std::vector <uchar> edgeCases;    // size (nx-1)*ny*nz
    std::vector <uchar> cubeCases;    // size (nx-1)*(ny-1)*(nz-1)

    std::vector <std::array<scalar_t, 3>> points;  //
    std::vector <std::array<scalar_t, 3>> normals; // The output
    std::vector <std::array<size_t, 3>> tris;     //

private:
    bool isCutEdge(size_t const &i, size_t const &j, size_t const &k) const;

    inline uchar
    calcCaseEdge(bool const &prevEdge, bool const &currEdge) const;

    inline uchar
    calcCubeCase(uchar const &ec0, uchar const &ec1,
                 uchar const &ec2, uchar const &ec3) const;

    inline void calcTrimValues(
            size_t &xl, size_t &xr, size_t const &j, size_t const &k) const;

    inline std::array<scalar_t, 3>
    interpolateOnCube(
            cube_t const &pts,
            scalarCube_t const &isovals,
            uchar const &edge) const;

    inline std::array<scalar_t, 3>
    interpolate(
            std::array<scalar_t, 3> const &a,
            std::array<scalar_t, 3> const &b,
            scalar_t const &weight) const;
};


#endif //CS729_FE_FLYINGEDGESSERIAL_H

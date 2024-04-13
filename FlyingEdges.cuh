//
// Created by Yanliang Li on 4/9/24.
//

#ifndef CS729_FE_FLYINGEDGES_CUH
#define CS729_FE_FLYINGEDGES_CUH

#include <iostream>
#include <string>
#include "/config/config.h"
#include "/utils/Image3D.h"
#include "/utils/marchingCubesTables.h"


class FlyingEdges {
    FlyingEdges(Image3D image, scalar_t const &isoval)
            : isoval(isoval),
              nx(image.getX()),
              ny(image.getY()),
              nz(image.getZ()),
              deallocated(false) {

        // Allocate device memory
        cudaError_t cuErr;

        cuErr = cudaMalloc(&pointValues, nx * ny * nz * sizeof(scalar_t));
        if (cuErr != cudaSuccess) {
            std::cout << "Error occured when allocating memory for pointValues" << std::endl;
            throw;
        }

        cuErr = cudaMalloc(&zero_pos, 3 * sizeof(scalar_t));
        if (cuErr != cudaSuccess) {
            std::cout << "Error occured when allocating memory for zero_pos" << std::endl;
            throw;
        }

        cuErr = cudaMalloc(&spacing, 3 * sizeof(scalar_t));
        if (cuErr != cudaSuccess) {
            std::cout << "Error occured when allocating memory for spacing" << std::endl;
            throw;
        }

        cuErr = cudaMalloc(&gridEdges, ny * nz * sizeof(gridEdge));
        if (cuErr != cudaSuccess) {
            std::cout << "Error occured when allocating memory for gridEdges" << std::endl;
            throw;
        }

        cuErr = cudaMalloc(&triCounter, (ny - 1) * (nz - 1) * sizeof(int));
        if (cuErr != cudaSuccess) {
            std::cout << "Error occured when allocating memory for triCounter" << std::endl;
            throw;
        }

        cuErr = cudaMalloc(&edgeCases, (nx - 1) * ny * nz * sizeof(uchar));
        if (cuErr != cudaSuccess) {
            std::cout << "Error occured when allocating memory for edgeCases" << std::endl;
            throw;
        }

        cuErr = cudaMalloc(&cubeCases, (nx - 1) * (ny - 1) * (nz - 1) * sizeof(uchar));
        if (cuErr != cudaSuccess) {
            std::cout << "Error occured when allocating memory for cubeCases" << std::endl;
            throw;
        }

        // Move data from host to device
        cudaMemcpy(
                pointValues,
                image.getDataPointer,
                nx * ny * nz * sizeof(scalar_t),
                cudaMemcpyHostToDevice);

        scalar_t *tmp = (scalar_t *) malloc(3 * sizeof(scalar_t));
        auto imageZp = image.getZeroPos();
        tmp[0] = imageZp[0];
        tmp[1] = imageZp[1];
        tmp[2] = imageZp[2];

        cudaMemcpy(
                zeroPos,
                tmp,
                3 * sizeof(scalar_t),
                cudaMemcpyHostToDevice);

        scalar_t *imageSp = image.getSpacing();
        tmp[0] = imageSp[0];
        tmp[1] = imageSp[1];
        tmp[2] = imageSp[2];

        cudaMemcpy(
                spacing,
                tmp,
                3 * sizeof(scalar_t),
                cudaMemcpyHostToDevice);
    }

    ~FlyingEdgesAlgorithm() {
        deallocate();
    }

    void loadMCtables(){

    }

    void deallocate() {
        if (!deallocated) {
            cudaFree(pointValues);
            cudaFree(zeroPos);
            cudaFree(spacing);
            cudaFree(gridEdges);
            cudaFree(triCounter);
            cudaFree(edgeCases);
            cudaFree(cubeCases);
            // Haven't done output
            // Process output later
            //  if(outputAllocated)
            //  {
            //      cudaFree(points);
            //      cudaFree(normals);
            //      cudaFree(tris);
            //  }
        }
        deallocated = true;
    }

    void pass1();

    void pass2();

    void pass3();

    void pass4();


public:
    struct gridEdge {
        gridEdge()
                : xl(0),
                  xr(0),
                  xstart(0),
                  ystart(0),
                  zstart(0) {}

        // trim values
        // set on pass 1
        int xl;
        int xr;

        // modified on pass 2
        // set on pass 3
        int xstart;
        int ystart;
        int zstart;
    };

private:
    // Important meta data to be set during construct
    scalar_t *pointValues;
    scalar_t *zero_pos;
    scalar_t *spacing;
    scalar_t const isoval;

    // true data to be processed during the whole algorithm
    gridEdge *gridEdges; // size of ny*nz
    int *triCounter;     // size of (ny-1)*(nz-1)

    uchar *edgeCases; // size of (nx-1)*ny*nz
    uchar *cubeCases; // size of (nx-1)*(ny-1)*(nz-1)

    // data for ouput
    scalar_t *points;     // the vertex processed with spacing
    scalar_t *normals;    // normal of each vertexes
    int *tris;            // triangle position for output
    int numPoints;        // number of vertexes
    int numTris;          // number of triangles

    // indexing
    int nx;
    int ny;
    int nz;

    // mark for deallocated or not
    bool deallocated;
};


#endif //CS729_FE_FLYINGEDGES_CUH

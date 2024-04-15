//
// Created by Yanliang Li on 4/9/24.
//

#include "FlyingEdges.cuh"
#include "utils/marchingCubesTables.h"
#include "config/config.h"

// TODO: We forget not free some of symbols for outputs. Try to fix bugs in future.

// Pass 1 of the algorithm labels each edge parallel to the x-axis as cut
// or not. In the process, each gridEdge is assigned an xl and xr.
// All edges before xl are uncut and all edges after xr are uncut.
// Subsequent passes of the algorithm don't look outside of [xl, xr).
// A gridEdge E_jk can be thought of as the row of edges parallel to the
// x-axis for some fixed j and k.

__device__
uchar calcCaseEdge(
        bool const &prevEdge,
        bool const &currEdge) {
    // o -- is greater than or equal to
    // case 0: (i-1) o-----o (i) | (_,j,k)
    // case 1: (i-1) x-----o (i) | (_,j+1,k)
    // case 2: (i-1) o-----x (i) | (_,j,k+1)
    // case 3: (i-1) x-----x (i) | (_,j+1,k+1)
    if (prevEdge && currEdge)
        return 0;
    if (!prevEdge && currEdge)
        return 1;
    if (prevEdge && !currEdge)
        return 2;
    else // !prevEdge && !currEdge
        return 3;
}

__global__
void assign_edgeCases(
        scalar_t *pointValues,
        scalar_t isoval,
        int nx, int ny,
        uchar *edgeCases) {
    // Each row has several blocks
    // Each thread is one point

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y;
    int k = blockIdx.z;

    // Why plus 1?
    // Because we need to make sure that we can deal with boundary case
    // To make compare for the last element of current block and the first element of next block
    // But we cannot straightly access across different blocks, so we do pre-load here.
    __shared__ bool isGE[FE_BLOCK_WIDTH_PLUS_ONE];

    if (i < nx)
        isGE[threadIdx.x] = pointValues[k * nx * ny + j * nx + i] >= isoval;

    // Maybe bug fixed?
    if (threadIdx.x == 0 && i < nx - 1 && (j * nx + i + blockDim.x) < nx * ny) {
        // Deal with the first element of next block
        isGE[blockDim.x] = pointValues[k * nx * ny + j * nx + i + blockDim.x] >= isoval;
    }

    __syncthreads();

    // write results back to edgeCases
    if (i < nx - 1) {
        uchar caseEdge = calcCaseEdge(isGE[threadIdx.x], isGE[threadIdx.x + 1]);
        edgeCases[k * (nx - 1) * ny + j * (nx - 1) + i] = caseEdge;
    }
}

__global__
void x_edge_trim(
        int nx, int ny, int nz,                    // input
        uchar *edgeCases,                          // input
        FlyingEdges::gridEdge *gridEdges) // output
{
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    int k = blockIdx.y * blockDim.y + threadIdx.y;

    if (j >= ny || k >= nz)
        return;

    size_t xl = nx;
    size_t xr = 0;

    uchar *curEdgeCases = edgeCases + k * (nx - 1) * ny + j * (nx - 1);

    for (int i = 0; i != nx - 1; ++i) {
        if (curEdgeCases[i] == 1 || curEdgeCases[i] == 2) {
            // get first left trim pos
            if (xl == nx)
                xl = i;
            // get last right trim pos
            xr = i + 1;
        }
    }

    gridEdges[k * ny + j].xl = xl;
    gridEdges[k * ny + j].xr = xr;
}

void FlyingEdges::pass1() {
    int tx = FE_BLOCK_WIDTH;
    uint3 gridDim = make_uint3(((nx - 1) + tx - 1) / tx, ny, nz);
    uint3 blockDim = make_uint3(tx, 1, 1);

    assign_edgeCases<<<gridDim, blockDim>>>(
            pointValues,
            isoval,
            nx, ny,
            edgeCases);

    int ty = FE_BLOCK_WIDTH_Y;
    int tz = FE_BLOCK_WIDTH_Z;
    gridDim = make_uint3((ny + ty - 1) / ty, (nz + tz - 1) / tz, 1);
    blockDim = make_uint3(ty, tz, 1);

    x_edge_trim<<<gridDim, blockDim>>>(
            nx, ny, nz,
            edgeCases,
            gridEdges);

    cudaDeviceSynchronize();

}

// Pass 2 of the algorithm determines the marching cubes case ID of each
// cube. This is determined fully from information obtained in pass1, so
// there is no need to access the input image. Each cube starts at (i,j,k)
// and extends to (i+1, j+1, k+1).
// In addition to determining case ID of each cell, pass 2 counts the
// number of cuts on incident to each gridEdge.

// We can get the xl and xr of a signle R-jk cell array
__device__
void calcTrimValues(
        int &xl, int &xr,
        FlyingEdges::gridEdge const &ge0,
        FlyingEdges::gridEdge const &ge1,
        FlyingEdges::gridEdge const &ge2,
        FlyingEdges::gridEdge const &ge3) {

    xl = min(ge0.xl, min(ge1.xl, min(ge2.xl, ge3.xl)));
    xr = max(ge0.xr, max(ge1.xr, max(ge2.xr, ge3.xr)));

    if (xl > xr)
        xl = xr;
}

// Get case id by vertex status
__device__
uchar calcCubeCase(
        uchar const &ec0, uchar const &ec1,
        uchar const &ec2, uchar const &ec3) {
    // o -- is greater than or equal to
    // case 0: (i-1) o-----o (i) | (_,j,k)
    // case 1: (i-1) x-----o (i) | (_,j+1,k)
    // case 2: (i-1) o-----x (i) | (_,j,k+1)
    // case 3: (i-1) x-----x (i) | (_,j+1,k+1)

    // Get cube cases by each vertex status
    uchar caseId = 0;
    if ((ec0 == 0) || (ec0 == 2)) // 0 | (i,j,k)
        caseId |= 1;
    if ((ec0 == 0) || (ec0 == 1)) // 1 | (i+1,j,k)
        caseId |= 2;
    if ((ec1 == 0) || (ec1 == 1)) // 2 | (i+1,j+1,k)
        caseId |= 4;
    if ((ec1 == 0) || (ec1 == 2)) // 3 | (i,j+1,k)
        caseId |= 8;
    if ((ec2 == 0) || (ec2 == 2)) // 4 | (i,j,k+1)
        caseId |= 16;
    if ((ec2 == 0) || (ec2 == 1)) // 5 | (i+1,j,k+1)
        caseId |= 32;
    if ((ec3 == 0) || (ec3 == 1)) // 6 | (i+1,j+1,k+1)
        caseId |= 64;
    if ((ec3 == 0) || (ec3 == 2)) // 7 | (i,j+1,k+1)
        caseId |= 128;
    return caseId;
}

__global__
void get_cubeCases(
        int nx, int ny, int nz,
        uchar *edgeCases,
        FlyingEdges::gridEdge *gridEdges,
        int *triCounter,
        uchar *cubeCases) {
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    int k = blockIdx.y * blockDim.y + threadIdx.y;

    if (j >= ny - 1 || k >= nz - 1)
        return;

    FlyingEdges::gridEdge &ge0 = gridEdges[k * ny + j];
    FlyingEdges::gridEdge &ge1 = gridEdges[k * ny + j + 1];
    FlyingEdges::gridEdge &ge2 = gridEdges[(k + 1) * ny + j];
    FlyingEdges::gridEdge &ge3 = gridEdges[(k + 1) * ny + j + 1];

    uchar *ec0 = edgeCases + k * ny * (nx - 1) + j * (nx - 1);
    uchar *ec1 = edgeCases + k * ny * (nx - 1) + (j + 1) * (nx - 1);
    uchar *ec2 = edgeCases + (k + 1) * ny * (nx - 1) + j * (nx - 1);
    uchar *ec3 = edgeCases + (k + 1) * ny * (nx - 1) + (j + 1) * (nx - 1);

    int xl, xr;
    calcTrimValues(xl, xr, ge0, ge1, ge2, ge3);

    int triCount = 0;
    // Locate the pointer to Cubes represented by current 4 edges
    uchar *curCubeCases = cubeCases + k * (nx - 1) * (ny - 1) + j * (nx - 1);

    int xstart = 0;
    int ystart = 0;
    int zstart = 0;

    const bool *isCutCur;
    for (int i = xl; i != xr; ++i) {
        uchar caseId = calcCubeCase(ec0[i], ec1[i], ec2[i], ec3[i]);

        curCubeCases[i] = caseId;

        // All in or all out
        if (caseId == 0 || caseId == 255) {
            continue;
        }

        // Sum number of triangles
        triCount += numTris[caseId];
        isCutCur = isCut[caseId];

        // edge0, edge3, edge8 represents cut or not in x, y, z axis
        xstart += isCutCur[0];
        ystart += isCutCur[3];
        zstart += isCutCur[8];
    }

    // Write sum of triangles back
    triCounter[k * (ny - 1) + j] = triCount;

    if (xr == nx - 1) {
        // Process the last one
        ystart += isCutCur[1];
        zstart += isCutCur[9];
    }

    ge0.xstart = xstart;
    ge0.ystart = ystart;
    ge0.zstart = zstart;
}

// Process the last line parallel to y axises
__global__
void getGhostXZ(
        int nx, int ny, int nz,
        uchar *edgeCases,
        FlyingEdges::gridEdge *gridEdges) {
    int k = blockIdx.x * blockDim.x + threadIdx.x;

    if (k >= nz) // This function will deal with gridEdge at (_, ny-1, nz-1)
        return;

    bool isCorner = k == nz - 1;

    int j = ny - 1;

    FlyingEdges::gridEdge &ge0 = gridEdges[k * ny + j];
    // If isCorner, this is just bogus.
    FlyingEdges::gridEdge &ge1 = gridEdges[(1 - isCorner) * (k + 1) * ny + j];

    uchar *ec0 = edgeCases + k * ny * (nx - 1) + j * (nx - 1);
    // If isCorner, this is just bogus
    uchar *ec1 = edgeCases + (1 - isCorner) * (k + 1) * ny * (nx - 1) + j * (nx - 1);

    int xl = min(ge0.xl, nx * isCorner + (1 - isCorner) * ge1.xl);
    int xr = max(ge0.xr, (1 - isCorner) * ge1.xr);

    int xstart = 0;
    int zstart = 0; // TODO don't set initial values in gridEdge Constructor;

    uchar c0;
    uchar c1;

    for (int i = xl; i != xr; ++i) {
        c0 = ec0[i];
        c1 = ec1[i];

        // see if the edges are cut
        xstart += (c0 == 1 || c0 == 2);

        // bogus if isCorner
        zstart += ((c0 == 0 && c1 == 1) || (c0 == 0 && c1 == 3) ||
                   (c0 == 1 && c1 == 2) || (c0 == 2 && c1 == 3));
    }

    if (xr == nx - 1) {
        // bogus if isCorner
        zstart += ((c0 == 0 && c1 == 2) || (c0 == 0 && c1 == 3) ||
                   (c0 == 1 && c1 == 2) || (c0 == 1 && c1 == 3));
    }

    ge0.xstart = xstart;
    ge0.ystart = 0;
    ge0.zstart = zstart * (1 - isCorner);
}

// Process last line parallel to z axises
__global__
void getGhostXY(
        int nx, int ny, int nz,
        uchar *edgeCases,
        FlyingEdges::gridEdge *gridEdges) {
    int j = blockIdx.x * blockDim.x + threadIdx.x;

    if (j >= ny - 1)
        return;

    int k = nz - 1;

    FlyingEdges::gridEdge &ge0 = gridEdges[k * ny + j];
    FlyingEdges::gridEdge &ge1 = gridEdges[k * ny + j + 1];

    uchar *ec0 = edgeCases + k * ny * (nx - 1) + j * (nx - 1);
    uchar *ec1 = edgeCases + k * ny * (nx - 1) + (j + 1) * (nx - 1);

    int xl = min(ge0.xl, ge1.xl);
    int xr = max(ge0.xr, ge1.xr);

    if (xl >= xr)
        return;

    int xstart = 0;
    int ystart = 0;

    uchar c0;
    uchar c1;

    for (int i = xl; i != xr; ++i) {
        c0 = ec0[i];
        c1 = ec1[i];

        // see if the edges are cut
        xstart += (c0 == 1 || c0 == 2);
        ystart += ((c0 == 0 && c1 == 1) || (c0 == 0 && c1 == 3) ||
                   (c0 == 1 && c1 == 2) || (c0 == 2 && c1 == 3));
    }

    if (xr == nx - 1) {
        ystart += ((c0 == 0 && c1 == 2) || (c0 == 0 && c1 == 3) ||
                   (c0 == 1 && c1 == 2) || (c0 == 1 && c1 == 3));
    }

    ge0.xstart = xstart;
    ge0.ystart = ystart;
    ge0.zstart = 0;
}

void FlyingEdges::pass2() {
    // pass2 calculates
    //   1) cubeCases for each block ray
    //   2) triCount for each block ray
    //   3) edgeRay count

    // 1st kernel: Calculate the 0, 1, 2 edge ray, cube cases, tricount
    // 2nd kernel: Calculate lost edges

    int ty = FE_BLOCK_WIDTH_Y;
    int tz = FE_BLOCK_WIDTH_Z;
    uint3 gridDim = make_uint3(((ny - 1) + ty - 1) / ty, ((nz - 1) + tz - 1) / tz, 1);
    uint3 blockDim = make_uint3(ty, tz, 1);

    get_cubeCases<<<gridDim, blockDim>>>(
            nx, ny, nz,
            edgeCases,
            gridEdges,   // modified
            triCounter,  // modified
            cubeCases);  // modified

    // POSSIBLE to do this here TODO
    // cudaFree(edgeCases);

    // TODO these can be launched and executed independently of each other
//    int bw = FE_BLOCK_WIDTH;

    // Making sure that the xz plane takes care of the (_, ny-1, nz-1) gridEdge
    // BE CAREFUL. xz takes care of corner. don't use (nz-1)
    // TODO: Check if CUDA launch parms here correct?
//    getGhostXZ<<<(nz + bw - 1) / bw, bw>>>(
//            nx, ny, nz,
//            edgeCases,
//            gridEdges);
//    getGhostXY<<<((ny - 1) + bw - 1) / bw, bw>>>(
//            nx, ny, nz,
//            edgeCases,
//            gridEdges);

    cudaDeviceSynchronize();

}


// Pass 3 of the algorithm uses information from pass 2 to determine how
// many triangles and points there are. It also sets up starting indices
// on each gridEdge. Once these sizes are determined, memory is allocated
// for storing triangles, points and normals.

// Prefix sum in blockwise approach
__global__
void blockAccum(
        int nx, int ny, int nz,
        int *triCounter,
        FlyingEdges::gridEdge *gridEdges,
        int *blockAccum) {
    int k = blockIdx.y * blockDim.y + threadIdx.y;

    // step 1: accumulate individual y thread
    // step 2: calc block sum
    // step 3: __syncthreads
    // step 4: add to individual y thread

    __shared__ int accum[4 * FE_BLOCK_WIDTH];

    if (k < nz) {
        int tmp;
        int accumX = 0;
        int accumY = 0;
        int accumZ = 0;
        int accumTri = 0;
        for (int j = 0; j != ny; ++j) {
            FlyingEdges::gridEdge &ge = gridEdges[k * ny + j];

            // Exclusive prefix sum approach: Every x/y/zstart stores the sum of previous info
            tmp = ge.xstart;
            ge.xstart = accumX;
            accumX += tmp;

            tmp = ge.ystart;
            ge.ystart = accumY;
            accumY += tmp;

            tmp = ge.zstart;
            ge.zstart = accumZ;
            accumZ += tmp;
        }

        if (k < nz - 1) {
            for (int j = 0; j != ny - 1; ++j) {
                // Why (-1)? Because triangles are recorded by cell-based way
                int &curTriCount = triCounter[k * (ny - 1) + j];

                tmp = curTriCount;
                curTriCount = accumTri;
                accumTri += tmp;
            }
        }

        accum[4 * threadIdx.y + 0] = accumX;
        accum[4 * threadIdx.y + 1] = accumY;
        accum[4 * threadIdx.y + 2] = accumZ;
        accum[4 * threadIdx.y + 3] = accumTri;
    }

    __syncthreads();

    // Prefix sum
    if (k < nz) {
        if (threadIdx.y == 0) // agh!
        {
            for (int idx = 1; idx != blockDim.y; ++idx) {
                accum[4 * idx + 0] += accum[4 * (idx - 1) + 0];
                accum[4 * idx + 1] += accum[4 * (idx - 1) + 1];
                accum[4 * idx + 2] += accum[4 * (idx - 1) + 2];
                accum[4 * idx + 3] += accum[4 * (idx - 1) + 3];
            }

            // answer for global accumulation
            blockAccum[4 * blockIdx.y + 0] = accum[4 * (blockDim.y - 1) + 0];
            blockAccum[4 * blockIdx.y + 1] = accum[4 * (blockDim.y - 1) + 1];
            blockAccum[4 * blockIdx.y + 2] = accum[4 * (blockDim.y - 1) + 2];
            blockAccum[4 * blockIdx.y + 3] = accum[4 * (blockDim.y - 1) + 3];
        }
    }
    __syncthreads();

    if (threadIdx.y == 0 || k >= nz)
        return;

    bool isEndK = k == nz - 1;
    for (int j = 0; j != ny - 1; ++j) {
        FlyingEdges::gridEdge &ge = gridEdges[k * ny + j];

        ge.xstart += accum[4 * (threadIdx.y - 1) + 0];
        ge.ystart += accum[4 * (threadIdx.y - 1) + 1];
        ge.zstart += accum[4 * (threadIdx.y - 1) + 2];

        // put z stuff here..
        if (!isEndK)
            triCounter[k * (ny - 1) + j] = accum[4 * (threadIdx.y - 1) + 3];
    }

    FlyingEdges::gridEdge &ge = gridEdges[k * ny + (ny - 1)];
    ge.xstart += accum[4 * (threadIdx.y - 1) + 0];
    ge.ystart += accum[4 * (threadIdx.y - 1) + 1];
    ge.zstart += accum[4 * (threadIdx.y - 1) + 2];
}

// Prefix sum in grid approach
__global__ // TODO can split up along j here easy enough.
void gridAccum(
        int nx, int ny, int nz, // which are needed TODO?
        int *triCounter,
        FlyingEdges::gridEdge *gridEdges,
        int *blockAccum) // used as input here, it is pre-fixed grid sum info
{
    // skip the first block
    int k = (blockIdx.z + 1) * blockDim.z + threadIdx.z;

    if (k >= nz)
        return;

    int addX = blockAccum[4 * blockIdx.z + 0];
    int addY = blockAccum[4 * blockIdx.z + 1];
    int addZ = blockAccum[4 * blockIdx.z + 2];
    int addTri = blockAccum[4 * blockIdx.z + 3];

    for (int j = 0; j != ny; ++j) {
        FlyingEdges::gridEdge &ge = gridEdges[k * ny + j];
        ge.xstart += addX;
        ge.ystart += addY;
        ge.zstart += addZ;
    }

    if (k >= nz - 1)
        return;

    for (int j = 0; j != ny - 1; ++j) {
        triCounter[k * (ny - 1) + j] += addTri;
    }
}

void FlyingEdges::pass3() {
    // Split the z axis
    // Kernel 1: calculate the accum values on block sync
    //           then accum individual values
    // Use that info accum each block (except the first one)
    // Kernel 2: just add values to individual threads
    int tz = FE_BLOCK_WIDTH;

    int numBlocks = (nz + tz - 1) / tz;

    // there are four because: xstart, ystart, zstart, triaccum
    int sizeBlocks = 4 * numBlocks * sizeof(int);

    uint3 gridDim = make_uint3(1, numBlocks, 1);
    uint3 blockDim = make_uint3(1, tz, 1);

    int *hostBlockAccum = (int *) malloc(sizeBlocks);
    for (int idx = 0; idx != 4 * numBlocks; ++idx) {
        hostBlockAccum[idx] = 0;
    }

    int *deviceBlockAccum;
    cudaMalloc(&deviceBlockAccum, sizeBlocks);

    cudaMemcpy(deviceBlockAccum, hostBlockAccum,
               sizeBlocks, cudaMemcpyHostToDevice);

    // Accumulate values locally

    blockAccum<<<gridDim, blockDim>>>(
            nx, ny, nz,
            triCounter,
            gridEdges,
            deviceBlockAccum);

    cudaMemcpy(hostBlockAccum, deviceBlockAccum,
               sizeBlocks, cudaMemcpyDeviceToHost);


    if (numBlocks != 1) {

        for (int i = 4; i != 4 * numBlocks; i += 4) {
            hostBlockAccum[i + 0] += hostBlockAccum[i - 4];
            hostBlockAccum[i + 1] += hostBlockAccum[i - 3];
            hostBlockAccum[i + 2] += hostBlockAccum[i - 2];
            hostBlockAccum[i + 3] += hostBlockAccum[i - 1];
        }
        // note: the last values in hostBlockAccum should contain total counts

        // The first block is done so it is ignored
        // and the last info in BlockAccum isn't needed (its the total counts)
        cudaMemcpy(deviceBlockAccum, hostBlockAccum,
                   sizeBlocks - 4 * sizeof(int), cudaMemcpyHostToDevice);

        // Accumulate values from other blocks
        gridDim = make_uint3(1, 1, numBlocks - 1);
        gridAccum<<<gridDim, blockDim>>>(
                nx, ny, nz,
                triCounter,
                gridEdges,
                deviceBlockAccum);
    }

    // Allocate memory for points, normals and tris
    // outputAllocated = true;
    numPoints = hostBlockAccum[4 * (numBlocks - 1) + 0] +
                hostBlockAccum[4 * (numBlocks - 1) + 1] +
                hostBlockAccum[4 * (numBlocks - 1) + 2];
    numTris = hostBlockAccum[4 * (numBlocks - 1) + 3];

    cudaMalloc(&points, 3 * sizeof(scalar_t) * numPoints);
    cudaMalloc(&normals, 3 * sizeof(scalar_t) * numPoints);
    cudaMalloc(&tris, 3 * sizeof(int) * numTris);

    // free memory used in this function
    free(hostBlockAccum);
    cudaFree(deviceBlockAccum);

    cudaDeviceSynchronize();
}

__device__
void computeGradient(
        int const &i, int const &j, int const &k,
        int const &nx, int const &ny, int const &nz,
        scalar_t *data,
        scalar_t *spacing,
        scalar_t *point) {
    // x0, x1, x2 denotes x, y, z
    scalar_t x0[2];
    scalar_t x1[2];
    scalar_t x2[2];
    scalar_t run[3];

    size_t dataIdx = k * nx * ny + j * nx + i;
    // The sequence is [0]>[1]
    if (i == 0) {
        x0[0] = data[dataIdx + 1];
        x0[1] = data[dataIdx];
        run[0] = spacing[0];
    } else if (i == (nx - 1)) {
        x0[0] = data[dataIdx];
        x0[1] = data[dataIdx - 1];
        run[0] = spacing[0];
    } else {
        x0[0] = data[dataIdx + 1];
        x0[1] = data[dataIdx - 1];
        run[0] = 2 * spacing[0];
    }

    if (j == 0) {
        x1[0] = data[dataIdx + nx];
        x1[1] = data[dataIdx];
        run[1] = spacing[1];
    } else if (j == (ny - 1)) {
        x1[0] = data[dataIdx];
        x1[1] = data[dataIdx - nx];
        run[1] = spacing[1];
    } else {
        x1[0] = data[dataIdx + nx];
        x1[1] = data[dataIdx - nx];
        run[1] = 2 * spacing[1];
    }

    if (k == 0) {
        x2[0] = data[dataIdx + nx * ny];
        x2[1] = data[dataIdx];
        run[2] = spacing[2];
    } else if (k == (nz - 1)) {
        x2[0] = data[dataIdx];
        x2[1] = data[dataIdx - nx * ny];
        run[2] = spacing[2];
    } else {
        x2[0] = data[dataIdx + nx * ny];
        x2[1] = data[dataIdx - nx * ny];
        run[2] = 2 * spacing[2];
    }

    point[0] = (x0[1] - x0[0]) / run[0];
    point[1] = (x1[1] - x1[0]) / run[1];
    point[2] = (x2[1] - x2[0]) / run[2];
}


__device__
void getCubeInfo(
        int i, int j, int k,
        int nx, int ny, int nz,
        scalar_t *pointValues, scalar_t *zeroPos, scalar_t *spacing,
        scalar_t *pointCube, scalar_t *isovalCube, scalar_t *gradCube) {
    isovalCube[0] = pointValues[k * ny * nx + j * nx + i];
    isovalCube[1] = pointValues[k * ny * nx + j * nx + i + 1];
    isovalCube[2] = pointValues[k * ny * nx + (j + 1) * nx + i + 1];
    isovalCube[3] = pointValues[k * ny * nx + (j + 1) * nx + i];
    isovalCube[4] = pointValues[(k + 1) * ny * nx + j * nx + i];
    isovalCube[5] = pointValues[(k + 1) * ny * nx + j * nx + i + 1];
    isovalCube[6] = pointValues[(k + 1) * ny * nx + (j + 1) * nx + (i + 1)];
    isovalCube[7] = pointValues[(k + 1) * ny * nx + (j + 1) * nx + i];

    // Get position info of 8 vertexes using spacing and zeroPos info
    scalar_t xpos = zeroPos[0] + i * spacing[0];
    scalar_t ypos = zeroPos[1] + j * spacing[1];
    scalar_t zpos = zeroPos[2] + k * spacing[2];

    pointCube[0 * 3 + 0] = xpos;
    pointCube[0 * 3 + 1] = ypos;
    pointCube[0 * 3 + 2] = zpos;

    pointCube[1 * 3 + 0] = xpos + spacing[0];
    pointCube[1 * 3 + 1] = ypos;
    pointCube[1 * 3 + 2] = zpos;

    pointCube[2 * 3 + 0] = xpos + spacing[0];
    pointCube[2 * 3 + 1] = ypos + spacing[1];
    pointCube[2 * 3 + 2] = zpos;

    pointCube[3 * 3 + 0] = xpos;
    pointCube[3 * 3 + 1] = ypos + spacing[1];
    pointCube[3 * 3 + 2] = zpos;

    pointCube[4 * 3 + 0] = xpos;
    pointCube[4 * 3 + 1] = ypos;
    pointCube[4 * 3 + 2] = zpos + spacing[2];

    pointCube[5 * 3 + 0] = xpos + spacing[0];
    pointCube[5 * 3 + 1] = ypos;
    pointCube[5 * 3 + 2] = zpos + spacing[2];

    pointCube[6 * 3 + 0] = xpos + spacing[0];
    pointCube[6 * 3 + 1] = ypos + spacing[1];
    pointCube[6 * 3 + 2] = zpos + spacing[2];

    pointCube[7 * 3 + 0] = xpos;
    pointCube[7 * 3 + 1] = ypos + spacing[1];
    pointCube[7 * 3 + 2] = zpos + spacing[2];

    computeGradient(i, j, k, nx, ny, nz, pointValues, spacing, gradCube + 3 * 0);
    computeGradient(i + 1, j, k, nx, ny, nz, pointValues, spacing, gradCube + 3 * 1);
    computeGradient(i + 1, j + 1, k, nx, ny, nz, pointValues, spacing, gradCube + 3 * 2);
    computeGradient(i, j + 1, k, nx, ny, nz, pointValues, spacing, gradCube + 3 * 3);
    computeGradient(i, j, k + 1, nx, ny, nz, pointValues, spacing, gradCube + 3 * 4);
    computeGradient(i + 1, j, k + 1, nx, ny, nz, pointValues, spacing, gradCube + 3 * 5);
    computeGradient(i + 1, j + 1, k + 1, nx, ny, nz, pointValues, spacing, gradCube + 3 * 6);
    computeGradient(i, j + 1, k + 1, nx, ny, nz, pointValues, spacing, gradCube + 3 * 7);
}

__device__
void interpolate(
        scalar_t const &weight,
        scalar_t *a,
        scalar_t *b,
        scalar_t *out) {
    out[0] = a[0] + (weight * (b[0] - a[0]));
    out[1] = a[1] + (weight * (b[1] - a[1]));
    out[2] = a[2] + (weight * (b[2] - a[2]));
}

__device__
void interpolateOnCube(
        uchar const &edge,
        scalar_t const &isoval,
        scalar_t *pts,
        scalar_t *isovals,
        scalar_t *out) {

    uchar i0 = edgeVertices[edge][0];
    uchar i1 = edgeVertices[edge][1];

    scalar_t weight = (isoval - isovals[i0]) / (isovals[i1] - isovals[i0]);
    interpolate(weight, pts + 3 * i0, pts + 3 * i1, out);
}

__global__
void getPointsAndNormals(
        int nx, int ny, int nz,
        scalar_t *pointValues, scalar_t *zeroPos, scalar_t *spacing,
        scalar_t isoval,
        FlyingEdges::gridEdge *gridEdges,
        int *triCounter,
        uchar *cubeCases,
        scalar_t *points, scalar_t *normals, int *tris) {
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    int k = blockIdx.y * blockDim.y + threadIdx.y;

    if (j >= ny - 1 || k >= nz - 1)
        return;

    FlyingEdges::gridEdge &ge0 = gridEdges[k * ny + j];
    FlyingEdges::gridEdge &ge1 = gridEdges[k * ny + j + 1];
    FlyingEdges::gridEdge &ge2 = gridEdges[(k + 1) * ny + j];
    FlyingEdges::gridEdge &ge3 = gridEdges[(k + 1) * ny + j + 1];

    int xl, xr;
    calcTrimValues(xl, xr, ge0, ge1, ge2, ge3);

    if (xl == xr)
        return;

    size_t triIdx = triCounter[k * (ny - 1) + j];
    // get the start cell
    uchar *curCubeCaseIds = cubeCases + (nx - 1) * (k * (ny - 1) + j);

    size_t x0counter = 0;
    size_t y0counter = 0;
    size_t z0counter = 0;
    size_t x1counter = 0;
    size_t y1counter = 0;
    size_t z1counter = 0;
    size_t x2counter = 0;
    size_t y2counter = 0;
    size_t z2counter = 0;
    size_t x3counter = 0;
    size_t y3counter = 0;
    size_t z3counter = 0;


    // Why (-2)?
    // Because we are using Central Difference Method to calculate gradients and normals
    // TODO: Since we have set boundary process for k = nz-1 and j = ny - 1, we must check if we have overlooked some boundary cases.
    bool isYEnd = (j == ny - 2);
    bool isZEnd = (k == nz - 2);

    scalar_t pointCube[8 * 3];
    scalar_t isovalCube[8];
    scalar_t gradCube[8 * 3];

    for (size_t i = xl; i != xr; ++i) {
        bool isXEnd = (i == nx - 2);

        uchar caseId = curCubeCaseIds[i];

        // All in or all out
        if (caseId == 0 || caseId == 255) {
            continue;
        }

        const bool *isCutCur = isCut[caseId]; // has 12 elements

        // fill out pointCube, isovalCube and gradCube
        getCubeInfo(i, j, k,
                    nx, ny, nz,
                    pointValues, zeroPos, spacing,
                    pointCube, isovalCube, gradCube);

        // Add Points and normals.
        // Calculate global indices for triangles
        int globalIdxs[12];
        if (isCutCur[0]) {
            int idx = ge0.xstart + x0counter;
            interpolateOnCube(0, isoval, pointCube, isovalCube, points + 3 * idx);
            interpolateOnCube(0, isoval, gradCube, isovalCube, normals + 3 * idx);
            globalIdxs[0] = idx;
            ++x0counter;
        }

        if (isCutCur[3]) {
            int idx = ge0.ystart + y0counter;
            interpolateOnCube(3, isoval, pointCube, isovalCube, points + 3 * idx);
            interpolateOnCube(3, isoval, gradCube, isovalCube, normals + 3 * idx);
            globalIdxs[3] = idx;
            ++y0counter;
        }

        if (isCutCur[8]) {
            int idx = ge0.zstart + z0counter;
            interpolateOnCube(8, isoval, pointCube, isovalCube, points + 3 * idx);
            interpolateOnCube(8, isoval, gradCube, isovalCube, normals + 3 * idx);
            globalIdxs[8] = idx;
            ++z0counter;
        }
        // We are going from xl to xr in x-incremental way. So the right 4 edges will be revisited in the next iteration.
        // Note:
        //   e1, e5, e9 and e11 will be visited in the next iteration
        //   when they are e3, e7, e8 and 10 respectively. So don't
        //   increment their counters. When the cube is an edge cube,
        //   their counters don't need to be incremented because they
        //   won't be used again.

        // Manage boundary cases if needed. Otherwise just update globalIdx.

        // e1 is mapped to y1counter
        if (isCutCur[1]) {
            int idx = ge0.ystart + y0counter;
            if (isXEnd) {

                interpolateOnCube(1, isoval, pointCube, isovalCube, points + 3 * idx);
                interpolateOnCube(1, isoval, gradCube, isovalCube, normals + 3 * idx);
                // y1counter counter doesn't need to be incremented
                // because it won't be used again.
            }
            globalIdxs[1] = idx;
        }

        // e9 is mapped to z1counter
        if (isCutCur[9]) {
            int idx = ge0.zstart + z0counter;
            if (isXEnd) {
                interpolateOnCube(9, isoval, pointCube, isovalCube, points + 3 * idx);
                interpolateOnCube(9, isoval, gradCube, isovalCube, normals + 3 * idx);
                // z1counter doesn't need to in incremented.
            }
            globalIdxs[9] = idx;
        }

        // e2 is mapped to x1counter
        if (isCutCur[2]) {
            int idx = ge1.xstart + x1counter;
            if (isYEnd) {
                interpolateOnCube(2, isoval, pointCube, isovalCube, points + 3 * idx);
                interpolateOnCube(2, isoval, gradCube, isovalCube, normals + 3 * idx);
            }
            globalIdxs[2] = idx;
            ++x1counter;
        }

        // e10 is mapped to z2counter
        if (isCutCur[10]) {
            int idx = ge1.zstart + z2counter;
            if (isYEnd) {
                interpolateOnCube(10, isoval, pointCube, isovalCube, points + 3 * idx);
                interpolateOnCube(10, isoval, gradCube, isovalCube, normals + 3 * idx);
            }
            globalIdxs[10] = idx;
            ++z2counter;
        }

        // e4 is mapped to x2counter
        if (isCutCur[4]) {
            int idx = ge2.xstart + x2counter;
            if (isZEnd) {
                interpolateOnCube(4, isoval, pointCube, isovalCube, points + 3 * idx);
                interpolateOnCube(4, isoval, gradCube, isovalCube, normals + 3 * idx);
            }
            globalIdxs[4] = idx;
            ++x2counter;
        }

        // e7 is mapped to y2counter
        if (isCutCur[7]) {
            int idx = ge2.ystart + y2counter;
            if (isZEnd) {
                interpolateOnCube(7, isoval, pointCube, isovalCube, points + 3 * idx);
                interpolateOnCube(7, isoval, gradCube, isovalCube, normals + 3 * idx);
            }
            globalIdxs[7] = idx;
            ++y2counter;
        }

        // e11 is mapped to z3counter
        if (isCutCur[11]) {
            int idx = ge1.zstart + z3counter;
            if (isXEnd and isYEnd) {
                interpolateOnCube(11, isoval, pointCube, isovalCube, points + 3 * idx);
                interpolateOnCube(11, isoval, gradCube, isovalCube, normals + 3 * idx);
                // z3counter does not need to be incremented.
            }
            globalIdxs[11] = idx;
        }

        // e5 is mapped to y3counter
        if (isCutCur[5]) {
            int idx = ge2.ystart + y3counter;
            if (isXEnd and isZEnd) {
                interpolateOnCube(5, isoval, pointCube, isovalCube, points + 3 * idx);
                interpolateOnCube(5, isoval, gradCube, isovalCube, normals + 3 * idx);
                // y3 counter does not need to be incremented.
            }
            globalIdxs[5] = idx;
        }

        // e6 is mapped to x3counter
        if (isCutCur[6]) {
            int idx = ge3.xstart + x3counter;
            if (isYEnd and isZEnd) {
                interpolateOnCube(6, isoval, pointCube, isovalCube, points + 3 * idx);
                interpolateOnCube(6, isoval, gradCube, isovalCube, normals + 3 * idx);
            }
            globalIdxs[6] = idx;
            ++x3counter;
        }

        // Add triangles
        const char *caseTri = caseTriangles[caseId];
        for (int idx = 0; caseTri[idx] != -1; idx += 3) {
//            tris[3 * triIdx + 0] = i;
//            tris[3 * triIdx + 1] = j;
//            tris[3 * triIdx + 2] = k;

            tris[3 * triIdx + 0] = globalIdxs[caseTri[idx]];
            tris[3 * triIdx + 1] = globalIdxs[caseTri[idx + 1]];
            tris[3 * triIdx + 2] = globalIdxs[caseTri[idx + 2]];
            ++triIdx;
        }
    }
}


void FlyingEdges::pass4() {
    // pass4 calculates points and normals
    //   1) points and normals

    // 1st kernel:           Calculate the main cube rays
    // 2nd and third kernel:

    int ty = FE_BLOCK_WIDTH_Y / 2;
    int tz = FE_BLOCK_WIDTH_Z / 2;
    uint3 gridDim = make_uint3(((ny - 1) + ty - 1) / ty, ((nz - 1) + tz - 1) / tz, 1);
    uint3 blockDim = make_uint3(ty, tz, 1);

    std::cout << gridDim.x << ", " << gridDim.y << ", " << gridDim.z << std::endl;
    std::cout << blockDim.x << ", " << blockDim.y << ", " << blockDim.z << std::endl;

    getPointsAndNormals<<<gridDim, blockDim>>>(
            nx, ny, nz,                                    // input
            pointValues, zero_pos, spacing,                 // input
            isoval,                                        // input
            gridEdges, triCounter, cubeCases,              // input
            points, normals, tris);                        // output

}

void FlyingEdges::moveOutput() {
    host_points = (scalar_t *) malloc(3 * sizeof(scalar_t) * numPoints);
    host_normals = (scalar_t *) malloc(3 * sizeof(scalar_t) * numPoints);
    host_tris = (scalar_t *) malloc(3 * sizeof(scalar_t) * numTris);
    cudaMemcpy(host_points, points, 3 * sizeof(scalar_t) * numPoints, cudaMemcpyDeviceToHost);
    cudaMemcpy(host_normals, normals, 3 * sizeof(scalar_t) * numPoints, cudaMemcpyDeviceToHost);
    cudaMemcpy(host_tris, tris, 3 * sizeof(scalar_t) * numPoints, cudaMemcpyDeviceToHost);

    for (int i = 0; i < numPoints; i++) {
        std::cout << "i: " << host_points[i * 3 + 0] << "j: " << host_points[i * 3 + 1] << "k: " << host_points[i * 3 + 2];
                  << std::endl;
    }
    for (int i = 0; i < numTris; i++) {
        std::cout << "i: " << host_tris[i * 3 + 0] << "j: " << host_tris[i * 3 + 1] << "k: " << host_tris[i * 3 + 2];
        << std::endl;
    }
}
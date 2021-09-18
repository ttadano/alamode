/*
 interpolation.h

 Copyright (c) 2021 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include <vector>
#include <cmath>
#include <iomanip>
#include "memory.h"
#include "mathfunctions.h"

namespace PHON_NS {
class TrilinearInterpolation {
 public:
    TrilinearInterpolation()
    = default;

    TrilinearInterpolation(const unsigned int ngrid_coarse_in[3],
                           const unsigned int ngrid_dense_in[3])
    {
        for (auto i = 0; i < 3; ++i) {
            grid_c[i] = ngrid_coarse_in[i];
            grid_f[i] = ngrid_dense_in[i];
        }
        ngrid_c = grid_c[0] * grid_c[1] * grid_c[2];
        ngrid_f = grid_f[0] * grid_f[1] * grid_f[2];

        if (cubes) deallocate(cubes);
        if (xf) deallocate(xf);
        if (xc) deallocate(xc);

        allocate(cubes, ngrid_f, 8);
        allocate(xf, ngrid_f, 3);
        allocate(xc, ngrid_c, 3);
    }

    void setup()
    {
        set_grid(grid_c, xc);
        set_grid(grid_f, xf);
        set_cubes();
    };

    template<typename T>
    void interpolate(const T *val_c, T *val_f, const bool regular_grid = true)
    {

        T v_cubes[8];
        T tx, ty, tz;
        T invdel_x, invdel_y, invdel_z;

        if (regular_grid) {
            invdel_x = static_cast<T>(grid_c[0]);
            invdel_y = static_cast<T>(grid_c[1]);
            invdel_z = static_cast<T>(grid_c[2]);

            for (auto i = 0; i < ngrid_f; ++i) {
                for (auto j = 0; j < 8; ++j) {
                    v_cubes[j] = val_c[cubes[i][j]];
                }
                tx = static_cast<T>(xf[i][0] - xc[cubes[i][0]][0]) * invdel_x;
                ty = static_cast<T>(xf[i][1] - xc[cubes[i][0]][1]) * invdel_y;
                tz = static_cast<T>(xf[i][2] - xc[cubes[i][0]][2]) * invdel_z;

                const auto c0 = BilinearInterpolation(tx, ty,
                                                      v_cubes[0], v_cubes[1],
                                                      v_cubes[2], v_cubes[3]);
                const auto c1 = BilinearInterpolation(tx, ty,
                                                      v_cubes[4], v_cubes[5],
                                                      v_cubes[6], v_cubes[7]);
                val_f[i] = (c1 - c0) * tz + c0;
            }
        }
    }
    ~TrilinearInterpolation()
    {
        if (cubes) deallocate(cubes);
        if (xf) deallocate(xf);
        if (xc) deallocate(xc);
    };

 private:
    unsigned int grid_c[3]{};
    unsigned int grid_f[3]{};
    unsigned int ngrid_f, ngrid_c;
    int **cubes = nullptr;
    double **xf = nullptr; // coordinate of fine grid
    double **xc = nullptr; // coordinate of coarse grid


    static void set_grid(const unsigned int ngrid_in[3],
                  double **x_out)
    {
        size_t ik;
        double invn[3];
        for (auto i = 0; i < 3; ++i) {
            invn[i] = 1.0 / static_cast<double>(ngrid_in[i]);
        }
        for (size_t ix = 0; ix < ngrid_in[0]; ++ix) {
            for (size_t iy = 0; iy < ngrid_in[1]; ++iy) {
                for (size_t iz = 0; iz < ngrid_in[2]; ++iz) {
                    ik = iz + iy * ngrid_in[2] + ix * ngrid_in[2] * ngrid_in[1];
                    x_out[ik][0] = static_cast<double>(ix) * invn[0];
                    x_out[ik][1] = static_cast<double>(iy) * invn[1];
                    x_out[ik][2] = static_cast<double>(iz) * invn[2];
                }
            }
        }
    }
    void set_cubes()
    {
        int iloc[2], jloc[2], kloc[2];
        double dn_c[3];
        double tmp[3];
        int n23 = static_cast<int>(grid_c[1] * grid_c[2]);
        int igrid[3];

        for (auto i = 0; i < 3; ++i) {
            dn_c[i] = static_cast<double>(grid_c[i]);
            igird[i] = static_cast<int>(grid_c[i]);
        }

        for (auto i = 0; i < ngrid_f; ++i) {
            for (auto j = 0; j < 3; ++j) tmp[j] = xf[i][j] * dn_c[j];
            iloc[0] = nint(std::floor(tmp[0]));
            iloc[1] = nint(std::ceil(tmp[0]));
            jloc[0] = nint(std::floor(tmp[1]));
            jloc[1] = nint(std::ceil(tmp[1]));
            kloc[0] = nint(std::floor(tmp[2]));
            kloc[1] = nint(std::ceil(tmp[2]));

            if (iloc[1] == iloc[0]) ++iloc[1];
            if (jloc[1] == jloc[0]) ++jloc[1];
            if (kloc[1] == kloc[0]) ++kloc[1];

            iloc[0] = iloc[0] % igrid[0];
            iloc[1] = iloc[1] % igrid[0];
            jloc[0] = jloc[0] % igrid[1];
            jloc[1] = jloc[1] % igrid[1];
            kloc[0] = kloc[0] % igrid[2];
            kloc[1] = kloc[1] % igrid[2];

            cubes[i][0] = kloc[0] + jloc[0] * igrid[2] + iloc[0] * n23; // index of c000
            cubes[i][1] = kloc[0] + jloc[0] * igrid[2] + iloc[1] * n23; // index of c100
            cubes[i][2] = kloc[0] + jloc[1] * igrid[2] + iloc[0] * n23; // index of c010
            cubes[i][3] = kloc[0] + jloc[1] * igrid[2] + iloc[1] * n23; // index of c110
            cubes[i][4] = kloc[1] + jloc[0] * igrid[2] + iloc[0] * n23; // index of c001
            cubes[i][5] = kloc[1] + jloc[0] * igrid[2] + iloc[1] * n23; // index of c101
            cubes[i][6] = kloc[1] + jloc[1] * igrid[2] + iloc[0] * n23; // index of c011
            cubes[i][7] = kloc[1] + jloc[1] * igrid[2] + iloc[1] * n23; // index of c111
        }
    };

    template<typename T>
    T BilinearInterpolation(const T tx, const T ty,
                            const T c00, const T c10,
                            const T c01, const T c11)
    {
        return c00 + (c10 - c00) * tx + (c01 - c00) * ty
              + (c11 - c01 - c10 + c00) * tx * ty;
    }
};

}

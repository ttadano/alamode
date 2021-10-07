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
class TriLinearInterpolator {
 public:
    TriLinearInterpolator()
    = default;

    TriLinearInterpolator(const unsigned int ngrid_coarse_in[3],
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
        // set_cubes(); we can get cubes on the go
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

                const auto c0 = BiLinearInterpolation(tx, ty,
                                                      v_cubes[0], v_cubes[1],
                                                      v_cubes[2], v_cubes[3]);
                const auto c1 = BiLinearInterpolation(tx, ty,
                                                      v_cubes[4], v_cubes[5],
                                                      v_cubes[6], v_cubes[7]);
                val_f[i] = (c1 - c0) * tz + c0;
            }
        }
    }

    template<typename T>
    void interpolate_new(const T *val_c, T *val_f, const bool regular_grid = true)
    {

        T v_cubes[8];
        int corner_index[8];
        double **corner_coord;

        allocate(corner_coord, 8, 3);

        if (regular_grid) {

            for (auto i = 0; i < ngrid_f; ++i) {

                get_corners(xf[i], corner_index, corner_coord)

                for (auto j = 0; j < 8; ++j) {
                    v_cubes[j] = val_c[corner_index[j]];
                }

                val_f[i] = TriLinearInterpolation(xf[i], corner_coord, v_cubes);

            }
        }

        deallocate(corner_coord);
    }

    template<typename T>
    void improved_interpolate(const T *val_c, T *val_f, const unsigned is, const bool regular_grid = true)
    {
        // take branch number as input

        T v_cubes[8];
        T tx, ty, tz;
        bool contain_gamma;

        if (regular_grid) {

            for (auto i = 0; i < ngrid_f; ++i) {

                if (i == 0 && is < 3) {
                    val_f[i] = eps;
                    continue;
                }
                
                contain_gamma = false;

                for (auto j = 0; j < 8; ++j) {
                    if (cubes[i][j] == 0) {
                        contain_gamma = true;
                        break;
                    }
                }

                if (contain_gamma && is < 3) {

                    T val_sum{}; 
                    int counter = 0;

                    int neighbor_corners[8];
                    double **neighbor_coord;
                    allocate(neighbor_coord, 8, 3);

                    for (auto j = 0; j < 8; ++j) {
                        if (cubes[i][j] == 0) continue;

                        // we find the 7 neighboring cubes
                        // and averaged the interpolated value

                        double tmp_coord[3];

                        for (auto k = 0; k < 3; ++k) {
                            if (xc[cubes[i][j]][k] < 0.5) {
                                tmp_coord[k] = xf[i][k] + xc[cubes[i][j]][k];
                            } else {
                                tmp_coord[k] = xf[i][k] - (1.0 - xc[cubes[i][j]][k]);
                            }
                        }

                        get_corners(tmp_coord, neighbor_corners, neighbor_coord);

                        val_sum += TriLinearInterpolation(i, neighbor_corners,  neighbor_coord, val_c);
                        counter += 1;
                    }

                    deallocate(neighbor_coord);

                    val_f[i] = val_sum / static_cast<T>(counter);

                } else {
                    int neighbor_corners[8];
                    double **neighbor_coord;
                    allocate(neighbor_coord, 8, 3);
                    
                    get_corners(xf[i], neighbor_corners, neighbor_coord);

                    val_f[i] = TriLinearInterpolation(i, neighbor_corners, neighbor_coord, val_c);
                    deallocate(neighbor_coord);

                }
            }

        } // regular grid
    }
    
    
    ~TriLinearInterpolator()
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
        // gamma point will always be index 1
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
        double **tmp;
        allocate(tmp, 8, 3);
        for (auto i = 0; i < ngrid_f; ++i) {

        }
        
        get_cornerindex(i, cubes[i]);
    };

    void get_cornerindex(unsigned i_f, int* corners) 
    {
        // wrap around get_corners
        double **tmp;
        allocate(tmp, 8, 3);
        get_corners(xf[i_f], corners, tmp);
        deallocate(tmp);
    }

    void get_corners(double* xk_i, int* corner_index, double** corner_coord) {
        
        int iloc[2], jloc[2], kloc[2];
        double dn_c[3];
        double tmp[3];
        int n23 = static_cast<int>(grid_c[1] * grid_c[2]);
        int igrid[3];

        for (auto i = 0; i < 3; ++i) {
            dn_c[i] = static_cast<double>(grid_c[i]);
            igrid[i] = static_cast<int>(grid_c[i]);
        }

        for (auto j = 0; j < 3; ++j) tmp[j] = xk_i[j] * dn_c[j];
        iloc[0] = nint(std::floor(tmp[0]));
        iloc[1] = nint(std::ceil(tmp[0]));
        jloc[0] = nint(std::floor(tmp[1]));
        jloc[1] = nint(std::ceil(tmp[1]));
        kloc[0] = nint(std::floor(tmp[2]));
        kloc[1] = nint(std::ceil(tmp[2]));

        if (iloc[1] == iloc[0]) ++iloc[1];
        if (jloc[1] == jloc[0]) ++jloc[1];
        if (kloc[1] == kloc[0]) ++kloc[1];

        corner_coord[0][0] = static_cast<double>(iloc[0]) / dn_c[0];
        corner_coord[0][1] = static_cast<double>(jloc[0]) / dn_c[1];
        corner_coord[0][2] = static_cast<double>(kloc[0]) / dn_c[2];

        corner_coord[1][0] = static_cast<double>(iloc[1]) / dn_c[0];
        corner_coord[1][1] = static_cast<double>(jloc[0]) / dn_c[1];
        corner_coord[1][2] = static_cast<double>(kloc[0]) / dn_c[2];

        corner_coord[2][0] = static_cast<double>(iloc[0]) / dn_c[0];
        corner_coord[2][1] = static_cast<double>(jloc[1]) / dn_c[1];
        corner_coord[2][2] = static_cast<double>(kloc[0]) / dn_c[2];

        corner_coord[3][0] = static_cast<double>(iloc[1]) / dn_c[0];
        corner_coord[3][1] = static_cast<double>(jloc[1]) / dn_c[1];
        corner_coord[3][2] = static_cast<double>(kloc[0]) / dn_c[2];

        corner_coord[4][0] = static_cast<double>(iloc[0]) / dn_c[0];
        corner_coord[4][1] = static_cast<double>(jloc[0]) / dn_c[1];
        corner_coord[4][2] = static_cast<double>(kloc[1]) / dn_c[2];
        
        corner_coord[5][0] = static_cast<double>(iloc[1]) / dn_c[0];
        corner_coord[5][1] = static_cast<double>(jloc[0]) / dn_c[1];
        corner_coord[5][2] = static_cast<double>(kloc[1]) / dn_c[2];
        
        corner_coord[6][0] = static_cast<double>(iloc[0]) / dn_c[0];
        corner_coord[6][1] = static_cast<double>(jloc[1]) / dn_c[1];
        corner_coord[6][2] = static_cast<double>(kloc[1]) / dn_c[2];
        
        corner_coord[7][0] = static_cast<double>(iloc[1]) / dn_c[0];
        corner_coord[7][1] = static_cast<double>(jloc[1]) / dn_c[1];
        corner_coord[7][2] = static_cast<double>(kloc[1]) / dn_c[2];

        iloc[0] = iloc[0] % igrid[0];
        iloc[1] = iloc[1] % igrid[0];
        jloc[0] = jloc[0] % igrid[1];
        jloc[1] = jloc[1] % igrid[1];
        kloc[0] = kloc[0] % igrid[2];
        kloc[1] = kloc[1] % igrid[2];

        corner_index[0] = kloc[0] + jloc[0] * igrid[2] + iloc[0] * n23; // index of c000
        corner_index[1] = kloc[0] + jloc[0] * igrid[2] + iloc[1] * n23; // index of c100
        corner_index[2] = kloc[0] + jloc[1] * igrid[2] + iloc[0] * n23; // index of c010
        corner_index[3] = kloc[0] + jloc[1] * igrid[2] + iloc[1] * n23; // index of c110
        corner_index[4] = kloc[1] + jloc[0] * igrid[2] + iloc[0] * n23; // index of c001
        corner_index[5] = kloc[1] + jloc[0] * igrid[2] + iloc[1] * n23; // index of c101
        corner_index[6] = kloc[1] + jloc[1] * igrid[2] + iloc[0] * n23; // index of c011
        corner_index[7] = kloc[1] + jloc[1] * igrid[2] + iloc[1] * n23; // index of c111

    }

    template<typename T>
    T TriLinearInterpolation(int i, int *corners, const T *val_c) 
    {
        T v_cubes[8];
        for (auto j = 0; j < 8; ++j) {
            v_cubes[j] = val_c[corners[j]];
        }
        T tx = static_cast<T>(xf[i][0] - xc[corners[0]][0]) * static_cast<T>(grid_c[0]);
        T ty = static_cast<T>(xf[i][1] - xc[corners[0]][1]) * static_cast<T>(grid_c[1]);
        T tz = static_cast<T>(xf[i][2] - xc[corners[0]][2]) * static_cast<T>(grid_c[2]);

        const auto c0 = BiLinearInterpolation(tx, ty,
                                              v_cubes[0], v_cubes[1],
                                              v_cubes[2], v_cubes[3]);
        const auto c1 = BiLinearInterpolation(tx, ty,
                                              v_cubes[4], v_cubes[5],
                                              v_cubes[6], v_cubes[7]);
        return (c1 - c0) * tz + c0;
    }

    template<typename T>
    T TriLinearInterpolation(int i, int *corners, double **corners_coord, const T *val_c) 
    {
        T v_cubes[8];
        for (auto j = 0; j < 8; ++j) {
            v_cubes[j] = val_c[corners[j]];
        }
        T tx = static_cast<T>(xf[i][0] - corners_coord[0][0]) * static_cast<T>(grid_c[0]);
        T ty = static_cast<T>(xf[i][1] - corners_coord[0][1]) * static_cast<T>(grid_c[1]);
        T tz = static_cast<T>(xf[i][2] - corners_coord[0][2]) * static_cast<T>(grid_c[2]);

        const auto c0 = BiLinearInterpolation(tx, ty,
                                              v_cubes[0], v_cubes[1],
                                              v_cubes[2], v_cubes[3]);
        const auto c1 = BiLinearInterpolation(tx, ty,
                                              v_cubes[4], v_cubes[5],
                                              v_cubes[6], v_cubes[7]);
        return (c1 - c0) * tz + c0;
    }

    template<typename T>
    T TriLinearInterpolation(double *center, double **corners_coord, const T *val_corner) 
    {
        T tx = static_cast<T>(center[0] - corners_coord[0][0]) * static_cast<T>(grid_c[0]);
        T ty = static_cast<T>(center[1] - corners_coord[0][1]) * static_cast<T>(grid_c[1]);
        T tz = static_cast<T>(center[2] - corners_coord[0][2]) * static_cast<T>(grid_c[2]);

        const auto c0 = BiLinearInterpolation(tx, ty,
                                              val_corner[0], val_corner[1],
                                              val_corner[2], val_corner[3]);
        const auto c1 = BiLinearInterpolation(tx, ty,
                                              val_corner[4], val_corner[5],
                                              val_corner[6], val_corner[7]);
        return (c1 - c0) * tz + c0;
    }

    template<typename T>
    T BiLinearInterpolation(const T tx, const T ty,
                            const T c00, const T c10,
                            const T c01, const T c11)
    {
        return c00 + (c10 - c00) * tx + (c01 - c00) * ty
              + (c11 - c01 - c10 + c00) * tx * ty;
    }

};

}

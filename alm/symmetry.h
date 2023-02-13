/*
 symmetry.h

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include <string>
#include <vector>
#include "system.h"
#include "timer.h"
#include <Eigen/Core>

namespace ALM_NS {
class SymmetryOperation {
public:
    Eigen::Matrix3i rotation;         // in lattice basis
    Eigen::Vector3d tran;             // in Cartesian basis
    Eigen::Matrix3d rotation_cart;  // in Cartesian basis
    bool compatible_with_lattice;
    bool compatible_with_cartesian;
    bool is_translation;

    SymmetryOperation();

    SymmetryOperation(const Eigen::Matrix3i &rot_in,
                      const Eigen::Vector3d &tran_in,
                      const Eigen::Matrix3d &rot_cart_in,
                      const bool compatibility_lat,
                      const bool compatibility_cart,
                      const bool is_trans_in)
    {
        rotation = rot_in;
        rotation_cart = rot_cart_in;
        tran = tran_in;
        compatible_with_lattice = compatibility_lat;
        compatible_with_cartesian = compatibility_cart;
        is_translation = is_trans_in;
    }

    // Operator definition to sort
    bool operator<(const SymmetryOperation &a) const
    {
        std::vector<double> v1, v2;
        for (auto i = 0; i < 3; ++i) {
            for (auto j = 0; j < 3; ++j) {
                v1.push_back(static_cast<double>(rotation(i, j)));
                v2.push_back(static_cast<double>(a.rotation(i, j)));
            }
        }
        for (auto i = 0; i < 3; ++i) {
            if (tran[i] < 0.0) {
                v1.push_back(1.0 + tran[i]);
            } else {
                v1.push_back(tran[i]);
            }
            if (a.tran[i] < 0.0) {
                v2.push_back(1.0 + a.tran[i]);
            } else {
                v2.push_back(a.tran[i]);
            }
        }
        return std::lexicographical_compare(v1.begin(), v1.end(),
                                            v2.begin(), v2.end());
    }
};

class RotationMatrix {
public:
    //int mat[3][3];
    Eigen::Matrix3i mat;

    RotationMatrix();

    RotationMatrix(const int rot[3][3])
    {
        for (auto i = 0; i < 3; ++i) {
            for (auto j = 0; j < 3; ++j) {
                mat(i, j) = rot[i][j];
            }
        }
    }

    RotationMatrix(const Eigen::Matrix3i &rot_in)
    {
        mat = rot_in;
    }
};

class Maps {
public:
    int atom_num;
    int tran_num;
};

class Symmetry {
public:
    Symmetry();

    ~Symmetry();

    void init(const std::unique_ptr<System> &system,
              const int verbosity,
              std::unique_ptr<Timer> &timer);

    double get_tolerance() const;

    void set_tolerance(const double);

    int get_print_symmetry() const;

    void set_print_symmetry(const int);

    const std::vector<Maps> &get_map_s2p() const;

    const std::vector<std::vector<int>> &get_map_p2s() const;

    const std::vector<SymmetryOperation> &get_SymmData() const;

    const std::vector<std::vector<int>> &get_map_sym() const;

    const std::vector<int> &get_symnum_tran() const;

    size_t get_nsym() const;

    size_t get_ntran() const;

    size_t get_nat_prim() const;

private:
    size_t nsym, ntran, nat_prim;
    std::vector<std::vector<int>> map_sym;   // [nat_base, nsym]
    std::vector<std::vector<int>> map_p2s;   // [nat_prim, ntran]
    std::vector<Maps> map_s2p;               // [nat_base]
    std::vector<SymmetryOperation> SymmData; // [nsym]
    std::vector<int> symnum_tran;            // [ntran]

    double tolerance;
    bool use_internal_symm_finder;
    int printsymmetry;

    void set_default_variables();

    void deallocate_variables();

    void setup_symmetry_operation(const Cell &,
                                  const int [3],
                                  const std::vector<std::vector<unsigned int>> &,
                                  const Spin &,
                                  const int);

    void gen_mapping_information(const Cell &,
                                 const std::vector<std::vector<unsigned int>> &);

    void findsym_alm(const Cell &,
                     const int [3],
                     const std::vector<std::vector<unsigned int>> &,
                     const Spin &);

    int findsym_spglib(const Cell &,
                       const std::vector<std::vector<unsigned int>> &,
                       const Spin &,
                       std::string &);

    bool is_translation(const int [3][3]) const;

    bool is_translation(const Eigen::Matrix3i &rot) const;

    bool is_proper(const double [3][3]) const;

    bool is_proper(const Eigen::Matrix3d &rot) const;

    void symop_in_cart(double [3][3],
                       const int [3][3],
                       const double [3][3],
                       const double [3][3]) const;

    void symop_in_cart(Eigen::Matrix3d &rot_cart,
                       const Eigen::Matrix3i &rot_lattice,
                       const Eigen::Matrix3d &lavec,
                       const Eigen::Matrix3d &rlavec) const;

    void print_symminfo_stdout() const;

    template<typename T>
    bool is_compatible(const T [3][3],
                       double tolerance_zero = 1.0e-5);

    template<typename T>
    bool is_compatible(const Eigen::MatrixBase<T> &mat,
                       double tolerance_zero = 1.0e-5);

    void find_lattice_symmetry(const Eigen::Matrix3d &aa,
                               std::vector<RotationMatrix> &) const;

    void find_crystal_symmetry(const Cell &,
                               const std::vector<std::vector<unsigned int>> &,
                               const int [3],
                               const Spin &,
                               const std::vector<RotationMatrix> &);

    void set_primitive_lattice(const double aa[3][3],
                               const size_t nat,
                               const int *kd,
                               double **x,
                               double aa_prim[3][3],
                               size_t &nat_prim,
                               int *kd_prim,
                               double **x_prim,
                               const double symprec) const;

    std::string file_sym;
};
}

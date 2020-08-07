/*
 symmetry_core.h

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include "pointers.h"
#include <string>
#include <vector>

namespace PHON_NS {
    class SymmetryOperation {
    public:
        int rot[3][3];
        double tran[3];

        SymmetryOperation();

        SymmetryOperation(const int rot_in[3][3],
                          const double tran_in[3])
        {
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    rot[i][j] = rot_in[i][j];
                }
            }
            for (int i = 0; i < 3; ++i) {
                tran[i] = tran_in[i];
            }
        }

        bool operator<(const SymmetryOperation &a) const
        {
            std::vector<double> v1, v2;
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    v1.push_back(static_cast<double>(rot[i][j]));
                    v2.push_back(static_cast<double>(a.rot[i][j]));
                }
            }
            for (int i = 0; i < 3; ++i) {
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
        int mat[3][3];

        RotationMatrix();

        RotationMatrix(const int rot[3][3])
        {
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    mat[i][j] = rot[i][j];
                }
            }
        }
    };

    class SymmetryOperationWithMapping {
    public:
        std::vector<double> rot;            // Rotation matrix in Cartesian basis
        std::vector<double> rot_real;       // Rotation matrix in fractional basis
        std::vector<double> rot_reciprocal; // Rotation matrix in reciprocal (fractional) basis
        std::vector<unsigned int> mapping;
        double shift[3]; // Translation vector in fractional basis

        SymmetryOperationWithMapping();

        SymmetryOperationWithMapping(const double S[3][3],
                                     const double T[3][3],
                                     const double R[3][3],
                                     unsigned int *mapping_info,
                                     const unsigned int n,
                                     const double shift_in[3])
        {
            unsigned int i, j;

            for (i = 0; i < 3; ++i) {
                for (j = 0; j < 3; ++j) {
                    rot.push_back(S[i][j]);
                    rot_real.push_back(T[i][j]);
                    rot_reciprocal.push_back(R[i][j]);
                }
            }
            for (i = 0; i < n; ++i) {
                mapping.push_back(mapping_info[i]);
            }
            for (i = 0; i < 3; ++i) {
                shift[i] = shift_in[i];
            }
        }
    };


    class Symmetry : protected Pointers {
    public:
        Symmetry(class PHON *);

        ~Symmetry();

        unsigned int nsym;
        bool symmetry_flag, time_reversal_sym;
        int trev_sym_mag;
        bool printsymmetry;
        double tolerance;
        std::vector<SymmetryOperation> SymmList;
        std::vector<SymmetryOperationWithMapping> SymmListWithMap;

        void setup_symmetry();

    private:

        std::string file_sym;

        void set_default_variables();

        void setup_symmetry_operation(int,
                                      unsigned int &,
                                      double [3][3],
                                      double [3][3],
                                      double **,
                                      unsigned int *);

        void findsym(int,
                     double [3][3],
                     double **,
                     std::vector<SymmetryOperation> &) const;

        void gensym_withmap(double **,
                            const unsigned int *);

        bool is_proper(double [3][3]) const;

        void find_lattice_symmetry(double [3][3],
                                   std::vector<RotationMatrix> &) const;

        void find_crystal_symmetry(
                int,
                std::vector<unsigned int> *,
                double **x,
                const std::vector<RotationMatrix> &,
                std::vector<SymmetryOperation> &) const;

        void broadcast_symmlist(std::vector<SymmetryOperation> &) const;
    };
}

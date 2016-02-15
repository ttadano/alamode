/*
 symmetry_core.h

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include "pointers.h"
#include <string>
#include <vector>
#include <fstream>

#ifdef _USE_EIGEN
#include <Eigen/Core>
#endif

namespace PHON_NS {

    class SymmetryOperation {
    public:
        int rot[3][3];
        double tran[3];

        SymmetryOperation();

        // Declaration construction

        SymmetryOperation(const int rot_in[3][3], const double tran_in[3])
        {
            for (int i = 0; i < 3; ++i){
                for (int j = 0; j < 3; ++j){
                    rot[i][j] = rot_in[i][j];
                }
            }
            for (int i = 0; i < 3; ++i){
                tran[i] = tran_in[i];
            }
        }
    };

    class RotationMatrix {
    public:
        int mat[3][3];

        RotationMatrix();
        RotationMatrix(const int rot[3][3]) {
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    mat[i][j] = rot[i][j];
                }
            }
        }
    };

    class SymmetryOperationWithMapping {
    public:
        std::vector<double> rot; // Rotation matrix in Cartesian basis
        std::vector<double> rot_real;  // Rotation matrix in fractional basis
        std::vector<double> rot_reciprocal; // Rotation matrix in reciprocal (fractional) basis
        std::vector<unsigned int> mapping; 
        double shift[3]; // Translation vector in fractional basis

        SymmetryOperationWithMapping();

        SymmetryOperationWithMapping(const SymmetryOperationWithMapping &a) 
        {
            for (std::vector<double>::const_iterator p = a.rot.begin(); p != a.rot.end(); ++p) {
                rot.push_back(*p);
            }
            for (std::vector<double>::const_iterator p = a.rot_real.begin(); p != a.rot_real.end(); ++p) {
                rot_real.push_back(*p);
            }
            for (std::vector<double>::const_iterator p = a.rot_reciprocal.begin(); p != a.rot_reciprocal.end(); ++p) {
                rot_reciprocal.push_back(*p);
            }
            for (std::vector<unsigned int>::const_iterator p = a.mapping.begin(); p != a.mapping.end(); ++p) {
                mapping.push_back(*p);
            }
            for (unsigned int i = 0; i < 3; ++i) {
                shift[i] = a.shift[i];
            }
        }

        SymmetryOperationWithMapping(const double S[3][3], const double T[3][3], const double R[3][3], 
            unsigned int *mapping_info, const unsigned int n, const double shift_in[3])
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

    inline bool operator<(const SymmetryOperationWithMapping a, const SymmetryOperationWithMapping b) {
        return std::lexicographical_compare(a.rot.begin(), a.rot.end(), b.rot.begin(), b.rot.end());
    }

    class Symmetry: protected Pointers {
    public:
        Symmetry(class PHON *);
        ~Symmetry();

        unsigned int nsym;
        bool symmetry_flag, time_reversal_sym;
        bool printsymmetry;
        double tolerance;
        std::vector<SymmetryOperation> SymmList;
        std::vector<SymmetryOperationWithMapping> SymmListWithMap;

        void setup_symmetry();

    private:

        std::string file_sym;

        void setup_symmetry_operation(int, unsigned int&, double[3][3], double[3][3], 
            double **, unsigned int *);
        void findsym(int, double [3][3], double **, std::vector<SymmetryOperation> &);
        void gensym_withmap(double **, unsigned int *);
        void find_lattice_symmetry(double [3][3], std::vector<RotationMatrix> &);
        void find_crystal_symmetry(int, int, std::vector<unsigned int> *, double **x, 
            std::vector<RotationMatrix>, std::vector<SymmetryOperation> &);
        void find_nnp_for_translation(unsigned int &, std::vector<SymmetryOperation>);
        void broadcast_symmlist(std::vector<SymmetryOperation> &);
    };
}

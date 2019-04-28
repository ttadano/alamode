/*
 writer.h

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include <string>
#include <vector>
#include "alm.h"

namespace ALM_NS
{
    class AtomProperty
    {
    public:
        double x, y, z;
        int kind;
        size_t atom, tran;

        AtomProperty() = default;;

        AtomProperty(const AtomProperty &other) = default;

        AtomProperty(const double *pos,
                     const int kind_in,
                     const int atom_in,
                     const int tran_in)
        {
            x = pos[0];
            y = pos[1];
            z = pos[2];
            kind = kind_in;
            atom = atom_in;
            tran = tran_in;
        }
    };

    class SystemInfo
    {
    public:
        double lattice_vector[3][3];
        std::vector<AtomProperty> atoms;
        size_t nat, natmin, ntran;
        size_t nspecies;

        SystemInfo() = default;;
    };

    class Writer
    {
    public:
        Writer();
        ~Writer();

        void writeall(ALM *) const;
        void write_input_vars(const ALM *) const;
        void write_displacement_pattern(ALM *) const;

    private:
        void write_force_constants(ALM *) const;
        void write_misc_xml(ALM *) const;
        void write_hessian(ALM *) const;
        void write_in_QEformat(ALM *) const;
        void write_fc3_thirdorderpy_format(ALM *) const;
        std::string easyvizint(int) const;

        std::string double2string(double,
                                  int nprec = 15) const;
    };
}

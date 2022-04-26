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
#include "system.h"
#include "symmetry.h"
#include "cluster.h"
#include "patterndisp.h"
#include "fcs.h"
#include "constraint.h"
#include "optimize.h"
#include "files.h"

//#include "alm.h"

namespace ALM_NS {
class AtomProperty {
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

class SystemInfo {
public:
    double lattice_vector[3][3];
    std::vector<AtomProperty> atoms;
    size_t nat, natmin, ntran;
    size_t nspecies;

    SystemInfo() = default;;
};

class Writer {
public:
    Writer();

    ~Writer();

    void writeall(const System *system,
                  const Symmetry *symmetry,
                  const Cluster *cluster,
                  const Constraint *constraint,
                  const Fcs *fcs,
                  const Optimize *optimize,
                  const Files *files,
                  const int verbosity) const;

    void write_input_vars(const System *system,
                          const Symmetry *symmetry,
                          const Cluster *cluster,
                          const Displace *displace,
                          const Fcs *fcs,
                          const Constraint *constraint,
                          const Optimize *optimize,
                          const Files *files,
                          const std::string run_mode) const;

    void write_displacement_pattern(const Cluster *cluster,
                                    const Displace *displace,
                                    const std::string prefix,
                                    const int verbosity) const;

private:
    void write_force_constants(const Cluster *cluster,
                               const Fcs *fcs,
                               const Symmetry *symmetry,
                               const double *fcs_vals,
                               const int verbosity,
                               const std::string fname_save) const;

    void write_misc_xml(const System *system,
                        const Symmetry *symmetry,
                        const Cluster *cluster,
                        const Fcs *fcs,
                        const Constraint *constraint,
                        const Files *files,
                        const double *fcs_vals,
                        const int verbosity) const;

    void write_hessian(const System *system,
                       const Symmetry *symmetry,
                       const Fcs *fcs,
                       const int verbosity,
                       const std::string fname_out) const;

    void write_in_QEformat(const System *system,
                           const Symmetry *symmetry,
                           const Fcs *fcs,
                           const std::string fname_out) const;

    void write_fc3_thirdorderpy_format(const System *system,
                                       const Symmetry *symmetry,
                                       const Cluster *cluster,
                                       const Constraint *constraint,
                                       const Fcs *fcs,
                                       const std::string fname_out) const;

    std::string easyvizint(int) const;

    std::string double2string(double,
                              int nprec = 15) const;
};
}

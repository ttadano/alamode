/*
 strain_coupling_types.h

 Copyright (c) 2026 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include <array>
#include <cstddef>
#include <string>
#include <vector>
#include "strain_reference_cell.h"

// Plain data of the strain-coupling inputs of the SCPH/QHA cell relaxation,
// shared by the text parsers (strain_file_parsers.h), the HDF5 container
// reader (strain_coupling_io.h) and the consumers in DerivativeIFC, so that
// both input routes feed the same processing code. Eigen and the standard
// library only: no MPI, no HDF5, no System.

namespace PHON_NS
{
namespace strain_coupling
{
// Where the strain couplings and the elastic constants come from: the legacy
// text files below STRAIN_IFC_DIR (ifc_dir, with a trailing slash) or the
// HDF5 container named by STRAINFILE (file).
struct StrainSource
{
    std::string ifc_dir;
    std::string file;

    bool use_file() const
    {
        return !file.empty();
    }
};

// One block of the strain-force coupling: the forces (eV/Angstrom) in the
// cell strained along `mode` by `smag`, one row of three components per atom.
struct StrainForceBlock
{
    std::string mode;
    double smag{0.0};
    double weight{0.0};
    std::vector<double> forces; // natom_rows * 3, row-major
};

struct StrainForceSet
{
    std::vector<StrainForceBlock> blocks;
    std::size_t natom_rows{0};
    bool has_cell{false};               // the rows follow the atoms of `cell` ...
    strain_parsers::ReferenceCell cell; // ... otherwise those of the current primitive cell
    std::string origin;                 // for messages: the file (and group) the data came from
    std::string cell_description;       // for the log: where `cell` was declared
    bool trailing_data{false};          // text route: tokens remained after the last block
};

// One strained supercell of the strain-harmonic-IFC coupling. The force
// constants themselves are loaded separately (Fcs_phonon); `label` names the
// file relative to STRAIN_IFC_DIR or the group inside the container.
struct StrainHarmonicEntry
{
    std::string mode;
    double smag{0.0};
    double weight{0.0};
    std::string label;
    bool has_cell{false};
    strain_parsers::ReferenceCell supercell; // the strained supercell, when the source records it
};

struct StrainHarmonicSet
{
    std::vector<StrainHarmonicEntry> entries;
    std::string origin;
    bool trailing_data{false};
};

// The /Elastic group of the container: reference stress and elastic
// constants in GPa (the layout of the text files, i = 3*mu + nu).
struct ElasticSet
{
    bool has_stress{false};
    bool has_c2c3{false};
    std::array<double, 9> stress_gpa{};
    std::vector<double> soec_gpa; // 81
    std::vector<double> toec_gpa; // 729
    std::string source;
};
} // namespace strain_coupling
} // namespace PHON_NS

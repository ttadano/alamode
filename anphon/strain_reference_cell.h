/*
 strain_reference_cell.h

 Copyright (c) 2026 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include <Eigen/Core>
#include <cstddef>
#include <istream>
#include <string>
#include <vector>

// The optional "&reference_cell ... /" header of strain_force.in, and the
// mapping from the cell it declares to the primitive cell of the current run.
//
// strain_force.in stores one force row per atom of the cell the strained DFT
// calculations were done for. When the &cell field of anphon selects a
// different (nested) cell, the rows must be tiled or averaged onto its atoms;
// the header makes the file self-describing so that anphon can do that itself:
//
//   &reference_cell
//     <scale>                 bohr per unit of the vectors below (as in &cell)
//     a1x a1y a1z             one lattice vector per line
//     a2x a2y a2z
//     a3x a3y a3z
//     <natom>
//     <symbol> x y z          fractional coordinates, natom lines
//   /
//
// No comments are allowed anywhere: the reader is a plain token stream.
// This translation unit depends on Eigen and the standard library only (no
// MPI, no System), so that it can be unit-tested; errors are reported by
// throwing std::runtime_error, which the caller turns into exit().

namespace PHON_NS
{
namespace strain_parsers
{
struct ReferenceCell
{
    // lattice(i, j) : i-th Cartesian component of the j-th lattice vector, bohr
    Eigen::Matrix3d lattice{Eigen::Matrix3d::Identity()};
    std::vector<std::string> symbols;
    Eigen::MatrixXd x_fractional; // natom x 3

    std::size_t natom() const
    {
        return symbols.size();
    }

    Eigen::MatrixXd x_cartesian() const
    {
        return x_fractional * lattice.transpose();
    }
};

// True if `token` is the header sentinel "&reference_cell" (case-insensitive).
bool is_reference_cell_tag(const std::string &token);

// Parse the header body that follows the sentinel, up to and including the
// terminating "/". Throws std::runtime_error when it is malformed.
ReferenceCell parse_reference_cell(std::istream &fin, const char *filename);

// Which reference atoms feed each atom of the current primitive cell.
struct AtomMatch
{
    double ratio{1.0};                 // V_current / V_reference
    std::vector<std::vector<int>> src; // src[i]: reference atoms mapped onto current atom i
};

// Match the atoms of the current primitive cell (lattice vectors as columns
// in bohr, Cartesian positions in bohr, element symbols) against the
// reference cell, folding the difference vectors through the smaller of the
// two lattices. The cells must be nested: one has to be an integer supercell
// of the other. Throws std::runtime_error on an incommensurate reference
// cell, an inconsistent atom count, or an atom without a counterpart.
AtomMatch match_atoms(const ReferenceCell &ref, const Eigen::Matrix3d &lattice_cur, const Eigen::MatrixXd &xc_cur,
                      const std::vector<std::string> &symbols_cur, const char *filename);

// Map one block of force rows from the reference atoms (rows_in, natom_ref x 3,
// row-major) onto the current atoms (rows_out, natom_cur x 3): a copy when the
// current cell is the larger one, the average of the translation images
// otherwise. Returns the largest spread among the averaged images (0 when
// nothing is averaged) so that the caller can warn about it.
double expand_atom_rows(const AtomMatch &match, const std::vector<double> &rows_in, std::vector<double> &rows_out);
} // namespace strain_parsers
} // namespace PHON_NS

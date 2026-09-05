/*
 strain_reference_cell.cpp

 Copyright (c) 2026 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include "strain_reference_cell.h"
#include <Eigen/Dense> // inverse() and determinant() are defined in Eigen/LU
#include <algorithm>
#include <cctype>
#include <cmath>
#include <sstream>
#include <stdexcept>

namespace PHON_NS
{
namespace strain_parsers
{
namespace
{
std::string lowercase(std::string s)
{
    for (auto &c: s) {
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }
    return s;
}

std::string malformed(const char *filename)
{
    return std::string("Could not parse the &reference_cell header of ") + filename +
           ". The expected content is a scale factor,\n"
           " three lattice vectors (one per line, as in the &cell field), the number of atoms, that many\n"
           " 'element x y z' lines in fractional coordinates, and a closing '/'. Comments are not allowed.";
}

bool is_integer_matrix(const Eigen::Matrix3d &m, const double tol)
{
    return (m - m.array().round().matrix()).cwiseAbs().maxCoeff() < tol;
}

std::string atom_label(const std::string &symbol, const Eigen::RowVector3d &xf)
{
    std::ostringstream os;
    os << symbol << ", fractional " << xf(0) << " " << xf(1) << " " << xf(2);
    return os.str();
}

// Fractional coordinates of Cartesian positions in the given lattice, reduced
// to [0, 1) - used only to label atoms in messages.
Eigen::RowVector3d fractional_in(const Eigen::Matrix3d &lattice, const Eigen::RowVector3d &xc)
{
    Eigen::Vector3d f = lattice.inverse() * xc.transpose();
    for (int k = 0; k < 3; ++k) f(k) -= std::floor(f(k));
    return f.transpose();
}
} // namespace

bool is_reference_cell_tag(const std::string &token)
{
    return lowercase(token) == "&reference_cell";
}

ReferenceCell parse_reference_cell(std::istream &fin, const char *filename)
{
    ReferenceCell ref;
    double scale;
    if (!(fin >> scale)) throw std::runtime_error(malformed(filename));

    for (int i = 0; i < 3; ++i) { // i-th input line = i-th lattice vector = i-th column
        for (int j = 0; j < 3; ++j) {
            double v;
            if (!(fin >> v)) throw std::runtime_error(malformed(filename));
            ref.lattice(j, i) = scale * v;
        }
    }
    if (std::fabs(ref.lattice.determinant()) < 1.0e-8) {
        throw std::runtime_error(std::string("The lattice vectors in the &reference_cell header of ") + filename +
                                 " are linearly dependent.");
    }

    long natom;
    if (!(fin >> natom) || natom < 0) throw std::runtime_error(malformed(filename));

    ref.symbols.resize(static_cast<std::size_t>(natom));
    ref.x_fractional.resize(natom, 3);
    for (long i = 0; i < natom; ++i) {
        double x, y, z;
        if (!(fin >> ref.symbols[i] >> x >> y >> z)) throw std::runtime_error(malformed(filename));
        ref.x_fractional.row(i) << x, y, z;
    }

    std::string tok;
    if (!(fin >> tok) || tok != "/") {
        throw std::runtime_error(std::string("The &reference_cell header of ") + filename +
                                 " is not terminated by '/' (check the number of atoms declared).");
    }
    return ref;
}

AtomMatch match_atoms(const ReferenceCell &ref, const Eigen::Matrix3d &lattice_cur, const Eigen::MatrixXd &xc_cur,
                      const std::vector<std::string> &symbols_cur, const char *filename)
{
    const auto n_ref = static_cast<int>(ref.natom());
    const auto n_cur = static_cast<int>(symbols_cur.size());
    const std::string fname(filename);

    if (n_ref == 0) {
        throw std::runtime_error("The &reference_cell header of " + fname +
                                 " declares no atoms (natom = 0).\n strain_force.in stores one force row per atom, "
                                 "so the header must contain the atom list.");
    }

    // Lattice relation: the current cell must be an integer supercell of the
    // reference cell, or the other way round (nested cells only).
    AtomMatch match;
    bool cur_is_larger;
    const Eigen::Matrix3d m = ref.lattice.inverse() * lattice_cur;
    if (is_integer_matrix(m, 1.0e-5 * std::max(1.0, m.cwiseAbs().maxCoeff()))) {
        match.ratio = std::fabs(m.determinant());
        cur_is_larger = true;
    } else {
        const Eigen::Matrix3d m2 = lattice_cur.inverse() * ref.lattice;
        if (!is_integer_matrix(m2, 1.0e-5 * std::max(1.0, m2.cwiseAbs().maxCoeff()))) {
            throw std::runtime_error(
                "The reference cell declared in the &reference_cell header of " + fname +
                " is incommensurate\n with the primitive cell of this run: neither lattice is an integer "
                "supercell of the other\n (rotated settings and non-nested commensurate cells are not "
                "supported). Please check the header and the &cell field.");
        }
        match.ratio = 1.0 / std::fabs(m2.determinant());
        cur_is_larger = false;
    }
    const auto nimage = static_cast<int>(std::lround(cur_is_larger ? match.ratio : 1.0 / match.ratio));
    if (nimage < 1 || std::fabs((cur_is_larger ? match.ratio : 1.0 / match.ratio) - nimage) > 1.0e-4) {
        std::ostringstream os;
        os << "The volume ratio V(current) / V(reference) = " << match.ratio << " obtained from the &reference_cell "
           << "header of\n " << fname << " is neither an integer nor the inverse of one. Please check the header "
           << "and the &cell field.";
        throw std::runtime_error(os.str());
    }

    const auto expected_cur = cur_is_larger ? n_ref * nimage : n_ref / nimage;
    if ((!cur_is_larger && n_ref % nimage != 0) || expected_cur != n_cur) {
        std::ostringstream os;
        os << "Inconsistent atom counts: the &reference_cell header of " << fname << " declares " << n_ref
           << " atoms and\n V(current) / V(reference) = " << match.ratio << ", so the primitive cell of this run "
           << "should contain " << (cur_is_larger ? n_ref * nimage : (n_ref + nimage - 1) / nimage)
           << " atoms; it contains " << n_cur << ".\n Please check the header and the &cell field.";
        throw std::runtime_error(os.str());
    }

    // Fold through the smaller lattice and match by Cartesian distance, the
    // criterion Fcs_phonon::replicate_force_constant uses (1.0e-3 bohr).
    const Eigen::Matrix3d b = cur_is_larger ? ref.lattice : lattice_cur;
    const Eigen::Matrix3d binv = b.inverse();
    const Eigen::MatrixXd xc_ref = ref.x_cartesian();
    std::vector<std::string> sym_ref_lower(n_ref), sym_cur_lower(n_cur);
    for (int j = 0; j < n_ref; ++j) sym_ref_lower[j] = lowercase(ref.symbols[j]);
    for (int i = 0; i < n_cur; ++i) sym_cur_lower[i] = lowercase(symbols_cur[i]);

    match.src.assign(n_cur, {});
    std::vector<int> used(n_ref, 0);
    for (int i = 0; i < n_cur; ++i) {
        for (int j = 0; j < n_ref; ++j) {
            if (sym_cur_lower[i] != sym_ref_lower[j]) continue;
            Eigen::Vector3d d = binv * (xc_cur.row(i) - xc_ref.row(j)).transpose();
            for (int k = 0; k < 3; ++k) d(k) -= std::round(d(k));
            if ((b * d).norm() < 1.0e-3) {
                match.src[i].push_back(j);
                ++used[j];
            }
        }
    }

    const int want_per_cur = cur_is_larger ? 1 : nimage;
    for (int i = 0; i < n_cur; ++i) {
        const auto found = static_cast<int>(match.src[i].size());
        if (found == 0) {
            throw std::runtime_error("Atom " + std::to_string(i + 1) + " (" +
                                     atom_label(symbols_cur[i], fractional_in(lattice_cur, xc_cur.row(i))) +
                                     ") of the primitive cell of this run is not\n a translation image of any atom "
                                     "declared in the &reference_cell header of " +
                                     fname +
                                     ".\n Both cells must describe the same crystal in the same Cartesian "
                                     "frame. Please check the header and the &cell field.");
        }
        if (found != want_per_cur) {
            throw std::runtime_error("Atom " + std::to_string(i + 1) + " (" +
                                     atom_label(symbols_cur[i], fractional_in(lattice_cur, xc_cur.row(i))) +
                                     ") of the primitive cell of this run matches " + std::to_string(found) +
                                     " atoms declared in the\n &reference_cell header of " + fname + " (expected " +
                                     std::to_string(want_per_cur) +
                                     "). The atoms in the header must be inequivalent under the\n reference "
                                     "lattice translations.");
        }
    }
    const int want_per_ref = cur_is_larger ? nimage : 1;
    for (int j = 0; j < n_ref; ++j) {
        if (used[j] != want_per_ref) {
            throw std::runtime_error("Atom " + std::to_string(j + 1) + " (" +
                                     atom_label(ref.symbols[j], ref.x_fractional.row(j)) +
                                     ") declared in the &reference_cell header of " + fname + " is used " +
                                     std::to_string(used[j]) + " time(s) instead of " + std::to_string(want_per_ref) +
                                     " by the atoms of the primitive cell of this run.\n Both cells must describe the "
                                     "same crystal. Please check the header and the &cell field.");
        }
    }
    return match;
}

double expand_atom_rows(const AtomMatch &match, const std::vector<double> &rows_in, std::vector<double> &rows_out)
{
    const auto n_cur = match.src.size();
    rows_out.assign(n_cur * 3, 0.0);
    double spread = 0.0;
    for (std::size_t i = 0; i < n_cur; ++i) {
        const auto &src = match.src[i];
        if (src.empty()) continue;
        for (int k = 0; k < 3; ++k) {
            double mean = 0.0;
            for (const auto j: src) mean += rows_in[static_cast<std::size_t>(j) * 3 + k];
            mean /= static_cast<double>(src.size());
            rows_out[i * 3 + k] = mean;
            for (const auto j: src) {
                spread = std::max(spread, std::fabs(rows_in[static_cast<std::size_t>(j) * 3 + k] - mean));
            }
        }
    }
    return spread;
}
} // namespace strain_parsers
} // namespace PHON_NS

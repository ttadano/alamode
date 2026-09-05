/*
 test_strain_parsers.cpp

 Unit test for anphon/strain_file_parsers.h: the token parsers of
 C1_array.in and elastic_constants.in. Covers the optional unit token
 (GPa / Ry, case-insensitive), legacy files without a token (the token after
 the label is the first value), unknown units, short reads, mismatched units
 between the SOEC and TOEC sections, trailing data, and the fact that the
 label spelling is not enforced.

 Built by the anphon CMake project as `test_strain_parsers`. Exits 0 on
 success and prints the first failing check otherwise.
*/

#include <Eigen/Core>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include "strain_file_parsers.h"
#include "strain_reference_cell.h"

using namespace PHON_NS::strain_parsers;

namespace
{
int nfail = 0;

void check(const bool ok, const char *what)
{
    if (!ok) {
        std::printf("FAIL: %s\n", what);
        ++nfail;
    }
}

bool near(const double a, const double b)
{
    return std::fabs(a - b) < 1.0e-12 * std::fmax(1.0, std::fabs(b));
}

// A section body: values start at `first` and increase by 1.
std::string body(const int n, const double first = 0.5)
{
    std::string s;
    for (int i = 0; i < n; ++i) s += std::to_string(first + i) + "\n";
    return s;
}

template <typename F>
bool throws(F f)
{
    try {
        f();
    } catch (const std::runtime_error &) {
        return true;
    }
    return false;
}
} // namespace

int main()
{
    // ---- parse_double_strict --------------------------------------------
    double v = 0.0;
    check(parse_double_strict("1e-3", v) && near(v, 1.0e-3), "strict: 1e-3");
    check(parse_double_strict("-2.0", v) && near(v, -2.0), "strict: -2.0");
    check(parse_double_strict("7", v) && near(v, 7.0), "strict: integer");
    check(!parse_double_strict("1.5x", v), "strict: trailing garbage rejected");
    check(!parse_double_strict("", v), "strict: empty rejected");
    check(!parse_double_strict("GPa", v), "strict: unit token rejected");
    check(!parse_double_strict("Ry", v), "strict: Ry rejected");

    // ---- C1_array.in ---------------------------------------------------
    {
        std::istringstream in("C1\n" + body(9));
        const auto d = parse_c1_array(in);
        check(d.unit == ElasticUnit::Legacy, "C1 legacy: unit");
        check(near(d.values[0], 0.5) && near(d.values[8], 8.5), "C1 legacy: first token is value 0, last is value 8");
        check(!has_trailing_data(in), "C1 legacy: no trailing data");
    }
    {
        std::istringstream in("C1 GPa\n" + body(9, -3.0));
        const auto d = parse_c1_array(in);
        check(d.unit == ElasticUnit::GPa, "C1 GPa: unit");
        check(near(d.values[0], -3.0) && near(d.values[8], 5.0), "C1 GPa: values");
    }
    {
        std::istringstream in("c1   ry\n" + body(9));
        check(parse_c1_array(in).unit == ElasticUnit::Ry, "C1 Ry: case-insensitive unit");
    }
    {
        std::istringstream in("whatever_label\n" + body(9));
        check(parse_c1_array(in).unit == ElasticUnit::Legacy, "C1: label spelling is not enforced");
    }
    {
        std::istringstream in("C1 MPa\n" + body(9));
        check(throws([&] { parse_c1_array(in); }), "C1: unknown unit throws");
    }
    {
        std::istringstream in("C1\n" + body(8));
        check(throws([&] { parse_c1_array(in); }), "C1 legacy: 8 values throws");
    }
    {
        std::istringstream in("C1 GPa\n" + body(8));
        check(throws([&] { parse_c1_array(in); }), "C1 GPa: 8 values throws");
    }
    {
        std::istringstream in("C1\n" + body(9) + "extra\n");
        parse_c1_array(in);
        check(has_trailing_data(in), "C1: trailing token detected");
    }
    {
        std::istringstream in("");
        check(throws([&] { parse_c1_array(in); }), "C1: empty stream throws");
    }

    // ---- elastic_constants.in ------------------------------------------
    {
        std::istringstream in("SOEC\n" + body(81) + "TOEC\n" + body(729, 1000.0));
        const auto d = parse_elastic_constants(in);
        check(d.unit_soec == ElasticUnit::Legacy && d.unit_toec == ElasticUnit::Legacy, "elastic legacy: units");
        check(d.c2.size() == 81 && d.c3.size() == 729, "elastic legacy: sizes");
        check(near(d.c2[0], 0.5) && near(d.c2[80], 80.5), "elastic legacy: SOEC values");
        check(near(d.c3[0], 1000.0) && near(d.c3[728], 1728.0), "elastic legacy: TOEC values");
        check(!has_trailing_data(in), "elastic legacy: no trailing data");
    }
    {
        std::istringstream in("SOEC GPa\n" + body(81) + "TOEC gpa\n" + body(729));
        const auto d = parse_elastic_constants(in);
        check(d.unit_soec == ElasticUnit::GPa && d.unit_toec == ElasticUnit::GPa,
              "elastic GPa: units (case-insensitive)");
    }
    {
        // The unit token is per section: a mixed file parses, and the caller
        // converts each section according to its own unit.
        std::istringstream in("SOEC GPa\n" + body(81) + "TOEC\n" + body(729, 1000.0));
        const auto d = parse_elastic_constants(in);
        check(d.unit_soec == ElasticUnit::GPa && d.unit_toec == ElasticUnit::Legacy, "elastic: per-section units");
        check(near(d.c2[0], 0.5) && near(d.c3[0], 1000.0), "elastic: per-section values");
    }
    {
        std::istringstream in("SOEC\n" + body(81) + "TOEC Ry\n" + body(729));
        const auto d = parse_elastic_constants(in);
        check(d.unit_soec == ElasticUnit::Legacy && d.unit_toec == ElasticUnit::Ry, "elastic: legacy SOEC + Ry TOEC");
    }
    {
        std::istringstream in("SOEC\n" + body(81) + "TOEC\n" + body(729) + "extra\n");
        parse_elastic_constants(in);
        check(has_trailing_data(in), "elastic: trailing token detected");
    }
    {
        std::istringstream in("SOEC\n" + body(81) + "TOEC\n" + body(728));
        check(throws([&] { parse_elastic_constants(in); }), "elastic: short TOEC throws");
    }
    {
        std::istringstream in("SOEC\n" + body(80));
        check(throws([&] { parse_elastic_constants(in); }), "elastic: short SOEC throws");
    }
    {
        std::istringstream in("SOEC\n" + body(81));
        check(throws([&] { parse_elastic_constants(in); }), "elastic: missing TOEC section throws");
    }
    {
        std::istringstream in("SOEC\n" + body(81) + "TOEC\n" + body(729) + "\n  \n");
        parse_elastic_constants(in);
        check(!has_trailing_data(in), "elastic: trailing whitespace is not data");
    }

    // ---- &reference_cell header of strain_force.in ---------------------
    check(is_reference_cell_tag("&reference_cell") && is_reference_cell_tag("&REFERENCE_CELL"),
          "header: sentinel is case-insensitive");
    check(!is_reference_cell_tag("xx") && !is_reference_cell_tag("&cell"), "header: other tokens are not the sentinel");
    const std::string hdr_body = "  2.0\n 5 0 0\n 0 5 0\n 0 0 5\n 2\n A 0 0 0\n B 0.5 0.5 0.5\n/\n";
    {
        std::istringstream in(hdr_body + "xx 0.005 1\n");
        const auto ref = parse_reference_cell(in, "strain_force.in");
        check(ref.natom() == 2 && ref.symbols[0] == "A" && ref.symbols[1] == "B", "header: atoms");
        check(near(ref.lattice(0, 0), 10.0) && near(ref.lattice(1, 1), 10.0) && near(ref.lattice(0, 1), 0.0),
              "header: lattice = scale * vectors, vectors as columns");
        check(near(ref.x_fractional(1, 2), 0.5), "header: fractional coordinates");
        check(near(ref.x_cartesian()(1, 0), 5.0), "header: Cartesian positions");
        std::string tok;
        in >> tok;
        check(tok == "xx", "header: the stream continues with the first block");
    }
    {
        std::istringstream in("  2.0\n 5 0 0\n 0 5 0\n 0 0 5\n 2\n A 0 0 0\n B 0.5 0.5 0.5\n");
        check(throws([&] { parse_reference_cell(in, "f"); }), "header: missing '/' throws");
    }
    {
        std::istringstream in("  2.0\n 5 0 0\n 0 5 0\n 0 0 5\n 3\n A 0 0 0\n B 0.5 0.5 0.5\n/\n");
        check(throws([&] { parse_reference_cell(in, "f"); }), "header: natom larger than the list throws");
    }
    {
        std::istringstream in("  2.0\n 5 0 0\n 0 5 0\n 0 0 5\n 1\n A 0 0 0\n B 0.5 0.5 0.5\n/\n");
        check(throws([&] { parse_reference_cell(in, "f"); }), "header: natom smaller than the list throws");
    }
    {
        std::istringstream in("  2.0\n 5 0 0\n 0 5 0\n 0 0 0\n 2\n A 0 0 0\n B 0.5 0.5 0.5\n/\n");
        check(throws([&] { parse_reference_cell(in, "f"); }), "header: singular lattice throws");
    }
    {
        std::istringstream in("  2.0\n 5 0 x\n");
        check(throws([&] { parse_reference_cell(in, "f"); }), "header: non-numeric lattice throws");
    }

    // ---- atom matching between nested cells ----------------------------
    // Reference: simple cubic, a = 10 bohr, A at the origin and B at the body
    // center. Current: the same crystal in a cell doubled along z, with the
    // atoms in anphon's first-occurrence order.
    ReferenceCell ref;
    ref.lattice = Eigen::Matrix3d::Identity() * 10.0;
    ref.symbols = {"A", "B"};
    ref.x_fractional.resize(2, 3);
    ref.x_fractional << 0.0, 0.0, 0.0, 0.5, 0.5, 0.5;

    Eigen::Matrix3d lat2 = Eigen::Matrix3d::Identity() * 10.0;
    lat2(2, 2) = 20.0;
    Eigen::MatrixXd xf2(4, 3);
    xf2 << 0.0, 0.0, 0.0, 0.5, 0.5, 0.25, 0.0, 0.0, 0.5, 0.5, 0.5, 0.75;
    const Eigen::MatrixXd xc2 = xf2 * lat2.transpose();
    const std::vector<std::string> sym2 = {"A", "B", "A", "B"};
    {
        const auto m = match_atoms(ref, lat2, xc2, sym2, "f");
        check(near(m.ratio, 2.0), "match: ratio 2 for the doubled cell");
        check(m.src.size() == 4 && m.src[0] == std::vector<int>{0} && m.src[1] == std::vector<int>{1} &&
                  m.src[2] == std::vector<int>{0} && m.src[3] == std::vector<int>{1},
              "match: each current atom maps to one reference atom (tiling)");
        std::vector<double> rows_in = {1, 2, 3, 4, 5, 6}, rows_out;
        const auto spread = expand_atom_rows(m, rows_in, rows_out);
        check(rows_out.size() == 12 && near(rows_out[6], 1.0) && near(rows_out[11], 6.0) && near(spread, 0.0),
              "expand: rows are copied onto the translation images");
    }
    {
        // Element symbols are compared case-insensitively.
        const std::vector<std::string> sym_lc = {"a", "b", "a", "b"};
        check(near(match_atoms(ref, lat2, xc2, sym_lc, "f").ratio, 2.0), "match: symbols are case-insensitive");
    }
    {
        // Non-trivial atom order in the current cell.
        Eigen::MatrixXd xfp(4, 3);
        xfp << 0.5, 0.5, 0.25, 0.0, 0.0, 0.5, 0.5, 0.5, 0.75, 0.0, 0.0, 0.0;
        const std::vector<std::string> symp = {"B", "A", "B", "A"};
        const auto m = match_atoms(ref, lat2, xfp * lat2.transpose(), symp, "f");
        check(m.src[0] == std::vector<int>{1} && m.src[1] == std::vector<int>{0} && m.src[3] == std::vector<int>{0},
              "match: permuted current atoms are mapped correctly");
    }
    {
        // Reverse direction: the reference is the larger cell, rows are averaged.
        ReferenceCell big;
        big.lattice = lat2;
        big.symbols = sym2;
        big.x_fractional = xf2;
        const Eigen::MatrixXd xc_small = ref.x_fractional * ref.lattice.transpose();
        const auto m = match_atoms(big, ref.lattice, xc_small, ref.symbols, "f");
        check(near(m.ratio, 0.5), "match: ratio 1/2 for the smaller current cell");
        check(m.src.size() == 2 && m.src[0] == std::vector<int>{0, 2} && m.src[1] == std::vector<int>{1, 3},
              "match: each current atom collects its translation images");
        std::vector<double> rows_in = {1, 0, 0, 10, 0, 0, 3, 0, 0, 10, 0, 0}, rows_out;
        const auto spread = expand_atom_rows(m, rows_in, rows_out);
        check(rows_out.size() == 6 && near(rows_out[0], 2.0) && near(rows_out[3], 10.0) && near(spread, 1.0),
              "expand: images are averaged and the spread is reported");
    }
    {
        Eigen::Matrix3d lat_bad = lat2;
        lat_bad(0, 0) = 15.0; // 1.5 x: not an integer supercell in either direction
        check(throws([&] { match_atoms(ref, lat_bad, xf2 * lat_bad.transpose(), sym2, "f"); }),
              "match: incommensurate cell throws");
    }
    {
        check(throws([&] { match_atoms(ref, lat2, xc2.topRows(3), {"A", "B", "A"}, "f"); }),
              "match: inconsistent atom count throws");
    }
    {
        Eigen::MatrixXd xf_off = xf2;
        xf_off(2, 0) += 0.05;
        check(throws([&] { match_atoms(ref, lat2, xf_off * lat2.transpose(), sym2, "f"); }),
              "match: displaced (unmatched) atom throws");
    }
    {
        check(throws([&] { match_atoms(ref, lat2, xc2, {"A", "B", "B", "A"}, "f"); }), "match: wrong species throws");
    }
    {
        ReferenceCell empty;
        empty.lattice = ref.lattice;
        check(throws([&] { match_atoms(empty, lat2, xc2, sym2, "f"); }), "match: header without atoms throws");
    }

    if (nfail == 0) {
        std::printf("test_strain_parsers: all checks passed\n");
        return EXIT_SUCCESS;
    }
    std::printf("test_strain_parsers: %d check(s) failed\n", nfail);
    return EXIT_FAILURE;
}

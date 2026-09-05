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

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <stdexcept>
#include <string>
#include "strain_file_parsers.h"

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

    if (nfail == 0) {
        std::printf("test_strain_parsers: all checks passed\n");
        return EXIT_SUCCESS;
    }
    std::printf("test_strain_parsers: %d check(s) failed\n", nfail);
    return EXIT_FAILURE;
}

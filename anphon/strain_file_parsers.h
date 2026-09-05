/*
 strain_file_parsers.h

 Copyright (c) 2026 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include <array>
#include <cstddef>
#include <istream>
#include <string>
#include <vector>

// Token-stream parsers of the elastic-constant input files of the SCPH/QHA
// structural optimization (C1_array.in and elastic_constants.in).
//
// This translation unit has no dependency on the rest of anphon (no MPI, no
// System) so that the parsing can be unit-tested; errors are reported by
// throwing std::runtime_error, which ElasticTensor turns into exit().
//
// File format: each section is a label token, an OPTIONAL unit token, and the
// values.  The label spelling is not checked (it never was).  Without a unit
// token, or with "Ry", the values are the legacy V*C in Ry for one specific
// cell; with "GPa" they are intensive densities that anphon multiplies by the
// volume of the current primitive cell, which makes the file independent of
// the &cell field.

namespace PHON_NS
{
namespace strain_parsers
{
enum class ElasticUnit
{
    Legacy,
    Ry,
    GPa
};

const char *unit_name(ElasticUnit unit);

// True if the whole token is a floating-point number (no trailing characters).
bool parse_double_strict(const std::string &token, double &value);

// Read "<label> [unit] value*nvalues" from fin. The token after the label is
// the unit if it is not a number, otherwise it is the first value.
// Throws std::runtime_error on a short read or an unknown unit.
ElasticUnit read_section(std::istream &fin, const char *filename, const char *what, std::size_t nvalues,
                         double *values);

struct C1Data
{
    ElasticUnit unit{ElasticUnit::Legacy};
    std::array<double, 9> values{};
};

struct ElasticData
{
    ElasticUnit unit_soec{ElasticUnit::Legacy};
    ElasticUnit unit_toec{ElasticUnit::Legacy};
    std::vector<double> c2; // 81 values, row-major i = 3*mu + nu
    std::vector<double> c3; // 729 values
};

// C1_array.in: one section of 9 values (the stress tensor).
C1Data parse_c1_array(std::istream &fin);

// elastic_constants.in: SOEC (81 values) followed by TOEC (729 values), each
// section with its own optional unit token.
ElasticData parse_elastic_constants(std::istream &fin);

// True if any token remains after the last value (the caller reports it).
bool has_trailing_data(std::istream &fin);
} // namespace strain_parsers
} // namespace PHON_NS

/*
 strain_file_parsers.cpp

 Copyright (c) 2026 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include "strain_file_parsers.h"
#include <cctype>
#include <cerrno>
#include <cstdlib>
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

std::runtime_error short_read(const char *filename, const char *what, const std::size_t nvalues)
{
    return std::runtime_error(std::string(filename) + " ended before all " + std::to_string(nvalues) +
                              " values of the " + what + " section were read.");
}
} // namespace

const char *unit_name(const ElasticUnit unit)
{
    switch (unit) {
    case ElasticUnit::GPa:
        return "GPa";
    case ElasticUnit::Ry:
        return "Ry";
    default:
        return "legacy (V*C in Ry, no unit token)";
    }
}

bool parse_double_strict(const std::string &token, double &value)
{
    if (token.empty()) return false;
    const char *begin = token.c_str();
    char *end = nullptr;
    errno = 0;
    const double v = std::strtod(begin, &end);
    if (end != begin + token.size() || errno == ERANGE) return false;
    value = v;
    return true;
}

ElasticUnit read_section(std::istream &fin, const char *filename, const char *what, const std::size_t nvalues,
                         double *values)
{
    std::string tok;

    // The section label. Its spelling is deliberately not checked: the readers
    // never did, and rejecting it now would break files that work today.
    if (!(fin >> tok)) {
        throw std::runtime_error(std::string(filename) + " ended before the " + what + " section.");
    }
    if (!(fin >> tok)) throw short_read(filename, what, nvalues);

    auto unit = ElasticUnit::Legacy;
    std::size_t i = 0;
    double v;
    if (parse_double_strict(tok, v)) {
        // Legacy file: the token after the label is already the first value.
        values[i++] = v;
    } else {
        const auto low = lowercase(tok);
        if (low == "gpa") {
            unit = ElasticUnit::GPa;
        } else if (low == "ry") {
            unit = ElasticUnit::Ry;
        } else {
            throw std::runtime_error("Unknown unit '" + tok + "' after the " + what + " label of " + filename +
                                     ".\n Accepted units: GPa (cell-independent) and Ry (V*C per cell).");
        }
    }
    for (; i < nvalues; ++i) {
        if (!(fin >> values[i])) throw short_read(filename, what, nvalues);
    }
    return unit;
}

C1Data parse_c1_array(std::istream &fin)
{
    C1Data data;
    data.unit = read_section(fin, "C1_array.in", "C1", 9, data.values.data());
    return data;
}

ElasticData parse_elastic_constants(std::istream &fin)
{
    ElasticData data;
    data.c2.assign(81, 0.0);
    data.c3.assign(729, 0.0);
    data.unit_soec = read_section(fin, "elastic_constants.in", "SOEC", 81, data.c2.data());
    data.unit_toec = read_section(fin, "elastic_constants.in", "TOEC", 729, data.c3.data());
    return data;
}

bool has_trailing_data(std::istream &fin)
{
    std::string tok;
    return static_cast<bool>(fin >> tok);
}
} // namespace strain_parsers
} // namespace PHON_NS

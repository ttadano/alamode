/*
 stage_timer.h

 Wall-clock stage timers for the SCPH/QHA setup and postprocess, printed as
 "  [timer] <label> <seconds> sec." on rank 0 when VERBOSITY >= 2.

 Copyright (c) 2026 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include <chrono>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

namespace PHON_NS
{
inline double stage_clock()
{
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

// newline_first: start on a fresh line when the caller left a progress message open ("... ").
inline void print_stage_line(const std::string &label, const double seconds, const int my_rank,
                             const unsigned int verbosity, const bool newline_first = false)
{
    if (my_rank != 0 || verbosity < 2) return;
    // format in a local stream so that the precision does not leak into std::cout
    std::ostringstream line;
    if (newline_first) line << '\n';
    line << "  [timer] " << std::left << std::setw(40) << label << std::right << std::fixed << std::setprecision(3)
         << seconds << " sec.\n";
    std::cout << line.str();
}
} // namespace PHON_NS

/*
 error.h

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#ifdef _WIN32
#include <mpi.h>
#else

#include "mpi.h"

#endif

#include <cstdlib>
#include <iostream>

namespace PHON_NS
{
inline void warn(const char *file, const char *message)
{
    std::cout << '\n' << " WARNING in " << file << "  MESSAGE: " << message << '\n';
}


// Terminate the whole MPI job, not only the calling process. A failure that
// is detected on one rank only (a file that cannot be opened, an exception
// while reading it, ...) would otherwise leave the other ranks waiting in the
// next collective forever. MPI_Abort is skipped when MPI is not (or no longer)
// initialized so that the helpers stay usable in serial tools.
inline void abort_all_ranks()
{
    std::cout.flush();
    int initialized = 0, finalized = 0;
    MPI_Initialized(&initialized);
    MPI_Finalized(&finalized);
    if (initialized && !finalized) MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    std::exit(EXIT_FAILURE);
}

inline void exit(const char *file, const char *message)
{
    std::cout << '\n' << " ERROR in " << file << "  MESSAGE: " << message << '\n';
    abort_all_ranks();
}

template <typename T>
void exit(const char *file, const char *message, const T info)
{
    std::cout << '\n' << " ERROR in " << file << "  MESSAGE: " << message << info << '\n';
    abort_all_ranks();
}

inline void exit(const char *file, const char *message, const char *info)
{
    std::cout << '\n' << " ERROR in " << file << "  MESSAGE: " << message << info << '\n';
    abort_all_ranks();
}

inline void exitall(const char *file, const char *message)
{
    MPI_Finalize();
    std::cout << '\n' << "ERROR in " << file << "  MESSAGE: " << message << '\n';
    std::exit(EXIT_FAILURE);
}
} // namespace PHON_NS

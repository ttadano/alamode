/*
 error.h

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include <iostream>
#include <cstdlib>

namespace PHON_NS {
inline void warn(const char *file,
                 const char *message)
{
    std::cout << std::endl << " WARNING in " << file << "  MESSAGE: " << message << std::endl;
}


inline void exit(const char *file,
                 const char *message)
{
    std::cout << std::endl << " ERROR in " << file << "  MESSAGE: " << message << std::endl;
    std::exit(EXIT_FAILURE);
}

template<typename T>
void exit(const char *file,
          const char *message,
          const T info)
{
    std::cout << std::endl << " ERROR in " << file << "  MESSAGE: " << message << info << std::endl;
    std::exit(EXIT_FAILURE);
}

inline void exit(const char *file,
                 const char *message,
                 const char *info)
{
    std::cout << std::endl << " ERROR in " << file << "  MESSAGE: " << message << info << std::endl;
    std::exit(EXIT_FAILURE);
}

inline void exitall(const char *file,
                    const char *message)
{
    MPI_Finalize();
    std::cout << std::endl << "ERROR in " << file << "  MESSAGE: " << message << std::endl;
    std::exit(EXIT_FAILURE);
}
}


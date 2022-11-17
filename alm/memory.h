/*
 memory.h

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include <iostream>

// memsize calculator

namespace ALM_NS {
inline size_t memsize_in_MB(const size_t size_of_one,
                            const size_t n1)
{
    const auto n = n1 * size_of_one;
    return n / 1000000;
}

inline size_t memsize_in_MB(const size_t size_of_one,
                            const size_t n1,
                            const size_t n2)
{
    const auto n = n1 * n2 * size_of_one;
    return n / 1000000;
}

inline size_t memsize_in_MB(const size_t size_of_one,
                            const size_t n1,
                            const size_t n2,
                            const size_t n3)
{
    const auto n = n1 * n2 * n3 * size_of_one;
    return n / 1000000;
}

inline size_t memsize_in_MB(const size_t size_of_one,
                            const size_t n1,
                            const size_t n2,
                            const size_t n3,
                            const size_t n4)
{
    const auto n = n1 * n2 * n3 * n4 * size_of_one;
    return n / 1000000;
}

// Declaration and definition must be located in the same file for template functions.

/* allocator */

template<typename T>
T *allocate(T *&arr,
            const size_t n1)
{
    try {
        arr = new T[n1];
    }
    catch (std::bad_alloc &ba) {
        std::cout << " Caught an exception when trying to allocate 1-dimensional array" << std::endl;
        std::cout << " " << ba.what() << " : Array shape = " << n1 << std::endl;
        std::cout << " " << ba.what() << " : Array size (MB) = " << memsize_in_MB(sizeof(T), n1) << std::endl;
        std::exit(EXIT_FAILURE);
    }
    return arr;
}

template<typename T>
T **allocate(T **&arr,
             const size_t n1,
             const size_t n2)
{
    try {
        arr = new T *[n1];
        arr[0] = new T[n1 * n2];
        for (size_t i = 1; i < n1; ++i) {
            arr[i] = arr[0] + i * n2;
        }
    }
    catch (std::bad_alloc &ba) {
        std::cout << " Caught an exception when trying to allocate 2-dimensional array" << std::endl;
        std::cout << " " << ba.what() << " : Array shape = " << n1 << "x" << n2 << std::endl;
        std::cout << " " << ba.what() << " : Array size (MB) = " << memsize_in_MB(sizeof(T), n1, n2) << std::endl;
        std::exit(EXIT_FAILURE);
    }
    return arr;
}

template<typename T>
T ***allocate(T ***&arr,
              const size_t n1,
              const size_t n2,
              const size_t n3)
{
    try {
        arr = new T **[n1];
        arr[0] = new T *[n1 * n2];
        arr[0][0] = new T[n1 * n2 * n3];
        for (size_t i = 0; i < n1; ++i) {
            arr[i] = arr[0] + i * n2;
            for (size_t j = 0; j < n2; ++j) {
                arr[i][j] = arr[0][0] + i * n2 * n3 + j * n3;
            }
        }
    }
    catch (std::bad_alloc &ba) {
        std::cout << " Caught an exception when trying to allocate 3-dimensional array" << std::endl;
        std::cout << " " << ba.what() << " : Array shape = " << n1 << "x" << n2 << "x" << n3 << std::endl;
        std::cout << " " << ba.what() << " : Array size (MB) = " << memsize_in_MB(sizeof(T), n1, n2, n3) << std::
        endl;
        std::exit(EXIT_FAILURE);
    }
    return arr;
}

template<typename T>
T ****allocate(T ****&arr,
               const size_t n1,
               const size_t n2,
               const size_t n3,
               const size_t n4)
{
    try {
        arr = new T ***[n1];
        arr[0] = new T **[n1 * n2];
        arr[0][0] = new T *[n1 * n2 * n3];
        arr[0][0][0] = new T[n1 * n2 * n3 * n4];

        for (size_t i = 0; i < n1; ++i) {
            arr[i] = arr[0] + i * n2;
            for (size_t j = 0; j < n2; ++j) {
                arr[i][j] = arr[0][0] + i * n2 * n3 + j * n3;
                for (size_t k = 0; k < n3; ++k) {
                    arr[i][j][k] = arr[0][0][0] + i * n2 * n3 * n4 + j * n3 * n4 + k * n4;
                }
            }
        }
    }
    catch (std::bad_alloc &ba) {
        std::cout << " Caught an exception when trying to allocate 4-dimensional array" << std::endl;
        std::cout << " " << ba.what() << " : Array shape = " << n1 << "x" << n2 << "x" << n3 << "x" << n4 << std::
        endl;
        std::cout << " " << ba.what() << " : Array size (MB) = " << memsize_in_MB(sizeof(T), n1, n2, n3, n4) << std
        ::
        endl;
        std::exit(EXIT_FAILURE);
    }
    return arr;
}


/* deallocator */


template<typename T>
void deallocate(T *&arr)
{
    delete[] arr;
    arr = nullptr;
}

template<typename T>
void deallocate(T **&arr)
{
    delete[] arr[0];
    delete[] arr;
    arr = nullptr;
}

template<typename T>
void deallocate(T ***&arr)
{
    delete[] arr[0][0];
    delete[] arr[0];
    delete[] arr;
    arr = nullptr;
}

template<typename T>
void deallocate(T ****&arr)
{
    delete[] arr[0][0][0];
    delete[] arr[0][0];
    delete[] arr[0];
    delete[] arr;
    arr = nullptr;
}
}

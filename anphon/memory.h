/*
 memory.h

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include <iostream>
#include <new>
#include <cstdlib>
#include "pointers.h"

namespace PHON_NS
{
    class Memory : protected Pointers
    {
    public:
        Memory(class PHON *);

        ~Memory();

        // allocator 

        template <typename T>
        T* allocate(T *&arr,
                    const unsigned int n1)
        {
#ifdef _SX
            arr = new T [n1];
#else
            try {
                arr = new T [n1];
            }
            catch (std::bad_alloc &ba) {
                std::cout << " Caught an exception when trying to allocate 1-dimensional array" << std::endl;
                std::cout << " " << ba.what() << " : Array shape = " << n1 << std::endl;
                std::cout << " " << ba.what() << " : Array size (MB) = " << memsize_in_MB(sizeof(T), n1) << std::endl;
                exit(EXIT_FAILURE);
            }
#endif
            return arr;
        }

        template <typename T>
        T** allocate(T **&arr,
                     const unsigned int n1,
                     const unsigned int n2)
        {
#ifdef _SX
            arr = new T *[n1];
            arr[0] = new T [n1 * n2];
            for (unsigned int i = 1; i < n1; ++i){
                arr[i] = arr[0] + i * n2;
            }
#else
            try {
                arr = new T *[n1];
                arr[0] = new T [n1 * n2];
                for (unsigned int i = 1; i < n1; ++i) {
                    arr[i] = arr[0] + i * n2;
                }
            }
            catch (std::bad_alloc &ba) {
                std::cout << " Caught an exception when trying to allocate 2-dimensional array" << std::endl;
                std::cout << " " << ba.what() << " : Array shape = " << n1 << "x" << n2 << std::endl;
                std::cout << " " << ba.what() << " : Array size (MB) = " << memsize_in_MB(sizeof(T), n1, n2) << std::
                    endl;
                exit(EXIT_FAILURE);
            }
#endif
            return arr;
        }

        template <typename T>
        T*** allocate(T ***&arr,
                      const unsigned int n1,
                      const unsigned int n2,
                      const unsigned int n3)
        {
#ifdef _SX
            arr = new T **[n1];
            arr[0] = new T *[n1 * n2];
            arr[0][0] = new T [n1 * n2 * n3];
            for (unsigned int i = 0; i < n1; ++i){
                arr[i] = arr[0] + i * n2;
                for (unsigned int j = 0; j < n2; ++j){
                    arr[i][j] = arr[0][0] + i * n2 * n3 + j * n3;
                }
            }
#else
            try {
                arr = new T **[n1];
                arr[0] = new T *[n1 * n2];
                arr[0][0] = new T [n1 * n2 * n3];
                for (unsigned int i = 0; i < n1; ++i) {
                    arr[i] = arr[0] + i * n2;
                    for (unsigned int j = 0; j < n2; ++j) {
                        arr[i][j] = arr[0][0] + i * n2 * n3 + j * n3;
                    }
                }
            }
            catch (std::bad_alloc &ba) {
                std::cout << " Caught an exception when trying to allocate 3-dimensional array" << std::endl;
                std::cout << " " << ba.what() << " : Array shape = " << n1 << "x" << n2 << "x" << n3 << std::endl;
                std::cout << " " << ba.what() << " : Array size (MB) = " << memsize_in_MB(sizeof(T), n1, n2, n3) << std
                    ::endl;
                exit(EXIT_FAILURE);
            }
#endif
            return arr;
        }

        template <typename T>
        T**** allocate(T ****&arr,
                       const unsigned int n1,
                       const unsigned int n2,
                       const unsigned int n3,
                       const unsigned int n4)
        {
#ifdef _SX
            arr = new T ***[n1];
            arr[0] = new T **[n1 * n2];
            arr[0][0] = new T *[n1 * n2 * n3];
            arr[0][0][0] = new T [n1 * n2 * n3 * n4];

            for (unsigned int i = 0; i < n1; ++i){
                arr[i] = arr[0] + i * n2;
                for (unsigned int j = 0; j < n2; ++j){
                    arr[i][j] = arr[0][0] + i * n2 * n3 + j * n3;
                    for (unsigned int k = 0; k < n3; ++k){
                        arr[i][j][k] = arr[0][0][0] + i * n2 * n3 * n4 + j * n3 * n4 + k * n4;
                    }
                }
            }
#else
            try {
                arr = new T ***[n1];
                arr[0] = new T **[n1 * n2];
                arr[0][0] = new T *[n1 * n2 * n3];
                arr[0][0][0] = new T [n1 * n2 * n3 * n4];

                for (unsigned int i = 0; i < n1; ++i) {
                    arr[i] = arr[0] + i * n2;
                    for (unsigned int j = 0; j < n2; ++j) {
                        arr[i][j] = arr[0][0] + i * n2 * n3 + j * n3;
                        for (unsigned int k = 0; k < n3; ++k) {
                            arr[i][j][k] = arr[0][0][0] + i * n2 * n3 * n4 + j * n3 * n4 + k * n4;
                        }
                    }
                }
            }
            catch (std::bad_alloc &ba) {
                std::cout << " Caught an exception when trying to allocate 4-dimensional array" << std::endl;
                std::cout << " " << ba.what() << " : Array shape = " << n1 << "x" << n2 << "x" << n3 << "x" << n4 << std
                    ::endl;
                std::cout << " " << ba.what() << " : Array size (MB) = " << memsize_in_MB(sizeof(T), n1, n2, n3, n4) <<
                    std::endl;
                exit(EXIT_FAILURE);
            }
#endif
            return arr;
        }

        // deallocator

        template <typename T>
        void deallocate(T *&arr)
        {
            delete [] arr;
        }

        template <typename T>
        void deallocate(T **&arr)
        {
            delete [] arr[0];
            delete [] arr;
        }

        template <typename T>
        void deallocate(T ***&arr)
        {
            delete [] arr[0][0];
            delete [] arr[0];
            delete [] arr;
        }

        template <typename T>
        void deallocate(T ****&arr)
        {
            delete [] arr[0][0][0];
            delete [] arr[0][0];
            delete [] arr[0];
            delete [] arr;
        }

        // memsize calculator

        unsigned long memsize_in_MB(const int size_of_one,
                                    const unsigned int n1)
        {
            unsigned long n = n1 * size_of_one;
            return n / 1000000;
        }

        unsigned long memsize_in_MB(const int size_of_one,
                                    const unsigned int n1,
                                    const unsigned int n2)
        {
            unsigned long n = n1 * n2 * size_of_one;
            return n / 1000000;
        }

        unsigned long memsize_in_MB(const int size_of_one,
                                    const unsigned int n1,
                                    const unsigned int n2,
                                    const unsigned int n3)
        {
            unsigned long n = n1 * n2 * n3 * size_of_one;
            return n / 1000000;
        }

        unsigned long memsize_in_MB(const int size_of_one,
                                    const unsigned int n1,
                                    const unsigned int n2,
                                    const unsigned int n3,
                                    const unsigned int n4)
        {
            unsigned long n = n1 * n2 * n3 * n4 * size_of_one;
            return n / 1000000;
        }
    };
}

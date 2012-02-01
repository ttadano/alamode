#ifndef ALM_MEMORY_HEADER
#define ALM_MEMORY_HEADER

#include <iostream>
#include <new>
#include "pointers.h"

namespace ALM_NS {
    class Memory: protected Pointers {
    public:
        Memory(class ALM *);

        // allocator (need to be improved)

        template <typename T>
        T *allocate(T *&arr, int n1){
            try{
                arr = new T [n1];
            }
            catch (std::bad_alloc &ba)
            {
                std::cout << "Caught an exception when trying to allocate 1-dimensional array" << std::endl;
                std::cout << ba.what() << " : Array size (MB) = " << memsize_in_MB(sizeof(T), n1) << std::endl;
                exit(EXIT_FAILURE);
            }
            return arr;
        }

        template <typename T>
        T **allocate(T **&arr, int n1, int n2){
            try{
                arr = new T *[n1];
                arr[0] = new T [n1 * n2];
                for (int i = 1; i < n1; i++){
                    arr[i] = arr[0] + i * n2;
                }
            }
            catch (std::bad_alloc &ba)
            {
                std::cout << "Caught an exception when trying to allocate 2-dimensional array" << std::endl;
                std::cout << ba.what() << " : Array size (MB) = " << memsize_in_MB(sizeof(T), n1, n2) << std::endl;
                exit(EXIT_FAILURE);
            }
            return arr;
        }

        template <typename T>
        T ***allocate(T ***&arr, int n1, int n2, int n3){
            try{
                arr = new T **[n1];
                arr[0] = new T *[n1 * n2];
                arr[0][0] = new T [n1 * n2 * n3];
                for (int i = 0; i < n1; i++){
                    arr[i] = arr[0] + i * n2;
                    for (int j = 0; j < n2; j++){
                        arr[i][j] = arr[0][0] + i * n2 * n3 + j * n3;
                    }
                }
            }
            catch(std::bad_alloc &ba)
            {
                std::cout << "Caught an exception when trying to allocate 3-dimensional array" << std::endl;
                std::cout << ba.what() << " : Array size (MB) = " << memsize_in_MB(sizeof(T), n1, n2, n3) << std::endl;
                exit(EXIT_FAILURE);
            }
            return arr;
        }

        // deallocator

        template <typename T>
        void deallocate(T *&arr){
            delete [] arr;
        }

        template <typename T>
        void deallocate(T **&arr){
            delete [] arr[0];
            delete [] arr;
        }

        template <typename T>
        void deallocate(T ***&arr){
            delete [] arr[0][0];
            delete [] arr[0];
            delete [] arr;
        }

        // memsize calculator

        int memsize_in_MB(int size_of_one, int n1){
            unsigned long n = n1 * size_of_one;
            return n / 1000000;
        }
        int memsize_in_MB(int size_of_one, int n1, int n2){
            unsigned long n = n1 * n2 * size_of_one;
            return n / 1000000;
        }
        int memsize_in_MB(int size_of_one, int n1, int n2, int n3){
            unsigned long n = n1 * n2 * n3 * size_of_one;
            return n / 1000000;
        }
    };
}

#endif
#ifndef ALM_MEMORY_HEADER
#define ALM_MEMORY_HEADER

#include "pointers.h"

namespace ALM_NS {
    class Memory: protected Pointers {
    public:
        Memory(class ALM *);

        //		double **create_d2_array(int, int, double **);
        //		void destroy_d2_array();
        //		int **create_i2_array(int, int);
        //	void destroy_i2_array();

        // 2d-array allocation
        template <typename T>
        T **allocate(T **&arr, int n1, int n2){
            arr = new T *[n1];
            arr[0] = new T [n1 * n2];
            for (int i = 1; i < n1; i++){
                arr[i] = arr[0] + i * n2;
            }
            return arr;
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

        template <typename T>
        T ***allocate(T ***&arr, int n1, int n2, int n3){
            arr = new T **[n1];
            arr[0] = new T *[n1 * n2];
            arr[0][0] = new T [n1 * n2 * n3];
            for (int i = 0; i < n1; i++){
                arr[i] = arr[0] + i * n2;
                for (int j = 0; j < n2; j++){
                    arr[i][j] = arr[0][0] + i * n2 * n3 + j * n3;
                }
            }
            return arr;
        }
        //		double **create(int n1, int n2){
        //	bigint nlen = ((bigint) sizeof(T)) * n1 * n2;
        //		double **arr = new double*[n1];
        //	for (int i = 0; i < n1; i++){
        //	arr[i] = new double[n2];
        //}
        //		return arr;
        //	}
    };
}

#endif
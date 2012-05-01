#pragma once

#include <string>
#include "pointers.h"


namespace ALM_NS {
    class Interaction: protected Pointers {
    public:
        Interaction(class ALM *);
        ~Interaction();
        void init();

        double **rcs;
        double **distlist;

        double distance(double *, double *);

        bool is_periodic[3];

        std::string *str_order;

        int nneib;
        int maxorder;
        double ***xcrd;

        int ***intpairs;
        int **ninter;
        double ****relvec;

        double ***minvec;

        bool is_incutoff(int, int *);

        template <typename T>
        T maxval(int n, T *arr)
        {
            T tmp;
            tmp = arr[0];

            for (int i = 0; i < n; i++) {
                tmp = std::max<T>(tmp, arr[i]);
            }
            return tmp;
        }

        template <typename T>
        T maxval(int n1, int n2, T **arr)
        {
            T tmp;
            tmp = arr[0][0];

            for (int i = 0; i < n1; i++) {
                for (int j = 0; j < n2; j++){
                    tmp = std::max<T>(tmp, arr[i][j]);
                } 
            }
            return tmp;
        }

        template <typename T>
        T maxval(int n1, int n2, int n3, T ***arr)
        {
            T tmp;
            tmp = arr[0][0][0];

            for (int i = 0; i < n1; i++) {
                for (int j = 0; j < n2; j++){
                    for (int k = 0; k < n3; k++){
                        tmp = std::max<T>(tmp, arr[i][j][k]);
                    } 
                }
            }
            return tmp;
        }

        template <typename T>
        void insort(int n, T *arr)
        {
            int i, j;
            T tmp;

            for (i = 1; i < n; ++i){
                tmp = arr[i];
                for (j = i - 1; j >= 0 && arr[j] > tmp; --j){
                    arr[j + 1] = arr[j];
                }
                arr[j + 1] = tmp;
            }
        }

    private:
        int nsize[3];
        void calc_distlist(int, double **);
        void search_interactions();
        void set_ordername();
        void calc_minvec();
    };
}

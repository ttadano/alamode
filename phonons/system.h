#pragma once

#include "pointers.h"
#include <string>

namespace PHON_NS {
    class System: protected Pointers {
    public:
        System(class PHON *);
        ~System();

        void setup();

        double lavec_s[3][3], rlavec_s[3][3];
        double lavec_p[3][3], rlavec_p[3][3];
        double **xr_p, **xr_s, **xc;
        double volume_p;

        unsigned int nat, natmin, ntran;
        unsigned int *kd, nkd;

        unsigned int **map_p2s;
        class Maps {
        public:
            unsigned int atom_num;
            unsigned int tran_num;
        };
        Maps *map_s2p;

        std::string *symbol_kd;
        double *mass_kd, *mass;

        double Tmin, Tmax, dT;

        void rotvec(double [3], double [3], double [3][3], char mode = 'N');
        double volume(double [3], double [3], double [3]);

    private:
        void load_system_info();
        void recips(double [3][3], double [3][3]);
    };
}
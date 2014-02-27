#pragma once

#include "pointers.h"
#include "kpoint.h"
#include <vector>
#include <complex>

namespace PHON_NS
{
    class Dos: protected Pointers{
    public:
        Dos(class PHON *);
        ~Dos();

        void setup();
        void calc_dos_all();

        bool flag_dos;
        bool projected_dos, two_phonon_dos;

        int n_energy;
        double emin, emax, delta_e;
        double *energy_dos;
        double *dos_phonon;
        double **pdos_phonon;
        double **dos2_phonon;

    private:
        unsigned int nk_irreducible;
        int *kmap_irreducible;
        std::vector<int> k_irreducible;

        void calc_dos(const unsigned int, int *, double **, const unsigned int, double *, 
            double *, const unsigned int, const int, std::vector<std::vector<KpointList> > &);
        void calc_atom_projected_dos(const unsigned int, double **, const unsigned int, double *,
            double **, const unsigned int, const unsigned int, const int, std::complex<double> ***);

        void calc_two_phonon_dos(const unsigned int, double *, double **, const int,
            std::vector<std::vector<KpointList> >);

    };
}

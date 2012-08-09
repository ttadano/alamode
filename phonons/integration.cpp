#include "integration.h"
#include "kpoint.h"
#include "memory.h"

using namespace PHON_NS;

Integration::Integration(PHON *phon): Pointers(phon) {}

Integration::~Integration(){};


void Integration::setup_integration()
{
    unsigned int nk = kpoint->nk;
    unsigned int nkx = kpoint->nkx;
    unsigned int nky = kpoint->nky;
    unsigned int nkz = kpoint->nkz;

    unsigned int ntetra = 6 * nk;
    memory->allocate(tetras, ntetra, 4);
    prepare_tetrahedron(nkx, nky, nkz); 
}

void Integration::finish_integration()
{
    memory->deallocate(tetras);
}

void Integration::prepare_tetrahedron(const int nk1, const int nk2, const int nk3)
{
    int i, j, k, ii, jj, kk;
    int n1, n2, n3, n4, n5, n6, n7, n8;
    int m;

    int nk23 = nk2 * nk3;

    for (i = 0; i < nk1; ++i){
        for (j = 0; j < nk2; ++j){
            for (k = 0; k < nk3; ++k){

                ii = (i + 1) % nk1;
                jj = (j + 1) % nk2;
                kk = (k + 1) % nk3;

                n1 =  k +  j*nk3 +  i*nk23;
                n2 =  k +  j*nk3 + ii*nk23;
                n3 =  k + jj*nk3 +  i*nk23;
                n4 =  k + jj*nk3 + ii*nk23;
                n5 = kk +  j*nk3 +  i*nk23;
                n6 = kk +  j*nk3 + ii*nk23;
                n7 = kk + jj*nk3 +  i*nk23;
                n8 = kk + jj*nk3 + ii*nk23;

                m = 6 * (k + j*nk3 + i*nk23);

                tetras[m][0] = n1;
                tetras[m][1] = n2;
                tetras[m][2] = n3;
                tetras[m][3] = n6;

                ++m;

                tetras[m][0] = n2;
                tetras[m][1] = n3;
                tetras[m][2] = n4;
                tetras[m][3] = n6;

                ++m;

                tetras[m][0] = n1;
                tetras[m][1] = n3;
                tetras[m][2] = n5;
                tetras[m][3] = n6;

                ++m;

                tetras[m][0] = n3;
                tetras[m][1] = n4;
                tetras[m][2] = n6;
                tetras[m][3] = n8;

                ++m;

                tetras[m][0] = n3;
                tetras[m][1] = n6;
                tetras[m][2] = n7;
                tetras[m][3] = n8;

                ++m;

                tetras[m][0] = n3;
                tetras[m][1] = n5;
                tetras[m][2] = n6;
                tetras[m][3] = n7;
            }
        }
    }
}

double Integration::do_tetrahedron(double *energy, double *f, const double e_ref)
{
    /*
     This function returns the summation of the given function f_{k}
     over the k-points which have the energy "e_ref" using the tetrahedron method.
  
     Ret(e_ref) = \int f(k) \delta(e_ref - energy(k))

    */

    double ret = 0.0;


}

double Integration::fij(const double ei, const double ej, const double e)
{
    return (e - ej) / (ei - ej);
}
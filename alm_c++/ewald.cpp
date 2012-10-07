#include "ewald.h"
#include "memory.h"
#include "error.h"
#include "system.h"
#include "constants.h"
#include <fstream>
#include <cmath>
#include <boost/math/special_functions/erf.hpp>
#include <vector>

using namespace ALM_NS;

Ewald::Ewald(ALM *alm): Pointers(alm) {
    Nmax = 1;
}

Ewald::~Ewald() {
    if (is_longrange) memory->deallocate(Q);
}

void Ewald::init()
{
    unsigned int i, j;

    if (is_longrange) {
        memory->allocate(Q, system->nat);
        load_charge();

        double tmp[3], d_tmp;
        Lmax = 10.0 / Bohr_in_Angstrom;

        for (i = 0; i < 3; ++i) {
            for (j = 0; j < 3; ++j) tmp[j] = 0.0;
            tmp[i] = 1.0;
            system->rotvec(tmp, tmp, system->lavec);
            d_tmp = 0.5 * std::sqrt(tmp[0]*tmp[0] + tmp[1]*tmp[1] + tmp[2]*tmp[2]);
            Lmax = std::min<double>(Lmax, d_tmp);
        }

        std::cout << "Ewald summation will be performed with the following parameters:" << std::endl;
        std::cout << "alpha = " << alpha << " , Gmax = " << Gmax << std::endl;
        std::cout << "Nmax = " << Nmax << " , Lmax = " << Lmax << std::endl;

        prepare_G(Gmax);

        double **f_tmp;
        memory->allocate(f_tmp, system->nat, 3);
        
        std::cout << "Short" << std::endl;
        ewald_force_short(system->nat, system->xcoord, f_tmp);

        for (i = 0; i < system->nat; ++i){
            for (j = 0; j < 3; ++j) {
                std::cout << std::setw(15) << f_tmp[i][0];
            }
            std::cout << std::endl;
        }
        std::cout << "Long" << std::endl;
        ewald_force_long(system->nat, system->xcoord, f_tmp);

        for (i = 0; i < system->nat; ++i){
            
            for (j = 0; j < 3; ++j) {
                std::cout << std::setw(15) << f_tmp[i][0];
            }
            std::cout << std::endl;
        }

        for (i = 0; i <= 100; ++i) {
            for (j = 8; j < 24; ++j) Q[j] = -1.0 + static_cast<double>(i) * 0.01;
            for (j = 24; j < 54; ++j) Q[j] = static_cast<double>(i) * 0.01 * 16.0 / 30.0;
            ewald_force_short(system->nat, system->xcoord, f_tmp);
            std::cout << "#Q_Ga = " << Q[8];
            std::cout << " f1_x, f9_x, f25_x = " << f_tmp[0][0] << " " << f_tmp[8][0] << " " << f_tmp[24][0];
            prepare_G(static_cast<double>(j));
            ewald_force_long(system->nat, system->xcoord, f_tmp);
            std::cout << " : " << f_tmp[0][0] << " " << f_tmp[8][0] << " " << f_tmp[24][0] << std::endl;
        }

/*
        for (i = 1; i < 500; ++i){
            for (j = 1; j < 31; ++j){
                alpha = static_cast<double>(i) * 0.01;
                ewald_force_short(system->nat, system->xcoord, f_tmp);
                std::cout << "#alpha = " << alpha << " Gmax = " << j;
                std::cout << " f1_x,f9_x = " << f_tmp[0][0] << " " << f_tmp[8][0];
                prepare_G(static_cast<double>(j));
                ewald_force_long(system->nat, system->xcoord, f_tmp);
                std::cout << " : " << f_tmp[0][0] << " " << f_tmp[8][0] << std::endl;
            }
        }
*/

        memory->deallocate(f_tmp);

        error->exit("hoge", "hoge");
    }    
}

void Ewald::load_charge()
{
    unsigned int i;
    std::ifstream ifs_long;

    ifs_long.open(file_longrange.c_str(), std::ios::in);
    if (!ifs_long) error->exit("load_charge", "cannot open file_longrange");

    ifs_long >> alpha >> Gmax;

    double Qsum = 0.0;
    for (i = 0; i < system->nat; ++i){
        ifs_long >> Q[i];
        Qsum += Q[i];
    }

    if (std::abs(Qsum) > eps8) error->exit("load_charge", "The system is not neutral.");

    ifs_long.close();
}

void Ewald::prepare_G(const double gmax)
{
    int i, j, k;
    int m;

    int ix, iy, iz;

    unsigned int nx_G = static_cast<unsigned int>(gmax);
    G_vector.clear();

    double x_tmp[3];

//    std::cout << "Generating G vectors for the ewald sum ...";

    for (i = 0; i < 2 * nx_G + 1; ++i){
        ix = i;
        if (ix > nx_G) ix -= 2 * nx_G + 1;
        for (j = 0; j < 2 * nx_G + 1; ++j){
            iy = j;
            if (iy > nx_G) iy -= 2 * nx_G + 1;
            for (k = 0; k < 2 * nx_G + 1; ++k){
                iz = k;
                if (iz > nx_G) iz -= 2 * nx_G + 1;

                if (ix == 0 && iy == 0 && iz == 0) continue;

                x_tmp[0] = static_cast<double>(ix);
                x_tmp[1] = static_cast<double>(iy);
                x_tmp[2] = static_cast<double>(iz);

                if (std::sqrt(std::pow(x_tmp[0], 2) + std::pow(x_tmp[1], 2) + std::pow(x_tmp[2], 2)) <= gmax)
                {
                    system->rotvec(x_tmp, x_tmp, system->rlavec, 'T');
                    G_vector.push_back(Gvecs(x_tmp));
                }            
            }
        }
    }

//    std::cout << "done." << std::endl;
//    std::cout << "The number of G vectors : " << G_vector.size() << std::endl;
}

void Ewald::ewald_force(const int N, double **r, double **f)
{
    // r: atomic coordinate in fractional representation of the supercell
    int i, j;

    for (i = 0; i < N; ++i) {
        for (j = 0; j < 3; ++j) {
            f[i][j] = 0.0;
        }
    }

    double **f_tmp;

    memory->allocate(f_tmp, N, 3);
    ewald_force_short(N, r, f_tmp);
    ewald_force_long(N, r, f_tmp);
}

void Ewald::ewald_force_short(const int N, double **r, double **f)
{
    // Assume the cutoff length in real-space as L/2.

    int i, j;
    int iat, jat;
    int icell, jcell, kcell;
    double rij, tmp;
    double r_tmp[3], f_tmp[3];

    for (i = 0; i < N; ++i){
        for (j = 0; j < 3; ++j){
            f[i][j] = 0.0;
        }
    }

    for (iat = 0; iat < N; ++iat){

        for (i = 0; i < 3; ++i) f_tmp[i]= 0.0;

        for (jat = 0; jat < N; ++jat){

            for (icell = -Nmax; icell <= Nmax; ++icell) { 
                for (jcell = -Nmax; jcell <= Nmax; ++jcell) { 
                    for (kcell = -Nmax; kcell <= Nmax; ++kcell) { 

                        if (jat == iat && icell == 0 && jcell == 0 && kcell == 0) continue;

                        r_tmp[0] = r[iat][0] - r[jat][0] + static_cast<double>(icell);
                        r_tmp[1] = r[iat][1] - r[jat][1] + static_cast<double>(jcell);
                        r_tmp[2] = r[iat][2] - r[jat][2] + static_cast<double>(kcell);

/*                        for (j = 0; j < 3; ++j) {
                            r_tmp[j] = r[iat][j] - r[jat][j];
                            r_tmp[j] = r_tmp[j] - static_cast<double>(nint(r_tmp[j]));
                        }
*/
                        system->rotvec(r_tmp, r_tmp, system->lavec);
                        rij = std::sqrt(std::pow(r_tmp[0], 2) + std::pow(r_tmp[1], 2) + std::pow(r_tmp[2], 2));
                        if (rij >= Lmax) continue;

//                        std::cout << "r[" << iat << "," << jat << "] = " << r_tmp[0] << " " << r_tmp[1] << " " << r_tmp[2] << std::endl;
//                        std::cout << "r[" << iat << "," << jat << "] = " << rij << std::endl;

                        tmp = Q[jat] * (2.0 * alpha * std::exp(-alpha * rij) / std::sqrt(pi) + boost::math::erfc(alpha * rij) / rij) / std::pow(rij, 2);
                        for (i = 0; i < 3; ++i) f_tmp[i] += tmp * r_tmp[i];
                    }
                }
            }

        }

        for (i = 0; i < 3; ++i) f[iat][i] = 4.0 * Q[iat] * f_tmp[i];
    }
}

void Ewald::ewald_force_long(const int N, double **r, double **f)
{
    unsigned int iat, jat;
    unsigned int i, j;
    double f_tmp[3], r_tmp[3];
    double g_tmp[3], gnorm;
    double tmp;

    for (i = 0; i < N; ++i){
        for (j = 0; j < 3; ++j){
            f[i][j] = 0.0;
        }
    }

    for (iat = 0; iat < N; ++iat){

        for (i = 0; i < 3; ++i ) f_tmp[i]= 0.0;

        for (std::vector<Gvecs>::iterator it = G_vector.begin(); it != G_vector.end(); ++it) {

            for (i = 0; i < 3; ++i) g_tmp[i] = (*it).vec[i];

            gnorm = std::pow(g_tmp[0], 2) + std::pow(g_tmp[1], 2) + std::pow(g_tmp[2], 2);

            tmp = 0.0;

            for (jat = 0; jat < N; ++jat){

                for (j = 0; j < 3; ++j) {
                    r_tmp[j] = r[iat][j] - r[jat][j];
//                    r_tmp[j] = r_tmp[j] - static_cast<double>(nint(r_tmp[j]));
                }
                system->rotvec(r_tmp, r_tmp, system->lavec);
                tmp += Q[jat] * std::sin(g_tmp[0] * r_tmp[0] + g_tmp[1] * r_tmp[1] + g_tmp[2] * r_tmp[2]);
            }

            tmp *= std::exp(-gnorm / (4.0 * std::pow(alpha, 2))) / gnorm;

            for (i = 0; i < 3; ++i) f_tmp[i] += tmp * g_tmp[i];
        }
        for (i = 0; i < 3; ++i) f[iat][i] = 8.0 * pi * Q[iat] * f_tmp[i] / system->cell_volume;
    }
}

double Ewald::dij(double r1[3], double r2[3])
{
    return std::abs(std::pow(r1[0] - r2[0], 2) + std::pow(r1[1] - r2[1], 2) + std::pow(r1[2] - r2[2], 2));
}

int Ewald::nint(const double x)
{
    return static_cast<int>(x + 0.5 - (x < 0.0));
}

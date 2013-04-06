#include "mpi_common.h"
#include <iostream>
#include <iomanip>
#include "dynamical.h"
#include "error.h"
#include "fcs_phonon.h"
#include "gruneisen.h"
#include "pointers.h"
#include "kpoint.h"
#include "memory.h"
#include "system.h"

using namespace PHON_NS;

Gruneisen::Gruneisen(PHON *phon): Pointers(phon){};
Gruneisen::~Gruneisen(){};

void Gruneisen::setup()
{

	MPI_Bcast(&delta_a, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if (mympi->my_rank == 0) {
    std::cout << std::endl; 
    std::cout << "Change in cell parameters : " << delta_a << std::endl;
	}

    memory->allocate(gruneisen, kpoint->nk, dynamical->neval);

    prepare_delta_fc2();
    prepare_newfc2();
}

void Gruneisen::calc_gruneisen()
{
    std::cout << "Calculating Gruneisen parameters ..." << std::endl;

    unsigned int nk = kpoint->nk;
    unsigned int ns = dynamical->neval;

    unsigned int i;
    unsigned int ik, is;

    double *eval_orig;
    double *eval_plus;
    double *eval_minus;

    double xk_tmp[3];

    std::complex<double> **evec_tmp;
    memory->allocate(evec_tmp, 1, 1); // dummy allocation
    
    memory->allocate(eval_orig, ns);
    memory->allocate(eval_plus, ns);
    memory->allocate(eval_minus, ns);


    for (ik = 0; ik < nk; ++ik){
        for (i = 0; i < 3; ++i) {
            xk_tmp[i] = kpoint->xk[ik][i];
        }

        dynamical->eval_k(xk_tmp, fcs_phonon->fc2, eval_orig, evec_tmp, false);
        dynamical->eval_k(xk_tmp, fc2_plus, eval_plus, evec_tmp, false);
        dynamical->eval_k(xk_tmp, fc2_minus, eval_minus, evec_tmp, false);

         for (is = 0; is < ns; ++is) {
             gruneisen[ik][is] = (eval_plus[is] - eval_minus[is]) / (2.0 * delta_a) / (-6.0 * eval_orig[is]);
         }
    }

    memory->deallocate(evec_tmp);
    memory->deallocate(eval_orig);
    memory->deallocate(eval_plus);
    memory->deallocate(eval_minus);

    std::cout << "done !" << std::endl;
}

void Gruneisen::finish_gruneisen()
{
    memory->deallocate(gruneisen);
    memory->deallocate(dfc2);
    memory->deallocate(fc2_plus);
    memory->deallocate(fc2_minus);
}

void Gruneisen::prepare_delta_fc2()
{
    unsigned int i;
    unsigned int iat, jat, kat;
    unsigned int icrd, jcrd, kcrd;
    double coord_tmp[3];

    unsigned int nalpha[3], ncell[3], ixyz[3];

    unsigned int natmin = system->natmin;
    unsigned int nat = system->nat;

    std::cout << "Preparing delta FC2 from cubic force constants ...";

    memory->allocate(dfc2, natmin, nat, 3, 3);

    for (iat = 0; iat < natmin; ++iat){
        for (jat = 0; jat < nat; ++jat){
            for (icrd = 0; icrd < 3; ++icrd){
                for (jcrd = 0; jcrd < 3; ++jcrd){
                    dfc2[iat][jat][icrd][jcrd] = 0.0;
                }
            }
        }
    }

    for (std::vector<FcsClass>::iterator it = fcs_phonon->force_constant[1].begin(); it != fcs_phonon->force_constant[1].end(); ++it) {

        FcsClass fc3_tmp = *it;

        iat = fc3_tmp.elems[0].cell * natmin + fc3_tmp.elems[0].atom;
        jat = fc3_tmp.elems[1].cell * natmin + fc3_tmp.elems[1].atom;
        kat = fc3_tmp.elems[2].cell * natmin + fc3_tmp.elems[2].atom;

        icrd = fc3_tmp.elems[0].xyz;
        jcrd = fc3_tmp.elems[1].xyz;
        kcrd = fc3_tmp.elems[2].xyz;

        for (i = 0; i < 3; ++i){
            nalpha[i] = fc3_tmp.elems[i].atom;
            ncell[i]  = fc3_tmp.elems[i].cell;
            ixyz[i] = fc3_tmp.elems[i].xyz;
        }

        iat = system->map_p2s[nalpha[0]][ncell[0]];
        jat = system->map_p2s[nalpha[1]][ncell[1]];
        kat = system->map_p2s[nalpha[2]][ncell[2]];

        for (i = 0; i < 3; ++i) {
            coord_tmp[i] = system->xr_s[kat][i] - system->xr_s[iat][i];
            coord_tmp[i] = dynamical->fold(coord_tmp[i]);
        }
        system->rotvec(coord_tmp, coord_tmp, system->lavec_s);

        dfc2[nalpha[0]][jat][icrd][jcrd] += fc3_tmp.fcs_val * coord_tmp[ixyz[2]];
    }

#ifdef _DEBUG
    for (iat = 0; iat < natmin; ++iat){
        for (jat = 0; jat < nat; ++jat){
            for (icrd = 0; icrd < 3; ++icrd){
                for (jcrd = 0; jcrd < 3; ++jcrd){
                    std::cout << std::setw(15) << dfc2[iat][jat][icrd][jcrd] << std::endl;
                }
            }
        }
    }
#endif

    std::cout << "done !" << std::endl;
}

void Gruneisen::prepare_newfc2()
{
    unsigned int iat, jat;
    unsigned int icrd, jcrd;

    unsigned int natmin = system->natmin;
    unsigned int nat = system->nat;

    memory->allocate(fc2_plus, system->natmin, system->nat, 3, 3);
    memory->allocate(fc2_minus, system->natmin, system->nat, 3, 3);


    for (iat = 0; iat < natmin; ++iat){
        for (jat = 0; jat < nat; ++jat){
            for (icrd = 0; icrd < 3; ++icrd){
                for (jcrd = 0; jcrd < 3; ++jcrd){
                    fc2_plus[iat][jat][icrd][jcrd]  = fcs_phonon->fc2[iat][jat][icrd][jcrd] + delta_a * dfc2[iat][jat][icrd][jcrd];
                    fc2_minus[iat][jat][icrd][jcrd] = fcs_phonon->fc2[iat][jat][icrd][jcrd] - delta_a * dfc2[iat][jat][icrd][jcrd];
                }
            }
        }
    }

#ifdef _DEBUG
    double fc2_tmp;

    for (iat = 0; iat < natmin; ++iat){

        for (icrd = 0; icrd < 3; ++icrd){
            for (jcrd = 0; jcrd < 3; ++jcrd){

                fc2_tmp = 0.0;
                for (jat = 0; jat < natmin; ++jat){
                    for (unsigned int itran = 0; itran < system->ntran; ++itran){

                        fc2_tmp += fc2_plus[iat][system->map_p2s[jat][itran]][icrd][jcrd];
                    }
                }
                std::cout << "fc2_tmp = " << fc2_tmp << std::endl;
            }
        }
    }
#endif
}

void Gruneisen::calc_gruneisen2()
{
    unsigned int is, ik;
    unsigned int i, j;
    unsigned int ns = dynamical->neval;
    unsigned int nk = kpoint->nk;

    double gamma_imag;

    dynamical->diagonalize_dynamical_all();

    memory->allocate(gruneisen, nk, ns);

    std::complex<double> **dfc2_reciprocal;

    memory->allocate(dfc2_reciprocal, ns, ns);

    std::cout << "Calculating Gruneisen parameters ..." << std::endl;

    for (ik = 0; ik < nk; ++ik){
        for (is = 0; is < ns; ++is){

            calc_dfc2_reciprocal(dfc2_reciprocal, kpoint->xk[ik]);

            gruneisen[ik][is] = std::complex<double>(0.0, 0.0);

            for (i = 0; i < ns; ++i){
                for (j = 0; j < ns; ++j){
                    gruneisen[ik][is] += std::conj(dynamical->evec_phonon[ik][is][i]) * dfc2_reciprocal[i][j] * dynamical->evec_phonon[ik][is][j];
                }
            }

            gamma_imag = gruneisen[ik][is].imag();
            if (std::abs(gamma_imag) > eps10) {
                error->warn("calc_gruneisen", "Gruneisen parameter is not real");
            }

            gruneisen[ik][is] /= - 6.0 * std::pow(dynamical->eval_phonon[ik][is], 2);
        }
    }
    memory->deallocate(dfc2_reciprocal);
}

void Gruneisen::calc_dfc2_reciprocal(std::complex<double> **dphi2, double *xk_in)
{
    unsigned int i, j;
    unsigned int icrd, jcrd, itran;
    unsigned int ns = dynamical->neval;
    unsigned int natmin = system->natmin;
    unsigned int ntran = system->ntran;

    unsigned int atm_p1, atm_p2, atm_s2;

    double phase;
    double vec[3];

    std::complex<double> exp_phase;
    std::complex<double> im(0.0, 1.0);
    std::complex<double> ctmp[3][3];


    for (i = 0; i < ns; ++i){
        for (j = 0; j < ns; ++j){
            dphi2[i][j] = std::complex<double>(0.0, 0.0);
        }
    }

    for (i = 0; i < natmin; ++i){

        atm_p1 = system->map_p2s[i][0];

        for (j = 0; j < natmin; ++j){

            atm_p2 = system->map_p2s[j][0];

            for (icrd = 0; icrd < 3; ++icrd){
                for (jcrd = 0; jcrd < 3; ++jcrd){
                    ctmp[icrd][jcrd] = std::complex<double>(0.0, 0.0);
                }
            }

            for (itran = 0; itran < ntran; ++itran){

                atm_s2 = system->map_p2s[j][itran];

                for (icrd = 0; icrd < 3; ++icrd){
                    if (system->cell_dimension[icrd] == 1) {
                        vec[icrd] = system->xr_s[atm_p1][icrd] - system->xr_s[atm_s2][icrd];
                        if (std::abs(vec[icrd]) < 0.5) {
                            vec[icrd] = 0.0;
                        } else {
                            if (system->xr_s[atm_p1][icrd] < 0.5) {
                                vec[icrd] = -1.0;
                            } else {
                                vec[icrd] = 1.0;
                            }
                        }
                    } else if (system->cell_dimension[icrd] == 2){
                        vec[icrd] = system->xr_s[atm_p2][icrd] - system->xr_s[atm_s2][icrd];
                        vec[icrd] = dynamical->fold(vec[icrd]);
                        if (std::abs(system->xr_s[atm_p1][icrd] - system->xr_s[atm_s2][icrd]) > 0.5) vec[icrd] *= -1.0;
                    } else {
                        //     vec[icrd] = system->xr_s[atm_p2][icrd] - system->xr_s[atm_s2][icrd];
                        vec[icrd] = system->xr_s[atm_p1][icrd] - system->xr_s[atm_s2][icrd];
                        vec[icrd] = dynamical->fold(vec[icrd]);
                    }
                }

                system->rotvec(vec, vec, system->lavec_s);
                system->rotvec(vec, vec, system->rlavec_p);

                phase = vec[0] * xk_in[0] + vec[1] * xk_in[1] + vec[2] * xk_in[2];
                exp_phase = std::exp(im * phase);

                if (std::abs(exp_phase.real()) < 1.0e-14) {
                    exp_phase = std::complex<double>(0.0, exp_phase.imag());
                }
                if (std::abs(exp_phase.imag()) < 1.0e-14) {
                    exp_phase = std::complex<double>(exp_phase.real(), 0.0);
                }


                for (icrd = 0; icrd < 3; ++icrd){
                    for (jcrd = 0; jcrd < 3; ++jcrd){
                        ctmp[icrd][jcrd] += dfc2[i][atm_s2][icrd][jcrd] * std::exp(im * phase);
                    }
                }
            }

#ifdef _DEBUG
            if (i != j) {
                std::cout << "i = " << i << " , j = " << j << std::endl;
                std::cout << "xk = " << xk_in[0] << " " << xk_in[1] << " " << xk_in[2] << std::endl;
                std::cout << "ctmp = " << ctmp[0][0] << std::endl;
            }
#endif

            for (icrd = 0; icrd < 3; ++icrd){
                for (jcrd = 0; jcrd < 3; ++jcrd){
                    dphi2[3 * i + icrd][3 * j + jcrd] = ctmp[icrd][jcrd] / std::sqrt(system->mass[atm_p1] * system->mass[atm_p2]);
                }
            }

        }
    }

}

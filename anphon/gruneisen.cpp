/*
gruneisen.cpp

Copyright (c) 2014 Terumasa Tadano

This file is distributed under the terms of the MIT license.
Please see the file 'LICENCE.txt' in the root directory 
or http://opensource.org/licenses/mit-license.php for information.
*/

#include "mpi_common.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include "dynamical.h"
#include "error.h"
#include "fcs_phonon.h"
#include "gruneisen.h"
#include "pointers.h"
#include "kpoint.h"
#include "memory.h"
#include "system.h"
#include "parsephon.h"
#include "write_phonons.h"
#include "mathfunctions.h"
#include "relaxation.h"

using namespace PHON_NS;

Gruneisen::Gruneisen(PHON *phon): Pointers(phon){};
Gruneisen::~Gruneisen(){};

void Gruneisen::setup()
{
    MPI_Bcast(&delta_a, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&print_newfcs, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD);

    if (print_gruneisen || print_newfcs) {
        prepare_delta_fc2();
        prepare_newfc2();
    }

    if (print_newfcs && relaxation->quartic_mode) prepare_newfc3();

    memory->allocate(gruneisen, kpoint->nk, dynamical->neval);

    //     if (mympi->my_rank == 0) {
    //         if (print_newfcs) {
    //             if (relaxation->quartic_mode) {
    //                 std::cout << " NEWFCS = 1 : Harmonic and cubic force constants of " << std::endl;
    //             }
    //             else {
    //                 std::cout << " NEWFCS = 1 : Harmonic force constants of " << std::endl;
    //             }
    //             std::cout << "              expanded/compressed systems will be estimated" << std::endl;
    //             std::cout << "              with DELTA_A = " << std::setw(5) << delta_a << std::endl;
    //         }
    //     }
}


void Gruneisen::calc_gruneisen()
{

    if (mympi->my_rank == 0) {
        std::cout << " GRUNEISEN = 1 : Calculating Gruneisen parameters ... ";
    }

    unsigned int nk = kpoint->nk;
    unsigned int ns = dynamical->neval;

    unsigned int i;
    unsigned int ik, is;

    double *eval_orig;
    double **eval_plus;
    double **eval_minus;

    double xk_tmp[3];

    std::complex<double> **evec_tmp;

    memory->allocate(evec_tmp, 1, 1); // dummy allocation
    memory->allocate(eval_orig, ns);
    memory->allocate(eval_plus, nk, ns);
    memory->allocate(eval_minus, nk, ns);

    for (ik = 0; ik < nk; ++ik){

        for (i = 0; i < 3; ++i) xk_tmp[i] = kpoint->xk[ik][i];

        dynamical->eval_k(xk_tmp, kpoint->kvec_na[ik], fcs_phonon->fc2_ext, eval_orig, evec_tmp, false);
        dynamical->eval_k(xk_tmp, kpoint->kvec_na[ik], fc2_plus_ext, eval_plus[ik], evec_tmp, false);
        dynamical->eval_k(xk_tmp, kpoint->kvec_na[ik], fc2_minus_ext, eval_minus[ik], evec_tmp, false);

        for (is = 0; is < ns; ++is) {
            gruneisen[ik][is] = (eval_plus[ik][is] - eval_minus[ik][is]) / (2.0 * delta_a) / (-6.0 * eval_orig[is]);
        }
    }

    //     for (ik = 0; ik < nk; ++ik) {
    //         for (is = 0; is < ns; ++is) {
    //             eval_plus[ik][is] = dynamical->freq(eval_plus[ik][is]);
    //             eval_minus[ik][is] = dynamical->freq(eval_minus[ik][is]);
    //         }
    //     }
    // 
    //     if (kpoint->kpoint_mode == 1) {
    //         std::string file_band_plus, file_band_minus;
    // 
    //         file_band_plus = input->job_title + ".band_+";
    //         file_band_minus = input->job_title + ".band_-";
    // 
    //         std::ofstream ofs_plus, ofs_minus;
    // 
    //         ofs_plus.open(file_band_plus.c_str(), std::ios::out);
    //         if (!ofs_plus) error->exit("calc_gruneisen", "Could not create band_plus file.");
    //         ofs_minus.open(file_band_minus.c_str(), std::ios::out);
    //         if (!ofs_minus) error->exit("calc_gruneisen", "Could not create band_minus file.");
    // 
    //         ofs_plus << "# Phonon energy (cm^-1) of the system expanded by " << std::setw(10) << delta_a * 100 << " %." << std::endl;
    //         ofs_minus << "# Phonon energy (cm^-1) of the system compressed by " << std::setw(10) << delta_a * 100 << " %." << std::endl;
    // 
    //         for (ik = 0; ik < nk; ++ik){
    //             ofs_plus << std::setw(8) << std::fixed << kpoint->kaxis[ik];
    //             ofs_minus << std::setw(8) << std::fixed << kpoint->kaxis[ik];
    // 
    //             for (is = 0; is < ns; ++is){
    //                 ofs_plus << std::setw(15) << std::scientific << writes->in_kayser(eval_plus[ik][is]);
    //                 ofs_minus << std::setw(15) << std::scientific << writes->in_kayser(eval_minus[ik][is]);
    //             }
    //             ofs_plus << std::endl;
    //             ofs_minus << std::endl;
    //         }
    // 
    //         ofs_plus.close();
    //         ofs_minus.close();
    //     }


    memory->deallocate(evec_tmp);
    memory->deallocate(eval_orig);
    memory->deallocate(eval_plus);
    memory->deallocate(eval_minus);

    std::cout << "done !" << std::endl;
}

void Gruneisen::finish_gruneisen()
{
    if (print_gruneisen) memory->deallocate(gruneisen);

    if (print_gruneisen || print_newfcs) {
        memory->deallocate(dfc2);

        fc2_plus_ext.clear();
        fc2_minus_ext.clear();

    }
}

void Gruneisen::prepare_delta_fc2()
{
    unsigned int i;
    unsigned int iat, jat, kat;
    unsigned int icrd, jcrd;
    double coord_tmp[3];

    unsigned int nalpha[3], ncell[3], ixyz[3];

    unsigned int natmin = system->natmin;
    unsigned int nat = system->nat;

    //   std::cout << "Preparing delta FC2 from cubic force constants ...";

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
        rotvec(coord_tmp, coord_tmp, system->lavec_s);

        dfc2[nalpha[0]][jat][ixyz[0]][ixyz[1]] += fc3_tmp.fcs_val * coord_tmp[ixyz[2]];
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

    //   std::cout << "done !" << std::endl;
}

void Gruneisen::prepare_newfc2()
{
    unsigned int iat, jat;
    unsigned int icrd, jcrd;

    unsigned int natmin = system->natmin;
    unsigned int nat = system->nat;

    int icell, jcell, kcell;
    int counter, ncell;

    double xdiff[3];
    int cell[3];
    double dfc2_tmp;
    FcsClassExtent dfc2_ext_tmp;

    //  if (fcs_phonon->is_fc2_ext) {
    fc2_plus_ext.clear();
    fc2_minus_ext.clear();

    for (std::vector<FcsClassExtent>::const_iterator it = fcs_phonon->fc2_ext.begin(); it != fcs_phonon->fc2_ext.end(); ++it) {
        fc2_plus_ext.push_back(*it);
        fc2_minus_ext.push_back(*it);
    }

    for (iat = 0; iat < natmin; ++iat) {
        for (jat = 0; jat < nat; ++jat) {

            for (icrd = 0; icrd < 3; ++icrd) {
                xdiff[icrd] = system->xr_s[jat][icrd] - system->xr_s[system->map_p2s[iat][0]][icrd];

                if (std::abs(xdiff[icrd]-0.5) < eps || std::abs(xdiff[icrd] + 0.5) < eps) {
                    //	error->exit("prepare_newfc2", "multiple interaction exist. This should not occur.");
                } else if (xdiff[icrd] > 0.5) {
                    cell[icrd] = -1;
                } else if (xdiff[icrd] < -0.5) {
                    cell[icrd] = 1;
                } else {
                    cell[icrd] = 0;
                }
            }

            counter = 0;
            ncell = 0;
            for (icell = -1; icell <= 1; ++icell) {
                for (jcell = -1; jcell <= 1; ++jcell) {
                    for (kcell = -1; kcell <= 1; ++kcell) {

                        if (icell == 0 && jcell == 0 && kcell == 0) continue;

                        ++counter;

                        if (icell == cell[0] && jcell == cell[1] && kcell == cell[2]) {
                            ncell = counter;
                            break;
                        }
                    }
                }
            }

            for (icrd = 0; icrd < 3; ++icrd) {
                for (jcrd = 0; jcrd < 3; ++jcrd) {
                    dfc2_tmp = dfc2[iat][jat][icrd][jcrd];

                    if (std::abs(dfc2_tmp) > eps) {
                        //	std::cout << iat << " " << jat << " " << ncell << " " <<  icrd << " " << jcrd << std::endl;
                        dfc2_ext_tmp.atm1 = iat;
                        dfc2_ext_tmp.atm2 = jat;
                        dfc2_ext_tmp.cell_s = ncell;
                        dfc2_ext_tmp.xyz1 = icrd;
                        dfc2_ext_tmp.xyz2 = jcrd;

                        dfc2_ext_tmp.fcs_val = delta_a * dfc2_tmp;
                        fc2_plus_ext.push_back(dfc2_ext_tmp);

                        dfc2_ext_tmp.fcs_val = -delta_a * dfc2_tmp;
                        fc2_minus_ext.push_back(dfc2_ext_tmp);
                    }
                }
            }
        }
    }

    //     } else {
    //         memory->allocate(fc2_plus, system->natmin, system->nat, 3, 3);
    //         memory->allocate(fc2_minus, system->natmin, system->nat, 3, 3);
    // 
    //         for (iat = 0; iat < natmin; ++iat){
    //             for (jat = 0; jat < nat; ++jat){
    //                 for (icrd = 0; icrd < 3; ++icrd){
    //                     for (jcrd = 0; jcrd < 3; ++jcrd){
    //                         fc2_plus[iat][jat][icrd][jcrd]  = fcs_phonon->fc2[iat][jat][icrd][jcrd] + delta_a * dfc2[iat][jat][icrd][jcrd];
    //                         fc2_minus[iat][jat][icrd][jcrd] = fcs_phonon->fc2[iat][jat][icrd][jcrd] - delta_a * dfc2[iat][jat][icrd][jcrd];
    //                     }
    //                 }
    //             }
    //         }
    //     }



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

void Gruneisen::prepare_newfc3()
{
    int i;
    int ind[4], arr_old[4], arr[4];
    int natmin = system->natmin;
    int nat = system->nat;
    unsigned int nalpha[4], ncell[4], ixyz[4];
    unsigned int atom_s[4];

    double coord_tmp[3];
    double fcs_tmp;

    std::vector<unsigned int> arr3, arr4;
    std::vector<FcsClassGru> fc3_copy, fc4_copy, dfc3;
    std::vector<FcsClassGru>::iterator it_lower;

    //    std::cout << " NEWFCS = 1Preparing new FC3 from quartic force constants ...";

    fc4_copy.clear();
    dfc3.clear();

    for (std::vector<FcsClass>::const_iterator it_fc4 = fcs_phonon->force_constant[2].begin(); it_fc4 != fcs_phonon->force_constant[2].end(); ++it_fc4) {
        FcsClass fc4_tmp = *it_fc4;

        for (i = 0; i < 4; ++i) {
            nalpha[i] = fc4_tmp.elems[i].atom;
            ncell[i] = fc4_tmp.elems[i].cell;
            ixyz[i] = fc4_tmp.elems[i].xyz;
            atom_s[i] = system->map_p2s[nalpha[i]][ncell[i]];
            ind[i] = 3 * atom_s[i] + ixyz[i];
        }

        fc4_copy.push_back(FcsClassGru(4, ind, fc4_tmp.fcs_val));
    }

    std::sort(fc4_copy.begin(), fc4_copy.end());


    for (i = 0; i < 3; ++i) arr_old[i] = -1;
    fcs_tmp = 0.0;

    for (std::vector<FcsClassGru>::const_iterator it_fc4 = fc4_copy.begin(); it_fc4 != fc4_copy.end(); ++it_fc4) {
        FcsClassGru fc4_tmp = *it_fc4;

        for (i = 0; i < 4; ++i) arr[i] = fc4_tmp.elems[i];

        if (arr[0] != arr_old[0] || arr[1] != arr_old[1] || arr[2] != arr_old[2]) {
            if (std::abs(fcs_tmp) > eps12) {
                dfc3.push_back(FcsClassGru(3, arr_old, fcs_tmp));
            }
            fcs_tmp = 0.0;

            for (i = 0; i < 3; ++i) arr_old[i] = arr[i];
        }

        for (i = 0; i < 3; ++i) {
            coord_tmp[i] = system->xr_s[arr[3]/3][i] - system->xr_s[arr[0]/3][i];
            coord_tmp[i] = dynamical->fold(coord_tmp[i]);
        }
        rotvec(coord_tmp, coord_tmp, system->lavec_s);

        fcs_tmp += fc4_tmp.fcs_val * coord_tmp[arr[3] % 3];
    }
    if (std::abs(fcs_tmp) > eps12) 	dfc3.push_back(FcsClassGru(3, arr_old, fcs_tmp));


    if (fcs_phonon->force_constant[1].size() < dfc3.size()) {
        error->exit("prepare_newfc3", "The number of DFC3 is greater than the number of FC3.\n \
                                      The cutoff radius for FC4 should not be greater than that of FC3.");
    }
    fc4_copy.clear();
    fc3_copy.clear();


    for (std::vector<FcsClass>::const_iterator it_fc3 = fcs_phonon->force_constant[1].begin(); it_fc3 != fcs_phonon->force_constant[1].end(); ++it_fc3) {
        FcsClass fc3_tmp = *it_fc3;

        for (i = 0; i < 3; ++i) {
            nalpha[i] = fc3_tmp.elems[i].atom;
            ncell[i] = fc3_tmp.elems[i].cell;
            ixyz[i] = fc3_tmp.elems[i].xyz;
            atom_s[i] = system->map_p2s[nalpha[i]][ncell[i]];
            ind[i] = 3 * atom_s[i] + ixyz[i];
        }
        fc3_copy.push_back(FcsClassGru(3, ind, fc3_tmp.fcs_val));
    }
    std::sort(fc3_copy.begin(), fc3_copy.end());

    fc3_plus.clear();
    fc3_minus.clear();

    for (std::vector<FcsClassGru>::const_iterator it = fc3_copy.begin(); it != fc3_copy.end(); ++it) {
        fc3_plus.push_back(*it);
        fc3_minus.push_back(*it);
    }


    for (std::vector<FcsClassGru>::const_iterator it = dfc3.begin(); it != dfc3.end(); ++it) {

        for (i = 0; i < 3; ++i) ind[i] = (*it).elems[i];

        it_lower = lower_bound(fc3_plus.begin(), fc3_plus.end(), FcsClassGru(3, ind, fcs_tmp));
        if (it_lower == fc3_plus.end()) {
            error->exit("prepare_newfc3", "The list of FC3 doesn't contain a force constant in DFC3");
        } else {
            (*it_lower).fcs_val += delta_a * (*it).fcs_val;
        }

        it_lower = lower_bound(fc3_minus.begin(), fc3_minus.end(), FcsClassGru(3, ind, fcs_tmp));
        if (it_lower == fc3_minus.end()) {
            error->exit("prepare_newfc3", "The list of FC3 doesn't contain a force constant in DFC3");
        } else {
            (*it_lower).fcs_val -= delta_a * (*it).fcs_val;
        }
    }

    //   std::cout << "done !" << std::endl;
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
                        vec[icrd] = system->xr_s[atm_p1][icrd] - system->xr_s[atm_s2][icrd];
                        vec[icrd] = dynamical->fold(vec[icrd]);
                        vec[icrd] += system->xr_s[atm_p2][icrd] - system->xr_s[atm_p1][icrd];
                    }
                }

                rotvec(vec, vec, system->lavec_s);
                rotvec(vec, vec, system->rlavec_p);

                phase = vec[0] * xk_in[0] + vec[1] * xk_in[1] + vec[2] * xk_in[2];
                exp_phase = std::exp(im * phase);

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


void Gruneisen::calc_gruneisen3()
{
    unsigned int i, j;
    unsigned int nk = kpoint->nk;
    unsigned int ns = dynamical->neval;

    unsigned int nalpha[3], ncell[3], ixyz[3];
    unsigned int iat, jat, kat;
    unsigned int iat2, jat2, kat2;

    unsigned int ik, is;

    double coord_tmp[3], x[3];
    double phase;
    std::complex<double> im(0.0, 1.0);

    memory->allocate(gruneisen, nk, ns);


    dynamical->diagonalize_dynamical_all();


    for (i = 0; i < nk; ++i) {
        for (j = 0; j < ns; ++j) {
            gruneisen[i][j] = std::complex<double>(0.0, 0.0);
        }
    }

    for (std::vector<FcsClass>::iterator it = fcs_phonon->force_constant[1].begin(); it != fcs_phonon->force_constant[1].end(); ++it) {

        FcsClass fc3_tmp = *it;

        for (i = 0; i < 3; ++i){
            nalpha[i] = fc3_tmp.elems[i].atom;
            ncell[i]  = fc3_tmp.elems[i].cell;
            ixyz[i] = fc3_tmp.elems[i].xyz;
            std::cout << std::setw(5) << nalpha[i];
            std::cout << std::setw(5) << ncell[i];
            std::cout << std::setw(5) << ixyz[i];
            std::cout << "     ";
        }
        std::cout << std::setw(15) << fc3_tmp.fcs_val;
        std::cout << std::endl;

        iat = system->map_p2s[nalpha[0]][ncell[0]];
        jat = system->map_p2s[nalpha[1]][ncell[1]];
        kat = system->map_p2s[nalpha[2]][ncell[2]];

        iat2 = system->map_p2s[nalpha[0]][0];

        jat2 = system->map_p2s[0][ncell[1]];
        kat2 = system->map_p2s[0][ncell[2]];

        for (i = 0; i < 3; ++i) {
            coord_tmp[i] = system->xr_s[kat2][i] - system->xr_s[jat2][i];
            //	coord_tmp[i] = dynamical->fold(coord_tmp[i]);
            x[i] = system->xr_s[iat2][i];
        }
        rotvec(coord_tmp, coord_tmp, system->lavec_s);
        rotvec(coord_tmp, coord_tmp, system->rlavec_p);

        // 		for (i = 0; i < 3; ++i) std::cout << std::setw(15) << coord_tmp[i] / (2.0 * pi);
        // 		std::cout << std::endl;

        rotvec(x, x, system->lavec_s);
        std::cout << "x = " << std::setw(15) << x[0] << std::setw(15) << x[1] << std::setw(15) << x[2] << std::endl;

        for (ik = 0; ik < nk; ++ik) {
            phase = coord_tmp[0] * kpoint->xk[ik][0] + coord_tmp[1] * kpoint->xk[ik][1] + coord_tmp[2] * kpoint->xk[ik][2];
            //	std::cout << "phase = " << phase << std::endl;

            for (is = 0; is < ns; ++is) {
                gruneisen[ik][is] += fc3_tmp.fcs_val * std::exp(im * phase) 
                    * x[ixyz[0]] * std::conj(dynamical->evec_phonon[ik][is][3 * nalpha[1] + ixyz[1]])
                    * dynamical->evec_phonon[ik][is][3 * nalpha[2] + ixyz[2]] / std::sqrt(system->mass[jat] * system->mass[kat]);
                std::cout << "ik = " << std::setw(5) << ik;
                std::cout << "is = " << std::setw(5) << is;
                std::cout << "gru = " << std::setw(15) << gruneisen[ik][is] << std::endl;
            }
        }
    }

    for (ik = 0; ik < nk; ++ik) {
        for (is = 0; is < ns; ++is) {
            gruneisen[ik][is] /= - 6.0 * std::pow(dynamical->eval_phonon[ik][is], 2);
        }
    }
}

void Gruneisen::write_newinfo_all()
{
    std::string file_info_plus, file_info_minus;
    std::ofstream ofs_plus, ofs_minus;
    std::ifstream ifs_orig;

    file_info_plus = input->job_title + "_+.info";
    file_info_minus = input->job_title + "_-.info";

    if (print_newfcs) {
        ofs_plus.open(file_info_plus.c_str(), std::ios::out);
        if (!ofs_plus) error->exit("write_newinfo", "Cannot create the file file_info_plus");
        ofs_minus.open(file_info_minus.c_str(), std::ios::out);
        if (!ofs_minus) error->exit("write_newinfo", "Cannot create the file file_info_minus");
        ifs_orig.open(fcs_phonon->file_fcs.c_str(), std::ios::in);
        if (!ifs_orig) error->exit("write_newinfo", "Cannot open the info file");

        write_newinfo(ifs_orig, ofs_plus, delta_a, fc2_plus, fc2_plus_ext, fc3_plus);
        ifs_orig.close();
        ifs_orig.open(fcs_phonon->file_fcs.c_str(), std::ios::in);
        write_newinfo(ifs_orig, ofs_minus, -delta_a, fc2_minus, fc2_minus_ext, fc3_minus);

        ofs_plus.close();
        ofs_minus.close();
        ifs_orig.close();
    }

}

void Gruneisen::write_newinfo(std::ifstream &ifs, std::ofstream &ofs, const double delta,
                              double ****fc2_new, std::vector<FcsClassExtent> fc2_new_ext, std::vector<FcsClassGru> fc3_new)
{

    int i, j, k;
    std::string line;

    double factor = 1.0 + delta;
    double aa[3][3];

    int nat = system->nat;
    int natmin = system->natmin;
    int nfc2;

    double dummy;

    int ind[3];
    int *pairs, *pairs2;
    double relvec[3];

    std::vector<int> elems;

    std::string str_pair[3];
    int len[3], xyz[3], atom[3];

    std::vector<FcsClassGru>::iterator it_lower;

    memory->allocate(pairs, natmin);
    memory->allocate(pairs2, natmin * 3);

    for (i = 0; i < 2; ++i) {
        std::getline(ifs, line);
        ofs << line << std::endl;
    }

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            ifs >> aa[i][j];
            aa[i][j] *= factor;
        }
    }

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            ofs << std::setw(25) << std::setprecision(16) << aa[i][j];
        }
        ofs << std::endl;
    }
    ifs.ignore();

    for (i = 0; i < nat + 7; ++i) {
        std::getline(ifs, line);
        ofs << line << std::endl;
    }

    ifs >> nfc2;
    ofs << nfc2 << std::endl;

    for (i = 0; i < nfc2; ++i) {
        ifs >> dummy >> ind[0] >> ind[1];

        if (!fcs_phonon->is_fc2_ext) {
            ofs << std::scientific << std::setprecision(16) << std::setw(25) << fc2_new[system->map_s2p[ind[0]/3].atom_num][ind[1]/3][ind[0]%3][ind[1]%3];
        }
        for (j = 0; j < 2; ++j) {
            ofs << std::setw(7) << ind[j];
        }
        ofs << std::endl;
    }
    ifs.ignore();

    for (i = 0; i < 2; ++i) {
        std::getline(ifs, line);
        ofs << line << std::endl;
    }

    while (std::getline(ifs,line)) {

        if (line == "##FORCE CONSTANTS") break;

        if (line == "#LIST_HARMONIC" || line == "#LIST_ANHARM3") {

            ofs << line << std::endl;

            for (i = 0; i < natmin; ++i) ifs >> pairs[i];
            for (i = 0; i < natmin; ++i) ofs << std::setw(6) << pairs[i];
            ofs << std::endl;

            for (i = 0; i < natmin; ++i) {
                for (j = 0; j < pairs[i]; ++j) {
                    ifs >> ind[0] >> ind[1] >> relvec[0] >> relvec[1] >> relvec[2];
                    ofs << std::setw(6) << ind[0] << std::setw(6) << ind[1];
                    for (k = 0; k < 3; ++k) {
                        relvec[k] *= factor;
                        ofs << std::scientific << std::setprecision(16) << std::setw(25) << relvec[k];
                    }
                    ofs << std::endl;
                    ;				}
            }
        }
    }

    ofs << line << std::endl;
    std::getline(ifs, line);
    ofs << line << std::endl;

    while (std::getline(ifs, line)) {

        if (line == "#FCS_HARMONIC") {
            ofs << line << std::endl;
            std::getline(ifs, line);
            ofs << line << std::endl;

            for (i = 0; i < 3 * natmin; ++i) ifs >> pairs2[i];
            for (i = 0; i < 3 * natmin; ++i) ofs << std::setw(6) << pairs2[i];
            ofs << std::endl;

            for (i = 0; i < 3 * natmin; ++i) {
                for (j = 0; j < pairs2[i]; ++j) {
                    ifs >> dummy;
                    ifs >> str_pair[0] >> str_pair[1];

                    for (k = 0; k < 2; ++k) len[k] = str_pair[k].length();

                    for (k = 0; k < 2; ++k) {
                        if (str_pair[k][len[k] - 1] == 'x') {
                            xyz[k] = 0;
                        } else if (str_pair[k][len[k] - 1] == 'y') {
                            xyz[k] = 1;
                        } else if (str_pair[k][len[k] - 1] == 'z') {
                            xyz[k] = 2;
                        } else {
                            error->exit("write_newinfo", "This cannot happen.");
                        }
                    }

                    for (k = 0; k < 2; ++k) {
                        atom[k] = std::atoi(str_pair[k].substr(0,len[k]-1).c_str()) - 1;
                    }
                    atom[0] = system->map_s2p[atom[0]].atom_num;
                    if (!fcs_phonon->is_fc2_ext) {
                        ofs << std::scientific << std::setprecision(16) << std::setw(25) << fc2_new[atom[0]][atom[1]][xyz[0]][xyz[1]] << std::endl;
                    } else {
                        ofs << 0.0 << std::endl;
                    }
                    ofs << std::setw(5) << str_pair[0] << std::setw(5) << str_pair[1] << std::endl;
                }
            }
        }

        if (line == "#FCS_ANHARM3") {
            ofs << line << std::endl;
            std::getline(ifs, line);
            ofs << line << std::endl;

            for (i = 0; i < 3 * natmin; ++i) ifs >> pairs2[i];
            for (i = 0; i < 3 * natmin; ++i) ofs << std::setw(6) << pairs2[i];
            ofs << std::endl;

            for (i = 0; i < 3 * natmin; ++i) {
                for (j = 0; j < pairs2[i]; ++j) {
                    ifs >> dummy;
                    ifs >> str_pair[0] >> str_pair[1] >> str_pair[2];

                    for (k = 0; k < 3; ++k) len[k] = str_pair[k].length();

                    for (k = 0; k < 3; ++k) {
                        if (str_pair[k][len[k] - 1] == 'x') {
                            xyz[k] = 0;
                        } else if (str_pair[k][len[k] - 1] == 'y') {
                            xyz[k] = 1;
                        } else if (str_pair[k][len[k] - 1] == 'z') {
                            xyz[k] = 2;
                        } else {
                            error->exit("write_newinfo", "This cannot happen.");
                        }
                    }

                    elems.clear();
                    for (k = 0; k < 3; ++k) {
                        atom[k] = std::atoi(str_pair[k].substr(0,len[k]-1).c_str()) - 1;
                        ind[k] = 3 * atom[k] + xyz[k];
                        elems.push_back(3 * atom[k] + xyz[k]);
                    }

                    it_lower = lower_bound(fc3_new.begin(), fc3_new.end(), FcsClassGru(3, ind, dummy));

                    if (it_lower == fc3_new.end()) {
                        error->exit("write_newinfo", "FC3 with the same index could not be found.");
                    } else {
                        ofs << std::scientific << std::setprecision(16) << std::setw(25) << (*it_lower).fcs_val << std::endl;
                        ofs << std::setw(5) << str_pair[0] << std::setw(5) << str_pair[1] << std::setw(5) << str_pair[2] << std::endl;
                    }
                }
            }

        }

        if (line == "#FCS_HARMONIC_EXT") {
            ofs << line << std::endl;


            for (i = 0; i < 3 * natmin; ++i) pairs2[i] = 0;
            for (std::vector<FcsClassExtent>::const_iterator it = fc2_new_ext.begin(); it != fc2_new_ext.end(); ++it) {
                pairs2[3 * (*it).atm1 + (*it).xyz1] += 1;
            }

            nfc2 = 0;
            for (i = 0; i < 3 * natmin; ++i) nfc2 += pairs2[i];

            ofs << std::setw(10) << nfc2 << std::endl;
            for (i = 0; i < 3 * natmin; ++i) ofs << std::setw(6) << pairs2[i];
            ofs << std::endl;

            for (std::vector<FcsClassExtent>::const_iterator it = fc2_new_ext.begin(); it != fc2_new_ext.end(); ++it) {
                ofs << std::setw(5) << (*it).atm1 << std::setw(5) << (*it).xyz1;
                ofs << std::setw(8) << (*it).atm2 << std::setw(5) << (*it).xyz2;
                ofs << std::setw(5) << (*it).cell_s;
                ofs << std::scientific << std::setprecision(16) << std::setw(25) << (*it).fcs_val << std::endl;
            }
        }
    }


}


// void Gruneisen::calc_pressure()
// {
//     // Test function for calculating P
// 
//     int i, j, k;
//     int icrd, jcrd;
//     int itran;
//     int natmin = system->natmin;
//     int nat = system->nat;
// 
//     int atom0, atom1;
// 
//     double R1[3], R2[3];
//     double volume;
//     double pressure;
//     double sum = 0.0;
// 
//     for (itran = 0; itran < system->ntran; ++itran) {
// 
//         for (i = 0; i < natmin; ++i) {
//             atom0 = system->map_p2s[i][0];
//             atom1 = system->map_p2s[i][itran];
// 
//             for (k = 0; k < 3; ++k) R1[k] = system->xr_s[atom1][k];
// 
//             rotvec(R1, R1, system->lavec_s);
// 
//             for (j = 0; j < nat; ++j) {
// 
//                 for (k = 0; k < 3; ++k) {
//                     R2[k] = system->xr_s[j][k];
//                     R2[k] += system->xr_s[atom1][k] - system->xr_s[atom0][k];
//                     if (R2[k] >= 1.0) R2[k] -= 1.0;
//                 }
//                 rotvec(R2, R2, system->lavec_s);
// 
//                 for (icrd = 0; icrd < 3; ++icrd) {
//                     for (jcrd = 0; jcrd < 3; ++jcrd) {
//                         sum += fcs_phonon->fc2[i][j][icrd][jcrd] * R1[icrd] * R2[jcrd];
//                     }
//                 }
//             }
// 
//         }
//     }
// 
//     volume = system->volume_p * std::pow(Bohr_in_Angstrom, 3) * 1.0e-30 * static_cast<double>(system->ntran);
//     sum *= Ryd;
//     pressure = - delta_a * sum / (3.0 * volume) * 1.0e-9;
// 
//     std::cout << "Pressure (GPa) = " << pressure << std::endl;
//     std::cout << "Bulk Modulus (GPa) = " << sum / (9.0 * volume) * 1.0e-9 << std::endl;
// }


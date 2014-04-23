/*
 write_phonons.cpp

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include "mpi_common.h"
#include <fstream>
#include <iomanip>
#include <sys/stat.h>
#include "constants.h"
#include "conductivity.h"
#include "dynamical.h"
#include "error.h"
#include "fcs_phonon.h"
#include "gruneisen.h"
#include "kpoint.h"
#include "memory.h"
#include "parsephon.h"
#include "phonon_dos.h"
#include "phonon_thermodynamics.h"
#include "phonon_velocity.h"
#include "relaxation.h"
#include "symmetry_core.h"
#include "system.h"
#include "write_phonons.h"
#include "mathfunctions.h"
#include "isotope.h"
#include "integration.h"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/lexical_cast.hpp>

using namespace PHON_NS;

Writes::Writes(PHON *phon): Pointers(phon){
    Ry_to_kayser = Hz_to_kayser / time_ry;
};

Writes::~Writes(){};

void Writes::write_input_vars()
{
    unsigned int i;

    std::cout << std::endl;
    std::cout << " Input variables :" << std::endl;
    std::cout << " ------------------------------------------------------------" << std::endl;
    std::cout << " General:" << std::endl;
    std::cout << "  PREFIX = " << input->job_title << std::endl;
    std::cout << "  MODE = " << phon->mode << std::endl;
    std::cout << "  FCSXML = " << fcs_phonon->file_fcs << std::endl;
    std::cout << std::endl;

    std::cout << "  NKD = " << system->nkd << "; KD = ";
    for (i = 0; i < system->nkd; ++i) {
        std::cout << std::setw(4) << system->symbol_kd[i];
    }
    std::cout << std::endl;
    std::cout << "  MASS = ";
    for (i = 0; i < system->nkd; ++i) {
        std::cout << std::setw(10) << system->mass_kd[i];
    }
    std::cout << std::endl;
    std::cout << "  NSYM = " << symmetry->nsym << "; TOLERANCE = " << symmetry->tolerance;
    std::cout << "; PRINTSYM = " << symmetry->printsymmetry << std::endl;
    // std::cout << "  TREVSYM = " << symmetry->time_reversal_sym << std::endl;
 
    std::cout << std::endl;

     std::cout << "  NONANALYTIC = " << dynamical->nonanalytic << std::endl;
     if (dynamical->nonanalytic) {
         std::cout << "  BORNINFO = " << dynamical->file_born << "; NA_SIGMA = " << dynamical->na_sigma << std::endl;
     }
    std::cout << std::endl;
    if (writes->nbands >= 0) {
        std::cout << "  NBANDS = " << writes->nbands << std::endl;
    }

    std::cout << "  TMIN = " << system->Tmin << "; TMAX = " << system->Tmax << "; DT = " << system->dT << std::endl;
    std::cout << "  EMIN = " << dos->emin << "; EMAX = " << dos->emax << "; DELTA_E = " << dos->delta_e << std::endl;
    std::cout << std::endl;

    std::cout << "  ISMEAR = " << integration->ismear << "; EPSILON = " << integration->epsilon << std::endl;
    std::cout << std::endl;

    if (phon->mode == "RTA") {
        std::cout << "  RESTART = " << phon->restart_flag << std::endl;
        std::cout << "  TRISYM = " << relaxation->use_triplet_symmetry << std::endl;
        std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << " Kpoint:" << std::endl;
    std::cout << "  KPMODE (1st entry for &kpoint) = " << kpoint->kpoint_mode << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;

    std::cout << " Analysis:" << std::endl;
    if (phon->mode == "PHONONS") {
        std::cout << "  PRINTVEL = " << phonon_velocity->print_velocity << std::endl;
        std::cout << "  PRINTVEC = " << dynamical->print_eigenvectors << std::endl;
        std::cout << "  PRINTXSF = " << writes->writeanime << std::endl;
        std::cout << std::endl;

        if (kpoint->kpoint_mode == 2) {
             std::cout << "  PDOS = " << dos->projected_dos << "; TDOS = " << dos->two_phonon_dos << std::endl;
             std::cout << "  PRINTMSD = " << writes->print_msd << std::endl;
             std::cout << std::endl;
        }
        std::cout << "  GRUNEISEN = " << gruneisen->print_gruneisen << std::endl;
        std::cout << "  NEWFCS = " << gruneisen->print_newfcs;
        if (gruneisen->print_newfcs) {
            std::cout << "; DELTA_A = " << gruneisen->delta_a << std::endl;
            std::cout << "  QUARTIC = " << relaxation->quartic_mode;
        }
        std::cout << std::endl;

    } else if (phon->mode == "RTA") {
        std::cout << "  ISOTOPE = " << isotope->include_isotope << std::endl;
        if (isotope->include_isotope) {
            std::cout << "  ISOFACT = ";
            for (i = 0; i < system->nkd; ++i) {
                std::cout << std::scientific << std::setw(10) << isotope->isotope_factor[i];
            }
            std::cout << std::endl;
        }
        std::cout << "  LCLASSICAL = " << conductivity->use_classical_Cv << std::endl;
        std::cout << std::endl;
        std::cout << "  KS_INPUT = " << relaxation->ks_input << std::endl;
        std::cout << "  QUARTIC = " << relaxation->quartic_mode << std::endl;
        std::cout << "  REALPART = " << relaxation->calc_realpart << std::endl;
        std::cout << "  ATOMPROJ = " << relaxation->atom_project_mode << std::endl;
        std::cout << "  FSTATE_W = " << relaxation->calc_fstate_omega << std::endl;
        std::cout << "  FSTATE_K = " << relaxation->calc_fstate_k << std::endl;
    } else {
        error->exit("write_input_vars", "This cannot happen");
    }
   
    std::cout << std::endl << std::endl;
    std::cout << " ------------------------------------------------------------" << std::endl;
    std::cout << std::endl;
}

void Writes::setup_result_io() 
{

    if (mympi->my_rank == 0) {

        if (phon->restart_flag) {
            // Restart
            fs_result.open(file_result.c_str(), std::ios::in | std::ios::out);
            if (!fs_result) {
                error->exit("setup_result_io", "Could not open file_result");
            }

            // Check the consistency

            std::string line_tmp, str_tmp;
            int natmin_tmp, nkd_tmp;
            int nk_tmp[3], nksym_tmp;
            int is_classical, ismear;
            double epsilon_tmp, T1, T2, delta_T;		


            bool found_tag;

            found_tag = false;
            while (fs_result >> line_tmp)
            {
                if (line_tmp == "#SYSTEM") {
                    found_tag = true;
                    break;
                }
            }
            if (!found_tag) error->exit("setup_result_io", "Could not find #SYSTEM tag");

            fs_result >> natmin_tmp >> nkd_tmp;

            if (!(natmin_tmp == system->natmin && nkd_tmp == system->nkd)) {
                error->exit("setup_result_io", "SYSTEM information is not consistent");
            }

            found_tag = false;
            while (fs_result >> line_tmp)
            {
                if (line_tmp == "#KPOINT") {
                    found_tag = true;
                    break;
                }
            }
            if (!found_tag) error->exit("setup_result_io", "Could not find #KPOINT tag");

            fs_result >> nk_tmp[0] >> nk_tmp[1] >> nk_tmp[2];
            fs_result >> nksym_tmp;

            if (!(kpoint->nkx == nk_tmp[0] && kpoint->nky == nk_tmp[1] && kpoint->nkz == nk_tmp[2]
            && kpoint->nk_reduced == nksym_tmp)) {
                error->exit("setup_result_io", "KPOINT information is not consistent");
            }


            found_tag = false;
            while (fs_result >> line_tmp)
            {
                if (line_tmp == "#CLASSICAL") {
                    found_tag = true;
                    break;
                }
            }
            if (!found_tag) error->exit("setup_result_io", "Could not find #CLASSICAL tag");

            fs_result >> is_classical;
            if (is_classical != conductivity->use_classical_Cv) {
                error->exit("setup_result_io", "CLASSICAL val is not consistent");
            }


            found_tag = false;
            while (fs_result >> line_tmp)
            {
                if (line_tmp == "#FCSXML") {
                    found_tag = true;
                    break;
                }
            }
            if (!found_tag) error->exit("setup_result_io", "Could not find #FCSXML tag");

            fs_result >> str_tmp;
            if (str_tmp != fcs_phonon->file_fcs) {
                error->warn("setup_result_io", "FCSXML is not consistent");
            }


            found_tag = false;
            while (fs_result >> line_tmp)
            {
                if (line_tmp == "#SMEARING") {
                    found_tag = true;
                    break;
                }
            }
            if (!found_tag) error->exit("setup_result_io", "Could not find #SMEARING tag");

            fs_result >> ismear;
            fs_result >> epsilon_tmp;

            if (ismear != integration->ismear) {
                error->exit("setup_result_io", "Smearing method is not consistent");
            }
            if (ismear != -1 && epsilon_tmp != integration->epsilon) {
                error->exit("setup_result_io", "Smearing width is not consistent");
            }


            found_tag = false;
            while (fs_result >> line_tmp)
            {
                if (line_tmp == "#TEMPERATURE") {
                    found_tag = true;
                    break;
                }
            }
            if (!found_tag) error->exit("setup_result_io", "Could not find #TEMPERATURE tag");

            fs_result >> T1 >> T2 >> delta_T;

            if (!(T1 == system->Tmin && T2 == system->Tmax && delta_T == system->dT)) {
                error->exit("setup_result_io", "Temperature information is not consistent");
            }

        } else {
            // From scratch
            fs_result.open(file_result.c_str(), std::ios::out);
            if (!fs_result) {
                error->exit("setup_result_io", "Could not open file_result");
            }

            fs_result << "## General information" << std::endl;
            fs_result << "#SYSTEM" << std::endl;
            fs_result << system->natmin << " " << system->nkd << std::endl;
            fs_result << system->volume_p << std::endl;
            fs_result << "#END SYSTEM" << std::endl;

            fs_result << "#KPOINT" << std::endl;
            fs_result << kpoint->nkx << " " << kpoint->nky << " " << kpoint->nkz << std::endl;
            fs_result << kpoint->nk_reduced << std::endl;

            for (int i = 0; i < kpoint->nk_reduced; ++i){
                fs_result << std::setw(6) << i + 1 << ":";
                for (int j = 0; j < 3; ++j){
                    fs_result << std::setw(15) << std::scientific << kpoint->kpoint_irred_all[i][0].kval[j];
                }
                fs_result << std::setw(12) << std::fixed << kpoint->weight_k[i] << std::endl;
            }
            fs_result.unsetf(std::ios::fixed);

            fs_result << "#END KPOINT" << std::endl;

            fs_result << "#CLASSICAL" << std::endl;
            fs_result << conductivity->use_classical_Cv << std::endl;
            fs_result << "#END CLASSICAL" << std::endl;

            fs_result << "#FCSXML" << std::endl;
            fs_result << fcs_phonon->file_fcs << std::endl;
            fs_result << "#END  FCSXML" << std::endl;

            fs_result << "#SMEARING" << std::endl;
            fs_result << integration->ismear << std::endl;
            fs_result << integration->epsilon << std::endl;
            fs_result << "#END SMEARING" << std::endl;

            fs_result << "#TEMPERATURE" << std::endl;
            fs_result << system->Tmin << " " << system->Tmax << " " << system->dT << std::endl;
            fs_result << "#END TEMPERATURE" << std::endl;

            fs_result << "##END General information" << std::endl;
        }
    }
}

void Writes::write_phonon_info()
{
//     if (nbands < 0 || nbands > 3 * system->natmin) {
//         std::cout << "nbands < 0 or nbands > 3 * natmin" << std::endl;
//         std::cout << "All modes will be printed." << std::endl;    
//         nbands =  3 * system->natmin;
//     }

    if (nbands < 0) {
        nbands = 3 * system->natmin;
    }
 
    std::cout << std::endl;
    std::cout << " ------------------------------------------------------------" << std::endl << std::endl;
    std::cout << " The following files are created: " << std::endl;

    
    if (kpoint->kpoint_mode == 1){
        write_phonon_bands();
    }

    if (phonon_velocity->print_velocity) {
        if (kpoint->kpoint_mode == 1) {
            write_phonon_vel();
        } else if (kpoint->kpoint_mode == 2) {
            write_phonon_vel_all();
        }
    }

    if (dos->flag_dos) {
        write_phonon_dos();

        if (dos->two_phonon_dos) {
            write_two_phonon_dos();
        }

        write_thermodynamics();
        if (print_msd) write_msd();
    }

    if(writeanime) {
        write_mode_anime();
    }

    if (dynamical->print_eigenvectors) {
        write_eigenvectors();
    }

    if (gruneisen->print_gruneisen) {
        write_gruneisen();
    }


}

void Writes::write_phonon_bands()
{
    std::ofstream ofs_bands;
    std::string file_bands = input->job_title + ".bands";

    ofs_bands.open(file_bands.c_str(), std::ios::out);
    if(!ofs_bands) error->exit("write_phonon_bands", "cannot open file_bands");

    unsigned int i, j;

    unsigned int nk = kpoint->nk;

    double *kaxis = kpoint->kaxis;
    double **eval = dynamical->eval_phonon;

    int kcount = 0;

    std::string str_tmp = "NONE";
    std::string str_kpath = "";
    std::string str_kval = "";

    for (i = 0; i < kpoint->kpInp.size(); ++i) {
        if (str_tmp != kpoint->kpInp[i].kpelem[0]) {
            str_tmp = kpoint->kpInp[i].kpelem[0];
            str_kpath += " " + str_tmp;
 
            std::ostringstream ss;
            ss << std::fixed << std::setprecision(6) << kpoint->kaxis[kcount];
            str_kval += " " + ss.str();
        }
        kcount += std::atoi(kpoint->kpInp[i].kpelem[8].c_str());

        if (str_tmp != kpoint->kpInp[i].kpelem[4]) {
            str_tmp = kpoint->kpInp[i].kpelem[4];
            str_kpath += " " + str_tmp;

            std::ostringstream ss;
            ss << std::fixed << std::setprecision(6) << kpoint->kaxis[kcount-1];
            str_kval += " " + ss.str();
        }
    }

    ofs_bands << "# " << str_kpath << std::endl;
    ofs_bands << "#" << str_kval << std::endl;
    ofs_bands << "# k-axis, Eigenvalues [cm^-1]" << std::endl;

    for (i = 0; i < nk; ++i){
        ofs_bands << std::setw(8) << std::fixed << kaxis[i];
        for (j = 0; j < nbands; ++j){
            ofs_bands << std::setw(15) << std::scientific << in_kayser(eval[i][j]);
        }
        ofs_bands << std::endl;
    }

    ofs_bands.close();

    std::cout << "  " <<  std::setw(input->job_title.length() + 12) << std::left << file_bands;
    std::cout << " : Phonon band structure" << std::endl;

}

void Writes::write_phonon_vel()
{
    std::ofstream ofs_vel;
    std::string file_vel = input->job_title + ".phvel";

    ofs_vel.open(file_vel.c_str(), std::ios::out);
    if(!ofs_vel) error->exit("write_phonon_vel", "cannot open file_vel");

    unsigned int i, j;
    unsigned int nk = kpoint->nk;

    double *kaxis = kpoint->kaxis;
    double Ry_to_SI_vel = Bohr_in_Angstrom*1.0e-10/time_ry;

    ofs_vel << "# k-axis, |Velocity| [m / sec]" << std::endl;
    ofs_vel.setf(std::ios::fixed);

    for (i = 0; i < nk; ++i){
        ofs_vel << std::setw(8) << kaxis[i];
        for (j = 0; j < nbands; ++j){
            ofs_vel << std::setw(15) << std::abs(phonon_velocity->phvel[i][j]*Ry_to_SI_vel);
        }
        ofs_vel << std::endl;
    }

    ofs_vel.close();

    std::cout << "  " <<  std::setw(input->job_title.length() + 12) << std::left << file_vel;
    std::cout << " : Phonon velocity along given k path" << std::endl;
}

void Writes::write_phonon_vel_all()
{
    std::ofstream ofs_vel;
    std::string file_vel = input->job_title + ".phvel_all";

    ofs_vel.open(file_vel.c_str(), std::ios::out);
    if(!ofs_vel) error->exit("write_phonon_vel_all", "cannot open file_vel_all");

    unsigned int i, j, k;
    unsigned int nk = kpoint->nk;
    unsigned int ns = dynamical->neval;

    double Ry_to_SI_vel = Bohr_in_Angstrom*1.0e-10/time_ry;
    double **eval = dynamical->eval_phonon;

    ofs_vel << "# Frequency [cm^-1], |Velocity| [m / sec]" << std::endl;
    ofs_vel.setf(std::ios::fixed);

    for (i = 0; i < nk; ++i){

        ofs_vel << "# ik = " << std::setw(8);
        for (j = 0; j < 3; ++j){
            ofs_vel << std::setw(15) << kpoint->xk[i][j];
        }
        ofs_vel << std::endl;

//         phonon_velocity->phonon_vel_k(kpoint->xk[i], vel);
// 
//         for (j = 0; j < ns; ++j){
//             rotvec(vel[j], vel[j], system->lavec_p, 'T');
//             for (k = 0; k < 3; ++k) vel[j][k] /= 2.0 * pi;
//         }

        for (j = 0; j < ns; ++j){
            ofs_vel << std::setw(5) << i;
            ofs_vel << std::setw(5) << j;
            ofs_vel << std::setw(15) << in_kayser(eval[i][j]);
            ofs_vel << std::setw(15) << phonon_velocity->phvel[i][j]*Ry_to_SI_vel;
            ofs_vel << std::endl;
        }
        ofs_vel << std::endl;
    }

    ofs_vel.close();

    std::cout << "  " << std::setw(input->job_title.length() + 12) << std::left << file_vel;
    std::cout << " : Phonon velocity at all k points" << std::endl;
}

void Writes::write_phonon_dos()
{
    int i, iat;
    std::ofstream ofs_dos;
    std::string file_dos = input->job_title + ".dos";

    ofs_dos.open(file_dos.c_str(), std::ios::out);
    if(!ofs_dos) error->exit("write_phonon_dos", "cannot open file_dos");

    ofs_dos << "#";
    for (i = 0; i < system->nkd; ++i) {
        ofs_dos << std::setw(5) << system->symbol_kd[i];
    }
    ofs_dos << std::endl;
    ofs_dos << "#";

    unsigned int *nat_each_kd;
    memory->allocate(nat_each_kd, system->nkd);
    for (i = 0; i < system->nkd; ++i) nat_each_kd[i] = 0;
    for (i = 0; i < system->natmin; ++i) {
        ++nat_each_kd[system->kd[system->map_p2s[i][0]]];
    }
    for (i = 0; i < system->nkd; ++i) {
        ofs_dos << std::setw(5) << nat_each_kd[i];
    }
    ofs_dos << std::endl;
    memory->deallocate(nat_each_kd);

    ofs_dos << "# Energy [cm^-1], TOTAL-DOS";
    if (dynamical->eigenvectors){
        ofs_dos << ", Atom Projected-DOS";   
    }
    ofs_dos << std::endl;
    ofs_dos.setf(std::ios::scientific);

    for (i = 0; i < dos->n_energy; ++i){
        ofs_dos << std::setw(15) << dos->energy_dos[i] << std::setw(15) << dos->dos_phonon[i];
        if(dynamical->eigenvectors) {
            for (iat = 0; iat < system->natmin; ++iat){
                ofs_dos << std::setw(15) << dos->pdos_phonon[iat][i];
            }
        }
        ofs_dos << std::endl;
    } 
    ofs_dos.close();


    std::cout << "  " << std::setw(input->job_title.length() + 12) << std::left << file_dos;

    if (dos->projected_dos) {
        std::cout << " : Phonon DOS and atom projected DOS" << std::endl;
    } else {
        std::cout << " : Phonon DOS" << std::endl;
    }
}

void Writes::write_two_phonon_dos()
{  
    int i, j;

    std::string file_tdos;
    std::ofstream ofs_tdos;

    file_tdos = input->job_title + ".dos2";
    ofs_tdos.open(file_tdos.c_str(), std::ios::out);

    ofs_tdos << "# Two-phonon DOS " << std::endl;
    ofs_tdos << "# delta(e+e1+e2), delta(e-e1-e2), delta(e-e1+e2), delta(e+e1-e2)" << std::endl;
    ofs_tdos << "# Energy [cm^-1], TDOS" << std::endl;

    int n = static_cast<int>((dos->emax * 2.0 - dos->emin) / dos->delta_e);

    for (i = 0; i < n; ++i) {
        ofs_tdos << std::setw(15) << dos->emin + dos->delta_e * static_cast<double>(i);

        for (j = 0; j < 4; ++j) ofs_tdos << std::setw(15) << dos->dos2_phonon[i][j];
        ofs_tdos << std::endl;
    }
    ofs_tdos.close();

    std::cout << "  " <<  std::setw(input->job_title.length() + 12) << std::left << file_tdos;
    std::cout << " : Two-phonon DOS" << std::endl;
}

void Writes::write_mode_anime()
{
    std::ofstream ofs_anime;
    std::string file_anime = input->job_title + ".axsf";

    ofs_anime.open(file_anime.c_str(), std::ios::out);
    if(!ofs_anime) error->exit("write_mode_anime", "cannot open file_anime");

    ofs_anime.setf(std::ios::scientific);

    unsigned int i, j, k;
    unsigned int natmin = system->natmin;
    unsigned int nk = kpoint->nk;

    double force_factor = 100.0;

    double **xmod;
    std::string *kd_tmp;

    memory->allocate(xmod, natmin, 3);
    memory->allocate(kd_tmp, natmin);

    ofs_anime << "ANIMSTEPS " << nbands * nk << std::endl;
    ofs_anime << "CRYSTAL" << std::endl;
    ofs_anime << "PRIMVEC" << std::endl;

    for (i = 0; i < 3; ++i){
        for (j = 0; j < 3; ++j){
            ofs_anime << std::setw(15) << system->lavec_p[j][i]*Bohr_in_Angstrom;
        }
        ofs_anime << std::endl;
    }

    for (i = 0; i < natmin; ++i){
        k = system->map_p2s[i][0];
        for (j = 0; j < 3; ++j){
            xmod[i][j] = system->xc[k][j];
        }

        for (j = 0; j < 3; ++j){
            xmod[i][j] *= Bohr_in_Angstrom;
        }
        kd_tmp[i] = system->symbol_kd[system->kd[k]];
    }

    unsigned int ik, imode;
    double norm;
    std::complex<double> evec_tmp;
    unsigned int m;
    i = 0;

    for (ik = 0; ik < nk; ++ik){
        for (imode = 0; imode < nbands; ++imode){
            ofs_anime << "PRIMCOORD " << std::setw(10) << i + 1 << std::endl;
            ofs_anime << std::setw(10) << natmin << std::setw(10) << 1 << std::endl;
            norm = 0.0;

            for (j = 0; j < 3 * natmin; ++j){
                evec_tmp = dynamical->evec_phonon[ik][imode][j];
                norm += std::pow(evec_tmp.real(), 2) + std::pow(evec_tmp.imag(), 2);
            }

            norm *= force_factor / static_cast<double>(natmin);

            for (j = 0; j < natmin; ++j){

                m = system->map_p2s[j][0];

                ofs_anime << std::setw(10) << kd_tmp[j];

                for (k = 0; k < 3; ++k){
                    ofs_anime << std::setw(15) << xmod[j][k];
                }
                for (k = 0; k < 3; ++k){
                    ofs_anime << std::setw(15) << dynamical->evec_phonon[ik][imode][3 * j + k].real() / (std::sqrt(system->mass[m]) * norm);
                }
                ofs_anime << std::endl;
            }

            ++i;
        }
    }

    memory->deallocate(xmod);
    memory->deallocate(kd_tmp);

    ofs_anime.close();
    std::cout << "  " <<  std::setw(input->job_title.length() + 12) << std::left << file_anime;
    std::cout << " : XcrysDen AXSF file to visualize phonon mode" << std::endl;
}

void Writes::write_eigenvectors()
{
    std::ofstream ofs_evec;
    std::string file_evec = input->job_title + ".evec";

    ofs_evec.open(file_evec.c_str(), std::ios::out);
    if(!ofs_evec) error->exit("write_eigenvectors", "cannot open file_evec");

    ofs_evec.setf(std::ios::scientific);
    unsigned int i, j, k;

    ofs_evec << "Lattice vectors of the primitive lattice" << std::endl;

    for (i = 0; i < 3; ++i){
        for (j = 0; j < 3; ++j){
            ofs_evec << std::setw(15) << system->lavec_p[j][i];
        }
        ofs_evec << std::endl;
    }

    ofs_evec << std::endl;

    for (i = 0; i < 3; ++i){
        for (j = 0; j < 3; ++j){
            ofs_evec << std::setw(15) << system->rlavec_p[i][j];
        }
        ofs_evec << std::endl;
    }

    unsigned int nk = kpoint->nk;
    unsigned int neval = dynamical->neval;

    ofs_evec << "Modes and k-points information below" << std::endl;
    ofs_evec << std::setw(10) << nbands;
    ofs_evec << std::setw(10) << nk << std::endl;

    for (i = 0; i < nk; ++i){
        ofs_evec << "#" << std::setw(10) << i + 1;
        for (j = 0; j < 3; ++j){
            ofs_evec << std::setw(15) << kpoint->xk[i][j];
        }
        ofs_evec << std::endl;
        for (j = 0; j < nbands; ++j){
            ofs_evec << std::setw(15) << dynamical->eval_phonon[i][j] << std::endl;

            for (k = 0; k < neval; ++k){
                ofs_evec << std::setw(15) << real(dynamical->evec_phonon[i][j][k]);
                ofs_evec << std::setw(15) << imag(dynamical->evec_phonon[i][j][k]) << std::endl;
            }
            ofs_evec << std::endl;
        }
        ofs_evec << std::endl;
    }
    ofs_evec.close();

    std::cout << "  " <<  std::setw(input->job_title.length() + 12) << std::left << file_evec;
    std::cout << " : Eigenvector of all k points" << std::endl;
}

double Writes::in_kayser(const double x)
{
    return x * Ry_to_kayser;
}

void Writes::write_thermodynamics()
{
    unsigned int i, NT;
    double Tmin = system->Tmin;
    double Tmax = system->Tmax;
    double dT = system->dT;

    double T, TD;
    std::string file_thermo;

    NT = static_cast<unsigned int>((Tmax - Tmin) / dT);

    std::ofstream ofs_thermo;
    file_thermo = input->job_title + ".thermo";
    ofs_thermo.open(file_thermo.c_str(), std::ios::out);
    if(!ofs_thermo) error->exit("write_thermodynamics", "cannot open file_cv");
    ofs_thermo << "# Temperature [K], Heat capacity / kB, Entropy / kB, Internal energy [Ry], Free energy [Ry]" << std::endl;
//  
//     for (i = 0; i <= NT; ++i) {
//    
//         T = Tmin + dT * static_cast<double>(i);
//         std::cout << " T = " << std::setw(15) << T;
//         TD = 1000.0;
//         phonon_thermodynamics->Debye_T(T, TD);
//         std::cout << "TD = " << TD << std::endl;
//     }
    

    for (i = 0; i <= NT; ++i){
        T = Tmin + dT * static_cast<double>(i);

        ofs_thermo << std::setw(16) << T;
        ofs_thermo << std::setw(18) << phonon_thermodynamics->Cv_tot(T) / k_Boltzmann;
        ofs_thermo << std::setw(18) << phonon_thermodynamics->vibrational_entropy(T) / k_Boltzmann;
        ofs_thermo << std::setw(18) << phonon_thermodynamics->internal_energy(T);
        ofs_thermo << std::setw(18) << phonon_thermodynamics->free_energy(T) << std::endl;
    }

    ofs_thermo.close();

    std::cout << "  " <<  std::setw(input->job_title.length() + 12) << std::left << file_thermo;
    std::cout << " : Thermodynamic quantities" << std::endl;
}

void Writes::write_gruneisen()
{

    if (kpoint->kpoint_mode == 1) {
        if (nbands < 0 || nbands > 3 * system->natmin) {
            nbands =  3 * system->natmin;
        }

        std::ofstream ofs_gruneisen;

        std::string file_gru = input->job_title + ".gruneisen";
        ofs_gruneisen.open(file_gru.c_str(), std::ios::out);
        if(!ofs_gruneisen) error->exit("write_gruneisen", "cannot open file_vel");

        unsigned int i, j;
        unsigned int nk = kpoint->nk;

        double *kaxis = kpoint->kaxis;

        ofs_gruneisen << "# k-axis, gamma" << std::endl;
        ofs_gruneisen.setf(std::ios::fixed);

        for (i = 0; i < nk; ++i){
            ofs_gruneisen << std::setw(8) << kaxis[i];
            for (j = 0; j < nbands; ++j){
                ofs_gruneisen << std::setw(15) << gruneisen->gruneisen[i][j].real();
            }
            ofs_gruneisen << std::endl;
        }

        ofs_gruneisen.close();

        std::cout << "  " <<  std::setw(input->job_title.length() + 12) << std::left << file_gru;
        std::cout << " : Gruneisen parameters along given k-path" << std::endl;

    } else {

        std::ofstream ofs_gruall;
        std::string file_gruall;
        file_gruall = input->job_title + ".gru_all";
        ofs_gruall.open(file_gruall.c_str(), std::ios::out);
        if (!ofs_gruall) error->exit("write_gruneisen", "cannot open file_gruall");

        unsigned int i, j, k;
        unsigned int nk = kpoint->nk;
        unsigned int ns = dynamical->neval;

        ofs_gruall << "# knum, snum, omega [cm^-1], gruneisen parameter" << std::endl;

        for (i = 0; i < nk; ++i){
            ofs_gruall << "# knum = " << i;
            for (k = 0; k < 3; ++k) {
                ofs_gruall << std::setw(15)<< kpoint->xk[i][k];
            }
            ofs_gruall << std::endl;

            for (j = 0; j < ns; ++j){
                ofs_gruall << std::setw(5) << i;
                ofs_gruall << std::setw(5) << j;
                ofs_gruall << std::setw(15) << in_kayser(dynamical->eval_phonon[i][j]);
                ofs_gruall << std::setw(15) << gruneisen->gruneisen[i][j].real();
                ofs_gruall << std::endl;
            }
        }

        ofs_gruall.close();


        std::cout << "  " <<  std::setw(input->job_title.length() + 12) << std::left << file_gruall;
        std::cout << " : Gruneisen parameters at all k points" << std::endl;
    }
}

void Writes::write_msd()
{
    // Write room mean square displacement of atoms

    std::string file_rmsd = input->job_title + ".msd";
    std::ofstream ofs_rmsd;

    unsigned int i, j;
    unsigned int NT;
    unsigned int ns = dynamical->neval;

    double Tmin = system->Tmin;
    double Tmax = system->Tmax;
    double dT = system->dT;
    double T, d2_tmp;

    ofs_rmsd.open(file_rmsd.c_str(), std::ios::out);
    if (!ofs_rmsd) error->exit("write_rmsd", "Could not open file_rmsd");

    ofs_rmsd << "# Mean Square Displacements at a function of temperature." << std::endl;
    ofs_rmsd << "# Temperature [K], <(u_{1}^{x})^{2}>, <(u_{1}^{y})^{2}>, <(u_{1}^{z})^{2}>, .... [Angstrom^2]" << std::endl;

    NT = static_cast<unsigned int>((Tmax - Tmin) / dT);

    for (i = 0; i < NT; ++i) {

        T = Tmin + static_cast<double>(i) * dT;
        ofs_rmsd << std::setw(15) << T;

        for (j = 0; j < ns; ++j){
            d2_tmp = phonon_thermodynamics->disp2_avg(T, j, j);
            ofs_rmsd << std::setw(15) << std::sqrt(d2_tmp)*Bohr_in_Angstrom;
        }
        ofs_rmsd << std::endl;
    }
    ofs_rmsd.close();

    std::cout << "  " <<  std::setw(input->job_title.length() + 12) << std::left << file_rmsd;
    std::cout << " : Mean-square-displacement (MSD)" << std::endl;

}

void Writes::write_kappa()
{
    // Write lattice thermal conductivity

    if (mympi->my_rank == 0) {
        int i, j, k;

        std::string file_kappa = input->job_title + ".kl";
        std::ofstream ofs_kl;

        ofs_kl.open(file_kappa.c_str(), std::ios::out);
        if (!ofs_kl) error->exit("write_kappa", "Could not open file_kappa");

        ofs_kl << "# Temperature [K], Thermal Conductivity (xx, xy, xz, yx, yy, yz, zx, zy, zz) [W/mK]" << std::endl;

        if (isotope->include_isotope) {
            ofs_kl << "# Isotope effects are included." << std::endl;
        }

        for (i = 0; i < conductivity->ntemp; ++i) {
            ofs_kl << std::setw(10) << std::right << std::fixed << std::setprecision(2)
                << conductivity->Temperature[i];
            for (j = 0; j < 3; ++j){
                for (k = 0; k < 3; ++k){
                    ofs_kl << std::setw(15) << std::fixed 
                        << std::setprecision(4) << conductivity->kappa[i][j][k];
                }
            }
            ofs_kl << std::endl;
        }
        ofs_kl.close();

        std::cout << std::endl;
        std::cout << " Lattice thermal conductivity is store in the file " << file_kappa << std::endl;
    }
}

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

using namespace PHON_NS;

Writes::Writes(PHON *phon): Pointers(phon){
    Ry_to_kayser = Hz_to_kayser / time_ry;
};

Writes::~Writes(){};

void Writes::write_input_vars()
{
    unsigned int i;

    std::cout << std::endl;
    std::cout << "Input variables below:" << std::endl;
    std::cout << "---------------------------------------------------" << std::endl;
    std::cout << "General:" << std::endl;
    std::cout << " PREFIX = " << input->job_title << std::endl;
    std::cout << " NSYM = " << symmetry->nsym << "; TOLERANCE = " << symmetry->tolerance << std::endl;
    std::cout << " PRINTSYMM = " << symmetry->printsymmetry << std::endl;
    std::cout << " TREVSYM = " << symmetry->time_reversal_sym << std::endl;
    std::cout << " CELLDIM = ";
    for (i = 0; i < 3; ++i) std::cout << std::setw(4) << system->cell_dimension[i];
    std::cout << std::endl << std::endl;

    std::cout << " MODE = " << phon->mode << std::endl;
    std::cout << " FCSINFO = " << fcs_phonon->file_fcs << std::endl;
    std::cout << std::endl;

    std::cout << " EIGENVECTOR = " << dynamical->eigenvectors << std::endl;
    std::cout << " PRINTXSF = " << writes->writeanime << "; NBANDS = " << writes->nbands << std::endl;
    std::cout << " TMIN = " << system->Tmin << "; TMAX = " << system->Tmax << "; DT = " << system->dT << std::endl;
    std::cout << " NONANALYTIC = " << dynamical->nonanalytic << "; BORNINFO = " << dynamical->file_born << "; NA_SIGMA = " << dynamical->na_sigma << std::endl;
    std::cout << " EMIN = " << dos->emin << "; EMAX = " << dos->emax << "; DELTA_E = " << dos->delta_e << std::endl;
    std::cout << std::endl;

    std::cout << " DELTA_A = " << gruneisen->delta_a << std::endl;
    std::cout << std::endl;

    std::cout << " RESTART = " << phon->restart_flag << std::endl;
    std::cout << " ISMEAR = " << relaxation->ksum_mode << "; EPSILON = " << relaxation->epsilon << std::endl;
    std::cout << " LCLASSICAL = " << conductivity->use_classical_Cv << std::endl;
    std::cout << " KS_INPUT = " << relaxation->ks_input << "; QUARTIC = " << relaxation->quartic_mode << std::endl;
    std::cout << " ATOMPROJ = " << relaxation->atom_project_mode << "; REALPART = " << relaxation->calc_realpart << std::endl;

    std::cout << std::endl << std::endl;

    std::cout << "Kpoint:" << std::endl;
    std::cout << " KPMODE (1st entry for &kpoint) = " << kpoint->kpoint_mode << std::endl;
    std::cout << std::endl;
    std::cout << "---------------------------------------------------" << std::endl;
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
            && kpoint->nk_equiv.size() == nksym_tmp)) {
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
                if (line_tmp == "#FCSINFO") {
                    found_tag = true;
                    break;
                }
            }
            if (!found_tag) error->exit("setup_result_io", "Could not find #FCSINFO tag");

            fs_result >> str_tmp;
            if (str_tmp != fcs_phonon->file_fcs) {
                error->warn("setup_result_io", "FCSINFO is not consistent");
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

            if (ismear != relaxation->ksum_mode) {
                error->exit("setup_result_io", "Smearing method is not consistent");
            }
            if (ismear != -1 && epsilon_tmp != relaxation->epsilon) {
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
            fs_result << kpoint->nk_equiv.size() << std::endl;

            int ik = 0;
            for (int i = 0; i < kpoint->nk_equiv.size(); ++i){
                fs_result << std::setw(6) << i + 1 << ":";
                for (int j = 0; j < 3; ++j){
                    fs_result << std::setw(15) << std::scientific << kpoint->kpIBZ[ik].kval[j];
                }
                fs_result << std::setw(12) << std::fixed << kpoint->weight_k[i] << std::endl;
                ik += kpoint->nk_equiv[i];
            }
            fs_result.unsetf(std::ios::fixed);

            fs_result << "#END KPOINT" << std::endl;

            fs_result << "#CLASSICAL" << std::endl;
            fs_result << conductivity->use_classical_Cv << std::endl;
            fs_result << "#END CLASSICAL" << std::endl;

            fs_result << "#FCSINFO" << std::endl;
            fs_result << fcs_phonon->file_fcs << std::endl;
            fs_result << "#END  FCSINFO" << std::endl;

            fs_result << "#SMEARING" << std::endl;
            fs_result << relaxation->ksum_mode << std::endl;
            fs_result << relaxation->epsilon << std::endl;
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
    if (nbands < 0 || nbands > 3 * system->natmin) {
        std::cout << "nbands < 0 or nbands > 3 * natmin" << std::endl;
        std::cout << "All modes will be printed." << std::endl;    
        nbands =  3 * system->natmin;
    }

    if(kpoint->kpoint_mode == 1){
        write_phonon_bands();
        write_phonon_vel();
    }

    if(dos->flag_dos) {
        write_phonon_dos();
        write_thermodynamics();
        write_phonon_vel_all();
    }

    if(writeanime) {
        write_mode_anime();
    }

    if(dynamical->eigenvectors) {
        write_eigenvectors();
        write_rmsd();
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

    ofs_bands << "# k-axis, Eigenvalues [cm^-1]" << std::endl;

    for (i = 0; i < nk; ++i){
        ofs_bands << std::setw(8) << std::fixed << kaxis[i];
        for (j = 0; j < nbands; ++j){
            ofs_bands << std::setw(15) << std::scientific << in_kayser(eval[i][j]);
        }
        ofs_bands << std::endl;
    }

    ofs_bands.close();
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
    double **vel;

    memory->allocate(vel, ns, 3);

    ofs_vel << "# Frequency [cm^-1], |Velocity| [m / sec]" << std::endl;
    ofs_vel.setf(std::ios::fixed);

    for (i = 0; i < nk; ++i){

        ofs_vel << "# ik = " << std::setw(8);
        for (j = 0; j < 3; ++j){
            ofs_vel << std::setw(15) << kpoint->xk[i][j];
        }
        ofs_vel << std::endl;

        phonon_velocity->phonon_vel_k(kpoint->xk[i], vel);

        for (j = 0; j < ns; ++j){
            rotvec(vel[j], vel[j], system->lavec_p, 'T');
            for (k = 0; k < 3; ++k) vel[j][k] /= 2.0 * pi;
        }

        for (j = 0; j < ns; ++j){
            ofs_vel << std::setw(5) << i;
            ofs_vel << std::setw(5) << j;
            ofs_vel << std::setw(15) << in_kayser(eval[i][j]);
            ofs_vel << std::setw(15) << std::sqrt(std::pow(vel[j][0], 2) + std::pow(vel[j][1], 2) + std::pow(vel[j][2], 2))*Ry_to_SI_vel;
            ofs_vel << std::endl;
        }
        ofs_vel << std::endl;
    }

    ofs_vel.close();

    memory->deallocate(vel);
}

void Writes::write_phonon_dos()
{
    int i, iat;
    std::ofstream ofs_dos;
    std::string file_bands = input->job_title + ".dos";

    ofs_dos.open(file_bands.c_str(), std::ios::out);
    if(!ofs_dos) error->exit("write_phonon_dos", "cannot open file_dos");

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

    std::cout << std::endl << "Total DOS ";
    if(dynamical->eigenvectors) {
        std::cout << "and atom projected-DOS ";
    }
    std::cout << "are printed in the file: " << file_bands << std::endl << std::endl;
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
    ofs_thermo << "# Temperature [K], Internal Energy [Ry], Heat Capacity / kB" << std::endl;

    TD = 1000.0;
    phonon_thermodynamics->Debye_T(Tmax, TD);
    std::cout << "TD = " << TD << std::endl;

    for (i = 0; i <= NT; ++i){
        T = Tmin + dT * static_cast<double>(i);

        ofs_thermo << std::setw(15) << T;
        ofs_thermo << std::setw(15) << phonon_thermodynamics->Internal_Energy(T);
        ofs_thermo << std::setw(15) << phonon_thermodynamics->Cv_tot(T) / k_Boltzmann << std::endl;
    }

    ofs_thermo.close();
}

void Writes::write_gruneisen()
{

    if (kpoint->kpoint_mode == 1) {
        if (nbands < 0 || nbands > 3 * system->natmin) {
            std::cout << "WARNING: nbands < 0 or nbands > 3 * natmin" << std::endl;
            std::cout << "All modes will be printed." << std::endl;    
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
    }
}

void Writes::write_rmsd()
{
    // Write room mean square displacement of atoms

    std::string file_rmsd = input->job_title + ".rmsd";
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

    ofs_rmsd << "# Root Mean Square Displacements at a function of temperature." << std::endl;
    ofs_rmsd << "# Temperature [K], u1_x, u1_y, u1_z, .... [Angstrom]" << std::endl;

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
            ofs_kl << std::setw(5) << conductivity->Temperature[i];
            for (j = 0; j < 3; ++j){
                for (k = 0; k < 3; ++k){
                    ofs_kl << std::setw(15) << conductivity->kappa[i][j][k];
                }
            }
            ofs_kl << std::endl;
        }
        ofs_kl.close();

        std::cout << std::endl;
        std::cout << "Lattice thermal conductivity is store in the file " << file_kappa << std::endl;
    }
}

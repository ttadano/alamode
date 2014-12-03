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
#include <boost/lexical_cast.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/version.hpp>

using namespace PHON_NS;

Gruneisen::Gruneisen(PHON *phon): Pointers(phon){};
Gruneisen::~Gruneisen(){};

void Gruneisen::setup()
{
    MPI_Bcast(&delta_a, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&print_newfcs, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD);

    if (print_gruneisen || print_newfcs) {
        int i, icell, ix, iy, iz;

        memory->allocate(xshift_s, 27, 3);

        for (i = 0; i < 3; ++i) xshift_s[0][i] = 0.0;

        icell =0;

        for (ix = -1; ix <= 1; ++ix) {
            for (iy = -1; iy <= 1; ++iy) {
                for (iz = -1; iz <= 1; ++iz) {
                    if (ix == 0 && iy == 0 && iz == 0) continue;

                    ++icell;

                    xshift_s[icell][0] = static_cast<double>(ix);
                    xshift_s[icell][1] = static_cast<double>(iy);
                    xshift_s[icell][2] = static_cast<double>(iz);
                }
            }
        }
        prepare_delta_fcs(fcs_phonon->force_constant_with_cell[1], delta_fc2);
    }

    if (print_newfcs && relaxation->quartic_mode > 0) {
        prepare_delta_fcs(fcs_phonon->force_constant_with_cell[2], delta_fc3);
    }
    if (print_gruneisen) {
        memory->allocate(gruneisen, kpoint->nk, dynamical->neval);
    }

    if (mympi->my_rank == 0) {
        if (print_newfcs) {
            std::cout << std::endl;
            if (relaxation->quartic_mode > 0) {
                std::cout << " NEWFCS = 1 : Harmonic and cubic force constants of " << std::endl;
            }
            else {
                std::cout << " NEWFCS = 1 : Harmonic force constants of " << std::endl;
            }
            std::cout << "              expanded/compressed systems will be estimated" << std::endl;
            std::cout << "              with DELTA_A = " << std::setw(5) << delta_a << std::endl;
        }
    }
}

void Gruneisen::finish_gruneisen()
{
    if (print_gruneisen) memory->deallocate(gruneisen);

    if (print_gruneisen || print_newfcs) {

        memory->deallocate(xshift_s);

        delta_fc2.clear();
        delta_fc3.clear();

    }
}

void Gruneisen::calc_gruneisen()
{
    unsigned int is, ik;
    unsigned int i, j;
    unsigned int ns = dynamical->neval;
    unsigned int nk = kpoint->nk;
    double gamma_imag;
    std::complex<double> **dfc2_reciprocal;

    memory->allocate(dfc2_reciprocal, ns, ns);

    if (mympi->my_rank == 0) {
        std::cout << std::endl;
        std::cout << " GRUNEISEN = 1 : Calculating Gruneisen parameters ... ";
    }

    for (ik = 0; ik < nk; ++ik){
        for (is = 0; is < ns; ++is){

            calc_dfc2_reciprocal(dfc2_reciprocal, kpoint->xk[ik]);

            gruneisen[ik][is] = std::complex<double>(0.0, 0.0);

            for (i = 0; i < ns; ++i){
                for (j = 0; j < ns; ++j){
                    gruneisen[ik][is] += std::conj(dynamical->evec_phonon[ik][is][i])
                                       * dfc2_reciprocal[i][j] 
                                       * dynamical->evec_phonon[ik][is][j];
                }
            }

            gamma_imag = gruneisen[ik][is].imag();
            if (std::abs(gamma_imag) > eps10) {
                error->warn("calc_gruneisen", "Gruneisen parameter is not real");
            }

            if (std::abs(dynamical->eval_phonon[ik][is]) < eps8) {
                gruneisen[ik][is] = 0.0;
            } else {
                gruneisen[ik][is] /= - 6.0 * std::pow(dynamical->eval_phonon[ik][is], 2);
            }
        }
    }
    memory->deallocate(dfc2_reciprocal);

    if (mympi->my_rank == 0) {
        std::cout << " done !" << std::endl;
    }
}

void Gruneisen::calc_dfc2_reciprocal(std::complex<double> **dphi2, double *xk_in) {

    unsigned int i, j;
    unsigned int ns = dynamical->neval;

    unsigned int atm1, atm2, xyz1, xyz2;
    unsigned int atm1_s, atm2_s;
    unsigned int tran, cell_s;

    double vec[3];
    double phase;

    std::complex<double> im(0.0, 1.0);


    for (i = 0; i < ns; ++i){
        for (j = 0; j < ns; ++j){
            dphi2[i][j] = std::complex<double>(0.0, 0.0);
        }
    }

    for (std::vector<FcsArrayWithCell>::const_iterator it = delta_fc2.begin(); it != delta_fc2.end(); ++it) {

        atm1 = (*it).pairs[0].index / 3;
        xyz1 = (*it).pairs[0].index % 3;
        atm2 = (*it).pairs[1].index / 3;
        xyz2 = (*it).pairs[1].index % 3;

        tran = (*it).pairs[1].tran;
        cell_s = (*it).pairs[1].cell_s;

        atm1_s = system->map_p2s[atm1][0];
        atm2_s = system->map_p2s[atm2][tran];

        for (i = 0; i < 3; ++i) {
            vec[i] = system->xr_s[atm2_s][i] + xshift_s[cell_s][i] 
            - system->xr_s[system->map_p2s[atm2][0]][i];
        }

        rotvec(vec, vec, system->lavec_s);
        rotvec(vec, vec, system->rlavec_p);

        phase = vec[0] * xk_in[0] + vec[1] * xk_in[1] + vec[2] * xk_in[2];

        dphi2[3 * atm1 + xyz1][3 * atm2 + xyz2] 
        += (*it).fcs_val * std::exp(-im * phase) / std::sqrt(system->mass[atm1_s] * system->mass[atm2_s]);
    }

}

void Gruneisen::prepare_delta_fcs(const std::vector<FcsArrayWithCell> fcs_in, std::vector<FcsArrayWithCell> &delta_fcs)
{
    unsigned int i;

    double vec[3];
    double fcs_tmp = 0.0;

    std::vector<FcsAlignedForGruneisen> fcs_aligned;

    std::vector<AtomCellSuper> pairs_vec;
    std::vector<int> arr_old, arr_tmp;

    unsigned int norder = fcs_in[0].pairs.size();
    unsigned int nelems;

    delta_fcs.clear();
    fcs_aligned.clear();

    for (std::vector<FcsArrayWithCell>::const_iterator it  = fcs_in.begin(); it != fcs_in.end(); ++it) {
            fcs_aligned.push_back(FcsAlignedForGruneisen((*it).fcs_val, (*it).pairs));
    }

    std::sort(fcs_aligned.begin(), fcs_aligned.end());

    delta_fcs.clear();
    arr_old.clear();

    nelems = 4 * (norder - 2) + 2;

    for (i = 0; i < nelems; ++i) arr_old.push_back(-1);

    for (std::vector<FcsAlignedForGruneisen>::const_iterator it = fcs_aligned.begin(); it != fcs_aligned.end(); ++it) {

        arr_tmp.clear();

        arr_tmp.push_back((*it).pairs[0].index / 3);
        arr_tmp.push_back((*it).pairs[0].index % 3);

        for (i = 1; i < norder - 1; ++i) {
            arr_tmp.push_back((*it).pairs[i].index / 3);
            arr_tmp.push_back((*it).pairs[i].tran);
            arr_tmp.push_back((*it).pairs[i].cell_s);
            arr_tmp.push_back((*it).pairs[i].index % 3);
        }
       

        if (arr_tmp != arr_old) {

            if (arr_old[0] != -1) { // Neglect the initial entry
                delta_fcs.push_back(FcsArrayWithCell(fcs_tmp, pairs_vec));
            }

            pairs_vec.clear();
            for (i = 0; i < norder - 1; ++i) {
                pairs_vec.push_back((*it).pairs[i]);
            }           

            fcs_tmp = 0.0;
            arr_old.clear();
            arr_old.reserve(arr_tmp.size());
            std::copy(arr_tmp.begin(), arr_tmp.end(), std::back_inserter(arr_old));
        }

        for (i = 0; i < 3; ++i) {
            vec[i] = system->xr_s[system->map_p2s[(*it).pairs[norder - 1].index / 3][(*it).pairs[norder - 1].tran]][i]
                   + xshift_s[(*it).pairs[norder - 1].cell_s][i];
        }

        rotvec(vec, vec, system->lavec_s);

        fcs_tmp += (*it).fcs_val * vec[(*it).pairs[norder - 1].index % 3];
    }

    delta_fcs.push_back(FcsArrayWithCell(fcs_tmp, pairs_vec));

    fcs_aligned.clear();
}

void Gruneisen::write_new_fcsxml_all()
{
    std::cout << std::endl;
    std::cout << " NEWFCS = 1 : Following XML files are created. " << std::endl;

    std::string file_xml;

    file_xml = input->job_title + "_+.xml";
    write_new_fcsxml(file_xml, delta_a);

    std::cout << "  " <<  std::setw(input->job_title.length() + 12) << std::left << file_xml;
    std::cout << " : Force constants of the system expanded by " 
        << std::fixed << std::setprecision(3) << delta_a * 100 << " %" << std::endl;

    file_xml = input->job_title + "_-.xml";
    write_new_fcsxml(file_xml, -delta_a);

    std::cout << "  " <<  std::setw(input->job_title.length() + 12) << std::left << file_xml;
    std::cout << " : Force constants of the system compressed by "
        << std::fixed << std::setprecision(3) << delta_a * 100 << " %" << std::endl;
}

void Gruneisen::write_new_fcsxml(const std::string filename_xml, const double change_ratio_of_a)
{
    int i, j;
    double lattice_vector[3][3];

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            lattice_vector[i][j] = (1.0 + change_ratio_of_a) * system->lavec_s[i][j];
        }
    }

    using boost::property_tree::ptree;

    ptree pt;
    std::string str_pos[3];

    pt.put("Data.ANPHON_version", "0.9.2");
    pt.put("Data.Description.OriginalXML", fcs_phonon->file_fcs);
    pt.put("Data.Description.Delta_A", double2string(change_ratio_of_a));

    pt.put("Data.Structure.NumberOfAtoms", system->nat);
    pt.put("Data.Structure.NumberOfElements", system->nkd);

    for (i = 0; i < system->nkd; ++i) { 
        ptree &child = pt.add("Data.Structure.AtomicElements.element", system->symbol_kd[i]);
        child.put("<xmlattr>.number", i + 1);
    }

    for (i = 0; i < 3; ++i) {
        str_pos[i].clear();
        for (j = 0; j < 3; ++j) {
            str_pos[i] += " " + double2string(lattice_vector[j][i]);
        }
    }
    pt.put("Data.Structure.LatticeVector", "");
    pt.put("Data.Structure.LatticeVector.a1", str_pos[0]);
    pt.put("Data.Structure.LatticeVector.a2", str_pos[1]);
    pt.put("Data.Structure.LatticeVector.a3", str_pos[2]);

    pt.put("Data.Structure.Position", "");
    std::string str_tmp;

    for (i = 0; i < system->nat; ++i) {
        str_tmp.clear();
        for (j = 0; j < 3; ++j) str_tmp += " " + double2string(system->xr_s[i][j]);
        ptree &child = pt.add("Data.Structure.Position.pos", str_tmp);
        child.put("<xmlattr>.index", i + 1);
        child.put("<xmlattr>.element", system->symbol_kd[system->kd[i]]);
    }

    pt.put("Data.Symmetry.NumberOfTranslations", system->ntran);
    for (i = 0; i < system->ntran; ++i) {
        for (j = 0; j < system->natmin; ++j) {
            ptree &child = pt.add("Data.Symmetry.Translations.map", system->map_p2s[j][i] + 1);
            child.put("<xmlattr>.tran", i + 1);
            child.put("<xmlattr>.atom", j + 1);
        }
    }

    pt.put("Data.ForceConstants", "");
    str_tmp.clear();

    for (std::vector<FcsArrayWithCell>::const_iterator it = fcs_phonon->force_constant_with_cell[0].begin();
        it != fcs_phonon->force_constant_with_cell[0].end(); ++it) {

            ptree &child = pt.add("Data.ForceConstants.HARMONIC.FC2", double2string((*it).fcs_val));

            child.put("<xmlattr>.pair1", boost::lexical_cast<std::string>((*it).pairs[0].index / 3 + 1)
                + " " + boost::lexical_cast<std::string>((*it).pairs[0].index %3 + 1));
            child.put("<xmlattr>.pair2", boost::lexical_cast<std::string>(system->map_p2s[(*it).pairs[1].index / 3][(*it).pairs[1].tran] + 1) 
                + " " + boost::lexical_cast<std::string>((*it).pairs[1].index % 3 + 1)
                + " " + boost::lexical_cast<std::string>((*it).pairs[1].cell_s + 1));
    }

    for (std::vector<FcsArrayWithCell>::const_iterator it = delta_fc2.begin(); it != delta_fc2.end(); ++it) {

        if (std::abs((*it).fcs_val) < eps12) continue;

        ptree &child = pt.add("Data.ForceConstants.HARMONIC.FC2", double2string(change_ratio_of_a * (*it).fcs_val));

        child.put("<xmlattr>.pair1", boost::lexical_cast<std::string>((*it).pairs[0].index / 3 + 1)
            + " " + boost::lexical_cast<std::string>((*it).pairs[0].index %3 + 1));
        child.put("<xmlattr>.pair2", boost::lexical_cast<std::string>(system->map_p2s[(*it).pairs[1].index / 3][(*it).pairs[1].tran] + 1) 
            + " " + boost::lexical_cast<std::string>((*it).pairs[1].index % 3 + 1)
            + " " + boost::lexical_cast<std::string>((*it).pairs[1].cell_s + 1));
    }

    if (relaxation->quartic_mode) {
        for (std::vector<FcsArrayWithCell>::const_iterator it = fcs_phonon->force_constant_with_cell[1].begin();
            it != fcs_phonon->force_constant_with_cell[1].end(); ++it) {

                if ((*it).pairs[1].index > (*it).pairs[2].index) continue;

                ptree &child = pt.add("Data.ForceConstants.ANHARM3.FC3", double2string((*it).fcs_val));

                child.put("<xmlattr>.pair1", boost::lexical_cast<std::string>((*it).pairs[0].index / 3 + 1)
                    + " " + boost::lexical_cast<std::string>((*it).pairs[0].index %3 + 1));
                child.put("<xmlattr>.pair2", boost::lexical_cast<std::string>(system->map_p2s[(*it).pairs[1].index / 3][(*it).pairs[1].tran] + 1) 
                    + " " + boost::lexical_cast<std::string>((*it).pairs[1].index % 3 + 1)
                    + " " + boost::lexical_cast<std::string>((*it).pairs[1].cell_s + 1));
                child.put("<xmlattr>.pair3", boost::lexical_cast<std::string>(system->map_p2s[(*it).pairs[2].index / 3][(*it).pairs[2].tran] + 1) 
                    + " " + boost::lexical_cast<std::string>((*it).pairs[2].index % 3 + 1)
                    + " " + boost::lexical_cast<std::string>((*it).pairs[2].cell_s + 1));
        }

        for (std::vector<FcsArrayWithCell>::const_iterator it = delta_fc3.begin(); it != delta_fc3.end(); ++it) {

                if (std::abs((*it).fcs_val) < eps12) continue;

                if ((*it).pairs[1].index > (*it).pairs[2].index) continue;

                ptree &child = pt.add("Data.ForceConstants.ANHARM3.FC3", double2string(change_ratio_of_a * (*it).fcs_val));

                child.put("<xmlattr>.pair1", boost::lexical_cast<std::string>((*it).pairs[0].index / 3 + 1)
                    + " " + boost::lexical_cast<std::string>((*it).pairs[0].index %3 + 1));
                child.put("<xmlattr>.pair2", boost::lexical_cast<std::string>(system->map_p2s[(*it).pairs[1].index / 3][(*it).pairs[1].tran] + 1) 
                    + " " + boost::lexical_cast<std::string>((*it).pairs[1].index % 3 + 1)
                    + " " + boost::lexical_cast<std::string>((*it).pairs[1].cell_s + 1));
                child.put("<xmlattr>.pair3", boost::lexical_cast<std::string>(system->map_p2s[(*it).pairs[2].index / 3][(*it).pairs[2].tran] + 1) 
                    + " " + boost::lexical_cast<std::string>((*it).pairs[2].index % 3 + 1)
                    + " " + boost::lexical_cast<std::string>((*it).pairs[2].cell_s + 1));
        }
    }

    using namespace boost::property_tree::xml_parser;
    const int indent = 2;

#if BOOST_VERSION >= 105600
    write_xml(filename_xml, pt, std::locale(),
	      xml_writer_make_settings<ptree::key_type>(' ', indent, widen<std::string>("utf-8")));
#else
    write_xml(filename_xml, pt, std::locale(),
	      xml_writer_make_settings(' ', indent, widen<char>("utf-8")));
#endif
}


std::string Gruneisen::double2string(const double d){

    std::string rt;
    std::stringstream ss;

    ss << std::scientific << std::setprecision(15) << d;
    ss >> rt;
    return rt;
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



// void Gruneisen::prepare_newfc2_mod()
// {
//     FcsClassExtent dfc2_ext_tmp;
// 
//     fc2_plus_ext.clear();
//     fc2_minus_ext.clear();
// 
//     for (std::vector<FcsClassExtent>::const_iterator it  = fcs_phonon->fc2_ext.begin(); 
//                                                      it != fcs_phonon->fc2_ext.end(); ++it) {
//         fc2_plus_ext.push_back(*it);
//         fc2_minus_ext.push_back(*it);
//     }
// 
//     for (std::vector<FcsArrayWithCell>::const_iterator it = delta_fc2.begin(); it != delta_fc2.end(); ++it) {
// 
//         dfc2_ext_tmp.atm1 = (*it).pairs[0].index / 3;
//         dfc2_ext_tmp.atm2 = system->map_p2s[(*it).pairs[1].index / 3][(*it).pairs[1].tran];
//         dfc2_ext_tmp.cell_s = (*it).pairs[1].cell_s;
//         dfc2_ext_tmp.xyz1 = (*it).pairs[0].index % 3;
//         dfc2_ext_tmp.xyz2 = (*it).pairs[1].index % 3;
// 
//         dfc2_ext_tmp.fcs_val = delta_a * (*it).fcs_val;
//         fc2_plus_ext.push_back(dfc2_ext_tmp);
// 
//         dfc2_ext_tmp.fcs_val = -delta_a * (*it).fcs_val;
//         fc2_minus_ext.push_back(dfc2_ext_tmp);
//     }
// }

// void Gruneisen::calc_gruneisen_old()
// {
// 
//     if (mympi->my_rank == 0) {
//         std::cout << " GRUNEISEN = 1 : Calculating Gruneisen parameters ... ";
//     }
// 
//     unsigned int nk = kpoint->nk;
//     unsigned int ns = dynamical->neval;
// 
//     unsigned int i;
//     unsigned int ik, is;
// 
//     double *eval_orig;
//     double **eval_plus;
//     double **eval_minus;
// 
//     double xk_tmp[3];
// 
//     std::complex<double> **evec_tmp;
// 
//     memory->allocate(evec_tmp, 1, 1); // dummy allocation
//     memory->allocate(eval_orig, ns);
//     memory->allocate(eval_plus, nk, ns);
//     memory->allocate(eval_minus, nk, ns);
// 
//     for (ik = 0; ik < nk; ++ik){
// 
//         for (i = 0; i < 3; ++i) xk_tmp[i] = kpoint->xk[ik][i];
// 
//         dynamical->eval_k(xk_tmp, kpoint->kvec_na[ik], fcs_phonon->fc2_ext, eval_orig, evec_tmp, false);
//         dynamical->eval_k(xk_tmp, kpoint->kvec_na[ik], fc2_plus_ext, eval_plus[ik], evec_tmp, false);
//         dynamical->eval_k(xk_tmp, kpoint->kvec_na[ik], fc2_minus_ext, eval_minus[ik], evec_tmp, false);
// 
//         for (is = 0; is < ns; ++is) {
//             gruneisen[ik][is] = (eval_plus[ik][is] - eval_minus[ik][is]) / (2.0 * delta_a) / (-6.0 * eval_orig[is]);
//         }
//     }
// 
//     memory->deallocate(evec_tmp);
//     memory->deallocate(eval_orig);
//     memory->deallocate(eval_plus);
//     memory->deallocate(eval_minus);
// 
//     std::cout << "done !" << std::endl;
// }

// void Gruneisen::prepare_delta_fc2()
// {
//     unsigned int i;
// 
//     double vec[3];
//     double fcs_tmp = 0.0;
// 
//     std::vector<FcsAlignedForGruneisen> fc3_aligned;
// 
//     std::vector<AtomCellSuper> pairs_vec;
//     std::vector<unsigned int> arr_old, arr_tmp;
// 
//     fc3_aligned.clear();
// 
//     for (std::vector<FcsArrayWithCell>::const_iterator it  = fcs_phonon->force_constant_with_cell[1].begin(); 
//                                                        it != fcs_phonon->force_constant_with_cell[1].end(); ++it) {
//         fc3_aligned.push_back(FcsAlignedForGruneisen((*it).fcs_val, (*it).pairs));
//     }
// 
//     std::sort(fc3_aligned.begin(), fc3_aligned.end());
// 
// 
//     delta_fc2.clear();
//     arr_old.clear();
// 
//     for (i = 0; i < 6; ++i) arr_old.push_back(-1);
// 
//     for (std::vector<FcsAlignedForGruneisen>::const_iterator it = fc3_aligned.begin(); it != fc3_aligned.end(); ++it) {
// 
//         arr_tmp.clear();
// 
//         arr_tmp.push_back((*it).pairs[0].index / 3);
//         arr_tmp.push_back((*it).pairs[0].index % 3);
// 
//         arr_tmp.push_back((*it).pairs[1].index / 3);
//         arr_tmp.push_back((*it).pairs[1].tran);
//         arr_tmp.push_back((*it).pairs[1].cell_s);
//         arr_tmp.push_back((*it).pairs[1].index % 3);
// 
//         if (arr_tmp != arr_old) {
// 
//             if (arr_old[0] != -1) { // Neglect the initial entry
//                 delta_fc2.push_back(FcsArrayWithCell(fcs_tmp, pairs_vec));
//             }
// 
//             pairs_vec.clear();
//             for (i = 0; i < 2; ++i) {
//                 pairs_vec.push_back((*it).pairs[i]);
//             }           
// 
//             fcs_tmp = 0.0;
//             arr_old.clear();
//             arr_old.reserve(arr_tmp.size());
//             std::copy(arr_tmp.begin(), arr_tmp.end(), std::back_inserter(arr_old));
//         }
// 
//         for (i = 0; i < 3; ++i) {
//             vec[i] = system->xr_s[system->map_p2s[(*it).pairs[2].index / 3][(*it).pairs[2].tran]][i] + xshift_s[(*it).pairs[2].cell_s][i];
//         }
//         
//         rotvec(vec, vec, system->lavec_s);
// 
//         fcs_tmp += (*it).fcs_val * vec[(*it).pairs[2].index % 3];
//     }
// 
//     delta_fc2.push_back(FcsArrayWithCell(fcs_tmp, pairs_vec));
// }


// void Gruneisen::write_newinfo_all()
// {
//     std::string file_info_plus, file_info_minus;
//     std::ofstream ofs_plus, ofs_minus;
//     std::ifstream ifs_orig;
// 
//     file_info_plus = input->job_title + "_+.info";
//     file_info_minus = input->job_title + "_-.info";
// 
//     if (print_newfcs) {
//         ofs_plus.open(file_info_plus.c_str(), std::ios::out);
//         if (!ofs_plus) error->exit("write_newinfo", "Cannot create the file file_info_plus");
//         ofs_minus.open(file_info_minus.c_str(), std::ios::out);
//         if (!ofs_minus) error->exit("write_newinfo", "Cannot create the file file_info_minus");
//         ifs_orig.open(fcs_phonon->file_fcs.c_str(), std::ios::in);
//         if (!ifs_orig) error->exit("write_newinfo", "Cannot open the info file");
// 
//         write_newinfo(ifs_orig, ofs_plus, delta_a, fc2_plus, fc2_plus_ext, fc3_plus);
//         ifs_orig.close();
//         ifs_orig.open(fcs_phonon->file_fcs.c_str(), std::ios::in);
//         write_newinfo(ifs_orig, ofs_minus, -delta_a, fc2_minus, fc2_minus_ext, fc3_minus);
// 
//         ofs_plus.close();
//         ofs_minus.close();
//         ifs_orig.close();
//     }
// 
// }

// void Gruneisen::write_newinfo(std::ifstream &ifs, std::ofstream &ofs, const double delta,
//                               double ****fc2_new, std::vector<FcsClassExtent> fc2_new_ext, std::vector<FcsClassGru> fc3_new)
// {
// 
//     int i, j, k;
//     std::string line;
// 
//     double factor = 1.0 + delta;
//     double aa[3][3];
// 
//     int nat = system->nat;
//     int natmin = system->natmin;
//     int nfc2;
// 
//     double dummy;
// 
//     int ind[3];
//     int *pairs, *pairs2;
//     double relvec[3];
// 
//     std::vector<int> elems;
// 
//     std::string str_pair[3];
//     int len[3], xyz[3], atom[3];
// 
//     std::vector<FcsClassGru>::iterator it_lower;
// 
//     memory->allocate(pairs, natmin);
//     memory->allocate(pairs2, natmin * 3);
// 
//     for (i = 0; i < 2; ++i) {
//         std::getline(ifs, line);
//         ofs << line << std::endl;
//     }
// 
//     for (i = 0; i < 3; ++i) {
//         for (j = 0; j < 3; ++j) {
//             ifs >> aa[i][j];
//             aa[i][j] *= factor;
//         }
//     }
// 
//     for (i = 0; i < 3; ++i) {
//         for (j = 0; j < 3; ++j) {
//             ofs << std::setw(25) << std::setprecision(16) << aa[i][j];
//         }
//         ofs << std::endl;
//     }
//     ifs.ignore();
// 
//     for (i = 0; i < nat + 7; ++i) {
//         std::getline(ifs, line);
//         ofs << line << std::endl;
//     }
// 
//     ifs >> nfc2;
//     ofs << nfc2 << std::endl;
// 
//     for (i = 0; i < nfc2; ++i) {
//         ifs >> dummy >> ind[0] >> ind[1];
// 
//         if (!fcs_phonon->is_fc2_ext) {
//             ofs << std::scientific << std::setprecision(16) << std::setw(25) << fc2_new[system->map_s2p[ind[0]/3].atom_num][ind[1]/3][ind[0]%3][ind[1]%3];
//         }
//         for (j = 0; j < 2; ++j) {
//             ofs << std::setw(7) << ind[j];
//         }
//         ofs << std::endl;
//     }
//     ifs.ignore();
// 
//     for (i = 0; i < 2; ++i) {
//         std::getline(ifs, line);
//         ofs << line << std::endl;
//     }
// 
//     while (std::getline(ifs,line)) {
// 
//         if (line == "##FORCE CONSTANTS") break;
// 
//         if (line == "#LIST_HARMONIC" || line == "#LIST_ANHARM3") {
// 
//             ofs << line << std::endl;
// 
//             for (i = 0; i < natmin; ++i) ifs >> pairs[i];
//             for (i = 0; i < natmin; ++i) ofs << std::setw(6) << pairs[i];
//             ofs << std::endl;
// 
//             for (i = 0; i < natmin; ++i) {
//                 for (j = 0; j < pairs[i]; ++j) {
//                     ifs >> ind[0] >> ind[1] >> relvec[0] >> relvec[1] >> relvec[2];
//                     ofs << std::setw(6) << ind[0] << std::setw(6) << ind[1];
//                     for (k = 0; k < 3; ++k) {
//                         relvec[k] *= factor;
//                         ofs << std::scientific << std::setprecision(16) << std::setw(25) << relvec[k];
//                     }
//                     ofs << std::endl;
//                     ;				}
//             }
//         }
//     }
// 
//     ofs << line << std::endl;
//     std::getline(ifs, line);
//     ofs << line << std::endl;
// 
//     while (std::getline(ifs, line)) {
// 
//         if (line == "#FCS_HARMONIC") {
//             ofs << line << std::endl;
//             std::getline(ifs, line);
//             ofs << line << std::endl;
// 
//             for (i = 0; i < 3 * natmin; ++i) ifs >> pairs2[i];
//             for (i = 0; i < 3 * natmin; ++i) ofs << std::setw(6) << pairs2[i];
//             ofs << std::endl;
// 
//             for (i = 0; i < 3 * natmin; ++i) {
//                 for (j = 0; j < pairs2[i]; ++j) {
//                     ifs >> dummy;
//                     ifs >> str_pair[0] >> str_pair[1];
// 
//                     for (k = 0; k < 2; ++k) len[k] = str_pair[k].length();
// 
//                     for (k = 0; k < 2; ++k) {
//                         if (str_pair[k][len[k] - 1] == 'x') {
//                             xyz[k] = 0;
//                         } else if (str_pair[k][len[k] - 1] == 'y') {
//                             xyz[k] = 1;
//                         } else if (str_pair[k][len[k] - 1] == 'z') {
//                             xyz[k] = 2;
//                         } else {
//                             error->exit("write_newinfo", "This cannot happen.");
//                         }
//                     }
// 
//                     for (k = 0; k < 2; ++k) {
//                         atom[k] = std::atoi(str_pair[k].substr(0,len[k]-1).c_str()) - 1;
//                     }
//                     atom[0] = system->map_s2p[atom[0]].atom_num;
//                     if (!fcs_phonon->is_fc2_ext) {
//                         ofs << std::scientific << std::setprecision(16) << std::setw(25) << fc2_new[atom[0]][atom[1]][xyz[0]][xyz[1]] << std::endl;
//                     } else {
//                         ofs << 0.0 << std::endl;
//                     }
//                     ofs << std::setw(5) << str_pair[0] << std::setw(5) << str_pair[1] << std::endl;
//                 }
//             }
//         }
// 
//         if (line == "#FCS_ANHARM3") {
//             ofs << line << std::endl;
//             std::getline(ifs, line);
//             ofs << line << std::endl;
// 
//             for (i = 0; i < 3 * natmin; ++i) ifs >> pairs2[i];
//             for (i = 0; i < 3 * natmin; ++i) ofs << std::setw(6) << pairs2[i];
//             ofs << std::endl;
// 
//             for (i = 0; i < 3 * natmin; ++i) {
//                 for (j = 0; j < pairs2[i]; ++j) {
//                     ifs >> dummy;
//                     ifs >> str_pair[0] >> str_pair[1] >> str_pair[2];
// 
//                     for (k = 0; k < 3; ++k) len[k] = str_pair[k].length();
// 
//                     for (k = 0; k < 3; ++k) {
//                         if (str_pair[k][len[k] - 1] == 'x') {
//                             xyz[k] = 0;
//                         } else if (str_pair[k][len[k] - 1] == 'y') {
//                             xyz[k] = 1;
//                         } else if (str_pair[k][len[k] - 1] == 'z') {
//                             xyz[k] = 2;
//                         } else {
//                             error->exit("write_newinfo", "This cannot happen.");
//                         }
//                     }
// 
//                     elems.clear();
//                     for (k = 0; k < 3; ++k) {
//                         atom[k] = std::atoi(str_pair[k].substr(0,len[k]-1).c_str()) - 1;
//                         ind[k] = 3 * atom[k] + xyz[k];
//                         elems.push_back(3 * atom[k] + xyz[k]);
//                     }
// 
//                     it_lower = lower_bound(fc3_new.begin(), fc3_new.end(), FcsClassGru(3, ind, dummy));
// 
//                     if (it_lower == fc3_new.end()) {
//                         error->exit("write_newinfo", "FC3 with the same index could not be found.");
//                     } else {
//                         ofs << std::scientific << std::setprecision(16) << std::setw(25) << (*it_lower).fcs_val << std::endl;
//                         ofs << std::setw(5) << str_pair[0] << std::setw(5) << str_pair[1] << std::setw(5) << str_pair[2] << std::endl;
//                     }
//                 }
//             }
// 
//         }
// 
//         if (line == "#FCS_HARMONIC_EXT") {
//             ofs << line << std::endl;
// 
// 
//             for (i = 0; i < 3 * natmin; ++i) pairs2[i] = 0;
//             for (std::vector<FcsClassExtent>::const_iterator it = fc2_new_ext.begin(); it != fc2_new_ext.end(); ++it) {
//                 pairs2[3 * (*it).atm1 + (*it).xyz1] += 1;
//             }
// 
//             nfc2 = 0;
//             for (i = 0; i < 3 * natmin; ++i) nfc2 += pairs2[i];
// 
//             ofs << std::setw(10) << nfc2 << std::endl;
//             for (i = 0; i < 3 * natmin; ++i) ofs << std::setw(6) << pairs2[i];
//             ofs << std::endl;
// 
//             for (std::vector<FcsClassExtent>::const_iterator it = fc2_new_ext.begin(); it != fc2_new_ext.end(); ++it) {
//                 ofs << std::setw(5) << (*it).atm1 << std::setw(5) << (*it).xyz1;
//                 ofs << std::setw(8) << (*it).atm2 << std::setw(5) << (*it).xyz2;
//                 ofs << std::setw(5) << (*it).cell_s;
//                 ofs << std::scientific << std::setprecision(16) << std::setw(25) << (*it).fcs_val << std::endl;
//             }
//         }
//     }
// }

// 
// void Gruneisen::prepare_newfc3()
// {
//     int i;
//     int ind[4], arr_old[4], arr[4];
//     int natmin = system->natmin;
//     int nat = system->nat;
//     unsigned int nalpha[4], ncell[4], ixyz[4];
//     unsigned int atom_s[4];
// 
//     double coord_tmp[3];
//     double fcs_tmp;
// 
//     std::vector<unsigned int> arr3, arr4;
//     std::vector<FcsClassGru> fc3_copy, fc4_copy, dfc3;
//     std::vector<FcsClassGru>::iterator it_lower;
// 
//     fc4_copy.clear();
//     dfc3.clear();
// 
//     for (std::vector<FcsClass>::const_iterator it_fc4 = fcs_phonon->force_constant[2].begin(); it_fc4 != fcs_phonon->force_constant[2].end(); ++it_fc4) {
//         FcsClass fc4_tmp = *it_fc4;
// 
//         for (i = 0; i < 4; ++i) {
//             nalpha[i] = fc4_tmp.elems[i].atom;
//             ncell[i] = fc4_tmp.elems[i].cell;
//             ixyz[i] = fc4_tmp.elems[i].xyz;
//             atom_s[i] = system->map_p2s[nalpha[i]][ncell[i]];
//             ind[i] = 3 * atom_s[i] + ixyz[i];
//         }
// 
//         fc4_copy.push_back(FcsClassGru(4, ind, fc4_tmp.fcs_val));
//     }
// 
//     std::sort(fc4_copy.begin(), fc4_copy.end());
// 
// 
//     for (i = 0; i < 3; ++i) arr_old[i] = -1;
//     fcs_tmp = 0.0;
// 
//     for (std::vector<FcsClassGru>::const_iterator it_fc4 = fc4_copy.begin(); it_fc4 != fc4_copy.end(); ++it_fc4) {
//         FcsClassGru fc4_tmp = *it_fc4;
// 
//         for (i = 0; i < 4; ++i) arr[i] = fc4_tmp.elems[i];
// 
//         if (arr[0] != arr_old[0] || arr[1] != arr_old[1] || arr[2] != arr_old[2]) {
//             if (std::abs(fcs_tmp) > eps12) {
//                 dfc3.push_back(FcsClassGru(3, arr_old, fcs_tmp));
//             }
//             fcs_tmp = 0.0;
// 
//             for (i = 0; i < 3; ++i) arr_old[i] = arr[i];
//         }
// 
//         for (i = 0; i < 3; ++i) {
//             coord_tmp[i] = system->xr_s[arr[3]/3][i] - system->xr_s[arr[0]/3][i];
//             coord_tmp[i] = dynamical->fold(coord_tmp[i]);
//         }
//         rotvec(coord_tmp, coord_tmp, system->lavec_s);
// 
//         fcs_tmp += fc4_tmp.fcs_val * coord_tmp[arr[3] % 3];
//     }
//     if (std::abs(fcs_tmp) > eps12) 	dfc3.push_back(FcsClassGru(3, arr_old, fcs_tmp));
// 
// 
//     if (fcs_phonon->force_constant[1].size() < dfc3.size()) {
//         error->exit("prepare_newfc3", "The number of DFC3 is greater than the number of FC3.\n \
//                                       The cutoff radius for FC4 should not be greater than that of FC3.");
//     }
//     fc4_copy.clear();
//     fc3_copy.clear();
// 
// 
//     for (std::vector<FcsClass>::const_iterator it_fc3 = fcs_phonon->force_constant[1].begin(); it_fc3 != fcs_phonon->force_constant[1].end(); ++it_fc3) {
//         FcsClass fc3_tmp = *it_fc3;
// 
//         for (i = 0; i < 3; ++i) {
//             nalpha[i] = fc3_tmp.elems[i].atom;
//             ncell[i] = fc3_tmp.elems[i].cell;
//             ixyz[i] = fc3_tmp.elems[i].xyz;
//             atom_s[i] = system->map_p2s[nalpha[i]][ncell[i]];
//             ind[i] = 3 * atom_s[i] + ixyz[i];
//         }
//         fc3_copy.push_back(FcsClassGru(3, ind, fc3_tmp.fcs_val));
//     }
//     std::sort(fc3_copy.begin(), fc3_copy.end());
// 
//     fc3_plus.clear();
//     fc3_minus.clear();
// 
//     for (std::vector<FcsClassGru>::const_iterator it = fc3_copy.begin(); it != fc3_copy.end(); ++it) {
//         fc3_plus.push_back(*it);
//         fc3_minus.push_back(*it);
//     }
// 
// 
//     for (std::vector<FcsClassGru>::const_iterator it = dfc3.begin(); it != dfc3.end(); ++it) {
// 
//         for (i = 0; i < 3; ++i) ind[i] = (*it).elems[i];
// 
//         it_lower = lower_bound(fc3_plus.begin(), fc3_plus.end(), FcsClassGru(3, ind, fcs_tmp));
//         if (it_lower == fc3_plus.end()) {
//             error->exit("prepare_newfc3", "The list of FC3 doesn't contain a force constant in DFC3");
//         } else {
//             (*it_lower).fcs_val += delta_a * (*it).fcs_val;
//         }
// 
//         it_lower = lower_bound(fc3_minus.begin(), fc3_minus.end(), FcsClassGru(3, ind, fcs_tmp));
//         if (it_lower == fc3_minus.end()) {
//             error->exit("prepare_newfc3", "The list of FC3 doesn't contain a force constant in DFC3");
//         } else {
//             (*it_lower).fcs_val -= delta_a * (*it).fcs_val;
//         }
//     }
// }
// 

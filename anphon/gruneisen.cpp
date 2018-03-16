/*
gruneisen.cpp

Copyright (c) 2014, 2015, 2016 Terumasa Tadano

This file is distributed under the terms of the MIT license.
Please see the file 'LICENCE.txt' in the root directory 
or http://opensource.org/licenses/mit-license.php for information.
*/

#include "mpi_common.h"
#include "gruneisen.h"
#include "constants.h"
#include "dynamical.h"
#include "error.h"
#include "fcs_phonon.h"
#include "kpoint.h"
#include "mathfunctions.h"
#include "memory.h"
#include "parsephon.h"
#include "pointers.h"
#include "system.h"
#include "anharmonic_core.h"
#include "version.h"
#include <iostream>
#include <iomanip>
#include <boost/lexical_cast.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/version.hpp>
#include <cmath>

using namespace PHON_NS;

Gruneisen::Gruneisen(PHON *phon): Pointers(phon)
{
    set_default_variables();
};

Gruneisen::~Gruneisen()
{
    deallocate_variables();
};

void Gruneisen::set_default_variables()
{
    delta_a = 0.01;
    print_gruneisen = false;
    print_newfcs = false;
    gruneisen = nullptr;
    xshift_s = nullptr;
}

void Gruneisen::deallocate_variables()
{
    if (gruneisen) {
        memory->deallocate(gruneisen);
    }
    if (xshift_s) {
        memory->deallocate(xshift_s);
    }
    delta_fc2.clear();
    delta_fc3.clear();
}


void Gruneisen::setup()
{
    MPI_Bcast(&delta_a, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&print_newfcs, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD);

    memory->allocate(xshift_s, 27, 3);

    for (int i = 0; i < 3; ++i) xshift_s[0][i] = 0.0;

    int icell = 0;

    for (int ix = -1; ix <= 1; ++ix) {
        for (int iy = -1; iy <= 1; ++iy) {
            for (int iz = -1; iz <= 1; ++iz) {
                if (ix == 0 && iy == 0 && iz == 0) continue;

                ++icell;

                xshift_s[icell][0] = static_cast<double>(ix);
                xshift_s[icell][1] = static_cast<double>(iy);
                xshift_s[icell][2] = static_cast<double>(iz);
            }
        }
    }

    if (print_gruneisen || print_newfcs) {
        prepare_delta_fcs(fcs_phonon->force_constant_with_cell[1], delta_fc2);
    }

    if (print_newfcs && anharmonic_core->quartic_mode > 0) {
        prepare_delta_fcs(fcs_phonon->force_constant_with_cell[2], delta_fc3);
    }
    if (print_gruneisen) {
        memory->allocate(gruneisen, kpoint->nk, dynamical->neval);
    }

    if (mympi->my_rank == 0) {
        if (print_newfcs) {
            std::cout << std::endl;
            if (anharmonic_core->quartic_mode > 0) {
                std::cout << " NEWFCS = 1 : Harmonic and cubic force constants of " << std::endl;
            } else {
                std::cout << " NEWFCS = 1 : Harmonic force constants of " << std::endl;
            }
            std::cout << "              expanded/compressed systems will be estimated" << std::endl;
            std::cout << "              with DELTA_A = " << std::setw(5) << delta_a << std::endl;
        }
    }
    //   print_stress_energy();
}


void Gruneisen::calc_gruneisen()
{
    unsigned int i, j;
    auto ns = dynamical->neval;
    auto nk = kpoint->nk;
    double gamma_imag;
    std::complex<double> **dfc2_reciprocal;

    memory->allocate(dfc2_reciprocal, ns, ns);

    if (mympi->my_rank == 0) {
        std::cout << std::endl;
        std::cout << " GRUNEISEN = 1 : Calculating Gruneisen parameters ... ";
    }

    for (auto ik = 0; ik < nk; ++ik) {

        calc_dfc2_reciprocal(dfc2_reciprocal, kpoint->xk[ik]);

        for (auto is = 0; is < ns; ++is) {

            gruneisen[ik][is] = std::complex<double>(0.0, 0.0);

            for (i = 0; i < ns; ++i) {
                for (j = 0; j < ns; ++j) {
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
        std::cout << "done!" << std::endl;
    }
}

void Gruneisen::calc_dfc2_reciprocal(std::complex<double> **dphi2,
                                     const double *xk_in)
{
    unsigned int i, j;
    unsigned int ns = dynamical->neval;

    unsigned int atm1, atm2, xyz1, xyz2;
    unsigned int atm1_s, atm2_s;
    unsigned int tran, cell_s;

    double vec[3];
    double phase;

    std::complex<double> im(0.0, 1.0);


    for (i = 0; i < ns; ++i) {
        for (j = 0; j < ns; ++j) {
            dphi2[i][j] = std::complex<double>(0.0, 0.0);
        }
    }

    for (const auto &it : delta_fc2) {

        atm1 = it.pairs[0].index / 3;
        xyz1 = it.pairs[0].index % 3;
        atm2 = it.pairs[1].index / 3;
        xyz2 = it.pairs[1].index % 3;

        tran = it.pairs[1].tran;
        cell_s = it.pairs[1].cell_s;

        atm1_s = system->map_p2s_anharm[atm1][0];
        atm2_s = system->map_p2s_anharm[atm2][tran];


        for (i = 0; i < 3; ++i) {
            vec[i] = system->xr_s_anharm[atm2_s][i] + xshift_s[cell_s][i]
                - system->xr_s_anharm[system->map_p2s_anharm[atm2][0]][i];
        }

        rotvec(vec, vec, system->lavec_s_anharm);
        rotvec(vec, vec, system->rlavec_p);

        phase = vec[0] * xk_in[0] + vec[1] * xk_in[1] + vec[2] * xk_in[2];

        dphi2[3 * atm1 + xyz1][3 * atm2 + xyz2]
            += it.fcs_val * std::exp(im * phase)
            / std::sqrt(system->mass_anharm[atm1_s] * system->mass_anharm[atm2_s]);

    }
}


void Gruneisen::prepare_delta_fcs(const std::vector<FcsArrayWithCell> &fcs_in,
                                  std::vector<FcsArrayWithCell> &delta_fcs)
{
    unsigned int i;
    double vec[3];
    double fcs_tmp = 0.0;

    std::vector<FcsAlignedForGruneisen> fcs_aligned;
    std::vector<AtomCellSuper> pairs_vec;
    std::vector<int> index_old, index_now;
    std::vector<int> index_with_cell;
    std::set<std::vector<int>> set_index_uniq;
    AtomCellSuper pairs_tmp;

    unsigned int norder = fcs_in[0].pairs.size();
    unsigned int nelems;
    unsigned int nmulti;

    delta_fcs.clear();
    fcs_aligned.clear();

    for (auto it = fcs_in.cbegin(); it != fcs_in.cend(); ++it) {
        fcs_aligned.emplace_back((*it).fcs_val, (*it).pairs);
    }
    std::sort(fcs_aligned.begin(), fcs_aligned.end());

    index_old.clear();
    nelems = 2 * (norder - 2) + 1;
    for (i = 0; i < nelems; ++i) index_old.push_back(-1);

    index_with_cell.clear();
    set_index_uniq.clear();

    for (auto it = fcs_aligned.cbegin(); it != fcs_aligned.cend(); ++it) {

        index_now.clear();
        index_with_cell.clear();

        index_now.push_back((*it).pairs[0].index);
        index_with_cell.push_back((*it).pairs[0].index);

        for (i = 1; i < norder - 1; ++i) {
            index_now.push_back((*it).pairs[i].index);
            index_now.push_back((*it).pairs[i].tran);

            index_with_cell.push_back((*it).pairs[i].index);
            index_with_cell.push_back((*it).pairs[i].tran);
            index_with_cell.push_back((*it).pairs[i].cell_s);
        }

        if (index_now != index_old) {

            if (index_old[0] != -1) {

                nmulti = set_index_uniq.size();
                fcs_tmp /= static_cast<double>(nmulti);

                if (std::abs(fcs_tmp) > eps15) {
                    for (const auto &it2 : set_index_uniq) {

                        pairs_vec.clear();

                        pairs_tmp.index = it2[0];
                        pairs_tmp.tran = 0;
                        pairs_tmp.cell_s = 0;
                        pairs_vec.push_back(pairs_tmp);
                        for (i = 1; i < norder - 1; ++i) {
                            pairs_tmp.index = it2[3 * i - 2];
                            pairs_tmp.tran = it2[3 * i - 1];
                            pairs_tmp.cell_s = it2[3 * i];
                            pairs_vec.push_back(pairs_tmp);
                        }
                        delta_fcs.emplace_back(fcs_tmp, pairs_vec);
                    }
                }
                set_index_uniq.clear();
            }

            fcs_tmp = 0.0;
            index_old.clear();
            index_old.reserve(index_now.size());
            std::copy(index_now.begin(), index_now.end(), std::back_inserter(index_old));
        }

        set_index_uniq.insert(index_with_cell);

        for (i = 0; i < 3; ++i) {
            vec[i] = system->xr_s_anharm[system->map_p2s_anharm[(*it).pairs[norder - 1].index / 3][(*it).pairs[norder -
                    1].tran]]
                [i]
                - system->xr_s_anharm[system->map_p2s_anharm[(*it).pairs[0].index / 3][0]][i]
                + xshift_s[(*it).pairs[norder - 1].cell_s][i];
        }

        rotvec(vec, vec, system->lavec_s_anharm);

        fcs_tmp += (*it).fcs_val * vec[(*it).pairs[norder - 1].index % 3];
    }

    nmulti = set_index_uniq.size();
    fcs_tmp /= static_cast<double>(nmulti);

    if (std::abs(fcs_tmp) > eps15) {
        for (const auto &it2 : set_index_uniq) {

            pairs_vec.clear();

            pairs_tmp.index = it2[0];
            pairs_tmp.tran = 0;
            pairs_tmp.cell_s = 0;
            pairs_vec.push_back(pairs_tmp);
            for (i = 1; i < norder - 1; ++i) {
                pairs_tmp.index = it2[3 * i - 2];
                pairs_tmp.tran = it2[3 * i - 1];
                pairs_tmp.cell_s = it2[3 * i];
                pairs_vec.push_back(pairs_tmp);
            }
            delta_fcs.emplace_back(fcs_tmp, pairs_vec);
        }
    }

    fcs_aligned.clear();
    set_index_uniq.clear();
}

void Gruneisen::write_new_fcsxml_all()
{
    std::cout << std::endl;

    if (fcs_phonon->update_fc2) {
        error->warn("write_new_fcsxml_all",
                    "NEWFCS = 1 cannot be combined with the FC2XML.");
    } else {
        std::cout << " NEWFCS = 1 : Following XML files are created. " << std::endl;

        std::string file_xml = input->job_title + "_+.xml";
        write_new_fcsxml(file_xml, delta_a);

        std::cout << "  " << std::setw(input->job_title.length() + 12) << std::left << file_xml;
        std::cout << " : Force constants of the system expanded by "
            << std::fixed << std::setprecision(3) << delta_a * 100 << " %" << std::endl;

        file_xml = input->job_title + "_-.xml";
        write_new_fcsxml(file_xml, -delta_a);

        std::cout << "  " << std::setw(input->job_title.length() + 12) << std::left << file_xml;
        std::cout << " : Force constants of the system compressed by "
            << std::fixed << std::setprecision(3) << delta_a * 100 << " %" << std::endl;
    }
}

void Gruneisen::write_new_fcsxml(const std::string filename_xml,
                                 const double change_ratio_of_a)
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

    pt.put("Data.ANPHON_version", ALAMODE_VERSION);
    pt.put("Data.Description.OriginalXML", fcs_phonon->file_fcs);
    pt.put("Data.Description.Delta_A", double2string(change_ratio_of_a));

    pt.put("Data.Structure.NumberOfAtoms", system->nat);
    pt.put("Data.Structure.NumberOfElements", system->nkd);

    for (i = 0; i < system->nkd; ++i) {
        ptree &child = pt.add("Data.Structure.AtomicElements.element",
                              system->symbol_kd[i]);
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
            ptree &child = pt.add("Data.Symmetry.Translations.map",
                                  system->map_p2s[j][i] + 1);
            child.put("<xmlattr>.tran", i + 1);
            child.put("<xmlattr>.atom", j + 1);
        }
    }

    pt.put("Data.ForceConstants", "");
    str_tmp.clear();

    for (const auto &it : fcs_phonon->force_constant_with_cell[0]) {

        ptree &child = pt.add("Data.ForceConstants.HARMONIC.FC2", double2string(it.fcs_val));

        child.put("<xmlattr>.pair1",
                  std::to_string(it.pairs[0].index / 3 + 1)
                  + " " + std::to_string(it.pairs[0].index % 3 + 1));
        child.put("<xmlattr>.pair2",
                  std::to_string(system->map_p2s[it.pairs[1].index / 3][it.pairs[1].tran] + 1)
                  + " " + std::to_string(it.pairs[1].index % 3 + 1)
                  + " " + std::to_string(it.pairs[1].cell_s + 1));
    }

    for (const auto &it : delta_fc2) {

        if (std::abs(it.fcs_val) < eps12) continue;

        ptree &child = pt.add("Data.ForceConstants.HARMONIC.FC2",
                              double2string(change_ratio_of_a * it.fcs_val));

        child.put("<xmlattr>.pair1",
                  std::to_string(it.pairs[0].index / 3 + 1)
                  + " " + std::to_string(it.pairs[0].index % 3 + 1));
        child.put("<xmlattr>.pair2",
                  std::to_string(system->map_p2s[it.pairs[1].index / 3][it.pairs[1].tran] + 1)
                  + " " + std::to_string(it.pairs[1].index % 3 + 1)
                  + " " + std::to_string(it.pairs[1].cell_s + 1));
    }

    if (anharmonic_core->quartic_mode) {
        for (const auto &it : fcs_phonon->force_constant_with_cell[1]) {

            if (it.pairs[1].index > it.pairs[2].index) continue;

            ptree &child = pt.add("Data.ForceConstants.ANHARM3.FC3",
                                  double2string(it.fcs_val));

            child.put("<xmlattr>.pair1",
                      std::to_string(it.pairs[0].index / 3 + 1)
                      + " " + std::to_string(it.pairs[0].index % 3 + 1));
            child.put("<xmlattr>.pair2",
                      std::to_string(system->map_p2s[it.pairs[1].index / 3][it.pairs[1].tran] + 1)
                      + " " + std::to_string(it.pairs[1].index % 3 + 1)
                      + " " + std::to_string(it.pairs[1].cell_s + 1));
            child.put("<xmlattr>.pair3",
                      std::to_string(system->map_p2s[it.pairs[2].index / 3][it.pairs[2].tran] + 1)
                      + " " + std::to_string(it.pairs[2].index % 3 + 1)
                      + " " + std::to_string(it.pairs[2].cell_s + 1));
        }

        for (const auto &it : delta_fc3) {

            if (std::abs(it.fcs_val) < eps12) continue;

            if (it.pairs[1].index > it.pairs[2].index) continue;

            ptree &child = pt.add("Data.ForceConstants.ANHARM3.FC3",
                                  double2string(change_ratio_of_a * it.fcs_val));

            child.put("<xmlattr>.pair1",
                      std::to_string(it.pairs[0].index / 3 + 1)
                      + " " + std::to_string(it.pairs[0].index % 3 + 1));
            child.put("<xmlattr>.pair2",
                      std::to_string(system->map_p2s[it.pairs[1].index / 3][it.pairs[1].tran] + 1)
                      + " " + std::to_string(it.pairs[1].index % 3 + 1)
                      + " " + std::to_string(it.pairs[1].cell_s + 1));
            child.put("<xmlattr>.pair3",
                      std::to_string(system->map_p2s[it.pairs[2].index / 3][it.pairs[2].tran] + 1)
                      + " " + std::to_string(it.pairs[2].index % 3 + 1)
                      + " " + std::to_string(it.pairs[2].cell_s + 1));
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


std::string Gruneisen::double2string(const double d)
{
    std::string rt;
    std::stringstream ss;

    ss << std::scientific << std::setprecision(15) << d;
    ss >> rt;
    return rt;
}


// double Gruneisen::calc_stress_energy2(const std::vector<FcsArrayWithCell> fcs_in)
// {
//     unsigned int i, j;
//     double ret = 0.0;
//     double **vec, **pos;
//     double tmp, tmp2;
//     double xshift[3];
//     unsigned int itran;
//     unsigned int norder = fcs_in[0].pairs.size();
// 
//     memory->allocate(vec, norder, 3);
//     memory->allocate(pos, norder, 3);
// 
//     for (std::vector<FcsArrayWithCell>::const_iterator it = fcs_in.begin(); it != fcs_in.end(); ++it) {
// 
//         for (i = 0; i < norder; ++i) {
//             for (j = 0; j < 3; ++j) {
//                 vec[i][j] = system->xr_s[system->map_p2s[(*it).pairs[i].index / 3][(*it).pairs[i].tran]][j]
//                 + xshift_s[(*it).pairs[i].cell_s][j];
// 
//                 pos[i][j] = system->xr_s[system->map_p2s[(*it).pairs[i].index / 3][0]][j];
//             //    vec[i][j] = system->xr_s[system->map_p2s[0][(*it).pairs[i].tran]][j] + xshift_s[(*it).pairs[i].cell_s][j];
//             }
//             rotvec(vec[i], vec[i], system->lavec_s);
//             rotvec(pos[i], pos[i], system->lavec_s);
//         } 
// 
//         
//         ret += (*it).fcs_val 
//             * (vec[1][(*it).pairs[0].index % 3] - pos[0][(*it).pairs[0].index % 3])
//             * (vec[1][(*it).pairs[1].index % 3] - pos[0][(*it).pairs[1].index % 3]);
//     }
// 
//     memory->deallocate(vec);
//     memory->deallocate(pos);
//     return ret;
// }
// 
// void Gruneisen::calc_stress_energy3(const std::vector<FcsArrayWithCell> fcs_in, double ****ret)
// {
//     unsigned int i, j, k, l;
//     double **vec, **pos;
//     double tmp, tmp2;
//     double xshift[3];
//     unsigned int itran;
//     unsigned int norder = fcs_in[0].pairs.size();
//     unsigned int crd[4];
// 
//     memory->allocate(vec, norder, 3);
//     memory->allocate(pos, norder, 3);
// 
//     for (i = 0; i < 3; ++i) {
//         for (j = 0; j < 3; ++j) {
//             for (k = 0; k < 3; ++k) {
//                 for (l = 0; l < 3; ++l) {
//                     ret[i][j][k][l] = 0.0;
//                 }
//             }
//         }
//     }
// 
//     for (std::vector<FcsArrayWithCell>::const_iterator it = fcs_in.begin(); it != fcs_in.end(); ++it) {
// 
//         for (i = 0; i < norder; ++i) {
//             for (j = 0; j < 3; ++j) {
//                 vec[i][j] = system->xr_s[system->map_p2s[(*it).pairs[i].index / 3][(*it).pairs[i].tran]][j]
//                 + xshift_s[(*it).pairs[i].cell_s][j];
// 
//                 pos[i][j] = system->xr_s[system->map_p2s[(*it).pairs[i].index / 3][0]][j];
//             }
//             rotvec(vec[i], vec[i], system->lavec_s);
//             rotvec(pos[i], pos[i], system->lavec_s);
//         }
// 
//         crd[0] = (*it).pairs[0].index % 3;
//         crd[1] = (*it).pairs[1].index % 3;
// 
//         for (k = 0; k < 3; ++k) {
//             
//             crd[2] = k;
// 
//             for (l = 0; l < 3; ++l) {
// 
//                 crd[3] = l;
// 
//                 ret[crd[0]][crd[1]][k][l] += (*it).fcs_val * (vec[1][k] - pos[0][k]) * (vec[1][l] - pos[0][l]);
//             }
//         }
//     }
// 
//     memory->deallocate(vec);
//     memory->deallocate(pos);
// 
//     for (i = 0; i < 3; ++i) {
//         for (j = 0; j < 3; ++j) {
//             for (k = 0; k < 3; ++k) {
//                 for (l = 0; l < 3; ++l) {
//                     ret[i][j][k][l] *= -0.5;
//                 }
//             }
//         }
//     }
// }
// 
// 
// void Gruneisen::print_stress_energy()
// {
// 
//     double volume = system->volume_p * std::pow(Bohr_in_Angstrom, 3) * 1.0e-30;
// 
// 
//     double ****A, ****C;
// 
//     memory->allocate(A, 3, 3, 3, 3);
//     memory->allocate(C, 3, 3, 3, 3);
// 
//     calc_stress_energy3(fcs_phonon->force_constant_with_cell[0], A);
// 
//     unsigned int i, j, k, l;
// 
//     std::cout << "# A [Ryd]" << std::endl;
// 
//     for (i = 0; i < 3; ++i) {
//         for (j = 0; j < 3; ++j) {
//             for (k = 0; k < 3; ++k) {
//                 for (l = 0; l < 3; ++l) {
//                     std::cout << std::setw(3) << i + 1;
//                     std::cout << std::setw(3) << j + 1;
//                     std::cout << std::setw(3) << k + 1;
//                     std::cout << std::setw(3) << l + 1;
//                     std::cout << std::setw(15) << std::fixed << A[i][j][k][l];
//                     std::cout << std::endl;
//                 }
//             }
//         }
//     }
// 
//     std::cout << std::endl;
//     std::cout << "# C [GPa]" << std::endl;
// 
//     for (i = 0; i < 3; ++i) {
//         for (j = 0; j < 3; ++j) {
//             for (k = 0; k < 3; ++k) {
//                 for (l = 0; l < 3; ++l) {
//                     C[i][j][k][l] = A[i][k][j][l] + A[j][k][i][l] - A[i][j][k][l];
//                     C[i][j][k][l] *= 1.0e-9 * Ryd / volume;
//                     std::cout << std::setw(3) << i + 1;
//                     std::cout << std::setw(3) << j + 1;
//                     std::cout << std::setw(3) << k + 1;
//                     std::cout << std::setw(3) << l + 1;
//                     std::cout << std::setw(15) << std::fixed << C[i][j][k][l];
//                     std::cout << std::endl;
// 
//                 }
//             }
//         }
//     }
// 
//     std::cout << "Bulk Modulus [GPa] = " << (C[0][0][0][0] + 2.0 * C[0][0][1][1]) / 3.0 << std::endl;
// 
//     memory->deallocate(A);
//     memory->deallocate(C);
// }

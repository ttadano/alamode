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
#include "phonon_dos.h"
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

Gruneisen::Gruneisen(PHON *phon) : Pointers(phon)
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
    gruneisen_bs = nullptr;
    gruneisen_dos = nullptr;
    xshift_s = nullptr;
}

void Gruneisen::deallocate_variables()
{
    if (gruneisen_bs) {
        deallocate(gruneisen_bs);
    }
    if (gruneisen_dos) {
        deallocate(gruneisen_dos);
    }
    if (xshift_s) {
        deallocate(xshift_s);
    }
    delta_fc2.clear();
    delta_fc3.clear();
}

void Gruneisen::setup()
{
    MPI_Bcast(&delta_a, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&print_newfcs, 1, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);

    allocate(xshift_s, 27, 3);

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
        prepare_delta_fcs(fcs_phonon->force_constant_with_cell[1], delta_fc2, 1);

        // impose_ASR_on_harmonic_IFC(delta_fc2, 0);
    }

    if (print_newfcs && anharmonic_core->quartic_mode > 0) {
        prepare_delta_fcs(fcs_phonon->force_constant_with_cell[2], delta_fc3, 1);
    }
    if (print_gruneisen) {
        if (kpoint->kpoint_bs) {
            allocate(gruneisen_bs, kpoint->kpoint_bs->nk, dynamical->neval);
        }
        if (dos->kmesh_dos) {
            allocate(gruneisen_dos, dos->kmesh_dos->nk, dynamical->neval);
        }
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
    const auto ns = dynamical->neval;
    std::complex<double> **dfc2_reciprocal;

    allocate(dfc2_reciprocal, ns, ns);

    if (mympi->my_rank == 0) {
        std::cout << std::endl;
        std::cout << " GRUNEISEN = 1 : Calculating Gruneisen parameters ... ";
    }

    if (kpoint->kpoint_bs) {
        const auto nk = kpoint->kpoint_bs->nk;
        const auto xk = kpoint->kpoint_bs->xk;
        const auto eval = dynamical->dymat_band->get_eigenvalues();
        const auto evec = dynamical->dymat_band->get_eigenvectors();

        for (auto ik = 0; ik < nk; ++ik) {

            calc_dfc2_reciprocal(dfc2_reciprocal, xk[ik]);

            for (auto is = 0; is < ns; ++is) {

                gruneisen_bs[ik][is] = std::complex<double>(0.0, 0.0);

                for (unsigned int i = 0; i < ns; ++i) {
                    for (unsigned int j = 0; j < ns; ++j) {
                        gruneisen_bs[ik][is] += std::conj(evec[ik][is][i])
                                                * dfc2_reciprocal[i][j]
                                                * evec[ik][is][j];
                    }
                }

                const auto gamma_imag = gruneisen_bs[ik][is].imag();
                if (std::abs(gamma_imag) > eps10) {
                    warn("calc_gruneisen", "Gruneisen parameter is not real");
                }

                if (std::abs(eval[ik][is]) < eps8) {
                    gruneisen_bs[ik][is] = 0.0;
                } else {
                    gruneisen_bs[ik][is] /= -6.0 * std::pow(eval[ik][is], 2);
                }
            }
        }
    }

    if (dos->kmesh_dos) {
        const auto nk = dos->kmesh_dos->nk;
        const auto xk = dos->kmesh_dos->xk;
        const auto eval = dos->dymat_dos->get_eigenvalues();
        const auto evec = dos->dymat_dos->get_eigenvectors();

        for (auto ik = 0; ik < nk; ++ik) {

            calc_dfc2_reciprocal(dfc2_reciprocal, xk[ik]);

            for (auto is = 0; is < ns; ++is) {

                gruneisen_dos[ik][is] = std::complex<double>(0.0, 0.0);

                for (unsigned int i = 0; i < ns; ++i) {
                    for (unsigned int j = 0; j < ns; ++j) {
                        gruneisen_dos[ik][is] += std::conj(evec[ik][is][i])
                                                 * dfc2_reciprocal[i][j]
                                                 * evec[ik][is][j];
                    }
                }

                const auto gamma_imag = gruneisen_dos[ik][is].imag();
                if (std::abs(gamma_imag) > eps10) {
                    warn("calc_gruneisen", "Gruneisen parameter is not real");
                }

                if (std::abs(eval[ik][is]) < eps8) {
                    gruneisen_dos[ik][is] = 0.0;
                } else {
                    gruneisen_dos[ik][is] /= -6.0 * std::pow(eval[ik][is], 2);
                }
            }
        }
    }

    deallocate(dfc2_reciprocal);

    if (mympi->my_rank == 0) {
        std::cout << "done!" << std::endl;
    }
}

void Gruneisen::calc_dfc2_reciprocal(std::complex<double> **dphi2,
                                     const double *xk_in)
{
    unsigned int i;
    const auto ns = dynamical->neval;

    double vec[3];

    const std::complex<double> im(0.0, 1.0);

    for (i = 0; i < ns; ++i) {
        for (unsigned int j = 0; j < ns; ++j) {
            dphi2[i][j] = std::complex<double>(0.0, 0.0);
        }
    }

    for (const auto &it: delta_fc2) {

        const auto atm1 = it.pairs[0].index / 3;
        const auto xyz1 = it.pairs[0].index % 3;
        const auto atm2 = it.pairs[1].index / 3;
        const auto xyz2 = it.pairs[1].index % 3;

        const auto tran = it.pairs[1].tran;
        const auto cell_s = it.pairs[1].cell_s;

        const auto atm1_s = system->map_p2s_anharm[atm1][0];
        const auto atm2_s = system->map_p2s_anharm[atm2][tran];

        for (i = 0; i < 3; ++i) {
            vec[i] = system->xr_s_anharm[atm2_s][i] + xshift_s[cell_s][i]
                     - system->xr_s_anharm[system->map_p2s_anharm[atm2][0]][i];
        }

        rotvec(vec, vec, system->lavec_s_anharm);
        rotvec(vec, vec, system->rlavec_p);

        double phase = vec[0] * xk_in[0] + vec[1] * xk_in[1] + vec[2] * xk_in[2];

        dphi2[3 * atm1 + xyz1][3 * atm2 + xyz2]
                += it.fcs_val * std::exp(im * phase)
                   / std::sqrt(system->mass_anharm[atm1_s] * system->mass_anharm[atm2_s]);

    }
}

void Gruneisen::prepare_delta_fcs(const std::vector<FcsArrayWithCell> &fcs_in,
                                  std::vector<FcsArrayWithCell> &delta_fcs,
                                  const int mirror_image_mode) const
{
    unsigned int i, j;
    double vec[3], vec_origin[3];
    double fcs_tmp = 0.0;

    std::vector<FcsAlignedForGruneisen> fcs_aligned;
    std::vector<AtomCellSuper> pairs_vec;
    std::vector<int> index_old, index_now;
    std::vector<int> index_with_cell, index_with_cell_old;
    std::set<std::vector<int>> set_index_uniq;
    AtomCellSuper pairs_tmp;

    const auto norder = fcs_in[0].pairs.size();
    unsigned int nmulti;

    delta_fcs.clear();
    fcs_aligned.clear();

    for (const auto &it: fcs_in) {
        fcs_aligned.emplace_back(it.fcs_val, it.pairs);
    }
    std::sort(fcs_aligned.begin(), fcs_aligned.end());

    if (mirror_image_mode == 0) {
        // original implementation
        // calculate IFC renormalization for the same atomic combination in the supercell
        // but with different mirror images at once and assign the average value for each
        // mirror images.
        index_old.clear();
        const auto nelems = 2 * (norder - 2) + 1;
        for (i = 0; i < nelems; ++i) index_old.push_back(-1);

        index_with_cell.clear();
        set_index_uniq.clear();

        for (const auto &it: fcs_aligned) {

            index_now.clear();
            index_with_cell.clear();

            index_now.push_back(it.pairs[0].index);
            index_with_cell.push_back(it.pairs[0].index);

            for (i = 1; i < norder - 1; ++i) {
                index_now.push_back(it.pairs[i].index);
                index_now.push_back(it.pairs[i].tran);

                index_with_cell.push_back(it.pairs[i].index);
                index_with_cell.push_back(it.pairs[i].tran);
                index_with_cell.push_back(it.pairs[i].cell_s);
            }

            if (index_now != index_old) {

                if (index_old[0] != -1) {

                    nmulti = set_index_uniq.size();
                    fcs_tmp /= static_cast<double>(nmulti);

                    if (std::abs(fcs_tmp) > eps15) {
                        for (const auto &it2: set_index_uniq) {

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

            for (i = 0; i < 3; i++) {
                vec_origin[i] = system->xr_s_anharm[system->map_p2s_anharm[it.pairs[0].index / 3][0]][i];
                for (j = 1; j < norder - 1; j++) {
                    vec_origin[i] +=
                            system->xr_s_anharm[system->map_p2s_anharm[it.pairs[j].index / 3][it.pairs[j].tran]][i]
                            + xshift_s[it.pairs[j].cell_s][i];
                }
                vec_origin[i] /= static_cast<double>(norder - 1);
            }

            for (i = 0; i < 3; ++i) {
                vec[i] = system->xr_s_anharm[system->map_p2s_anharm[it.pairs[norder - 1].index / 3][it.pairs[norder -
                                                                                                             1].tran]][i]
                         // - system->xr_s_anharm[system->map_p2s_anharm[it.pairs[0].index / 3][0]][i]
                         - vec_origin[i]
                         + xshift_s[it.pairs[norder - 1].cell_s][i];
            }

            rotvec(vec, vec, system->lavec_s_anharm);

            fcs_tmp += it.fcs_val * vec[it.pairs[norder - 1].index % 3];
        }
        nmulti = set_index_uniq.size();

        fcs_tmp /= static_cast<double>(nmulti);

        if (std::abs(fcs_tmp) > eps15) {
            for (const auto &it2: set_index_uniq) {

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
    } else { // if(mirror_image_mode != 0)
        // new implementation
        // calculate IFC renormalization separately for each mirror image combinations.
        index_with_cell_old.clear();
        const auto nelems = 3 * (norder - 2) + 1;
        for (i = 0; i < nelems; ++i) index_with_cell_old.push_back(-1);

        index_with_cell.clear();

        for (const auto &it: fcs_aligned) {

            index_with_cell.clear();

            index_with_cell.push_back(it.pairs[0].index);

            for (i = 1; i < norder - 1; ++i) {
                index_with_cell.push_back(it.pairs[i].index);
                index_with_cell.push_back(it.pairs[i].tran);
                index_with_cell.push_back(it.pairs[i].cell_s);
            }

            if (index_with_cell != index_with_cell_old) {

                if (index_with_cell_old[0] != -1) {

                    if (std::abs(fcs_tmp) > eps15) {

                        pairs_vec.clear();

                        pairs_tmp.index = index_with_cell_old[0];
                        pairs_tmp.tran = 0;
                        pairs_tmp.cell_s = 0;
                        pairs_vec.push_back(pairs_tmp);
                        for (i = 1; i < norder - 1; ++i) {
                            pairs_tmp.index = index_with_cell_old[3 * i - 2];
                            pairs_tmp.tran = index_with_cell_old[3 * i - 1];
                            pairs_tmp.cell_s = index_with_cell_old[3 * i];
                            pairs_vec.push_back(pairs_tmp);
                        }
                        delta_fcs.emplace_back(fcs_tmp, pairs_vec);
                    }
                }

                fcs_tmp = 0.0;
                index_with_cell_old.clear();
                index_with_cell_old.reserve(index_with_cell.size());
                std::copy(index_with_cell.begin(), index_with_cell.end(), std::back_inserter(index_with_cell_old));
            }

            for (i = 0; i < 3; i++) {
                vec_origin[i] = system->xr_s_anharm[system->map_p2s_anharm[it.pairs[0].index / 3][0]][i];
                for (j = 1; j < norder - 1; j++) {
                    vec_origin[i] +=
                            system->xr_s_anharm[system->map_p2s_anharm[it.pairs[j].index / 3][it.pairs[j].tran]][i]
                            + xshift_s[it.pairs[j].cell_s][i];
                }
                vec_origin[i] /= static_cast<double>(norder - 1);
            }

            for (i = 0; i < 3; ++i) {
                vec[i] = system->xr_s_anharm[system->map_p2s_anharm[it.pairs[norder - 1].index / 3][it.pairs[norder -
                                                                                                             1].tran]][i]
                         // - system->xr_s_anharm[system->map_p2s_anharm[it.pairs[0].index / 3][0]][i]
                         - vec_origin[i]
                         + xshift_s[it.pairs[norder - 1].cell_s][i];
            }

            rotvec(vec, vec, system->lavec_s_anharm);

            fcs_tmp += it.fcs_val * vec[it.pairs[norder - 1].index % 3];
        }

        if (std::abs(fcs_tmp) > eps15) {

            pairs_vec.clear();

            pairs_tmp.index = index_with_cell[0];
            pairs_tmp.tran = 0;
            pairs_tmp.cell_s = 0;
            pairs_vec.push_back(pairs_tmp);
            for (i = 1; i < norder - 1; ++i) {
                pairs_tmp.index = index_with_cell[3 * i - 2];
                pairs_tmp.tran = index_with_cell[3 * i - 1];
                pairs_tmp.cell_s = index_with_cell[3 * i];
                pairs_vec.push_back(pairs_tmp);
            }
            delta_fcs.emplace_back(fcs_tmp, pairs_vec);
        }
    }

    fcs_aligned.clear();
    set_index_uniq.clear();
}

// only mode == 0 is implemented
// This process breaks the permutation symmetry of the IFC.
// void Gruneisen::impose_ASR_on_harmonic_IFC(std::vector<FcsArrayWithCell> &delta_fcs,
//                            int mode)
// {
//     // should check if delta_fcs consists of harmonic IFCs?
// 
//     // change self-interaction
//     if(mode == 0){
// 
//         int nat = system->nat;
// 
//         int i_ind1, i_tran1, i_ind2, i_tran2, iat1, xyz1, xyz2;
//         double **phi;
//         allocate(phi, nat, 9);
// 
//         // calculate correction
//         for(auto &it: delta_fcs){
//             i_ind1 = it.pairs[0].index;
//             i_tran1 = it.pairs[0].tran;
//             i_ind2 = it.pairs[1].index;
//             iat1 = system->map_p2s_anharm[i_ind1/3][i_tran1];
//             xyz1 = i_ind1%3;
//             xyz2 = i_ind2%3;
// 
//             phi[iat1][xyz1*3 + xyz2] += it.fcs_val;
//         }
// 
//         // apply correction
//         for(auto &it: delta_fcs){
//             // if self-interaction
//             if((it.pairs[0].index/3 == it.pairs[1].index/3) && 
//                (it.pairs[0].tran == it.pairs[1].tran))
//             {
//                 i_ind1 = it.pairs[0].index;
//                 i_tran1 = it.pairs[0].tran;
//                 iat1 = system->map_p2s_anharm[i_ind1/3][i_tran1];
//                 xyz1 = i_ind1%3;
//                 xyz2 = i_ind2%3;
//                 it.fcs_val -= phi[iat1][xyz1*3 + xyz2];
//             }
//         }
// 
//     }
// 
//     return;
// }


void Gruneisen::write_new_fcsxml_all()
{
    std::cout << std::endl;

    if (fcs_phonon->update_fc2) {
        warn("write_new_fcsxml_all",
             "NEWFCS = 1 cannot be combined with the FC2XML.");
    } else {
        std::cout << " NEWFCS = 1 : Following XML files are created. " << std::endl;

        auto file_xml = input->job_title + "_+.xml";
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

void Gruneisen::write_new_fcsxml(const std::string &filename_xml,
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
        auto &child = pt.add("Data.Structure.Position.pos", str_tmp);
        child.put("<xmlattr>.index", i + 1);
        child.put("<xmlattr>.element", system->symbol_kd[system->kd[i]]);
    }

    pt.put("Data.Symmetry.NumberOfTranslations", system->ntran);
    for (i = 0; i < system->ntran; ++i) {
        for (j = 0; j < system->natmin; ++j) {
            auto &child = pt.add("Data.Symmetry.Translations.map",
                                 system->map_p2s[j][i] + 1);
            child.put("<xmlattr>.tran", i + 1);
            child.put("<xmlattr>.atom", j + 1);
        }
    }

    pt.put("Data.ForceConstants", "");
    str_tmp.clear();

    for (const auto &it: fcs_phonon->force_constant_with_cell[0]) {

        auto &child = pt.add("Data.ForceConstants.HARMONIC.FC2", double2string(it.fcs_val));

        child.put("<xmlattr>.pair1",
                  std::to_string(it.pairs[0].index / 3 + 1)
                  + " " + std::to_string(it.pairs[0].index % 3 + 1));
        child.put("<xmlattr>.pair2",
                  std::to_string(system->map_p2s[it.pairs[1].index / 3][it.pairs[1].tran] + 1)
                  + " " + std::to_string(it.pairs[1].index % 3 + 1)
                  + " " + std::to_string(it.pairs[1].cell_s + 1));
    }

    for (const auto &it: delta_fc2) {

        if (std::abs(it.fcs_val) < eps12) continue;

        auto &child = pt.add("Data.ForceConstants.HARMONIC.FC2",
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
        for (const auto &it: fcs_phonon->force_constant_with_cell[1]) {

            if (it.pairs[1].index > it.pairs[2].index) continue;

            auto &child = pt.add("Data.ForceConstants.ANHARM3.FC3",
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

        for (const auto &it: delta_fc3) {

            if (std::abs(it.fcs_val) < eps12) continue;

            if (it.pairs[1].index > it.pairs[2].index) continue;

            auto &child = pt.add("Data.ForceConstants.ANHARM3.FC3",
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
    const auto indent = 2;

#if BOOST_VERSION >= 105600
    write_xml(filename_xml, pt, std::locale(),
              xml_writer_make_settings<ptree::key_type>(' ', indent, widen<std::string>("utf-8")));
#else
    write_xml(filename_xml, pt, std::locale(),
        xml_writer_make_settings(' ', indent, widen<char>("utf-8")));
#endif
}

std::string Gruneisen::double2string(const double d) const
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
//     allocate(vec, norder, 3);
//     allocate(pos, norder, 3);
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
//     deallocate(vec);
//     deallocate(pos);
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
//     allocate(vec, norder, 3);
//     allocate(pos, norder, 3);
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
//     deallocate(vec);
//     deallocate(pos);
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
//     allocate(A, 3, 3, 3, 3);
//     allocate(C, 3, 3, 3, 3);
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
//     deallocate(A);
//     deallocate(C);
// }

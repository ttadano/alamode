/*
 writer.cpp

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include "writer.h"
#include "constraint.h"
#include "error.h"
#include "fcs.h"
#include "files.h"
#include "optimize.h"
#include "cluster.h"
#include "memory.h"
#include "patterndisp.h"
#include "symmetry.h"
#include "system.h"
#include "timer.h"
#include "version.h"
#include <iostream>
#include <fstream>
#include <map>
#include <unordered_map>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/version.hpp>

using namespace ALM_NS;

Writer::Writer() : output_maxorder(5)
{
    save_format_flags["alamode"] = 1;
    save_format_flags["shengbte"] = 0;
    save_format_flags["shengbte4"] = 0;
    save_format_flags["qefc"] = 0;
    save_format_flags["hessian"] = 0;
    filename_fcs = "";
}

Writer::~Writer() = default;

void Writer::write_input_vars(const System *system,
                              const Symmetry *symmetry,
                              const Cluster *cluster,
                              const Displace *displace,
                              const Fcs *fcs,
                              const Constraint *constraint,
                              const Optimize *optimize,
                              const Files *files,
                              const std::string run_mode) const
{
    size_t i;

    const auto nat = system->get_supercell().number_of_atoms;
    const auto nkd = system->get_supercell().number_of_elems;

    //alm->timer->start_clock("writer");

    std::cout << '\n';
    std::cout << " Input variables:\n";
    std::cout << " -------------------------------------------------------------------" << '\n';
    std::cout << " General:\n";
    std::cout << "  PREFIX = " << files->get_prefix() << '\n';
    std::cout << "  MODE = " << run_mode << '\n';
    std::cout << "  NAT = " << nat << "; NKD = " << nkd << '\n';
    std::cout << "  PRINTSYM = " << symmetry->get_print_symmetry()
              << "; TOLERANCE = " << symmetry->get_tolerance() << '\n';
    std::cout << "  KD = ";
    for (i = 0; i < nkd; ++i) std::cout << std::setw(4) << system->get_kdname()[i];
    std::cout << '\n';
    std::cout << "  FCSYM_BASIS = " << fcs->get_forceconstant_basis() << '\n';
    std::cout << "  PERIODIC = ";
    for (i = 0; i < 3; ++i) std::cout << std::setw(3) << system->get_periodicity()[i];
    std::cout << '\n';
    std::cout << "  MAGMOM = " << system->get_str_magmom() << '\n';
    std::cout << "  FCS_ALAMODE = " << save_format_flags.at("alamode") << ';';
    std::cout << "  NMAXSAVE = " << get_output_maxorder() << '\n';
    std::cout << "  FC3_SHENGBTE = " << save_format_flags.at("shengbte") << '\n';
    std::cout << "  FC4_SHENGBTE = " << save_format_flags.at("shengbte4") << '\n';
    std::cout << "  FC2_QEFC = " << save_format_flags.at("qefc") << '\n';
    std::cout << "  HESSIAN = " << save_format_flags.at("hessian") << '\n';
    std::cout << "  FC_ZERO_THR = " << fcs->get_fc_zero_threshold() << '\n';
    std::cout << '\n';

    std::cout << " Interaction:\n";
    std::cout << "  NORDER = " << cluster->get_maxorder() << '\n';
    std::cout << "  NBODY = ";
    for (auto m = 0; m < cluster->get_maxorder(); ++m) {
        std::cout << std::setw(3) << cluster->get_nbody_include()[m];
    }
    std::cout << "\n\n";

    if (run_mode == "suggest") {
        std::cout << "  DBASIS = " << displace->get_disp_basis() << "\n\n";

    } else if (run_mode == "optimize") {
        const auto optctrl = optimize->get_optimizer_control();
        std::vector<std::string> str_linearmodel{"least-squares", "elastic-net", "adaptive-lasso"};
        std::cout << " Optimize:\n";
        std::cout << "  LMODEL = "
                  << str_linearmodel[optctrl.linear_model - 1] << '\n';
        std::cout << "  DFSET = " << files->get_datfile_train().filename << '\n';
        std::cout << "  NDATA = " << files->get_datfile_train().ndata
                  << "; NSTART = " << files->get_datfile_train().nstart
                  << "; NEND = " << files->get_datfile_train().nend;
        if (files->get_datfile_train().skip_s
            < files->get_datfile_train().skip_e) {
            std::cout << "   SKIP = " << files->get_datfile_train().skip_s
                      << "-" << files->get_datfile_train().skip_e - 1 << "\n\n";
        } else {
            std::cout << "   SKIP = \n\n";
        }

        std::cout << "  ICONST = " << constraint->get_constraint_mode() << '\n';
        std::cout << "  ROTAXIS = " << constraint->get_rotation_axis() << '\n';
        std::cout << "  FC2XML = " << constraint->get_fc_file(2) << '\n';
        std::cout << "  FC3XML = " << constraint->get_fc_file(3) << "\n\n";
        std::cout << "  SPARSE = " << optctrl.use_sparse_solver << '\n';
        std::cout << "  SPARSESOLVER = " << optctrl.sparsesolver << '\n';
        std::cout << "  CONV_TOL = " << optctrl.tolerance_iteration << '\n';
        std::cout << "  MAXITER = " << optctrl.maxnum_iteration << "\n";
        std::cout << "  PERIODIC_IMAGE_CONV = " << optctrl.periodic_image_conv << "\n\n";
        if (optctrl.linear_model == 2) {
            std::cout << " Elastic-net related variables:\n";
            std::cout << "  CV = " << std::setw(5) << optctrl.cross_validation << '\n';
            std::cout << "  DFSET_CV = " << files->get_datfile_validation().filename << '\n';
            std::cout << "  NDATA_CV = " << files->get_datfile_validation().ndata
                      << "; NSTART_CV = " << files->get_datfile_validation().nstart
                      << "; NEND_CV = " << files->get_datfile_validation().nend << "\n\n";
            std::cout << "  L1_RATIO = " << optctrl.l1_ratio << '\n';
            std::cout << "  L1_ALPHA = " << optctrl.l1_alpha << '\n';
            std::cout << "  CV_MINALPHA = " << optctrl.l1_alpha_min
                      << "; CV_MAXALPHA = " << optctrl.l1_alpha_max
                      << ";  CV_NALPHA = " << optctrl.num_l1_alpha << '\n';
            std::cout << "  STANDARDIZE = " << optctrl.standardize << '\n';
            std::cout << "  ENET_DNORM = " << optctrl.displacement_normalization_factor << '\n';
            std::cout << "  NWRITE = " << std::setw(5) << optctrl.output_frequency << '\n';
            std::cout << "  DEBIAS_OLS = " << optctrl.debiase_after_l1opt << '\n';
            std::cout << '\n';
        }
    }
    std::cout << " -------------------------------------------------------------------\n\n";
    std::cout << std::flush;
    //alm->timer->stop_clock("writer");
}

void Writer::save_fcs_with_specific_format(const std::string fcs_format,
                                           const System *system,
                                           const Symmetry *symmetry,
                                           const Cluster *cluster,
                                           const Constraint *constraint,
                                           const Fcs *fcs,
                                           const Optimize *optimize,
                                           const Files *files,
                                           const int verbosity) const
{
    if (fcs_format == "alamode") {
        auto fname_save = get_filename_fcs();
        if (fname_save.empty()) {
            fname_save = files->get_prefix() + ".xml";
        }

        save_fcs_alamode_oldformat(system,
                                   symmetry,
                                   cluster,
                                   fcs,
                                   constraint,
                                   optimize->get_params(),
                                   files->get_datfile_train().filename,
                                   fname_save,
                                   verbosity);

    } else if (fcs_format == "shengbte") {

        if (cluster->get_maxorder() > 1) {
            auto fname_save = get_filename_fcs();
            if (fname_save.empty()) {
                fname_save = files->get_prefix() + ".FORCE_CONSTANT_3RD";;
            }

            save_fc3_shengbte_format(system,
                                     symmetry,
                                     cluster,
                                     constraint,
                                     fcs,
                                     fname_save,
                                     verbosity);
        }
    } else if (fcs_format == "shengbte4") {
        if (cluster->get_maxorder() > 2) {
            auto fname_save = get_filename_fcs();
            if (fname_save.empty()) {
                fname_save = files->get_prefix() + ".FORCE_CONSTANT_4TH";;
            }

            save_fc4_shengbte_format(system,
                                     symmetry,
                                     cluster,
                                     constraint,
                                     fcs,
                                     fname_save,
                                     verbosity);
        }
    } else if (fcs_format == "qefc") {
        auto fname_save = get_filename_fcs();
        if (fname_save.empty()) {
            fname_save = files->get_prefix() + ".fc";;
        }
        save_fc2_QEfc_format(system,
                             symmetry,
                             fcs,
                             fname_save,
                             verbosity);

    } else if (fcs_format == "hessian") {
        auto fname_save = get_filename_fcs();
        if (fname_save.empty()) {
            fname_save = files->get_prefix() + ".hessian";;
        }

        write_hessian(system,
                      symmetry,
                      fcs,
                      fname_save,
                      verbosity);
    }

}


void Writer::writeall(const System *system,
                      const Symmetry *symmetry,
                      const Cluster *cluster,
                      const Constraint *constraint,
                      const Fcs *fcs,
                      const Optimize *optimize,
                      const Files *files,
                      const int verbosity) const
{
//    alm->timer->start_clock("writer");

    if (verbosity > 0) {
        std::cout << " The following files are created:" << std::endl << std::endl;
    }

    const auto fname_save = files->get_prefix() + ".fcs";
    write_force_constants(cluster,
                          fcs,
                          symmetry,
                          optimize->get_params(),
                          verbosity,
                          fname_save);

    for (const auto &save_format_flag: save_format_flags) {
        if (save_format_flag.second) {
            save_fcs_with_specific_format(save_format_flag.first,
                                          system,
                                          symmetry,
                                          cluster,
                                          constraint,
                                          fcs,
                                          optimize,
                                          files,
                                          verbosity);
        }
    }

//    alm->timer->stop_clock("writer");
}

void Writer::write_force_constants(const Cluster *cluster,
                                   const Fcs *fcs,
                                   const Symmetry *symmetry,
                                   const double *fcs_vals,
                                   const int verbosity,
                                   const std::string fname_save) const
{
    int order, j, l;
    std::string *str_fcs;
    std::ofstream ofs_fcs;
    std::vector<int> atom_tmp;
    std::vector<std::vector<int>> cell_dummy;

    const auto maxorder = cluster->get_maxorder();

    ofs_fcs.open(fname_save.c_str(), std::ios::out);
    if (!ofs_fcs) exit("write_force_constants", "cannot open fcs file");

    ofs_fcs << " *********************** Force Constants (FCs) ***********************" << std::endl;
    ofs_fcs << " *        Force constants are printed in Rydberg atomic units.       *" << std::endl;
    ofs_fcs << " *        FC2: Ry/a0^2     FC3: Ry/a0^3     FC4: Ry/a0^4   etc.      *" << std::endl;
    ofs_fcs << " *        FC?: Ry/a0^?     a0 = Bohr radius                          *" << std::endl;
    ofs_fcs << " *                                                                   *" << std::endl;
    ofs_fcs << " *        The value shown in the last column is the distance         *" << std::endl;
    ofs_fcs << " *        between the most distant atomic pairs.                     *" << std::endl;
    ofs_fcs << " *********************************************************************" << std::endl;
    ofs_fcs << std::endl;
    ofs_fcs << " ----------------------------------------------------------------------" << std::endl;
    ofs_fcs << "      Index              FCs         P        Pairs     Distance [Bohr]" << std::endl;
    ofs_fcs << " (Global, Local)              (Multiplicity)                           " << std::endl;
    ofs_fcs << " ----------------------------------------------------------------------" << std::endl;

    allocate(str_fcs, maxorder);

    for (order = 0; order < maxorder; ++order) {
        str_fcs[order] = "*FC" + std::to_string(order + 2);
    }

    size_t k = 0;

    for (order = 0; order < maxorder; ++order) {

        size_t m = 0;

        if (!fcs->get_nequiv()[order].empty()) {

            ofs_fcs << std::endl << std::setw(6) << str_fcs[order] << std::endl;

            for (unsigned int ui = 0; ui < fcs->get_nequiv()[order].size(); ++ui) {

                ofs_fcs << std::setw(8) << k + 1 << std::setw(8) << ui + 1
                        << std::setw(18) << std::setprecision(7)
                        << std::scientific << fcs_vals[k];

                atom_tmp.clear();
                for (l = 1; l < order + 2; ++l) {
                    atom_tmp.push_back(fcs->get_fc_table()[order][m].elems[l] / 3);
                }
                j = symmetry->get_map_s2p()[fcs->get_fc_table()[order][m].elems[0] / 3].atom_num;
                std::sort(atom_tmp.begin(), atom_tmp.end());

                const auto iter_cluster
                        = cluster->get_interaction_cluster(order, j).
                                find(InteractionCluster(atom_tmp, cell_dummy));

                if (iter_cluster == cluster->get_interaction_cluster(order, j).end()) {
                    std::cout << std::setw(5) << j;
                    for (l = 0; l < order + 1; ++l) {
                        std::cout << std::setw(5) << atom_tmp[l];
                    }
                    std::cout << std::endl;
                    exit("write_force_constants",
                         "This cannot happen.");
                }

                const auto multiplicity = (*iter_cluster).cell.size();
                const auto distmax = (*iter_cluster).distmax;
                ofs_fcs << std::setw(4) << multiplicity;

                for (l = 0; l < order + 2; ++l) {
                    ofs_fcs << std::setw(7)
                            << easyvizint(fcs->get_fc_table()[order][m].elems[l]);
                }
                ofs_fcs << std::setw(12) << std::setprecision(3)
                        << std::fixed << distmax << std::endl;

                m += fcs->get_nequiv()[order][ui];
                ++k;
            }
        }
    }

    ofs_fcs << std::endl;

    ofs_fcs.unsetf(std::ios::showpos);

    for (order = 0; order < maxorder; ++order) {
        str_fcs[order] = "**FC" + std::to_string(order + 2);
    }

    ofs_fcs << std::endl << std::endl;
    ofs_fcs << " ------------------------ All FCs below ------------------------" << std::endl;

    auto ip = 0;

    for (order = 0; order < maxorder; ++order) {

        auto id = 0;

        if (!fcs->get_nequiv()[order].empty()) {
            ofs_fcs << std::endl << std::setw(6) << str_fcs[order] << std::endl;

            for (unsigned int iuniq = 0; iuniq < fcs->get_nequiv()[order].size(); ++iuniq) {

                auto str_tmp = "  # FC" + std::to_string(order + 2) + "_";
                str_tmp += std::to_string(iuniq + 1);

                ofs_fcs << str_tmp << std::setw(5) << fcs->get_nequiv()[order][iuniq]
                        << std::setw(16) << std::scientific
                        << std::setprecision(7) << fcs_vals[ip] << std::endl;

                for (j = 0; j < fcs->get_nequiv()[order][iuniq]; ++j) {
                    ofs_fcs << std::setw(5) << j + 1 << std::setw(12)
                            << std::setprecision(5) << std::fixed << fcs->get_fc_table()[order][id].sign;
                    for (k = 0; k < order + 2; ++k) {
                        ofs_fcs << std::setw(6)
                                << easyvizint(fcs->get_fc_table()[order][id].elems[k]);
                    }
                    ofs_fcs << std::endl;
                    ++id;
                }
                ofs_fcs << std::endl;
                ++ip;
            }
        }
    }
    deallocate(str_fcs);
    ofs_fcs.close();

    if (verbosity > 0) {
        std::cout << " Force constants in a human-readable format : "
                  << fname_save << std::endl;
    }
}

void Writer::write_displacement_pattern(const Cluster *cluster,
                                        const Displace *displace,
                                        const std::string prefix,
                                        const int verbosity) const
{
    const auto maxorder = cluster->get_maxorder();

    std::ofstream ofs_pattern;
    std::string file_disp_pattern;

    if (verbosity > 0) {
        std::cout << " Suggested displacement patterns are printed in the following files: " << std::endl;
    }

    for (auto order = 0; order < maxorder; ++order) {

        if (order == 0) {
            file_disp_pattern = prefix + ".pattern_HARMONIC";
        } else {
            file_disp_pattern = prefix + ".pattern_ANHARM"
                                + std::to_string(order + 2);
        }

        ofs_pattern.open(file_disp_pattern.c_str(), std::ios::out);
        if (!ofs_pattern) {
            exit("write_displacement_pattern",
                 "Cannot open file_disp_pattern");
        }

        auto counter = 0;

        ofs_pattern << "Basis : " << displace->get_disp_basis()[0] << std::endl;

        for (auto entry: displace->get_pattern_all(order)) {
            ++counter;

            ofs_pattern << std::setw(5) << counter << ":"
                        << std::setw(5) << entry.atoms.size() << std::endl;
            for (size_t i = 0; i < entry.atoms.size(); ++i) {
                ofs_pattern << std::setw(7) << entry.atoms[i] + 1;
                for (auto j = 0; j < 3; ++j) {
                    ofs_pattern << std::setw(15) << entry.directions[3 * i + j];
                }
                ofs_pattern << std::endl;
            }
        }

        ofs_pattern.close();

        if (verbosity > 0) {
            std::cout << "  " << cluster->get_ordername(order)
                      << " : " << file_disp_pattern << std::endl;
        }

    }
    if (verbosity > 0) std::cout << std::endl;
}


void Writer::save_fcs_alamode_oldformat(const System *system,
                                        const Symmetry *symmetry,
                                        const Cluster *cluster,
                                        const Fcs *fcs,
                                        const Constraint *constraint,
                                        const double *fcs_vals,
                                        const std::string fname_dfset,
                                        const std::string fname_fcs,
                                        const int verbosity) const
{
    SystemInfo system_structure;

    size_t i, j;

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            system_structure.lattice_vector[i][j]
                    = system->get_supercell().lattice_vector[i][j];
        }
    }

    system_structure.nat = system->get_supercell().number_of_atoms;
    system_structure.natmin = symmetry->get_nat_prim();
    system_structure.ntran = symmetry->get_ntran();
    system_structure.nspecies = system->get_supercell().number_of_elems;

    AtomProperty prop_tmp{};

    for (i = 0; i < system->get_supercell().number_of_atoms; ++i) {
        prop_tmp.x = system->get_supercell().x_fractional[i][0];
        prop_tmp.y = system->get_supercell().x_fractional[i][1];
        prop_tmp.z = system->get_supercell().x_fractional[i][2];
        prop_tmp.kind = system->get_supercell().kind[i];
        prop_tmp.atom = symmetry->get_map_s2p()[i].atom_num + 1;
        prop_tmp.tran = symmetry->get_map_s2p()[i].tran_num + 1;

        system_structure.atoms.emplace_back(AtomProperty(prop_tmp));
    }

    using boost::property_tree::ptree;

    ptree pt;
    std::string str_pos[3];

    pt.put("Data.ALM_version", ALAMODE_VERSION);
//    pt.put("Data.Optimize.DFSET", files->get_datfile_train().filename);
    pt.put("Data.Optimize.DFSET", fname_dfset);
    pt.put("Data.Optimize.Constraint", constraint->get_constraint_mode());

    pt.put("Data.Structure.NumberOfAtoms", system_structure.nat);
    pt.put("Data.Structure.NumberOfElements", system_structure.nspecies);

    for (i = 0; i < system_structure.nspecies; ++i) {
        auto &child = pt.add("Data.Structure.AtomicElements.element",
                             system->get_kdname()[i]);
        child.put("<xmlattr>.number", i + 1);
    }

    for (i = 0; i < 3; ++i) {
        str_pos[i].clear();
        for (j = 0; j < 3; ++j) {
            str_pos[i] += " " + double2string(system_structure.lattice_vector[j][i]);
        }
    }
    pt.put("Data.Structure.LatticeVector", "");
    pt.put("Data.Structure.LatticeVector.a1", str_pos[0]);
    pt.put("Data.Structure.LatticeVector.a2", str_pos[1]);
    pt.put("Data.Structure.LatticeVector.a3", str_pos[2]);

    std::stringstream ss;
    ss << system->get_periodicity()[0] << " "
       << system->get_periodicity()[1] << " "
       << system->get_periodicity()[2];
    pt.put("Data.Structure.Periodicity", ss.str());

    pt.put("Data.Structure.Position", "");
    std::string str_tmp;

    for (i = 0; i < system_structure.nat; ++i) {
        str_tmp.clear();
        for (j = 0; j < 3; ++j) str_tmp += " " + double2string(system->get_supercell().x_fractional[i][j]);
        auto &child = pt.add("Data.Structure.Position.pos", str_tmp);
        child.put("<xmlattr>.index", i + 1);
        child.put("<xmlattr>.element", system->get_kdname()[system->get_supercell().kind[i] - 1]);
    }

    pt.put("Data.Symmetry.NumberOfTranslations", symmetry->get_ntran());
    for (i = 0; i < system_structure.ntran; ++i) {
        for (j = 0; j < system_structure.natmin; ++j) {
            auto &child = pt.add("Data.Symmetry.Translations.map",
                                 symmetry->get_map_p2s()[j][i] + 1);
            child.put("<xmlattr>.tran", i + 1);
            child.put("<xmlattr>.atom", j + 1);
        }
    }

    if (system->get_spin().lspin) {
        pt.put("Data.MagneticMoments", "");
        pt.put("Data.MagneticMoments.Noncollinear", system->get_spin().noncollinear);
        pt.put("Data.MagneticMoments.TimeReversalSymmetry", system->get_spin().time_reversal_symm);
        for (i = 0; i < system_structure.nat; ++i) {
            str_tmp.clear();
            for (j = 0; j < 3; ++j) str_tmp += " " + double2string(system->get_spin().magmom[i][j], 5);
            auto &child = pt.add("Data.MagneticMoments.mag", str_tmp);
            child.put("<xmlattr>.index", i + 1);
        }
    }

    pt.put("Data.ForceConstants", "");
    str_tmp.clear();

    pt.put("Data.ForceConstants.HarmonicUnique.NFC2", fcs->get_nequiv()[0].size());
    pt.put("Data.ForceConstants.HarmonicUnique.Basis", fcs->get_forceconstant_basis());
    size_t ihead = 0;
    size_t k = 0;
    const auto nelem = cluster->get_maxorder() + 1;
    int *pair_tmp;
    std::vector<int> atom_tmp;
    std::vector<std::vector<int>> cell_dummy;
    std::set<InteractionCluster>::iterator iter_cluster;
    size_t multiplicity;


    allocate(pair_tmp, nelem);

    for (unsigned int ui = 0; ui < fcs->get_nequiv()[0].size(); ++ui) {

        for (i = 0; i < 2; ++i) {
            pair_tmp[i] = fcs->get_fc_table()[0][ihead].elems[i] / 3;
        }
        j = symmetry->get_map_s2p()[pair_tmp[0]].atom_num;

        atom_tmp.clear();
        atom_tmp.push_back(pair_tmp[1]);

        iter_cluster = cluster->get_interaction_cluster(0, j).find(
                InteractionCluster(atom_tmp, cell_dummy));
        if (iter_cluster == cluster->get_interaction_cluster(0, j).end()) {
            exit("load_reference_system_xml",
                 "Cubic force constant is not found.");
        }

        multiplicity = (*iter_cluster).cell.size();

        auto &child = pt.add("Data.ForceConstants.HarmonicUnique.FC2",
                             double2string(fcs_vals[k] * fcs->get_fc_table()[0][ihead].sign));
        child.put("<xmlattr>.pairs",
                  std::to_string(fcs->get_fc_table()[0][ihead].elems[0])
                  + " " + std::to_string(fcs->get_fc_table()[0][ihead].elems[1]));
        child.put("<xmlattr>.multiplicity", multiplicity);
        ihead += fcs->get_nequiv()[0][ui];
        ++k;
    }
    ihead = 0;


    if (cluster->get_maxorder() > 1) {

        pt.put("Data.ForceConstants.CubicUnique.NFC3", fcs->get_nequiv()[1].size());
        pt.put("Data.ForceConstants.CubicUnique.Basis", fcs->get_forceconstant_basis());

        for (unsigned int ui = 0; ui < fcs->get_nequiv()[1].size(); ++ui) {
            for (i = 0; i < 3; ++i) {
                pair_tmp[i] = fcs->get_fc_table()[1][ihead].elems[i] / 3;
            }
            j = symmetry->get_map_s2p()[pair_tmp[0]].atom_num;

            atom_tmp.clear();
            for (i = 1; i < 3; ++i) {
                atom_tmp.push_back(pair_tmp[i]);
            }
            std::sort(atom_tmp.begin(), atom_tmp.end());

            iter_cluster = cluster->get_interaction_cluster(1, j).find(
                    InteractionCluster(atom_tmp, cell_dummy));
            if (iter_cluster == cluster->get_interaction_cluster(1, j).end()) {
                exit("load_reference_system_xml",
                     "Cubic force constant is not found.");
            }
            multiplicity = (*iter_cluster).cell.size();


            auto &child = pt.add("Data.ForceConstants.CubicUnique.FC3",
                                 double2string(fcs_vals[k] * fcs->get_fc_table()[1][ihead].sign));
            child.put("<xmlattr>.pairs",
                      std::to_string(fcs->get_fc_table()[1][ihead].elems[0])
                      + " " + std::to_string(fcs->get_fc_table()[1][ihead].elems[1])
                      + " " + std::to_string(fcs->get_fc_table()[1][ihead].elems[2]));
            child.put("<xmlattr>.multiplicity", multiplicity);
            ihead += fcs->get_nequiv()[1][ui];
            ++k;
        }
    }

    int imult;
    std::string elementname = "Data.ForceConstants.HARMONIC.FC2";

    auto fc_cart_harmonic = fcs->get_fc_cart()[0];

    std::sort(fc_cart_harmonic.begin(),
              fc_cart_harmonic.end());

    for (const auto &it: fc_cart_harmonic) {

        for (k = 0; k < 2; ++k) {
            pair_tmp[k] = it.atoms[k];
        }

        j = symmetry->get_map_s2p()[pair_tmp[0]].atom_num;

        atom_tmp.clear();
        atom_tmp.push_back(pair_tmp[1]);

        iter_cluster = cluster->get_interaction_cluster(0, j).find(
                InteractionCluster(atom_tmp, cell_dummy));

        if (iter_cluster != cluster->get_interaction_cluster(0, j).end()) {
            multiplicity = (*iter_cluster).cell.size();

            for (imult = 0; imult < multiplicity; ++imult) {
                auto cell_now = (*iter_cluster).cell[imult];

                auto &child = pt.add(elementname,
                                     double2string(it.fc_value / static_cast<double>(multiplicity)));

                child.put("<xmlattr>.pair1", std::to_string(j + 1)
                                             + " " + std::to_string(it.coords[0] + 1));

                for (k = 1; k < 2; ++k) {
                    child.put("<xmlattr>.pair" + std::to_string(k + 1),
                              std::to_string(pair_tmp[k] + 1)
                              + " " + std::to_string(it.coords[k] + 1)
                              + " " + std::to_string(cell_now[k - 1] + 1));
                }
            }
        } else {
            exit("save_fcs_alamode_oldformat", "This cannot happen.");
        }
    }

    auto ishift = fcs->get_nequiv()[0].size();

    // Print anharmonic force constants to the xml file.

    for (auto order = 1; order < cluster->get_maxorder(); ++order) {

        if (order >= get_output_maxorder()) break;

        auto fc_cart_anharm = fcs->get_fc_cart()[order];

        std::sort(fc_cart_anharm.begin(),
                  fc_cart_anharm.end());

        for (const auto &it: fc_cart_anharm) {

            // Print force constants only when the coefficient is nonzero
            // and the last (order + 1) elements are sorted in ascending order.

            if (!it.is_ascending_order) continue;

            for (k = 0; k < order + 2; ++k) {
                pair_tmp[k] = it.atoms[k];
            }
            j = symmetry->get_map_s2p()[pair_tmp[0]].atom_num;

            atom_tmp.clear();

            for (k = 1; k < order + 2; ++k) {
                atom_tmp.push_back(pair_tmp[k]);
            }
            std::sort(atom_tmp.begin(), atom_tmp.end());

            elementname = "Data.ForceConstants.ANHARM"
                          + std::to_string(order + 2)
                          + ".FC" + std::to_string(order + 2);

            iter_cluster = cluster->get_interaction_cluster(order, j).find(
                    InteractionCluster(atom_tmp, cell_dummy));

            if (iter_cluster != cluster->get_interaction_cluster(order, j).end()) {
                multiplicity = (*iter_cluster).cell.size();

                for (imult = 0; imult < multiplicity; ++imult) {
                    auto cell_now = (*iter_cluster).cell[imult];

                    auto &child = pt.add(elementname,
                                         double2string(it.fc_value / static_cast<double>(multiplicity)));

                    child.put("<xmlattr>.pair1", std::to_string(j + 1)
                                                 + " " + std::to_string(it.coords[0] + 1));

                    for (k = 1; k < order + 2; ++k) {
                        child.put("<xmlattr>.pair" + std::to_string(k + 1),
                                  std::to_string(pair_tmp[k] + 1)
                                  + " " + std::to_string(it.coords[k] + 1)
                                  + " " + std::to_string(cell_now[k - 1] + 1));
                    }
                }
            } else {
                exit("save_fcs_alamode_oldformat", "This cannot happen.");
            }
        }
        ishift += fcs->get_nequiv()[order].size();
    }

    using namespace boost::property_tree::xml_parser;
    const auto indent = 2;

    //const auto file_xml = files->get_prefix() + ".xml";

#if BOOST_VERSION >= 105600
    write_xml(fname_fcs, pt, std::locale(),
              xml_writer_make_settings<ptree::key_type>(' ', indent,
                                                        widen<std::string>("utf-8")));
#else
    write_xml(fname_fcs, pt, std::locale(),
              xml_writer_make_settings(' ', indent, widen<char>("utf-8")));
#endif

    deallocate(pair_tmp);

    if (verbosity > 0) {
        std::cout << " Input data for the phonon code ANPHON      : " << fname_fcs << std::endl;
    }
}

void Writer::write_hessian(const System *system,
                           const Symmetry *symmetry,
                           const Fcs *fcs,
                           const std::string fname_out,
                           const int verbosity) const
{
    size_t i, j;
    int pair_tmp[2];
    int pair_tran[2];
    std::ofstream ofs_hes;
    double **hessian;

    const auto nat3 = 3 * system->get_supercell().number_of_atoms;

    allocate(hessian, nat3, nat3);

    for (i = 0; i < nat3; ++i) {
        for (j = 0; j < nat3; ++j) {
            hessian[i][j] = 0.0;
        }
    }

    for (const auto &it: fcs->get_fc_cart()[0]) {

        for (i = 0; i < 2; ++i) pair_tmp[i] = it.atoms[i];

        for (size_t itran = 0; itran < symmetry->get_ntran(); ++itran) {
            for (i = 0; i < 2; ++i) {
                pair_tran[i] = symmetry->get_map_sym()[pair_tmp[i]][symmetry->get_symnum_tran()[itran]];
            }
            hessian[3 * pair_tran[0] + it.coords[0]][3 * pair_tran[1] + it.coords[1]] = it.fc_value;
        }
    }

    ofs_hes.open(fname_out.c_str(), std::ios::out);
    if (!ofs_hes) exit("write_hessian", "cannot create hessian file");

    ofs_hes << "# atom1, xyz1, atom2, xyz2, FC2 (Ryd/Bohr^2)" << std::endl;
    for (i = 0; i < nat3; ++i) {
        for (j = 0; j < nat3; ++j) {
            ofs_hes << std::setw(5) << i / 3 + 1;
            ofs_hes << std::setw(5) << i % 3 + 1;
            ofs_hes << std::setw(6) << j / 3 + 1;
            ofs_hes << std::setw(5) << j % 3 + 1;
            ofs_hes << std::setw(25) << std::setprecision(15)
                    << std::scientific << hessian[i][j];
            ofs_hes << std::endl;
        }
    }
    ofs_hes.close();
    deallocate(hessian);

    if (verbosity) {
        std::cout << " Complete Hessian matrix                    : " << fname_out << std::endl;
    }
}

std::string Writer::double2string(const double d,
                                  const int nprec) const
{
    std::string rt;
    std::stringstream ss;

    ss << std::scientific << std::setprecision(nprec) << d;
    ss >> rt;
    return rt;
}

void Writer::save_fc2_QEfc_format(const System *system,
                                  const Symmetry *symmetry,
                                  const Fcs *fcs,
                                  const std::string fname_out,
                                  const int verbosity) const
{
    size_t i, j;
    int pair_tmp[2];
    int pair_tran[2];
    std::ofstream ofs_hes;
    double **hessian;
    const auto nat3 = 3 * system->get_supercell().number_of_atoms;

    allocate(hessian, nat3, nat3);

    for (i = 0; i < nat3; ++i) {
        for (j = 0; j < nat3; ++j) {
            hessian[i][j] = 0.0;
        }
    }
    for (const auto &it: fcs->get_fc_cart()[0]) {

        for (i = 0; i < 2; ++i) pair_tmp[i] = it.atoms[i];
        for (size_t itran = 0; itran < symmetry->get_ntran(); ++itran) {
            for (i = 0; i < 2; ++i) {
                pair_tran[i] = symmetry->get_map_sym()[pair_tmp[i]][symmetry->get_symnum_tran()[itran]];
            }
            hessian[3 * pair_tran[0] + it.coords[0]][3 * pair_tran[1] + it.coords[1]] = it.fc_value;
        }
    }

//    auto file_fc = files->get_prefix() + ".fc";

    ofs_hes.open(fname_out.c_str(), std::ios::out);
    if (!ofs_hes) exit("save_fc2_QEfc_format", "cannot create fc file");

    ofs_hes << "  1  1  1" << std::endl;
    for (auto icrd = 0; icrd < 3; ++icrd) {
        for (auto jcrd = 0; jcrd < 3; ++jcrd) {
            for (i = 0; i < system->get_supercell().number_of_atoms; ++i) {
                for (j = 0; j < system->get_supercell().number_of_atoms; ++j) {
                    ofs_hes << std::setw(6) << icrd + 1;
                    ofs_hes << std::setw(6) << jcrd + 1;
                    ofs_hes << std::setw(6) << i + 1;
                    ofs_hes << std::setw(6) << j + 1;
                    ofs_hes << std::endl;
                    ofs_hes << "  1  1  1 " << std::setw(20) << std::setprecision(13)
                            << std::scientific << hessian[3 * j + jcrd][3 * i + icrd];
                    ofs_hes << std::endl;
                }
            }
        }
    }
    ofs_hes.close();
    deallocate(hessian);

    if (verbosity) {
        std::cout << " Harmonic force constants in QE fc format     : " << fname_out << std::endl;
    }
}

void Writer::save_fc3_shengbte_format(const System *system,
                                      const Symmetry *symmetry,
                                      const Cluster *cluster,
                                      const Constraint *constraint,
                                      const Fcs *fcs,
                                      const std::string fname_out,
                                      const int verbosity) const
{
    size_t i, j, k;
    int pair_tmp[3], coord_tmp[3];
    std::ofstream ofs_fc3;
    size_t nelems = 0;
    const auto nat3 = 3 * system->get_supercell().number_of_atoms;
    const auto natmin = symmetry->get_nat_prim();
    const auto nat = system->get_supercell().number_of_atoms;
    const auto ntran = symmetry->get_ntran();

    std::vector<int> atom_tmp, flatten_array;
    std::vector<std::vector<int>> cell_dummy;
    std::set<InteractionCluster>::iterator iter_cluster;
    atom_tmp.resize(2);
    flatten_array.resize(2);
    cell_dummy.resize(2);

    std::unordered_map<std::tuple<size_t, size_t, size_t>, size_t> has_element;
    std::unordered_map<std::tuple<size_t, size_t, size_t>, double> fc3;

    double ***x_image = system->get_x_image();

    for (const auto &it: fcs->get_fc_cart()[1]) {

        if (!it.is_ascending_order) continue;

        for (i = 0; i < 3; ++i) {
            pair_tmp[i] = it.atoms[i];
            coord_tmp[i] = it.coords[i];
        }

        for (i = 0; i < 2; ++i) {
            atom_tmp[i] = pair_tmp[i + 1];
            flatten_array[i] = it.flattenarray[i + 1];
        }

        j = symmetry->get_map_s2p()[pair_tmp[0]].atom_num;
        iter_cluster = cluster->get_interaction_cluster(1, j).find(InteractionCluster(atom_tmp, cell_dummy));

        do {
            std::tuple<int, int, int> atom_tuple(j, flatten_array[0] / 3, flatten_array[1] / 3);
            if (has_element.find(atom_tuple) == has_element.end()) {
                nelems += (*iter_cluster).cell.size();
                has_element[atom_tuple] = 1;
            }
            fc3[std::make_tuple(3 * j + coord_tmp[0], flatten_array[0], flatten_array[1])] = it.fc_value;
        } while (std::next_permutation(flatten_array.begin(), flatten_array.end()));
    }

    ofs_fc3.open(fname_out.c_str(), std::ios::out);
    if (!ofs_fc3) exit("save_fc3_shengbte_format", "cannot create the file");
    ofs_fc3 << nelems << std::endl;

    double vecs[2][3];
    std::vector<int> atom_tmp_orig(2), index_sort_atom(2);
    std::vector<int> index_found(2);
    auto ielem = 0;
    const auto factor = Ryd / 1.6021766208e-19 / std::pow(Bohr_in_Angstrom, 3);

    for (i = 0; i < natmin; ++i) {
        for (auto jtran = 0; jtran < ntran; ++jtran) {
            for (j = 0; j < natmin; ++j) {
                for (auto ktran = 0; ktran < ntran; ++ktran) {
                    for (k = 0; k < natmin; ++k) {

                        const auto jat = symmetry->get_map_p2s()[j][jtran];
                        const auto kat = symmetry->get_map_p2s()[k][ktran];

                        if (has_element.find(std::make_tuple(i, jat, kat)) == has_element.end()) continue;

                        atom_tmp[0] = jat;
                        atom_tmp[1] = kat;
                        for (auto m = 0; m < 2; ++m) {
                            atom_tmp_orig[m] = atom_tmp[m];
                        }
                        std::sort(atom_tmp.begin(), atom_tmp.end());
                        for (auto m = 0; m < 2; ++m) index_found[m] = 0;
                        for (auto m = 0; m < 2; ++m) {
                            for (auto mm = 0; mm < 2; ++mm) {
                                if (index_found[mm] == 0) {
                                    if (atom_tmp_orig[m] == atom_tmp[mm]) {
                                        index_sort_atom[m] = mm;
                                        index_found[mm] = 1;
                                        break;
                                    }
                                }
                            }
                        }

                        iter_cluster = cluster->get_interaction_cluster(1, i).find(
                                InteractionCluster(atom_tmp, cell_dummy));
                        if (iter_cluster == cluster->get_interaction_cluster(1, i).end()) {
                            exit("save_fc3_shengbte_format", "This cannot happen.");
                        }

                        const auto multiplicity = (*iter_cluster).cell.size();

                        const auto jat0 = symmetry->get_map_p2s()[symmetry->get_map_s2p()[atom_tmp[0]].atom_num][0];
                        const auto kat0 = symmetry->get_map_p2s()[symmetry->get_map_s2p()[atom_tmp[1]].atom_num][0];

                        for (size_t imult = 0; imult < multiplicity; ++imult) {
                            auto cell_now = (*iter_cluster).cell[imult];

                            for (auto m = 0; m < 3; ++m) {
                                vecs[0][m] = (x_image[0][atom_tmp[0]][m]
                                              - x_image[0][jat0][m]
                                              + x_image[cell_now[0]][0][m]
                                              - x_image[0][0][m]) * Bohr_in_Angstrom;
                                vecs[1][m] = (x_image[0][atom_tmp[1]][m]
                                              - x_image[0][kat0][m]
                                              + x_image[cell_now[1]][0][m]
                                              - x_image[0][0][m]) * Bohr_in_Angstrom;
                            }

                            ++ielem;
                            ofs_fc3 << std::endl;
                            ofs_fc3 << ielem << std::endl;
                            ofs_fc3 << std::scientific;
                            ofs_fc3 << std::setprecision(10);

                            for (auto m = 0; m < 2; ++m) {
                                const auto m_map = index_sort_atom[m];
                                for (auto mm = 0; mm < 3; ++mm) {
                                    ofs_fc3 << std::setw(20) << vecs[m_map][mm];
                                }
                                ofs_fc3 << std::endl;
                            }

                            ofs_fc3 << std::setw(5) << i + 1;
                            ofs_fc3 << std::setw(5) << j + 1;
                            ofs_fc3 << std::setw(5) << k + 1 << std::endl;

                            for (auto ii = 0; ii < 3; ++ii) {
                                for (auto jj = 0; jj < 3; ++jj) {
                                    for (auto kk = 0; kk < 3; ++kk) {
                                        ofs_fc3 << std::setw(2) << ii + 1;
                                        ofs_fc3 << std::setw(3) << jj + 1;
                                        ofs_fc3 << std::setw(3) << kk + 1;
                                        if (fc3.find(std::make_tuple(3 * i + ii, 3 * jat + jj, 3 * kat + kk)) !=
                                            fc3.end()) {
                                            ofs_fc3 << std::setw(20)
                                                    << fc3[std::make_tuple(3 * i + ii, 3 * jat + jj, 3 * kat + kk)]
                                                       * factor / static_cast<double>(multiplicity) << std::endl;
                                        } else {
                                            ofs_fc3 << std::setw(20) << 0.0 << std::endl;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    ofs_fc3.close();

    if (verbosity) {
        std::cout << " Third-order FCs in ShengBTE format         : " << fname_out << std::endl;
    }
}


void Writer::save_fc4_shengbte_format(const System *system,
                                      const Symmetry *symmetry,
                                      const Cluster *cluster,
                                      const Constraint *constraint,
                                      const Fcs *fcs,
                                      const std::string fname_out,
                                      const int verbosity) const
{
    size_t i, j, k, l;
    int pair_tmp[4], coord_tmp[4];
    std::ofstream ofs_fc4;
    size_t nelems = 0;
    const auto nat3 = 3 * system->get_supercell().number_of_atoms;
    const auto natmin = symmetry->get_nat_prim();
    const auto nat = system->get_supercell().number_of_atoms;
    const auto ntran = symmetry->get_ntran();

    std::vector<int> atom_tmp, flatten_array;
    std::vector<std::vector<int>> cell_dummy;
    std::set<InteractionCluster>::iterator iter_cluster;
    atom_tmp.resize(3);
    flatten_array.resize(3);
    cell_dummy.resize(3);

    double ***x_image = system->get_x_image();

    if (cluster->get_maxorder() < 3) {
        exit("save_fc4_shengbte_format",
             "No fourth-order terms when NORDER < 3.");
    }

    std::unordered_map<std::tuple<size_t, size_t, size_t, size_t>, size_t> has_element;
    std::unordered_map<std::tuple<size_t, size_t, size_t, size_t>, double> fc4;

    for (const auto &it: fcs->get_fc_cart()[2]) {

        if (!it.is_ascending_order) continue;

        for (i = 0; i < 4; ++i) {
            pair_tmp[i] = it.atoms[i];
            coord_tmp[i] = it.coords[i];
        }

        j = symmetry->get_map_s2p()[pair_tmp[0]].atom_num;

        for (i = 0; i < 3; ++i) {
            atom_tmp[i] = pair_tmp[i + 1];
            flatten_array[i] = it.flattenarray[i + 1];
        }

        iter_cluster = cluster->get_interaction_cluster(2, j).find(InteractionCluster(atom_tmp, cell_dummy));

        do {
            std::tuple<size_t, size_t, size_t, size_t> atom_tuple(j,
                                                                  flatten_array[0] / 3,
                                                                  flatten_array[1] / 3,
                                                                  flatten_array[2] / 3);
            if (has_element.find(atom_tuple) == has_element.end()) {
                nelems += (*iter_cluster).cell.size();
                has_element[atom_tuple] = 1;
            }
            fc4[std::make_tuple(3 * j + coord_tmp[0],
                                flatten_array[0],
                                flatten_array[1],
                                flatten_array[2])] = it.fc_value;
        } while (std::next_permutation(flatten_array.begin(), flatten_array.end()));
    }

    ofs_fc4.open(fname_out.c_str(), std::ios::out);
    if (!ofs_fc4) exit("save_fc4_shengbte_format", "cannot create the file");
    ofs_fc4 << nelems << std::endl;

    double vecs[3][3];
    std::vector<int> atom_tmp_orig(3), index_sort_atom(3);
    std::vector<int> index_found(3);
    auto ielem = 0;
    const auto factor = Ryd / 1.6021766208e-19 / std::pow(Bohr_in_Angstrom, 4);

    for (i = 0; i < natmin; ++i) {
        for (auto jtran = 0; jtran < ntran; ++jtran) {
            for (j = 0; j < natmin; ++j) {
                for (auto ktran = 0; ktran < ntran; ++ktran) {
                    for (k = 0; k < natmin; ++k) {
                        for (auto ltran = 0; ltran < ntran; ++ltran) {
                            for (l = 0; l < natmin; ++l) {
                                const auto jat = symmetry->get_map_p2s()[j][jtran];
                                const auto kat = symmetry->get_map_p2s()[k][ktran];
                                const auto lat = symmetry->get_map_p2s()[l][ltran];

                                if (has_element.find(std::make_tuple(i, jat, kat, lat)) == has_element.end()) continue;

                                atom_tmp[0] = jat;
                                atom_tmp[1] = kat;
                                atom_tmp[2] = lat;
                                for (auto m = 0; m < 3; ++m) {
                                    atom_tmp_orig[m] = atom_tmp[m];
                                }
                                std::sort(atom_tmp.begin(), atom_tmp.end());
                                for (auto m = 0; m < 3; ++m) index_found[m] = 0;
                                for (auto m = 0; m < 3; ++m) {
                                    for (auto mm = 0; mm < 3; ++mm) {
                                        if (index_found[mm] == 0) {
                                            if (atom_tmp_orig[m] == atom_tmp[mm]) {
                                                index_sort_atom[m] = mm;
                                                index_found[mm] = 1;
                                                break;
                                            }
                                        }
                                    }
                                }

                                iter_cluster = cluster->get_interaction_cluster(2, i).find(
                                        InteractionCluster(atom_tmp, cell_dummy));
                                if (iter_cluster == cluster->get_interaction_cluster(2, i).end()) {
                                    exit("save_fc4_shengbte_format", "This cannot happen.");
                                }

                                const auto multiplicity = (*iter_cluster).cell.size();

                                const auto jat0 = symmetry->get_map_p2s()[symmetry->get_map_s2p()[atom_tmp[0]].atom_num][0];
                                const auto kat0 = symmetry->get_map_p2s()[symmetry->get_map_s2p()[atom_tmp[1]].atom_num][0];
                                const auto lat0 = symmetry->get_map_p2s()[symmetry->get_map_s2p()[atom_tmp[2]].atom_num][0];

                                for (size_t imult = 0; imult < multiplicity; ++imult) {
                                    auto cell_now = (*iter_cluster).cell[imult];

                                    for (auto m = 0; m < 3; ++m) {
                                        vecs[0][m] = (x_image[0][atom_tmp[0]][m]
                                                      - x_image[0][jat0][m]
                                                      + x_image[cell_now[0]][0][m]
                                                      - x_image[0][0][m]) * Bohr_in_Angstrom;
                                        vecs[1][m] = (x_image[0][atom_tmp[1]][m]
                                                      - x_image[0][kat0][m]
                                                      + x_image[cell_now[1]][0][m]
                                                      - x_image[0][0][m]) * Bohr_in_Angstrom;
                                        vecs[2][m] = (x_image[0][atom_tmp[2]][m]
                                                      - x_image[0][lat0][m]
                                                      + x_image[cell_now[2]][0][m]
                                                      - x_image[0][0][m]) * Bohr_in_Angstrom;
                                    }

                                    ++ielem;
                                    ofs_fc4 << std::endl;
                                    ofs_fc4 << ielem << std::endl;
                                    ofs_fc4 << std::scientific;
                                    ofs_fc4 << std::setprecision(10);

                                    for (auto m = 0; m < 3; ++m) {
                                        const auto m_map = index_sort_atom[m];
                                        for (auto mm = 0; mm < 3; ++mm) {
                                            ofs_fc4 << std::setw(20) << vecs[m_map][mm];
                                        }
                                        ofs_fc4 << std::endl;
                                    }

                                    ofs_fc4 << std::setw(5) << i + 1;
                                    ofs_fc4 << std::setw(5) << j + 1;
                                    ofs_fc4 << std::setw(5) << k + 1;
                                    ofs_fc4 << std::setw(5) << l + 1 << std::endl;

                                    for (auto ii = 0; ii < 3; ++ii) {
                                        for (auto jj = 0; jj < 3; ++jj) {
                                            for (auto kk = 0; kk < 3; ++kk) {
                                                for (auto ll = 0; ll < 3; ++ll) {
                                                    ofs_fc4 << std::setw(2) << ii + 1;
                                                    ofs_fc4 << std::setw(3) << jj + 1;
                                                    ofs_fc4 << std::setw(3) << kk + 1;
                                                    ofs_fc4 << std::setw(3) << ll + 1;
                                                    if (fc4.find(std::make_tuple(3 * i + ii,
                                                                                 3 * jat + jj,
                                                                                 3 * kat + kk,
                                                                                 3 * lat + ll)) != fc4.end()) {
                                                        ofs_fc4 << std::setw(20)
                                                                << fc4[std::make_tuple(3 * i + ii,
                                                                                       3 * jat + jj,
                                                                                       3 * kat + kk,
                                                                                       3 * lat + ll)]
                                                                   * factor / static_cast<double>(multiplicity)
                                                                << std::endl;
                                                    } else {
                                                        ofs_fc4 << std::setw(20) << 0.0 << std::endl;
                                                    }

                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    ofs_fc4.close();

    if (verbosity) {
        std::cout << " Fourth-order FCs in ShengBTE format        : " << fname_out << std::endl;
    }
}

std::string Writer::easyvizint(const int n) const
{
    const auto atmn = n / 3 + 1;
    const auto crdn = n % 3;
    std::string str_crd[3] = {"x", "y", "z"};
    auto str_tmp = std::to_string(atmn);
    str_tmp += str_crd[crdn];

    return str_tmp;
}

void Writer::set_fcs_save_flag(const std::string key_str,
                               const int val)
{
    if (save_format_flags.find(key_str) != save_format_flags.end()) {
        save_format_flags[key_str] = val;
    }
}

int Writer::get_fcs_save_flag(const std::string key_str)
{
    if (save_format_flags.find(key_str) != save_format_flags.end()) {
        return save_format_flags[key_str];
    }
    return -1;
}

void Writer::set_output_maxorder(const int maxorder)
{
    output_maxorder = maxorder;
}

int Writer::get_output_maxorder() const
{
    return output_maxorder;
}

void Writer::set_filename_fcs(const std::string filename_in)
{
    filename_fcs = filename_in;
}

std::string Writer::get_filename_fcs() const
{
    return filename_fcs;
}


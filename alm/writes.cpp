/*
 writes.cpp

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include <fstream>
#include <boost/lexical_cast.hpp>
#include "writes.h"
#include "system.h"
#include "interaction.h"
#include "memory.h"
#include "symmetry.h"
#include "error.h"
#include "files.h"
#include "fcs.h"
#include "fitting.h"
#include "constraint.h"
#include "input.h"
#include "timer.h"
#include "patterndisp.h"
#include "version.h"
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/version.hpp>

using namespace ALM_NS;

Writes::Writes(ALM *alm): Pointers(alm)
{
}

Writes::~Writes()
{
}

void Writes::write_input_vars()
{
    unsigned int i;

    std::cout << std::endl;
    std::cout << " Input variables:" << std::endl;
    std::cout << " -------------------------------------------------------------------" << std::endl;
    std::cout << " General:" << std::endl;
    std::cout << "  PREFIX = " << files->job_title << std::endl;
    std::cout << "  MODE = " << alm->mode << std::endl;
    std::cout << "  NAT = " << system->nat << "; NKD = " << system->nkd << std::endl;
    std::cout << "  NSYM = " << symmetry->nsym << "; PRINTSYM = " << symmetry->is_printsymmetry
        << "; TOLERANCE = " << symmetry->tolerance << std::endl;
    std::cout << "  KD = ";
    for (i = 0; i < system->nkd; ++i) std::cout << std::setw(4) << system->kdname[i];
    std::cout << std::endl;
    std::cout << "  PERIODIC = ";
    for (i = 0; i < 3; ++i) std::cout << std::setw(3) << interaction->is_periodic[i];
    std::cout << std::endl;
    std::cout << "  MAGMOM = " << input->str_magmom << std::endl;
    std::cout << "  HESSIAN = " << writes->print_hessian << std::endl;
    std::cout << std::endl;


    std::cout << " Interaction:" << std::endl;
    std::cout << "  NORDER = " << interaction->maxorder << std::endl;
    std::cout << "  NBODY = ";
    for (i = 0; i < interaction->maxorder; ++i)
        std::cout << std::setw(3) << interaction->nbody_include[i];

    std::cout << std::endl << std::endl;


    if (alm->mode == "suggest") {
        std::cout << "  DBASIS = " << displace->disp_basis << std::endl;
        std::cout << std::endl;

    } else if (alm->mode == "fitting") {
        std::cout << " Fitting:" << std::endl;
        std::cout << "  DFILE = " << files->file_disp << std::endl;
        std::cout << "  FFILE = " << files->file_force << std::endl;
        std::cout << "  NDATA = " << system->ndata << "; NSTART = " << system->nstart
            << "; NEND = " << system->nend << "; NSKIP = " << system->nskip << std::endl;
        std::cout << "  NBOOT = " << fitting->nboot << std::endl;
        std::cout << "  MULTDAT = " << symmetry->multiply_data << std::endl;
        std::cout << "  ICONST = " << constraint->constraint_mode << std::endl;
        std::cout << "  ROTAXIS = " << constraint->rotation_axis << std::endl;
        std::cout << "  FC2XML = " << constraint->fc2_file << std::endl;
        std::cout << "  FC3XML = " << constraint->fc3_file << std::endl;
        std::cout << std::endl;
    }
    std::cout << " -------------------------------------------------------------------" << std::endl;
    std::cout << std::endl;
}

void Writes::writeall()
{
    std::cout << " The following files are created:" << std::endl << std::endl;
    write_force_constants();
    write_misc_xml();
    if (print_hessian) write_hessian();
    std::cout << std::endl;
}

void Writes::write_force_constants()
{
    int order, j, k, l, m;
    unsigned int ui;
    int multiplicity;
    int maxorder = interaction->maxorder;
    double distmax;
    std::string *str_fcs;
    std::string str_tmp;
    std::ofstream ofs_fcs;
    std::vector<int> atom_tmp;
    std::vector<std::vector<int>> cell_dummy;
    std::set<MinimumDistanceCluster>::iterator iter_cluster;

    ofs_fcs.open(files->file_fcs.c_str(), std::ios::out);
    if (!ofs_fcs) error->exit("openfiles", "cannot open fcs file");

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

    memory->allocate(str_fcs, maxorder);

    for (order = 0; order < maxorder; ++order) {
        str_fcs[order] = "*FC" + boost::lexical_cast<std::string>(order + 2);
    }

    k = 0;

    for (order = 0; order < maxorder; ++order) {

        m = 0;

        if (fcs->ndup[order].size() > 0) {

            ofs_fcs << std::endl << std::setw(6) << str_fcs[order] << std::endl;

            for (ui = 0; ui < fcs->ndup[order].size(); ++ui) {

                ofs_fcs << std::setw(8) << k + 1 << std::setw(8) << ui + 1
                    << std::setw(18) << std::setprecision(7)
                    << std::scientific << fitting->params[k];

                atom_tmp.clear();
                for (l = 1; l < order + 2; ++l) {
                    atom_tmp.push_back(fcs->fc_set[order][m].elems[l] / 3);
                }
                j = symmetry->map_s2p[fcs->fc_set[order][m].elems[0] / 3].atom_num;
                std::sort(atom_tmp.begin(), atom_tmp.end());

                iter_cluster = interaction->mindist_cluster[order][j].find(
                    MinimumDistanceCluster(atom_tmp, cell_dummy));

                if (iter_cluster != interaction->mindist_cluster[order][j].end()) {
                    multiplicity = (*iter_cluster).cell.size();
                    distmax = (*iter_cluster).distmax;
                } else {
                    std::cout << std::setw(5) << j;
                    for (l = 0; l < order + 1; ++l) {
                        std::cout << std::setw(5) << atom_tmp[l];
                    }
                    std::cout << std::endl;
                    error->exit("write_force_constants",
                                "This cannot happen.");
                }
                ofs_fcs << std::setw(4) << multiplicity;

                for (l = 0; l < order + 2; ++l) {
                    ofs_fcs << std::setw(7)
                        << fcs->easyvizint(fcs->fc_set[order][m].elems[l]);
                }
                ofs_fcs << std::setw(12) << std::setprecision(3)
                    << std::fixed << distmax << std::endl;

                m += fcs->ndup[order][ui];
                ++k;
            }
        }
    }

    ofs_fcs << std::endl;

    if (constraint->extra_constraint_from_symmetry) {

        ofs_fcs << " -------------- Constraints from crystal symmetry --------------" << std::endl << std::endl;;
        for (order = 0; order < maxorder; ++order) {
            int nparam = fcs->ndup[order].size();


            for (std::vector<ConstraintClass>::iterator p = constraint->const_symmetry[order].begin();
                 p != constraint->const_symmetry[order].end();
                 ++p) {
                ofs_fcs << "   0 = " << std::scientific << std::setprecision(6);
                ConstraintClass const_pointer = *p;
                for (j = 0; j < nparam; ++j) {
                    if (std::abs(const_pointer.w_const[j]) > eps8) {
                        str_tmp = " * (FC" + boost::lexical_cast<std::string>(order + 2)
                            + "_" + boost::lexical_cast<std::string>(j + 1) + ")";
                        ofs_fcs << std::setw(10) << std::right
                            << std::showpos << const_pointer.w_const[j];
                        ofs_fcs << std::setw(12) << std::left << str_tmp;
                    }
                }
                ofs_fcs << std::endl;
            }
            ofs_fcs << std::endl;
        }
        ofs_fcs << std::endl;
    }

    ofs_fcs.unsetf(std::ios::showpos);

    for (order = 0; order < maxorder; ++order) {
        str_fcs[order] = "**FC" + boost::lexical_cast<std::string>(order + 2);
    }

    ofs_fcs << std::endl << std::endl;
    ofs_fcs << " ------------------------ All FCs below ------------------------" << std::endl;

    int ip = 0;
    int id;

    for (order = 0; order < maxorder; ++order) {

        id = 0;

        if (fcs->ndup[order].size() > 0) {
            ofs_fcs << std::endl << std::setw(6) << str_fcs[order] << std::endl;

            for (unsigned int iuniq = 0; iuniq < fcs->ndup[order].size(); ++iuniq) {

                str_tmp = "  # FC" + boost::lexical_cast<std::string>(order + 2) + "_";
                str_tmp += boost::lexical_cast<std::string>(iuniq + 1);

                ofs_fcs << str_tmp << std::setw(5) << fcs->ndup[order][iuniq]
                    << std::setw(16) << std::scientific
                    << std::setprecision(7) << fitting->params[ip] << std::endl;

                for (j = 0; j < fcs->ndup[order][iuniq]; ++j) {
                    ofs_fcs << std::setw(5) << j + 1 << std::setw(12)
                        << std::setprecision(5) << std::fixed << fcs->fc_set[order][id].coef;
                    for (k = 0; k < order + 2; ++k) {
                        ofs_fcs << std::setw(6)
                            << fcs->easyvizint(fcs->fc_set[order][id].elems[k]);
                    }
                    ofs_fcs << std::endl;
                    ++id;
                }
                ofs_fcs << std::endl;
                ++ip;
            }
        }
    }
    memory->deallocate(str_fcs);
    ofs_fcs.close();

    std::cout << " Force constants in a human-readable format : "
        << files->file_fcs << std::endl;
}

void Writes::write_displacement_pattern()
{
    int i, j;
    int order;
    int maxorder = interaction->maxorder;
    int counter;

    std::ofstream ofs_pattern;

    std::cout << " Suggested displacement patterns are printed in the following files: " << std::endl;

    for (order = 0; order < maxorder; ++order) {
        ofs_pattern.open(files->file_disp_pattern[order].c_str(), std::ios::out);
        if (!ofs_pattern)
            error->exit("write_displacement_pattern",
                        "Cannot open file_disp_pattern");

        counter = 0;

        ofs_pattern << "Basis : " << displace->disp_basis[0] << std::endl;

        for (std::vector<AtomWithDirection>::iterator it = displace->pattern_all[order].begin();
             it != displace->pattern_all[order].end(); ++it) {
            AtomWithDirection entry = *it;

            ++counter;

            ofs_pattern << std::setw(5) << counter << ":"
                << std::setw(5) << entry.atoms.size() << std::endl;
            for (i = 0; i < entry.atoms.size(); ++i) {
                ofs_pattern << std::setw(7) << entry.atoms[i] + 1;
                for (j = 0; j < 3; ++j) {
                    ofs_pattern << std::setw(15) << entry.directions[3 * i + j];
                }
                ofs_pattern << std::endl;
            }
        }

        ofs_pattern.close();

        std::cout << "  " << interaction->str_order[order]
            << " : " << files->file_disp_pattern[order] << std::endl;
    }
    std::cout << std::endl;
}

void Writes::write_misc_xml()
{
    SystemInfo system_structure;

    int i, j;

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            system_structure.lattice_vector[i][j] = system->lavec[i][j];
        }
    }

    system_structure.nat = system->nat;
    system_structure.natmin = symmetry->natmin;
    system_structure.ntran = symmetry->ntran;
    system_structure.nspecies = system->nkd;

    AtomProperty prop_tmp;

    for (i = 0; i < system->nat; ++i) {
        prop_tmp.x = system->xcoord[i][0];
        prop_tmp.y = system->xcoord[i][1];
        prop_tmp.z = system->xcoord[i][2];
        prop_tmp.kind = system->kd[i];
        prop_tmp.atom = symmetry->map_s2p[i].atom_num + 1;
        prop_tmp.tran = symmetry->map_s2p[i].tran_num + 1;

        system_structure.atoms.push_back(AtomProperty(prop_tmp));
    }

    using boost::property_tree::ptree;

    ptree pt;
    std::string str_pos[3];

    pt.put("Data.ALM_version", ALAMODE_VERSION);
    pt.put("Data.Fitting.DisplaceFile", files->file_disp);
    pt.put("Data.Fitting.ForceFile", files->file_force);
    pt.put("Data.Fitting.Constraint", constraint->constraint_mode);

    pt.put("Data.Structure.NumberOfAtoms", system_structure.nat);
    pt.put("Data.Structure.NumberOfElements", system_structure.nspecies);

    for (i = 0; i < system_structure.nspecies; ++i) {
        ptree &child = pt.add("Data.Structure.AtomicElements.element",
                              system->kdname[i]);
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
    ss << interaction->is_periodic[0] << " "
        << interaction->is_periodic[1] << " "
        << interaction->is_periodic[2];
    pt.put("Data.Structure.Periodicity", ss.str());

    pt.put("Data.Structure.Position", "");
    std::string str_tmp;

    for (i = 0; i < system_structure.nat; ++i) {
        str_tmp.clear();
        for (j = 0; j < 3; ++j) str_tmp += " " + double2string(system->xcoord[i][j]);
        ptree &child = pt.add("Data.Structure.Position.pos", str_tmp);
        child.put("<xmlattr>.index", i + 1);
        child.put("<xmlattr>.element", system->kdname[system->kd[i] - 1]);
    }

    pt.put("Data.Symmetry.NumberOfTranslations", symmetry->ntran);
    for (i = 0; i < system_structure.ntran; ++i) {
        for (j = 0; j < system_structure.natmin; ++j) {
            ptree &child = pt.add("Data.Symmetry.Translations.map",
                                  symmetry->map_p2s[j][i] + 1);
            child.put("<xmlattr>.tran", i + 1);
            child.put("<xmlattr>.atom", j + 1);
        }
    }

    if (system->lspin) {
        pt.put("Data.MagneticMoments", "");
        pt.put("Data.MagneticMoments.Noncollinear", system->noncollinear);
        pt.put("Data.MagneticMoments.TimeReversalSymmetry", symmetry->trev_sym_mag);
        for (i = 0; i < system_structure.nat; ++i) {
            str_tmp.clear();
            for (j = 0; j < 3; ++j) str_tmp += " " + double2string(system->magmom[i][j], 5);
            ptree &child = pt.add("Data.MagneticMoments.mag", str_tmp);
            child.put("<xmlattr>.index", i + 1);
        }
    }

    pt.put("Data.ForceConstants", "");
    str_tmp.clear();

    pt.put("Data.ForceConstants.HarmonicUnique.NFC2", fcs->ndup[0].size());

    int ihead = 0;
    int k = 0;
    int nelem = interaction->maxorder + 1;
    int *pair_tmp;

    memory->allocate(pair_tmp, nelem);

    for (unsigned int ui = 0; ui < fcs->ndup[0].size(); ++ui) {

        for (i = 0; i < 2; ++i) {
            pair_tmp[i] = fcs->fc_set[0][ihead].elems[i] / 3;
        }
        j = symmetry->map_s2p[pair_tmp[0]].atom_num;

        ptree &child = pt.add("Data.ForceConstants.HarmonicUnique.FC2",
                              double2string(fitting->params[k]));
        child.put("<xmlattr>.pairs",
                  boost::lexical_cast<std::string>(fcs->fc_set[0][ihead].elems[0])
                  + " " + boost::lexical_cast<std::string>(fcs->fc_set[0][ihead].elems[1]));
        child.put("<xmlattr>.multiplicity",
                  interaction->mindist_pairs[pair_tmp[0]][pair_tmp[1]].size());
        ihead += fcs->ndup[0][ui];
        ++k;
    }
    ihead = 0;

    std::vector<int> atom_tmp;
    std::vector<std::vector<int>> cell_dummy;
    std::set<MinimumDistanceCluster>::iterator iter_cluster;
    int multiplicity;

    if (interaction->maxorder > 1) {

        pt.put("Data.ForceConstants.CubicUnique.NFC3", fcs->ndup[1].size());

        for (unsigned int ui = 0; ui < fcs->ndup[1].size(); ++ui) {
            for (i = 0; i < 3; ++i) {
                pair_tmp[i] = fcs->fc_set[1][ihead].elems[i] / 3;
            }
            j = symmetry->map_s2p[pair_tmp[0]].atom_num;

            atom_tmp.clear();
            for (i = 1; i < 3; ++i) {
                atom_tmp.push_back(pair_tmp[i]);
            }
            std::sort(atom_tmp.begin(), atom_tmp.end());

            iter_cluster = interaction->mindist_cluster[1][j].find(
                MinimumDistanceCluster(atom_tmp, cell_dummy));
            if (iter_cluster == interaction->mindist_cluster[1][j].end()) {
                error->exit("load_reference_system_xml",
                            "Cubic force constant is not found.");
            } else {
                multiplicity = (*iter_cluster).cell.size();
            }

            ptree &child = pt.add("Data.ForceConstants.CubicUnique.FC3",
                                  double2string(fitting->params[k]));
            child.put("<xmlattr>.pairs",
                      boost::lexical_cast<std::string>(fcs->fc_set[1][ihead].elems[0])
                      + " " + boost::lexical_cast<std::string>(fcs->fc_set[1][ihead].elems[1])
                      + " " + boost::lexical_cast<std::string>(fcs->fc_set[1][ihead].elems[2]));
            child.put("<xmlattr>.multiplicity", multiplicity);
            ihead += fcs->ndup[1][ui];
            ++k;
        }
    }

    int ip, ishift;

    std::sort(fcs->fc_set[0].begin(), fcs->fc_set[0].end());

    for (std::vector<FcProperty>::iterator it = fcs->fc_set[0].begin();
         it != fcs->fc_set[0].end(); ++it) {
        FcProperty fctmp = *it;
        ip = fctmp.mother;

        for (k = 0; k < 2; ++k) {
            pair_tmp[k] = fctmp.elems[k] / 3;
        }
        j = symmetry->map_s2p[pair_tmp[0]].atom_num;
        for (std::vector<DistInfo>::iterator it2 = interaction->mindist_pairs[pair_tmp[0]][pair_tmp[1]].begin();
             it2 != interaction->mindist_pairs[pair_tmp[0]][pair_tmp[1]].end(); ++it2) {
            ptree &child = pt.add("Data.ForceConstants.HARMONIC.FC2",
                                  double2string(fitting->params[ip] * fctmp.coef
                                      / static_cast<double>(interaction->mindist_pairs[pair_tmp[0]][pair_tmp[1]].size())));

            child.put("<xmlattr>.pair1", boost::lexical_cast<std::string>(j + 1)
                      + " " + boost::lexical_cast<std::string>(fctmp.elems[0] % 3 + 1));
            child.put("<xmlattr>.pair2", boost::lexical_cast<std::string>(pair_tmp[1] + 1)
                      + " " + boost::lexical_cast<std::string>(fctmp.elems[1] % 3 + 1)
                      + " " + boost::lexical_cast<std::string>((*it2).cell + 1));
        }
    }

    ishift = fcs->ndup[0].size();

    // Print anharmonic force constants to the xml file.

    int imult;

    int order;
    std::string elementname;
    for (order = 1; order < interaction->maxorder; ++order) {

        std::sort(fcs->fc_set[order].begin(), fcs->fc_set[order].end());

        for (std::vector<FcProperty>::iterator it = fcs->fc_set[order].begin();
             it != fcs->fc_set[order].end(); ++it) {
            FcProperty fctmp = *it;
            ip = fctmp.mother + ishift;

            for (k = 0; k < order + 2; ++k) {
                pair_tmp[k] = fctmp.elems[k] / 3;
            }
            j = symmetry->map_s2p[pair_tmp[0]].atom_num;

            atom_tmp.clear();

            for (k = 1; k < order + 2; ++k) {
                atom_tmp.push_back(pair_tmp[k]);
            }
            std::sort(atom_tmp.begin(), atom_tmp.end());

            elementname = "Data.ForceConstants.ANHARM"
                + boost::lexical_cast<std::string>(order + 2)
                + ".FC" + boost::lexical_cast<std::string>(order + 2);


            iter_cluster = interaction->mindist_cluster[order][j].find(
                MinimumDistanceCluster(atom_tmp, cell_dummy));

            if (iter_cluster != interaction->mindist_cluster[order][j].end()) {
                multiplicity = (*iter_cluster).cell.size();

                for (imult = 0; imult < multiplicity; ++imult) {
                    std::vector<int> cell_now = (*iter_cluster).cell[imult];

                    ptree &child = pt.add(elementname,
                                          double2string(fitting->params[ip] * fctmp.coef
                                              / static_cast<double>(multiplicity)));

                    child.put("<xmlattr>.pair1", boost::lexical_cast<std::string>(j + 1)
                              + " " + boost::lexical_cast<std::string>(fctmp.elems[0] % 3 + 1));

                    for (k = 1; k < order + 2; ++k) {
                        child.put("<xmlattr>.pair" + boost::lexical_cast<std::string>(k + 1),
                                  boost::lexical_cast<std::string>(pair_tmp[k] + 1)
                                  + " " + boost::lexical_cast<std::string>(fctmp.elems[k] % 3 + 1)
                                  + " " + boost::lexical_cast<std::string>(cell_now[k - 1] + 1));
                    }
                }
            } else {
                error->exit("write_misc_xml", "This cannot happen.");
            }
        }
        ishift += fcs->ndup[order].size();
    }

    using namespace boost::property_tree::xml_parser;
    const int indent = 2;

    std::string file_xml = files->job_title + ".xml";

#if BOOST_VERSION >= 105600
    write_xml(file_xml, pt, std::locale(),
              xml_writer_make_settings<ptree::key_type>(' ', indent,
                                                        widen<std::string>("utf-8")));
#else
    write_xml(file_xml, pt, std::locale(),
        xml_writer_make_settings(' ', indent, widen<char>("utf-8")));
#endif

    memory->deallocate(pair_tmp);

    std::cout << " Input data for the phonon code ANPHON      : " << file_xml << std::endl;
}

void Writes::write_hessian()
{
    int i, j, itran, ip;
    int pair_tmp[2];
    int pair_tran[2];
    std::ofstream ofs_hes;
    double **hessian;
    int nat3 = 3 * system->nat;

    memory->allocate(hessian, nat3, nat3);

    for (i = 0; i < nat3; ++i) {
        for (j = 0; j < nat3; ++j) {
            hessian[i][j] = 0.0;
        }
    }

    for (std::vector<FcProperty>::iterator it = fcs->fc_set[0].begin();
         it != fcs->fc_set[0].end(); ++it) {
        FcProperty fctmp = *it;
        ip = fctmp.mother;

        for (i = 0; i < 2; ++i) pair_tmp[i] = fctmp.elems[i] / 3;
        for (itran = 0; itran < symmetry->ntran; ++itran) {
            for (i = 0; i < 2; ++i) {
                pair_tran[i] = symmetry->map_sym[pair_tmp[i]][symmetry->symnum_tran[itran]];
            }
            hessian[3 * pair_tran[0] + fctmp.elems[0] % 3][3 * pair_tran[1] + fctmp.elems[1] % 3]
                = fitting->params[ip] * fctmp.coef;
        }
    }

    ofs_hes.open(files->file_hes.c_str(), std::ios::out);
    if (!ofs_hes) error->exit("write_hessian", "cannot create hessian file");

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
    memory->deallocate(hessian);

    std::cout << " Complete Hessian matrix                    : " << files->file_hes << std::endl;

    /*
      std::string file_fc2 = files->job_title + ".fc2";
      std::ofstream ofs_fc2;
      ofs_fc2.open(file_fc2.c_str(), std::ios::out);
      ofs_fc2 << " # iat, icrd, jat, icrd, icell, relvec, fc2" << std::endl;
      double vec[3];
      for (std::vector<FcProperty>::iterator it = fcs->fc_set[0].begin();
           it != fcs->fc_set[0].end(); ++it) {
          FcProperty fctmp = *it;
          ip = fctmp.mother;
  
          for (i = 0; i < 2; ++i) pair_tmp[i] = fctmp.elems[i] / 3;
          for (itran = 0; itran < symmetry->ntran; ++itran) {
              for (i = 0; i < 2; ++i) {
                  pair_tran[i] = symmetry->map_sym[pair_tmp[i]][symmetry->symnum_tran[itran]];
              }
              for (std::vector<DistInfo>::iterator
                   it2  = interaction->mindist_pairs[pair_tran[0]][pair_tran[1]].begin();
                   it2 != interaction->mindist_pairs[pair_tran[0]][pair_tran[1]].end(); ++it2) {
                    int multiplicity = interaction->mindist_pairs[pair_tran[0]][pair_tran[1]].size();
                    for (i = 0; i < 3; ++i) {
                      vec[i] = interaction->x_image[(*it2).cell][pair_tran[1]][i]
                             - interaction->x_image[0][pair_tran[0]][i];
                    }
                    ofs_fc2 << std::setw(5) << pair_tran[0] + 1 << std::setw(5) << fctmp.elems[0] % 3 + 1;
                    ofs_fc2 << std::setw(5) << pair_tran[1] + 1 << std::setw(5) << fctmp.elems[1] % 3 + 1;
                    ofs_fc2 << std::setw(5) << (*it2).cell + 1;
                    ofs_fc2 << std::setw(15) << vec[0];
                    ofs_fc2 << std::setw(15) << vec[1];
                    ofs_fc2 << std::setw(15) << vec[2];
  
                    ofs_fc2 << std::setw(15)
                    << fitting->params[ip] * fctmp.coef / static_cast<double>(multiplicity);
                    ofs_fc2 << std::endl;
              }
          }
      }
      ofs_fc2.close();
     */
}

std::string Writes::double2string(const double d, const int nprec)
{
    std::string rt;
    std::stringstream ss;

    ss << std::scientific << std::setprecision(nprec) << d;
    ss >> rt;
    return rt;
}

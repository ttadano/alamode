/*
system.cpp

Copyright (c) 2014, 2015, 2016 Terumasa Tadano

This file is distributed under the terms of the MIT license.
Please see the file 'LICENCE.txt' in the root directory 
or http://opensource.org/licenses/mit-license.php for information.
*/

#include "mpi_common.h"
#include "system.h"
#include "fcs_phonon.h"
#include "error.h"
#include "memory.h"
#include "constants.h"
#include "symmetry_core.h"
#include <string>
#include <iostream>
#include <iomanip>
#include "mathfunctions.h"
#include "xml_parser.h"
#include <sstream>
#include <map>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>

using namespace PHON_NS;

System::System(PHON *phon) : Pointers(phon)
{
    set_default_variables();
}

System::~System()
{
    deallocate_variables();
}

void System::set_default_variables()
{
    xr_p = nullptr;
    xr_s = nullptr;
    xc = nullptr;
    xr_s_anharm = nullptr;
    kd = nullptr;
    kd_anharm = nullptr;
    mass_kd = nullptr;
    mass = nullptr;
    mass_anharm = nullptr;
    symbol_kd = nullptr;
    map_p2s = nullptr;
    map_p2s_anharm = nullptr;
    map_p2s_anharm_orig = nullptr;
    map_s2p = nullptr;
    map_s2p_anharm = nullptr;
    magmom = nullptr;
    atomlist_class = nullptr;
}

void System::deallocate_variables()
{
    if (xr_p) {
        memory->deallocate(xr_p);
    }
    if (xr_s) {
        memory->deallocate(xr_s);
    }
    if (xc) {
        memory->deallocate(xc);
    }
    if (xr_s_anharm) {
        memory->deallocate(xr_s_anharm);
    }
    if (kd) {
        memory->deallocate(kd);
    }
    if (kd_anharm) {
        memory->deallocate(kd_anharm);
    }
    if (mass_kd) {
        memory->deallocate(mass_kd);
    }
    if (mass) {
        memory->deallocate(mass);
    }
    if (mass_anharm) {
        memory->deallocate(mass_anharm);
    }
    if (symbol_kd) {
        memory->deallocate(symbol_kd);
    }
    if (map_p2s) {
        memory->deallocate(map_p2s);
    }
    if (map_p2s_anharm) {
        memory->deallocate(map_p2s_anharm);
    }
    if (map_s2p) {
        memory->deallocate(map_s2p);
    }
    if (map_s2p_anharm) {
        memory->deallocate(map_s2p_anharm);
    }
    if (magmom) {
        memory->deallocate(magmom);
    }
    if (map_p2s_anharm_orig) {
        memory->deallocate(map_p2s_anharm_orig);
    }
    if (atomlist_class) {
        memory->deallocate(atomlist_class);
    }
}


void System::setup()
{
    using namespace std;

    unsigned int i, j;
    double vec_tmp[3][3];
    unsigned int *kd_prim;
    double **xtmp;

    if (mympi->my_rank == 0) {

        if (!mass_kd) {
            memory->allocate(mass_kd, nkd);
            set_mass_elem_from_database(nkd, symbol_kd, mass_kd);
        }
    }
    load_system_info_from_XML();

    recips(lavec_s, rlavec_s);
    recips(lavec_s_anharm, rlavec_s_anharm);
    recips(lavec_p, rlavec_p);

    memory->allocate(xr_p, nat, 3);
    memory->allocate(xc, nat, 3);

    for (i = 0; i < nat; ++i) {
        rotvec(xc[i], xr_s[i], lavec_s);
        rotvec(xr_p[i], xc[i], rlavec_p);
        for (j = 0; j < 3; ++j) {
            xr_p[i][j] /= 2.0 * pi;
        }
    }

    if (mympi->my_rank == 0) {
        cout << " -----------------------------------------------------------------" << endl;
        cout << endl;
        cout << " Crystal structure" << endl;
        cout << " =================" << endl << endl;
        cout << " Lattice Vectors:" << endl << endl;
        cout.setf(ios::scientific);

        cout << " * Supercell (from " << fcs_phonon->file_fcs << " )" << endl << endl;

        cout << setw(16) << lavec_s_anharm[0][0];
        cout << setw(15) << lavec_s_anharm[1][0];
        cout << setw(15) << lavec_s_anharm[2][0];
        cout << " : a1" << endl;

        cout << setw(16) << lavec_s_anharm[0][1];
        cout << setw(15) << lavec_s_anharm[1][1];
        cout << setw(15) << lavec_s_anharm[2][1];
        cout << " : a2" << endl;

        cout << setw(16) << lavec_s_anharm[0][2];
        cout << setw(15) << lavec_s_anharm[1][2];
        cout << setw(15) << lavec_s_anharm[2][2];
        cout << " : a3" << endl;
        cout << endl;

        cout << setw(16) << rlavec_s_anharm[0][0];
        cout << setw(15) << rlavec_s_anharm[0][1];
        cout << setw(15) << rlavec_s_anharm[0][2];
        cout << " : b1" << endl;

        cout << setw(16) << rlavec_s_anharm[1][0];
        cout << setw(15) << rlavec_s_anharm[1][1];
        cout << setw(15) << rlavec_s_anharm[1][2];
        cout << " : b2" << endl;

        cout << setw(16) << rlavec_s_anharm[2][0];
        cout << setw(15) << rlavec_s_anharm[2][1];
        cout << setw(15) << rlavec_s_anharm[2][2];
        cout << " : b3" << endl;
        cout << endl;

        cout << " * Primitive cell " << endl << endl;
        cout << setw(16) << lavec_p[0][0];
        cout << setw(15) << lavec_p[1][0];
        cout << setw(15) << lavec_p[2][0];
        cout << " : a1" << endl;

        cout << setw(16) << lavec_p[0][1];
        cout << setw(15) << lavec_p[1][1];
        cout << setw(15) << lavec_p[2][1];
        cout << " : a2" << endl;

        cout << setw(16) << lavec_p[0][2];
        cout << setw(15) << lavec_p[1][2];
        cout << setw(15) << lavec_p[2][2];
        cout << " : a3" << endl;
        cout << endl;

        cout << setw(16) << rlavec_p[0][0];
        cout << setw(15) << rlavec_p[0][1];
        cout << setw(15) << rlavec_p[0][2];
        cout << " : b1" << endl;

        cout << setw(16) << rlavec_p[1][0];
        cout << setw(15) << rlavec_p[1][1];
        cout << setw(15) << rlavec_p[1][2];
        cout << " : b2" << endl;

        cout << setw(16) << rlavec_p[2][0];
        cout << setw(15) << rlavec_p[2][1];
        cout << setw(15) << rlavec_p[2][2];
        cout << " : b3" << endl;
        cout << endl << endl;


        for (i = 0; i < 3; ++i) {
            for (j = 0; j < 3; ++j) {
                vec_tmp[i][j] = lavec_p[j][i];
            }
        }
        volume_p = volume(vec_tmp[0], vec_tmp[1], vec_tmp[2]);

        cout << "  Volume of the primitive cell : "
             << volume_p << " (a.u.)^3" << endl << endl;
        cout << "  Number of atoms in the supercell     : "
             << nat_anharm << endl;
        cout << "  Number of atoms in the primitive cell: "
             << natmin << endl << endl;

        if (fcs_phonon->update_fc2) {
            cout << endl;
            cout << "  FC2XML is given: Harmonic IFCs will be replaced by the values in "
                 << fcs_phonon->file_fc2 << endl;
            cout << endl;

            cout << " * Supercell for HARMONIC (from "
                 << fcs_phonon->file_fc2 << " )" << endl << endl;

            cout << setw(16) << lavec_s[0][0];
            cout << setw(15) << lavec_s[1][0];
            cout << setw(15) << lavec_s[2][0];
            cout << " : a1" << endl;

            cout << setw(16) << lavec_s[0][1];
            cout << setw(15) << lavec_s[1][1];
            cout << setw(15) << lavec_s[2][1];
            cout << " : a2" << endl;

            cout << setw(16) << lavec_s[0][2];
            cout << setw(15) << lavec_s[1][2];
            cout << setw(15) << lavec_s[2][2];
            cout << " : a3" << endl;
            cout << endl;

            cout << setw(16) << rlavec_s[0][0];
            cout << setw(15) << rlavec_s[0][1];
            cout << setw(15) << rlavec_s[0][2];
            cout << " : b1" << endl;

            cout << setw(16) << rlavec_s[1][0];
            cout << setw(15) << rlavec_s[1][1];
            cout << setw(15) << rlavec_s[1][2];
            cout << " : b2" << endl;

            cout << setw(16) << rlavec_s[2][0];
            cout << setw(15) << rlavec_s[2][1];
            cout << setw(15) << rlavec_s[2][2];
            cout << " : b3" << endl;
            cout << endl;

            cout << "  Number of atoms in the supercell (HARMONIC)   : " << nat << endl;
            cout << endl;
        }

        memory->allocate(xtmp, natmin, 3);

        for (i = 0; i < natmin; ++i) {
            rotvec(xtmp[i], xr_s[map_p2s[i][0]], lavec_s);
            rotvec(xtmp[i], xtmp[i], rlavec_p);
            for (j = 0; j < 3; ++j) xtmp[i][j] /= 2.0 * pi;
        }

        cout << "  Atomic positions in the primitive cell (fractional):" << endl;
        for (i = 0; i < natmin; ++i) {
            cout << setw(4) << i + 1 << ":";
            for (j = 0; j < 3; ++j) {
                cout << setw(15) << xtmp[i][j];
            }
            cout << setw(4) << symbol_kd[kd[map_p2s[i][0]]] << endl;
        }
        cout << endl;

        memory->deallocate(xtmp);

        if (lspin) {
            cout << "  MagneticMoments entry found in the XML file. " << endl;
            cout << "  Magnetic moment in Cartesian coordinates: " << endl;
            for (i = 0; i < natmin; ++i) {
                cout << setw(4) << i + 1 << ":";
                for (j = 0; j < 3; ++j) {
                    cout << setw(15) << magmom[i][j];
                }
                cout << endl;
            }
            cout << endl;
            if (noncollinear == 0) {
                cout << "  Collinear calculation: magnetic moments are considered as scalar variables." << endl;
            } else if (noncollinear == 1) {
                cout << "  Noncollinear calculation: magnetic moments are considered as vector variables." << endl;
                if (symmetry->trev_sym_mag) {
                    cout << "  Time-reversal symmetry will be considered for generating magnetic space group" << endl;
                } else {
                    cout << "  Time-reversal symmetry will NOT be considered for generating magnetic space group" <<
                         endl;
                }
            }
            cout << endl;
        }

        cout << "  Mass of atomic species (u):" << endl;
        for (i = 0; i < nkd; ++i) {
            cout << setw(4) << symbol_kd[i] << ":";
            cout << fixed << setw(12) << mass_kd[i] << endl;
        }
        cout << endl << endl;
    }

    // Check the consistency of FCSXML and FC2XML
    MPI_Bcast(&fcs_phonon->update_fc2, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD);
    if (fcs_phonon->update_fc2) {
        memory->allocate(map_p2s_anharm_orig, natmin, ntran_anharm);
        check_consistency_primitive_lattice();
    }
    // Atomic masses in Rydberg unit

    memory->allocate(mass, nat);
    memory->allocate(mass_anharm, nat_anharm);
    for (i = 0; i < nat; ++i) {
        mass[i] = mass_kd[kd[i]] * amu_ry;
    }
    for (i = 0; i < nat_anharm; ++i) {
        mass_anharm[i] = mass_kd[kd_anharm[i]] * amu_ry;
    }
    MPI_Bcast(&Tmin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Tmax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&dT, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&volume_p, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    memory->allocate(kd_prim, natmin);

    for (i = 0; i < natmin; ++i) {
        kd_prim[i] = kd[map_p2s[i][0]];
    }
    MPI_Bcast(&lspin, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD);
    if (mympi->my_rank > 0) {
        memory->allocate(magmom, natmin, 3);
    }
    MPI_Bcast(&magmom[0][0], 3 * natmin, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&noncollinear, 1, MPI_INT, 0, MPI_COMM_WORLD);

    setup_atomic_class(natmin, kd_prim, magmom);

    memory->deallocate(kd_prim);
}

void System::load_system_info_from_XML()
{
    if (mympi->my_rank == 0) {

        int i;
        using namespace boost::property_tree;
        ptree pt;

        std::map<std::string, int> dict_atomic_kind;

        try {
            read_xml(fcs_phonon->file_fcs, pt);
        }
        catch (std::exception &e) {
            std::string str_error = "Cannot open file FCSXML ( "
                                    + fcs_phonon->file_fcs + " )";
            error->exit("load_system_info_from_XML",
                        str_error.c_str());
        }

        // Parse nat and ntran

        nat = boost::lexical_cast<unsigned int>(
                get_value_from_xml(pt,
                                   "Data.Structure.NumberOfAtoms"));
        int nkd_tmp = boost::lexical_cast<unsigned int>(
                get_value_from_xml(pt,
                                   "Data.Structure.NumberOfElements"));

        if (nkd != nkd_tmp)
            error->exit("load_system_info_from_XML",
                        "NKD in the FCSXML file is not consistent with that given in the input file.");

        ntran = boost::lexical_cast<unsigned int>(
                get_value_from_xml(pt,
                                   "Data.Symmetry.NumberOfTranslations"));

        natmin = nat / ntran;

        // Parse lattice vectors

        std::stringstream ss;

        for (i = 0; i < 3; ++i) {
            ss.str("");
            ss.clear();
            ss << get_value_from_xml(pt,
                                     "Data.Structure.LatticeVector.a"
                                     + std::to_string(i + 1));
            ss >> lavec_s[0][i] >> lavec_s[1][i] >> lavec_s[2][i];
        }

        // Parse atomic elements and coordinates

        memory->allocate(xr_s, nat, 3);
        memory->allocate(kd, nat);

        BOOST_FOREACH (const ptree::value_type &child_, pt.get_child("Data.Structure.AtomicElements")) {
                        const auto &child = child_.second;
                        const auto icount_kd = child.get<unsigned int>("<xmlattr>.number");
                        dict_atomic_kind[boost::lexical_cast<std::string>(child_.second.data())] = icount_kd - 1;
                    }

        unsigned int index;

        BOOST_FOREACH (const ptree::value_type &child_, pt.get_child("Data.Structure.Position")) {
                        const auto &child = child_.second;
                        const auto str_index = child.get<std::string>("<xmlattr>.index");
                        const auto str_element = child.get<std::string>("<xmlattr>.element");

                        ss.str("");
                        ss.clear();
                        ss << child.data();

                        index = boost::lexical_cast<unsigned int>(str_index) - 1;

                        if (index >= nat)
                            error->exit("load_system_info_xml",
                                        "index is out of range");

                        kd[index] = dict_atomic_kind[str_element];
                        ss >> xr_s[index][0] >> xr_s[index][1] >> xr_s[index][2];
                    }

        dict_atomic_kind.clear();

        // Parse mapping information

        memory->allocate(map_p2s, natmin, ntran);
        memory->allocate(map_s2p, nat);

        BOOST_FOREACH (const ptree::value_type &child_, pt.get_child("Data.Symmetry.Translations")) {
                        const auto &child = child_.second;
                        const auto str_tran = child.get<std::string>("<xmlattr>.tran");
                        const auto str_atom = child.get<std::string>("<xmlattr>.atom");

                        const auto tran = boost::lexical_cast<unsigned int>(str_tran) - 1;
                        const auto atom_p = boost::lexical_cast<unsigned int>(str_atom) - 1;
                        const auto atom_s = boost::lexical_cast<unsigned int>(child.data()) - 1;

                        if (tran >= ntran || atom_p >= natmin || atom_s >= nat) {
                            error->exit("load_system_info_xml",
                                        "index is out of range");
                        }

                        map_p2s[atom_p][tran] = atom_s;
                        map_s2p[atom_s].atom_num = atom_p;
                        map_s2p[atom_s].tran_num = tran;
                    }

        // Parse magnetic moments

        double **magmom_tmp;
        memory->allocate(magmom_tmp, nat, 3);
        memory->allocate(magmom, natmin, 3);

        lspin = true;
        try {
            BOOST_FOREACH(const ptree::value_type &child_, pt.get_child("Data.MagneticMoments")) {
                            if (child_.first == "mag") {
                                const auto &child = child_.second;
                                const auto str_index = child.get<std::string>("<xmlattr>.index");

                                ss.str("");
                                ss.clear();
                                ss << child.data();

                                index = boost::lexical_cast<unsigned int>(str_index) - 1;

                                if (index >= nat)
                                    error->exit("load_system_info_xml",
                                                "index is out of range");

                                ss >> magmom_tmp[index][0]
                                   >> magmom_tmp[index][1]
                                   >> magmom_tmp[index][2];
                            }
                        }

        }
        catch (...) {
            lspin = false;
        }

        if (lspin) {
            for (i = 0; i < natmin; ++i) {
                for (int j = 0; j < 3; ++j) {
                    magmom[i][j] = magmom_tmp[map_p2s[i][0]][j];
                }
            }

            try {
                noncollinear = boost::lexical_cast<int>(
                        get_value_from_xml(pt,
                                           "Data.MagneticMoments.Noncollinear"));
            }
            catch (...) {
                noncollinear = 0;
            }

            try {
                symmetry->trev_sym_mag = boost::lexical_cast<int>(
                        get_value_from_xml(pt,
                                           "Data.MagneticMoments.TimeReversalSymmetry"));
            }
            catch (...) {
                symmetry->trev_sym_mag = true;
            }
        } else {
            for (i = 0; i < natmin; ++i) {
                for (int j = 0; j < 3; ++j) {
                    magmom[i][j] = 0.0;
                }
            }
            noncollinear = 0;
            symmetry->trev_sym_mag = true;
        }
        memory->deallocate(magmom_tmp);

        // Now, replicate the information for anharmonic terms.

        int j;
        nat_anharm = nat;
        ntran_anharm = ntran;
        memory->allocate(xr_s_anharm, nat_anharm, 3);
        memory->allocate(kd_anharm, nat_anharm);
        memory->allocate(map_p2s_anharm, natmin, ntran_anharm);
        memory->allocate(map_s2p_anharm, nat_anharm);

        for (i = 0; i < 3; ++i) {
            for (j = 0; j < 3; ++j) lavec_s_anharm[i][j] = lavec_s[i][j];
        }
        for (i = 0; i < nat_anharm; ++i) {
            for (j = 0; j < 3; ++j) xr_s_anharm[i][j] = xr_s[i][j];
            kd_anharm[i] = kd[i];
            map_s2p_anharm[i] = map_s2p[i];
        }
        for (i = 0; i < natmin; ++i) {
            for (j = 0; j < ntran_anharm; ++j) {
                map_p2s_anharm[i][j] = map_p2s[i][j];
            }
        }

        if (fcs_phonon->update_fc2) {

            // When FC2XML is given, structural information is updated only for harmonic terms.

            try {
                read_xml(fcs_phonon->file_fc2, pt);
            }
            catch (std::exception &e) {
                auto str_error = "Cannot open file FC2XML ( "
                                 + fcs_phonon->file_fc2 + " )";
                error->exit("load_system_info_from_XML",
                            str_error.c_str());
            }

            // Parse nat and ntran

            nat = boost::lexical_cast<unsigned int>(
                    get_value_from_xml(pt,
                                       "Data.Structure.NumberOfAtoms"));
            nkd_tmp = boost::lexical_cast<unsigned int>(
                    get_value_from_xml(pt,
                                       "Data.Structure.NumberOfElements"));

            if (nkd != nkd_tmp)
                error->exit("load_system_info_from_XML",
                            "NKD in the FC2XML file is not consistent with that given in the input file.");

            ntran = boost::lexical_cast<unsigned int>(
                    get_value_from_xml(pt,
                                       "Data.Symmetry.NumberOfTranslations"));

            const int natmin_tmp = nat / ntran;

            if (natmin_tmp != natmin)
                error->exit("load_system_info_from_XML",
                            "Number of atoms in a primitive cell is different in FCSXML and FC2XML.");

            memory->deallocate(xr_s);
            memory->deallocate(kd);
            memory->deallocate(map_p2s);
            memory->deallocate(map_s2p);


            // Parse lattice vectors

            std::stringstream ss;

            for (i = 0; i < 3; ++i) {
                ss.str("");
                ss.clear();
                ss << get_value_from_xml(pt,
                                         "Data.Structure.LatticeVector.a"
                                         + std::to_string(i + 1));
                ss >> lavec_s[0][i] >> lavec_s[1][i] >> lavec_s[2][i];
            }

            // Parse atomic elements and coordinates

            memory->allocate(xr_s, nat, 3);
            memory->allocate(kd, nat);

            BOOST_FOREACH (const ptree::value_type &child_, pt.get_child("Data.Structure.AtomicElements")) {
                            const auto &child = child_.second;
                            const auto icount_kd = child.get<unsigned int>("<xmlattr>.number");
                            dict_atomic_kind[boost::lexical_cast<std::string>(child_.second.data())] = icount_kd - 1;
                        }

            BOOST_FOREACH (const ptree::value_type &child_, pt.get_child("Data.Structure.Position")) {
                            const auto &child = child_.second;
                            const auto str_index = child.get<std::string>("<xmlattr>.index");
                            const auto str_element = child.get<std::string>("<xmlattr>.element");

                            ss.str("");
                            ss.clear();
                            ss << child.data();

                            auto index_kd = boost::lexical_cast<unsigned int>(str_index) - 1;

                            if (index_kd >= nat)
                                error->exit("load_system_info_xml",
                                            "index is out of range");

                            kd[index_kd] = dict_atomic_kind[str_element];
                            ss >> xr_s[index_kd][0] >> xr_s[index_kd][1] >> xr_s[index_kd][2];
                        }

            dict_atomic_kind.clear();

            // Parse mapping information

            memory->allocate(map_p2s, natmin, ntran);
            memory->allocate(map_s2p, nat);

            BOOST_FOREACH (const ptree::value_type &child_, pt.get_child("Data.Symmetry.Translations")) {
                            const auto &child = child_.second;
                            const auto str_tran = child.get<std::string>("<xmlattr>.tran");
                            const auto str_atom = child.get<std::string>("<xmlattr>.atom");

                            const auto tran = boost::lexical_cast<unsigned int>(str_tran) - 1;
                            const auto atom_p = boost::lexical_cast<unsigned int>(str_atom) - 1;
                            const auto atom_s = boost::lexical_cast<unsigned int>(child.data()) - 1;

                            if (tran >= ntran || atom_p >= natmin || atom_s >= nat) {
                                error->exit("load_system_info_xml", "index is out of range");
                            }

                            map_p2s[atom_p][tran] = atom_s;
                            map_s2p[atom_s].atom_num = atom_p;
                            map_s2p[atom_s].tran_num = tran;
                        }
        }

    }

    MPI_Bcast(&lavec_s[0][0], 9, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&lavec_p[0][0], 9, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&lavec_s_anharm[0][0], 9, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Bcast(&nkd, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nat, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nat_anharm, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(&natmin, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ntran, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ntran_anharm, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(&lspin, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD);

    if (mympi->my_rank > 0) {
        memory->allocate(mass_kd, nkd);
        memory->allocate(xr_s, nat, 3);
        memory->allocate(xr_s_anharm, nat_anharm, 3);
        memory->allocate(kd, nat);
        memory->allocate(kd_anharm, nat_anharm);
        memory->allocate(map_p2s, natmin, ntran);
        memory->allocate(map_p2s_anharm, natmin, ntran_anharm);
        memory->allocate(map_s2p, nat);
        memory->allocate(map_s2p_anharm, nat_anharm);
        if (lspin) {
            memory->allocate(magmom, natmin, 3);
        }
    }

    MPI_Bcast(&mass_kd[0], nkd, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&xr_s[0][0], 3 * nat, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&xr_s_anharm[0][0], 3 * nat_anharm, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&kd[0], nat, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(&kd_anharm[0], nat_anharm, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(&map_p2s[0][0], natmin * ntran, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(&map_p2s_anharm[0][0], natmin * ntran_anharm, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(&map_s2p[0], nat * sizeof map_s2p[0], MPI_BYTE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&map_s2p_anharm[0], nat_anharm * sizeof map_s2p_anharm[0], MPI_BYTE, 0, MPI_COMM_WORLD);
    if (lspin) MPI_Bcast(&magmom[0][0], 3 * natmin, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}


void System::recips(double vec[3][3],
                    double inverse[3][3]) const
{
    const auto det = vec[0][0] * vec[1][1] * vec[2][2]
                     + vec[1][0] * vec[2][1] * vec[0][2]
                     + vec[2][0] * vec[0][1] * vec[1][2]
                     - vec[0][0] * vec[2][1] * vec[1][2]
                     - vec[2][0] * vec[1][1] * vec[0][2]
                     - vec[1][0] * vec[0][1] * vec[2][2];

    if (std::abs(det) < eps12) {
        error->exit("recips", "Lattice Vector is singular");
    }

    const auto factor = 2.0 * pi / det;

    inverse[0][0] = (vec[1][1] * vec[2][2] - vec[1][2] * vec[2][1]) * factor;
    inverse[0][1] = (vec[0][2] * vec[2][1] - vec[0][1] * vec[2][2]) * factor;
    inverse[0][2] = (vec[0][1] * vec[1][2] - vec[0][2] * vec[1][1]) * factor;

    inverse[1][0] = (vec[1][2] * vec[2][0] - vec[1][0] * vec[2][2]) * factor;
    inverse[1][1] = (vec[0][0] * vec[2][2] - vec[0][2] * vec[2][0]) * factor;
    inverse[1][2] = (vec[0][2] * vec[1][0] - vec[0][0] * vec[1][2]) * factor;

    inverse[2][0] = (vec[1][0] * vec[2][1] - vec[1][1] * vec[2][0]) * factor;
    inverse[2][1] = (vec[0][1] * vec[2][0] - vec[0][0] * vec[2][1]) * factor;
    inverse[2][2] = (vec[0][0] * vec[1][1] - vec[0][1] * vec[1][0]) * factor;
}

double System::volume(const double vec1[3],
                      const double vec2[3],
                      const double vec3[3]) const
{
    const auto vol = std::abs(vec1[0] * (vec2[1] * vec3[2] - vec2[2] * vec3[1])
                              + vec1[1] * (vec2[2] * vec3[0] - vec2[0] * vec3[2])
                              + vec1[2] * (vec2[0] * vec3[1] - vec2[1] * vec3[0]));

    return vol;
}


void System::setup_atomic_class(const unsigned int N,
                                const unsigned int *kd,
                                double **magmom_in)
{
    // In the case of collinear calculation, spin moments are considered as scalar
    // variables. Therefore, the same elements with different magnetic moments are
    // considered as different types. In noncollinear calculations, 
    // magnetic moments are not considered in this stage. They will be treated
    // separately in symmetry.cpp where spin moments will be rotated and flipped 
    // using time-reversal symmetry.

    unsigned int i;
    AtomType type_tmp;
    std::set<AtomType> set_type;
    set_type.clear();

    for (i = 0; i < N; ++i) {
        type_tmp.element = kd[i];

        if (noncollinear == 0) {
            type_tmp.magmom = magmom_in[i][2];
        } else {
            type_tmp.magmom = 0.0;
        }
        set_type.insert(type_tmp);
    }

    nclassatom = set_type.size();

    memory->allocate(atomlist_class, nclassatom);

    for (i = 0; i < N; ++i) {
        int count = 0;
        for (auto it : set_type) {
            if (noncollinear) {
                if (kd[i] == it.element) {
                    atomlist_class[count].push_back(i);
                }
            } else {
                if (kd[i] == it.element &&
                    std::abs(magmom[i][2] - it.magmom) < eps6) {
                    atomlist_class[count].push_back(i);
                }
            }
            ++count;
        }
    }
    set_type.clear();
}

void System::check_consistency_primitive_lattice() const
{
    // Check if the ordering of atoms in the primitive cells derived 
    // from FCSXML and FC2XML are same or not. If not, the ordering for the
    // FCSXML (anharmonic terms) will be changed so that it becomes equivalent
    // to that of FC2XML. 
    // This operation is necessary for obtaining correct computational results.

    int i, j, k;
    double xdiff[3];
    double **x_harm, **x_anharm;

    memory->allocate(x_harm, natmin, 3);
    memory->allocate(x_anharm, natmin, 3);

    std::vector<int> map_anh2harm;
    map_anh2harm.resize(natmin);

    for (i = 0; i < natmin; ++i) {
        rotvec(x_harm[i], xr_s[map_p2s[i][0]], lavec_s);
        rotvec(x_harm[i], x_harm[i], rlavec_p);
        for (j = 0; j < 3; ++j) x_harm[i][j] /= 2.0 * pi;
    }

    for (i = 0; i < natmin; ++i) {
        rotvec(x_anharm[i], xr_s_anharm[map_p2s_anharm[i][0]], lavec_s_anharm);
        rotvec(x_anharm[i], x_anharm[i], rlavec_p);
        for (j = 0; j < 3; ++j) x_anharm[i][j] /= 2.0 * pi;
    }

    for (i = 0; i < natmin; ++i) {

        int iloc = -1;

        for (j = 0; j < natmin; ++j) {

            for (k = 0; k < 3; ++k) {
                xdiff[k] = x_anharm[i][k] - x_harm[j][k];
                xdiff[k] = xdiff[k] - static_cast<double>(nint(xdiff[k]));
            }

            const auto norm = xdiff[0] * xdiff[0] + xdiff[1] * xdiff[1] + xdiff[2] * xdiff[2];
            if (norm < eps4 && kd[map_p2s[j][0]] == kd_anharm[map_p2s_anharm[i][0]]) {
                iloc = j;
                break;
            }
        }

        if (iloc == -1) {
            error->exit("check_consistency_primitive",
                        "Could not find equivalent atom. Probably, the crystal structure is different.");
        }

        map_anh2harm[i] = iloc;
    }

    memory->deallocate(x_harm);
    memory->deallocate(x_anharm);

    // Rebuild the mapping information for anharmonic terms.

    unsigned int **map_p2s_tmp;

    memory->allocate(map_p2s_tmp, natmin, ntran_anharm);

    for (i = 0; i < natmin; ++i) {
        for (j = 0; j < ntran_anharm; ++j) {
            map_p2s_anharm_orig[i][j] = map_p2s_anharm[i][j];
        }
    }

    for (i = 0; i < ntran_anharm; ++i) {
        for (j = 0; j < natmin; ++j) {
            map_p2s_tmp[j][i] = map_p2s_anharm[j][i];
        }
    }

    for (i = 0; i < ntran_anharm; ++i) {
        for (j = 0; j < natmin; ++j) {
            map_p2s_anharm[map_anh2harm[j]][i] = map_p2s_tmp[j][i];
        }
    }

    for (i = 0; i < ntran_anharm; ++i) {
        for (j = 0; j < natmin; ++j) {
            k = map_p2s_anharm[j][i];
            map_s2p_anharm[k].atom_num = j;
            map_s2p_anharm[k].tran_num = i;
        }
    }

    memory->deallocate(map_p2s_tmp);
    map_anh2harm.clear();
}

int System::get_atomic_number_by_name(const std::string &kdname_in)
{
    auto kdname_copy = kdname_in;
    kdname_copy[0] = toupper(kdname_copy[0]);

    int ret = -1;

    for (auto i = 0; i < element_names.size(); ++i) {
        if (kdname_copy == element_names[i]) {
            ret = i;
            break;
        }
    }
    return ret;
}

void System::set_mass_elem_from_database(const int nkd,
                                         const std::string *symbol_in,
                                         double *mass_kd_out)
{
    for (int i = 0; i < nkd; ++i) {
        const auto atom_number = get_atomic_number_by_name(symbol_in[i]);
        if (atom_number >= element_names.size() || atom_number == -1) {
            error->exit("set_mass_elem_from_database",
                        "Atomic mass for the given element doesn't exist in the database.\nTherefore, please input MASS manually.");
        }
        const auto mass_tmp = atomic_masses[atom_number];
        if (mass_tmp < 0.0) {
            error->exit("set_mass_elem_from_database",
                        "One of the elements in the KD-tag is unstable. \nTherefore, please input MASS manually.");
        }
        mass_kd_out[i] = mass_tmp;
    }
}

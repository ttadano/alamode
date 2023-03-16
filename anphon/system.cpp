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
#include <Eigen/LU>
#include <Eigen/Geometry>

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
//    xr_p = nullptr;
//    xr_s = nullptr;
//    xc = nullptr;
//    xr_s_anharm = nullptr;
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
//    if (xr_p) {
//        deallocate(xr_p);
//    }
//    if (xr_s) {
//        deallocate(xr_s);
//    }
//    if (xc) {
//        deallocate(xc);
//    }
//    if (xr_s_anharm) {
//        deallocate(xr_s_anharm);
//    }
    if (kd) {
        deallocate(kd);
    }
    if (kd_anharm) {
        deallocate(kd_anharm);
    }
    if (mass_kd) {
        deallocate(mass_kd);
    }
    if (mass) {
        deallocate(mass);
    }
    if (mass_anharm) {
        deallocate(mass_anharm);
    }
    if (symbol_kd) {
        deallocate(symbol_kd);
    }
    if (map_p2s) {
        deallocate(map_p2s);
    }
    if (map_p2s_anharm) {
        deallocate(map_p2s_anharm);
    }
    if (map_s2p) {
        deallocate(map_s2p);
    }
    if (map_s2p_anharm) {
        deallocate(map_s2p_anharm);
    }
    if (magmom) {
        deallocate(magmom);
    }
    if (map_p2s_anharm_orig) {
        deallocate(map_p2s_anharm_orig);
    }
    if (atomlist_class) {
        deallocate(atomlist_class);
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
            allocate(mass_kd, nkd);
            set_mass_elem_from_database(nkd, symbol_kd, mass_kd);
        }
    }
    load_system_info_from_file();
    load_system_info_from_XML();

    recips(lavec_s, rlavec_s);
    recips(lavec_s_anharm, rlavec_s_anharm);
    recips(lavec_p, rlavec_p);

    xr_p.resize(nat, 3);
    xc.resize(nat, 3);

    xc = xr_s * lavec_s.transpose();
    xr_p = xc * lavec_p.inverse().transpose();

//    for (i = 0; i < nat; ++i) {
//        rotvec(xc[i], xr_s[i], lavec_s);
//        rotvec(xr_p[i], xc[i], rlavec_p);
//        for (j = 0; j < 3; ++j) {
//            xr_p[i][j] /= 2.0 * pi;
//        }
//    }

    if (mympi->my_rank == 0) {
        cout << " -----------------------------------------------------------------" << endl;
        cout << endl;
        cout << " ===================\n";
        cout << "  Crystal Structure \n";
        cout << " ===================\n\n";
        cout << " Lattice Vectors:\n\n";
        cout.setf(ios::scientific);

        cout << " * Supercell (from " << fcs_phonon->file_fcs << " )" << endl << endl;
        cout << setw(16) << supercell_base.lattice_vector(0, 0);
        cout << setw(15) << supercell_base.lattice_vector(1, 0);
        cout << setw(15) << supercell_base.lattice_vector(2, 0);
        cout << " : a1" << endl;

        cout << setw(16) << supercell_base.lattice_vector(0, 1);
        cout << setw(15) << supercell_base.lattice_vector(1, 1);
        cout << setw(15) << supercell_base.lattice_vector(2, 1);
        cout << " : a2" << endl;

        cout << setw(16) << supercell_base.lattice_vector(0, 2);
        cout << setw(15) << supercell_base.lattice_vector(1, 2);
        cout << setw(15) << supercell_base.lattice_vector(2, 2);
        cout << " : a3" << endl;
        cout << endl;

//        cout << setw(16) << rlavec_s_anharm[0][0];
//        cout << setw(15) << rlavec_s_anharm[0][1];
//        cout << setw(15) << rlavec_s_anharm[0][2];
//        cout << " : b1" << endl;
//
//        cout << setw(16) << rlavec_s_anharm[1][0];
//        cout << setw(15) << rlavec_s_anharm[1][1];
//        cout << setw(15) << rlavec_s_anharm[1][2];
//        cout << " : b2" << endl;
//
//        cout << setw(16) << rlavec_s_anharm[2][0];
//        cout << setw(15) << rlavec_s_anharm[2][1];
//        cout << setw(15) << rlavec_s_anharm[2][2];
//        cout << " : b3" << endl;
//        cout << endl;

        cout << " * Primitive cell " << endl << endl;
        cout << setw(16) << lavec_p(0, 0);
        cout << setw(15) << lavec_p(1, 0);
        cout << setw(15) << lavec_p(2, 0);
        cout << " : a1" << endl;

        cout << setw(16) << lavec_p(0, 1);
        cout << setw(15) << lavec_p(1, 1);
        cout << setw(15) << lavec_p(2, 1);
        cout << " : a2" << endl;

        cout << setw(16) << lavec_p(0, 2);
        cout << setw(15) << lavec_p(1, 2);
        cout << setw(15) << lavec_p(2, 2);
        cout << " : a3" << endl;
        cout << endl;

        cout << setw(16) << rlavec_p(0, 0);
        cout << setw(15) << rlavec_p(0, 1);
        cout << setw(15) << rlavec_p(0, 2);
        cout << " : b1" << endl;

        cout << setw(16) << rlavec_p(1, 0);
        cout << setw(15) << rlavec_p(1, 1);
        cout << setw(15) << rlavec_p(1, 2);
        cout << " : b2" << endl;

        cout << setw(16) << rlavec_p(2, 0);
        cout << setw(15) << rlavec_p(2, 1);
        cout << setw(15) << rlavec_p(2, 2);
        cout << " : b3" << endl;
        cout << endl << endl;

        for (i = 0; i < 3; ++i) {
            for (j = 0; j < 3; ++j) {
                vec_tmp[i][j] = lavec_p(j, i);
            }
        }
        volume_p = volume(vec_tmp[0], vec_tmp[1], vec_tmp[2]);

        cout << "  Volume of the primitive cell : "
             << volume_p << " (a.u.)^3" << endl << endl;
        cout << "  Number of atoms in the supercell     : "
             << supercell_base.number_of_atoms << endl;
        cout << "  Number of atoms in the primitive cell: "
             << natmin << endl << endl;

        if (fcs_phonon->update_fc2) {
            cout << endl;
            cout << "  FC2XML is given: Harmonic IFCs will be replaced by the values in "
                 << fcs_phonon->file_fc2 << endl;
            cout << endl;

            cout << " * Supercell for HARMONIC (from "
                 << fcs_phonon->file_fc2 << " )" << endl << endl;

            cout << setw(16) << supercell_fc2.lattice_vector(0, 0);
            cout << setw(15) << supercell_fc2.lattice_vector(1, 0);
            cout << setw(15) << supercell_fc2.lattice_vector(2, 0);
            cout << " : a1" << endl;

            cout << setw(16) << supercell_fc2.lattice_vector(0, 1);
            cout << setw(15) << supercell_fc2.lattice_vector(1, 1);
            cout << setw(15) << supercell_fc2.lattice_vector(2, 1);
            cout << " : a2" << endl;

            cout << setw(16) << supercell_fc2.lattice_vector(0, 2);
            cout << setw(15) << supercell_fc2.lattice_vector(1, 2);
            cout << setw(15) << supercell_fc2.lattice_vector(2, 2);
            cout << " : a3" << endl;
            cout << endl;

//            cout << setw(16) << rlavec_s[0][0];
//            cout << setw(15) << rlavec_s[0][1];
//            cout << setw(15) << rlavec_s[0][2];
//            cout << " : b1" << endl;
//
//            cout << setw(16) << rlavec_s[1][0];
//            cout << setw(15) << rlavec_s[1][1];
//            cout << setw(15) << rlavec_s[1][2];
//            cout << " : b2" << endl;
//
//            cout << setw(16) << rlavec_s[2][0];
//            cout << setw(15) << rlavec_s[2][1];
//            cout << setw(15) << rlavec_s[2][2];
//            cout << " : b3" << endl;
//            cout << endl;

            cout << "  Number of atoms in the supercell (HARMONIC)   : " << supercell_fc2.number_of_atoms << endl;
            cout << endl;
        }

        Eigen::MatrixXd xtmp(natmin, 3);

//        allocate(xtmp, natmin, 3);

        for (i = 0; i < natmin; ++i) {
            for (j = 0; j < 3; ++j) {
                xtmp(i, j) = xr_s(map_p2s[i][0], j);
            }
        }
        xtmp = xtmp * lavec_s.transpose();
        xtmp = xtmp * lavec_p.transpose().inverse();

//        for (i = 0; i < natmin; ++i) {
//            rotvec(xtmp[i], xr_s[map_p2s[i][0]], lavec_s);
//            rotvec(xtmp[i], xtmp[i], rlavec_p);
//            for (j = 0; j < 3; ++j) xtmp[i][j] /= 2.0 * pi;
//        }

        cout << "  Atomic positions in the primitive cell (fractional):" << endl;
        for (i = 0; i < natmin; ++i) {
            cout << setw(4) << i + 1 << ":";
            for (j = 0; j < 3; ++j) {
                cout << setw(15) << xtmp(i, j);
            }
            cout << setw(4) << symbol_kd[kd[map_p2s[i][0]]] << endl;
        }
        cout << endl;

//        deallocate(xtmp);

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
                if (symmetry->time_reversal_sym_from_alm) {
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
    MPI_Bcast(&fcs_phonon->update_fc2, 1, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);
    if (fcs_phonon->update_fc2) {
        allocate(map_p2s_anharm_orig, natmin, ntran_anharm);
        check_consistency_primitive_lattice();
    }
    // Atomic masses in Rydberg unit

    allocate(mass, nat);
    allocate(mass_anharm, nat_anharm);
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

    allocate(kd_prim, natmin);

    for (i = 0; i < natmin; ++i) {
        kd_prim[i] = kd[map_p2s[i][0]];
    }
    MPI_Bcast(&lspin, 1, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);
    if (mympi->my_rank > 0) {
        allocate(magmom, natmin, 3);
    }
    MPI_Bcast(&magmom[0][0], 3 * natmin, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&noncollinear, 1, MPI_INT, 0, MPI_COMM_WORLD);

    setup_atomic_class(natmin, kd_prim, magmom);

    deallocate(kd_prim);
}


void System::load_system_info_from_file()
{
    // Parse structure information either from the XML file or h5 file
    int filetype[4]; // filetype[0] for FCSFILE, filetype[X-1] for FCXFILE (X=2, 3, 4)
    // -1: FC?FILE not given, 0: FC?FILE in XML format, 1: FC?FLIE in HDF5 format
    std::vector<std::string> filename_list{fcs_phonon->file_fcs, fcs_phonon->file_fc2,
                                           fcs_phonon->file_fc3, fcs_phonon->file_fc4};

    if (mympi->my_rank == 0) {
        for (auto i = 0; i < filename_list.size(); ++i) {
            const auto filename = filename_list[i];
            if (!filename.empty()) {
                const auto file_extension = filename.substr(filename.find_last_of('.') + 1);
                if (file_extension == "xml" || file_extension == "XML") {
                    filetype[i] = 0;
                } else if (file_extension == "h5" || file_extension == "hdf5") {
                    filetype[i] = 1;
                }
            } else {
                filetype[i] = -1;
            }
        }
    } else {
        filename_list.resize(4);
    }
    MPI_Bcast(&filetype[0], 4, MPI_INT, 0, MPI_COMM_WORLD);
    for (auto i = 0; i < 4; ++i) {
        mympi->MPI_Bcast_string(filename_list[i], 0, MPI_COMM_WORLD);
    }

    Spin spin_fc2;

    for (auto i = 0; i < 4; ++i) {

        Cell scell, pcell;
        Spin spin_s, spin_p;
        MappingTable map_s, map_p;

        if (filetype[i] == -1) continue;

        if (filetype[i] == 0) {
            get_structure_and_mapping_table_xml(filename_list[i], scell, pcell,
                                                spin_s, spin_p,
                                                map_s, map_p);
        } else if (filetype[i] == 1) {

        } else {
            exit("load_system_info_from_file", "This cannot happen.");
        }

        if (i == 0) {
            supercell_base = scell;
            primcell_base = pcell;
            map_scell_base = map_s;
            map_pcell_base = map_p;
            spin_base = spin_s;

        } else if (i == 1) {
            supercell_fc2 = scell;
            primcell_fc2 = pcell;
            map_scell_fc2 = map_s;
            map_pcell_fc2 = map_p;
            spin_fc2 = spin_s;

        } else if (i == 2) {
            supercell_fc3 = scell;
            primcell_fc3 = pcell;
            map_scell_fc3 = map_s;
            map_pcell_fc3 = map_p;

        } else if (i == 3) {
            supercell_fc4 = scell;
            primcell_fc4 = pcell;
            map_scell_fc4 = map_s;
            map_pcell_fc4 = map_p;
        }
    }

    if (filetype[0] == -1) {
        supercell_base = supercell_fc2;
        primcell_base = primcell_fc2;
        map_scell_base = map_scell_fc2;
        map_pcell_base = map_pcell_fc2;
        spin_base = spin_fc2;
    }

}

void System::get_structure_and_mapping_table_xml(const std::string &filename,
                                                 Cell &scell_out,
                                                 Cell &pcell_out,
                                                 Spin &spin_super_out,
                                                 Spin &spin_prim_out,
                                                 MappingTable &map_super_out,
                                                 MappingTable &map_prim_out) const
{
    unsigned int nat_tmp, nkd_tmp, ntran_tmp, natmin_tmp;
    double lavec_s_tmp[3][3];
    double **xr_s_tmp = nullptr;
    int *kd_tmp = nullptr;
    unsigned int **map_p2s_tmp = nullptr;
    Maps *map_s2p_tmp = nullptr;
    double **magmom_tmp = nullptr;
    int lspin_tmp, noncollinear_tmp, time_reversal_symmetry_tmp;

    if (mympi->my_rank == 0) {

        using namespace boost::property_tree;
        ptree pt;

        std::map<std::string, int> dict_atomic_kind;

        try {
            read_xml(filename, pt);
        }
        catch (std::exception &e) {
            std::string str_error = "Cannot open file FCSFILE ( "
                                    + filename + " )";
            exit("load_system_info_from_XML",
                 str_error.c_str());
        }

        // Parse nat_base and ntran_super
        nat_tmp = boost::lexical_cast<unsigned int>(
                get_value_from_xml(pt,
                                   "Data.Structure.NumberOfAtoms"));
        nkd_tmp = boost::lexical_cast<unsigned int>(
                get_value_from_xml(pt,
                                   "Data.Structure.NumberOfElements"));
        ntran_tmp = boost::lexical_cast<unsigned int>(
                get_value_from_xml(pt,
                                   "Data.Symmetry.NumberOfTranslations"));

        natmin_tmp = nat_tmp / ntran_tmp;

        // Parse lattice vectors
        std::stringstream ss;

        for (auto i = 0; i < 3; ++i) {
            ss.str("");
            ss.clear();
            ss << get_value_from_xml(pt,
                                     "Data.Structure.LatticeVector.a"
                                     + std::to_string(i + 1));
            ss >> lavec_s_tmp[0][i] >> lavec_s_tmp[1][i] >> lavec_s_tmp[2][i];
        }

        // Parse atomic elements and coordinates

        allocate(xr_s_tmp, nat_tmp, 3);
        allocate(kd_tmp, nat_tmp);

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

                        if (index >= nat_tmp)
                            exit("load_system_info_xml",
                                 "index is out of range");

                        kd_tmp[index] = dict_atomic_kind[str_element];
                        ss >> xr_s_tmp[index][0] >> xr_s_tmp[index][1] >> xr_s_tmp[index][2];
                    }

        dict_atomic_kind.clear();

        // Parse mapping information

        allocate(map_p2s_tmp, natmin_tmp, ntran_tmp);
        allocate(map_s2p_tmp, nat_tmp);

        BOOST_FOREACH (const ptree::value_type &child_, pt.get_child("Data.Symmetry.Translations")) {
                        const auto &child = child_.second;
                        const auto str_tran = child.get<std::string>("<xmlattr>.tran");
                        const auto str_atom = child.get<std::string>("<xmlattr>.atom");

                        const auto tran = boost::lexical_cast<unsigned int>(str_tran) - 1;
                        const auto atom_p = boost::lexical_cast<unsigned int>(str_atom) - 1;
                        const auto atom_s = boost::lexical_cast<unsigned int>(child.data()) - 1;

                        if (tran >= ntran_tmp || atom_p >= natmin_tmp || atom_s >= nat_tmp) {
                            exit("load_system_info_xml",
                                 "index is out of range");
                        }

                        map_p2s_tmp[atom_p][tran] = atom_s;
                        map_s2p_tmp[atom_s].atom_num = atom_p;
                        map_s2p_tmp[atom_s].tran_num = tran;
                    }

        // Parse magnetic moments

        double **magmom_tmp2;
        allocate(magmom_tmp2, nat_tmp, 3);
        allocate(magmom_tmp, natmin_tmp, 3);

        lspin_tmp = true;
        try {
            BOOST_FOREACH(const ptree::value_type &child_, pt.get_child("Data.MagneticMoments")) {
                            if (child_.first == "mag") {
                                const auto &child = child_.second;
                                const auto str_index = child.get<std::string>("<xmlattr>.index");

                                ss.str("");
                                ss.clear();
                                ss << child.data();

                                index = boost::lexical_cast<unsigned int>(str_index) - 1;

                                if (index >= nat_tmp)
                                    exit("load_system_info_xml",
                                         "index is out of range");

                                ss >> magmom_tmp2[index][0]
                                   >> magmom_tmp2[index][1]
                                   >> magmom_tmp2[index][2];
                            }
                        }

        }
        catch (...) {
            lspin_tmp = false;
        }

        if (lspin_tmp) {
            for (auto i = 0; i < natmin_tmp; ++i) {
                for (int j = 0; j < 3; ++j) {
                    magmom_tmp[i][j] = magmom_tmp2[map_p2s_tmp[i][0]][j];
                }
            }

            try {
                noncollinear_tmp = boost::lexical_cast<int>(
                        get_value_from_xml(pt,
                                           "Data.MagneticMoments.Noncollinear"));
            }
            catch (...) {
                noncollinear_tmp = 0;
            }

            try {

                time_reversal_symmetry_tmp = boost::lexical_cast<int>(
                        get_value_from_xml(pt,
                                           "Data.MagneticMoments.TimeReversalSymmetry"));
            }
            catch (...) {
                time_reversal_symmetry_tmp = 1;
            }
        } else {
            for (auto i = 0; i < natmin_tmp; ++i) {
                for (int j = 0; j < 3; ++j) {
                    magmom_tmp[i][j] = 0.0;
                }
            }
            noncollinear_tmp = 0;
            time_reversal_symmetry_tmp = 1;
        }
        deallocate(magmom_tmp2);
    }

    MPI_Bcast(&lavec_s_tmp[0][0], 9, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nkd_tmp, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nat_tmp, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(&natmin_tmp, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ntran_tmp, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(&lspin_tmp, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&noncollinear_tmp, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&time_reversal_symmetry_tmp, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (mympi->my_rank > 0) {
        allocate(xr_s_tmp, nat_tmp, 3);
        allocate(kd_tmp, nat_tmp);
        allocate(map_p2s_tmp, natmin_tmp, ntran_tmp);
        allocate(map_s2p_tmp, nat_tmp);

        if (lspin_tmp) {
            allocate(magmom_tmp, natmin_tmp, 3);
        }
    }

    MPI_Bcast(&xr_s_tmp[0][0], 3 * nat_tmp, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&kd_tmp[0], nat_tmp, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(&map_p2s_tmp[0][0], natmin_tmp * ntran_tmp, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(&map_s2p_tmp[0], nat_tmp * sizeof map_s2p_tmp[0], MPI_BYTE, 0, MPI_COMM_WORLD);
    if (lspin_tmp) MPI_Bcast(&magmom_tmp[0][0], 3 * natmin_tmp, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    Eigen::Matrix3d lattice_vector, reciprocal_lattice_vector;
    Eigen::MatrixXd xf(nat_tmp, 3);
    std::vector<int> kind(nat_tmp);

    // Structure information
    scell_out.number_of_atoms = nat_tmp;
    scell_out.number_of_elems = nkd_tmp;
    for (auto i = 0; i < 3; ++i) {
        for (auto j = 0; j < 3; ++j) {
            lattice_vector(i, j) = lavec_s_tmp[i][j];
        }
    }
    reciprocal_lattice_vector = lattice_vector.inverse();
    scell_out.lattice_vector = lattice_vector;
    scell_out.reciprocal_lattice_vector = reciprocal_lattice_vector;
    for (auto i = 0; i < nat_tmp; ++i) {
        for (auto j = 0; j < 3; ++j) {
            xf(i, j) = xr_s_tmp[i][j];
        }
        kind[i] = kd_tmp[i];
    }
    scell_out.x_fractional = xf;
    scell_out.x_cartesian = xf * lattice_vector.transpose();
    scell_out.kind = kind;
    scell_out.volume = volume(lattice_vector, Direct);

    // Magnetism information
    spin_super_out.lspin = lspin_tmp;
    spin_super_out.time_reversal_symm = time_reversal_symmetry_tmp;
    spin_super_out.noncollinear = noncollinear_tmp;

    if (lspin_tmp) {
        spin_super_out.magmom.resize(nat_tmp, std::vector<double>(3));
        for (auto i = 0; i < nat_tmp; ++i) {
            for (auto j = 0; j < 3; ++j) {
                spin_super_out.magmom[i][j] = magmom_tmp[i][j];
            }
        }
    }

    // Mapping table
    map_super_out.to_true_primitive.resize(nat_tmp);
    map_super_out.from_true_primitive.resize(natmin_tmp, std::vector<size_t>(ntran_tmp));

    for (auto i = 0; i < nat_tmp; ++i) {
        map_super_out.to_true_primitive[i] = map_s2p_tmp[i];
    }

    for (auto i = 0; i < natmin_tmp; ++i) {
        for (auto j = 0; j < ntran_tmp; ++j) {
            map_super_out.from_true_primitive[i][j] = map_p2s_tmp[i][j];
        }
    }

    if (xr_s_tmp) deallocate(xr_s_tmp);
    if (kd_tmp) deallocate(kd_tmp);
    if (magmom_tmp) deallocate(magmom_tmp);
    if (map_p2s_tmp) deallocate(map_p2s_tmp);
    if (map_s2p_tmp) deallocate(map_s2p_tmp);
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
            exit("load_system_info_from_XML",
                 str_error.c_str());
        }

        // Parse nat_base and ntran_super

        nat = boost::lexical_cast<unsigned int>(
                get_value_from_xml(pt,
                                   "Data.Structure.NumberOfAtoms"));
        int nkd_tmp = boost::lexical_cast<unsigned int>(
                get_value_from_xml(pt,
                                   "Data.Structure.NumberOfElements"));

        if (nkd != nkd_tmp)
            exit("load_system_info_from_XML",
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
            ss >> lavec_s(0, i) >> lavec_s(1, i) >> lavec_s(2, i);
        }

        // Parse atomic elements and coordinates

        //allocate(xr_s, nat, 3);
        xr_s.resize(nat, 3);
        allocate(kd, nat);

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
                            exit("load_system_info_xml",
                                 "index is out of range");

                        kd[index] = dict_atomic_kind[str_element];
                        ss >> xr_s(index, 0) >> xr_s(index, 1) >> xr_s(index, 2);
                    }

        dict_atomic_kind.clear();

        // Parse mapping information

        allocate(map_p2s, natmin, ntran);
        allocate(map_s2p, nat);

        BOOST_FOREACH (const ptree::value_type &child_, pt.get_child("Data.Symmetry.Translations")) {
                        const auto &child = child_.second;
                        const auto str_tran = child.get<std::string>("<xmlattr>.tran");
                        const auto str_atom = child.get<std::string>("<xmlattr>.atom");

                        const auto tran = boost::lexical_cast<unsigned int>(str_tran) - 1;
                        const auto atom_p = boost::lexical_cast<unsigned int>(str_atom) - 1;
                        const auto atom_s = boost::lexical_cast<unsigned int>(child.data()) - 1;

                        if (tran >= ntran || atom_p >= natmin || atom_s >= nat) {
                            exit("load_system_info_xml",
                                 "index is out of range");
                        }

                        map_p2s[atom_p][tran] = atom_s;
                        map_s2p[atom_s].atom_num = atom_p;
                        map_s2p[atom_s].tran_num = tran;
                    }

        // Parse magnetic moments

        double **magmom_tmp;
        allocate(magmom_tmp, nat, 3);
        allocate(magmom, natmin, 3);

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
                                    exit("load_system_info_xml",
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
                symmetry->time_reversal_sym_from_alm = boost::lexical_cast<int>(
                        get_value_from_xml(pt,
                                           "Data.MagneticMoments.TimeReversalSymmetry"));
            }
            catch (...) {
                symmetry->time_reversal_sym_from_alm = true;
            }
        } else {
            for (i = 0; i < natmin; ++i) {
                for (int j = 0; j < 3; ++j) {
                    magmom[i][j] = 0.0;
                }
            }
            noncollinear = 0;
            symmetry->time_reversal_sym_from_alm = true;
        }
        deallocate(magmom_tmp);

        // Now, replicate the information for anharmonic terms.

        int j;
        nat_anharm = nat;
        ntran_anharm = ntran;
        //allocate(xr_s_anharm, nat_anharm, 3);
        xr_s_anharm.resize(nat_anharm, 3);
        allocate(kd_anharm, nat_anharm);
        allocate(map_p2s_anharm, natmin, ntran_anharm);
        allocate(map_s2p_anharm, nat_anharm);

        lavec_s_anharm = lavec_s;
//        for (i = 0; i < 3; ++i) {
//            for (j = 0; j < 3; ++j) lavec_s_anharm[i][j] = lavec_s[i][j];
//        }
        xr_s_anharm = xr_s;
        for (i = 0; i < nat_anharm; ++i) {
//            for (j = 0; j < 3; ++j) xr_s_anharm[i][j] = xr_s[i][j];
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
                exit("load_system_info_from_XML",
                     str_error.c_str());
            }

            // Parse nat_base and ntran_super

            nat = boost::lexical_cast<unsigned int>(
                    get_value_from_xml(pt,
                                       "Data.Structure.NumberOfAtoms"));
            nkd_tmp = boost::lexical_cast<unsigned int>(
                    get_value_from_xml(pt,
                                       "Data.Structure.NumberOfElements"));

            if (nkd != nkd_tmp)
                exit("load_system_info_from_XML",
                     "NKD in the FC2XML file is not consistent with that given in the input file.");

            ntran = boost::lexical_cast<unsigned int>(
                    get_value_from_xml(pt,
                                       "Data.Symmetry.NumberOfTranslations"));

            const int natmin_tmp = nat / ntran;

            if (natmin_tmp != natmin)
                exit("load_system_info_from_XML",
                     "Number of atoms in a primitive cell is different in FCSXML and FC2XML.");

            //deallocate(xr_s);
            deallocate(kd);
            deallocate(map_p2s);
            deallocate(map_s2p);


            // Parse lattice vectors

            std::stringstream ss;

            for (i = 0; i < 3; ++i) {
                ss.str("");
                ss.clear();
                ss << get_value_from_xml(pt,
                                         "Data.Structure.LatticeVector.a"
                                         + std::to_string(i + 1));
                ss >> lavec_s(0, i) >> lavec_s(1, i) >> lavec_s(2, i);
            }

            // Parse atomic elements and coordinates

            //allocate(xr_s, nat, 3);
            xr_s.resize(nat, 3);
            allocate(kd, nat);

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
                                exit("load_system_info_xml",
                                     "index is out of range");

                            kd[index_kd] = dict_atomic_kind[str_element];
                            ss >> xr_s(index_kd, 0) >> xr_s(index_kd, 1) >> xr_s(index_kd, 2);
                        }

            dict_atomic_kind.clear();

            // Parse mapping information

            allocate(map_p2s, natmin, ntran);
            allocate(map_s2p, nat);

            BOOST_FOREACH (const ptree::value_type &child_, pt.get_child("Data.Symmetry.Translations")) {
                            const auto &child = child_.second;
                            const auto str_tran = child.get<std::string>("<xmlattr>.tran");
                            const auto str_atom = child.get<std::string>("<xmlattr>.atom");

                            const auto tran = boost::lexical_cast<unsigned int>(str_tran) - 1;
                            const auto atom_p = boost::lexical_cast<unsigned int>(str_atom) - 1;
                            const auto atom_s = boost::lexical_cast<unsigned int>(child.data()) - 1;

                            if (tran >= ntran || atom_p >= natmin || atom_s >= nat) {
                                exit("load_system_info_xml", "index is out of range");
                            }

                            map_p2s[atom_p][tran] = atom_s;
                            map_s2p[atom_s].atom_num = atom_p;
                            map_s2p[atom_s].tran_num = tran;
                        }
        }

    }

    MPI_Bcast(lavec_s.data(), 9, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(lavec_p.data(), 9, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(lavec_s_anharm.data(), 9, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Bcast(&nkd, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nat, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nat_anharm, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(&natmin, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ntran, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ntran_anharm, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(&lspin, 1, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);

    if (mympi->my_rank > 0) {
        allocate(mass_kd, nkd);
        //allocate(xr_s, nat, 3);
        xr_s.resize(nat, 3);
        //allocate(xr_s_anharm, nat_anharm, 3);
        xr_s_anharm.resize(nat_anharm, 3);
        allocate(kd, nat);
        allocate(kd_anharm, nat_anharm);
        allocate(map_p2s, natmin, ntran);
        allocate(map_p2s_anharm, natmin, ntran_anharm);
        allocate(map_s2p, nat);
        allocate(map_s2p_anharm, nat_anharm);
        if (lspin) {
            allocate(magmom, natmin, 3);
        }
    }

    MPI_Bcast(&mass_kd[0], nkd, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(xr_s.data(), 3 * nat, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(xr_s_anharm.data(), 3 * nat_anharm, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&kd[0], nat, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(&kd_anharm[0], nat_anharm, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(&map_p2s[0][0], natmin * ntran, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(&map_p2s_anharm[0][0], natmin * ntran_anharm, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(&map_s2p[0], nat * sizeof map_s2p[0], MPI_BYTE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&map_s2p_anharm[0], nat_anharm * sizeof map_s2p_anharm[0], MPI_BYTE, 0, MPI_COMM_WORLD);
    if (lspin) MPI_Bcast(&magmom[0][0], 3 * natmin, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void System::recips(const Eigen::Matrix3d &mat_in,
                    Eigen::Matrix3d &rmat_out) const
{
    const auto det = mat_in.determinant();

    if (std::abs(det) < eps12) {
        exit("recips", "Lattice Vector is singular");
    }
    rmat_out = tpi * mat_in.inverse();
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
        exit("recips", "Lattice Vector is singular");
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

    allocate(atomlist_class, nclassatom);

    for (i = 0; i < N; ++i) {
        int count = 0;
        for (auto it: set_type) {
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
    Eigen::Vector3d xdiff;
    Eigen::MatrixXd x_harm(natmin, 3), x_anharm(natmin, 3);

    std::vector<int> map_anh2harm;
    map_anh2harm.resize(natmin);

    for (i = 0; i < natmin; ++i) {
        for (j = 0; j < 3; ++j) {
            x_harm(i, j) = xr_s(map_p2s[i][0], j);
            x_anharm(i, j) = xr_s_anharm(map_p2s_anharm[i][0], j);
        }
//        rotvec(x_harm[i], xr_s[map_p2s[i][0]], lavec_s);
//        rotvec(x_harm[i], x_harm[i], rlavec_p);
//        for (j = 0; j < 3; ++j) x_harm[i][j] /= 2.0 * pi;
    }

    x_harm = x_harm * lavec_s.transpose() * lavec_p.inverse().transpose();
    x_anharm = x_anharm * lavec_s_anharm.transpose() * lavec_p.inverse().transpose();
//
//    for (i = 0; i < natmin; ++i) {
//        rotvec(x_anharm[i], xr_s_anharm[map_p2s_anharm[i][0]], lavec_s_anharm);
//        rotvec(x_anharm[i], x_anharm[i], rlavec_p);
//        for (j = 0; j < 3; ++j) x_anharm[i][j] /= 2.0 * pi;
//    }

    for (i = 0; i < natmin; ++i) {

        int iloc = -1;

        for (j = 0; j < natmin; ++j) {

            xdiff = (x_anharm.row(i) - x_harm.row(j)).unaryExpr(
                    [](const double x) { return x - static_cast<double>(nint(x)); });

//            for (k = 0; k < 3; ++k) {
//                xdiff[k] = x_anharm[i][k] - x_harm[j][k];
//                xdiff[k] = xdiff[k] - static_cast<double>(nint(xdiff[k]));
//            }

            //const auto norm = xdiff[0] * xdiff[0] + xdiff[1] * xdiff[1] + xdiff[2] * xdiff[2];
            const auto norm = xdiff.squaredNorm();
            if (norm < eps4 && kd[map_p2s[j][0]] == kd_anharm[map_p2s_anharm[i][0]]) {
                iloc = j;
                break;
            }
        }

        if (iloc == -1) {
            exit("check_consistency_primitive",
                 "Could not find equivalent atom. Probably, the crystal structure is different.");
        }

        map_anh2harm[i] = iloc;
    }

    //deallocate(x_harm);
    //deallocate(x_anharm);

    // Rebuild the mapping information for anharmonic terms.

    unsigned int **map_p2s_tmp;

    allocate(map_p2s_tmp, natmin, ntran_anharm);

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

    deallocate(map_p2s_tmp);
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
            exit("set_mass_elem_from_database",
                 "Atomic mass for the given element doesn't exist in the database.\nTherefore, please input MASS manually.");
        }
        const auto mass_tmp = atomic_masses[atom_number];
        if (mass_tmp < 0.0) {
            exit("set_mass_elem_from_database",
                 "One of the elements in the KD-tag is unstable. \nTherefore, please input MASS manually.");
        }
        mass_kd_out[i] = mass_tmp;
    }
}

double System::volume(const Eigen::Matrix3d &mat_in,
                      const LatticeType latttype_in) const
{
    Eigen::Matrix3d mat;
    Eigen::Vector3d v1, v2, v3;

    if (latttype_in == Direct) {
        mat = mat_in.transpose();
    } else if (latttype_in == Reciprocal) {
        mat = mat_in;
    } else {
        exit("volume", "Invalid LatticeType is given");
    }

    v1 = mat.row(0);
    v2 = mat.row(1);
    v3 = mat.row(2);

    const auto vol = std::abs(v1.dot(v2.cross(v3)));
    return vol;
}

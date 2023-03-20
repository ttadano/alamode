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
#include "hdf5_parser.h"
#include <sstream>
#include <map>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <Eigen/LU>
#include <Eigen/Geometry>

#include <highfive/H5Easy.hpp>

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
    mass_kd = nullptr;
    symbol_kd = nullptr;
    map_p2s_anharm = nullptr;
    map_p2s_anharm_orig = nullptr;
    map_s2p = nullptr;
    map_s2p_anharm = nullptr;
    load_primitive_from_file = 0;
}

void System::deallocate_variables()
{
    if (mass_kd) {
        deallocate(mass_kd);
    }
    if (symbol_kd) {
        deallocate(symbol_kd);
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
    if (map_p2s_anharm_orig) {
        deallocate(map_p2s_anharm_orig);
    }
}

void System::setup()
{
    using namespace std;

    unsigned int i, j;
    double vec_tmp[3][3];

    MPI_Bcast(&load_primitive_from_file, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(lavec_p_input.data(), 9, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (mympi->my_rank == 0) {

        if (!mass_kd) {
            allocate(mass_kd, nkd);
            set_mass_elem_from_database(nkd, symbol_kd, mass_kd);
        }
    }
    load_system_info_from_file();
    update_primitive_lattice();
    generate_mapping_tables();

    //load_system_info_from_XML();

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

        cout << " * Primitive cell " << endl << endl;
        cout << setw(16) << primcell_base.lattice_vector(0, 0);
        cout << setw(15) << primcell_base.lattice_vector(1, 0);
        cout << setw(15) << primcell_base.lattice_vector(2, 0);
        cout << " : a1" << endl;

        cout << setw(16) << primcell_base.lattice_vector(0, 1);
        cout << setw(15) << primcell_base.lattice_vector(1, 1);
        cout << setw(15) << primcell_base.lattice_vector(2, 1);
        cout << " : a2" << endl;

        cout << setw(16) << primcell_base.lattice_vector(0, 2);
        cout << setw(15) << primcell_base.lattice_vector(1, 2);
        cout << setw(15) << primcell_base.lattice_vector(2, 2);
        cout << " : a3" << endl;
        cout << endl;

        cout << setw(16) << primcell_base.reciprocal_lattice_vector(0, 0);
        cout << setw(15) << primcell_base.reciprocal_lattice_vector(0, 1);
        cout << setw(15) << primcell_base.reciprocal_lattice_vector(0, 2);
        cout << " : b1" << endl;

        cout << setw(16) << primcell_base.reciprocal_lattice_vector(1, 0);
        cout << setw(15) << primcell_base.reciprocal_lattice_vector(1, 1);
        cout << setw(15) << primcell_base.reciprocal_lattice_vector(1, 2);
        cout << " : b2" << endl;

        cout << setw(16) << primcell_base.reciprocal_lattice_vector(2, 0);
        cout << setw(15) << primcell_base.reciprocal_lattice_vector(2, 1);
        cout << setw(15) << primcell_base.reciprocal_lattice_vector(2, 2);
        cout << " : b3" << endl;
        cout << endl << endl;

        volume_p = primcell_base.volume;

        cout << "  Volume of the primitive cell : "
             << primcell_base.volume << " (a.u.)^3" << endl << endl;
        cout << "  Number of atoms in the supercell     : "
             << supercell_base.number_of_atoms << endl;
        cout << "  Number of atoms in the primitive cell: "
             << primcell_base.number_of_atoms << endl << endl;

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

            cout << "  Number of atoms in the supercell (HARMONIC)   : "
                 << supercell_fc2.number_of_atoms << endl;
            cout << endl;
        }

        cout << "  Atomic positions in the primitive cell (fractional):" << endl;
        for (i = 0; i < primcell_base.number_of_atoms; ++i) {
            cout << setw(4) << i + 1 << ":";
            for (j = 0; j < 3; ++j) {
                cout << setw(15) << primcell_base.x_fractional(i, j);
            }
            cout << setw(4) << symbol_kd[primcell_base.kind[i]] << endl;
        }
        cout << endl;

        if (spin_prim_base.lspin) {
            cout << "  MagneticMoments entry found in the XML file. " << endl;
            cout << "  Magnetic moment in Cartesian coordinates: " << endl;
            for (i = 0; i < primcell_base.number_of_atoms; ++i) {
                cout << setw(4) << i + 1 << ":";
                for (j = 0; j < 3; ++j) {
                    cout << setw(15) << spin_prim_base.magmom[i][j];
                }
                cout << endl;
            }
            cout << endl;
            if (spin_prim_base.noncollinear == 0) {
                cout << "  Collinear calculation: magnetic moments are considered as scalar variables." << endl;
            } else if (spin_prim_base.noncollinear == 1) {
                cout << "  Noncollinear calculation: magnetic moments are considered as vector variables." << endl;
                if (spin_prim_base.time_reversal_symm) {
                    cout << "  Time-reversal symmetry will be considered for generating magnetic space group" << endl;
                } else {
                    cout << "  Time-reversal symmetry will NOT be considered for generating magnetic space group" <<
                         endl;
                }
            }
            cout << endl;
        }

        cout << "  Mass of atomic species (u):" << endl;
        for (i = 0; i < primcell_base.number_of_elems; ++i) {
            cout << setw(4) << symbol_kd[i] << ":";
            cout << fixed << setw(12) << mass_kd[i] << endl;
        }
        cout << endl << endl;
    }

    // Check the consistency of FCSXML and FC2XML
    MPI_Bcast(&fcs_phonon->update_fc2, 1, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);
//    if (fcs_phonon->update_fc2) {
//        allocate(map_p2s_anharm_orig, natmin, ntran_anharm);
//        check_consistency_primitive_lattice();
//    }
    // Atomic masses in Rydberg unit

    mass_super.resize(supercell_base.number_of_atoms);
    mass_prim.resize(primcell_base.number_of_atoms);
    invsqrt_mass_p.resize(primcell_base.number_of_atoms);

    for (i = 0; i < supercell_base.number_of_atoms; ++i) {
        mass_super[i] = mass_kd[supercell_base.kind[i]] * amu_ry;
    }
    for (i = 0; i < primcell_base.number_of_atoms; ++i) {
        mass_prim[i] = mass_kd[primcell_base.kind[i]] * amu_ry;
    }
    for (i = 0; i < primcell_base.number_of_atoms; ++i) {
        invsqrt_mass_p[i] = 1.0 / std::sqrt(mass_prim[i]);
    }

    MPI_Bcast(&Tmin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Tmax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&dT, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&volume_p, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}


void System::load_system_info_from_file()
{
    // Parse structure information either from the XML file or h5 file
    int filetype[4]; // filetype[0] for FCSFILE, filetype[X-1] for FCXFILE (X=2, 3, 4)
    // -1: FC?FILE not given, 0: FC?FILE in XML format, 1: FC?FLIE in HDF5 format
    std::vector<std::string> filename_list{fcs_phonon->file_fcs,
                                           fcs_phonon->file_fc2,
                                           fcs_phonon->file_fc3,
                                           fcs_phonon->file_fc4};

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

    Spin spin_super_fc2, spin_prim_fc2;
    std::vector<std::string> elements_base, elements_fc2, elements_tmp;

    for (auto i = 0; i < 4; ++i) {

        Cell scell, pcell;
        Spin spin_s, spin_p;
        MappingTable map_s, map_p;

        if (filetype[i] == -1) continue;

        if (filetype[i] == 0) {
            get_structure_and_mapping_table_xml(filename_list[i],
                                                scell, pcell,
                                                spin_s, spin_p,
                                                map_s, map_p,
                                                elements_tmp);
        } else if (filetype[i] == 1) {
            get_structure_and_mapping_table_h5(filename_list[i],
                                               scell, pcell,
                                               spin_s, spin_p,
                                               map_s, map_p,
                                               elements_tmp);
        } else {
            exit("load_system_info_from_file", "This cannot happen.");
        }

        if (i == 0) {
            supercell_base = scell;
            primcell_base = pcell;
            map_scell_base = map_s;
            map_pcell_base = map_p;
            spin_super_base = spin_s;
            spin_prim_base = spin_p;
            elements_base = elements_tmp;

        } else if (i == 1) {
            supercell_fc2 = scell;
            primcell_fc2 = pcell;
            map_scell_fc2 = map_s;
            map_pcell_fc2 = map_p;
            spin_super_fc2 = spin_s;
            spin_prim_fc2 = spin_p;
            elements_fc2 = elements_tmp;

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

    // If FCSFILE is not defined in input, set FC2FILE information as base.
    if (filetype[0] == -1) {
        supercell_base = supercell_fc2;
        primcell_base = primcell_fc2;
        map_scell_base = map_scell_fc2;
        map_pcell_base = map_pcell_fc2;
        spin_super_base = spin_super_fc2;
        spin_prim_base = spin_prim_fc2;
        elements_base = elements_fc2;
    }

    if (filetype[1] == -1) {
        supercell_fc2 = supercell_base;
        primcell_fc2 = primcell_base;
        map_scell_fc2 = map_scell_base;
        map_pcell_fc2 = map_pcell_base;
        spin_super_fc2 = spin_super_base;
        spin_prim_fc2 = spin_prim_base;
    }

    // Copy data if FC3FILE is not given.
    if (filetype[2] == -1) {
        supercell_fc3 = supercell_base;
        primcell_fc2 = primcell_base;
        map_scell_fc3 = map_scell_base;
        map_pcell_fc3 = map_pcell_base;
    }

    // Copy data if FC4FILE is not given.
    if (filetype[3] == -1) {
        supercell_fc4 = supercell_base;
        primcell_fc4 = primcell_base;
        map_scell_fc4 = map_scell_base;
        map_pcell_fc4 = map_pcell_base;
    }
}

void System::get_structure_and_mapping_table_xml(const std::string &filename,
                                                 Cell &scell_out,
                                                 Cell &pcell_out,
                                                 Spin &spin_super_out,
                                                 Spin &spin_prim_out,
                                                 MappingTable &map_super_out,
                                                 MappingTable &map_prim_out,
                                                 std::vector<std::string> &elements) const
{
    unsigned int nat_tmp, nkd_tmp, ntran_tmp, natmin_tmp;
    double lavec_s_tmp[3][3];
    double **xr_s_tmp = nullptr;
    int *kd_tmp = nullptr;
    unsigned int **map_p2s_tmp = nullptr;
    Maps *map_s2p_tmp = nullptr;
    double **magmom_tmp = nullptr;
    int lspin_tmp, noncollinear_tmp, time_reversal_symmetry_tmp;

    elements.clear();

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
                        elements.emplace_back(child_.second.data());
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
        elements.resize(nkd_tmp);
    }

    MPI_Bcast(&xr_s_tmp[0][0], 3 * nat_tmp, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&kd_tmp[0], nat_tmp, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    for (auto i = 0; i < nkd_tmp; ++i) {
        mympi->MPI_Bcast_string(elements[i], 0, MPI_COMM_WORLD);
    }

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
    reciprocal_lattice_vector = tpi * lattice_vector.inverse();
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
    map_super_out.from_true_primitive.resize(natmin_tmp, std::vector<unsigned int>(ntran_tmp));

    for (auto i = 0; i < nat_tmp; ++i) {
        map_super_out.to_true_primitive[i] = map_s2p_tmp[i];
    }

    for (auto i = 0; i < natmin_tmp; ++i) {
        for (auto j = 0; j < ntran_tmp; ++j) {
            map_super_out.from_true_primitive[i][j] = map_p2s_tmp[i][j];
        }
    }

    scell_out.has_entry = 1;

    if (xr_s_tmp) deallocate(xr_s_tmp);
    if (kd_tmp) deallocate(kd_tmp);
    if (magmom_tmp) deallocate(magmom_tmp);
    if (map_p2s_tmp) deallocate(map_p2s_tmp);
    if (map_s2p_tmp) deallocate(map_s2p_tmp);
}

void System::get_structure_and_mapping_table_h5(const std::string &filename,
                                                Cell &scell_out,
                                                Cell &pcell_out,
                                                Spin &spin_super_out,
                                                Spin &spin_prim_out,
                                                MappingTable &map_super_out,
                                                MappingTable &map_prim_out,
                                                std::vector<std::string> &elements) const
{
    using namespace H5Easy;

    int natmin_tmp, ntran_tmp;

    if (mympi->my_rank == 0) {
        File file(filename, File::ReadOnly);

        const std::string celltype_s = "SuperCell";
        const std::string celltype_p = "PrimitiveCell";

        get_structures_from_h5(file,
                               celltype_s,
                               scell_out.lattice_vector,
                               scell_out.x_fractional,
                               scell_out.kind,
                               elements);

        get_structures_from_h5(file,
                               celltype_p,
                               pcell_out.lattice_vector,
                               pcell_out.x_fractional,
                               pcell_out.kind,
                               elements);

        get_magnetism_from_h5(file,
                              celltype_s,
                              spin_super_out.lspin,
                              spin_super_out.magmom,
                              spin_super_out.noncollinear,
                              spin_super_out.time_reversal_symm);

        get_magnetism_from_h5(file,
                              celltype_p,
                              spin_prim_out.lspin,
                              spin_prim_out.magmom,
                              spin_prim_out.noncollinear,
                              spin_prim_out.time_reversal_symm);

        std::vector<std::vector<int>> mapping_table;
        get_mapping_table_from_h5(file, celltype_s, mapping_table);

        natmin_tmp = mapping_table.size();
        ntran_tmp = mapping_table[0].size();

        map_super_out.from_true_primitive.resize(natmin_tmp, std::vector<unsigned int>(ntran_tmp));
        map_super_out.to_true_primitive.resize(natmin_tmp * ntran_tmp);

        size_t atom_s;

        for (auto i = 0; i < natmin_tmp; ++i) {
            for (auto j = 0; j < ntran_tmp; ++j) {
                atom_s = mapping_table[i][j];
                map_super_out.from_true_primitive[i][j] = atom_s;
                map_super_out.to_true_primitive[atom_s].atom_num = i;
                map_super_out.to_true_primitive[atom_s].tran_num = j;
            }
        }

        mapping_table.clear();
        get_mapping_table_from_h5(file, celltype_p, mapping_table);

        natmin_tmp = mapping_table.size();
        ntran_tmp = mapping_table[0].size();

        map_prim_out.from_true_primitive.resize(natmin_tmp, std::vector<unsigned int>(ntran_tmp));
        map_prim_out.to_true_primitive.resize(natmin_tmp * ntran_tmp);

        for (auto i = 0; i < natmin_tmp; ++i) {
            for (auto j = 0; j < ntran_tmp; ++j) {
                atom_s = mapping_table[i][j];
                map_prim_out.from_true_primitive[i][j] = atom_s;
                map_prim_out.to_true_primitive[atom_s].atom_num = i;
                map_prim_out.to_true_primitive[atom_s].tran_num = j;
            }
        }
    }

    // Broadcast data
    mympi->MPI_Bcast_CellClass(scell_out, 0, MPI_COMM_WORLD);
    mympi->MPI_Bcast_CellClass(pcell_out, 0, MPI_COMM_WORLD);
    mympi->MPI_Bcast_SpinClass(spin_super_out, 0, MPI_COMM_WORLD);
    mympi->MPI_Bcast_SpinClass(spin_prim_out, 0, MPI_COMM_WORLD);
    mympi->MPI_Bcast_MappingTable(map_super_out, 0, MPI_COMM_WORLD);
    mympi->MPI_Bcast_MappingTable(map_prim_out, 0, MPI_COMM_WORLD);

    int nkd_tmp;
    if (mympi->my_rank == 0) {
        nkd_tmp = elements.size();
    }
    MPI_Bcast(&nkd_tmp, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (mympi->my_rank != 0) {
        elements.resize(nkd_tmp);
    }
    for (auto i = 0; i < nkd_tmp; ++i) {
        mympi->MPI_Bcast_string(elements[i], 0, MPI_COMM_WORLD);
    }

    // Fill in missing data
    scell_out.number_of_atoms = scell_out.x_fractional.rows();
    scell_out.number_of_elems = elements.size();
    pcell_out.number_of_atoms = pcell_out.x_fractional.rows();
    pcell_out.number_of_elems = elements.size();
    scell_out.reciprocal_lattice_vector = tpi * scell_out.lattice_vector.inverse();
    pcell_out.reciprocal_lattice_vector = tpi * pcell_out.lattice_vector.inverse();
    scell_out.x_cartesian = scell_out.x_fractional * scell_out.lattice_vector.transpose();
    pcell_out.x_cartesian = pcell_out.x_fractional * pcell_out.lattice_vector.transpose();

    scell_out.has_entry = 1;
    pcell_out.has_entry = 1;
}

//
//void System::load_system_info_from_XML()
//{
//    if (mympi->my_rank == 0) {
//
//        int i;
//        using namespace boost::property_tree;
//        ptree pt;
//
//        std::map<std::string, int> dict_atomic_kind;
//
//        try {
//            read_xml(fcs_phonon->file_fcs, pt);
//        }
//        catch (std::exception &e) {
//            std::string str_error = "Cannot open file FCSXML ( "
//                                    + fcs_phonon->file_fcs + " )";
//            exit("load_system_info_from_XML",
//                 str_error.c_str());
//        }
//
//        // Parse nat_base and ntran_super
//
//        auto nat_tmp = boost::lexical_cast<unsigned int>(
//                get_value_from_xml(pt,
//                                   "Data.Structure.NumberOfAtoms"));
//        int nkd_tmp = boost::lexical_cast<unsigned int>(
//                get_value_from_xml(pt,
//                                   "Data.Structure.NumberOfElements"));
//
//        if (nkd != nkd_tmp)
//            exit("load_system_info_from_XML",
//                 "NKD in the FCSXML file is not consistent with that given in the input file.");
//
//        ntran = boost::lexical_cast<unsigned int>(
//                get_value_from_xml(pt,
//                                   "Data.Symmetry.NumberOfTranslations"));
//
//        natmin = nat / ntran;
//
//        // Parse lattice vectors
//
//        std::stringstream ss;
//
//        for (i = 0; i < 3; ++i) {
//            ss.str("");
//            ss.clear();
//            ss << get_value_from_xml(pt,
//                                     "Data.Structure.LatticeVector.a"
//                                     + std::to_string(i + 1));
////            ss >> lavec_s(0, i) >> lavec_s(1, i) >> lavec_s(2, i);
//        }
//
//        // Parse atomic elements and coordinates
//
//        //xr_s.resize(nat, 3);
//        allocate(kd, nat);
//
//        BOOST_FOREACH (const ptree::value_type &child_, pt.get_child("Data.Structure.AtomicElements")) {
//                        const auto &child = child_.second;
//                        const auto icount_kd = child.get<unsigned int>("<xmlattr>.number");
//                        dict_atomic_kind[boost::lexical_cast<std::string>(child_.second.data())] = icount_kd - 1;
//                    }
//
//        unsigned int index;
//
//        BOOST_FOREACH (const ptree::value_type &child_, pt.get_child("Data.Structure.Position")) {
//                        const auto &child = child_.second;
//                        const auto str_index = child.get<std::string>("<xmlattr>.index");
//                        const auto str_element = child.get<std::string>("<xmlattr>.element");
//
//                        ss.str("");
//                        ss.clear();
//                        ss << child.data();
//
//                        index = boost::lexical_cast<unsigned int>(str_index) - 1;
//
//                        if (index >= nat)
//                            exit("load_system_info_xml",
//                                 "index is out of range");
//
//                        kd[index] = dict_atomic_kind[str_element];
////                        ss >> xr_s(index, 0) >> xr_s(index, 1) >> xr_s(index, 2);
//                    }
//
//        dict_atomic_kind.clear();
//
//        // Parse mapping information
//
//        allocate(map_p2s, natmin, ntran);
//        allocate(map_s2p, nat);
//
//        BOOST_FOREACH (const ptree::value_type &child_, pt.get_child("Data.Symmetry.Translations")) {
//                        const auto &child = child_.second;
//                        const auto str_tran = child.get<std::string>("<xmlattr>.tran");
//                        const auto str_atom = child.get<std::string>("<xmlattr>.atom");
//
//                        const auto tran = boost::lexical_cast<unsigned int>(str_tran) - 1;
//                        const auto atom_p = boost::lexical_cast<unsigned int>(str_atom) - 1;
//                        const auto atom_s = boost::lexical_cast<unsigned int>(child.data()) - 1;
//
//                        if (tran >= ntran || atom_p >= natmin || atom_s >= nat) {
//                            exit("load_system_info_xml",
//                                 "index is out of range");
//                        }
//
//                        map_p2s[atom_p][tran] = atom_s;
//                        map_s2p[atom_s].atom_num = atom_p;
//                        map_s2p[atom_s].tran_num = tran;
//                    }
//
//        // Parse magnetic moments
//
//        double **magmom_tmp;
//        allocate(magmom_tmp, nat, 3);
//
//        //lspin = true;
//        try {
//            BOOST_FOREACH(const ptree::value_type &child_, pt.get_child("Data.MagneticMoments")) {
//                            if (child_.first == "mag") {
//                                const auto &child = child_.second;
//                                const auto str_index = child.get<std::string>("<xmlattr>.index");
//
//                                ss.str("");
//                                ss.clear();
//                                ss << child.data();
//
//                                index = boost::lexical_cast<unsigned int>(str_index) - 1;
//
//                                if (index >= nat)
//                                    exit("load_system_info_xml",
//                                         "index is out of range");
//
//                                ss >> magmom_tmp[index][0]
//                                   >> magmom_tmp[index][1]
//                                   >> magmom_tmp[index][2];
//                            }
//                        }
//
//        }
//        catch (...) {
//            //    lspin = false;
//        }
//
//        deallocate(magmom_tmp);
//
//        // Now, replicate the information for anharmonic terms.
//
//        int j;
//        nat_anharm = nat;
//        ntran_anharm = ntran;
//        //xr_s_anharm.resize(nat_anharm, 3);
//        allocate(kd_anharm, nat_anharm);
//        allocate(map_p2s_anharm, natmin, ntran_anharm);
//        allocate(map_s2p_anharm, nat_anharm);
//
//        //lavec_s_anharm = lavec_s;
//        //xr_s_anharm = xr_s;
//        for (i = 0; i < nat_anharm; ++i) {
//            kd_anharm[i] = kd[i];
//            map_s2p_anharm[i] = map_s2p[i];
//        }
//        for (i = 0; i < natmin; ++i) {
//            for (j = 0; j < ntran_anharm; ++j) {
//                map_p2s_anharm[i][j] = map_p2s[i][j];
//            }
//        }
//
//        if (fcs_phonon->update_fc2) {
//
//            // When FC2XML is given, structural information is updated only for harmonic terms.
//
//            try {
//                read_xml(fcs_phonon->file_fc2, pt);
//            }
//            catch (std::exception &e) {
//                auto str_error = "Cannot open file FC2XML ( "
//                                 + fcs_phonon->file_fc2 + " )";
//                exit("load_system_info_from_XML",
//                     str_error.c_str());
//            }
//
//            // Parse nat_base and ntran_super
//
//            nat_tmp = boost::lexical_cast<unsigned int>(
//                    get_value_from_xml(pt,
//                                       "Data.Structure.NumberOfAtoms"));
//            nkd_tmp = boost::lexical_cast<unsigned int>(
//                    get_value_from_xml(pt,
//                                       "Data.Structure.NumberOfElements"));
//
//            if (nkd != nkd_tmp)
//                exit("load_system_info_from_XML",
//                     "NKD in the FC2XML file is not consistent with that given in the input file.");
//
//            ntran = boost::lexical_cast<unsigned int>(
//                    get_value_from_xml(pt,
//                                       "Data.Symmetry.NumberOfTranslations"));
//
//            const int natmin_tmp = nat / ntran;
//
//            if (natmin_tmp != natmin)
//                exit("load_system_info_from_XML",
//                     "Number of atoms in a primitive cell is different in FCSXML and FC2XML.");
//
//            deallocate(kd);
//            deallocate(map_p2s);
//            deallocate(map_s2p);
//
//
//            // Parse lattice vectors
//
//            std::stringstream ss;
//
//            for (i = 0; i < 3; ++i) {
//                ss.str("");
//                ss.clear();
//                ss << get_value_from_xml(pt,
//                                         "Data.Structure.LatticeVector.a"
//                                         + std::to_string(i + 1));
////                ss >> lavec_s(0, i) >> lavec_s(1, i) >> lavec_s(2, i);
//            }
//
//            // Parse atomic elements and coordinates
//
//            //xr_s.resize(nat, 3);
//            allocate(kd, nat);
//
//            BOOST_FOREACH (const ptree::value_type &child_, pt.get_child("Data.Structure.AtomicElements")) {
//                            const auto &child = child_.second;
//                            const auto icount_kd = child.get<unsigned int>("<xmlattr>.number");
//                            dict_atomic_kind[boost::lexical_cast<std::string>(child_.second.data())] = icount_kd - 1;
//                        }
//
//            BOOST_FOREACH (const ptree::value_type &child_, pt.get_child("Data.Structure.Position")) {
//                            const auto &child = child_.second;
//                            const auto str_index = child.get<std::string>("<xmlattr>.index");
//                            const auto str_element = child.get<std::string>("<xmlattr>.element");
//
//                            ss.str("");
//                            ss.clear();
//                            ss << child.data();
//
//                            auto index_kd = boost::lexical_cast<unsigned int>(str_index) - 1;
//
//                            if (index_kd >= nat)
//                                exit("load_system_info_xml",
//                                     "index is out of range");
//
//                            kd[index_kd] = dict_atomic_kind[str_element];
////                            ss >> xr_s(index_kd, 0) >> xr_s(index_kd, 1) >> xr_s(index_kd, 2);
//                        }
//
//            dict_atomic_kind.clear();
//
//            // Parse mapping information
//
//            allocate(map_p2s, natmin, ntran);
//            allocate(map_s2p, nat);
//
//            BOOST_FOREACH (const ptree::value_type &child_, pt.get_child("Data.Symmetry.Translations")) {
//                            const auto &child = child_.second;
//                            const auto str_tran = child.get<std::string>("<xmlattr>.tran");
//                            const auto str_atom = child.get<std::string>("<xmlattr>.atom");
//
//                            const auto tran = boost::lexical_cast<unsigned int>(str_tran) - 1;
//                            const auto atom_p = boost::lexical_cast<unsigned int>(str_atom) - 1;
//                            const auto atom_s = boost::lexical_cast<unsigned int>(child.data()) - 1;
//
//                            if (tran >= ntran || atom_p >= natmin || atom_s >= nat) {
//                                exit("load_system_info_xml", "index is out of range");
//                            }
//
//                            map_p2s[atom_p][tran] = atom_s;
//                            map_s2p[atom_s].atom_num = atom_p;
//                            map_s2p[atom_s].tran_num = tran;
//                        }
//        }
//
//    }
//
////    MPI_Bcast(lavec_s.data(), 9, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//    MPI_Bcast(lavec_p_input.data(), 9, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//    //MPI_Bcast(lavec_s_anharm.data(), 9, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//
//    MPI_Bcast(&nkd, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
//    MPI_Bcast(&nat, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
//    MPI_Bcast(&nat_anharm, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
//    MPI_Bcast(&natmin, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
//    MPI_Bcast(&ntran, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
//    MPI_Bcast(&ntran_anharm, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
//
//    if (mympi->my_rank > 0) {
//        allocate(mass_kd, nkd);
////        xr_s.resize(nat, 3);
////        xr_s_anharm.resize(nat_anharm, 3);
//        allocate(kd, nat);
//        allocate(kd_anharm, nat_anharm);
//        allocate(map_p2s, natmin, ntran);
//        allocate(map_p2s_anharm, natmin, ntran_anharm);
//        allocate(map_s2p, nat);
//        allocate(map_s2p_anharm, nat_anharm);
//    }
//
//    MPI_Bcast(&mass_kd[0], nkd, MPI_DOUBLE, 0, MPI_COMM_WORLD);
////    MPI_Bcast(xr_s.data(), 3 * nat, MPI_DOUBLE, 0, MPI_COMM_WORLD);
////    MPI_Bcast(xr_s_anharm.data(), 3 * nat_anharm, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//    MPI_Bcast(&kd[0], nat, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
//    MPI_Bcast(&kd_anharm[0], nat_anharm, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
//    MPI_Bcast(&map_p2s[0][0], natmin * ntran, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
//    MPI_Bcast(&map_p2s_anharm[0][0], natmin * ntran_anharm, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
//    MPI_Bcast(&map_s2p[0], nat * sizeof map_s2p[0], MPI_BYTE, 0, MPI_COMM_WORLD);
//    MPI_Bcast(&map_s2p_anharm[0], nat_anharm * sizeof map_s2p_anharm[0], MPI_BYTE, 0, MPI_COMM_WORLD);
//}

void System::update_primitive_lattice()
{
    // Update primcell_base and spin_prim_base.

    if (!primcell_base.has_entry && (load_primitive_from_file == 0)) {
        exit("System::update_primitive_lattice()",
             "Primitive lattice vectors must be given in the &cell field, \n"
             "because the corresponding data does not exist in FCSFILE.");
    }

    MPI_Bcast(lavec_p_input.data(), 9, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (load_primitive_from_file == 1) {
        // Update the information of primcell_base using the information of supercell_base
        // and the given lattice vector.

        Eigen::Matrix3d transmat_to_prim = supercell_base.lattice_vector.inverse() * lavec_p_input;

        primcell_base.lattice_vector = lavec_p_input;
        recips(primcell_base.lattice_vector, primcell_base.reciprocal_lattice_vector);
        primcell_base.volume = volume(primcell_base.lattice_vector, Direct);

        const auto ndiv = nint(1.0 / transmat_to_prim.determinant());
        if (supercell_base.number_of_atoms % ndiv != 0) {
            exit("update_primitive_lattice",
                 "The input primitive cell lattice vector is incommensurate \n "
                 " with the supercell information parsed from FCSFILE.");
        }

        primcell_base.number_of_atoms = supercell_base.number_of_atoms / ndiv;
        primcell_base.number_of_elems = supercell_base.number_of_elems;

        spin_prim_base.lspin = spin_super_base.lspin;
        spin_prim_base.noncollinear = spin_super_base.noncollinear;
        spin_prim_base.time_reversal_symm = spin_super_base.time_reversal_symm;
        spin_prim_base.magmom.clear();
        spin_prim_base.magmom.shrink_to_fit();

        // Convert the basis of coordinates from inputcell to primitive fractional
        // (a_in, b_in, c_in) * xf_in = xc_in
        // (a_p, b_p, c_p) * xf_p = xc_in
        // xf_p = (a_p, b_p, c_p)^{-1} * (a_in, b_in, c_in) * xf_in
        //      = [Mat(inp->p)]^{-1} * xf_in
        Eigen::Matrix3d conversion_mat = transmat_to_prim.inverse().transpose();
        Eigen::MatrixXd xf_prim_all = supercell_base.x_fractional * conversion_mat;

        std::vector<std::vector<double>> xf_unique, magmom_unique;
        std::vector<double> xf_tmp_vec(3), magmom_tmp_vec(3);
        std::vector<int> kind_unique;
        Eigen::VectorXd xf_tmp(3), xf_tmp2(3), xf_diff(3);

        for (auto i = 0; i < supercell_base.number_of_atoms; ++i) {
            xf_tmp = xf_prim_all.row(i);
            xf_tmp = xf_tmp.unaryExpr([](const double x) { return std::fmod(x, 1.0); });
            for (auto j = 0; j < 3; ++j) {
                if (xf_tmp[j] < -eps6) xf_tmp[j] += 1.0;
            }

            bool is_duplicate = false;
            // Just linear search for simplicity. Should be OK for relatively small number of inputs.
            for (auto k = 0; k < xf_unique.size(); ++k) {
                for (auto kk = 0; kk < 3; ++kk) {
                    xf_tmp2[kk] = xf_unique[k][kk];
                }
                xf_diff = (xf_tmp - xf_tmp2).unaryExpr([](const double x) { return std::fmod(x, 1.0); });
                for (auto j = 0; j < 3; ++j) {
                    if (xf_diff[j] < -0.5) xf_diff[j] += 1.0;
                    if (xf_diff[j] >= 0.5) xf_diff[j] -= 1.0;
                }
                if (xf_diff.norm() < eps6) {
                    is_duplicate = true;

                    if (kind_unique[k] != supercell_base.kind[i]) {
                        exit("update_primitive_lattice",
                             "Different atoms with different element types occupy the same atomic site.\n"
                             "This is strange. Please check the PRIMCELL and input structure carefully.");
                    }

                    if (spin_super_base.lspin) {
                        double norm_magmom = 0.0;
                        for (auto kk = 0; kk < 3; ++kk) {
                            norm_magmom += std::pow(spin_super_base.magmom[i][kk] - magmom_unique[k][kk], 2);
                        }
                        if (std::sqrt(norm_magmom) > eps6) {
                            exit("update_primitive_lattice",
                                 "Different atoms with different MAGMOM entries occupy the same atomic site.\n"
                                 "This is strange. Please check the PRIMCELL, MAGMOM, and input structure carefully.");
                        }
                    }
                }
            }

            if (!is_duplicate) {
                for (auto j = 0; j < 3; ++j) xf_tmp_vec[j] = xf_tmp[j];
                xf_unique.emplace_back(xf_tmp_vec);
                kind_unique.emplace_back(supercell_base.kind[i]);
                if (spin_super_base.lspin) {
                    for (auto j = 0; j < 3; ++j) magmom_tmp_vec[j] = spin_super_base.magmom[i][j];
                    magmom_unique.emplace_back(magmom_tmp_vec);
                }
            }
        }

        if (xf_unique.size() != primcell_base.number_of_atoms) {
            std::cout << "primcell.number_of_atoms = " << primcell_base.number_of_atoms << '\n';
            std::cout << "xf_unique.size() = " << xf_unique.size() << '\n';
            exit("update_primitive_lattice",
                 "Mapping to the primitive cell failed. "
                 "Please check the lattice constant in the &cell field.");
        }

        primcell_base.x_fractional.resize(primcell_base.number_of_atoms, 3);
        primcell_base.kind.resize(primcell_base.number_of_atoms);

        for (auto i = 0; i < primcell_base.number_of_atoms; ++i) {
            for (auto j = 0; j < 3; ++j) {
                primcell_base.x_fractional(i, j) = xf_unique[i][j];
            }
            primcell_base.kind[i] = kind_unique[i];
        }

        primcell_base.x_cartesian = primcell_base.x_fractional * primcell_base.lattice_vector.transpose();

        if (spin_prim_base.lspin) {
            std::copy(magmom_unique.begin(),
                      magmom_unique.end(),
                      std::back_inserter(spin_prim_base.magmom));
        }
    }
}


void System::generate_mapping_tables()
{
    std::vector<std::vector<unsigned int>> map_p2s_tmp;
    std::vector<Maps> map_s2p_tmp;

    map_s2p_new.resize(4);
    map_p2s_new.resize(4);

    generate_mapping_primitive_super(primcell_base, supercell_base,
                                     map_p2s_new[0], map_s2p_new[0]);

    generate_mapping_primitive_super(primcell_base, supercell_fc2,
                                     map_p2s_new[1], map_s2p_new[1]);

    generate_mapping_primitive_super(primcell_base, supercell_fc3,
                                     map_p2s_new[2], map_s2p_new[2]);

    generate_mapping_primitive_super(primcell_base, supercell_fc4,
                                     map_p2s_new[3], map_s2p_new[3]);
}


void System::generate_mapping_primitive_super(const Cell &pcell,
                                                 const Cell &scell,
                                                 std::vector<std::vector<unsigned int>> &map_p2s_out,
                                                 std::vector<Maps> &map_s2p_out) const
{

    Eigen::Vector3d x1, x2, x3, xdiff, tran_d;
    Eigen::Vector3i tran;

    const Eigen::MatrixXd x_super_in_primitive_frac = scell.x_fractional * scell.lattice_vector.transpose()
                                                * pcell.lattice_vector.inverse().transpose();
    const Eigen::Matrix3d transform_basis_primitive_to_super
            = scell.lattice_vector.inverse() * pcell.lattice_vector;

    const auto ntran_tmp = scell.number_of_atoms / pcell.number_of_atoms;

    if (scell.number_of_atoms != pcell.number_of_atoms * ntran_tmp) {
        exit("generate_mapping_primitive_super",
             "The number of atoms in the supercell is not divisible "
             "by the number of atoms in the primitive cell. Something is wrong.");
    }

    map_p2s_out.resize(pcell.number_of_atoms, std::vector<unsigned int>(ntran_tmp));

    std::vector<int> flag_found(scell.number_of_atoms, 0);
    std::map<int, int> map_index;
    std::vector<std::map<int, int>> map_index_lists;
    std::vector<int> atom_num_prim;
    std::vector<std::vector<int>> trans_vecs;

    atom_num_prim.clear();
    trans_vecs.clear();
    map_index_lists.clear();

    for (auto iat = 0; iat < scell.number_of_atoms; ++iat) {

        if (flag_found[iat]) continue;

        x1 = x_super_in_primitive_frac.row(iat);

        auto iloc = -1;
        for (auto jat = 0; jat < pcell.number_of_atoms; ++jat) {
            x2 = pcell.x_fractional.row(jat);
            xdiff = (x1 - x2).unaryExpr([](const double x) { return x - static_cast<double>(nint(x)); });

            if (xdiff.norm() < eps6 && (scell.kind[iat] == pcell.kind[jat])) {
                tran_d = (x1 - x2).unaryExpr([](const double x) { return static_cast<double>(nint(x)); });
                // First move back to the fractional coordinate of the supercell and
                // make sure that the shift vectors are in the 0<=x<1 region in that basis.
                tran_d = transform_basis_primitive_to_super * tran_d;
                tran_d = tran_d.unaryExpr([](const double x) { return x - static_cast<double>(nint(x)); });
                // Then, transform it back to the components in the primitive cell basis.
                // All components should be integer.
                tran_d = transform_basis_primitive_to_super.inverse() * tran_d;
                tran = tran_d.unaryExpr([](const double x) { return nint(x); });
                iloc = jat;
                break;
            }
        }

        if (iloc == -1) {
            exit("generate_mapping_primitive_super",
                 "An equivalent atom not found.");
        } else {
            atom_num_prim.emplace_back(iloc);
            std::vector<int> vtmp(&tran[0], tran.data() + tran.cols() * tran.rows());
            trans_vecs.emplace_back(vtmp);
        }

        tran_d = tran.unaryExpr([](const int x) { return static_cast<double>(x); });

        map_index.clear();

        // Now, we search all atoms in the supercell that can be mapped by the found trans value.
        for (auto j = 0; j < pcell.number_of_atoms; ++j) {
            x2 = pcell.x_fractional.row(j);
            x2 += tran_d;
            x2 = transform_basis_primitive_to_super * x2;

            for (auto k = 0; k < scell.number_of_atoms; ++k) {
                if (flag_found[k]) continue;

                x3 = scell.x_fractional.row(k);
                xdiff = (x3 - x2).unaryExpr([](const double x) { return x - static_cast<double>(nint(x)); });

                if (xdiff.norm() < eps6) {
                    flag_found[k] = 1;
                    map_index.insert({j, k});
                }
            }
        }
        map_index_lists.emplace_back(map_index);
    }

    std::set<std::vector<int>> unique_shifts_set;
    std::vector<std::vector<int>> unique_shifts_vec;
    std::vector<int> unique_itran_indices;

    for (auto itran = 0; itran < trans_vecs.size(); ++ itran) {
        if (unique_shifts_set.find(trans_vecs[itran]) == unique_shifts_set.end()) {

            std::cout << " itran = " << itran << ' ';
            for (const auto &it : trans_vecs[itran]) {
                std::cout << std::setw(4) << it;
            }
            std::cout << '\n';
            unique_shifts_set.insert(trans_vecs[itran]);
            unique_shifts_vec.emplace_back(trans_vecs[itran]);
            unique_itran_indices.emplace_back(itran);
        }
    }

    map_s2p_out.resize(scell.number_of_atoms);

    for (auto iat = 0; iat < pcell.number_of_atoms; ++iat) {
        for (auto itran = 0; itran < unique_itran_indices.size(); ++itran) {
            auto jat = map_index_lists[unique_itran_indices[itran]][iat];
            map_p2s_out[iat][itran] = jat;
            map_s2p_out[jat].atom_num = iat;
            map_s2p_out[jat].tran_num = itran;
        }
    }
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

//void System::check_consistency_primitive_lattice() const
//{
//    // Check if the ordering of atoms in the primitive cells derived
//    // from FCSXML and FC2XML are same or not. If not, the ordering for the
//    // FCSXML (anharmonic terms) will be changed so that it becomes equivalent
//    // to that of FC2XML.
//    // This operation is necessary for obtaining correct computational results.
//
//    int i, j, k;
//    Eigen::Vector3d xdiff;
//    Eigen::MatrixXd x_harm(natmin, 3), x_anharm(natmin, 3);
//
//    std::vector<int> map_anh2harm;
//    map_anh2harm.resize(natmin);
//
//    for (i = 0; i < natmin; ++i) {
//        for (j = 0; j < 3; ++j) {
////            x_harm(i, j) = xr_s(map_p2s[i][0], j);
////            x_anharm(i, j) = xr_s_anharm(map_p2s_anharm[i][0], j);
//        }
//    }
//
////    x_harm = x_harm * lavec_s.transpose() * lavec_p_input.inverse().transpose();
//    //x_anharm = x_anharm * lavec_s_anharm.transpose() * lavec_p_input.inverse().transpose();
//
//    for (i = 0; i < natmin; ++i) {
//
//        int iloc = -1;
//
//        for (j = 0; j < natmin; ++j) {
//
//            xdiff = (x_anharm.row(i) - x_harm.row(j)).unaryExpr(
//                    [](const double x) { return x - static_cast<double>(nint(x)); });
//
//            //const auto norm = xdiff[0] * xdiff[0] + xdiff[1] * xdiff[1] + xdiff[2] * xdiff[2];
//            const auto norm = xdiff.squaredNorm();
//            if (norm < eps4 && kd[map_p2s[j][0]] == kd_anharm[map_p2s_anharm[i][0]]) {
//                iloc = j;
//                break;
//            }
//        }
//
//        if (iloc == -1) {
//            exit("check_consistency_primitive",
//                 "Could not find equivalent atom. Probably, the crystal structure is different.");
//        }
//
//        map_anh2harm[i] = iloc;
//    }
//
//    // Rebuild the mapping information for anharmonic terms.
//
//    unsigned int **map_p2s_tmp;
//
//    allocate(map_p2s_tmp, natmin, ntran_anharm);
//
//    for (i = 0; i < natmin; ++i) {
//        for (j = 0; j < ntran_anharm; ++j) {
//            map_p2s_anharm_orig[i][j] = map_p2s_anharm[i][j];
//        }
//    }
//
//    for (i = 0; i < ntran_anharm; ++i) {
//        for (j = 0; j < natmin; ++j) {
//            map_p2s_tmp[j][i] = map_p2s_anharm[j][i];
//        }
//    }
//
//    for (i = 0; i < ntran_anharm; ++i) {
//        for (j = 0; j < natmin; ++j) {
//            map_p2s_anharm[map_anh2harm[j]][i] = map_p2s_tmp[j][i];
//        }
//    }
//
//    for (i = 0; i < ntran_anharm; ++i) {
//        for (j = 0; j < natmin; ++j) {
//            k = map_p2s_anharm[j][i];
//            map_s2p_anharm[k].atom_num = j;
//            map_s2p_anharm[k].tran_num = i;
//        }
//    }
//
//    deallocate(map_p2s_tmp);
//    map_anh2harm.clear();
//}

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

const Cell &System::get_primcell() const
{
    return primcell_base;
}


const Cell &System::get_cell(const std::string celltype,
                             const std::string filetype) const
{
    // This function returns Cell object of the given celltype and filetype.
    // It is very general but the overhead is high.
    // So, this function should NOT be called from the inside the loop.
    std::string celltype_tmp, filetype_tmp;
    if (celltype == "SuperCell" || celltype == "Super" || celltype == "super") {
        celltype_tmp = "super";
    } else if (celltype == "PrimCell" || celltype == "Prim" || celltype == "prim") {
        celltype_tmp = "prim";
    } else {
        exit("get_cell", "invalid celltype");
    }

    if (filetype == "base" || filetype == "default" || filetype == "Base") {
        filetype_tmp = "base";
    } else if (filetype == "fc2" || filetype == "FC2") {
        filetype_tmp = "fc2";
    } else if (filetype == "fc3" || filetype == "FC3") {
        filetype_tmp = "fc3";
    } else if (filetype == "fc4" || filetype == "FC4") {
        filetype_tmp = "fc4";
    } else {
        exit("get_cell", "invalid filetype");
    }

    if (celltype_tmp == "super" && filetype_tmp == "base") {
        return supercell_base;
    } else if (celltype_tmp == "prim" && filetype_tmp == "base") {
        return primcell_base;
    } else if (celltype_tmp == "super" && filetype_tmp == "fc2") {
        return supercell_fc2;
    } else if (celltype_tmp == "prim" && filetype_tmp == "fc2") {
        return primcell_fc2;
    } else if (celltype_tmp == "super" && filetype_tmp == "fc3") {
        return supercell_fc3;
    } else if (celltype_tmp == "prim" && filetype_tmp == "fc3") {
        return primcell_fc3;
    } else if (celltype_tmp == "super" && filetype_tmp == "fc4") {
        return supercell_fc4;
    } else if (celltype_tmp == "prim" && filetype_tmp == "fc4") {
        return primcell_fc4;
    }
    return supercell_base; // dummy for supressing compiler warning
}

const Spin &System::get_spin(const std::string celltype) const
{
    std::string celltype_tmp;
    if (celltype == "SuperCell" || celltype == "Super" || celltype == "super") {
        celltype_tmp = "super";
    } else if (celltype == "PrimCell" || celltype == "Prim" || celltype == "prim") {
        celltype_tmp = "prim";
    } else {
        exit("get_spin", "invalid celltype");
    }

    if (celltype_tmp == "super") {
        return spin_super_base;
    } else if (celltype_tmp == "prim") {
        return spin_prim_base;
    }
    return spin_super_base; // dummy for supressing compiler warning
}

const MappingTable &System::get_mapping_table(const std::string celltype,
                                              const std::string filetype) const
{
    std::string celltype_tmp, filetype_tmp;
    if (celltype == "SuperCell" || celltype == "Super" || celltype == "super") {
        celltype_tmp = "super";
    } else if (celltype == "PrimCell" || celltype == "Prim" || celltype == "prim") {
        celltype_tmp = "prim";
    } else {
        exit("get_mapping_table", "invalid celltype");
    }

    if (filetype == "base" || filetype == "default" || filetype == "Base") {
        filetype_tmp = "base";
    } else if (filetype == "fc2" || filetype == "FC2") {
        filetype_tmp = "fc2";
    } else if (filetype == "fc3" || filetype == "FC3") {
        filetype_tmp = "fc3";
    } else if (filetype == "fc4" || filetype == "FC4") {
        filetype_tmp = "fc4";
    } else {
        exit("get_mapping_table", "invalid filetype");
    }

    if (celltype_tmp == "super" && filetype_tmp == "base") {
        return map_scell_base;
    } else if (celltype_tmp == "prim" && filetype_tmp == "base") {
        return map_pcell_base;
    } else if (celltype_tmp == "super" && filetype_tmp == "fc2") {
        return map_scell_fc2;
    } else if (celltype_tmp == "prim" && filetype_tmp == "fc2") {
        return map_pcell_fc2;
    } else if (celltype_tmp == "super" && filetype_tmp == "fc3") {
        return map_scell_fc3;
    } else if (celltype_tmp == "prim" && filetype_tmp == "fc3") {
        return map_pcell_fc3;
    } else if (celltype_tmp == "super" && filetype_tmp == "fc4") {
        return map_scell_fc4;
    } else if (celltype_tmp == "prim" && filetype_tmp == "fc4") {
        return map_pcell_fc4;
    }

    return map_scell_base; // dummy for supressing compiler warning
}

const std::vector<double> &System::get_mass_prim() const
{
    return mass_prim;
}

const std::vector<double> &System::get_mass_super() const
{
    return mass_super;
}

const std::vector<double> &System::get_invsqrt_mass() const
{
    return invsqrt_mass_p;
}

const std::vector<Maps> &System::get_map_s2p(const int index) const
{
    return map_s2p_new[index];
}

const std::vector<std::vector<unsigned int>> &System::get_map_p2s(const int index) const
{
    return map_p2s_new[index];
}


/*
 system.cpp

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include "system.h"
#include "constants.h"
#include "mathfunctions.h"
#include "timer.h"
#include "memory.h"
#include "error.h"
#include "constraint.h"
#include "fcs.h"
#include "symmetry.h"
#include "fitting.h"
#include "xml_parser.h"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>
#include <boost/optional.hpp>
#include <boost/lexical_cast.hpp>

using namespace ALM_NS;

System::System(ALM *alm): Pointers(alm)
{
}

System::~System()
{
    memory->deallocate(x_cartesian);
    memory->deallocate(atomlist_class);
    memory->deallocate(magmom);
}

void System::init()
{
    using namespace std;

    int i, j;

    cout << " SYSTEM" << endl;
    cout << " ======" << endl << endl;

    recips(lavec, rlavec);

    cout.setf(ios::scientific);

    cout << "  Lattice Vector" << endl;
    cout << setw(16) << lavec[0][0];
    cout << setw(15) << lavec[1][0];
    cout << setw(15) << lavec[2][0];
    cout << " : a1" << endl;

    cout << setw(16) << lavec[0][1];
    cout << setw(15) << lavec[1][1];
    cout << setw(15) << lavec[2][1];
    cout << " : a2" << endl;

    cout << setw(16) << lavec[0][2];
    cout << setw(15) << lavec[1][2];
    cout << setw(15) << lavec[2][2];
    cout << " : a3" << endl;
    cout << endl;

    double vec_tmp[3][3];
    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            vec_tmp[i][j] = lavec[j][i];
        }
    }

    cell_volume = volume(vec_tmp[0], vec_tmp[1], vec_tmp[2]);
    cout << "  Cell volume = " << cell_volume << " (a.u)^3"
        << endl << endl;

    cout << "  Reciprocal Lattice Vector" << std::endl;
    cout << setw(16) << rlavec[0][0];
    cout << setw(15) << rlavec[0][1];
    cout << setw(15) << rlavec[0][2];
    cout << " : b1" << endl;

    cout << setw(16) << rlavec[1][0];
    cout << setw(15) << rlavec[1][1];
    cout << setw(15) << rlavec[1][2];
    cout << " : b2" << endl;

    cout << setw(16) << rlavec[2][0];
    cout << setw(15) << rlavec[2][1];
    cout << setw(15) << rlavec[2][2];
    cout << " : b3" << endl;
    cout << endl;

    cout << "  Atomic species:" << endl;
    for (i = 0; i < nkd; ++i) {
        cout << setw(6) << i + 1 << setw(5) << kdname[i] << endl;
    }
    cout << endl;

    cout << "  Atomic positions in fractional basis and atomic species" << endl;
    for (i = 0; i < nat; ++i) {
        cout << setw(6) << i + 1;
        cout << setw(15) << xcoord[i][0];
        cout << setw(15) << xcoord[i][1];
        cout << setw(15) << xcoord[i][2];
        cout << setw(5) << kd[i] << endl;
    }
    cout << endl << endl;
    cout.unsetf(ios::scientific);

    // Generate Cartesian coordinate

    memory->allocate(x_cartesian, nat, 3);

    for (i = 0; i < nat; ++i) {
        for (j = 0; j < 3; ++j) {
            x_cartesian[i][j] = xcoord[i][j];
        }
    }
    frac2cart(x_cartesian);
    setup_atomic_class(kd);

    if (lspin) {
        cout << "  MAGMOM is given. The magnetic moments of each atom are as follows:" << endl;
        for (i = 0; i < nat; ++i) {
            cout << setw(6) << i + 1;
            cout << setw(5) << magmom[i][0];
            cout << setw(5) << magmom[i][1];
            cout << setw(5) << magmom[i][2];
            cout << endl;
        }
        cout << endl;
        if (noncollinear == 0) {
            cout << "  NONCOLLINEAR = 0: magnetic moments are considered as scalar variables." << endl;
        } else if (noncollinear == 1) {
            cout << "  NONCOLLINEAR = 1: magnetic moments are considered as vector variables." << endl;
            if (symmetry->trev_sym_mag) {
                cout << "  TREVSYM = 1: Time-reversal symmetry will be considered for generating magnetic space group" << endl;
            } else {
                cout << "  TREVSYM = 0: Time-reversal symmetry will NOT be considered for generating magnetic space group" << endl;
            }
        }
        cout << endl << endl;
    }

    timer->print_elapsed();
    cout << " -------------------------------------------------------------------" << endl;
    cout << endl;
}

void System::recips(double aa[3][3], double bb[3][3])
{
    /*
    Calculate Reciprocal Lattice Vectors

    Here, BB is just the inverse matrix of AA (multiplied by factor 2 Pi)

    BB = 2 Pi AA^{-1},
    = t(b1, b2, b3)

    (b11 b12 b13)
    = (b21 b22 b23)
    (b31 b32 b33),

    b1 = t(b11, b12, b13) etc.
    */

    double det;
    det = aa[0][0] * aa[1][1] * aa[2][2]
        + aa[1][0] * aa[2][1] * aa[0][2]
        + aa[2][0] * aa[0][1] * aa[1][2]
        - aa[0][0] * aa[2][1] * aa[1][2]
        - aa[2][0] * aa[1][1] * aa[0][2]
        - aa[1][0] * aa[0][1] * aa[2][2];

    if (std::abs(det) < eps12) {
        error->exit("recips", "Lattice Vector is singular");
    }

    double factor = 2.0 * pi / det;

    bb[0][0] = (aa[1][1] * aa[2][2] - aa[1][2] * aa[2][1]) * factor;
    bb[0][1] = (aa[0][2] * aa[2][1] - aa[0][1] * aa[2][2]) * factor;
    bb[0][2] = (aa[0][1] * aa[1][2] - aa[0][2] * aa[1][1]) * factor;

    bb[1][0] = (aa[1][2] * aa[2][0] - aa[1][0] * aa[2][2]) * factor;
    bb[1][1] = (aa[0][0] * aa[2][2] - aa[0][2] * aa[2][0]) * factor;
    bb[1][2] = (aa[0][2] * aa[1][0] - aa[0][0] * aa[1][2]) * factor;

    bb[2][0] = (aa[1][0] * aa[2][1] - aa[1][1] * aa[2][0]) * factor;
    bb[2][1] = (aa[0][1] * aa[2][0] - aa[0][0] * aa[2][1]) * factor;
    bb[2][2] = (aa[0][0] * aa[1][1] - aa[0][1] * aa[1][0]) * factor;
}

void System::frac2cart(double **xf)
{
    // x_cartesian = A x_fractional

    int i, j;

    double *x_tmp;
    memory->allocate(x_tmp, 3);

    for (i = 0; i < nat; ++i) {

        rotvec(x_tmp, xf[i], lavec);

        for (j = 0; j < 3; ++j) {
            xf[i][j] = x_tmp[j];
        }
    }
    memory->deallocate(x_tmp);
}

void System::load_reference_system_xml(std::string file_reference_fcs,
                                       const int order_fcs,
                                       double *const_out)
{
    using namespace boost::property_tree;
    ptree pt;

    int nat_ref, natmin_ref, ntran_ref;
    int **intpair_ref;
    std::string str_error;
    double *fcs_ref;
    int nfcs_ref;

    try {
        read_xml(file_reference_fcs, pt);
    }
    catch (std::exception &e) {
        if (order_fcs == 0) {
            str_error = "Cannot open file FC2XML ( " + file_reference_fcs + " )";
        } else if (order_fcs == 1) {
            str_error = "Cannot open file FC3XML ( " + file_reference_fcs + " )";
        }
        error->exit("load_reference_system_xml", str_error.c_str());
    }

    nat_ref = boost::lexical_cast<int>(
        get_value_from_xml(pt, "Data.Structure.NumberOfAtoms"));
    ntran_ref = boost::lexical_cast<int>(
        get_value_from_xml(pt, "Data.Symmetry.NumberOfTranslations"));
    natmin_ref = nat_ref / ntran_ref;

    if (natmin_ref != symmetry->natmin) {
        error->exit("load_reference_system_xml",
                    "The number of atoms in the primitive cell is not consistent.");
    }

    if (order_fcs == 0) {
        nfcs_ref = boost::lexical_cast<int>(
            get_value_from_xml(pt, "Data.ForceConstants.HarmonicUnique.NFC2"));

        if (nfcs_ref != fcs->ndup[0].size()) {
            error->exit("load_reference_system_xml",
                        "The number of harmonic force constants is not the same.");
        }

    } else if (order_fcs == 1) {
        nfcs_ref = boost::lexical_cast<int>(
            get_value_from_xml(pt, "Data.ForceConstants.CubicUnique.NFC3"));

        if (nfcs_ref != fcs->ndup[1].size()) {
            error->exit("load_reference_system_xml",
                        "The number of cubic force constants is not the same.");
        }
    }
    memory->allocate(fcs_ref, nfcs_ref);
    memory->allocate(intpair_ref, nfcs_ref, 3);

    int counter = 0;

    if (order_fcs == 0) {
        BOOST_FOREACH (const ptree::value_type& child_, pt.get_child("Data.ForceConstants.HarmonicUnique")) {
                if (child_.first == "FC2") {
                    const ptree &child = child_.second;
                    const std::string str_intpair = child.get<std::string>("<xmlattr>.pairs");
                    const std::string str_multiplicity = child.get<std::string>("<xmlattr>.multiplicity");

                    std::istringstream is(str_intpair);
                    is >> intpair_ref[counter][0] >> intpair_ref[counter][1];
                    fcs_ref[counter] = boost::lexical_cast<double>(child.data());
                    ++counter;
                }
            }
    } else if (order_fcs == 1) {
        BOOST_FOREACH (const ptree::value_type& child_, pt.get_child("Data.ForceConstants.CubicUnique")) {
                if (child_.first == "FC3") {
                    const ptree &child = child_.second;
                    const std::string str_intpair = child.get<std::string>("<xmlattr>.pairs");
                    const std::string str_multiplicity = child.get<std::string>("<xmlattr>.multiplicity");

                    std::istringstream is(str_intpair);
                    is >> intpair_ref[counter][0] >> intpair_ref[counter][1] >> intpair_ref[counter][2];
                    fcs_ref[counter] = boost::lexical_cast<double>(child.data());
                    ++counter;
                }
            }
    }

    int i;
    std::set<FcProperty> list_found;
    std::set<FcProperty>::iterator iter_found;
    int *ind;
    int nterms = order_fcs + 2;
    memory->allocate(ind, nterms);

    list_found.clear();

    for (std::vector<FcProperty>::iterator p = fcs->fc_set[order_fcs].begin();
         p != fcs->fc_set[order_fcs].end(); ++p) {
        FcProperty list_tmp = *p; // Using copy constructor
        for (i = 0; i < nterms; ++i) {
            ind[i] = list_tmp.elems[i];
        }
        list_found.insert(FcProperty(nterms, list_tmp.coef,
                                     ind, list_tmp.mother));
    }

    for (i = 0; i < nfcs_ref; ++i) {
        iter_found = list_found.find(FcProperty(nterms, 1.0,
                                                intpair_ref[i], 1));
        if (iter_found == list_found.end()) {
            error->exit("load_reference_system",
                        "Cannot find equivalent force constant, number: ",
                        i + 1);
        }
        FcProperty arrtmp = *iter_found;
        const_out[arrtmp.mother] = fcs_ref[i];
    }

    memory->deallocate(intpair_ref);
    memory->deallocate(fcs_ref);
    memory->deallocate(ind);
    list_found.clear();
}

void System::load_reference_system()
{
    int i;
    int iat, jat;
    int icrd;

    unsigned int nat_s, nkd_s;
    unsigned int natmin_ref, ntran_ref;
    double lavec_s[3][3];
    int *kd_s;
    double **xcoord_s;
    int *map_ref;
    int **map_p2s_s;
    Symmetry::Maps *map_s2p_s;
    std::ifstream ifs_fc2;

    ifs_fc2.open(constraint->fc2_file.c_str(), std::ios::in);
    if (!ifs_fc2)
        error->exit("load_reference_system",
                    "cannot open file fc2_file");

    bool is_found_system = false;

    int nparam_harmonic_ref;
    int nparam_harmonic = fcs->ndup[0].size();

    std::string str_tmp;

    while (!ifs_fc2.eof() && !is_found_system) {
        std::getline(ifs_fc2, str_tmp);
        if (str_tmp == "##SYSTEM INFO") {

            is_found_system = true;

            std::getline(ifs_fc2, str_tmp);
            for (i = 0; i < 3; ++i) {
                ifs_fc2 >> lavec_s[0][i] >> lavec_s[1][i] >> lavec_s[2][i];
            }
            ifs_fc2.ignore();
            std::getline(ifs_fc2, str_tmp);
            ifs_fc2 >> nkd_s;
            ifs_fc2.ignore();
            std::getline(ifs_fc2, str_tmp);
            std::getline(ifs_fc2, str_tmp);

            ifs_fc2 >> nat_s >> natmin_ref >> ntran_ref;

            if (natmin_ref != symmetry->natmin) {
                error->exit("load_reference_system",
                            "The number of atoms in the primitive cell is not consistent");
            }

            if (nat_s != nat) {
                std::cout << "The number of atoms in the reference system differs from input." << std::endl;
                std::cout << "Trying to map the related force constants (^o^)" << std::endl << std::endl;
            }

            memory->allocate(xcoord_s, nat_s, 3);
            memory->allocate(kd_s, nat_s);
            memory->allocate(map_p2s_s, natmin_ref, ntran_ref);
            memory->allocate(map_s2p_s, nat_s);

            unsigned int ikd, itran, icell;
            std::getline(ifs_fc2, str_tmp);
            std::getline(ifs_fc2, str_tmp);
            for (i = 0; i < nat_s; ++i) {
                ifs_fc2 >> str_tmp >> ikd
                    >> xcoord_s[i][0] >> xcoord_s[i][1] >> xcoord_s[i][2] >> itran >> icell;
                kd_s[i] = ikd;
                map_p2s_s[icell - 1][itran - 1] = i;
                map_s2p_s[i].atom_num = icell - 1;
                map_s2p_s[i].tran_num = itran - 1;
            }
        }
    }
    if (!is_found_system) {
        error->exit("load_reference_system",
                    "SYSTEM INFO flag not found in the fc2_file");
    }

    //
    // Generate Mapping Information (big supercell -> small supercell)
    //

    double *xtmp;
    double *xdiff;
    int **intpair_tmp;

    memory->allocate(xtmp, 3);
    memory->allocate(xdiff, 3);
    memory->allocate(map_ref, nat_s);

    bool map_found;
    double dist;

    for (iat = 0; iat < nat_s; ++iat) {
        map_found = false;

        rotvec(xtmp, xcoord_s[iat], lavec_s);
        rotvec(xtmp, xtmp, rlavec);

        for (icrd = 0; icrd < 3; ++icrd) xtmp[icrd] /= 2.0 * pi;

        for (jat = 0; jat < nat; ++jat) {
            for (icrd = 0; icrd < 3; ++icrd) {
                xdiff[icrd] = xtmp[icrd] - xcoord[jat][icrd];
                xdiff[icrd] = std::fmod(xdiff[icrd], 1.0);
            }
            dist = xdiff[0] * xdiff[0] + xdiff[1] * xdiff[1] + xdiff[2] * xdiff[2];

            if (dist < eps12 && kd_s[iat] == kd[jat]) {
                map_ref[iat] = jat;
                map_found = true;
                break;
            }
        }
        if (!map_found) {
            error->exit("load_reference_system",
                        "Could not find an equivalent atom for atom ",
                        iat + 1);
        }
    }

    memory->deallocate(xtmp);
    memory->deallocate(xdiff);
    memory->deallocate(xcoord_s);
    memory->deallocate(kd_s);

    ifs_fc2.clear();
    ifs_fc2.seekg(0, std::ios_base::beg);

    double *fc2_ref;

    bool is_found_fc2 = false;

    while (!ifs_fc2.eof() && !is_found_fc2) {
        std::getline(ifs_fc2, str_tmp);
        if (str_tmp == "##HARMONIC FORCE CONSTANTS") {
            ifs_fc2 >> nparam_harmonic_ref;
            if (nparam_harmonic_ref < nparam_harmonic) {
                error->exit("load_reference_system",
                            "Reference file doesn't contain necessary fc2. (too few)");
            } else if (nparam_harmonic_ref > nparam_harmonic) {
                error->exit("load_reference_system",
                            "Reference file contains extra force constants.");
            }

            is_found_fc2 = true;

            memory->allocate(fc2_ref, nparam_harmonic);
            memory->allocate(intpair_tmp, nparam_harmonic, 2);

            for (i = 0; i < nparam_harmonic; ++i) {
                ifs_fc2 >> fc2_ref[i] >> intpair_tmp[i][0] >> intpair_tmp[i][1];
            }

            std::set<FcProperty> list_found;
            std::set<FcProperty>::iterator iter_found;
            int *ind;
            memory->allocate(ind, 2);

            list_found.clear();
            for (std::vector<FcProperty>::iterator p = fcs->fc_set[0].begin();
                 p != fcs->fc_set[0].end(); ++p) {
                for (i = 0; i < 2; ++i) ind[i] = (*p).elems[i];
                list_found.insert(FcProperty(2, (*p).coef, ind, (*p).mother));
            }

            for (i = 0; i < nparam_harmonic; ++i) {
                constraint->const_mat[i][i] = 1.0;
            }

            for (i = 0; i < nparam_harmonic; ++i) {

                iter_found = list_found.find(FcProperty(2, 1.0, intpair_tmp[i], 1));
                if (iter_found == list_found.end()) {
                    error->exit("load_reference_system",
                                "Cannot find equivalent force constant, number: ",
                                i + 1);
                }
                constraint->const_rhs[(*iter_found).mother] = fc2_ref[i];
            }

            memory->deallocate(intpair_tmp);
            memory->deallocate(ind);
            memory->deallocate(fc2_ref);
            list_found.clear();
        }
    }

    if (!is_found_fc2) {
        error->exit("load_reference_system",
                    "HARMONIC FORCE CONSTANTS flag not found in the fc2_file");
    }
    ifs_fc2.close();
}

double System::volume(double vec1[3], double vec2[3], double vec3[3])
{
    double vol;

    vol = std::abs(vec1[0] * (vec2[1] * vec3[2] - vec2[2] * vec3[1])
        + vec1[1] * (vec2[2] * vec3[0] - vec2[0] * vec3[2])
        + vec1[2] * (vec2[0] * vec3[1] - vec2[1] * vec3[0]));

    return vol;
}

void System::setup_atomic_class(int *kd)
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

    for (i = 0; i < nat; ++i) {
        type_tmp.element = kd[i];

        if (noncollinear == 0) {
            type_tmp.magmom = magmom[i][2];
        } else {
            type_tmp.magmom = 0.0;
        }
        set_type.insert(type_tmp);
    }

    nclassatom = set_type.size();

    memory->allocate(atomlist_class, nclassatom);

    for (i = 0; i < nat; ++i) {
        int count = 0;
        for (std::set<AtomType>::iterator it = set_type.begin();
             it != set_type.end(); ++it) {
            if (noncollinear) {
                if (kd[i] == (*it).element) {
                    atomlist_class[count].push_back(i);
                }
            } else {
                if ((kd[i] == (*it).element)
                    && (std::abs(magmom[i][2] - (*it).magmom) < eps6)) {
                    atomlist_class[count].push_back(i);
                }
            }
            ++count;
        }
    }
    set_type.clear();
}

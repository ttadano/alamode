/*
 system.cpp

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include "system.h"
#include "constants.h"
#include "error.h"
#include "mathfunctions.h"
#include "memory.h"
#include "timer.h"
#include "smith.h"
#include <iostream>
#include <iomanip>
#include <set>
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Geometry>

using namespace ALM_NS;

System::System()
{
    set_default_variables();
}

System::~System()
{
    deallocate_variables();
}

void System::init(const int verbosity,
                  std::unique_ptr<Timer> &timer)
{
    timer->start_clock("system");

    // Compute primitive cell and supercell first.
    build_cells();
    const auto nat = supercell.number_of_atoms;

    // Set atomic types (kind + magmom)
    set_atomtype_group();

    const auto nneib = 27;
    if (x_image) {
        deallocate(x_image);
    }
    allocate(x_image, nneib, nat, 3);

    if (exist_image) {
        deallocate(exist_image);
    }
    allocate(exist_image, nneib);
    generate_coordinate_of_periodic_images();

    if (verbosity > 0) {
        print_structure_stdout(verbosity);
        if (spin_super.lspin) print_magmom_stdout();
        timer->print_elapsed();
        std::cout << " -------------------------------------------------------------------" << std::endl;
        std::cout << std::endl;
    }

    timer->stop_clock("system");
}

void System::set_basecell(const double lavec_in[3][3],
                          const size_t nat_in,
                          const int *kind_in,
                          const double xf_in[][3])
{
    size_t i, j;
    std::vector<int> unique_nums(nat_in);
    auto wrong_number = false;
    bool in_unique_nums;

    for (i = 0; i < nat_in; i++) {
        unique_nums[i] = 0;
    }

    size_t nkd = 0;
    for (i = 0; i < nat_in; i++) {
        in_unique_nums = false;
        for (j = 0; j < nkd; j++) {
            if (unique_nums[j] == kind_in[i]) {
                in_unique_nums = true;
                break;
            }
        }
        if (!in_unique_nums) {
            unique_nums[nkd] = kind_in[i];
            nkd++;
        }
    }

    for (i = 0; i < nkd; i++) {
        if (static_cast<size_t>(unique_nums[i]) > nkd) {
            std::cout << " WARNING : integers assigned to atoms are wrong. "
                      << " The numbers will be resorted." << std::endl;
            wrong_number = true;
            break;
        }
    }

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            inputcell.lattice_vector(i, j) = lavec_in[i][j];
        }
    }
    set_reciprocal_latt(inputcell.lattice_vector,
                        inputcell.reciprocal_lattice_vector);

    inputcell.volume = volume(inputcell.lattice_vector, Direct);
    inputcell.number_of_atoms = nat_in;
    inputcell.number_of_elems = nkd;
    inputcell.kind.clear();
    inputcell.kind.shrink_to_fit();
    inputcell.x_fractional.resize(nat_in, 3);

    std::vector<double> xtmp;

    if (!wrong_number) {
        for (i = 0; i < nat_in; ++i) {
            inputcell.kind.push_back(kind_in[i]);
        }
    } else {
        for (i = 0; i < nat_in; ++i) {
            for (j = 0; j < nkd; j++) {
                if (kind_in[i] == unique_nums[j]) {
                    inputcell.kind.push_back(static_cast<int>(j + 1));
                }
            }
        }
    }

    xtmp.resize(3);
    for (i = 0; i < nat_in; ++i) {
        for (j = 0; j < 3; ++j) {
            xtmp[j] = xf_in[i][j];
        }
        // The fractional coordinate should be in the range of 0<=xf<1.
        for (j = 0; j < 3; ++j) {
//            while (xtmp[j] >= 1.0) {
//                xtmp[j] -= 1.0;
//            }
//            while (xtmp[j] < 0.0) {
//                xtmp[j] += 1.0;
//            }
            inputcell.x_fractional(i, j) = xtmp[j];
        }
    }

    inputcell.x_cartesian = inputcell.x_fractional * inputcell.lattice_vector.transpose();

    // This is needed to avoid segmentation fault.
    spin_input.magmom.clear();
    std::vector<double> vec(3);
    for (i = 0; i < nat_in; ++i) {
        for (j = 0; j < 3; ++j) {
            vec[j] = 0;
        }
        spin_input.magmom.push_back(vec);
    }
}

void System::build_cells()
{
    build_primcell();
    build_supercell();
}


void System::build_primcell()
{
    // Generate primitive cell information from the inputcell and PRIMCELL value
    // The symmetry detection is not performed here as it will be done
    // when the Symmetry::init() is called.

    const auto ndiv = nint(1.0 / transmat_to_prim.determinant());
    if (inputcell.number_of_atoms % ndiv != 0) {
        exit("build_primcell",
             "The determinant of PRIMCELL is not a divisor of NAT of the input.");
    }

    primcell.number_of_atoms = inputcell.number_of_atoms / ndiv;
    primcell.number_of_elems = inputcell.number_of_elems;

    spin_prim.lspin = spin_input.lspin;
    spin_prim.noncollinear = spin_input.noncollinear;
    spin_prim.time_reversal_symm = spin_input.time_reversal_symm;
    spin_prim.magmom.clear();
    spin_prim.magmom.shrink_to_fit();

    // (a_p, b_p, c_p) = (a_in, b_in, c_in) * Mat(inp->p)
    primcell.lattice_vector = inputcell.lattice_vector * transmat_to_prim;
    set_reciprocal_latt(primcell.lattice_vector,
                        primcell.reciprocal_lattice_vector);
    primcell.volume = volume(primcell.lattice_vector, Direct);

    // Convert the basis of coordinates from inputcell to primitive fractional
    // (a_in, b_in, c_in) * xf_in = xc_in
    // (a_p, b_p, c_p) * xf_p = xc_in
    // xf_p = (a_p, b_p, c_p)^{-1} * (a_in, b_in, c_in) * xf_in
    //      = [Mat(inp->p)]^{-1} * xf_in
    Eigen::Matrix3d conversion_mat = transmat_to_prim.inverse().transpose();
    Eigen::MatrixXd xf_prim_all = inputcell.x_fractional * conversion_mat;

    std::vector<std::vector<double>> xf_unique, magmom_unique;
    std::vector<double> xf_tmp_vec(3), magmom_tmp_vec(3);
    std::vector<int> kind_unique;
    Eigen::VectorXd xf_tmp(3), xf_tmp2(3);

    for (auto i = 0; i < inputcell.number_of_atoms; ++i) {
        xf_tmp = xf_prim_all.row(i);
        xf_tmp = xf_tmp.unaryExpr([](const double x) { return std::fmod(x, 1.0); });
        for (auto j = 0; j < 3; ++j) {
            if (xf_tmp[j] < 0.0) xf_tmp[j] += 1.0;
        }

        bool is_duplicate = false;
        // Just linear search for simplicity. Should be OK for relatively small number of inputs.
        for (auto k = 0; k < xf_unique.size(); ++k) {
            for (auto kk = 0; kk < 3; ++kk) {
                xf_tmp2[kk] = xf_unique[k][kk];
            }
            if ((xf_tmp - xf_tmp2).norm() < eps6) {
                is_duplicate = true;

                if (kind_unique[k] != inputcell.kind[i]) {
                    exit("build_primcell",
                         "Different atoms with different element types occupy the same atomic site.\n"
                         "This is strange. Please check the PRIMCELL and input structure carefully.");
                }

                if (spin_input.lspin) {
                    double norm_magmom = 0.0;
                    for (auto kk = 0; kk < 3; ++kk) {
                        norm_magmom += std::pow(spin_input.magmom[i][kk] - magmom_unique[k][kk], 2);
                    }
                    if (std::sqrt(norm_magmom) > eps6) {
                        exit("build_primcell",
                             "Different atoms with different MAGMOM entries occupy the same atomic site.\n"
                             "This is strange. Please check the PRIMCELL, MAGMOM, and input structure carefully.");
                    }
                }
            }
        }

        if (!is_duplicate) {
            for (auto j = 0; j < 3; ++j) xf_tmp_vec[j] = xf_tmp[j];
            xf_unique.emplace_back(xf_tmp_vec);
            kind_unique.emplace_back(inputcell.kind[i]);
            if (spin_input.lspin) {
                for (auto j = 0; j < 3; ++j) magmom_tmp_vec[j] = spin_input.magmom[i][j];
                magmom_unique.emplace_back(magmom_tmp_vec);
            }
        }
    }

    if (xf_unique.size() != primcell.number_of_atoms) {
        exit("build_primcell",
             "Mapping to the primitive cell failed. Please check inputs (PRIMCELL and input structure).");
    }

    primcell.x_fractional.resize(primcell.number_of_atoms, 3);
    primcell.kind.resize(primcell.number_of_atoms);

    for (auto i = 0; i < primcell.number_of_atoms; ++i) {
        for (auto j = 0; j < 3; ++j) {
            primcell.x_fractional(i, j) = xf_unique[i][j];
        }
        primcell.kind[i] = kind_unique[i];
    }

    primcell.x_cartesian = primcell.x_fractional * primcell.lattice_vector.transpose();

    if (spin_prim.lspin) {
        std::copy(magmom_unique.begin(),
                  magmom_unique.end(),
                  std::back_inserter(spin_prim.magmom));
    }
}


void System::build_supercell()
{
    Eigen::MatrixXi transmat_int(3, 3);
    Eigen::MatrixXi D, R, L;

    for (auto i = 0; i < 3; ++i) {
        for (auto j = 0; j < 3; ++j) {
            transmat_int(i, j) = nint(transmat_to_super(i, j));
        }
    }
    smith_decomposition(transmat_int, D, L, R);

    Eigen::MatrixXd D_double = D.cast<double>();
    Eigen::MatrixXd R_double = R.cast<double>();
    Eigen::MatrixXd L_double = L.cast<double>();

    supercell.lattice_vector = inputcell.lattice_vector * transmat_to_super;
    set_reciprocal_latt(supercell.lattice_vector, supercell.reciprocal_lattice_vector);
    supercell.volume = volume(supercell.lattice_vector, Direct);

    supercell.number_of_atoms = inputcell.number_of_atoms * nint(transmat_to_super.determinant());
    supercell.number_of_elems = inputcell.number_of_elems;
    supercell.x_fractional.resize(supercell.number_of_atoms, 3);

    Eigen::Matrix3d inv_transmat_to_super = transmat_to_super.inverse();

    Eigen::VectorXd DD(3), DD_mod(3);
    Eigen::MatrixXd DD_inverse = D_double.inverse();
    Eigen::MatrixXd transmat_newprim_to_origsuper = R_double * DD_inverse;

    supercell.kind.clear();

    std::vector<std::vector<double>> magmom_tmp;

    int counter = 0;
    for (auto iat = 0; iat < inputcell.number_of_atoms; ++iat) {
        Eigen::VectorXd xp_f = inputcell.x_fractional.row(iat);
        Eigen::VectorXd xs_f = inv_transmat_to_super * xp_f;
        int kind_now = inputcell.kind[iat];
        for (auto i = 0; i < D(0, 0); ++i) {
            for (auto j = 0; j < D(1, 1); ++j) {
                for (auto k = 0; k < D(2, 2); ++k) {
                    DD << static_cast<double>(i), static_cast<double>(j), static_cast<double>(k);
                    DD = transmat_newprim_to_origsuper * DD + xs_f;
                    DD_mod = DD.unaryExpr([](const double x) { return std::fmod(x, 1.0); });
                    supercell.x_fractional.row(counter) = DD_mod;
                    supercell.kind.emplace_back(kind_now);
                    if (spin_input.lspin) magmom_tmp.emplace_back(spin_input.magmom[iat]);
                    ++counter;
                }
            }
        }
    }

    supercell.x_cartesian = supercell.x_fractional * supercell.lattice_vector.transpose();

    spin_super.lspin = spin_input.lspin;
    spin_super.noncollinear = spin_input.noncollinear;
    spin_super.time_reversal_symm = spin_input.time_reversal_symm;
    spin_super.magmom.clear();
    spin_super.magmom.shrink_to_fit();
    if (spin_super.lspin) {
        std::copy(magmom_tmp.begin(),
                  magmom_tmp.end(),
                  std::back_inserter(spin_super.magmom));
    }
}

const Cell &System::get_supercell() const
{
    return supercell;
}

const Cell &System::get_primcell() const
{
    return primcell;
}

const Cell &System::get_inputcell() const
{
    return inputcell;
}

double ***System::get_x_image() const
{
    return x_image;
}

int *System::get_exist_image() const
{
    return exist_image;
}

void System::set_periodicity(const int is_periodic_in[3])
{
    if (!is_periodic) {
        // This should be already allocated though.
        allocate(is_periodic, 3);
    }
    for (unsigned int i = 0; i < 3; i++) {
        is_periodic[i] = is_periodic_in[i];
    }
}

int *System::get_periodicity() const
{
    return is_periodic;
}

void System::set_kdname(const std::string *kdname_in)
{
    //TODO: modify below
    const auto nkd = inputcell.number_of_elems;

    if (kdname) {
        deallocate(kdname);
    }
    allocate(kdname, nkd);
    for (size_t i = 0; i < nkd; ++i) {
        kdname[i] = kdname_in[i];
    }
}

std::string *System::get_kdname() const
{
    return kdname;
}

void System::set_reciprocal_latt(const Eigen::Matrix3d &lavec_in,
                                 Eigen::Matrix3d &rlavec_out) const
{
    // Compute reciprocal lattice vectors
    const auto det = lavec_in.determinant();

    if (std::abs(det) < eps12) {
        exit("set_reciprocal_latt", "Lattice Vector is singular");
    }
    const auto factor = 2.0 * pi;
    rlavec_out = factor * lavec_in.inverse();
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


void System::set_default_variables()
{
    kdname = nullptr;

    supercell.number_of_atoms = 0;
    supercell.number_of_elems = 0;

    allocate(is_periodic, 3);
    is_periodic[0] = 1;
    is_periodic[1] = 1;
    is_periodic[2] = 1;

    x_image = nullptr;
    exist_image = nullptr;
    str_magmom = "";

    spin_input.lspin = false;
    spin_input.noncollinear = 0;
    spin_input.time_reversal_symm = 1;
    transmat_to_super = Eigen::Matrix3d::Identity();
    transmat_to_prim = Eigen::Matrix3d::Identity();
}

void System::deallocate_variables()
{
    if (kdname) {
        deallocate(kdname);
    }
    if (x_image) {
        deallocate(x_image);
    }
    if (is_periodic) {
        deallocate(is_periodic);
    }
    if (exist_image) {
        deallocate(exist_image);
    }
}

void System::set_spin_variables(const size_t nat_in,
                                const bool lspin_in,
                                const int noncol_in,
                                const int trev_sym_in,
                                const double (*magmom_in)[3])
{
    spin_input.lspin = lspin_in;
    spin_input.noncollinear = noncol_in;
    spin_input.time_reversal_symm = trev_sym_in;
    spin_input.magmom.clear();

    std::vector<double> vec(3);
    for (size_t i = 0; i < nat_in; ++i) {
        for (auto j = 0; j < 3; ++j) {
            vec[j] = magmom_in[i][j];
        }
        spin_input.magmom.push_back(vec);
    }
}

const Spin &System::get_spin() const
{
    return spin_input;
}

void System::set_str_magmom(std::string str_magmom_in)
{
    str_magmom = str_magmom_in;
}

const std::string &System::get_str_magmom() const
{
    return str_magmom;
}

const std::vector<std::vector<unsigned int>> &System::get_atomtype_group() const
{
    return atomtype_group;
}


void System::set_atomtype_group()
{
    // In the case of collinear calculation, spin moments are considered as scalar
    // variables. Therefore, the same elements with different magnetic moments are
    // considered as different types. In noncollinear calculations,
    // magnetic moments are not considered in this stage. They will be treated
    // separately in symmetry.cpp where spin moments will be rotated and flipped
    // using time-reversal symmetry.

    unsigned int i;
    AtomType type_tmp{};
    std::set<AtomType> set_type;
    set_type.clear();

    for (i = 0; i < supercell.number_of_atoms; ++i) {
        type_tmp.element = supercell.kind[i];

        if (spin_input.noncollinear == 0) {
            type_tmp.magmom = spin_input.magmom[i][2];
        } else {
            type_tmp.magmom = 0.0;
        }
        set_type.insert(type_tmp);
    }

    const auto natomtypes = set_type.size();
    atomtype_group.resize(natomtypes);

    for (i = 0; i < supercell.number_of_atoms; ++i) {
        int count = 0;
        for (auto it: set_type) {
            if (spin_input.noncollinear) {
                if (supercell.kind[i] == it.element) {
                    atomtype_group[count].push_back(i);
                }
            } else {
                if ((supercell.kind[i] == it.element)
                    && (std::abs(spin_input.magmom[i][2] - it.magmom) < eps6)) {
                    atomtype_group[count].push_back(i);
                }
            }
            ++count;
        }
    }
    set_type.clear();
}

void System::set_transformation_matrices(const double transmat_to_super_in[3][3],
                                         const double transmat_to_prim_in[3][3])
{
    for (auto i = 0; i < 3; ++i) {
        for (auto j = 0; j < 3; ++j) {
            transmat_to_super(i, j) = transmat_to_super_in[i][j];
            transmat_to_prim(i, j) = transmat_to_prim_in[i][j];
        }
    }
}

void System::generate_coordinate_of_periodic_images()
{
    // Generate Cartesian coordinates of atoms in the neighboring 27 supercells

    unsigned int i;
    int ia, ja, ka;

    const auto nat = supercell.number_of_atoms;
    const Eigen::MatrixXd xf_in = supercell.x_fractional;

    Eigen::MatrixXd x_tmp(nat, 3);

    auto icell = 0;

    x_tmp = xf_in * supercell.lattice_vector.transpose();

    for (i = 0; i < nat; ++i) {
        for (unsigned int j = 0; j < 3; ++j) {
            x_image[0][i][j] = x_tmp(i, j);
        }
    }

    // Convert to Cartesian coordinate
    for (ia = -1; ia <= 1; ++ia) {
        for (ja = -1; ja <= 1; ++ja) {
            for (ka = -1; ka <= 1; ++ka) {
                if (ia == 0 && ja == 0 && ka == 0) continue;
                ++icell;
                x_tmp = supercell.x_fractional;
                for (i = 0; i < nat; ++i) {
                    x_tmp(i, 0) += static_cast<double>(ia);
                    x_tmp(i, 1) += static_cast<double>(ja);
                    x_tmp(i, 2) += static_cast<double>(ka);
                }
                x_tmp = x_tmp * supercell.lattice_vector.transpose();
                for (i = 0; i < nat; ++i) {
                    x_image[icell][i][0] = x_tmp(i, 0);
                    x_image[icell][i][1] = x_tmp(i, 1);
                    x_image[icell][i][2] = x_tmp(i, 2);
                }
            }
        }
    }

    icell = 0;
    exist_image[0] = 1;

    for (ia = -1; ia <= 1; ++ia) {
        for (ja = -1; ja <= 1; ++ja) {
            for (ka = -1; ka <= 1; ++ka) {
                if (ia == 0 && ja == 0 && ka == 0) continue;
                ++icell;
                // When periodic flag is zero along an axis,
                // periodic images along that axis cannot be considered.
                if (((std::abs(ia) == 1) && (is_periodic[0] == 0)) ||
                    ((std::abs(ja) == 1) && (is_periodic[1] == 0)) ||
                    ((std::abs(ka) == 1) && (is_periodic[2] == 0))) {
                    exist_image[icell] = 0;
                } else {
                    exist_image[icell] = 1;
                }
            }
        }
    }
}


void System::print_structure_stdout(const int verbosity)
{
    using namespace std;
    size_t i;

    cout << " ===================\n";
    cout << "  CRYSTAL STRUCTURE \n";
    cout << " ===================\n\n";

    cout.setf(ios::scientific);

    auto cell = get_inputcell();

    cout << "  ++++++++++++\n";
    cout << "   Input Cell \n";
    cout << "  ++++++++++++\n\n";

    cout << "   Lattice Vector (bohr)" << endl;
    cout << setw(16) << cell.lattice_vector(0, 0);
    cout << setw(15) << cell.lattice_vector(1, 0);
    cout << setw(15) << cell.lattice_vector(2, 0);
    cout << " : a1" << endl;

    cout << setw(16) << cell.lattice_vector(0, 1);
    cout << setw(15) << cell.lattice_vector(1, 1);
    cout << setw(15) << cell.lattice_vector(2, 1);
    cout << " : a2" << endl;

    cout << setw(16) << cell.lattice_vector(0, 2);
    cout << setw(15) << cell.lattice_vector(1, 2);
    cout << setw(15) << cell.lattice_vector(2, 2);
    cout << " : a3" << endl;
    cout << endl;

    cout << "   Number of atoms : " << cell.number_of_atoms << "\n\n";

    if (verbosity > 1) {
        cout << "   Cell volume = " << cell.volume << " (bohr^3)" << endl << endl;

        cout << "   Reciprocal Lattice Vector (1/bohr)" << std::endl;
        cout << setw(16) << cell.reciprocal_lattice_vector(0, 0);
        cout << setw(15) << cell.reciprocal_lattice_vector(0, 1);
        cout << setw(15) << cell.reciprocal_lattice_vector(0, 2);
        cout << " : b1" << endl;

        cout << setw(16) << cell.reciprocal_lattice_vector(1, 0);
        cout << setw(15) << cell.reciprocal_lattice_vector(1, 1);
        cout << setw(15) << cell.reciprocal_lattice_vector(1, 2);
        cout << " : b2" << endl;

        cout << setw(16) << cell.reciprocal_lattice_vector(2, 0);
        cout << setw(15) << cell.reciprocal_lattice_vector(2, 1);
        cout << setw(15) << cell.reciprocal_lattice_vector(2, 2);
        cout << " : b3" << endl;
        cout << endl;
    }

    cout << "   Atomic species:" << endl;
    for (i = 0; i < cell.number_of_elems; ++i) {
        cout << setw(6) << i + 1 << setw(5) << kdname[i] << endl;
    }

    if (verbosity > 1) {
        cout << '\n';
        cout << "   Atomic positions in fractional basis and atomic species" << endl;
        for (i = 0; i < cell.number_of_atoms; ++i) {
            cout << setw(6) << i + 1;
            cout << setw(15) << cell.x_fractional(i, 0);
            cout << setw(15) << cell.x_fractional(i, 1);
            cout << setw(15) << cell.x_fractional(i, 2);
            cout << setw(5) << kdname[cell.kind[i] - 1] << endl;
        }
    }

    cout << endl << endl;

    cell = get_primcell();

    cout << "  ++++++++++++++++\n";
    cout << "   Primitive Cell \n";
    cout << "  ++++++++++++++++\n\n";

    cout << "   Lattice Vector (bohr)" << endl;
    cout << setw(16) << cell.lattice_vector(0, 0);
    cout << setw(15) << cell.lattice_vector(1, 0);
    cout << setw(15) << cell.lattice_vector(2, 0);
    cout << " : a1" << endl;

    cout << setw(16) << cell.lattice_vector(0, 1);
    cout << setw(15) << cell.lattice_vector(1, 1);
    cout << setw(15) << cell.lattice_vector(2, 1);
    cout << " : a2" << endl;

    cout << setw(16) << cell.lattice_vector(0, 2);
    cout << setw(15) << cell.lattice_vector(1, 2);
    cout << setw(15) << cell.lattice_vector(2, 2);
    cout << " : a3" << endl;
    cout << endl;

    const auto nat_prim = cell.number_of_atoms;

    cout << "   Number of atoms : " << cell.number_of_atoms << "\n\n";

    if (verbosity > 1) {
        cout << "   Cell volume = " << cell.volume << " (bohr^3)" << endl << endl;

        cout << "   Reciprocal Lattice Vector (1/bohr)" << std::endl;
        cout << setw(16) << cell.reciprocal_lattice_vector(0, 0);
        cout << setw(15) << cell.reciprocal_lattice_vector(0, 1);
        cout << setw(15) << cell.reciprocal_lattice_vector(0, 2);
        cout << " : b1" << endl;

        cout << setw(16) << cell.reciprocal_lattice_vector(1, 0);
        cout << setw(15) << cell.reciprocal_lattice_vector(1, 1);
        cout << setw(15) << cell.reciprocal_lattice_vector(1, 2);
        cout << " : b2" << endl;

        cout << setw(16) << cell.reciprocal_lattice_vector(2, 0);
        cout << setw(15) << cell.reciprocal_lattice_vector(2, 1);
        cout << setw(15) << cell.reciprocal_lattice_vector(2, 2);
        cout << " : b3" << endl;
        cout << endl;
    }

    cout << "   Atomic positions in fractional basis and atomic species" << endl;
    for (i = 0; i < cell.number_of_atoms; ++i) {
        cout << setw(6) << i + 1;
        cout << setw(15) << cell.x_fractional(i, 0);
        cout << setw(15) << cell.x_fractional(i, 1);
        cout << setw(15) << cell.x_fractional(i, 2);
        cout << setw(5) << kdname[cell.kind[i] - 1] << endl;
    }
    cout << endl << endl;

    cell = get_supercell();
    cout << "  ++++++++++++\n";
    cout << "   Super Cell \n";
    cout << "  ++++++++++++\n\n";

    cout << "   Lattice Vector (bohr)" << endl;
    cout << setw(16) << cell.lattice_vector(0, 0);
    cout << setw(15) << cell.lattice_vector(1, 0);
    cout << setw(15) << cell.lattice_vector(2, 0);
    cout << " : a1" << endl;

    cout << setw(16) << cell.lattice_vector(0, 1);
    cout << setw(15) << cell.lattice_vector(1, 1);
    cout << setw(15) << cell.lattice_vector(2, 1);
    cout << " : a2" << endl;

    cout << setw(16) << cell.lattice_vector(0, 2);
    cout << setw(15) << cell.lattice_vector(1, 2);
    cout << setw(15) << cell.lattice_vector(2, 2);
    cout << " : a3" << endl;
    cout << endl;

    cout << "   Number of atoms : " << cell.number_of_atoms << "\n\n";
    cout << "   Supercell contains " << std::setw(5) << cell.number_of_atoms / nat_prim
         << " primitive cells\n\n";

    if (verbosity > 1) {
        cout << "   Cell volume = " << cell.volume << " (bohr^3)" << endl << endl;

        cout << "   Reciprocal Lattice Vector (1/bohr)" << std::endl;
        cout << setw(16) << cell.reciprocal_lattice_vector(0, 0);
        cout << setw(15) << cell.reciprocal_lattice_vector(0, 1);
        cout << setw(15) << cell.reciprocal_lattice_vector(0, 2);
        cout << " : b1" << endl;

        cout << setw(16) << cell.reciprocal_lattice_vector(1, 0);
        cout << setw(15) << cell.reciprocal_lattice_vector(1, 1);
        cout << setw(15) << cell.reciprocal_lattice_vector(1, 2);
        cout << " : b2" << endl;

        cout << setw(16) << cell.reciprocal_lattice_vector(2, 0);
        cout << setw(15) << cell.reciprocal_lattice_vector(2, 1);
        cout << setw(15) << cell.reciprocal_lattice_vector(2, 2);
        cout << " : b3" << endl;
        cout << endl;

        cout << "   Atomic positions in fractional basis and atomic species" << endl;
        for (i = 0; i < cell.number_of_atoms; ++i) {
            cout << setw(6) << i + 1;
            cout << setw(15) << cell.x_fractional(i, 0);
            cout << setw(15) << cell.x_fractional(i, 1);
            cout << setw(15) << cell.x_fractional(i, 2);
            cout << setw(5) << kdname[cell.kind[i] - 1] << endl;
        }
    }
    cout << endl << endl;

    cout.unsetf(ios::scientific);
}


void System::print_magmom_stdout() const
{
    using namespace std;

    cout << " ====================\n";
    cout << "  MAGNETIC STRUCTURE \n";
    cout << " ====================\n\n";

    cout << "  MAGMOM is given.\n"
            "  The magnetic moments of each atom in the primitive cell are as follows:" << endl;
    for (size_t i = 0; i < primcell.number_of_atoms; ++i) {
        cout << setw(6) << i + 1;
        cout << setw(5) << spin_prim.magmom[i][0];
        cout << setw(5) << spin_prim.magmom[i][1];
        cout << setw(5) << spin_prim.magmom[i][2];
        cout << endl;
    }
    cout << endl;
    if (spin_prim.noncollinear == 0) {
        cout << "  NONCOLLINEAR = 0: magnetic moments are considered as scalar variables." << endl;
    } else if (spin_prim.noncollinear == 1) {
        cout << "  NONCOLLINEAR = 1: magnetic moments are considered as vector variables." << endl;
        if (spin_prim.time_reversal_symm) {
            cout << "  TREVSYM = 1: Time-reversal symmetry will be considered for generating magnetic space group"
                 << endl;
        } else {
            cout <<
                 "  TREVSYM = 0: Time-reversal symmetry will NOT be considered for generating magnetic space group"
                 << endl;
        }
    }
    cout << endl << endl;
}

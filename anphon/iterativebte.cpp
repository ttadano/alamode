#include "mpi_common.h"
#include "conductivity.h"
#include "iterativebte.h"
#include "constants.h"
#include "dynamical.h"
#include "error.h"
#include "integration.h"
#include "parsephon.h"
#include "isotope.h"
#include "kpoint.h"
#include "mathfunctions.h"
#include "memory.h"
#include "phonon_dos.h"
#include "thermodynamics.h"
#include "phonon_velocity.h"
#include "anharmonic_core.h"
#include "system.h"
#include "write_phonons.h"
#include "symmetry_core.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <iterator>

//mpic++ -o3 -std=c++11 -I../include -I/Users/wenhao/mylib/include -I/Users/wenhao/mylib/spg/include -I/Users/wenhao/mylib/fftw/3.3.9/include -c iterativebte.cpp
// 
// DONE: test tetrahedron method
// DONE: test with isotope
// TODO: test with more grid
// TODO: calculation from a restart
// TODO: write .result file

using namespace PHON_NS;

Iterativebte::Iterativebte(PHON *phon) : Pointers(phon)
{
    set_default_variables();
}

Iterativebte::~Iterativebte()
{
    deallocate_variables();
}

void Iterativebte::set_default_variables()
{
    // public
    do_iterative = true;
    direct_solution = false;
    Temperature = nullptr;
    ntemp = 0;
    max_cycle = 20;
    convergence_criteria = 0.005;
    kappa = nullptr;
    use_triplet_symmetry = true;
    sym_permutation = true;

    // private
    vel = nullptr;
    dFold = nullptr;
    dFnew = nullptr;
    L_absorb = nullptr;
    L_emitt = nullptr;
    damping4 = nullptr;
}

void Iterativebte::deallocate_variables()
{
    if (Temperature) { deallocate(Temperature); }
    if (kappa) { deallocate(kappa); }
    if (vel) { deallocate(vel); }
    if (dFold) { deallocate(dFold); }
    if (dFnew) { deallocate(dFnew); }
    if (L_absorb) { deallocate(L_absorb); }
    if (L_emitt) { deallocate(L_emitt); }
    if (damping4) { deallocate(damping4); }
}

void Iterativebte::setup_iterative()
{
    nktot = dos->kmesh_dos->nk;
    ns = dynamical->neval;
    ns2 = ns * ns;

    MPI_Bcast(&max_cycle, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&convergence_criteria, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    sym_permutation = false;
    use_triplet_symmetry = true;

    // Temperature in K
    ntemp = static_cast<unsigned int>((system->Tmax - system->Tmin) / system->dT) + 1;
    allocate(Temperature, ntemp);

    for (auto i = 0; i < ntemp; ++i) {
        Temperature[i] = system->Tmin + static_cast<double>(i) * system->dT;
    }

    // calculate vel

    allocate(vel, nktot, ns, 3);            // velocity of the full grid
    phonon_velocity->get_phonon_group_velocity_mesh_mpi(*dos->kmesh_dos,
                                                        system->lavec_p,
                                                        fcs_phonon->fc2_ext,
                                                        vel); //this will gather to rank0 process

    MPI_Bcast(&vel[0][0][0], nktot * ns * 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    //std::cout << mympi->my_rank << " " << typeid(vel[0][0][0]).name() << std::endl;
    /*
    if (mympi->my_rank == 1) {
        std::cout << "here" << std::endl;
        for (auto ik = 0; ik < nktot; ++ik) {
            for (auto is = 0; is < ns; ++is) {
                std::cout << std::setw(12) << dynamical->eval_phonon[ik][is] * 1.1e5 << " ";
            }
            std::cout << std::endl;
        } 
    }*/

    allocate(kappa, ntemp, 3, 3);

    // distribute q point among the processors
    auto nk_ir = dos->kmesh_dos->nk_irred;

    nk_l.clear();
    for (auto i = 0; i < nk_ir; ++i) {
        if (i % mympi->nprocs == mympi->my_rank) nk_l.push_back(i);
    }

    nklocal = nk_l.size();

    if (mympi->my_rank == 0) {
        std::cout << std::endl;
        std::cout << " distribute q point ... " << std::endl;
        std::cout << " number of q point each process: " << std::setw(5) << nklocal << std::endl;
        std::cout << std::endl;
        std::cout << " Generating all k pairs ... " << std::endl;
    }

    get_triplets();

    if (mympi->my_rank == 0) {
        std::cout << "     DONE ! " << std::endl;
    }
    //setup_kpindex();

    write_result();
}

void Iterativebte::get_triplets()
{
    localnk_triplets_emitt.clear();  // pairs k3 = k1 - k2 ( -k1 + k2 + k3 = G )
    localnk_triplets_absorb.clear(); // pairs k3 = k1 + k2 (  k1 + k2 + k3 = G )

    int counter = 0;
    int counter2 = 0;

    for (unsigned int i = 0; i < nklocal; ++i) {

        auto ik = nk_l[i];
        std::vector<KsListGroup> triplet;
        std::vector<KsListGroup> triplet2;

        // k3 = k1 - k2
        dos->kmesh_dos->get_unique_triplet_k(ik, symmetry->SymmList,
                                             use_triplet_symmetry,
                                             sym_permutation, triplet);

        // k3 = - (k1 + k2)
        dos->kmesh_dos->get_unique_triplet_k(ik, symmetry->SymmList,
                                             use_triplet_symmetry,
                                             sym_permutation, triplet2, 1);

        counter += triplet.size();
        counter2 += triplet2.size();

        localnk_triplets_emitt.push_back(triplet);
        localnk_triplets_absorb.push_back(triplet2);
    }

    kplength_emitt = counter;   // remember number of unique pairs
    kplength_absorb = counter2;

}

void Iterativebte::do_iterativebte()
{
    if (direct_solution) {
        direct_solver();
    } else {
        setup_L();
        iterative_solver();
    }
}

void Iterativebte::setup_L()
{
    if (mympi->my_rank == 0) {
        std::cout << std::endl;
        std::cout << " Calculate once for the transition probability L(absorb) and L(emitt) ..." << std::endl;
        std::cout << "     size of L (MB) (approx.) = " << memsize_in_MB(sizeof(double), kplength_absorb, ns, ns2)
                  << std::endl;
        if (integration->ismear == 0 || integration->ismear == 1) {
            std::cout << "     smearing will be used for the delta function" << std::endl;
        } else if (integration->ismear == 2) {
            std::cout << "     adaptive methods will be used for the delta function" << std::endl;
        } else if (integration->ismear == -1) {
            std::cout << "     tetrahedron method will be used for the delta function" << std::endl;
        } else {
            exit("calc_L", "ismear other than -1, 0, 1, 2 are not supported");
        }
    }
    if (integration->ismear >= 0) {
        setup_L_smear();
    } else if (integration->ismear == -1) {
        setup_L_tetra();
    }
    if (mympi->my_rank == 0) {
        std::cout << "     DONE !" << std::endl;
    }
}

void Iterativebte::setup_L_smear()
{
    // we calculate V for all pairs L+(local_nk*eachpair,ns,ns2) and L-

    allocate(L_absorb, kplength_absorb, ns, ns2);
    allocate(L_emitt, kplength_emitt, ns, ns2);

    unsigned int arr[3];
    int k1, k2, k3, k1_minus;
    int s1, s2, s3;
    int ib;
    double omega1, omega2, omega3;

    double v3_tmp;

    unsigned int counter;
    double delta = 0;

    auto epsilon = integration->epsilon;

    // we generate the counters

    const auto omega_tmp = dos->dymat_dos->get_eigenvalues();
    const auto evec_tmp = dos->dymat_dos->get_eigenvectors();

    // emitt
    counter = 0;
    for (auto ik = 0; ik < nklocal; ++ik) {

        auto tmpk = nk_l[ik];
        k1 = dos->kmesh_dos->kpoint_irred_all[tmpk][0].knum;    // k index in full grid

        // emitt k1 -> k2 + k3 
        // V(-q1, q2, q3) delta(w1 - w2 - w3)
        for (auto j = 0; j < localnk_triplets_emitt[ik].size(); ++j) {

            auto pair = localnk_triplets_emitt[ik][j];

            k2 = pair.group[0].ks[0];
            k3 = pair.group[0].ks[1];
            //counter = ikp_emitt[ik][j];

            for (s1 = 0; s1 < ns; ++s1) {
                arr[0] = dos->kmesh_dos->kindex_minus_xk[k1] * ns + s1;
                omega1 = omega_tmp[k1][s1];

                for (ib = 0; ib < ns2; ++ib) {
                    s2 = ib / ns;
                    s3 = ib % ns;

                    arr[1] = k2 * ns + s2;
                    arr[2] = k3 * ns + s3;
                    omega2 = omega_tmp[k2][s2];
                    omega3 = omega_tmp[k3][s3];

                    if (integration->ismear == 0) {
                        delta = delta_lorentz(omega1 - omega2 - omega3, epsilon);
                    } else if (integration->ismear == 1) {
                        delta = delta_gauss(omega1 - omega2 - omega3, epsilon);
                    } else if (integration->ismear == 2) {
                        double epsilon2[2];
                        integration->adaptive_sigma->get_sigma(k2, s2, k3, s3, epsilon2);
                        //integration->adaptive_smearing(k2, s2, k3, s3, epsilon2);
                        //sum_smear += epsilon2[0] + epsilon2[1];
                        delta = delta_gauss(omega1 - omega2 - omega3, epsilon2[0]);
                    }

                    v3_tmp = std::norm(anharmonic_core->V3(arr,
                                                           dos->kmesh_dos->xk,
                                                           omega_tmp,
                                                           evec_tmp));

                    L_emitt[counter][s1][ib] = (pi / 4.0) * v3_tmp * delta / static_cast<double>(nktot);
                    //sum1 += L_emitt[counter][s1][ib];
                }
            }
            counter += 1;
        }
    }
    if (counter != kplength_emitt) {
        exit("setup_L", "Emitt: pair length not equal!");
    }

    counter = 0;
    for (auto ik = 0; ik < nklocal; ++ik) {

        auto tmpk = nk_l[ik];
        k1 = dos->kmesh_dos->kpoint_irred_all[tmpk][0].knum;    // k index in full grid

        // absorption k1 + k2 -> -k3
        // V(q1, q2, q3) since k3 = - (k1 + k2)
        for (auto j = 0; j < localnk_triplets_absorb[ik].size(); ++j) {

            auto pair = localnk_triplets_absorb[ik][j];
            //counter = ikp_absorb[ik][j];

            k2 = pair.group[0].ks[0];
            k3 = pair.group[0].ks[1];

            for (s1 = 0; s1 < ns; ++s1) {
                arr[0] = k1 * ns + s1;
                omega1 = omega_tmp[k1][s1];

                for (ib = 0; ib < ns2; ++ib) {
                    s2 = ib / ns;
                    s3 = ib % ns;

                    arr[1] = k2 * ns + s2;
                    arr[2] = k3 * ns + s3;
                    omega2 = omega_tmp[k2][s2];
                    omega3 = omega_tmp[k3][s3];

                    if (integration->ismear == 0) {
                        delta = delta_lorentz(omega1 + omega2 - omega3, epsilon);
                    } else if (integration->ismear == 1) {
                        delta = delta_gauss(omega1 + omega2 - omega3, epsilon);
                    } else if (integration->ismear == 2) {
                        double epsilon2[2];
                        integration->adaptive_sigma->get_sigma(k2, s2, k3, s3, epsilon2);
                        //integration->adaptive_smearing(k2, s2, k3, s3, epsilon2);
                        //sum_smear += epsilon2[0] + epsilon2[1];
                        delta = delta_gauss(omega1 + omega2 - omega3, epsilon2[1]);
                    }

                    v3_tmp = std::norm(anharmonic_core->V3(arr,
                                                           dos->kmesh_dos->xk,
                                                           omega_tmp,
                                                           evec_tmp));

                    L_absorb[counter][s1][ib] = (pi / 4.0) * v3_tmp * delta / static_cast<double>(nktot);
                    //sum2 += L_absorb[counter][s1][ib];
                }
            }
            counter += 1;
        }
    }

    if (counter != kplength_absorb) {
        exit("setup_L", "absorb: pair length not equal!");
    }

    //if (mympi->my_rank == 0) {
    //    std::cout << "  DONE !" << std::endl;
    //}
}

void Iterativebte::setup_L_tetra()
{
    // we calculate V for all pairs L+(local_nk*eachpair,ns,ns2) and L-

    allocate(L_absorb, kplength_absorb, ns, ns2);
    allocate(L_emitt, kplength_emitt, ns, ns2);

    unsigned int arr[3];
    int k1, k2, k3, k1_minus;
    int s1, s2, s3;
    int ib;
    double omega1, omega2, omega3;

    double v3_tmp;
    double xk_tmp[3];

    unsigned int counter;
    double delta = 0;

    auto epsilon = integration->epsilon;

    unsigned int *kmap_identity;
    allocate(kmap_identity, nktot);
    for (auto i = 0; i < nktot; ++i) kmap_identity[i] = i;

    double *energy_tmp;
    double *weight_tetra;
    allocate(energy_tmp, nktot);
    allocate(weight_tetra, nktot);

    // emitt
    std::vector<std::vector<int>> ikp_emitt;
    ikp_emitt.clear();
    int cnt = 0;
    for (auto ik = 0; ik < nklocal; ++ik) {
        std::vector<int> counterk;
        counterk.clear();
        for (auto j = 0; j < localnk_triplets_emitt[ik].size(); ++j) {
            counterk.push_back(cnt);
            cnt += 1;
        }
        ikp_emitt.push_back(counterk);
    }
    // absorb
    std::vector<std::vector<int>> ikp_absorb;
    ikp_absorb.clear();
    cnt = 0;
    for (auto ik = 0; ik < nklocal; ++ik) {
        std::vector<int> counterk;
        counterk.clear();
        for (auto j = 0; j < localnk_triplets_absorb[ik].size(); ++j) {
            counterk.push_back(cnt);
            cnt += 1;
        }
        ikp_absorb.push_back(counterk);
    }

    const auto omega_tmp = dos->dymat_dos->get_eigenvalues();
    const auto evec_tmp = dos->dymat_dos->get_eigenvectors();

    for (auto ik = 0; ik < nklocal; ++ik) {

        auto tmpk = nk_l[ik];
        k1 = dos->kmesh_dos->kpoint_irred_all[tmpk][0].knum;    // k index in full grid

        // emission: calc delta function w1 - w2 - w3, with k3 = k1 - k2
        for (s1 = 0; s1 < ns; ++s1) {

            omega1 = omega_tmp[k1][s1];

            for (ib = 0; ib < ns2; ++ib) {
                s2 = ib / ns;
                s3 = ib % ns;

                for (k2 = 0; k2 < nktot; k2++) {
                    // k3 = k1 - k2
                    for (auto i = 0; i < 3; ++i) {
                        xk_tmp[i] = dos->kmesh_dos->xk[k1][i] - dos->kmesh_dos->xk[k2][i];
                    }

                    k3 = dos->kmesh_dos->get_knum(xk_tmp);

                    omega2 = omega_tmp[k2][s2];
                    omega3 = omega_tmp[k3][s3];

                    energy_tmp[k2] = omega2 + omega3;
                }
                integration->calc_weight_tetrahedron(nktot,
                                                     kmap_identity,
                                                     energy_tmp,
                                                     omega1,
                                                     dos->tetra_nodes_dos->get_ntetra(),
                                                     dos->tetra_nodes_dos->get_tetras(),
                                                     weight_tetra);

                // emitt k1 -> k2 + k3 
                // V(-q1, q2, q3) delta(w1 - w2 - w3)
                for (auto j = 0; j < localnk_triplets_emitt[ik].size(); ++j) {

                    auto pair = localnk_triplets_emitt[ik][j];
                    auto counter = ikp_emitt[ik][j];

                    k2 = pair.group[0].ks[0];
                    k3 = pair.group[0].ks[1];

                    arr[0] = dos->kmesh_dos->kindex_minus_xk[k1] * ns + s1;
                    arr[1] = k2 * ns + s2;
                    arr[2] = k3 * ns + s3;
                    delta = weight_tetra[k2];
                    v3_tmp = std::norm(anharmonic_core->V3(arr,
                                                           dos->kmesh_dos->xk,
                                                           omega_tmp,
                                                           evec_tmp));

                    L_emitt[counter][s1][ib] = (pi / 4.0) * v3_tmp * delta;
                }

                // absorb
                for (k2 = 0; k2 < nktot; k2++) {

                    for (auto i = 0; i < 3; ++i) {
                        xk_tmp[i] = dos->kmesh_dos->xk[k1][i] + dos->kmesh_dos->xk[k2][i];
                    }
                    k3 = dos->kmesh_dos->get_knum(xk_tmp);

                    omega2 = omega_tmp[k2][s2];
                    omega3 = omega_tmp[k3][s3];

                    energy_tmp[k2] = -omega2 + omega3;
                }
                integration->calc_weight_tetrahedron(nktot,
                                                     kmap_identity,
                                                     energy_tmp,
                                                     omega1,
                                                     dos->tetra_nodes_dos->get_ntetra(),
                                                     dos->tetra_nodes_dos->get_tetras(),
                                                     weight_tetra);

                // absorption k1 + k2 -> -k3
                // V(q1, q2, q3) since k3 = - (k1 + k2)
                for (auto j = 0; j < localnk_triplets_absorb[ik].size(); ++j) {

                    auto pair = localnk_triplets_absorb[ik][j];
                    auto counter = ikp_absorb[ik][j];

                    k2 = pair.group[0].ks[0];
                    k3 = pair.group[0].ks[1];

                    arr[0] = k1 * ns + s1;
                    arr[1] = k2 * ns + s2;
                    arr[2] = k3 * ns + s3;
                    delta = weight_tetra[k2];
                    v3_tmp = std::norm(anharmonic_core->V3(arr,
                                                           dos->kmesh_dos->xk,
                                                           omega_tmp,
                                                           evec_tmp));

                    L_absorb[counter][s1][ib] = (pi / 4.0) * v3_tmp * delta;
                }
            } // ib
        } // s1    
    } // ik

    //if (mympi->my_rank == 0) {
    //    std::cout << "  DONE !" << std::endl;
    //}
    deallocate(kmap_identity);
    deallocate(energy_tmp);
    deallocate(weight_tetra);
}

void Iterativebte::calc_Q_from_L(double **&n, double **&q1)
{
    int s1, s2, s3;
    double ph1;
    double n1, n2, n3;

    double **Qemit;
    double **Qabsorb;
    allocate(Qemit, nklocal, ns);
    allocate(Qabsorb, nklocal, ns);

    for (auto ik = 0; ik < nklocal; ++ik) {
        for (s1 = 0; s1 < ns; ++s1) {
            Qemit[ik][s1] = 0.0;
            Qabsorb[ik][s1] = 0.0;
        }
    }

    unsigned int counter;

    // emit
    counter = 0;
    for (auto ik = 0; ik < nklocal; ++ik) {

        auto tmpk = nk_l[ik];
        const int k1 = dos->kmesh_dos->kpoint_irred_all[tmpk][0].knum;

        for (auto j = 0; j < localnk_triplets_emitt[ik].size(); ++j) {

            auto pair = localnk_triplets_emitt[ik][j];
            auto multi = static_cast<double>(pair.group.size());
            auto k2 = pair.group[0].ks[0];
            auto k3 = pair.group[0].ks[1];

            for (s1 = 0; s1 < ns; ++s1) {
                n1 = n[k1][s1];

                for (int ib = 0; ib < ns2; ++ib) {
                    s2 = ib / ns;
                    s3 = ib % ns;
                    n2 = n[k2][s2];
                    n3 = n[k3][s3];
                    Qemit[ik][s1] += 0.5 * (n1 * (n2 + 1.0) * (n3 + 1.0)) * L_emitt[counter][s1][ib] * multi;
                }
            }
            counter += 1;
        }

    }
    if (counter != kplength_emitt) {
        exit("setup_L", "Emitt: pair length not equal!");
    }

    // absorb k1 + k2 -> -k3
    counter = 0;
    for (auto ik = 0; ik < nklocal; ++ik) {

        auto tmpk = nk_l[ik];
        const int k1 = dos->kmesh_dos->kpoint_irred_all[tmpk][0].knum;

        for (auto j = 0; j < localnk_triplets_absorb[ik].size(); ++j) {

            auto pair = localnk_triplets_absorb[ik][j];
            auto multi = static_cast<double>(pair.group.size());
            auto k2 = pair.group[0].ks[0];
            auto k3 = pair.group[0].ks[1];

            for (s1 = 0; s1 < ns; ++s1) {
                n1 = n[k1][s1];

                for (int ib = 0; ib < ns2; ++ib) {
                    s2 = ib / ns;
                    s3 = ib % ns;
                    n2 = n[k2][s2];
                    n3 = n[k3][s3];
                    Qabsorb[ik][s1] += (n1 * n2 * (n3 + 1.0)) * L_absorb[counter][s1][ib] * multi;
                }
            }
            counter += 1;

        }
        //MPI_Barrier(MPI_COMM_WORLD);
        //exit("1","2");
    }

    if (counter != kplength_absorb) {
        exit("setup_L", "absorb: pair length not equal!");
    }

    for (auto ik = 0; ik < nklocal; ++ik) {
        for (s1 = 0; s1 < ns; ++s1) {
            q1[ik][s1] = Qemit[ik][s1] + Qabsorb[ik][s1];
            //if (mympi->my_rank == 1) {
            //    std::cout << ik << " " << s1 << "emit:" << std::setw(12) << std::scientific
            //             << std::setprecision(2) << Qemit[ik][s1] << std::endl;  
            //    std::cout << ik << " " << s1 << "abso:" << std::setw(12) << std::scientific
            //             << std::setprecision(2) << Qabsorb[ik][s1] << std::endl;  
            //}
        }
    }
    deallocate(Qemit);
    deallocate(Qabsorb);
}

void Iterativebte::calc_anharmonic_imagself4()
{
    // TODO: merge this duplicate function to conductivity class
    unsigned int i;
    unsigned int *nks_thread;

    // Distribute (k,s) to individual MPI threads
    allocate(damping4, ntemp, nklocal, ns);

    double *damping4_loc;
    allocate(damping4_loc, ntemp);

    if (mympi->my_rank == 0) {
        std::cout << std::endl;
        std::cout << " Computing 4-phonon scattering amplitude ... " << std::endl;
        std::cout << " WARNING: This is very expensive!! Please be patient." << std::endl;
    }

    for (auto ik = 0; ik < nklocal; ++ik) {
        auto tmpk = nk_l[ik];
        auto k1 = dos->kmesh_dos->kpoint_irred_all[tmpk][0].knum;
        for (auto s1 = 0; s1 < ns; ++s1) {
            auto omega = dos->dymat_dos->get_eigenvalues()[k1][s1];

            anharmonic_core->calc_damping4_smearing(ntemp,
                                                    Temperature,
                                                    omega,
                                                    tmpk,
                                                    s1,
                                                    dos->kmesh_dos,
                                                    dos->dymat_dos->get_eigenvalues(),
                                                    dos->dymat_dos->get_eigenvectors(),
                                                    damping4_loc);

            for (auto itemp = 0; itemp < ntemp; ++itemp) {
                damping4[itemp][ik][s1] = damping4_loc[itemp];
            }

        }
    }

    deallocate(damping4_loc);
}

void Iterativebte::iterative_solver()
{
    // calculate the Vs
    double mixing_factor = 0.75;

    double **Q;
    double **kappa_new;
    double **kappa_old;

    //bool converged;

    allocate(kappa_new, 3, 3);
    allocate(kappa_old, 3, 3);
    allocate(Q, nklocal, ns);

    allocate(dFold, nktot, ns, 3);
    allocate(dFnew, nktot, ns, 3);

    double **fb;
    double **dndt;
    allocate(dndt, nklocal, ns);
    allocate(fb, nktot, ns);

    double **isotope_damping_loc;
    if (isotope->include_isotope) {
        double **isotope_damping;
        allocate(isotope_damping, dos->kmesh_dos->nk_irred, ns);

        if (mympi->my_rank == 0) {
            for (auto ik = 0; ik < dos->kmesh_dos->nk_irred; ik++) {
                for (auto is = 0; is < ns; is++) {
                    isotope_damping[ik][is] = isotope->gamma_isotope[ik][is];
                }
            }
        }

        MPI_Bcast(&isotope_damping[0][0], dos->kmesh_dos->nk_irred * ns, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        allocate(isotope_damping_loc, nklocal, ns); // this is for reducing some memory usage
        for (auto ik = 0; ik < nklocal; ik++) {
            auto tmpk = nk_l[ik];
            for (auto is = 0; is < ns; is++) {
                isotope_damping_loc[ik][is] = isotope_damping[tmpk][is];
            }
        }

        deallocate(isotope_damping);
    }

    if (conductivity->fph_rta > 0) {
        calc_anharmonic_imagself4();
    }

    if (mympi->my_rank == 0) {
        std::cout << std::endl << " Iteration starts ..." << std::endl << std::endl;
    }


    // we solve iteratively for each temperature
    int ik, is, ix, iy;
    double n1, n2, n3;

    // generate index for, emitt
    std::vector<std::vector<int>> ikp_emitt;
    ikp_emitt.clear();
    int cnt = 0;
    for (ik = 0; ik < nklocal; ++ik) {
        std::vector<int> counterk;
        counterk.clear();
        for (auto j = 0; j < localnk_triplets_emitt[ik].size(); ++j) {
            counterk.push_back(cnt);
            cnt += 1;
        }
        ikp_emitt.push_back(counterk);
    }
    // absorb
    std::vector<std::vector<int>> ikp_absorb;
    ikp_absorb.clear();
    cnt = 0;
    for (ik = 0; ik < nklocal; ++ik) {
        std::vector<int> counterk;
        counterk.clear();
        for (auto j = 0; j < localnk_triplets_absorb[ik].size(); ++j) {
            counterk.push_back(cnt);
            cnt += 1;
        }
        ikp_absorb.push_back(counterk);
    }

    // start iteration
    double **Wks;
    allocate(Wks, ns, 3);

    for (auto itemp = 0; itemp < ntemp; ++itemp) {

        double beta = 1.0 / (thermodynamics->T_to_Ryd * Temperature[itemp]);

        calc_boson(itemp, fb, dndt);

        if (mympi->my_rank == 0) {
            std::cout << " Temperature step ..." << std::setw(10) << std::right
                      << std::fixed << std::setprecision(2) << Temperature[itemp] << " K" <<
                      "    -----------------------------" << std::endl;
            std::cout << "      Kappa [W/mK]        xx          xy          xz" <<
                      "          yx          yy          yz" <<
                      "          zx          zy          zz" << std::endl;
        }

        calc_Q_from_L(fb, Q);
        average_Q(Q);

        for (ik = 0; ik < nktot; ++ik) {
            for (is = 0; is < ns; ++is) {
                for (ix = 0; ix < 3; ++ix) {
                    dFold[ik][is][ix] = 0.0;
                }
            }
        }

        int s1, s2, s3;
        int k1, k2, k3, k3_minus;

        int generating_sym;
        int nsym = symmetry->SymmList.size();
        int isym;

        for (auto itr = 0; itr < max_cycle; ++itr) {

            if (mympi->my_rank == 0) {
                std::cout << "   -> iter " << std::setw(3) << itr << ": ";
            }

            // zero dFnew because we will do MPI_allreduce
            for (ik = 0; ik < nktot; ++ik) {
                for (is = 0; is < ns; ++is) {
                    for (ix = 0; ix < 3; ++ix) {
                        dFnew[ik][is][ix] = 0.0;
                    }
                }
            }

            for (ik = 0; ik < nklocal; ++ik) {

                auto tmpk = nk_l[ik];
                auto num_equivalent = dos->kmesh_dos->kpoint_irred_all[tmpk].size();
                auto kref = dos->kmesh_dos->kpoint_irred_all[tmpk][0].knum;

                for (auto ieq = 0; ieq < num_equivalent; ++ieq) {

                    k1 = dos->kmesh_dos->kpoint_irred_all[tmpk][ieq].knum;      // k1 will go through all points

                    generating_sym = -1;
                    for (isym = 0; isym < nsym; ++isym) {
                        auto krot = dos->kmesh_dos->knum_sym(kref, symmetry->SymmList[isym].rot);
                        if (k1 == krot) generating_sym = isym;
                    }
                    if (generating_sym == -1) {
                        exit("iterative solution", "cannot find all equivalent k");
                    }

                    // calculate W here
                    for (s1 = 0; s1 < ns; ++s1) {

                        for (ix = 0; ix < 3; ++ix) {
                            Wks[s1][ix] = 0.0;
                        }

                        // emitt k1 -> k2 + k3
                        for (auto j = 0; j < localnk_triplets_emitt[ik].size(); ++j) {

                            auto pair = localnk_triplets_emitt[ik][j];
                            int kp_index = ikp_emitt[ik][j];

                            for (auto ig = 0; ig < pair.group.size(); ig++) {

                                k2 = dos->kmesh_dos->knum_sym(pair.group[ig].ks[0],
                                                              symmetry->SymmList[generating_sym].rot);
                                k3 = dos->kmesh_dos->knum_sym(pair.group[ig].ks[1],
                                                              symmetry->SymmList[generating_sym].rot);

                                for (int ib = 0; ib < ns2; ++ib) {
                                    s2 = ib / ns;
                                    s3 = ib % ns;

                                    n1 = fb[k1][s1];
                                    n2 = fb[k2][s2];
                                    n3 = fb[k3][s3];
                                    for (ix = 0; ix < 3; ++ix) {
                                        Wks[s1][ix] -= 0.5 * (dFold[k2][s2][ix] + dFold[k3][s3][ix]) * n1 * (n2 + 1.0)
                                                       * (n3 + 1.0) * L_emitt[kp_index][s1][ib];
                                    }
                                }
                            }
                        }

                        // absorb k1 + k2 -> -k3
                        for (auto j = 0; j < localnk_triplets_absorb[ik].size(); ++j) {

                            auto pair = localnk_triplets_absorb[ik][j];
                            int kp_index = ikp_absorb[ik][j];

                            for (auto ig = 0; ig < pair.group.size(); ig++) {

                                k2 = dos->kmesh_dos->knum_sym(pair.group[ig].ks[0],
                                                              symmetry->SymmList[generating_sym].rot);
                                k3 = dos->kmesh_dos->knum_sym(pair.group[ig].ks[1],
                                                              symmetry->SymmList[generating_sym].rot);
                                k3_minus = dos->kmesh_dos->kindex_minus_xk[k3];

                                for (int ib = 0; ib < ns2; ++ib) {
                                    s2 = ib / ns;
                                    s3 = ib % ns;

                                    n1 = fb[k1][s1];
                                    n2 = fb[k2][s2];
                                    n3 = fb[k3][s3];
                                    for (ix = 0; ix < 3; ++ix) {
                                        Wks[s1][ix] +=
                                                (dFold[k2][s2][ix] - dFold[k3_minus][s3][ix]) * n1 * n2 * (n3 + 1.0)
                                                * L_absorb[kp_index][s1][ib];
                                    }
                                }
                            }
                        }

                    } // s1

                    average_W_at_k(k1, Wks);

                    for (s1 = 0; s1 < ns; ++s1) {

                        double Q_final = Q[ik][s1];
                        if (isotope->include_isotope) {
                            Q_final += fb[k1][s1] * (fb[k1][s1] + 1.0) * 2.0 * isotope_damping_loc[ik][s1];
                        }

                        if (conductivity->fph_rta > 0) {
                            Q_final += fb[k1][s1] * (fb[k1][s1] + 1.0) * 2.0 * damping4[itemp][ik][s1];
                        }

                        if (Q_final < 1.0e-50 || dos->dymat_dos->get_eigenvalues()[k1][s1] < eps8) {
                            for (ix = 0; ix < 3; ix++) dFnew[k1][s1][ix] = 0.0;
                        } else {
                            for (ix = 0; ix < 3; ix++) {
                                dFnew[k1][s1][ix] = (-vel[k1][s1][ix] * dndt[ik][s1] / beta - Wks[s1][ix]) / Q_final;
                            }
                        }
                        if (itr > 0) {
                            for (ix = 0; ix < 3; ix++) {
                                // a mixing factor of 0.75
                                // unstability in convergence happens sometime at low temperature.
                                dFnew[k1][s1][ix] =
                                        dFnew[k1][s1][ix] * mixing_factor + dFold[k1][s1][ix] * (1.0 - mixing_factor);
                            }
                        }

                    }

                } // ieq
            } // ik

            // reduce dFnew to dFold
            MPI_Allreduce(&dFnew[0][0][0],
                          &dFold[0][0][0],
                          nktot * ns * 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

            //average_dF(dFold);
            calc_kappa(itemp, dFold, kappa_new);

            //print kappa
            if (mympi->my_rank == 0) {
                for (ix = 0; ix < 3; ++ix) {
                    for (iy = 0; iy < 3; ++iy) {
                        std::cout << std::setw(12) << std::scientific
                                  << std::setprecision(2) << kappa_new[ix][iy];
                    }
                }
                std::cout << std::endl;
            }

            // check_convergence here
            auto converged = false;
            if (itr > 0) converged = check_convergence(kappa_old, kappa_new);

            if (converged) {
                // write to final kappa
                for (ix = 0; ix < 3; ++ix) {
                    for (iy = 0; iy < 3; ++iy) {
                        kappa[itemp][ix][iy] = kappa_new[ix][iy];
                    }
                }
                if (mympi->my_rank == 0) {
                    std::cout << "   -> iter   Kappa converged " << std::endl;
                }
                break;
            } else {
                // continue next loop
                for (ix = 0; ix < 3; ++ix) {
                    for (iy = 0; iy < 3; ++iy) {
                        kappa_old[ix][iy] = kappa_new[ix][iy];
                    }
                }
            }

            if (itr == (max_cycle - 1)) {
                // update kappa even if not converged
                for (ix = 0; ix < 3; ++ix) {
                    for (iy = 0; iy < 3; ++iy) {
                        kappa[itemp][ix][iy] = kappa_new[ix][iy];
                    }
                }
                if (mympi->my_rank == 0) {
                    std::cout << "   -> iter     Warning !! max cycle reached but kappa not converged " << std::endl;
                }
            }

        } // iter
        write_Q_dF(itemp, Q, dFold);

    } // itemp

    deallocate(Q);
    deallocate(dndt);

    deallocate(kappa_new);
    deallocate(kappa_old);
    deallocate(fb);
    deallocate(Wks);
    if (isotope->include_isotope) {
        deallocate(isotope_damping_loc);
        //deallocate(isotope_damping);
    }
    if (mympi->my_rank == 0) {
        fs_result.close();
    }
}

void Iterativebte::direct_solver()
{
    return;
}

void Iterativebte::calc_boson(int itemp, double **&b_out, double **&dndt_out)
{
    auto etemp = Temperature[itemp];
    double omega;
    for (auto ik = 0; ik < nktot; ++ik) {
        for (auto is = 0; is < ns; ++is) {
            omega = dos->dymat_dos->get_eigenvalues()[ik][is];
            b_out[ik][is] = thermodynamics->fB(omega, etemp);
        }
    }

    const double t_to_ryd = thermodynamics->T_to_Ryd;

    for (auto ik = 0; ik < nklocal; ++ik) {
        auto ikr = nk_l[ik];
        auto k1 = dos->kmesh_dos->kpoint_irred_all[ikr][0].knum;
        for (auto is = 0; is < ns; ++is) {
            omega = dos->dymat_dos->get_eigenvalues()[k1][is];
            auto x = omega / (t_to_ryd * etemp);
            dndt_out[ik][is] = std::pow(1.0 / (2.0 * sinh(0.5 * x)), 2) * x / etemp;
        }
    }
}

void Iterativebte::average_Q(double **&q1)
{
    double *q;
    double *tmp_omega;
    allocate(q, ns);
    allocate(tmp_omega, ns);

    const auto tol_omega = 1.0e-7; // Approximately equal to 0.01 cm^{-1}
    int s1, s2;

    for (auto ik = 0; ik < nklocal; ik++) {

        auto tmpk = nk_l[ik];

        const int k1 = dos->kmesh_dos->kpoint_irred_all[tmpk][0].knum;  // k index in full grid

        for (s1 = 0; s1 < ns; ++s1) tmp_omega[s1] = dos->dymat_dos->get_eigenvalues()[k1][s1];

        for (s1 = 0; s1 < ns; ++s1) {
            double sumq1 = 0.0;
            int countq1 = 0;
            for (s2 = 0; s2 < ns; ++s2) {
                if (std::abs(tmp_omega[s2] - tmp_omega[s1]) < tol_omega) {
                    sumq1 += q1[ik][s2];
                    countq1 += 1;
                }
            }
            q[s1] = sumq1 / static_cast<double>(countq1);
        }

        for (s1 = 0; s1 < ns; ++s1) q1[ik][s1] = q[s1];
    }

    deallocate(q);
    deallocate(tmp_omega);
}

void Iterativebte::average_dF(double ***&dF)
{
    double *tmp_omega;
    allocate(tmp_omega, ns);
    double **tmp_df;
    allocate(tmp_df, ns, 3);

    const auto tol_omega = 1.0e-7; // Approximately equal to 0.01 cm^{-1}
    int s1, s2;

    for (auto k1 = 0; k1 < nktot; k1++) {

        for (s1 = 0; s1 < ns; ++s1) {
            tmp_omega[s1] = dos->dymat_dos->get_eigenvalues()[k1][s1];
        }

        for (s1 = 0; s1 < ns; ++s1) {
            double sumdf[3];
            for (auto i = 0; i < 3; ++i) sumdf[i] = 0.0;
            int countq1 = 0;
            for (s2 = 0; s2 < ns; ++s2) {
                if (std::abs(tmp_omega[s2] - tmp_omega[s1]) < tol_omega) {
                    for (auto i = 0; i < 3; ++i) {
                        sumdf[i] += dF[k1][s2][i];
                    }
                    countq1 += 1;
                }
            }
            for (auto i = 0; i < 3; ++i) {
                tmp_df[s1][i] = sumdf[i] / static_cast<double>(countq1);
            }
        }

        for (s1 = 0; s1 < ns; ++s1) {
            for (auto i = 0; i < 3; ++i) {
                dF[k1][s1][i] = tmp_df[s1][i];
            }
        }
    }

    deallocate(tmp_df);
    deallocate(tmp_omega);
}

void Iterativebte::average_W_at_k(int k1, double **&W)
{
    double *tmp_omega;
    allocate(tmp_omega, ns);
    double **tmp_W;
    allocate(tmp_W, ns, 3);

    const auto tol_omega = 1.0e-7; // Approximately equal to 0.01 cm^{-1}
    int s1, s2;

    for (s1 = 0; s1 < ns; ++s1) {
        tmp_omega[s1] = dos->dymat_dos->get_eigenvalues()[k1][s1];
    }

    for (s1 = 0; s1 < ns; ++s1) {
        double sumdf[3];
        for (auto i = 0; i < 3; ++i) sumdf[i] = 0.0;
        int countq1 = 0;
        for (s2 = 0; s2 < ns; ++s2) {
            if (std::abs(tmp_omega[s2] - tmp_omega[s1]) < tol_omega) {
                for (auto i = 0; i < 3; ++i) {
                    sumdf[i] += W[s2][i];
                }
                countq1 += 1;
            }
        }
        for (auto i = 0; i < 3; ++i) {
            tmp_W[s1][i] = sumdf[i] / static_cast<double>(countq1);
        }
    }

    for (s1 = 0; s1 < ns; ++s1) {
        for (auto i = 0; i < 3; ++i) {
            W[s1][i] = tmp_W[s1][i];
        }
    }

    deallocate(tmp_W);
    deallocate(tmp_omega);
}

void Iterativebte::calc_kappa(int itemp, double ***&df, double **&kappa_out)
{
    auto etemp = Temperature[itemp];

    double omega;

    double beta = 1.0 / (thermodynamics->T_to_Ryd * etemp);
    double **tmpkappa;
    allocate(tmpkappa, 3, 3);

    for (auto ix = 0; ix < 3; ++ix) {
        for (auto iy = 0; iy < 3; ++iy) {
            tmpkappa[ix][iy] = 0.0;
        }
    }

    const double factor = Ryd / (time_ry * Bohr_in_Angstrom * 1.0e-10 * nktot * system->volume_p);

    for (auto k1 = 0; k1 < nktot; ++k1) {
        for (auto s1 = 0; s1 < ns; ++s1) {

            //if (mympi->my_rank == 1) {
            //    std::cout << df[k1][s1][0] << std::endl;
            //}

            omega = dos->dymat_dos->get_eigenvalues()[k1][s1];  // in Ry
            // df in unit bohr/K
            // vel, omega in atomic unit
            // kappa in W/mK
            double n1 = thermodynamics->fB(omega, etemp);

            for (auto ix = 0; ix < 3; ++ix) {
                for (auto iy = 0; iy < 3; ++iy) {
                    tmpkappa[ix][iy] += -omega * vel[k1][s1][ix] * beta * n1 * (n1 + 1.0) * df[k1][s1][iy];
                }
            }
        }
    }

    for (auto ix = 0; ix < 3; ++ix) {
        for (auto iy = 0; iy < 3; ++iy) {
            kappa_out[ix][iy] = tmpkappa[ix][iy] * factor;
            //std::cout << mympi->my_rank << " " << std::scientific << factor << std::endl;
        }
    }

    deallocate(tmpkappa);
}

bool Iterativebte::check_convergence(double **&k_old, double **&k_new)
{
    // it's better to compair three eigenvalues, which should be more robost
    double max_diff = -100;
    double diff;
    for (auto ix = 0; ix < 3; ++ix) {
        for (auto iy = 0; iy < 3; ++iy) {
            diff = std::abs(k_new[ix][iy] - k_old[ix][iy]) / std::abs(k_old[ix][iy]);
            if (diff < 0.0) exit("iterative solution", "negative kappa! ");
            if (diff > max_diff) max_diff = diff;
        }
    }
    if (max_diff < convergence_criteria) {
        return true;
    } else {
        return false;
    }
}

void Iterativebte::write_kappa()
{
    // taken from write_phonon
    if (mympi->my_rank == 0) {

        auto file_kappa = input->job_title + ".kl";

        std::ofstream ofs_kl;

        ofs_kl.open(file_kappa.c_str(), std::ios::out);
        if (!ofs_kl) exit("write_kappa", "Could not open file_kappa");

        ofs_kl << "# Iterative result." << std::endl;
        ofs_kl << "# Temperature [K], Thermal Conductivity (xx, xy, xz, yx, yy, yz, zx, zy, zz) [W/mK]" << std::endl;

        for (auto itemp = 0; itemp < ntemp; ++itemp) {
            ofs_kl << std::setw(10) << std::right << std::fixed << std::setprecision(2)
                   << Temperature[itemp];
            for (auto ix = 0; ix < 3; ++ix) {
                for (auto iy = 0; iy < 3; ++iy) {
                    ofs_kl << std::setw(15) << std::scientific
                           << std::setprecision(4) << kappa[itemp][ix][iy];
                }
            }
            ofs_kl << std::endl;
        }
        ofs_kl.close();
        std::cout << std::endl;
        std::cout << " -----------------------------------------------------------------" << std::endl << std::endl;
        std::cout << " Lattice thermal conductivity is stored in the file " << file_kappa << std::endl;
    }

}

void Iterativebte::write_result()
{
    // write Q and W for all phonon, only phonon in irreducible BZ is written
    int i;
    int nk_ir = dos->kmesh_dos->nk_irred;
    double Ry_to_kayser = Hz_to_kayser / time_ry;

    if (mympi->my_rank == 0) {
        std::cout << " Writing Q and W to file ..." << std::endl;

        fs_result.open(conductivity->get_filename_results(3).c_str(), std::ios::out);

        if (!fs_result) {
            exit("setup_result_io",
                 "Could not open file_result3");
        }

        fs_result << "## General information" << std::endl;
        fs_result << "#SYSTEM" << std::endl;
        fs_result << system->natmin << " " << system->nkd << std::endl;
        fs_result << system->volume_p << std::endl;
        fs_result << "#END SYSTEM" << std::endl;

        fs_result << "#KPOINT" << std::endl;
        fs_result << dos->kmesh_dos->nk_i[0] << " " << dos->kmesh_dos->nk_i[1] << " " << dos->kmesh_dos->nk_i[2]
                  << std::endl;
        fs_result << dos->kmesh_dos->nk_irred << std::endl;

        for (int i = 0; i < dos->kmesh_dos->nk_irred; ++i) {
            fs_result << std::setw(6) << i + 1 << ":";
            for (int j = 0; j < 3; ++j) {
                fs_result << std::setw(15)
                          << std::scientific << dos->kmesh_dos->kpoint_irred_all[i][0].kval[j];
            }
            fs_result << std::setw(12)
                      << std::fixed << dos->kmesh_dos->weight_k[i] << std::endl;
        }
        fs_result.unsetf(std::ios::fixed);

        fs_result << "#END KPOINT" << std::endl;

        fs_result << "#CLASSICAL" << std::endl;
        fs_result << thermodynamics->classical << std::endl;
        fs_result << "#END CLASSICAL" << std::endl;

        fs_result << "#FCSXML" << std::endl;
        fs_result << fcs_phonon->file_fcs << std::endl;
        fs_result << "#END  FCSXML" << std::endl;

        fs_result << "#SMEARING" << std::endl;
        fs_result << integration->ismear << std::endl;
        fs_result << integration->epsilon * Ry_to_kayser << std::endl;
        fs_result << "#END SMEARING" << std::endl;

        fs_result << "#TEMPERATURE" << std::endl;
        fs_result << system->Tmin << " " << system->Tmax << " " << system->dT << std::endl;
        fs_result << "#END TEMPERATURE" << std::endl;

        fs_result << "##END General information" << std::endl;

        fs_result << "##Phonon Frequency" << std::endl;
        fs_result << "#K-point (irreducible), Branch, Omega (cm^-1), Group velocity (m/s)" << std::endl;

        double factor = Bohr_in_Angstrom * 1.0e-10 / time_ry;
        for (i = 0; i < dos->kmesh_dos->nk_irred; ++i) {
            const int ik = dos->kmesh_dos->kpoint_irred_all[i][0].knum;
            for (auto is = 0; is < dynamical->neval; ++is) {
                fs_result << std::setw(6) << i + 1 << std::setw(6) << is + 1;
                fs_result << std::setw(15) << writes->in_kayser(dos->dymat_dos->get_eigenvalues()[ik][is]);
                fs_result << std::setw(15) << vel[ik][is][0] * factor
                          << std::setw(15) << vel[ik][is][1] * factor
                          << std::setw(15) << vel[ik][is][2] * factor << std::endl;
            }
        }

        fs_result << "##END Phonon Frequency" << std::endl << std::endl;
        fs_result << "##Q and W at each temperature" << std::endl;
    }
}

void Iterativebte::write_Q_dF(int itemp, double **&q, double ***&df)
{
    auto etemp = Temperature[itemp];

    auto nk_ir = dos->kmesh_dos->nk_irred;
    double **Q_tmp;
    double **Q_all;
    allocate(Q_all, nk_ir, ns);
    allocate(Q_tmp, nk_ir, ns);
    for (auto ik = 0; ik < nk_ir; ++ik) {
        for (auto is = 0; is < ns; ++is) {
            Q_all[ik][is] = 0.0;
            Q_tmp[ik][is] = 0.0;
        }
    }
    for (auto ik = 0; ik < nklocal; ++ik) {
        auto tmpk = nk_l[ik];
        for (auto is = 0; is < ns; ++is) {
            Q_tmp[tmpk][is] = q[ik][is];
        }
    }
    MPI_Allreduce(&Q_tmp[0][0], &Q_all[0][0],
                  nk_ir * ns, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    deallocate(Q_tmp);

    // now we have Q
    if (mympi->my_rank == 0) {
        fs_result << std::setw(10) << etemp << std::endl;

        for (auto ik = 0; ik < nk_ir; ++ik) {
            for (auto is = 0; is < ns; ++is) {
                auto k1 = dos->kmesh_dos->kpoint_irred_all[ik][0].knum;
                fs_result << std::setw(6) << ik + 1 << std::setw(6) << is + 1 << std::endl;
                fs_result
                        << std::setw(15) << std::scientific << std::setprecision(5) << Q_all[ik][is]
                        << std::setw(15) << std::scientific << std::setprecision(5) << df[k1][is][0]
                        << std::setw(15) << std::scientific << std::setprecision(5) << df[k1][is][1]
                        << std::setw(15) << std::scientific << std::setprecision(5) << df[k1][is][2] << std::endl;
            }
        }
        fs_result << std::endl;
    }
    deallocate(Q_all);
}

/*
void Iterativebte::setup_kpindex()
{
    // emitt
    int cnt;
    
    std::cout << "OK0" << std::endl;
    ikp_emitt.clear();
    cnt = 0;
    std::cout << "OK1" << std::endl;
    for (auto ik = 0; ik < nklocal; ++ik) {
        std::vector<int> counterk;
        counterk.clear();
        for (auto j = 0; j < localnk_triplets_emitt[ik].size(); ++j) {
            counterk.push_back(cnt);
            cnt += 1;
        }
        ikp_emitt.push_back(counterk);
    }
    std::cout << "OK2" << std::endl;
    // absorb
    //ikp_absorb.clear();
    cnt = 0;
    for (auto ik = 0; ik < nklocal; ++ik) {
        std::vector<int> counterk;
        counterk.clear();
        for (auto j = 0; j < localnk_triplets_absorb[ik].size(); ++j) {
            counterk.push_back(cnt);
            cnt += 1;
        }
        ikp_absorb.push_back(counterk);
    }
} */
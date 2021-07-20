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
}

void Iterativebte::deallocate_variables()
{
    if (Temperature) {memory->deallocate(Temperature);}
    if (kappa) {memory->deallocate(kappa);}
    if (vel) {memory->deallocate(vel);}
    if (dFold) {memory->deallocate(dFold);}
    if (dFnew) {memory->deallocate(dFnew);}  
    if (L_absorb) {memory->deallocate(L_absorb);}
    if (L_emitt) {memory->deallocate(L_emitt);}
}


void Iterativebte::setup_iterative()
{ 
    nktot = kpoint->nk;
    ns = dynamical->neval;
    ns2 = ns * ns;

    MPI_Bcast(&max_cycle, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&convergence_criteria, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
   
    sym_permutation = false;
    use_triplet_symmetry = true;

    // Temperature in K
    ntemp = static_cast<unsigned int>((system->Tmax - system->Tmin) / system->dT) + 1;
    memory->allocate(Temperature, ntemp);

    for (auto i = 0; i < ntemp; ++i) {
        Temperature[i] = system->Tmin + static_cast<double>(i) * system->dT;
    }

    // calculate vel

    memory->allocate(vel,nktot,ns,3);            // velocity of the full grid
    phonon_velocity->calc_phonon_vel_mesh(vel);  //this will gather to rank0 process

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

    memory->allocate(kappa,ntemp, 3, 3);

    // distribute q point among the processors
    auto nk_ir = kpoint->nk_irred;

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
    
    for (unsigned int i = 0; i < nklocal; ++i){

        auto ik = nk_l[i]; 
        std::vector<KsListGroup> triplet;
        std::vector<KsListGroup> triplet2;

        // k3 = k1 - k2
        kpoint->get_unique_triplet_k(ik, use_triplet_symmetry, 
                                        sym_permutation, triplet);

        // k3 = - (k1 + k2)
        kpoint->get_unique_triplet_k(ik, use_triplet_symmetry, 
                                        sym_permutation, triplet2, 1);
        
        counter += triplet.size();
        counter2+= triplet2.size();

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
        std::cout << "     size of L (MB) (approx.) = " << memory->memsize_in_MB(sizeof(double),kplength_absorb, ns, ns2) << std::endl;
        if (integration->ismear == 0 || integration->ismear == 1) {
            std::cout << "     smearing will be used for the delta function" << std::endl;
        } else if (integration->ismear == 2) {
            std::cout << "     adaptive methods will be used for the delta function" << std::endl;
        } else if (integration->ismear == -1) {
            std::cout << "     tetrahedron method will be used for the delta function" << std::endl;
        } else {
            error->exit("calc_L","ismear other than -1, 0, 1, 2 are not supported");
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
    
    memory->allocate(L_absorb, kplength_absorb, ns, ns2);
    memory->allocate(L_emitt,  kplength_emitt, ns, ns2);

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


    // emitt
    counter = 0;
    for (auto ik = 0; ik < nklocal; ++ik) {

        auto tmpk = nk_l[ik];
        k1 = kpoint->kpoint_irred_all[tmpk][0].knum;    // k index in full grid

        // emitt k1 -> k2 + k3 
        // V(-q1, q2, q3) delta(w1 - w2 - w3)
        for (auto j = 0; j < localnk_triplets_emitt[ik].size(); ++j) {
            
            auto pair = localnk_triplets_emitt[ik][j];

            k2 = pair.group[0].ks[0];
            k3 = pair.group[0].ks[1];
            //counter = ikp_emitt[ik][j];

            for (s1 = 0; s1 < ns; ++s1) {
                arr[0] = kpoint->knum_minus[k1] * ns + s1;
                omega1 = dynamical->eval_phonon[k1][s1];
                
                for (ib = 0; ib < ns2; ++ib) {
                    s2 = ib / ns;
                    s3 = ib % ns;

                    arr[1] = k2 * ns + s2;
                    arr[2] = k3 * ns + s3;
                    omega2 = dynamical->eval_phonon[k2][s2];
                    omega3 = dynamical->eval_phonon[k3][s3];
                    
                    if (integration->ismear == 0) {
                        delta = delta_lorentz(omega1 - omega2 - omega3, epsilon);
                    } else if (integration->ismear == 1) {
                        delta = delta_gauss(omega1 - omega2 - omega3, epsilon);
                    } else if (integration->ismear == 2) {
                        double epsilon2[2];
                        integration->adaptive_smearing(k2,s2,k3,s3,epsilon2);
                        //sum_smear += epsilon2[0] + epsilon2[1];
                        delta = delta_gauss(omega1 - omega2 - omega3, epsilon2[0]);
                    }
                    
                    v3_tmp =  std::norm(anharmonic_core->V3(arr,
                                        dynamical->eval_phonon,
                                        dynamical->evec_phonon));

                    L_emitt[counter][s1][ib] = (pi / 4.0) * v3_tmp * delta / static_cast<double>(nktot);
                    //sum1 += L_emitt[counter][s1][ib];
                }
            }
            counter += 1;
        }
    }
    if (counter != kplength_emitt) {
        error->exit("setup_L","Emitt: pair length not equal!");
    }

    counter = 0;
    for (auto ik = 0; ik < nklocal; ++ik) {

        auto tmpk = nk_l[ik];
        k1 = kpoint->kpoint_irred_all[tmpk][0].knum;    // k index in full grid

        // absorption k1 + k2 -> -k3
        // V(q1, q2, q3) since k3 = - (k1 + k2)
        for (auto j = 0; j < localnk_triplets_absorb[ik].size(); ++j) {
            
            auto pair = localnk_triplets_absorb[ik][j];
            //counter = ikp_absorb[ik][j];

            k2 = pair.group[0].ks[0];
            k3 = pair.group[0].ks[1];

            for (s1 = 0; s1 < ns; ++s1) {
                arr[0] = k1 * ns + s1;
                omega1 = dynamical->eval_phonon[k1][s1];
                
                for (ib = 0; ib < ns2; ++ib) {
                    s2 = ib / ns;
                    s3 = ib % ns;

                    arr[1] = k2 * ns + s2;
                    arr[2] = k3 * ns + s3;
                    omega2 = dynamical->eval_phonon[k2][s2];
                    omega3 = dynamical->eval_phonon[k3][s3];
                    
                    if (integration->ismear == 0) {
                        delta = delta_lorentz(omega1 + omega2 - omega3, epsilon);
                    } else if (integration->ismear == 1) {
                        delta = delta_gauss(omega1 + omega2 - omega3, epsilon);
                    } else if (integration->ismear == 2) {
                        double epsilon2[2];
                        integration->adaptive_smearing(k2,s2,k3,s3,epsilon2);
                        //sum_smear += epsilon2[0] + epsilon2[1];
                        delta = delta_gauss(omega1 + omega2 - omega3, epsilon2[1]);
                    }

                    v3_tmp =  std::norm(anharmonic_core->V3(arr,
                                        dynamical->eval_phonon,
                                        dynamical->evec_phonon));

                    L_absorb[counter][s1][ib] = (pi / 4.0) * v3_tmp * delta / static_cast<double>(nktot);
                    //sum2 += L_absorb[counter][s1][ib];
                }
            }
            counter += 1;
        }        
    }
    
    if (counter != kplength_absorb) {
        error->exit("setup_L","absorb: pair length not equal!");
    }

    //if (mympi->my_rank == 0) {
    //    std::cout << "  DONE !" << std::endl;
    //}
}

void Iterativebte::setup_L_tetra()
{
    // we calculate V for all pairs L+(local_nk*eachpair,ns,ns2) and L-
    
    memory->allocate(L_absorb, kplength_absorb, ns, ns2);
    memory->allocate(L_emitt,  kplength_emitt, ns, ns2);

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

    int *kmap_identity;
    memory->allocate(kmap_identity, nktot);
    for (auto i = 0; i < nktot; ++i) kmap_identity[i] = i;

    double *energy_tmp;
    double *weight_tetra;
    memory->allocate(energy_tmp, nktot);
    memory->allocate(weight_tetra, nktot);

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
    
    for (auto ik = 0; ik < nklocal; ++ik) {

        auto tmpk = nk_l[ik];
        k1 = kpoint->kpoint_irred_all[tmpk][0].knum;    // k index in full grid

        // emission: calc delta function w1 - w2 - w3, with k3 = k1 - k2
        for (s1 = 0; s1 < ns; ++s1) {

            omega1 = dynamical->eval_phonon[k1][s1];
            
            for (ib = 0; ib < ns2; ++ib) {
                s2 = ib / ns;
                s3 = ib % ns;

                for (k2 = 0; k2 < nktot; k2++) {
                    // k3 = k1 - k2
                    for (auto i = 0; i < 3; ++i) {
                        xk_tmp[i] = kpoint->xk[k1][i] - kpoint->xk[k2][i];
                    }
                    k3 = kpoint->get_knum(xk_tmp[0], xk_tmp[1], xk_tmp[2]);

                    omega2 = dynamical->eval_phonon[k2][s2];
                    omega3 = dynamical->eval_phonon[k3][s3];

                    energy_tmp[k2] = omega2 + omega3;
                } 
                integration->calc_weight_tetrahedron(nktot,
                                                     kmap_identity,
                                                     weight_tetra,
                                                     energy_tmp,
                                                     omega1);
                
                // emitt k1 -> k2 + k3 
                // V(-q1, q2, q3) delta(w1 - w2 - w3)
                for (auto j = 0; j < localnk_triplets_emitt[ik].size(); ++j) {
                    
                    auto pair = localnk_triplets_emitt[ik][j];
                    auto counter = ikp_emitt[ik][j];

                    k2 = pair.group[0].ks[0];
                    k3 = pair.group[0].ks[1];

                    arr[0] = kpoint->knum_minus[k1] * ns + s1;
                    arr[1] = k2 * ns + s2;
                    arr[2] = k3 * ns + s3;
                    delta = weight_tetra[k2];
                    v3_tmp =  std::norm(anharmonic_core->V3(arr,
                                        dynamical->eval_phonon,
                                        dynamical->evec_phonon));

                    L_emitt[counter][s1][ib] = (pi / 4.0) * v3_tmp * delta ;
                }

                // absorb
                for (k2 = 0; k2 < nktot; k2++) {

                    for (auto i = 0; i < 3; ++i) {
                        xk_tmp[i] = kpoint->xk[k1][i] + kpoint->xk[k2][i];
                    }
                    k3 = kpoint->get_knum(xk_tmp[0], xk_tmp[1], xk_tmp[2]);

                    omega2 = dynamical->eval_phonon[k2][s2];
                    omega3 = dynamical->eval_phonon[k3][s3];

                    energy_tmp[k2] = -omega2 + omega3;
                } 
                integration->calc_weight_tetrahedron(nktot,
                                                     kmap_identity,
                                                     weight_tetra,
                                                     energy_tmp,
                                                     omega1);
                
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
                    v3_tmp =  std::norm(anharmonic_core->V3(arr,
                                        dynamical->eval_phonon,
                                        dynamical->evec_phonon));

                    L_absorb[counter][s1][ib] = (pi / 4.0) * v3_tmp * delta ;
                }
            } // ib
        } // s1    
    } // ik
    
    //if (mympi->my_rank == 0) {
    //    std::cout << "  DONE !" << std::endl;
    //}
    memory->deallocate(kmap_identity);
    memory->deallocate(energy_tmp);
    memory->deallocate(weight_tetra);
}


void Iterativebte::calc_Q_from_L(double** &n, double** &q1)
{ 
    int s1, s2, s3;
    double ph1;
    double n1, n2, n3;


    double **Qemit;
    double **Qabsorb;
    memory->allocate(Qemit,nklocal,ns);
    memory->allocate(Qabsorb,nklocal,ns);

    for (auto ik = 0; ik < nklocal; ++ik) {
        for (s1 = 0; s1 < ns; ++s1) { 
            Qemit[ik][s1] = 0.0;
            Qabsorb[ik][s1] = 0.0;
        }
    }

    unsigned int counter ;

    // emit
    counter = 0;
    for (auto ik = 0; ik < nklocal; ++ik) {

        auto tmpk = nk_l[ik];
        const int k1 = kpoint->kpoint_irred_all[tmpk][0].knum; 

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
                    Qemit[ik][s1] += 0.5 * ( n1 * (n2 + 1.0) * (n3 + 1.0) ) * L_emitt[counter][s1][ib] * multi;
                }
            }
            counter += 1;
        }

    }
    if (counter != kplength_emitt) {
        error->exit("setup_L","Emitt: pair length not equal!");
    }

    // absorb k1 + k2 -> -k3
    counter = 0;
    for (auto ik = 0; ik < nklocal; ++ik) {

        auto tmpk = nk_l[ik];
        const int k1 = kpoint->kpoint_irred_all[tmpk][0].knum; 

        for (auto j = 0; j < localnk_triplets_absorb[ik].size(); ++j) {
                
            auto pair = localnk_triplets_absorb[ik][j];
            auto multi = static_cast<double>(pair.group.size());
            auto k2 = pair.group[0].ks[0];
            auto k3 = pair.group[0].ks[1];
            
            for (s1 = 0; s1 < ns; ++s1){
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
            //error->exit("1","2");
    }
    
    if (counter != kplength_absorb) {
        error->exit("setup_L","absorb: pair length not equal!");
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
    memory->deallocate(Qemit);
    memory->deallocate(Qabsorb);
}


void Iterativebte::iterative_solver()
{
    // calculate the Vs
    double mixing_factor = 0.75;
    
    double **Q;
    double **kappa_new;
    double **kappa_old;

    //bool converged;

    memory->allocate(kappa_new,3,3);
    memory->allocate(kappa_old,3,3);
    memory->allocate(Q,nklocal,ns);

    memory->allocate(dFold,nktot,ns,3);
    memory->allocate(dFnew,nktot,ns,3);

    double **fb;
    double **dndt;
    memory->allocate(dndt,nklocal,ns);
    memory->allocate(fb, nktot, ns);


    double **isotope_damping_loc;
    if (isotope->include_isotope) {
        double **isotope_damping;
        memory->allocate(isotope_damping, kpoint->nk_irred, ns);

        if (mympi->my_rank == 0) {
            for (auto ik = 0; ik < kpoint->nk_irred; ik++) {
                for (auto is = 0; is < ns; is++) {
                    isotope_damping[ik][is] = isotope->gamma_isotope[ik][is];
                }
            }
        }

        MPI_Bcast(&isotope_damping[0][0], kpoint->nk_irred * ns, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        memory->allocate(isotope_damping_loc, nklocal, ns); // this is for reducing some memory usage
        for (auto ik = 0; ik < nklocal; ik++) {
            auto tmpk = nk_l[ik];
            for (auto is = 0; is < ns; is++) {
                isotope_damping_loc[ik][is] = isotope_damping[tmpk][is];
            }
        }

        memory->deallocate(isotope_damping); 
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
    double** Wks;
    memory->allocate(Wks,ns,3);

    for (auto itemp = 0; itemp < ntemp; ++itemp){
        
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

        calc_Q_from_L(fb,Q);
        average_Q(Q);

        for (ik = 0; ik < nktot; ++ik) {
            for (is = 0; is < ns; ++is){
                for (ix = 0; ix < 3; ++ix){
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
                for (is = 0; is < ns; ++is){
                    for (ix = 0; ix < 3; ++ix){
                        dFnew[ik][is][ix] = 0.0;
                    }
                }
            }

            for (ik = 0; ik < nklocal; ++ik){

                auto tmpk = nk_l[ik];
                auto num_equivalent = kpoint->kpoint_irred_all[tmpk].size();
                auto kref = kpoint->kpoint_irred_all[tmpk][0].knum;    

                for (auto ieq = 0; ieq < num_equivalent; ++ieq){

                    k1 = kpoint->kpoint_irred_all[tmpk][ieq].knum;      // k1 will go through all points

                    generating_sym = -1;
                    for (isym = 0; isym < nsym; ++isym){
                        auto krot = kpoint->knum_sym(kref,isym);
                        if (k1 == krot) generating_sym = isym;      
                    }
                    if (generating_sym == -1) {
                        error->exit("iterative solution","cannot find all equivalent k");
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

                            for (auto ig = 0; ig < pair.group.size(); ig++){

                                k2 = kpoint->knum_sym(pair.group[ig].ks[0],generating_sym);
                                k3 = kpoint->knum_sym(pair.group[ig].ks[1],generating_sym);
                                
                                for (int ib = 0; ib < ns2; ++ib) {
                                    s2 = ib / ns;
                                    s3 = ib % ns;
                                    
                                    n1 = fb[k1][s1];
                                    n2 = fb[k2][s2];
                                    n3 = fb[k3][s3];
                                    for (ix = 0; ix < 3; ++ix) {
                                        Wks[s1][ix] -= 0.5 * (dFold[k2][s2][ix] + dFold[k3][s3][ix]) * n1 * (n2+1.0) * (n3+1.0) * L_emitt[kp_index][s1][ib];
                                    }
                                }
                            }
                        }
                        
                        // absorb k1 + k2 -> -k3
                        for (auto j = 0; j < localnk_triplets_absorb[ik].size(); ++j) {
            
                            auto pair = localnk_triplets_absorb[ik][j];
                            int kp_index = ikp_absorb[ik][j];

                            for (auto ig = 0; ig < pair.group.size(); ig++){

                                k2 = kpoint->knum_sym(pair.group[ig].ks[0],generating_sym);
                                k3 = kpoint->knum_sym(pair.group[ig].ks[1],generating_sym);
                                k3_minus = kpoint->knum_minus[k3];
                                
                                for (int ib = 0; ib < ns2; ++ib) {
                                    s2 = ib / ns;
                                    s3 = ib % ns;
                                    
                                    n1 = fb[k1][s1];
                                    n2 = fb[k2][s2];
                                    n3 = fb[k3][s3];
                                    for (ix = 0; ix < 3; ++ix) {
                                        Wks[s1][ix] += (dFold[k2][s2][ix] - dFold[k3_minus][s3][ix]) * n1 * n2 * (n3+1.0) * L_absorb[kp_index][s1][ib];
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

                        if (Q_final < 1.0e-50 || dynamical->eval_phonon[k1][s1] < eps8) {
                            for (ix = 0; ix < 3; ix++) dFnew[k1][s1][ix] = 0.0;
                        } else {
                            for (ix = 0; ix < 3; ix++) {
                                dFnew[k1][s1][ix] = ( - vel[k1][s1][ix] * dndt[ik][s1] / beta - Wks[s1][ix] ) / Q_final;
                            }
                        }
                        if (itr > 0) {
                            for (ix = 0; ix < 3; ix++) {
                                // a mixing factor of 0.75
                                // unstability in convergence happens sometime at low temperature.
                                dFnew[k1][s1][ix] = dFnew[k1][s1][ix] * mixing_factor + dFold[k1][s1][ix] * (1.0 - mixing_factor); 
                            }
                        }

                    }

                } // ieq
            } // ik

            // reduce dFnew to dFold
            MPI_Allreduce(  &dFnew[0][0][0], 
                            &dFold[0][0][0], 
                            nktot * ns * 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            
            //average_dF(dFold);
            calc_kappa(itemp, dFold, kappa_new);

            //print kappa
            if (mympi->my_rank == 0) {
                for (ix = 0; ix < 3; ++ix) {
                    for (iy = 0; iy < 3; ++ iy) {
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
                    for (iy = 0; iy < 3; ++ iy) {
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
                    for (iy = 0; iy < 3; ++ iy) {
                        kappa_old[ix][iy] = kappa_new[ix][iy];
                    }
                }
            }

            if (itr == (max_cycle-1)) {
                // update kappa even if not converged
                for (ix = 0; ix < 3; ++ix) {
                    for (iy = 0; iy < 3; ++ iy) {
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

    memory->deallocate(Q);
    memory->deallocate(dndt);

    memory->deallocate(kappa_new);
    memory->deallocate(kappa_old);
    memory->deallocate(fb);
    memory->deallocate(Wks);
    if (isotope->include_isotope) {
        memory->deallocate(isotope_damping_loc);
        //memory->deallocate(isotope_damping);
    }
    if (mympi->my_rank == 0) {
        writes->fs_result.close();
    }
}


void Iterativebte::direct_solver()
{
    return;
}

void Iterativebte::calc_boson(int itemp, double** &b_out, double** &dndt_out) 
{
    auto etemp = Temperature[itemp];
    double omega;
    for (auto ik = 0; ik < nktot; ++ik) {
        for (auto is = 0; is < ns; ++is) {
            omega = dynamical->eval_phonon[ik][is];
            b_out[ik][is] = thermodynamics->fB(omega,etemp);
        }
    }

    const double t_to_ryd = thermodynamics->T_to_Ryd;

    for (auto ik = 0; ik < nklocal; ++ik) {
        auto ikr = nk_l[ik];
        auto k1 = kpoint->kpoint_irred_all[ikr][0].knum; 
        for (auto is = 0; is < ns; ++is) {
            omega = dynamical->eval_phonon[k1][is];
            auto x = omega / (t_to_ryd * etemp);
            dndt_out[ik][is] = std::pow(1.0/(2.0 * sinh(0.5 * x)),2) * x / etemp;
        }
    }
}


void Iterativebte::average_Q(double** &q1)
{
    double *q;
    double *tmp_omega;
    memory->allocate(q, ns);
    memory->allocate(tmp_omega, ns);

    const auto tol_omega = 1.0e-7; // Approximately equal to 0.01 cm^{-1}
    int s1, s2;

    for (auto ik = 0; ik < nklocal; ik++) {

        auto tmpk = nk_l[ik];
        
        const int k1 = kpoint->kpoint_irred_all[tmpk][0].knum;  // k index in full grid

        for (s1 = 0; s1 < ns; ++s1) tmp_omega[s1] = dynamical->eval_phonon[k1][s1];

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

    memory->deallocate(q);
    memory->deallocate(tmp_omega);
}


void Iterativebte::average_dF(double*** &dF)
{
    double *tmp_omega;
    memory->allocate(tmp_omega, ns);
    double **tmp_df;
    memory->allocate(tmp_df, ns, 3);

    const auto tol_omega = 1.0e-7; // Approximately equal to 0.01 cm^{-1}
    int s1, s2;

    for (auto k1 = 0; k1 < nktot; k1++) {

        for (s1 = 0; s1 < ns; ++s1) { 
            tmp_omega[s1] = dynamical->eval_phonon[k1][s1];
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

    memory->deallocate(tmp_df);
    memory->deallocate(tmp_omega);
}


void Iterativebte::average_W_at_k(int k1, double** &W)
{
    double *tmp_omega;
    memory->allocate(tmp_omega, ns);
    double **tmp_W;
    memory->allocate(tmp_W, ns, 3);

    const auto tol_omega = 1.0e-7; // Approximately equal to 0.01 cm^{-1}
    int s1, s2;


    for (s1 = 0; s1 < ns; ++s1) { 
        tmp_omega[s1] = dynamical->eval_phonon[k1][s1];
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

    memory->deallocate(tmp_W);
    memory->deallocate(tmp_omega);
}


void Iterativebte::calc_kappa(int itemp, double*** &df, double** &kappa_out)
{
    auto etemp = Temperature[itemp];

    double omega;

    double beta = 1.0 / (thermodynamics->T_to_Ryd * etemp);
    double **tmpkappa;
    memory->allocate(tmpkappa,3,3);

    for (auto ix = 0; ix < 3; ++ix) {
        for (auto iy = 0; iy < 3; ++iy) {
            tmpkappa[ix][iy] = 0.0;
        }
    }
    
    const double factor = Ryd / ( time_ry * Bohr_in_Angstrom * 1.0e-10 * nktot * system->volume_p);

    for (auto k1 = 0; k1 < nktot; ++k1) {
        for (auto s1 = 0; s1 < ns; ++s1){

            //if (mympi->my_rank == 1) {
            //    std::cout << df[k1][s1][0] << std::endl;
            //}

            omega = dynamical->eval_phonon[k1][s1];  // in Ry
            // df in unit bohr/K
            // vel, omega in atomic unit
            // kappa in W/mK
            double n1 = thermodynamics->fB(omega,etemp);

            for (auto ix = 0; ix < 3; ++ix) {
                for (auto iy = 0; iy < 3; ++iy) {
                    tmpkappa[ix][iy] += - omega * vel[k1][s1][ix] * beta * n1 * (n1+1.0) * df[k1][s1][iy];
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

    memory->deallocate(tmpkappa);
}


bool Iterativebte::check_convergence(double** &k_old, double** &k_new)
{
    // it's better to compair three eigenvalues, which should be more robost
    double max_diff = -100;
    double diff;
    for (auto ix = 0; ix < 3; ++ix) {
        for (auto iy = 0; iy < 3; ++iy) {
            diff = std::abs( k_new[ix][iy] - k_old[ix][iy] ) / std::abs(k_old[ix][iy]);
            if (diff < 0.0) error->exit("iterative solution","negative kappa! ");
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
        if (!ofs_kl) error->exit("write_kappa", "Could not open file_kappa");
        
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
    int nk_ir = kpoint->nk_irred;
    double Ry_to_kayser = Hz_to_kayser / time_ry;
    
    if (mympi->my_rank == 0) {
        std::cout << " Writing Q and W to file ..." << std::endl;
        
        writes->fs_result.open(writes->file_result.c_str(), std::ios::out);

        if (!writes->fs_result) {
            error->exit("setup_result_io",
                            "Could not open file_result");
        }

        writes->fs_result << "## General information" << std::endl;
        writes->fs_result << "#SYSTEM" << std::endl;
        writes->fs_result << system->natmin << " " << system->nkd << std::endl;
        writes->fs_result << system->volume_p << std::endl;
        writes->fs_result << "#END SYSTEM" << std::endl;

        writes->fs_result << "#KPOINT" << std::endl;
        writes->fs_result << kpoint->nkx << " " << kpoint->nky << " " << kpoint->nkz << std::endl;
        writes->fs_result << kpoint->nk_irred << std::endl;

        for (int i = 0; i < kpoint->nk_irred; ++i) {
            writes->fs_result << std::setw(6) << i + 1 << ":";
            for (int j = 0; j < 3; ++j) {
                writes->fs_result << std::setw(15)
                              << std::scientific << kpoint->kpoint_irred_all[i][0].kval[j];
            }
            writes->fs_result << std::setw(12)
                          << std::fixed << kpoint->weight_k[i] << std::endl;
        }
        writes->fs_result.unsetf(std::ios::fixed);

        writes->fs_result << "#END KPOINT" << std::endl;

        writes->fs_result << "#CLASSICAL" << std::endl;
        writes->fs_result << thermodynamics->classical << std::endl;
        writes->fs_result << "#END CLASSICAL" << std::endl;

        writes->fs_result << "#FCSXML" << std::endl;
        writes->fs_result << fcs_phonon->file_fcs << std::endl;
        writes->fs_result << "#END  FCSXML" << std::endl;

        writes->fs_result << "#SMEARING" << std::endl;
        writes->fs_result << integration->ismear << std::endl;
        writes->fs_result << integration->epsilon * Ry_to_kayser << std::endl;
        writes->fs_result << "#END SMEARING" << std::endl;

        writes->fs_result << "#TEMPERATURE" << std::endl;
        writes->fs_result << system->Tmin << " " << system->Tmax << " " << system->dT << std::endl;
        writes->fs_result << "#END TEMPERATURE" << std::endl;

        writes->fs_result << "##END General information" << std::endl;

        writes->fs_result << "##Phonon Frequency" << std::endl;
        writes->fs_result << "#K-point (irreducible), Branch, Omega (cm^-1), Group velocity (m/s)" << std::endl;

        double factor = Bohr_in_Angstrom * 1.0e-10 / time_ry;
        for (i = 0; i < kpoint->nk_irred; ++i) {
            const int ik = kpoint->kpoint_irred_all[i][0].knum;
            for (auto is = 0; is < dynamical->neval; ++is) {
                writes->fs_result << std::setw(6) << i + 1 << std::setw(6) << is + 1;
                writes->fs_result << std::setw(15) << writes->in_kayser(dynamical->eval_phonon[ik][is]); 
                writes->fs_result << std::setw(15) << vel[ik][is][0] * factor 
                                  << std::setw(15) << vel[ik][is][1] * factor 
                                  << std::setw(15) << vel[ik][is][2] * factor << std::endl;
            }
        }

        writes->fs_result << "##END Phonon Frequency" << std::endl << std::endl;
        writes->fs_result << "##Q and W at each temperature" << std::endl;
    }
}

void Iterativebte::write_Q_dF(int itemp, double** &q, double*** &df)
{
    auto etemp = Temperature[itemp];

    auto nk_ir = kpoint->nk_irred;
    double** Q_tmp;
    double** Q_all;
    memory->allocate(Q_all, nk_ir, ns);
    memory->allocate(Q_tmp, nk_ir, ns);
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
                  nk_ir * ns , MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    memory->deallocate(Q_tmp);

    // now we have Q
    if (mympi->my_rank == 0) {
        writes->fs_result << std::setw(10) << etemp << std::endl;
        
        for (auto ik = 0; ik < nk_ir; ++ik) {
            for (auto is = 0; is < ns; ++is) {
                auto k1 = kpoint->kpoint_irred_all[ik][0].knum;
                writes->fs_result << std::setw(6) << ik + 1 << std::setw(6) << is + 1 << std::endl;
                writes->fs_result 
                    << std::setw(15) << std::scientific << std::setprecision(5) << Q_all[ik][is]
                    << std::setw(15) << std::scientific << std::setprecision(5) << df[k1][is][0]
                    << std::setw(15) << std::scientific << std::setprecision(5) << df[k1][is][1]
                    << std::setw(15) << std::scientific << std::setprecision(5) << df[k1][is][2] << std::endl;
            }
        }
        writes->fs_result << std::endl;
    }
    memory->deallocate(Q_all);
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
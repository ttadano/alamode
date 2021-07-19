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

// perhaps I need to separate delta(w1 + w2 - w3) and delta(w1 + w3 - w2) 
// for the tetrahedron method to work

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
    ntemp = 0;
    max_cycle = 20;
    convergence_criteria = 0.005; 
    use_triplet_symmetry = true;
    sym_permutation = true;
    stable_version = true;

    Temperature = nullptr;
    kappa = nullptr;
    v3 = nullptr;
    vel = nullptr;
    dFold = nullptr;
    dFnew = nullptr;
    delta1 = nullptr;
    delta2 = nullptr;
    delta3 = nullptr;
}


void Iterativebte::deallocate_variables()
{
    if (Temperature) {memory->deallocate(Temperature);}
    if (kappa) {memory->deallocate(kappa);}
    if (v3) {memory->deallocate(v3);}
    if (delta1) {memory->deallocate(delta1);}
    if (delta2) {memory->deallocate(delta2);}
    if (delta3) {memory->deallocate(delta3);}
    if (vel) {memory->deallocate(vel);}
    if (dFold) {memory->deallocate(dFold);}
    if (dFnew) {memory->deallocate(dFnew);}  
}


void Iterativebte::setup_iterative()
{ 
    nktot = kpoint->nk;
    ns = dynamical->neval;
    ns2 = ns * ns;

    MPI_Bcast(&max_cycle, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&convergence_criteria, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
   
    sym_permutation = true;
    use_triplet_symmetry = true;

    // Temperature in K
    ntemp = static_cast<unsigned int>((system->Tmax - system->Tmin) / system->dT) + 1;
    memory->allocate(Temperature, ntemp);

    for (auto i = 0; i < ntemp; ++i) {
        Temperature[i] = system->Tmin + static_cast<double>(i) * system->dT;
    }

    // calculate vel
    memory->allocate(vel,nktot,ns,3);
    phonon_velocity->calc_phonon_vel_mesh(vel);  //this will gather to rank0 process

    const auto factor = Bohr_in_Angstrom * 1.0e-10 / time_ry;
    MPI_Bcast(&vel[0][0][0], nktot * ns * 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);

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

    write_result();
}


void Iterativebte::get_triplets()
{
    num_unipair.clear();
    start_unipair.clear();
    localnk_triplets.clear();

    int counter = 0;
    
    for (unsigned int i = 0; i < nklocal; ++i){

        auto ik = nk_l[i]; 
        std::vector<KsListGroup> triplet;

        kpoint->get_unique_triplet_k(ik, use_triplet_symmetry, 
                                        sym_permutation, triplet);
        
        auto size = triplet.size();

        start_unipair.push_back(counter);
        counter += size;

        num_unipair.push_back(size);

        for (unsigned int j = 0; j < size; ++j){
            localnk_triplets.push_back(triplet[j]);
        }
        
    }

    kplength = localnk_triplets.size();  
}


void Iterativebte::do_iterativebte()
{
    setup_v3();
    iterative_solver();
    //write_kappa();
}


void Iterativebte::setup_v3()
{
    // we calculate V for all pairs V(local_nk*eachpair,ns,ns2)
    // -> 2pi/hbar |V3|^2
    
    memory->allocate(v3, kplength, ns, ns2);
    if (mympi->my_rank == 0) {
        std::cout << " we calculate once for the matrix elements v3 ..." << std::endl;
        std::cout << " size of v3 (MB) = " << memory->memsize_in_MB(sizeof(double),kplength, ns, ns2) << std::endl;
    }

    unsigned int arr[3];
    int k1, k2, k3, k1_minus;
    int s1, s2, s3;
    int ib;
    
    for (auto ik = 0; ik < nklocal; ++ik) {

        auto tmpk = nk_l[ik];
        
        k1 = kpoint->kpoint_irred_all[tmpk][0].knum;    // k index in full grid
        k1_minus = kpoint->knum_minus[k1];            // -k index in full grid 

        for (auto j = 0; j < num_unipair[ik]; ++j) {
            
            auto kp = start_unipair[ik]+j;
            auto pair = localnk_triplets[kp];

            k2 = pair.group[0].ks[0];
            k3 = pair.group[0].ks[1];

            for (s1 = 0; s1 < ns; ++s1) {
                arr[0] = k1_minus * ns + s1;
                
                for (ib = 0; ib < ns2; ++ib) {
                    s2 = ib / ns;
                    s3 = ib % ns;

                    arr[1] = k2 * ns + s2;
                    arr[2] = k3 * ns + s3;
                    
                    v3[kp][s1][ib] = (pi / 4.0) * std::norm(anharmonic_core->V3(arr,
                                                            dynamical->eval_phonon,
                                                            dynamical->evec_phonon)) / static_cast<double>(nktot);
                    
                }
            }

        }
    }
    if (mympi->my_rank == 0) {
        std::cout << "  DONE !" << std::endl;
    }
}


void Iterativebte::iterative_solver()
{
    // calculate the Vs
    double mixing_factor = 0.75;
    
    double **Q;
    double **kappa_new;
    double **kappa_old;

    //bool converged;

    //double ***delta1;
    //double ***delta2;
    //double ***delta3;
    
    memory->allocate(delta1, kplength, ns, ns2);
    memory->allocate(delta2, kplength, ns, ns2);
    memory->allocate(delta3, kplength, ns, ns2);

    if (integration->ismear == 0 || integration->ismear == 1) {
        calc_delta_smear(delta1, delta2, delta3);
    } else if (integration->ismear == -1) {
        calc_delta_tetra(delta1, delta2, delta3);
    }

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

    // we solve iteratively for each temperature
    int ik, is, ix, iy;
    double n1, n2, n3;
     
    if (mympi->my_rank == 0) {
        std::cout << std::endl << " we start iteration ..." << std::endl << std::endl;
    }
    
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

    
        //calc_phi_smear(itemp, phi1, phi2, phi3);

        calc_Q_from_phi1(fb,Q);

        average_Q(Q);

        for (ik = 0; ik < nktot; ++ik) {
            for (is = 0; is < ns; ++is){
                for (ix = 0; ix < 3; ++ix){
                    dFold[ik][is][ix] = 0.0;
                }
            }
        }

        int s1, s2, s3;
        int k1, k2, k3, k2_minus, k3_minus;
        int  generating_sym;
        int nsym = symmetry->SymmList.size();
        int isym;
        double d1, d2, d3;

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

            // calculate W
            // W = \sum_{q2,q3} 2pi/hbar |v3|^2 [
            // (f(-q2) - f(q3)) * n1 * n2 * (n3 + 1) * delta(w1 + w2 - w3) 
            //  -0.5 * (f(q2) + f(q3)) * n1 * (n2 + 1) * (n3 + 1) * delta(w1 - w2 - w3) ]
            for (ik = 0; ik < nklocal; ++ik){

                auto tmpk = nk_l[ik];
                auto num_equivalent = kpoint->kpoint_irred_all[tmpk].size();
                auto kref = kpoint->kpoint_irred_all[tmpk][0].knum;     // the first of the equivalent q point

                for (auto ieq = 0; ieq < num_equivalent; ++ieq){

                    k1 = kpoint->kpoint_irred_all[tmpk][ieq].knum;      // one of the equivalent q point

                    generating_sym = -1;
                    for (isym = 0; isym < nsym; ++isym){
                        auto krot = kpoint->knum_sym(kref,isym);
                        if (k1 == krot) generating_sym = isym;          // which symmetry generate k1 from kref
                    }
                    if (generating_sym == -1) {
                        error->exit("iterative solution","cannot find all equivalent k");
                    }

                    for (s1 = 0; s1 < ns; ++s1) {

                        for (ix = 0; ix < 3; ++ix) {
                            Wks[s1][ix] = 0.0;
                        }

                        for (auto j = 0; j < num_unipair[ik]; ++j) {

                            auto kp = start_unipair[ik]+j;
                            auto pair = localnk_triplets[kp];

                            //we go through each pairs
                            for (auto ig = 0; ig < pair.group.size(); ig++){

                                k2 = kpoint->knum_sym(pair.group[ig].ks[0],generating_sym);
                                k3 = kpoint->knum_sym(pair.group[ig].ks[1],generating_sym);
                                k2_minus = kpoint->knum_minus[k2];
                                k3_minus = kpoint->knum_minus[k3];
                                
                                for (int ib = 0; ib < ns2; ++ib) {
                                    s2 = ib / ns;
                                    s3 = ib % ns;
                                    
                                    n1 = fb[k1][s1];
                                    n2 = fb[k2][s2];
                                    n3 = fb[k3][s3];
                                    d1 = delta1[kp][s1][ib];
                                    d2 = delta2[kp][s1][ib];
                                    d3 = delta3[kp][s1][ib];

                                    for (ix = 0; ix < 3; ++ix) {
                                        double tmp_prod = 0.5 * (dFold[k2_minus][s2][ix] - dFold[k3][s3][ix]) * n1 * n2 * (n3 + 1.0) * d1
                                                        + 0.5 * (dFold[k3_minus][s3][ix] - dFold[k2][s2][ix]) * n1 * n3 * (n2 + 1.0) * d3
                                                        - 0.5 * (dFold[k2][s2][ix] + dFold[k3][s3][ix]) * n1 * (n2 + 1.0) * (n3 + 1.0) * d2 ;
                                        Wks[s1][ix] += v3[kp][s1][ib] * tmp_prod;
                                        //                * (dFold[k2_minus][s2][ix] - dFold[k3][s3][ix]) * n1 * n2 * (n3 + 1.0) * d1 ;
                                        //Wks[s1][ix] -= v3[kp][s1][ib] * 0.5 * (dFold[k2][s2][ix] + dFold[k3][s3][ix]) * n1 * (n2 + 1.0) * (n3 + 1.0) * d2 ;
                                        //f (mympi->my_rank == 1) std::cout << phi2[kp][s1][ib] << std::endl;
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
                            for (ix = 0; ix < 3; ++ix){
                                dFnew[k1][s1][ix] = dFnew[k1][s1][ix] * mixing_factor + dFold[k1][s1][ix] * (1.0 - mixing_factor); 
                            } 
                        }
                    }

                }
            } // ik

            MPI_Allreduce(  &dFnew[0][0][0], 
                            &dFold[0][0][0], 
                            nktot * ns * 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            
            calc_kappa(itemp, dFold, kappa_new);
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
            //
            
        } // iter
        write_Q_dF(itemp, Q, dFold);

    } // itemp
    if (mympi->my_rank == 0) {
        writes->fs_result.close();
    }

    memory->deallocate(Q);
    memory->deallocate(dndt);

    memory->deallocate(kappa_new);
    memory->deallocate(kappa_old);
    memory->deallocate(fb);
    if (isotope->include_isotope) {
        memory->deallocate(isotope_damping_loc);
        //memory->deallocate(isotope_damping);
    }
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

void Iterativebte::calc_delta_smear(double*** &d1, double*** &d2, double*** &d3)
{
    //calculate once for the three delta functions that contain delta
    // d1 = delta(omega1 + omega2 - omega3)
    // d2 = delta(omega1 - omega2 - omega3)
    // d3 = delta(omega1 - omega2 + omega3)
    
    int s1, s2, s3;
    double omega1, omega2, omega3;
    auto epsilon = integration->epsilon;
    
    for (auto ik = 0; ik < nklocal; ++ik) {

        auto tmpk = nk_l[ik];
        
        const int k1 = kpoint->kpoint_irred_all[tmpk][0].knum;  // k index in full grid

        for (auto j = 0; j < num_unipair[ik]; ++j) {
            
            auto kp = start_unipair[ik]+j;
            auto pair = localnk_triplets[kp];

            auto k2 = pair.group[0].ks[0];
            auto k3 = pair.group[0].ks[1];

            for (s1 = 0; s1 < ns; ++s1) {

                omega1 = dynamical->eval_phonon[k1][s1];

                for (int ib = 0; ib < ns2; ++ib) {
                    s2 = ib / ns;
                    s3 = ib % ns;

                    omega2 = dynamical->eval_phonon[k2][s2];
                    omega3 = dynamical->eval_phonon[k3][s3];
                    
                    if (integration->ismear == 0) {
                        d1[kp][s1][ib] = delta_lorentz(omega1 + omega2 - omega3, epsilon);
                        d2[kp][s1][ib] = delta_lorentz(omega1 - omega2 - omega3, epsilon);
                        d3[kp][s1][ib] = delta_lorentz(omega1 - omega2 + omega3, epsilon);
                    } else if (integration->ismear == 1) {
                        d1[kp][s1][ib] = delta_gauss(omega1 + omega2 - omega3, epsilon);
                        d2[kp][s1][ib] = delta_gauss(omega1 - omega2 - omega3, epsilon);
                        d3[kp][s1][ib] = delta_gauss(omega1 - omega2 + omega3, epsilon);
                    } else {
                        error->exit("phi","ismear");
                    }
                }
            }
        }
    }
}

void Iterativebte::calc_delta_tetra(double*** &d1, double*** &d2, double*** &d3)
{
    //calculate once for the three delta functions that contain delta
    // d1 = delta(omega1 + omega2 - omega3)
    // d2 = delta(omega1 - omega2 - omega3)
    
    int s1, s2, s3;
    int k1, k2, k3;
    double omega1, omega2, omega3;
    // auto epsilon = integration->epsilon;
    int i;

    int *kmap_identity;
    memory->allocate(kmap_identity, nktot);
    for (i = 0; i < nktot; ++i) kmap_identity[i] = i;

    double **energy_tmp;
    double **weight_tetra;
    memory->allocate(energy_tmp, 3, nktot);
    memory->allocate(weight_tetra, 3, nktot);

    double xk_tmp[3];

    for (auto ik = 0; ik < nklocal; ++ik) {

        auto tmpk = nk_l[ik];
        
        k1 = kpoint->kpoint_irred_all[tmpk][0].knum;  // k index in full grid

        for (s1 = 0; s1 < ns; ++s1) {

            omega1 = dynamical->eval_phonon[k1][s1];

            for (int ib = 0; ib < ns2; ++ib) {

                s2 = ib / ns;
                s3 = ib % ns;

                for (k2 = 0; k2 < nktot; k2++) {

                    for (i = 0; i < 3; ++i) xk_tmp[i] = kpoint->xk[k1][i] - kpoint->xk[k2][i];

                    k3 = kpoint->get_knum(xk_tmp[0], xk_tmp[1], xk_tmp[2]);
                    
                    omega2 = dynamical->eval_phonon[k2][s2];
                    omega3 = dynamical->eval_phonon[k3][s3];

                    energy_tmp[0][k2] = -omega2 + omega3;
                    energy_tmp[1][k2] =  omega2 + omega3;
                    energy_tmp[2][k2] =  omega2 - omega3;
                } 

                for (i = 0; i < 3; ++i) {
                    integration->calc_weight_tetrahedron(nktot,
                                                     kmap_identity,
                                                     weight_tetra[i],
                                                     energy_tmp[i],
                                                     omega1);
                }

                for (auto j = 0; j < num_unipair[ik]; ++j) {
                    auto kp = start_unipair[ik]+j;
                    auto pair = localnk_triplets[kp];
                    k2 = pair.group[0].ks[0];

                    d1[kp][s1][ib] = weight_tetra[0][k2] * static_cast<double>(nktot);
                    d2[kp][s1][ib] = weight_tetra[1][k2] * static_cast<double>(nktot);
                    d3[kp][s1][ib] = weight_tetra[2][k2] * static_cast<double>(nktot);
                }
            } // ib
        } // s1
    }

    memory->deallocate(energy_tmp);
    memory->deallocate(weight_tetra);
    memory->deallocate(kmap_identity);
}


void Iterativebte::calc_Q_from_phi1(double** &n, double** &q1)
{ 
    // Q = sum{q2,q3} [ 2pi/hbar*|v3|^2*{ n1 * n2 * (n3+1) * delta(w1+w2-w3) + 0.5 * n1 * (n2+1) * (n3+1) * delta(w1-w2-w3) }
    int s1, s2, s3;
    double n1, n2, n3;
    double ph1;

    for (auto ik = 0; ik < nklocal; ++ik) {
        for (s1 = 0; s1 < ns; ++s1) { 
            q1[ik][s1] = 0.0;
        }
    }

    for (auto ik = 0; ik < nklocal; ++ik) {

        auto tmpk = nk_l[ik];
        const int k1 = kpoint->kpoint_irred_all[tmpk][0].knum; 

        for (s1 = 0; s1 < ns; ++s1) { 

            for (auto j = 0; j < num_unipair[ik]; ++j) {

                auto kp = start_unipair[ik]+j;
                auto pair = localnk_triplets[kp];
                auto k2 = pair.group[0].ks[0];
                auto k3 = pair.group[0].ks[1];

                auto multi = static_cast<double>(pair.group.size());
                
                for (int ib = 0; ib < ns2; ++ib) {
                    auto s2 = ib / ns;
                    auto s3 = ib % ns;
                    n1 = n[k1][s1];
                    n2 = n[k2][s2];
                    n3 = n[k3][s3];
                    double phi1 = 0.5 * n1 * n2 * (n3 + 1.0) * delta1[kp][s1][ib] 
                                + 0.5 * n1 * n3 * (n2 + 1.0) * delta3[kp][s1][ib]
                                + 0.5 * n1 * (n2 + 1.0) * (n3 + 1.0) * delta2[kp][s1][ib];

                    q1[ik][s1] += v3[kp][s1][ib] * phi1 * multi;
                    
                }
            }

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
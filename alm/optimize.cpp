/*
 optimize.cpp

 Copyright (c) 2014-2018 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include "optimize.h"
#include "files.h"
#include "constants.h"
#include "constraint.h"
#include "error.h"
#include "fcs.h"
#include "input_parser.h"
#include "mathfunctions.h"
#include "memory.h"
#include "symmetry.h"
#include "timer.h"
#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/SparseQR>
#include <Eigen/SparseCholesky>
#include <Eigen/IterativeLinearSolvers>

using namespace ALM_NS;

Optimize::Optimize()
{
    set_default_variables();
}

Optimize::~Optimize()
{
    deallocate_variables();
}

void Optimize::set_default_variables()
{
    params = nullptr;
    cv_l1_alpha = 0.0;
}

void Optimize::deallocate_variables()
{
    if (params) {
        deallocate(params);
    }
}

int Optimize::optimize_main(const Symmetry *symmetry,
                            Constraint *constraint,
                            Fcs *fcs,
                            const int maxorder,
                            const std::string file_prefix,
                            const std::vector<std::string> &str_order,
                            const int verbosity,
                            const DispForceFile &filedata_train,
                            const DispForceFile &filedata_validation,
                            Timer *timer)
{
    timer->start_clock("optimize");

    const auto natmin = symmetry->get_nat_prim();

    const auto ndata_used = filedata_train.nend - filedata_train.nstart + 1
        - filedata_train.skip_e + filedata_train.skip_s;
    const auto ndata_used_validation = filedata_validation.nend - filedata_validation.nstart + 1;
    const auto ntran = symmetry->get_ntran();
    auto info_fitting = 0;
    const auto M = get_number_of_rows_sensing_matrix();
    const auto M_validation = 3 * natmin * ndata_used_validation * ntran;
    size_t N = 0;
    size_t N_new = 0;
    for (auto i = 0; i < maxorder; ++i) {
        N += fcs->get_nequiv()[i].size();
    }

    if (constraint->get_constraint_algebraic()) {
        for (auto i = 0; i < maxorder; ++i) {
            N_new += constraint->get_index_bimap(i).size();
        }
    }

    if (verbosity > 0) {
        std::vector<std::string> str_linearmodel{"least-squares", "elastic-net"};
        std::cout << " OPTIMIZATION\n";
        std::cout << " ============\n\n";
        std::cout << "  LMODEL = " << str_linearmodel[optcontrol.linear_model - 1] << "\n\n";
        if (filedata_train.filename_second.empty()) {
            std::cout << "  Training data file (DFSET) : " << filedata_train.filename << "\n\n";
        } else {
            std::cout << "  Training data file (DFILE) : " << filedata_train.filename << "\n";
            std::cout << "  Training data file (FFILE) : " << filedata_train.filename_second << "\n\n";
        }
        std::cout << "  NSTART = " << filedata_train.nstart << "; NEND = " << filedata_train.nend << '\n';
        if (filedata_train.skip_s < filedata_train.skip_e)
            std::cout << ": SKIP = " << filedata_train.skip_s << "-" <<
                filedata_train.skip_e - 1 << '\n';
        std::cout << "  " << ndata_used
            << " entries will be used for training.\n\n";

        if (optcontrol.cross_validation == -1) {
            std::cout << "  CV = -1 : Manual cross-validation mode is selected\n";
            std::cout << "  Validation data file (DFSET_CV) : " << filedata_validation.filename << "\n\n";
            std::cout << "  NSTART_CV = " << filedata_validation.nstart << "; NEND_CV = " << filedata_validation.nend <<
                std::endl;
            std::cout << "  " << ndata_used_validation
                << " entries will be used for validation." << std::endl << std::endl;
        }

        std::cout << "  Total Number of Parameters : " << N << '\n';
        if (constraint->get_constraint_algebraic()) {
            std::cout << "  Total Number of Free Parameters : " << N_new << '\n';
        }
        std::cout << '\n';
    }

    // Run optimization and obtain force constants

    std::vector<double> fcs_tmp(N, 0.0);
    if (optcontrol.linear_model == 1) {
        // Use ordinary least-squares


        info_fitting = least_squares(maxorder,
                                     N,
                                     N_new,
                                     M,
                                     verbosity,
                                     symmetry,
                                     fcs,
                                     constraint,
                                     fcs_tmp);

    } else if (optcontrol.linear_model == 2) {
        // Use elastic net

        if (!constraint->get_constraint_algebraic()) {
            exit("optimize_main",
                 "Sorry, ICONST = 10 or ICONST = 11 must be used when using elastic net.");
        }

        info_fitting = elastic_net(file_prefix,
                                   maxorder,
                                   N_new,
                                   M,
                                   symmetry,
                                   str_order,
                                   fcs,
                                   constraint,
                                   verbosity,
                                   fcs_tmp);
    }

    if (info_fitting == 0) {
        // I should copy fcs_tmp to parameters in the Fcs class?
        // Copy force constants to public variable "params"
        if (params) {
            deallocate(params);
        }
        allocate(params, N);
        for (auto i = 0; i < N; ++i) params[i] = fcs_tmp[i];

        // Set calculated force constants in FCS class
        fcs->set_forceconstant_cartesian(maxorder,
                                         params);
    }

    fcs_tmp.clear();
    fcs_tmp.shrink_to_fit();

    if (verbosity > 0) {
        std::cout << std::endl;
        timer->print_elapsed();
        std::cout << " -------------------------------------------------------------------" << std::endl;
        std::cout << std::endl;
    }

    timer->stop_clock("optimize");

    return info_fitting;
}

int Optimize::least_squares(const int maxorder,
                            const size_t N,
                            const size_t N_new,
                            const size_t M,
                            const int verbosity,
                            const Symmetry *symmetry,
                            const Fcs *fcs,
                            const Constraint *constraint,
                            std::vector<double> &param_out)
{
    auto info_fitting = 0;

    std::vector<double> amat;
    std::vector<double> bvec;

    if (constraint->get_constraint_algebraic()) {

        // Apply constraints algebraically. (ICONST = 2, 3 is not supported.)
        // SPARSE = 1 is used only when the constraints are considered algebraically.

        // Calculate matrix elements for fitting

        double fnorm;
        const auto nrows = get_number_of_rows_sensing_matrix();
        const unsigned long ncols = static_cast<long>(N_new);

        if (optcontrol.use_sparse_solver) {

            // Use a solver for sparse matrix 
            // (Requires less memory for sparse inputs.)

            SpMat sp_amat(nrows, ncols);
            Eigen::VectorXd sp_bvec(nrows);

            get_matrix_elements_in_sparse_form(maxorder,
                                               sp_amat,
                                               sp_bvec,
                                               u_train,
                                               f_train,
                                               fnorm,
                                               symmetry,
                                               fcs,
                                               constraint);
            if (verbosity > 0) {
                std::cout << "  Now, start fitting ..." << std::endl;
            }

            info_fitting = run_eigen_sparse_solver(sp_amat,
                                                   sp_bvec,
                                                   param_out,
                                                   fnorm,
                                                   maxorder,
                                                   fcs,
                                                   constraint,
                                                   optcontrol.sparsesolver,
                                                   verbosity);

        } else {

            // Use a direct solver for a dense matrix

            amat.resize(nrows * ncols, 0.0);
            bvec.resize(nrows, 0.0);

            get_matrix_elements_algebraic_constraint(maxorder,
                                                     amat,
                                                     bvec,
                                                     u_train,
                                                     f_train,
                                                     fnorm,
                                                     symmetry,
                                                     fcs,
                                                     constraint);

            // Perform fitting with SVD

            info_fitting
                = fit_algebraic_constraints(N_new,
                                            M,
                                            &amat[0],
                                            &bvec[0],
                                            param_out,
                                            fnorm,
                                            maxorder,
                                            fcs,
                                            constraint,
                                            verbosity);
        }

    } else {

        // Apply constraints numerically (ICONST=2 is supported)

        if (optcontrol.use_sparse_solver && verbosity > 0) {
            std::cout << "  WARNING: SPARSE = 1 works only with ICONST = 10 or ICONST = 11." << std::endl;
            std::cout << "  Use a solver for dense matrix." << std::endl;
        }

        get_matrix_elements(maxorder,
                            amat,
                            bvec,
                            u_train,
                            f_train,
                            symmetry,
                            fcs);

        // Perform fitting with SVD or QRD

        assert(!amat.empty());
        assert(!bvec.empty());

        if (constraint->get_exist_constraint()) {
            info_fitting
                = fit_with_constraints(N,
                                       M,
                                       constraint->get_number_of_constraints(),
                                       &amat[0],
                                       &bvec[0],
                                       &param_out[0],
                                       constraint->get_const_mat(),
                                       constraint->get_const_rhs(),
                                       verbosity);
        } else {
            info_fitting
                = fit_without_constraints(N,
                                          M,
                                          &amat[0],
                                          &bvec[0],
                                          &param_out[0],
                                          verbosity);
        }
    }

    return info_fitting;
}


int Optimize::elastic_net(const std::string job_prefix,
                          const int maxorder,
                          const size_t N_new,
                          const size_t M,
                          const Symmetry *symmetry,
                          const std::vector<std::string> &str_order,
                          const Fcs *fcs,
                          Constraint *constraint,
                          const int verbosity,
                          std::vector<double> &param_out)
{
    int info_fitting;

    std::vector<double> param_tmp(N_new, 0.0);

    // Scale displacements if DNORM is not 1 and the data is not standardized.
    const int scale_displacement
        = std::abs(optcontrol.displacement_normalization_factor - 1.0) > eps
        && optcontrol.standardize == 0;

    if (optcontrol.cross_validation == 0) {
        // Calculate force
        // constants at given L1-alpha

        if (scale_displacement) {
            apply_scalers(maxorder, constraint);
        }

        // Optimize with a given L1 coefficient (l1_alpha)
        run_elastic_net_optimization(maxorder,
                                     M,
                                     N_new,
                                     fcs,
                                     symmetry,
                                     constraint,
                                     verbosity,
                                     param_tmp);

        if (verbosity > 0) {
            size_t iparam = 0;
            std::vector<int> nzero_lasso(maxorder);

            for (auto i = 0; i < maxorder; ++i) {
                nzero_lasso[i] = 0;
                for (const auto &it : constraint->get_index_bimap(i)) {
                    const auto inew = it.left + iparam;
                    if (std::abs(param_tmp[inew]) < eps) ++nzero_lasso[i];
                }
                iparam += constraint->get_index_bimap(i).size();
            }

            for (auto order = 0; order < maxorder; ++order) {
                std::cout << "  Number of non-zero " << std::setw(9) << str_order[order] << " FCs : "
                    << constraint->get_index_bimap(order).size() - nzero_lasso[order] << std::endl;
            }
            std::cout << std::endl;
        }

        // Scale back force constants

        if (scale_displacement) {
            apply_scaler_force_constants(maxorder,
                                         optcontrol.displacement_normalization_factor,
                                         constraint,
                                         param_tmp);
            finalize_scalers(maxorder, constraint);
        }

        recover_original_forceconstants(maxorder,
                                        param_tmp,
                                        param_out,
                                        fcs->get_nequiv(),
                                        constraint);
        info_fitting = 0;

    } else {
        // Run cross validation (manually or automatically) to
        // get L1 alpha to give minimum CV score
        if (scale_displacement) {
            apply_scalers(maxorder, constraint);
        }

        // cv_l1_alpha is a private variable of Optimize class.
        cv_l1_alpha = run_elastic_net_crossvalidation(job_prefix,
                                                      maxorder,
                                                      fcs,
                                                      symmetry,
                                                      constraint,
                                                      verbosity);
        if (scale_displacement) {
            finalize_scalers(maxorder, constraint);
        }

        info_fitting = 1;
    }

    return info_fitting;
}

double Optimize::run_elastic_net_crossvalidation(const std::string job_prefix,
                                                 const int maxorder,
                                                 const Fcs *fcs,
                                                 const Symmetry *symmetry,
                                                 const Constraint *constraint,
                                                 const int verbosity)
{
    // Cross-validation mode:
    // Returns alpha giving minimum CV score


    if (verbosity > 0) {
        std::cout << "  Elastic-net cross-validation with the following parameters:" << std::endl;
        std::cout << "   L1_RATIO = " << optcontrol.l1_ratio << std::endl;
        std::cout << "   CV = " << std::setw(15) << optcontrol.cross_validation << std::endl;
        if (optcontrol.l1_alpha_min > 0) {
            std::cout << "   CV_MINALPHA = " << std::setw(15) << optcontrol.l1_alpha_min;
        } else {
            std::cout << "   CV_MINALPHA = CV_MAXALPHA*1e-6 ";
        }
        if (optcontrol.l1_alpha_max > 0) {
            std::cout << "  CV_MAXALPHA = " << std::setw(15) << optcontrol.l1_alpha_max << std::endl;
        } else {
            std::cout << " CV_MAXALPHA = (Use recommended value)" << std::endl;
        }
        std::cout << "   CV_NALPHA = " << std::setw(5) << optcontrol.num_l1_alpha << std::endl;
        std::cout << "   CONV_TOL = " << std::setw(15) << optcontrol.tolerance_iteration << std::endl;
        std::cout << "   MAXITER = " << std::setw(5) << optcontrol.maxnum_iteration << std::endl;
        std::cout << "   ENET_DNORM = " << std::setw(15) << optcontrol.displacement_normalization_factor << std::endl;
        std::cout << std::endl;

        if (optcontrol.standardize) {
            std::cout << "  STANDARDIZE = 1 : Standardization will be performed for matrix A and vector b." << std::
                endl;
            std::cout << "                    The ENET_DNORM-tag will be neglected." << std::endl;
        } else {
            std::cout << "  STANDARDIZE = 0 : No standardization of matrix A and vector b." << std::endl;
            std::cout << "                    Columns of matrix A will be scaled by the ENET_DNORM value." << std::
                endl;
        }
        std::cout << std::endl;

        if (optcontrol.cross_validation == -1) {
            std::cout << "  CV = -1: Manual CV mode." << std::endl;
            std::cout << "           Validation data is read from DFSET_CV" << std::endl;
        } else if (optcontrol.cross_validation > 0) {
            std::cout << "  CV > 0: Automatic CV mode." << std::endl;
        } else {
            exit("run_elastic_net_crossvalidation",
                 "This cannot happen.");
        }
        std::cout << std::endl;
    }


    // Returns alpha at minimum CV
    if (optcontrol.cross_validation == -1) {
        return run_enetcv_manual(job_prefix,
                                 maxorder,
                                 fcs,
                                 symmetry,
                                 constraint,
                                 verbosity);
    } else {
        return run_enetcv_auto(job_prefix,
                               maxorder,
                               fcs,
                               symmetry,
                               constraint,
                               verbosity);
    }
}

double Optimize::run_enetcv_manual(const std::string job_prefix,
                                   const int maxorder,
                                   const Fcs *fcs,
                                   const Symmetry *symmetry,
                                   const Constraint *constraint,
                                   const int verbosity)
{
    // Manual CV mode where the test data is read from the user-defined file.
    // Indeed, the test data is already read in the input_parser and stored in u_validation and f_validation.

    std::vector<double> amat_1D, amat_1D_validation;
    std::vector<double> bvec, bvec_validation;
    std::vector<double> alphas, training_error, validation_error;
    std::vector<std::vector<int>> nonzeros;
    double fnorm, fnorm_validation;

    size_t N_new = 0;
    if (constraint->get_constraint_algebraic()) {
        for (auto i = 0; i < maxorder; ++i) {
            N_new += constraint->get_index_bimap(i).size();
        }
    }

    get_matrix_elements_algebraic_constraint(maxorder,
                                             amat_1D,
                                             bvec,
                                             u_train,
                                             f_train,
                                             fnorm,
                                             symmetry,
                                             fcs,
                                             constraint);

    get_matrix_elements_algebraic_constraint(maxorder,
                                             amat_1D_validation,
                                             bvec_validation,
                                             u_validation,
                                             f_validation,
                                             fnorm_validation,
                                             symmetry,
                                             fcs,
                                             constraint);

    Eigen::MatrixXd A = Eigen::Map<Eigen::MatrixXd>(&amat_1D[0], amat_1D.size() / N_new, N_new);
    Eigen::VectorXd b = Eigen::Map<Eigen::VectorXd>(&bvec[0], bvec.size());
    Eigen::MatrixXd A_validation = Eigen::Map<Eigen::MatrixXd>(&amat_1D_validation[0],
                                                               amat_1D_validation.size() / N_new, N_new);
    Eigen::VectorXd b_validation = Eigen::Map<Eigen::VectorXd>(&bvec_validation[0], bvec_validation.size());

    const auto estimated_max_alpha = get_estimated_max_alpha(A, b);

    if (verbosity > 0) {
        std::cout << "  Recommended CV_MAXALPHA = "
            << estimated_max_alpha
            << std::endl << std::endl;
    }

    const auto file_coef = job_prefix + ".solution_path";
    const auto file_cv = job_prefix + ".enet_cv";

    compute_alphas(optcontrol.l1_alpha_max,
                   optcontrol.l1_alpha_min,
                   optcontrol.num_l1_alpha,
                   alphas);

    run_enet_solution_path(maxorder, A, b, A_validation, b_validation,
                           fnorm, fnorm_validation,
                           file_coef, verbosity,
                           constraint,
                           alphas,
                           training_error, validation_error, nonzeros);

    write_cvresult_to_file(file_cv,
                           alphas,
                           training_error,
                           validation_error,
                           nonzeros);

    const auto ialpha = get_ialpha_at_minimum_validation_error(validation_error);

    if (verbosity > 0) {
        std::cout << "  The manual CV has been done." << std::endl;
        std::cout << "  Minimum validation error at alpha = "
            << alphas[ialpha] << std::endl;
        std::cout << "  The CV result is saved in " << file_cv << std::endl;
    }

    return alphas[ialpha];
}

double Optimize::run_enetcv_auto(const std::string job_prefix,
                                 const int maxorder,
                                 const Fcs *fcs,
                                 const Symmetry *symmetry,
                                 const Constraint *constraint,
                                 const int verbosity)
{
    // Automatic CV mode.

    size_t N_new = 0;
    if (constraint->get_constraint_algebraic()) {
        for (auto i = 0; i < maxorder; ++i) {
            N_new += constraint->get_index_bimap(i).size();
        }
    }

    const auto nstructures = u_train.size();
    const auto nsets = optcontrol.cross_validation;

    if (nsets > nstructures) {
        exit("run_elastic_net_crossvalidation",
             "The input CV is larger than the total number of training data.");
    }

    std::vector<int> ndata_block(nsets, nstructures / nsets);
    for (auto iset = 0; iset < nsets; ++iset) {
        if (nstructures - nsets * (nstructures / nsets) > iset) {
            ++ndata_block[iset];
        }
    }

    std::vector<std::vector<double>> u_train_tmp, u_validation_tmp;
    std::vector<std::vector<double>> f_train_tmp, f_validation_tmp;

    std::vector<double> amat_1D, amat_1D_validation;
    std::vector<double> bvec, bvec_validation;
    std::vector<double> alphas, training_error, validation_error;
    std::vector<std::vector<int>> nonzeros;
    std::vector<std::vector<double>> training_error_accum, validation_error_accum;
    double fnorm, fnorm_validation, estimated_max_alpha;

    auto ishift = 0;

    if (verbosity > 0) {
        std::cout << "  Start " << nsets << "-fold CV with "
            << u_train.size() << " Datasets" << std::endl;
        std::cout << std::endl;
    }

    if (!(optcontrol.l1_alpha_max > 0)) {
        estimated_max_alpha = 0;
        for (auto iset = 0; iset < nsets; ++iset) {
            const auto istart_validation = ishift;
            const auto iend_validation = istart_validation + ndata_block[iset];

            u_train_tmp.clear();
            f_train_tmp.clear();
            u_validation_tmp.clear();
            f_validation_tmp.clear();

            for (auto idata = 0; idata < nstructures; ++idata) {
                if (idata >= istart_validation && idata < iend_validation) {
                    u_validation_tmp.emplace_back(u_train[idata]);
                    f_validation_tmp.emplace_back(f_train[idata]);
                } else {
                    u_train_tmp.emplace_back(u_train[idata]);
                    f_train_tmp.emplace_back(f_train[idata]);
                }
            }
            ishift += ndata_block[iset];

            get_matrix_elements_algebraic_constraint(maxorder,
                                                     amat_1D,
                                                     bvec,
                                                     u_train_tmp,
                                                     f_train_tmp,
                                                     fnorm,
                                                     symmetry,
                                                     fcs,
                                                     constraint);

            Eigen::MatrixXd A = Eigen::Map<Eigen::MatrixXd>(&amat_1D[0], amat_1D.size() / N_new, N_new);
            Eigen::VectorXd b = Eigen::Map<Eigen::VectorXd>(&bvec[0], bvec.size());
            const auto this_estimated_max_alpha = get_estimated_max_alpha(A, b);

            if (verbosity > 0) {
                std::cout << "  Recommended CV_MAXALPHA (" << std::setw(3)
                    << iset + 1 << ") = "
                    << this_estimated_max_alpha << std::endl;
            }

            if (this_estimated_max_alpha > estimated_max_alpha) {
                estimated_max_alpha = this_estimated_max_alpha;
            }
        }
        ishift = 0;
    }

    for (auto iset = 0; iset < nsets; ++iset) {

        if (verbosity > 0) {
            std::cout << std::endl;
            std::cout << "  SET : " << std::setw(3) << iset + 1 << std::endl;
        }
        const auto istart_validation = ishift;
        const auto iend_validation = istart_validation + ndata_block[iset];

        u_train_tmp.clear();
        f_train_tmp.clear();
        u_validation_tmp.clear();
        f_validation_tmp.clear();

        for (auto idata = 0; idata < nstructures; ++idata) {
            if (idata >= istart_validation && idata < iend_validation) {
                u_validation_tmp.emplace_back(u_train[idata]);
                f_validation_tmp.emplace_back(f_train[idata]);
            } else {
                u_train_tmp.emplace_back(u_train[idata]);
                f_train_tmp.emplace_back(f_train[idata]);
            }
        }
        ishift += ndata_block[iset];

        get_matrix_elements_algebraic_constraint(maxorder,
                                                 amat_1D,
                                                 bvec,
                                                 u_train_tmp,
                                                 f_train_tmp,
                                                 fnorm,
                                                 symmetry,
                                                 fcs,
                                                 constraint);

        get_matrix_elements_algebraic_constraint(maxorder,
                                                 amat_1D_validation,
                                                 bvec_validation,
                                                 u_validation_tmp,
                                                 f_validation_tmp,
                                                 fnorm_validation,
                                                 symmetry,
                                                 fcs,
                                                 constraint);

        Eigen::MatrixXd A = Eigen::Map<Eigen::MatrixXd>(&amat_1D[0], amat_1D.size() / N_new, N_new);
        Eigen::VectorXd b = Eigen::Map<Eigen::VectorXd>(&bvec[0], bvec.size());
        Eigen::MatrixXd A_validation = Eigen::Map<Eigen::MatrixXd>(&amat_1D_validation[0],
                                                                   amat_1D_validation.size() / N_new, N_new);
        Eigen::VectorXd b_validation = Eigen::Map<Eigen::VectorXd>(&bvec_validation[0], bvec_validation.size());

        if (verbosity > 0) {
            std::cout << "  Recommended CV_MAXALPHA = "
                << get_estimated_max_alpha(A, b)
                << std::endl << std::endl;
        }

        const auto file_coef = job_prefix + ".solution_path" + std::to_string(iset + 1);
        const auto file_cv = job_prefix + ".enet_cvset" + std::to_string(iset + 1);

        if (optcontrol.l1_alpha_max > 0) {
            compute_alphas(optcontrol.l1_alpha_max,
                           optcontrol.l1_alpha_min,
                           optcontrol.num_l1_alpha,
                           alphas);
        } else {
            if (optcontrol.l1_alpha_max > 0) {
                compute_alphas(estimated_max_alpha,
                               optcontrol.l1_alpha_min,
                               optcontrol.num_l1_alpha,
                               alphas);
            } else {
                compute_alphas(estimated_max_alpha,
                               estimated_max_alpha * 1e-6,
                               optcontrol.num_l1_alpha,
                               alphas);
            }
        }

        run_enet_solution_path(maxorder, A, b, A_validation, b_validation,
                               fnorm, fnorm_validation,
                               file_coef, verbosity,
                               constraint,
                               alphas,
                               training_error, validation_error, nonzeros);

        if (job_prefix != "") {
            write_cvresult_to_file(file_cv,
                                   alphas,
                                   training_error,
                                   validation_error,
                                   nonzeros);
        }

        if (verbosity > 0) {
            std::cout << "  SET " << std::setw(3) << iset + 1 << " has been finished." << std::endl;
            std::cout << "  Minimum validation error at alpha = "
                << alphas[get_ialpha_at_minimum_validation_error(validation_error)] << std::endl;
            if (job_prefix != "") {
                std::cout << "  The CV result is saved in " << file_cv << std::endl << std::endl;
            }
            std::cout << "  ---------------------------------------------------" << std::endl;
        }

        training_error_accum.emplace_back(training_error);
        validation_error_accum.emplace_back(validation_error);
    }

    std::vector<double> terr_mean, terr_std;
    std::vector<double> verr_mean, verr_std;

    const auto nalphas = alphas.size();

    terr_mean.resize(nalphas);
    terr_std.resize(nalphas);
    verr_mean.resize(nalphas);
    verr_std.resize(nalphas);

    set_errors_of_cvscore(terr_mean, terr_std, verr_mean, verr_std,
                          training_error_accum, validation_error_accum);
    const auto ialpha_minimum = get_ialpha_at_minimum_validation_error(verr_mean);

    if (job_prefix != "") {
        const auto file_cvscore = job_prefix + ".cvscore";
        write_cvscore_to_file(file_cvscore,
                              alphas,
                              terr_mean,
                              terr_std,
                              verr_mean,
                              verr_std,
                              ialpha_minimum,
                              nsets);

        if (verbosity > 0) {
            std::cout << " Average and standard deviation of the CV error are" << std::endl;
            std::cout << " saved in " << file_cvscore << std::endl;
            std::cout << " Minimum CVSCORE at alpha = " << alphas[ialpha_minimum] << std::endl;
            std::cout << std::endl;
        }
    }

    return alphas[ialpha_minimum];
}

void Optimize::write_cvresult_to_file(const std::string file_out,
                                      const std::vector<double> &alphas,
                                      const std::vector<double> &training_error,
                                      const std::vector<double> &validation_error,
                                      const std::vector<std::vector<int>> &nonzeros) const
{
    std::ofstream ofs_cv;
    ofs_cv.open(file_out.c_str(), std::ios::out);
    ofs_cv << "# Algorithm : Coordinate descent" << std::endl;
    ofs_cv << "# L1_RATIO = " << optcontrol.l1_ratio << std::endl;
    ofs_cv << "# ENET_DNORM = " << std::setw(15) << optcontrol.displacement_normalization_factor << std::endl;
    ofs_cv << "# STANDARDIZE = " << optcontrol.standardize << std::endl;
    ofs_cv << "# CONV_TOL = " << std::setw(15) << optcontrol.tolerance_iteration << std::endl;
    ofs_cv << "# L1 ALPHA, Fitting error, Validation error, Num. zero IFCs (2nd, 3rd, ...) " << std::endl;

    const auto maxorder = nonzeros[0].size();
    for (auto ialpha = 0; ialpha < alphas.size(); ++ialpha) {
        ofs_cv << std::setw(15) << alphas[ialpha];
        ofs_cv << std::setw(15) << training_error[ialpha];
        ofs_cv << std::setw(15) << validation_error[ialpha];
        for (auto i = 0; i < maxorder; ++i) {
            ofs_cv << std::setw(10) << nonzeros[ialpha][i];
        }
        ofs_cv << std::endl;
    }
    ofs_cv.close();
}

void Optimize::write_cvscore_to_file(const std::string file_out,
                                     const std::vector<double> &alphas,
                                     const std::vector<double> &terr_mean,
                                     const std::vector<double> &terr_std,
                                     const std::vector<double> &verr_mean,
                                     const std::vector<double> &verr_std,
                                     const int ialpha_minimum,
                                     const size_t nsets) const
{
    const auto nalphas = alphas.size();

    std::ofstream ofs_cv;
    ofs_cv.open(file_out.c_str(), std::ios::out);
    ofs_cv << "# Algorithm : Coordinate descent" << std::endl;
    ofs_cv << "# L1_RATIO = " << optcontrol.l1_ratio << std::endl;
    ofs_cv << "# ENET_DNORM = " << std::setw(15) << optcontrol.displacement_normalization_factor << std::endl;
    ofs_cv << "# STANDARDIZE = " << optcontrol.standardize << std::endl;
    ofs_cv << "# CONV_TOL = " << std::setw(15) << optcontrol.tolerance_iteration << std::endl;
    ofs_cv << "# " << nsets << "-fold cross-validation scores" << std::endl;
    ofs_cv << "# L1 ALPHA, Fitting error (mean, std), Validation error (mean, std) " << std::endl;

    for (size_t ialpha = 0; ialpha < nalphas; ++ialpha) {
        ofs_cv << std::setw(15) << alphas[ialpha];
        ofs_cv << std::setw(15) << terr_mean[ialpha];
        ofs_cv << std::setw(15) << terr_std[ialpha];
        ofs_cv << std::setw(15) << verr_mean[ialpha];
        ofs_cv << std::setw(15) << verr_std[ialpha];
        ofs_cv << std::endl;
    }

    ofs_cv << "# Minimum CVSCORE at alpha = " << alphas[ialpha_minimum] << std::endl;
    ofs_cv.close();
}

void Optimize::set_errors_of_cvscore(std::vector<double> &terr_mean,
                                     std::vector<double> &terr_std,
                                     std::vector<double> &verr_mean,
                                     std::vector<double> &verr_std,
                                     const std::vector<std::vector<double>> &training_error_accum,
                                     const std::vector<std::vector<double>> &validation_error_accum) const
{
    const auto nsets = training_error_accum.size();
    const auto nalphas = terr_mean.size();

    double sum_t, sum2_t;
    double sum_v, sum2_v;
    const auto factor = 1.0 / static_cast<double>(nsets);

    for (size_t ialpha = 0; ialpha < nalphas; ++ialpha) {
        sum_t = 0.0;
        sum2_t = 0.0;
        sum_v = 0.0;
        sum2_v = 0.0;
        for (size_t iset = 0; iset < nsets; ++iset) {
            sum_t += training_error_accum[iset][ialpha];
            sum2_t += training_error_accum[iset][ialpha]
                * training_error_accum[iset][ialpha];
            sum_v += validation_error_accum[iset][ialpha];
            sum2_v += validation_error_accum[iset][ialpha]
                * validation_error_accum[iset][ialpha];
        }
        sum_t *= factor;
        sum2_t *= factor;
        sum_v *= factor;
        sum2_v *= factor;
        terr_mean[ialpha] = sum_t;
        terr_std[ialpha] = std::sqrt(sum2_t - sum_t * sum_t);
        verr_mean[ialpha] = sum_v;
        verr_std[ialpha] = std::sqrt(sum2_v - sum_v * sum_v);
    }
}

int Optimize::get_ialpha_at_minimum_validation_error(const std::vector<double> &validation_error) const
{
    return std::distance(validation_error.begin(),
                         std::min_element(validation_error.begin(),
                                          validation_error.end()));
}

void Optimize::run_enet_solution_path(const int maxorder,
                                      Eigen::MatrixXd &A,
                                      Eigen::VectorXd &b,
                                      Eigen::MatrixXd &A_validation,
                                      Eigen::VectorXd &b_validation,
                                      const double fnorm,
                                      const double fnorm_validation,
                                      const std::string file_coef,
                                      const int verbosity,
                                      const Constraint *constraint,
                                      const std::vector<double> &alphas,
                                      std::vector<double> &training_error,
                                      std::vector<double> &validation_error,
                                      std::vector<std::vector<int>> &nonzeros) const
{
    int initialize_mode;

    std::ofstream ofs_coef;

    std::vector<double> params_tmp;
    std::vector<int> nzero_lasso(maxorder);

    bool *has_prod;

    Eigen::MatrixXd Prod;
    Eigen::VectorXd grad0, grad, x;
    Eigen::VectorXd scale_beta, scale_beta_enet;
    Eigen::VectorXd factor_std;
    Eigen::VectorXd fdiff, fdiff_validation;
    Eigen::VectorXd mean, dev;

    size_t N_new = A.cols();
    size_t M = A.rows();
    size_t M_validation = A_validation.rows();

    Prod.setZero(N_new, N_new);
    grad0.resize(N_new);
    grad.resize(N_new);
    x.setZero(N_new);
    scale_beta.resize(N_new);
    scale_beta_enet.resize(N_new);
    factor_std.resize(N_new);
    fdiff.resize(M);
    fdiff_validation.resize(M_validation);

    allocate(has_prod, N_new);

    for (size_t i = 0; i < N_new; ++i) {
        has_prod[i] = false;
    }

    if (optcontrol.save_solution_path) {
        ofs_coef.open(file_coef.c_str(), std::ios::out);
        ofs_coef << "# L1 ALPHA, coefficients" << std::endl;
        params_tmp.resize(N_new);
    }

    if (optcontrol.standardize) {
        get_standardizer(A, mean, dev, factor_std, scale_beta);
        apply_standardizer(A, mean, dev);
        apply_standardizer(A_validation, mean, dev);
    } else {
        get_standardizer(A, mean, dev, factor_std, scale_beta);
    }

    training_error.clear();
    validation_error.clear();
    nonzeros.clear();

    // Start iteration

    grad0 = A.transpose() * b;
    grad = grad0;

    if (verbosity == 1) std::cout << std::setw(3);

    for (size_t ialpha = 0; ialpha < alphas.size(); ++ialpha) {

        const auto l1_alpha = alphas[ialpha];

        if (ialpha == 0) {
            initialize_mode = 0;
        } else {
            initialize_mode = 1;
        }

#pragma omp parallel for
        for (auto i = 0; i < N_new; ++i) {
            scale_beta_enet(i) = 1.0 / (1.0 / scale_beta(i) + (1.0 - optcontrol.l1_ratio) * l1_alpha);
        }

        coordinate_descent(M, N_new, l1_alpha,
                           initialize_mode,
                           x, A, b, grad0, has_prod, Prod, grad, fnorm,
                           scale_beta_enet,
                           verbosity);

        double correction_intercept = 0.0;
        for (size_t i = 0; i < N_new; ++i) {
            correction_intercept += x(i) * mean(i) * factor_std(i);
        }
        fdiff = A * x - b + correction_intercept * Eigen::VectorXd::Ones(M);
        fdiff_validation = A_validation * x - b_validation + correction_intercept * Eigen::VectorXd::Ones(M_validation);
        const auto res1 = fdiff.dot(fdiff) / (fnorm * fnorm);
        const auto res2 = fdiff_validation.dot(fdiff_validation) / (fnorm_validation * fnorm_validation);

        get_number_of_zero_coefs(maxorder,
                                 constraint,
                                 x,
                                 nzero_lasso);

        training_error.push_back(std::sqrt(res1));
        validation_error.push_back(std::sqrt(res2));
        nonzeros.push_back(nzero_lasso);

        if (optcontrol.save_solution_path) {
            ofs_coef << std::setw(15) << l1_alpha;

            for (auto i = 0; i < N_new; ++i) params_tmp[i] = x[i];

            apply_scaler_force_constants(maxorder,
                                         optcontrol.displacement_normalization_factor,
                                         constraint,
                                         params_tmp);

            for (auto i = 0; i < N_new; ++i) {
                ofs_coef << std::setw(15) << params_tmp[i];
            }
            ofs_coef << std::endl;
        }

        if (verbosity == 1) {

            std::cout << '.' << std::flush;
            if (ialpha % 25 == 24) {
                std::cout << std::endl;
                std::cout << std::setw(3);
            }
        }
    }

    if (verbosity == 1) std::cout << std::endl;

    if (optcontrol.save_solution_path) {
        ofs_coef.close();
        params_tmp.clear();
        params_tmp.shrink_to_fit();
    }
    deallocate(has_prod);
}

void Optimize::compute_alphas(const double l1_alpha_max,
                              const double l1_alpha_min,
                              const int num_l1_alpha,
                              std::vector<double> &alphas) const
{
    alphas.resize(num_l1_alpha);

    for (auto ialpha = 0; ialpha < num_l1_alpha; ++ialpha) {

        const auto l1_alpha = l1_alpha_min
            * std::pow(l1_alpha_max / l1_alpha_min,
                       static_cast<double>(num_l1_alpha - ialpha - 1) /
                       static_cast<double>(num_l1_alpha));

        alphas[ialpha] = l1_alpha;
    }
}

void Optimize::run_elastic_net_optimization(const int maxorder,
                                            const size_t M,
                                            const size_t N_new,
                                            const Fcs *fcs,
                                            const Symmetry *symmetry,
                                            const Constraint *constraint,
                                            const int verbosity,
                                            std::vector<double> &param_out) const
{
    // Start Lasso optimization
    int i;
    bool *has_prod;
    double fnorm;

    Eigen::MatrixXd A, Prod;
    Eigen::VectorXd b, grad0, grad, x;
    Eigen::VectorXd scale_beta, factor_std;
    Eigen::VectorXd fdiff;
    Eigen::VectorXd mean, dev;

    std::vector<double> amat_1D;
    std::vector<double> bvec;

    get_matrix_elements_algebraic_constraint(maxorder,
                                             amat_1D,
                                             bvec,
                                             u_train,
                                             f_train,
                                             fnorm,
                                             symmetry,
                                             fcs,
                                             constraint);

    // Coordinate descent

    A = Eigen::Map<Eigen::MatrixXd>(&amat_1D[0], M, N_new);
    b = Eigen::Map<Eigen::VectorXd>(&bvec[0], M);

    Prod.setZero(N_new, N_new);
    grad0.resize(N_new);
    grad.resize(N_new);
    x.setZero(N_new);
    scale_beta.resize(N_new);
    factor_std.resize(N_new);
    fdiff.resize(M);

    allocate(has_prod, N_new);

    for (i = 0; i < N_new; ++i) {
        has_prod[i] = false;
    }

    if (verbosity > 0) {
        std::cout << "  Elastic-net minimization with the following parameters:" << std::endl;
        std::cout << "   L1_RATIO = " << optcontrol.l1_ratio << std::endl;
        std::cout << "   L1_ALPHA = " << std::setw(15) << optcontrol.l1_alpha << std::endl;
        std::cout << "   CONV_TOL = " << std::setw(15) << optcontrol.tolerance_iteration << std::endl;
        std::cout << "   MAXITER = " << std::setw(5) << optcontrol.maxnum_iteration << std::endl;
        std::cout << "   ENET_DNORM = " << std::setw(15) << optcontrol.displacement_normalization_factor << std::endl;

        std::cout << std::endl;
        if (optcontrol.standardize) {
            std::cout << " STANDARDIZE = 1 : Standardization will be performed for matrix A and vector b." << std::endl;
            std::cout << "                   The ENET_DNORM-tag will be neglected." << std::endl;
        } else {
            std::cout << " STANDARDIZE = 0 : No standardization of matrix A and vector b." << std::endl;
            std::cout << "                   Columns of matrix A will be scaled by the ENET_DNORM value." << std::endl;
        }
    }

    // Standardize if necessary

    if (optcontrol.standardize) {
        get_standardizer(A, mean, dev, factor_std, scale_beta);
        apply_standardizer(A, mean, dev);
    } else {
        get_standardizer(A, mean, dev, factor_std, scale_beta);
    }

    grad0 = A.transpose() * b;
    grad = grad0;

#pragma omp parallel for
    for (i = 0; i < N_new; ++i) {
        scale_beta(i) = 1.0 / (1.0 / scale_beta(i) + (1.0 - optcontrol.l1_ratio) * optcontrol.l1_alpha);
    }

    // Coordinate Descent Method
    coordinate_descent(M, N_new, optcontrol.l1_alpha,
                       0,
                       x, A, b, grad0, has_prod, Prod, grad, fnorm,
                       scale_beta,
                       verbosity);

    for (i = 0; i < N_new; ++i) {
        param_out[i] = x[i] * factor_std[i];
    }

    if (verbosity > 0) {
        double correction_intercept = 0.0;
        for (size_t i = 0; i < N_new; ++i) {
            correction_intercept += x(i) * mean(i) * factor_std(i);
        }
        fdiff = A * x - b + correction_intercept * Eigen::VectorXd::Ones(M);
        const auto res1 = fdiff.dot(fdiff) / (fnorm * fnorm);
        std::cout << "  RESIDUAL (%): " << std::sqrt(res1) * 100.0 << std::endl;
    }

    deallocate(has_prod);

    if (optcontrol.debiase_after_l1opt) {
        run_least_squares_with_nonzero_coefs(A, b,
                                             factor_std,
                                             param_out,
                                             verbosity);
    }
}

void Optimize::run_least_squares_with_nonzero_coefs(const Eigen::MatrixXd &A_in,
                                                    const Eigen::VectorXd &b_in,
                                                    const Eigen::VectorXd &factor_std,
                                                    std::vector<double> &params_inout,
                                                    const int verbosity) const
{
    // Perform OLS fitting to the features selected by LASSO for reducing the bias.

    if (verbosity > 0) {
        std::cout << " DEBIAS_OLS = 1: Attempt to reduce the bias of LASSO by performing OLS fitting" << std::endl;
        std::cout << "                 with features selected by LASSO." << std::endl;
    }

    const auto N_new = A_in.cols();
    const auto M = A_in.rows();

    std::vector<int> nonzero_index, zero_index;

    for (auto i = 0; i < N_new; ++i) {
        if (std::abs(params_inout[i]) >= eps) {
            nonzero_index.push_back(i);
        } else {
            zero_index.push_back(i);
        }
    }

    const auto N_nonzero = nonzero_index.size();
    Eigen::MatrixXd A_nonzero(M, N_nonzero);

    for (auto i = 0; i < N_nonzero; ++i) {
        A_nonzero.col(i) = A_in.col(nonzero_index[i]);
    }
    Eigen::VectorXd x_nonzero = A_nonzero.colPivHouseholderQr().solve(b_in);

    for (auto i = 0; i < N_new; ++i) params_inout[i] = 0.0;
    for (auto i = 0; i < N_nonzero; ++i) {
        params_inout[nonzero_index[i]] = x_nonzero[i] * factor_std[nonzero_index[i]];
    }
}

void Optimize::get_number_of_zero_coefs(const int maxorder,
                                        const Constraint *constraint,
                                        const Eigen::VectorXd &x,
                                        std::vector<int> &nzeros) const
{
    // Count the number of zero parameters
    size_t iparam = 0;
    nzeros.resize(maxorder);
    for (auto i = 0; i < maxorder; ++i) {
        nzeros[i] = 0;
        for (const auto &it : constraint->get_index_bimap(i)) {
            const auto inew = it.left + iparam;
            if (std::abs(x[inew]) < eps) ++nzeros[i];

        }
        iparam += constraint->get_index_bimap(i).size();
    }
}


void Optimize::get_standardizer(const Eigen::MatrixXd &Amat,
                                Eigen::VectorXd &mean,
                                Eigen::VectorXd &dev,
                                Eigen::VectorXd &factor_std,
                                Eigen::VectorXd &scale_beta) const
{
    const auto nrows = Amat.rows();
    const auto ncols = Amat.cols();

    if (mean.size() != ncols) mean.resize(ncols);
    if (dev.size() != ncols) dev.resize(ncols);
    if (factor_std.size() != ncols) factor_std.resize(ncols);
    if (scale_beta.size() != ncols) scale_beta.resize(ncols);

    const auto inv_nrows = 1.0 / static_cast<double>(nrows);
    double sum1, sum2;

    if (optcontrol.standardize) {
        for (auto j = 0; j < ncols; ++j) {
            sum1 = Amat.col(j).sum() * inv_nrows;
            sum2 = Amat.col(j).dot(Amat.col(j)) * inv_nrows;
            mean(j) = sum1;
            dev(j) = std::sqrt(sum2 - sum1 * sum1);
            factor_std(j) = 1.0 / dev(j);
            scale_beta(j) = 1.0;
        }
    } else {
        for (auto j = 0; j < ncols; ++j) {
            sum2 = Amat.col(j).dot(Amat.col(j)) * inv_nrows;
            mean(j) = 0.0;
            dev(j) = 1.0;
            factor_std(j) = 1.0;
            scale_beta(j) = 1.0 / sum2;
        }
    }
}

void Optimize::apply_standardizer(Eigen::MatrixXd &Amat,
                                  const Eigen::VectorXd &mean,
                                  const Eigen::VectorXd &dev) const
{
    const auto ncols = Amat.cols();
    const auto nrows = Amat.rows();
    if (mean.size() != ncols || dev.size() != ncols) {
        exit("apply_standardizer", "The number of colums is inconsistent.");
    }

    for (auto i = 0; i < nrows; ++i) {
        for (auto j = 0; j < ncols; ++j) {
            Amat(i, j) = (Amat(i, j) - mean(j)) / dev(j);
        }
    }
}

double Optimize::get_estimated_max_alpha(const Eigen::MatrixXd &Amat,
                                         const Eigen::VectorXd &bvec) const
{
    const auto ncols = Amat.cols();
    const auto nrows = Amat.rows();
    Eigen::MatrixXd C = Amat;

    Eigen::VectorXd mean = Eigen::VectorXd::Zero(Amat.cols());
    Eigen::VectorXd dev = Eigen::VectorXd::Ones(Amat.cols());

    if (optcontrol.standardize) {
        Eigen::VectorXd factor_std, scale_beta;
        factor_std.resize(Amat.cols());
        scale_beta.resize(Amat.cols());
        get_standardizer(Amat, mean, dev, factor_std, scale_beta);
    }

    for (auto i = 0; i < nrows; ++i) {
        for (auto j = 0; j < ncols; ++j) {
            C(i, j) = (C(i, j) - mean(j)) / dev(j);
        }
    }

    C = C.transpose() * bvec;
    auto lambda_max = 0.0;

    for (auto i = 0; i < ncols; ++i) {
        lambda_max = std::max<double>(lambda_max, std::abs(C(i)));
    }
    lambda_max /= static_cast<double>(nrows);

    return lambda_max;
}

void Optimize::apply_scaler_displacement(std::vector<std::vector<double>> &u_inout,
                                         const double normalization_factor,
                                         const bool scale_back) const
{
    const auto nrows = u_inout.size();
    const auto ncols = u_inout[0].size();

    if (scale_back) {
        for (auto i = 0; i < nrows; ++i) {
            for (auto j = 0; j < ncols; ++j) {
                u_inout[i][j] *= normalization_factor;
            }
        }
    } else {
        const auto inv_scale_factor = 1.0 / normalization_factor;
        for (auto i = 0; i < nrows; ++i) {
            for (auto j = 0; j < ncols; ++j) {
                u_inout[i][j] *= inv_scale_factor;
            }
        }
    }
}

void Optimize::apply_scaler_constraint(const int maxorder,
                                       const double normalization_factor,
                                       Constraint *constraint,
                                       const bool scale_back) const
{
    if (scale_back) {
        for (auto i = 0; i < maxorder; ++i) {
            const auto scale_factor = 1.0 / std::pow(normalization_factor, i + 1);
            for (auto j = 0; j < constraint->get_const_fix(i).size(); ++j) {
                const auto scaled_val = constraint->get_const_fix(i)[j].val_to_fix * scale_factor;
                constraint->set_const_fix_val_to_fix(i, j, scaled_val);
            }
        }
    } else {
        for (auto i = 0; i < maxorder; ++i) {
            const auto scale_factor = std::pow(normalization_factor, i + 1);
            for (auto j = 0; j < constraint->get_const_fix(i).size(); ++j) {
                const auto scaled_val = constraint->get_const_fix(i)[j].val_to_fix * scale_factor;
                constraint->set_const_fix_val_to_fix(i, j, scaled_val);
            }
        }
    }
}

void Optimize::apply_scaler_force_constants(const int maxorder,
                                            const double normalization_factor,
                                            const Constraint *constraint,
                                            std::vector<double> &param_inout) const
{
    auto k = 0;
    for (auto i = 0; i < maxorder; ++i) {
        const auto scale_factor = 1.0 / std::pow(normalization_factor, i + 1);

        for (auto j = 0; j < constraint->get_index_bimap(i).size(); ++j) {
            param_inout[k] *= scale_factor;
            ++k;
        }
    }
}

void Optimize::apply_scalers(const int maxorder,
                             Constraint *constraint)
{
    apply_scaler_displacement(u_train,
                              optcontrol.displacement_normalization_factor);
    apply_scaler_constraint(maxorder,
                            optcontrol.displacement_normalization_factor,
                            constraint);

    if (optcontrol.cross_validation == -1) {
        apply_scaler_displacement(u_validation,
                                  optcontrol.displacement_normalization_factor);

    }
}

void Optimize::finalize_scalers(const int maxorder,
                                Constraint *constraint)
{
    apply_scaler_displacement(u_train,
                              optcontrol.displacement_normalization_factor,
                              true);
    apply_scaler_constraint(maxorder,
                            optcontrol.displacement_normalization_factor,
                            constraint,
                            true);
    if (optcontrol.cross_validation == -1) {
        apply_scaler_displacement(u_validation,
                                  optcontrol.displacement_normalization_factor,
                                  true);
    }
}

void Optimize::apply_basis_converter(std::vector<std::vector<double>> &u_multi,
                                     Eigen::Matrix3d cmat) const
{
    // Convert the basis of displacements from Cartesian to fractional
    const auto nrows = u_multi.size();
    const auto ncols = u_multi[0].size();
    size_t i, j;
    Eigen::Vector3d vec_tmp;

    const auto nat = ncols / 3;
    for (i = 0; i < nrows; ++i) {
        for (j = 0; j < nat; ++j) {
            for (int k = 0; k < 3; ++k) {
                vec_tmp(k) = u_multi[i][3 * j + k];
            }
            vec_tmp = cmat * vec_tmp;
            for (int k = 0; k < 3; ++k) {
                u_multi[i][3 * j + k] = vec_tmp(k);
            }
        }
    }
}

void Optimize::apply_basis_converter_amat(const int natmin3,
                                          const int ncols,
                                          double **amat_orig_tmp,
                                          Eigen::Matrix3d cmat) const
{
    const auto natmin = natmin3 / 3;
    Eigen::Vector3d vec_tmp;
    const Eigen::Matrix3d cmat_t = cmat.transpose();

    for (auto icol = 0; icol < ncols; ++ icol) {
        for (auto iat = 0; iat < natmin; ++iat) {
            for (auto i = 0; i < 3; ++i) {
                vec_tmp(i) = amat_orig_tmp[3 * iat + i][icol];
            }
            vec_tmp = cmat_t * vec_tmp;
            for (auto i = 0; i < 3; ++i) {
                amat_orig_tmp[3 * iat + i][icol] = vec_tmp(i);
            }
        }
    }
}


void Optimize::set_training_data(const std::vector<std::vector<double>> &u_train_in,
                                 const std::vector<std::vector<double>> &f_train_in)
{
    u_train.clear();
    f_train.clear();
    u_train = u_train_in;
    f_train = f_train_in;
    u_train.shrink_to_fit();
    f_train.shrink_to_fit();
}

void Optimize::set_validation_data(const std::vector<std::vector<double>> &u_validation_in,
                                   const std::vector<std::vector<double>> &f_validation_in)
{
    u_validation.clear();
    f_validation.clear();
    u_validation = u_validation_in;
    f_validation = f_validation_in;
    u_validation.shrink_to_fit();
    f_validation.shrink_to_fit();
}

std::vector<std::vector<double>> Optimize::get_u_train() const
{
    return u_train;
}

std::vector<std::vector<double>> Optimize::get_f_train() const
{
    return f_train;
}


void Optimize::set_fcs_values(const int maxorder,
                              double *fc_in,
                              std::vector<size_t> *nequiv,
                              const Constraint *constraint)
{
    // fc_in: irreducible set of force constants
    // fc_length: dimension of params (can differ from that of fc_in)

    size_t i;

    size_t N = 0;
    size_t Nirred = 0;
    for (i = 0; i < maxorder; ++i) {
        N += nequiv[i].size();
        Nirred += constraint->get_index_bimap(i).size();
    }

    std::vector<double> param_in(Nirred, 0.0);
    std::vector<double> param_out(N, 0.0);

    for (i = 0; i < Nirred; ++i) {
        param_in[i] = fc_in[i];
    }
    recover_original_forceconstants(maxorder,
                                    param_in,
                                    param_out,
                                    nequiv,
                                    constraint);
    if (params) {
        deallocate(params);
    }
    allocate(params, N);
    for (i = 0; i < N; ++i) {
        params[i] = param_out[i];
    }
}

size_t Optimize::get_number_of_rows_sensing_matrix() const
{
    return u_train.size() * u_train[0].size();
}


int Optimize::fit_without_constraints(const size_t N,
                                      const size_t M,
                                      double *amat,
                                      const double *bvec,
                                      double *param_out,
                                      const int verbosity) const
{
    int i;
    int nrhs = 1, nrank, INFO, M_tmp, N_tmp;
    auto rcond = -1.0;
    auto f_square = 0.0;
    double *WORK, *S, *fsum2;


    const auto LMIN = std::min<int>(M, N);
    auto LMAX = std::max<int>(M, N);

    auto LWORK = 3 * LMIN + std::max<int>(2 * LMIN, LMAX);
    LWORK = 2 * LWORK;

    if (verbosity > 0) {
        std::cout << "  Entering fitting routine: SVD without constraints" << std::endl;
    }


    allocate(WORK, LWORK);
    allocate(S, LMIN);
    allocate(fsum2, LMAX);

    for (i = 0; i < M; ++i) {
        fsum2[i] = bvec[i];
        f_square += std::pow(bvec[i], 2);
    }
    for (i = M; i < LMAX; ++i) fsum2[i] = 0.0;

    if (verbosity > 0) std::cout << "  SVD has started ... ";

    // Fitting with singular value decomposition
    // M_tmp and N_tmp are prepared to cast N and M to (non-const) int.
    M_tmp = M;
    N_tmp = N;
    dgelss_(&M_tmp, &N_tmp, &nrhs, amat, &M_tmp, fsum2, &LMAX,
            S, &rcond, &nrank, WORK, &LWORK, &INFO);

    if (verbosity > 0) {
        std::cout << "finished !" << std::endl << std::endl;
        std::cout << "  RANK of the matrix = " << nrank << std::endl;
    }

    if (nrank < N)
        warn("fit_without_constraints",
             "Matrix is rank-deficient. Force constants could not be determined uniquely :(");

    if (nrank == N && verbosity > 0) {
        auto f_residual = 0.0;
        for (i = N; i < M; ++i) {
            f_residual += std::pow(fsum2[i], 2);
        }
        std::cout << std::endl << "  Residual sum of squares for the solution: "
            << sqrt(f_residual) << std::endl;
        std::cout << "  Fitting error (%) : "
            << sqrt(f_residual / f_square) * 100.0 << std::endl;
    }

    for (i = 0; i < N; ++i) {
        param_out[i] = fsum2[i];
    }

    deallocate(WORK);
    deallocate(S);
    deallocate(fsum2);

    return INFO;
}

int Optimize::fit_with_constraints(const size_t N,
                                   const size_t M,
                                   const size_t P,
                                   double *amat,
                                   const double *bvec,
                                   double *param_out,
                                   const double * const *cmat,
                                   double *dvec,
                                   const int verbosity) const
{
    size_t i, j;
    int N_tmp, M_tmp, P_tmp;
    double *fsum2;
    double *mat_tmp;

    if (verbosity > 0) {
        std::cout << "  Entering fitting routine: QRD with constraints" << std::endl;
    }

    allocate(fsum2, M);
    allocate(mat_tmp, (M + P) * N);

    size_t k = 0;
    size_t l = 0;

    // Concatenate two matrices as 1D array
    for (j = 0; j < N; ++j) {
        for (i = 0; i < M; ++i) {
            mat_tmp[k++] = amat[l++];
        }
        for (i = 0; i < P; ++i) {
            mat_tmp[k++] = cmat[i][j];
        }
    }

    const auto nrank = rankQRD((M + P), N, mat_tmp, eps12);
    deallocate(mat_tmp);

    if (nrank != N) {
        std::cout << std::endl;
        std::cout << " **************************************************************************" << std::endl;
        std::cout << "  WARNING : rank deficient.                                                " << std::endl;
        std::cout << "  rank ( (A) ) ! = N            A: Fitting matrix     B: Constraint matrix " << std::endl;
        std::cout << "       ( (B) )                  N: The number of parameters                " << std::endl;
        std::cout << "  rank = " << nrank << " N = " << N << std::endl << std::endl;
        std::cout << "  This can cause a difficulty in solving the fitting problem properly      " << std::endl;
        std::cout << "  with DGGLSE, especially when the difference is large. Please check if    " << std::endl;
        std::cout << "  you obtain reliable force constants in the .fcs file.                    " << std::endl << std::
            endl;
        std::cout << "  You may need to reduce the cutoff radii and/or increase NDATA            " << std::endl;
        std::cout << "  by giving linearly-independent displacement patterns.                    " << std::endl;
        std::cout << " **************************************************************************" << std::endl;
        std::cout << std::endl;
    }

    auto f_square = 0.0;
    for (i = 0; i < M; ++i) {
        fsum2[i] = bvec[i];
        f_square += std::pow(bvec[i], 2);
    }
    if (verbosity > 0) std::cout << "  QR-Decomposition has started ...";

    double *cmat_mod;
    allocate(cmat_mod, P * N);

    k = 0;
    for (j = 0; j < N; ++j) {
        for (i = 0; i < P; ++i) {
            cmat_mod[k++] = cmat[i][j];
        }
    }

    // Fitting

    auto LWORK = static_cast<int>(P) + std::min<int>(M, N) + 10 * std::max<int>(M, N);
    int INFO;
    double *WORK, *x;
    allocate(WORK, LWORK);
    allocate(x, N);

    // M_tmp, N_tmp, P_tmp are prepared to cast N, M, P to (non-const)
    // int.
    M_tmp = M;
    N_tmp = N;
    P_tmp = P;
    dgglse_(&M_tmp, &N_tmp, &P_tmp, amat, &M_tmp, cmat_mod, &P_tmp,
            fsum2, dvec, x, WORK, &LWORK, &INFO);

    if (verbosity > 0) std::cout << " finished. " << std::endl;

    auto f_residual = 0.0;
    for (i = N - P; i < M; ++i) {
        f_residual += std::pow(fsum2[i], 2);
    }

    if (verbosity > 0) {
        std::cout << std::endl << "  Residual sum of squares for the solution: "
            << sqrt(f_residual) << std::endl;
        std::cout << "  Fitting error (%) : "
            << std::sqrt(f_residual / f_square) * 100.0 << std::endl;
    }

    // copy fcs to bvec
    for (i = 0; i < N; ++i) {
        param_out[i] = x[i];
    }

    deallocate(cmat_mod);
    deallocate(WORK);
    deallocate(x);
    deallocate(fsum2);

    return INFO;
}

int Optimize::fit_algebraic_constraints(const size_t N,
                                        const size_t M,
                                        double *amat,
                                        const double *bvec,
                                        std::vector<double> &param_out,
                                        const double fnorm,
                                        const int maxorder,
                                        const Fcs *fcs,
                                        const Constraint *constraint,
                                        const int verbosity) const
{
    int i;
    int nrhs = 1, nrank, INFO, M_tmp, N_tmp;
    auto rcond = -1.0;
    double *WORK, *S, *fsum2;

    if (verbosity > 0) {
        std::cout << "  Entering fitting routine: SVD with constraints considered algebraically." << std::endl;
    }

    auto LMIN = std::min<int>(M, N);
    auto LMAX = std::max<int>(M, N);

    auto LWORK = 3 * LMIN + std::max<int>(2 * LMIN, LMAX);
    LWORK = 2 * LWORK;

    allocate(WORK, LWORK);
    allocate(S, LMIN);
    allocate(fsum2, LMAX);

    for (i = 0; i < M; ++i) {
        fsum2[i] = bvec[i];
    }
    for (i = M; i < LMAX; ++i) fsum2[i] = 0.0;

    if (verbosity > 0) std::cout << "  SVD has started ... ";

    // Fitting with singular value decomposition
    // M_tmp and N_tmp are prepared to cast N and M to (non-const) int.
    M_tmp = M;
    N_tmp = N;
    dgelss_(&M_tmp, &N_tmp, &nrhs, amat, &M_tmp, fsum2, &LMAX,
            S, &rcond, &nrank, WORK, &LWORK, &INFO);

    deallocate(WORK);
    deallocate(S);

    if (verbosity > 0) {
        std::cout << "finished !" << std::endl << std::endl;
        std::cout << "  RANK of the matrix = " << nrank << std::endl;
    }

    if (nrank < N) {
        std::cout << " **************************************************************************\n";
        std::cout << "  WARNING : Rank deficient                                                 \n\n";
        std::cout << "  Force constants could not be determined uniquely because                 \n";
        std::cout << "  the sensing matrix is not full rank.                                     \n";
        std::cout << "  You may need to reduce the cutoff radii and/or increase the number of    \n";
        std::cout << "  training datasets.                                                       \n";
        std::cout << " **************************************************************************\n";
    }

    if (nrank == N && verbosity > 0) {
        auto f_residual = 0.0;
        for (i = N; i < M; ++i) {
            f_residual += std::pow(fsum2[i], 2);
        }
        std::cout << std::endl;
        std::cout << "  Residual sum of squares for the solution: "
            << sqrt(f_residual) << std::endl;
        std::cout << "  Fitting error (%) : "
            << sqrt(f_residual / (fnorm * fnorm)) * 100.0 << std::endl;
    }

    if (INFO == 0) {
        std::vector<double> param_irred(N, 0.0);
        for (i = 0; i < LMIN; ++i) param_irred[i] = fsum2[i];
        deallocate(fsum2);

        // Recover reducible set of force constants

        recover_original_forceconstants(maxorder,
                                        param_irred,
                                        param_out,
                                        fcs->get_nequiv(),
                                        constraint);
    }

    return INFO;
}


void Optimize::get_matrix_elements(const int maxorder,
                                   std::vector<double> &amat,
                                   std::vector<double> &bvec,
                                   const std::vector<std::vector<double>> &u_in,
                                   const std::vector<std::vector<double>> &f_in,
                                   const Symmetry *symmetry,
                                   const Fcs *fcs) const
{
    size_t i, j;
    long irow;
    const auto natmin = symmetry->get_nat_prim();
    const auto natmin3 = 3 * natmin;
    std::vector<std::vector<double>> u_multi, f_multi;

    if (u_in.size() != f_in.size()) {
        exit("get_matrix_elements",
             "The lengths of displacement array and force array are diferent.");
    }

    const auto ndata_fit = u_in.size();
    const auto ncycle = ndata_fit * symmetry->get_ntran();
    const auto nrows = ndata_fit * u_in[0].size();
    size_t ncols = 0;
    for (i = 0; i < maxorder; ++i) {
        ncols += fcs->get_nequiv()[i].size();
    }

    if (amat.size() != nrows * ncols) {
        amat.resize(nrows * ncols, 0.0);
    }
    if (bvec.size() != nrows) {
        bvec.resize(nrows, 0.0);
    }

    data_multiplier(u_in, u_multi, symmetry);
    data_multiplier(f_in, f_multi, symmetry);

    if (fcs->get_forceconstant_basis() == "Lattice") {
        apply_basis_converter(u_multi,
                              fcs->get_basis_conversion_matrix());
    }


#ifdef _OPENMP
#pragma omp parallel private(irow, i, j)
#endif
    {
        int *ind;
        int mm, order, iat, k;
        size_t im, iparam;
        size_t idata;
        double amat_tmp;
        double **amat_orig_tmp;

        allocate(ind, maxorder + 1);
        allocate(amat_orig_tmp, natmin3, ncols);

#ifdef _OPENMP
#pragma omp for schedule(guided)
#endif
        for (irow = 0; irow < ncycle; ++irow) {

            // generate r.h.s vector B
            for (i = 0; i < natmin; ++i) {
                iat = symmetry->get_map_p2s()[i][0];
                for (j = 0; j < 3; ++j) {
                    im = 3 * i + j + natmin3 * irow;
                    bvec[im] = f_multi[irow][3 * iat + j];
                }
            }

            for (i = 0; i < natmin3; ++i) {
                for (j = 0; j < ncols; ++j) {
                    amat_orig_tmp[i][j] = 0.0;
                }
            }

            // generate l.h.s. matrix A

            idata = natmin3 * irow;
            iparam = 0;

            for (order = 0; order < maxorder; ++order) {

                mm = 0;

                for (const auto &iter : fcs->get_nequiv()[order]) {
                    for (i = 0; i < iter; ++i) {
                        ind[0] = fcs->get_fc_table()[order][mm].elems[0];
                        k = inprim_index(ind[0], symmetry);
                        amat_tmp = 1.0;
                        for (j = 1; j < order + 2; ++j) {
                            ind[j] = fcs->get_fc_table()[order][mm].elems[j];
                            amat_tmp *= u_multi[irow][fcs->get_fc_table()[order][mm].elems[j]];
                        }
                        amat_orig_tmp[k][iparam] -= gamma(order + 2, ind) * fcs->get_fc_table()[order][mm].sign *
                            amat_tmp;
                        ++mm;
                    }
                    ++iparam;
                }
            }

            // When the force constants are defined in the fractional coordinate,
            // we need to multiply the basis_conversion_matrix to obtain atomic forces
            // in the Cartesian coordinate.
            if (fcs->get_forceconstant_basis() == "Lattice") {
                apply_basis_converter_amat(natmin3,
                                           ncols,
                                           amat_orig_tmp,
                                           fcs->get_basis_conversion_matrix());
            }

            for (i = 0; i < natmin3; ++i) {
                for (j = 0; j < ncols; ++j) {
                    // Transpose here for later use of lapack without transpose
                    amat[natmin3 * ncycle * j + i + idata] = amat_orig_tmp[i][j];
                }
            }
        }

        deallocate(ind);
        deallocate(amat_orig_tmp);
    }

    u_multi.clear();
    f_multi.clear();
}


void Optimize::get_matrix_elements_algebraic_constraint(const int maxorder,
                                                        std::vector<double> &amat,
                                                        std::vector<double> &bvec,
                                                        const std::vector<std::vector<double>> &u_in,
                                                        const std::vector<std::vector<double>> &f_in,
                                                        double &fnorm,
                                                        const Symmetry *symmetry,
                                                        const Fcs *fcs,
                                                        const Constraint *constraint) const
{
    size_t i, j;
    long irow;

    if (u_in.size() != f_in.size()) {
        exit("get_matrix_elements",
             "The lengths of displacement array and force array are diferent.");
    }

    const auto ndata_fit = u_in.size();
    const auto natmin = symmetry->get_nat_prim();
    const auto natmin3 = 3 * natmin;
    const auto nrows = u_in.size() * u_in[0].size();
    size_t ncols = 0;
    size_t ncols_new = 0;

    for (i = 0; i < maxorder; ++i) {
        ncols += fcs->get_nequiv()[i].size();
        ncols_new += constraint->get_index_bimap(i).size();
    }

    const auto ncycle = ndata_fit * symmetry->get_ntran();

    if (amat.size() != nrows * ncols_new) {
        amat.resize(nrows * ncols_new, 0.0);
    }
    if (bvec.size() != nrows) {
        bvec.resize(nrows, 0.0);
    }

    std::vector<double> bvec_orig(nrows, 0.0);
    std::vector<std::vector<double>> u_multi, f_multi;

    data_multiplier(u_in, u_multi, symmetry);
    data_multiplier(f_in, f_multi, symmetry);

    if (fcs->get_forceconstant_basis() == "Lattice") {
        apply_basis_converter(u_multi,
                              fcs->get_basis_conversion_matrix());
    }

#ifdef _OPENMP
#pragma omp parallel private(irow, i, j)
#endif
    {
        int *ind;
        int mm, order, iat, k;
        size_t im;
        size_t idata;
        size_t ishift, iparam;
        size_t iold, inew;
        double amat_tmp;
        double **amat_orig_tmp;
        double **amat_mod_tmp;

        allocate(ind, maxorder + 1);
        allocate(amat_orig_tmp, natmin3, ncols);
        allocate(amat_mod_tmp, natmin3, ncols_new);

#ifdef _OPENMP
#pragma omp for schedule(guided)
#endif
        for (irow = 0; irow < ncycle; ++irow) {

            // generate r.h.s vector B
            for (i = 0; i < natmin; ++i) {
                iat = symmetry->get_map_p2s()[i][0];
                for (j = 0; j < 3; ++j) {
                    im = 3 * i + j + natmin3 * irow;
                    bvec[im] = f_multi[irow][3 * iat + j];
                    bvec_orig[im] = f_multi[irow][3 * iat + j];
                }
            }

            for (i = 0; i < natmin3; ++i) {
                for (j = 0; j < ncols; ++j) {
                    amat_orig_tmp[i][j] = 0.0;
                }
                for (j = 0; j < ncols_new; ++j) {
                    amat_mod_tmp[i][j] = 0.0;
                }
            }

            // generate l.h.s. matrix A

            idata = natmin3 * irow;
            iparam = 0;

            for (order = 0; order < maxorder; ++order) {

                mm = 0;

                for (const auto &iter : fcs->get_nequiv()[order]) {
                    for (i = 0; i < iter; ++i) {
                        ind[0] = fcs->get_fc_table()[order][mm].elems[0];
                        k = inprim_index(ind[0], symmetry);

                        amat_tmp = 1.0;
                        for (j = 1; j < order + 2; ++j) {
                            ind[j] = fcs->get_fc_table()[order][mm].elems[j];
                            amat_tmp *= u_multi[irow][fcs->get_fc_table()[order][mm].elems[j]];
                        }
                        amat_orig_tmp[k][iparam] -= gamma(order + 2, ind)
                            * fcs->get_fc_table()[order][mm].sign * amat_tmp;
                        ++mm;
                    }
                    ++iparam;
                }
            }

            // When the force constants are defined in the fractional coordinate,
            // we need to multiply the basis_conversion_matrix to obtain atomic forces
            // in the Cartesian coordinate.
            if (fcs->get_forceconstant_basis() == "Lattice") {
                apply_basis_converter_amat(natmin3,
                                           ncols,
                                           amat_orig_tmp,
                                           fcs->get_basis_conversion_matrix());
            }

            // Convert the full matrix and vector into a smaller irreducible form
            // by using constraint information.

            ishift = 0;
            iparam = 0;

            for (order = 0; order < maxorder; ++order) {

                for (i = 0; i < constraint->get_const_fix(order).size(); ++i) {

                    for (j = 0; j < natmin3; ++j) {
                        bvec[j + idata] -= constraint->get_const_fix(order)[i].val_to_fix
                            * amat_orig_tmp[j][ishift + constraint->get_const_fix(order)[i].p_index_target];
                    }
                }

                for (const auto &it : constraint->get_index_bimap(order)) {
                    inew = it.left + iparam;
                    iold = it.right + ishift;

                    for (j = 0; j < natmin3; ++j) {
                        amat_mod_tmp[j][inew] = amat_orig_tmp[j][iold];
                    }
                }

                for (i = 0; i < constraint->get_const_relate(order).size(); ++i) {

                    iold = constraint->get_const_relate(order)[i].p_index_target + ishift;

                    for (j = 0; j < constraint->get_const_relate(order)[i].alpha.size(); ++j) {

                        inew = constraint->get_index_bimap(order).right.at(
                                constraint->get_const_relate(order)[i].p_index_orig[j]) +
                            iparam;

                        for (k = 0; k < natmin3; ++k) {
                            amat_mod_tmp[k][inew] -= amat_orig_tmp[k][iold]
                                * constraint->get_const_relate(order)[i].alpha[j];
                        }
                    }
                }

                ishift += fcs->get_nequiv()[order].size();
                iparam += constraint->get_index_bimap(order).size();
            }

            for (i = 0; i < natmin3; ++i) {
                for (j = 0; j < ncols_new; ++j) {
                    // Transpose here for later use of lapack without transpose
                    amat[natmin3 * ncycle * j + i + idata] = amat_mod_tmp[i][j];
                }
            }
        }

        deallocate(ind);
        deallocate(amat_orig_tmp);
        deallocate(amat_mod_tmp);
    }

    fnorm = 0.0;
    for (i = 0; i < bvec_orig.size(); ++i) {
        fnorm += bvec_orig[i] * bvec_orig[i];
    }
    fnorm = std::sqrt(fnorm);

    u_multi.clear();
    f_multi.clear();
}

void Optimize::get_matrix_elements_in_sparse_form(const int maxorder,
                                                  SpMat &sp_amat,
                                                  Eigen::VectorXd &sp_bvec,
                                                  const std::vector<std::vector<double>> &u_in,
                                                  const std::vector<std::vector<double>> &f_in,
                                                  double &fnorm,
                                                  const Symmetry *symmetry,
                                                  const Fcs *fcs,
                                                  const Constraint *constraint) const
{
    size_t i, j;
    long irow;
    typedef Eigen::Triplet<double, size_t> T;
    std::vector<T> nonzero_entries;

    if (u_in.size() != f_in.size()) {
        exit("get_matrix_elements",
             "The lengths of displacement array and force array are diferent.");
    }

    const auto ndata_fit = u_in.size();
    const auto natmin = symmetry->get_nat_prim();
    const auto natmin3 = 3 * natmin;
    const auto nrows = u_in.size() * u_in[0].size();

    size_t ncols = 0;
    size_t ncols_new = 0;

    for (i = 0; i < maxorder; ++i) {
        ncols += fcs->get_nequiv()[i].size();
        ncols_new += constraint->get_index_bimap(i).size();
    }

    const auto ncycle = ndata_fit * symmetry->get_ntran();

    std::vector<double> bvec_orig(nrows, 0.0);
    std::vector<std::vector<double>> u_multi, f_multi;

    data_multiplier(u_in, u_multi, symmetry);
    data_multiplier(f_in, f_multi, symmetry);

    if (fcs->get_forceconstant_basis() == "Lattice") {
        apply_basis_converter(u_multi,
                              fcs->get_basis_conversion_matrix());
    }

#ifdef _OPENMP
#pragma omp parallel private(irow, i, j)
#endif
    {
        int *ind;
        int mm, order, iat, k;
        size_t im, iparam;
        size_t idata;
        size_t ishift;
        size_t iold, inew;
        double amat_tmp;
        double **amat_orig_tmp;
        double **amat_mod_tmp;

        std::vector<T> nonzero_omp;

        allocate(ind, maxorder + 1);
        allocate(amat_orig_tmp, natmin3, ncols);
        allocate(amat_mod_tmp, natmin3, ncols_new);

#ifdef _OPENMP
#pragma omp for schedule(guided)
#endif
        for (irow = 0; irow < ncycle; ++irow) {

            // generate r.h.s vector B
            for (i = 0; i < natmin; ++i) {
                iat = symmetry->get_map_p2s()[i][0];
                for (j = 0; j < 3; ++j) {
                    im = 3 * i + j + natmin3 * irow;
                    sp_bvec(im) = f_multi[irow][3 * iat + j];
                    bvec_orig[im] = f_multi[irow][3 * iat + j];
                }
            }

            for (i = 0; i < natmin3; ++i) {
                for (j = 0; j < ncols; ++j) {
                    amat_orig_tmp[i][j] = 0.0;
                }
                for (j = 0; j < ncols_new; ++j) {
                    amat_mod_tmp[i][j] = 0.0;
                }
            }

            // generate l.h.s. matrix A

            idata = natmin3 * irow;
            iparam = 0;

            for (order = 0; order < maxorder; ++order) {

                mm = 0;

                for (const auto &iter : fcs->get_nequiv()[order]) {
                    for (i = 0; i < iter; ++i) {
                        ind[0] = fcs->get_fc_table()[order][mm].elems[0];
                        k = inprim_index(ind[0], symmetry);

                        amat_tmp = 1.0;
                        for (j = 1; j < order + 2; ++j) {
                            ind[j] = fcs->get_fc_table()[order][mm].elems[j];
                            amat_tmp *= u_multi[irow][fcs->get_fc_table()[order][mm].elems[j]];
                        }
                        amat_orig_tmp[k][iparam] -= gamma(order + 2, ind)
                            * fcs->get_fc_table()[order][mm].sign * amat_tmp;
                        ++mm;
                    }
                    ++iparam;
                }
            }

            // When the force constants are defined in the fractional coordinate,
            // we need to multiply the basis_conversion_matrix to obtain atomic forces
            // in the Cartesian coordinate.
            if (fcs->get_forceconstant_basis() == "Lattice") {
                apply_basis_converter_amat(natmin3,
                                           ncols,
                                           amat_orig_tmp,
                                           fcs->get_basis_conversion_matrix());
            }

            // Convert the full matrix and vector into a smaller irreducible form
            // by using constraint information.

            ishift = 0;
            iparam = 0;

            for (order = 0; order < maxorder; ++order) {

                for (i = 0; i < constraint->get_const_fix(order).size(); ++i) {

                    for (j = 0; j < natmin3; ++j) {
                        sp_bvec(j + idata) -= constraint->get_const_fix(order)[i].val_to_fix
                            * amat_orig_tmp[j][ishift + constraint->get_const_fix(order)[i].p_index_target];
                    }
                }

                for (const auto &it : constraint->get_index_bimap(order)) {
                    inew = it.left + iparam;
                    iold = it.right + ishift;

                    for (j = 0; j < natmin3; ++j) {
                        amat_mod_tmp[j][inew] = amat_orig_tmp[j][iold];
                    }
                }

                for (i = 0; i < constraint->get_const_relate(order).size(); ++i) {

                    iold = constraint->get_const_relate(order)[i].p_index_target + ishift;

                    for (j = 0; j < constraint->get_const_relate(order)[i].alpha.size(); ++j) {

                        inew = constraint->get_index_bimap(order).right.at(
                                constraint->get_const_relate(order)[i].p_index_orig[j]) +
                            iparam;

                        for (k = 0; k < natmin3; ++k) {
                            amat_mod_tmp[k][inew] -= amat_orig_tmp[k][iold]
                                * constraint->get_const_relate(order)[i].alpha[j];
                        }
                    }
                }

                ishift += fcs->get_nequiv()[order].size();
                iparam += constraint->get_index_bimap(order).size();
            }

            for (i = 0; i < natmin3; ++i) {
                for (j = 0; j < ncols_new; ++j) {
                    if (std::abs(amat_mod_tmp[i][j]) > eps) {
                        nonzero_omp.emplace_back(T(idata + i, j, amat_mod_tmp[i][j]));
                    }
                }
            }
        }

        deallocate(ind);
        deallocate(amat_orig_tmp);
        deallocate(amat_mod_tmp);

#pragma omp critical
        {
            for (const auto &it : nonzero_omp) {
                nonzero_entries.emplace_back(it);
            }
        }
    }

    fnorm = 0.0;
    for (i = 0; i < bvec_orig.size(); ++i) {
        fnorm += bvec_orig[i] * bvec_orig[i];
    }
    fnorm = std::sqrt(fnorm);
    sp_amat.setFromTriplets(nonzero_entries.begin(), nonzero_entries.end());
    sp_amat.makeCompressed();
}


void Optimize::recover_original_forceconstants(const int maxorder,
                                               const std::vector<double> &param_in,
                                               std::vector<double> &param_out,
                                               const std::vector<size_t> *nequiv,
                                               const Constraint *constraint) const
{
    // Expand the given force constants into the larger sets
    // by using the constraint matrix.

    size_t i, j, k;
    size_t ishift = 0;
    size_t iparam = 0;
    double tmp;
    size_t inew, iold;

    size_t nparams = 0;

    for (i = 0; i < maxorder; ++i) nparams += nequiv[i].size();

    param_out.resize(nparams, 0.0);

    for (i = 0; i < maxorder; ++i) {
        for (j = 0; j < constraint->get_const_fix(i).size(); ++j) {
            param_out[constraint->get_const_fix(i)[j].p_index_target + ishift]
                = constraint->get_const_fix(i)[j].val_to_fix;
        }

        for (const auto &it : constraint->get_index_bimap(i)) {
            inew = it.left + iparam;
            iold = it.right + ishift;

            param_out[iold] = param_in[inew];
        }

        for (j = 0; j < constraint->get_const_relate(i).size(); ++j) {
            tmp = 0.0;

            for (k = 0; k < constraint->get_const_relate(i)[j].alpha.size(); ++k) {
                tmp += constraint->get_const_relate(i)[j].alpha[k]
                    * param_out[constraint->get_const_relate(i)[j].p_index_orig[k] + ishift];
            }
            param_out[constraint->get_const_relate(i)[j].p_index_target + ishift] = -tmp;
        }

        ishift += nequiv[i].size();
        iparam += constraint->get_index_bimap(i).size();
    }
}


void Optimize::data_multiplier(const std::vector<std::vector<double>> &data_in,
                               std::vector<std::vector<double>> &data_out,
                               const Symmetry *symmetry) const
{
    const auto nat = symmetry->get_nat_prim() * symmetry->get_ntran();
    const auto ndata_used = data_in.size();
    const auto ntran = symmetry->get_ntran();

    data_out.resize(ntran * ndata_used, std::vector<double>(3 * nat));

    auto idata = 0;
    for (auto i = 0; i < ndata_used; ++i) {

        for (auto itran = 0; itran < symmetry->get_ntran(); ++itran) {
            for (auto j = 0; j < nat; ++j) {
                const auto n_mapped = symmetry->get_map_sym()[j][symmetry->get_symnum_tran()[itran]];
                for (auto k = 0; k < 3; ++k) {
                    data_out[idata][3 * n_mapped + k] = data_in[i][3 * j + k];
                }
            }
            ++idata;
        }
    }
}

int Optimize::inprim_index(const int n,
                           const Symmetry *symmetry) const
{
    auto in = -1;
    const auto atmn = n / 3;
    const auto crdn = n % 3;

    for (size_t i = 0; i < symmetry->get_nat_prim(); ++i) {
        if (symmetry->get_map_p2s()[i][0] == atmn) {
            in = 3 * i + crdn;
            break;
        }
    }
    return in;
}

double Optimize::gamma(const int n,
                       const int *arr) const
{
    int *arr_tmp, *nsame;
    int i;

    allocate(arr_tmp, n);
    allocate(nsame, n);

    for (i = 0; i < n; ++i) {
        arr_tmp[i] = arr[i];
        nsame[i] = 0;
    }

    const auto ind_front = arr[0];
    auto nsame_to_front = 1;

    insort(n, arr_tmp);

    auto nuniq = 1;
    auto iuniq = 0;

    nsame[0] = 1;

    for (i = 1; i < n; ++i) {
        if (arr_tmp[i] == arr_tmp[i - 1]) {
            ++nsame[iuniq];
        } else {
            ++nsame[++iuniq];
            ++nuniq;
        }

        if (arr[i] == ind_front) ++nsame_to_front;
    }

    auto denom = 1;

    for (i = 0; i < nuniq; ++i) {
        denom *= factorial(nsame[i]);
    }

    deallocate(arr_tmp);
    deallocate(nsame);

    return static_cast<double>(nsame_to_front) / static_cast<double>(denom);
}


double* Optimize::get_params() const
{
    return params;
}

int Optimize::factorial(const int n) const
{
    if (n == 1 || n == 0) {
        return 1;
    }
    return n * factorial(n - 1);
}


int Optimize::rankQRD(const size_t m,
                      const size_t n,
                      double *mat,
                      const double tolerance) const
{
    // Return the rank of matrix mat revealed by the column pivoting QR decomposition
    // The matrix mat is destroyed.

    auto m_ = static_cast<int>(m);
    auto n_ = static_cast<int>(n);

    auto LDA = m_;

    auto LWORK = 10 * n_;
    int INFO;
    int *JPVT;
    double *WORK, *TAU;

    const auto nmin = std::min<int>(m_, n_);

    allocate(JPVT, n_);
    allocate(WORK, LWORK);
    allocate(TAU, nmin);

    for (auto i = 0; i < n_; ++i) JPVT[i] = 0;

    dgeqp3_(&m_, &n_, mat, &LDA, JPVT, TAU, WORK, &LWORK, &INFO);

    deallocate(JPVT);
    deallocate(WORK);
    deallocate(TAU);

    if (std::abs(mat[0]) < eps) return 0;

    double **mat_tmp;
    allocate(mat_tmp, m_, n_);

    unsigned long k = 0;

    for (auto j = 0; j < n_; ++j) {
        for (auto i = 0; i < m_; ++i) {
            mat_tmp[i][j] = mat[k++];
        }
    }

    auto nrank = 0;
    for (auto i = 0; i < nmin; ++i) {
        if (std::abs(mat_tmp[i][i]) > tolerance * std::abs(mat[0])) ++nrank;
    }

    deallocate(mat_tmp);

    return nrank;
}


int Optimize::run_eigen_sparse_solver(const SpMat &sp_mat,
                                      const Eigen::VectorXd &sp_bvec,
                                      std::vector<double> &param_out,
                                      const double fnorm,
                                      const int maxorder,
                                      const Fcs *fcs,
                                      const Constraint *constraint,
                                      const std::string solver_type,
                                      const int verbosity) const
{
    const auto solver_type_lower = boost::algorithm::to_lower_copy(solver_type);
    Eigen::VectorXd x;

    if (verbosity > 0) {
        std::cout << "  Solve least-squares problem by Eigen " + solver_type + ".\n";
    }

    if (solver_type_lower == "simplicialldlt") {
        SpMat AtA = sp_mat.transpose() * sp_mat;
        Eigen::VectorXd AtB = sp_mat.transpose() * sp_bvec;

        Eigen::SimplicialLDLT<SpMat> ldlt(AtA);
        x = ldlt.solve(AtB);

        if (ldlt.info() != Eigen::Success) {
            std::cerr << "  Fitting by " + solver_type + " failed." << std::endl;
            std::cerr << ldlt.info() << std::endl;
            return 1;
        }

    } else if (solver_type_lower == "sparseqr") {

        Eigen::SparseQR<SpMat, Eigen::COLAMDOrdering<int>> qr(sp_mat);
        x = qr.solve(sp_bvec);

        if (qr.info() != Eigen::Success) {
            std::cerr << "  Fitting by " + solver_type + " failed." << std::endl;
            std::cerr << qr.info() << std::endl;
            return 1;
        }

    } else if (solver_type_lower == "conjugategradient") {
        SpMat AtA = sp_mat.transpose() * sp_mat;
        Eigen::VectorXd AtB = sp_mat.transpose() * sp_bvec;

        Eigen::ConjugateGradient<SpMat> cg(AtA);
        cg.setTolerance(optcontrol.tolerance_iteration);
        cg.setMaxIterations(optcontrol.maxnum_iteration);
        x.setZero();
        x = cg.solve(AtB);

        if (cg.info() != Eigen::Success) {
            std::cerr << "  Fitting by " + solver_type + " failed." << std::endl;
            std::cerr << cg.info() << std::endl;
            return 1;
        }

    } else if (solver_type_lower == "leastsquaresconjugategradient") {

#if EIGEN_VERSION_AT_LEAST(3, 3, 0)
        Eigen::LeastSquaresConjugateGradient<SpMat> lscg(sp_mat);
        lscg.setTolerance(optcontrol.tolerance_iteration);
        lscg.setMaxIterations(optcontrol.maxnum_iteration);
        x.setZero();
        x = lscg.solve(sp_bvec);

        if (lscg.info() != Eigen::Success) {
            std::cerr << "  Fitting by " + solver_type + " failed." << std::endl;
            std::cerr << lscg.info() << std::endl;
            return 1;
        }

#else
        std::cerr << "The linked Eigen version is too old\n";
        std::cerr << solver_type + " is available as of 3.3.0\n";
        return 1;
#endif

    } else if (solver_type_lower == "bicgstab") {
        SpMat AtA = sp_mat.transpose() * sp_mat;
        Eigen::VectorXd AtB = sp_mat.transpose() * sp_bvec;

        Eigen::BiCGSTAB<SpMat> bicg(AtA);
        bicg.setTolerance(optcontrol.tolerance_iteration);
        bicg.setMaxIterations(optcontrol.maxnum_iteration);
        x.setZero();
        x = bicg.solve(AtB);

        if (bicg.info() != Eigen::Success) {
            std::cerr << "  Fitting by " + solver_type + " failed." << std::endl;
            std::cerr << bicg.info() << std::endl;
            return 1;
        }
    }

    auto res = sp_bvec - sp_mat * x;
    const auto res2norm = res.squaredNorm();
    const auto nparams = x.size();
    std::vector<double> param_irred(nparams);

    for (auto i = 0; i < nparams; ++i) {
        param_irred[i] = x(i);
    }

    // Recover reducible set of force constants

    recover_original_forceconstants(maxorder,
                                    param_irred,
                                    param_out,
                                    fcs->get_nequiv(),
                                    constraint);

    if (verbosity > 0) {
        std::cout << "  Residual sum of squares for the solution: "
            << sqrt(res2norm) << std::endl;
        std::cout << "  Fitting error (%) : "
            << sqrt(res2norm / (fnorm * fnorm)) * 100.0 << std::endl;
    }

    return 0;
}


void Optimize::set_optimizer_control(const OptimizerControl &optcontrol_in)
{
    // Check the validity of the options before copying it.

    if (optcontrol_in.cross_validation < -1) {
        exit("set_optimizer_control", "cross_validation must be -1, 0, or larger");
    }
    if (optcontrol_in.linear_model == 2) {
        if (optcontrol_in.l1_ratio <= eps || optcontrol_in.l1_ratio > 1.0) {
            exit("set_optimizer_control", "L1_RATIO must be 0 < L1_RATIO <= 1.");
        }

        if (optcontrol_in.cross_validation >= 1 || optcontrol_in.cross_validation == -1) {
            if (optcontrol_in.l1_alpha_max > 0) {
                if (optcontrol_in.l1_alpha_min >= optcontrol_in.l1_alpha_max) {
                    exit("set_optimizer_control", "L1_ALPHA_MIN must be smaller than L1_ALPHA_MAX.");
                }
            }
        }
    }

    optcontrol = optcontrol_in;
}

OptimizerControl Optimize::get_optimizer_control() const
{
    return optcontrol;
}

double Optimize::get_cv_l1_alpha() const
{
    return cv_l1_alpha;
}

void Optimize::coordinate_descent(const int M,
                                  const int N,
                                  const double alpha,
                                  const int warm_start,
                                  Eigen::VectorXd &x,
                                  const Eigen::MatrixXd &A,
                                  const Eigen::VectorXd &b,
                                  const Eigen::VectorXd &grad0,
                                  bool *has_prod,
                                  Eigen::MatrixXd &Prod,
                                  Eigen::VectorXd &grad,
                                  const double fnorm,
                                  const Eigen::VectorXd &scale_beta,
                                  const int verbosity) const
{
    int i, j;
    double diff{0.0};
    Eigen::VectorXd beta(N), delta(N);
    Eigen::VectorXd res(N);
    bool do_print_log;

    if (warm_start) {
        for (i = 0; i < N; ++i) beta(i) = x(i);
    } else {
        for (i = 0; i < N; ++i) beta(i) = 0.0;
        grad = grad0;
    }

    if (verbosity > 1) {
        std::cout << "-----------------------------------------------------------------" << std::endl;
        std::cout << "  L1_ALPHA = " << std::setw(15) << alpha << std::endl;
    }

    const auto Minv = 1.0 / static_cast<double>(M);
    const auto alphlambda = alpha * optcontrol.l1_ratio;

    auto iloop = 0;

    if (optcontrol.standardize) {
        while (iloop < optcontrol.maxnum_iteration) {
            do_print_log = !((iloop + 1) % optcontrol.output_frequency) && (verbosity > 1);

            if (do_print_log) {
                std::cout << "   Coordinate Descent : " << std::setw(5) << iloop + 1 << std::endl;
            }
            delta = beta;
            for (i = 0; i < N; ++i) {
                beta(i) = shrink(Minv * grad(i) + beta(i), alphlambda);
                delta(i) -= beta(i);
                if (std::abs(delta(i)) > 0.0) {
                    if (!has_prod[i]) {
#pragma omp parallel for
                        for (j = 0; j < N; ++j) {
                            Prod(j, i) = A.col(j).dot(A.col(i));
                        }
                        has_prod[i] = true;
                    }
                    grad = grad + Prod.col(i) * delta(i);
                }
            }
            ++iloop;
            diff = 0.0;
#pragma omp parallel for reduction(+:diff)
            for (i = 0; i < N; ++i) {
                diff += delta(i) * delta(i);
            }
            diff = std::sqrt(diff / static_cast<double>(N));
            //diff = std::sqrt(delta.dot(delta) / static_cast<double>(N));

            if (diff < optcontrol.tolerance_iteration) break;

            if (do_print_log) {
                std::cout << "    1: ||u_{k}-u_{k-1}||_2     = " << std::setw(15) << diff
                    << std::setw(15) << diff * std::sqrt(static_cast<double>(N) / beta.dot(beta)) << std::endl;
                auto tmp = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+:tmp)
#endif
                for (i = 0; i < N; ++i) {
                    tmp += std::abs(beta(i));
                }
                std::cout << "    2: ||u_{k}||_1             = " << std::setw(15) << tmp << std::endl;
                res = A * beta - b;
                tmp = res.dot(res);
                std::cout << "    3: ||Au_{k}-f||_2          = " << std::setw(15) << std::sqrt(tmp)
                    << std::setw(15) << std::sqrt(tmp / (fnorm * fnorm)) << std::endl;
                std::cout << std::endl;
            }
        }
    } else {
        // Non-standardized version. Needs additional operations
        while (iloop < optcontrol.maxnum_iteration) {
            do_print_log = !((iloop + 1) % optcontrol.output_frequency) && (verbosity > 1);

            if (do_print_log) {
                std::cout << "   Coordinate Descent : " << std::setw(5) << iloop + 1 << std::endl;
            }
            delta = beta;
            for (i = 0; i < N; ++i) {
                beta(i) = shrink(Minv * grad(i) + beta(i) / scale_beta(i), alphlambda) * scale_beta(i);
                delta(i) -= beta(i);
                if (std::abs(delta(i)) > 0.0) {
                    if (!has_prod[i]) {
                        for (j = 0; j < N; ++j) {
                            Prod(j, i) = A.col(j).dot(A.col(i));
                        }
                        has_prod[i] = true;
                    }
                    grad = grad + Prod.col(i) * delta(i);
                }
            }
            ++iloop;
            diff = std::sqrt(delta.dot(delta) / static_cast<double>(N));

            if (diff < optcontrol.tolerance_iteration) break;

            if (do_print_log) {
                std::cout << "    1: ||u_{k}-u_{k-1}||_2     = " << std::setw(15) << diff
                    << std::setw(15) << diff * std::sqrt(static_cast<double>(N) / beta.dot(beta)) << std::endl;
                auto tmp = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+:tmp)
#endif
                for (i = 0; i < N; ++i) {
                    tmp += std::abs(beta(i));
                }
                std::cout << "    2: ||u_{k}||_1             = " << std::setw(15) << tmp << std::endl;
                res = A * beta - b;
                tmp = res.dot(res);
                std::cout << "    3: ||Au_{k}-f||_2          = " << std::setw(15) << std::sqrt(tmp)
                    << std::setw(15) << std::sqrt(tmp / (fnorm * fnorm)) << std::endl;
                std::cout << std::endl;
            }
        }
    }

    if (verbosity > 1) {
        if (iloop >= optcontrol.maxnum_iteration) {
            std::cout << "WARNING: Convergence NOT achieved within " << optcontrol.maxnum_iteration
                << " coordinate descent iterations." << std::endl;
        } else {
            std::cout << "  Convergence achieved in " << iloop << " iterations." << std::endl;
        }
        const auto param2norm = beta.dot(beta);
        if (std::abs(param2norm) < eps) {
            std::cout << "    1': ||u_{k}-u_{k-1}||_2     = " << std::setw(15) << 0.0
                << std::setw(15) << 0.0 << std::endl;
        } else {
            std::cout << "    1': ||u_{k}-u_{k-1}||_2     = " << std::setw(15) << diff
                << std::setw(15) << diff * std::sqrt(static_cast<double>(N) / param2norm) << std::endl;
        }
        double tmp = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+:tmp)
#endif
        for (i = 0; i < N; ++i) {
            tmp += std::abs(beta(i));
        }
        std::cout << "    2': ||u_{k}||_1             = " << std::setw(15) << tmp << std::endl;
        res = A * beta - b;
        tmp = res.dot(res);
        std::cout << "    3': ||Au_{k}-f||_2          = " << std::setw(15) << std::sqrt(tmp)
            << std::setw(15) << std::sqrt(tmp / (fnorm * fnorm)) << std::endl;
        std::cout << std::endl;
    }

    for (i = 0; i < N; ++i) x[i] = beta(i);
}

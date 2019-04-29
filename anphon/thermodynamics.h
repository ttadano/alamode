/*
 phonon_thermodynamics.h

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include "pointers.h"
#include <complex>
#include <vector>
#include "kpoint.h"

namespace PHON_NS
{
    class temperature_information
    {
    public:
        double tmin{}, tmax{}, delta_t{};
        size_t number_of_grids{};
        std::vector<double> temperature_grid;

        temperature_information() = default;
        ~temperature_information() = default;

        temperature_information(const double tmin_,
                                const double tmax_,
                                const double delta_t_)
        {
            tmin = tmin_;
            tmax = tmax_;
            delta_t = delta_t_;
            number_of_grids = static_cast<unsigned int>((tmax - tmin) / delta_t) + 1;
            temperature_grid.resize(number_of_grids);

            for (auto i = 0; i < number_of_grids; ++i) {
                temperature_grid[i] = tmin + static_cast<double>(i) * delta_t;
            }
        }
    };

    class Thermodynamics : protected Pointers
    {
    public:
        Thermodynamics(class PHON *);

        ~Thermodynamics();

        double T_to_Ryd;
        bool classical;
        bool calc_FE_bubble;
        double *FE_bubble;

        void setup();
        void set_temperature_info(const double tmin_in,
                                  const double tmax_in,
                                  const double dt_in);

        temperature_information get_temperature_info() const;

        double Cv(const double omega,
                  const double temp_in) const;

        double Cv_classical(const double omega,
                            const double temp_in) const;

        double fB(const double omega,
                  const double temp_in) const;

        double fC(const double omega,
                  const double temp_in) const;

        double Cv_tot(const double temp_in,
                      const unsigned int nk_irred,
                      const unsigned int ns,
                      const std::vector<std::vector<KpointList>> &kp_irred,
                      double *weight_k_irred,
                      double **eval_in) const;

        double internal_energy(const double temp_in,
                               const unsigned int nk_irred,
                               const unsigned int ns,
                               const std::vector<std::vector<KpointList>> &kp_irred,
                               double *weight_k_irred,
                               double **eval_in) const;

        double vibrational_entropy(const double temp_in,
                                   const unsigned int nk_irred,
                                   const unsigned int ns,
                                   const std::vector<std::vector<KpointList>> &kp_irred,
                                   double *weight_k_irred,
                                   double **eval_in) const;

        double free_energy_QHA(const double temp_in,
                               const unsigned int nk_irred,
                               const unsigned int ns,
                               const std::vector<std::vector<KpointList>> &kp_irred,
                               double *weight_k_irred,
                               double **eval_in) const;


        double disp2_avg(double,
                         unsigned int,
                         unsigned int) const;

        double disp_corrfunc(const double T_in,
                             const unsigned int ncrd1,
                             const unsigned int ncrd2,
                             const double cell_shift[3],
                             const unsigned int nk,
                             const unsigned int ns,
                             double **xk_in,
                             double **eval_in,
                             std::complex<double> ***evec_in) const;

        double coth_T(double,
                      double) const;

        void compute_free_energy_bubble();


        void compute_FE_bubble(double **,
                               std::complex<double> ***,
                               double *) const;

        double compute_FE_bubble_SCPH(double,
                                      double **,
                                      std::complex<double> ***) const;

        double FE_scph_correction(unsigned int,
                                  double **,
                                  std::complex<double> ***) const;

    private:
        temperature_information tempinfo;
    };
}

/*
 constants.h

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include <cfloat>
#include <cmath>

static const double pi = 4.0*atan(1.0);
static const double amu = 1.660538782e-27;
static const double electron_mass = 9.10938215e-31;
static const double amu_ry = amu / electron_mass / 2.0;
static const double c_light = 299792458;
static const double h_planck = 6.62606896e-34;
static const double Ryd = 4.35974394e-18 / 2.0;
static const double time_ry = h_planck / (2.0 * pi) / Ryd;
static const double Hz_to_kayser = 1.0e-2 / (2.0 * pi * c_light);
static const double Bohr_in_Angstrom = 0.52917721092;
static const double k_Boltzmann = 1.3806488e-23; // J/K
static const double eps = DBL_EPSILON;
static const double eps15 = 1.0e-15;
static const double eps12 = 1.0e-12;
static const double eps10 = 1.0e-10;
static const double eps8 = 1.0e-8;
static const double eps6 = 1.0e-6;
static const double inverse_pi = 1.0 / pi;

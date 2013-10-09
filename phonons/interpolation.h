#pragma once

#include "pointers.h"
#include <vector>
#include <set>
#include <complex>

namespace PHON_NS {

	class Interpolation: protected Pointers {
	public:
		Interpolation(class PHON *);
		~Interpolation();

		void prepare_interpolation();

		int nk1, nk2, nk3;

	private:
		std::complex<double> ***mat_k;
		std::complex<double> ***mat_r;

		double ***damp;

		void setup_damping();

		void r2q(double *, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, std::complex<double> ***, std::complex<double> **);

	};	
}

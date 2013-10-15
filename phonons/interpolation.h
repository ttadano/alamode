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
		void parse_self_energy();

		bool is_interpolate_mode;

		int nk1, nk2, nk3;

	private:

		int nksym_ref;

		std::complex<double> ***mat_k;
		std::complex<double> ***mat_r;

		double ***damp;
		std::complex<double> ***matrix_k;

		double **xk_interpolate;
		unsigned int **k_reduced_interpolate;
		unsigned int nk_reduced_interpolate;

		void setup_damping();
		void setup_polarization_matrix();
		void create_matrix_k(unsigned int);

		void r2q(double *, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, std::complex<double> ***, std::complex<double> **);

	};	

	extern "C" {

		void zgemm_(const char *transa, const char *transb, int *m, int *n,
			int *k, std::complex<double> *alpha, std::complex<double> *a, int *lda,
			std::complex<double> *b, int *ldb, std::complex<double> *beta,
			std::complex<double> *c, int *ldc);

		void zgeev_(const char *jobvl, const char *jobvr, int *n, std::complex<double> *a, int *lda,
				std::complex<double> *w, std::complex<double> *vl, int *ldvl, std::complex<double> *vr, int *ldvr,
				std::complex<double> *work, int *lwork, double *rwork, int *info);
	}
}

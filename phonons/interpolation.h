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
		void exec_interpolation();
		void finish_interpolation();

		bool is_interpolate_mode;

		unsigned int nk1, nk2, nk3;
		unsigned int nk_ref;
		unsigned int nk_reduced_ref, nequiv_max_ref;
		int **k_reduced_ref;


	private:

		int nksym_ref;
		unsigned int ntemp_ref;
		std::vector<unsigned int> nk_equiv_ref;


		double **kvec_na_interpolate;

		std::complex<double> ***mat_r;
		std::complex<double> **self_energy;
		std::complex<double> ***self_energy_extend;

		double ***damp;

		double **eigval;
		std::complex<double> ***eigvec;

		double **xk_interpolate;
		unsigned int nk_reduced_interpolate;

		void setup_damping();
		void parse_self_energy();
		void prepare_dymat_for_interpolation();
		void prepare_self_energy_extend();
		void diagonalize_interpolated_matrix(std::complex<double> **, std::complex<double> *);

		void setup_polarization_matrix();
		void create_matrix_k(unsigned int, std::complex<double> ***);

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

	inline bool compare_real(const std::complex<double> &a, const std::complex<double> &b) 
	{
		return a.real() < b.real();
	}
}

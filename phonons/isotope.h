#include "mpi_common.h"
#include "pointers.h"


namespace PHON_NS
{
	class Isotope: protected Pointers {
	public:

		Isotope(class PHON *);
		~Isotope();

		bool include_isotope;
		double *isotope_factor;

		double **gamma_isotope;

		void setup_isotope_scattering();
		void calc_isotope_selfenergy_all();

	private:
		void calc_isotope_selfenergy(int, int, double, double &);
		void calc_isotope_selfenergy_tetra(int, int, double, double &);


	};
}
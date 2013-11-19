#pragma once

#include "pointers.h"
#include <string>
#include <fstream>
#include <vector>

#ifdef _USE_EIGEN
#include <Eigen/Core>
#endif

namespace ALM_NS {

	class SymmetryOperation {
	public:
		std::vector<int> symop;

		SymmetryOperation();

		// Declaration construction

		SymmetryOperation(const SymmetryOperation &a)
		{
			for(std::vector<int>::const_iterator p = a.symop.begin(); p != a.symop.end(); ++p){
				symop.push_back(*p);
			}
		}

		SymmetryOperation(const int rot[3][3], const int trans[3])
		{
			for (int i = 0; i < 3; ++i){
				for (int j = 0; j < 3; ++j){
					symop.push_back(rot[i][j]);
				}
			}
			for (int i = 0; i < 3; ++i){
				symop.push_back(trans[i]);
			}
		}
	};

	inline bool operator<(const SymmetryOperation a, const SymmetryOperation b){
		return std::lexicographical_compare(a.symop.begin(), a.symop.end(), b.symop.begin(), b.symop.end());
	}

	class SymmetryOperationTransFloat {
	public:
	    int rot[3][3];
		double tran[3];

		SymmetryOperationTransFloat();

		// Declaration construction


		SymmetryOperationTransFloat(const int rot_in[3][3], const double tran_in[3])
		{
			for (int i = 0; i < 3; ++i){
				for (int j = 0; j < 3; ++j){
					rot[i][j] = rot_in[i][j];
				}
			}
			for (int i = 0; i < 3; ++i){
				tran[i] = tran_in[i];
			}
		}
	};

	class RotationMatrix {
	public:
		int mat[3][3];

		RotationMatrix();
		RotationMatrix(const int rot[3][3]) {
			for (int i = 0; i < 3; ++i) {
				for (int j = 0; j < 3; ++j) {
					mat[i][j] = rot[i][j];
				}
			}
		}
	};

	class Symmetry: protected Pointers {
	public:
		Symmetry(class ALM *);
		~Symmetry();

		void init();
		void setup_symmetry_operation(int, unsigned int&, unsigned int&, double[3][3], double[3][3], 
			double **, int *);
// 		void findsym(int, int, int *, double[3][3], double[3][3],
// 			double **);

		void findsym(int, double [3][3], double **, std::vector<SymmetryOperation> &);

		unsigned int nsym, nnp;
		int ntran, natmin;
		int nsym_s, ntran_s, natmin_s; // for reference system (supercell?)
		int ntran_ref;

		int is_printsymmetry;
		int multiply_data;

		int ***symrel_int;
		int *symnum_tran;
		double tolerance;
		double ***symrel;
		double **tnons;

		int **map_sym;
		int **map_p2s;
		int **map_p2s_s;
		class Maps {
		public:
			int atom_num;
			int tran_num;
		};
		Maps *map_s2p, *map_s2p_s;

		void genmaps(int, double **, int **, int **, class Symmetry::Maps *);

		bool *sym_available;

		std::string file_sym, refsys_file;
		std::ofstream ofs_sym;
		std::ifstream ifs_sym;

		int numsymop(int, double **x, double);

	private:
// #ifdef _USE_EIGEN
// 		bool is_ortho(Eigen::Matrix3d, Eigen::Matrix3d, Eigen::Matrix3d);
// 		bool is_invariant(Eigen::Matrix3d, int, int*, double **, int[3], int);
// #endif
// 		bool is_ortho(double [3][3], double [3][3], double [3][3]);
// 		bool is_invariant(double [3][3], int, int *, double **, int [3], int);
		void matmul3(double [3][3], const double [3][3], const double [3][3]);
		void transpose3(double [3][3], const double [3][3]);
		void symop_in_cart(double [3][3], double[3][3]);
		void pure_translations();
		void data_multiplier(int, int, int);

		void print_symmetrized_coordinate(double **);
		void symop_availability_check(double ***, bool *, const int, int &);

		void find_lattice_symmetry(double [3][3], std::vector<RotationMatrix> &);
		void find_crystal_symmetry(int, int, std::vector<unsigned int> *, double **x, 
			std::vector<RotationMatrix>, std::vector<SymmetryOperationTransFloat> &);
		void find_nnp_for_translation(unsigned int &, std::vector<SymmetryOperationTransFloat>);

		std::vector<SymmetryOperation> SymmList;
	};


}

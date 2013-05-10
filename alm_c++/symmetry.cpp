#include <cmath>
#include <iostream>
#include <iomanip>
#include "symmetry.h"
#include "system.h"
#include "memory.h"
#include "constants.h"
#include "timer.h"
#include "error.h"
#include "files.h"
#include <cmath>
#include <vector>
#include <algorithm>

#ifdef _USE_EIGEN
#include <Eigen/Core>
#endif

using namespace ALM_NS;

Symmetry::Symmetry(ALM *alm) : Pointers(alm) 
{
	file_sym = "SYMM_INFO";
}

Symmetry::~Symmetry() {
	memory->deallocate(tnons);
	memory->deallocate(symrel_int);
	memory->deallocate(symrel);
	memory->deallocate(map_sym);
	memory->deallocate(map_p2s);
	memory->deallocate(map_s2p);
	memory->deallocate(symnum_tran);
}

void Symmetry::init()
{
	int i, j, k;
	int nat = system->nat;

	SymmList.clear();
	gensym(nat, nsym, nnp, system->lavec, system->rlavec, system->xcoord, system->kd);

	memory->allocate(tnons, nsym, 3);
	memory->allocate(symrel_int, nsym, 3, 3);

	int isym = 0;
	for (std::vector<SymmetryOperation>::iterator iter = SymmList.begin(); iter != SymmList.end(); ++iter){
		SymmetryOperation symop_tmp = *iter;
		for (i = 0; i < 3; ++i){
			for (j = 0; j < 3; ++j){
				symrel_int[isym][i][j] = symop_tmp.symop[3 * i + j];
			}
		}
		for (i = 0; i < 3; ++i){
			tnons[isym][i] = static_cast<double>(symop_tmp.symop[i + 9]) / static_cast<double>(nnp);
		}
		++isym;
	}
	SymmList.clear();

	std::cout << std::endl << "Number of symmetry operations = " << nsym << std::endl;

	timer->print_elapsed();

	memory->allocate(symrel, nsym, 3, 3);
	symop_in_cart(system->lavec, system->rlavec);

	memory->allocate(sym_available, nsym);
	int nsym_fc;
	symop_availability_check(symrel, sym_available, nsym, nsym_fc);

	std::cout << "Among " << nsym << " symmetries, " << nsym_fc << " symmetries will be used to reduce the number of parameters." << std::endl;
	std::cout << "Used symmetry operations (rotational part in Cartesian coordinate only) are printed below:" << std::endl;

	for (i = 0; i < nsym; ++i){
		if (sym_available[i]) {
			std::cout << std::setw(5) << i + 1;
			for (j = 0; j < 3; ++j){
				for (k = 0; k < 3; ++k){
					std::cout << std::setw(4) << symrel[i][j][k];
				}
			}
			std::cout << std::endl;
		}
	}
	std::cout << std::endl;

	pure_translations();

	memory->allocate(map_sym, nat, nsym);
	memory->allocate(map_p2s, natmin, ntran);
	memory->allocate(map_s2p, nat);

	genmaps(nat, system->xcoord, map_sym, map_p2s, map_s2p);

	std::cout << std::endl;
	std::cout << "**Cell-Atom Correspondens Below**" << std::endl;
	std::cout << std::setw(6) << "CELL" << " | " << std::setw(5) << "ATOM" << std::endl;

	for (int i = 0; i < ntran; ++i){
		std::cout << std::setw(6) << i + 1 << " | ";
		for (int j = 0; j < natmin; ++j)  {
			std::cout << std::setw(5) << map_p2s[j][i] + 1;
			if((j + 1)%5 == 0) {
				std::cout << std::endl << "       | ";
			}
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;

	if(multiply_data) data_multiplier(nat, system->ndata, multiply_data);

	timer->print_elapsed();
}

void Symmetry::gensym(int nat, int &nsym, int nnp, double aa[3][3], double bb[3][3], double **x, int *kd)
{
	int i, j;

	if(nsym == 0) {

		// Automatically find symmetries.

		std::cout << "Generating Symmetry Operations: This can take a while." << std::endl << std::endl;

		// findsym(nat, nnp, kd, aa, bb, x, nsym, rot, tran_int);
		findsym(nat, nnp, kd, aa, bb, x);
		nsym = SymmList.size();

		ofs_sym.open(file_sym.c_str(), std::ios::out);
		ofs_sym << nsym << std::endl;
		ofs_sym << nnp << std::endl;

		for(std::vector<SymmetryOperation>::iterator p = SymmList.begin(); p != SymmList.end(); ++p){
			SymmetryOperation sym_tmp = *p;
			for (i = 0; i < 9; ++i){
				ofs_sym << std::setw(4) << sym_tmp.symop[i];
			}
			ofs_sym << std::setw(7) << " ";
			for(i = 9; i < 12; ++i){
				ofs_sym << sym_tmp.symop[i] << std::setw(4);
			}
			ofs_sym << std::endl;
		}

		ofs_sym.close();
	} 
	else if(nsym == 1) {

		// Identity operation only !

		int rot_tmp[3][3], tran_tmp[3];

		for (i = 0; i < 3; ++i){
			for (j = 0; j < 3; ++j){
				if(i == j) {
					rot_tmp[i][j] = 1;
				} else {
					rot_tmp[i][j] = 0;
				}
			}
			tran_tmp[i] = 0;
		}

		SymmList.push_back(SymmetryOperation(rot_tmp, tran_tmp));
	} 
	else {
		int nsym2;
		int rot_tmp[3][3], tran_tmp[3];

		ifs_sym.open(file_sym.c_str(), std::ios::in);
		ifs_sym >> nsym2 >> nnp;

		if(nsym != nsym2) error->exit("gensym", "nsym in the given file and the input file are not consistent.");

		for (i = 0; i < nsym; ++i) {
			ifs_sym >> rot_tmp[0][0] >> rot_tmp[0][1] >> rot_tmp[0][2]
			>> rot_tmp[1][0] >> rot_tmp[1][1] >> rot_tmp[1][2] 
			>> rot_tmp[2][0] >> rot_tmp[2][1] >> rot_tmp[2][2]
			>> tran_tmp[0] >> tran_tmp[1] >> tran_tmp[2];

			SymmList.push_back(SymmetryOperation(rot_tmp, tran_tmp));
		}
		ifs_sym.close();
	}

#ifdef _DEBUG
	print_symmetrized_coordinate(x);
#endif
}

void Symmetry::findsym(int nat, int nnp, int *kd, double aa[3][3], double bb[3][3], double **x)
{
	int i, j;

	int m11, m12, m13, m21, m22, m23, m31, m32, m33;
	int det, np1, np2, np3;

#ifdef _USE_EIGEN
	Eigen::Matrix3d amat, bmat;
	Eigen::Matrix3d rot2;
#else
	double rot2[3][3];
#endif

	int rot_tmp[3][3], rot_reciprocal[3][3];
	int tran_tmp[3];

	for (i = 0; i < 3; ++i){
		for (j = 0; j < 3; ++j){

#ifdef _USE_EIGEN
			amat(i,j) = aa[i][j];
			bmat(i,j) = bb[i][j];
#endif

			if(i == j) {
				rot_tmp[i][j] = 1;
			} else {
				rot_tmp[i][j] = 0;
			}
		}
		tran_tmp[i] = 0;
	}

	// Add the identity operation to the list

	SymmList.push_back(SymmetryOperation(rot_tmp, tran_tmp));

	int nnps = nnp * nnp * nnp;
	int **arr_trans;

	memory->allocate(arr_trans, nnps, 3);

	int itran = 0;

	for (np1 = 0; np1 < nnp; ++np1){
		for (np2 = 0; np2 < nnp; ++np2){
			for (np3 = 0; np3 < nnp; ++np3){
				arr_trans[itran][0] = np1;
				arr_trans[itran][1] = np2;
				arr_trans[itran][2] = np3;
				++itran;
			}
		}
	}   

	for (m11 = -1; m11 <= 1; ++m11){
		for (m12 = -1; m12 <= 1; ++m12) {
			for (m13 = -1; m13 <= 1; ++m13){
				for (m21 = -1; m21 <= 1; ++m21){
					for (m22 = -1; m22 <= 1; ++m22){
						for (m23 = -1; m23 <= 1; ++m23){
							for (m31 = -1; m31 <= 1; ++m31){
								for (m32 = -1; m32 <= 1; ++m32){
									for (m33 = -1; m33 <= 1; ++m33){

										det = m11 * (m22 * m33 - m32 * m23)
											- m21 * (m12 * m33 - m32 * m13)
											+ m31 * (m12 * m23 - m22 * m13);

										if (det != 1 && det != -1) continue;

										rot_tmp[0][0] = m11;
										rot_tmp[0][1] = m12;
										rot_tmp[0][2] = m13;
										rot_tmp[1][0] = m21;
										rot_tmp[1][1] = m22;
										rot_tmp[1][2] = m23;
										rot_tmp[2][0] = m31;
										rot_tmp[2][1] = m32;
										rot_tmp[2][2] = m33;

										det = 1 / det;

										rot_reciprocal[0][0] = (m22 * m33 - m23 * m32) * det ;
										rot_reciprocal[0][1] = (m23 * m31 - m21 * m33) * det ;
										rot_reciprocal[0][2] = (m21 * m32 - m22 * m31) * det ;
										rot_reciprocal[1][0] = (m32 * m13 - m33 * m12) * det ;
										rot_reciprocal[1][1] = (m33 * m11 - m31 * m13) * det ;
										rot_reciprocal[1][2] = (m31 * m12 - m32 * m11) * det ;
										rot_reciprocal[2][0] = (m12 * m23 - m13 * m22) * det ;
										rot_reciprocal[2][1] = (m13 * m21 - m11 * m23) * det ;
										rot_reciprocal[2][2] = (m11 * m22 - m12 * m21) * det ;

										for (i = 0; i < 3; ++i) {
											for (j = 0; j < 3; ++j){
#ifdef _USE_EIGEN
												rot2(i,j) = static_cast<double>(rot_reciprocal[i][j]);
#else
												rot2[i][j] = static_cast<double>(rot_reciprocal[i][j]);
#endif
											}
										}

#ifdef _USE_EIGEN
										if(!is_ortho(rot2, amat, bmat)) continue;
#else
										if(!is_ortho(rot2, aa, bb)) continue;
#endif

#pragma omp parallel for private(np1, np2, np3)
										for (itran = 0; itran < nnps; ++itran){

											np1 = arr_trans[itran][0];
											np2 = arr_trans[itran][1];
											np3 = arr_trans[itran][2];

											if(m11 == 1 && m12 == 0 && m13 ==0 &&
												m21 == 0 && m22 == 1 && m23 == 0 &&
												m31 == 0 && m32 == 0 && m33 == 1 &&
												np1 == 0 && np2 == 0 && np3 == 0) continue;

											if(!is_invariant(rot2, nat, kd, x, arr_trans[itran], nnp)) continue;

											// STL containers are not thread-safe
#pragma omp critical
											SymmList.push_back(SymmetryOperation(rot_tmp, arr_trans[itran]));
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}

	std::sort(SymmList.begin() + 1, SymmList.end());
	memory->deallocate(arr_trans);
}

#ifdef _USE_EIGEN
bool Symmetry::is_ortho(Eigen::Matrix3d rot, Eigen::Matrix3d amat, Eigen::Matrix3d bmat)
{
	double pi2 = 2.0 * pi;

	Eigen::Matrix3d sat, unit;

	double tmp;

	sat = rot * amat.transpose();
	unit = (sat.transpose() * (bmat * (bmat.transpose() * sat)));
	unit /= pow(pi2, 2);

	tmp = pow((unit(0,0) - 1.0), 2) + pow((unit(1,1) - 1.0), 2) + pow((unit(2,2) - 1.0), 2)
		+ pow(unit(0,1), 2) + pow(unit(0,2), 2)
		+ pow(unit(1,0), 2) + pow(unit(1,2), 2)
		+ pow(unit(2,0), 2) + pow(unit(2,1), 2);

	if(tmp > eps) {
		return false;
	} else {
		return true;
	}
}

bool Symmetry::is_invariant(Eigen::Matrix3d rot, int nat, int *kd, double **x, int tran[3], int nnp)
{

	int i, j, k, l;
	Eigen::Vector3d wsi, usi, vsi, tmp;

	bool value = true;

	for (i = 0; i < nat; ++i){

		for (j = 0; j < 3; ++j){   
			wsi(j) = x[i][j] - static_cast<double>(tran[j]) / static_cast<double>(nnp);
		}

		usi = rot.transpose() * wsi;

		l = -1;

		for (j = 0; j < nat; ++j){

			if(kd[j] == kd[i]) {

				for (k = 0; k < 3; ++k) { 
					vsi(k) = x[j][k];
					tmp(k) = fmod(std::abs(usi(k) - vsi(k)), 1.0); 
					// need "std" to specify floating point operation
					// especially for intel compiler (there was no problem in MSVC)
					tmp(k) = std::min<double>(tmp(k), 1.0 - tmp(k)) ;
				}
				double diff = tmp.dot(tmp);
				if (diff < eps12) l = j;
			}
		}

		if(l == -1) value = false;

	}
	return value;
}
#endif

bool Symmetry::is_ortho(double rot[3][3], double amat[3][3], double bmat[3][3]) 
{
	int i, j;
	double pi2 = 2.0 * pi;
	double res;

	double sat[3][3], unit[3][3], tmp[3][3], tmp2[3][3];
	double amat_t[3][3], bmat_t[3][3], sat_t[3][3];

	transpose3(amat_t, amat);
	transpose3(bmat_t, bmat);

	matmul3(sat, rot, amat_t);
	transpose3(sat_t, sat);

	matmul3(tmp, bmat_t, sat);
	matmul3(tmp2, bmat, tmp);
	matmul3(unit, sat_t, tmp2);

	for (i = 0; i < 3; ++i) {
		for (j = 0; j < 3; ++j) {
			unit[i][j] /= pow(pi2, 2);
		}
	}

	res = pow((unit[0][0] - 1.0), 2) + pow((unit[1][1] - 1.0), 2) + pow((unit[2][2] - 1.0), 2)
		+ pow(unit[0][1], 2) + pow(unit[0][2], 2)
		+ pow(unit[1][0], 2) + pow(unit[1][2], 2)
		+ pow(unit[2][0], 2) + pow(unit[2][1], 2);

	if(res > eps) {
		return false;
	} else {
		return true;
	}
}

void Symmetry::matmul3(double ret[3][3], const double amat[3][3], const double bmat[3][3]) 
{
	int i, j, k;

	for (i = 0; i < 3; ++i) {
		for (j = 0; j < 3; ++j) {
			ret[i][j] = 0.0;
			for (k = 0; k < 3; ++k) ret[i][j] += amat[i][k] * bmat[k][j]; 	        
		}
	}
}

void Symmetry::transpose3(double ret[3][3], const double mat[3][3]) 
{
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			ret[i][j] = mat[j][i];
		}
	}
}


bool Symmetry::is_invariant(double rot[3][3], int nat, int *kd, double **x, int tran[3], int nnp) {

	int i, j, k, l;
	double wsi[3], usi[3], vsi[3], tmp[3];
	double diff;

	bool value = true;

	for (i = 0; i < nat; ++i) {

		for (j = 0; j < 3; ++j) {
			wsi[j] = x[i][j] - static_cast<double>(tran[j]) / static_cast<double>(nnp);
		}

		system->rotvec(usi, wsi, rot, 'T');

		l = -1;

		for (j = 0; j < nat; ++j) {

			if (kd[j] == kd[i]) {

				for (k = 0; k < 3; ++k) {
					vsi[k] = x[j][k];
					tmp[k] = std::fmod(std::abs(usi[k] - vsi[k]), 1.0);
					tmp[k] = std::min<double>(tmp[k], 1.0 - tmp[k]);
				}
				diff = tmp[0] * tmp[0] + tmp[1] * tmp[1] + tmp[2] * tmp[2];
				if (diff < eps12) l = j;
			}
		}

		if (l == -1) value = false;


	}
	return value;
}

void Symmetry::symop_in_cart(double lavec[3][3], double rlavec[3][3])
{
	int i, j;

#ifdef _USE_EIGEN
	Eigen::Matrix3d aa, bb, sym_tmp;
	Eigen::Matrix3d sym_crt;

	for (i = 0; i < 3; ++i){
		for (j = 0; j < 3; ++j){
			aa(i,j) = lavec[i][j];
			bb(i,j) = rlavec[i][j];
		}
	}

#else 
	double sym_tmp[3][3], sym_crt[3][3];
	double tmp[3][3];
#endif



	for (int isym = 0; isym < nsym; ++isym) {

		for (i = 0; i < 3; ++i){
			for (j = 0; j < 3; ++j){
#ifdef _USE_EIGEN
				sym_tmp(i,j) = static_cast<double>(symrel_int[isym][i][j]);
#else
				sym_tmp[i][j] = static_cast<double>(symrel_int[isym][i][j]);
#endif
			}
		}
#ifdef _USE_EIGEN
		sym_crt = (aa * (sym_tmp * bb)) / (2.0 * pi);
		for (i = 0; i < 3; ++i){
			for (j = 0; j < 3; ++j){
				symrel[isym][i][j] = sym_crt(i,j);
			}
		}
#else
		matmul3(tmp, sym_tmp, rlavec);
		matmul3(sym_crt, lavec, tmp);
		for (i = 0; i < 3; ++i){
			for (j = 0; j < 3; ++j){
				symrel[isym][i][j] = sym_crt[i][j] / (2.0 * pi);
			}
		}
#endif


	}

#ifdef _DEBUG

	std::cout << "Symmetry Operations in Cartesian Coordinate" << std::endl;
	for (int isym = 0; isym < nsym; ++isym){
		for (i = 0; i < 3; ++i){
			for (j = 0; j < 3; ++j){
				std::cout << std::setw(8) << symrel[isym][i][j];    
			}
		}
		std::cout << std::endl;
	}
#endif
}

void Symmetry::pure_translations()
{
	int i;

	ntran = 0;
	for(i = 0; i < nsym; ++i){
		if(symrel_int[i][0][0] == 1 && symrel_int[i][1][1] == 1 && symrel_int[i][2][2] == 1) {
			++ntran;
		}
	}

	natmin = system->nat / ntran;

	if(ntran > 1) {
		std::cout << "Given system is not primitive cell;" << std::endl;
		std::cout << std::setw(8) <<  ntran << " translation operations exist." << std::endl;
	} else {
		std::cout << "Given system is a primitive cell." << std::endl;
	}
	std::cout << "Each cell contains " << natmin << " atoms" << std::endl;

	memory->allocate(symnum_tran, ntran);

	int isym = 0;

	for (i = 0; i < nsym; ++i){
		if(symrel_int[i][0][0] == 1 && symrel_int[i][1][1] == 1 && symrel_int[i][2][2] == 1) {
			symnum_tran[isym++] = i; 
		}
	}
}

void Symmetry::genmaps(int nat, double **x, int **map_sym, int **map_p2s, Maps *map_s2p)
{
	double **xnew;
	int isym, iat, jat;
	int i;
	double tmp[3], dist; 

	memory->allocate(xnew, nat, 3);

	for(iat = 0; iat < nat; ++iat){
		for(isym = 0; isym < nsym; ++isym){
			map_sym[iat][isym] = -1;
		}
	}

	for (isym = 0; isym < nsym; ++isym){

		for (iat = 0; iat < nat; ++iat){

			for (i = 0; i < 3; ++i){
				xnew[iat][i] = static_cast<double>(symrel_int[isym][i][0]) * x[iat][0] 
				+ static_cast<double>(symrel_int[isym][i][1]) * x[iat][1] 
				+ static_cast<double>(symrel_int[isym][i][2]) * x[iat][2] 
				+ tnons[isym][i];
			}

			for (jat = 0; jat < nat; ++jat){

				for (i = 0; i < 3; ++i){
					tmp[i] = fmod(std::abs(xnew[iat][i] - x[jat][i]), 1.0);
					tmp[i] = std::min<double>(tmp[i], 1.0 - tmp[i]);
				}

				dist = tmp[0] * tmp[0] + tmp[1] * tmp[1] + tmp[2] * tmp[2];
				if(dist < eps10) {
					map_sym[iat][isym] = jat;
					break;
				}
			}
			if (map_sym[iat][isym] == -1) error->exit("genmaps", "cannot find symmetry for operation # ", isym + 1);
		}
	}
	memory->deallocate(xnew);    

	bool *is_checked;
	memory->allocate(is_checked, nat);

	for (i = 0; i < nat; ++i) is_checked[i] = false;

	jat = 0;
	int atomnum_translated;
	for (iat = 0; iat < nat; ++iat){

		if(is_checked[iat]) continue;
		for (i = 0; i < ntran; ++i){
			atomnum_translated = map_sym[iat][symnum_tran[i]];
			map_p2s[jat][i] = atomnum_translated;
			is_checked[atomnum_translated] = true;
		}
		++jat;
	}

	memory->deallocate(is_checked);

	for (iat = 0; iat < natmin; ++iat){
		for (i  = 0; i < ntran; ++i){
			atomnum_translated = map_p2s[iat][i];
			map_s2p[atomnum_translated].atom_num = iat;
			map_s2p[atomnum_translated].tran_num = i;
		}
	}
}

void Symmetry::data_multiplier(int nat, int ndata, int multiply_data)
{
	int i, j, k, itran;
	int n_map;
	double **u, **f;
	double ***u_sym, ***f_sym;

	files->ofs_disp_sym.open(files->file_disp_sym.c_str(), std::ios::out | std::ios::binary);
	files->ofs_force_sym.open(files->file_force_sym.c_str(), std::ios::out | std::ios::binary);

	if(!files->ofs_disp_sym)  error->exit("data_multiplier", "cannot open file_disp"); 
	if(!files->ofs_force_sym) error->exit("data_multiplier", "cannot open file_force");

	if(multiply_data == 3) {

		std::cout << "**Displacement-force data will be expanded to the bigger supercell**" << std::endl << std::endl;

		int nsym_ref = 0;
		unsigned int nnp_ref;

		// Read from reference file

		std::ifstream ifs_refsys;
		ifs_refsys.open(refsys_file.c_str(), std::ios::in);
		if(!ifs_refsys) error->exit("data_multiplier", "cannot open refsys_file");

		double lavec_ref[3][3], rlavec_ref[3][3];
		int nat_ref;

		for (i = 0; i < 3; ++i){
			for (j = 0; j < 3; ++j){
				ifs_refsys >> lavec_ref[j][i];
			}
		}
		system->recips(lavec_ref, rlavec_ref);
		ifs_refsys >> nat_ref >> nnp_ref;

		int *kd_ref, *map_ref;
		double **x_ref;

		memory->allocate(kd_ref, nat_ref);
		memory->allocate(x_ref, nat_ref, 3);

		for (i = 0; i < nat_ref; ++i){
			ifs_refsys >> kd_ref[i];
			for (j = 0; j < 3; ++j){
				ifs_refsys >> x_ref[i][j];
			}
		}
		ifs_refsys.close();

		// Generate symmetry operations of the reference system

		findsym(nat_ref, nnp_ref, kd_ref, lavec_ref, rlavec_ref, x_ref);

		nsym_ref = SymmList.size();

		double **tnons_ref;
		int ***symrel_int_ref;

		memory->allocate(tnons_ref, nsym_ref, 3);
		memory->allocate(symrel_int_ref, nsym_ref, 3, 3);

		int isym = 0;

		ntran_ref = 0;

		for (std::vector<SymmetryOperation>::iterator iter = SymmList.begin(); iter != SymmList.end(); ++iter){
			SymmetryOperation symop_tmp = *iter;
			for (i = 0; i < 3; ++i){
				for (j = 0; j < 3; ++j){
					symrel_int_ref[isym][i][j] = symop_tmp.symop[3 * i + j];
				}
			}
			for (i = 0; i < 3; ++i){
				tnons_ref[isym][i] = static_cast<double>(symop_tmp.symop[i + 9]) / static_cast<double>(nnp_ref);
			}

			if (symrel_int_ref[isym][0][0] == 1 && symrel_int_ref[isym][1][1] == 1 && symrel_int_ref[isym][2][2] == 1) {
				++ntran_ref;
			}

			++isym;
		}
		SymmList.clear();

		double **xnew, tmp[3];
		int **map_sym_ref;
		int *symnum_tran_ref;

		bool map_found;
		int iat, jat, icrd, jcrd;
		double dist;


		memory->allocate(symnum_tran_ref, ntran_ref);
		itran = 0;

		for (isym = 0; isym < nsym_ref; ++isym) {
			if (symrel_int_ref[isym][0][0] == 1 && symrel_int_ref[isym][1][1] == 1 && symrel_int_ref[isym][2][2] == 1) { 
				symnum_tran_ref[itran++] = isym;
			}
		}

		memory->allocate(xnew, nat_ref, 3);
		memory->allocate(map_sym_ref, nat_ref, nsym_ref);

		for(iat = 0; iat < nat_ref; ++iat){
			for(isym = 0; isym < nsym_ref; ++isym){
				map_sym_ref[iat][isym] = -1;
			}
		}

		for (isym = 0; isym < nsym_ref; ++isym){
			for (iat = 0; iat < nat_ref; ++iat){

				for (i = 0; i < 3; ++i){
					xnew[iat][i] = static_cast<double>(symrel_int_ref[isym][i][0]) * x_ref[iat][0] 
					+ static_cast<double>(symrel_int_ref[isym][i][1]) * x_ref[iat][1] 
					+ static_cast<double>(symrel_int_ref[isym][i][2]) * x_ref[iat][2] 
					+ tnons_ref[isym][i];
				}

				for (jat = 0; jat < nat_ref; ++jat){
					for (i = 0; i < 3; ++i){
						tmp[i] = fmod(std::abs(xnew[iat][i] - x_ref[jat][i]), 1.0);
						tmp[i] = std::min<double>(tmp[i], 1.0 - tmp[i]);
					}

					dist = tmp[0] * tmp[0] + tmp[1] * tmp[1] + tmp[2] * tmp[2];
					if(dist < eps10) {
						map_sym_ref[iat][isym] = jat;
						break;
					}
				}
				if (map_sym_ref[iat][isym] == -1) error->exit("data_multiplier", "cannot find symmetry for operation # ", isym + 1);
			}
		}
		memory->deallocate(xnew); 


		// Generate mapping information of larger supercell to smaller one.

		int **map_large_to_small;
		double xtmp[3], xdiff[3], xshift[3];

		memory->allocate(map_large_to_small, ntran_ref, nat);

		for (itran = 0; itran < ntran_ref; ++itran){

			for (icrd = 0; icrd < 3; ++icrd) {
				xshift[icrd] = x_ref[0][icrd] - x_ref[map_sym_ref[0][symnum_tran_ref[itran]]][icrd];
			}

			for (iat = 0; iat < nat; ++iat) {

				map_found = false;

				system->rotvec(xtmp, system->xcoord[iat], system->lavec);
				system->rotvec(xtmp, xtmp, rlavec_ref);

				for (icrd = 0; icrd < 3; ++icrd) {
					xtmp[icrd] /= 2.0 * pi;
				}

				for (jat = 0; jat < nat_ref; ++jat) {

					for (jcrd = 0; jcrd < 3; ++jcrd) {
						xdiff[jcrd] = xtmp[jcrd] - x_ref[jat][jcrd] - xshift[jcrd];
						xdiff[jcrd] = std::fmod(xdiff[jcrd], 1.0);
					}
					dist = xdiff[0]*xdiff[0] + xdiff[1]*xdiff[1] + xdiff[2]*xdiff[2];

					if (dist < eps12 && kd_ref[jat] == system->kd[iat]) {
						map_large_to_small[itran][iat] = jat;
						map_found = true;
						break;
					}
				}

				if (!map_found) error->exit("data_multiplier", "cannot find equivalent atom");
			}

			std::cout << "itran = " << itran << std::endl;
			for (iat = 0; iat < nat; ++iat) {
				std::cout << "iat = " << iat << " mapped = " << map_large_to_small[itran][iat] << std::endl;
			}
		}

		// Write mapped displacement-force data set

		memory->allocate(u_sym, ntran_ref, nat, 3);
		memory->allocate(f_sym, ntran_ref, nat, 3);

		memory->allocate(u, nat_ref, 3);
		memory->allocate(f, nat_ref, 3);

		for (i = 0; i < ndata; ++i){
			for (j = 0; j < nat_ref; ++j){
				files->ifs_disp  >> u[j][0] >> u[j][1] >> u[j][2];
				files->ifs_force >> f[j][0] >> f[j][1] >> f[j][2];
			}

			for (itran = 0; itran < ntran_ref; ++itran){
				for (j = 0; j < nat; ++j){
					for (k = 0; k < 3; ++k){
						n_map = map_large_to_small[itran][j];
						u_sym[itran][j][k] = u[n_map][k];
						f_sym[itran][j][k] = f[n_map][k];
					}
				}
			}

			for (itran = 0; itran < ntran_ref; ++itran){
				for (j = 0; j < nat; ++j){
					files->ofs_disp_sym.write((char *) &u_sym[itran][j][0], sizeof(double));
					files->ofs_disp_sym.write((char *) &u_sym[itran][j][1], sizeof(double));
					files->ofs_disp_sym.write((char *) &u_sym[itran][j][2], sizeof(double));
					files->ofs_force_sym.write((char *) &f_sym[itran][j][0], sizeof(double));
					files->ofs_force_sym.write((char *) &f_sym[itran][j][1], sizeof(double));
					files->ofs_force_sym.write((char *) &f_sym[itran][j][2], sizeof(double));
				}
			} 

		}
	//	error->exit("hoge", "hoge");

	} else if (multiply_data == 2) {

		int isym;

		memory->allocate(u, nat, 3);
		memory->allocate(f, nat, 3);

		memory->allocate(u_sym, nsym, nat, 3);
		memory->allocate(f_sym, nsym, nat, 3);

		for (i = 0; i < ndata; ++i) {
			for (j = 0; j < nat; ++j) {
				files->ifs_disp >> u[j][0] >> u[j][1] >> u[j][2];
				files->ifs_disp >> f[j][0] >> f[j][1] >> f[j][2];
			}

			for (isym = 0; isym < nsym; ++isym) {
				for (j = 0; j < nat; ++j) {
					n_map = map_sym[j][isym];

					for (k = 0; k < 3; ++k) {
						u_sym[isym][n_map][k] = u[j][k];
						f_sym[isym][n_map][k] = f[j][k];
					}
				}
			}

			for (isym = 0; isym < nsym; ++isym) {
				for (j = 0; j < nat; ++j) {
					files->ofs_disp_sym.write((char *) &u_sym[isym][j][0], sizeof(double));
					files->ofs_disp_sym.write((char *) &u_sym[isym][j][1], sizeof(double));
					files->ofs_disp_sym.write((char *) &u_sym[isym][j][2], sizeof(double));
					files->ofs_force_sym.write((char *) &f_sym[isym][j][0], sizeof(double));
					files->ofs_force_sym.write((char *) &f_sym[isym][j][1], sizeof(double));
					files->ofs_force_sym.write((char *) &f_sym[isym][j][2], sizeof(double));
				}
			}
		}

		std::cout << "Symmetrically equivalent displacements and forces data are" << std::endl;
		std::cout << "stored in files: " << files->file_disp_sym << " " << files->file_force_sym << std::endl;

	} else if (multiply_data == 1) {

		memory->allocate(u, nat, 3);
		memory->allocate(f, nat, 3);

		memory->allocate(u_sym, ntran, nat, 3);
		memory->allocate(f_sym, ntran, nat, 3);

		for (i = 0; i < ndata; ++i){
			for (j = 0; j < nat; ++j){
				files->ifs_disp >> u[j][0] >> u[j][1] >> u[j][2];
				files->ifs_force >> f[j][0] >> f[j][1] >> f[j][2];
			}


			for (itran = 0; itran < ntran; ++itran){
				for (j = 0; j < nat; ++j){
					n_map = map_sym[j][symnum_tran[itran]];

					for (k = 0; k < 3; ++k){
						u_sym[itran][n_map][k] = u[j][k];
						f_sym[itran][n_map][k] = f[j][k];
					}
				}
			}

			for (itran = 0; itran < ntran; ++itran){
				for (j = 0; j < nat; ++j){
					files->ofs_disp_sym.write((char *) &u_sym[itran][j][0], sizeof(double));
					files->ofs_disp_sym.write((char *) &u_sym[itran][j][1], sizeof(double));
					files->ofs_disp_sym.write((char *) &u_sym[itran][j][2], sizeof(double));
					files->ofs_force_sym.write((char *) &f_sym[itran][j][0], sizeof(double));
					files->ofs_force_sym.write((char *) &f_sym[itran][j][1], sizeof(double));
					files->ofs_force_sym.write((char *) &f_sym[itran][j][2], sizeof(double));
				}
			}
		}
		std::cout << "Symmetrically (Only translational part) equivalent displacements and forces data are" << std::endl;
		std::cout << "stored in files: " << files->file_disp_sym << " " << files->file_force_sym << std::endl;

	} else if (multiply_data == 0) {

		memory->allocate(u, nat, 3);
		memory->allocate(f, nat, 3);

		for (i = 0; i < ndata; ++i) {
			for (j = 0; j < nat; ++j) {
				files->ifs_disp >> u[j][0] >> u[j][1] >> u[j][2];
				files->ifs_force >> f[j][0] >> f[j][1] >> f[j][2];
			}
			for (j = 0; j < nat; ++j){
				files->ofs_disp_sym.write((char *) &u[j][0], sizeof(double));
				files->ofs_disp_sym.write((char *) &u[j][1], sizeof(double));
				files->ofs_disp_sym.write((char *) &u[j][2], sizeof(double));
				files->ofs_force_sym.write((char *) &f[j][0], sizeof(double));
				files->ofs_force_sym.write((char *) &f[j][1], sizeof(double));
				files->ofs_force_sym.write((char *) &f[j][2], sizeof(double));
			}
		}

		ntran = 1;
	} else if (multiply_data == -1) {
		ntran = 1;
	} 

	files->ofs_disp_sym.close();
	files->ofs_force_sym.close();

	memory->deallocate(u);
	memory->deallocate(f);
	memory->deallocate(u_sym);
	memory->deallocate(f_sym);
}

void Symmetry::print_symmetrized_coordinate(double **x)
{
	int i, j, k, l;
	int isym = 0;
	double **x_symm, **x_avg;
	int nat = system->nat;
	int m11, m12, m13, m21, m22, m23, m31, m32, m33;
	int det;

#ifdef _USE_EIGEN
	Eigen::Matrix3d rot;
	Eigen::Vector3d wsi, usi, vsi, tmp;
#else 
	double rot[3][3];
	double wsi[3], usi[3], vsi[3], tmp[3];
#endif

	int tran[3];

	memory->allocate(x_symm, nat, 3);
	memory->allocate(x_avg, nat, 3);

	for (i = 0; i < nat; ++i){
		for (j = 0; j < 3; ++j){
			x_avg[i][j] = 0.0;
		}
	}

	for (std::vector<SymmetryOperation>::iterator p = SymmList.begin(); p != SymmList.end(); ++p){
		SymmetryOperation symm_tmp = *p;

		++isym;
		std::cout << "Symmetry No. : " << std::setw(5) << isym << std::endl;

		m11 = symm_tmp.symop[0];
		m12 = symm_tmp.symop[1];
		m13 = symm_tmp.symop[2];
		m21 = symm_tmp.symop[3];
		m22 = symm_tmp.symop[4];
		m23 = symm_tmp.symop[5];
		m31 = symm_tmp.symop[6];
		m32 = symm_tmp.symop[7];
		m33 = symm_tmp.symop[8];

		det = m11 * (m22 * m33 - m32 * m23)
			- m21 * (m12 * m33 - m32 * m13)
			+ m31 * (m12 * m23 - m22 * m13);

#ifdef _USE_EIGEN
		rot(0,0) = static_cast<double>((m22 * m33 - m23 * m32) * det);
		rot(0,1) = static_cast<double>((m23 * m31 - m21 * m33) * det);
		rot(0,2) = static_cast<double>((m21 * m32 - m22 * m31) * det);
		rot(1,0) = static_cast<double>((m32 * m13 - m33 * m12) * det);
		rot(1,1) = static_cast<double>((m33 * m11 - m31 * m13) * det);
		rot(1,2) = static_cast<double>((m31 * m12 - m32 * m11) * det);
		rot(2,0) = static_cast<double>((m12 * m23 - m13 * m22) * det);
		rot(2,1) = static_cast<double>((m13 * m21 - m11 * m23) * det);
		rot(2,2) = static_cast<double>((m11 * m22 - m12 * m21) * det);
#else
		rot[0][0] = static_cast<double>((m22 * m33 - m23 * m32) * det);
		rot[0][1] = static_cast<double>((m23 * m31 - m21 * m33) * det);
		rot[0][2] = static_cast<double>((m21 * m32 - m22 * m31) * det);
		rot[1][0] = static_cast<double>((m32 * m13 - m33 * m12) * det);
		rot[1][1] = static_cast<double>((m33 * m11 - m31 * m13) * det);
		rot[1][2] = static_cast<double>((m31 * m12 - m32 * m11) * det);
		rot[2][0] = static_cast<double>((m12 * m23 - m13 * m22) * det);
		rot[2][1] = static_cast<double>((m13 * m21 - m11 * m23) * det);
		rot[2][2] = static_cast<double>((m11 * m22 - m12 * m21) * det);
#endif

		for (i = 9; i < 12; ++i){
			tran[i - 9] = symm_tmp.symop[i];
		}

		for (i = 0; i < nat; ++i){

			for (j = 0; j < 3; ++j){   
#ifdef _USE_EIGEN
				wsi(j) = x[i][j] - static_cast<double>(tran[j]) / static_cast<double>(nnp);
#else 
				wsi[j] = x[i][j] - static_cast<double>(tran[j]) / static_cast<double>(nnp);
#endif
			}

#ifdef _USE_EIGEN
			usi = rot * wsi;
#else
			system->rotvec(usi, wsi, rot);
#endif

			l = -1;

			for (j = 0; j < nat; ++j){
				for (k = 0; k < 3; ++k) {
#ifdef _USE_EIGEN
					vsi(k) = x[j][k];
					tmp(k) = fmod(std::abs(usi(k) - vsi(k)), 1.0); 
					// need "std" to specify floating point operation
					// especially for intel compiler (there was no problem in MSVC)
					tmp(k) = std::min<double>(tmp(k), 1.0 - tmp(k)) ;
#else
					vsi[k] = x[j][k];
					tmp[k] = std::fmod(std::abs(usi[k] - vsi[k]), 1.0);
					tmp[k] = std::min<double>(tmp[k], 1.0 - tmp[k]);
#endif
				}
#ifdef _USE_EIGEN
				double diff = tmp.dot(tmp);
#else
				double diff = tmp[0]*tmp[0] + tmp[1]*tmp[1] + tmp[2]*tmp[2];
#endif
				if (diff < eps12) l = j;
			}

			for (j = 0; j < 3; ++j){
#ifdef _USE_EIGEN
				x_symm[l][j] = usi(j);
#else 
				x_symm[l][j] = usi[j];
#endif
				do {
					if (x_symm[l][j] < 0.0) {
						x_symm[l][j] += 1.0;
					} else if (x_symm[l][j] > 1.0){
						x_symm[l][j] -= 1.0;
					}
				} while(x_symm[l][j] < 0.0 || x_symm[l][j] > 1.0);
			}

		}

		for (i = 0; i < nat; ++i){
			for (j = 0; j < 3; ++j){
				std::cout << std::setw(15) << std::fixed << x_symm[i][j];
			}
			std::cout << " ( ";
			for (j = 0; j < 3; ++j){
				std::cout << std::setw(15) << std::scientific << x_symm[i][j]-x[i][j];
			}
			std::cout << " )" << std::endl;

			for (j = 0; j < 3; ++j){
				x_avg[i][j] += x_symm[i][j];
			}
		}

	}

	for (i = 0; i < nat; ++i){
		for (j = 0; j < 3; ++j){
			x_avg[i][j] /= static_cast<double>(SymmList.size());
		}
	}

	std::cout << "Symmetrically Averaged Coordinate" << std::endl;
	for (i = 0; i < nat; ++i){
		for (j = 0; j < 3; ++j){
			std::cout << std::setw(15) << std::fixed << std::setprecision(9) << x_avg[i][j];
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;

	std::cout.setf(std::ios::floatfield);

	memory->deallocate(x_symm);
	memory->deallocate(x_avg);
}

void Symmetry::symop_availability_check(double ***rot, bool *flag, const int n, int &nsym_fc)
{
	int i, j, k;
	int nfinite;

	nsym_fc = 0;

	for (i = 0; i < nsym; ++i){

		nfinite = 0;
		for (j = 0; j < 3; ++j){
			for (k = 0; k < 3; ++k){
				if(std::abs(rot[i][j][k]) > eps12) ++nfinite;
			}
		}

		if (nfinite == 3){
			++nsym_fc;
			flag[i] = true;
		} else {
			flag[i] = false;
		}
	}
}

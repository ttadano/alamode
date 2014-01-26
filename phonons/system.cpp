#include "mpi_common.h"
#include "system.h"
#include "fcs_phonon.h"
#include "error.h"
#include "memory.h"
#include "constants.h"
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "mathfunctions.h"

using namespace PHON_NS;

System::System(PHON *phon): Pointers(phon) {}

System::~System() {}

void System::setup()
{
	unsigned int i, j;

	double **xtmp;

	load_system_info();

	recips(lavec_s, rlavec_s);
	recips(lavec_p, rlavec_p);

	memory->allocate(xr_p, nat, 3);
	memory->allocate(xc, nat, 3);

	for (i = 0; i < nat; ++i){
		rotvec(xc[i], xr_s[i], lavec_s);
		rotvec(xr_p[i], xc[i], rlavec_p);
		for(j = 0; j < 3; ++j){
			xr_p[i][j] /=  2.0 * pi;
		}
	}

	if (mympi->my_rank == 0) {
		std::cout << "**Lattice Vectors**" << std::endl;
		std::cout.setf(std::ios::scientific);

		std::cout << " *Super Cell* " << std::endl << std::endl;
		std::cout << " " << lavec_s[0][0] << " " << lavec_s[1][0] << " " << lavec_s[2][0] << " : a1" << std::endl;
		std::cout << " " << lavec_s[0][1] << " " << lavec_s[1][1] << " " << lavec_s[2][1] << " : a2" << std::endl;
		std::cout << " " << lavec_s[0][2] << " " << lavec_s[1][2] << " " << lavec_s[2][2] << " : a3" << std::endl;
		std::cout << std::endl;

		std::cout << " " << rlavec_s[0][0] << " " << rlavec_s[0][1] << " " << rlavec_s[0][2] << " : b1" << std::endl;
		std::cout << " " << rlavec_s[1][0] << " " << rlavec_s[1][1] << " " << rlavec_s[1][2] << " : b2" << std::endl;
		std::cout << " " << rlavec_s[2][0] << " " << rlavec_s[2][1] << " " << rlavec_s[2][2] << " : b3" << std::endl;
		std::cout << std::endl << std::endl;

		std::cout << " *Primitive Cell* " << std::endl << std::endl;
		std::cout << " " << lavec_p[0][0] << " " << lavec_p[1][0] << " " << lavec_p[2][0] << " : a1" << std::endl;
		std::cout << " " << lavec_p[0][1] << " " << lavec_p[1][1] << " " << lavec_p[2][1] << " : a2" << std::endl;
		std::cout << " " << lavec_p[0][2] << " " << lavec_p[1][2] << " " << lavec_p[2][2] << " : a3" << std::endl;
		std::cout << std::endl;

		std::cout << " " << rlavec_p[0][0] << " " << rlavec_p[0][1] << " " << rlavec_p[0][2] << " : b1" << std::endl;
		std::cout << " " << rlavec_p[1][0] << " " << rlavec_p[1][1] << " " << rlavec_p[1][2] << " : b2" << std::endl;
		std::cout << " " << rlavec_p[2][0] << " " << rlavec_p[2][1] << " " << rlavec_p[2][2] << " : b3" << std::endl;
		std::cout << std::endl << std::endl;

		double vec_tmp[3][3];

		for (i = 0; i < 3; ++i){
			for (j = 0; j < 3; ++j){
				vec_tmp[i][j] = lavec_p[j][i];
			}
		}
		volume_p = volume(vec_tmp[0], vec_tmp[1], vec_tmp[2]);

		std::cout << " Unit cell volume = " << volume_p << " (a.u)^3" << std::endl;
		std::cout << " Number of Atoms in the supercell: " << nat << std::endl << std::endl;

		memory->allocate(xtmp, natmin, 3);

		for (i = 0; i < natmin; ++i) {
			rotvec(xtmp[i], xr_s[map_p2s[i][0]], lavec_s);
			rotvec(xtmp[i], xtmp[i], rlavec_p);
			for (j = 0; j < 3; ++j) xtmp[i][j] /= 2.0 * pi;
		}

		std::cout << "Atoms in the primitive cell" << std::endl;
		for (i = 0; i < natmin; ++i){
			std::cout << std::setw(6) << i + 1 << ":";
			for (j = 0; j < 3; ++j) {
				std::cout << std::setw(15) << xtmp[i][j];
			}
			std::cout << std::setw(5) << symbol_kd[kd[map_p2s[i][0]]] << std::endl;
		}
		std::cout << std::endl << std::endl;

		memory->deallocate(xtmp);

		std::cout << "Mass of atomic species:" << std::endl;
		for (i = 0; i < nkd; ++i) {
			std::cout << std::setw(4) << symbol_kd[i] << ":";
			std::cout << std::fixed << std::setw(12) << mass_kd[i] << std::endl;
		}
		std::cout << std::endl << std::endl;
	}

	// Atomic masses in Rydberg unit

	memory->allocate(mass, nat);
	for (i = 0; i < nat; ++i){
		mass[i] = mass_kd[kd[i]]*amu_ry;
	}
	MPI_Bcast(&Tmin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&Tmax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&dT, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&cell_dimension[0], 3, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
	MPI_Bcast(&volume_p, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);


	memory->allocate(kd_prim, natmin);

	for (i = 0; i < natmin; ++i) {
		kd_prim[i] = kd[map_p2s[i][0]];
	}
	setup_atomic_class(natmin, kd_prim);

	memory->deallocate(kd_prim);
}

void System::load_system_info()
{
	unsigned int i;
	int nkd_tmp;

	if (mympi->my_rank == 0) {
		std::string file_fcs = fcs_phonon->file_fcs;
		std::ifstream ifs_fcs;

		ifs_fcs.open(file_fcs.c_str(), std::ios::in);
		if(!ifs_fcs) error->exit("load_system_info", "cannot open file file_fcs");

		std::string str_tmp;

		while(!ifs_fcs.eof())
		{
			std::getline(ifs_fcs, str_tmp);
			if(str_tmp == "##SYSTEM INFO"){

				std::getline(ifs_fcs, str_tmp);

				for (i = 0; i < 3; ++i){
					ifs_fcs >> lavec_s[0][i] >> lavec_s[1][i] >> lavec_s[2][i];
				}
				ifs_fcs.ignore();
				std::getline(ifs_fcs, str_tmp);
				ifs_fcs >> nkd_tmp;

				if (nkd != nkd_tmp) error->exit("load_system_info", "NKD in the info file is not consistent with that given in the input file.");
				// 				memory->allocate(symbol_kd, nkd);
				// 				memory->allocate(mass_kd, nkd);

				// 				for (i = 0; i < nkd; ++i){
				// 					ifs_fcs >> str_tmp >> symbol_kd[i] >> mass_kd[i];
				// 				}

				ifs_fcs.ignore();
				std::getline(ifs_fcs, str_tmp);
				std::getline(ifs_fcs, str_tmp);
				ifs_fcs >> nat >> natmin >> ntran;

				memory->allocate(xr_s, nat, 3);
				memory->allocate(kd, nat);
				memory->allocate(map_p2s, natmin, ntran);
				memory->allocate(map_s2p, nat);

				unsigned int ikd, itran, icell;
				std::getline(ifs_fcs, str_tmp);
				std::getline(ifs_fcs, str_tmp);

				for (i = 0; i < nat; ++i){
					ifs_fcs >> str_tmp >> ikd >> xr_s[i][0] >> xr_s[i][1] >> xr_s[i][2] >> itran >> icell;
					kd[i] = ikd - 1;
					map_p2s[icell - 1][itran - 1] = i;
					map_s2p[i].atom_num = icell - 1;
					map_s2p[i].tran_num = itran - 1;
				}
			}
		}
		ifs_fcs.close();
	}

	MPI_Bcast(&lavec_s[0][0], 9, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&lavec_p[0][0], 9, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	MPI_Bcast(&nkd, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
	MPI_Bcast(&nat, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
	MPI_Bcast(&natmin, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
	MPI_Bcast(&ntran, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

	if (mympi->my_rank > 0){
		memory->allocate(mass_kd, nkd);
		memory->allocate(xr_s, nat, 3);
		memory->allocate(kd, nat);
		memory->allocate(map_p2s, natmin, ntran);
		memory->allocate(map_s2p, nat);
	}

	MPI_Bcast(&mass_kd[0], nkd, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&xr_s[0][0], 3*nat, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&kd[0], nat, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
	MPI_Bcast(&map_p2s[0][0], natmin*ntran, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

	// need to define MPI_MYTYPE for the exchange of user defined class ?

	for (i = 0; i < nat; ++i) {
		MPI_Bcast(&map_s2p[i], sizeof(map_s2p[i]), MPI_BYTE, 0, MPI_COMM_WORLD);
	}

}

void System::recips(double vec[3][3], double inverse[3][3])
{
	double det;
	det = vec[0][0] * vec[1][1] * vec[2][2] 
	+ vec[1][0] * vec[2][1] * vec[0][2] 
	+ vec[2][0] * vec[0][1] * vec[1][2]
	- vec[0][0] * vec[2][1] * vec[1][2] 
	- vec[2][0] * vec[1][1] * vec[0][2]
	- vec[1][0] * vec[0][1] * vec[2][2];

	if(std::abs(det) < eps12) {
		error->exit("recips", "Lattice Vector is singular");
	}

	double factor = 2.0 * pi / det;

	inverse[0][0] = (vec[1][1] * vec[2][2] - vec[1][2] * vec[2][1]) * factor;
	inverse[0][1] = (vec[0][2] * vec[2][1] - vec[0][1] * vec[2][2]) * factor;
	inverse[0][2] = (vec[0][1] * vec[1][2] - vec[0][2] * vec[1][1]) * factor;

	inverse[1][0] = (vec[1][2] * vec[2][0] - vec[1][0] * vec[2][2]) * factor;
	inverse[1][1] = (vec[0][0] * vec[2][2] - vec[0][2] * vec[2][0]) * factor;
	inverse[1][2] = (vec[0][2] * vec[1][0] - vec[0][0] * vec[1][2]) * factor;

	inverse[2][0] = (vec[1][0] * vec[2][1] - vec[1][1] * vec[2][0]) * factor;
	inverse[2][1] = (vec[0][1] * vec[2][0] - vec[0][0] * vec[2][1]) * factor;
	inverse[2][2] = (vec[0][0] * vec[1][1] - vec[0][1] * vec[1][0]) * factor;
}

double System::volume(double vec1[3], double vec2[3], double vec3[3])
{
	double vol;

	vol = std::abs(vec1[0]*(vec2[1]*vec3[2] - vec2[2]*vec3[1]) 
		+ vec1[1]*(vec2[2]*vec3[0] - vec2[0]*vec3[2]) 
		+ vec1[2]*(vec2[0]*vec3[1] - vec2[1]*vec3[0]));

	return vol;
}


void System::setup_atomic_class(unsigned int N, unsigned int *kd) {

	// This function can be modified when one needs to 
	// compute symmetry operations of spin polarized systems.

	unsigned int i;
	std::set<unsigned int> kd_uniq;
	kd_uniq.clear();

	for (i = 0; i < N; ++i) {
		kd_uniq.insert(kd[i]);
	}
	nclassatom = kd_uniq.size();

	memory->allocate(atomlist_class, nclassatom);

	for (i = 0; i < N; ++i) {
		int count = 0;
		for (std::set<unsigned int>::iterator it = kd_uniq.begin(); it != kd_uniq.end(); ++it)  {
			if (kd[i] == (*it)) {
				atomlist_class[count].push_back(i);
			}
			++count;
		}
	}

	kd_uniq.clear();
}

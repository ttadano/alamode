#include "mpi_common.h"
#include "dynamical.h"
#include "kpoint.h"
#include "fcs_phonon.h"
#include "system.h"
#include "memory.h"
#include "error.h"
#include "phonons.h"
#include "../alm_c++/constants.h"
#include <string>
#include <iomanip>
#include <fstream>
#include <algorithm>

#ifdef _USE_BOOST
#include <boost/lexical_cast.hpp>
#endif

using namespace PHON_NS;

Fcs_phonon::Fcs_phonon(PHON *phon): Pointers(phon) {
    maxorder = 2;
}

Fcs_phonon::~Fcs_phonon(){}

void Fcs_phonon::setup(std::string mode)
{
    unsigned int nat = system->nat;
    unsigned int natmin = system->natmin;
    memory->allocate(fc2, natmin, nat, 3, 3);

    unsigned int i, j, icrd, jcrd;
    for (i = 0; i < natmin; ++i){
        for (j = 0; j < nat; ++j){
            for (icrd = 0; icrd < 3; ++icrd){
                for (jcrd = 0; jcrd < 3; ++jcrd){
                    fc2[i][j][icrd][jcrd] = 0.0;
                }
            }
        }
    }

    if (mympi->my_rank == 0) load_fc2();
    // This is not necessary
    MPI_Bcast(&fc2[0][0][0][0], 9*natmin*nat, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (mode == "boltzmann" || mode == "gruneisen"){
        memory->allocate(force_constant, 2);

        if (mympi->my_rank == 0) load_fcs();

        MPI_Bcast_fc_class(2);

        if (mympi->my_rank == 0) {
            for (i = 0; i < maxorder; ++i){
                std::cout << "Number of non-zero IFCs for " << i + 2 << " order: ";
                std::cout << force_constant[i].size() << std::endl;
            }
        }
    }
}

void Fcs_phonon::load_fc2()
{
    unsigned int i, icrd, jcrd;
    unsigned int iat, jat;
    unsigned int inatmin;

    unsigned int len1, len2;

    std::string str_tmp, str1, str2;
    std::ifstream ifs_fcs;
    unsigned int nfc2;
    double fc2_tmp;

    ifs_fcs.open(file_fcs.c_str(), std::ios::in);

    while(!ifs_fcs.eof())
    {
        std::getline(ifs_fcs, str_tmp);
        if(str_tmp == "#FCS_HARMONIC"){
            ifs_fcs >> nfc2;
            ifs_fcs.ignore();
            std::getline(ifs_fcs, str_tmp);

            for (i = 0; i < nfc2; ++i){

                ifs_fcs >> fc2_tmp;

                ifs_fcs >> str1 >> str2;
                len1 = str1.size();
                len2 = str2.size();

                iat = atoi(str1.substr(0, len1 - 1).c_str()) - 1;
                jat = atoi(str2.substr(0, len2 - 1).c_str()) - 1;
                icrd = coordinate_index(str1[len1 - 1]);
                jcrd = coordinate_index(str2[len2 - 1]);

                inatmin = system->map_s2p[iat].atom_num;
                fc2[inatmin][jat][icrd][jcrd] = fc2_tmp;
            }
        }
    }
    ifs_fcs.close();
}

void Fcs_phonon::load_fcs()
{
    unsigned int i, iorder;
    unsigned int ifcs, nfcs;

    double val;

    unsigned int atmn, xyz;
    std::string str_int;
    unsigned int len;
    unsigned int *ind;

    Triplet tri_tmp;
    std::vector<unsigned int> ivec;
    std::vector<Triplet> tri_vec;

    bool flag_found;
    std::ifstream ifs_fcs;
    std::string str_tmp;
    std::string *flag_str;

    std::cout << "Reading harmonic and anharmonic force constants from the info file..." << std::endl;

    memory->allocate(flag_str, maxorder);
    ifs_fcs.open(file_fcs.c_str(), std::ios::in);

    memory->allocate(ind, maxorder + 1);

    for (iorder = 0; iorder < maxorder; ++iorder){

        ifs_fcs.clear();
        ifs_fcs.seekg(0, std::ios_base::beg);

        if (iorder == 0) {
            flag_str[iorder] = "#FCS_HARMONIC";
        } else {
#ifdef _USE_BOOST
            flag_str[iorder] = "#FCS_ANHARM" + boost::lexical_cast<std::string>(iorder + 2);
#else
			std::stringstream ss_tmp;
			ss_tmp << iorder + 2;
			flag_str[iorder] = "#FCS_ANHARM" + ss_tmp.str();
#endif
        }

        flag_found = false;

        while(!ifs_fcs.eof() && !flag_found)
        {
            std::getline(ifs_fcs, str_tmp);
            if (str_tmp == flag_str[iorder]) {

                flag_found = true;
                ifs_fcs >> nfcs;
                ifs_fcs.ignore();
                std::getline(ifs_fcs,str_tmp);

                for (ifcs = 0; ifcs < nfcs; ++ifcs){

                    ifs_fcs >> val;
                    ivec.clear();

                    for (i = 0; i < iorder + 2; ++i){
                        ifs_fcs >> str_int;
                        len = str_int.size();
                        atmn = atoi(str_int.substr(0, len - 1).c_str()) - 1;
                        xyz = coordinate_index(str_int[len - 1]);
                        ivec.push_back(3 * atmn + xyz);
                    }

                    if (std::abs(val) > eps) {
                        do {
                            tri_vec.clear();

                            for (i = 0; i < iorder + 2; ++i){
                                tri_tmp.atom = system->map_s2p[ivec[i] / 3].atom_num;
                                tri_tmp.cell = system->map_s2p[ivec[i] / 3].tran_num;
                                tri_tmp.xyz  = ivec[i] % 3;

                                tri_vec.push_back(tri_tmp);
                            }

                            force_constant[iorder].push_back(FcsClass(val, tri_vec));

                        } while (std::next_permutation(ivec.begin() + 1, ivec.end()));            
                    }
                }
            }
        }

        if (!flag_found){
            str_tmp = flag_str[iorder] + " flag not found in the info file";
            error->exit("load_fcs", str_tmp.c_str());
        }
    }

    ifs_fcs.close();

    memory->deallocate(ind);
    std::cout << "Done !" << std::endl;

}

unsigned int Fcs_phonon::coordinate_index(const char char_coord)
{
    unsigned int m;
    if(char_coord == 'x'){
        m = 0;
    } else if (char_coord == 'y'){
        m = 1;
    } else if (char_coord == 'z'){
        m = 2;
    } else {
        error->exit("coordinate_index", "invalid char_coord", char_coord);
    }
    return m;
}

void Fcs_phonon::MPI_Bcast_fc_class(const unsigned int N)
{
    unsigned int i;
    int j, k;
    int nelem;
    double *fcs_tmp;
    unsigned int ***ind;

    Triplet tri_tmp;
    std::vector<Triplet> tri_vec;

    for (i = 0; i < N; ++i){
        int len = force_constant[i].size();
        nelem = i + 2;

        MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);

        memory->allocate(fcs_tmp, len);
        memory->allocate(ind, len, nelem, 3);

        if (mympi->my_rank == 0) {
            for (j = 0; j < len; ++j){
                fcs_tmp[j] = force_constant[i][j].fcs_val;
                for (k = 0; k < nelem; ++k){
                    ind[j][k][0] = force_constant[i][j].elems[k].atom;
                    ind[j][k][1] = force_constant[i][j].elems[k].cell;
                    ind[j][k][2] = force_constant[i][j].elems[k].xyz;
                }
            }
        }

        MPI_Bcast(&fcs_tmp[0], len, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&ind[0][0][0], 3*nelem*len, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

        if (mympi->my_rank > 0) {
            force_constant[i].clear();

            for (j = 0; j < len; ++j){

                tri_vec.clear();

                for (k = 0; k < nelem; ++k){
                    tri_tmp.atom = ind[j][k][0];
                    tri_tmp.cell = ind[j][k][1];
                    tri_tmp.xyz  = ind[j][k][2];

                    tri_vec.push_back(tri_tmp);
                }
                force_constant[i].push_back(FcsClass(fcs_tmp[j], tri_vec));
            }
        }

        memory->deallocate(fcs_tmp);
        memory->deallocate(ind);
    }
}

#include "writes.h"
#include "system.h"
#include "interaction.h"
#include "memory.h"
#include "symmetry.h"
#include "error.h"
#include "files.h"
#include "fcs.h"
#include "fitting.h"
#include <boost/lexical_cast.hpp>
#include <fstream>


using namespace ALM_NS;

Writes::Writes(ALM *alm): Pointers(alm){}

Writes::~Writes() {}

void Writes::writeall()
{
    wrtfcs();


    ofs_info.open(files->file_info.c_str(), std::ios::out);
    if(!ofs_info) error->exit("writeall", "cannot open file_info");

    wrtmisc();

    ofs_info.close();
}

void Writes::wrtfcs()
{
    int i, j, k, l, m;

    int maxorder = interaction->maxorder;
    std::string *str_fcs;

    memory->allocate(str_fcs, maxorder);

    std::string str_tmp;

    std::ofstream ofs_fcs;
    ofs_fcs.open(files->file_fcs.c_str(), std::ios::out);
    if(!ofs_fcs) error->exit("openfiles", "cannot open fcs file");

    for (i = 0; i < maxorder; ++i){
        str_fcs[i] = "*FC" + boost::lexical_cast<std::string>(i + 2);
    }

    ofs_fcs <<  "********************Force Constants (FCs)********************" << std::endl;
    ofs_fcs <<  "!     Force Constants will be printed in atomic unit        !" << std::endl;
    ofs_fcs <<  "!     FC2: Ry/a0^2     FC3: Ry/a0^3     FC4: Ry/a0^4   etc. !" << std::endl;
    ofs_fcs <<  "!     FC?: Ry/a0^?                                          !" << std::endl;
    ofs_fcs <<  "!     a0= Bohr radius                                       !" << std::endl;
    ofs_fcs << "*************************************************************"  << std::endl << std::endl;
    ofs_fcs << "---------------Symmetrically Independent FCs---------------" << std::endl;
    ofs_fcs << " Global No." << "  Local No." << "            FCs" << "            Pairs"  << "        Distance (for IFC2)"<< std::endl;

    k = 0;

    ofs_fcs.setf(std::ios::scientific);

    for (i = 0; i < maxorder; ++i){

        m = 0;

        if(fcs->ndup[i].size() > 0) {

            ofs_fcs << std::endl << std::setw(6) << str_fcs[i] << std::endl;

            for (j = 0; j < fcs->ndup[i].size(); ++j){

                ofs_fcs << std::setw(6) << k + 1 << std::setw(6) << j + 1 << std::setw(16) <<  fitting->params[k];
                for (l = 0; l < i + 2; ++l){
                    ofs_fcs << std::setw(7) << fcs->easyvizint(fcs->fc_set[i][m].elems[l]);    
                }
                if(i==0) {
                ofs_fcs << std::setw(15) << interaction->distlist[fcs->fc_set[i][m].elems[0]/3][fcs->fc_set[i][m].elems[1]/3];
                }
                ofs_fcs << std::endl;
                m += fcs->ndup[i][j];
                ++k;
            }
        }
    }

    for (i = 0; i < maxorder; ++i){
        str_fcs[i] = "**FC" + boost::lexical_cast<std::string>(i + 2);
    }

    ofs_fcs << std::endl << std::endl
        << "---------------All FCs below---------------" << std::endl;

    int ip = 0;
    int id;

    for (i = 0; i < maxorder; ++i){

        id = 0;

        if(fcs->ndup[i].size() > 0){
            ofs_fcs << std::endl << std::setw(6) << str_fcs[i] << std::endl;

            for (unsigned int iuniq = 0; iuniq < fcs->ndup[i].size(); ++iuniq){

                str_tmp = "# FC" + boost::lexical_cast<std::string>(i + 2) + "_";
                str_tmp += boost::lexical_cast<std::string>(iuniq + 1);

                ofs_fcs << str_tmp << std::setw(6) << fcs->ndup[i][iuniq] << std::setw(16) << fitting->params[ip] << std::endl;

                for (j = 0; j < fcs->ndup[i][iuniq]; ++j){
                    ofs_fcs << std::setw(5) << j + 1 << std::setw(16) << fcs->fc_set[i][id].coef;
                    for (k = 0; k < i + 2; ++k){
                        ofs_fcs << std::setw(6) << fcs->easyvizint(fcs->fc_set[i][id].elems[k]);
                    }
                    ofs_fcs << std::endl;
                    ++id;
                }
                ofs_fcs << std::endl;
                ++ip;
            }

        }
    }

    memory->deallocate(str_fcs);
    ofs_fcs.close();

    std::cout << std::endl << "Force Constants are written to file: " << files->file_fcs << std::endl;
}

void Writes::wrtmisc(){

    // write miscellaneous information to file_info 
    // for subsequent calculations (phonons, md, alm)
    unsigned int i, j;



    ofs_info << "##SYSTEM INFO" << std::endl;
    ofs_info << "Lattice Vector (in Bohr unit)" << std::endl;
    for (j = 0; j < 3; ++j){
        for(i = 0; i < 3; ++i){
            ofs_info <<  std::setw(15) << system->lavec[j][i];
        }
        ofs_info << std::endl;
    }
    ofs_info << "Atomic Species" << std::endl;
    ofs_info << std::setw(6) << system->nkd << std::endl;

    for(i = 0; i < system->nkd; ++i){
        ofs_info << std::setw(6) << i + 1 << std::setw(5) << system->kdname[i] << std::setw(20) << system->mass_kd[i] << std::endl;
    }
    ofs_info << "Translational Symmetry Information" << std::endl;
    ofs_info << std::setw(6) << system->nat << std::setw(6) << symmetry->natmin << std::setw(6) << symmetry->ntran << std::endl;
    ofs_info << std::setw(11) << "'Atoms'" << std::setw(11) << "'Species'" 
        << std::setw(75) <<  "'Atomic Coordinates (Fractional)'                      " 
        << std::setw(15) << "'TRANSLATION'" << std::setw(15) << "'IN THE CELL'" << std::endl;
    for(i = 0; i < system->nat; ++i){
        ofs_info << std::setw(11) << i + 1 << std::setw(11) << system->kd[i];
        for(j = 0; j < 3; ++j){
            ofs_info << std::scientific << std::setprecision(16) << std::setw(25) << system->xcoord[i][j];
        }
        ofs_info << std::setw(15) << symmetry->map_s2p[i].tran_num + 1 
            << std::setw(15) << symmetry->map_s2p[i].atom_num + 1 << std::endl;
    }


    ofs_info << "##HARMONIC FORCE CONSTANTS" << std::endl;
    ofs_info << fcs->ndup[0].size() << std::endl;

    for(i = 0; i < fcs->ndup[0].size(); ++i){
        ofs_info << std::scientific << std::setprecision(16) << std::setw(25) << fitting->params[i] << std::endl;
    }

    int k, iat;

    ofs_info << "##INTERACTION LISTS" << std::endl;
    ofs_info << "Interaction List and Reference Vectors(Cartesian) for each order" << std::endl;

    for(int order = 0; order < interaction->maxorder; ++order){
        ofs_info << "#LIST_" + interaction->str_order[order] << std::endl;

        for(k = 0; k < symmetry->natmin; ++k) ofs_info << std::setw(6) << interaction->ninter[k][order];
        ofs_info << std::endl;

        for(k = 0; k < symmetry->natmin; ++k){
            iat = symmetry->map_p2s[k][0];
            for(int m = 0; m < interaction->ninter[k][order]; ++m){
                ofs_info << std::setw(6) << iat + 1 << std::setw(6) << interaction->intpairs[k][order][m] + 1;
                for(i = 0; i < 3; ++i){
                    ofs_info << std::scientific << std::setprecision(16) << std::setw(25) << interaction->relvec[k][order][m][i];
                }
                ofs_info << std::endl;
            }
        }
    }

    int *ncount;
    int ind_tmp;
    int id;

    int ip = 0;
    int ishift = 0;

    memory->allocate(ncount, 3*symmetry->natmin);

    ofs_info << "##FORCE CONSTANTS" << std::endl;
    ofs_info << "All force constants and interaction info" << std::endl;
    for(int order = 0; order < interaction->maxorder; ++order){
        ofs_info << "#FCS_" + interaction->str_order[order] << std::endl;

        int nelem = 0;
        for(std::vector<int>::iterator it = fcs->ndup[order].begin(); it != fcs->ndup[order].end(); ++it){
            nelem += *it;
        }
        ofs_info << std::setw(10) << nelem << std::endl;

        for(i = 0; i < 3*symmetry->natmin; ++i) ncount[i] = 0;

        id = 0;

        for(i = 0; i < fcs->ndup[order].size(); ++i){
            for(j = 0; j < fcs->ndup[order][i]; ++j){
                ind_tmp = fcs->fc_set[order][id].elems[0];
                for(k = 0; k < symmetry->natmin; ++k){
                    if(ind_tmp / 3 == symmetry->map_p2s[k][0]) {
                        ++ncount[3 * k + ind_tmp % 3];
                        break;
                    }
                }
                ++id;
            }
        }

        for(i = 0; i < 3*symmetry->natmin; ++i){
            ofs_info << std::setw(6) << ncount[i];
        }
        ofs_info << std::endl;

      
        //id = 0;
        //for(i = 0; i < fcs->ndup[order].size(); ++i){
        //    for(j = 0; j < fcs->ndup[order][i]; ++j){
        //        ofs_info << std::scientific << std::setprecision(16) << std::setw(25) << fitting->params[ip]*fcs->fc_set[order][id].coef << std::endl;
        //        for(k = 0; k < order + 2; ++k){
        //            ofs_info << std::setw(5) << fcs->easyvizint(fcs->fc_set[order][id].elems[k]);
        //        }
        //        ofs_info << std::endl;
        //        ++id;
        //    }
        //    ++ip;
        //}

        // this sorting is necessary for linking to molecular dynamics program.
        std::sort(fcs->fc_set[order].begin(), fcs->fc_set[order].end());

        for(std::vector<FcProperty>::iterator it = fcs->fc_set[order].begin(); it != fcs->fc_set[order].end(); ++it){
            FcProperty fctmp = *it;
            ip = fctmp.mother + ishift;
            ofs_info << std::scientific << std::setprecision(16) << std::setw(25) << fitting->params[ip]*fctmp.coef << std::endl;
            for(k = 0; k < order + 2; ++k){
                ofs_info << std::setw(5) << fcs->easyvizint(fctmp.elems[k]);
            }
            ofs_info << std::endl;
        }

        ishift += fcs->ndup[order].size();
    }
    memory->deallocate(ncount);

    std::cout << std::endl << "Miscellaneous information needed for post-process was stored to file: " << files->file_info << std::endl;
}

#include "fcs_phonon.h"
#include "system.h"
#include "memory.h"
#include "error.h"
#include <string>
#include <fstream>

using namespace PHON_NS;

Fcs_phonon::Fcs_phonon(PHON *phon): Pointers(phon) {}

Fcs_phonon::~Fcs_phonon(){}

void Fcs_phonon::setup()
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

    load_fc2();
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
		//                fc2[jat][inatmin][jcrd][icrd] = fc2_tmp;

            }
        }
    }
    ifs_fcs.close();
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

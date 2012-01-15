#include "symmetry.h"
#include "system.h"
#include "memory.h"
#include "constants.h"

using namespace ALM_NS;

Symmetry::Symmetry(ALM *alm) : Pointers(alm) 
{
    maxsym = 10000;
}

Symmetry::~Symmetry() {}

void Symmetry::init()
{
    int nat = system->nat;

    memory->allocate(tnons, 3, maxsym);
    memory->allocate(symrel, 3, 3, maxsym);

    gensym(nat, nsym, nnp, system->lavec, system->rlavec, system->xcoord,
        system->kd, symrel_int, tnons);
}

void Symmetry::gensym(int nat, int nsym, int nnp,
    double aa[3][3], double bb[3][3], double **x, int *kd, int ***rot, double **tran)
{
    if(nsym == 0) {
        // Automatically find symmetries.

        int **tran_int;
        memory->allocate(tran_int, 3, maxsym);
        findsym(nat, nnp, kd, aa, bb, x, nsym, rot, tran_int);
    } 
    else if(nsym == 1) {
        // Just identity operation.

        for (int i = 0; i < 3; i++){
            for (int j = 0; j < 3; j++){
                if(i == j) {
                    rot[i][j][0] = 1;
                } else {
                    rot[i][j][0] = 0;
                }
            }
            tran[i][0] = 0.0;
        }
    } 
    else {
    }
}

void Symmetry::findsym(int nat, int nnp, int *kd, double aa[3][3], double bb[3][3],
    double **x, int nsym, int ***rot, int **tran_int)
{
    // Symmetry Finder (originally from TAPP code)

    int m11, m12, m13, m21, m22, m23, m31, m32, m33;
    int det, np1, np2, np3;

    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
            if(i == j) {
                rot[i][j][0] = 1;
            } else {
                rot[i][j][0] = 0;
            }
        }
        tran_int[i][0] = 0;
    }

    nsym = 1;

    int rot_tmp[3][3], rot_reciprocal[3][3];
    int tran_tmp[3];

    for (m11 = -1; m11 <= 1; m11++){
        for (m12 = -1; m12 <= 1; m12++) {
            for (m13 = -1; m13 <= 1; m13++){
                for (m21 = -1; m21 <= 1; m21++){
                    for (m22 = -1; m22 <= 1; m22++){
                        for (m23 = -1; m23 <= 1; m23++){
                            for (m31 = -1; m31 <= 1; m31++){
                                for (m32 = -1; m32 <= 1; m32++){
                                    for (m33 = -1; m33 <= 1; m33++){

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

                                        rot_reciprocal[0][0] = (m22 * m33 - m23 * m32) / det ;
                                        rot_reciprocal[0][1] = (m23 * m31 - m21 * m33) / det ;
                                        rot_reciprocal[0][2] = (m21 * m32 - m22 * m33) / det ;
                                        rot_reciprocal[1][0] = (m32 * m13 - m33 * m12) / det ;
                                        rot_reciprocal[1][1] = (m33 * m11 - m31 * m13) / det ;
                                        rot_reciprocal[1][2] = (m31 * m12 - m32 * m11) / det ;
                                        rot_reciprocal[2][0] = (m12 * m23 - m13 * m22) / det ;
                                        rot_reciprocal[2][1] = (m13 * m21 - m11 * m23) / det ;
                                        rot_reciprocal[2][2] = (m11 * m22 - m12 * m21) / det ;

                                        if(!is_ortho(rot_reciprocal, aa, bb)) continue;

                                        for (np1 = 0; np1 < nnp; np1++){
                                            for (np2 = 0; np2 < nnp; np2++){
                                                for (np3 = 0; np3 < nnp; np3++){

                                                    if(m11 == 1 && m12 == 0 && m13 ==0 &&
                                                        m21 == 0 && m22 == 1 && m23 == 0 &&
                                                        m31 == 0 && m32 == 0 && m33 == 1 &&
                                                        np1 == 0 && np2 == 0 && np3 == 0) continue;

                                                    tran_tmp[0] = np1;
                                                    tran_tmp[1] = np2;
                                                    tran_tmp[2] = np3;

                                                    if(!is_invariant(rot_reciprocal, nat, kd, x, tran_tmp, nnp)) continue;

                                                    nsym++;

                                                    for (int i = 0; i < 3; i++){
                                                        for (int j = 0; j < 3; j++){
                                                            rot[i][j][nsym] = rot_tmp[i][j];
                                                        }
                                                        tran_int[i][nsym] = tran_tmp[i];
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
            }
        }
    }
}

bool Symmetry::is_ortho(int rot[3][3], double aa[3][3], double bb[3][3])
{
    double pi2 = 2.0 * pi;

    int i, j;
    double rot2[3][3];
    double sat[3][3], unit[3][3];

    // upcasting
    for (i = 0; i < 3; i++){
        for (j = 0; j < 3; j++){
            rot2[i][j] = static_cast<double>(rot[i][j]);
        }
    }


}
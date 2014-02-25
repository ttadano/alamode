#pragma once

#include "pointers.h"
#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include "constants.h"

namespace PHON_NS {

    class KpointList {
    public:
        std::vector<double> kval;
        unsigned int knum;

        KpointList(){};
        KpointList(const KpointList &obj){
            knum = obj.knum;
            for (std::vector<double>::const_iterator it = obj.kval.begin(); it != obj.kval.end(); ++it){
                kval.push_back(*it);
            }
        }

        KpointList(const unsigned int knum_in, const std::vector<double> vec)
        {
            knum = knum_in;
            for (std::vector<double>::const_iterator it = vec.begin(); it != vec.end(); ++it){
                kval.push_back(*it);
            }
        }

    };

    class KpointInp {
    public:
        std::vector<std::string> kpelem;

        KpointInp(){};

        KpointInp(const std::vector<std::string> &obj) {
            for (std::vector<std::string>::const_iterator it = obj.begin(); it != obj.end(); ++it) {
                kpelem.push_back(*it);
            }
        }
    };

    class KpointPlane {
    public:
        double k[3];
        int n[2];

        KpointPlane(){};

        KpointPlane(double *xk_in, int *n_in) {
            for (int i = 0; i < 3; ++i) k[i] = xk_in[i];
            for (int i = 0; i < 2; ++i) n[i] = n_in[i];
        }
    };

    class Kpoint: protected Pointers {
    public:
        Kpoint(class PHON *);
        ~Kpoint();

        void kpoint_setups(std::string);

        int kpoint_mode;
        unsigned int nkx, nky, nkz;
        unsigned int nk;
        unsigned int *knum_minus;

        double **xk;
        double *kaxis;
        double **kvec_na;

        std::vector<KpointInp> kpInp;
        std::vector<double> weight_k;
        std::vector<std::vector<KpointList> > kpoint_irred_all;

        unsigned int nplanes;
        std::vector<KpointPlane> *kp_planes; 
        unsigned int nk_reduced;
        std::map<int, int> kmap_to_irreducible;


        int get_knum(const double, const double, const double);

        void generate_irreducible_kmap(int *, unsigned int &, std::vector<int> &,
            const unsigned int, const unsigned int, const unsigned int, 
            double **, const int, int ***);
        void gen_kmesh(const bool, const unsigned int [3], double **, std::vector<std::vector<KpointList> > &);


    private:
        void setup_kpoint_given(std::vector<KpointInp> &, unsigned int &, double **&, double **&);
        void setup_kpoint_band(std::vector<KpointInp> &, unsigned int &, double **&, double **&, double *&);
        void setup_kpoint_mesh(std::vector<KpointInp> &, unsigned int &, unsigned int &, unsigned int &, unsigned int &,
            double **&, double **&, const bool, std::vector<std::vector<KpointList> > &);
        void setup_kpoint_plane(std::vector<KpointInp> &, unsigned int &, std::vector<KpointPlane> *&);

        void reduce_kpoints(const unsigned int, double **, const unsigned int [3], std::vector<std::vector<KpointList> > &);
        void gen_nkminus(const unsigned int, unsigned int *, double **);
        void gen_kpoints_plane(std::vector<KpointInp>, std::vector<KpointPlane> *);
        bool in_first_BZ(double *);

        void mpi_broadcast_kpoint_vector(std::vector<std::vector<KpointList> > &);
        void mpi_broadcast_kplane_vector(const unsigned int, std::vector<KpointPlane> *&);
    };
}

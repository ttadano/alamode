#ifndef ALM_FILES_HEADER
#define ALM_FILES_HEADER

#include <string>
#include <fstream>
#include "pointers.h"

namespace ALM_NS {
    class Files : protected Pointers {
    public:
        Files(class ALM *);
        ~Files();

        void init();
        void openfile(std::string, std::ifstream);
        void openfile(std::string, std::ofstream);

        std::string job_title;
        std::string file_fcs, file_premd;

        std::string file_disp, file_force;
        std::string file_disp_sym, file_force_sym;

        std::string file_int;

        std::ofstream ofs_log;
        std::ofstream ofs_fcs, ofs_premd;

        std::ofstream ofs_int;
        std::ifstream ifs_int;

        std::ifstream ifs_disp, ifs_force;
        std::ofstream ofs_disp_sym, ofs_force_sym;
        std::ifstream ifs_disp_sym, ifs_force_sym;

    private:
        void openfiles();
        void closefiles();
        void setfilenames();
    };
}
#endif
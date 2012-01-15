#ifndef ALM_FILES_HEADER
#define ALM_FILES_HEADER

#include "pointers.h"
#include <string>
#include <fstream>

namespace ALM_NS {
    class Files : protected Pointers {
    public:
        Files(class ALM *);
        ~Files();

        void init();
        void openfiles();
        void closefiles();
        void setfilenames();

  //  private:
        std::string job_title;
        std::string file_log, file_fcs, file_premd;

        std::string file_disp, file_force;
        std::string file_disp_sym, file_force_sym;

        std::ofstream ofs_log;
        std::ofstream ofs_fcs, ofs_premd;

        std::ifstream ifs_disp, ifs_force;
        std::ofstream ofs_disp_sym, ofs_force_sym;

    };
}
#endif
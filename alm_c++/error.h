#ifndef ALM_ERROR_HEADER
#define ALM_ERROR_HEADER

#include "pointers.h"

namespace ALM_NS {
    class Error : protected Pointers{
    public:
        Error(class ALM *);
        ~Error();

        void exit(const char *, const char *);
        void warn(const char *, const char *);
        void exit(const char *, const char *, int);
    };
}

#endif
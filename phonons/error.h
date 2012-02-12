#ifndef ALM_ERROR_HEADER
#define ALM_ERROR_HEADER

#include <string>
#include "pointers.h"

namespace PHON_NS {
    class Error : protected Pointers{
    public:
        Error(class PHON *);
        ~Error();

        void exit(const char *, const char *);
        void warn(const char *, const char *);
        void exit(const char *, const char *, int);
    };
}

#endif
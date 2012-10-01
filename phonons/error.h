#pragma once

#include <string>
#include "pointers.h"

namespace PHON_NS {
    class Error : protected Pointers{
    public:
        Error(class PHON *);
        ~Error();

        void warn(const char *, const char *);
        void exit(const char *, const char *);        
        void exit(const char *, const char *, int);
        void exitall(const char *, const char *);        
    };
}

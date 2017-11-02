/*
 error.h

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include <string>
#include "pointers.h"

namespace PHON_NS
{
    class Error : protected Pointers
    {
    public:
        Error(class PHON *);
        ~Error();

        void warn(const char *, const char *);
        void exit(const char *, const char *);
        void exit(const char *, const char *, int);
        void exitall(const char *, const char *);
        void exit(const char *, const char *, const char *);
    };
}

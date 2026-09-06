/*
 main.cpp

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include "phonon_cui.h"
#include "error.h"
#include "ndarray.h"

using namespace PHON_NS;

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    // an allocation failure on any rank aborts the whole job (see include/ndarray.h)
    ndarray_detail::ndarray_out_of_memory_hook() = PHON_NS::abort_all_ranks;

    const auto phon_cui = new PhononCUI();

    phon_cui->run(argc, argv, MPI_COMM_WORLD);

    delete phon_cui;

    MPI_Finalize();

    return 0;
}

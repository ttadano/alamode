/* Declaration of pointers used in the whole program. */

#include "mpi.h"

namespace ALM_NS {

	class ALM {
	public:
		class Memory *memory; // Memory Allocation
		class Error *error;

	ALM(int, char **, MPI_Comm);
	~ALM();
	};
}
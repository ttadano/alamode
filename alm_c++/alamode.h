#include "mpi.h"

namespace ALM_NS {
	class ALM {
	public:
		class Memory *memory;
		ALM(int, char **, MPI_Comm);
		};
}

#include "constructors.h"
#include "file_work.h"
#include "potentials.h"
#include "supporting.h"
#include "constants.h"
#include "dynamics.h"


//#pragma comment(lib, "libname.lib")



int main(int argc, char* argv[]) 
{

	int rank, size;
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	
	moleculardynamic(500, 1000, 1.0);

	MPI_Finalize();
	return 0;

}


#include "constructors.h"
#include "file_work.h"
#include "montecarlo.h"
#include "potentials.h"
#include "supporting.h"
#include "constants.h"
#include "dynamics.h"


//#pragma comment(lib, "libname.lib")



int main(int argc, char* argv[]) // code, Gamma_numerator, Gamma_denominator, n_particles, n_iterations
								// code 0 -- ewald, code 1 -- sqpherical
{

	int rank, size;
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	//argv[0] = argc
	//argv[1] = code
	//argv[2] = Gamma_numerator
	//argv[3] = Gamma_denominator
	//argv[4] = n_particles
	//argv[5] = n_iterations

	//if (argc > 1)
	//{
		/*double GAMMA = atoi(argv[2]) * 1.0 / atoi(argv[3]);

		if (atoi(argv[1]) == 0)
		{
			montecarlo(atoi(argv[4]), atoi(argv[5]), GAMMA);
		}

		if (atoi(argv[1]) == 1)
		{
			montecarlo2(atoi(argv[4]), atoi(argv[5]), GAMMA);
		}*/

		/*double L = 0.0;
		std::vector<std::vector<double>> ions;
		L = random_state(ions, 100);
		print(ions);
		dump_step("dump.txt", ions, 15, L);*/

	moleculardynamic(500, 1000, 1.0);
		
	//}
	





	MPI_Finalize();
	return 0;

}
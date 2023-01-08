#include "constructors.h"

double simple_cubic(std::vector<std::vector<double>>& ions, int elementary_max)
{
	ions.clear();
	int n_particles = std::pow(elementary_max, 3);
	double L_calc_cell = std::pow(4 * PI * n_particles / 3.0, 1.0 / 3.0);
	double cube_length = L_calc_cell / elementary_max;

	for (auto i = 0u; i < elementary_max; i++)
	{
		for (auto j = 0u; j < elementary_max; j++)
		{
			for (auto k = 0u; k < elementary_max; k++)
			{
				ions.push_back({ i * cube_length, j * cube_length, k * cube_length });
			}
		}
	}
	return L_calc_cell;
}

double base_cubic(std::vector<std::vector<double>>& ions, int elementary_max)
{
	ions.clear();
	int n_particles = 2 * std::pow(elementary_max, 3);
	double L_calc_cell = std::pow(4 * PI * n_particles / 3.0, 1.0 / 3.0);
	double cube_length = L_calc_cell / elementary_max;
	for (auto i = 0u; i < elementary_max; i++)
	{
		for (auto j = 0u; j < elementary_max; j++)
		{
			for (auto k = 0u; k < elementary_max; k++)
			{
				ions.push_back({ i * cube_length, j * cube_length, k * cube_length });
				ions.push_back({ (i + 0.5) * cube_length, (j + 0.5) * cube_length, (k + 0.5) * cube_length });

			}
		}
	}
	return L_calc_cell;
}

double face_cubic(std::vector<std::vector<double>>& ions, int elementary_max)
{
	ions.clear();
	int n_particles = 4 * std::pow(elementary_max, 3);
	double L_calc_cell = std::pow(4 * PI * n_particles / 3.0, 1.0 / 3.0);
	double cube_length = L_calc_cell / elementary_max;

	for (auto i = 0u; i < elementary_max; i++)
	{
		for (auto j = 0u; j < elementary_max; j++)
		{
			for (auto k = 0u; k < elementary_max; k++)
			{
				ions.push_back({ i * cube_length, j * cube_length, k * cube_length });

				ions.push_back({ i * cube_length , (j + 0.5) * cube_length, (k + 0.5) * cube_length });

				ions.push_back({ (i + 0.5) * cube_length, (j + 0.5) * cube_length, k * cube_length });

				ions.push_back({ (i + 0.5) * cube_length , j * cube_length, (k + 0.5) * cube_length });

			}
		}
	}

	return L_calc_cell;
}

double random_state(std::vector<std::vector<double>>& ions, int n_particles)
{
	ions.clear();
	double L_calc_cell = std::pow(4 * PI * n_particles / 3.0, 1.0 / 3.0);
	
	
	size_t seed = 32;
	std::random_device rd{};
	std::mt19937 engine{ seed };
	std::uniform_real_distribution<double> dist{ 0.0, L_calc_cell };
	
	

	for (auto i = 0u; i < n_particles; i++)
	{
		double x = dist(engine);
		double y = dist(engine);
		double z = dist(engine);

		MPI_Bcast(&x, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&y, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&z, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		ions.push_back({ x ,y , z });
		
	}
	//int rank, size;
	//MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	//MPI_Comm_size(MPI_COMM_WORLD, &size);

	//for (int i = 0; i < size; i++)
	//{
	//	MPI_Barrier(MPI_COMM_WORLD);
	//	if (i == rank) {
	//		std::cout << "rank: " << rank << std::endl;
	//		print(ions);
	//		//MPI_Barrier(MPI_COMM_WORLD);
	//	}
	//	MPI_Barrier(MPI_COMM_WORLD);
	//}
	MPI_Barrier(MPI_COMM_WORLD);
	return L_calc_cell;
}
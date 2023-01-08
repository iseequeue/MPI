#include "dynamics.h"

void moleculardynamic(const int& n_particles, const int& n_iterations, const double& GAMMA)
{

	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	std::vector<std::vector<double>> ions;
	double L_calc_cell = random_state(ions, n_particles);


	Timer<std::chrono::seconds> t1("Molecular Dynamics");

	std::random_device rd{};
	std::mt19937 engine{ 5 };

	std::normal_distribution<> norm_dist{ 0,1 };

	std::vector<std::vector<double>> v_ions;
	std::vector<std::vector<double>> a_ions;

	//test example
	//double L_calc_cell = 20.0;
	//ions = { {15,10,10},{5,10,10} };
	//v_ions = { {-1.0,0.0,0.0}, {1.0, 0.0, 0.0} };
	//a_ions = { {0.0,0.0,0.0}, {0.0, 0.0, 0.0} };

	for (auto i = 0u; i < n_particles; i++)
	{
		double vx = norm_dist(engine);
		double vy = norm_dist(engine);
		double vz = norm_dist(engine);

		MPI_Bcast(&vx, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&vy, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(&vz, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		v_ions.push_back({ vx, vy, vz });
		a_ions.push_back({ 0.0, 0.0, 0.0 });
	}
	MPI_Barrier(MPI_COMM_WORLD);

	std::vector<std::vector<double>> v_dubbing_ions = v_ions;

	double r_max = L_calc_cell * std::pow(3.0 / (4.0 * PI), 1.0 / 3.0);
	std::vector<double> energy;
	std::vector<double> kinetic_energy;
	double potential = 0.0;
	double kinetic = 0.0;


	double T = 0.0;
	double dT = 0.01;
	double f = 0.0;



	for (auto step = 0u; step < n_iterations; step++)
	{

		for (auto i = 0u; i < ions.size(); i++)
		{
			for (auto j = 0u; j < 3; j++)
			{
				v_dubbing_ions[i][j] += v_ions[i][j] * a_ions[i][j] * dT / 2;
				ions[i][j] += v_dubbing_ions[i][j] * dT;
			}
		}

		potential = 0.0;
		std::vector<double> bias{ 0.0, 0.0, 0.0 };
		std::vector<double> distance_ij{ 0.0, 0.0, 0.0 };
		double f = 0.0;

		for (auto i = 0u; i < ions.size(); i++)
		{
			if (i % size == rank)
			{
				distance_ij = { 0.0, 0.0, 0.0 };
				double r_ij = 0.0;

				bias = { 0.5 * L_calc_cell - ions[i][0] , 0.5 * L_calc_cell - ions[i][1] , 0.5 * L_calc_cell - ions[i][2] };
				for (auto j = 0u; j < ions.size(); j++)
				{
					if (i == j) continue;

					{
						distance_ij = { pbc(ions[i][0],bias[0], L_calc_cell) - pbc(ions[j][0], bias[0], L_calc_cell),
										pbc(ions[i][1],bias[1], L_calc_cell) - pbc(ions[j][1], bias[1], L_calc_cell),
										pbc(ions[i][2],bias[2], L_calc_cell) - pbc(ions[j][2], bias[2], L_calc_cell) };
						r_ij = length(distance_ij);


						if (r_ij < r_max)
						{
							potential += spher_av_pot(r_ij, r_max);
							f = 1 / std::pow(r_ij, 3) - 1 / std::pow(r_max, 3);
							a_ions[i][0] = -distance_ij[0] * f;
							a_ions[i][1] = -distance_ij[1] * f;
							a_ions[i][2] = -distance_ij[2] * f;


							if (r_ij > L_calc_cell - r_max)
							{
								for (auto q = 0u; q < UNITS.size(); q++)
								{
									r_ij = length({ distance_ij[0] + L_calc_cell * UNITS[q][0], distance_ij[1] + L_calc_cell * UNITS[q][1],
											distance_ij[2] + L_calc_cell * UNITS[q][2] });

									if (r_ij < r_max)
									{
										potential += spher_av_pot(r_ij, r_max);
										f = 1 / std::pow(r_ij, 3) - 1 / std::pow(r_max, 3);
										a_ions[i][0] = -distance_ij[0] * f;
										a_ions[i][1] = -distance_ij[1] * f;
										a_ions[i][2] = -distance_ij[2] * f;
									}
								}
							}
						}
					}
				}
				MPI_Bcast(&a_ions[i], 3, MPI_DOUBLE, rank, MPI_COMM_WORLD);
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
		double total_energy;
		MPI_Allreduce(&potential, &total_energy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		#potential = total_energy;
		kinetic = 0.0;

		for (auto i = 0u; i < ions.size(); i++)
		{
			for (auto j = 0u; j < 3; j++)
			{
				v_ions[i][j] = v_dubbing_ions[i][j] + a_ions[i][j] * dT / 2;
				kinetic += v_ions[i][j] * v_ions[i][j];
				if (ions[i][j] > L_calc_cell) ions[i][j] -= L_calc_cell;
				if (ions[i][j] < 0.0) ions[i][j] += L_calc_cell;
			}
		}

		energy.push_back((kinetic + total_energy) / 2);
		kinetic_energy.push_back(kinetic / 2);

		if (rank == 0)
		{
			dump_step("x.txt", ions, step, L_calc_cell);
			dump_step("v.txt", v_ions, step, L_calc_cell);
			dump_step("a.txt", a_ions, step, L_calc_cell);

		}
		if (rank == 0 && (step % 100 == 0)) std::cout << step << std::endl;
	}
	if (rank == 0)
	{
		write_vector("energy.txt", energy);
		write_vector("kinetic.txt", kinetic_energy);
	}
}

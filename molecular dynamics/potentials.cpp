#include "potentials.h"

double v1(const std::vector<double>& vector_ij, const double& L)
{
	double distance_ij = length(vector_ij);

	return erfc(SQRTPI * distance_ij / L) / distance_ij - 1 / L;
}

double v2(const std::vector<double>& vector_ij, const double& L_calc_cell)
{
	double energy = 0.0;
	double d = 0.0; //helPIng variable |r/n + L|
	int n2 = 0; //helPIng variable n^2
	double current_energy = 0.0;


	for (auto nx = -EWALD_MAX; nx <= EWALD_MAX; nx++)
	{
		for (auto ny = -EWALD_MAX; ny <= EWALD_MAX; ny++)
		{
			for (auto nz = -EWALD_MAX; nz <= EWALD_MAX; nz++)
			{

				n2 = nx * nx + ny * ny + nz * nz;
				d = length({ vector_ij[0] / L_calc_cell + nx, vector_ij[1] / L_calc_cell + ny, vector_ij[2] / L_calc_cell + nz });
				if (n2 == 0) continue;

				current_energy = erfc(SQRTPI * d) / d + exp(-PI * n2) / (PI * n2)
					* cos(2 * PI * scalar(vector_ij, { nx * 1.0, ny * 1.0, nz * 1.0 }) / L_calc_cell);
				energy += current_energy;


			}
		}
	}
	return energy / L_calc_cell;
}

double V0(const std::vector<std::vector<double>>& ions, const double &GAMMA)
{
	return -0.5 * GAMMA * MSC * std::pow(ions.size(), 2.0 / 3.0);
}

double u(const int& i, const std::vector<std::vector<double>>& ions, const double& L_calc_cell)
{
	double energy = 0.0;
	for (auto j = 0u; j < ions.size(); j++)
	{
		if (i == j) continue;
		energy += v1({ ions[i][0] - ions[j][0], ions[i][1] - ions[j][1], ions[i][2] - ions[j][2] }, L_calc_cell)
			+ v2({ ions[i][0] - ions[j][0], ions[i][1] - ions[j][1], ions[i][2] - ions[j][2] }, L_calc_cell);
	}

	return energy;
}



double full_energy(const std::vector<std::vector<double>>& ions, const double& L_calc_cell, const double& GAMMA)
{
	double energy = 0.0;

	for (auto i = 0u; i < ions.size(); i++)
	{

		energy += u(i, ions, L_calc_cell);

	}
	return V0(ions, GAMMA) + GAMMA * energy / 2;
}



double tricky_energy(const std::vector<std::vector<double>>& ions_origin, const double& L_calc_cell, const double& GAMMA)
{
	double energy = 0.0;
	std::vector<double> bias{ 0.0, 0.0, 0.0 };

	for (auto i = 0u; i < ions_origin.size(); i++)
	{
		std::vector<std::vector<double>> ions = ions_origin;

		bias = { 0.5 * L_calc_cell - ions[i][0] , 0.5 * L_calc_cell - ions[i][1] , 0.5 * L_calc_cell - ions[i][2] };

		for (auto k = 0u; k < ions.size(); k++)
		{
			for (auto p = 0u; p < 3; p++)
			{
				ions[k][p] = pbc(ions[k][p], bias[p], L_calc_cell);
			}
		}



		for (auto j = 0u; j < ions.size(); j++)
		{
			if (i == j) continue;
			energy += (v1({ ions[i][0] - ions[j][0], ions[i][1] - ions[j][1], ions[i][2] - ions[j][2] }, L_calc_cell)
				+ v2({ ions[i][0] - ions[j][0], ions[i][1] - ions[j][1], ions[i][2] - ions[j][2] }, L_calc_cell));
		}
	}

	return 0.5 * energy * GAMMA + V0(ions_origin, GAMMA);
}



double spher_av_pot(const double& r_ij, const double& r_max)
{
	double x = r_ij / r_max;
	return (1 + 0.5 * x * (x * x - 3.0)) / r_ij;
}

std::vector<double> u_spher(const int& i, const std::vector<std::vector<double>>& ions, const double& L_calc_cell)
{
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);	


	double energy = 0.0;
	double virial = 0.0;
	double r_max = L_calc_cell * std::pow(3.0 / (4.0 * PI), 1.0 / 3.0);
	std::vector<double> bias{ 0.0, 0.0, 0.0 };

	std::vector<double> distance_ij{ 0.0, 0.0, 0.0 };
	double r_ij = 0.0;
	
	bias = { 0.5 * L_calc_cell - ions[i][0] , 0.5 * L_calc_cell - ions[i][1] , 0.5 * L_calc_cell - ions[i][2] };

	for (auto j = 0u; j < ions.size(); j++)
	{
		if (i == j) continue;

		if (j % size == rank)
		{
			distance_ij = { pbc(ions[i][0],bias[0], L_calc_cell) - pbc(ions[j][0], bias[0], L_calc_cell),
							pbc(ions[i][1],bias[1], L_calc_cell) - pbc(ions[j][1], bias[1], L_calc_cell),
							pbc(ions[i][2],bias[2], L_calc_cell) - pbc(ions[j][2], bias[2], L_calc_cell) };
			r_ij = length(distance_ij);


			if (r_ij < r_max)
			{
				energy += spher_av_pot(r_ij, r_max);
				virial += 1 / r_ij - r_ij * r_ij / std::pow(r_max, 3);

				if (r_ij > L_calc_cell - r_max)
				{
					for (auto q = 0u; q < UNITS.size(); q++)
					{
						r_ij = length({ distance_ij[0] + L_calc_cell * UNITS[q][0], distance_ij[1] + L_calc_cell * UNITS[q][1],
								distance_ij[2] + L_calc_cell * UNITS[q][2] });

						if (r_ij < r_max)
						{
							energy += spher_av_pot(r_ij, r_max);

							virial += 1 / r_ij - r_ij * r_ij / std::pow(r_max, 3);
						}
					}
				}
			}
		}		
	}
	MPI_Barrier(MPI_COMM_WORLD);
	double total_energy;
	MPI_Allreduce(&energy, &total_energy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
	energy = total_energy;

	double total_w;
	MPI_Allreduce(&virial, &total_w, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	virial = total_w;
	
	return { energy, virial };
}

std::vector<double> spherical_energy(const std::vector<std::vector<double>>& ions_origin, const double& L_calc_cell, const double& GAMMA)
{

	int Ns = 0;

	std::vector<double> ew = { 0.0, 0.0 };
	std::vector<double> ew_current = { 0.0, 0.0 };

	for (auto i = 0u; i < ions_origin.size(); i++)
	{
		ew_current = u_spher(i, ions_origin, L_calc_cell);
		ew[0] += ew_current[0];
		ew[1] += ew_current[1];
	}

	return { -0.15 * GAMMA * std::pow(ions_origin.size(), 2.0 / 3) * (ions_origin.size() + 5.0) + 0.5 * GAMMA * ew[0] , ew[1]};
}



std::vector<double> ewald_i(const int& i, const std::vector<double> &coord_i, const std::vector<std::vector<double>>& ions, const double& L_calc_cell)
{
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);


	double energy = 0.0;
	double virial = 0.0;


	for (auto j = 0u; j < ions.size(); j++)
	{
		if (i == j) continue;
		if (j % size == rank)
		{
			energy += v1({ coord_i[0] - ions[j][0], coord_i[1] - ions[j][1], coord_i[2] - ions[j][2] }, L_calc_cell)
				+ v2({ coord_i[0] - ions[j][0], coord_i[1] - ions[j][1], coord_i[2] - ions[j][2] }, L_calc_cell);

			virial += w_ewald({ coord_i[0] - ions[j][0], coord_i[1] - ions[j][1], coord_i[2] - ions[j][2] }, L_calc_cell);
		}
	}


	MPI_Barrier(MPI_COMM_WORLD);
	double total_energy;
	MPI_Allreduce(&energy, &total_energy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	energy = total_energy;

	double total_w;
	MPI_Allreduce(&virial, &total_w, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	virial = total_w;

	return { energy, virial };
}

std::vector<double> ewald_energy(const std::vector<std::vector<double>>& ions_origin, const double& L_calc_cell, const double& GAMMA)
{

	std::vector<double> ew = { 0.0, 0.0 };
	std::vector<double> ew_current = { 0.0, 0.0 };

	for (auto i = 0u; i < ions_origin.size(); i++)
	{
		ew_current = ewald_i(i, ions_origin[i], ions_origin, L_calc_cell);
		ew[0] += ew_current[0];
		ew[1] += ew_current[1];
	}

	return { V0(ions_origin, GAMMA) + 0.5 * GAMMA * ew[0] , ew[1] };
}







double w_ewald(const std::vector<double>& r_ij, const double& L)
{
	double r = length(r_ij);

	double angular = 0.0;

	std::vector<double> grad = { 0.0, 0.0, 0.0 };

	for (auto q = 0u; q < grad.size(); q++)
	{
		grad[q] = (-2 * r / L * exp(-PI * r * r / (L * L)) - erfc(SQRTPI * r / L)) * r_ij[q] / (r * r * r);

		for (auto nx = -EWALD_MAX; nx <= EWALD_MAX; nx++)
		{
			for (auto ny = -EWALD_MAX; ny <= EWALD_MAX; ny++)
			{
				for (auto nz = -EWALD_MAX; nz <= EWALD_MAX; nz++)
				{
					if (nx * nx + ny * ny + nz * nz == 0) continue;

					angular = 0.0;

					std::vector<double> N = { nx * 1.0, ny * 1.0, nz * 1.0 };
					double n = length(N);

					double d = length({ r_ij[0] / L + N[0], r_ij[1] / L + N[1], r_ij[2] / L + N[2] });

					angular += (-2 * d * exp(-PI * d * d) - erfc(SQRTPI * d)) * (r_ij[q] / L + N[q]) / (L * d * d * d);

					angular += -exp(-PI * n * n) / (PI * n * n) * sin(2 * PI / L * scalar(r_ij, N)) * 2 * PI / L * N[q];

					grad[q] += angular / L;

				}
			}
		}

	}

	return -scalar(r_ij, grad) / 3;

}

double w_i_ewald(const int& i, const std::vector<std::vector<double >>& ions, const double& L)
{
	double w = 0.0;
	for (auto j = 0u; j < ions.size(); j++)
	{
		if (i == j) continue;
		w += w_ewald({ ions[i][0] - ions[j][0], ions[i][1] - ions[j][1], ions[i][2] - ions[j][2] }, L);
	}
	return w;
}

double W_ewald(const std::vector<std::vector<double>>& ions, const double& L)
{
	double w = 0;
	for (auto i = 0u; i < ions.size(); i++)
	{
		w += w_i_ewald(i, ions, L);
	}
	return w / 2;
}
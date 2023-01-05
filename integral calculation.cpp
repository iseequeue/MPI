#include <cstdio>
#include "mpi.h"
#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>
#include <iomanip> 


double f(const double& x)
{
	return 4/(1+x*x);
}

int main(int argc, char** argv)
{
	int n = atoi(argv[1]); // number of splits

	int rank, size;
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	double sum = 0.0, a = 0.0, b = 1.0;
	double h = (b - a) / n;
	double x = a + h / 2.0 + rank * h * n / size;

	double time_work = 0;
	if (rank == 0)
	{
		time_work = MPI_Wtime();
	}

	for (int i = 1; i <= n / size; i++)   //summation at all processes
	{
		sum = sum + f(x) * h;
		x = x + h;
	}

	if (rank == 0)
	{
		for (int i = 1; i < size; i++)  // collect data
		{
			MPI_Recv(&a, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &status);
			sum = sum + a;
		}
		time_work = MPI_Wtime() - time_work;
		std::cout << "number of iterations = " << n << std::endl;
		std::cout << " The Program is RUN on " << size << " CPU" << std::endl;
		std::cout << "Integral = " << sum << std::endl;
		std::cout << "Time     = " << time_work << std::endl;
	}
	else
	{
		MPI_Send(&sum, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
	}

	MPI_Finalize();
	return 0;

}
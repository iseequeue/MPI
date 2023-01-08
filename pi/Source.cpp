#include <cstdio>
#include "mpi.h"
#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>
#include <iomanip> 




// Temporary array for slave process

double factorial(int N)
{
	double f = 1.0;
	if (N <= 1) return f;
	for (int i = 1; i <= N ; i++)
		f *= i;
	return f;
}


int main(int argc, char* argv[])
{

	const int n = 10000;
	double a[n];  //1/factorial array
	double a2[10000]; //local arrays
	

	int pid, np,
		elements_per_process,
		n_elements_recieved, start;
	// np -> no. of processes
	// pid -> process id

	MPI_Status status;

	// Creation of parallel processes
	MPI_Init(&argc, &argv);

	// find out process ID,
	// and how many processes were started
	MPI_Comm_rank(MPI_COMM_WORLD, &pid);
	MPI_Comm_size(MPI_COMM_WORLD, &np);

	double time_work = 0;
	// master process
	if (pid == 0)
	{
		time_work = MPI_Wtime();
		int index, i;
		elements_per_process = n / np;

		// check if more than 1 processes are run
		if (np > 1)
		{
			// distributes the portion of array
			// to child processes to calculate
			// their partial sums
			for (i = 1; i < np - 1; i++)
			{
				index = i * elements_per_process;

				MPI_Send(&elements_per_process, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
				MPI_Send(&index, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
			}

			// last process adds remaining elements
			index = i * elements_per_process;
			int elements_left = n - index;

			MPI_Send(&elements_left, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
			MPI_Send(&index, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
		}

		// master process add its own sub array
		for (i = 0; i < elements_per_process; i++)
			a[i]=1/factorial(i);

		// collects subarrays from other processes
		
		for (i = 1; i < np; i++)
		{
			MPI_Recv(&a[i * elements_per_process], elements_per_process, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
			int sender = status.MPI_SOURCE;
		}

		// prints the final sum of array
		std::cout << std::endl;
		std::cout << " The Program is RUN on " << np << " CPU" << std::endl;
	}

	// slave processes
	else
	{
		int index2;
		MPI_Recv(&n_elements_recieved, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);

		MPI_Recv(&index2, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);

		// calculates its subarr
		double subarr[n];
		for (int i = 0; i < n_elements_recieved; i++)
			subarr[i] = 1/factorial(i+index2);

		// sends the partial sum to the root process
		MPI_Send(&subarr[0], n_elements_recieved, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if (pid == 0)
	{
		int index, i;
		elements_per_process = n / np;

		// check if more than 1 processes are run
		if (np > 1)
		{
			// distributes the portion of array
			// to child processes to calculate
			// their partial sums
			for (i = 1; i < np - 1; i++)
			{
				index = i * elements_per_process;

				MPI_Send(&elements_per_process, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
				MPI_Send(&a[index], elements_per_process, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
			}

			// last process adds remaining elements
			index = i * elements_per_process;
			int elements_left = n - index;

			MPI_Send(&elements_left, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
			MPI_Send(&a[index], elements_left, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
		}

		// master process add its own sub array
		double sum = 0.0;
		for (i = 0; i < elements_per_process; i++)
			sum += a[i];
		// collects partial sums from other processes
		double tmp=0.0;
		for (i = 1; i < np; i++)
		{
			MPI_Recv(&tmp, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
			int sender = status.MPI_SOURCE;
			sum += tmp;
		}

		// prints the final sum of array
		std::cout << "exp = " << sum << std::endl;
		std::cout << " The Program is RUN on " << np << " CPU" << std::endl;
		time_work = MPI_Wtime() - time_work;
		std::cout << "Time     = " << time_work << std::endl;
	}

	// slave processes
	else
	{
		MPI_Recv(&n_elements_recieved, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);

		// stores the received array segment
		// in local array a2
		MPI_Recv(&a2, n_elements_recieved, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);

		// calculates its partial sum
		double partial_sum = 0.0;
		for (int i = 0; i < n_elements_recieved; i++)
			partial_sum += a2[i];

		// sends the partial sum to the root process
		MPI_Send(&partial_sum, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}

	// cleans up all MPI state before exit of process
	MPI_Finalize();


}




#include <stdio.h>
#include <stdlib.h>
#include <random>
#include <fstream>
#include <string>
#include <mpi.h>


//checks if a given point is inside the circle R1
bool inside_circle_1(double x, double y, const double R1){
	return x*x+y*y <= R1*R1;
}
//checks if a given point is inside the circle R2
bool inside_circle_2(double x, double y, const double x2, const double R2){
	return (x-x2)*(x-x2) + y*y <= R2*R2;
}
/*
for a given number of samples to throw this function returns the number of sambples (pts_inside) that
land inside the intersection of the two circles
*/
unsigned overlapMC(const double x2, const double R1, const double R2, size_t n, int rank, int procs)
{
	unsigned pts_inside = 0;

	const double x_min = std::min(-1*R1, x2-R2);
	const double x_max = std::max(R1, x2+R2);
	const double y_max = std::max(R1,R2);
	const double y_min = -1*y_max;

	std::default_random_engine g(42 + rank);  // random generator with seed 42
	std::uniform_real_distribution<double> u_x(x_min, x_max);
	std::uniform_real_distribution<double> u_y(y_min, y_max);

	// TODO_b: split the amount of work as equally as possible for each process.
	size_t n_local_size = (n+procs-1)/procs;
	size_t n_start = n_local_size * rank;
	size_t n_end   = std::min(n_local_size*(rank+1), n);

	for (size_t i = n_start; i < n_end; ++i)
	{
		// TODO_a: implement the MC integration part here!
		double x = u_x(g);
		double y = u_y(g);
		if(inside_circle_1(x,y,R1) && inside_circle_2(x,y,x2,R2))
			++pts_inside;
	}

	return  pts_inside;
}




int main(int argc, char *argv[])
{
	// TODO_b: Start-up the MPI environment and determine this process' rank ID as
    // well as the total number of processes (=ranks) involved in the
    // communicator
    int rank, procs = -1;
    rank = 0;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &procs);


	const double R1 = 5.0;		// Radius of first circle
	const double R2 = 10.0;		// Radius of second circle
	const double x2 = 12.0;		// x2 coordinate of second circle center

	// TODO_a: calculate the rectangle area for which you uniformly sample x & y
	const double x_min = std::min(-R1, x2-R2);
	const double x_max = std::max(R1, x2+R2);
	const double y_max = std::max(R1,R2);
	const double area_rectangle = 2*y_max * (std::abs(x_min) + x_max);

	size_t n = 1e9 + 1;// default number of MC samples

	double ref = 17.0097776; // reference solution

	double t0 = MPI_Wtime(); 

	unsigned local_sum = overlapMC(x2, R1, R2, n, rank, procs);

	unsigned global_sum = (procs == -1) ? local_sum : 0;

	// TODO_b: Sum up and average all the local_sums to the master ranks global_sum using MPI_Reduce
	MPI_Reduce(&local_sum, &global_sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	double area = global_sum / double(n) * area_rectangle;

	double t1 = MPI_Wtime();

	if(rank == 0){
		double error = std::abs(area - ref);
		if(error > 1e-2){
			printf("ERROR: you should get pi, but got: %.10f\n", area);
		}
		else{
			printf("result:  %.10f\nref: %.10f\ntime: %.10f\n", area, ref, t1 - t0);

			std::string file_name = "out/";
			file_name += std::to_string(procs);
			file_name += ".txt";
			//output time in a file
			std::ofstream myfile;
			myfile.open (file_name);
			myfile << procs << " " << t1 - t0 << "\n";
			myfile.close();
		}
	}

	// TODO_b: Shutdown the MPI environment
	MPI_Finalize();

	return 0;
}

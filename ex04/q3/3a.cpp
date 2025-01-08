#include <mpi.h>
#include <ostream>
#include <iostream>
#include <fstream>


int main(int argc, char** argv)
{

    MPI_Init(&argc,&argv);
    const int N = 8;
    double* result = new double[N];
    int rank, size;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &size );

    //assume N is divisible by size
    if (N%size !=0)
    {
        std::cerr << "N = " << N << " must be divisible by the number of ranks.\n";
        MPI_Abort(MPI_COMM_WORLD,1);
    }

    //size of each local array
    int batch_size = N/size;
    double* local_result = new double[batch_size];


    MPI_Scatter(result, batch_size, MPI_DOUBLE, local_result, batch_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // do a very computationally expensive calculation
    //...
    for(int i=0; i<batch_size; i++)
    {
        local_result[i] = rank;
    }
    // write the result to a file
    MPI_Gather(local_result, batch_size, MPI_DOUBLE, result, batch_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if(0 == rank)
    {
        std::ofstream file("result.txt");
        for(int i=0; i <= N; ++i)
            file << result[i] << std::endl;
    }


    delete[] result;
    MPI_Finalize();

    return 0;
}


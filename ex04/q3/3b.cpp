int rank, size, bufsize;
MPI_Init(&argc, &argv);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &size);

// only 2 ranks: 0, 1
double important_value;
// obtain the important value
//...
// exchange the value

MPI_Pack_size(1, MPI_INT, MPI_COMM_WORLD, &bufsize);
bufsize += MPI_BSEND_OVERHEAD;
char *buffer = new char[bufsize];
MPI_Buffer_attach(buffer, bufsize);

if(rank == 0)
    MPI_Bsend(&important_value, 1, MPI_DOUBLE, 1, 123, MPI_COMM_WORLD);
else
    MPI_Bsend(&important_value, 1, MPI_DOUBLE, 0, 123, ,MPI_COMM_WORLD);

MPI_Recv(&important_value , 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,  MPI_STATUS_IGNORE);

MPI_Buffer_detach(buffer, &bufsize);
delete[] buffer;

 // do other work

MPI_Finalize();



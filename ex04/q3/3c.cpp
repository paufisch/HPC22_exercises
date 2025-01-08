

MPI_Init(&argc, &argv);

int rank, size;
MPI_Init(&argc, &argv);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &size);

int bval;
bval=rank;
MPI_Bcast(&bval, 1, MPI_INT, 0, MPI_COMM_WORLD);

MPI_Status stat;
MPI_Recv(&bval, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);

cout << "[" << rank << "]" << endl;

MPI_Finalize();
return 0;


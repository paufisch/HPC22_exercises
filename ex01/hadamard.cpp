#include<chrono>
#include<stdio.h>
#include<iostream>
#include<cmath>

void hadamard_row_wise(double* A, double* B, double* C, int N){
    for (size_t i = 0; i < N*N; i++){
        C[i] = 0;

    // TODO: Hadamard product row-wise
        int j = i / N;
        int k = i % N;
        C[i] = A[j * N + k] * B[j * N + k];
    }
}

void hadamard_col_wise(double* A, double* B, double* C, int N){
    for (size_t i = 0; i < N*N; i++){
        C[i] = 0;


    // TODO: Hadamard product col-wise
        int j = i / N;
        int k = i % N;
        C[i] = A[k * N + j] * B[k * N + j];
    }
}

int main(int argc, char *argv[]){
    if(argc < 2){
        printf("Error: give me matrix size\n");
        return 1;
    }

    const int N = atoi(argv[1]);

    printf("\n================================\n");
    printf("N = %u\n", N);


    // TODO: Randomly initialize the matrices A,B
    double* A = (double*) calloc(N*N, sizeof(double));
    double* B = (double*) calloc(N*N, sizeof(double));
    double* C = (double*) calloc(N*N, sizeof(double)); //calloc = malloc but set to zeros
    

    const int reps = 1e3;

    printf("start running Hadamard C = A°B row-wise benchmark...\n");
    // TODO: benchmark your hadamard_row_wise function using std::chrono and with reps=1e3. 
    // Save the time in the variable "time_row". Measure in seconds!
    auto start = std::chrono::system_clock::now();
    hadamard_row_wise(A,B,C,N);
    auto end = std::chrono::system_clock::now();
    auto elapsed = end - start;
    double time_row = elapsed.count();

    printf("Hadamard row-wise time %9.4fs\n", time_row);

    double* C_temp = (double*) calloc(N*N, sizeof(double));
    for(int i = 0; i < N*N; i++){
        C_temp[i] = C[i];
    }

    printf("start running Hadamard C = A°B col-wise benchmark...\n");
    // TODO: benchmark your hadamard_col_wise function using std::chrono and with reps=1e3.
    // Save the time in the variable "time_col". Measure in seconds!¨
    auto start2 = std::chrono::system_clock::now();
    hadamard_col_wise(A,B,C,N);
    auto end2 = std::chrono::system_clock::now();
    auto elapsed2 = end2 - start2;
    double time_col = elapsed2.count();


    printf("Hadamard row-wise time %9.4fs\n", time_col);

    double error = 0.0;
    for(int i = 0; i < N*N; i++){
        error += std::abs(C_temp[i] - C[i]);
    }

    if(error > 1e-10){
        printf("======================================================\n");
        printf("ERROR: row-wise and col-wise are not similar! error = %.40f\n", error);
        printf("======================================================\n");
        return 1;
    }

    printf("T(col-wise)/T(row-wise) = %9.4f\n", time_col/time_row);

    //free allocated memory
    free(A);
    free(B);
    free(C);

    return 0;
}

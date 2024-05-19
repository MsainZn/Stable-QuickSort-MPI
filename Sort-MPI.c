#include "SortUtilities.h"

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int nprocs, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int n = 1007;
    double *array = NULL;
    double sort_time;
    int i;
    if (rank == 0) {
        array = (double *)malloc(n * sizeof(double));
        for (i = 0; i < n; i++) {
            array[i] = 0.0 + (double)rand() / RAND_MAX * (1000.0 - 0.0);
        }
    }

    QuickSortMPI(&array, n, nprocs, rank, &sort_time, 0);

    if (rank == 0) {
        printf("\nElapsed time: %10.6lf s\n", sort_time);
        PrintArray(array, n);
        free(array);
    }

    MPI_Finalize();
    return 0;
}

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <math.h>
#include <float.h>
#include <string.h>

// ProtoType
void MergeSort(double array[], int count);
void MergeSortRecursion (double array[], int left_index, int right_index);
void MergeSortReorder(double array[], int left_index, int mid_index, int right_index);

void QuickSort(double array[], int count);
void QuickSortRecursive(double array[], int left_index, int right_index);
int QuickSortMPI(double** array, int size, int nprocs, int rank, double* sort_time, int use_tree);
void SwapElements(double* a, double* b);
int QuickSortPartitionerBase(double array[], int left_index, int right_index);
void MergeTwoSortedArrays(double *arrA, int arrA_count, double *arrB, int arrB_count, double *merged_arr);
int IsSortingDoublesCorrect(double array[], int arr_count, int* j);
void PrintArray(double vector[], int length);
void RunTest01(void);

// Implementation
void MergeSort(double array[], int count){
    MergeSortRecursion(array, 0, count-1);
}

void MergeSortRecursion (double array[], int left_index, int right_index){
    if (left_index<right_index){
        int mid_index = left_index + (right_index - left_index)/2;
        MergeSortRecursion(array, left_index, mid_index);
        MergeSortRecursion(array, mid_index+1, right_index);
        MergeSortReorder(array, left_index, mid_index, right_index);
    }
}

void MergeSortReorder (double array[], int left_index, int mid_index, int right_index){
    int lcount = mid_index - left_index + 1;
    int rcount = right_index - mid_index;
    
    double* ltemp = (double*) malloc(lcount * sizeof(double));
    double* rtemp = (double*) malloc(rcount * sizeof(double));
	
    int i, j, k;

    for (i = 0; i < lcount; i++){
        ltemp[i] = array[left_index+i];
    }
    
    for (j = 0; j < rcount; j++){
        rtemp[j] = array[mid_index+1+j];
    }

     k = left_index; i = 0; j = 0;
     while (i < lcount && j < rcount)
     {
        if (ltemp[i] > rtemp[j])
        {
            array[k] = rtemp[j];
            j++;
        }
        else{
            array[k] = ltemp[i];
            i++;
        }
        k++;
     }

     while (i<lcount)
     {
        array[k] = ltemp[i];
        k++;
        i++;
     }
     
     while (j<rcount)
     {
        array[k] = rtemp[j];
        k++;
        j++;
     }

    free(ltemp); ltemp = NULL;
    free(rtemp); rtemp = NULL;

}

void QuickSort(double array[], int count){
    QuickSortRecursive(array, 0, count-1);
}

void QuickSortRecursive(double array[], int left_index, int right_index){
    
    if (left_index < right_index){
        // int pivot_index = QuickPartitioner(array, left_index, right_index);
        int pivot_index = QuickSortPartitionerBase(array, left_index, right_index);
        QuickSortRecursive(array, left_index, pivot_index-1);
        QuickSortRecursive(array, pivot_index+1, right_index);
    }
}

void SwapElements(double* a, double* b){
    double temp = *a;
    *a = *b;
    *b = temp;
}

int QuickSortPartitionerBase(double array[], int left_index, int right_index){

    // Select Element On Random and then push it at array right-end
    int pivot_index = left_index + rand()%(right_index-left_index+1);
    SwapElements(&array[pivot_index], &array[right_index]);
    pivot_index = right_index;
    double pivot_element = array[pivot_index];
    // Sort based on pivot
    int i, j;
    for (i = left_index, j = left_index; i < right_index; i++){
        if (array[i] < pivot_element){
            SwapElements(&array[j], &array[i]);
            j++;
        }
    }
    SwapElements(&array[pivot_index], &array[j]);
    
    return j;
}

void MergeTwoSortedArrays(double *arrA, int arrA_count, double *arrB, int arrB_count, double *merged_arr)
{
    int i = 0, j = 0, index = 0;

    while (i < arrA_count && j < arrB_count)
    {
        if (arrA[i] <= arrB[j])
        {
            merged_arr[index] = arrA[i];
            i++;
        }
        else
        {
            merged_arr[index] = arrB[j];
            j++;
        }
        index++;
    }

    while (i < arrA_count)
    {
        merged_arr[index] = arrA[i];
        i++;
        index++;
    }

    while (j < arrB_count)
    {
        merged_arr[index] = arrB[j];
        j++;
        index++;
    }
}

int IsSortingDoublesCorrect(double array[], int arr_count, int* j)
{
    int i;
    *j = -1;
    for (i = 0; i < arr_count - 1; i++)
    {
        if (array[i] > array[i + 1])
        {
            *j = i;
            return 0;
        }
    }
    return 1;
}

void PrintArray(double array [], int arr_count)
{
    int i;
    for (i = 0; i < arr_count; i++){
        printf("arr[%d]: %15.10lf\n", i,array[i]);
    }

}

int QuickSortMPI(double** array, int size, int nprocs, int rank, double* sort_time, int use_tree) {

    int i;
    double *data_sub = NULL;
    int *send_counts = NULL, *send_displacements = NULL;
    double time_init, time_end;

    if (rank == 0) {
        send_counts = (int *)malloc(nprocs * sizeof(int));
        send_displacements = (int *)malloc(nprocs * sizeof(int));

        // Determine the send counts and displacements
        int remainder = size % nprocs;
        int sum = 0;
        for (i = 0; i < nprocs; i++) {
            send_counts[i] = size / nprocs + (i < remainder ? 1 : 0);
            send_displacements[i] = sum;
            sum += send_counts[i];
        }

        time_init = MPI_Wtime();
    }

    // Broadcast the size of each chunk to all processes
    int elements_per_proc;
    MPI_Scatter(send_counts, 1, MPI_INT, &elements_per_proc, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Allocate space for each process's sub-array
    data_sub = (double *)malloc(elements_per_proc * sizeof(double));

    // Scatter the data using MPI_Scatterv
    MPI_Scatterv(*array, send_counts, send_displacements, MPI_DOUBLE, data_sub,
                 elements_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Quick sort in serial
    QuickSort(data_sub, elements_per_proc);
    
    // Merge Algorithm Tree-based
    if (use_tree){
        // Tree-based merge
        if (rank == 0) {
            free(send_counts);  send_counts = NULL;
            free(send_displacements); send_displacements = NULL;
            free(*array); *array = NULL;
        }
        int step = 1;
        MPI_Barrier(MPI_COMM_WORLD);
        while (step < nprocs) {
            if (rank % (2 * step) == 0) {
                if (rank + step < nprocs) {
                    int recv_count = -1;
                    MPI_Recv(&recv_count, 1, MPI_INT, rank + step, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    double* buffer_recv = (double*) malloc(recv_count * sizeof(double));
                    MPI_Recv(buffer_recv, recv_count, MPI_DOUBLE, rank + step, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    double* merger = (double*) malloc((recv_count + elements_per_proc) * sizeof(double));
                    MergeTwoSortedArrays(buffer_recv, recv_count, data_sub, elements_per_proc, merger);
                    elements_per_proc += recv_count;
                    free(buffer_recv); buffer_recv = NULL;
                    free(data_sub);    data_sub = NULL;
                    data_sub = merger;
                }
            } else {
                int rank_recver = rank - step;
                MPI_Send(&elements_per_proc, 1, MPI_INT, rank_recver, 0, MPI_COMM_WORLD);
                MPI_Send(data_sub, elements_per_proc, MPI_DOUBLE, rank_recver, 0, MPI_COMM_WORLD);
                break;
            }
            step *= 2;
        }
        if (rank == 0) {
            // Evaluate Validity
            int err_idx;
            time_end = MPI_Wtime() - time_init;
            if (!IsSortingDoublesCorrect(data_sub, size, &err_idx)) {
                return -1;
            }
            // Re-asemble
            *array = data_sub; 
            *sort_time = time_end;
        }
    // Merge Algorithm linear-based (ONLY PROCESS 0)
    } else {
        MPI_Gatherv(data_sub, elements_per_proc, MPI_DOUBLE, *array, send_counts, send_displacements, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        // Perform linear merge at root process directly into data
        if (rank == 0) {
            int *current_indices = (int *)calloc(nprocs, sizeof(int));
            int sorted_index = 0;

            while (sorted_index < size) {
                double min_value = INFINITY;
                int min_index = -1;

                for (i = 0; i < nprocs; i++) {
                    int idx = send_displacements[i] + current_indices[i];
                    if (current_indices[i] < send_counts[i] && (*array)[idx] < min_value) {
                        min_value = (*array)[idx];
                        min_index = i;
                    }
                }
                (*array)[sorted_index++] = min_value;
                current_indices[min_index]++;
            }
            free(current_indices);
            free(send_counts);
            free(send_displacements);
        }
        free(data_sub);

        if (rank == 0) {
            // Evaluate Validity
            int err_idx;
            time_end = MPI_Wtime() - time_init;
            if (!IsSortingDoublesCorrect(*array, size, &err_idx)) {
                return -1;
            }
            // Re-asemble
            *sort_time = time_end;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    return 0;
}

void RunTest01(void)
{
    double array[] = {10.6, 7.7, 8., 94.1, 1.46, 56.2};
    int array_size = sizeof(array)/sizeof(array[0]);
    int i;
    printf("Given array is \n");
    for (i = 0; i < array_size; i++)
        printf("%lf ", array[i]);
    printf("\n");

    QuickSort(array, array_size);
    // MergeSort(array, array_size);

    printf("Sorted array is \n");
    for (int i = 0; i < array_size; i++)
        printf("%lf ", array[i]);
    printf("\n");
}

#include <stdio.h>
#include <time.h>
#include <stdlib.h>


/**
 * Function for swapping elements
 * @param a Element to be swapped
 * @param b Element to be swapped
 */
void swapMpi(double* a, double* b){
    double temp = *a;
    *a = *b;
    *b = temp;
}

/**
 * Finds the right position of the pivot
 * @param arr
 * @param left lowest element of the array
 * @param right highest element of the array
 * @return pivot's index
 */
int partitionMpi(double arr[], int left, int right){
    double pivot = arr[right];

    int i = left;
    int j = left;

    while (j < right){
        if(arr[j] < pivot){
            swapMpi(&arr[i], &arr[j]);
            i++;
        }
        j++;
    }

    swapMpi(&arr[i], &arr[right]);

    return i;
}


/**
 * Picks a random pivot element between left and right and partition arr[l...r] around the randomly picked element using partition()
 * @param arr
 * @param left
 * @param right
 * @return the right position of randomly picked pivot
 */
int randomPartitionMpi(double arr[], int left, int right){
    srand(time(NULL));
    int length = right - left + 1;
    int pivot = rand() % length;

    swapMpi(&arr[left + pivot], &arr[right]);

    return partitionMpi(arr, left, right);
}

/**
 * Finds median
 * @param arr
 * @param left
 * @param right
 * @param middle
 * @param a
 * @param b
 */
void MedianMpiUtil(double arr[], int left, int right, int middle, double *a, double *b){

    if (left <= right) {

        // Find the partition index
        int partitionIndex = randomPartitionMpi(arr, left, right);


        if (partitionIndex == middle) {
            *b = arr[partitionIndex];
            if (*a != -1)
                return;
        }


        if (partitionIndex == middle - 1) {
            *a = arr[partitionIndex];
            if (*b != -1)
                return;
        }

        // If partitionIndex >= k then
        // find the index in first half
        // of the arr[]
        if (partitionIndex >= middle){
            return MedianMpiUtil(arr, left, partitionIndex - 1, middle, a, b);
        }
            // If partitionIndex <= k then
            // find the index in second half
            // of the arr[]
        else{
            return MedianMpiUtil(arr, partitionIndex + 1, right, middle, a, b);
        }
    }

}

/**
 * Finds the median
 * @param arr
 * @param length
 * @return median
 */
double findMpiMedian(double *arr, int length){
    double result;

    double a = -1.0, b = -1.0;

    // If n is odd
    if (length % 2 == 1) {
        MedianMpiUtil(arr, 0, length - 1, length / 2, &a, &b);
        result = b;
    }
        // If n is even
    else {
        MedianMpiUtil(arr, 0, length - 1, length / 2, &a, &b);
        result = (a + b) / 2;
    }

    // Print the Median of arr[]
//    printf("The median is: %f\n", result);

    return result;

}


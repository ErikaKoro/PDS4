#include "mpi_tree.h"

#include "quick_select.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>


/**
 * Function for swapping elements
 * @param a Element to be swapped
 * @param b Element to be swapped
 */
void swapPerProc(double* a, double* b){
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
int64_t partitionPerProc(double **points, double arr[], int64_t left, int64_t right){
    double pivot = arr[right];

    int64_t i = left;
    int64_t j = left;

    while (j < right){
        if(arr[j] < pivot){
            swapPerProc(&arr[i], &arr[j]);
            double *tmp = points[i];
            points[i] = points[j];
            points[j] = tmp;
            i++;
        }
        j++;
    }

    swapPerProc(&arr[i], &arr[right]);
    double *temp = points[i];
    points[i] = points[right];
    points[right] = temp;

    return i;
}


/**
 * Picks a random pivot element between left and right and partition arr[l...r] around the randomly picked element using partition()
 * @param arr
 * @param left
 * @param right
 * @return the right position of randomly picked pivot
 */
int64_t randomPartitionPerProc(double **points, double arr[], int64_t left, int64_t right){
    srand(time(NULL));
    int length = right - left + 1;
    int pivot = rand() % length;

    swapPerProc(&arr[left + pivot], &arr[right]);
    double *temp = points[left + pivot];
    points[left + pivot] = points[right];
    points[right] = temp;

    return partitionPerProc(points, arr, left, right);
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
void MedianUtilPerProc(double **points, double arr[], int64_t left, int64_t right, int64_t middle, double *a, double *b){

    if (left <= right) {

        // Find the partition index
        int64_t partitionIndex = randomPartitionPerProc(points, arr, left, right);


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
            return MedianUtilPerProc(points, arr, left, partitionIndex - 1, middle, a, b);
        }
            // If partitionIndex <= k then
            // find the index in second half
            // of the arr[]
        else{
            return MedianUtilPerProc(points, arr, partitionIndex + 1, right, middle, a, b);
        }
    }

}


/**
 * Finds the median
 * @param arr
 * @param length
 * @return median
 */
double findMedianPerProc(double **points, double *arr, int64_t length){
    double result;

    double a = -1.0, b = -1.0;

    // If n is odd
    if (length % 2 == 1) {
        MedianUtilPerProc(points, arr, 0, length - 1, length / 2, &a, &b);
        result = b;
    }
        // If n is even
    else {
        MedianUtilPerProc(points, arr, 0, length - 1, length / 2, &a, &b);
        result = (a + b) / 2;
    }

    // Print the Median of arr[]
//    printf("The median is: %f\n", result);

    return result;

}

void findDistance(double *dist, double **points, int64_t dimension, const double *pivot, int64_t numberOfPoints){

    for (int i = 0; i < numberOfPoints; ++i) {
        for (int j = 0; j < dimension; ++j) {
            dist[i] += (points[i][j] - pivot[j]) * (points[i][j] - pivot[j]);
        }
       // printf("\nThe distance is %.10f\n", dist[i]);
    }
}

void buildVPTree(vptree *parentTree, double **points, double *distances, int64_t dimension, int64_t numberOfPoints) {
    // Condition that terminates the recursion
    if (numberOfPoints == 1)
        return;

    // Find the median distance, using "quick select" algorithm
    // In findMedian we sort the array of points according to the sort of distances.
    parentTree->median = findMedianPerProc(points + parentTree->start, distances + parentTree->start, numberOfPoints);

    // Initialize inner tree
    parentTree->inner = (vptree *) malloc(sizeof(vptree));      // FREEMEM
    parentTree->inner->inner = NULL;
    parentTree->inner->outer = NULL;
    parentTree->inner->start = parentTree->start;
    parentTree->inner->stop = numberOfPoints / 2 - 1 + parentTree->start;
//    parentTree->inner->vpPoint = parentTree->vpPoint;

    buildVPTree(
            parentTree->inner,
            points,
            distances,
            dimension,
            parentTree->inner->stop - parentTree->inner->start + 1
    );


    // Initialize outer tree
    parentTree->outer = (vptree *) malloc(sizeof(vptree));      // FREEMEM
    parentTree->outer->inner = NULL;
    parentTree->outer->outer = NULL;
    parentTree->outer->start = numberOfPoints / 2 + parentTree->start;
    parentTree->outer->stop = parentTree->stop;
//    parentTree->outer->vpPoint = parentTree->vpPoint;

    buildVPTree(
            parentTree->outer,
            points,
            distances,
            dimension,
            parentTree->outer->stop - parentTree->outer->start + 1
    );
}


void freeMpiMemory(vptree *tree) {

    if (tree->inner == NULL && tree->outer == NULL) {        // If the tree's leaves are empty
        free(tree);
        return;
    }

    if (tree->inner == NULL) {
        freeMpiMemory(tree->outer);
        tree->outer = NULL;
    }
    else if (tree->outer == NULL){
        freeMpiMemory(tree->inner);
        tree->inner = NULL;
    }
    else {
        freeMpiMemory(tree->inner);
        tree->inner = NULL;
        freeMpiMemory(tree->outer);
        tree->outer = NULL;
    }

    if(tree->inner == NULL && tree->outer == NULL) {        // If the tree's leaves are empty
        free(tree);
    }

}

/**
* This function calculates the k nearest neighbours of a point.
*
* @param k The number of the k nearest neighbours
* @param nearest The array that holds the k nearest points
* @param points The array with all the points but sorted
* @param numberOfPoints
*/
void knn_mpi_search(int k, double **nearest, double **points, int64_t numberOfPoints, int64_t dimension){
    // The array points is sorted according to the sort of the array distances. So, when it is
    // given as an argument the first k points of it are the nearest to the chosen point
    if(k < numberOfPoints) {
        printf("\nThe nearest are: \n");
        for (int i = 0; i < k; ++i) {
            for (int j = 0; j < dimension; ++j) {
                nearest[i][j] = points[i + 1][j];
                printf("%.1f ", nearest[i][j]);
            }
            printf("\n");
        }
    }
    else if(k == numberOfPoints || k > numberOfPoints){
        printf("\nThe nearest are: \n");
        for (int i = 0; i < k; ++i) {
            for (int j = 0; j < dimension; ++j) {
                if( i == numberOfPoints - 1)        // The points array has no more elements
                    break;
                nearest[i][j] = points[i + 1][j];
                printf("%.1f ", nearest[i][j]);
            }
            printf("\n");
        }
    }
}

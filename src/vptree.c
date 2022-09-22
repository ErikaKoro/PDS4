#include "vptree.h"
#include "quick_select.h"
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>


/**
 * Finds the distance of each element, per process, from the pivot
 *
 * @param dist array to hold the results
 * @param points points per process with their coordinates
 * @param dimension the dimension of our space (e.g 3D)
 * @param pivot array with the coordinates of the pivot
 * @param pointsPerProc
 */
void findDistance(double *dist, double **points, int64_t dimension, const double *pivot, int64_t numberOfPoints){

    for (int i = 0; i < numberOfPoints; ++i) {
        for (int j = 0; j < dimension; ++j) {
            dist[i] += (points[i][j] - pivot[j]) * (points[i][j] - pivot[j]);
        }
//        printf("The distance is %.10f\n", dist[i]);
    }
}


/**
 * Builds the vantage point tree according to the median distance from the vantage point
 *
 * @param parentTree The initial tree of the vantage point
 * @param points The array that holds all the points' coordinates
 * @param distances The array that holds the distances calculated once at the beggining of the program
 * @param dimension The points' dimension
 * @param numberOfPoints The number of points
 */
void buildVPTree(vptree *parentTree, double **points, double *distances, int64_t dimension, int64_t numberOfPoints){
    // Condition that terminates the recursion
    if(numberOfPoints == 1)
        return ;

    // Find the median distance, using "quick select" algorithm
    // In findMedian we sort the array of points according to the sort of distances.
    parentTree->median = findMedian(points + parentTree->start, distances + parentTree->start, numberOfPoints);

    // Initialize inner tree
    parentTree->inner = (vptree *)malloc(sizeof (vptree));      // FREEMEM
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
    parentTree->outer = (vptree *)malloc(sizeof (vptree));      // FREEMEM
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


/**
 * Function that is called in order to free the memory in a tree
 *
 * @param tree the memory to be de-allocated
 */
void freeMemory(vptree *tree) {

    if (tree->inner == NULL && tree->outer == NULL) {        // If the tree's leaves are empty
        free(tree);
        return;
    }

    if (tree->inner == NULL) {
        freeMemory(tree->outer);
        tree->outer = NULL;
    }
    else if (tree->outer == NULL){
        freeMemory(tree->inner);
        tree->inner = NULL;
    }
    else {
        freeMemory(tree->inner);
        tree->inner = NULL;
        freeMemory(tree->outer);
        tree->outer = NULL;
    }

    if(tree->inner == NULL && tree->outer == NULL) {        // If the tree's leaves are empty
        free(tree);
    }

}


int main(int argc, char **argv){

    if(argc != 2){
        printf("gimme args you stupid\n");
        return -1;
    }

    /* Read the file */
    int64_t numberOfPoints;
    int64_t dimension;
    double **holdThePoints;

    FILE *fh = fopen (argv[1], "rb");
    if (fh != NULL) {
        fread(&dimension, sizeof(int64_t), 1, fh);
        fread(&numberOfPoints, sizeof(int64_t), 1, fh);

        holdThePoints = (double **)malloc(sizeof(double *) * numberOfPoints);       // FREEMEM

        for(int i = 0; i < numberOfPoints; i++){
            holdThePoints[i] = (double *) malloc(sizeof (double) * dimension);      // FREEMEM
        }

        for (int i = 0; i < numberOfPoints; i++) {
            fread(*(holdThePoints + i), sizeof(double), dimension, fh);
        }

        fclose (fh);
    }

    printf("The array with the points is:\n\n");
    for (int i = 0; i < numberOfPoints; ++i) {
        for (int j = 0; j < dimension; ++j) {
            printf("%.2f ", holdThePoints[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    vptree initial;
    initial.inner = NULL;
    initial.outer = NULL;
    initial.start = 0;
    initial.stop = numberOfPoints - 1;
    initial.vpIndex = numberOfPoints - 1;

    // Allocate memory for the vantage point and fill it with the last point's coordinates(random choice)
    initial.vpPoint = (double *)malloc(dimension * sizeof(double));     // FREEMEM
    for (int i = 0; i < dimension; ++i) {
        initial.vpPoint[i] = holdThePoints[numberOfPoints - 1][i];
    }

    // Allocate the array that will hold the distances from the vantage point to the others
    double *distances = (double *)calloc(numberOfPoints, sizeof(double ));     // FREEMEM

    // Calculate the distances from the chosen pivot
    findDistance(distances, holdThePoints, dimension, initial.vpPoint, numberOfPoints);

    buildVPTree(&initial, holdThePoints, distances, dimension, numberOfPoints);

    for (int i = 0; i < numberOfPoints; ++i) {
        printf("The distance is %.10f\n", distances[i]);
    }

    // Free memory SO AS NOT TO HAVE MEMORY LEAKS(I DON'T DO SUCH THINGS, A FRIEND TOLD ME)
    for (int i = 0; i < numberOfPoints; ++i) {
        free(holdThePoints[i]);
    }
    free(holdThePoints);
    free(distances);
    free(initial.vpPoint);
    freeMemory(initial.inner);
    freeMemory(initial.outer);

}
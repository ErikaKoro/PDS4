#include "sequential_vptree.h"
#include "quick_select.h"
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <sys/time.h>

typedef struct knn{
    double **nearest;
}KNN;

/**
 * Calculates the elapse time
 *
 * @param begin the starting timestamp
 * @param end the ending timestamp
 * @return elapsed time in seconds
 */
double measureTime(struct timeval begin, struct timeval end) {
    long seconds;
    long microseconds;
    double elapsed;

    seconds = end.tv_sec - begin.tv_sec;
    microseconds = end.tv_usec - begin.tv_usec;
    elapsed = seconds + microseconds * 1e-6;

    return elapsed;
}


/**
 * Finds the distance of each element, per process, from the pivot
 *
 * @param dist array to hold the results
 * @param points points per process with their coordinates
 * @param dimension the dimension of our space (e.g 3D)
 * @param pivot array with the coordinates of the pivot
 * @param pointsPerProc
 */
double *findDistance(double *dist, double **points, int64_t dimension, const double *pivot, int64_t numberOfPoints){

    for (int i = 0; i < numberOfPoints; ++i) {
        for (int j = 0; j < dimension; ++j) {
            dist[i] += (points[i][j] - pivot[j]) * (points[i][j] - pivot[j]);
        }
        // printf("The distance is %.10f\n", dist[i]);
    }
    return dist;
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
        free(tree);     // free the tree
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

/**
 * This function testes whether the sort of the distances is right after the call of the buildVPTree function
 * @param distances
 * @param numberOfPoints
 */
void testFunction(double const *distances, int64_t numberOfPoints){
    for (int i = 1; i < numberOfPoints; ++i) {
        if(distances[i - 1] > distances[i]){
            printf("\nTest Failed\n");
            return;
        }
    }
    // printf("Test passed\n\n");
}

/**
 * This function copies one array to another
 * @param points The array that need to be copied
 * @param copied The array where the other arrays is copied to
 * @param numberOfPoints The number of array's rows
 * @param dimension The number of array's cols
 */
void deepCopyArray(double **points, double **copied, int64_t numberOfPoints, int64_t dimension){

    for (int i = 0; i < numberOfPoints; ++i) {
        for (int j = 0; j < dimension; ++j) {
            copied[i][j] = points[i][j];
        }
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
void knn_search(int k, double **nearest, double **points, int64_t numberOfPoints, int64_t dimension){
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

/**
 * This function calculates for each point its distances from the others and sorts the distances array
 * according to the median distance and finds its vantage tree
 * @param dimension
 * @param numberOfPoints
 * @param k
 * @param points
 */
void calculateKNN(int64_t dimension, int64_t numberOfPoints, int k, double **points){

    KNN *total = (KNN *)malloc(sizeof (KNN) * numberOfPoints);      // Create an object for each point      // FREEMEM
    double **copied = (double **)malloc(sizeof (double *) * numberOfPoints);        // FREEMEM
    for (int i = 0; i < numberOfPoints; ++i) {
        copied[i] = (double *) malloc(sizeof (double ) * dimension);        // FREEMEM
    }

    // Copy all the points here(so as we can choose the correct vpPoint after the call of VPTree
    deepCopyArray(points, copied, numberOfPoints, dimension);

    // Now, for each point allocate memory for its array with the neighbours
    for (int i = 0; i < numberOfPoints; ++i) {      // for all points, their k neighbours need to be calculated

        vptree *initial = (vptree *) malloc(sizeof (vptree) * numberOfPoints);      // FREEMEM
        initial->inner = NULL;
        initial->outer = NULL;
        initial->start = 0;
        initial->stop = numberOfPoints - 1;

        // ALLOCATE
        total[i].nearest = (double **) malloc(sizeof(double) * k);      // for each point, each nearest array has k points
        for (int j = 0; j < k; ++j) {
            total[i].nearest[j] = (double *)malloc(sizeof (double ) * dimension);
        }

        // Allocate memory for the vantage point and fill it with the last point's coordinates(random choice)
        initial->vpPoint = (double *)malloc(dimension * sizeof(double));     // FREEMEM
        for (int j = 0; j < dimension; ++j) {
            initial->vpPoint[j] = copied[i][j];
        }

        // Allocate the array that will hold the distances from the vantage point to the others
        double *distances = (double *)calloc(numberOfPoints, sizeof(double ));     // FREEMEM

        // Calculate the distances from the chosen pivot
        distances = findDistance(distances, points, dimension, initial->vpPoint, numberOfPoints);


        buildVPTree(initial, points, distances, dimension, numberOfPoints);


        testFunction(distances, numberOfPoints);

//        printf("The sorted distances are: \n");
//        for (int j = 0; j < numberOfPoints; ++j) {
//            printf("%.1f ", distances[j]);
//        }
//        printf("\n\n");
//
//        printf("The sorted points are: \n");
//        for (int j = 0; j < numberOfPoints; ++j) {
//            for (int l = 0; l < dimension; ++l) {
//                printf("%.1f ", points[j][l]);
//            }
//            printf("\n");
//
//        }

        knn_search(k, total[i].nearest, points, numberOfPoints, dimension);

        free(distances);
        free(initial->vpPoint);
        freeMemory(initial);
    }


    // FREE ME
    free(total);
    for (int i = 0; i < numberOfPoints; ++i) {
        free(copied[i]);
    }
    free(copied);
}


int main(int argc, char **argv){

    if(argc != 3){
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
    printf("The dimension is %ld ", dimension);

//    printf("\n\nThe array with the points is:\n\n");
//    for (int i = 0; i < numberOfPoints; ++i) {
//        for (int j = 0; j < dimension; ++j) {
//            printf("%.2f ", holdThePoints[i][j]);
//        }
//        printf("\n");
//    }
//    printf("\n");

    // Vars needed for execution time measurement
    struct timeval begin, end;
    int k = atoi(argv[2]);

    gettimeofday(&begin, 0);
    calculateKNN(dimension, numberOfPoints, k, holdThePoints);
    gettimeofday(&end, 0);
    printf("Time for tree construction: %.5f seconds.\n", measureTime(begin, end));


     // Free memory SO AS NOT TO HAVE MEMORY LEAKS(I DON'T DO SUCH THINGS, A FRIEND TOLD ME)
    for (int i = 0; i < numberOfPoints; ++i) {
        free(holdThePoints[i]);
    }
    free(holdThePoints);
}
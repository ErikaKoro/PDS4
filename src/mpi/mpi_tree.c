#include "mpi_tree.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "../quick_select.h"


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
    parentTree->median = findMedian(points + parentTree->start, distances + parentTree->start, numberOfPoints);

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

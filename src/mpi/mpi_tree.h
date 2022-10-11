#ifndef PDS4_MPI_TREE_H
#define PDS4_MPI_TREE_H

#include <inttypes.h>

typedef struct vptree{
    int64_t start;      // The index of the tree's start point
    int64_t stop;
    double *vpPoint;   // vantage point
    double median;     // the euclidean median distance from the vantage point
    struct vptree *inner;   // vantage point subtrees
    struct vptree *outer;
}vptree;

void findDistance(double *dist, double **points, int64_t dimension, const double *pivot, int64_t numberOfPoints);

void buildVPTree(vptree *parentTree, double **points, double *distances, int64_t dimension, int64_t numberOfPoints);

void knn_mpi_search(int k, double **nearest, double **points, int64_t numberOfPoints, int64_t dimension);

void freeMpiMemory(vptree *tree);

#endif //PDS4_MPI_TREE_H

#include <inttypes.h>
#ifndef VPTREE_H
#define VPTREE_H

typedef struct vptree{
    int64_t start;      // The index of the tree's start point
    int64_t stop;
    double *vpPoint;   // vantage point
    double median;     // the euclidean median distance from the vantage point
    struct vptree *inner;   // vantage point subtrees
    struct vptree *outer;
}vptree;

void buildVPTree(vptree *parentTree, double **points, double *distances, int64_t dimension, int64_t numberOfPoints);

void testFunction(double const *distances, int64_t numberOfPoints);

#endif

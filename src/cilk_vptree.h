#include <inttypes.h>
#include <cilk/cilk.h>
#ifndef PDS4_MPI_VPTREE_H
#define PDS4_MPI_VPTREE_H

typedef struct vptree{
    int64_t start;      // The index of the tree's start point
    int64_t stop;
    double *vpPoint;   // vantage point
    double median;     // the euclidean median distance from the vantage point
    int64_t vpIndex;   // the index of the vantage point tree in the original set
    struct vptree *inner;   // vantage point subtrees
    struct vptree *outer;
}vptree;

#endif //PDS4_MPI_VPTREE_H

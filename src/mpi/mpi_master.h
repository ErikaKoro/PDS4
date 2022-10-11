#ifndef PDS4_MPI_MASTER_H
#define PDS4_MPI_MASTER_H
#include <inttypes.h>

typedef struct mpiKNN{
    double **nearest;
}KNN;

void calculateMasterKNN(KNN perPivot,
                       const double *pivot,
                       int rank,
                       int size,
                       int64_t pointsPerProc,
                       double **holdPoints,
                       int64_t dimension,
                       int k_neighbours
);

#endif //PDS4_MPI_MASTER_H

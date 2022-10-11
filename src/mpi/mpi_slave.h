#ifndef PDS4_MPI_SLAVE_H
#define PDS4_MPI_SLAVE_H
#include <inttypes.h>

void calculateSlaveKNN(const double *pivot, int rank, int size, int64_t pointsPerProc, double **holdPoints, int64_t dimension, int k_neighbours);

#endif //PDS4_MPI_SLAVE_H

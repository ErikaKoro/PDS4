#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi_exchanges.h"
#include "mpi_quick_select.h"
#include "distribute_by_median.h"


/**
 * Finds the distance of each element, per process, from the pivot
 *
 * @param dist array to hold the results
 * @param points points per process with their coordinates
 * @param dimension the dimension of our space (e.g 3D)
 * @param pivot array with the coordinates of the pivot
 * @param pointsPerProc
 */
void findMPIDistances(int rank, double *dist, double **points, int dimension, const double *pivot, int pointsPerProc){

    for (int i = 0; i < pointsPerProc; ++i) {
        for (int j = 0; j < dimension; ++j) {
            dist[i] += (points[i][j] - pivot[j]) * (points[i][j] - pivot[j]);
        }
//        printf("Rank: %d    The distance is %.10f\n", rank, dist[i]);
    }
}


/**
 * Sorts the array with the points per process, holding the ones that each process wants to keep in the first indexes of the array
 * and the other part of the array has the points that it wants to give.
 *
 * @param worldSize
 * @param rank
 * @param holdThePoints
 * @param pointsPerProc
 * @param distance
 * @param median
 * @return the counter of the points that will remain in the process
 */
int partitionByMedian(int worldSize, int rank, double **holdThePoints, int pointsPerProc, double *distance, double median){
    int j = 0;
    int i = 0;
    while(j < pointsPerProc){
        if(distance[j] < median && rank < worldSize/2){
            double temp;
            double *tmp;
            temp = distance[j];
            distance[j] = distance[i];
            distance[i] = temp;

            // In the first indexes of the array we want to keep the points that will remain in the process and in the last the elements that the process wants to give
            tmp = holdThePoints[j];  // so swap the elements
            holdThePoints[j] = holdThePoints[i];
            holdThePoints[i] = tmp;
            i++;
            j++;

        } else if(distance[j] > median && rank >= worldSize/2){
            double temp;
            double *tmp;
            temp = distance[j];
            distance[j] = distance[i];
            distance[i] = temp;

            // In the last indexes of the array we want to keep the points that will be exchanged with other points from other processes
            tmp = holdThePoints[j];
            holdThePoints[j] = holdThePoints[i];
            holdThePoints[i] = tmp;
            j++;
            i++;
        } else {
            j++;
        }
    }

    return i;  // At the end of this function the i index points to the first element of the second part of the array with the points of the process that will be exchanged.
}


/**
 * Calculates the distances, gathers them to the master, which finds their median, then broadcasts it again to the
 * other processes in order to calculate how many points do they want to exchange and sort their array with their points based on their distances from the median. After that,
 * the master gathers from each process the number of points each one of them wants to exchange. The master based on the array with the counters finds which process should communicate with
 * with which process in order to the processes with rank < worldSize/2 have distances(from the pivot) < median and the processes with rank >= worldSize/2 have
 * distances(from the pivot) > median. Continuously, a buffer info is sent to each process with which one should communicate aiming at the completion of the exchanges.
 * Finally, the processes communicate and exchange their points.
 *
 * @param master master rank
 * @param rank process id
 * @param dimension
 * @param holdPoints the array with the points and their coordinates
 * @param pointsPerProc
 * @param worldSize the number of processes
 * @param communicator "holds" the processes with which the function will be called with, recursively
 */
void distributeByMedian(double *pivot,int master, int rank, int dimension, double **holdPoints, int pointsPerProc, int worldSize, MPI_Comm communicator) {
    double *distance = (double *) calloc(pointsPerProc,sizeof(double)); // for each point hold the squares of the subtractions  // FREEMEM(CHECK)

    // Find the distance of each point of the process from the pivot point chosen by the master
    findMPIDistances(rank, distance, holdPoints, dimension, pivot, pointsPerProc);


    //Allocate the master's buffer to hold the distances
    double *receiver = NULL;
    if (rank == master) {
        receiver = (double *) malloc(sizeof(double) * worldSize * pointsPerProc);   // FREEMEM(CHECK)
    }


    //The master gathers from all the processes their points' distances from the pivot
    MPI_Gather(distance, (int) pointsPerProc, MPI_DOUBLE, receiver, (int) pointsPerProc, MPI_DOUBLE, master,
               communicator);

    double median;
    if (rank == master) {
        // The master, having the array with all the distances finds their median using the "quickselect" algorithm
        median = findMpiMedian(receiver, worldSize * pointsPerProc);
//        printf("The median is %.1f\n", median);

    }


    //Broadcast the median to the processes
    MPI_Bcast(&median, 1, MPI_DOUBLE, master, communicator);


    // Call the function which calculates how many points in each process will remain
    // Sort the array holdPoints by holding in the first indexes the points that will remain in the process and in the last indexes the ones that will be exchanged
    int counter = partitionByMedian(worldSize, rank, holdPoints, pointsPerProc, distance, median);


//    printf("The counter with rank %d is %d\n", rank, counter);

    int pointsToGive = pointsPerProc - counter;  // The number of points each process wants to give
    //printf("The points to give with rank %d are %d\n", rank, pointsToGive);

    int *counterReceiver = NULL;
    if(rank == master) {
        // Allocate a buffer in which the master will store the pointsToGive from each process
        counterReceiver = (int *) malloc(worldSize * sizeof(int));  // FREEMEM
    }
    // Give the number of points each process wants to exchange to the master
    MPI_Gather(&pointsToGive, 1, MPI_INT, counterReceiver, 1, MPI_INT, master, communicator);


    ///--------------------------------------------NUMBER OF EXCHANGES AND THE INFO ARRAY--------------------------------------------
    int numberOfExchanges;  // This variable after calling the findExchanges function will contain the number of communications per process
    // and the number of how many points will send to each one of them.So diving this variable by two it will give us the rank communications
    int *exchangePerProcess;  // Holds the array with which ranks will the current process communicate and how many exchanges will it have with each one of them

    // The function has different returns whether it is called by the master process or by the others
    // That's why we call it with this condition
    exchangePerProcess = findExchanges(rank, counterReceiver, worldSize, master, communicator,&numberOfExchanges);



    /// ----------------------------------------------------------EXCHANGES------------------------------------------------------
    // For each process allocate memory for its requests
    MPI_Request *requests = (MPI_Request *) malloc((numberOfExchanges / 2) * sizeof(MPI_Request));  // FREEMEM(CHECK)

    double *sendPoints = (double *) malloc(pointsToGive * dimension * sizeof(double));  // In this array will be copied the part of the holdPoints, // FREEMEM(CHECK)
    // after it is sorted with the partitionByMedian,
    //  that contains the points that each process wants to exchange
    double *recvBuffer = (double *) malloc(pointsToGive * dimension * sizeof(double));  // In this array will be received the info from the sendPoints buffer // FREEMEM(CHECK)

    for (int i = counter; i < pointsPerProc; i++) {
        memcpy(sendPoints + (i - counter) * dimension, holdPoints[i], dimension * sizeof(double));  // copy the exchange points from holdPoints to the sendPoints
    }

    int offset = 0;  // use this offset to know until which index has info been sent

    for (int i = 0; i < numberOfExchanges / 2; i++) {
        MPI_Isend(sendPoints + offset, dimension * exchangePerProcess[2 * i + 1], MPI_DOUBLE, exchangePerProcess[2 * i],
                  0, communicator, &requests[i]);  // current rank sends a request to the rank that wants to exchange with, which is stored
        // in the exchangePerProcess[2 * i]. In the exchangePerProcess[2 * i + 1] is stored the number of points
        // it wants to send, each one of them being represented by as many coordinates as it is the space dimension


        MPI_Irecv(recvBuffer + offset, dimension * exchangePerProcess[2 * i + 1], MPI_DOUBLE, exchangePerProcess[2 * i],
                  0, communicator, &requests[i]);
        offset += dimension * exchangePerProcess[2 * i + 1];  // move offset as many elements as have been sent plus the new indexes that were sent in this loop

    }


    // Wait for the process to finish the communications
    for (int j = 0; j < numberOfExchanges / 2; j++) {
        MPI_Wait(&requests[j], MPI_STATUS_IGNORE);
    }


    //Update the new holdPoints array
    for (int k = counter; k < pointsPerProc; k++) {
        for (int l = 0; l < dimension; l++) {
            holdPoints[k][l] = recvBuffer[l + (k - counter) * dimension];
        }
    }

    if(rank == master){
        free(receiver);
        free(counterReceiver);
    }
    free(distance);
    free(exchangePerProcess);
    free(sendPoints);
    free(recvBuffer);
    free(requests);
    /// --------------------------------------------- Recursive call ---------------------------------------------


    // Termination condition
    if(worldSize == 2) {
        return;
    }

    MPI_Comm smallDistComm; // Define the new communicator after splitting the current worldSize in two for the recursive call
    MPI_Comm bigDistComm;  // The new communicator for the processes with their distances > median

    if (rank < worldSize / 2) {

        MPI_Comm_split(communicator, 1, rank,
                       &smallDistComm);  // Split the communicator in two parts for the processes with ranks > worldSize/2
        // that have distances < median

        MPI_Comm_rank(smallDistComm, &rank);  // Get the current rank after the split
        MPI_Comm_size(smallDistComm, &worldSize);  // Get the current worldSize after the split

        //        printf("New world size %d small comm\n", worldSize);
        //        printf("New world rank %d small comm\n", rank);

        // Call recursively for the new rank, size and communicator
        distributeByMedian(pivot, 0, rank, dimension, holdPoints, pointsPerProc, worldSize, smallDistComm);

        MPI_Comm_free(&smallDistComm);

    } else {
        // recursive call for the ranks >= worldSize/2 that have distances > median
        MPI_Comm_split(communicator, 2, rank, &bigDistComm);
        MPI_Comm_rank(bigDistComm, &rank);
        MPI_Comm_size(bigDistComm, &worldSize);

        //        printf("New world size %d big comm\n", worldSize);
        //        printf("New world rank %d big comm\n", rank);
        //printf("My rank is: %d ", rank);
        distributeByMedian(pivot, 0, rank, dimension, holdPoints, pointsPerProc, worldSize, bigDistComm);

        MPI_Comm_free(&bigDistComm);
    }

}
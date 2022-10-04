#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include "mpi_quick_select.h"
#include "mpi_knn.h"
#include "sequential/sequential_vptree.h"


/**
 *
 * This function finds the exchanges between the processes according to the following algorithm:
 * counterReceiver: The array in which the master has gathered the number of points each process wants to exchange
 * The main idea is to compare one by one the processes with ranks smaller than worldSize/2 to the processes with ranks bigger than worldSize/2,
 * based on how many points each one wants to exchange. For example, let's consider the array with processes' ids
 * and the counterReceiver array to be:
 * ranks = [0 1 2 3 ]
 * counterReceiver = [2 3 5 0]
 * The aim is to give the points of the first half of the array to the second half. So, in the beginning we compare the counterReceiver[0] = 2 to
 * counterReceiver[0 + worldSize/2] = counterReceiver[0 + 4/2] = counterReceiver[2] = 5.
 * We see that 2 < 5 which means process with rank 0 should give to process with rank 2 2 points and vice versa.
 * The process with rank 0 has completed its exchanges in contrast to the process with rank 2. The last one wants to give 5 - 2 = 3 more points to another process.
 * So, now we have to compare the (i + 1)th element of the counterReceiver to the left points of the incomplete . We observe that counterReceiver[1] = 3
 * is equal to the number of points that has remained in the process with rank 2. So, the process with rank 1 should give to process with rank 2
 * 3 points and take from it 3 points as well.
 * Then, both processes have completed their exchanges.This algorithm continues till all processes have exchanged as many elements as they wanted.
 *
 * @param rank
 * @param counterReceiver The array with the number of points, each process wants to give, that has bees sent to the master
 * @param worldSize
 * @param master
 * @param communicator the cluster of processes that calls the function
 * @param numberOfExchanges points to the length of the array with the infos for each rank
 * @return the array with the infos per process
 *
 */
int *findExchanges(int rank, int *counterReceiver, int worldSize, int master, MPI_Comm communicator, int *numberOfExchanges){
    int *helperIndex;

    /// IN THE EVEN INDEXES WILL BE STORED THE NUMBER OF POINTS THAT WILL SEND IN THE RANK THAT LIES IN THE ODD INDEXES
    int **info;

    if(rank == master) {
        // Allocate memory for the array that will contain with which process should the current process communicate and how many points will they exchange
        info = (int **) malloc(worldSize * sizeof(int *));      // FREEMEM (CHECK)
        for (int i = 0; i < worldSize; i++) {
            info[i] = (int *) malloc(2 * worldSize * sizeof(int));
        }

        // This array holds the index for its line in the info array in which we should write the next element
        helperIndex = (int *) calloc(worldSize, sizeof(int));       // FREEMEM (CHECK)

        int nextIndex = worldSize / 2;  // Indicates the index of the second half of the array
        for (int i = 0; i < worldSize / 2; i++) {

            while (counterReceiver[i] > 0) {  // Repeat until each process gets the info for its exchanges. When counterReceiver[i] = 0 the rank doesn't want to exchange anymore or not at all

                if (counterReceiver[i] == counterReceiver[nextIndex]) {  // Compare elements of the first half of array to the second half's elements

                    info[i][helperIndex[i]] = nextIndex;  // Store for the rank i that it should give points to the nextIndex(it lies in the second half of the array as it is initialized with the value worldSize/2

                    helperIndex[i]++;  // Increase the index which indicates the last index we wrote sth about the i rank

                    info[i][helperIndex[i]] = counterReceiver[i];  // In the new index store how many points will the i rank send to the newIndex rank

                    helperIndex[i]++;

                    info[nextIndex][helperIndex[nextIndex]] = i;  // Store the info exchanges for the nextIndex rank as well
                    // Send to the rank i

                    helperIndex[nextIndex]++;  // Increase the index which indicates the last index we wrote sth about the nextIndex rank

                    info[nextIndex][helperIndex[nextIndex]] = counterReceiver[i];  // Send to rank i as many elements as it sent to rank nextIndex

                    helperIndex[nextIndex]++;

                    counterReceiver[i] = 0; // The i rank has gotten the info for its exchanges

                    counterReceiver[nextIndex] = 0;  // The nextIndex has gotten the info for its exchanges too

                    nextIndex++;  // Point to the new index we' re going to write
                } else if (counterReceiver[i] < counterReceiver[nextIndex]) {  // If the current rank wants to give less elements than the nextIndex requires
                    info[i][helperIndex[i]] = nextIndex;

                    helperIndex[i]++;

                    info[i][helperIndex[i]] = counterReceiver[i];

                    helperIndex[i]++;

                    info[nextIndex][helperIndex[nextIndex]] = i;

                    helperIndex[nextIndex]++;

                    info[nextIndex][helperIndex[nextIndex]] = counterReceiver[i];

                    helperIndex[nextIndex]++;

                    counterReceiver[nextIndex] -= counterReceiver[i]; // the points remained that the nextIndex rank should exchange to finish its exchanges

                    counterReceiver[i] = 0; // The i rank has gotten all the info
                } else if (counterReceiver[i] > counterReceiver[nextIndex] && counterReceiver[nextIndex] != 0) {
                    info[i][helperIndex[i]] = nextIndex;

                    helperIndex[i]++;

                    info[i][helperIndex[i]] = counterReceiver[nextIndex];

                    helperIndex[i]++;

                    info[nextIndex][helperIndex[nextIndex]] = i;

                    helperIndex[nextIndex]++;

                    info[nextIndex][helperIndex[nextIndex]] = counterReceiver[nextIndex];

                    helperIndex[nextIndex]++;

                    counterReceiver[i] -= counterReceiver[nextIndex];  // The points left in i rank for which we should find with which next process should they be exchanged

                    counterReceiver[nextIndex] = 0;  // The nextIndex has gotten all its info

                    nextIndex++;
                } else {
                    nextIndex++;  // In case we have counterReceiver[i] > counterReceiver[nextIndex] && counterReceiver[nextIndex] == 0
                }
            }
        }

//        printf("The array with exchanges is: \n");
//        for (int i = 0; i < worldSize; i++) {
//            printf("The i is %d ", i);
//            for (int j = 0; j < helperIndex[i]; ++j) {
//                printf("%d ", info[i][j]);
//            }
//            printf("\n");
//        }
    }

    // Allocate memory for each process for the info it will receive
    int *infoPerProc;

    // If you are the master rank send the infos for the exchanges to the other processes
    if(rank == master) {
        // Allocate memory for the request objects
        MPI_Request *requests = (MPI_Request *) malloc((worldSize - 1) * sizeof (MPI_Request));     // FREEMEM (CHECK)

        // For each rank in the worldSize send a request from the master to the other processes
        for (int i = 1; i < worldSize; i++) {
            MPI_Isend(info[i], helperIndex[i], MPI_INT, i, 10, communicator, &requests[i - 1]);
        }

//        printf("\nThe rank %d has receive buffer: ", rank);
//        for(int i = 0; i < helperIndex[0]; i++){
//            printf("%d ", info[0][i]);
//        }
//        printf("\n");

        for(int i = 0; i < worldSize - 1; i++){
            MPI_Wait(&requests[i], MPI_STATUS_IGNORE);  // Wait for the requests to finish
        }
        *numberOfExchanges = helperIndex[0];  // The length of the array for the infos of the rank master
        return info[0];

    } else {
        int elements;  // The number of elements that have been requested to be sent
        MPI_Status status;
        MPI_Probe(master, 10, communicator, &status);  // Ask MPI to give you the size of the message

        MPI_Get_count(&status, MPI_INT, &elements);


        infoPerProc = (int *) malloc(elements * sizeof (int));  // Allocate memory for the receiver buffer of the processes that will get infos from the master     // FREEMEM (CHECK)

        MPI_Recv(infoPerProc, elements, MPI_INT, master, 10, communicator, MPI_STATUS_IGNORE);  // Receive the array per process!


//        printf("\nThe rank %d has receive buffer: ", rank);
//        for(int i = 1; i < elements + 1; i++){
//            printf("%d ", infoPerProc[i]);
//        }
//        printf("\n");
        *numberOfExchanges = elements;

        free(helperIndex);
        for (int i = 0; i < worldSize; ++i) {
            free(info[i]);
        }
        free(info);
        return infoPerProc;
    }
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
void findDistance(double *dist, double **points, int dimension, const double *pivot, int pointsPerProc){

    for (int i = 0; i < pointsPerProc; ++i) {
        for (int j = 0; j < dimension; ++j) {
            dist[i] += (points[i][j] - pivot[j]) * (points[i][j] - pivot[j]);
        }
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
 * check if the processes with ranks between 0 and worldSize(exclusive) have distances bigger than the biggest distance they received from the previous process(sorted processes)
 * so as to know if the processes in the end have the right points
 *
 * @param points
 * @param pointsPerProc
 * @param dimension
 * @param pivot
 */
void testFunction(double **points, int pointsPerProc, int dimension, double *pivot) {

    int rank;
    int worldSize;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);

    double *dist = (double *) malloc(pointsPerProc * sizeof (double));

    findDistance(dist, points, dimension, pivot, pointsPerProc);

    double maxDist = dist[0];
    double minDist = dist[0];
    for (int i = 0; i < pointsPerProc; ++i) {
        if(dist[i] > maxDist){
            maxDist = dist[i];
        }

        if(dist[i] < minDist) {
            minDist = dist[i];
        }
    }

    if (rank == 0){
        MPI_Send(&maxDist, 1, MPI_DOUBLE, rank + 1, 50, MPI_COMM_WORLD);
    }
    else if (rank == worldSize - 1){
        double prevMax;

        MPI_Recv(&prevMax, 1, MPI_DOUBLE, rank - 1, 50, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        if (prevMax < maxDist) {
            printf("Rank: %d SUCCESS\n", rank - 1);
            printf("Rank: %d SUCCESS\n", rank);
        } else {
            printf("Rank: %d FAILURE\n", rank - 1);

        }


    } else {
        double prevMax;
        MPI_Request request_1;
        MPI_Request request_2;

        MPI_Irecv(&prevMax, 1, MPI_DOUBLE, rank - 1, 50, MPI_COMM_WORLD, &request_1);
        MPI_Isend(&maxDist, 1, MPI_DOUBLE, rank + 1, 50, MPI_COMM_WORLD, &request_2);

        MPI_Wait(&request_1, NULL);
        MPI_Wait(&request_2, NULL);

        if (prevMax < maxDist) {
            printf("Rank: %d SUCCESS\n", rank - 1);
        } else {
            printf("Rank: %d FAILURE\n", rank - 1);

        }
    }

}

bool isPowerOfTwo(int64_t number){
    return (ceil(log2((double)number)) == floor(log2((double)number)));
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

    double *distance = NULL;
    if(rank != master) {
        distance = (double *) calloc(pointsPerProc,sizeof(double)); // for each point hold the squares of the subtractions      // FREEMEM

        // Find the distance of each point of the process from the pivot point chosen by the master
        findDistance(distance, holdPoints, dimension, pivot, pointsPerProc);
    }


    //Allocate the master's buffer to hold the distances
    double *receiver = NULL;
    if (rank == master) {

        // for the master process the array that will hold the distances is the receiver buffer
        receiver = (double *) malloc(sizeof(double) * worldSize * pointsPerProc);       // FREEMEM (CHECK)

        // Find the distance of each point of the process from the pivot point chosen by the master
        findDistance(receiver, holdPoints, dimension, pivot, pointsPerProc);
    }


    //The master gathers from all the processes their points' distances from the pivot
    MPI_Gather(
            distance,
            (int) pointsPerProc,
            MPI_DOUBLE,
            receiver + pointsPerProc,
            (int) pointsPerProc,
            MPI_DOUBLE,
            master,
            communicator
    );

    double median;
    if (rank == master) {
        // The master, having the array with all the distances finds their median using the "quickselect" algorithm
        median = findMedian(receiver, worldSize * pointsPerProc);

    }


    //Broadcast the median to the processes
    MPI_Bcast(&median, 1, MPI_DOUBLE, master, communicator);


    // Call the function which calculates how many points in each process will remain
    // Sort the array holdPoints by holding in the first indexes the points that will remain in the process and in the last indexes the ones that will be exchanged
    int counter = partitionByMedian(worldSize, rank, holdPoints, pointsPerProc, receiver, median);
    int pointsToGive = pointsPerProc - counter;  // The number of points each process wants to give


    int *counterReceiver = (int *) malloc(worldSize * sizeof(int));  // Allocate a buffer in which the master will store the pointsToGive from each process     // FREEMEM

    // Give the number of points each process wants to exchange to the master
    MPI_Gather(&pointsToGive, 1, MPI_INT, counterReceiver, 1, MPI_INT, master, communicator);

//    // DEALLOCATE
//    if(rank == master)
//        free(receiver);
//    else
//        free(distance);


    ///--------------------------------------------NUMBER OF EXCHANGES AND THE INFO ARRAY--------------------------------------------
    int numberOfExchanges;  // This variable after calling the findExchanges function will contain the number of communications per process
    // and the number of how many points will send to each one of them.So dividing this variable by two it will give us the rank communications
    int *exchangePerProcess;  // Holds the array with which ranks will the current process communicate and how many exchanges will it have with each one of them

    // The function has different returns whether it is called by the master process or by the others
    // That's why we call it with this condition
    exchangePerProcess = findExchanges(rank, counterReceiver, worldSize, master, communicator,&numberOfExchanges);



    /// ----------------------------------------------------------EXCHANGES------------------------------------------------------
    // For each process allocate memory for its requests
    MPI_Request *requests = (MPI_Request *) malloc((numberOfExchanges / 2) * sizeof(MPI_Request));      // FREEMEM (CHECK)

    double *sendPoints = (double *) malloc(pointsToGive * dimension * sizeof(double));  // In this array will be copied the part of the holdPoints,     // FREEMEM
                                                                                            // after it is sorted with the partitionByMedian,
                                                                                            //  that contains the points that each process wants to exchange

    double *recvBuffer = (double *) malloc(pointsToGive * dimension * sizeof(double));  // In this array will be received the info from the sendPoints buffer   // FREEMEM

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


    /// --------------------------------------------- Recursive call ---------------------------------------------
    // Termination condition
    if(worldSize == 2){

        // REACHING THIS POINT MEANS THAT EACH PROCESS WITH RANK < WORLD SIZE / 2 HAS POINTS WITH DISTANCES FROM THE PIVOT
        // SMALLER THAN MEDIAN AND PROCESSES WITH RANK > WORLD SIZE / 2 HAVE POINTS WITH DISTANCES > MEDIAN
        // BUT, THIS DOESN'T MEAN THAT EACH PROCESS HAS ITS POINTS SORTED
        /// SO, IN THIS PHASE THE SEQUENTIAL FUNCTION BUILD VP TREE IS CALLED SO AS TO SORT THE POINTS AND CREATE
        /// THE BINARY TREE CORRECTLY

        vptree initial;  // Each process creates its struct of vptree
        initial.inner = NULL;
        initial.outer = NULL;
        initial.start = 0;
        initial.stop = pointsPerProc - 1;

        if(rank != master) {
            buildVPTree(&initial, holdPoints, distance, dimension, pointsPerProc);
        }
        else{
            buildVPTree(&initial, holdPoints, receiver, dimension, pointsPerProc);
        }

        // Given that each process has created its tree, the master process concatenates them to create the binary tree
        // each process sends its tree to the master

        return;
    }



    //else
    MPI_Comm smallDistComm; // Define the new communicator after splitting the current worldSize in two for the recursive call
    MPI_Comm bigDistComm;  // The new communicator for the processes with their distances > median

    if(rank < worldSize/2){
        MPI_Comm_split(communicator, 1, rank, &smallDistComm);  // Split the communicator in two parts for the processes with ranks > worldSize/2
        // that have distances < median

        MPI_Comm_rank(smallDistComm, &rank);  // Get the current rank after the split
        MPI_Comm_size(smallDistComm, &worldSize);  // Get the current worldSize after the split

//        printf("New world size %d small comm\n", worldSize);
//        printf("New world rank %d small comm\n", rank);

        free(requests);
        free(exchangePerProcess);
        // Call recursively for the new rank, size and communicator
        distributeByMedian(pivot, 0, rank, dimension, holdPoints, pointsPerProc, worldSize, smallDistComm);
    }
    else{
        // recursive call for the ranks >= worldSize/2 that have distances > median
        MPI_Comm_split(communicator, 2, rank, &bigDistComm);
        MPI_Comm_rank(bigDistComm, &rank);
        MPI_Comm_size(bigDistComm, &worldSize);

//        printf("New world size %d big comm\n", worldSize);
//        printf("New world rank %d big comm\n", rank);

        free(requests);
        free(exchangePerProcess);


        distributeByMedian(pivot, 0, rank, dimension, holdPoints, pointsPerProc, worldSize, bigDistComm);

    }
}


int main(int argc, char **argv) {
    int size;
    int rank;
    vptree initial;


    MPI_Init(&argc, &argv);
    double start, end;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* Read the binary file with our data back in */
    // The first two 64-bit numbers from the file give us the number of our points and the space's dimension
    int64_t numberOfPoints;
    int64_t dimension;
    long int pointsPerProc;
    double **holdThePoints;

    FILE *fh = fopen(argv[1], "rb");
    if (fh != NULL) {
        fread(&dimension, sizeof(int64_t), 1, fh);
        fread(&numberOfPoints, sizeof(int64_t), 1, fh);

        // hold number of points that is a power of 2
        if(!isPowerOfTwo(numberOfPoints)){
            int64_t power = 1;
            while(power < numberOfPoints)
                power*=2;

            numberOfPoints = power / 2;
            if(rank == 0)
                printf("The number of points is: %ld\n", numberOfPoints);
        }
        else{
            if(rank == 0)
                printf("The number of points is: %ld\n", numberOfPoints);
        }


        pointsPerProc = numberOfPoints / size;
        printf("Rank: %d, Points per proc: %ld\n", rank, pointsPerProc);

        // hold the points with their coordinates read from the file
        holdThePoints = (double **) malloc(pointsPerProc * sizeof(double *));       // FREEMEM (CHECK)
        for (int i = 0; i < pointsPerProc; i++) {
            holdThePoints[i] = (double *) malloc(dimension * sizeof(double));
        }

        // Seek where is the pointer in our file
        fseek(fh, rank * dimension * pointsPerProc * sizeof(double), SEEK_CUR);

        for (int i = 0; i < pointsPerProc; i++) {
            double *x2 = (double *) malloc(sizeof(double) * dimension); // represents only one point        // FREEMEM (CHECK)

            fread(x2, sizeof(double), dimension, fh);
            holdThePoints[i] = x2;
        }
        fclose(fh);

        initial.inner = NULL;
        initial.outer = NULL;
        initial.start = 0;
        initial.stop = numberOfPoints - 1;
        initial.vpPoint = (double *) calloc(dimension, sizeof(double));


        // The master process chooses one pivot and broadcasts it to the other processes
        if(rank == 0) {
            //printf("The pivot is: ");
            for (int j = 0; j < dimension; j++) {
                initial.vpPoint[j] = holdThePoints[0][j];
                //printf("%.10f ", pivot[j]);
            }
            //putchar('\n');
        }

//        MPI_Barrier(MPI_COMM_WORLD);

        if(rank == 0){
            start = MPI_Wtime();
        }
        //Broadcast the pivot to the processes
        MPI_Bcast(initial.vpPoint, (int)dimension, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        distributeByMedian(initial.vpPoint, 0, rank, (int)dimension, holdThePoints, (int)pointsPerProc, size, MPI_COMM_WORLD);


        if(rank == 0){
            end = MPI_Wtime();
            printf("The time is: %.4f\n", end - start);
        }

        testFunction(holdThePoints, (int)pointsPerProc, (int)dimension, initial.vpPoint);

        free(initial.vpPoint);
    }


    for (int i = 0; i < pointsPerProc; ++i) {
        free(holdThePoints[i]);
    }
    free(holdThePoints);
    MPI_Finalize();

    return 0;
}

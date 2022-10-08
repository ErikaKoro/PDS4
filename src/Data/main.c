#include <stdint.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include "timer.h"

int main(int argc, char **argv){
    if(argc != 4){
        printf("gimme args you stupid\n");
        return -1;
    }
    // printf("%s\n", argv[3]);     // DEBUG COMMENT
    srand(time(NULL));
    /* Create the file */
    int64_t dimension= atoi(argv[1]);
    int64_t numberOfPoints = atoi(argv[2]);

    FILE *fh = fopen (argv[3], "wb");
    // printf("%p\n", fh);      // DEBUG COMMENT
    if (fh != NULL) {
        fwrite (&dimension, sizeof (int64_t), 1, fh);
        fwrite (&numberOfPoints, sizeof (int64_t), 1, fh);
        double range = 200.00;
        double div = RAND_MAX / range;
        for (int j = 0; j < numberOfPoints; j++){
            for(int i = 0; i < dimension; i++){

                double coordinates = -100 + (rand() / div);
                // printf("%.10f ", coordinates);
                fwrite(&coordinates, sizeof(double), 1, fh);
            }
            //printf("\n");
        }
        fclose (fh);
    }
    else{
        printf("Hi\n");
    }
    printf("\n");
    printf("\n");
    printf("\n");

    return 0;
}

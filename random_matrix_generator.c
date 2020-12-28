//
//  random_matrix_generator.c
//  
//
//  Created by Loris Cino on 27/11/2020.
//


#include "stdio.h"
#include "stdlib.h"
#include <sys/time.h>
#include <stdint.h>

int print_matrix(int** matrix, int size);
int** malloc_matrix(int size);

int main(int argc, char **argv)
{
    int size = atoi(argv[1]);
    int** random_matrix =malloc_matrix(size);
    int i, o;
    
    struct timeval t;
    gettimeofday(&t, 0);

    srand(t.tv_usec);
    for(o = 0; o<size; o++)
        for(i = 0; i<size; i++)
            random_matrix[o][i] = (rand() % size) + 1;
    
    print_matrix(random_matrix, size);
    return 0;
}


 int** malloc_matrix(int size){
     int** matrix =malloc(size*sizeof(int*));
     for(int i=0;i<size;++i)
         matrix[i]=malloc(size*sizeof(int));

     return matrix;
 }

int print_matrix(int** matrix, int size){
    for(int i = 0; i < size; i++) {
        for(int j = 0; j < size; j++)
                printf("%d\t", matrix[i][j]);
    printf("\n");
    }
    return 1;
}

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "mpi.h"
#include "math.h"
#include <sys/time.h>
#include <stdint.h>
#include <ctype.h>
#include <unistd.h>

#define PART_PER_PROCESS 32
#define TERMINATE 31000
#define TERMINATE_NO_SUCC 31001

typedef struct {
    short int value;
    char state;
    //Il valore e':
    //'b' if the cell is deleted
    //'w' if the cell cannot be deleted
    //'u' if we don't know
}block;

//size of the matrix
int size;
block** global_matrix;
//count the node examinated
int node_count = 0;
int rank;
//Total number of process
int n_process = 0;
MPI_Datatype block_datatype;
MPI_Datatype col_datatype;
MPI_Datatype ncol_type;
int* sendcounts;
int* displs;
int* sendcountsCol;
int* displsCol;
int my_size;


int read_matrix_from_file(block*** matrix, char* file);
int read_matrix_from_stdin(block*** m);
int print_matrix(block** matrix);
int print_matrix_result(block** matrix);
int solve_hitori(block** matrix,int i,int unknown);
int check_row_and_column(block*** m, int* unknown);
int find_connection_of_white(int row,int col,block** matrix);
int check_adjacent_rules(block*** m, int unknown);
int check_neighborhood(int row,int col,block** matrix, int** flag_vector);
//FUNZIONI USILIARIE
void free_matrix(block** arr2D);
block** malloc_matrix();
int copy_matrix(block*** input_matrix, block*** output_matrix);
int apply_rule(block*** matrix, int unknown);
double logbase(int base, int x);
int isNumber(char number[]);
//PARALLEL FUNCTION
int scatterElements(block*** matrix, block*** my_elements,int my_size);
int scatterElementsColumn(block*** matrix, block*** my_elements,int my_size);
int print_sub_matrix(block** my_elements,int my_size);
int print_once(int prev,int next, block** mypartition);
int check_row_and_column_par(block*** m, int unknown);
int gatherUnknown(int* unknowns,int unknown);

int main(int argc, char **argv){
    char* input_file;

    //Initialization of the MPI data structures
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &n_process);

    //Controllo sull'input
    if(argc<3 || argc>3){
        printf("\nThe program needs two arguments: \n1)The file containing the puzzle\n2) The dimension o f the puzzle");
        return -1;
    }
    if( isNumber(argv[2]) ) size = atoi(argv[2]);
    else{
      printf("Il secondo argomento deve essere un intero.\n");
      return -1;
    }

    input_file = argv[1];

    if(!strcmp(input_file ,"stdin"))
      	read_matrix_from_stdin(&global_matrix);
    else if( read_matrix_from_file(&global_matrix,input_file) == -1) {
      	MPI_Finalize();
        return -1;
    }

    my_size = ((rank+1) * size)/n_process - (rank * size)/n_process;
    int unknown = size*size;
    //Calcolo la dimensione della sottomatrice di ogni processo
    sendcounts = (int *) malloc(n_process*sizeof(int));
    //Calcolo lo spiazzamento di ogni processo
    displs = (int *) malloc(n_process*sizeof(int));
    //Stessa cosa ma per le colonne
    sendcountsCol = (int *) malloc(n_process*sizeof(int));
    displsCol = (int *) malloc(n_process*sizeof(int));
    for(int i = 0; i < n_process; i++){
      sendcounts[i] = (((i+1) * size)/n_process - (i * size)/n_process)*size;
      sendcountsCol[i] = (((i+1) * size)/n_process - (i * size)/n_process);
      //if(!rank) printf(" sendcounts[%d] = %d \n",i , sendcounts[i]);
      displs[i] =( i * size)/n_process * size;
      displsCol[i] =( i * size)/n_process;
      //if(!rank) printf(" displs[%d] = %d \n",i , displs[i]);
    }
    //Making the datatype
    MPI_Aint displ[2];
    int bl[2] = {1, 1};
    MPI_Datatype dt[2] = {MPI_INT, MPI_CHAR};

    MPI_Get_address(&(global_matrix[0][0].value), &displ[0]);
    MPI_Get_address(&(global_matrix[0][0].state), &displ[1]);

    //printf("%d)Creo il tipo di dato. \n", rank);

    displ[1] -= displ[0];
    displ[0] = 0;

    if (MPI_Type_create_struct(2, bl, displ, dt, &block_datatype) == -1){
      	printf("Creation of the datatype failed.\n");
      	return -1;
    }
    MPI_Type_commit(&block_datatype);

    //'Single_Col' Data Type
    MPI_Type_vector(size, 1, size, block_datatype, &col_datatype);
    MPI_Type_commit(&col_datatype);

    MPI_Type_create_resized(col_datatype, 0 , sizeof(block),  &ncol_type);
    MPI_Type_commit(&ncol_type);

    //if (!rank) printf("Solving\n");
		//Per misurare il tempo
		MPI_Barrier(MPI_COMM_WORLD);
		double time = - MPI_Wtime();

    if(rank == 0){
      //print_matrix(backup_matrix);
      if(!solve_hitori(global_matrix,0,unknown)){
          printf("Solution found.\n");
          print_matrix(global_matrix);
          fprintf(stderr,"Examinated %d node by process %d.\n", node_count, rank);
      }else{
          printf("No solution.\n");
          fprintf(stderr,"Examinated %d node by process %d.\n", node_count, rank);
      }
    }else{
      int message = 1;
      MPI_Status recv_status;

      //Il numero di righe/colonne per ogni processo
      int my_size = ((rank+1) * size)/n_process - (rank * size)/n_process;

      //Prima verrà utilizzato per delle righe e poi per delle colonne
      block** my_elements = calloc (size, sizeof(block*));
      my_elements[0] = calloc (size * my_size, sizeof(block));

      for (int i = 1; i < my_size; i++) my_elements[i] = my_elements[i-1] + size;

      while(message){
        unknown = 0;

        scatterElements(&my_elements, &my_elements,my_size);

        if( check_row_and_column(&my_elements, &unknown) ) unknown = 1;

        //printf("%d)Il mio unknown è %d \n", rank, unknown);
        scatterElementsColumn(&my_elements,&my_elements,my_size);

        print_once((rank + n_process - 1)%n_process, (rank + 1)%n_process, my_elements);

        if( check_row_and_column(&my_elements, &unknown) ) unknown = 1;

        gatherUnknown(NULL, unknown);

        sleep(100);

        if ( MPI_Recv(&message, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &recv_status) == -1){
            printf("MPI_Recv failed.\n");
            MPI_Abort(MPI_COMM_WORLD, -1);
            return -1;
        }

        printf("%d) message = %d. \n", rank, message);
      }

      free_matrix(my_elements);
    }

    //printf("%d) Terminate.\n", rank);
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);
    time += MPI_Wtime();

    if (rank == 0){
      fprintf(stderr, "dim = %d)Execution time: %f\n",size, time);
    }

    MPI_Finalize();
    free(sendcounts);
    free(displs);
    return 0;
}

int solve_hitori(block** matrix,int i,int unknown){
    //My block's row and colum
    int r = i / size;
    int c = i % size;
    //If the white ipotesis is not correct
    block** backup_matrix = malloc_matrix();
    int unknown_copy;

    //If I have 0 unknown I have solved the puzzle :)
    if(unknown == 0){
        //printf("Returnung to main.\n" );
        global_matrix = matrix;
		    free_matrix(backup_matrix);
        return 0;
    }


    //Find the next unknown state block
    while(matrix[r][c].state != 'u'){
        i++;
        r = i / size;
        c = i % size;
    }

    //printf("Now there are %d unkowns \n", unknown);
    //tring to set current block to white
    //printf("Setting %d element to white. \n", i+1);
    node_count++;
    matrix[r][c].state = 'w';
    //Make a safe copy
    copy_matrix(&matrix, &backup_matrix);
    //anyway I will decrease variable unknown or fail
    unknown = unknown - 1;
    unknown_copy = unknown;

    //If the function apply rule return -1 means that it failed
    unknown = apply_rule(&backup_matrix, unknown);
    if( unknown != -1){
      //printf("Applied succesfully rules. \n");
  		int return_code = solve_hitori(backup_matrix, i + 1, unknown);
  		if( return_code == 0) return 0;
    }

    //failed white trying black
    //printf("Failed apply to rules. \n");
    //printf("Setting %d element to black. \n", i+1);
    node_count++;
    matrix[r][c].state = 'b';
    unknown = unknown_copy;

    //erase wrong part
    free_matrix(backup_matrix);

    unknown = apply_rule(&matrix, unknown);
	if( unknown != -1){
    //printf("Applied succesfully rules. \n");
		int return_code = solve_hitori(matrix, i + 1, unknown);
		if( return_code == 0) return 0;
  }

  return -1;
}

int check_row_and_column(block*** m, int* u){
    block** matrix = *m;
    int unknown = *u;

    //printf("%d) unknown = %d\n", rank, unknown);
    //printf("%d) displs[rank] = %d\n", rank, displs[rank]);

    for(int i = displs[rank]; i < sendcounts[rank]; i++){
        //My current block's row and colum
        int r = i / size;
        int c = i % size;

        //if the current block's state is unknown i can't do anything
        if(matrix[r][c].state != 'u'){
            //Let's do the same thing with for the rows
            for(int cnew = 0; cnew < size; cnew++ ){
                if(cnew != c)
                    if(matrix[r][c].value == matrix[r][cnew].value){
                        if(matrix[r][c].state == 'w' && matrix[r][cnew].state == 'w') return 1;
                        else if(matrix[r][c].state == 'w' && matrix[r][cnew].state == 'u'){
                                    unknown = unknown - 1;
                                    matrix[r][cnew].state = 'b';
                        }
                    }
              }

        }
    }
    *u = unknown;
    return 0;
}

int check_adjacent_rules(block*** m, int unknown){
    //Find if there are two black cell adjacent
    block** matrix = *m;
    int first_white = 1;

    for(int i = 0; i < size * size; i++){
            //My current block's row and colum
            int r = i / size;
            int c = i % size;

            //if the state of the current block is black
            if(matrix[r][c].state ==  'b'){
                //For all the adjacent cell apply the second rule, but before I check that the cell exist
                if (r > 0 ){
                    if (matrix[r-1][c].state ==  'b')
                        return -1;
                    else if(matrix[r-1][c].state ==  'u'){
                            matrix[r-1][c].state =  'w';
                            unknown = unknown - 1;
                        }
                }if(c > 0){
                    if (matrix[r][c-1].state ==  'b')
                        return -1;
                    else if(matrix[r][c-1].state ==  'u'){
                            matrix[r][c-1].state =  'w';
                            unknown = unknown - 1;
                        }
                }if(r + 1 < size){
                    if (matrix[r+1][c].state ==  'b')
                        return -1;
                    else if(matrix[r+1][c].state ==  'u'){
                        matrix[r+1][c].state =  'w';
                        unknown = unknown - 1;
                    }
                }if(c + 1 < size){
                        if (matrix[r][c+1].state ==  'b')
                            return -1;
                        else if(matrix[r][c+1].state ==  'u'){
                            matrix[r][c+1].state =  'w';
                            unknown = unknown - 1;
                        }
                }

            }else if (first_white) {
                first_white--;
                if(find_connection_of_white(r, c, matrix))
                    return -1;
            }
    }

    return unknown;
}

int find_connection_of_white(int row,int col,block** matrix){
    int* flag_vector = (int *) calloc(size*size,sizeof(int));
    check_neighborhood(row, col, matrix, &flag_vector);
    for(int j = 0; j < size*size; j++){
        int r = j / size;
        int c = j % size;
        if(flag_vector[j] == 0 && matrix[r][c].state == 'w'){
			free(flag_vector);
            return -1;
        }
    }
	free(flag_vector);
    return 0;
}

int check_neighborhood(int row,int col,block** matrix, int** flag_vector){
    int* f_vector = *flag_vector;
    //setting current block as visited
    f_vector[row*size+col] = 1;
    //if the block exist AND if it is not black AND it isnt visited yet
    if( row > 0 && matrix[row-1][col].state != 'b' && f_vector[(row-1)*size+col] != 1){
        check_neighborhood(row - 1, col, matrix, flag_vector);
    }if( col > 0 && matrix[row][col-1].state != 'b' && f_vector[row*size+col -1] != 1){
        check_neighborhood(row, col - 1, matrix, flag_vector);
    }if( row + 1 < size && matrix[row + 1][col].state != 'b' && f_vector[(row+1)*size+col] != 1){
        check_neighborhood(row+1, col, matrix, flag_vector);
    }if( col + 1 < size && matrix[row][col + 1].state != 'b' && f_vector[row*size+col+1] != 1){
        check_neighborhood(row, col+1, matrix, flag_vector);
    }
    return 0;
}

int check_row_and_column_par(block*** matrix,int unknown){
  int* unknowns = (int *) malloc(n_process*sizeof(int));
  block** my_elements = calloc (size, sizeof(block*));
  my_elements[0] = calloc (size * my_size, sizeof(block));

  for (int i = 1; i < my_size; i++) my_elements[i] = my_elements[i-1] + size;

  //printf("Avvio scatter. \n");
  //print_matrix(*matrix);

  scatterElements(matrix,&my_elements, my_size);

  if( check_row_and_column(&my_elements, &unknown) ) {
    gatherUnknown(unknowns, unknown );
    free(unknowns);
    return -1;
  }

  scatterElementsColumn(matrix,&my_elements, my_size);

  print_once(n_process-1, 1, my_elements);

  if( check_row_and_column(&my_elements, &unknown) ) {
    //La gather serve per sbloccare gli altri processi dall'attesa
    gatherUnknown(unknowns, unknown );
    free(unknowns);
    return -1;
  }

  gatherUnknown(unknowns, unknown );

  //Aggiorno il nuovo valore degli unknown
  //Poichè i numeri sono negativi equivale ad una sottrazione
  for(int i = 1; i < n_process; i++){
    printf("%d) unknown = %d\n", i, unknown );
    //1 vuol dire che le regole non sono rispettate
    if(unknowns[i] == 1) {
      free(unknowns);
      return -1;
    }
    unknown = unknown + unknowns[i];
  }

  sleep(100);

  free(unknowns);
  return unknown;
}

int apply_rule(block*** matrix,int unknown){
    int prev_u;
    do{
        prev_u = unknown;
        unknown = check_row_and_column_par(matrix, unknown);
        if (unknown < 0)
            return -1;
		    unknown = check_adjacent_rules(matrix,unknown);
        if (unknown < 0)
            return -1;
    }while(prev_u != unknown);
    //the do while cycle continue until appling the rules makes no change
  return unknown;
}

/*******************************************************************************
                                 AUXILIARY FUNCTIONS
 ********************************************************************************/

void free_matrix(block** arr2D)
{
	free(arr2D[0]);
    free(arr2D);
}

block** malloc_matrix(){
    //Alloco le matrici
    block** matrix = calloc (size, sizeof(block*));
    matrix[0] = calloc (size * size, sizeof(block));

    for (int i = 1; i < size; i++) {
        matrix[i] = matrix[i-1] + size;
    }

    return matrix;
}

int copy_matrix(block*** input_matrix, block*** output_matrix){
    //The matrixs should be already allocated
    block** in = *input_matrix;
    block** out = *output_matrix;

    memcpy(out[0],in[0], size*size*sizeof(block));

    return 0;
}

int read_matrix_from_file(block*** m, char* file_name){
	int i;
    int j;
    block** matrix = *m;

    matrix = malloc_matrix();

	if(!rank){
		FILE *file;
		file=fopen(file_name, "r");
		if(file == NULL){
			printf("\nError while opening the file. \n");
      MPI_Abort(MPI_COMM_WORLD, -1);
			return -1;
		}

		for(i = 0; i < size; i++){
			for(j = 0; j < size; j++) {
				if (!fscanf(file, "%hi", &(matrix[i][j].value)))
					break;// mat[i][j] -= '0';
				matrix[i][j].state = 'u';
			}

		  }

	}
  *m = matrix;
  return 0;
}

int read_matrix_from_stdin(block*** m){
    int i;
    int j;
    block** matrix = *m;

    matrix = malloc_matrix();

	if(!rank){
		for(i = 0; i < size; i++){
			for(j = 0; j < size; j++) {
				if (!scanf("%hi", &(matrix[i][j].value)))
					break;// mat[i][j] -= '0';
				matrix[i][j].state = 'u';
			}

		  }
	}
  *m = matrix;
  return 0;
}

int print_matrix(block** matrix){
    //printf("The matrix is: \n");
    for(int i = 0; i < size; i++) {
        for(int j = 0; j < size; j++) {
            if(matrix[i][j].state == 'w')
                printf("%d\t", matrix[i][j].value);
            else if(matrix[i][j].state == 'b')
                printf("x\t");
            else
                printf("%d,%c\t", matrix[i][j].value, matrix[i][j].state);
        }
    printf("\n");
    }
    return 1;
}

int print_matrix_result(block** matrix){
    for(int i = 0; i < size; i++) {
        for(int j = 0; j < size; j++) {
             if ( matrix[i][j].state == 'b' ) printf("x\t");
			       else printf("%i\t", matrix[i][j].value);
        }
    printf("\n");
    }
    return 1;
}

double logbase(int base, int x) {
    return log((double) x) / log((double)base);
}

int isNumber(char number[])
{
    int i = 0;

    //checking for negative numbers
    if (number[0] == '-')
        i = 1;
    for (; number[i] != 0; i++)
    {
        //if (number[i] > '9' || number[i] < '0')
        if (!isdigit(number[i]))
            return 0;
    }
    return 1;
}

/*******************************************************************************
                                PARALLEL FUNCTIONS
********************************************************************************/

int scatterElements(block*** matrix, block*** my_elements,int my_size){
  //printf("%d) %d %d %d \n",rank, my_size,sendcounts[rank], displs[rank] );
  MPI_Scatterv(**matrix, sendcounts, displs, block_datatype, **my_elements, size*my_size,block_datatype, 0, MPI_COMM_WORLD);
  //printf("%s\n", "Dopo scatter");
  return 0;
}

int scatterElementsColumn(block*** matrix, block*** my_elements,int my_size){
  printf("%d) %d %d %d \n",rank, my_size,sendcountsCol[rank], displsCol[rank] );
  MPI_Scatterv(**matrix, sendcountsCol, displsCol, ncol_type, **my_elements, my_size , ncol_type, 0, MPI_COMM_WORLD);
  //printf("%s\n", "Dopo scatterCol");
  return 0;
}

int gatherUnknown(int* unknowns,int unknown){
  //printf("%d) %s\n", rank, "Prima gatherUnknown");
  return MPI_Gather(&unknown, 1, MPI_INT,unknowns, 1, MPI_INT, 0, MPI_COMM_WORLD);
}

int print_sub_matrix(block** matrix,int my_size){
  for(int i = 0; i < my_size; i++) {
      for(int j = 0; j < size; j++) {
          if(matrix[i][j].state == 'w')
              printf("%d\t", matrix[i][j].value);
          else if(matrix[i][j].state == 'b')
              printf("x\t");
          else
              printf("%d,%c\t", matrix[i][j].value, matrix[i][j].state);
      }
  printf("\n");
  }

  printf("\n \n \n");
  return 0;
}

int print_once(int prev,int next, block** mypartition){
	int stampare = 0;
	MPI_Status recv_status;

	if(rank == n_process - 1 ){
		if ( MPI_Recv(&stampare, 1, MPI_INT, prev, 0, MPI_COMM_WORLD, &recv_status) == -1){
			printf("MPI_Recv failed.\n");
			return -1;
		}
		print_sub_matrix(mypartition, sendcounts[rank]/size);
	}else if(rank == 0){
		print_sub_matrix(mypartition, sendcounts[rank]/size);
		stampare = 1;

		if ( MPI_Send(&stampare, 1, MPI_INT, next, 0, MPI_COMM_WORLD) == -1){
            printf("MPI_Send failed.\n");
            return -1;
        }
	}else{
		if ( MPI_Recv(&stampare, 1, MPI_INT, prev, 0, MPI_COMM_WORLD, &recv_status) == -1){
			printf("MPI_Recv failed.\n");
			return -1;
		}

		print_sub_matrix(mypartition, sendcounts[rank]/size);

		if ( MPI_Send(&stampare, 1, MPI_INT, next, 0, MPI_COMM_WORLD) == -1){
            printf("MPI_Send failed.\n");
            return -1;
        }
	}

	return 0;
}

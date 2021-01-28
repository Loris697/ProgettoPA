#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "mpi.h"
#include "math.h"
#include <sys/time.h>
#include <stdint.h>
#include <ctype.h>
#include <unistd.h>

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
MPI_Datatype block_datatype, block_n_element_col_r;
MPI_Datatype col_datatype;
MPI_Datatype ncol_type;
MPI_Datatype col_datatype_block, col_datatype_block1, col_datatype_block_r;
MPI_Datatype block_n_element_col;
MPI_Datatype block_n1_element_col;
int* sendcounts;
int* displs;
int* sendcountsCol;
int* displsCol;
int my_size, my_n_col, my_n_row;


int read_matrix_from_file(block*** matrix, char* file);
int read_matrix_from_stdin(block*** m);
int print_matrix(block** matrix);
int print_matrix_result(block** matrix);
int solve_hitori(block** matrix,int i,int unknown);
int check_row_and_column(block*** m, int* unknown);
int find_connection_of_white(int row,int col,block** matrix);
int check_adjacent_rules(block*** m, int *unknown);
int check_neighborhood(int row,int col,block** matrix, int** flag_vector);
//FUNZIONI USILIARIE
void free_matrix(block** arr2D);
block** malloc_matrix();
int copy_matrix(block*** input_matrix, block*** output_matrix);
int apply_rule(block*** matrix, int unknown);
double logbase(int base, int x);
int isNumber(char number[]);
int check_input(char** input_file, int argc, char** argv);
//PARALLEL FUNCTION
int scatterElements(block*** matrix, block*** my_elements,int my_size);
int scatterElementsColumn(block*** matrix, block*** my_elements,int my_size);
int check_row_and_column_par(block*** m, int unknown);
int gatherUnknown(int* unknowns,int unknown);
int gatherElements(block*** matrix, block*** my_elements);
int gatherElementsCol(block*** matrix, block*** my_elements);
int Bcast(int message);
void init_input_data();
void init_datatype();
void halo_exchange(block **my_block);
int gatherBlocks(block ***matrix, block ***my_block);
int checkUnknown(int *unknowns, int* unknown);
//SEQ FUNCTION
int check_adjacent_rules_seq(block*** m, int unknown);
int check_row_and_column_seq(block*** m, int unknown);

int main(int argc, char **argv){
    //Initialization of the MPI data structures
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &n_process);

    int radp = sqrt(n_process);
    int rank_y = rank / radp;
    int rank_x = rank % radp;
    //Solo per rimuovere il warning
    char* input_file = NULL;

    //Controllo sull'input
    if ( check_input(&input_file, argc, argv) == -1) {
      MPI_Finalize();
      return 0;
    }

    my_size = ((rank+1) * size)/n_process - (rank * size)/n_process;
    my_n_col = ((rank_x+1) * size)/radp - (rank_x * size)/radp + 2;
    my_n_row = ((rank_y+1) * size)/radp - (rank_y * size)/radp + 2;
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
      displs[i] =( i * size)/n_process * size;
      displsCol[i] =( i * size)/n_process;
    }

    init_datatype();

		//Per misurare il tempo
		MPI_Barrier(MPI_COMM_WORLD);
		double time = - MPI_Wtime();

    if(rank == 0){
      if(!solve_hitori(global_matrix,0,unknown)){
          printf("Solution found.\n");
          print_matrix(global_matrix);
          fprintf(stderr,"Examinated %d node by process %d.\n", node_count, rank);
      }else{
          printf("No solution.\n");
          fprintf(stderr,"Examinated %d node by process %d.\n", node_count, rank);
      }
      Bcast(2);
    }else{
      int message = 1;
      int failed = 0;

      //Vettore per per righe e per le colonne
      block** my_elements = calloc (size, sizeof(block*));
      my_elements[0] = calloc (size * my_size, sizeof(block));

      for (int i = 1; i < my_size; i++) my_elements[i] = my_elements[i-1] + size;

      while(1){
        unknown = 0;

        //Il processo 0 mi dice se devo continuare o meno
        if ( MPI_Bcast(&message,1, MPI_INT,0,MPI_COMM_WORLD) == -1){
            printf("MPI_Bcast failed.\n");
            MPI_Abort(MPI_COMM_WORLD, -1);
            return -1;
        }

        if( message == 0 ){
            //Per la verifica della regola su righe e colonne
            scatterElements(&my_elements, &my_elements,my_size);

            if( check_row_and_column(&my_elements, &unknown) ){
              printf("%d)Regole righe violate\n", rank);
              failed = 1;
            }

            //Mando indietro i dati al processo zero
            gatherElements(&my_elements, &my_elements);

            scatterElementsColumn(&my_elements,&my_elements,my_size);

            if( check_row_and_column(&my_elements, &unknown) ) {
              printf("%d)Regole colonne violate\n", rank);
              failed = 1;
            }

            gatherElementsCol(&my_elements, &my_elements);

            //Invio il valore degli unknown restanti al processo zero
            if (failed) gatherUnknown(NULL, 1);
            else gatherUnknown(NULL, unknown);
        }else if ( message == 1 ) {
          //Per la verifica delle regole di adiacenza
          block** my_block;
          init_input_data(NULL, &my_block, &my_n_row, &my_n_col);

          halo_exchange(my_block);

          if( check_adjacent_rules(&my_block,&unknown) ) unknown = 1;

          gatherBlocks( NULL, &my_block);
          free(my_block);

          gatherUnknown(NULL, unknown);
        }
        else break;
      }
      free_matrix(my_elements);
    }

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
  		int return_code = solve_hitori(backup_matrix, i + 1, unknown);
  		if( return_code == 0) return 0;
    }

    node_count++;
    matrix[r][c].state = 'b';
    unknown = unknown_copy;

    //erase wrong part
    free_matrix(backup_matrix);

    unknown = apply_rule(&matrix, unknown);
	if( unknown != -1){
		int return_code = solve_hitori(matrix, i + 1, unknown);
		if( return_code == 0) return 0;
  }

  return -1;
}

int check_row_and_column(block*** m, int* u){
    block** matrix = *m;
    int unknown = *u;

    for(int i = 0; i < sendcounts[rank]; i++){
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

int check_adjacent_rules(block*** m, int *u){
    //Find if there are two black cell adjacent
    block** matrix = *m;
    int unknown = *u;

    for(int i = 0; i < my_n_col * my_n_row; i++){
            //My current block's row and colum
            int r = i / my_n_col;
            int c = i % my_n_col;

            //Nei bordi non devo fare niente
            if(r == 0 || r == my_n_row - 1 || c == 0 || c == my_n_col - 1 ) continue;

            if (r > 0 ){
                if (matrix[r-1][c].state ==  'b' && matrix[r][c].state ==  'b') return 1;
                else if(matrix[r-1][c].state ==  'b' && matrix[r][c].state ==  'u'){
                    matrix[r][c].state =  'w';
                    unknown = unknown - 1;
                }
            }if(c > 0){
                if (matrix[r][c-1].state ==  'b' && matrix[r][c].state ==  'b' ) return 1;
                else if(matrix[r][c-1].state ==  'b' && matrix[r][c].state ==  'u'){
                    matrix[r][c].state =  'w';
                    unknown = unknown - 1;
                }
            }if(r + 1 < my_n_row){
                if (matrix[r+1][c].state ==  'b' && matrix[r][c].state ==  'b' ) return 1;
                else if(matrix[r+1][c].state ==  'b' && matrix[r][c].state ==  'u'){
                    matrix[r][c].state =  'w';
                    unknown = unknown - 1;
                }
            }if(c + 1 < my_n_row){
                if (matrix[r][c+1].state ==  'b' && matrix[r][c].state ==  'b') return 1;
                else if(matrix[r][c+1].state ==  'b' && matrix[r][c].state ==  'u'){
                    matrix[r][c].state =  'w';
                    unknown = unknown - 1;
                }
          }
    }

    *u = unknown;
    return 0;
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
  if( n_process == 1 ) return check_row_and_column_seq(matrix,unknown);
  int* unknowns = (int *) malloc(n_process*sizeof(int));
  block** my_elements = calloc (size, sizeof(block*));
  my_elements[0] = calloc (size * my_size, sizeof(block));

  for (int i = 1; i < my_size; i++) my_elements[i] = my_elements[i-1] + size;

  scatterElements(matrix,&my_elements, my_size);

  if( check_row_and_column(&my_elements, &unknown) ) {
    printf("%s\n", "Regole non rispettate");
    unknown = -1;
  }

  gatherElements(matrix, &my_elements);

  scatterElementsColumn(matrix,&my_elements, my_size);

  if( check_row_and_column(&my_elements, &unknown) ) {
      //La gather serve per sbloccare gli altri processi dall'attesa
      printf("%s\n", "Regole non rispettate");
      unknown = -1;
  }

  gatherElementsCol(matrix, &my_elements);

  gatherUnknown(unknowns, unknown );

  free(my_elements);
  return checkUnknown(unknowns, &unknown);
}

int check_adjacent_rules_par(block*** matrix,int unknown ){
  if( n_process == 1 ) return check_adjacent_rules_seq(matrix,unknown);
  block** my_block;
  int* unknowns = (int *) malloc(n_process*sizeof(int));
  int r, c;

  init_input_data(matrix, &my_block, &my_n_col, &my_n_row);

  halo_exchange(my_block);

  if( check_adjacent_rules(&my_block,&unknown) ){
    printf("%s\n", "Regole non rispettate");
    unknown = -1;
  }

  //copio nel nuovo vettore i dati
  for (int i = 1; i < my_n_row - 1; i++) {
      memcpy((*matrix)[i-1], &my_block[i][1], (my_n_col - 2)*sizeof(block));
  }

  gatherBlocks(matrix, &my_block);

  gatherUnknown(unknowns, unknown );

  free(my_block);
  for(int i = 0; i < size * size; i++){
    r = i / size;
    c = i % size;
    if (((*matrix)[r][c]).state == 'w') break;
  }

  if(find_connection_of_white(r, c, *matrix)) return -1;
  return checkUnknown(unknowns, &unknown);
}

int apply_rule(block*** matrix,int unknown){
    int prev_u;
    do{
        prev_u = unknown;
        //int unknown_real = 0;
        Bcast(0);
        //printf("%d)Inizio( verifica di regole righe/colonne.\n", rank);
        unknown = check_row_and_column_par(matrix, unknown);
        if (unknown < 0) return -1;
        Bcast(1);
        //printf("%d)Inizio verifica di regole adiacenza.\n", rank);
		    unknown = check_adjacent_rules_par(matrix,unknown);
        if (unknown < 0) return -1;
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
    //printf("The matrix is: %d\n", size);
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
*******************************************************************************/
void init_datatype(){
  //Making the block datatype
  MPI_Aint displ[2];
  int bl[2] = {1, 1};
  MPI_Datatype dt[2] = {MPI_INT, MPI_CHAR};

  MPI_Get_address(&(global_matrix[0][0].value), &displ[0]);
  MPI_Get_address(&(global_matrix[0][0].state), &displ[1]);

  //printf("%d)Creo il tipo di dato. \n", rank);

  displ[1] -= displ[0];
  displ[0] = 0;

  //Questo tipo di dato serve per mandare un singolo blocco
  if (MPI_Type_create_struct(2, bl, displ, dt, &block_datatype) == -1){
      printf("Creation of the datatype failed.\n");
  }
  MPI_Type_commit(&block_datatype);

  //Questo tipo di dato permette di mandare solo una colonna alla volta
  //Serve solo per costruire il tipo ncol_type
  MPI_Type_vector(size, 1, size, block_datatype, &col_datatype);
  MPI_Type_commit(&col_datatype);

  //Stesso ragionamento di prima ma questi du servono per lo scambio degli halo e dei blocchi
  MPI_Type_vector(my_n_row-2, 1, size, block_datatype, &col_datatype_block);
  MPI_Type_commit(&col_datatype_block);

  MPI_Type_vector(my_n_row-1, 1, size, block_datatype, &col_datatype_block1);
  MPI_Type_commit(&col_datatype_block1);

  //Questo servirà per la ricezione invece
  MPI_Type_vector(my_n_row-2, 1, my_n_col, block_datatype, &col_datatype_block_r);
  MPI_Type_commit(&col_datatype_block_r);


  //Per lo scambio efffettivo delle colonne
  MPI_Type_create_resized(col_datatype, 0 , (MPI_Aint) sizeof(block),  &ncol_type);
  MPI_Type_commit(&ncol_type);

	MPI_Type_create_resized(col_datatype_block, 0 , (MPI_Aint) sizeof(block),  &block_n_element_col);
  MPI_Type_create_resized(col_datatype_block1, 0 , (MPI_Aint) sizeof(block),  &block_n1_element_col);
  MPI_Type_create_resized(col_datatype_block_r, 0 , (MPI_Aint) sizeof(block),  &block_n_element_col_r);

	MPI_Type_commit(&block_n_element_col);
  MPI_Type_commit(&block_n1_element_col);
  MPI_Type_commit(&block_n_element_col_r);
}

int scatterElements(block*** matrix, block*** my_elements,int my_size){
  MPI_Scatterv(**matrix, sendcounts, displs, block_datatype, **my_elements, size*my_size,block_datatype, 0, MPI_COMM_WORLD);
  return 0;
}

int scatterElementsColumn(block*** matrix, block*** my_elements,int my_size){
  MPI_Scatterv(**matrix, sendcountsCol, displsCol, ncol_type, **my_elements, size*my_size , block_datatype, 0, MPI_COMM_WORLD);
  return 0;
}

int gatherUnknown(int* unknowns,int unknown){
  MPI_Gather(&unknown, 1, MPI_INT,unknowns, 1, MPI_INT, 0, MPI_COMM_WORLD);
  return 0;
}

int gatherElements(block*** matrix, block*** my_elements){
  MPI_Gatherv(**my_elements, sendcounts[rank], block_datatype, **matrix, sendcounts, displs, block_datatype, 0, MPI_COMM_WORLD);
  return 0;
}

int gatherElementsCol(block*** matrix, block*** my_elements){
  MPI_Gatherv(**my_elements, sendcounts[rank], block_datatype, **matrix, sendcountsCol, displsCol , ncol_type, 0, MPI_COMM_WORLD);
  return 0;
}

int gatherBlocks(block*** matrix, block*** my_block){
  MPI_Status st;
  int radp = sqrt(n_process);
  int rank_x, rank_y;

  rank_y = rank / radp;
  rank_x = rank % radp;

  //each process computes dim rows
  int row = ((rank_y+1) * size)/radp - (rank_y * size)/radp + 2;
  int col = ((rank_x+1) * size)/radp - (rank_x * size)/radp + 2;
  int* nrow = &row;
  int* ncol = &col;

  if (n_process == 1) *nrow = *ncol = size;

  if (!rank){
      for (int i = 1; i < n_process; i++) {
        int rank_x, rank_y;
        int start_row, start_col;
        int num_row, num_col;

        rank_y = i / radp;
        rank_x = i % radp;

        start_row = (rank_y * size)/radp;
        start_col = (rank_x * size)/radp;

        num_row = ((rank_y + 1) * (size))/radp - start_row;
        num_col = ((rank_x + 1) * (size))/radp - start_col;

        if ((num_row) + 2 == *nrow)
          MPI_Recv((*matrix)[start_row] + start_col, num_col, block_n_element_col, i, 0, MPI_COMM_WORLD, &st);
        else
          MPI_Recv((*matrix)[start_row] + start_col, num_col, block_n1_element_col, i, 0, MPI_COMM_WORLD, &st);
      }
  }
  else{
      MPI_Send(&((*my_block)[1][1]), (*ncol)-2, block_n_element_col_r, 0, 0, MPI_COMM_WORLD);
  }
  return 0;
}

int Bcast(int message){
  if ( MPI_Bcast(&message,1, MPI_INT,0,MPI_COMM_WORLD) == -1){
      printf("MPI_Bcast failed.\n");
      MPI_Abort(MPI_COMM_WORLD, -1);
      return -1;
  }
  return 0;
}

void init_input_data(block*** matrix, block*** my_block, int* nrow, int* ncol)
{
    MPI_Status st;
    int radp = sqrt(n_process);
    int rank_x, rank_y;

    rank_y = rank / radp;
    rank_x = rank % radp;

    //each process computes dim rows
    *nrow = ((rank_y+1) * size)/radp - (rank_y * size)/radp + 2;
    *ncol = ((rank_x+1) * size)/radp - (rank_x * size)/radp + 2;

    if (n_process == 1) *nrow = *ncol = size;

    *my_block = calloc ((*nrow ) + 2, sizeof(block*));
    (*my_block)[0] = calloc ((((*nrow ) + 2) * (*ncol ) + 2), sizeof(block));

    for (int i = 1; i < (*nrow ); i++) {
        (*my_block)[i] = (*my_block)[i-1] + (*ncol );
    }

    if (!rank){
        //copio nel nuovo vettore i dati
        for (int i = 1; i < (*nrow ) - 1; i++) {
            memcpy(&(*my_block)[i][1], (*matrix)[i-1], ((*ncol ) - 2)*sizeof(block));
        }

        for (int i = 1; i < n_process; i++) {
        	int rank_x, rank_y;
        	int start_row, start_col;
          int num_row, num_col;

          rank_y = i / radp;
          rank_x = i % radp;

          start_row = (rank_y * size)/radp;
          start_col = (rank_x * size)/radp;

          num_row = ((rank_y + 1) * (size))/radp - start_row;
          num_col = ((rank_x + 1) * (size))/radp - start_col;

          if ((num_row  + 2) == *nrow)
          	MPI_Send((*matrix)[start_row] + start_col, num_col, block_n_element_col, i, 0, MPI_COMM_WORLD);
          else
          	MPI_Send((*matrix)[start_row] + start_col, num_col, block_n1_element_col, i, 0, MPI_COMM_WORLD);
        }
    }
    else{
        MPI_Recv(&((*my_block)[1][1]), (*ncol-2), block_n_element_col_r, 0, 0, MPI_COMM_WORLD, &st);
    }
}

void halo_exchange(block **my_block)
{
    int size;
    int rank;
    MPI_Status st;
    int rank_x, rank_y;
    int radp;
    int rank_east, rank_west, rank_north, rank_south;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

  	radp = sqrt(n_process);
  	rank_y = rank / radp;
  	rank_x = rank % radp;

  	rank_east = rank_y * radp + rank_x + 1;
  	rank_west = rank_y * radp + rank_x - 1;
  	rank_north = (rank_y - 1) * radp + rank_x;
  	rank_south = (rank_y + 1) * radp + rank_x;

    if (size == 1) return;

    if (rank_x != 0) //Sendrecv towards west
        MPI_Sendrecv(my_block[1]+1, 1, col_datatype_block_r, rank_west, 0, my_block[1], 1, col_datatype_block_r, rank_west, 0, MPI_COMM_WORLD, &st);

    if (rank_x != radp-1) //Sendrecv towards east
        MPI_Sendrecv(my_block[1]+my_n_col-2, 1, col_datatype_block_r, rank_east, 0, my_block[1]+my_n_col-1, 1, col_datatype_block_r, rank_east, 0, MPI_COMM_WORLD, &st);

    if (rank_y != 0) //Sendrecv towards north
        MPI_Sendrecv(my_block[1]+1, my_n_col-2, block_datatype, rank_north, 0, my_block[0]+1, my_n_col-2, block_datatype, rank_north, 0, MPI_COMM_WORLD, &st);

    if (rank_y != radp-1) //Sendrecv towards south
        MPI_Sendrecv(my_block[my_n_row-2]+1, my_n_col-2, block_datatype, rank_south, 0, my_block[my_n_row-1]+1, my_n_col-2, block_datatype, rank_south, 0, MPI_COMM_WORLD, &st);

}

int checkUnknown(int* unknowns, int* unknown){
  //Funzione richiamata dal processo uno per controllare e aggiornare il numero delle celle sconosciute
  for(int i = 1; i < n_process; i++){
    if(unknowns[i] == 1) {
      free(unknowns);
      printf("Returning -1\n");
      return -1;
    }
    *unknown = *unknown + unknowns[i];
  }
  free(unknowns);
  return *unknown;
}

/******************************************************************************
                          SEQ FUNCTION
********************************************************************************/
int check_row_and_column_seq(block*** m, int unknown){
	block** matrix = *m;

	for(int i = 0; i < size * size; i++){
		//My current block's row and colum
		int r = i / size;
		int c = i % size;

		//if the current block's state is unknown i can't do anything
		if(matrix[r][c].state != 'u'){
			for(int rnew = 0; rnew < size; rnew++ )
				//I have to not check the block with itself
				//It's varing 'r' so i'm verifing the column rule
				if(rnew != r){
					if(matrix[r][c].value == matrix[rnew][c].value){
						//If the number is the same
						//The same white state means wrong solution
						if(matrix[r][c].state == 'w' && matrix[rnew][c].state == 'w'){
							return -1;
					}
					}
				}

			//Let's do the same thing with for the rows
			for(int cnew = 0; cnew < size; cnew++ ){
				if(cnew != c)
					if(matrix[r][c].value == matrix[r][cnew].value){
						if(matrix[r][c].state == 'w' && matrix[r][cnew].state == 'w'){
							return -1;
						}
				}

			}
		}
	}
	return unknown;
}

int check_adjacent_rules_seq(block*** m, int unknown){
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
				}if(c > 0){
					if (matrix[r][c-1].state ==  'b')
						return -1;
				}if(r + 1 < size){
					if (matrix[r+1][c].state ==  'b')
						return -1;
				}if(c + 1 < size){
                    if (matrix[r][c+1].state ==  'b')
                        return -1;
				}

			}else if (first_white) {
				first_white--;
				if(find_connection_of_white(r, c, matrix))
					return -1;
			}
	}

	return unknown;
}

int check_input(char** input_file, int argc, char** argv){
  int radp = sqrt(n_process);
  size = 0;
  int c;

  opterr = 0;

  while ((c = getopt (argc, argv, "f:s:")) != -1)
   switch (c)
     {
     case 'f':
       *input_file = (char*) malloc(strlen(optarg)*sizeof(char));
       strcpy(*input_file,optarg);
       break;
     case 's':
       if( isNumber(optarg )) size = atol(optarg);
       else{
         if(!rank) printf("Il parametro passato con -s è la dimensione del puzzle, deve essere un intero.\n");
         return -1;
       }
       break;
     default:
       if(!rank) printf("Paramentro sconosciuto.\n");
     }


  if( *input_file == NULL){
      //Setto il valore di defult a stdin
      if(!rank) printf("Il parametro passato con -f è il file da dove prendere la matrice. Di default verrà aspettata la matrice da stdin\n");
      *input_file = (char*) malloc(6*sizeof(char));
      strcpy(*input_file,"stdin");
  }

  if (size < 1 ){
    if(!rank) printf("Il parametro passato con -s è la dimensione del puzzle, è obbligatorio e deve essere un intero maggiore di 0.\n");
    return -1;
  }

  if(!strcmp(*input_file ,"stdin"))
      read_matrix_from_stdin(&global_matrix);
  else if( read_matrix_from_file(&global_matrix,*input_file) == -1) {
      if (!rank) printf("Impossibile leggere dal file.\n");
      return -1;
  }

  if ((float) radp != sqrt(n_process)){
    if (!rank) printf("Il numero di processi deve essere un quadrato perfetto.\n");
    return -1;
  }

  return 0;
}

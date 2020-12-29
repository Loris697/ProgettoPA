#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "mpi.h"
#include "math.h"
#include <sys/time.h>
#include <stdint.h>
#include <ctype.h>

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

struct message_status {
    int count;
    int color;
};

//size of the matrix
int size;
block** global_matrix;
//count the node examinated
int node_count = 0;
int rank;
//Total number of process
int n_process = 0;
MPI_Datatype message_datatype;
MPI_Datatype block_datatype;
struct message_status status;
//Init the value of each process value
int message_count = 0;
//Usata dal processo 0 per contare quanti thread hanno finito
int no_success_process= 0;
//Usata dal processo per segnare quando lui ha finito;
int i_terminated = 0;
int return_code = 0;
int n_cell_assigned = 0;
int termination_started = 0;

int read_matrix_from_file(block*** matrix, char* file);
int read_matrix_from_stdin(block*** m);
int print_matrix(block** matrix);
int print_matrix_result(block** matrix);
int solve_hitori(block** matrix,int i,int unknown);
int check_row_and_column(block*** m, int unknown);
int find_connection_of_white(int row,int col,block** matrix);
int check_adjacent_rules(block*** m, int unknown);
int check_neighborhood(int row,int col,block** matrix, int** flag_vector);
void free_matrix(block** arr2D);
block** malloc_matrix();
int copy_matrix(block*** input_matrix, block*** output_matrix);
int apply_rule(block*** matrix, int unknown);
double logbase(int base, int x);
//COMMUNICATION FUNCTIONS
int read_message(int prev, int* color);
int send_message(int next, int* color);
int check_for_termination(int next,int prev);
int check_for_termination_waiting(int next,int prev);
int ring(int next,int prev, int* color);
int ring_waiting(int next,int prev, int* color);
int start_termination(int next,int prev);
int increasing_no_success_proc(int idproc);
int print_once(int prev,int next);
void shuffle(int *array, size_t n);
int isNumber(char number[]);

int main(int argc, char **argv)
    {
        char* input_file;
        //count how many block I've to find state
        int unknown;

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

        //Initialization of the MPI data structures
        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &n_process);

		//Making the logical ring to avoid bad termination due to the call of MPI_finalize
        //with message in the receive buffer
        int next = (rank + 1)%n_process;
        int prev = (rank + n_process - 1)%n_process;
        int color = 0; // 0 = white 1 = black

        //Making the datatype
        if( MPI_Type_contiguous(2, MPI_INT, &message_datatype) == -1){
            printf("Creation of the datatype failed.\n");
            return -1;
        }
        MPI_Type_commit(&message_datatype);

        input_file = argv[1];
        unknown = size*size;
        block** backup_matrix = malloc_matrix();

        MPI_Aint displ[2];
      	int bl[2] = {1, 1};
      	MPI_Datatype dt[2] = {MPI_INT, MPI_CHAR};

      	MPI_Get_address(&(backup_matrix[0][0].value), &displ[0]);
      	MPI_Get_address(&(backup_matrix[0][0].state), &displ[1]);

      	displ[1] -= displ[0];
      	displ[0] = 0;

      	if (MPI_Type_create_struct(2, bl, displ, dt, &block_datatype) == -1){
      		printf("Creation of the datatype failed.\n");
      		return -1;
      	}
      	MPI_Type_commit(&block_datatype);

        if(!strcmp(input_file ,"stdin"))
      	read_matrix_from_stdin(&global_matrix);
        else if( read_matrix_from_file(&global_matrix,input_file) == -1) {
      		MPI_Finalize();
      		return -1;
      	}

        //if (!rank) printf("Solving\n");
		//Per misurare il tempo
		MPI_Barrier(MPI_COMM_WORLD);
		double time = - MPI_Wtime();
    n_cell_assigned = (int) ceil(logbase((double)2,(double) n_process * PART_PER_PROCESS));
  	int parts = (int) pow(2, n_cell_assigned);
  	int ppp_effective = parts / n_process;
  	if (rank < parts % n_process) ppp_effective++;
  	int i;
  	int return_code = 1;
  	int k = 0;
  	//printf("%d) PPP = %d\n", rank, ppp_effective);
  	int *mypart = (int *) malloc(sizeof(int)*(ppp_effective));
  	for (i = rank; i < parts; i = i + n_process){
  		mypart[k] = i;
  		k++;
  	}

  	//if(rank == 0) printf("%d) P Total = %d\n", rank, parts);

  	//shuffle(mypart, ppp_effective);

  	for (k = 0; k < ppp_effective; k++){
  		int j = 0;
  		//printf("%d) k = %d, mypart[k] = %d\n", rank, k, mypart[k]);
  		copy_matrix(&global_matrix, &backup_matrix);

      int i = mypart[k];
      for(j = 0; j < n_cell_assigned; j++){
          int r = j / size;
          int c = j % size;
          if( i & (int) pow(2, j)){
            backup_matrix[r][c].state = 'b';
            //Se ci sono due neri adiacenti
            if( i & (int) pow(2, j + 1)) break;
          }else backup_matrix[r][c].state = 'w';
        }

        //Se il ciclo non si è concluso
        if(j != n_cell_assigned) continue;

        //print_matrix(backup_matrix);

        return_code = solve_hitori(backup_matrix,0,unknown-n_cell_assigned);
        if(!return_code){
          //printf("Solution found.\n");
          i_terminated = 1;
          fprintf(stderr,"Examinated %d node by process %d.\n", node_count, rank);
          break;
        }else if(return_code == -2) break;

      }
      if(rank){
        int terminate = TERMINATE_NO_SUCC;
        if ( return_code != -2 ){
          //Se sono uscito dal ciclo perchè ho trocato una soluzione
          if(k < ppp_effective) send_message(0, &color);
          //Se sono uscito dal ciclo perchè ho finito le mie parti
          else  send_message(0, &terminate);
          ring_waiting(next, prev, &color);
        }
      }else{
        if ( return_code != -2 ){
          //Se sono uscito dal ciclo perchè ho trocato una soluzione
          if(k < ppp_effective) start_termination(next,prev);
          //Se sono uscito dal ciclo perchè ho finito le mie parti
          else{
            increasing_no_success_proc(rank);
            check_for_termination_waiting(next,prev);
          }
        }
      }

      //printf("%d) Terminate.\n", rank);
      fflush(stdout);
      MPI_Barrier(MPI_COMM_WORLD);
      time += MPI_Wtime();

      //per stampare una sola matrice
      if(n_process != 1) print_once(prev, next);
      else print_matrix_result(global_matrix);

      if (rank == 0) fprintf(stderr, "dim = %d)Execution time: %f\n",size, time);
      free_matrix(backup_matrix);
      MPI_Finalize();
      return 0;
}

int solve_hitori(block** matrix,int i,int unknown){
    //My block's row and colum
    int r = i / size;
    int c = i % size;
    //If the white ipotesis is not correct
    block** backup_matrix = malloc_matrix();
    int unknown_copy;
    int next = (rank + 1)%n_process;
    int prev = (rank + n_process - 1)%n_process;
    int color = 0; // 0 = white 1 = black

	//Prima di iniziare i calcoli il processo zero controlla se ci sono messaggi di terminazione
	//-2 significa terminazione immediata
    if(rank == 0){
        if( check_for_termination(next, prev)) return -2;
	  }else{
        if( ring(next, prev, &color)) return -2;
    }

    //If I have 0 unknown I have solved the puzzle :)
    if(unknown == 0){
        //printf("Returnung to main.\n" );
        global_matrix = matrix;
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
  		else if(return_code == -2) return -2;
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
		else if(return_code == -2) return -2;
    }
    return -1;
}

int check_row_and_column(block*** m, int unknown){
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
                        }else if(matrix[r][c].state == 'w' && matrix[rnew][c].state == 'u'){
                            unknown = unknown - 1;
                            matrix[rnew][c].state = 'b';
                        }
                    }
                }

                    //Let's do the same thing with for the rows
                    for(int cnew = 0; cnew < size; cnew++ ){
                        if(cnew != c)
                            if(matrix[r][c].value == matrix[r][cnew].value){
                                if(matrix[r][c].state == 'w' && matrix[r][cnew].state == 'w'){
                                    return -1;
                                }else if(matrix[r][c].state == 'w' && matrix[r][cnew].state == 'u'){
                                    unknown = unknown - 1;
                                    matrix[r][cnew].state = 'b';
                                }
                            }
                    }

        }
    }
    return unknown;
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
            return -1;
        }
    }
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

int apply_rule(block*** matrix,int unknown){
    int prev_u;
    do{
        prev_u = unknown;
        unknown = check_row_and_column(matrix, unknown);
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

	//Se non sono solo condivido la matrice trasmetto(ricevo) i dati
	if(n_process != 1) MPI_Bcast(matrix[0],size*size, block_datatype,0,MPI_COMM_WORLD);
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

	//Se non sono solo condivido la matrice trasmetto(ricevo) i dati
	if(n_process != 1) MPI_Bcast(matrix[0],size*size, block_datatype,0,MPI_COMM_WORLD);
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


void shuffle(int *array, size_t n)
{
    struct timeval t;
    gettimeofday(&t, 0);

    srand(t.tv_usec);

    if (n > 1)
    {
        size_t i;
        for (i = 0; i < n - 1; i++)
        {
          size_t j = i + rand() / (RAND_MAX / (n - i) + 1);
          int t = array[j];
          array[j] = array[i];
          array[i] = t;
        }
    }
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
                                COMMUNICATION FUNCTIONS
********************************************************************************/

int check_for_termination(int next,int prev){
    //Flag to check incoming messages
    int not_recv = 0;
    MPI_Status st;
    MPI_Status recv_status;
    int return_code = 0;

	//Vedo se ci sono messaggi sul tag 0
    if( MPI_Iprobe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD,&not_recv, &st) == -1){
        printf("MPI_Probe failed.\n");
        return -1;
    }


    if(not_recv){
        //printf("RECEIVED SOLUTION FOUND\n");
        //Reading the message
        if ( MPI_Recv(&status, 1, message_datatype, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &recv_status) == -1){
            printf("MPI_Recv failed.\n");
            return -1;
        }
	termination_started = 1;
        start_termination(next, prev);
        return 1;
    }

    //Sul tag 1 si mandano i messaggio della terminazione senza successo, leggo tutti i messaggi
    while(1){
        if( MPI_Iprobe(MPI_ANY_SOURCE, 1, MPI_COMM_WORLD,&not_recv, &st) == -1){
            printf("MPI_Probe failed.\n");
            return -1;
        }

        if (not_recv) {
            //cprintf("Process %d received not success message.\n", rank);
            //Reading the message
            //printf("RECEIVED NO SOLUTION\n");
            if ( MPI_Recv(&status, 1, message_datatype, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &recv_status) == -1){
                printf("MPI_Recv failed.\n");
                return -1;
            }

            return_code = increasing_no_success_proc(status.color);
        }else break;
    }

    return return_code;
}

int check_for_termination_waiting(int next,int prev){
    //Flag to check incoming messages
    MPI_Status recv_status;

	while(1){
		if ( MPI_Recv(&status, 1, message_datatype, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &recv_status) == -1){
				printf("MPI_Recv failed.\n");
				return -1;
		}

		if( recv_status.MPI_TAG == 0 ){
			start_termination(next, prev);
			termination_started = 1;
			return 1;
		}else if( recv_status.MPI_TAG == 1 && !termination_started){
			//Sul tag 1 si mandano i messaggio della terminazione senza successo
			if(increasing_no_success_proc(status.color))
				return TERMINATE;
		}
	}

    return 0;
}

int start_termination(int next,int prev){
    int color = 0;
    int message_count = 0;

	//Se il programma è eseguito in sequenziale
	if(n_process == 1) return 1;

    while(1){
        send_message(next, &color);

        //Chi manda il messaggio non deve cambiare colore
        color = 0;

        read_message(prev, &color);

        if(color == 0 && status.color == 0 && message_count + status.count == 0){
            status.color = TERMINATE;
            return send_message(next, &color);
        }
    }
    return 0;
}

int read_message(int prev, int* color){
    MPI_Status recv_status;

    if ( MPI_Recv(&status, 1, message_datatype, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &recv_status) == -1){
        printf("MPI_Recv failed.\n");
        return -1;
    }

    if( status.color == TERMINATE )  return TERMINATE;

    status.count += message_count;

    if ( *color == 1)
        status.color = 1;

    *color = 0;

    return 0;
}

int send_message(int next, int* color){
    //Messaggio di terminazione senza successo
    if(*color == TERMINATE_NO_SUCC){
        //Imposto il colore uguale al mio rank per far capire al processo 0 chi ha terminato
        status.color = rank;
	//printf("Process %d send a message (%d) to %d\n", rank,status.color, next);
        if ( MPI_Send(&status, 1, message_datatype, next, 1, MPI_COMM_WORLD) == -1){
            printf("MPI_Send failed.\n");
            return -1;
        }
    }else{
    //Messaggio su tag normale
        if ( MPI_Send(&status, 1, message_datatype, next, 0, MPI_COMM_WORLD) == -1){
            printf("MPI_Send failed.\n");
            return -1;
        }
	//printf("Process %d send a message (%d) to %d\n", rank,status.color, next);
        if(status.color == TERMINATE) return TERMINATE;
    }

    return 0;
}

int ring(int next,int prev, int* color){
    int not_recv;
    MPI_Status st;
    int terminate;

    if( MPI_Iprobe(prev, 0, MPI_COMM_WORLD,&not_recv, &st) == -1){
        printf("MPI_Probe failed.\n");
        return -1;
    }

    if(not_recv)
        do{
            terminate = read_message(prev, color);
            if (rank != n_process - 1 || terminate != TERMINATE) send_message(next, color);
            if (terminate == TERMINATE) return TERMINATE;
        }while(1);

    return 0;
}

int ring_waiting(int next,int prev, int* color){
    int terminate;

    do{
        terminate = read_message(prev, color);

        if (rank != n_process - 1 || terminate != TERMINATE) send_message(next, color);
        if (terminate == TERMINATE) return TERMINATE;
    }while(1);

    return 0;
}

int increasing_no_success_proc(int idproc){
    no_success_process++;
    //Rispristino il valore zero al colore (bianco)
    status.color = 0;
    if (no_success_process == n_process){
        printf("NO SOLUTION FOUND.\n");
        start_termination(1, n_process - 1);
        return 1;
    }
    return 0;
}

int print_once(int prev,int next){
	int stampare = 0;
	MPI_Status recv_status;

  //printf("%d) terminated = %d\n", rank, i_terminated);

	if(rank == n_process - 1 ){
		if ( MPI_Recv(&stampare, 1, MPI_INT, prev, 0, MPI_COMM_WORLD, &recv_status) == -1){
			printf("MPI_Recv failed.\n");
			return -1;
		}
		if( stampare == 1 && i_terminated == 1) print_matrix_result(global_matrix);
	}else if(rank == 0){
		if (i_terminated) print_matrix_result(global_matrix);
		else stampare = 1;

		if ( MPI_Send(&stampare, 1, MPI_INT, next, 0, MPI_COMM_WORLD) == -1){
            printf("MPI_Send failed.\n");
            return -1;
        }
	}else{
		if ( MPI_Recv(&stampare, 1, MPI_INT, prev, 0, MPI_COMM_WORLD, &recv_status) == -1){
			printf("MPI_Recv failed.\n");
			return -1;
		}
		if( stampare == 1 && i_terminated == 1){
			print_matrix_result(global_matrix);
			stampare = 0;
		}
		if ( MPI_Send(&stampare, 1, MPI_INT, next, 0, MPI_COMM_WORLD) == -1){
            printf("MPI_Send failed.\n");
            return -1;
        }
	}

	return 0;
}

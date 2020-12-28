#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "mpi.h"
#include "math.h"

#define PART_PER_PROCESS 16
#define TERMINATE 31000
#define TERMINATE_NO_SUCC 31001

typedef struct {
    int value;
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
struct message_status status;
//Init the value of each process value
int message_count = 0;
//Usata dal processo 0 per contare quanti thread hanno finito
int no_success_process= 0;
//Usata dal processo per segnare quando lui ha finito;
int i_terminated = 0;
int return_code = 0;

int termination_started = 0;

int read_matrix_from_file(block*** matrix, char* file);
int read_matrix_from_stdin(block*** m);
int print_matrix(block** matrix);
int solve_hitori(block** matrix,int i,int unknown);
int check_row_and_column(block*** m, int unknown);
int find_connection_of_white(int row,int col,block** matrix);
int check_adjacent_rules(block*** m, int unknown);
int check_neighborhood(int row,int col,block** matrix, int** flag_vector);
void free_matrix(block** arr2D);
block** malloc_matrix();
int copy_matrix(block*** input_matrix, block*** output_matrix);
int apply_rule(block*** matrix, int unknown);
//COMMUNICATION FUNCTIONS
int read_message(int prev, int* color);
int send_message(int next, int* color);
int check_for_termination(int next,int prev);
int check_for_termination_waiting(int next,int prev);
int ring(int next,int prev, int* color);
int ring_waiting(int next,int prev, int* color);
int start_termination(int next,int prev);
int increasing_no_success_proc(int idproc);

int main(int argc, char **argv)
    {
        char* input_file;
        //count how many block I've to find state
        int unknown, i;

        if(argc<3 || argc>3){
            printf("\nThe program needs two arguments: \n1)The file containing the puzzle\n2) The dimension o f the puzzle");
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
        size = atoi(argv[2]);
        unknown = size*size;
        block** backup_matrix = malloc_matrix();

        if(!strcmp(input_file ,"stdin"))
            read_matrix_from_stdin(&global_matrix);
        else
            read_matrix_from_file(&global_matrix,input_file);
    
        if (!rank) printf("Solving\n");
		//Per misurare il tempo
		MPI_Barrier(MPI_COMM_WORLD);
		double time = - MPI_Wtime();
        
        //Dividing the work at least 8 part per process to better balance the workload
        //Discover how many part there are
        int n_cell_assigned = (int) ceil(log2((double) n_process * PART_PER_PROCESS));
        int parts = (int) pow(2, n_cell_assigned);
        //Il ciclo continua ad essere eseguito anche quando non ci sono più elementi da elaborare per attendere il mesaggio di terminazione sicura
        for (i = rank; i < parts; i = i + n_process){
            int j = 0;
            copy_matrix(&global_matrix, &backup_matrix);
            //Adesso rimane da codificare le possibile soluzioni, utilizzando i "parts" meno significativi della variabile "i". 0 corrisponde ad una cella bianca e 1 ad una cella nera
            for(j = 0; j < n_cell_assigned; j++){
                int r = j / size;
                int c = j % size;
                if( i & (int) pow(2, j)){
                    backup_matrix[r][c].state = 'b';
                    //Se ci sono due neri adiacenti
                    if( i & (int) pow(2, j + 1))
                        break;
                }else backup_matrix[r][c].state = 'w';
            }
            //Se il ciclo non si è concluso
            if(j != n_cell_assigned)
                continue;
				
			return_code = solve_hitori(backup_matrix,n_cell_assigned ,unknown - n_cell_assigned);
			if(!return_code){
				printf("Solution found.\n");
				//print_matrix(global_matrix);
				//Avverto il processo zero di iniziare la terminazione sicura
				fprintf(stderr,"Examinated %d node by process %d.\n", node_count, rank);
				break;
			}else if(return_code == -2){
				//printf("Exit from main loop\n");
				break;
			}
        }
		
		//printf("%d) Fuori dal ciclo. \n", rank);
        
        if(rank){
            int terminate = TERMINATE_NO_SUCC;
			if ( return_code != -2 ){
				//Se sono uscito dal ciclo perchè ho trocato una soluzione
				if(i < parts) send_message(0, &color);
				//Se sono uscito dal ciclo perchè ho finito le mie parti
				else send_message(0, &terminate);
				ring_waiting(next, prev, &color);
			}
        }else{
			if ( return_code != -2 ){
				//Se sono uscito dal ciclo perchè ho trocato una soluzione
				if(i < parts) start_termination(next,prev);
				//Se sono uscito dal ciclo perchè ho finito le mie parti
				else{
					increasing_no_success_proc(rank);
					check_for_termination_waiting(next,prev);
				}
			}
		}
		
		//printf("%d) Terminazione. \n", rank);
		
		MPI_Barrier(MPI_COMM_WORLD);
		time += MPI_Wtime();
        if (rank == 0) fprintf(stderr, "Execution time: %f\n", time);
        free_matrix(backup_matrix);
        MPI_Finalize();
        return -1;
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
        global_matrix = matrix;
        return 0;
    }


    //Find the next unknown state block
    while(matrix[r][c].state != 'u'){
        i++;
        r = i / size;
        c = i % size;
    }

    //tring to set current block to white
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
		else if(return_code == -2) return -2;
    }
    
    //failed white trying black
    node_count++;
    matrix[r][c].state = 'b';
    unknown = unknown_copy;
    
    //erase wrong part
    free_matrix(backup_matrix);

    unknown = apply_rule(&matrix, unknown);
	if( unknown != -1){
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
	//Making the datatype
	MPI_Datatype block_datatype;
	MPI_Aint displ[2];
	int bl[2] = {1, 1};
	MPI_Datatype dt[2] = {MPI_INT, MPI_CHAR};
			
	MPI_Get_address(&(matrix[0][0].value), &displ[0]);
	MPI_Get_address(&(matrix[0][0].state), &displ[1]);
			
	displ[1] -= displ[0];
	displ[0] = 0;
			
	if ( MPI_Type_struct(2, bl, displ, dt, &block_datatype) == -1){
		printf("Creation of the datatype failed.\n");
		return -1;
	}
	MPI_Type_commit(&block_datatype);

	if(!rank){
		FILE *file;
		file=fopen(file_name, "r");
		if(file == NULL)
			printf("\nError while opening the file. \n");

		for(i = 0; i < size; i++){
			for(j = 0; j < size; j++) {
				if (!fscanf(file, "%i", &(matrix[i][j].value)))
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
	//Making the datatype
	MPI_Datatype block_datatype;
	MPI_Aint displ[2];
	int bl[2] = {1, 1};
	MPI_Datatype dt[2] = {MPI_INT, MPI_CHAR};
			
	MPI_Get_address(&(matrix[0][0].value), &displ[0]);
	MPI_Get_address(&(matrix[0][0].state), &displ[1]);
			
	displ[1] -= displ[0];
	displ[0] = 0;
			
	if ( MPI_Type_struct(2, bl, displ, dt, &block_datatype) == -1){
		printf("Creation of the datatype failed.\n");
		return -1;
	}
	MPI_Type_commit(&block_datatype);

	if(!rank){
		for(i = 0; i < size; i++){
			for(j = 0; j < size; j++) {
				if (!scanf("%i", &(matrix[i][j].value)))
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
    printf("The matrix is: \n");
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

/*******************************************************************************
                                COMMUNICATION FUNCTIONS
********************************************************************************/

int check_for_termination(int next,int prev){
    //Flag to check incoming messages
    int not_recv = 0;
    MPI_Status st;
    MPI_Status recv_status;
    
    if( MPI_Iprobe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD,&not_recv, &st) == -1){
        printf("MPI_Probe failed.\n");
        return -1;
    }
    
    if(not_recv){
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
    int return_code = 0;
    while(1){
        if( MPI_Iprobe(MPI_ANY_SOURCE, 1, MPI_COMM_WORLD,&not_recv, &st) == -1){
            printf("MPI_Probe failed.\n");
            return -1;
        }
        
        if (not_recv) {
            //cprintf("Process %d received not success message.\n", rank);
            //Reading the message
            if ( MPI_Recv(&status, 1, message_datatype, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &recv_status) == -1){
                printf("MPI_Recv failed.\n");
                return -1;
            }
            
            return_code = increasing_no_success_proc(status.color);
        }else return return_code;
    }
    
    return 0;
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

    if ( MPI_Recv(&status, 1, message_datatype, prev, 0, MPI_COMM_WORLD, &recv_status) == -1){
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
        if ( MPI_Send(&status, 1, message_datatype, next, 1, MPI_COMM_WORLD) == -1){
            printf("MPI_Send failed.\n");
            return -1;
        }
		//printf("Process %d send a message (%d) to %d\n", rank,status.color, next);
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
#include "stdio.h"
#include "stdlib.h"
#include "string.h"

typedef struct {
	int value;
	char state;
	//Il valore e':
	//'b' if the cell is deleted
	//'w' if the cell cannot be deleted
	//'u' if we don't know
}block;

int size;
block** global_matrix;
//count the node examinated
int node_count = 0;

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

int main(int argc, char **argv)
	{
		char* input_file;
		//count how many block I've to find state
		int unknown;

		if(argc<3 || argc>3){
			printf("\nThe program needs two arguments: \n1)The file containing the puzzle\n2) The dimension o f the puzzle");
		}

		input_file = argv[1];
		size = atoi(argv[2]);
		unknown = size*size;

		if(!strcmp(input_file ,"stdin"))
            read_matrix_from_stdin(&global_matrix);
        else
            read_matrix_from_file(&global_matrix,input_file);
		//printf("unknown = %d\n", unknown);

		print_matrix(global_matrix);

		printf("Solving\n");

		if(!solve_hitori(global_matrix,0,unknown)){
			printf("Solution found.\n");
            printf("Examinated %d node.\n", node_count);
			print_matrix(global_matrix);
            return 0;
		}else
			printf("No solution found.\n");

		return -1;
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
    printf("Setting %d element to white. \n", i+1);
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
        printf("Applied succesfully rules. \n");
        if( solve_hitori(backup_matrix, i + 1, unknown) == 0)
            return 0;
    }
    
    //failed white trying black
    printf("Failed apply to rules. \n");
    printf("Setting %d element to black. \n", i+1);
    node_count++;
    matrix[r][c].state = 'b';
    unknown = unknown_copy;
    
    //erase wrong part
    free_matrix(backup_matrix);

    //printf("Now there are %d unkowns \n", unknown);

    unknown = apply_rule(&matrix, unknown);
    if(unknown != -1){
        printf("Applied succesfully rules. \n");
        if( solve_hitori(matrix, i + 1, unknown) == 0)
            return 0;
    }else{
        //Failed to find solution
        return -1;
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
							printf("%d) rnew = %d c = %d -returning -1\n", r, rnew, c);
							return -1;
						}else if(matrix[r][c].state == 'w' && matrix[rnew][c].state == 'u'){
							//printf("%d) rnew = %d c = %d -SET BLACK\n", r, rnew, c);
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
									printf("%d) cnew = %d r = %d -returning -1\n", c, cnew, r);
									return -1;
								}else if(matrix[r][c].state == 'w' && matrix[r][cnew].state == 'u'){
									//printf("%d) cnew = %d r = %d -SET BLACK\n", c, cnew, r);
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
				//printf("r = %d c = %d -IS BLACK\n", r, c);
				//For all the adjacent cell apply the second rule, but before I check that the cell exist
				if (r > 0 ){
					if (matrix[r-1][c].state ==  'b')
						return -1;
					else if(matrix[r-1][c].state ==  'u'){
						  //printf("r = %d c = %d -SET WHITE\n", r-1, c);
							matrix[r-1][c].state =  'w';
							unknown = unknown - 1;
						}
				}if(c > 0){
					if (matrix[r][c-1].state ==  'b')
						return -1;
					else if(matrix[r][c-1].state ==  'u'){
						  //printf("r = %d c = %d -SET WHITE\n", r, c-1);
							matrix[r][c-1].state =  'w';
							unknown = unknown - 1;
						}
				}if(r + 1 < size){
					if (matrix[r+1][c].state ==  'b')
						return -1;
					else if(matrix[r+1][c].state ==  'u'){
						//printf("r = %d c = %d -SET WHITE\n", r+1, c);
						matrix[r+1][c].state =  'w';
						unknown = unknown - 1;
					}
				}if(c + 1 < size){
						if (matrix[r][c+1].state ==  'b')
							return -1;
						else if(matrix[r][c+1].state ==  'u'){
							//printf("r = %d c = %d -SET WHITE\n", r, c+1);
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
	//printf("Function find_connection_of_white called by element %d, %d \n", row, col );
	int* flag_vector = (int *) calloc(size*size,sizeof(int));
	check_neighborhood(row, col, matrix, &flag_vector);
	for(int j = 0; j < size*size; j++){
		int r = j / size;
		int c = j % size;
		//printf("(%d, %d) - %d position of falgs' vector is %d (state = %c)\n", r, c, j, flag_vector[j], matrix[r][c].state);
		if(flag_vector[j] == 0 && matrix[r][c].state == 'w'){
			printf("Error while checking condition on white, (%d, %d) is separeted\n", r, c);
            //print_matrix(matrix);
			return -1;
		}
	}
	//printf("White condition OK\n");
	return 0;
}

int check_neighborhood(int row,int col,block** matrix, int** flag_vector){
	int* f_vector = *flag_vector;
	//printf("Travel through (%d, %d ) block. \n", row , col);
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
		//printf("Check columns and rows returned %d \n", unknown);
		//print_matrix(*matrix);
		if (unknown < 0)
			return -1;
	  unknown = check_adjacent_rules(matrix,unknown);
		//printf("Check adjacent rules returned %d \n", unknown);
		//print_matrix(*matrix);
		if (unknown < 0)
			return -1;
	}while(prev_u != unknown);
	//the do while cycle continue until appling the rules makes no change
  return unknown;
}

/*
                                 AUXILIARY FUNCTION
 ********************************************************************************/

void free_matrix(block** arr2D)
{
    for(int i=0;i<size;i++)
        free(arr2D[i]);

    free(arr2D);
}

block** malloc_matrix(){
    block** matrix =malloc(size*sizeof(block*));
    for(int i=0;i<size;++i)
        matrix[i]=malloc(size*sizeof(block));

    return matrix;
}

int copy_matrix(block*** input_matrix, block*** output_matrix){
    //The matrixs should be already allocated
    block** in = *input_matrix;
    block** out = *output_matrix;

    for(int i = 0; i < size; i++ )
        memcpy(out[i],in[i], size*sizeof(block));

    return 0;
}

int read_matrix_from_file(block*** m, char* file_name){
    int i;
    int j;
    block** matrix = *m;

    printf("Reading the matrix... \n");

    matrix = malloc_matrix();

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
    printf("Read ended. \n");
    fclose(file);
    *m = matrix;
    return 0;
}

int read_matrix_from_stdin(block*** m){
    int i;
    int j;
    block** matrix = *m;

    printf("Reading the matrix from stdin... \n");

    matrix = malloc_matrix();

    for(i = 0; i < size; i++){
        for(j = 0; j < size; j++) {
            if (!scanf("%i", &(matrix[i][j].value)))
                break;// mat[i][j] -= '0';
            matrix[i][j].state = 'u';
        }

      }
    printf("Read ended. \n");
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

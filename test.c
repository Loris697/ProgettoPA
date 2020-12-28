#include "stdio.h"
#include "stdlib.h"

typedef struct {
	int value;
	char state;
	//Il valore e':
	//'b' if the cell is deleted
	//'w' if the cell cannot be deleted
	//'u' if we don't know
}block;

int size;

int read_matrix_from_file(block*** matrix, char* file);
int print_matrix(block** matrix);
int find_connection_of_white(int row,int col,block** matrix);
int check_adjacent_rules(block*** m, int* unknown);
int check_neighborhood(int row,int col,block** matrix, int** flag_vector);

int main(int argc, char **argv)
{
	block** matrix;
	char* input_file;
	//count how many block I've to find state
	int unknown;

	int counter;
	if(argc<3 || argc>3){
		printf("\nThe program needs two arguments: \n1)The file containing the puzzle\n2) The dimension o f the puzzle");
	}

	input_file = argv[1];
	size = atoi(argv[2]);
	unknown = size*size;

	read_matrix_from_file(&matrix,input_file);

	//Valori di test
	matrix[0][1].state = 'b';
	unknown--;
	matrix[1][2].state = 'b';
	unknown--;
	printf("unknown = %d\n", unknown);

	print_matrix(matrix);

	printf("Solving\n");

	if(!check_adjacent_rules(&matrix, &unknown)){
		printf("Solution found.\n");
		printf("unknown = %d\n", unknown);
		print_matrix(matrix);
	}else
		printf("No solution found.\n");

	return 0;
}

int read_matrix_from_file(block*** m, char* file_name){
	int i;
	int j;
	block** matrix = *m;

	printf("Reading the matrix... \n");

	matrix =malloc(size*sizeof(int*));
	for(i=0;i<size;++i)
		matrix[i]=malloc(size*sizeof(int));

	FILE *file;
	file=fopen(file_name, "r");
	if(file == NULL)
		printf("\nError while opening the file. \n");

	for(i = 0; i < size; i++){
		for(j = 0; j < size; j++) {
			if (!fscanf(file, "%i", &(matrix[i][j].value)))
				break;
			matrix[i][j].state = 'u';
		}

	  }
	printf("Read ended. \n");
	fclose(file);
	*m = matrix;
	return 0;
}

int print_matrix(block** matrix){
	printf("The matrix is: \n");
	for(int i = 0; i < size; i++) {
		for(int j = 0; j < size; j++) {
			printf("%d,%c ", matrix[i][j].value, matrix[i][j].state);
		}
	printf("\n");
	}
	return 1;
}

int check_adjacent_rules(block*** m, int* unknown){
	//Find if there are two black cell adjacent
	block** matrix = *m;
	int first_white = 1;

	for(int i = 0; i < size * size; i++){
			//My current block's row and colum
			int r = i / size;
			int c = i % size;

			//if the state of the current block is black
			if(matrix[r][c].state ==  'b'){
				printf("r = %d c = %d -IS BLACK\n", r, c);
				//For all the adjacent cell apply the second rule, but before I check that the cell exist
				if (r > 0 ){
					if (matrix[r-1][c].state ==  'b')
						return -1;
					else if(matrix[r-1][c].state ==  'u'){
						  //printf("r = %d c = %d -SET WHITE\n", r-1, c);
							matrix[r-1][c].state =  'w';
							*unknown = *unknown - 1;
						}
				}if(c > 0){
					if (matrix[r][c-1].state ==  'b')
						return -1;
					else if(matrix[r][c-1].state ==  'u'){
						  //printf("r = %d c = %d -SET WHITE\n", r, c-1);
							matrix[r][c-1].state =  'w';
							*unknown = *unknown - 1;
						}
				}if(r + 1 < size){
					if (matrix[r+1][c].state ==  'b')
						return -1;
					else if(matrix[r+1][c].state ==  'u'){
						//printf("r = %d c = %d -SET WHITE\n", r+1, c);
						matrix[r+1][c].state =  'w';
						*unknown = *unknown - 1;
					}
				}if(c + 1 < size){
						if (matrix[r][c+1].state ==  'b')
							return -1;
						else if(matrix[r][c+1].state ==  'u'){
							//printf("r = %d c = %d -SET WHITE\n", r, c+1);
							matrix[r][c+1].state =  'w';
							*unknown = *unknown - 1;
						}
				}

			}else if (first_white) {
				first_white--;
				if(find_connection_of_white(r, c, matrix))
					return -1;
			}
	}

	return 0;
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
			return -1;
		}
	}
	printf("White condition OK\n");
	return 0;
}

int check_neighborhood(int row,int col,block** matrix, int** flag_vector){
	int* f_vector = *flag_vector;
	printf("Travel through (%d, %d ) block. \n", row , col);
	f_vector[row*size+col] = 1;
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

Il progetto è comspo da tre eseguibili:

random_matrix_generator
Genera una matrice allo stdout di dimensione pari al primo paramentro.

./random_matrix_generator 20

hitori_generator_par
Che genera una puzzle che ammette una soluzione.Può funzionare in due modi prendendo la matrice di partenza dal "stdin" o da un file, nel 
primo caso bisogna specificare come primo parametro "stdin", mentre nel secondo caso bisogna il file dove è presente la matrice. Il secondo paramentro è sempre la dimensione della matrice:

./random_matrix_generator 20 | mpirun -n 4 ./hitori_generator_par stdin 20

mpirun -n 4 ./hitori_generator_par input_file 20

hitori_par.c(ATTENZIONE: il numero di processi deve essere un quadrato perfetto)
Semplicemente si occupa di risolvere un puzzle. Il funzionamento è molto simile al generatore. Pure questo può funzionare in due modalità.

./random_matrix_generator 20 | mpirun -n 4 ./hitori_par stdin 20

mpirun -n 4 ./hitori_par input_file 20

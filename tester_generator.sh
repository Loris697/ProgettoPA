#!/bin/bash
GOAL=10
FOLDER=$1
SIZE=$2
TEMP="temp.txt"

for SIZE in 5 10 15 20 25 30 35
do	
	SOLUTION_FOUND=0
	for (( ; SOLUTION_FOUND < GOAL ; ))
	do
		./random_matrix_generator $SIZE > $TEMP
		for PROCESS in 1 2 4 6 8 12 16
		do
			echo "CERCO SOLUZIONE PER SOLUTION = "$SOLUTION_FOUND", SIZE = "$SIZE", PROCESS ="$PROCESS
			mpirun --oversubscribe -n $PROCESS ./hitori_generator_par $TEMP $SIZE 1> "hitori"$SIZE"by"$SIZE"/("$SIZE")solution"$SOLUTION_FOUND".in"
		done
		SOLUTION_FOUND=$((SOLUTION_FOUND+1))
	done
done


#!/bin/bash
GOAL=10
FOLDER=$1
SIZE=$2

for SIZE in 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90
do	
	for PROCESS in 1 4 16
	do
		SOLUTION_FOUND=0
		for (( ; SOLUTION_FOUND < GOAL ; ))
		do
			echo "CERCO SOLUZIONE PER SOLUTION = "$SOLUTION_FOUND", SIZE = "$SIZE", PROCESS ="$PROCESS
			mpirun -n $PROCESS --oversubscribe ./hitori_par "hitori"$SIZE"by"$SIZE"/("$SIZE")solution"$SOLUTION_FOUND".in" $SIZE 1> /dev/null
			SOLUTION_FOUND=$((SOLUTION_FOUND+1))
		done
	done
done

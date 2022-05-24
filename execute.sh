#!/bin/bash

# Command for compile c mat with openMP + MPI
# bash execute.sh compile run
compile='compile' # compile file 
run='run' # run file 
#if you use both can compile and run
args=("$@")

exec(){
	for trys in {1..2}
	do
        echo "try 1"
		for n in 10000 100000 1000000 10000000 100000000
        do
            for wnodes in 4
            do
                for numthreads in 4
                do
                    export NUM__OF_THREADS=$numthreads
                    mpirun -np $wnodes -machinefile mfile ./run $n
                done
            done
            export OMP_NUM_THREADS=4
            mpirun -np 4 -machinefile mfile ./run $n
		done
	done
}

echo " ${args[0]} ${args[1]} "

if [ "${args[0]}" == "$compile" ] &&  [ "${args[1]}" == "$run" ];
then
	mpicc -fopenmp needle_MPI.c -o run
	exec
    echo "success"
else
	echo "revisa los comandos de entrada"

fi#!/bin/bash

# Command for compile c mat with openMP + MPI
# bash execute.sh compile run
compile='compile' # compile file 
run='run' # run file 
#if you use both can compile and run
args=("$@")

exec(){
	for trys in {1..2}
	do
        echo "try 1"
		for n in 10000 100000 1000000 10000000 100000000
        do
            for wnodes in 4
            do
                for numthreads in 4
                do
                    export OMP_NUM_THREADS=$numthreads
                    mpirun -np $wnodes -machinefile mfile ./run $n
                done
            done
            export NUM_OF_THREADS=4
            mpirun -np 4 -machinefile mfile ./run $n
		done
	done
}

echo " ${args[0]} ${args[1]} "

if [ "${args[0]}" == "$compile" ] &&  [ "${args[1]}" == "$run" ];
then
	mpicc -fopenmp dartboard_MPI.c -o run
	exec
    echo "success"
else
	echo "revisa los comandos de entrada"

fi
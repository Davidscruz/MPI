#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/time.h>
#include <time.h>
#include <sys/wait.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>
#define NUM_OF_THREADS 4
//g++ -o dartboard_omp -fopenmp dartboard_omp.cpp

using namespace std; 

void lanzar(double n, double *p, double *x, double *y, int *hilos, int numranks, int rank, double* scatter, double* gather){
	
	long i; 
	long num_per_thr;
	long start_index; 
	long final_index;

	#pragma omp parallel firstprivate(p,x,y), private (i, num_per_thr, start_index, final_index), num_threads(NUM_OF_THREADS)
	{
	num_per_thr = n / omp_get_num_threads();
	start_index = omp_get_thread_num() * num_per_thr; 
	*nhilos = NUM_OF_THREADS;
	final_index = start_index + num_per_thr;

	for(i = start_index; i < final_index; i++){
        if (x[i]*x[i] + y[i]*y[i] < 1.0){
            p[omp_get_thread_num()]++;                                    
        }
	}
	}
}

void gendata(double n, double factor, double *x, double *y) {
	for(long i = 0; i < n; i++){
	x[i] = (double)rand()*factor;       // random number (0 - 1)
	y[i] = (double)rand()*factor;       // random number (0 - 1)
	}
}

int main(int argc, char* argv[]){

	double factor; // limit of rand function
    double n; // amount of trials
    double pi; // output
	double aux; //Aciertos de todos los threads
	double *p; //Aciertos 
	double *x; //Posición del dardo en x
	double *y; //Posición del dardo el y 
	
	srand(time(NULL));

	n = (double)atoi(argv[1]);
	l = 1; //Tamaño de la aguja 

	int nhilos;

	int numranks, rank, len;

    //variables to measure time
    double startTime;
    double endTime;
    double tiempo;

	//Reserva de memoria    
	p = new double [NUM_OF_THREADS];
	x = new double [(long)n];
	y = new double [(long)n];

	//Generar un conjunto de datos para los lanzamientos 
	gendata(n, factor, x, y);

    //start MPI zone
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numranks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	//Empezar tiempo
    startTime = MPI_Wtime();

    //enviar la posicion del dardo en y
    MPI_Bcast(y, n*n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    //enviar la posicion del dardo en x
    MPI_Scatter(&x[(n*n/numranks)*rank], n*n/numranks, MPI_DOUBLE, scatterNed, n*n/numranks, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    
    //MPI kernel 
    lanzar(n, p, x, y, &nhilos, numranks, rank,scatterDar, gatherDar);

    //Juntas los resultados obtenidos en cada nodo
    MPI_Gather(gatherDar, n*n/numranks, MPI_DOUBLE, &result[(n*n/numranks)*rank], n*n/numranks, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    //Parar tiempo
    endTime = MPI_Wtime();

    //Escribir el tiempo y resultados
    if(rank == 0){	
        tiempo = endTime - startTime;
        writeTime(numranks, nhilos, n, tiempo);
	}

    //Finalizar la zona MPI
    MPI_Finalize();


	//Acumular la cantidad de acierdos de cada thread
	for(long i = 0; i <NUM_OF_THREADS; i++) {
		gatherDar+= p[i];
	}

	pi = 4*aux/n;


    return 0;
}

void writeTime(int wnodos, int nhilos, int tam, double tiempo){
    FILE *f = fopen("timesOpenMPI.txt","a+");
    fprintf(f,"%i;%i;%i;%.6lf\n", wnodos, nhilos, tam, tiempo);
    fclose(f);
}
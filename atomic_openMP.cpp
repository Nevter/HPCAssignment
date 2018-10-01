/**
 * OpenMP Implementation of the atomic modeling thing
 */
#include <omp.h>  //include OMP library
#include <stdio.h>

int main(int argc, char *argv[]) {

#pragma omp parallel
	{
		printf("Hello World from thread = %d\n", omp_get_thread_num());
	}
}


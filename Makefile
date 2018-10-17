CC=g++
flags=-std=c++11

all: atomic_serial atomic_MPI atomic_openMP

#Serial Implementation
atomic_serial: atomic_serial.cpp lib/dcdplugin.c
	$(CC) -g atomic_serial.cpp -o bin/atomic_serial.o  $(flags) -O2

#OpenMP Implementation
atomic_openMP: atomic_openMP.cpp
	$(CC) -fopenmp atomic_openMP.cpp -o bin/atomic_openMP.o $(flags) -O2

#MPI Implementation
atomic_MPI: atomic_MPI.cpp
	mpic++ atomic_MPI.cpp -o bin/atomic_MPI.o -lm -O2

.cpp.o:
	$(CC) -c $< $(flags)

clean: 
	rm ./bin/*.o
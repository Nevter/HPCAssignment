CC=g++
flags=-std=c++11

all: atomic_serial atomic_MPI atomic_openMP

#Serial Implementation
atomic_serial: atomic_serial.cpp
	$(CC) -g atomic_serial.cpp -o bin/atomic_serial.o  $(flags)

#OpenMP Implementation
atomic_openMP: atomic_openMP.cpp
	$(CC) -fopenmp atomic_openMP.cpp -o bin/atomic_openMP.o $(flags)

#MPI Implementation
atomic_MPI: atomic_MPI.cpp
	$(CC) -g -o bin/atomic_MPI.o atomic_MPI.cpp $(flags)

.cpp.o:
	$(CC) -c $< $(flags)

clean: 
	rm ./bin/*.o
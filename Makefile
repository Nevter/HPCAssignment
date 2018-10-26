CC=g++
flags=-std=c++11

all: atomic_serialBF atomic_serialDC atomic_MPI atomic_openMP

#Serial BF Implementation
atomic_serialDC: atomic_serialDC.cpp lib/dcdplugin.c
	$(CC) -g -fopenmp atomic_serialDC.cpp -o bin/atomic_serialDC.o  $(flags) -O2

#Serial BF Implementation
atomic_serialBF: atomic_serialBF.cpp lib/dcdplugin.c
	$(CC) -g -fopenmp atomic_serialBF.cpp -o bin/atomic_serialBF.o  $(flags) -O2

#OpenMP Implementation
atomic_openMP: atomic_openMP.cpp
	$(CC) -fopenmp atomic_openMP.cpp -o bin/atomic_openMP.o $(flags) -O2

#MPI Implementation
atomic_MPI: atomic_MPI.cpp
	mpic++ -fopenmp atomic_MPI.cpp -o bin/atomic_MPI.o -lm -O2

.cpp.o:
	$(CC) -c $< $(flags)

clean: 
	rm ./bin/*.o
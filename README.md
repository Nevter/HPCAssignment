NVLLUK001

Compile the programs by running >make

*Note: the latest version of g++ and openMPI/mpic++ is required. 

Each program can be run using the provided scripts. 

serial:
./serialBF -i <inputfile> -o <outputfile>

openMP:
./openMP -i <inputfile> -o <outputfile>

MPI:
./MPI <numberOfProcesses> -i <inputfile> -o <outputfile>

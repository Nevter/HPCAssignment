/**
 * MPI Implementation of the atomic modeling thing
 */

//Includes 
#include "mpi.h"
#include <stdio.h>
#include <iostream>
#include <vector>
#include <string> 
#include <fstream>
#include <math.h>       /* sqrt */
#include <queue>
#include <cstdlib>

#include "lib/dcdplugin.c"


//Forward declerations
void initInputFileParameters(std::string inputFileName);
std::vector<int> getParticleSet(std::string strParticleLine);
void printToFile(std::string content, std::string outputFile);

//Structs
struct inputFileParameters {
    std::string dcdInputFile;
    int kCutOff;
    std::vector<int> particleSetA;
    std::vector<int> particleSetB;
} ifParams;

struct atomPair {
    int timeStep;
    int atomAIndex;
    int atomBIndex;
    double distance;

    std::string toString(){
        std::string asString = "";
        asString += std::to_string(timeStep) + ",";
        asString += std::to_string(atomAIndex) + ",";
        asString += std::to_string(atomBIndex) + ",";
        asString += std::to_string(distance);
        return asString;
    }
};

/**
 * 
 * 
 * 
 */
int main(int argc, char *argv[])  {
    
    int threadNum, numThreads;
    MPI_Init(&argc, &argv);                     //Start MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &threadNum);   //get rank of node's process
    MPI_Comm_size(MPI_COMM_WORLD, &numThreads);  //get the number of processes
 

    // Get the input and output file names
    std::string inputFileName = "";
    std::string outputFilename = "output.csv";
    if (argc >= 3){
        for (int i = 0; i < argc; i++){
            if (std::string(argv[i]) == "-i"){
                inputFileName = std::string(argv[i+1]);
            }
            if (std::string(argv[i]) == "-o"){
                outputFilename = std::string(argv[i+1]);
            }
        }
    }
    else {
        std::cout << "Incorrect program flags" << std::endl;
        std::cout << "-i <inputFileName>" << std::endl;
        std::cout << "[-o <outputFileName>]" << std::endl;
        exit(1);
    }

    // Read the input file parameters
    initInputFileParameters(inputFileName); 

    // Get the input file name
    //TODO: Fix this.
    //std::string fp = "/home/luke/Honours/HPC/assignment/data/"+ifParams.dcdInputFile;
    std::string fp = "/home/luke/Honours/HPC/assignment/data/example_pn3_10RU_751frames.dcd";
    const char* dcdFileName = fp.c_str();

    //Work out how many frames each process needs to process
    // Each process gets its own dcd reader
    // instance of a new DCD reading object
    int numAtoms = 0;
    void *raw_data = open_dcd_read(dcdFileName, "dcd", &numAtoms);
    dcdhandle *dcd = (dcdhandle *) raw_data;
    
    //hold info of each frame. 
    molfile_timestep_t timestep;
    timestep.coords = (float *) malloc(3 * sizeof(float) * numAtoms);
    
    //read the number of frames
    int numFrames = dcd->nsets;

    //work out which frames this thread needs to process
    int numFramesPerThread = numFrames/numThreads;
    int startingFrame = threadNum * numFramesPerThread;
    int endingFrame;
    if (threadNum==numThreads-1) endingFrame = numFrames;
    else endingFrame = startingFrame + numFramesPerThread;



    std::string threadOutput = "";

    //accelerate the read_next_timestep to the right frame
    for(int i=0; i<startingFrame;i++){
        // Read the next frame
        int rc = read_next_timestep(raw_data, numAtoms, &timestep);
    }
    
    for(int frame=startingFrame; frame<endingFrame; frame++){
    
        // Read the next frame
        int rc = read_next_timestep(raw_data, numAtoms, &timestep);

        // Holds the smallest k atomPairs    
        auto cmp = []( atomPair& lhs, atomPair& rhs) { return lhs.distance < rhs.distance; };   
        std::priority_queue<atomPair, std::vector<atomPair>, decltype(cmp)> smallestSet(cmp);

        // For each atom in set A 
        for (int a : ifParams.particleSetA){
            
            // Get each atom in set B
            for (int b : ifParams.particleSetB){
                
                //get the atoms
                float *atomA = timestep.coords + 3 * a;
                float *atomB = timestep.coords + 3 * b;

                // Get the distance between the two atoms
                double distance = std::sqrt(std::pow(atomB[0]-atomA[0],2.0) + std::pow(atomB[1]-atomA[1],2.0) + std::pow(atomB[2]-atomA[2],2.0));   
                
                atomPair ap;
                ap.timeStep = frame;
                ap.atomAIndex = a;
                ap.atomBIndex = b;
                ap.distance = distance;

                //If the distance is smaller than any in the smallest set, save it        
                if (smallestSet.size() < ifParams.kCutOff) {
                    smallestSet.push(ap);
                }
                else {
                    if (smallestSet.top().distance > distance){
                        smallestSet.pop();
                        smallestSet.push(ap);
                    }
                }  
            }
        }
        
        std::string frameOutput = "";
        // Put the smallest k atom pairs in an output String
        while (!smallestSet.empty()){
            atomPair ap = smallestSet.top();
            frameOutput = ap.toString() + "\n" + frameOutput;
            smallestSet.pop();
        }
        
        threadOutput += frameOutput;
    }   

    //close the dcd reader
    free(timestep.coords);    
    close_file_read(raw_data);

    //threadOutput = threadOutput.substr(0, threadOutput.find("\n"));
    //threadOutput += "\n";
    //MASTER THREAD
    if (threadNum == 0){
        
        std::vector<std::string> outputVector(numThreads);
        outputVector[0] = threadOutput;

        for (int i=1; i<numThreads; i++){
            //get the answers from the other threads
            int len;
            MPI_Status status;
            MPI_Recv(&len,              //address of data to store
                     1,                         //count - number of items incoming
                     MPI_INT,                   //datatype 
                     i,            //source (as rank)
                     1,                         //tag
                     MPI_COMM_WORLD,            //comm group
                     &status);        //status

            int receivedThread = status.MPI_SOURCE;
            char receivedString[len];
            //if (i == 1)
            //std::cout << "Receiving from: " << i << ", with len " << len << std::endl;
            
            MPI_Recv(receivedString,              //address of data to store
                     len,                         //count - number of items incoming
                     MPI_CHAR,                   //datatype 
                     i,            //source (as rank)
                     1,                         //tag
                     MPI_COMM_WORLD,            //comm group
                     MPI_STATUS_IGNORE);        //status
            
            
            std::string output(receivedString);
            //if (i == 1)
            //std::cout << "Receiving from: " << i << ", content: " << output;

            outputVector[i] = output;


        }
        
        //write the output to a file  
        std::string finalOutput = "";
        for (int i=0; i<outputVector.size(); i++){
            finalOutput += outputVector[i];
        }

        printToFile(finalOutput, outputFilename);

        //std::cout << "\n\n" << finalOutput << std::endl;


    }
    // SLAVE THREADS
    else {
        

        int len = threadOutput.length();
        //std::cout << "Thread " << threadNum << ", sending(" << len << "): " << threadOutput;
        char outputArr[len+1];
        strcpy(outputArr, threadOutput.c_str());

        //send length of string to main
        MPI_Send(&len,                     //address of data to send
                 1,                     //count - number of items
                 MPI_INT,               //datatype
                 0,                     //destination (as rank)
                 1,                     //tag
                 MPI_COMM_WORLD);       //comm group

        //send answer to main
        MPI_Send(outputArr,                     //address of data to send
                 len,                     //count - number of items
                 MPI_CHAR,               //datatype
                 0,                     //destination (as rank)
                 1,                     //tag
                 MPI_COMM_WORLD);       //comm group
        
    }


    //std::cout << "Thread " << threadNum << ", start: " << startingFrame << ", end: " << endingFrame << std::endl;
    //std::cout << "Thread " << threadNum << " finishing." << std::endl;
    MPI_Finalize();                             //End MPI

    return 0;
}


void printToFile(std::string content, std::string outputFile){
    std::ofstream outFile;
    outFile.open (outputFile);
    outFile << content;
    outFile.close();
}


/**
 * Read input file parameters from an input file as specified 
 * in the assignment.
 */
void initInputFileParameters(std::string inputFileName){
    // Read the input file into an array of each line;
    std::string inputFileLines[4];
    std::string line;
    std::ifstream inputFile (inputFileName);
    int lineCounter = 0;
    if (inputFile.is_open()){
        while ( getline (inputFile,line)){
            inputFileLines[lineCounter] = line;
            lineCounter++;
        }
        inputFile.close();
    }
    else {
        std::cout << "Unable to open file:" << inputFileName << std::endl;
        exit(1); 
    }

    // Set the inputDCD and k value in the input file parameters
    ifParams.dcdInputFile = inputFileLines[0];
    ifParams.kCutOff = stoi(inputFileLines[1]);

    // Set the particle sets
    ifParams.particleSetA = getParticleSet(inputFileLines[2]);
    ifParams.particleSetB = getParticleSet(inputFileLines[3]);
}

/**
 * Get a list of the particles in a particle set. 
 * The input is a string containing atom indices 
 * comprising particle set A (a range or itemized 
 * list of comma separated values or mix of two).
 */
std::vector<int> getParticleSet(std::string strParticleLine){
    std::vector<std::string> strParticleSet;
    std::vector<int> particleSet;
    
    // Split the strParticle Line by commas, separating individual values and value ranges
    std::string s = strParticleLine;
    std::string delimiter = ",";
    size_t pos = 0;
    std::string token;
    while ((pos = s.find(delimiter)) != std::string::npos) {
        token = s.substr(0, pos);
        strParticleSet.push_back(token);
        s.erase(0, pos + delimiter.length());
    }
    strParticleSet.push_back(s);
    
    // Get each particle, including unfolding ranges
    delimiter = "-";
    for(std::string entry : strParticleSet) {
        // range
        if (entry.find(delimiter) != std::string::npos){
            pos = entry.find(delimiter);
            std::string first = entry.substr(0,pos);
            std::string last = entry.substr(pos+1,entry.size());
            for (int i = std::stoi(first); i <= std::stoi(last); i++){
                particleSet.push_back(i);
            }
        }
        // single value
        else {
            particleSet.push_back(std::stoi(entry));
        }
    }
    return particleSet;
}

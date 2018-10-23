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
#include <sys/time.h>
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
    std::string outputFilename = "output";
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

    struct timeval start, end;
    if (threadNum == 0){
        //start a timer
        gettimeofday(&start, NULL);
    }

    // Read the input file parameters
    initInputFileParameters(inputFileName); 

    // Get the input file name
    std::string fp = ifParams.dcdInputFile;
    fp = fp.substr(0, fp.find(".dcd")+4);
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
                     MPI_ANY_SOURCE,            //source (as rank)
                     1,                         //tag
                     MPI_COMM_WORLD,            //comm group
                     &status);        //status

            int receivedThread = status.MPI_SOURCE;
            char receivedString[len];
            
            MPI_Recv(receivedString,              //address of data to store
                     len,                         //count - number of items incoming
                     MPI_CHAR,                   //datatype 
                     receivedThread,            //source (as rank)
                     1,                         //tag
                     MPI_COMM_WORLD,            //comm group
                     MPI_STATUS_IGNORE);        //status
            
            
            std::string output(receivedString);
            outputVector[receivedThread] = output;

        }
        
        //collate the answers
        std::string finalOutput = "";
        for (int i=0; i<outputVector.size(); i++){
            std::string thOut = outputVector[i];
            thOut = thOut.substr(0,thOut.find("!"));
            finalOutput += thOut;
        }

        //end timer
        gettimeofday(&end, NULL);
        float delta = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1.e6;
        
        //write the output to a file 
        printToFile(finalOutput, outputFilename);

        //print out time 
        std::string progtime = "Time: " + std::to_string(delta);
        printToFile(progtime, outputFilename+"Time");

    }
    // SLAVE THREADS
    else {
        
        threadOutput = threadOutput + "!";
        int len = threadOutput.length();
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


    MPI_Finalize();                             //End MPI

    return 0;
}


void printToFile(std::string content, std::string outputFile){
    outputFile += ".txt";
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

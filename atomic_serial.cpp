/**
 * Serial Implementation of the atomic modeling thing
 */

// Includes
#include <stdio.h>
#include <iostream>
#include <vector>
#include <string> 
#include <fstream>
#include <math.h>       /* sqrt */
#include <queue>

#include <cstdlib>

#include "lib/dcdplugin.c"


// Forward declerations
void initInputFileParameters(std::string inputFileName);
std::vector<int> getParticleSet(std::string strParticleLine);
void printToFile(std::string content, std::string outputFile);

// Structs
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
    std::cout << " ~~~ Atomic Serial ~~~ " << std::endl;
    std::cout << "Config Input file: " << inputFileName << std::endl;
    std::cout << "Output file: " << outputFilename << std::endl;

    // Read the input file parameters
    initInputFileParameters(inputFileName);

    // Open the input DCD file
    //TODO: Fix this.
    //std::string fp = "/home/luke/Honours/HPC/assignment/data/"+ifParams.dcdInputFile;
    std::string fp = "/home/luke/Honours/HPC/assignment/data/example_pn3_10RU_751frames.dcd";
    const char* dcdFileName = fp.c_str();

    // instance of a new DCD reading object
    int numAtoms = 0;
    
    // read the header
    void *raw_data = open_dcd_read(dcdFileName, "dcd", &numAtoms);
    dcdhandle *dcd = (dcdhandle *) raw_data;
    
    //hold info of each frame. 
    molfile_timestep_t timestep;
    timestep.coords = (float *) malloc(3 * sizeof(float) * numAtoms);
    
    //read and print out the atom 1 from the first 5 frames
    int numFrames = dcd->nsets;
    std::cout << "Number of frames: " << numFrames << std::endl;
    std::cout << "Number of atoms: " << numAtoms << std::endl;

    //DEBUG: change number of frames so not looking at all frames
    //numFrames = 1;

    std::string output = "";

    for(int frame=0; frame<numFrames; frame++){
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
        output += frameOutput;
        
    }   

    //close the dcd reader
    free(timestep.coords);    
    close_file_read(raw_data);

    //write the output to a file  
    printToFile(output, outputFilename);
    
    //return a success
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
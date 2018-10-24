/**
 * Serial Implementation of the atomic modeling thing
 */

// Includes
#include <omp.h>  //include OMP library
#include <stdio.h>
#include <iostream>
#include <vector>
#include <string> 
#include <fstream>
#include <math.h>       /* sqrt */
#include <queue>
#include <sys/time.h>
#include <cstdlib>
#include <algorithm>    // std::sort
#include "lib/dcdplugin.c"

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

struct atom {
    float x;
    float y;
    float z;
    int id;
};

// Forward declerations
void initInputFileParameters(std::string inputFileName);
std::vector<int> getParticleSet(std::string strParticleLine);
void printToFile(std::string content, std::string outputFile);
std::vector<atomPair> findNearestSet(int k, std::vector<atom> setA,  std::vector<atom> setB, float lowerBound, float upperBound, int frame);
std::vector<atomPair> getSmallestSet(int k, std::vector<atom> setA,  std::vector<atom> setB, float lowerBound, float upperBound, int frame, int numSegments);

/**
 * 
 * 
 * 
 */
int main(int argc, char *argv[])  {
    
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
    std::cout << " ~~~ Atomic Serial ~~~ " << std::endl;
    std::cout << "Config Input file: " << inputFileName << std::endl;
    std::cout << "Output file: " << outputFilename << std::endl;

    // Read the input file parameters
    initInputFileParameters(inputFileName);

    //start a timer
    double start = omp_get_wtime();

    // Open the input DCD file
    std::string fp = ifParams.dcdInputFile;
    fp = fp.substr(0, fp.find(".dcd")+4);
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

    std::string output = "";

    for(int frame=0; frame<numFrames; frame++){
        // Read the next frame
        int rc = read_next_timestep(raw_data, numAtoms, &timestep);

        // Holds the smallest k atomPairs    
        auto cmp = []( atomPair& lhs, atomPair& rhs) { return lhs.distance < rhs.distance; };   
        std::priority_queue<atomPair, std::vector<atomPair>, decltype(cmp)> smallestSet(cmp);
        
        //create a set of atoms
        std::vector<atom> setA;
        std::vector<atom> setB;
        // Get each atom in set A 
        for (int a : ifParams.particleSetA){
            //get the atom
            float *atomA = timestep.coords + 3 * a;
            atom at;
            at.x = atomA[0];
            at.y = atomA[1];
            at.z = atomA[2];
            at.id = a;
            setA.push_back(at);
        }
        // Get each atom in set B
        for (int b : ifParams.particleSetB){
            //get the atoms
            float *atomB = timestep.coords + 3 * b;
            atom at;
            at.x = atomB[0];
            at.y = atomB[1];
            at.z = atomB[2];
            at.id = b;
            setB.push_back(at);
        }
        //sort set A to be in accending x order 
        auto compAtom = [](atom a1, atom a2){ return a1.x < a2.x; };
        std::sort(setA.begin(), setA.end(), compAtom);
        //sort set B to in accending x order 
        std::sort(setB.begin(), setB.end(), compAtom);
        
        float lowerBound = std::min(setA.front().x, setB.front().x);
        float upperBound = std::min(setA.back().x, setB.back().x);
        

        std::vector<atomPair> frameSet = getSmallestSet(ifParams.kCutOff, setA, setB, lowerBound, upperBound, frame, 4);

        std::string frameOutput = "";
        for (int i=0; i < ifParams.kCutOff; i++){
            frameOutput += frameSet[i].toString() + "\n";
        }

        output += frameOutput;
 
    }   

    //end timer
    double end=omp_get_wtime();
    std::string progtime = "Time: " + std::to_string(end-start);
    
    //close the dcd reader
    free(timestep.coords);    
    close_file_read(raw_data);

    //write the output to a file  
    printToFile(output, outputFilename);
    
    //print out time 
    std::string info = "~Serial~ \nInput file: " + inputFileName + "\nNumber of threads: 1" + '\n' + progtime + '\n';
    printToFile(info, outputFilename+"Info");
    std::cout << info << std::endl;

    //return a success
    return 0;
}

std::vector<atomPair> getSmallestSet(int k, std::vector<atom> setA,  std::vector<atom> setB, float lowerBound, float upperBound, int frame, int numSegments){

    if (numSegments == 0){
        return findNearestSet(k, setA, setB, lowerBound, upperBound, frame);
    }
    

    float midPoint =  upperBound - lowerBound;

    std::vector<atomPair> lowerHalf = getSmallestSet(ifParams.kCutOff, setA, setB, lowerBound, midPoint, frame, numSegments-1);
    std::vector<atomPair> upperHalf = getSmallestSet(ifParams.kCutOff, setA, setB, midPoint, upperBound, frame, numSegments-1);
    std::vector<atomPair> middleSection = findNearestSet(ifParams.kCutOff, setA, setB, midPoint-1, midPoint+1, frame);

    std::vector<atomPair> combined;
    combined.reserve(lowerHalf.size() + upperHalf.size() + middleSection.size() );
    combined.insert(combined.end(), lowerHalf.begin(), lowerHalf.end());
    combined.insert(combined.end(), upperHalf.begin(), upperHalf.end());
    combined.insert(combined.end(), middleSection.begin(), middleSection.end());

    //sort by smallest 
    auto compAtomPair = [](atomPair a1, atomPair a2){ return a1.distance < a2.distance; };
    std::sort(combined.begin(), combined.end(), compAtomPair);
    combined.resize(k);
    
    return combined;

}


std::vector<atomPair> findNearestSet(int k, std::vector<atom> setA,  std::vector<atom> setB, float lowerBound, float upperBound, int frame){
    // Holds the smallest k atomPairs    
    auto cmp = []( atomPair& lhs, atomPair& rhs) { return lhs.distance < rhs.distance; };   
    std::priority_queue<atomPair, std::vector<atomPair>, decltype(cmp)> smallestSet(cmp);

    int count = 0;

    for (atom a : setA){
        if (a.x >= lowerBound && a.x <= upperBound){
            for (atom b : setB){
                if (b.x >= lowerBound && a.x <= upperBound){
                    double distance = std::sqrt(std::pow(b.x-a.x,2.0) + std::pow(b.y-a.y,2.0) + std::pow(b.z-a.z,2.0));

                    atomPair ap;
                    ap.timeStep = frame;
                    ap.atomAIndex = a.id;
                    ap.atomBIndex = b.id;
                    ap.distance = distance;

                    //If the distance is smaller than any in the smallest set, save it        
                    if (smallestSet.size() < k) {
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
        }
    }

    std::vector<atomPair> ksmallest;
    // Put the smallest k atom pairs in an output String
    while (!smallestSet.empty()){
        atomPair ap = smallestSet.top();
        ksmallest.push_back(ap);
        smallestSet.pop();
    }

    return ksmallest;

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
/**
 * Serial Implementation of the atomic modeling thing
 */

// Includes
#include <stdio.h>
#include <iostream>
#include <vector>
#include <string> 
#include <fstream>
#include <algorithm>
#include <math.h>       /* sqrt */
#include <set> 
#include <queue>
#include <functional>

#include <cstdlib>

#include "lib/array_tools.hpp"
#include "lib/dcd_r.hpp"
#include "lib/dcd_r.cpp"
#include "lib/dcd.cpp"


// Forward declerations
void initInputFileParameters(std::string inputFileName);
std::vector<int> getParticleSet(std::string strParticleLine);
double calculate3DDistance(float x0, float y0, float z0, float x1, float y1, float z1);

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

    // instance of a new object DCD_R attached to a dcd file 
    DCD_R dcdf(dcdFileName);
    
    // read the header
    dcdf.read_header();
    
    //floats to hold the x,y,z coords of each atom in each frame. 
    const float *x,*y,*z;
    
    //read and print out the atom 1 from the first 5 frames
    int numFrames = dcdf.getNFILE();
    int numAtoms = dcdf.getNATOM();
    std::cout << "Number of frames: " << numFrames << std::endl;
    std::cout << "Number of atoms: " << numAtoms << std::endl;

    //DEBUG: change number of frames so not looking at all frames
    //numFrames = 1;

    std::vector<atomPair> outputV;

    for(int frame=0; frame<numFrames; frame++){
        // Read the next frame
        dcdf.read_oneFrame();

        // Get the atom positions
        x=dcdf.getX();
        y=dcdf.getY();
        z=dcdf.getZ();

        // Holds the smallest k atomPairs    
        auto cmp = []( atomPair& lhs, atomPair& rhs) { return lhs.distance < rhs.distance; };   
        std::priority_queue<atomPair, std::vector<atomPair>, decltype(cmp)> smallestSet(cmp);

        // For each atom in set A 
        for (int a : ifParams.particleSetA){
            
            // Get each atom in set B
            for (int b : ifParams.particleSetB){
                
                // Ensure that we aren't looking at the same atom 
                if (a==b) continue;
                
                // Get the distance between the two atoms
                double distance = calculate3DDistance(x[a],y[a],z[a],x[b],y[b],z[b]);
                

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
                    if (smallestSet.top().distance > ap.distance){
                        smallestSet.pop();
                        smallestSet.push(ap);
                    }
                }
                
            }
        }
    
        // Put the smallest k atom pairs in a master vector
        std::vector<atomPair> smallestK;
        for (int i =0; i < ifParams.kCutOff; i++){
            smallestK.push_back(smallestSet.top()); smallestSet.pop();
        }
        std::reverse(smallestK.begin(),smallestK.end());
        outputV.insert(outputV.end(), smallestK.begin(), smallestK.end());
    }   

    //print the output vector 
    /*
    for (atomPair ap : outputV){
        std::cout << ap.toString() << std::endl;
    }
    */
    
    //return a success
    return 0;
}

double calculate3DDistance(float xa, float ya, float za, float xb, float yb, float zb){

    return std::sqrt(std::pow(xb-xa,2.0) + std::pow(yb-ya,2.0) + std::pow(zb-za,2.0));
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
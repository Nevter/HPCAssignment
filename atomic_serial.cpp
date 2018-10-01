/**
 * Serial Implementation of the atomic modeling thing
 */

// Includes
#include <stdio.h>
#include <iostream>
#include <vector>
#include <string> 
#include <fstream>

// Forward declerations
void initInputFileParameters(std::string inputFileName);
std::vector<int> getParticleSet(std::string strParticleLine);

//Structs
struct inputFileParameters {
    std::string dcdInputFile;
    std::string kCutOff;
    std::vector<int> particleSetA;
    std::vector<int> particleSetB;
} ifParams;

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
    std::cout << "Input file: " << inputFileName << std::endl;
    std::cout << "Output file: " << outputFilename << std::endl;

    // Read the input file parameters
    initInputFileParameters(inputFileName);

    
    
    return 0;
}

void initInputFileParameters(std::string inputFileName){
    //Read the input file into an array of each line;
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
        std::cout << "Unable to open file" << std::endl; 
    }

    //Set the inputDCD and k value in the input file parameters
    ifParams.dcdInputFile = inputFileLines[0];
    ifParams.kCutOff = inputFileLines[1];

    //Set the particle sets
    ifParams.particleSetA = getParticleSet(inputFileLines[2]);
    ifParams.particleSetB = getParticleSet(inputFileLines[3]);

}

std::vector<int> getParticleSet(std::string strParticleLine){
    std::vector<std::string> strParticleSetA;
    
    std::string s = strParticleLine;
    std::string delimiter = ",";
    size_t pos = 0;
    std::string token;
    while ((pos = s.find(delimiter)) != std::string::npos) {
        token = s.substr(0, pos);
        strParticleSetA.push_back(token);
        s.erase(0, pos + delimiter.length());
    }
    strParticleSetA.push_back(s);
    
    std::vector<int> particleSetA;

    delimiter = "-";
    for(std::string entry : strParticleSetA) {
        //a range
        if (entry.find(delimiter) != std::string::npos){
            pos = entry.find(delimiter);
            std::string first = entry.substr(0,pos);
            std::string last = entry.substr(pos+1,entry.size());
            for (int i = std::stoi(first); i <= std::stoi(last); i++){
                particleSetA.push_back(i);
            }
        }
        //a single value
        else {
            particleSetA.push_back(std::stoi(entry));
        }
    }
    return particleSetA;

}
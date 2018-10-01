/**
 * Serial Implementation of the atomic modeling thing
 */
#include <stdio.h>
#include <iostream>

int main(int argc, char *argv[])  {
    
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
    }


    printf("Hello Serial World \n");

    return 0;
}

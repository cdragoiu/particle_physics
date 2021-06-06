#include "HCALCalibration/HFCalibration/bin/TreeWrapper.h"
#include "HCALCalibration/HFCalibration/bin/Calibration.h"
#include "TFile.h"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char **argv) {
    
    // program arguments
    if (argc < 3) {
        std::cout << std::endl;
        std::cout << "usage: " << argv[0] << " txtFile outFile" << std::endl;
        std::cout << "       txtFile = file with list of input files" << std::endl;
        std::cout << "       outFile = output root file" << std::endl;
        std::cout << std::endl;
        return 0;
    }
    std::string txtFileName(argv[1]), outFileName(argv[2]);
    
    // load input files
    std::ifstream txtFile(txtFileName.c_str());
    std::vector<std::string> inFileNames;
    while (!txtFile.eof()) {
        std::string name;
        txtFile >> name;
        inFileNames.push_back(name);
    }
    
    // setup analysis
    TreeWrapper *tree = new TreeWrapper(inFileNames);
    TFile *outFile = new TFile(outFileName.c_str(), "RECREATE");
    Calibration *calib = new Calibration();
    
    // loop over events
    for (unsigned int e = 0; e != tree->Size(); ++e) {
        tree->Read(e);
        calib->ProcessEvent(tree->Event());
    }
    
    // wrap it up
    std::cout << "\n" << "Total number of events processed: " << tree->Size() << std::endl;
    calib->Write();
    delete calib;
    outFile->Close();
    delete outFile;
    delete tree;
    return 0;
    
}

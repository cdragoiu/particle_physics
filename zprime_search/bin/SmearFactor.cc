#include "TFile.h"
#include "TRandom3.h"
#include "TTree.h"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char **argv) {
    
    // program arguments
    if (argc < 2) {
        std::cout << "usage: " << argv[0] << " txtFile" << std::endl;
        std::cout << "       txtFile    = file with list of input files" << std::endl;
        return 0;
    }
    std::string txtFileName(argv[1]);
    
    // load input files
    std::ifstream txtFile(txtFileName.c_str());
    std::vector<std::string> inFileNames;
    while (!txtFile.eof()) {
        std::string name;
        txtFile >> name;
        inFileNames.push_back(name);
    }
    
    // process files
    TRandom3 rnd = TRandom3(0);
    for (unsigned int i = 0; i != inFileNames.size(); ++i) {
        std::cout << "--> " << inFileNames.at(i) << std::endl;
        TFile *file = TFile::Open(inFileNames.at(i).c_str(), "update");
        TTree *tree = (TTree*)file->Get("Events");
        float hfEle_smrf[99];
        TBranch *branch = tree->Branch("hfEle_smrf", hfEle_smrf, "hfEle_smrf[hfEle_size]/F");
        unsigned int hfEle_size;
        tree->SetBranchAddress("hfEle_size", &hfEle_size);
        unsigned int N = tree->GetEntries();
        for (unsigned int e = 0; e != N; ++e) {
            tree->GetEntry(e);
            for (unsigned int h = 0; h != hfEle_size && h != 99; ++h) {
                hfEle_smrf[h] = rnd.Gaus(90.9, 4.8) / 91.2;
            }
            branch->Fill();
        }
        tree->Write("Events", 6);
        file->Close();
    }
    
    return 0;
}

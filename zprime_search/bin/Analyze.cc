#include "MyAnalyses/CommonTools/interface/Utils.h"
#include "MyAnalyses/ZprimeSearch/bin/Corrections.h"
#include "MyAnalyses/ZprimeSearch/bin/TreeWrapper.h"
#include "MyAnalyses/ZprimeSearch/bin/Zprime.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char **argv) {
    
    // program arguments
    if (argc < 5) {
        std::cout<<"usage: "<<argv[0]<<" txtFile outFile classType runType [iterations]"<<std::endl;
        std::cout<<"       txtFile    = file with list of input files"<<std::endl;
        std::cout<<"       outFile    = output root file"<<std::endl;
        std::cout<<"       classType  = data, mc"<<std::endl;
        std::cout<<"       runType    = qcd, unf, ele, hf"<<std::endl;
        std::cout<<"                    ebeb, ebee, eeeb, eeee, ebhf, eehf"<<std::endl;
        std::cout<<"                    sysEleErg, sysHfErg, sysPu"<<std::endl;
        std::cout<<"                    sysTrg8, sysTrg17, sysTrg27"<<std::endl;
        std::cout<<"                    sysEleId, sysHfId, sysHfEta"<<std::endl;
        std::cout<<"                    sysFsr, sysCt10pdf, sysNnpdf"<<std::endl;
        std::cout<<"       iterations = number of systematic variations"<<std::endl;
        return 0;
    }
    std::string txtFileName(argv[1]), outFileName(argv[2]);
    std::string classType(argv[3]), runType(argv[4]);
    unsigned int N(0);
    if (argc == 6) N = atoi(argv[5]);
    
    // load input files
    std::ifstream txtFile(txtFileName.c_str());
    std::vector<std::string> inFileNames;
    while (!txtFile.eof()) {
        std::string name;
        txtFile >> name;
        inFileNames.push_back(name);
    }
    
    // setup analysis
    Corrections corr;
    TreeWrapper *tree = new TreeWrapper(inFileNames);
    TFile *outFile = new TFile(outFileName.c_str(), "RECREATE");
    std::vector<std::string> dirs;
    std::vector<int> ords;
    for (unsigned int i = 0; i != N; ++i) {
        if (runType.find("pdf") != std::string::npos) {
            dirs.push_back(runType+"_"+ToStr(i+1));
            ords.push_back(i+1);
        }
        else {
            dirs.push_back(runType+"_max_"+ToStr(i+1));
            ords.push_back(i+1);
            dirs.push_back(runType+"_min_"+ToStr(i+1));
            ords.push_back(-(i+1));
        }
    }
    std::vector<Zprime*> runs;
    for (unsigned int i = 0; i != dirs.size(); ++i) {
        outFile->mkdir(dirs.at(i).c_str());
        outFile->cd(dirs.at(i).c_str());
        runs.push_back(new Zprime(classType+"_"+runType));
    }
    if (N == 0) runs.push_back(new Zprime(classType+"_"+runType));
    
    // loop over events
    int percent(-1);
    for (unsigned int e = 0; e != tree->Size(); ++e) {
        int percent_tmp = 100.0 * (e + 1) / tree->Size();
        if (percent_tmp != percent) {
            percent = percent_tmp;
            ShowPercentage(percent);
        }
        tree->Read(e);
        std::vector<TLorentzVector> electrons;
        for (unsigned int v = 0; v != tree->Event().electrons.size(); ++v)
            electrons.push_back(tree->Event().electrons.at(v).p4);
        std::vector<TLorentzVector> hfElectrons;
        for (unsigned int v = 0; v != tree->Event().hfElectrons.size(); ++v)
            hfElectrons.push_back(tree->Event().hfElectrons.at(v).p4);
        for (unsigned int i = 0; i != runs.size(); ++i) {
            for (unsigned int v = 0; v != electrons.size(); ++v) {
                if (runType.find("sysEleErg") != std::string::npos) {
                    double sf(1.0);
                    if (fabs(electrons.at(v).Eta()) < 1.442) sf = 1.0 + 0.004*ords.at(i);
                    if (fabs(electrons.at(v).Eta()) > 1.566) sf = 1.0 + 0.005*ords.at(i);
                    tree->Event().electrons.at(v).p4.SetPxPyPzE(sf*electrons.at(v).Px(),
                                                                sf*electrons.at(v).Py(),
                                                                sf*electrons.at(v).Pz(),
                                                                sf*electrons.at(v).E());
                }
            }
            for (unsigned int v = 0; v != hfElectrons.size(); ++v) {
                double sf(1.0);
                if (classType == "mc") {
                    if (runType.find("sysHfErg") != std::string::npos)
                         sf = corr.GetHfEleEnergyCorrMC(hfElectrons.at(v).Et(), ords.at(i));
                    else sf = corr.GetHfEleEnergyCorrMC(hfElectrons.at(v).Et());
                    sf *= tree->Event().hfElectrons.at(v).smrf;
                }
                else {
                    if (runType.find("sysHfErg") != std::string::npos)
                         sf = corr.GetHfEleEnergyCorrData(hfElectrons.at(v).Et(), ords.at(i));
                    else sf = corr.GetHfEleEnergyCorrData(hfElectrons.at(v).Et());
                }
                tree->Event().hfElectrons.at(v).p4.SetPxPyPzE(sf*hfElectrons.at(v).Px(),
                                                              sf*hfElectrons.at(v).Py(),
                                                              sf*hfElectrons.at(v).Pz(),
                                                              sf*hfElectrons.at(v).E());
            }
            if (ords.size()) runs.at(i)->ProcessEvent(tree->Event(), corr, ords.at(i));
            else runs.at(i)->ProcessEvent(tree->Event(), corr);
        }
    }
    
    // wrap it up
    std::cout << "\n" << "Total number of events processed: " << tree->Size() << std::endl;
    for (unsigned int i = 0; i != runs.size(); ++i) {
        if (dirs.size()) outFile->cd(dirs.at(i).c_str());
        runs.at(i)->Write();
        delete runs.at(i);
    }
    outFile->Close();
    delete outFile;
    delete tree;
    return 0;
    
}

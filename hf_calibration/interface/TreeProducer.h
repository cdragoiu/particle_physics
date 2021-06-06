#ifndef TreeProducer_GUARD
#define TreeProducer_GUARD

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "TFile.h"
#include "TTree.h"
#include <string>

class TreeProducer : public edm::EDAnalyzer {
    
public:
    
    // constructor
	TreeProducer(const edm::ParameterSet &config);
    
    // destructor
	~TreeProducer();
    
private:
    
    // analyze event
	void analyze(const edm::Event &event, const edm::EventSetup &setup);
	
	// input variables
	std::string outFileName;
    edm::InputTag vertexTag;
    edm::InputTag electronTag;
    edm::InputTag conversionTag;
    edm::InputTag beamSpotTag;
    edm::InputTag hfElectronTag;
    edm::InputTag hfClusterMapTag;
    
    // output variables
	TFile *outFile;
	TTree *tree;
    // event information
    unsigned int ev_run, ev_lumi, ev_nr, ev_bunch;
    // primary vertex information
    unsigned int pv_size, pv_nTracks[99];
    float pv_chi2[99], pv_normalizedChi2[99], pv_ndof[99], pv_x[99], pv_y[99], pv_z[99], pv_rho[99];
    // electron information
    unsigned int ele_size;
    int ele_charge[99];
    float ele_pt[99], ele_phi[99], ele_eta[99], ele_energy[99], ele_etaSC[99];
    float ele_dEtaIn[99], ele_dPhiIn[99], ele_sigmaIEtaIEta[99], ele_full5x5SigmaIEtaIEta[99];
    float ele_HoE[99], ele_d0[99], ele_dZ[99], ele_IoEmIoP[99];
    float ele_iso[99];
    bool ele_hasConversion[99];
    unsigned int ele_missingHits[99];
    // HF electron information
    unsigned int hfEle_size;
    float hfEle_pt[99], hfEle_phi[99], hfEle_eta[99], hfEle_energy[99];
    float hfEle_e1x1[99], hfEle_e3x3[99], hfEle_e5x5[99];
    float hfEle_eLong1x1[99], hfEle_eLong3x3[99], hfEle_eLong5x5[99];
    float hfEle_eShort1x1[99], hfEle_eShort3x3[99], hfEle_eShort5x5[99];
    float hfEle_e9e25[99], hfEle_eCOREe9[99], hfEle_eCore[99], hfEle_eSeL[99];
    
};

#endif

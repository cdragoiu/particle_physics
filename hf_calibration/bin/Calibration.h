#ifndef Calibration_GUARD
#define Calibration_GUARD

#include "HCALCalibration/HFCalibration/bin/TreeWrapper.h"
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include <string>
#include <vector>

class Calibration {
    
public:
    
    // constructor
    Calibration();
    
    // destructor
    ~Calibration();
    
    // process event
    void ProcessEvent(const EventInfo &event);
    
    // write histograms
    void Write();
    
private:
    
    // select electrons
    void SelectElectrons(const std::vector<EventInfo::Electron> &electrons,
                         const std::vector<EventInfo::HfElectron> &hfElectrons);
    
    // fill histograms
    void FillHistograms();
    
    // book histograms
    void BookHistograms();
    
    // test electron kinematics
    bool IsGoodEleP4(const TLorentzVector &p4) const;
    
    // test electron ID
    bool IsGoodEleID(const EventInfo::Electron &ele) const;
    
    // test HF electron kinematics
    bool IsGoodHfEleP4(const TLorentzVector &p4) const;
    
    // test HF electron ID
    bool IsGoodHfEleID(const EventInfo::HfElectron &hfEle) const;
    
    // selected electrons
    const EventInfo::Electron *electron;
    const EventInfo::HfElectron *hfElectron;
    bool goodElectrons;
    
    // histograms
    TH1D *ele_pt, *ele_eta, *ele_phi;
    TH1D *hfEle_et, *hfEle_eta, *hfEle_phi;
    TH1D *hfEle_eLong3x3, *hfEle_eLong5x5, *hfEle_eLCeL9, *hfEle_eS9eL9;
    TH1D *hfEle_eL9eL25, *hfEle_cut2D;
    TH1D *ll_pt, *ll_mass, *ll_rapidity, *ll_phi;
    TH1D *ele_Mcut_pt, *ele_Mcut_eta, *ele_Mcut_phi;
    TH1D *hfEle_Mcut_et, *hfEle_Mcut_eta, *hfEle_Mcut_phi;
    TH2D *hfEle_etaX_corrY;
    
};

#endif

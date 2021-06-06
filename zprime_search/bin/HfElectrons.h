#ifndef HfElectrons_GUARD
#define HfElectrons_GUARD

#include "MyAnalyses/CommonTools/interface/P4Hists.h"
#include "MyAnalyses/ZprimeSearch/bin/TreeWrapper.h"
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include <string>
#include <vector>

class HfElectrons {
    
public:
    
    // constructor
    HfElectrons(bool _isMC, bool _doQCD, bool _doHist);
    
    // destructor
    ~HfElectrons();
    
    // get electron
    const EventInfo::Electron& GetEle() const;
    
    // get HF electron
    const EventInfo::HfElectron& GetHfEle() const;
    
    // test electron kinematics
    bool IsGoodEleP4(const TLorentzVector &p4, double etaMin, double etaMax) const;
    
    // test HF electron kinematics
    bool IsGoodHfEleP4(const TLorentzVector &p4) const;
    
    // get status
    bool IsGood() const;
    
    // reset information
    void Reset();
    
    // select electrons
    void SelectElectrons(const std::vector<EventInfo::Electron> &electrons,
                         double etaMin, double etaMax,
                         const std::vector<EventInfo::HfElectron> &hfElectrons);
    
    // fill eventCount histogram
    void FillEventCount(std::string tag = "");
    
    // fill trigger histograms
    void FillTriggerHistograms(const EventInfo &event, double weight);
    
    // fill raw histograms
    void FillRawHistograms(const EventInfo &event, double weight);
    
    // fill raw MC only histograms
    void FillRawMcHistograms(const EventInfo &event, const EventInfo::Particle &ele,
                             const EventInfo::Particle &pos, double weight);
    
    // fill MC only histograms
    void FillMcHistograms(const EventInfo::Particle &ele,
                          const EventInfo::Particle &pos, double weight);
    
    // fill histograms
    void FillHistograms(double weight);
    
    // write histograms
    void Write();
    
private:
    
    // book histograms
    void BookHistograms();
    
    // test electron ID
    bool IsGoodEleID(const EventInfo::Electron &ele) const;
    
    // test HF electron ID
    bool IsGoodHfEleID(const EventInfo::HfElectron &hfEle) const;
    
    // private objects
    bool isMC, doQCD, doHist;
    
    // selected electrons
    const EventInfo::Electron *electron;
    const EventInfo::HfElectron *hfElectron;
    unsigned int eleSize, hfEleSize;
    bool goodElectronsNoID, goodEleID, goodHfEleID;
    bool goodElectrons;
    
    // histograms
    TH1D *eventCount;
    TH1D *raw_ele_size, *raw_hfEle_size;
    TH2D *raw_eleSizeX_pvSizeY, *raw_hfEleSizeX_pvSizeY;
    TH1D *raw_ele_dEtaIn, *raw_ele_dPhiIn, *raw_ele_sigmaIEtaIEta, *raw_ele_HoE;
    TH1D *raw_ele_d0, *raw_ele_dZ, *raw_ele_IoEmIoP, *raw_ele_iso;
    TH1D *raw_hfEle_eCore, *raw_hfEle_eLong1x1, *raw_hfEle_eLong3x3, *raw_hfEle_eLong5x5;
    TH1D *raw_hfEle_eShort1x1, *raw_hfEle_eShort3x3, *raw_hfEle_eShort5x5;
    TH1D *raw_hfEle_eLCeL9, *raw_hfEle_eL1eL9, *raw_hfEle_eL9eL25;
    TH1D *raw_hfEle_eS1eS9, *raw_hfEle_eS9eS25;
    TH1D *raw_hfEle_eS1eL1, *raw_hfEle_eS9eL9, *raw_hfEle_eS25eL25;
    TH1D *raw_hfEle_cut2D;
    TH2D *raw_hfEleEtaX_elePtY;
    TH2D *raw_hfEle_etX_etaY, *raw_hfEle_etX_etaY_ID;
    P4Hists *ele_p4, *hfEle_p4;
    TH1D *ele_dPhi, *ele_dEta, *ele_dPt, *ele_dE;
    TH2D *ele_etaX_phiY, *hfEle_etaX_phiY;
    TH2D *eleEtaX_hfEleEtaY, *elePtX_hfElePtY;
    TH2D *hfEleEnergyX_llMassY;
    P4Hists *ele_Mcut_p4, *hfEle_Mcut_p4;
    TH2D *hfEle_etX_corrY;
    P4Hists *eleGD_p4, *eleBD_p4, *hfEleGD_p4, *hfEleBD_p4;
    TH1D *eleGD_dR, *eleBD_dR, *hfEleGD_dR, *hfEleBD_dR;
    TH1D *eleGD_dPt, *eleBD_dPt, *hfEleGD_dPt, *hfEleBD_dPt;
    TH1D *eleGD_dE, *eleBD_dE, *hfEleGD_dE, *hfEleBD_dE;
    TH2D *hfEleGD_etX_corrY;
    TH2D *hfEleGD_etX_etaY, *hfEleGD_etX_etaY_ID;
    TH2D *hltEle27_ptX_etaY_T, *hltEle27_ptX_etaY_TnP;
    
};

#endif

#ifndef Electrons_GUARD
#define Electrons_GUARD

#include "MyAnalyses/CommonTools/interface/P4Hists.h"
#include "MyAnalyses/ZprimeSearch/bin/TreeWrapper.h"
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include <string>
#include <vector>

class Electrons {
    
public:
    
    // constructor
    Electrons(bool _isMC, bool _doQCD, bool _doHist);
    
    // destructor
    ~Electrons();
    
    // get first electron
    const EventInfo::Electron& GetEle1() const;
    
    // get second electron
    const EventInfo::Electron& GetEle2() const;
    
    // test kinematics
    bool IsGoodP4(const TLorentzVector &p4, double etaMin, double etaMax) const;
    
    // get status
    bool IsGood() const;
    
    // reset information
    void Reset();
    
    // select electrons
    void SelectElectrons(const std::vector<EventInfo::Electron> &electrons,
                         double eta1Min, double eta1Max, double eta2Min, double eta2Max);
    
    // fill eventCount histogram
    void FillEventCount(std::string tag = "");
    
    // fill trigger histograms
    void FillTriggerHistograms(const EventInfo &event, double weight);
    
    // fill raw histograms
    void FillRawHistograms(const EventInfo &event, double weight);
    
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
    bool IsGoodID(const EventInfo::Electron &ele) const;
    
    // private objects
    bool isMC, doQCD, doHist;
    
    // selected electrons
    const EventInfo::Electron *electron[2];
    unsigned int eleSize;
    bool goodElectronsNoID, goodElectrons;
    
    // histograms
    TH1D *eventCount;
    TH1D *raw_ele_size;
    TH2D *raw_eleSizeX_pvSizeY;
    TH1D *raw_ele_dEtaIn[2], *raw_ele_dPhiIn[2], *raw_ele_sigmaIEtaIEta[2], *raw_ele_HoE[2];
    TH1D *raw_ele_d0[2], *raw_ele_dZ[2], *raw_ele_IoEmIoP[2], *raw_ele_iso[2];
    P4Hists *ele_p4[2];
    TH1D *ele_dPhi, *ele_dEta, *ele_dPt, *ele_dE;
    TH2D *ele1_etaX_phiY, *ele2_etaX_phiY;
    TH2D *ele_eta1X_eta2Y, *ele_pt1X_pt2Y;
    P4Hists *ele1_Mcut_p4, *ele2_Mcut_p4;
    P4Hists *eleGD_p4[2], *eleBD_p4[2];
    TH1D *eleGD_dR[2], *eleBD_dR[2], *eleGD_dPt[2], *eleBD_dPt[2], *eleGD_dE[2], *eleBD_dE[2];
    TH2D *hltEle27_ptX_etaY_T, *hltEle8_ptX_etaY_TnP, *hltEle17_ptX_etaY_TnP;
    
};

#endif

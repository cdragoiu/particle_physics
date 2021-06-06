#ifndef Zprime_GUARD
#define Zprime_GUARD

#include "MyAnalyses/CommonTools/interface/P4Hists.h"
#include "MyAnalyses/ZprimeSearch/bin/Corrections.h"
#include "MyAnalyses/ZprimeSearch/bin/Electrons.h"
#include "MyAnalyses/ZprimeSearch/bin/GenInfo.h"
#include "MyAnalyses/ZprimeSearch/bin/HfElectrons.h"
#include "MyAnalyses/ZprimeSearch/bin/TreeWrapper.h"
#include "RooUnfold/src/RooUnfoldResponse.h"
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include <string>
#include <vector>

class Zprime {
    
public:
    
    // constructor
    Zprime(std::string _runType);
    
    // destructor
    ~Zprime();
    
    // process event
    void ProcessEvent(const EventInfo &event, const Corrections &corr, int ord=0);
    
    // write histograms and responses
    void Write();
    
    // compute cos(theta)
    static double GetCosTheta(const TLorentzVector &ele, const TLorentzVector &pos);
    
    // AFB mass binning
    static double massBins_ele[17], massBins_hf[10];
    
private:
    
    // test trigger conditions
    bool IsGoodTrigger(const std::vector<bool> &triggers) const;
    
    // test for a good primary vertex
    bool IsGoodPV(const std::vector<EventInfo::PrimVertex> &primVertices) const;
    
    // reconstruct the ll object
    void GetLL(const TLorentzVector &p1, const TLorentzVector &p2);
    
    // fill raw histograms
    void FillRawHistograms(const EventInfo &event, double weight);
    
    // fill MC only histograms
    void FillMcHistograms(const EventInfo::Particle &par, double weight);
    
    // fill histograms
    void FillHistograms(const EventInfo &event, double weight);
    
    // fill responses
    void FillResponses(double weight);
    
    // book histograms
    void BookHistograms();
    
    // book responses
    void BookResponses();
    
    // private objects
    std::string runType;
    bool isMC, isHF, doUNF, doSYS;
    
    // selected objects
    GenInfo *genInfo;
    Electrons *electrons;
    HfElectrons *hfElectrons;
    TLorentzVector ll;
    
    // histograms
    TH1D *raw_pv_size_true, *raw_pv_size;
    P4Hists *llGD_p4, *llBD_p4;
    TH1D *llGD_dR, *llBD_dR, *llGD_dPt, *llBD_dPt, *llGD_dE, *llBD_dE;
    P4Hists *ll_p4;
    TH2D *ll_massX_cosThetaY_0y1, *ll_massX_cosThetaY_1y1p25, *ll_massX_cosThetaY_1p25y1p5;
    TH2D *ll_massX_cosThetaY_1p5y2p4, *ll_massX_cosThetaY_2p4y5;
    TH2D *llMassX_eleEtaY, *llMassX_posEtaY;
    TH2D *llRapidityX_pvSizeY, *llMassX_pvSizeY;
    TH2D *llRapidityX_pfJetSizeY;
    TH2D *llRapidityX_pfMetY, *llRapidityX_pfMetOverSumEtY;
    TH2D *llRapidityX_eleDphiY, *llMassX_eleDphiY, *llRapidityX_eleDetaY, *llRapidityX_eleEtaProdY;
    TH2D *llRapidityX_llMassY;
    
    // responses
    RooUnfoldResponse *response_0y1, *response_1y1p25, *response_1p25y1p5;
    RooUnfoldResponse *response_1p5y2p4, *response_2p4y5;
    
};

#endif

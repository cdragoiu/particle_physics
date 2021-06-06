#ifndef GenInfo_GUARD
#define GenInfo_GUARD

#include "MyAnalyses/CommonTools/interface/P4Hists.h"
#include "MyAnalyses/ZprimeSearch/bin/TreeWrapper.h"
#include "TH1.h"
#include "TH2.h"
#include <vector>

class GenInfo {
    
public:
    
    // constructor
    GenInfo(bool _doHist);
    
    // destructor
    ~GenInfo();
    
    // get electron
    const EventInfo::Particle& GetEle() const;
    
    // get positron
    const EventInfo::Particle& GetPos() const;
    
    // get Z boson
    const EventInfo::Particle& GetZ() const;
    
    // get status
    bool IsGood() const;
    
    // reset information
    void Reset();
    
    // select particles
    void SelectParticles(const std::vector<EventInfo::Particle> &particles);
    
    // fill raw histograms
    void FillRawHistograms(double weight);
    
    // fill histograms
    void FillHistograms(double weight);
    
    // write histograms
    void Write();
    
private:
    
    // book histograms
    void BookHistograms();
    
    // private objects
    bool doHist;
    
    // selected particles
    const EventInfo::Particle *ele, *pos, *Z;
    bool goodParticles;
    
    // histograms
    P4Hists *gen_ele_p4, *gen_pos_p4, *gen_Z_p4;
    TH1D *gen_dPhi, *gen_dEta, *gen_dPt, *gen_dE;
    TH2D *gen_ele_etaX_phiY, *gen_pos_etaX_phiY;
    TH2D *gen_ll_massX_cosThetaY_0y1, *gen_ll_massX_cosThetaY_1y1p25;
    TH2D *gen_ll_massX_cosThetaY_1p25y1p5, *gen_ll_massX_cosThetaY_1p5y2p4;
    TH2D *gen_ll_massX_cosThetaY_2p4y5;
    
};

#endif

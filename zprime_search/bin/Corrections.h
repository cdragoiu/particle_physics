#ifndef Corrections_GUARD
#define Corrections_GUARD

#include "TLorentzVector.h"
#include <vector>

// scale factor structure
struct ScaleFactor {
    ScaleFactor(double _ptMin, double _ptMax, double _etaMin, double _etaMax,
                double _corr, double _err):
                ptMin(_ptMin), ptMax(_ptMax), etaMin(_etaMin), etaMax(_etaMax),
                corr(_corr), err(_err) {}
    double ptMin, ptMax, etaMin, etaMax, corr, err;
};

class Corrections {
    
public:
    
    // constructor
    Corrections();
    
    // destructor
    ~Corrections();
    
    // get the Ele8 leg of Ele17_Ele8 scale factor
    double GetTrgEle8SF(const TLorentzVector &p4, int ord=0) const;
    
    // get the Ele17 leg of Ele17_Ele8 scale factor
    double GetTrgEle17SF(const TLorentzVector &p4, int ord=0) const;
    
    // get the Ele27_WP80 scale factor
    double GetTrgEle27SF(const TLorentzVector &p4, int ord=0) const;
    
    // get the tight electron ID scale factor
    double GetTightEleIdSF(const TLorentzVector &p4, int ord=0) const;
    
    // get the HF electron ID scale factor
    double GetHfEleIdSF(const TLorentzVector &p4, int ord=0) const;
    
    // get the HF electron Eta scale factor
    double GetHfEleEtaSF(const TLorentzVector &p4, int ord=0) const;
    
    // get the HF electron energy correction (Data)
    double GetHfEleEnergyCorrData(double et, int ord=0) const;
    
    // get the HF electron energy correction (MC)
    double GetHfEleEnergyCorrMC(double et, int ord=0);
    
private:
    
    // Ele8 leg of Ele17_Ele8 scale factors initialization
    void InitTrgEle8SF();
    
    // Ele17 leg of Ele17_Ele8 scale factors initialization
    void InitTrgEle17SF();
    
    // Ele27_WP80 scale factors initialization
    void InitTrgEle27SF();
    
    // tight electron ID scale factors initialization
    void InitTightEleIdSF();
    
    // HF electron ID scale factors initialization
    void InitHfEleIdSF();
    
    // HF electron Eta scale factors initialization
    void InitHfEleEtaSF();
    
    // scale factors
    std::vector<ScaleFactor> trgEle8SF, trgEle17SF, trgEle27SF;
    std::vector<ScaleFactor> tightEleIdSF, hfEleIdSF, hfEleEtaSF;
    
};

#endif

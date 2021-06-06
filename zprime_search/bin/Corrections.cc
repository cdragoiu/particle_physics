#include "MyAnalyses/ZprimeSearch/bin/Corrections.h"
#include <cmath>
#include <cstdlib>

// constructor -------------------------------------------------------------------------------------
Corrections::Corrections() {
    InitTrgEle8SF();
    InitTrgEle17SF();
    InitTrgEle27SF();
    InitTightEleIdSF();
    InitHfEleIdSF();
    InitHfEleEtaSF();
}

// destructor --------------------------------------------------------------------------------------
Corrections::~Corrections() {
}

// get the Ele8 leg of Ele17_Ele8 scale factor -----------------------------------------------------
double Corrections::GetTrgEle8SF(const TLorentzVector &p4, int ord) const {
    for (unsigned int i = 0; i != trgEle8SF.size(); ++i) {
        if (p4.Pt() < trgEle8SF.at(i).ptMin) continue;
        if (p4.Pt() > trgEle8SF.at(i).ptMax) continue;
        if (fabs(p4.Eta()) < trgEle8SF.at(i).etaMin) continue;
        if (fabs(p4.Eta()) > trgEle8SF.at(i).etaMax) continue;
        int c(0);
        if (ord) c = (abs(ord)==i+1)*ord/abs(ord);
        return trgEle8SF.at(i).corr + c*trgEle8SF.at(i).err;
    }
    return 1.0;
}

// get the Ele17 leg of Ele17_Ele8 scale factor ----------------------------------------------------
double Corrections::GetTrgEle17SF(const TLorentzVector &p4, int ord) const {
    for (unsigned int i = 0; i != trgEle17SF.size(); ++i) {
        if (p4.Pt() < trgEle17SF.at(i).ptMin) continue;
        if (p4.Pt() > trgEle17SF.at(i).ptMax) continue;
        if (fabs(p4.Eta()) < trgEle17SF.at(i).etaMin) continue;
        if (fabs(p4.Eta()) > trgEle17SF.at(i).etaMax) continue;
        int c(0);
        if (ord) c = (abs(ord)==i+1)*ord/abs(ord);
        return trgEle17SF.at(i).corr + c*trgEle17SF.at(i).err;
    }
    return 1.0;
}

// get the Ele27_WP80 scale factor -----------------------------------------------------------------
double Corrections::GetTrgEle27SF(const TLorentzVector &p4, int ord) const {
    for (unsigned int i = 0; i != trgEle27SF.size(); ++i) {
        if (p4.Pt() < trgEle27SF.at(i).ptMin) continue;
        if (p4.Pt() > trgEle27SF.at(i).ptMax) continue;
        if (fabs(p4.Eta()) < trgEle27SF.at(i).etaMin) continue;
        if (fabs(p4.Eta()) > trgEle27SF.at(i).etaMax) continue;
        int c(0);
        if (ord) c = (abs(ord)==i+1)*ord/abs(ord);
        return trgEle27SF.at(i).corr + c*trgEle27SF.at(i).err;
    }
    return 1.0;
}

// get the tight electron ID scale factor ----------------------------------------------------------
double Corrections::GetTightEleIdSF(const TLorentzVector &p4, int ord) const {
    for (unsigned int i = 0; i != tightEleIdSF.size(); ++i) {
        if (p4.Pt() < tightEleIdSF.at(i).ptMin) continue;
        if (p4.Pt() > tightEleIdSF.at(i).ptMax) continue;
        if (fabs(p4.Eta()) < tightEleIdSF.at(i).etaMin) continue;
        if (fabs(p4.Eta()) > tightEleIdSF.at(i).etaMax) continue;
        int c(0);
        if (ord) c = (abs(ord)==i+1)*ord/abs(ord);
        return tightEleIdSF.at(i).corr + c*tightEleIdSF.at(i).err;
    }
    return 1.0;
}

// get the HF electron ID scale factor -------------------------------------------------------------
double Corrections::GetHfEleIdSF(const TLorentzVector &p4, int ord) const {
    for (unsigned int i = 0; i != hfEleIdSF.size(); ++i) {
        if (p4.Et() < hfEleIdSF.at(i).ptMin) continue;
        if (p4.Et() > hfEleIdSF.at(i).ptMax) continue;
        if (fabs(p4.Eta()) < hfEleIdSF.at(i).etaMin) continue;
        if (fabs(p4.Eta()) > hfEleIdSF.at(i).etaMax) continue;
        int c(0);
        if (ord) c = (abs(ord)==i+1)*ord/abs(ord);
        double err = sqrt(pow(0.023,2)+pow(hfEleIdSF.at(i).err,2));
        return hfEleIdSF.at(i).corr + c*err;
    }
    return 1.0;
}

// get the HF electron Eta scale factor ------------------------------------------------------------
double Corrections::GetHfEleEtaSF(const TLorentzVector &p4, int ord) const {
    for (unsigned int i = 0; i != hfEleEtaSF.size(); ++i) {
        if (p4.Eta() < hfEleEtaSF.at(i).etaMin) continue;
        if (p4.Eta() > hfEleEtaSF.at(i).etaMax) continue;
        int c(0);
        if (ord) c = (abs(ord)==i+1)*ord/abs(ord);
        return hfEleEtaSF.at(i).corr + c*hfEleEtaSF.at(i).err;
    }
    return 1.0;
}

// get the HF electron energy correction (Data) ----------------------------------------------------
double Corrections::GetHfEleEnergyCorrData(double et, int ord) const {
    if (et < 20.0) return 1.0;
    if (et > 120.0) return 1.0;
    double p[6] = {0.919, 0.001477, 15.54, 0.9255, 0.02056, 0.210};
    double e[6] = {0.011, 0.000032,  0.25, 0.0092, 0.00031, 0.011};
    for (unsigned int i = 0; i != 6; ++i) {
        int c(0);
        if (ord) c = (abs(ord)==i+1)*ord/abs(ord);
        p[i] += c*e[i];
    }
    double corr = p[0]*exp(-p[1]*pow(et-p[2],2))+p[3]*erf(p[4]*et-p[5]);
    return corr;
}

// get the HF electron energy correction (MC) ------------------------------------------------------
double Corrections::GetHfEleEnergyCorrMC(double et, int ord) {
    if (et < 20.0) return 1.0;
    if (et > 120.0) return 1.0;
    double p[7] = {1.3538, 0.2721, 0.001142, 26.28, 0.4085, 0.001104, 20.47};
    double e[7] = {0.0046, 0.0047, 0.000081,  0.66, 0.0077, 0.000030,  0.32};
    for (unsigned int i = 0; i != 7; ++i) {
        int c(0);
        if (ord) c = (abs(ord)==i+1)*ord/abs(ord);
        p[i] += c*e[i];
    }
    double corr = p[0]-p[1]*exp(-p[2]*pow(et-p[3],2))+p[4]*erf(-p[5]*pow(et-p[6],2));
    return corr;
}

// Ele8 leg of Ele17_Ele8 scale factors initialization ---------------------------------------------
void Corrections::InitTrgEle8SF() {
    trgEle8SF.clear();
    trgEle8SF.push_back(ScaleFactor(20.0, 30.0, 0.0, 0.8, 0.9961, 0.0004));
    trgEle8SF.push_back(ScaleFactor(20.0, 30.0, 0.8, 1.442, 0.9934, 0.0005));
    trgEle8SF.push_back(ScaleFactor(20.0, 30.0, 1.556, 2.0, 0.9973, 0.0006));
    trgEle8SF.push_back(ScaleFactor(20.0, 30.0, 2.0, 2.5, 0.9975, 0.0008));
    trgEle8SF.push_back(ScaleFactor(30.0, 40.0, 0.0, 0.8, 0.9962, 0.0002));
    trgEle8SF.push_back(ScaleFactor(30.0, 40.0, 0.8, 1.442, 0.9943, 0.0002));
    trgEle8SF.push_back(ScaleFactor(30.0, 40.0, 1.556, 2.0, 0.9981, 0.0003));
    trgEle8SF.push_back(ScaleFactor(30.0, 40.0, 2.0, 2.5, 0.9988, 0.0005));
    trgEle8SF.push_back(ScaleFactor(40.0, 50.0, 0.0, 0.8, 0.9948, 0.0002));
    trgEle8SF.push_back(ScaleFactor(40.0, 50.0, 0.8, 1.442, 0.9952, 0.0002));
    trgEle8SF.push_back(ScaleFactor(40.0, 50.0, 1.556, 2.0, 0.9982, 0.0003));
    trgEle8SF.push_back(ScaleFactor(40.0, 50.0, 2.0, 2.5, 0.9986, 0.0005));
    trgEle8SF.push_back(ScaleFactor(50.0, 70.0, 0.0, 0.8, 0.9963, 0.0011));
    trgEle8SF.push_back(ScaleFactor(50.0, 70.0, 0.8, 1.442, 0.9959, 0.0013));
    trgEle8SF.push_back(ScaleFactor(50.0, 70.0, 1.556, 2.0, 0.9986, 0.0017));
    trgEle8SF.push_back(ScaleFactor(50.0, 70.0, 2.0, 2.5, 1.0041, 0.0026));
    trgEle8SF.push_back(ScaleFactor(70.0, 250.0, 0.0, 0.8, 0.9938, 0.0018));
    trgEle8SF.push_back(ScaleFactor(70.0, 250.0, 0.8, 1.442, 0.9976, 0.0019));
    trgEle8SF.push_back(ScaleFactor(70.0, 250.0, 1.556, 2.0, 1.0020, 0.0035));
    trgEle8SF.push_back(ScaleFactor(70.0, 250.0, 2.0, 2.5, 0.9906, 0.0042));
    return;
}

// Ele17 leg of Ele17_Ele8 scale factors initialization --------------------------------------------
void Corrections::InitTrgEle17SF() {
    trgEle17SF.clear();
    trgEle17SF.push_back(ScaleFactor(20.0, 30.0, 0.0, 0.8, 0.9935, 0.0004));
    trgEle17SF.push_back(ScaleFactor(20.0, 30.0, 0.8, 1.442, 0.9892, 0.0005));
    trgEle17SF.push_back(ScaleFactor(20.0, 30.0, 1.556, 2.0, 0.9973, 0.0005));
    trgEle17SF.push_back(ScaleFactor(20.0, 30.0, 2.0, 2.5, 0.9955, 0.0008));
    trgEle17SF.push_back(ScaleFactor(30.0, 40.0, 0.0, 0.8, 0.9953, 0.0002));
    trgEle17SF.push_back(ScaleFactor(30.0, 40.0, 0.8, 1.442, 0.9938, 0.0002));
    trgEle17SF.push_back(ScaleFactor(30.0, 40.0, 1.556, 2.0, 0.9978, 0.0002));
    trgEle17SF.push_back(ScaleFactor(30.0, 40.0, 2.0, 2.5, 0.9985, 0.0004));
    trgEle17SF.push_back(ScaleFactor(40.0, 50.0, 0.0, 0.8, 0.9947, 0.0002));
    trgEle17SF.push_back(ScaleFactor(40.0, 50.0, 0.8, 1.442, 0.9944, 0.0002));
    trgEle17SF.push_back(ScaleFactor(40.0, 50.0, 1.556, 2.0, 0.9977, 0.0003));
    trgEle17SF.push_back(ScaleFactor(40.0, 50.0, 2.0, 2.5, 0.9985, 0.0005));
    trgEle17SF.push_back(ScaleFactor(50.0, 70.0, 0.0, 0.8, 0.9943, 0.0010));
    trgEle17SF.push_back(ScaleFactor(50.0, 70.0, 0.8, 1.442, 0.9959, 0.0012));
    trgEle17SF.push_back(ScaleFactor(50.0, 70.0, 1.556, 2.0, 0.9968, 0.0014));
    trgEle17SF.push_back(ScaleFactor(50.0, 70.0, 2.0, 2.5, 1.0015, 0.0024));
    trgEle17SF.push_back(ScaleFactor(70.0, 250.0, 0.0, 0.8, 0.9933, 0.0017));
    trgEle17SF.push_back(ScaleFactor(70.0, 250.0, 0.8, 1.442, 0.9968, 0.0012));
    trgEle17SF.push_back(ScaleFactor(70.0, 250.0, 1.556, 2.0, 0.9971, 0.0024));
    trgEle17SF.push_back(ScaleFactor(70.0, 250.0, 2.0, 2.5, 0.9906, 0.0042));
    return;
}

// Ele27_WP80 scale factors initialization ---------------------------------------------------------
void Corrections::InitTrgEle27SF() {
    trgEle27SF.clear();
    trgEle27SF.push_back(ScaleFactor(30.0, 40.0, 0.0, 0.8, 0.9740, 0.0004));
    trgEle27SF.push_back(ScaleFactor(30.0, 40.0, 0.8, 1.442, 0.9595, 0.0006));
    trgEle27SF.push_back(ScaleFactor(30.0, 40.0, 1.556, 2.0, 1.0144, 0.0019));
    trgEle27SF.push_back(ScaleFactor(30.0, 40.0, 2.0, 2.5, 1.0022, 0.0025));
    trgEle27SF.push_back(ScaleFactor(40.0, 50.0, 0.0, 0.8, 0.9827, 0.0004));
    trgEle27SF.push_back(ScaleFactor(40.0, 50.0, 0.8, 1.442, 0.9818, 0.0005));
    trgEle27SF.push_back(ScaleFactor(40.0, 50.0, 1.556, 2.0, 1.0136, 0.0020));
    trgEle27SF.push_back(ScaleFactor(40.0, 50.0, 2.0, 2.5, 1.0084, 0.0029));
    trgEle27SF.push_back(ScaleFactor(50.0, 70.0, 0.0, 0.8, 0.9825, 0.0021));
    trgEle27SF.push_back(ScaleFactor(50.0, 70.0, 0.8, 1.442, 0.9871, 0.0027));
    trgEle27SF.push_back(ScaleFactor(50.0, 70.0, 1.556, 2.0, 1.0009, 0.0104));
    trgEle27SF.push_back(ScaleFactor(50.0, 70.0, 2.0, 2.5, 1.0120, 0.0155));
    trgEle27SF.push_back(ScaleFactor(70.0, 250.0, 0.0, 0.8, 0.9834, 0.0035));
    trgEle27SF.push_back(ScaleFactor(70.0, 250.0, 0.8, 1.442, 0.9897, 0.0042));
    trgEle27SF.push_back(ScaleFactor(70.0, 250.0, 1.556, 2.0, 0.9711, 0.0230));
    trgEle27SF.push_back(ScaleFactor(70.0, 250.0, 2.0, 2.5, 1.0374, 0.0410));
    return;
}

// tight electron ID scale factors initialization --------------------------------------------------
void Corrections::InitTightEleIdSF() {
    tightEleIdSF.clear();
    tightEleIdSF.push_back(ScaleFactor(10.0, 15.0, 0.0, 0.8, 0.827, 0.021));
    tightEleIdSF.push_back(ScaleFactor(10.0, 15.0, 0.8, 1.442, 0.948, 0.024));
    tightEleIdSF.push_back(ScaleFactor(10.0, 15.0, 1.556, 2.0, 0.854, 0.048));
    tightEleIdSF.push_back(ScaleFactor(10.0, 15.0, 2.0, 2.5, 1.007, 0.047));
    tightEleIdSF.push_back(ScaleFactor(15.0, 20.0, 0.0, 0.8, 0.924, 0.01));
    tightEleIdSF.push_back(ScaleFactor(15.0, 20.0, 0.8, 1.442, 0.932, 0.012));
    tightEleIdSF.push_back(ScaleFactor(15.0, 20.0, 1.556, 2.0, 0.853, 0.022));
    tightEleIdSF.push_back(ScaleFactor(15.0, 20.0, 2.0, 2.5, 0.903, 0.029));
    tightEleIdSF.push_back(ScaleFactor(20.0, 30.0, 0.0, 0.8, 0.96, 0.003));
    tightEleIdSF.push_back(ScaleFactor(20.0, 30.0, 0.8, 1.442, 0.936, 0.004));
    tightEleIdSF.push_back(ScaleFactor(20.0, 30.0, 1.556, 2.0, 0.879, 0.007));
    tightEleIdSF.push_back(ScaleFactor(20.0, 30.0, 2.0, 2.5, 0.974, 0.004));
    tightEleIdSF.push_back(ScaleFactor(30.0, 40.0, 0.0, 0.8, 0.978, 0.001));
    tightEleIdSF.push_back(ScaleFactor(30.0, 40.0, 0.8, 1.442, 0.958, 0.002));
    tightEleIdSF.push_back(ScaleFactor(30.0, 40.0, 1.556, 2.0, 0.909, 0.003));
    tightEleIdSF.push_back(ScaleFactor(30.0, 40.0, 2.0, 2.5, 0.987, 0.004));
    tightEleIdSF.push_back(ScaleFactor(40.0, 50.0, 0.0, 0.8, 0.981, 0.001));
    tightEleIdSF.push_back(ScaleFactor(40.0, 50.0, 0.8, 1.442, 0.969, 0.001));
    tightEleIdSF.push_back(ScaleFactor(40.0, 50.0, 1.556, 2.0, 0.942, 0.002));
    tightEleIdSF.push_back(ScaleFactor(40.0, 50.0, 2.0, 2.5, 0.991, 0.003));
    tightEleIdSF.push_back(ScaleFactor(50.0, 200.0, 0.0, 0.8, 0.982, 0.002));
    tightEleIdSF.push_back(ScaleFactor(50.0, 200.0, 0.8, 1.442, 0.969, 0.002));
    tightEleIdSF.push_back(ScaleFactor(50.0, 200.0, 1.556, 2.0, 0.957, 0.004));
    tightEleIdSF.push_back(ScaleFactor(50.0, 200.0, 2.0, 2.5, 0.999, 0.005));
    return;
}

// HF electron ID scale factors initialization -----------------------------------------------------
void Corrections::InitHfEleIdSF() {
    hfEleIdSF.clear();
    hfEleIdSF.push_back(ScaleFactor(20.0, 30.0, 3.0, 3.5, 0.911, 0.010));
    hfEleIdSF.push_back(ScaleFactor(20.0, 30.0, 3.5, 4.0, 0.872, 0.023));
    hfEleIdSF.push_back(ScaleFactor(20.0, 30.0, 4.0, 5.0, 0.815, 0.041));
    hfEleIdSF.push_back(ScaleFactor(30.0, 40.0, 3.0, 3.5, 0.912, 0.005));
    hfEleIdSF.push_back(ScaleFactor(30.0, 40.0, 3.5, 4.0, 0.884, 0.013));
    hfEleIdSF.push_back(ScaleFactor(30.0, 40.0, 4.0, 5.0, 0.834, 0.057));
    hfEleIdSF.push_back(ScaleFactor(40.0, 50.0, 3.0, 3.5, 0.897, 0.010));
    hfEleIdSF.push_back(ScaleFactor(40.0, 50.0, 3.5, 5.0, 0.864, 0.041));
    hfEleIdSF.push_back(ScaleFactor(50.0, 70.0, 3.0, 5.0, 0.846, 0.031));
    hfEleIdSF.push_back(ScaleFactor(70.0, 150.0, 3.0, 5.0, 0.906, 0.071));
    return;
}

// HF electron Eta scale factors initialization ----------------------------------------------------
void Corrections::InitHfEleEtaSF() {
    hfEleEtaSF.clear();
    hfEleEtaSF.push_back(ScaleFactor(0.0, 0.0, -5.0, -4.0, 1.324, 0.095));
    hfEleEtaSF.push_back(ScaleFactor(0.0, 0.0, -4.0, -3.8, 1.394, 0.070));
    hfEleEtaSF.push_back(ScaleFactor(0.0, 0.0, -3.8, -3.6, 1.150, 0.039));
    hfEleEtaSF.push_back(ScaleFactor(0.0, 0.0, -3.6, -3.4, 1.005, 0.023));
    hfEleEtaSF.push_back(ScaleFactor(0.0, 0.0, -3.4, -3.2, 1.170, 0.017));
    hfEleEtaSF.push_back(ScaleFactor(0.0, 0.0, -3.2, -3.0, 0.754, 0.012));
    hfEleEtaSF.push_back(ScaleFactor(0.0, 0.0, 3.0, 3.2, 0.752, 0.011));
    hfEleEtaSF.push_back(ScaleFactor(0.0, 0.0, 3.2, 3.4, 1.245, 0.018));
    hfEleEtaSF.push_back(ScaleFactor(0.0, 0.0, 3.4, 3.6, 1.104, 0.025));
    hfEleEtaSF.push_back(ScaleFactor(0.0, 0.0, 3.6, 3.8, 1.212, 0.042));
    hfEleEtaSF.push_back(ScaleFactor(0.0, 0.0, 3.8, 4.0, 1.364, 0.072));
    hfEleEtaSF.push_back(ScaleFactor(0.0, 0.0, 4.0, 5.0, 1.368, 0.098));
    return;
}

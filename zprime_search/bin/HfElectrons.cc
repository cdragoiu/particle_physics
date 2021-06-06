#include "MyAnalyses/CommonTools/interface/Utils.h"
#include "MyAnalyses/ZprimeSearch/bin/HfElectrons.h"
#include "TMath.h"
#include <algorithm>
#include <cmath>

// constructor -------------------------------------------------------------------------------------
HfElectrons::HfElectrons(bool _isMC, bool _doQCD, bool _doHist) :
                         isMC(_isMC), doQCD(_doQCD), doHist(_doHist), electron(0), hfElectron(0),
                         eleSize(0), hfEleSize(0), goodElectronsNoID(false),
                         goodEleID(false), goodHfEleID(false), goodElectrons(false) {
    if (doHist) BookHistograms();
}

// destructor --------------------------------------------------------------------------------------
HfElectrons::~HfElectrons() {
    if (doHist) {
        delete eventCount;
        delete raw_ele_size;
        delete raw_hfEle_size;
        delete raw_eleSizeX_pvSizeY;
        delete raw_hfEleSizeX_pvSizeY;
        delete raw_ele_dEtaIn;
        delete raw_ele_dPhiIn;
        delete raw_ele_sigmaIEtaIEta;
        delete raw_ele_HoE;
        delete raw_ele_d0;
        delete raw_ele_dZ;
        delete raw_ele_IoEmIoP;
        delete raw_ele_iso;
        delete raw_hfEle_eCore;
        delete raw_hfEle_eLong1x1;
        delete raw_hfEle_eLong3x3;
        delete raw_hfEle_eLong5x5;
        delete raw_hfEle_eShort1x1;
        delete raw_hfEle_eShort3x3;
        delete raw_hfEle_eShort5x5;
        delete raw_hfEle_eLCeL9;
        delete raw_hfEle_eL1eL9;
        delete raw_hfEle_eL9eL25;
        delete raw_hfEle_eS1eS9;
        delete raw_hfEle_eS9eS25;
        delete raw_hfEle_eS1eL1;
        delete raw_hfEle_eS9eL9;
        delete raw_hfEle_eS25eL25;
        delete raw_hfEle_cut2D;
        delete raw_hfEleEtaX_elePtY;
        delete raw_hfEle_etX_etaY;
        delete raw_hfEle_etX_etaY_ID;
        delete ele_p4;
        delete hfEle_p4;
        delete ele_dPhi;
        delete ele_dEta;
        delete ele_dPt;
        delete ele_dE;
        delete ele_etaX_phiY;
        delete hfEle_etaX_phiY;
        delete eleEtaX_hfEleEtaY;
        delete elePtX_hfElePtY;
        delete hfEleEnergyX_llMassY;
        delete ele_Mcut_p4;
        delete hfEle_Mcut_p4;
        delete hfEle_etX_corrY;
        if (isMC) {
            delete eleGD_p4;
            delete eleBD_p4;
            delete hfEleGD_p4;
            delete hfEleBD_p4;
            delete eleGD_dR;
            delete eleBD_dR;
            delete hfEleGD_dR;
            delete hfEleBD_dR;
            delete eleGD_dPt;
            delete eleBD_dPt;
            delete hfEleGD_dPt;
            delete hfEleBD_dPt;
            delete eleGD_dE;
            delete eleBD_dE;
            delete hfEleGD_dE;
            delete hfEleBD_dE;
            delete hfEleGD_etX_corrY;
            delete hfEleGD_etX_etaY;
            delete hfEleGD_etX_etaY_ID;
        }
        delete hltEle27_ptX_etaY_T;
        delete hltEle27_ptX_etaY_TnP;
    }
}

// get electron ------------------------------------------------------------------------------------
const EventInfo::Electron& HfElectrons::GetEle() const {
    return *electron;
}

// get HF electron ---------------------------------------------------------------------------------
const EventInfo::HfElectron& HfElectrons::GetHfEle() const {
    return *hfElectron;
}

// test electron kinematics ------------------------------------------------------------------------
bool HfElectrons::IsGoodEleP4(const TLorentzVector &p4, double etaMin, double etaMax) const {
    if (fabs(p4.Eta()) < etaMin || fabs(p4.Eta()) > etaMax) return false;
    if (fabs(p4.Eta()) > 1.442 && fabs(p4.Eta()) < 1.556) return false;
    if (p4.Pt() < 30.0) return false;
    return true;
}

// test HF electron kinematics ---------------------------------------------------------------------
bool HfElectrons::IsGoodHfEleP4(const TLorentzVector &p4) const {
    if (fabs(p4.Eta()) < 3.0 || fabs(p4.Eta()) > 5.0) return false;
    if (p4.Et() < 20.0) return false;
    return true;
}

// get status --------------------------------------------------------------------------------------
bool HfElectrons::IsGood() const {
    return goodElectrons;
}

// reset information -------------------------------------------------------------------------------
void HfElectrons::Reset() {
    electron = 0;
    hfElectron = 0;
    eleSize = 0;
    hfEleSize = 0;
    goodElectronsNoID = false;
    goodEleID = false;
    goodHfEleID = false;
    goodElectrons = false;
    return;
}

// select electrons --------------------------------------------------------------------------------
void HfElectrons::SelectElectrons(const std::vector<EventInfo::Electron> &electrons,
                                  double etaMin, double etaMax,
                                  const std::vector<EventInfo::HfElectron> &hfElectrons) {
    Reset();
    for (unsigned int i = 0; i != electrons.size(); ++i) {
        if (!IsGoodEleP4(electrons.at(i).p4, 0.0, 2.4)) continue;
        ++eleSize;
    }
    for (unsigned int i = 0; i != hfElectrons.size(); ++i) {
        if (!IsGoodHfEleP4(hfElectrons.at(i).p4)) continue;
        ++hfEleSize;
    }
    if (eleSize != 1 || hfEleSize != 1) return;
    if (doHist) FillEventCount("goodEleSize");
    if (electrons.at(0).p4.Eta() * hfElectrons.at(0).p4.Eta() < 0.0) return;
    if (doHist) FillEventCount("goodEleEta");
    if (fabs(electrons.at(0).p4.DeltaPhi(hfElectrons.at(0).p4)) < 2.0/3.0*TMath::Pi()) return;
    if (doHist) FillEventCount("goodElePhi");
    if (!IsGoodEleP4(electrons.at(0).p4, etaMin, etaMax)) return;
    if (!IsGoodHfEleP4(hfElectrons.at(0).p4)) return;
    if (doHist) FillEventCount("goodEleNoID");
    goodElectronsNoID = true;
    if (IsGoodEleID(electrons.at(0))) {
        if (doHist) FillEventCount("goodEleID");
        goodEleID = true;
    }
    if (IsGoodHfEleID(hfElectrons.at(0))) {
        if (doHist) FillEventCount("goodHfEleID");
        goodHfEleID = true;
    }
    if (!goodEleID || !goodHfEleID) return;
    if (doHist) FillEventCount("goodEle");
    electron = &electrons.at(0);
    hfElectron = &hfElectrons.at(0);
    goodElectrons = true;
    return;
}

// fill eventCount histogram -----------------------------------------------------------------------
void HfElectrons::FillEventCount(std::string tag) {
    if (tag == "goodTrigger") eventCount->Fill(1.5);
    else if (tag == "goodPV") eventCount->Fill(2.5);
    else if (tag == "goodEleSize") eventCount->Fill(3.5);
    else if (tag == "goodEleEta") eventCount->Fill(4.5);
    else if (tag == "goodElePhi") eventCount->Fill(5.5);
    else if (tag == "goodEleNoID") eventCount->Fill(6.5);
    else if (tag == "goodEleID") eventCount->Fill(7.5);
    else if (tag == "goodHfEleID") eventCount->Fill(8.5);
    else if (tag == "goodEle") eventCount->Fill(9.5);
    else eventCount->Fill(0.5);
    return;
}

// fill trigger histograms -------------------------------------------------------------------------
void HfElectrons::FillTriggerHistograms(const EventInfo &event, double weight) {
    if (event.electrons.size() < 2) return;
    if (!IsGoodEleID(event.electrons.at(0))) return;
    if (!IsGoodEleID(event.electrons.at(1))) return;
    if (event.electrons.at(0).charge == event.electrons.at(1).charge) return;
    TLorentzVector p4 = event.electrons.at(0).p4 + event.electrons.at(1).p4;
    if (p4.M() < 80.0 || p4.M() > 100.0) return;
    if (event.electrons.at(0).ele27wp80) {
        hltEle27_ptX_etaY_T->Fill(event.electrons.at(1).p4.Pt(),
                                  fabs(event.electrons.at(1).p4.Eta()), weight);
        if (event.electrons.at(1).ele27wp80) {
            hltEle27_ptX_etaY_TnP->Fill(event.electrons.at(1).p4.Pt(),
                                        fabs(event.electrons.at(1).p4.Eta()), weight);
        }
    }
    return;
}

// fill raw histograms -----------------------------------------------------------------------------
void HfElectrons::FillRawHistograms(const EventInfo &event, double weight) {
    raw_ele_size->Fill(eleSize, weight);
    raw_eleSizeX_pvSizeY->Fill(eleSize, event.primVertices.size(), weight);
    raw_hfEle_size->Fill(hfEleSize, weight);
    raw_hfEleSizeX_pvSizeY->Fill(hfEleSize, event.primVertices.size(), weight);
    if (!goodElectronsNoID) return;
    TLorentzVector p4 = event.electrons.at(0).p4 + event.hfElectrons.at(0).p4;
    if (p4.M() < 80.0 || p4.M() > 100.0) return;
    raw_ele_dEtaIn->Fill(event.electrons.at(0).dEtaIn, weight);
    raw_ele_dPhiIn->Fill(event.electrons.at(0).dPhiIn, weight);
    raw_ele_sigmaIEtaIEta->Fill(event.electrons.at(0).sigmaIEtaIEta, weight);
    raw_ele_HoE->Fill(event.electrons.at(0).HoE, weight);
    raw_ele_d0->Fill(event.electrons.at(0).d0, weight);
    raw_ele_dZ->Fill(event.electrons.at(0).dZ, weight);
    raw_ele_IoEmIoP->Fill(event.electrons.at(0).IoEmIoP, weight);
    raw_ele_iso->Fill(event.electrons.at(0).iso, weight);
    if (!goodEleID) return;
    raw_hfEle_eCore->Fill(event.hfElectrons.at(0).eCore, weight);
    raw_hfEle_eLong1x1->Fill(event.hfElectrons.at(0).eLong1x1, weight);
    raw_hfEle_eLong3x3->Fill(event.hfElectrons.at(0).eLong3x3, weight);
    raw_hfEle_eLong5x5->Fill(event.hfElectrons.at(0).eLong5x5, weight);
    raw_hfEle_eShort1x1->Fill(event.hfElectrons.at(0).eShort1x1, weight);
    raw_hfEle_eShort3x3->Fill(event.hfElectrons.at(0).eShort3x3, weight);
    raw_hfEle_eShort5x5->Fill(event.hfElectrons.at(0).eShort5x5, weight);
    double eLCeL9 = event.hfElectrons.at(0).eCore / event.hfElectrons.at(0).eLong3x3;
    raw_hfEle_eLCeL9->Fill(eLCeL9, weight);
    double eL1eL9 = event.hfElectrons.at(0).eLong1x1 / event.hfElectrons.at(0).eLong3x3;
    raw_hfEle_eL1eL9->Fill(eL1eL9, weight);
    double eL9eL25 = event.hfElectrons.at(0).eLong3x3 / event.hfElectrons.at(0).eLong5x5;
    raw_hfEle_eL9eL25->Fill(eL9eL25, weight);
    double eS1eS9 = event.hfElectrons.at(0).eShort1x1 / event.hfElectrons.at(0).eShort3x3;
    raw_hfEle_eS1eS9->Fill(eS1eS9, weight);
    double eS9eS25 = event.hfElectrons.at(0).eShort3x3 / event.hfElectrons.at(0).eShort5x5;
    raw_hfEle_eS9eS25->Fill(eS9eS25, weight);
    double eS1eL1 = event.hfElectrons.at(0).eShort1x1 / event.hfElectrons.at(0).eLong1x1;
    raw_hfEle_eS1eL1->Fill(eS1eL1, weight);
    double eS9eL9 = event.hfElectrons.at(0).eShort3x3 / event.hfElectrons.at(0).eLong3x3;
    raw_hfEle_eS9eL9->Fill(eS9eL9, weight);
    double eS25eL25 = event.hfElectrons.at(0).eShort5x5 / event.hfElectrons.at(0).eLong5x5;
    raw_hfEle_eS25eL25->Fill(eS25eL25, weight);
    double cut2D = event.hfElectrons.at(0).eCOREe9 - 1.125 * event.hfElectrons.at(0).eSeL;
    raw_hfEle_cut2D->Fill(cut2D, weight);
    raw_hfEle_etX_etaY->Fill(event.hfElectrons.at(0).p4.Et(),
                             fabs(event.hfElectrons.at(0).p4.Eta()), weight);
    if (!goodHfEleID) return;
    raw_hfEle_etX_etaY_ID->Fill(event.hfElectrons.at(0).p4.Et(),
                                fabs(event.hfElectrons.at(0).p4.Eta()), weight);
    raw_hfEleEtaX_elePtY->Fill(event.hfElectrons.at(0).p4.Eta(),
                               event.electrons.at(0).p4.Pt(), weight);
    return;
}

// fill raw MC only histograms ---------------------------------------------------------------------
void HfElectrons::FillRawMcHistograms(const EventInfo &event, const EventInfo::Particle &ele,
                                      const EventInfo::Particle &pos, double weight) {
    for (unsigned int i = 0; i != event.hfElectrons.size(); ++i) {
        double dR_ele = event.hfElectrons.at(i).p4.DeltaR(ele.p4);
        double dR_pos = event.hfElectrons.at(i).p4.DeltaR(pos.p4);
        if (std::min(dR_ele, dR_pos) > 0.3) continue;
        if (!IsGoodHfEleP4(event.hfElectrons.at(i).p4)) return;
        hfEleGD_etX_etaY->Fill(event.hfElectrons.at(i).p4.Et(),
                               fabs(event.hfElectrons.at(i).p4.Eta()), weight);
        if (!IsGoodHfEleID(event.hfElectrons.at(i))) return;
        hfEleGD_etX_etaY_ID->Fill(event.hfElectrons.at(i).p4.Et(),
                                  fabs(event.hfElectrons.at(i).p4.Eta()), weight);
        return;
    }
    return;
}

// fill MC only histograms -------------------------------------------------------------------------
void HfElectrons::FillMcHistograms(const EventInfo::Particle &ele,
                                   const EventInfo::Particle &pos, double weight) {
    double dR_ele = electron->p4.DeltaR(ele.p4);
    double dR_pos = electron->p4.DeltaR(pos.p4);
    const EventInfo::Particle *gen = dR_ele < dR_pos ? &ele : &pos;
    double dR = std::min(dR_ele, dR_pos);
    if (dR < 0.3) {
        eleGD_p4->Fill(electron->p4, weight);
        eleGD_dR->Fill(dR, weight);
        eleGD_dPt->Fill(electron->p4.Pt() / gen->p4.Pt(), weight);
        eleGD_dE->Fill(electron->p4.Energy() / gen->p4.Energy(), weight);
    }
    else {
        eleBD_p4->Fill(electron->p4, weight);
        eleBD_dR->Fill(dR, weight);
        eleBD_dPt->Fill(electron->p4.Pt() / gen->p4.Pt(), weight);
        eleBD_dE->Fill(electron->p4.Energy() / gen->p4.Energy(), weight);
    }
    dR_ele = hfElectron->p4.DeltaR(ele.p4);
    dR_pos = hfElectron->p4.DeltaR(pos.p4);
    gen = dR_ele < dR_pos ? &ele : &pos;
    dR = std::min(dR_ele, dR_pos);
    if (dR < 0.3) {
        hfEleGD_p4->Fill(hfElectron->p4, weight);
        hfEleGD_dR->Fill(dR, weight);
        hfEleGD_dPt->Fill(hfElectron->p4.Pt() / gen->p4.Pt(), weight);
        hfEleGD_dE->Fill(hfElectron->p4.Energy() / gen->p4.Energy(), weight);
    }
    else {
        hfEleBD_p4->Fill(hfElectron->p4, weight);
        hfEleBD_dR->Fill(dR, weight);
        hfEleBD_dPt->Fill(hfElectron->p4.Pt() / gen->p4.Pt(), weight);
        hfEleBD_dE->Fill(hfElectron->p4.Energy() / gen->p4.Energy(), weight);
    }
    TLorentzVector p4 = electron->p4 + hfElectron->p4;
    if (p4.M() > 70.0 && p4.M() < 110.0 && dR < 0.3) {
        hfEleGD_etX_corrY->Fill(hfElectron->p4.Et(),
                                gen->p4.Energy() / hfElectron->p4.Energy(), weight);
    }
    return;
}

// fill histograms ---------------------------------------------------------------------------------
void HfElectrons::FillHistograms(double weight) {
    ele_p4->Fill(electron->p4, weight);
    hfEle_p4->Fill(hfElectron->p4, weight);
    ele_dPhi->Fill(electron->p4.DeltaPhi(hfElectron->p4), weight);
    ele_dEta->Fill(electron->p4.Eta() - hfElectron->p4.Eta(), weight);
    ele_dPt->Fill(electron->p4.Pt() / hfElectron->p4.Pt(), weight);
    ele_dE->Fill(electron->p4.Energy() / hfElectron->p4.Energy(), weight);
    ele_etaX_phiY->Fill(electron->p4.Eta(), electron->p4.Phi(), weight);
    hfEle_etaX_phiY->Fill(hfElectron->p4.Eta(), hfElectron->p4.Phi(), weight);
    eleEtaX_hfEleEtaY->Fill(electron->p4.Eta(), hfElectron->p4.Eta(), weight);
    elePtX_hfElePtY->Fill(electron->p4.Pt(), hfElectron->p4.Pt(), weight);
    TLorentzVector p4 = electron->p4 + hfElectron->p4;
    hfEleEnergyX_llMassY->Fill(hfElectron->p4.Energy(), p4.M(), weight);
    if (p4.M() < 70.0 || p4.M() > 110.0) return;
    double M = 91.1876;
    double emEta = electron->p4.Eta();
    double hfEta = hfElectron->p4.Eta();
    double emPhi = electron->p4.Phi();
    double hfPhi = hfElectron->p4.Phi();
    double emE = electron->p4.Energy();
    double predE = (pow(M,2)*cosh(emEta)*cosh(hfEta))/(2*emE*(cosh(emEta-hfEta)-cos(emPhi-hfPhi)));
    hfEle_etX_corrY->Fill(hfElectron->p4.Et(), predE/hfElectron->p4.Energy(), weight);
    if (p4.M() < 80.0 || p4.M() > 100.0) return;
    ele_Mcut_p4->Fill(electron->p4, weight);
    hfEle_Mcut_p4->Fill(hfElectron->p4, weight);
    return;
}

// write histograms --------------------------------------------------------------------------------
void HfElectrons::Write() {
    eventCount->Write();
    raw_ele_size->Write();
    raw_hfEle_size->Write();
    raw_eleSizeX_pvSizeY->Write();
    raw_hfEleSizeX_pvSizeY->Write();
    raw_ele_dEtaIn->Write();
    raw_ele_dPhiIn->Write();
    raw_ele_sigmaIEtaIEta->Write();
    raw_ele_HoE->Write();
    raw_ele_d0->Write();
    raw_ele_dZ->Write();
    raw_ele_IoEmIoP->Write();
    raw_ele_iso->Write();
    raw_hfEle_eCore->Write();
    raw_hfEle_eLong1x1->Write();
    raw_hfEle_eLong3x3->Write();
    raw_hfEle_eLong5x5->Write();
    raw_hfEle_eShort1x1->Write();
    raw_hfEle_eShort3x3->Write();
    raw_hfEle_eShort5x5->Write();
    raw_hfEle_eLCeL9->Write();
    raw_hfEle_eL1eL9->Write();
    raw_hfEle_eL9eL25->Write();
    raw_hfEle_eS1eS9->Write();
    raw_hfEle_eS9eS25->Write();
    raw_hfEle_eS1eL1->Write();
    raw_hfEle_eS9eL9->Write();
    raw_hfEle_eS25eL25->Write();
    raw_hfEle_cut2D->Write();
    raw_hfEleEtaX_elePtY->Write();
    raw_hfEle_etX_etaY->Write();
    raw_hfEle_etX_etaY_ID->Write();
    ele_p4->Write();
    hfEle_p4->Write();
    ele_dPhi->Write();
    ele_dEta->Write();
    ele_dPt->Write();
    ele_dE->Write();
    ele_etaX_phiY->Write();
    hfEle_etaX_phiY->Write();
    eleEtaX_hfEleEtaY->Write();
    elePtX_hfElePtY->Write();
    hfEleEnergyX_llMassY->Write();
    ele_Mcut_p4->Write();
    hfEle_Mcut_p4->Write();
    hfEle_etX_corrY->Write();
    if (isMC) {
        eleGD_p4->Write();
        eleBD_p4->Write();
        hfEleGD_p4->Write();
        hfEleBD_p4->Write();
        eleGD_dR->Write();
        eleBD_dR->Write();
        hfEleGD_dR->Write();
        hfEleBD_dR->Write();
        eleGD_dPt->Write();
        eleBD_dPt->Write();
        hfEleGD_dPt->Write();
        hfEleBD_dPt->Write();
        eleGD_dE->Write();
        eleBD_dE->Write();
        hfEleGD_dE->Write();
        hfEleBD_dE->Write();
        hfEleGD_etX_corrY->Write();
        hfEleGD_etX_etaY->Write();
        hfEleGD_etX_etaY_ID->Write();
    }
    hltEle27_ptX_etaY_T->Write();
    hltEle27_ptX_etaY_TnP->Write();
    return;
}

// book histograms ---------------------------------------------------------------------------------
void HfElectrons::BookHistograms() {
    eventCount = Book1D("eventCount", 10, 0.0, 10.0);
    raw_ele_size = Book1D("raw_ele_size", 20, -0.5, 19.5);
    raw_hfEle_size = Book1D("raw_hfEle_size", 20, -0.5, 19.5);
    raw_eleSizeX_pvSizeY = Book2D("raw_eleSizeX_pvSizeY", 20, -0.5, 19.5, 80, -0.5, 79.5);
    raw_hfEleSizeX_pvSizeY = Book2D("raw_hfEleSizeX_pvSizeY", 20, -0.5, 19.5, 80, -0.5, 79.5);
    raw_ele_dEtaIn = Book1D("raw_ele_dEtaIn", 1000, -0.5, 0.5);
    raw_ele_dPhiIn = Book1D("raw_ele_dPhiIn", 1000, -0.5, 0.5);
    raw_ele_sigmaIEtaIEta = Book1D("raw_ele_sigmaIEtaIEta", 1000, 0.0, 0.1);
    raw_ele_HoE = Book1D("raw_ele_HoE", 1000, 0.0, 1.0);
    raw_ele_d0 = Book1D("raw_ele_d0", 1000, -10.0, 10.0);
    raw_ele_dZ = Book1D("raw_ele_dZ", 1000, -40.0, 40.0);
    raw_ele_IoEmIoP = Book1D("raw_ele_IoEmIoP", 1000, -1.0, 1.0);
    raw_ele_iso = Book1D("raw_ele_iso", 5000, 0.0, 5.0);
    raw_hfEle_eCore = Book1D("raw_hfEle_eCore", 2000, 0.0, 2000.0);
    raw_hfEle_eLong1x1 = Book1D("raw_hfEle_eLong1x1", 2000, 0.0, 2000.0);
    raw_hfEle_eLong3x3 = Book1D("raw_hfEle_eLong3x3", 2000, 0.0, 2000.0);
    raw_hfEle_eLong5x5 = Book1D("raw_hfEle_eLong5x5", 2000, 0.0, 2000.0);
    raw_hfEle_eShort1x1 = Book1D("raw_hfEle_eShort1x1", 2000, 0.0, 2000.0);
    raw_hfEle_eShort3x3 = Book1D("raw_hfEle_eShort3x3", 2000, 0.0, 2000.0);
    raw_hfEle_eShort5x5 = Book1D("raw_hfEle_eShort5x5", 2000, 0.0, 2000.0);
    raw_hfEle_eLCeL9 = Book1D("raw_hfEle_eLCeL9", 1100, 0.0, 1.1);
    raw_hfEle_eL1eL9 = Book1D("raw_hfEle_eL1eL9", 1100, 0.0, 1.1);
    raw_hfEle_eL9eL25 = Book1D("raw_hfEle_eL9eL25", 1100, 0.0, 1.1);
    raw_hfEle_eS1eS9 = Book1D("raw_hfEle_eS1eS9", 1100, 0.0, 1.1);
    raw_hfEle_eS9eS25 = Book1D("raw_hfEle_eS9eS25", 1100, 0.0, 1.1);
    raw_hfEle_eS1eL1 = Book1D("raw_hfEle_eS1eL1", 1100, 0.0, 1.1);
    raw_hfEle_eS9eL9 = Book1D("raw_hfEle_eS9eL9", 1100, 0.0, 1.1);
    raw_hfEle_eS25eL25 = Book1D("raw_hfEle_eS25eL25", 1100, 0.0, 1.1);
    raw_hfEle_cut2D = Book1D("raw_hfEle_cut2D", 1100, 0.0, 1.1);
    raw_hfEleEtaX_elePtY = Book2D("raw_hfEleEtaX_elePtY", 200, -5.0, 5.0, 2000, 0.0, 2000.0);
    raw_hfEle_etX_etaY = Book2D("raw_hfEle_etX_etaY", 500, 0.0, 500.0, 240, 2.8, 5.2);
    raw_hfEle_etX_etaY_ID = Book2D("raw_hfEle_etX_etaY_ID", 500, 0.0, 500.0, 240, 2.8, 5.2);
    ele_p4 = new P4Hists("ele");
    hfEle_p4 = new P4Hists("hfEle");
    ele_dPhi = Book1D("ele_dPhi", 400, -TMath::Pi(), TMath::Pi());
    ele_dEta = Book1D("ele_dEta", 400, -10.0, 10.0);
    ele_dPt = Book1D("ele_dPt", 400, 0.0, 4.0);
    ele_dE = Book1D("ele_dE", 400, 0.0, 4.0);
    ele_etaX_phiY = Book2D("ele_etaX_phiY", 200, -5.0, 5.0, 200, -TMath::Pi(), TMath::Pi());
    hfEle_etaX_phiY = Book2D("hfEle_etaX_phiY", 200, -5.0, 5.0, 200, -TMath::Pi(), TMath::Pi());
    eleEtaX_hfEleEtaY = Book2D("eleEtaX_hfEleEtaY", 200, -5.0, 5.0, 200, -5.0, 5.0);
    elePtX_hfElePtY = Book2D("elePtX_hfElePtY", 2000, 0.0, 2000.0, 2000, 0.0, 2000.0);
    hfEleEnergyX_llMassY = Book2D("hfEleEnergyX_llMassY", 2000, 0.0, 2000.0, 1000, 0.0, 1000.0);
    ele_Mcut_p4 = new P4Hists("ele_Mcut");
    hfEle_Mcut_p4 = new P4Hists("hfEle_Mcut");
    hfEle_etX_corrY = Book2D("hfEle_etX_corrY", 200, 0.0, 200.0, 200, 0.0, 2.0);
    if (isMC) {
        eleGD_p4 = new P4Hists("eleGD");
        eleBD_p4 = new P4Hists("eleBD");
        hfEleGD_p4 = new P4Hists("hfEleGD");
        hfEleBD_p4 = new P4Hists("hfEleBD");
        eleGD_dR = Book1D("eleGD_dR", 400, 0.0, 4.0);
        eleBD_dR = Book1D("eleBD_dR", 400, 0.0, 4.0);
        hfEleGD_dR = Book1D("hfEleGD_dR", 400, 0.0, 4.0);
        hfEleBD_dR = Book1D("hfEleBD_dR", 400, 0.0, 4.0);
        eleGD_dPt = Book1D("eleGD_dPt", 400, 0.0, 4.0);
        eleBD_dPt = Book1D("eleBD_dPt", 400, 0.0, 4.0);
        hfEleGD_dPt = Book1D("hfEleGD_dPt", 400, 0.0, 4.0);
        hfEleBD_dPt = Book1D("hfEleBD_dPt", 400, 0.0, 4.0);
        eleGD_dE = Book1D("eleGD_dE", 400, 0.0, 4.0);
        eleBD_dE = Book1D("eleBD_dE", 400, 0.0, 4.0);
        hfEleGD_dE = Book1D("hfEleGD_dE", 400, 0.0, 4.0);
        hfEleBD_dE = Book1D("hfEleBD_dE", 400, 0.0, 4.0);
        hfEleGD_etX_corrY = Book2D("hfEleGD_etX_corrY", 200, 0.0, 200.0, 200, 0.0, 2.0);
        hfEleGD_etX_etaY = Book2D("hfEleGD_etX_etaY", 500, 0.0, 500.0, 240, 2.8, 5.2);
        hfEleGD_etX_etaY_ID = Book2D("hfEleGD_etX_etaY_ID", 500, 0.0, 500.0, 240, 2.8, 5.2);
    }
    hltEle27_ptX_etaY_T = Book2D("hltEle27_ptX_etaY_T", 1000, 0.0, 1000.0, 5000, 0.0, 5.0);
    hltEle27_ptX_etaY_TnP = Book2D("hltEle27_ptX_etaY_TnP", 1000, 0.0, 1000.0, 5000, 0.0, 5.0);
    return;
}

// test electron ID --------------------------------------------------------------------------------
bool HfElectrons::IsGoodEleID(const EventInfo::Electron &ele) const {
    if (doQCD) {
        bool goodTightId = (
                            (fabs(ele.etaSC) < 1.442)     &&
                            (fabs(ele.dEtaIn) < 0.004)    &&
                            (fabs(ele.dPhiIn) < 0.03)     &&
                            (ele.sigmaIEtaIEta < 0.01)    &&
                            (ele.HoE < 0.12)              &&
                            (fabs(ele.d0) < 0.02)         &&
                            (fabs(ele.dZ) < 0.1)          &&
                            (fabs(ele.IoEmIoP) < 0.05)    &&
                            (ele.iso > 0.2)               && // <--
                            (!ele.hasConversion)          &&
                            (ele.missingHits <= 0)
                           )                              ||
                           (
                            (fabs(ele.etaSC) > 1.566)     &&
                            (fabs(ele.dEtaIn) < 0.005)    &&
                            (fabs(ele.dPhiIn) < 0.02)     &&
                            (ele.sigmaIEtaIEta < 0.03)    &&
                            (ele.HoE < 0.10)              &&
                            (fabs(ele.d0) < 0.02)         &&
                            (fabs(ele.dZ) < 0.1)          &&
                            (fabs(ele.IoEmIoP) < 0.05)    &&
                            (ele.iso > 0.2)               && // <--
                            (!ele.hasConversion)          &&
                            (ele.missingHits <= 0)
                           );
        if (!goodTightId) return false;
    }
    else {
        if (!ele.goodTightId) return false;
    }
    return true;
}

// test HF electron ID -----------------------------------------------------------------------------
bool HfElectrons::IsGoodHfEleID(const EventInfo::HfElectron &hfEle) const {
    if ((hfEle.eLong3x3 / hfEle.eLong5x5) < 0.96) return false;
    if ((hfEle.eCOREe9 - 1.125 * hfEle.eSeL) < 0.4) return false;
    return true;
}

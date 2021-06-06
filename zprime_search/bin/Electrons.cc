#include "MyAnalyses/CommonTools/interface/Utils.h"
#include "MyAnalyses/ZprimeSearch/bin/Electrons.h"
#include "TMath.h"
#include <algorithm>
#include <cmath>

// constructor -------------------------------------------------------------------------------------
Electrons::Electrons(bool _isMC, bool _doQCD, bool _doHist) :
                     isMC(_isMC), doQCD(_doQCD), doHist(_doHist), electron(), eleSize(0),
                     goodElectronsNoID(false), goodElectrons(false) {
    if (doHist) BookHistograms();
}

// destructor --------------------------------------------------------------------------------------
Electrons::~Electrons() {
    if (doHist) {
        delete eventCount;
        delete raw_ele_size;
        delete raw_eleSizeX_pvSizeY;
        for (unsigned int i = 0; i != 2; ++i) {
            delete raw_ele_dEtaIn[i];
            delete raw_ele_dPhiIn[i];
            delete raw_ele_sigmaIEtaIEta[i];
            delete raw_ele_HoE[i];
            delete raw_ele_d0[i];
            delete raw_ele_dZ[i];
            delete raw_ele_IoEmIoP[i];
            delete raw_ele_iso[i];
            delete ele_p4[i];
        }
        delete ele_dPhi;
        delete ele_dEta;
        delete ele_dPt;
        delete ele_dE;
        delete ele1_etaX_phiY;
        delete ele2_etaX_phiY;
        delete ele_eta1X_eta2Y;
        delete ele_pt1X_pt2Y;
        delete ele1_Mcut_p4;
        delete ele2_Mcut_p4;
        if (isMC) {
            for (unsigned int i = 0; i != 2; ++i) {
                delete eleGD_p4[i];
                delete eleBD_p4[i];
                delete eleGD_dR[i];
                delete eleBD_dR[i];
                delete eleGD_dPt[i];
                delete eleBD_dPt[i];
                delete eleGD_dE[i];
                delete eleBD_dE[i];
            }
        }
        delete hltEle27_ptX_etaY_T;
        delete hltEle8_ptX_etaY_TnP;
        delete hltEle17_ptX_etaY_TnP;
    }
}

// get first electron ------------------------------------------------------------------------------
const EventInfo::Electron& Electrons::GetEle1() const {
    return *electron[0];
}

// get second electron -----------------------------------------------------------------------------
const EventInfo::Electron& Electrons::GetEle2() const {
    return *electron[1];
}

// test kinematics ---------------------------------------------------------------------------------
bool Electrons::IsGoodP4(const TLorentzVector &p4, double etaMin, double etaMax) const {
    if (fabs(p4.Eta()) < etaMin || fabs(p4.Eta()) > etaMax) return false;
    if (fabs(p4.Eta()) > 1.442 && fabs(p4.Eta()) < 1.556) return false;
    if (p4.Pt() < 20.0) return false;
    return true;
}

// get status --------------------------------------------------------------------------------------
bool Electrons::IsGood() const {
    return goodElectrons;
}

// reset information -------------------------------------------------------------------------------
void Electrons::Reset() {
    electron[0] = 0;
    electron[1] = 0;
    eleSize = 0;
    goodElectronsNoID = false;
    goodElectrons = false;
    return;
}

// select electrons --------------------------------------------------------------------------------
void Electrons::SelectElectrons(const std::vector<EventInfo::Electron> &electrons,
                                double eta1Min, double eta1Max, double eta2Min, double eta2Max) {
    Reset();
    for (unsigned int i = 0; i != electrons.size(); ++i) {
        if (!IsGoodP4(electrons.at(i).p4, 0.0, 2.4)) continue;
        ++eleSize;
    }
    if (eleSize != 2) return;
    if (doHist) FillEventCount("goodEleSize");
    if (electrons.at(0).charge == electrons.at(1).charge) return;
    if (doHist) FillEventCount("goodEleCharge");
    if (!IsGoodP4(electrons.at(0).p4, eta1Min, eta1Max)) return;
    if (!IsGoodP4(electrons.at(1).p4, eta2Min, eta2Max)) return;
    if (doHist) FillEventCount("goodEleNoID");
    goodElectronsNoID = true;
    if (!IsGoodID(electrons.at(0))) return;
    if (!IsGoodID(electrons.at(1))) return;
    if (doHist) FillEventCount("goodEle");
    electron[0] = &electrons.at(0);
    electron[1] = &electrons.at(1);
    goodElectrons = true;
    return;
}

// fill eventCount histogram -----------------------------------------------------------------------
void Electrons::FillEventCount(std::string tag) {
    if (tag == "goodTrigger") eventCount->Fill(1.5);
    else if (tag == "goodPV") eventCount->Fill(2.5);
    else if (tag == "goodEleSize") eventCount->Fill(3.5);
    else if (tag == "goodEleCharge") eventCount->Fill(4.5);
    else if (tag == "goodEleNoID") eventCount->Fill(5.5);
    else if (tag == "goodEle") eventCount->Fill(6.5);
    else eventCount->Fill(0.5);
    return;
}

// fill trigger histograms -------------------------------------------------------------------------
void Electrons::FillTriggerHistograms(const EventInfo &event, double weight) {
    if (event.electrons.size() < 2) return;
    if (!IsGoodID(event.electrons.at(0))) return;
    if (!IsGoodID(event.electrons.at(1))) return;
    if (event.electrons.at(0).charge == event.electrons.at(1).charge) return;
    TLorentzVector p4 = event.electrons.at(0).p4 + event.electrons.at(1).p4;
    if (p4.M() < 80.0 || p4.M() > 100.0) return;
    if (event.electrons.at(0).ele27wp80) {
        hltEle27_ptX_etaY_T->Fill(event.electrons.at(1).p4.Pt(),
                                  fabs(event.electrons.at(1).p4.Eta()), weight);
        if (event.electrons.at(1).ele17ele8_f2) {
            hltEle8_ptX_etaY_TnP->Fill(event.electrons.at(1).p4.Pt(),
                                       fabs(event.electrons.at(1).p4.Eta()), weight);
        }
        if (event.electrons.at(1).ele17ele8_f1) {
            hltEle17_ptX_etaY_TnP->Fill(event.electrons.at(1).p4.Pt(),
                                        fabs(event.electrons.at(1).p4.Eta()), weight);
        }
    }
    return;
}

// fill raw histograms -----------------------------------------------------------------------------
void Electrons::FillRawHistograms(const EventInfo &event, double weight) {
    raw_ele_size->Fill(eleSize, weight);
    raw_eleSizeX_pvSizeY->Fill(eleSize, event.primVertices.size(), weight);
    if (!goodElectronsNoID) return;
    TLorentzVector p4 = event.electrons.at(0).p4 + event.electrons.at(1).p4;
    if (p4.M() < 80.0 || p4.M() > 100.0) return;
    for (unsigned int i = 0; i != 2; ++i) {
        raw_ele_dEtaIn[i]->Fill(event.electrons.at(i).dEtaIn, weight);
        raw_ele_dPhiIn[i]->Fill(event.electrons.at(i).dPhiIn, weight);
        raw_ele_sigmaIEtaIEta[i]->Fill(event.electrons.at(i).sigmaIEtaIEta, weight);
        raw_ele_HoE[i]->Fill(event.electrons.at(i).HoE, weight);
        raw_ele_d0[i]->Fill(event.electrons.at(i).d0, weight);
        raw_ele_dZ[i]->Fill(event.electrons.at(i).dZ, weight);
        raw_ele_IoEmIoP[i]->Fill(event.electrons.at(i).IoEmIoP, weight);
        raw_ele_iso[i]->Fill(event.electrons.at(i).iso, weight);
    }
    return;
}

// fill MC only histograms -------------------------------------------------------------------------
void Electrons::FillMcHistograms(const EventInfo::Particle &ele,
                                 const EventInfo::Particle &pos, double weight) {
    for (unsigned int i = 0; i != 2; ++i) {
        double dR_ele = electron[i]->p4.DeltaR(ele.p4);
        double dR_pos = electron[i]->p4.DeltaR(pos.p4);
        double dR = std::min(dR_ele, dR_pos);
        double pt = dR_ele < dR_pos ? ele.p4.Pt() : pos.p4.Pt();
        double energy = dR_ele < dR_pos ? ele.p4.Energy() : pos.p4.Energy();
        if (dR < 0.3) {
            eleGD_p4[i]->Fill(electron[i]->p4, weight);
            eleGD_dR[i]->Fill(dR, weight);
            eleGD_dPt[i]->Fill(electron[i]->p4.Pt() / pt, weight);
            eleGD_dE[i]->Fill(electron[i]->p4.Energy() / energy, weight);
        }
        else {
            eleBD_p4[i]->Fill(electron[i]->p4, weight);
            eleBD_dR[i]->Fill(dR, weight);
            eleBD_dPt[i]->Fill(electron[i]->p4.Pt() / pt, weight);
            eleBD_dE[i]->Fill(electron[i]->p4.Energy() / energy, weight);
        }
    }
    return;
}

// fill histograms ---------------------------------------------------------------------------------
void Electrons::FillHistograms(double weight) {
    ele_p4[0]->Fill(electron[0]->p4, weight);
    ele_p4[1]->Fill(electron[1]->p4, weight);
    ele_dPhi->Fill(electron[0]->p4.DeltaPhi(electron[1]->p4), weight);
    ele_dEta->Fill(electron[0]->p4.Eta() - electron[1]->p4.Eta(), weight);
    ele_dPt->Fill(electron[0]->p4.Pt() / electron[1]->p4.Pt(), weight);
    ele_dE->Fill(electron[0]->p4.Energy() / electron[1]->p4.Energy(), weight);
    ele1_etaX_phiY->Fill(electron[0]->p4.Eta(), electron[0]->p4.Phi(), weight);
    ele2_etaX_phiY->Fill(electron[1]->p4.Eta(), electron[1]->p4.Phi(), weight);
    ele_eta1X_eta2Y->Fill(electron[0]->p4.Eta(), electron[1]->p4.Eta(), weight);
    ele_pt1X_pt2Y->Fill(electron[0]->p4.Pt(), electron[1]->p4.Pt(), weight);
    TLorentzVector p4 = electron[0]->p4 + electron[1]->p4;
    if (p4.M() < 80.0 || p4.M() > 100.0) return;
    ele1_Mcut_p4->Fill(electron[0]->p4, weight);
    ele2_Mcut_p4->Fill(electron[1]->p4, weight);
    return;
}

// write histograms --------------------------------------------------------------------------------
void Electrons::Write() {
    eventCount->Write();
    raw_ele_size->Write();
    raw_eleSizeX_pvSizeY->Write();
    for (unsigned int i = 0; i != 2; ++i) {
        raw_ele_dEtaIn[i]->Write();
        raw_ele_dPhiIn[i]->Write();
        raw_ele_sigmaIEtaIEta[i]->Write();
        raw_ele_HoE[i]->Write();
        raw_ele_d0[i]->Write();
        raw_ele_dZ[i]->Write();
        raw_ele_IoEmIoP[i]->Write();
        raw_ele_iso[i]->Write();
        ele_p4[i]->Write();
    }
    ele_dPhi->Write();
    ele_dEta->Write();
    ele_dPt->Write();
    ele_dE->Write();
    ele1_etaX_phiY->Write();
    ele2_etaX_phiY->Write();
    ele_eta1X_eta2Y->Write();
    ele_pt1X_pt2Y->Write();
    ele1_Mcut_p4->Write();
    ele2_Mcut_p4->Write();
    if (isMC) {
        for (unsigned int i = 0; i != 2; ++i) {
            eleGD_p4[i]->Write();
            eleBD_p4[i]->Write();
            eleGD_dR[i]->Write();
            eleBD_dR[i]->Write();
            eleGD_dPt[i]->Write();
            eleBD_dPt[i]->Write();
            eleGD_dE[i]->Write();
            eleBD_dE[i]->Write();
        }
    }
    hltEle27_ptX_etaY_T->Write();
    hltEle8_ptX_etaY_TnP->Write();
    hltEle17_ptX_etaY_TnP->Write();
    return;
}

// book histograms ---------------------------------------------------------------------------------
void Electrons::BookHistograms() {
    eventCount = Book1D("eventCount", 7, 0.0, 7.0);
    raw_ele_size = Book1D("raw_ele_size", 20, -0.5, 19.5);
    raw_eleSizeX_pvSizeY = Book2D("raw_eleSizeX_pvSizeY", 20, -0.5, 19.5, 80, -0.5, 79.5);
    for (unsigned int i = 0; i != 2; ++i) {
        std::string tag = ToStr(i+1);
        raw_ele_dEtaIn[i] = Book1D("raw_ele" + tag + "_dEtaIn", 1000, -0.5, 0.5);
        raw_ele_dPhiIn[i] = Book1D("raw_ele" + tag + "_dPhiIn", 1000, -0.5, 0.5);
        raw_ele_sigmaIEtaIEta[i] = Book1D("raw_ele" + tag + "_sigmaIEtaIEta", 1000, 0.0, 0.1);
        raw_ele_HoE[i] = Book1D("raw_ele" + tag + "_HoE", 1000, 0.0, 1.0);
        raw_ele_d0[i] = Book1D("raw_ele" + tag + "_d0", 1000, -10.0, 10.0);
        raw_ele_dZ[i] = Book1D("raw_ele" + tag + "_dZ", 1000, -40.0, 40.0);
        raw_ele_IoEmIoP[i] = Book1D("raw_ele" + tag + "_IoEmIoP", 1000, -1.0, 1.0);
        raw_ele_iso[i] = Book1D("raw_ele" + tag + "_iso", 5000, 0.0, 5.0);
        ele_p4[i] = new P4Hists("ele" + tag);
    }
    ele_dPhi = Book1D("ele_dPhi", 400, -TMath::Pi(), TMath::Pi());
    ele_dEta = Book1D("ele_dEta", 400, -10.0, 10.0);
    ele_dPt = Book1D("ele_dPt", 400, 0.0, 4.0);
    ele_dE = Book1D("ele_dE", 400, 0.0, 4.0);
    ele1_etaX_phiY = Book2D("ele1_etaX_phiY", 200, -5.0, 5.0, 200, -TMath::Pi(), TMath::Pi());
    ele2_etaX_phiY = Book2D("ele2_etaX_phiY", 200, -5.0, 5.0, 200, -TMath::Pi(), TMath::Pi());
    ele_eta1X_eta2Y = Book2D("ele_eta1X_eta2Y", 200, -5.0, 5.0, 200, -5.0, 5.0);
    ele_pt1X_pt2Y = Book2D("ele_pt1X_pt2Y", 2000, 0.0, 2000.0, 2000, 0.0, 2000.0);
    ele1_Mcut_p4 = new P4Hists("ele1_Mcut");
    ele2_Mcut_p4 = new P4Hists("ele2_Mcut");
    if (isMC) {
        for (unsigned int i = 0; i != 2; ++i) {
            std::string tag = ToStr(i+1);
            eleGD_p4[i] = new P4Hists("ele" + tag + "GD");
            eleBD_p4[i] = new P4Hists("ele" + tag + "BD");
            eleGD_dR[i] = Book1D("ele" + tag + "GD_dR", 400, 0.0, 4.0);
            eleBD_dR[i] = Book1D("ele" + tag + "BD_dR", 400, 0.0, 4.0);
            eleGD_dPt[i] = Book1D("ele" + tag + "GD_dPt", 400, 0.0, 4.0);
            eleBD_dPt[i] = Book1D("ele" + tag + "BD_dPt", 400, 0.0, 4.0);
            eleGD_dE[i] = Book1D("ele" + tag + "GD_dE", 400, 0.0, 4.0);
            eleBD_dE[i] = Book1D("ele" + tag + "BD_dE", 400, 0.0, 4.0);
        }
    }
    hltEle27_ptX_etaY_T = Book2D("hltEle27_ptX_etaY_T", 1000, 0.0, 1000.0, 5000, 0.0, 5.0);
    hltEle8_ptX_etaY_TnP = Book2D("hltEle8_ptX_etaY_TnP", 1000, 0.0, 1000.0, 5000, 0.0, 5.0);
    hltEle17_ptX_etaY_TnP = Book2D("hltEle17_ptX_etaY_TnP", 1000, 0.0, 1000.0, 5000, 0.0, 5.0);
    return;
}

// test electron ID --------------------------------------------------------------------------------
bool Electrons::IsGoodID(const EventInfo::Electron &ele) const {
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

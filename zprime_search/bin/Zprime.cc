#include "MyAnalyses/CommonTools/interface/Utils.h"
#include "MyAnalyses/ZprimeSearch/bin/Zprime.h"
#include "TMath.h"
#include <algorithm>
#include <cmath>

// constructor -------------------------------------------------------------------------------------
Zprime::Zprime(std::string _runType) : runType(_runType) {
    isMC = runType.find("mc") != std::string::npos;
    isHF = runType.find("hf") != std::string::npos;
    doUNF = runType.find("unf") != std::string::npos;
    doSYS = runType.find("sys") != std::string::npos;
    if (isMC) genInfo = new GenInfo(!doSYS);
    bool doQCD = runType.find("qcd") != std::string::npos;
    if (isHF) hfElectrons = new HfElectrons(isMC, doQCD, !doSYS);
    else electrons = new Electrons(isMC, doQCD, !doSYS);
    BookHistograms();
    if (doUNF) BookResponses();
}

// destructor --------------------------------------------------------------------------------------
Zprime::~Zprime() {
    if (isMC) delete genInfo;
    if (isHF) delete hfElectrons;
    else delete electrons;
    if (!doSYS) {
        delete raw_pv_size_true;
        delete raw_pv_size;
        if (isMC) {
            delete llGD_p4;
            delete llBD_p4;
            delete llGD_dR;
            delete llBD_dR;
            delete llGD_dPt;
            delete llBD_dPt;
            delete llGD_dE;
            delete llBD_dE;
        }
        delete ll_p4;
    }
    delete ll_massX_cosThetaY_0y1;
    delete ll_massX_cosThetaY_1y1p25;
    delete ll_massX_cosThetaY_1p25y1p5;
    delete ll_massX_cosThetaY_1p5y2p4;
    delete ll_massX_cosThetaY_2p4y5;
    if (!doSYS) {
        delete llMassX_eleEtaY;
        delete llMassX_posEtaY;
        delete llRapidityX_pvSizeY;
        delete llMassX_pvSizeY;
        delete llRapidityX_pfJetSizeY;
        delete llRapidityX_pfMetY;
        delete llRapidityX_pfMetOverSumEtY;
        delete llRapidityX_eleDphiY;
        delete llMassX_eleDphiY;
        delete llRapidityX_eleDetaY;
        delete llRapidityX_eleEtaProdY;
        delete llRapidityX_llMassY;
    }
    if (doUNF) {
        delete response_0y1;
        delete response_1y1p25;
        delete response_1p25y1p5;
        delete response_1p5y2p4;
        delete response_2p4y5;
    }
}

// process event -----------------------------------------------------------------------------------
void Zprime::ProcessEvent(const EventInfo &event, const Corrections &corr, int ord) {
    double weight(1.0), rawWeight(1.0);
    if (isMC) {
        weight = event.auxInfo.weightMC;
        rawWeight = event.auxInfo.weightMC * event.auxInfo.weightPU;
        if (runType.find("sysFsr") != std::string::npos) weight *= event.auxInfo.weightFSR;
        if (runType.find("sysCt10pdf") != std::string::npos)
            weight *= event.auxInfo.weightPDF_ct10.at(ord-1);
        if (runType.find("sysNnpdf") != std::string::npos)
            weight *= event.auxInfo.weightPDF_nnpdf.at(ord-1);
    }
    if (isHF) {
        hfElectrons->Reset();
        if (!doSYS) hfElectrons->FillEventCount();
        double etaMin(0.0), etaMax(2.4);
        if (runType.find("ebhf") != std::string::npos) {
            etaMin = 0.0; etaMax = 1.442;
        }
        else if (runType.find("eehf") != std::string::npos) {
            etaMin = 1.556; etaMax = 2.4;
        }
        if (isMC) {
            genInfo->SelectParticles(event.particles);
            if (!doSYS && genInfo->IsGood()) genInfo->FillRawHistograms(weight);
        }
        if (!doSYS) hfElectrons->FillTriggerHistograms(event, weight);
        if (IsGoodTrigger(event.ele27wp80)) {
            if (!doSYS) hfElectrons->FillEventCount("goodTrigger");
            if (IsGoodPV(event.primVertices)) {
                if (!doSYS) hfElectrons->FillEventCount("goodPV");
                hfElectrons->SelectElectrons(event.electrons, etaMin, etaMax, event.hfElectrons);
                if (!doSYS && isMC && genInfo->IsGood())
                    hfElectrons->FillRawMcHistograms(event, genInfo->GetEle(),
                                                     genInfo->GetPos(), weight);
                if (!doSYS) {
                    hfElectrons->FillRawHistograms(event, rawWeight);
                    FillRawHistograms(event, rawWeight);
                }
                if (hfElectrons->IsGood()) {
                    GetLL(hfElectrons->GetEle().p4, hfElectrons->GetHfEle().p4);
                    if (isMC) {
                        if (!doSYS && genInfo->IsGood()) {
                            genInfo->FillHistograms(weight);
                            hfElectrons->FillMcHistograms(genInfo->GetEle(),
                                                          genInfo->GetPos(), weight);
                            FillMcHistograms(genInfo->GetZ(), weight);
                        }
                        if (runType.find("sysPu") != std::string::npos)
                             weight *= event.auxInfo.weightPU + 0.05*ord;
                        else weight *= event.auxInfo.weightPU;
                        if (runType.find("sysTrg27") != std::string::npos)
                             weight *= corr.GetTrgEle27SF(hfElectrons->GetEle().p4, ord);
                        else weight *= corr.GetTrgEle27SF(hfElectrons->GetEle().p4);
                        if (runType.find("sysEleId") != std::string::npos)
                             weight *= corr.GetTightEleIdSF(hfElectrons->GetEle().p4, ord);
                        else weight *= corr.GetTightEleIdSF(hfElectrons->GetEle().p4);
                        if (runType.find("sysHfId") != std::string::npos)
                             weight *= corr.GetHfEleIdSF(hfElectrons->GetHfEle().p4, ord);
                        else weight *= corr.GetHfEleIdSF(hfElectrons->GetHfEle().p4);
                        if (runType.find("sysHfEta") != std::string::npos)
                             weight *= corr.GetHfEleEtaSF(hfElectrons->GetHfEle().p4, ord);
                        else weight *= corr.GetHfEleEtaSF(hfElectrons->GetHfEle().p4);
                    }
                    if (!doSYS) hfElectrons->FillHistograms(weight);
                    FillHistograms(event, weight);
                }
            }
        }
        if (doUNF) FillResponses(weight);
    }
    else {
        electrons->Reset();
        if (!doSYS) electrons->FillEventCount();
        double eta1Min(0.0), eta1Max(2.4), eta2Min(0.0), eta2Max(2.4);
        if (runType.find("ebeb") != std::string::npos) {
            eta1Min = 0.0; eta1Max = 1.442;
            eta2Min = 0.0; eta2Max = 1.442;
        }
        else if (runType.find("ebee") != std::string::npos) {
            eta1Min = 0.0; eta1Max = 1.442;
            eta2Min = 1.556; eta2Max = 2.4;
        }
        else if (runType.find("eeeb") != std::string::npos) {
            eta1Min = 1.556; eta1Max = 2.4;
            eta2Min = 0.0; eta2Max = 1.442;
        }
        else if (runType.find("eeee") != std::string::npos) {
            eta1Min = 1.556; eta1Max = 2.4;
            eta2Min = 1.556; eta2Max = 2.4;
        }
        if (isMC) {
            genInfo->SelectParticles(event.particles);
            if (!doSYS && genInfo->IsGood()) genInfo->FillRawHistograms(weight);
        }
        if (!doSYS) electrons->FillTriggerHistograms(event, weight);
        if (IsGoodTrigger(event.ele17ele8)) {
            if (!doSYS) electrons->FillEventCount("goodTrigger");
            if (IsGoodPV(event.primVertices)) {
                if (!doSYS) electrons->FillEventCount("goodPV");
                electrons->SelectElectrons(event.electrons, eta1Min, eta1Max, eta2Min, eta2Max);
                if (!doSYS) {
                    electrons->FillRawHistograms(event, rawWeight);
                    FillRawHistograms(event, rawWeight);
                }
                if (electrons->IsGood()) {
                    GetLL(electrons->GetEle1().p4, electrons->GetEle2().p4);
                    if (isMC) {
                        if (!doSYS && genInfo->IsGood()) {
                            genInfo->FillHistograms(weight);
                            electrons->FillMcHistograms(genInfo->GetEle(),genInfo->GetPos(),weight);
                            FillMcHistograms(genInfo->GetZ(), weight);
                        }
                        if (runType.find("sysPu") != std::string::npos)
                             weight *= event.auxInfo.weightPU + 0.05*ord;
                        else weight *= event.auxInfo.weightPU;
                        if (runType.find("sysTrg17") != std::string::npos)
                             weight *= corr.GetTrgEle17SF(electrons->GetEle1().p4, ord);
                        else weight *= corr.GetTrgEle17SF(electrons->GetEle1().p4);
                        if (runType.find("sysTrg8") != std::string::npos)
                             weight *= corr.GetTrgEle8SF(electrons->GetEle2().p4, ord);
                        else weight *= corr.GetTrgEle8SF(electrons->GetEle2().p4);
                        if (runType.find("sysEleId") != std::string::npos) {
                            weight *= corr.GetTightEleIdSF(electrons->GetEle1().p4, ord);
                            weight *= corr.GetTightEleIdSF(electrons->GetEle2().p4, ord);
                        }
                        else {
                            weight *= corr.GetTightEleIdSF(electrons->GetEle1().p4);
                            weight *= corr.GetTightEleIdSF(electrons->GetEle2().p4);
                        }
                    }
                    if (!doSYS) electrons->FillHistograms(weight);
                    FillHistograms(event, weight);
                }
            }
        }
        if (doUNF) FillResponses(weight);
    }
    return;
}

// write histograms and responses ------------------------------------------------------------------
void Zprime::Write() {
    if (!doSYS) {
        if (isMC) genInfo->Write();
        if (isHF) hfElectrons->Write();
        else electrons->Write();
        raw_pv_size_true->Write();
        raw_pv_size->Write();
        if (isMC) {
            llGD_p4->Write();
            llBD_p4->Write();
            llGD_dR->Write();
            llBD_dR->Write();
            llGD_dPt->Write();
            llBD_dPt->Write();
            llGD_dE->Write();
            llBD_dE->Write();
        }
        ll_p4->Write();
    }
    ll_massX_cosThetaY_0y1->Write();
    ll_massX_cosThetaY_1y1p25->Write();
    ll_massX_cosThetaY_1p25y1p5->Write();
    ll_massX_cosThetaY_1p5y2p4->Write();
    ll_massX_cosThetaY_2p4y5->Write();
    if (!doSYS) {
        llMassX_eleEtaY->Write();
        llMassX_posEtaY->Write();
        llRapidityX_pvSizeY->Write();
        llMassX_pvSizeY->Write();
        llRapidityX_pfJetSizeY->Write();
        llRapidityX_pfMetY->Write();
        llRapidityX_pfMetOverSumEtY->Write();
        llRapidityX_eleDphiY->Write();
        llMassX_eleDphiY->Write();
        llRapidityX_eleDetaY->Write();
        llRapidityX_eleEtaProdY->Write();
        llRapidityX_llMassY->Write();
    }
    if (doUNF) {
        response_0y1->Write();
        response_1y1p25->Write();
        response_1p25y1p5->Write();
        response_1p5y2p4->Write();
        response_2p4y5->Write();
    }
    return;
}

// compute cos(theta) ------------------------------------------------------------------------------
double Zprime::GetCosTheta(const TLorentzVector &ele, const TLorentzVector &pos) {
    double eleP = (ele.E() + ele.Pz()) / sqrt(2);
    double eleM = (ele.E() - ele.Pz()) / sqrt(2);
    double posP = (pos.E() + pos.Pz()) / sqrt(2);
    double posM = (pos.E() - pos.Pz()) / sqrt(2);
    TLorentzVector p4 = ele + pos;
    double cosTheta = 2 * (eleP * posM - eleM * posP) /
    sqrt(pow(p4.M(),2) * (pow(p4.M(),2) + pow(p4.Pt(),2)));
    if (p4.Pz() < 0.0) cosTheta = -cosTheta;
    return cosTheta;
}

// AFB mass binning --------------------------------------------------------------------------------
double Zprime::massBins_ele[17] = {0.0, 40.0, 50.0, 60.0, 76.0, 86.0, 96.0, 106.0, 120.0,
                                   133.0, 150.0, 171.0, 200.0, 320.0, 500.0, 2000.0, 8000.0};
double Zprime::massBins_hf[10] = {0.0, 40.0, 76.0, 86.0, 96.0, 106.0, 120.0, 150.0, 320.0, 8000.0};

// test trigger conditions -------------------------------------------------------------------------
bool Zprime::IsGoodTrigger(const std::vector<bool> &triggers) const {
    for (unsigned int i = 0; i != triggers.size(); ++i) if (triggers.at(i)) return true;
    return false;
}

// test for a good primary vertex ------------------------------------------------------------------
bool Zprime::IsGoodPV(const std::vector<EventInfo::PrimVertex> &primVertices) const {
    for (unsigned int i = 0; i != primVertices.size(); ++i) {
        if ((primVertices.at(i).ndof > 4.0) &&
            (primVertices.at(i).rho < 2.0) &&
            (fabs(primVertices.at(i).v3.z()) < 24.0)) return true;
    }
    return false;
}

// reconstruct the ll object -----------------------------------------------------------------------
void Zprime::GetLL(const TLorentzVector &p1, const TLorentzVector &p2) {
    ll = p1 + p2;
    return;
}

// fill raw histograms -----------------------------------------------------------------------------
void Zprime::FillRawHistograms(const EventInfo &event, double weight) {
    raw_pv_size_true->Fill(event.primVertices.size(), event.auxInfo.weightMC);
    raw_pv_size->Fill(event.primVertices.size(), weight);
    return;
}

// fill MC only histograms -------------------------------------------------------------------------
void Zprime::FillMcHistograms(const EventInfo::Particle &par, double weight) {
    double dR = ll.DeltaR(par.p4);
    if (dR < 0.3) {
        llGD_p4->Fill(ll, weight);
        llGD_dR->Fill(dR, weight);
        llGD_dPt->Fill(ll.Pt() / par.p4.Pt(), weight);
        llGD_dE->Fill(ll.Energy() / par.p4.Energy(), weight);
    }
    else {
        llBD_p4->Fill(ll, weight);
        llBD_dR->Fill(dR, weight);
        llBD_dPt->Fill(ll.Pt() / par.p4.Pt(), weight);
        llBD_dE->Fill(ll.Energy() / par.p4.Energy(), weight);
    }
    return;
}

// fill histograms ---------------------------------------------------------------------------------
void Zprime::FillHistograms(const EventInfo &event, double weight) {
    if (!doSYS) ll_p4->Fill(ll, weight);
    const TLorentzVector *ele(0), *pos(0);
    if (isHF) {
        if (hfElectrons->GetEle().charge < 0) {
            ele = &hfElectrons->GetEle().p4;
            pos = &hfElectrons->GetHfEle().p4;
        }
        else {
            ele = &hfElectrons->GetHfEle().p4;
            pos = &hfElectrons->GetEle().p4;
        }
    }
    else {
        if (electrons->GetEle1().charge < 0) {
            ele = &electrons->GetEle1().p4;
            pos = &electrons->GetEle2().p4;
        }
        else {
            ele = &electrons->GetEle2().p4;
            pos = &electrons->GetEle1().p4;
        }
    }
    double cosTheta = GetCosTheta(*ele, *pos);
    if (fabs(ll.Rapidity()) < 1.0) {
        ll_massX_cosThetaY_0y1->Fill(ll.M(), cosTheta, weight);
    }
    else if (fabs(ll.Rapidity()) < 1.25) {
        ll_massX_cosThetaY_1y1p25->Fill(ll.M(), cosTheta, weight);
    }
    else if (fabs(ll.Rapidity()) < 1.5) {
        ll_massX_cosThetaY_1p25y1p5->Fill(ll.M(), cosTheta, weight);
    }
    else if (fabs(ll.Rapidity()) < 2.4) {
        ll_massX_cosThetaY_1p5y2p4->Fill(ll.M(), cosTheta, weight);
    }
    else if (fabs(ll.Rapidity()) < 5.0) {
        ll_massX_cosThetaY_2p4y5->Fill(ll.M(), cosTheta, weight);
    }
    if (!doSYS) {
        llMassX_eleEtaY->Fill(ll.M(), ele->Eta(), weight);
        llMassX_posEtaY->Fill(ll.M(), pos->Eta(), weight);
        llRapidityX_pvSizeY->Fill(ll.Rapidity(), event.primVertices.size(), weight);
        llMassX_pvSizeY->Fill(ll.M(), event.primVertices.size(), weight);
        unsigned int pfJetSize(0);
        for (unsigned int i = 0; i != event.pfJets.size(); ++i) {
            if (event.pfJets.at(i).p4.Pt() < 20.0) continue;
            if (ele->DeltaR(event.pfJets.at(i).p4) < 0.3) continue;
            if (pos->DeltaR(event.pfJets.at(i).p4) < 0.3) continue;
            ++pfJetSize;
        }
        llRapidityX_pfJetSizeY->Fill(ll.Rapidity(), pfJetSize, weight);
        llRapidityX_pfMetY->Fill(ll.Rapidity(), event.pfMet.et, weight);
        llRapidityX_pfMetOverSumEtY->Fill(ll.Rapidity(), event.pfMet.et/event.pfMet.sumEt, weight);
        llRapidityX_eleDphiY->Fill(ll.Rapidity(), fabs(ele->DeltaPhi(*pos)), weight);
        llMassX_eleDphiY->Fill(ll.M(), fabs(ele->DeltaPhi(*pos)), weight);
        llRapidityX_eleDetaY->Fill(ll.Rapidity(), fabs(ele->Eta()-pos->Eta()), weight);
        llRapidityX_eleEtaProdY->Fill(ll.Rapidity(), ele->Eta()*pos->Eta(), weight);
        llRapidityX_llMassY->Fill(ll.Rapidity(), ll.M(), weight);
    }
    return;
}

// fill responses ----------------------------------------------------------------------------------
void Zprime::FillResponses(double weight) {
    double recoY(9.9), recoM(0.0), recoCT(9.9);
    if (isHF && hfElectrons->IsGood()) {
        recoY = ll.Rapidity();
        recoM = ll.M();
        if (hfElectrons->GetEle().charge < 0) {
            recoCT = GetCosTheta(hfElectrons->GetEle().p4, hfElectrons->GetHfEle().p4);
        }
        else {
            recoCT = GetCosTheta(hfElectrons->GetHfEle().p4, hfElectrons->GetEle().p4);
        }
    }
    else if (!isHF && electrons->IsGood()) {
        recoY = ll.Rapidity();
        recoM = ll.M();
        if (electrons->GetEle1().charge < 0) {
            recoCT = GetCosTheta(electrons->GetEle1().p4, electrons->GetEle2().p4);
        }
        else {
            recoCT = GetCosTheta(electrons->GetEle2().p4, electrons->GetEle1().p4);
        }
    }
    double genY(9.9), genM(0.0), genCT(9.9);
    if (genInfo->IsGood()) {
        TLorentzVector p4 = genInfo->GetEle().p4 + genInfo->GetPos().p4;
        genY = p4.Rapidity();
        genM = p4.M();
        genCT = GetCosTheta(genInfo->GetEle().p4, genInfo->GetPos().p4);
    }
    if (fabs(genY) > 0.0 && fabs(genY) < 1.0 && fabs(recoY) > 0.0 && fabs(recoY) < 1.0) {
        response_0y1->Fill(recoM, recoCT, genM, genCT, weight);
    }
    else {
        if (fabs(genY) > 0.0 && fabs(genY) < 1.0) response_0y1->Miss(genM, genCT, weight);
        if (fabs(recoY) > 0.0 && fabs(recoY) < 1.0) response_0y1->Fake(recoM, recoCT, weight);
    }
    if (fabs(genY) > 1.0 && fabs(genY) < 1.25 && fabs(recoY) > 1.0 && fabs(recoY) < 1.25) {
        response_1y1p25->Fill(recoM, recoCT, genM, genCT, weight);
    }
    else {
        if (fabs(genY) > 1.0 && fabs(genY) < 1.25) response_1y1p25->Miss(genM, genCT, weight);
        if (fabs(recoY) > 1.0 && fabs(recoY) < 1.25) response_1y1p25->Fake(recoM, recoCT, weight);
    }
    if (fabs(genY) > 1.25 && fabs(genY) < 1.5 && fabs(recoY) > 1.25 && fabs(recoY) < 1.5) {
        response_1p25y1p5->Fill(recoM, recoCT, genM, genCT, weight);
    }
    else {
        if (fabs(genY) > 1.25 && fabs(genY) < 1.5) response_1p25y1p5->Miss(genM, genCT, weight);
        if (fabs(recoY) > 1.25 && fabs(recoY) < 1.5) response_1p25y1p5->Fake(recoM, recoCT, weight);
    }
    if (fabs(genY) > 1.5 && fabs(genY) < 2.4 && fabs(recoY) > 1.5 && fabs(recoY) < 2.4) {
        response_1p5y2p4->Fill(recoM, recoCT, genM, genCT, weight);
    }
    else {
        if (fabs(genY) > 1.5 && fabs(genY) < 2.4) response_1p5y2p4->Miss(genM, genCT, weight);
        if (fabs(recoY) > 1.5 && fabs(recoY) < 2.4) response_1p5y2p4->Fake(recoM, recoCT, weight);
    }
    if (fabs(genY) > 2.4 && fabs(genY) < 5.0 && fabs(recoY) > 2.4 && fabs(recoY) < 5.0) {
        response_2p4y5->Fill(recoM, recoCT, genM, genCT, weight);
    }
    else {
        if (fabs(genY) > 2.4 && fabs(genY) < 5.0) response_2p4y5->Miss(genM, genCT, weight);
        if (fabs(recoY) > 2.4 && fabs(recoY) < 5.0) response_2p4y5->Fake(recoM, recoCT, weight);
    }
    return;
}

// book histograms ---------------------------------------------------------------------------------
void Zprime::BookHistograms() {
    if (!doSYS) {
        raw_pv_size_true = Book1D("raw_pv_size_true", 80, -0.5, 79.5);
        raw_pv_size = Book1D("raw_pv_size", 80, -0.5, 79.5);
        if (isMC) {
            llGD_p4 = new P4Hists("llGD");
            llBD_p4 = new P4Hists("llBD");
            llGD_dR = Book1D("llGD_dR", 400, 0.0, 4.0);
            llBD_dR = Book1D("llBD_dR", 400, 0.0, 4.0);
            llGD_dPt = Book1D("llGD_dPt", 400, 0.0, 4.0);
            llBD_dPt = Book1D("llBD_dPt", 400, 0.0, 4.0);
            llGD_dE = Book1D("llGD_dE", 400, 0.0, 4.0);
            llBD_dE = Book1D("llBD_dE", 400, 0.0, 4.0);
        }
        ll_p4 = new P4Hists("ll");
    }
    ll_massX_cosThetaY_0y1 = new TH2D("ll_massX_cosThetaY_0y1", "",
                                      16, massBins_ele, 200, -1.0, 1.0);
    ll_massX_cosThetaY_0y1->Sumw2();
    ll_massX_cosThetaY_1y1p25 = new TH2D("ll_massX_cosThetaY_1y1p25", "",
                                         16, massBins_ele, 200, -1.0, 1.0);
    ll_massX_cosThetaY_1y1p25->Sumw2();
    ll_massX_cosThetaY_1p25y1p5 = new TH2D("ll_massX_cosThetaY_1p25y1p5", "",
                                           16, massBins_ele, 200, -1.0, 1.0);
    ll_massX_cosThetaY_1p25y1p5->Sumw2();
    ll_massX_cosThetaY_1p5y2p4 = new TH2D("ll_massX_cosThetaY_1p5y2p4", "",
                                          16, massBins_ele, 200, -1.0, 1.0);
    ll_massX_cosThetaY_1p5y2p4->Sumw2();
    ll_massX_cosThetaY_2p4y5 = new TH2D("ll_massX_cosThetaY_2p4y5", "",
                                        9, massBins_hf, 200, -1.0, 1.0);
    ll_massX_cosThetaY_2p4y5->Sumw2();
    if (!doSYS) {
        llMassX_eleEtaY = Book2D("llMassX_eleEtaY", 4000, 0.0, 4000.0, 200, -5.0, 5.0);
        llMassX_posEtaY = Book2D("llMassX_posEtaY", 4000, 0.0, 4000.0, 200, -5.0, 5.0);
        llRapidityX_pvSizeY = Book2D("llRapidityX_pvSizeY", 200, -5.0, 5.0, 80, -0.5, 79.5);
        llMassX_pvSizeY = Book2D("llMassX_pvSizeY", 4000, 0.0, 4000.0, 80, -0.5, 79.5);
        llRapidityX_pfJetSizeY = Book2D("llRapidityX_pfJetSizeY", 200, -5.0, 5.0, 10, -0.5, 9.5);
        llRapidityX_pfMetY = Book2D("llRapidityX_pfMetY", 200, -5.0, 5.0, 2000, 0.0, 2000.0);
        llRapidityX_pfMetOverSumEtY = Book2D("llRapidityX_pfMetOverSumEtY",
                                             200, -5.0, 5.0, 400, 0.0, 1.0);
        llRapidityX_eleDphiY = Book2D("llRapidityX_eleDphiY", 200, -5.0, 5.0, 200, 0.0, TMath::Pi());
        llMassX_eleDphiY = Book2D("llMassX_eleDphiY", 4000, 0.0, 4000.0, 200, 0.0, TMath::Pi());
        llRapidityX_eleDetaY = Book2D("llRapidityX_eleDetaY", 200, -5.0, 5.0, 200, 0.0, 10.0);
        llRapidityX_eleEtaProdY = Book2D("llRapidityX_eleEtaProdY",
                                         200, -5.0, 5.0, 480, -12.0, 12.0);
        llRapidityX_llMassY = Book2D("llRapidityX_llMassY", 200, -5.0, 5.0, 4000, 0.0, 4000.0);
    }
    return;
}

// book responses ----------------------------------------------------------------------------------
void Zprime::BookResponses() {
    TH2D hist_ele("hist_ele", "", 16, massBins_ele, 2, -1.0, 1.0);
    response_0y1 = new RooUnfoldResponse(&hist_ele, &hist_ele, "response_0y1");
    response_1y1p25 = new RooUnfoldResponse(&hist_ele, &hist_ele, "response_1y1p25");
    response_1p25y1p5 = new RooUnfoldResponse(&hist_ele, &hist_ele, "response_1p25y1p5");
    response_1p5y2p4 = new RooUnfoldResponse(&hist_ele, &hist_ele, "response_1p5y2p4");
    TH2D hist_hf("hist_hf", "", 9, massBins_hf, 2, -1.0, 1.0);
    response_2p4y5 = new RooUnfoldResponse(&hist_hf, &hist_hf, "response_2p4y5");
    return;
}

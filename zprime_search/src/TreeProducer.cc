#include "MyAnalyses/CommonTools/interface/Utils.h"
#include "MyAnalyses/ZprimeSearch/interface/TreeProducer.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaReco/interface/HFEMClusterShape.h"
#include "DataFormats/EgammaReco/interface/HFEMClusterShapeAssociation.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "EgammaAnalysis/ElectronTools/interface/ElectronEffectiveArea.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>

// constructor -------------------------------------------------------------------------------------
TreeProducer::TreeProducer(const edm::ParameterSet &config) {
    isMC = config.getParameter<bool>("isMC");
    outFileName = config.getParameter<std::string>("outFileName");
    if (isMC) {
        fsrWeightTag = config.getParameter<edm::InputTag>("fsrWeightTag");
        pdfWeightTag_ct10 = config.getParameter<edm::InputTag>("pdfWeightTag_ct10");
        pdfWeightTag_mstw = config.getParameter<edm::InputTag>("pdfWeightTag_mstw");
        pdfWeightTag_nnpdf = config.getParameter<edm::InputTag>("pdfWeightTag_nnpdf");
        genParticleTag = config.getParameter<edm::InputTag>("genParticleTag");
    }
    vertexTag = config.getParameter<edm::InputTag>("vertexTag");
    triggerResultsTag = config.getParameter<edm::InputTag>("triggerResultsTag");
    triggerSummaryTag = config.getParameter<edm::InputTag>("triggerSummaryTag");
    electronTag = config.getParameter<edm::InputTag>("electronTag");
    isoRhoTag = config.getParameter<edm::InputTag>("isoRhoTag");
    conversionTag = config.getParameter<edm::InputTag>("conversionTag");
    beamSpotTag = config.getParameter<edm::InputTag>("beamSpotTag");
    hfElectronTag = config.getParameter<edm::InputTag>("hfElectronTag");
    hfClusterMapTag = config.getParameter<edm::InputTag>("hfClusterMapTag");
    pfJetTag = config.getParameter<edm::InputTag>("pfJetTag");
    pfMetTag = config.getParameter<edm::InputTag>("pfMetTag");
	outFile = new TFile(outFileName.c_str(), "RECREATE");
	tree = new TTree("Events", "");
	tree->Branch("ev_run", &ev_run, "ev_run/i");
    tree->Branch("ev_lumi", &ev_lumi, "ev_lumi/i");
    tree->Branch("ev_nr", &ev_nr, "ev_nr/i");
    tree->Branch("ev_bunch", &ev_bunch, "ev_bunch/i");
    if (isMC) {
        tree->Branch("ev_weightPU", &ev_weightPU, "ev_weightPU/F");
        tree->Branch("ev_weightFSR", &ev_weightFSR, "ev_weightFSR/F");
        tree->Branch("ev_weightPDF_ct10", ev_weightPDF_ct10, "ev_weightPDF_ct10[53]/F");
        tree->Branch("ev_weightPDF_mstw", ev_weightPDF_mstw, "ev_weightPDF_mstw[41]/F");
        tree->Branch("ev_weightPDF_nnpdf", ev_weightPDF_nnpdf, "ev_weightPDF_nnpdf[101]/F");
        tree->Branch("gen_size", &gen_size, "gen_size/i");
        tree->Branch("gen_numberOfDaughters", gen_numberOfDaughters,
                     "gen_numberOfDaughters[gen_size]/i");
        tree->Branch("gen_charge", gen_charge, "gen_charge[gen_size]/I");
        tree->Branch("gen_status", gen_status, "gen_status[gen_size]/I");
        tree->Branch("gen_pdgId", gen_pdgId, "gen_pdgId[gen_size]/I");
        tree->Branch("gen_pdgIdMother", gen_pdgIdMother, "gen_pdgIdMother[gen_size]/I");
        tree->Branch("gen_pdgIdGrandMother", gen_pdgIdGrandMother,
                     "gen_pdgIdGrandMother[gen_size]/I");
        tree->Branch("gen_pt", gen_pt, "gen_pt[gen_size]/F");
        tree->Branch("gen_phi", gen_phi, "gen_phi[gen_size]/F");
        tree->Branch("gen_eta", gen_eta, "gen_eta[gen_size]/F");
        tree->Branch("gen_mass", gen_mass, "gen_mass[gen_size]/F");
        tree->Branch("gen_vx", gen_vx, "gen_vx[gen_size]/F");
        tree->Branch("gen_vy", gen_vy, "gen_vy[gen_size]/F");
        tree->Branch("gen_vz", gen_vz, "gen_vz[gen_size]/F");
    }
    tree->Branch("pv_size", &pv_size, "pv_size/i");
    tree->Branch("pv_nTracks", pv_nTracks, "pv_nTracks[pv_size]/i");
    tree->Branch("pv_chi2", pv_chi2, "pv_chi2[pv_size]/F");
    tree->Branch("pv_normalizedChi2", pv_normalizedChi2, "pv_normalizedChi2[pv_size]/F");
    tree->Branch("pv_ndof", pv_ndof, "pv_ndof[pv_size]/F");
    tree->Branch("pv_x", pv_x, "pv_x[pv_size]/F");
    tree->Branch("pv_y", pv_y, "pv_y[pv_size]/F");
    tree->Branch("pv_z", pv_z, "pv_z[pv_size]/F");
    tree->Branch("pv_rho", pv_rho, "pv_rho[pv_size]/F");
    tree->Branch("hlt_ele8", hlt_ele8, "hlt_ele8[7]/O");
    tree->Branch("hlt_ele17", hlt_ele17, "hlt_ele17[7]/O");
    tree->Branch("hlt_ele17ele8", hlt_ele17ele8, "hlt_ele17ele8[9]/O");
    tree->Branch("hlt_ele23hft30", hlt_ele23hft30, "hlt_ele23hft30[9]/O");
    tree->Branch("hlt_ele27hft15", hlt_ele27hft15, "hlt_ele27hft15[9]/O");
    tree->Branch("hlt_ele32sc17m50", hlt_ele32sc17m50, "hlt_ele32sc17m50[7]/O");
    tree->Branch("hlt_ele27", hlt_ele27, "hlt_ele27[8]/O");
    tree->Branch("hlt_ele27wp80", hlt_ele27wp80, "hlt_ele27wp80[9]/O");
    tree->Branch("ele_size", &ele_size, "ele_size/i");
    tree->Branch("ele_ele8", ele_ele8, "ele_ele8[ele_size]/O");
    tree->Branch("ele_ele17", ele_ele17, "ele_ele17[ele_size]/O");
    tree->Branch("ele_ele17ele8_f1", ele_ele17ele8_f1, "ele_ele17ele8_f1[ele_size]/O");
    tree->Branch("ele_ele17ele8_f2", ele_ele17ele8_f2, "ele_ele17ele8_f2[ele_size]/O");
    tree->Branch("ele_ele23hft30", ele_ele23hft30, "ele_ele23hft30[ele_size]/O");
    tree->Branch("ele_ele27hft15", ele_ele27hft15, "ele_ele27hft15[ele_size]/O");
    tree->Branch("ele_ele32sc17m50_f1", ele_ele32sc17m50_f1, "ele_ele32sc17m50_f1[ele_size]/O");
    tree->Branch("ele_ele32sc17m50_f2", ele_ele32sc17m50_f2, "ele_ele32sc17m50_f2[ele_size]/O");
    tree->Branch("ele_ele32sc17m50_f3", ele_ele32sc17m50_f3, "ele_ele32sc17m50_f3[ele_size]/O");
    tree->Branch("ele_ele27", ele_ele27, "ele_ele27[ele_size]/O");
    tree->Branch("ele_ele27wp80", ele_ele27wp80, "ele_ele27wp80[ele_size]/O");
    tree->Branch("ele_charge", ele_charge, "ele_charge[ele_size]/I");
    tree->Branch("ele_pt", ele_pt, "ele_pt[ele_size]/F");
    tree->Branch("ele_phi", ele_phi, "ele_phi[ele_size]/F");
    tree->Branch("ele_eta", ele_eta, "ele_eta[ele_size]/F");
    tree->Branch("ele_energy", ele_energy, "ele_energy[ele_size]/F");
    tree->Branch("ele_calibError", ele_calibError, "ele_calibError[ele_size]/F");
    tree->Branch("ele_etaSC", ele_etaSC, "ele_etaSC[ele_size]/F");
    tree->Branch("ele_dEtaIn", ele_dEtaIn, "ele_dEtaIn[ele_size]/F");
    tree->Branch("ele_dPhiIn", ele_dPhiIn, "ele_dPhiIn[ele_size]/F");
    tree->Branch("ele_sigmaIEtaIEta", ele_sigmaIEtaIEta, "ele_sigmaIEtaIEta[ele_size]/F");
    tree->Branch("ele_HoE", ele_HoE, "ele_HoE[ele_size]/F");
    tree->Branch("ele_d0", ele_d0, "ele_d0[ele_size]/F");
    tree->Branch("ele_dZ", ele_dZ, "ele_dZ[ele_size]/F");
    tree->Branch("ele_IoEmIoP", ele_IoEmIoP, "ele_IoEmIoP[ele_size]/F");
    tree->Branch("ele_e1x5", ele_e1x5, "ele_e1x5[ele_size]/F");
    tree->Branch("ele_e2x5", ele_e2x5, "ele_e2x5[ele_size]/F");
    tree->Branch("ele_e5x5", ele_e5x5, "ele_e5x5[ele_size]/F");
    tree->Branch("ele_isoCH", ele_isoCH, "ele_isoCH[ele_size]/F");
    tree->Branch("ele_isoNH", ele_isoNH, "ele_isoNH[ele_size]/F");
    tree->Branch("ele_isoEM", ele_isoEM, "ele_isoEM[ele_size]/F");
    tree->Branch("ele_isoNeutralNoPU", ele_isoNeutralNoPU, "ele_isoNeutralNoPU[ele_size]/F");
    tree->Branch("ele_iso", ele_iso, "ele_iso[ele_size]/F");
    tree->Branch("ele_isoDetEM", ele_isoDetEM, "ele_isoDetEM[ele_size]/F");
    tree->Branch("ele_isoDetHad1", ele_isoDetHad1, "ele_isoDetHad1[ele_size]/F");
    tree->Branch("ele_isoDetTrk", ele_isoDetTrk, "ele_isoDetTrk[ele_size]/F");
    tree->Branch("ele_hasConversion", ele_hasConversion, "ele_hasConversion[ele_size]/O");
    tree->Branch("ele_isEcalDriven", ele_isEcalDriven, "ele_isEcalDriven[ele_size]/O");
    tree->Branch("ele_missingHits", ele_missingHits, "ele_missingHits[ele_size]/i");
    tree->Branch("ele_goodVetoId", ele_goodVetoId, "ele_goodVetoId[ele_size]/O");
    tree->Branch("ele_goodLooseId", ele_goodLooseId, "ele_goodLooseId[ele_size]/O");
    tree->Branch("ele_goodMediumId", ele_goodMediumId, "ele_goodMediumId[ele_size]/O");
    tree->Branch("ele_goodTightId", ele_goodTightId, "ele_goodTightId[ele_size]/O");
    tree->Branch("ele_goodHeepId", ele_goodHeepId, "ele_goodHeepId[ele_size]/O");
    tree->Branch("ele_mvaTrig", ele_mvaTrig, "ele_mvaTrig[ele_size]/F");
    tree->Branch("ele_mvaNonTrig", ele_mvaNonTrig, "ele_mvaNonTrig[ele_size]/F");
    tree->Branch("ele_vx", ele_vx, "ele_vx[ele_size]/F");
    tree->Branch("ele_vy", ele_vy, "ele_vy[ele_size]/F");
    tree->Branch("ele_vz", ele_vz, "ele_vz[ele_size]/F");
    tree->Branch("hfEle_size", &hfEle_size, "hfEle_size/i");
    tree->Branch("hfEle_ele23hft30", hfEle_ele23hft30, "hfEle_ele23hft30[hfEle_size]/O");
    tree->Branch("hfEle_ele27hft15", hfEle_ele27hft15, "hfEle_ele27hft15[hfEle_size]/O");
    tree->Branch("hfEle_pt", hfEle_pt, "hfEle_pt[hfEle_size]/F");
    tree->Branch("hfEle_phi", hfEle_phi, "hfEle_phi[hfEle_size]/F");
    tree->Branch("hfEle_eta", hfEle_eta, "hfEle_eta[hfEle_size]/F");
    tree->Branch("hfEle_energy", hfEle_energy, "hfEle_energy[hfEle_size]/F");
    tree->Branch("hfEle_e1x1", hfEle_e1x1, "hfEle_e1x1[hfEle_size]/F");
    tree->Branch("hfEle_e3x3", hfEle_e3x3, "hfEle_e3x3[hfEle_size]/F");
    tree->Branch("hfEle_e5x5", hfEle_e5x5, "hfEle_e5x5[hfEle_size]/F");
    tree->Branch("hfEle_eLong1x1", hfEle_eLong1x1, "hfEle_eLong1x1[hfEle_size]/F");
    tree->Branch("hfEle_eLong3x3", hfEle_eLong3x3, "hfEle_eLong3x3[hfEle_size]/F");
    tree->Branch("hfEle_eLong5x5", hfEle_eLong5x5, "hfEle_eLong5x5[hfEle_size]/F");
    tree->Branch("hfEle_eShort1x1", hfEle_eShort1x1, "hfEle_eShort1x1[hfEle_size]/F");
    tree->Branch("hfEle_eShort3x3", hfEle_eShort3x3, "hfEle_eShort3x3[hfEle_size]/F");
    tree->Branch("hfEle_eShort5x5", hfEle_eShort5x5, "hfEle_eShort5x5[hfEle_size]/F");
    tree->Branch("hfEle_e9e25", hfEle_e9e25, "hfEle_e9e25[hfEle_size]/F");
    tree->Branch("hfEle_eCOREe9", hfEle_eCOREe9, "hfEle_eCOREe9[hfEle_size]/F");
    tree->Branch("hfEle_eCore", hfEle_eCore, "hfEle_eCore[hfEle_size]/F");
    tree->Branch("hfEle_eSeL", hfEle_eSeL, "hfEle_eSeL[hfEle_size]/F");
    tree->Branch("hfEle_vx", hfEle_vx, "hfEle_vx[hfEle_size]/F");
    tree->Branch("hfEle_vy", hfEle_vy, "hfEle_vy[hfEle_size]/F");
    tree->Branch("hfEle_vz", hfEle_vz, "hfEle_vz[hfEle_size]/F");
    tree->Branch("pfJet_size", &pfJet_size, "pfJet_size/i");
    tree->Branch("pfJet_pt", pfJet_pt, "pfJet_pt[pfJet_size]/F");
    tree->Branch("pfJet_phi", pfJet_phi, "pfJet_phi[pfJet_size]/F");
    tree->Branch("pfJet_eta", pfJet_eta, "pfJet_eta[pfJet_size]/F");
    tree->Branch("pfJet_energy", pfJet_energy, "pfJet_energy[pfJet_size]/F");
    tree->Branch("pfJet_nConstituents", pfJet_nConstituents, "pfJet_nConstituents[pfJet_size]/i");
    tree->Branch("pfJet_chargedMultiplicity", pfJet_chargedMultiplicity,
                 "pfJet_chargedMultiplicity[pfJet_size]/i");
    tree->Branch("pfJet_chargedEmEnergyFraction", pfJet_chargedEmEnergyFraction,
                 "pfJet_chargedEmEnergyFraction[pfJet_size]/F");
    tree->Branch("pfJet_chargedHadronEnergyFraction", pfJet_chargedHadronEnergyFraction,
                 "pfJet_chargedHadronEnergyFraction[pfJet_size]/F");
    tree->Branch("pfJet_neutralEmEnergyFraction", pfJet_neutralEmEnergyFraction,
                 "pfJet_neutralEmEnergyFraction[pfJet_size]/F");
    tree->Branch("pfJet_neutralHadronEnergyFraction", pfJet_neutralHadronEnergyFraction,
                 "pfJet_neutralHadronEnergyFraction[pfJet_size]/F");
    tree->Branch("pfJet_goodLooseId", pfJet_goodLooseId, "pfJet_goodLooseId[pfJet_size]/O");
    tree->Branch("pfJet_goodMediumId", pfJet_goodMediumId, "pfJet_goodMediumId[pfJet_size]/O");
    tree->Branch("pfJet_goodTightId", pfJet_goodTightId, "pfJet_goodTightId[pfJet_size]/O");
    tree->Branch("pfJet_jetArea", pfJet_jetArea, "pfJet_jetArea[pfJet_size]/F");
    tree->Branch("pfJet_vx", pfJet_vx, "pfJet_vx[pfJet_size]/F");
    tree->Branch("pfJet_vy", pfJet_vy, "pfJet_vy[pfJet_size]/F");
    tree->Branch("pfJet_vz", pfJet_vz, "pfJet_vz[pfJet_size]/F");
    tree->Branch("pfMet_et", &pfMet_et, "pfMet_et/F");
    tree->Branch("pfMet_sumEt", &pfMet_sumEt, "pfMet_sumEt/F");
    tree->Branch("pfMet_phi", &pfMet_phi, "pfMet_phi/F");
}

// destructor --------------------------------------------------------------------------------------
TreeProducer::~TreeProducer() {
	outFile->cd();
	tree->Write();
	tree->Print();
	delete tree;
	delete outFile;
}

// high-level trigger matching ---------------------------------------------------------------------
bool TreeProducer::FindMatch(const trigger::TriggerEvent &triggerEvent, edm::InputTag &filterTag,
                             TLorentzVector &p4) {
    trigger::size_type index = triggerEvent.filterIndex(filterTag);
    if (index == triggerEvent.sizeFilters()) return false;
    const std::vector<trigger::size_type> keys = triggerEvent.filterKeys(index);
    const std::vector<trigger::TriggerObject> triggerObjects = triggerEvent.getObjects();
    for (unsigned int i = 0; i != keys.size(); ++i) {
        trigger::TriggerObject trgObject = triggerObjects.at(keys.at(i));
        TLorentzVector p4_trg(0.0, 0.0, 0.0, 0.0);
        p4_trg.SetPtEtaPhiM(trgObject.pt(), trgObject.eta(), trgObject.phi(), trgObject.mass());
        double dR = p4.DeltaR(p4_trg);
        if (dR < 0.2) return true;
    }
    return false;
}

// analyze event -----------------------------------------------------------------------------------
void TreeProducer::analyze(const edm::Event &event, const edm::EventSetup &setup) {
    
    // event information
    ev_run = event.id().run();
    ev_lumi = event.luminosityBlock();
	ev_nr = event.id().event();
    ev_bunch = event.bunchCrossing();
    if (isMC) {
        edm::Handle<std::vector<PileupSummaryInfo> > puObj;
        event.getByLabel(edm::InputTag("addPileupInfo"), puObj);
        for (unsigned int i = 0; i != puObj->size(); ++i) {
            if (puObj->at(i).getBunchCrossing() != 0) continue;
            int trueN = puObj->at(i).getTrueNumInteractions();
            edm::LumiReWeighting lumiReWeighting("pileup_MC.root", "pileup_DATA.root",
                                                 "pileup", "pileup");
            ev_weightPU = lumiReWeighting.weight(trueN);
            break;
        }
        edm::Handle<double> fsrObj;
        event.getByLabel(fsrWeightTag, fsrObj);
        ev_weightFSR = *fsrObj;
        edm::Handle<std::vector<double> > ct10Obj;
        event.getByLabel(pdfWeightTag_ct10, ct10Obj);
        for (unsigned int i = 0; i != 53; ++i) {
            ev_weightPDF_ct10[i] = ct10Obj->at(i);
        }
        edm::Handle<std::vector<double> > mstwObj;
        event.getByLabel(pdfWeightTag_mstw, mstwObj);
        for (unsigned int i = 0; i != 41; ++i) {
            ev_weightPDF_mstw[i] = mstwObj->at(i);
        }
        edm::Handle<std::vector<double> > nnpdfObj;
        event.getByLabel(pdfWeightTag_nnpdf, nnpdfObj);
        for (unsigned int i = 0; i != 101; ++i) {
            ev_weightPDF_nnpdf[i] = nnpdfObj->at(i);
        }
    }
    
    // generated particles (e, Z, W, Z')
    if (isMC) {
        edm::Handle<std::vector<reco::GenParticle> > genObj;
        event.getByLabel(genParticleTag, genObj);
        unsigned int genCount(0);
        for (unsigned int i = 0; i != genObj->size() && genCount != 99; ++i) {
            if ((genObj->at(i).status() != 1 && genObj->at(i).status() != 3) ||
                (abs(genObj->at(i).pdgId()) != 11 && abs(genObj->at(i).pdgId()) != 23 &&
                 abs(genObj->at(i).pdgId()) != 24 && abs(genObj->at(i).pdgId()) != 32)) continue;
            gen_numberOfDaughters[genCount] = genObj->at(i).numberOfDaughters();
            gen_charge[genCount] = genObj->at(i).charge();
            gen_status[genCount] = genObj->at(i).status();
            gen_pdgId[genCount] = genObj->at(i).pdgId();
            gen_pdgIdMother[genCount] = 0;
            gen_pdgIdGrandMother[genCount] = 0;
            if (genObj->at(i).mother() != 0) {
                gen_pdgIdMother[genCount] = genObj->at(i).mother()->pdgId();
                if (genObj->at(i).mother()->mother() != 0)
                gen_pdgIdGrandMother[genCount] = genObj->at(i).mother()->mother()->pdgId();
            }
            gen_pt[genCount] = genObj->at(i).pt();
            gen_phi[genCount] = genObj->at(i).phi();
            gen_eta[genCount] = genObj->at(i).eta();
            gen_mass[genCount] = genObj->at(i).mass();
            gen_vx[genCount] = genObj->at(i).vx();
            gen_vy[genCount] = genObj->at(i).vy();
            gen_vz[genCount] = genObj->at(i).vz();
            ++genCount;
        }
        gen_size = genCount;
    }
    
    // primary vertex information
    edm::Handle<std::vector<reco::Vertex> > pvObj;
    event.getByLabel(vertexTag, pvObj);
    pv_size = pvObj->size();
    for (unsigned int i = 0; i != pvObj->size() && i != 99; ++i) {
        pv_nTracks[i] = pvObj->at(i).nTracks();
        pv_chi2[i] = pvObj->at(i).chi2();
        pv_normalizedChi2[i] = pvObj->at(i).normalizedChi2();
        pv_ndof[i] = pvObj->at(i).ndof();
        pv_x[i] = pvObj->at(i).x();
        pv_y[i] = pvObj->at(i).y();
        pv_z[i] = pvObj->at(i).z();
        pv_rho[i] = pvObj->at(i).position().rho();
    }
    
    // trigger information
    for (unsigned int i = 0; i != 7; ++i) hlt_ele8[i] = false;
    for (unsigned int i = 0; i != 7; ++i) hlt_ele17[i] = false;
    for (unsigned int i = 0; i != 9; ++i) hlt_ele17ele8[i] = false;
    for (unsigned int i = 0; i != 9; ++i) hlt_ele23hft30[i] = false;
    for (unsigned int i = 0; i != 9; ++i) hlt_ele27hft15[i] = false;
    for (unsigned int i = 0; i != 7; ++i) hlt_ele32sc17m50[i] = false;
    for (unsigned int i = 0; i != 8; ++i) hlt_ele27[i] = false;
    for (unsigned int i = 0; i != 9; ++i) hlt_ele27wp80[i] = false;
    edm::Handle<edm::TriggerResults> trgResObj;
    event.getByLabel(triggerResultsTag, trgResObj);
    edm::TriggerNames hltNames(event.triggerNames(*trgResObj));
    for (unsigned int i = 0; i != trgResObj->size(); ++i) {
        std::string name = "HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v";
        for (unsigned int j = 0; j != 7; ++j) {
            if (hltNames.triggerName(i).compare(name + ToStr(10 + j)) == 0)
                hlt_ele8[j] = trgResObj->accept(i);
        }
        name = "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v";
        for (unsigned int j = 0; j != 7; ++j) {
            if (hltNames.triggerName(i).compare(name + ToStr(1 + j)) == 0)
                hlt_ele17[j] = trgResObj->accept(i);
        }
        name = std::string("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_") +
                           "Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v";
        for (unsigned int j = 0; j != 9; ++j) {
            if (hltNames.triggerName(i).compare(name + ToStr(12 + j)) == 0)
                hlt_ele17ele8[j] = trgResObj->accept(i);
        }
        name = "HLT_Ele23_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_HFT30_v";
        for (unsigned int j = 0; j != 9; ++j) {
            if (hltNames.triggerName(i).compare(name + ToStr(1 + j)) == 0)
                hlt_ele23hft30[j] = trgResObj->accept(i);
        }
        name = "HLT_Ele27_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_HFT15_v";
        for (unsigned int j = 0; j != 9; ++j) {
            if (hltNames.triggerName(i).compare(name + ToStr(1 + j)) == 0)
                hlt_ele27hft15[j] = trgResObj->accept(i);
        }
        name = "HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_v";
        for (unsigned int j = 0; j != 7; ++j) {
            if (hltNames.triggerName(i).compare(name + ToStr(1 + j)) == 0)
                hlt_ele32sc17m50[j] = trgResObj->accept(i);
        }
        name = "HLT_Ele27_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_v";
        for (unsigned int j = 0; j != 8; ++j) {
            if (hltNames.triggerName(i).compare(name + ToStr(5 + j)) == 0)
                hlt_ele27[j] = trgResObj->accept(i);
        }
        name = "HLT_Ele27_WP80_v";
        for (unsigned int j = 0; j != 9; ++j) {
            if (hltNames.triggerName(i).compare(name + ToStr(5 + j)) == 0)
                hlt_ele27wp80[j] = trgResObj->accept(i);
        }
    }
    
    // electron information
    edm::Handle<trigger::TriggerEvent> trgSumObj;
    event.getByLabel(triggerSummaryTag, trgSumObj);
    edm::Handle<std::vector<pat::Electron> > eleObj;
    event.getByLabel(electronTag, eleObj);
    edm::Handle<double> isoRhoObj;
    event.getByLabel(isoRhoTag, isoRhoObj);
    edm::Handle<std::vector<reco::Conversion> > convObj;
    event.getByLabel(conversionTag, convObj);
    edm::Handle<reco::BeamSpot> beamObj;
    event.getByLabel(beamSpotTag, beamObj);
    ele_size = eleObj->size();
    for (unsigned int i = 0; i != eleObj->size() && i != 99; ++i) {
        TLorentzVector p4(0.0, 0.0, 0.0, 0.0);
        p4.SetPtEtaPhiE(eleObj->at(i).pt(), eleObj->at(i).eta(),
                        eleObj->at(i).phi(), eleObj->at(i).energy());
        edm::InputTag filterTag("hltEle8TightIdLooseIsoTrackIsolFilter::HLT");
        ele_ele8[i] = FindMatch(*trgSumObj, filterTag, p4);
        if (!ele_ele8[i]) {
            filterTag = edm::InputTag("hltEle8TightIdLooseIsoTrackIsoFilter::HLT");
            ele_ele8[i] = FindMatch(*trgSumObj, filterTag, p4);
        }
        filterTag = edm::InputTag("hltEle17CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter::HLT");
        ele_ele17[i] = FindMatch(*trgSumObj, filterTag, p4);
        filterTag = edm::InputTag(std::string("hltEle17TightIdLooseIso") +
                                              "Ele8TightIdLooseIsoTrackIsolFilter::HLT");
        ele_ele17ele8_f1[i] = FindMatch(*trgSumObj, filterTag, p4);
        if (!ele_ele17ele8_f1[i]) {
            filterTag = edm::InputTag(std::string("hltEle17TightIdLooseIso") +
                                                  "Ele8TightIdLooseIsoTrackIsoFilter::HLT");
            ele_ele17ele8_f1[i] = FindMatch(*trgSumObj, filterTag, p4);
        }
        filterTag = edm::InputTag(std::string("hltEle17TightIdLooseIso") +
                                              "Ele8TightIdLooseIsoTrackIsolDoubleFilter::HLT");
        ele_ele17ele8_f2[i] = FindMatch(*trgSumObj, filterTag, p4);
        if (!ele_ele17ele8_f2[i]) {
            filterTag = edm::InputTag(std::string("hltEle17TightIdLooseIso") +
                                                  "Ele8TightIdLooseIsoTrackIsoDoubleFilter::HLT");
            ele_ele17ele8_f2[i] = FindMatch(*trgSumObj, filterTag, p4);
        }
        filterTag = edm::InputTag("hltEle23TightIdLooseIsoTrackIsolFilter::HLT");
        ele_ele23hft30[i] = FindMatch(*trgSumObj, filterTag, p4);
        if (!ele_ele23hft30[i]) {
            filterTag = edm::InputTag("hltEle23TightIdLooseIsoTrackIsoFilter::HLT");
            ele_ele23hft30[i] = FindMatch(*trgSumObj, filterTag, p4);
        }
        filterTag = edm::InputTag("hltEle27TightIdLooseIsoTrackIsolFilter::HLT");
        ele_ele27hft15[i] = FindMatch(*trgSumObj, filterTag, p4);
        if (!ele_ele27hft15[i]) {
            filterTag = edm::InputTag("hltEle27TightIdLooseIsoTrackIsoFilter::HLT");
            ele_ele27hft15[i] = FindMatch(*trgSumObj, filterTag, p4);
        }
        filterTag = edm::InputTag("hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsolFilter::HLT");
        ele_ele32sc17m50_f1[i] = FindMatch(*trgSumObj, filterTag, p4);
        if (!ele_ele32sc17m50_f1[i]) {
            filterTag = edm::InputTag("hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsoFilter::HLT");
            ele_ele32sc17m50_f1[i] = FindMatch(*trgSumObj, filterTag, p4);
        }
        filterTag = edm::InputTag("hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17HEDoubleFilter::HLT");
        ele_ele32sc17m50_f2[i] = FindMatch(*trgSumObj, filterTag, p4);
        filterTag = edm::InputTag("hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PMMassFilter::HLT");
        ele_ele32sc17m50_f3[i] = FindMatch(*trgSumObj, filterTag, p4);
        filterTag = edm::InputTag("hltEle27CaloIdLCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter::HLT");
        ele_ele27[i] = FindMatch(*trgSumObj, filterTag, p4);
        filterTag = edm::InputTag("hltEle27WP80TrackIsoFilter::HLT");
        ele_ele27wp80[i] = FindMatch(*trgSumObj, filterTag, p4);
        ele_charge[i] = eleObj->at(i).charge();
        ele_pt[i] = eleObj->at(i).pt();
        ele_phi[i] = eleObj->at(i).phi();
        ele_eta[i] = eleObj->at(i).eta();
        ele_energy[i] = eleObj->at(i).energy();
        ele_calibError[i] = eleObj->at(i).ecalRegressionError();
        ele_etaSC[i] = eleObj->at(i).superCluster()->eta();
        ele_dEtaIn[i] = eleObj->at(i).deltaEtaSuperClusterTrackAtVtx();
        ele_dPhiIn[i] = eleObj->at(i).deltaPhiSuperClusterTrackAtVtx();
        ele_sigmaIEtaIEta[i] = eleObj->at(i).sigmaIetaIeta();
        ele_HoE[i] = eleObj->at(i).hadronicOverEm();
        if (pvObj->size() > 0) {
            ele_d0[i] = eleObj->at(i).gsfTrack()->dxy(pvObj->at(0).position());
            ele_dZ[i] = eleObj->at(i).gsfTrack()->dz(pvObj->at(0).position());
        }
        else {
            ele_d0[i] = eleObj->at(i).gsfTrack()->dxy();
            ele_dZ[i] = eleObj->at(i).gsfTrack()->dz();
        }
        ele_IoEmIoP[i] = 1.0 / eleObj->at(i).ecalEnergy() -
                         eleObj->at(i).eSuperClusterOverP() / eleObj->at(i).ecalEnergy();
        ele_e1x5[i] = eleObj->at(i).e1x5();
        ele_e2x5[i] = eleObj->at(i).e2x5Max();
        ele_e5x5[i] = eleObj->at(i).e5x5();
        ele_isoCH[i] = eleObj->at(i).chargedHadronIso();
        ele_isoNH[i] = eleObj->at(i).neutralHadronIso();
        ele_isoEM[i] = eleObj->at(i).photonIso();
        double eleEA = ElectronEffectiveArea::GetElectronEffectiveArea(
                       ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03, ele_etaSC[i],
                       ElectronEffectiveArea::kEleEAData2012);
        ele_isoNeutralNoPU[i] = std::max(0.0, ele_isoNH[i] + ele_isoEM[i] -
                                         std::max(0.0, *isoRhoObj) * eleEA);
        ele_iso[i] = (ele_isoCH[i] + ele_isoNeutralNoPU[i]) / ele_pt[i];
        ele_isoDetEM[i] = eleObj->at(i).dr03EcalRecHitSumEt();
        ele_isoDetHad1[i] = eleObj->at(i).dr03HcalDepth1TowerSumEt();
        ele_isoDetTrk[i] = eleObj->at(i).dr03TkSumPt();
        ele_hasConversion[i] = ConversionTools::hasMatchedConversion(eleObj->at(i), convObj,
                                                                     beamObj->position());
        ele_isEcalDriven[i] = eleObj->at(i).ecalDrivenSeed();
        ele_missingHits[i] = eleObj->at(i).gsfTrack()->trackerExpectedHitsInner().numberOfHits();
        ele_goodVetoId[i]   = (
                               (fabs(ele_etaSC[i]) < 1.442)     &&
                               (fabs(ele_dEtaIn[i]) < 0.007)    &&
                               (fabs(ele_dPhiIn[i]) < 0.8)      &&
                               (ele_sigmaIEtaIEta[i] < 0.01)    &&
                               (ele_HoE[i] < 0.15)              &&
                               (fabs(ele_d0[i]) < 0.04)         &&
                               (fabs(ele_dZ[i]) < 0.2)          &&
                               (ele_iso[i] < 0.15)
                              )                                 ||
                              (
                               (fabs(ele_etaSC[i]) > 1.566)     &&
                               (fabs(ele_dEtaIn[i]) < 0.01)     &&
                               (fabs(ele_dPhiIn[i]) < 0.7)      &&
                               (ele_sigmaIEtaIEta[i] < 0.03)    &&
                               (fabs(ele_d0[i]) < 0.04)         &&
                               (fabs(ele_dZ[i]) < 0.2)          &&
                               (ele_iso[i] < 0.15)
                              );
        ele_goodLooseId[i]  = (
                               (fabs(ele_etaSC[i]) < 1.442)     &&
                               (fabs(ele_dEtaIn[i]) < 0.007)    &&
                               (fabs(ele_dPhiIn[i]) < 0.15)     &&
                               (ele_sigmaIEtaIEta[i] < 0.01)    &&
                               (ele_HoE[i] < 0.12)              &&
                               (fabs(ele_d0[i]) < 0.02)         &&
                               (fabs(ele_dZ[i]) < 0.2)          &&
                               (fabs(ele_IoEmIoP[i]) < 0.05)    &&
                               (ele_iso[i] < 0.15)              &&
                               (!ele_hasConversion[i])          &&
                               (ele_missingHits[i] <= 1)
                              )                                 ||
                              (
                               (fabs(ele_etaSC[i]) > 1.566)     &&
                               (fabs(ele_dEtaIn[i]) < 0.009)    &&
                               (fabs(ele_dPhiIn[i]) < 0.10)     &&
                               (ele_sigmaIEtaIEta[i] < 0.03)    &&
                               (ele_HoE[i] < 0.10)              &&
                               (fabs(ele_d0[i]) < 0.02)         &&
                               (fabs(ele_dZ[i]) < 0.2)          &&
                               (fabs(ele_IoEmIoP[i]) < 0.05)    &&
                               (ele_pt[i] > 20.0 ? (ele_iso[i] < 0.15) : (ele_iso[i] < 0.10)) &&
                               (!ele_hasConversion[i])          &&
                               (ele_missingHits[i] <= 1)
                              );
        ele_goodMediumId[i] = (
                               (fabs(ele_etaSC[i]) < 1.442)     &&
                               (fabs(ele_dEtaIn[i]) < 0.004)    &&
                               (fabs(ele_dPhiIn[i]) < 0.06)     &&
                               (ele_sigmaIEtaIEta[i] < 0.01)    &&
                               (ele_HoE[i] < 0.12)              &&
                               (fabs(ele_d0[i]) < 0.02)         &&
                               (fabs(ele_dZ[i]) < 0.1)          &&
                               (fabs(ele_IoEmIoP[i]) < 0.05)    &&
                               (ele_iso[i] < 0.15)              &&
                               (!ele_hasConversion[i])          &&
                               (ele_missingHits[i] <= 1)
                              )                                 ||
                              (
                               (fabs(ele_etaSC[i]) > 1.566)     &&
                               (fabs(ele_dEtaIn[i]) < 0.007)    &&
                               (fabs(ele_dPhiIn[i]) < 0.03)     &&
                               (ele_sigmaIEtaIEta[i] < 0.03)    &&
                               (ele_HoE[i] < 0.10)              &&
                               (fabs(ele_d0[i]) < 0.02)         &&
                               (fabs(ele_dZ[i]) < 0.1)          &&
                               (fabs(ele_IoEmIoP[i]) < 0.05)    &&
                               (ele_pt[i] > 20.0 ? (ele_iso[i] < 0.15) : (ele_iso[i] < 0.10)) &&
                               (!ele_hasConversion[i])          &&
                               (ele_missingHits[i] <= 1)
                              );
        ele_goodTightId[i]  = (
                               (fabs(ele_etaSC[i]) < 1.442)     &&
                               (fabs(ele_dEtaIn[i]) < 0.004)    &&
                               (fabs(ele_dPhiIn[i]) < 0.03)     &&
                               (ele_sigmaIEtaIEta[i] < 0.01)    &&
                               (ele_HoE[i] < 0.12)              &&
                               (fabs(ele_d0[i]) < 0.02)         &&
                               (fabs(ele_dZ[i]) < 0.1)          &&
                               (fabs(ele_IoEmIoP[i]) < 0.05)    &&
                               (ele_iso[i] < 0.10)              &&
                               (!ele_hasConversion[i])          &&
                               (ele_missingHits[i] <= 0)
                              )                                 ||
                              (
                               (fabs(ele_etaSC[i]) > 1.566)     &&
                               (fabs(ele_dEtaIn[i]) < 0.005)    &&
                               (fabs(ele_dPhiIn[i]) < 0.02)     &&
                               (ele_sigmaIEtaIEta[i] < 0.03)    &&
                               (ele_HoE[i] < 0.10)              &&
                               (fabs(ele_d0[i]) < 0.02)         &&
                               (fabs(ele_dZ[i]) < 0.1)          &&
                               (fabs(ele_IoEmIoP[i]) < 0.05)    &&
                               (ele_pt[i] > 20.0 ? (ele_iso[i] < 0.10) : (ele_iso[i] < 0.07)) &&
                               (!ele_hasConversion[i])          &&
                               (ele_missingHits[i] <= 0)
                              );
        double et = eleObj->at(i).et();
        double heep_iso = ele_isoDetEM[i] + ele_isoDetHad1[i] - 0.28 * (*isoRhoObj);
        ele_goodHeepId[i]   = (
                               (fabs(ele_etaSC[i]) < 1.442)     &&
                               (et > 35.0)                      &&
                               (ele_isEcalDriven[i])            &&
                               (fabs(ele_dEtaIn[i]) < 0.005)    &&
                               (fabs(ele_dPhiIn[i]) < 0.06)     &&
                               (ele_HoE[i] < 0.05)              &&
                               ((ele_e1x5[i]/ele_e5x5[i]>0.83)||(ele_e2x5[i]/ele_e5x5[i]>0.94)) &&
                               (heep_iso < 2.0 + 0.03 * et)     &&
                               (ele_isoDetTrk[i] < 5.0)         &&
                               (ele_missingHits[i] <= 1)        &&
                               (fabs(ele_d0[i]) < 0.02)
                              )                                 ||
                              (
                               (fabs(ele_etaSC[i]) > 1.566)     &&
                               (et > 35.0)                      &&
                               (ele_isEcalDriven[i])            &&
                               (fabs(ele_dEtaIn[i]) < 0.007)    &&
                               (fabs(ele_dPhiIn[i]) < 0.06)     &&
                               (ele_HoE[i] < 0.05)              &&
                               (ele_sigmaIEtaIEta[i] < 0.03)    &&
                               (et < 50.0 ? (heep_iso < 2.5) : (heep_iso < 2.5+0.03*(et-50.0))) &&
                               (ele_isoDetTrk[i] < 5.0)         &&
                               (ele_missingHits[i] <= 1)        &&
                               (fabs(ele_d0[i]) < 0.05)
                              );
        ele_mvaTrig[i] = eleObj->at(i).electronID("mvaTrigV0");
        ele_mvaNonTrig[i] = eleObj->at(i).electronID("mvaNonTrigV0");
        ele_vx[i] = eleObj->at(i).vx();
        ele_vy[i] = eleObj->at(i).vy();
        ele_vz[i] = eleObj->at(i).vz();
    }
    
    // HF electron information
    edm::Handle<std::vector<reco::RecoEcalCandidate> > hfEleObj;
    event.getByLabel(hfElectronTag, hfEleObj);
    edm::Handle<reco::HFEMClusterShapeAssociationCollection> hfCluMapObj;
    event.getByLabel(hfClusterMapTag, hfCluMapObj);
    hfEle_size = hfEleObj->size();
    for (unsigned int i = 0; i != hfEleObj->size() && i != 99; ++i) {
        TLorentzVector p4(0.0, 0.0, 0.0, 0.0);
        p4.SetPtEtaPhiE(hfEleObj->at(i).pt(), hfEleObj->at(i).eta(),
                        hfEleObj->at(i).phi(), hfEleObj->at(i).energy());
        edm::InputTag filterTag("hltHFEMPt30TightFilter");
        hfEle_ele23hft30[i] = FindMatch(*trgSumObj, filterTag, p4);
        filterTag = edm::InputTag("hltHFEMTightFilter");
        hfEle_ele27hft15[i] = FindMatch(*trgSumObj, filterTag, p4);
        hfEle_pt[i] = hfEleObj->at(i).pt();
        hfEle_phi[i] = hfEleObj->at(i).phi();
        hfEle_eta[i] = hfEleObj->at(i).eta();
        hfEle_energy[i] = hfEleObj->at(i).energy();
        edm::Ref<std::vector<reco::HFEMClusterShape> > hfCluster;
        hfCluster = hfCluMapObj->find(hfEleObj->at(i).superCluster())->val;
        hfEle_e1x1[i] = hfCluster->e1x1();
        hfEle_e3x3[i] = hfCluster->e3x3();
        hfEle_e5x5[i] = hfCluster->e5x5();
        hfEle_eLong1x1[i] = hfCluster->eLong1x1();
        hfEle_eLong3x3[i] = hfCluster->eLong3x3();
        hfEle_eLong5x5[i] = hfCluster->eLong5x5();
        hfEle_eShort1x1[i] = hfCluster->eShort1x1();
        hfEle_eShort3x3[i] = hfCluster->eShort3x3();
        hfEle_eShort5x5[i] = hfCluster->eShort5x5();
        hfEle_e9e25[i] = hfCluster->e9e25();
        hfEle_eCOREe9[i] = hfCluster->eCOREe9();
        hfEle_eCore[i] = hfCluster->eCore();
        hfEle_eSeL[i] = hfCluster->eSeL();
        hfEle_vx[i] = hfEleObj->at(i).vx();
        hfEle_vy[i] = hfEleObj->at(i).vy();
        hfEle_vz[i] = hfEleObj->at(i).vz();
    }
    
    // PF jet information
    edm::Handle<std::vector<pat::Jet> > pfJetObj;
    event.getByLabel(pfJetTag, pfJetObj);
    unsigned int pfCount(0);
    for (unsigned int i = 0; i != pfJetObj->size() && pfCount != 99; ++i) {
        if (pfJetObj->at(i).pt() < 10.0) continue;
        pfJet_pt[pfCount] = pfJetObj->at(i).pt();
        pfJet_phi[pfCount] = pfJetObj->at(i).phi();
        pfJet_eta[pfCount] = pfJetObj->at(i).eta();
        pfJet_energy[pfCount] = pfJetObj->at(i).energy();
        pfJet_nConstituents[pfCount] = pfJetObj->at(i).nConstituents();
        pfJet_chargedMultiplicity[pfCount] = pfJetObj->at(i).chargedMultiplicity();
        pfJet_chargedEmEnergyFraction[pfCount] = pfJetObj->at(i).chargedEmEnergyFraction();
        pfJet_chargedHadronEnergyFraction[pfCount] = pfJetObj->at(i).chargedHadronEnergyFraction();
        pfJet_neutralEmEnergyFraction[pfCount] = pfJetObj->at(i).neutralEmEnergyFraction();
        pfJet_neutralHadronEnergyFraction[pfCount] = pfJetObj->at(i).neutralHadronEnergyFraction();
        bool goodLowEta = (pfJet_chargedHadronEnergyFraction[i] > 0.0) &&
                          (pfJet_chargedEmEnergyFraction[i] < 0.99)    &&
                          (pfJet_chargedMultiplicity[i] > 0);
        pfJet_goodLooseId[pfCount]  = (pfJet_neutralHadronEnergyFraction[i] < 0.99) &&
                                      (pfJet_neutralEmEnergyFraction[i] < 0.99)     &&
                                      (pfJet_nConstituents[i] > 1)                  &&
                                      (fabs(pfJet_eta[i]) < 2.4 ? goodLowEta : true);
        pfJet_goodMediumId[pfCount] = (pfJet_neutralHadronEnergyFraction[i] < 0.95) &&
                                      (pfJet_neutralEmEnergyFraction[i] < 0.95)     &&
                                      (pfJet_nConstituents[i] > 1)                  &&
                                      (fabs(pfJet_eta[i]) < 2.4 ? goodLowEta : true);
        pfJet_goodTightId[pfCount]  = (pfJet_neutralHadronEnergyFraction[i] < 0.90) &&
                                      (pfJet_neutralEmEnergyFraction[i] < 0.90)     &&
                                      (pfJet_nConstituents[i] > 1)                  &&
                                      (fabs(pfJet_eta[i]) < 2.4 ? goodLowEta : true);
        pfJet_jetArea[pfCount] = pfJetObj->at(i).jetArea();
        pfJet_vx[pfCount] = pfJetObj->at(i).vx();
        pfJet_vy[pfCount] = pfJetObj->at(i).vy();
        pfJet_vz[pfCount] = pfJetObj->at(i).vz();
        ++pfCount;
    }
    pfJet_size = pfCount;
    
    // PF MET information
    edm::Handle<std::vector<pat::MET> > pfMetObj;
    event.getByLabel(pfMetTag, pfMetObj);
    if (pfMetObj->size() != 0) {
        pfMet_et = pfMetObj->at(0).et();
        pfMet_sumEt = pfMetObj->at(0).sumEt();
        pfMet_phi = pfMetObj->at(0).phi();
    }
    else {
        pfMet_et = -9.9;
        pfMet_sumEt = -9.9;
        pfMet_phi = -9.9;
    }
    
	tree->Fill();
    return;
}

DEFINE_FWK_MODULE(TreeProducer);

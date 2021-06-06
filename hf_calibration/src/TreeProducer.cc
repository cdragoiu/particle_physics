#include "HCALCalibration/HFCalibration/interface/TreeProducer.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaReco/interface/HFEMClusterShape.h"
#include "DataFormats/EgammaReco/interface/HFEMClusterShapeAssociation.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>

// constructor -------------------------------------------------------------------------------------
TreeProducer::TreeProducer(const edm::ParameterSet &config) {
    outFileName = config.getParameter<std::string>("outFileName");
    vertexTag = config.getParameter<edm::InputTag>("vertexTag");
    electronTag = config.getParameter<edm::InputTag>("electronTag");
    conversionTag = config.getParameter<edm::InputTag>("conversionTag");
    beamSpotTag = config.getParameter<edm::InputTag>("beamSpotTag");
    hfElectronTag = config.getParameter<edm::InputTag>("hfElectronTag");
    hfClusterMapTag = config.getParameter<edm::InputTag>("hfClusterMapTag");
	outFile = new TFile(outFileName.c_str(), "RECREATE");
	tree = new TTree("Events", "");
	tree->Branch("ev_run", &ev_run, "ev_run/i");
    tree->Branch("ev_lumi", &ev_lumi, "ev_lumi/i");
    tree->Branch("ev_nr", &ev_nr, "ev_nr/i");
    tree->Branch("ev_bunch", &ev_bunch, "ev_bunch/i");
    tree->Branch("pv_size", &pv_size, "pv_size/i");
    tree->Branch("pv_nTracks", pv_nTracks, "pv_nTracks[pv_size]/i");
    tree->Branch("pv_chi2", pv_chi2, "pv_chi2[pv_size]/F");
    tree->Branch("pv_normalizedChi2", pv_normalizedChi2, "pv_normalizedChi2[pv_size]/F");
    tree->Branch("pv_ndof", pv_ndof, "pv_ndof[pv_size]/F");
    tree->Branch("pv_x", pv_x, "pv_x[pv_size]/F");
    tree->Branch("pv_y", pv_y, "pv_y[pv_size]/F");
    tree->Branch("pv_z", pv_z, "pv_z[pv_size]/F");
    tree->Branch("pv_rho", pv_rho, "pv_rho[pv_size]/F");
    tree->Branch("ele_size", &ele_size, "ele_size/i");
    tree->Branch("ele_charge", ele_charge, "ele_charge[ele_size]/I");
    tree->Branch("ele_pt", ele_pt, "ele_pt[ele_size]/F");
    tree->Branch("ele_phi", ele_phi, "ele_phi[ele_size]/F");
    tree->Branch("ele_eta", ele_eta, "ele_eta[ele_size]/F");
    tree->Branch("ele_energy", ele_energy, "ele_energy[ele_size]/F");
    tree->Branch("ele_etaSC", ele_etaSC, "ele_etaSC[ele_size]/F");
    tree->Branch("ele_dEtaIn", ele_dEtaIn, "ele_dEtaIn[ele_size]/F");
    tree->Branch("ele_dPhiIn", ele_dPhiIn, "ele_dPhiIn[ele_size]/F");
    tree->Branch("ele_sigmaIEtaIEta", ele_sigmaIEtaIEta, "ele_sigmaIEtaIEta[ele_size]/F");
    tree->Branch("ele_full5x5SigmaIEtaIEta", ele_full5x5SigmaIEtaIEta,
                 "ele_full5x5SigmaIEtaIEta[ele_size]/F");
    tree->Branch("ele_HoE", ele_HoE, "ele_HoE[ele_size]/F");
    tree->Branch("ele_d0", ele_d0, "ele_d0[ele_size]/F");
    tree->Branch("ele_dZ", ele_dZ, "ele_dZ[ele_size]/F");
    tree->Branch("ele_IoEmIoP", ele_IoEmIoP, "ele_IoEmIoP[ele_size]/F");
    tree->Branch("ele_iso", ele_iso, "ele_iso[ele_size]/F");
    tree->Branch("ele_hasConversion", ele_hasConversion, "ele_hasConversion[ele_size]/O");
    tree->Branch("ele_missingHits", ele_missingHits, "ele_missingHits[ele_size]/i");
    tree->Branch("hfEle_size", &hfEle_size, "hfEle_size/i");
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
}

// destructor --------------------------------------------------------------------------------------
TreeProducer::~TreeProducer() {
	outFile->cd();
	tree->Write();
	tree->Print();
	delete tree;
	delete outFile;
}

// analyze event -----------------------------------------------------------------------------------
void TreeProducer::analyze(const edm::Event &event, const edm::EventSetup &setup) {
    
    // event information
    ev_run = event.id().run();
    ev_lumi = event.luminosityBlock();
	ev_nr = event.id().event();
    ev_bunch = event.bunchCrossing();
    
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
    
    // electron information
    edm::Handle<std::vector<reco::GsfElectron> > eleObj;
    event.getByLabel(electronTag, eleObj);
    edm::Handle<std::vector<reco::Conversion> > convObj;
    event.getByLabel(conversionTag, convObj);
    edm::Handle<reco::BeamSpot> beamObj;
    event.getByLabel(beamSpotTag, beamObj);
    ele_size = eleObj->size();
    for (unsigned int i = 0; i != eleObj->size() && i != 99; ++i) {
        ele_charge[i] = eleObj->at(i).charge();
        ele_pt[i] = eleObj->at(i).pt();
        ele_phi[i] = eleObj->at(i).phi();
        ele_eta[i] = eleObj->at(i).eta();
        ele_energy[i] = eleObj->at(i).energy();
        ele_etaSC[i] = eleObj->at(i).superCluster()->eta();
        ele_dEtaIn[i] = eleObj->at(i).deltaEtaSuperClusterTrackAtVtx();
        ele_dPhiIn[i] = eleObj->at(i).deltaPhiSuperClusterTrackAtVtx();
        ele_sigmaIEtaIEta[i] = eleObj->at(i).sigmaIetaIeta();
        ele_full5x5SigmaIEtaIEta[i] = eleObj->at(i).full5x5_sigmaIetaIeta();
        ele_HoE[i] = eleObj->at(i).hcalOverEcal();
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
        double sumCH = eleObj->at(i).pfIsolationVariables().sumChargedHadronPt;
        double sumNH = eleObj->at(i).pfIsolationVariables().sumNeutralHadronEt;
        double sumEM = eleObj->at(i).pfIsolationVariables().sumPhotonEt;
        double sumPU = eleObj->at(i).pfIsolationVariables().sumPUPt;
        ele_iso[i] = (sumCH + std::max(0.0, sumNH + sumEM - 0.5 * sumPU)) / ele_pt[i];
        ele_hasConversion[i] = ConversionTools::hasMatchedConversion(eleObj->at(i), convObj,
                                                                     beamObj->position());
        ele_missingHits[i] = eleObj->at(i).gsfTrack()->numberOfLostHits();
    }
    
    // HF electron information
    edm::Handle<std::vector<reco::RecoEcalCandidate> > hfEleObj;
    event.getByLabel(hfElectronTag, hfEleObj);
    edm::Handle<reco::HFEMClusterShapeAssociationCollection> hfCluMapObj;
    event.getByLabel(hfClusterMapTag, hfCluMapObj);
    hfEle_size = hfEleObj->size();
    for (unsigned int i = 0; i != hfEleObj->size() && i != 99; ++i) {
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
    }
    
	tree->Fill();
    return;
}

DEFINE_FWK_MODULE(TreeProducer);

#ifndef RecoParticleFlow_PFPatProducer_METScan_
#define RecoParticleFlow_PFPatProducer_METScan_

// system include files
#include <memory>
#include <string>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PFParticle.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"

#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgoRcd.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"

#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/METReco/interface/HcalNoiseSummary.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/PFClusterMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/PFClusterMETCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/METReco/interface/HcalNoiseSummary.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Framework/interface/EventSetup.h"




#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"




#include <math.h>

#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TVector3.h>
#include <fstream>

using namespace std;
using namespace edm;
using namespace reco;

struct BadRuns { 
        unsigned long Run_Bad;
        unsigned long LS_Bad;
        unsigned long event_Bad;
};

class METScanningNtupleMaker : public edm::EDAnalyzer {
	public:

		explicit METScanningNtupleMaker(const edm::ParameterSet&);

		~METScanningNtupleMaker();

		virtual void analyze(const edm::Event&, const edm::EventSetup&);

		virtual void beginRun(const edm::Run & r, const edm::EventSetup & c);

	private:
		edm::EDGetTokenT<reco::MuonCollection> Muon_token;
		edm::EDGetTokenT<reco::PFCandidateCollection> PfCandidates_token;
		edm::EDGetTokenT<reco::PFJetCollection> PfJets_token;
		edm::EDGetTokenT<reco::CaloMETCollection> CaloMET_token;
		edm::EDGetTokenT<reco::PFMETCollection> PFCaloMET_token;
		edm::EDGetTokenT<reco::PFClusterMETCollection> PFClusterMET_token;
		edm::EDGetTokenT<reco::PFMETCollection> PFMET_token;
		edm::EDGetTokenT<reco::PFClusterCollection> EcalPFClusters_token;
		edm::EDGetTokenT<reco::PFClusterCollection> HcalPFClusters_token;
		edm::EDGetTokenT<reco::PFClusterCollection> HBHEPFClusters_token;
		edm::EDGetTokenT<reco::PFClusterCollection> HOPFClusters_token;
		edm::EDGetTokenT<reco::PFClusterCollection> HFPFClusters_token;
		edm::EDGetTokenT<reco::TrackCollection> Tracks_token;
		edm::EDGetTokenT<bool> TrackingLETMC_token;
		edm::EDGetTokenT<bool> TrackingLETMS_token;
		edm::EDGetTokenT<bool> TrackingMSC_token;
		edm::EDGetTokenT<bool> TrackingTMSC_token;
		//edm::EDGetTokenT<bool> CSC2015_token;
		edm::EDGetTokenT<bool> GlobalTightHalo2016_token;
		edm::EDGetTokenT<bool> GlobalSuperTightHalo2016_token;
		edm::EDGetTokenT<bool> HcalStripHalo_token;
		edm::EDGetTokenT<bool> HBHER1_token;
		edm::EDGetTokenT<bool> HBHER2L_token;
		edm::EDGetTokenT<bool> HBHER2T_token;
		edm::EDGetTokenT<bool> ECALTP_token;
		edm::EDGetTokenT<bool> ECALSC_token;
		edm::EDGetTokenT<EcalRecHitCollection> RecHitsEB_token;
		edm::EDGetTokenT<EcalRecHitCollection> RecHitsEE_token;
		edm::EDGetTokenT<EcalRecHitCollection> RecHitsES_token;
		edm::EDGetTokenT<HcalNoiseSummary> hSummary_token;
		edm::EDGetTokenT<bool> BadChCandF_token;

		edm::EDGetTokenT<bool> BadPFMuon_token;

		//edm::EDGetTokenT<bool> BadChCandFOld_token;

		//edm::EDGetTokenT<bool> BadPFMuonOld_token;

		edm::EDGetTokenT<bool> EcalBadCalib_token;

		edm::EDGetTokenT<vector<reco::Vertex> >  vertex_token;

		bool isReco_token;

		size_t run,event,lumiBlock,time;
		//bool filtercsc2015
		//bool filterbadChCandidateOld, filterbadPFMuonOld
		bool filterglobaltighthalo2016,filterglobalsupertighthalo2016, filterhcalstriphalo, filterhbher1, filterhbher2l, filterhbher2t, filterhbher1nozeros, filterhbheiso, filterecaltp, filterecalsc; 
		bool filtertrackingletmc, filtertrackingletms, filtertrackingmsc, filtertrackingtmsc, filterbadChCandidate, filterbadPFMuon, filterEcalBadCalib  ;
		edm::RunNumber_t irun;
		edm::EventNumber_t ievent;
		edm::LuminosityBlockNumber_t ilumiBlock;
		edm::Timestamp itime;

		size_t nVtx;

		std::vector<bool>   muon_PF;
		std::vector<float>  muon_pt;
		std::vector<float>  muon_eta;
		std::vector<float>  muon_phi;
		std::vector<float>  muon_ptError;
		std::vector<float>  imuon_pt;
		std::vector<float>  imuon_eta;
		std::vector<float>  imuon_phi;
		std::vector<float>  imuon_ptError;
		std::vector<float>  muon_SC;

		std::vector<float>  pfLepton_pt;
		std::vector<float>  pfLepton_eta;
		std::vector<float>  pfLepton_phi;
		std::vector<float>  pfLepton_pdgId;

		std::vector<float>  pfHadron_pt;
		std::vector<float>  pfHadron_eta;
		std::vector<float>  pfHadron_phi;
		std::vector<float>  pfHadron_pdgId;  


		std::vector<float>  pfJet_pt;
		std::vector<float>  pfJet_eta;
		std::vector<float>  pfJet_phi;
		std::vector<float>  pfJet_looseId;
		std::vector<float>  pfJet_tightId;
		std::vector<float>  pfJet_tlvId;
		int pfJet_hpfl;
		int pfJet_hpft;
		int pfJet_hpfv;

		float caloMETPt;
		float caloMETPhi;
		float caloMETSumEt;

		float pfCaloMETPt;
		float pfCaloMETPhi;
		float pfCaloMETSumEt;

		float pfClusterMETPt;
		float pfClusterMETPhi;
		float pfClusterMETSumEt;

		float pfMETPt;
		float pfMETPhi;
		float pfMETSumEt;

		std::vector<float>  pfClusterEcal_energy;
		std::vector<float>  pfClusterEcal_time;
		std::vector<float>  pfClusterEcal_eta;
		std::vector<float>  pfClusterEcal_phi;
		std::vector<float>  pfClusterEcal_status13;
		std::vector<float>  pfClusterEcal_status14;

		std::vector<float>  pfClusterHcal_energy;
		std::vector<float>  pfClusterHcal_time;
		std::vector<float>  pfClusterHcal_eta;
		std::vector<float>  pfClusterHcal_phi;

		std::vector<float>  pfClusterHBHE_energy;
		std::vector<float>  pfClusterHBHE_time;
		std::vector<float>  pfClusterHBHE_eta;
		std::vector<float>  pfClusterHBHE_phi;

		std::vector<float>  pfClusterHO_energy;
		std::vector<float>  pfClusterHO_time;
		std::vector<float>  pfClusterHO_eta;
		std::vector<float>  pfClusterHO_phi;

		std::vector<float>  pfClusterHF_energy;
		std::vector<float>  pfClusterHF_time;
		std::vector<float>  pfClusterHF_eta;
		std::vector<float>  pfClusterHF_phi;

		std::vector<float> track_ptError;
		std::vector<float> track_pt;
		std::vector<float> track_eta;
		std::vector<float> track_phi;
		std::vector<float> track_d0;
		std::vector<float> track_d0Error;
		std::vector<float> track_dz;
		std::vector<float> track_dzError;

		vector<double> _muonEnergy;
		vector<double> _muonInnerTrackPt;
		vector<double> _muonInnerTrackPtError;
		vector<double> _muonBestTrackPt;
		vector<double> _muonBestTrackPtError;
		vector<bool> _muonGlobal;
		vector<int> _muonalgo;
		vector<double> _muonsegComp;
		vector<bool> _muonHighpurity;
		vector<bool> _muonPF;
		vector<int> _muonInnerTrackCharge;
		vector<int> _muonBestTrackCharge;

		vector<double> _muonGlobalTrackPt;
		vector<double> _muonGlobalTrackPtError;
		vector<int> _muonGlobalTrackCharge;

		vector<double> _muonTPFMSTrackPt;
		vector<double> _muonTPFMSTrackPtError;
		vector<int> _muonTPFMSTrackCharge;

		vector<double> _muonPickyTrackPt;
		vector<double> _muonPickyTrackPtError;
		vector<int> _muonPickyTrackCharge;


		vector<double> _muonDYTTrackPt;
		vector<double> _muonDYTTrackPtError;
		vector<int> _muonDYTTrackCharge;

		vector<double> _muonTunePMuonBestTrackPt;
		vector<double> _muonTunePMuonBestTrackPtError;
		vector<int> _muonTunePMuonBestTrackCharge;

		vector<double> _muonTrackIsoSumPt;

		vector<int> _muonValidPixelHits;
		vector<int> _muonValidMuonHits;
		vector<int> _muonTrackerLayersHits;
		vector<int> _muonMatchedStations;
		vector<double> _muonDxy;
		vector<double> _muonDz;

		vector<int> _muonHighPtTrackerId;
		vector<int> _muonHighPtId;


		vector<int> _muonLoose;
		vector<int> _muonMedium;
		vector<int> _muonTight;
		vector<int> _muonSoft;
		vector<int> _muonTracker;
		vector<int> _muonMatchedNStations;
		vector<int> _muonStationMask;
		vector<int> _muonRPCLayers;
		vector<int> _muonBestValidMuonHits;
		//tree stuff
		std::string outputfile_;
		TFile* tf1;
		TTree* s;

		string file_INPUT;
		std::ifstream badrunFile;
		std::vector<BadRuns> BadRunsList;

		std::vector<std::string> split(const std::string &s, char delim);

};

#endif

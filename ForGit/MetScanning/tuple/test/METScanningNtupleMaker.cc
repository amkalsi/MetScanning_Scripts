


#include "MetScanning/tuple/test/METScanningNtupleMaker.h"
#include <iostream>

//User

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

METScanningNtupleMaker::METScanningNtupleMaker(const edm::ParameterSet& iConfig) {


	//the input tags
	isReco_token         = iConfig.getParameter<bool>("isReco");
	Muon_token           = consumes<reco::MuonCollection>(iConfig.getParameter<edm::InputTag>("muonCandidates"           ));
	PfCandidates_token   = consumes<reco::PFCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCandidates"           ));
	PfJets_token         = consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("pfJets"                 ));
	CaloMET_token        = consumes<reco::CaloMETCollection>(iConfig.getParameter<edm::InputTag>("caloMET"));
	PFCaloMET_token      = consumes<reco::PFMETCollection>(iConfig.getParameter<edm::InputTag>("pfCaloMET"));
	PFClusterMET_token   = consumes<reco::PFClusterMETCollection>(iConfig.getParameter<edm::InputTag>("pfClusterMET"));
	PFMET_token          = consumes<reco::PFMETCollection>(iConfig.getParameter<edm::InputTag>("pfMET"                  ));
	EcalPFClusters_token = consumes<reco::PFClusterCollection>(iConfig.getParameter<edm::InputTag>("EcalPFClusterCollection"));
	HcalPFClusters_token = consumes<reco::PFClusterCollection>(iConfig.getParameter<edm::InputTag>("HcalPFClusterCollection"));
	HBHEPFClusters_token = consumes<reco::PFClusterCollection>(iConfig.getParameter<edm::InputTag>("HBHEPFClusterCollection"));
	HOPFClusters_token   = consumes<reco::PFClusterCollection>(iConfig.getParameter<edm::InputTag>("HOPFClusterCollection"));
	HFPFClusters_token   = consumes<reco::PFClusterCollection>(iConfig.getParameter<edm::InputTag>("HFPFClusterCollection"));
	Tracks_token         = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracksCollection"));
	TrackingLETMC_token  = consumes<bool>(iConfig.getParameter<edm::InputTag>("TRKfilterLETMC"         ));
	TrackingLETMS_token  = consumes<bool>(iConfig.getParameter<edm::InputTag>("TRKfilterLETMS"         ));
	TrackingMSC_token    = consumes<bool>(iConfig.getParameter<edm::InputTag>("TRKfilterMSC"           ));
	TrackingTMSC_token   = consumes<bool>(iConfig.getParameter<edm::InputTag>("TRKfilterTMSC"          ));
	//CSC2015_token = consumes<bool>(iConfig.getParameter<edm::InputTag>("CSC2015filter"));
	GlobalTightHalo2016_token = consumes<bool>(iConfig.getParameter<edm::InputTag>("GlobalHalofilterTight"));
	GlobalSuperTightHalo2016_token = consumes<bool>(iConfig.getParameter<edm::InputTag>("GlobalHalofilterSuperTight"));
	HcalStripHalo_token = consumes<bool>(iConfig.getParameter<edm::InputTag>("HcalStripHaloFilter"));
	HBHER1_token = consumes<bool>(iConfig.getParameter<edm::InputTag>("HBHEfilterR1"));
	HBHER2L_token = consumes<bool>(iConfig.getParameter<edm::InputTag>("HBHEfilterR2L"));
	HBHER2T_token = consumes<bool>(iConfig.getParameter<edm::InputTag>("HBHEfilterR2T"));
	ECALTP_token = consumes<bool>(iConfig.getParameter<edm::InputTag>("ECALTPfilter"));
	ECALSC_token = consumes<bool>(iConfig.getParameter<edm::InputTag>("ECALSCfilter"));
	RecHitsEB_token = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("EBRecHits"));
	RecHitsEE_token = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("EERecHits"));
	RecHitsES_token = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("ESRecHits"));
	hSummary_token = consumes<HcalNoiseSummary>(iConfig.getParameter<edm::InputTag>("HcalNoise"));
	BadChCandF_token = consumes<bool>(iConfig.getParameter<edm::InputTag>("BadChCandFilter"));
	BadPFMuon_token = consumes<bool>(iConfig.getParameter<edm::InputTag>("BadPFMuon"));

	//BadChCandFOld_token = consumes<bool>(iConfig.getParameter<edm::InputTag>("BadChCandSummer16Filter"));
	//BadPFMuonOld_token = consumes<bool>(iConfig.getParameter<edm::InputTag>("BadPFMuonOld"));
	EcalBadCalib_token = consumes<bool>(iConfig.getParameter<edm::InputTag>("EcalBadCalibFilter"));
	vertex_token = consumes<vector<reco::Vertex> >(iConfig.getParameter<edm::InputTag>("OfflinePrimaryVertices"));

	file_INPUT = iConfig.getParameter<string>("file_INPUT");

	// The root tuple
	outputfile_ = iConfig.getParameter<std::string>("rootOutputFile"); 
	tf1 = new TFile(outputfile_.c_str(), "RECREATE");  
	s = new TTree("tree","tree");



	//basic informations ==========================
	s->Branch("run",&run,"run/l");  
	s->Branch("lumi",&lumiBlock,"lumi/l");  
	s->Branch("event",&event,"event/l");  
	s->Branch("time",&time,"time/l");

	s->Branch("nVtx", &nVtx, "nVtx/l");
	//s->Branch("filter_csc2015",&filtercsc2015,"filter_csc2015/O");
	s->Branch("filter_globaltighthalo2016",&filterglobaltighthalo2016,"filter_globaltighthalo2016/O");
	s->Branch("filter_globalsupertighthalo2016",&filterglobalsupertighthalo2016,"filter_globalsupertighthalo2016/O");
	s->Branch("filter_hcalstriphalo",&filterhcalstriphalo,"filter_hcalstriphalo/O");
	s->Branch("filter_hbher1",&filterhbher1,"filter_hbher1/O");
	s->Branch("filter_hbher2l",&filterhbher2l,"filter_hbher2l/O");
	s->Branch("filter_hbher2t",&filterhbher2t,"filter_hbher2t/O");
	s->Branch("filter_hbheiso",&filterhbheiso,"filter_hbheiso/O");
	s->Branch("filter_ecaltp",&filterecaltp,"filter_ecaltp/O");
	s->Branch("filter_ecalsc",&filterecalsc,"filter_ecalsc/O");
	s->Branch("filter_badChCand",&filterbadChCandidate,"filter_badChCand/O");
	s->Branch("filter_badPFMuon",&filterbadPFMuon,"filter_badPFMuon/O");

	//s->Branch("filter_badChCandOld",&filterbadChCandidateOld,"filter_badChCandOld/O");
	//s->Branch("filter_badPFMuonOld",&filterbadPFMuonOld,"filter_badPFMuonOld/O");

	s->Branch("filter_EcalBadCalib",&filterEcalBadCalib,"filter_EcalBadCalib/O");

	//pfLeptons =====================================
	s->Branch("pfLepton_pt"             , &pfLepton_pt   );  
	s->Branch("pfLepton_eta"            , &pfLepton_eta  ); 
	s->Branch("pfLepton_phi"            , &pfLepton_phi  );  
	s->Branch("pfLepton_pdgId"          , &pfLepton_pdgId);


	//pfHadrons ==================================

	s->Branch("pfHadron_pt"             , &pfHadron_pt   );
	s->Branch("pfHadron_eta"            , &pfHadron_eta  );
	s->Branch("pfHadron_phi"            , &pfHadron_phi  );
	s->Branch("pfHadron_pdgId"          , &pfHadron_pdgId);

	//Muons  ======================================


	s->Branch("muon_PF"             , &muon_PF);
	s->Branch("muon_pt"             , &muon_pt   );
	s->Branch("muon_eta"             , &muon_eta   );
	s->Branch("muon_phi"             , &muon_phi   );
	s->Branch("muon_ptError"        , &muon_ptError   );
	s->Branch("imuon_pt"             , &imuon_pt   );
	s->Branch("imuon_eta"             , &imuon_eta   );
	s->Branch("imuon_phi"             , &imuon_phi   );
	s->Branch("imuon_ptError"        , &imuon_ptError   );
	s->Branch("muon_SC"              , &muon_SC    );


	//s->Branch("pfLepton_d0"             , &pfLepton_d0   );
	//Jets ========================================
	s->Branch("pfJet_pt"                , &pfJet_pt      );  
	s->Branch("pfJet_eta"               , &pfJet_eta     ); 
	s->Branch("pfJet_phi"               , &pfJet_phi     );  
	s->Branch("pfJet_looseId"           , &pfJet_looseId );  
	s->Branch("pfJet_tightId"           , &pfJet_tightId );  
	s->Branch("pfJet_tightLepVetoId"    , &pfJet_tlvId   );  
	s->Branch("pfJet_highestPtFailLoose", &pfJet_hpfl    );  
	s->Branch("pfJet_highestPtFailTight", &pfJet_hpft    );  
	s->Branch("pfJet_highestPtFailTLV"  , &pfJet_hpfv    );  

	//METs ========================================
	s->Branch("caloMETPt",&caloMETPt,"caloMETPt/F");  
	s->Branch("caloMETPhi",&caloMETPhi,"caloMETPhi/F"); 
	s->Branch("caloMETSumEt",&caloMETSumEt,"caloMETSumEt/F");  

	if(isReco_token){
		s->Branch("pfCaloMETPt",&pfCaloMETPt,"pfCaloMETPt/F");  
		s->Branch("pfCaloMETPhi",&pfCaloMETPhi,"pfCaloMETPhi/F"); 
		s->Branch("pfCaloMETSumEt",&pfCaloMETSumEt,"pfCaloMETSumEt/F");  

		s->Branch("pfClusterMETPt",&pfClusterMETPt,"pfClusterMETPt/F");  
		s->Branch("pfClusterMETPhi",&pfClusterMETPhi,"pfClusterMETPhi/F"); 
		s->Branch("pfClusterMETSumEt",&pfClusterMETSumEt,"pfClusterMETSumEt/F"); 

		//clusters ==========================================
		s->Branch("pfClusterEcal_energy",&pfClusterEcal_energy);
		s->Branch("pfClusterEcal_time",&pfClusterEcal_time);
		s->Branch("pfClusterEcal_eta",&pfClusterEcal_eta);
		s->Branch("pfClusterEcal_phi",&pfClusterEcal_phi);
		s->Branch("pfClusterEcal_status13",&pfClusterEcal_status13);
		s->Branch("pfClusterEcal_status14",&pfClusterEcal_status14);

		s->Branch("pfClusterHcal_energy",&pfClusterHcal_energy);
		s->Branch("pfClusterHcal_time",&pfClusterHcal_time);
		s->Branch("pfClusterHcal_eta",&pfClusterHcal_eta);
		s->Branch("pfClusterHcal_phi",&pfClusterHcal_phi);

		s->Branch("pfClusterHBHE_energy",&pfClusterHBHE_energy);
		s->Branch("pfClusterHBHE_time",&pfClusterHBHE_time);
		s->Branch("pfClusterHBHE_eta",&pfClusterHBHE_eta);
		s->Branch("pfClusterHBHE_phi",&pfClusterHBHE_phi);

		s->Branch("pfClusterHO_energy",&pfClusterHO_energy);
		s->Branch("pfClusterHO_time",&pfClusterHO_time);
		s->Branch("pfClusterHO_eta",&pfClusterHO_eta);
		s->Branch("pfClusterHO_phi",&pfClusterHO_phi);

		s->Branch("pfClusterHF_energy",&pfClusterHF_energy);
		s->Branch("pfClusterHF_time",&pfClusterHF_time);
		s->Branch("pfClusterHF_eta",&pfClusterHF_eta);
		s->Branch("pfClusterHF_phi",&pfClusterHF_phi);

	}

	s->Branch("pfMETPt",&pfMETPt,"pfMETPt/F");  
	s->Branch("pfMETPhi",&pfMETPhi,"pfMETPhi/F"); 
	s->Branch("pfMETSumEt",&pfMETSumEt,"pfMETSumEt/F"); 


	//tracks ============================================
	s->Branch("track_pt",&track_pt);
	s->Branch("track_eta",&track_eta);
	s->Branch("track_phi",&track_phi);
	s->Branch("track_d0",&track_d0);
	s->Branch("track_d0Error" ,&track_d0Error);
	s->Branch("track_dz",&track_dz);
	s->Branch("track_dzError" ,&track_dzError);
	s->Branch("track_ptError",&track_ptError);


	s->Branch("_muonEnergy", &_muonEnergy);
	s->Branch("_muonInnerTrackPt", & _muonInnerTrackPt);
	s->Branch("_muonInnerTrackPtError", &_muonInnerTrackPtError);
	s->Branch("_muonBestTrackPt", & _muonBestTrackPt);
	s->Branch("_muonBestTrackPtError",&_muonBestTrackPtError);
	s->Branch("_muonInnerTrackCharge", &_muonInnerTrackCharge);
	s->Branch("_muonBestTrackCharge", & _muonBestTrackCharge);

	s->Branch("_muonGlobal", &_muonGlobal);
	s->Branch("_muonalgo", &_muonalgo);
	s->Branch("_muonsegComp", &_muonsegComp);
	s->Branch("_muonHighpurity", &_muonHighpurity);
	s->Branch("_muonPF", &_muonPF);

	s->Branch("_muonGlobalTrackPt", & _muonGlobalTrackPt);
	s->Branch("_muonTPFMSTrackPt", & _muonTPFMSTrackPt);
	s->Branch("_muonPickyTrackPt", & _muonPickyTrackPt);
	s->Branch("_muonDYTTrackPt", & _muonDYTTrackPt);
	s->Branch("_muonTunePMuonBestTrackPt", & _muonTunePMuonBestTrackPt);


	s->Branch("_muonGlobalTrackPtError", & _muonGlobalTrackPtError);
	s->Branch("_muonTPFMSTrackPtError", & _muonTPFMSTrackPtError);
	s->Branch("_muonPickyTrackPtError", & _muonPickyTrackPtError);
	s->Branch("_muonDYTTrackPtError", & _muonDYTTrackPtError);
	s->Branch("_muonTunePMuonBestTrackPtError", & _muonTunePMuonBestTrackPtError);

	s->Branch("_muonGlobalTrackCharge", & _muonGlobalTrackCharge);
	s->Branch("_muonTPFMSTrackCharge", & _muonTPFMSTrackCharge);
	s->Branch("_muonPickyTrackCharge", & _muonPickyTrackCharge);
	s->Branch("_muonDYTTrackCharge", & _muonDYTTrackCharge);
	s->Branch("_muonTunePMuonBestTrackCharge", & _muonTunePMuonBestTrackCharge);
	s->Branch("_muonTrackIsoSumPt", &_muonTrackIsoSumPt);
	s->Branch("_muonValidPixelHits", &_muonValidPixelHits);



	s->Branch("_muonValidMuonHits", &_muonValidMuonHits);
	s->Branch("_muonTrackerLayersHits", &_muonTrackerLayersHits);
	s->Branch("_muonMatchedStations", &_muonMatchedStations);

	s->Branch("_muonDxy", &_muonDxy);
	s->Branch("_muonDz", &_muonDz);

	s->Branch("_muonHighPtTrackerId", &_muonHighPtTrackerId);
	s->Branch("_muonHighPtId", & _muonHighPtId);


	s->Branch("_muonLoose", &_muonLoose);
	s->Branch("_muonMedium", &_muonMedium);
	s->Branch("_muonTight",&_muonTight);
	s->Branch("_muonSoft",&_muonSoft);
	s->Branch("_muonTracker",&_muonTracker);
	s->Branch("_muonMatchedNStations",&_muonMatchedNStations);
	s->Branch("_muonStationMask",&_muonStationMask);
	s->Branch("_muonRPCLayers",&_muonRPCLayers);
	s->Branch("_muonBestValidMuonHits",&_muonBestValidMuonHits);

	std::string line;
	BadRunsList.clear();


	badrunFile.open(file_INPUT);

	while (getline(badrunFile, line))
	{
		std::vector<std::string> columns = split(line,':');
		BadRuns bad_run;

		bad_run.Run_Bad     = std::stoll(columns[0]);
		bad_run.LS_Bad     = std::stoll(columns[1]);
		bad_run.event_Bad   = std::stoll(columns[2]);

		BadRunsList.push_back(bad_run);
	}

	badrunFile.close();

}


METScanningNtupleMaker::~METScanningNtupleMaker() { 

	tf1->cd();
	s->Write();
	tf1->Write();
	tf1->Close();  
}



void 
METScanningNtupleMaker::beginRun(const edm::Run& run, 
		const edm::EventSetup & es) { }


void 
METScanningNtupleMaker::analyze(const Event& iEvent, 
		const EventSetup& iSetup) {


	edm::ESHandle<EcalSeverityLevelAlgo> sevlv;
	iSetup.get<EcalSeverityLevelAlgoRcd>().get(sevlv);

	irun  = iEvent.id().run();
	ievent  = iEvent.id().event();
	ilumiBlock = iEvent.id().luminosityBlock();
	itime = iEvent.time();

	run = (size_t)irun;
	event = (size_t)ievent;
	lumiBlock = (size_t)ilumiBlock;
	time = (size_t)((iEvent.time().value())>>32);

	bool run_found;
	run_found=false;
	for (const BadRuns & m: BadRunsList)
	{
		if(irun == m.Run_Bad && ilumiBlock == m.LS_Bad && ievent == m.event_Bad ) { run_found=true; std::cout<<"is matched:"<< irun << endl;}

	}

	if(!run_found) return;

	Handle<vector<reco::Vertex> > offlineVtx;
	iEvent.getByToken(vertex_token, offlineVtx);
	nVtx = (size_t) offlineVtx->size();
	//	std::cout<<" nVtx = "<<offlineVtx->size()
	//		<<" " <<nVtx<<std::endl;




	//get filters
	/*
	   Handle<bool> ifiltertrackingletmc;
	   iEvent.getByToken(TrackingLETMC_token, ifiltertrackingletmc);
	   filtertrackingletmc = *ifiltertrackingletmc;

	   Handle<bool> ifiltertrackingletms;
	   iEvent.getByToken(TrackingLETMS_token, ifiltertrackingletms);
	   filtertrackingletms = *ifiltertrackingletms;

	   Handle<bool> ifiltertrackingmsc;
	   iEvent.getByToken(TrackingMSC_token, ifiltertrackingmsc);
	   filtertrackingmsc = *ifiltertrackingmsc;

	   Handle<bool> ifiltertrackingtmsc;
	   iEvent.getByToken(TrackingTMSC_token, ifiltertrackingtmsc);
	   filtertrackingtmsc = *ifiltertrackingtmsc;
	   */

	//Drop CSC2015 Filter ~ Isabell Suggestion
	//Handle<bool> ifiltercsc2015;
	//iEvent.getByToken(CSC2015_token, ifiltercsc2015);
	//filtercsc2015 = *ifiltercsc2015;

	Handle<bool> ifilterglobaltighthalo2016;
	iEvent.getByToken(GlobalTightHalo2016_token, ifilterglobaltighthalo2016);
	filterglobaltighthalo2016 = *ifilterglobaltighthalo2016;

	Handle<bool> ifilterglobalsupertighthalo2016;
	iEvent.getByToken(GlobalSuperTightHalo2016_token, ifilterglobalsupertighthalo2016);
	filterglobalsupertighthalo2016 = *ifilterglobalsupertighthalo2016;

	Handle<bool> ifilterhcalstriphalo;
	iEvent.getByToken(HcalStripHalo_token, ifilterhcalstriphalo);
	filterhcalstriphalo = *ifilterhcalstriphalo;

	Handle<bool> ifilterhbher1;
	iEvent.getByToken(HBHER1_token, ifilterhbher1);
	filterhbher1 = *ifilterhbher1;

	Handle<bool> ifilterhbher2l;
	iEvent.getByToken(HBHER2L_token, ifilterhbher2l);
	filterhbher2l = *ifilterhbher2l;

	Handle<bool> ifilterhbher2t;
	iEvent.getByToken(HBHER2T_token, ifilterhbher2t);
	filterhbher2t = *ifilterhbher2t;

	//Handle<bool> ifilterhbheiso;
	//iEvent.getByToken(HBHEISO_token, ifilterhbheiso);
	//filterhbheiso = *ifilterhbheiso;

	Handle<bool> ifilterecaltp;
	iEvent.getByToken(ECALTP_token, ifilterecaltp);
	filterecaltp = *ifilterecaltp;

	Handle<bool> ifilterecalsc;
	iEvent.getByToken(ECALSC_token, ifilterecalsc);
	filterecalsc = *ifilterecalsc;

	edm::Handle<HcalNoiseSummary> hSummary;
	iEvent.getByToken(hSummary_token, hSummary);

	filterhbheiso = true;
	if( hSummary->numIsolatedNoiseChannels() >= 10 ) filterhbheiso = false;
	if( hSummary->isolatedNoiseSumE()        >= 50 ) filterhbheiso = false;
	if( hSummary->isolatedNoiseSumEt()       >= 25 ) filterhbheiso = false;

	filterhbher1nozeros = true;
	if( hSummary->maxHPDHits()               >= 17                           ) filterhbher1nozeros = false;
	if( hSummary->maxHPDNoOtherHits()        >= 10                           ) filterhbher1nozeros = false;
	if( hSummary->HasBadRBXTS4TS5() && !hSummary->goodJetFoundInLowBVRegion()) filterhbher1nozeros = false;

	Handle<bool> ifilterbadChCand;
	iEvent.getByToken(BadChCandF_token, ifilterbadChCand);
	filterbadChCandidate = *ifilterbadChCand;

	Handle<bool> ifilterbadPFMuon;
	iEvent.getByToken(BadPFMuon_token, ifilterbadPFMuon);
	filterbadPFMuon = *ifilterbadPFMuon;

	//Handle<bool> ifilterbadChCandOld;
	//iEvent.getByToken(BadChCandFOld_token, ifilterbadChCandOld);
	//filterbadChCandidateOld = *ifilterbadChCandOld;

	//Handle<bool> ifilterbadPFMuonOld;
	//iEvent.getByToken(BadPFMuonOld_token, ifilterbadPFMuonOld);
	//filterbadPFMuonOld = *ifilterbadPFMuonOld;

	Handle<bool> ifilterEcalBadCalib;
	iEvent.getByToken(EcalBadCalib_token,ifilterEcalBadCalib);
	filterEcalBadCalib = *ifilterEcalBadCalib;

	// get Leptons
	Handle<reco::PFCandidateCollection> pfCandidates;
	iEvent.getByToken(PfCandidates_token,pfCandidates);


	pfLepton_pt     .clear();
	pfLepton_eta    .clear();
	pfLepton_phi    .clear();
	pfLepton_pdgId  .clear();

	// get charged Hadrons (pdgId=211)

	pfHadron_pt     .clear();
	pfHadron_eta    .clear();
	pfHadron_phi    .clear();
	pfHadron_pdgId  .clear();


	// get muons

	//  typedef View<reco::MuonCollection> MuonCollectionView;
	Handle<MuonCollection> muonCandidates;
	iEvent.getByToken(Muon_token,muonCandidates);


	muon_PF         .clear();
	muon_pt         .clear();
	muon_eta         .clear();
	muon_phi         .clear();
	muon_ptError    .clear();
	imuon_pt         .clear();
	imuon_eta         .clear();
	imuon_phi         .clear();
	imuon_ptError    .clear();  
	muon_SC          .clear();

	// get Jets
	Handle<reco::PFJetCollection> pfJets;
	iEvent.getByToken(PfJets_token,pfJets);

	pfJet_pt     .clear();
	pfJet_eta    .clear();
	pfJet_phi    .clear();
	pfJet_looseId.clear();
	pfJet_tightId.clear();
	pfJet_tlvId  .clear();
	pfJet_hpfl = -1;
	pfJet_hpft = -1;
	pfJet_hpfv = -1;

	// get METs
	Handle<reco::CaloMETCollection> caloMET;
	iEvent.getByToken(CaloMET_token, caloMET);

	if (isReco_token){
		Handle<reco::PFMETCollection> pfCaloMET;
		iEvent.getByToken(PFCaloMET_token, pfCaloMET);

		pfCaloMETPt = pfCaloMET->begin()->pt();
		pfCaloMETPhi = pfCaloMET->begin()->phi();
		pfCaloMETSumEt = pfCaloMET->begin()->sumEt();

		Handle<reco::PFClusterMETCollection> pfClusterMET;
		iEvent.getByToken(PFClusterMET_token, pfClusterMET);

		pfClusterMETPt = pfClusterMET->begin()->pt();
		pfClusterMETPhi = pfClusterMET->begin()->phi();
		pfClusterMETSumEt = pfClusterMET->begin()->sumEt();


		//get PFClusters
		Handle<reco::PFClusterCollection> pfClustersEcal;
		iEvent.getByToken(EcalPFClusters_token,pfClustersEcal);

		Handle<reco::PFClusterCollection> pfClustersHcal;
		iEvent.getByToken(HcalPFClusters_token,pfClustersHcal);

		Handle<reco::PFClusterCollection> pfClustersHBHE;
		iEvent.getByToken(HBHEPFClusters_token,pfClustersHBHE);

		Handle<reco::PFClusterCollection> pfClustersHO;
		iEvent.getByToken(HOPFClusters_token,pfClustersHO);

		Handle<reco::PFClusterCollection> pfClustersHF;
		iEvent.getByToken(HFPFClusters_token,pfClustersHF);

		pfClusterEcal_energy.clear();
		pfClusterEcal_eta.clear();
		pfClusterEcal_phi.clear();
		pfClusterEcal_time.clear();
		pfClusterEcal_status13.clear();
		pfClusterEcal_status14.clear();

		pfClusterHcal_energy.clear();
		pfClusterHcal_eta.clear();
		pfClusterHcal_phi.clear();
		pfClusterHcal_time.clear();

		pfClusterHBHE_energy.clear();
		pfClusterHBHE_eta.clear();
		pfClusterHBHE_phi.clear();
		pfClusterHBHE_time.clear();

		pfClusterHO_energy.clear();
		pfClusterHO_eta.clear();
		pfClusterHO_phi.clear();
		pfClusterHO_time.clear();

		pfClusterHF_energy.clear();
		pfClusterHF_eta.clear();
		pfClusterHF_phi.clear();
		pfClusterHF_time.clear();

		// Barrel
		edm::Handle< EcalRecHitCollection > ebRecHits_h_;
		iEvent.getByToken( RecHitsEB_token, ebRecHits_h_ );

		// Endcaps
		edm::Handle< EcalRecHitCollection > eeRecHits_h_;
		iEvent.getByToken( RecHitsEE_token, eeRecHits_h_ );


		// Preshower
		edm::Handle< EcalRecHitCollection > esRecHits_h_;
		iEvent.getByToken( RecHitsES_token, esRecHits_h_ );


		//ECAL clusters
		for( size_t ibc=0; ibc<pfClustersEcal->size(); ++ibc ) {
			reco::PFClusterRef bcRef( pfClustersEcal, ibc );
			pfClusterEcal_energy.push_back( bcRef->energy() );
			pfClusterEcal_time.push_back( bcRef->time() );
			pfClusterEcal_eta.push_back( bcRef->eta() );
			pfClusterEcal_phi.push_back( bcRef->phi() );

			// retrieve the id in the list of rechits and get the severity level
			bool status13 = false;
			bool status14 = false;
			vector<std::pair<DetId, float> > detId_v = bcRef->hitsAndFractions();
			unsigned nhits = detId_v.size();
			for ( size_t ihit=0; ihit<nhits; ihit++ ) {
				int sev = 0;
				if( (detId_v[ihit].first).subdetId() == EcalBarrel) {
					EBDetId id( (detId_v[ihit].first).rawId() ); 
					sev =  (Int_t) sevlv->severityLevel( id, *ebRecHits_h_);
					//sev =  (Int_t) sevlv->severityLevel( id, *(ebRecHits_h_.product()) );
				}
				else if( (detId_v[ihit].first).subdetId() == EcalEndcap) {
					EEDetId id( (detId_v[ihit].first).rawId() );
					sev =  (Int_t) sevlv->severityLevel( id, *eeRecHits_h_ );
					//sev =  (Int_t) sevlv->severityLevel( id, *(eeRecHits_h_.product()) );
				}
				else if( (detId_v[ihit].first).subdetId() == EcalPreshower) {
					ESDetId id( (detId_v[ihit].first).rawId() );
					sev =  (Int_t) sevlv->severityLevel( id, *esRecHits_h_ );
					//sev =  (Int_t) sevlv->severityLevel( id, *(esRecHits_h_.product()) );
				}
				if     (sev == 13) status13 = true;
				else if(sev == 14) status14 = true;
			}

			pfClusterEcal_status13.push_back( status13 );
			pfClusterEcal_status14.push_back( status14 );
		}

		//HCAL clusters
		for( size_t ibc=0; ibc<pfClustersHcal->size(); ++ibc ) {
			reco::PFClusterRef bcRef( pfClustersHcal, ibc );
			pfClusterHcal_energy.push_back( bcRef->energy() );
			pfClusterHcal_time.push_back( bcRef->time() );
			pfClusterHcal_eta.push_back( bcRef->eta() );
			pfClusterHcal_phi.push_back( bcRef->phi() );
		}

		//HBHE clusters
		for( size_t ibc=0; ibc<pfClustersHBHE->size(); ++ibc ) {
			reco::PFClusterRef bcRef( pfClustersHBHE, ibc );
			pfClusterHBHE_energy.push_back( bcRef->energy() );
			pfClusterHBHE_time.push_back( bcRef->time() );
			pfClusterHBHE_eta.push_back( bcRef->eta() );
			pfClusterHBHE_phi.push_back( bcRef->phi() );
		}

		//HO clusters
		for( size_t ibc=0; ibc<pfClustersHO->size(); ++ibc ) {
			reco::PFClusterRef bcRef( pfClustersHO, ibc );
			pfClusterHO_energy.push_back( bcRef->energy() );
			pfClusterHO_time.push_back( bcRef->time() );
			pfClusterHO_eta.push_back( bcRef->eta() );
			pfClusterHO_phi.push_back( bcRef->phi() );
		}

		//HF clusters
		for( size_t ibc=0; ibc<pfClustersHF->size(); ++ibc ) {
			reco::PFClusterRef bcRef( pfClustersHF, ibc );
			pfClusterHF_energy.push_back( bcRef->energy() );
			pfClusterHF_time.push_back( bcRef->time() );
			pfClusterHF_eta.push_back( bcRef->eta() );
			pfClusterHF_phi.push_back( bcRef->phi() );
		}


	}


	Handle<reco::PFMETCollection> pfMET;
	iEvent.getByToken(PFMET_token, pfMET);


	//get tracks
	Handle<reco::TrackCollection> tracks;
	iEvent.getByToken(Tracks_token,tracks);

	//get Ecal RecHits
	//not super sure it is useful, keep commeted for now
	//edm::ESHandle<EcalChannelStatus> chStatus_;
	//iSetup.get<EcalChannelStatusRcd>().get(chStatus_);

	track_ptError.clear();
	track_pt.clear();
	track_eta.clear();
	track_phi.clear();
	track_d0.clear();
	track_d0Error.clear();
	track_dz.clear();
	track_dzError.clear();
	//================================================================

	_muonEnergy.clear();
	_muonInnerTrackPt.clear();
	_muonInnerTrackPtError.clear();
	_muonBestTrackPt.clear();
	_muonBestTrackPtError.clear();
	_muonInnerTrackCharge.clear();
	_muonBestTrackCharge.clear();
	_muonGlobal.clear();
	_muonalgo.clear();
	_muonsegComp.clear();
	_muonHighpurity.clear();
	_muonPF.clear();
	_muonGlobalTrackPt.clear();
	_muonPickyTrackPt.clear();
	_muonTPFMSTrackPt.clear();
	_muonDYTTrackPt.clear();

	_muonGlobalTrackPtError.clear();
	_muonPickyTrackPtError.clear();
	_muonTPFMSTrackPtError.clear();
	_muonDYTTrackPtError.clear();

	_muonGlobalTrackCharge.clear();
	_muonPickyTrackCharge.clear();
	_muonTPFMSTrackCharge.clear();
	_muonDYTTrackCharge.clear();

	_muonTrackIsoSumPt.clear();
	_muonTunePMuonBestTrackCharge.clear();
	_muonTunePMuonBestTrackPt.clear();
	_muonTunePMuonBestTrackPtError.clear();

	_muonValidPixelHits.clear();
	_muonValidMuonHits.clear();
	_muonTrackerLayersHits.clear();
	_muonMatchedStations.clear();
	_muonDxy.clear();
	_muonDz.clear();

	_muonHighPtTrackerId.clear();
	_muonHighPtId.clear();

	_muonLoose.clear();
	_muonMedium.clear();
	_muonTight.clear();
	_muonSoft.clear();
	_muonTracker.clear();
	_muonMatchedNStations.clear();
	_muonStationMask.clear();
	_muonRPCLayers.clear();
	_muonBestValidMuonHits.clear();


	//pfCandidates
	for( size_t ibc=0; ibc<pfCandidates->size(); ++ibc ) {

		int pdgId = pfCandidates->at(ibc).pdgId();
		if(std::abs(pdgId) < 11 || std::abs(pdgId) > 16) continue;


		pfLepton_pt   .push_back( pfCandidates->at(ibc).pt()  );
		pfLepton_eta  .push_back( pfCandidates->at(ibc).eta() );
		pfLepton_phi  .push_back( pfCandidates->at(ibc).phi() );
		pfLepton_pdgId.push_back( pdgId                   );
		//pfLepton_d0   .push_back( d0                      );

	}

	//pfHadrons

	for( size_t ibc=0; ibc<pfCandidates->size(); ++ibc ) {
		int pdgId = pfCandidates->at(ibc).pdgId();
		if(std::abs(pdgId) != 211 ) continue;

		pfHadron_pt   .push_back( pfCandidates->at(ibc).pt() );
		pfHadron_eta  .push_back( pfCandidates->at(ibc).eta() );
		pfHadron_phi  .push_back( pfCandidates->at(ibc).phi() );
		pfHadron_pdgId.push_back( pdgId                   );
	}


	//Muons
	for( size_t ibc=0; ibc<muonCandidates->size(); ++ibc ) {
		//int pdgId = muonCandidates->at(ibc).pdgId();

		reco::TrackRef innerMuonTrack = muonCandidates->at(ibc).innerTrack();
		// reco::TrackRef globalMuonTrack = muon.globalTrack();

		// const reco::MuonCollection & muon = (*muons)[ibc];
		reco::TrackRef bestMuonTrack = muonCandidates->at(ibc).muonBestTrack();

		//    if( innerMuonTrack.isNull() ) continue;

		//   if( bestMuonTrack.isNull()  ) continue;

		muon_PF    .push_back(muonCandidates->at(ibc).isPFMuon());
		muon_pt    .push_back( bestMuonTrack.isNonnull() ? bestMuonTrack->pt() : -1  );
		muon_eta    .push_back( bestMuonTrack.isNonnull() ? bestMuonTrack->eta() :-9  );
		muon_phi    .push_back( bestMuonTrack.isNonnull() ? bestMuonTrack->phi() :-9  );
		muon_ptError    .push_back( bestMuonTrack.isNonnull() ?  bestMuonTrack->ptError() :-1  );
		imuon_pt    .push_back( innerMuonTrack.isNonnull() ? innerMuonTrack->pt() : -1  );
		imuon_eta    .push_back( innerMuonTrack.isNonnull() ? innerMuonTrack->eta() : -9  );
		imuon_phi    .push_back( innerMuonTrack.isNonnull() ? innerMuonTrack->phi() :-9  );

		imuon_ptError    .push_back( innerMuonTrack.isNonnull() ? innerMuonTrack->ptError() : -1  );
		muon_SC          .push_back(muon::segmentCompatibility(muonCandidates->at(ibc)));

		_muonEnergy.push_back(muonCandidates->at(ibc).energy());
		_muonInnerTrackPt.push_back(muonCandidates->at(ibc).innerTrack().isNonnull() ? muonCandidates->at(ibc).innerTrack()->pt() : -1);
		_muonBestTrackPt.push_back(muonCandidates->at(ibc).muonBestTrack().isNonnull() ? muonCandidates->at(ibc).muonBestTrack()->pt() : -1);
		_muonInnerTrackPtError.push_back(muonCandidates->at(ibc).innerTrack().isNonnull() ? muonCandidates->at(ibc).innerTrack()->ptError() : -1);
		_muonBestTrackPtError.push_back(muonCandidates->at(ibc).muonBestTrack().isNonnull() ? muonCandidates->at(ibc).muonBestTrack()->ptError() : -1);
		_muonInnerTrackCharge.push_back(muonCandidates->at(ibc).innerTrack().isNonnull() ? muonCandidates->at(ibc).innerTrack()->charge() : -9);
		_muonBestTrackCharge.push_back(muonCandidates->at(ibc).muonBestTrack().isNonnull() ? muonCandidates->at(ibc).muonBestTrack()->charge() : -9);
		_muonGlobal.push_back(muonCandidates->at(ibc).isGlobalMuon() );
		_muonalgo.push_back(muonCandidates->at(ibc).innerTrack().isNonnull() ? muonCandidates->at(ibc).innerTrack()->originalAlgo() : -1);
		_muonsegComp.push_back(muon::segmentCompatibility(muonCandidates->at(ibc)));
		_muonHighpurity.push_back(muonCandidates->at(ibc).innerTrack().isNonnull() ? muonCandidates->at(ibc).innerTrack()->quality(reco::TrackBase::highPurity): 0);
		_muonPF.push_back(muonCandidates->at(ibc).isPFMuon());

		_muonGlobalTrackPt.push_back(muonCandidates->at(ibc).globalTrack().isNonnull() ? muonCandidates->at(ibc).globalTrack()->pt() : -1);
		_muonGlobalTrackPtError.push_back(muonCandidates->at(ibc).globalTrack().isNonnull() ? muonCandidates->at(ibc).globalTrack()->ptError() : -1);
		_muonGlobalTrackCharge.push_back(muonCandidates->at(ibc).globalTrack().isNonnull() ? muonCandidates->at(ibc).globalTrack()->charge() : -9);
		_muonTPFMSTrackPt.push_back(muonCandidates->at(ibc).tpfmsTrack().isAvailable() ? muonCandidates->at(ibc).tpfmsTrack()->pt() : -1);
		_muonTPFMSTrackPtError.push_back(muonCandidates->at(ibc).tpfmsTrack().isAvailable() ? muonCandidates->at(ibc).tpfmsTrack()->ptError() : -1);
		_muonTPFMSTrackCharge.push_back(muonCandidates->at(ibc).tpfmsTrack().isAvailable() ? muonCandidates->at(ibc).tpfmsTrack()->charge() : -9);

		_muonPickyTrackPt.push_back(muonCandidates->at(ibc).pickyTrack().isAvailable() ? muonCandidates->at(ibc).pickyTrack()->pt() : -1);
		_muonPickyTrackPtError.push_back(muonCandidates->at(ibc).pickyTrack().isAvailable() ? muonCandidates->at(ibc).pickyTrack()->ptError() : -1);
		_muonPickyTrackCharge.push_back(muonCandidates->at(ibc).pickyTrack().isAvailable() ? muonCandidates->at(ibc).pickyTrack()->charge() : -9);
		_muonDYTTrackPt.push_back(muonCandidates->at(ibc).dytTrack().isAvailable() ? muonCandidates->at(ibc).dytTrack()->pt() : -1);
		_muonDYTTrackPtError.push_back(muonCandidates->at(ibc).dytTrack().isAvailable() ? muonCandidates->at(ibc).dytTrack()->ptError() : -1);
		_muonDYTTrackCharge.push_back(muonCandidates->at(ibc).dytTrack().isAvailable() ? muonCandidates->at(ibc).dytTrack()->charge() : -9);

		_muonTunePMuonBestTrackPt.push_back(muonCandidates->at(ibc).tunePMuonBestTrack().isAvailable() ? muonCandidates->at(ibc).tunePMuonBestTrack()->pt() : -1);
		_muonTunePMuonBestTrackPtError.push_back(muonCandidates->at(ibc).tunePMuonBestTrack().isAvailable() ? muonCandidates->at(ibc).tunePMuonBestTrack()->ptError() : -1);
		_muonTunePMuonBestTrackCharge.push_back(muonCandidates->at(ibc).tunePMuonBestTrack().isAvailable() ? muonCandidates->at(ibc).tunePMuonBestTrack()->charge() : -9);

		_muonTrackIsoSumPt.push_back(muonCandidates->at(ibc).isolationR03().sumPt);
		_muonValidPixelHits.push_back(muonCandidates->at(ibc).innerTrack().isNonnull() ? muonCandidates->at(ibc).innerTrack()->hitPattern().numberOfValidPixelHits() : -1);
		_muonValidMuonHits.push_back(muonCandidates->at(ibc).globalTrack().isNonnull() ? muonCandidates->at(ibc).globalTrack()->hitPattern().numberOfValidMuonHits() : -1);
		_muonMatchedStations.push_back(muonCandidates->at(ibc).numberOfMatchedStations());
		_muonTrackerLayersHits.push_back(muonCandidates->at(ibc).innerTrack().isNonnull() ? muonCandidates->at(ibc).innerTrack()->hitPattern().trackerLayersWithMeasurement() : -1);
		_muonDxy.push_back(muonCandidates->at(ibc).muonBestTrack().isNonnull() ? fabs(muonCandidates->at(ibc).muonBestTrack()->dxy(offlineVtx->begin()->position())) : -999);
		_muonDz.push_back(muonCandidates->at(ibc).muonBestTrack().isNonnull() ? fabs(muonCandidates->at(ibc).muonBestTrack()->dz(offlineVtx->begin()->position())) : -999);
		_muonHighPtId.push_back(int(muon::isHighPtMuon(muonCandidates->at(ibc), *offlineVtx->begin())));
		_muonHighPtTrackerId.push_back(int(muon::isTrackerHighPtMuon(muonCandidates->at(ibc), *offlineVtx->begin())));
		_muonLoose.push_back(int(muon::isLooseMuon(muonCandidates->at(ibc))));
		_muonMedium.push_back(int(muon::isMediumMuon(muonCandidates->at(ibc))));
		_muonTight.push_back(int(muon::isTightMuon(muonCandidates->at(ibc), *offlineVtx->begin())));
		_muonSoft.push_back(int(muon::isSoftMuon(muonCandidates->at(ibc), *offlineVtx->begin())));
		_muonTracker.push_back(muonCandidates->at(ibc).isTrackerMuon() );
		_muonMatchedNStations.push_back(muonCandidates->at(ibc).expectedNnumberOfMatchedStations());
		_muonStationMask.push_back(muonCandidates->at(ibc).stationMask());
		_muonRPCLayers.push_back(muonCandidates->at(ibc).numberOfMatchedRPCLayers());
		_muonBestValidMuonHits.push_back(muonCandidates->at(ibc).tunePMuonBestTrack().isAvailable() ? muonCandidates->at(ibc).tunePMuonBestTrack()->hitPattern().numberOfValidMuonHits(): -1);




	}

	//======================================================================
	//  for ( unsigned i=0; i < muons->size(); ++i ) { // loop over all muons                                                                              

	//const reco::Muon & muon = (*muons)[i];

	//    reco::TrackRef bestMuonTrack = muon.muonBestTrack();
	//muon_pt    .push_back( bestMuonTrack->pt()   ); 
	//}





	//================================================================
	//pfJets
	int maxl = -1;
	int maxt = -1;
	int maxv = -1;
	for( size_t ibc=0; ibc<pfJets->size(); ++ibc ) {

		bool looseId = false, tightId = false, tlvId = false;

		//looseJetID = (NHF<0.99 && NEMF<0.99 && NumConst>1 && MUF<0.8) && ((abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || abs(eta)>2.4) 
		looseId = ((pfJets->at(ibc).neutralHadronEnergyFraction() < 0.99 && pfJets->at(ibc).neutralEmEnergyFraction() < 0.99 
					&& pfJets->at(ibc).nConstituents() > 1
					&& pfJets->at(ibc).muonEnergyFraction() < 0.8)
				&& ((std::abs(pfJets->at(ibc).eta()) <= 2.4 && pfJets->at(ibc).chargedHadronEnergyFraction() > 0
						&& pfJets->at(ibc).chargedHadronMultiplicity() > 0
						&& pfJets->at(ibc).chargedEmEnergyFraction() < 0.9)
					|| std::abs(pfJets->at(ibc).eta()) > 2.4));

		//tightJetID = (NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8) && ((abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.90) || abs(eta)>2.4) 
		tightId = ((pfJets->at(ibc).neutralHadronEnergyFraction() < 0.90 && pfJets->at(ibc).neutralEmEnergyFraction() < 0.90 
					&& pfJets->at(ibc).nConstituents() > 1)
				&& ((std::abs(pfJets->at(ibc).eta()) <= 2.4 && pfJets->at(ibc).chargedHadronEnergyFraction() > 0
						&& pfJets->at(ibc).chargedHadronMultiplicity() > 0
						&& pfJets->at(ibc).chargedEmEnergyFraction() < 0.9)
					|| std::abs(pfJets->at(ibc).eta()) > 2.4));

		//tightLepVetoJetID = (NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8) && ((abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.90) || abs(eta)>2.4) 
		tlvId = ((pfJets->at(ibc).neutralHadronEnergyFraction() < 0.90 && pfJets->at(ibc).neutralEmEnergyFraction() < 0.90 
					&& pfJets->at(ibc).nConstituents() > 1
					&& pfJets->at(ibc).muonEnergyFraction() < 0.8)
				&& ((std::abs(pfJets->at(ibc).eta()) <= 2.4 && pfJets->at(ibc).chargedHadronEnergyFraction() > 0
						&& pfJets->at(ibc).chargedHadronMultiplicity() > 0
						&& pfJets->at(ibc).chargedEmEnergyFraction() < 0.9)
					|| std::abs(pfJets->at(ibc).eta()) > 2.4));


		//indices of highest pT jets failing loose and tight Id
		if((maxl == -1 || (maxl > 0 && pfJets->at(maxl).pt() < pfJets->at(ibc).pt())) && looseId == false) maxl = ibc;
		if((maxt == -1 || (maxt > 0 && pfJets->at(maxt).pt() < pfJets->at(ibc).pt())) && tightId == false) maxt = ibc;
		if((maxv == -1 || (maxv > 0 && pfJets->at(maxv).pt() < pfJets->at(ibc).pt())) && tlvId   == false) maxv = ibc;

		pfJet_pt     .push_back( pfJets->at(ibc).pt()   );
		pfJet_eta    .push_back( pfJets->at(ibc).eta()  );
		pfJet_phi    .push_back( pfJets->at(ibc).phi()  );
		pfJet_looseId.push_back( looseId );
		pfJet_tightId.push_back( tightId );
		pfJet_tlvId  .push_back( tlvId );
	}

	pfJet_hpfl = maxl;
	pfJet_hpft = maxt;
	pfJet_hpfv = maxv;  


	//METs =======================
	caloMETPt = caloMET->begin()->pt();
	caloMETPhi = caloMET->begin()->phi();
	caloMETSumEt = caloMET->begin()->sumEt();

	pfMETPt = pfMET->begin()->pt();
	pfMETPhi = pfMET->begin()->phi();
	pfMETSumEt = pfMET->begin()->sumEt();



	//tracks
	for( size_t ibc=0; ibc<tracks->size(); ++ibc ) {
		reco::TrackRef trkRef( tracks, ibc );
		track_ptError.push_back( trkRef->ptError() );
		track_pt.push_back( trkRef->pt() );
		track_eta.push_back( trkRef->eta() );
		track_phi.push_back( trkRef->phi() );
		track_d0.push_back( trkRef->d0() );
		track_d0Error.push_back( trkRef->d0Error() );
		track_dz.push_back( trkRef->dz() );
		track_dzError.push_back( trkRef->dzError() );
	}


	//tree filling ===========================
	s->Fill();

}

std::vector<std::string> METScanningNtupleMaker::split(const std::string &s, char delim = ':'){

	std::vector<std::string> elems;

	std::stringstream ss(s);
	std::string item;
	while (std::getline(ss, item, delim))
	{
		elems.push_back(item);
	}
	return elems;

}

DEFINE_FWK_MODULE(METScanningNtupleMaker);

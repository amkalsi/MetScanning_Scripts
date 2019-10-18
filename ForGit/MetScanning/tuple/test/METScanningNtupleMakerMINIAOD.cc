// Original Author:  Laurent Thomas
//         Created:  Fri, 26 Apr 2019 12:51:46 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TString.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include <fstream>
#include <iostream>

const int  N_METFilters=17;
enum METFilterIndex{
	idx_Flag_goodVertices,
	idx_Flag_globalTightHalo2016Filter,
	idx_Flag_globalSuperTightHalo2016Filter,
	idx_Flag_HBHENoiseFilter,
	idx_Flag_HBHENoiseIsoFilter,
	idx_Flag_EcalDeadCellTriggerPrimitiveFilter,
	idx_Flag_BadPFMuonFilter,
	idx_Flag_BadChargedCandidateFilter,
	idx_Flag_eeBadScFilter,
	idx_Flag_ecalBadCalibFilter,
	idx_Flag_ecalLaserCorrFilter,
	idx_Flag_EcalDeadCellBoundaryEnergyFilter,
	idx_PassecalBadCalibFilter_Update,
	idx_PassecalLaserCorrFilter_Update,
	idx_PassEcalDeadCellBoundaryEnergyFilter_Update,
	idx_PassBadChargedCandidateFilter_Update,
        idx_PassBadPFMuonFilter_Update

};

struct BadRuns {
	unsigned long Run_Bad; 
	unsigned long LS_Bad;
	unsigned long event_Bad;
};


//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.
using namespace edm;
using namespace std;
using namespace reco;
using namespace pat;


class METScanningNtupleMakerMINIAOD : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
	public:
		explicit METScanningNtupleMakerMINIAOD(const edm::ParameterSet&);
		~METScanningNtupleMakerMINIAOD();

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


	private:
		virtual void beginJob() override;
		virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
		virtual void endJob() override;
		virtual bool GetMETFilterDecision(const edm::Event& iEvent, edm::Handle<TriggerResults> METFilterResults, TString studiedfilter);
		virtual bool GetIdxFilterDecision(int it);
		virtual TString GetIdxFilterName(int it);
		virtual void InitandClearStuff();
		std::vector<std::string> split(const std::string &s, char delim);

		// ----------member data ---------------------------
		edm::EDGetTokenT<TriggerResults> metfilterspatToken_; 
		edm::EDGetTokenT<TriggerResults> metfiltersrecoToken_; 
		edm::EDGetTokenT<bool> ecalBadCalibFilterUpdateToken_;
		edm::EDGetTokenT<bool> ecalLaserCorrFilterUpdateToken_;  
		edm::EDGetTokenT<bool> ecalDeadCellBoundaryEnergyFilterUpdateToken_;
		edm::EDGetTokenT<bool> badChargedCandidateFilterUpdateToken_;
                edm::EDGetTokenT<bool> badPFMuonFilterUpdateToken_;

		edm::EDGetTokenT<std::vector<Vertex> > verticesToken_; 

		edm::EDGetTokenT<std::vector< pat::Jet> > jetToken_;
		edm::EDGetTokenT<pat::PackedCandidateCollection> pfcandsToken_;

		edm::EDGetTokenT<std::vector< pat::MET> > metToken_;
		edm::EDGetTokenT<std::vector< pat::MET> > puppimetToken_;
		edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
		edm::EDGetTokenT<pat::MuonCollection> muonToken_;


		//  edm::Service<TFileService> fs;


		//Some histos to be saved for simple checks 
		TH1F *h_PFMet, *h_PuppiMet, *h_nvtx, *h_leadjetpt;
		TH1F *h_PFMet_num[N_METFilters], *h_PuppiMet_num[N_METFilters], *h_nvtx_num[N_METFilters], *h_leadjetpt_num[N_METFilters];
		TH2F *h_jet200etavsphi_fail[N_METFilters];
		//The output TTree
		TTree* outputTree;

		//Variables associated to leaves of the TTree

		unsigned long _eventNb;
		unsigned long _runNb;
		unsigned long _lumiBlock;
		unsigned long _bx;

		//Nb of primary vertices
		int _n_PV;
		bool run_found;

		//MINIAOD original MET filters decisions
		bool Flag_goodVertices;
		bool Flag_globalTightHalo2016Filter;
		bool Flag_globalSuperTightHalo2016Filter;
		bool Flag_HBHENoiseFilter;
		bool Flag_HBHENoiseIsoFilter;
		bool Flag_EcalDeadCellTriggerPrimitiveFilter;
		bool Flag_BadPFMuonFilter;
		bool Flag_BadChargedCandidateFilter;
		bool Flag_eeBadScFilter;
		bool Flag_ecalBadCalibFilter;
		bool Flag_ecalLaserCorrFilter; 
		bool Flag_EcalDeadCellBoundaryEnergyFilter;

		//Decision obtained rerunning the filters on top of MINIAOD
		bool PassecalBadCalibFilter_Update;
		bool PassecalLaserCorrFilter_Update;  
		bool PassEcalDeadCellBoundaryEnergyFilter_Update;
		bool PassBadChargedCandidateFilter_Update;
		bool PassBadPFMuonFilter_Update;
		//Jets 
		vector<double> _jetEta;
		vector<double>  _jetPhi;
		vector<double>  _jetPt;
		vector<double>  _jetRawPt;
		vector<double>  _jet_CHEF;
		vector<double>  _jet_NHEF;
		vector<double>  _jet_NEEF;
		vector<double>  _jet_CEEF;
		vector<double>  _jet_MUEF;
		vector <int>  _jet_CHM;
		vector <int>  _jet_NHM;
		vector <int>  _jet_PHM;
		vector <int>  _jet_NM;

		// Electron  candidate

		vector<double>_elePt;
		vector<double>_eleEta;
		vector<double>_elePhi;
		vector<double>_eleEnergy;
		vector<double>_eleEt;

		vector <int>  gsf_VID_cutBasedElectronID_Fall17_94X_V1_loose;
		vector <int>  gsf_VID_cutBasedElectronID_Fall17_94X_V1_medium;
		vector <int>  gsf_VID_cutBasedElectronID_Fall17_94X_V1_tight;
		vector <int>  gsf_VID_cutBasedElectronID_Fall17_94X_V1_veto;

		//PF candidates
		vector <double> _PFcand_pt;
		vector <double> _PFcand_eta;
		vector <double> _PFcand_phi;
		vector <int> _PFcand_pdgId;
		vector <int> _PFcand_fromPV;

		//MET
		double _met;
		double _met_phi;
		double _puppimet;
		double _puppimet_phi;

		// string
		string file_INPUT;
		std::ifstream badrunFile;
		std::vector<BadRuns> BadRunsList;

		vector<double> _muonPt;
		vector<double> _muonEta;
		vector<double> _muonPhi;
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
		vector<double> _muonInnerDxy;
		vector<double> _muonInnerDz;

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


		vector<double>   _muonGlobTrackNChiS;
		vector<double> _muonCombQualityChiSLP;
		vector<double> _muonCombQualitytrkKink;
		vector<double> _muonInnerTrackVFrac;

		vector<int> _muonPixelLayers;
		vector<int> _muonGoodMuon;



};


//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
METScanningNtupleMakerMINIAOD::METScanningNtupleMakerMINIAOD(const edm::ParameterSet& iConfig)
	:
		metfilterspatToken_(consumes<TriggerResults>(iConfig.getParameter<edm::InputTag>("METFiltersPAT"))),
		metfiltersrecoToken_(consumes<TriggerResults>(iConfig.getParameter<edm::InputTag>("METFiltersRECO"))),
		ecalBadCalibFilterUpdateToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("ECALBadCalibFilterUpdate"))),
		ecalLaserCorrFilterUpdateToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("ECALLaserCorrFilterUpdate"))),
		ecalDeadCellBoundaryEnergyFilterUpdateToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("ECALDeadCellBoundaryEnergyFilterUpdate"))),
		badChargedCandidateFilterUpdateToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("BadChargedCandidateFilterUpdate"))),
                badPFMuonFilterUpdateToken_(consumes<bool>(iConfig.getParameter<edm::InputTag>("BadPFMuonFilterUpdate"))),
		verticesToken_(consumes<std::vector<Vertex> > (iConfig.getParameter<edm::InputTag>("Vertices"))),
		jetToken_(consumes< std::vector< pat::Jet> >(iConfig.getParameter<edm::InputTag>("Jets"))),
		pfcandsToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("PFCandCollection"))),
		metToken_(consumes<std::vector<pat::MET> > (iConfig.getParameter<edm::InputTag>("PFMet"))),
		puppimetToken_(consumes<std::vector<pat::MET> > (iConfig.getParameter<edm::InputTag>("PuppiMet"))),
		electronToken_(consumes<pat::ElectronCollection> (iConfig.getParameter<edm::InputTag>("electronCollection"))),
		muonToken_(consumes<pat::MuonCollection> (iConfig.getParameter<edm::InputTag>("muonCollection"))),
		file_INPUT(iConfig.getParameter<string>("file_INPUT"))
{
	//now do what ever initialization is needed
	edm::Service<TFileService> fs; 
	h_nvtx  = fs->make<TH1F>("h_nvtx" , "Number of reco vertices (MET>200);N_{vtx};Events"  ,    100, 0., 100.);
	h_PFMet  = fs->make<TH1F>("h_PFMet" , "Type 1 PFMET (GeV);Type 1 PFMET (GeV);Events"  ,    1000, 0., 5000.);
	h_PuppiMet  = fs->make<TH1F>("h_PuppiMet" , "PUPPI MET (GeV);PUPPI MET (GeV);Events"  ,    1000, 0., 5000.);
	h_leadjetpt  = fs->make<TH1F>("h_leadjetpt" , "Leading jet p_T (GeV);p_{T} (leading jet) (GeV) (MET<100);Events"  ,    1000, 0., 5000.);

	for(int i =0; i< N_METFilters;i++){
		TString filtername = GetIdxFilterName(i);
		h_nvtx_num[i]  = fs->make<TH1F>("h_nvtx_" +filtername , "Number of reco vertices, events passing "+filtername+";N_{vtx};Events", 100, 0., 100.);
		h_PFMet_num[i]  = fs->make<TH1F>("h_PFMet_" +filtername , "Type 1 PFMET (GeV), events passing "+filtername+";Type 1 PFMET (GeV);Events"  ,    1000, 0., 5000.);
		h_PuppiMet_num[i]  = fs->make<TH1F>("h_PuppiMet_" +filtername , "PUPPI MET (GeV), events passing "+filtername+";PUPPI MET (GeV);Events"  ,    1000, 0., 5000.);
		h_leadjetpt_num[i]  = fs->make<TH1F>("h_leadjetpt_" +filtername , "Leading jet p_T (GeV), events passing "+filtername+";p_{T} (leading jet) (GeV) (MET<100);Events"  ,    1000, 0., 5000.);
		h_jet200etavsphi_fail[i] = fs->make<TH2F>("h_jet200etavsphi_fail_" +filtername , "Jet (pt>200) occupancy map, events failing "+filtername+";#eta(jet);#phi(jet);Events"  ,    200, -5,5,200,-3.1416,3.1416);
	}
	outputTree = fs->make<TTree>("tree","tree");

}


METScanningNtupleMakerMINIAOD::~METScanningNtupleMakerMINIAOD()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
	void
METScanningNtupleMakerMINIAOD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

	InitandClearStuff();


	_runNb = iEvent.id().run();
	_eventNb = iEvent.id().event();
	_lumiBlock = iEvent.luminosityBlock();
	_bx=iEvent.bunchCrossing();

//	bool run_found;
	run_found=false;
	for (const BadRuns & m: BadRunsList) 
	{
		if(_runNb == m.Run_Bad && _lumiBlock == m.LS_Bad && _eventNb == m.event_Bad ) { run_found=true;}
	}
//	std::cout<<run_found;
//	if(!run_found) return;
	//Vertices
	edm::Handle<std::vector<Vertex> > theVertices;
	iEvent.getByToken(verticesToken_,theVertices) ;
	_n_PV = theVertices->size();


	//MET filters are stored in TriggerResults::RECO or TriggerResults::PAT . Should take the latter if it exists
	edm::Handle<TriggerResults> METFilterResults;
	iEvent.getByToken(metfilterspatToken_, METFilterResults);
	if(!(METFilterResults.isValid())) iEvent.getByToken(metfiltersrecoToken_, METFilterResults);

	//std::cout<<"outside worl"<<std::endl;

	const edm::TriggerNames &names = iEvent.triggerNames(*METFilterResults);
	for(unsigned int i =0 ;  i < METFilterResults->size() ; i++) {
		//std::cout<<"name:"<< names.triggerName(i)<< std::endl;

	}
	//std::cout<<"outside worl11"<<std::endl;

	Flag_goodVertices= GetMETFilterDecision(iEvent,METFilterResults,"Flag_goodVertices");
	Flag_globalTightHalo2016Filter= GetMETFilterDecision(iEvent,METFilterResults,"Flag_globalTightHalo2016Filter");
	Flag_globalSuperTightHalo2016Filter= GetMETFilterDecision(iEvent,METFilterResults,"Flag_globalSuperTightHalo2016Filter");
	Flag_HBHENoiseFilter= GetMETFilterDecision(iEvent,METFilterResults,"Flag_HBHENoiseFilter");
	Flag_HBHENoiseIsoFilter= GetMETFilterDecision(iEvent,METFilterResults,"Flag_HBHENoiseIsoFilter");
	Flag_EcalDeadCellTriggerPrimitiveFilter= GetMETFilterDecision(iEvent,METFilterResults,"Flag_EcalDeadCellTriggerPrimitiveFilter");
	Flag_BadPFMuonFilter= GetMETFilterDecision(iEvent,METFilterResults,"Flag_BadPFMuonFilter");
	Flag_BadChargedCandidateFilter= GetMETFilterDecision(iEvent,METFilterResults,"Flag_BadChargedCandidateFilter");
	Flag_eeBadScFilter= GetMETFilterDecision(iEvent,METFilterResults,"Flag_eeBadScFilter");
	Flag_ecalBadCalibFilter= GetMETFilterDecision(iEvent,METFilterResults,"Flag_ecalBadCalibFilter");
	Flag_EcalDeadCellBoundaryEnergyFilter= GetMETFilterDecision(iEvent,METFilterResults,"Flag_EcalDeadCellBoundaryEnergyFilter");
	Flag_ecalLaserCorrFilter= GetMETFilterDecision(iEvent,METFilterResults,"Flag_ecalLaserCorrFilter");


	//Now accessing the decisions of some filters that we reran on top of MINIAOD
	edm::Handle<bool> handle_PassecalBadCalibFilter_Update ;
	iEvent.getByToken(ecalBadCalibFilterUpdateToken_,handle_PassecalBadCalibFilter_Update);
	if(handle_PassecalBadCalibFilter_Update.isValid()) PassecalBadCalibFilter_Update =  (*handle_PassecalBadCalibFilter_Update );
	else PassecalBadCalibFilter_Update = true;

	edm::Handle<bool> handle_PassecalLaserCorrFilter_Update ;
	iEvent.getByToken(ecalLaserCorrFilterUpdateToken_,handle_PassecalLaserCorrFilter_Update);
	if(handle_PassecalLaserCorrFilter_Update.isValid())PassecalLaserCorrFilter_Update =  (*handle_PassecalLaserCorrFilter_Update );
	else PassecalLaserCorrFilter_Update = true;

	edm::Handle<bool> handle_PassEcalDeadCellBoundaryEnergyFilter_Update;
	iEvent.getByToken(ecalDeadCellBoundaryEnergyFilterUpdateToken_,handle_PassEcalDeadCellBoundaryEnergyFilter_Update);
	if(handle_PassEcalDeadCellBoundaryEnergyFilter_Update.isValid())PassEcalDeadCellBoundaryEnergyFilter_Update =  (*handle_PassEcalDeadCellBoundaryEnergyFilter_Update );
	else{  std::cout <<"handle_PassEcalDeadCellBoundaryEnergyFilter_Update.isValid =false" <<endl; PassEcalDeadCellBoundaryEnergyFilter_Update = true;}

	edm::Handle<bool> handle_PassBadChargedCandidateFilter_Update;
	iEvent.getByToken(badChargedCandidateFilterUpdateToken_,handle_PassBadChargedCandidateFilter_Update);
	if(handle_PassBadChargedCandidateFilter_Update.isValid())PassBadChargedCandidateFilter_Update =  (*handle_PassBadChargedCandidateFilter_Update );
	else{  std::cout <<"handle_PassBadChargedCandidateFilter_Update.isValid =false" <<endl; PassBadChargedCandidateFilter_Update = true;}
 
        edm::Handle<bool> handle_PassBadPFMuonFilterUpdate;
        iEvent.getByToken(badPFMuonFilterUpdateToken_,handle_PassBadPFMuonFilterUpdate);
        if(handle_PassBadPFMuonFilterUpdate.isValid())PassBadPFMuonFilter_Update =  (*handle_PassBadPFMuonFilterUpdate );
        else{  std::cout <<"handle_PassBadPFMuonFilter_Update.isValid =false" <<endl; PassBadPFMuonFilter_Update = true;}


	//Jets
	edm::Handle< std::vector< pat::Jet> > theJets;
	iEvent.getByToken(jetToken_,theJets );
	// muons
	edm::Handle<pat::MuonCollection>     muonCollection_;
	iEvent.getByToken( muonToken_, muonCollection_);
	for( pat::MuonCollection::const_iterator imuon = (*muonCollection_).begin() ; imuon != (*muonCollection_).end() ; imuon++) {

		_muonPt.push_back(imuon->pt());
		_muonEta.push_back(imuon->eta());
		_muonPhi.push_back(imuon->phi());
		_muonEnergy.push_back(imuon->energy());
		/*		std::cout<<"imuon->innerTrack().isNonnull():"<< imuon->innerTrack().isNonnull()<< std::endl;
				std::cout<<"imuon->muonBestTrack().isNonnull():"<< imuon->muonBestTrack().isNonnull()<< std::endl;
		//std::cout<<"imuon->innerTrack()->originalAlgo():"<<imuon->innerTrack()->originalAlgo()<< 
		std::cout<<"imuon->globalTrack().isNonnull():"<< imuon->globalTrack().isNonnull()<< std::endl;
		std::cout<<"imuon->tpfmsTrack().isNonnull():"<< imuon->tpfmsTrack().isNonnull()<< std::endl;
		std::cout<<"imuon->pickyTrack().isNonnull():"<<imuon->pickyTrack().isNonnull()<< std::endl;
		std::cout<<"imuon->dytTrack().isNonnull():"<<imuon->dytTrack().isNonnull()<<std::endl;
		std::cout<<"imuon->tunePMuonBestTrack().isNonnull():"<<imuon->tunePMuonBestTrack().isNonnull()<< std::endl;
		*/
		_muonInnerTrackPt.push_back(imuon->innerTrack().isNonnull() ? imuon->innerTrack()->pt() : -1);
		_muonBestTrackPt.push_back(imuon->muonBestTrack().isNonnull() ? imuon->muonBestTrack()->pt() : -1);
		_muonInnerTrackPtError.push_back(imuon->innerTrack().isNonnull() ? imuon->innerTrack()->ptError() : -1);
		_muonBestTrackPtError.push_back(imuon->muonBestTrack().isNonnull() ? imuon->muonBestTrack()->ptError() : -1);
		//std::cout<<"debug 1"<< std::endl;

		_muonInnerTrackCharge.push_back(imuon->innerTrack().isNonnull() ? imuon->innerTrack()->charge() : -9);
		_muonBestTrackCharge.push_back(imuon->muonBestTrack().isNonnull() ? imuon->muonBestTrack()->charge() : -9);
		//std::cout<<"debug 2"<< std::endl;

		_muonGlobal.push_back(imuon->isGlobalMuon() );
		_muonalgo.push_back(imuon->innerTrack().isNonnull() ? imuon->innerTrack()->originalAlgo() : -1);
		_muonsegComp.push_back(muon::segmentCompatibility(*imuon));
		_muonHighpurity.push_back(imuon->innerTrack().isNonnull() ? imuon->innerTrack()->quality(reco::TrackBase::highPurity): 0);
		_muonPF.push_back(imuon->isPFMuon());
		//std::cout<<"debug 3"<< std::endl;

		_muonGlobalTrackPt.push_back(imuon->globalTrack().isNonnull() ? imuon->globalTrack()->pt() : -1);
		_muonGlobalTrackPtError.push_back(imuon->globalTrack().isNonnull() ? imuon->globalTrack()->ptError() : -1);
		_muonGlobalTrackCharge.push_back(imuon->globalTrack().isNonnull() ? imuon->globalTrack()->charge() : -9);
		//std::cout<<"debug 4"<< std::endl;

		//std::cout<<" imuon->tpfmsTrack()->pt() :"<<  imuon->tpfmsTrack()->pt()<< std::endl; 
		_muonTPFMSTrackPt.push_back(imuon->tpfmsTrack().isAvailable() ? imuon->tpfmsTrack()->pt() : -1);
		_muonTPFMSTrackPtError.push_back(imuon->tpfmsTrack().isAvailable() ? imuon->tpfmsTrack()->ptError() : -1);
		_muonTPFMSTrackCharge.push_back(imuon->tpfmsTrack().isAvailable() ? imuon->tpfmsTrack()->charge() : -9);
		//std::cout<<"debug 5"<< std::endl;


		_muonPickyTrackPt.push_back(imuon->pickyTrack().isAvailable() ? imuon->pickyTrack()->pt() : -1);
		_muonPickyTrackPtError.push_back(imuon->pickyTrack().isAvailable() ? imuon->pickyTrack()->ptError() : -1);
		_muonPickyTrackCharge.push_back(imuon->pickyTrack().isAvailable() ? imuon->pickyTrack()->charge() : -9);

		//std::cout<<"debug 6"<< std::endl;

		_muonDYTTrackPt.push_back(imuon->dytTrack().isAvailable() ? imuon->dytTrack()->pt() : -1);
		_muonDYTTrackPtError.push_back(imuon->dytTrack().isAvailable() ? imuon->dytTrack()->ptError() : -1);
		_muonDYTTrackCharge.push_back(imuon->dytTrack().isAvailable() ? imuon->dytTrack()->charge() : -9);

		//std::cout<<"debug 7"<< std::endl;

		_muonTunePMuonBestTrackPt.push_back(imuon->tunePMuonBestTrack().isAvailable() ? imuon->tunePMuonBestTrack()->pt() : -1);
		_muonTunePMuonBestTrackPtError.push_back(imuon->tunePMuonBestTrack().isAvailable() ? imuon->tunePMuonBestTrack()->ptError() : -1);
		_muonTunePMuonBestTrackCharge.push_back(imuon->tunePMuonBestTrack().isAvailable() ? imuon->tunePMuonBestTrack()->charge() : -9);
		//std::cout<<"debug 8"<< std::endl;

		_muonTrackIsoSumPt.push_back(imuon->isolationR03().sumPt);
		_muonValidPixelHits.push_back(imuon->innerTrack().isNonnull() ? imuon->innerTrack()->hitPattern().numberOfValidPixelHits() : -1);
		_muonValidMuonHits.push_back(imuon->globalTrack().isNonnull() ? imuon->globalTrack()->hitPattern().numberOfValidMuonHits() : -1);
		//std::cout<<"debug 9"<< std::endl;

		_muonMatchedStations.push_back(imuon->numberOfMatchedStations());
		_muonTrackerLayersHits.push_back(imuon->innerTrack().isNonnull() ? imuon->innerTrack()->hitPattern().trackerLayersWithMeasurement() : -1);
		_muonDxy.push_back(imuon->muonBestTrack().isNonnull() ? fabs(imuon->muonBestTrack()->dxy(theVertices->begin()->position())) : -999);
		_muonDz.push_back(imuon->muonBestTrack().isNonnull() ? fabs(imuon->muonBestTrack()->dz(theVertices->begin()->position())) : -999);
		_muonInnerDxy.push_back(imuon->innerTrack().isNonnull() ? fabs(imuon->innerTrack()->dxy(theVertices->begin()->position())) : -999);
		_muonInnerDz.push_back(imuon->innerTrack().isNonnull() ? fabs(imuon->innerTrack()->dz(theVertices->begin()->position())) : -999);


		//std::cout<<"debug 10"<< std::endl;
		_muonHighPtId.push_back(int(muon::isHighPtMuon(*imuon, *theVertices->begin())));
		_muonHighPtTrackerId.push_back(int(muon::isTrackerHighPtMuon(*imuon, *theVertices->begin())));
		_muonLoose.push_back(int(muon::isLooseMuon(*imuon)));
		_muonMedium.push_back(int(muon::isMediumMuon(*imuon)));
		_muonTight.push_back(int(muon::isTightMuon(*imuon, *theVertices->begin())));
		_muonSoft.push_back(int(muon::isSoftMuon(*imuon, *theVertices->begin())));
		_muonTracker.push_back(imuon->isTrackerMuon() );
		_muonMatchedNStations.push_back(imuon->expectedNnumberOfMatchedStations());
		_muonStationMask.push_back(imuon->stationMask());
		_muonRPCLayers.push_back(imuon->numberOfMatchedRPCLayers());
		_muonBestValidMuonHits.push_back(imuon->tunePMuonBestTrack().isAvailable() ? imuon->tunePMuonBestTrack()->hitPattern().numberOfValidMuonHits(): -1);

		_muonGlobTrackNChiS.push_back(imuon->globalTrack().isNonnull() ? imuon->globalTrack()->normalizedChi2(): -999);
		_muonCombQualityChiSLP.push_back(imuon->combinedQuality().chi2LocalPosition);
		_muonCombQualitytrkKink.push_back(imuon->combinedQuality().trkKink);
		_muonInnerTrackVFrac.push_back(imuon->innerTrack().isNonnull() ? imuon->innerTrack()->validFraction(): -999) ;

		_muonPixelLayers.push_back(imuon->innerTrack().isNonnull() ? imuon->innerTrack()->hitPattern().pixelLayersWithMeasurement() : -1) ;
		_muonGoodMuon.push_back( int(muon::isGoodMuon(*imuon, muon::TMOneStationTight)) );


	}

	// Electrons
	//
	edm::Handle<pat::ElectronCollection>     electronCollection_ ;
	iEvent.getByToken( electronToken_ , electronCollection_) ;


	for( pat::ElectronCollection::const_iterator gsfiter = (*electronCollection_).begin() ; gsfiter != (*electronCollection_).end() ; gsfiter++ ) {

		_elePt.push_back(gsfiter->pt() );
		_eleEta.push_back(gsfiter->eta() );
		_elePhi.push_back(gsfiter->phi() );
		_eleEnergy.push_back(gsfiter->energy() );
		_eleEt.push_back( gsfiter->et() ) ;

		gsf_VID_cutBasedElectronID_Fall17_94X_V1_loose.push_back( int(gsfiter->electronID("cutBasedElectronID-Fall17-94X-V1-loose")));
		gsf_VID_cutBasedElectronID_Fall17_94X_V1_medium.push_back( int(gsfiter->electronID("cutBasedElectronID-Fall17-94X-V1-medium")));
		gsf_VID_cutBasedElectronID_Fall17_94X_V1_tight.push_back( int(gsfiter->electronID("cutBasedElectronID-Fall17-94X-V1-tight")));
		gsf_VID_cutBasedElectronID_Fall17_94X_V1_veto.push_back( int(gsfiter->electronID("cutBasedElectronID-Fall17-94X-V1-veto")));

	}


	double leadjetpt (0.);
	for( std::vector<pat::Jet>::const_iterator jet = (*theJets).begin(); jet != (*theJets).end(); jet++ ) {
		if((&*jet)->pt() >leadjetpt) leadjetpt = (&*jet)->pt();
		if((&*jet)->pt()<200) continue;
		_jetEta.push_back((&*jet)->eta());
		_jetPhi.push_back((&*jet)->phi());
		_jetPt.push_back((&*jet)->pt());
		_jetRawPt.push_back( (&*jet)->correctedP4("Uncorrected").Pt() );
		_jet_CHEF.push_back((&*jet)->chargedHadronEnergyFraction());
		_jet_NHEF.push_back((&*jet)->neutralHadronEnergyFraction() );
		_jet_NEEF.push_back((&*jet)->neutralEmEnergyFraction() );
		_jet_CEEF.push_back((&*jet)->chargedEmEnergyFraction() );
		_jet_MUEF.push_back((&*jet)->muonEnergyFraction() );
		_jet_CHM.push_back((&*jet)->chargedMultiplicity());
		_jet_NHM.push_back((&*jet)->neutralHadronMultiplicity());
		_jet_PHM.push_back((&*jet)->photonMultiplicity());
		_jet_NM.push_back((&*jet)->neutralMultiplicity());

		for(int i = 0; i< N_METFilters ;i++){
			if(!GetIdxFilterDecision(i) ){
				h_jet200etavsphi_fail[i]->Fill((&*jet)->eta(),(&*jet)->phi());
			}
		}
	}



	//Type 1 PFMET
	edm::Handle< vector<pat::MET> > ThePFMET;
	iEvent.getByToken(metToken_, ThePFMET);
	const vector<pat::MET> *pfmetcol = ThePFMET.product();
	const pat::MET *pfmet;
	pfmet = &(pfmetcol->front());
	_met = pfmet->pt();
	_met_phi = pfmet->phi();

	//PUPPI MET
	edm::Handle< vector<pat::MET> > ThePUPPIMET;
	iEvent.getByToken(puppimetToken_, ThePUPPIMET);
	const vector<pat::MET> *puppimetcol = ThePUPPIMET.product();
	const pat::MET *puppimet;
	puppimet = &(puppimetcol->front());
	_puppimet = puppimet->pt();
	_puppimet_phi = puppimet->phi();

	//PF candidates
	edm::Handle<pat::PackedCandidateCollection> pfcands;
	iEvent.getByToken(pfcandsToken_ ,pfcands);
	for(pat::PackedCandidateCollection::const_reverse_iterator p = pfcands->rbegin() ; p != pfcands->rend() ; p++ ) {
		if(p->pt()<200)continue;
		_PFcand_pt.push_back(p->pt());
		_PFcand_eta.push_back(p->eta());
		_PFcand_phi.push_back(p->phi());
		_PFcand_pdgId.push_back(p->pdgId());
		_PFcand_fromPV.push_back(p->fromPV(0));//See https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2017#Packed_ParticleFlow_Candidates
	}



	//Filling trees and histos   
	outputTree->Fill();
	//For effcy vs MET: 
	h_PFMet->Fill(_met);
	h_PuppiMet->Fill(_puppimet);
	for(int i = 0; i< N_METFilters ;i++){
		if( GetIdxFilterDecision(i) ){
			h_PFMet_num[i]->Fill(_met);
			h_PuppiMet_num[i]->Fill(_puppimet);
		}
	}
	//For accept rate vs leading jet (events with low MET)
	//This is to measure the mistag rate
	if(_met<100) {
		h_leadjetpt->Fill(leadjetpt);
		for(int i = 0; i< N_METFilters ;i++){
			if( GetIdxFilterDecision(i) )    h_leadjetpt_num[i]->Fill(leadjetpt);
		}
	}
	//For accept rate vs Nvtx 
	//This is to measure the efficiency 
	if(_met>200){
		h_nvtx->Fill(_n_PV);
		for(int i = 0; i< N_METFilters ;i++){
			if( GetIdxFilterDecision(i) ) h_nvtx_num[i]->Fill(_n_PV);
		}
	}


}


// ------------ method called once each job just before starting event loop  ------------
	void
METScanningNtupleMakerMINIAOD::beginJob()
{


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


	outputTree->Branch("_eventNb",   &_eventNb,   "_eventNb/l");
	outputTree->Branch("_runNb",     &_runNb,     "_runNb/l");
	outputTree->Branch("_lumiBlock", &_lumiBlock, "_lumiBlock/l");
	outputTree->Branch("_bx", &_bx, "_bx/l");
	outputTree->Branch("_n_PV", &_n_PV, "_n_PV/I");

	outputTree->Branch("Flag_goodVertices",&Flag_goodVertices,"Flag_goodVertices/O");
	outputTree->Branch("Flag_globalTightHalo2016Filter",&Flag_globalTightHalo2016Filter,"Flag_globalTightHalo2016Filter/O");
	outputTree->Branch("Flag_globalSuperTightHalo2016Filter",&Flag_globalSuperTightHalo2016Filter,"Flag_globalSuperTightHalo2016Filter/O");
	outputTree->Branch("Flag_HBHENoiseFilter",&Flag_HBHENoiseFilter,"Flag_HBHENoiseFilter/O");
	outputTree->Branch("Flag_HBHENoiseIsoFilter",&Flag_HBHENoiseIsoFilter,"Flag_HBHENoiseIsoFilter/O");
	outputTree->Branch("Flag_EcalDeadCellTriggerPrimitiveFilter",&Flag_EcalDeadCellTriggerPrimitiveFilter,"Flag_EcalDeadCellTriggerPrimitiveFilter/O");
	outputTree->Branch("Flag_BadPFMuonFilter",&Flag_BadPFMuonFilter,"Flag_BadPFMuonFilter/O");
	outputTree->Branch("Flag_BadChargedCandidateFilter",&Flag_BadChargedCandidateFilter,"Flag_BadChargedCandidateFilter/O");
	outputTree->Branch("Flag_eeBadScFilter",&Flag_eeBadScFilter,"Flag_eeBadScFilter/O");
	outputTree->Branch("Flag_ecalBadCalibFilter",&Flag_ecalBadCalibFilter,"Flag_ecalBadCalibFilter/O");
	outputTree->Branch("Flag_ecalLaserCorrFilter",&Flag_ecalLaserCorrFilter,"Flag_ecalLaserCorrFilter/O");
	outputTree->Branch("Flag_EcalDeadCellBoundaryEnergyFilter",&Flag_EcalDeadCellBoundaryEnergyFilter,"Flag_EcalDeadCellBoundaryEnergyFilter/O");

	outputTree->Branch("PassecalBadCalibFilter_Update",&PassecalBadCalibFilter_Update,"PassecalBadCalibFilter_Update/O");
	outputTree->Branch("PassecalLaserCorrFilter_Update",&PassecalLaserCorrFilter_Update,"PassecalLaserCorrFilter_Update/O");
	outputTree->Branch("PassEcalDeadCellBoundaryEnergyFilter_Update",&PassEcalDeadCellBoundaryEnergyFilter_Update,"PassEcalDeadCellBoundaryEnergyFilter_Update/O");
	outputTree->Branch("PassBadChargedCandidateFilter_Update",&PassBadChargedCandidateFilter_Update,"PassBadChargedCandidateFilter_Update/O");
	outputTree->Branch("PassBadPFMuonFilter_Update", &PassBadPFMuonFilter_Update, "PassBadPFMuonFilter_Update/O");

	outputTree->Branch("_jetEta",&_jetEta);
	outputTree->Branch("_jetPhi",&_jetPhi);
	outputTree->Branch("_jetPt",&_jetPt);
	outputTree->Branch("_jetRawPt",&_jetRawPt);
	outputTree->Branch("_jet_CHEF",&_jet_CHEF);
	outputTree->Branch("_jet_NHEF",&_jet_NHEF);
	outputTree->Branch("_jet_NEEF",&_jet_NEEF);
	outputTree->Branch("_jet_CEEF",&_jet_CEEF);
	outputTree->Branch("_jet_MUEF",&_jet_MUEF);
	outputTree->Branch("_jet_CHM",&_jet_CHM);
	outputTree->Branch("_jet_NHM",&_jet_NHM);
	outputTree->Branch("_jet_PHM",&_jet_PHM);
	outputTree->Branch("_jet_NM",&_jet_NM);

	outputTree->Branch("_PFcand_pt",&_PFcand_pt);
	outputTree->Branch("_PFcand_eta",&_PFcand_eta);
	outputTree->Branch("_PFcand_phi",&_PFcand_phi);
	outputTree->Branch("_PFcand_pdgId",&_PFcand_pdgId);
	outputTree->Branch("_PFcand_fromPV",&_PFcand_fromPV);



	outputTree->Branch("_met", &_met, "_met/D");
	outputTree->Branch("_met_phi", &_met_phi, "_met_phi/D");
	outputTree->Branch("_puppimet", &_puppimet, "_puppimet/D");
	outputTree->Branch("_puppimet_phi", &_puppimet_phi, "_puppimet_phi/D");


	outputTree->Branch("_elePt",      &_elePt);
	outputTree->Branch("_eleEta",     &_eleEta);
	outputTree->Branch("_elePhi",     &_elePhi);
	outputTree->Branch("_eleEnergy",  &_eleEnergy);
	outputTree->Branch("_eleEt",      &_eleEt);

	outputTree->Branch("gsf_VID_cutBasedElectronID_Fall17_94X_V1_loose",  &gsf_VID_cutBasedElectronID_Fall17_94X_V1_loose);
	outputTree->Branch("gsf_VID_cutBasedElectronID_Fall17_94X_V1_medium",  & gsf_VID_cutBasedElectronID_Fall17_94X_V1_medium);
	outputTree->Branch("gsf_VID_cutBasedElectronID_Fall17_94X_V1_tight", &gsf_VID_cutBasedElectronID_Fall17_94X_V1_tight);
	outputTree->Branch("gsf_VID_cutBasedElectronID_Fall17_94X_V1_veto", &gsf_VID_cutBasedElectronID_Fall17_94X_V1_veto);

	outputTree->Branch("_muonPt", &_muonPt);
	outputTree->Branch("_muonEta", &_muonEta);
	outputTree->Branch("_muonPhi", &_muonPhi);
	outputTree->Branch("_muonEnergy", &_muonEnergy);
	outputTree->Branch("_muonInnerTrackPt", & _muonInnerTrackPt);
	outputTree->Branch("_muonInnerTrackPtError", &_muonInnerTrackPtError);
	outputTree->Branch("_muonBestTrackPt", & _muonBestTrackPt);
	outputTree->Branch("_muonBestTrackPtError",&_muonBestTrackPtError);
	outputTree->Branch("_muonInnerTrackCharge", &_muonInnerTrackCharge);
	outputTree->Branch("_muonBestTrackCharge", & _muonBestTrackCharge);

	outputTree->Branch("_muonGlobal", &_muonGlobal);
	outputTree->Branch("_muonalgo", &_muonalgo);
	outputTree->Branch("_muonsegComp", &_muonsegComp);
	outputTree->Branch("_muonHighpurity", &_muonHighpurity);
	outputTree->Branch("_muonPF", &_muonPF);

	outputTree->Branch("_muonGlobalTrackPt", & _muonGlobalTrackPt);
	outputTree->Branch("_muonTPFMSTrackPt", & _muonTPFMSTrackPt);
	outputTree->Branch("_muonPickyTrackPt", & _muonPickyTrackPt);
	outputTree->Branch("_muonDYTTrackPt", & _muonDYTTrackPt);
	outputTree->Branch("_muonTunePMuonBestTrackPt", & _muonTunePMuonBestTrackPt);


	outputTree->Branch("_muonGlobalTrackPtError", & _muonGlobalTrackPtError);
	outputTree->Branch("_muonTPFMSTrackPtError", & _muonTPFMSTrackPtError);
	outputTree->Branch("_muonPickyTrackPtError", & _muonPickyTrackPtError);
	outputTree->Branch("_muonDYTTrackPtError", & _muonDYTTrackPtError);
	outputTree->Branch("_muonTunePMuonBestTrackPtError", & _muonTunePMuonBestTrackPtError);

	outputTree->Branch("_muonGlobalTrackCharge", & _muonGlobalTrackCharge);
	outputTree->Branch("_muonTPFMSTrackCharge", & _muonTPFMSTrackCharge);
	outputTree->Branch("_muonPickyTrackCharge", & _muonPickyTrackCharge);
	outputTree->Branch("_muonDYTTrackCharge", & _muonDYTTrackCharge);
	outputTree->Branch("_muonTunePMuonBestTrackCharge", & _muonTunePMuonBestTrackCharge);
	outputTree->Branch("_muonTrackIsoSumPt", &_muonTrackIsoSumPt);
	outputTree->Branch("_muonValidPixelHits", &_muonValidPixelHits);
	outputTree->Branch("_muonValidMuonHits", &_muonValidMuonHits);
	outputTree->Branch("_muonTrackerLayersHits", &_muonTrackerLayersHits);
	outputTree->Branch("_muonMatchedStations", &_muonMatchedStations);

	outputTree->Branch("_muonDxy", &_muonDxy);
	outputTree->Branch("_muonDz", &_muonDz);
	outputTree->Branch("_muonInnerDxy", &_muonInnerDxy);
	outputTree->Branch("_muonInnerDz", &_muonInnerDz);

	outputTree->Branch("_muonHighPtTrackerId", &_muonHighPtTrackerId);
	outputTree->Branch("_muonHighPtId", & _muonHighPtId);


	outputTree->Branch("_muonLoose", &_muonLoose);
	outputTree->Branch("_muonMedium", &_muonMedium);
	outputTree->Branch("_muonTight",&_muonTight);
	outputTree->Branch("_muonSoft",&_muonSoft);
	outputTree->Branch("_muonTracker",&_muonTracker);
	outputTree->Branch("_muonMatchedNStations",&_muonMatchedNStations);
	outputTree->Branch("_muonStationMask",&_muonStationMask);
	outputTree->Branch("_muonRPCLayers",&_muonRPCLayers);
	outputTree->Branch("_muonBestValidMuonHits",&_muonBestValidMuonHits);



	outputTree->Branch("_muonGlobTrackNChiS",&_muonGlobTrackNChiS);
	outputTree->Branch("_muonCombQualityChiSLP", &_muonCombQualityChiSLP);
	outputTree->Branch("_muonCombQualitytrkKink",&_muonCombQualitytrkKink);
	outputTree->Branch("_muonInnerTrackVFrac",& _muonInnerTrackVFrac);

	outputTree->Branch("_muonPixelLayers",&_muonPixelLayers);
	outputTree->Branch("_muonGoodMuon", &_muonGoodMuon);


}

// ------------ method called once each job just after ending the event loop  ------------
	void
METScanningNtupleMakerMINIAOD::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
METScanningNtupleMakerMINIAOD::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);

	//Specify that only 'tracks' is allowed
	//To use, remove the default given above and uncomment below
	//ParameterSetDescription desc;
	//desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
	//descriptions.addDefault(desc);
}

bool METScanningNtupleMakerMINIAOD::GetMETFilterDecision(const edm::Event& iEvent,edm::Handle<TriggerResults> METFilterResults, TString studiedfilter){

	if( !METFilterResults.failedToGet() ) {
		int N_MetFilters = METFilterResults->size();
		const edm::TriggerNames & metfilterName = iEvent.triggerNames(*METFilterResults);
		for( int i_Metfilter = 0; i_Metfilter < N_MetFilters; ++i_Metfilter ) {
			TString MetfilterPath =metfilterName.triggerName(i_Metfilter);
			//			std::cout<<"MetfilterPath:"<<MetfilterPath<< "\t index \t" << MetfilterPath.Index(studiedfilter) <<  std::endl;
			if(MetfilterPath.Index(studiedfilter) >=0)  return METFilterResults.product()->accept(i_Metfilter);

		}
	}
	return true; 
}


bool METScanningNtupleMakerMINIAOD::GetIdxFilterDecision(int it){
	if(it== idx_Flag_goodVertices)return  Flag_goodVertices;
	else if(it==  idx_Flag_globalTightHalo2016Filter)return   Flag_globalTightHalo2016Filter;
	else if(it==  idx_Flag_globalSuperTightHalo2016Filter)return   Flag_globalSuperTightHalo2016Filter;
	else if(it==  idx_Flag_HBHENoiseFilter)return   Flag_HBHENoiseFilter;
	else if(it==  idx_Flag_HBHENoiseIsoFilter)return   Flag_HBHENoiseIsoFilter;
	else if(it==  idx_Flag_EcalDeadCellTriggerPrimitiveFilter)return   Flag_EcalDeadCellTriggerPrimitiveFilter;
	else if(it==  idx_Flag_BadPFMuonFilter)return   Flag_BadPFMuonFilter;
	else if(it==  idx_Flag_BadChargedCandidateFilter)return   Flag_BadChargedCandidateFilter;
	else if(it==  idx_Flag_eeBadScFilter)return   Flag_eeBadScFilter;
	else if(it==  idx_Flag_ecalBadCalibFilter)return   Flag_ecalBadCalibFilter;
	else if(it==  idx_Flag_ecalLaserCorrFilter)return   Flag_ecalLaserCorrFilter;
	else if(it==  idx_Flag_EcalDeadCellBoundaryEnergyFilter)return   Flag_EcalDeadCellBoundaryEnergyFilter;
	else if(it==  idx_PassecalBadCalibFilter_Update)return   PassecalBadCalibFilter_Update;
	else if(it==  idx_PassecalLaserCorrFilter_Update)return   PassecalLaserCorrFilter_Update;
	else if(it==  idx_PassEcalDeadCellBoundaryEnergyFilter_Update)return   PassEcalDeadCellBoundaryEnergyFilter_Update;
	else if(it==  idx_PassBadChargedCandidateFilter_Update)return   PassBadChargedCandidateFilter_Update;
        else if(it==  idx_PassBadPFMuonFilter_Update) return PassBadPFMuonFilter_Update;
	else return false;
}

TString METScanningNtupleMakerMINIAOD::GetIdxFilterName(int it){
	if(it==idx_Flag_goodVertices)return "Flag_goodVertices";
	else if(it== idx_Flag_globalTightHalo2016Filter)return "Flag_globalTightHalo2016Filter";
	else if(it== idx_Flag_globalSuperTightHalo2016Filter)return "Flag_globalSuperTightHalo2016Filter";
	else if(it== idx_Flag_HBHENoiseFilter)return "Flag_HBHENoiseFilter";
	else if(it== idx_Flag_HBHENoiseIsoFilter)return "Flag_HBHENoiseIsoFilter";
	else if(it== idx_Flag_EcalDeadCellTriggerPrimitiveFilter)return "Flag_EcalDeadCellTriggerPrimitiveFilter";
	else if(it== idx_Flag_BadPFMuonFilter)return "Flag_BadPFMuonFilter";
	else if(it== idx_Flag_BadChargedCandidateFilter)return "Flag_BadChargedCandidateFilter";
	else if(it== idx_Flag_eeBadScFilter)return "Flag_eeBadScFilter";
	else if(it== idx_Flag_ecalBadCalibFilter)return "Flag_ecalBadCalibFilter";
	else if(it== idx_Flag_ecalLaserCorrFilter)return "Flag_ecalLaserCorrFilter";
	else if(it== idx_Flag_EcalDeadCellBoundaryEnergyFilter)return "Flag_EcalDeadCellBoundaryEnergyFilter";
	else if(it== idx_PassecalBadCalibFilter_Update)return "PassecalBadCalibFilter_Update";
	else if(it== idx_PassecalLaserCorrFilter_Update)return "PassecalLaserCorrFilter_Update";
	else if(it== idx_PassEcalDeadCellBoundaryEnergyFilter_Update )return "PassEcalDeadCellBoundaryEnergyFilter_Update";
	else if(it== idx_PassBadChargedCandidateFilter_Update )return "PassBadChargedCandidateFilter_Update";
        else if(it==  idx_PassBadPFMuonFilter_Update) return "PassBadPFMuonFilter_Update";
	else return "";
}


void METScanningNtupleMakerMINIAOD::InitandClearStuff(){

	Flag_goodVertices=false;
	Flag_globalTightHalo2016Filter=false;
	Flag_globalSuperTightHalo2016Filter=false;
	Flag_HBHENoiseFilter=false;
	Flag_HBHENoiseIsoFilter=false;
	Flag_EcalDeadCellTriggerPrimitiveFilter=false;
	Flag_BadPFMuonFilter=false;
	Flag_BadChargedCandidateFilter=false;
	Flag_eeBadScFilter=false;
	Flag_ecalBadCalibFilter=false;
	Flag_ecalLaserCorrFilter=false;
	Flag_EcalDeadCellBoundaryEnergyFilter=false;
	PassecalBadCalibFilter_Update=false;
	PassecalLaserCorrFilter_Update=false;
	PassEcalDeadCellBoundaryEnergyFilter_Update=false;
	PassBadChargedCandidateFilter_Update=false;
	PassBadPFMuonFilter_Update=false;

	_jetEta.clear();
	_jetPhi.clear();
	_jetPt.clear();
	_jetRawPt.clear();
	_jet_CHEF.clear();
	_jet_NHEF.clear();
	_jet_NEEF.clear();
	_jet_CEEF.clear();
	_jet_MUEF.clear();
	_jet_CHM.clear();
	_jet_NHM.clear();
	_jet_PHM.clear();
	_jet_NM.clear();

	_PFcand_pt.clear();
	_PFcand_eta.clear();
	_PFcand_phi.clear();
	_PFcand_pdgId.clear();
	_PFcand_fromPV.clear();

	_elePt.clear();
	_eleEta.clear();
	_elePhi.clear();
	_eleEnergy.clear();
	_eleEt.clear();

	gsf_VID_cutBasedElectronID_Fall17_94X_V1_loose.clear();
	gsf_VID_cutBasedElectronID_Fall17_94X_V1_medium.clear();
	gsf_VID_cutBasedElectronID_Fall17_94X_V1_tight.clear();
	gsf_VID_cutBasedElectronID_Fall17_94X_V1_veto.clear();

	_muonPt.clear();
	_muonEta.clear();
	_muonPhi.clear();
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
	_muonInnerDxy.clear();
	_muonInnerDz.clear();

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


	_muonGlobTrackNChiS.clear();
	_muonCombQualityChiSLP.clear();
	_muonCombQualitytrkKink.clear();
	_muonInnerTrackVFrac.clear();

	_muonPixelLayers.clear();
	_muonGoodMuon.clear();


}

std::vector<std::string> METScanningNtupleMakerMINIAOD::split(const std::string &s, char delim = ':'){

	std::vector<std::string> elems;

	std::stringstream ss(s);
	std::string item;
	while (std::getline(ss, item, delim))
	{
		elems.push_back(item);
	}
	return elems;

}

//define this as a plug-in
DEFINE_FWK_MODULE(METScanningNtupleMakerMINIAOD);

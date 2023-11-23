// -*- C++ -*-
//
// Package:    PromptAna/PromptAnalyzer
// Class:      PromptAnalyzer
//
/**\class PromptAnalyzer PromptAnalyzer.cc PromptAna/PromptAnalyzer/plugins/PromptAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Robert Ciesielski
//         Created:  Wed, 27 Jun 2018 16:18:44 GMT
//        Modified:  Luiz Emediato - 24 Nov 2020 
//
// system include files
#include <memory>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <sstream>
#include <istream>
#include <fstream>
#include <iomanip>
#include <map>

//...reading files ...Luiz
#include "FWCore/ParameterSet/interface/FileInPath.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
// dEdx
#include "DataFormats/TrackReco/interface/DeDxHit.h"
#include "DataFormats/TrackReco/interface/DeDxHitInfo.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TLorentzVector.h"
//#include "DataFormats/Math/interface/deltaR.h"

#include "TCutG.h"
#include <TROOT.h>

// TTree memmory 1TB ...Luiz
//#include "TTree.h"

// pixel clusters
//#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
//#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
//#include "Geometry/Records/interface/TrackerTopologyRcd.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"

//PPS
#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"

#include "DataFormats/CTPPSDetId/interface/TotemRPDetId.h"
#include "DataFormats/CTPPSReco/interface/CTPPSLocalTrackLite.h"

// PFCandidates
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
// V0Producer: PFCandidate.h is enough ...Luiz 
//#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
//#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"

// Muon Includes
//#include "DataFormats/MuonReco/interface/MuonFwd.h"
//#include "DataFormats/MuonReco/interface/Muon.h"
//#include "DataFormats/MuonReco/interface/MuonSelectors.h"

// Ferenc's PID ...Luiz
#include "UserCode/EnergyLossPID/interface/ParticleType.h"
#include "UserCode/EnergyLossPID/interface/MostProbable.h"

// Forward Protons ...Luiz
//#include "DataFormats/ProtonReco/interface/ForwardProton.h"
//#include "DataFormats/ProtonReco/interface/ForwardProtonFwd.h"

#define M_LN10 2.30258509299404568402
#define Sqr(x) ((x) * (x))

//-----------------------------

// optics
const double v_x_R_1_F = -2.24791387053766;  const double L_x_R_1_F = 0.125396407127792E3; 
const double v_y_R_1_F = +0.025781593410852; const double L_y_R_1_F = 238.517247191010E3;
const double v_x_R_2_F = -1.92610996810677;  const double L_x_R_2_F = -3.00655323980445E3; 
const double v_y_R_2_F = -0.000000021508565; const double L_y_R_2_F = 271.511335947517E3;

const double v_x_L_1_F = -2.24791387053766;  const double L_x_L_1_F = 0.125396407127792E3;
const double v_y_L_1_F = +0.025781593410852; const double L_y_L_1_F = 238.517247191010E3;
const double v_x_L_2_F = -1.92610996810677;  const double L_x_L_2_F = -3.00655323980445E3; 
const double v_y_L_2_F = -0.000000021508565; const double L_y_L_2_F = 271.511335947517E3;


//double m_pi=0.13957;
double m_k = 0.493677;
//...Luiz
double m_k0 = 0.497611;
//double m_mu = 0.1056583715;
//double m_p =0.93827;
double m_rho = 0.77;
double m_phi = 1.019461;
//
double m_x1070 = 1.072; //PDG
double m_x1235 = 1.235; //from fit

//...Luiz
// new PDG
double m_pi = 0.13957061;
// new PDG
double m_mu = 0.1056583745;
double m_e = 0.0005109989461;
double m_p = 0.9382720813;
double m_kstar=0.89555;  // neutral

//...using Ferenc's PID: pion is 2 now, kaon is 3 ...Luiz
//enum EPID { pidUnknown, pidProton, pidKaon, pidPion }; //obsolete
//
//.............0...........1............2........3........4.......
//enum EPID { pidUnknown, pidElectron, pidPion, pidKaon, pidProton };

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


//using reco::TrackCollection;

using namespace edm;
using namespace reco;
using namespace std;

class PromptAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit PromptAnalyzer(const edm::ParameterSet&);
      ~PromptAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
  
      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);

      bool jsonLocal(int r, int ls);

      // ----------member data ---------------------------
  //      edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
  edm::EDGetTokenT<reco::TrackCollection> trkToken_;
  edm::EDGetTokenT<vector<CTPPSLocalTrackLite> > RPtrkToken_;
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  edm::EDGetTokenT<reco::BeamSpot> beamspotToken_;
  edm::EDGetTokenT<edm::TriggerResults>  trigToken_;
  // V0 ...Luiz
  edm::EDGetTokenT<reco::VertexCompositeCandidateCollection> kshortsToken_;
  edm::EDGetTokenT<reco::VertexCompositeCandidateCollection> lambdasToken_;
  edm::EDGetTokenT<reco::DeDxDataValueMap> dedxsToken_;
  edm::EDGetTokenT<reco::DeDxDataValueMap> dedxPIXsToken_;
  //edm::EDGetTokenT<reco::DeDxDataValueMap> dedxpixelsToken_;
  //
  //edm::EDGetTokenT<reco::ForwardProtonCollection> RecoProtonsSingleRPToken_;
  //edm::EDGetTokenT<reco::ForwardProtonCollection> RecoProtonsMultiRPToken_;
  
  HLTConfigProvider hltConfig_;

  map<string,TH1F*> histosTH1F;
  map<string,TH2F*> histosTH2F;

  //...py effic ...Luiz
  std::vector<std::vector<double>> efvec;

};
  
//                                                                                                             
// constants, enums and typedefs                                                                               
//                                                                                                             

//                                                                                                             
// static data member definitions                                                                              
//                                                                                                             

//                                                                                                             
// constructors and destructor ...including kshorts and lambdas ...Luiz                                     
//
// TAGs from python configuration : TOKEN <-> TAG
PromptAnalyzer::PromptAnalyzer(const edm::ParameterSet& iConfig)
  :
  trkToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks")))
  ,RPtrkToken_(consumes<vector<CTPPSLocalTrackLite> >(iConfig.getParameter<edm::InputTag>("RPtracks")))
  ,vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices")))
  ,beamspotToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamspot")))
  ,trigToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggers")))
  // V0 ...Luiz                                                                                                
  ,kshortsToken_(consumes<reco::VertexCompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("kshorts")))
  ,lambdasToken_(consumes<reco::VertexCompositeCandidateCollection>(iConfig.getParameter<edm::InputTag>("lambdas")))
  // dedx ...Luiz
  ,dedxsToken_(consumes<reco::DeDxDataValueMap>(iConfig.getParameter<edm::InputTag>("dedxs")))
  ,dedxPIXsToken_(consumes<reco::DeDxDataValueMap>(iConfig.getParameter<edm::InputTag>("dedxPIXs")))
  ////,dedxpixelsToken_(consumes<reco::DeDxDataValueMap>(iConfig.getParameter<edm::InputTag>("dedxpixels")))
  // reco RP protons ...Luiz
  //,RecoProtonsSingleRPToken_(consumes<reco::ForwardProtonCollection>(iConfig.getParameter<edm::InputTag>("tagRecoProtonsSingleRP")))
  //,RecoProtonsMultiRPToken_(consumes<reco::ForwardProtonCollection>(iConfig.getParameter<edm::InputTag>("tagRecoProtonsMultiRP")))
{
  
  //now do what ever initialization is needed

  //usesResource("TFileService");
  //edm::Service<TFileService> fs;

  /*

  //----------------------------------

  //--------------------------------------

  */

  
  /* 
  const std::vector< const edm::ValueMap<reco::DeDxData> *>& v_dEdx;
  const edm::RefToBase<reco::Track>& trackref;
  for (unsigned int i=0; i<v_dEdx.size(); i++) {
     const edm::ValueMap<reco::DeDxData>& dEdxTrack = *(v_dEdx[i]);
     const reco::DeDxData& dedx = dEdxTrack[trackref];
     histograms.h_dedx_estim[count][i].fill(dedx.dEdx());
  }
  */
  
}


PromptAnalyzer::~PromptAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

////edm::Service<TFileService> fs;

//
// member functions
//

// ------------ method called for each event  ------------
void
PromptAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  //   using namespace edm;

   //--------------------------
   
  edm::Handle<TrackCollection> tracks;
  edm::Handle<vector<CTPPSLocalTrackLite> > RPtracks;
  edm::Handle<VertexCollection> vertices;
  edm::Handle<reco::BeamSpot> beamspot;
  edm::Handle<edm::TriggerResults> triggers;
  // kshorts and lambdas ...Luiz
  edm::Handle<reco::VertexCompositeCandidateCollection> kshorts;
  edm::Handle<reco::VertexCompositeCandidateCollection> lambdas;
  edm::Handle<reco::DeDxDataValueMap> dedxs;
  edm::Handle<reco::DeDxDataValueMap> dedxPIXs;
  ////edm::Handle<reco::DeDxDataValueMap> dedxpixels;
  
  // reco RP protons ...Luiz
  //edm::Handle<reco::ForwardProtonCollection> ProtonsSingleRP;
  //edm::Handle<reco::ForwardProtonCollection> ProtonsMultiRP;

  
  iEvent.getByToken(trkToken_,tracks);
  iEvent.getByToken(RPtrkToken_,RPtracks);
  iEvent.getByToken(vtxToken_,vertices);
  iEvent.getByToken(beamspotToken_,beamspot);
  iEvent.getByToken(trigToken_,triggers); 
  // kshorts and lambdas ...Luiz
  iEvent.getByToken(kshortsToken_,kshorts);
  iEvent.getByToken(lambdasToken_,lambdas);
  iEvent.getByToken(dedxsToken_,dedxs);
  iEvent.getByToken(dedxPIXsToken_,dedxPIXs);
  ////iEvent.getByToken(dedxpixelsToken_,dedxpixels);
  // reco RP protons ...Luiz
  //iEvent.getByToken(RecoProtonsSingleRPToken_,ProtonsSingleRP);
  //iEvent.getByToken(RecoProtonsMultiRPToken_,ProtonsMultiRP);

  /*
  const std::vector< const edm::ValueMap<reco::DeDxData> *>& v_dEdx;
  const edm::RefToBase<reco::Track>& trackref;
  
  for (unsigned int i=0; i<v_dEdx.size(); i++) {
     const edm::ValueMap<reco::DeDxData>& dEdxTrack = *(v_dEdx[i]);
     const reco::DeDxData& dedx = dEdxTrack[trackref];
     histograms.h_dedx_estim[count][i].fill(dedx.dEdx());
  }
  */
  
  //------------------------------------------------------------------  
  /*
  std::cout<<"check triggers"<<std::endl;
  
  if (triggers.isValid()) {
    edm::TriggerNames triggerNames = iEvent.triggerNames(*triggers);
    unsigned int size = triggers->size();
    
    for(unsigned ij = 0; ij<size; ++ij) {
      std::string name = triggerNames.triggerName(ij);
      //      const char* variab1 = name.c_str(); 
      
      if( triggers->accept(ij) ){
	std::cout<<ij<<") "<<name<<" accept="<<triggers->accept(ij)<<std::endl;
      }
    }
  }
  */
  //------------------------------------------------------------------  
  //TB/BT or TT/BB are now separately in TOTEM2 or TOTEM4
  // 0 = for TB in TOTEM2 or TT in TOTEM4
  // 1 = for BT in TOTEM2 or BB in TOTEM4
  
  int runnr =  iEvent.id().run();
  long int event =  iEvent.id().event();

  histosTH1F["hrun"]->Fill(runnr);

  //std::cout << "run# = " << runnr  << " event = " << event << std::endl;
  
  //  int eventnr =  iEvent.id().event();
  int LS = iEvent.luminosityBlock();  

  //LS histos were here

  int ntrk0 = 0;
  int ntrk  = 0;

  int ndedx = 0;
	      
  int totcharge=0;

  //tracks in 2track-events (npixelhits>0)
  //TLorentzVector pi1(0.,0.,0.,0.);
  //TLorentzVector pi2(0.,0.,0.,0.);
  TLorentzVector pipiRec(0.,0.,0.,0.);
  //...Luiz
  TLorentzVector Pi1(0.,0.,0.,0.);
  TLorentzVector Pi2(0.,0.,0.,0.);

  //tracks in 4track-events (npixelhits>0)
  TLorentzVector pi4pos1(0.,0.,0.,0.);
  TLorentzVector pi4pos2(0.,0.,0.,0.);
  TLorentzVector pi4neg1(0.,0.,0.,0.);
  TLorentzVector pi4neg2(0.,0.,0.,0.);
  TLorentzVector pi4Rec(0.,0.,0.,0.);
  
  int ntrk4pos=0;
  int ntrk4neg=0;

  double xBS,yBS,zBS;
  if (beamspot.isValid()){
    xBS = beamspot->x0();
    yBS = beamspot->y0();
    zBS = beamspot->z0();
  }else{
    std::cout<<"sorry, no beamspot"<<std::endl;
    xBS=0;
    yBS=0;
    zBS=0;
  }

  // my cuts ...Luiz
  bool fiducialRegion4 = false;
  bool fiducialRegionPt4 = false;

  ////bool fiducialRegionK4 = false;
  ////bool fiducialRegionPtK4 = false;

  //  double etaCut = 2.5; .....................checking |eta|<3.0 ...it improves V0 reco
  //double etaCut = 2.5;
  double etaCut = 3.0;
  double ptCut = 0.0;

  int pid = 0;
  int pidV0 = 0;
  int pidno = 0;
  int pidV0no = 0;

  double mydedx = 0.0;
  double mydedxPIX = 0.0;
  double mydedxPIX4 = 0.0;
  double mydedxPIXno = 0.0;


  int npiA = 0;
  int npiB = 0;
  int npiC = 0;
  int npiD = 0;
  int nkA  = 0;
  int nkB  = 0;
  int nkC  = 0;
  int nkD  = 0;

  
  double mrecKstar=0.0;
  double mrecKstarbar=0.0;

  
  // my tracks in 4track-events (npixelhits>0)
  // pions & kaons ...Luiz
  TLorentzVector pi1(0.,0.,0.,0.);
  TLorentzVector pi2(0.,0.,0.,0.);
  TLorentzVector pi3(0.,0.,0.,0.);
  TLorentzVector pi4(0.,0.,0.,0.);
  //
  TLorentzVector piA(0.,0.,0.,0.);
  TLorentzVector piB(0.,0.,0.,0.);
  TLorentzVector piC(0.,0.,0.,0.);
  TLorentzVector piD(0.,0.,0.,0.);
  //
  TLorentzVector k1(0.,0.,0.,0.);
  TLorentzVector k2(0.,0.,0.,0.);
  TLorentzVector k3(0.,0.,0.,0.);
  TLorentzVector k4(0.,0.,0.,0.);
  //
  TLorentzVector kA(0.,0.,0.,0.);
  TLorentzVector kB(0.,0.,0.,0.);
  TLorentzVector kC(0.,0.,0.,0.);
  TLorentzVector kD(0.,0.,0.,0.);
  //
  TLorentzVector pi1pi2Rec(0.,0.,0.,0.);
  TLorentzVector pi3pi4Rec(0.,0.,0.,0.);
  TLorentzVector pi1pi3Rec(0.,0.,0.,0.);
  TLorentzVector pi2pi4Rec(0.,0.,0.,0.);
  TLorentzVector pi1pi4Rec(0.,0.,0.,0.);
  TLorentzVector pi2pi3Rec(0.,0.,0.,0.);
  // 
  TLorentzVector k1k2Rec(0.,0.,0.,0.);
  TLorentzVector k3k4Rec(0.,0.,0.,0.);
  TLorentzVector k1k3Rec(0.,0.,0.,0.);
  TLorentzVector k2k4Rec(0.,0.,0.,0.);
  TLorentzVector k1k4Rec(0.,0.,0.,0.);
  TLorentzVector k2k3Rec(0.,0.,0.,0.);

  // checking
  TLorentzVector pi1234Rec(0.,0.,0.,0.);
  TLorentzVector pi1324Rec(0.,0.,0.,0.);
  TLorentzVector pi1423Rec(0.,0.,0.,0.);

  //...combining pions and kaons for the event selection type = 11 (one primary & one Vee)
   /* 
  ...first combining, then select the Q_pairs=0
  pi1pi2 pi3k4
  pi1pi3 pi2k4
  pi2pi3 pi1k4

  pi1pi2 k3pi4
  pi1pi4 k3pi2
  pi2pi4 k3pi1

  pi1k2 pi3pi4
  pi3k2 pi1pi4
  pi4k2 pi1pi3

  k1pi2 pi3pi4
  k1pi3 pi2pi4
  k1pi4 pi2pi3
   */

  //
  TLorentzVector pi3k4Rec(0.,0.,0.,0.);
  TLorentzVector pi2k4Rec(0.,0.,0.,0.);
  TLorentzVector pi1k4Rec(0.,0.,0.,0.);
  //
  TLorentzVector k3pi4Rec(0.,0.,0.,0.);
  TLorentzVector k3pi2Rec(0.,0.,0.,0.);
  TLorentzVector k3pi1Rec(0.,0.,0.,0.);
  //
  TLorentzVector pi1k2Rec(0.,0.,0.,0.);
  TLorentzVector pi3k2Rec(0.,0.,0.,0.);
  TLorentzVector pi4k2Rec(0.,0.,0.,0.);
  //
  TLorentzVector k1pi2Rec(0.,0.,0.,0.);
  TLorentzVector k1pi3Rec(0.,0.,0.,0.);
  TLorentzVector k1pi4Rec(0.,0.,0.,0.);
  //
  TLorentzVector k1pi2pi3pi4Rec(0.,0.,0.,0.);
  TLorentzVector pi1k2pi3pi4Rec(0.,0.,0.,0.);
  TLorentzVector pi1pi2k3pi4Rec(0.,0.,0.,0.);
  TLorentzVector pi1pi2pi3k4Rec(0.,0.,0.,0.);
  //
	   //...K*Kbar* channel
	   //1234
            TLorentzVector k1pi2k3pi4Rec(0.,0.,0.,0.);
            TLorentzVector k1pi2pi3k4Rec(0.,0.,0.,0.);
            TLorentzVector pi1k2k3pi4Rec(0.,0.,0.,0.);
            TLorentzVector pi1k2pi3k4Rec(0.,0.,0.,0.);
	   //1324
            TLorentzVector k1pi3k2pi4Rec(0.,0.,0.,0.);
            TLorentzVector k1pi3pi2k4Rec(0.,0.,0.,0.);
            TLorentzVector pi1k3k2pi4Rec(0.,0.,0.,0.);
            TLorentzVector pi1k3pi2k4Rec(0.,0.,0.,0.);
	   //1423
            TLorentzVector k1pi4k2pi3Rec(0.,0.,0.,0.);
            TLorentzVector k1pi4pi2k3Rec(0.,0.,0.,0.);
            TLorentzVector pi1k4k2pi3Rec(0.,0.,0.,0.);
            TLorentzVector pi1k4pi2k3Rec(0.,0.,0.,0.);
  // !
  TLorentzVector k2pi4Rec(0.,0.,0.,0.);
  TLorentzVector k2pi3Rec(0.,0.,0.,0.);
  TLorentzVector pi2k3Rec(0.,0.,0.,0.);
  //
  TLorentzVector pi1k3Rec(0.,0.,0.,0.);
  //

  //
  TLorentzVector pipipipiRec(0.,0.,0.,0.);
  TLorentzVector kkkkRec(0.,0.,0.,0.);
  //
  
  int charray[4]={0,0,0,0};
  double chi2array[4]={0.,0.,0.,0.};
  double d0array[4]={0.,0.,0.,0.};
  double dzarray[4]={0.,0.,0.,0.};
  int pidarray[4]={0,0,0,0};
  int pidv0array[4]={0,0,0,0};
  double vtxdxyarray[4]={0.,0.,0.,0.};
  double vtxdzarray[4]={0.,0.,0.,0.};
  double phiarray[4]={0.,0.,0.,0.};

  double trkpx[4]={0.,0.,0.,0.};
  double trkpy[4]={0.,0.,0.,0.};
  double trkpz[4]={0.,0.,0.,0.};
  double atrkpx[4]={0.,0.,0.,0.};
  double atrkpy[4]={0.,0.,0.,0.};
  double atrkpz[4]={0.,0.,0.,0.};
  
  //...ordering
  int arraych[4]={0,0,0,0};
  double arraychi2[4]={0.,0.,0.,0.};
  double arrayd0[4]={0.,0.,0.,0.};
  double arraydz[4]={0.,0.,0.,0.};
  int arraypid[4]={0,0,0,0};
  int arraypidv0[4]={0,0,0,0};
  double arrayvtxdxy[4]={0.,0.,0.,0.};
  double arrayvtxdz[4]={0.,0.,0.,0.};
  double arrayphi[4]={0.,0.,0.,0.};
  
       // ...transverse impact parameter distribution 
       double d01 = 0.0; 
       double d02 = 0.0;
       double d03 = 0.0;
       double d04 = 0.0;
       // ...longitudinal impact parameter distribution 
       double dz1 = 0.0; 
       double dz2 = 0.0;
       double dz3 = 0.0;
       double dz4 = 0.0;
       // ...transverse impact parameter distribution 
       double vtxdxy1 = 0.0; 
       double vtxdxy2 = 0.0;
       double vtxdxy3 = 0.0;
       double vtxdxy4 = 0.0;
       // ...longitudinal impact parameter distribution 
       double vtxdz1 = 0.0; 
       double vtxdz2 = 0.0;
       double vtxdz3 = 0.0;
       double vtxdz4 = 0.0;

       
  //...Luiz
  // Get primary vertex PV
  //if (vertices->empty()) return;  // skip the event if no PV found
  //const reco::Vertex& pv = vertices->begin()->position();

    //...Luiz
    //...to avoid 'not declared in this scope'
    //math::XYZPoint ppv(0.,0.,0.);
    //math::XYZPoint pv(0.,0.,0.);
    double pvx = 0.;
    double pvy = 0.;
    double pvz = 0.;
    //int visfake = vertices->begin()->isFake();

    //...what does it mean vertices->empty() ???

    //...notice that *vertex* IS NOT *primary vertex*
    //   this difference is clear in the software
    //   vertex is the interaction point (0,0,0)
    //   track referencePoint() is the PCA w.r.t. vertex or (0,0,0)
    //   track referencePoint(pv) is the PCA w.r.t. primary vertex
    
    //...there are no vertices->empty events 
    //...I am disabling the IF statement

    //...vertices are primary vertices from the collection and we do have events with no primaries
    
    //if ( !( vertices->empty() ) ){
    //math::XYZPoint pvtx = itTrack->referencePoint();  // PCA w.r.t. the center of CMS or beamspot (0,0,0) 
    //math::XYZPoint ppv = (*vertices)[0].position();   // w.r.t. primary vertex
    //math::XYZPoint pv = vertices->begin()->position();// w.r.t. primary vertex
    //ppv = (*vertices)[0].position();
       
    pvx =  vertices->begin()->x();
    pvy =  vertices->begin()->y();
    pvz =  vertices->begin()->z();

    /*
    math::XYZPoint mypv(pvx,pvy,pvz); //...oh!

    pv = mypv; //...  :-(
    */

    math::XYZPoint pv(pvx,pvy,pvz); //...scrumdiddlyumptious! :-)
    
    //} //..IF above

    /* I want to see it! ...it works flawlessly!
    std::cout << " --- primary vertex position ----------------" << std::endl;
    std::cout << "pv  = " << pv  << std::endl;
    std::cout << "pvx = " << pvx << std::endl;
    std::cout << "pvy = " << pvy << std::endl;
    std::cout << "pvz = " << pvz << std::endl;
    */

    //...checking!
    if ( vertices->empty() ){
      std::cout << " +++ vertices +++ empty +++++++++++++++++  " << std::endl;
    }


    
    
    // ...loop over dedxPIXs
    //const edm::RefToBase<reco::Track>& trackref ;
    //const std::vector< const edm::ValueMap<reco::DeDxData> *>& v_dEdx;
    //
    //reco::DeDxHitInfoRef dedxHitsRef = dedxCollH->get(track.key());
    //
    //const edm::RefToBase<reco::Track>& trackref;
    //const std::vector< const edm::ValueMap<reco::DeDxData> *>& v_dEdx ;

    /*
    //...loop through dedx
    for(DeDxDataValueMap::const_iterator itDeDx = dedxPIXs->begin();itDeDx != dedxPIXs->end();++itDeDx) {

      int ndedxPIXsize = itDeDx.size();
      std::cout << " ndedxPIXsize = " << ndedxPIXsize << std::endl;
    
    }
    */

    
    
    //...instantiating for PID  ...Luiz
    //
    int mypass = -99;
    MostProbable* mymostProbable = new MostProbable();
    ParticleType* myparticleType = new ParticleType(mymostProbable, mypass);
    //
    //...getting dEdx
    const edm::ValueMap<reco::DeDxData> dedxPIXTrack = *dedxPIXs.product();
    

    
  for(TrackCollection::const_iterator itTrack = tracks->begin();itTrack != tracks->end();++itTrack) {

    int looper = itTrack->isLooper();  
    double p = itTrack->p();          // track momentum
    double pt = itTrack->pt();  
    double pz = itTrack->pz();  
    double eta = itTrack->eta();  
    double phi = itTrack->phi();  
    int charge = itTrack->charge();
    int npixelhits = itTrack->hitPattern().numberOfValidPixelHits();
    //    int nstriphits = itTrack->hitPattern().numberOfValidStripHits();
    int algo = itTrack->algo();
    double chi2 = itTrack->normalizedChi2();
    double d0 = itTrack->d0();
    double dz = itTrack->dz();

    double px = itTrack->px();  
    double py = itTrack->py();  
    //
    //double ptt = TMath::Sqrt( px*px + py*py );
    //std::cout << " ------" << std::endl;
    //std::cout << "ptt = " << ptt << std::endl;
    //std::cout << "pt = " << pt << std::endl;
     
    // Bene: One can use TransientTracks to estimate track impact parameters with
    // respect to the beam line or primary vertex, taking into account the curvature of the track.

    //...Luiz
    //math::XYZPoint pv = itTrack->referencePoint(); // PCA w.r.t. the center of CMS or beamspot (0,0,0) 
    //double pvx =  itTrack->vx();
    //double pvy =  itTrack->vy();
    //double pvz =  itTrack->vz();
    double vdxy = itTrack->dxy(pv); // vtxdxy : transverse impact parameter
    double vdz =  itTrack->dz(pv);  // vtxdz  : longitudinal impact parameter
    
    //...Bene: be sure the pt of the tracks is bigger than 0.5 GeV (?)

    /*
    std::cout.precision(10);
    std::cout << " --- primary vertex ----------------" << std::endl;
    std::cout << "pv  = " << pv  << std::endl;
    std::cout << "pvx = " << pvx << std::endl;
    std::cout << "pvy = " << pvy << std::endl;
    std::cout << "pvz = " << pvz << std::endl;
    */
    
    /*
    std::cout << " *** transverse i.p.   *** " << std::endl;
    std::cout << "vdxy = " << vdxy << std::endl;
    std::cout << " --- longitudinal i.p. --- " << std::endl;
    std::cout << "vdz  = " << vdz << std::endl;
    */
    
    /*
    //...considering the track curvature
    //
    // with respect to any specified vertex, such as primary vertex
    GlobalPoint vert(pv.x(), pv.y(), pv.z());
    TrajectoryStateClosestToPoint  traj = itTrack->trajectoryStateClosestToPoint(vert );
    double dd0 = traj.perigeeParameters().transverseImpactParameter();
    double zz0 = traj.perigeeParameters().longitudinalImpactParameter();
    */

   
    //...attention here!
    //// comment this line out only for the M(K0) window technique
    if(npixelhits>0){

      histosTH1F["hpt"]->Fill(pt);
      histosTH1F["heta"]->Fill(eta);
      histosTH1F["hphi"]->Fill(phi);
      histosTH1F["halgo"]->Fill(algo);
      histosTH1F["hnhits"]->Fill(npixelhits);
      
      histosTH1F["hlooper"]->Fill(looper);
      histosTH1F["hchi2"]->Fill(chi2);
      histosTH1F["hd0"]->Fill(d0);
      histosTH1F["hdz"]->Fill(dz);
      

      math::XYZPoint point(xBS,yBS,zBS);
      double d0BS = itTrack->dxy(point);
      double dzBS = itTrack->dz(point);

      histosTH1F["hd0BS"]->Fill(d0BS);
      histosTH1F["hdzBS"]->Fill(dzBS);

  
      totcharge += charge;

      // pions 
      double ene=TMath::Sqrt(pt*pt+pz*pz+m_pi*m_pi);
      TLorentzVector trk_lorentz(itTrack->px(),itTrack->py(),itTrack->pz(),ene);
     
      //--------------------------------------
      // 2 trk
      //double ene=TMath::Sqrt(pt*pt+pz*pz+m_pi*m_pi);
      //TLorentzVector trk_lorentz(itTrack->px(),itTrack->py(),itTrack->pz(),ene);
      pipiRec += trk_lorentz; 
      
      //if(ntrk==0) pi1 = trk_lorentz;
      //if(ntrk==1) pi2 = trk_lorentz;
      //...Luiz
      if(ntrk==0) Pi1 = trk_lorentz;
      if(ntrk==1) Pi2 = trk_lorentz;

      //--------------------------------------
      // 4trk

      pi4Rec += trk_lorentz; 
	   
      if(charge>0){
	if(ntrk4pos==0) pi4pos1 = trk_lorentz;
	if(ntrk4pos==1) pi4pos2 = trk_lorentz;
	     
	ntrk4pos++;
      }
	   
      if(charge<0){
	if(ntrk4neg==0) pi4neg1 = trk_lorentz;
	if(ntrk4neg==1) pi4neg2 = trk_lorentz;
	
	ntrk4neg++;
      }
      //--------------------------------------
      
      //--------------------------------------
      // my 4trk definitions ...Luiz

	   pipipipiRec += trk_lorentz; 

	   //...first, tagging by track number
	   if(ntrk==0) piA = trk_lorentz;
	   if(ntrk==1) piB = trk_lorentz;
	   if(ntrk==2) piC = trk_lorentz;
	   if(ntrk==3) piD = trk_lorentz;

	   
	   
           // PID ...Luiz
	   //
	   // pid = particleType->guess(p,values[all].first);
	   //
	   //const reco::TrackRef t = itTrack-> ;
	   //const reco::DeDxData& dedxPIXObj = dedxPIXTrack[t];
	   //mydedxPIX = dedxPIXObj.dEdx();
	   //
	   reco::TrackRef aTrackRef0 = reco::TrackRef(tracks, ntrk);
           const reco::DeDxData& dedxPIXObj0 = dedxPIXTrack[aTrackRef0];	   
           //const reco::DeDxData& dedxObj0 = dedxTrack[aTrackRef0];	   
	   mydedxPIX = dedxPIXObj0.dEdx();
	   //
	   //...
	   pid = myparticleType->guess( p , mydedxPIX );
	   //...V0
	   pidV0 = myparticleType->sure( p , mydedxPIX );
	   //
	   //...must this *object* be deleted?
	   //
	   ////std::cout << " mydedxPIX = " << mydedxPIX << std::endl;
	   ////std::cout << " pid = " << pid << std::endl;
	   ////std::cout << " pidV0 = " << pidV0 << std::endl;

	   //...dEdx histogram   ...ntrk=any
	   // histosTH2F["hdedx"]->Fill(itTrack->p,itTrack->harmonic2_dEdx);
	   histosTH2F["hdedxPIX"]->Fill( p , mydedxPIX );
	   
	   
	   //...beware of the index here: ntrk==4 for 4-track events
	   if(ntrk==0){
	     charray[0]=charge;	
	     chi2array[0]=chi2;
	     d0array[0]=d0;
	     dzarray[0]=dz;
	     pidarray[0]=pid; pidv0array[0]=pidV0;
	     vtxdxyarray[0]=vdxy;
	     vtxdzarray[0]=vdz;
	     phiarray[0]=phi;
	     trkpx[0]=px;
	     trkpy[0]=py;
	     trkpz[0]=pz;
	   }
	   if(ntrk==1){
	     charray[1]=charge;	
	     chi2array[1]=chi2;
	     d0array[1]=d0;
	     dzarray[1]=dz;
	     pidarray[1]=pid; pidv0array[1]=pidV0;	     
	     vtxdxyarray[1]=vdxy;
	     vtxdzarray[1]=vdz;
	     phiarray[1]=phi;
	     trkpx[1]=px;
	     trkpy[1]=py;
	     trkpz[1]=pz;
	   }
	   if(ntrk==2){
	     charray[2]=charge;	
	     chi2array[2]=chi2;
	     d0array[2]=d0;
	     dzarray[2]=dz;
	     pidarray[2]=pid; pidv0array[2]=pidV0;
	     vtxdxyarray[2]=vdxy;
	     vtxdzarray[2]=vdz;
	     phiarray[2]=phi;
	     trkpx[2]=px;
	     trkpy[2]=py;
	     trkpz[2]=pz;
	   }
	   if(ntrk==3){
	     charray[3]=charge;	
	     chi2array[3]=chi2;
	     d0array[3]=d0;
	     dzarray[3]=dz;
	     pidarray[3]=pid; pidv0array[3]=pidV0;
	     vtxdxyarray[3]=vdxy;
	     vtxdzarray[3]=vdz;
	     phiarray[3]=phi;
	     trkpx[3]=px;
	     trkpy[3]=py;
	     trkpz[3]=pz;
	   }

       //...ordering pions and kaons by momentum, index=1 is the highest Pt
       //...we need to include kaons for the selection 11 : one primary and one Vee
	   
       vector<Double_t> piVec = { piA.Pt(), piB.Pt(), piC.Pt(), piD.Pt() };
       //
       sort(piVec.begin(), piVec.end());
       
       //...ordering by Pt and connecting the charges & PID's to the particles...tricky!

       //if(piVec[3]!=0.0 && piVec[3]==piA.Pt()){ pi1 = piA ; npiA=1;
       if(piVec[3]!=0.0 && piVec[3]==piA.Pt()){ pi1 = piA ;
	 atrkpx[0]=trkpx[0];
	 atrkpy[0]=trkpy[0];
	 atrkpz[0]=trkpz[0];
	 arraych[0]=charray[0];
       arraychi2[0]=chi2array[0];
	 arrayd0[0]=d0array[0];
	 arraydz[0]=dzarray[0];
        arraypid[0]=pidarray[0];
        arraypidv0[0]=pidv0array[0];
	arrayphi[0]=phiarray[0];
	arrayvtxdxy[0]=vtxdxyarray[0];
	arrayvtxdz[0]=vtxdzarray[0];}

       //if(piVec[3]!=0.0 && piVec[3]==piB.Pt()){ pi1 = piB ; npiB=1;
       if(piVec[3]!=0.0 && piVec[3]==piB.Pt()){ pi1 = piB ;
	 atrkpx[0]=trkpx[1];
	 atrkpy[0]=trkpy[1];
	 atrkpz[0]=trkpz[1];
	 arraych[0]=charray[1];
       arraychi2[0]=chi2array[1];
	 arrayd0[0]=d0array[1];
	 arraydz[0]=dzarray[1];
        arraypid[0]=pidarray[1];
        arraypidv0[0]=pidv0array[1];
	arrayphi[0]=phiarray[1];
	arrayvtxdxy[0]=vtxdxyarray[1];
	arrayvtxdz[0]=vtxdzarray[1];}
 
       //if(piVec[3]!=0.0 && piVec[3]==piC.Pt()){ pi1 = piC ; npiC=1; 
       if(piVec[3]!=0.0 && piVec[3]==piC.Pt()){ pi1 = piC ; 
	 atrkpx[0]=trkpx[2];
	 atrkpy[0]=trkpy[2];
	 atrkpz[0]=trkpz[2];
	 arraych[0]=charray[2];
       arraychi2[0]=chi2array[2];
	 arrayd0[0]=d0array[2];
	 arraydz[0]=dzarray[2];
        arraypid[0]=pidarray[2];
        arraypidv0[0]=pidv0array[2];
	arrayphi[0]=phiarray[2];
       	arrayvtxdxy[0]=vtxdxyarray[2];
	arrayvtxdz[0]=vtxdzarray[2];}
 
       //if(piVec[3]!=0.0 && piVec[3]==piD.Pt()){ pi1 = piD ; npiD=1;
       if(piVec[3]!=0.0 && piVec[3]==piD.Pt()){ pi1 = piD ;
	 atrkpx[0]=trkpx[3];
	 atrkpy[0]=trkpy[3];
	 atrkpz[0]=trkpz[3];
	 arraych[0]=charray[3];
       arraychi2[0]=chi2array[3];
	 arrayd0[0]=d0array[3];
	 arraydz[0]=dzarray[3];
	arraypid[0]=pidarray[3];
	arraypidv0[0]=pidv0array[3];
	arrayphi[0]=phiarray[3];
	arrayvtxdxy[0]=vtxdxyarray[3];
	arrayvtxdz[0]=vtxdzarray[3];}
 
       //	 
       //if(piVec[2]!=0.0 && piVec[2]==piA.Pt()){ pi2 = piA ; npiA=2;
       if(piVec[2]!=0.0 && piVec[2]==piA.Pt()){ pi2 = piA ;
	 atrkpx[1]=trkpx[0];
	 atrkpy[1]=trkpy[0];
	 atrkpz[1]=trkpz[0];
	 arraych[1]=charray[0];
       arraychi2[1]=chi2array[0];
	 arrayd0[1]=d0array[0];
	 arraydz[1]=dzarray[0];
        arraypid[1]=pidarray[0];
        arraypidv0[1]=pidv0array[0];
	arrayphi[1]=phiarray[0];
       	arrayvtxdxy[1]=vtxdxyarray[0];
	arrayvtxdz[1]=vtxdzarray[0];}

       //if(piVec[2]!=0.0 && piVec[2]==piB.Pt()){ pi2 = piB ; npiB=2; 
       if(piVec[2]!=0.0 && piVec[2]==piB.Pt()){ pi2 = piB ; 
	 atrkpx[1]=trkpx[1];
	 atrkpy[1]=trkpy[1];
	 atrkpz[1]=trkpz[1];
	 arraych[1]=charray[1];
       arraychi2[1]=chi2array[1];
	 arrayd0[1]=d0array[1];
	 arraydz[1]=dzarray[1];
	arraypid[1]=pidarray[1];
	arraypidv0[1]=pidv0array[1];
	arrayphi[1]=phiarray[1];
       	arrayvtxdxy[1]=vtxdxyarray[1];
	arrayvtxdz[1]=vtxdzarray[1];}
	
       //if(piVec[2]!=0.0 && piVec[2]==piC.Pt()){ pi2 = piC ; npiC=2; 
       if(piVec[2]!=0.0 && piVec[2]==piC.Pt()){ pi2 = piC ; 
	 atrkpx[1]=trkpx[2];
	 atrkpy[1]=trkpy[2];
	 atrkpz[1]=trkpz[2];
	 arraych[1]=charray[2];
       arraychi2[1]=chi2array[2];
	 arrayd0[1]=d0array[2];
	 arraydz[1]=dzarray[2];
	arraypid[1]=pidarray[2];
	arraypidv0[1]=pidv0array[2];
	arrayphi[1]=phiarray[2];
       	arrayvtxdxy[1]=vtxdxyarray[2];
	arrayvtxdz[1]=vtxdzarray[2];}
	
       //if(piVec[2]!=0.0 && piVec[2]==piD.Pt()){ pi2 = piD ; npiD=2; 
       if(piVec[2]!=0.0 && piVec[2]==piD.Pt()){ pi2 = piD ; 
	 atrkpx[1]=trkpx[3];
	 atrkpy[1]=trkpy[3];
	 atrkpz[1]=trkpz[3];
	 arraych[1]=charray[3];
       arraychi2[1]=chi2array[3];
	 arrayd0[1]=d0array[3];
	 arraydz[1]=dzarray[3];
	arraypid[1]=pidarray[3];
	arraypidv0[1]=pidv0array[3];
	arrayphi[1]=phiarray[3];
       	arrayvtxdxy[1]=vtxdxyarray[3];
	arrayvtxdz[1]=vtxdzarray[3];}
       
       //	 
       //if(piVec[1]!=0.0 && piVec[1]==piA.Pt()){ pi3 = piA ; npiA=3; 
       if(piVec[1]!=0.0 && piVec[1]==piA.Pt()){ pi3 = piA ; 
	 atrkpx[2]=trkpx[0];
	 atrkpy[2]=trkpy[0];
	 atrkpz[2]=trkpz[0];
	 arraych[2]=charray[0];
       arraychi2[2]=chi2array[0];
	 arrayd0[2]=d0array[0];
	 arraydz[2]=dzarray[0];
        arraypid[2]=pidarray[0];
        arraypidv0[2]=pidv0array[0];
	arrayphi[2]=phiarray[0];
       	arrayvtxdxy[2]=vtxdxyarray[0];
	arrayvtxdz[2]=vtxdzarray[0];}
	
       //if(piVec[1]!=0.0 && piVec[1]==piB.Pt()){ pi3 = piB ; npiB=3;  
       if(piVec[1]!=0.0 && piVec[1]==piB.Pt()){ pi3 = piB ;  
	 atrkpx[2]=trkpx[1];
	 atrkpy[2]=trkpy[1];
	 atrkpz[2]=trkpz[1];
	 arraych[2]=charray[1];
       arraychi2[2]=chi2array[1];
	 arrayd0[2]=d0array[1];
	 arraydz[2]=dzarray[1];
	arraypid[2]=pidarray[1];
	arraypidv0[2]=pidv0array[1];
	arrayphi[2]=phiarray[1];
       	arrayvtxdxy[2]=vtxdxyarray[1];
	arrayvtxdz[2]=vtxdzarray[1];}

       //if(piVec[1]!=0.0 && piVec[1]==piC.Pt()){ pi3 = piC ; npiC=3; 
       if(piVec[1]!=0.0 && piVec[1]==piC.Pt()){ pi3 = piC ; 
	 atrkpx[2]=trkpx[2];
	 atrkpy[2]=trkpy[2];
	 atrkpz[2]=trkpz[2];
	 arraych[2]=charray[2];
       arraychi2[2]=chi2array[2];
	 arrayd0[2]=d0array[2];
	 arraydz[2]=dzarray[2];
	arraypid[2]=pidarray[2];
	arraypidv0[2]=pidv0array[2];
	arrayphi[2]=phiarray[2];
       	arrayvtxdxy[2]=vtxdxyarray[2];
	arrayvtxdz[2]=vtxdzarray[2];}
	
       //if(piVec[1]!=0.0 && piVec[1]==piD.Pt()){ pi3 = piD ; npiD=3;  
       if(piVec[1]!=0.0 && piVec[1]==piD.Pt()){ pi3 = piD ;  
	 atrkpx[2]=trkpx[3];
	 atrkpy[2]=trkpy[3];
	 atrkpz[2]=trkpz[3];
	 arraych[2]=charray[3];
       arraychi2[2]=chi2array[3];
	 arrayd0[2]=d0array[3];
	 arraydz[2]=dzarray[3];
	arraypid[2]=pidarray[3];
	arraypidv0[2]=pidv0array[3];
	arrayphi[2]=phiarray[3];
       	arrayvtxdxy[2]=vtxdxyarray[3];
	arrayvtxdz[2]=vtxdzarray[3];}
	
       //
       //if(piVec[0]!=0.0 && piVec[0]==piA.Pt()){ pi4 = piA ; npiA=4;  
       if(piVec[0]!=0.0 && piVec[0]==piA.Pt()){ pi4 = piA ;  
	 atrkpx[3]=trkpx[0];
	 atrkpy[3]=trkpy[0];
	 atrkpz[3]=trkpz[0];
	 arraych[3]=charray[0];
       arraychi2[3]=chi2array[0];
	 arrayd0[3]=d0array[0];
	 arraydz[3]=dzarray[0];
	arraypid[3]=pidarray[0];
	arraypidv0[3]=pidv0array[0];
	arrayphi[3]=phiarray[0];
       	arrayvtxdxy[3]=vtxdxyarray[0];
	arrayvtxdz[3]=vtxdzarray[0];}
	
       //if(piVec[0]!=0.0 && piVec[0]==piB.Pt()){ pi4 = piB ; npiB=4; 
       if(piVec[0]!=0.0 && piVec[0]==piB.Pt()){ pi4 = piB ; 
	 atrkpx[3]=trkpx[1];
	 atrkpy[3]=trkpy[1];
	 atrkpz[3]=trkpz[1];
	 arraych[3]=charray[1];
       arraychi2[3]=chi2array[1];
	 arrayd0[3]=d0array[1];
	 arraydz[3]=dzarray[1];
	arraypid[3]=pidarray[1];
	arraypidv0[3]=pidv0array[1];
        arrayphi[3]=phiarray[1];
       	arrayvtxdxy[3]=vtxdxyarray[1];
	arrayvtxdz[3]=vtxdzarray[1];}
	
       //if(piVec[0]!=0.0 && piVec[0]==piC.Pt()){ pi4 = piC ; npiC=4; 
       if(piVec[0]!=0.0 && piVec[0]==piC.Pt()){ pi4 = piC ; 
	 atrkpx[3]=trkpx[2];
	 atrkpy[3]=trkpy[2];
	 atrkpz[3]=trkpz[2];
	 arraych[3]=charray[2];
       arraychi2[3]=chi2array[2];
	 arrayd0[3]=d0array[2];
	 arraydz[3]=dzarray[2];
	arraypid[3]=pidarray[2];
	arraypidv0[3]=pidv0array[2];
	arrayphi[3]=phiarray[2];
       	arrayvtxdxy[3]=vtxdxyarray[2];
	arrayvtxdz[3]=vtxdzarray[2];}
	
       //if(piVec[0]!=0.0 && piVec[0]==piD.Pt()){ pi4 = piD ; npiD=4; 
       if(piVec[0]!=0.0 && piVec[0]==piD.Pt()){ pi4 = piD ; 
	 atrkpx[3]=trkpx[3];
	 atrkpy[3]=trkpy[3];
	 atrkpz[3]=trkpz[3];
	 arraych[3]=charray[3];       
       arraychi2[3]=chi2array[3];       
	 arrayd0[3]=d0array[3];       
	 arraydz[3]=dzarray[3];       
	arraypid[3]=pidarray[3];
	arraypidv0[3]=pidv0array[3];
	arrayphi[3]=phiarray[3];
       	arrayvtxdxy[3]=vtxdxyarray[3];
	arrayvtxdz[3]=vtxdzarray[3];}


       //-----------------------

       // kaons
       double eneK=TMath::Sqrt(pt*pt+pz*pz+m_k*m_k);
       TLorentzVector trk_lorentzK(itTrack->px(),itTrack->py(),itTrack->pz(),eneK);
       kkkkRec += trk_lorentzK; 

	   if(ntrk==0) kA = trk_lorentzK;
	   if(ntrk==1) kB = trk_lorentzK;
	   if(ntrk==2) kC = trk_lorentzK;
	   if(ntrk==3) kD = trk_lorentzK;

       vector<Double_t> kVec = { kA.Pt(), kB.Pt(), kC.Pt(), kD.Pt() };
       //
       sort(kVec.begin(), kVec.end());
        //
       //...checking the consistency between pion- and kaon-track pt orders
       //if(kVec[3]!=0.0 && kVec[3]==kA.Pt()){ k1 = kA ; nkA=1;}
       //if(kVec[3]!=0.0 && kVec[3]==kB.Pt()){ k1 = kB ; nkB=1;}
       //if(kVec[3]!=0.0 && kVec[3]==kC.Pt()){ k1 = kC ; nkC=1;}
       //if(kVec[3]!=0.0 && kVec[3]==kD.Pt()){ k1 = kD ; nkD=1;}
       //	 
       //if(kVec[2]!=0.0 && kVec[2]==kA.Pt()){ k2 = kA ; nkA=2;}
       //if(kVec[2]!=0.0 && kVec[2]==kB.Pt()){ k2 = kB ; nkB=2;}
       //if(kVec[2]!=0.0 && kVec[2]==kC.Pt()){ k2 = kC ; nkC=2;}
       //if(kVec[2]!=0.0 && kVec[2]==kD.Pt()){ k2 = kD ; nkD=2;}
       //	 
       //if(kVec[1]!=0.0 && kVec[1]==kA.Pt()){ k3 = kA ; nkA=3;}
       //if(kVec[1]!=0.0 && kVec[1]==kB.Pt()){ k3 = kB ; nkB=3;}
       //if(kVec[1]!=0.0 && kVec[1]==kC.Pt()){ k3 = kC ; nkC=3;}
       //if(kVec[1]!=0.0 && kVec[1]==kD.Pt()){ k3 = kD ; nkD=3;}
       //	  
       //if(kVec[0]!=0.0 && kVec[0]==kA.Pt()){ k4 = kA ; nkA=4;}
       //if(kVec[0]!=0.0 && kVec[0]==kB.Pt()){ k4 = kB ; nkB=4;}
       //if(kVec[0]!=0.0 && kVec[0]==kC.Pt()){ k4 = kC ; nkC=4;}
       //if(kVec[0]!=0.0 && kVec[0]==kD.Pt()){ k4 = kD ; nkD=4;}       

       if(kVec[3]!=0.0 && kVec[3]==kA.Pt()){ k1 = kA ;}
       if(kVec[3]!=0.0 && kVec[3]==kB.Pt()){ k1 = kB ;}
       if(kVec[3]!=0.0 && kVec[3]==kC.Pt()){ k1 = kC ;}
       if(kVec[3]!=0.0 && kVec[3]==kD.Pt()){ k1 = kD ;}
       //	 
       if(kVec[2]!=0.0 && kVec[2]==kA.Pt()){ k2 = kA ;}
       if(kVec[2]!=0.0 && kVec[2]==kB.Pt()){ k2 = kB ;}
       if(kVec[2]!=0.0 && kVec[2]==kC.Pt()){ k2 = kC ;}
       if(kVec[2]!=0.0 && kVec[2]==kD.Pt()){ k2 = kD ;}
       //	 
       if(kVec[1]!=0.0 && kVec[1]==kA.Pt()){ k3 = kA ;}
       if(kVec[1]!=0.0 && kVec[1]==kB.Pt()){ k3 = kB ;}
       if(kVec[1]!=0.0 && kVec[1]==kC.Pt()){ k3 = kC ;}
       if(kVec[1]!=0.0 && kVec[1]==kD.Pt()){ k3 = kD ;}
       //	  
       if(kVec[0]!=0.0 && kVec[0]==kA.Pt()){ k4 = kA ;}
       if(kVec[0]!=0.0 && kVec[0]==kB.Pt()){ k4 = kB ;}
       if(kVec[0]!=0.0 && kVec[0]==kC.Pt()){ k4 = kC ;}
       if(kVec[0]!=0.0 && kVec[0]==kD.Pt()){ k4 = kD ;}        
       
      //--------------------------------------
      //--------------------------------------
      // end of my 4trk definitions

      ntrk++; //pixel tracks
     } //...end of npixelhits>0

    ntrk0++; //all tracks

  } //...end of tracks

      //...checking the consistency between pion- and kaon-track pt orders
      //if(ntrk==4 && totcharge==0){
      //std::cout << "**********************************************************" << std::endl;
      //std::cout << "npiA=" << npiA << " npiB=" << npiB << " npiC=" << npiC << " npiD=" << npiD << std::endl;
      //std::cout << " nkA=" <<  nkA << "  nkB=" <<  nkB << "  nkC=" <<  nkC << "  nkD=" <<  nkD << std::endl;
      //}
  
  histosTH1F["hntrk0"]->Fill(ntrk0);
  histosTH1F["hntrk"]->Fill(ntrk);

       if(ntrk==4 && totcharge==0){
	 histosTH1F["hntrk4q0"]->Fill(ntrk);
       }    

       
       double pi1pt = pi1.Pt();
       double pi2pt = pi2.Pt();
       double pi3pt = pi3.Pt();
       double pi4pt = pi4.Pt();

       double pi1eta = pi1.Eta();
       double pi2eta = pi2.Eta();
       double pi3eta = pi3.Eta();
       double pi4eta = pi4.Eta();

       histosTH1F["heta4pi"]->Fill(pi1eta);
       histosTH1F["heta4pi"]->Fill(pi2eta);
       histosTH1F["heta4pi"]->Fill(pi3eta);
       histosTH1F["heta4pi"]->Fill(pi4eta);

           // my 4-vectors
           pi1pi2Rec = pi1 + pi2;
           pi3pi4Rec = pi3 + pi4;
           pi1pi3Rec = pi1 + pi3;
           pi2pi4Rec = pi2 + pi4;      
	   pi1pi4Rec = pi1 + pi4;
	   pi2pi3Rec = pi2 + pi3;
           //
           k1k2Rec = k1 + k2;
           k3k4Rec = k3 + k4;
           k1k3Rec = k1 + k3;
           k2k4Rec = k2 + k4;      
           k1k4Rec = k1 + k4;
           k2k3Rec = k2 + k3;      

	   // ...checking out
	   pi1234Rec = pi1pi2Rec +  pi3pi4Rec;
	   pi1324Rec = pi1pi3Rec +  pi2pi4Rec;
	   pi1423Rec = pi1pi4Rec +  pi2pi3Rec;


  
	   // Armenteros-Podolanski method
	   //............1234
	   double vee12px = pi1pi2Rec.Px();
	   double vee12py = pi1pi2Rec.Py();
	   double vee12pz = pi1pi2Rec.Pz();
	     //
	   double vee34px = pi3pi4Rec.Px();
	   double vee34py = pi3pi4Rec.Py();
           double vee34pz = pi3pi4Rec.Pz();
	   //............1324
	   double vee13px = pi1pi3Rec.Px();
	   double vee13py = pi1pi3Rec.Py();
	   double vee13pz = pi1pi3Rec.Pz();
	     //
	   double vee24px = pi2pi4Rec.Px();
	   double vee24py = pi2pi4Rec.Py();
           double vee24pz = pi2pi4Rec.Pz();
	   //............1423
	   double vee14px = pi1pi4Rec.Px();
	   double vee14py = pi1pi4Rec.Py();
	   double vee14pz = pi1pi4Rec.Pz();
	     //
	   double vee23px = pi2pi3Rec.Px();
	   double vee23py = pi2pi3Rec.Py();
           double vee23pz = pi2pi3Rec.Pz();

	   
	   // impact parameter
	   //double vee12x = pi1pi2Rec.vx();
	   //double vee12y = pi1pi2Rec.vy();
	   //double vee12z = pi1pi2Rec.vz();
	   //
	   //double rimpacvee12 = TMath::Sqrt( vee12x*vee12x + vee12y*vee12y + vee12z*vee12z );
	   
	   //
	   //std::cout << " ooooooooooooooooo " << std::endl;
	   //std::cout << " vee12x = " << vee12x << std::endl;
	   //std::cout << " vee12y = " << vee12y << std::endl;
	   //std::cout << " vee12z = " << vee12z << std::endl;
	   //std::cout << " vee12 = " << vee12 << std::endl;

	   
	   //.........................................
	   double pi1px = pi1.Px();
	   double pi2px = pi2.Px();
	   double pi3px = pi3.Px();
	   double pi4px = pi4.Px();
	   //
	   double pi1py = pi1.Py();
	   double pi2py = pi2.Py();
	   double pi3py = pi3.Py();
	   double pi4py = pi4.Py();
	   //
	   double pi1pz = pi1.Pz();
	   double pi2pz = pi2.Pz();
	   double pi3pz = pi3.Pz();
	   double pi4pz = pi4.Pz();
	   //.........................................

	   
	     double Ppi1 =  TMath::Sqrt( pi1px*pi1px + pi1py*pi1py + pi1pz*pi1pz );
	     double Ppi2 =  TMath::Sqrt( pi2px*pi2px + pi2py*pi2py + pi2pz*pi2pz );
	     double Ppi3 =  TMath::Sqrt( pi3px*pi3px + pi3py*pi3py + pi3pz*pi3pz );
	     double Ppi4 =  TMath::Sqrt( pi4px*pi4px + pi4py*pi4py + pi4pz*pi4pz );

	     
	     //...V0 ...1234
	     double Pvee12 = TMath::Sqrt( vee12px*vee12px + vee12py*vee12py + vee12pz*vee12pz );
	     double Pvee34 = TMath::Sqrt( vee34px*vee34px + vee34py*vee34py + vee34pz*vee34pz );

	     //...V0 ...1324
	     double Pvee13 = TMath::Sqrt( vee13px*vee13px + vee13py*vee13py + vee13pz*vee13pz );
	     double Pvee24 = TMath::Sqrt( vee24px*vee24px + vee24py*vee24py + vee24pz*vee24pz );

	     //...V0 ...1423
	     double Pvee14 = TMath::Sqrt( vee14px*vee14px + vee14py*vee14py + vee14pz*vee14pz );
	     double Pvee23 = TMath::Sqrt( vee23px*vee23px + vee23py*vee23py + vee23pz*vee23pz );
	     
	     

     //...combining pions and kaons for the event selection type = 11 (one primary & one Vee)
           //
	   /* 
     ...first combining, then select the Q_pairs=0
     pi1pi2 pi3k4
     pi1pi3 pi2k4
     pi2pi3 pi1k4

     pi1pi2 k3pi4
     pi1pi4 k3pi2
     pi2pi4 k3pi1

     pi1k2 pi3pi4
     pi3k2 pi1pi4
     pi4k2 pi1pi3

     k1pi2 pi3pi4
     k1pi3 pi2pi4
     k1pi4 pi2pi3
 	   */
	   
	   //...commented out means already defined
	   //
	   //pi1pi2Rec
	   pi3k4Rec  = pi3 + k4;
	   //pi1pi3Rec
	   pi2k4Rec  = pi2 + k4; 
	   //pi2pi3Rec
	   pi1k4Rec  = pi1 + k4;
	   //
	   //
	   //pi1pi2Rec
	   k3pi4Rec  = k3  + pi4;
	   //pi1pi4Rec
	   k3pi2Rec  = k3  + pi2;
	   //pi2pi4Rec
	   k3pi1Rec  = k3  + pi1;
	   //
	   //
	   pi1k2Rec  = pi1 + k2;
	   //pi3pi4Rec
	   pi3k2Rec  = pi3 + k2;
	   //pi1pi4Rec
	   pi4k2Rec  = pi4 + k2;
	   //pi1pi3Rec	   
	   //
	   //
	   k1pi2Rec  = k1  + pi2;
	   //pi3pi4Rec
	   k1pi3Rec  = k1  + pi3;
	   //pi2pi4Rec
	   k1pi4Rec  = k1  + pi4;
	   //pi2pi3Rec

	   /*
	   pipiRec = pi1pi2Rec +
	             pi3pi4Rec +
	             pi1pi3Rec +
	             pi2pi4Rec +
	             pi1pi4Rec +
	             pi2pi3Rec ;
	   */

	   //...type:11  ...important for (K0sK+pi-) or (K0sK-pi+)
           k1pi2pi3pi4Rec  = k1  + pi2 + pi3 + pi4;
           pi1k2pi3pi4Rec  = pi1 + k2  + pi3 + pi4;
           pi1pi2k3pi4Rec  = pi1 + pi2 + k3  + pi4;
           pi1pi2pi3k4Rec  = pi1 + pi2 + pi3 + k4 ;


	   
	   //...K*Kbar* channel
	   //1234
           k1pi2k3pi4Rec  = k1  + pi2 + k3 + pi4;
           k1pi2pi3k4Rec  = k1  + pi2 + pi3 + k4;
           pi1k2k3pi4Rec  = pi1  + k2 + k3 + pi4;
           pi1k2pi3k4Rec  = pi1  + k2 + pi3 + k4;
	   //1324
           k1pi3k2pi4Rec  = k1  + pi3 + k2 + pi4;
           k1pi3pi2k4Rec  = k1  + pi3 + pi2 + k4;
           pi1k3k2pi4Rec  = pi1  + k3 + k2 + pi4;
           pi1k3pi2k4Rec  = pi1  + k3 + pi2 + k4;
	   //1423
           k1pi4k2pi3Rec  = k1  + pi4 + k2 + pi3;
           k1pi4pi2k3Rec  = k1  + pi4 + pi2 + k3;
           pi1k4k2pi3Rec  = pi1  + k4 + k2 + pi3;
           pi1k4pi2k3Rec  = pi1  + k4 + pi2 + k3;
 	   
	 
	   //1234
	   //           k1pi2k3pi4Rec  = k1  + pi2 + k3 + pi4;
	   //...k1pi2Rec  = k1  + pi2;
           //...k3pi4Rec  = k3  + pi4;
	   //           k1pi2pi3k4Rec  = k1  + pi2 + pi3 + k4;
	   // k1pi2Rec  = k1  + pi2;
           //...pi3k4Rec  = pi3 + k4;
	   //           pi1k2k3pi4Rec  = pi1  + k2 + k3 + pi4;
           //...pi1k2Rec  = pi1 + k2;
           // k3pi4Rec  = k3  + pi4;
	   //           pi1k2pi3k4Rec  = pi1  + k2 + pi3 + k4;
           // pi1k2Rec  = pi1 + k2;
           // pi3k4Rec  = pi3 + k4;


	   //1324
	   //           k1pi3k2pi4Rec  = k1  + pi3 + k2 + pi4;
           //...k1pi3Rec  = k1  + pi3;
           k2pi4Rec  = k2  + pi4;   // * * * 
	   //           k1pi3pi2k4Rec  = k1  + pi3 + pi2 + k4;
           // k1pi3Rec  = k1  + pi3;
           //...pi2k4Rec  = pi2 + k4;
	   //           pi1k3k2pi4Rec  = pi1  + k3 + k2 + pi4;
           //...pi1k3Rec  = pi1 + k3;
           // k2pi4Rec  = k2  + pi4;
	   //           pi1k3pi2k4Rec  = pi1  + k3 + pi2 + k4;
           // pi1k3Rec  = pi1 + k3;
           // pi2k4Rec  = pi2 + k4;


	   //1423
	   //           k1pi4k2pi3Rec  = k1  + pi4 + k2 + pi3;
           //...k1pi4Rec  = k1  + pi4;
           k2pi3Rec  = k2  + pi3;   // * * * 
	   //           k1pi4pi2k3Rec  = k1  + pi4 + pi2 + k3;
           // k1pi4Rec  = k1  + pi4;
           pi2k3Rec  = pi2 + k3;    // * * *
	   //           pi1k4k2pi3Rec  = pi1  + k4 + k2 + pi3;
           //...pi1k4Rec  = pi1 + k4;
           // k2pi3Rec  = k2  + pi3;
	   //           pi1k4pi2k3Rec  = pi1  + k4 + pi2 + k3;
           // pi1k4Rec  = pi1 + k4;
           // pi2k3Rec  = pi2 + k3;

	   
	   //...K*Kbar*
	   pi1k3Rec  = pi1  + k3;

	   
	   //...reseting to the original definition with new order

         charray[0] =   arraych[0] ;// charge;	
       chi2array[0] = arraychi2[0] ;// chi2;
         d0array[0] =   arrayd0[0] ;// d0;
         dzarray[0] =   arraydz[0] ;// dz;
        pidarray[0] =  arraypid[0] ;// pid;
      pidv0array[0] =  arraypidv0[0];//pidV0;
	phiarray[0] =  arrayphi[0] ;// phi;
     vtxdxyarray[0] = arrayvtxdxy[0];//vtxdxy;
      vtxdzarray[0] =  arrayvtxdz[0];//vtxdz;
	//
	 charray[1] =   arraych[1] ;// charge;	
       chi2array[1] = arraychi2[1] ;// chi2;
         d0array[1] =   arrayd0[1] ;// d0;
         dzarray[1] =   arraydz[1] ;// dz;
        pidarray[1] =  arraypid[1] ;// pid;
      pidv0array[1] =  arraypidv0[1];//pidV0;
	phiarray[1] =  arrayphi[1] ;// phi;
     vtxdxyarray[1] = arrayvtxdxy[1];//vtxdxy;
      vtxdzarray[1] =  arrayvtxdz[1];//vtxdz;
	//
         charray[2] =   arraych[2] ;// charge;	
       chi2array[2] = arraychi2[2] ;// chi2;
	 d0array[2] =   arrayd0[2] ;// d0;
	 dzarray[2] =   arraydz[2] ;// dz;
	pidarray[2] =  arraypid[2] ;// pid;
      pidv0array[2] =  arraypidv0[2];//pidV0;
	phiarray[2] =  arrayphi[2] ;// phi;
     vtxdxyarray[2] = arrayvtxdxy[2];//vtxdxy;
      vtxdzarray[2] =  arrayvtxdz[2];//vtxdz;
	//
         charray[3] =   arraych[3] ;// charge;	
       chi2array[3] = arraychi2[3] ;// chi2;
	 d0array[3] =   arrayd0[3] ;// d0;
	 dzarray[3] =   arraydz[3] ;// dz;
	pidarray[3] =  arraypid[3] ;// pid;
      pidv0array[3] =  arraypidv0[3];//pidV0;
	phiarray[3] =  arrayphi[3] ;// phi;
     vtxdxyarray[3] = arrayvtxdxy[3];//vtxdxy;
      vtxdzarray[3] =  arrayvtxdz[3];//vtxdz;
	   
      trkpx[0] = atrkpx[0];
      trkpy[0] = atrkpy[0];
      trkpz[0] = atrkpz[0];
      //
      trkpx[1] = atrkpx[1];
      trkpy[1] = atrkpy[1];
      trkpz[1] = atrkpz[1];
      //
      trkpx[2] = atrkpx[2];
      trkpy[2] = atrkpy[2];
      trkpz[2] = atrkpz[2];
      //
      trkpx[3] = atrkpx[3];
      trkpy[3] = atrkpy[3];
      trkpz[3] = atrkpz[3];

      
      /* do we need a similar one here ?
       if(ntrk==0){
	 int nclusters=  sipixelcluster_coll->size();
	 int nclusters2= sistripcluster_coll->nStripClusters;
	 
	 histosTH1F["hnclusters"]->Fill(nclusters);
	 histosTH1F["hnclusters2"]->Fill(nclusters2);
       }
      */

      //...finding the tracks connected to the primary vertex using impact parameter...type:11 only
      bool isTrack1 = false ;
      bool isTrack2 = false ;
      bool isTrack3 = false ;
      bool isTrack4 = false ;

      /*
        vector<Double_t> vdxyVec = { TMath::Abs(vtxdxyarray[0]), TMath::Abs(vtxdxyarray[1]),
				     TMath::Abs(vtxdxyarray[2]), TMath::Abs(vtxdxyarray[3]) };
        // ...better to use 3D impact parameter...REVIEW!
        sort(vdxyVec.begin(), vdxyVec.end());
        //	  
        if(vdxyVec[0]!=0.0 && vdxyVec[0]==TMath::Abs(vtxdxyarray[0])){ isTrack1 = true ; }
        if(vdxyVec[0]!=0.0 && vdxyVec[0]==TMath::Abs(vtxdxyarray[1])){ isTrack2 = true ; }
        if(vdxyVec[0]!=0.0 && vdxyVec[0]==TMath::Abs(vtxdxyarray[2])){ isTrack3 = true ; }
        if(vdxyVec[0]!=0.0 && vdxyVec[0]==TMath::Abs(vtxdxyarray[3])){ isTrack4 = true ; }       
	//	
	if(vdxyVec[1]!=0.0 && vdxyVec[1]==TMath::Abs(vtxdxyarray[0])){ isTrack1 = true ; }
        if(vdxyVec[1]!=0.0 && vdxyVec[1]==TMath::Abs(vtxdxyarray[1])){ isTrack2 = true ; }
        if(vdxyVec[1]!=0.0 && vdxyVec[1]==TMath::Abs(vtxdxyarray[2])){ isTrack3 = true ; }
        if(vdxyVec[1]!=0.0 && vdxyVec[1]==TMath::Abs(vtxdxyarray[3])){ isTrack4 = true ; }
      */      

       // ...transverse impact parameter distribution d0 (dxy)
       d01 = d0array[0];
       d02 = d0array[1];
       d03 = d0array[2];
       d04 = d0array[3];
       // ...longitudinal impact parameter distribution dz
       dz1 = dzarray[0];
       dz2 = dzarray[1];
       dz3 = dzarray[2];
       dz4 = dzarray[3];
       
       if(ntrk==4 && totcharge==0){
	 histosTH1F["hd01"]->Fill(d01);
	 histosTH1F["hd02"]->Fill(d02);
	 histosTH1F["hd03"]->Fill(d03);
	 histosTH1F["hd04"]->Fill(d04);
	 //
	 histosTH1F["hdz1"]->Fill(dz1);
	 histosTH1F["hdz2"]->Fill(dz2);
	 histosTH1F["hdz3"]->Fill(dz3);
	 histosTH1F["hdz4"]->Fill(dz4);
       }    
       // ...transverse impact parameter distribution vtxdxy
       vtxdxy1 = vtxdxyarray[0];
       vtxdxy2 = vtxdxyarray[1];
       vtxdxy3 = vtxdxyarray[2];
       vtxdxy4 = vtxdxyarray[3];
       // ...longitudinal impact parameter distribution vtxdz
       vtxdz1 = vtxdzarray[0];
       vtxdz2 = vtxdzarray[1];
       vtxdz3 = vtxdzarray[2];
       vtxdz4 = vtxdzarray[3];

       if(ntrk==4 && totcharge==0){
	 histosTH1F["hvtxdxy1"]->Fill(vtxdxy1);
	 histosTH1F["hvtxdxy2"]->Fill(vtxdxy2);
	 histosTH1F["hvtxdxy3"]->Fill(vtxdxy3);
	 histosTH1F["hvtxdxy4"]->Fill(vtxdxy4);
         //	 
	 histosTH1F["hvtxdz1"]->Fill(vtxdz1);
	 histosTH1F["hvtxdz2"]->Fill(vtxdz2);
	 histosTH1F["hvtxdz3"]->Fill(vtxdz3);
	 histosTH1F["hvtxdz4"]->Fill(vtxdz4);
       }    

  // ...Luiz
  //...3D impact parameter
  double rimpac1 = TMath::Sqrt( vtxdxy1*vtxdxy1 + vtxdz1*vtxdz1 );
  double rimpac2 = TMath::Sqrt( vtxdxy2*vtxdxy2 + vtxdz2*vtxdz2 );
  double rimpac3 = TMath::Sqrt( vtxdxy3*vtxdxy3 + vtxdz3*vtxdz3 );
  double rimpac4 = TMath::Sqrt( vtxdxy4*vtxdxy4 + vtxdz4*vtxdz4 );

        vector<Double_t> rimpacVec = { rimpac1, rimpac2, rimpac3, rimpac4 };

        sort(rimpacVec.begin(), rimpacVec.end());
        //
	//...type:11 
	//
	//...prompt tracks have smaller impact parameters w.r.t. pv 
        if(rimpacVec[0]!=0.0 && rimpacVec[0]==rimpac1){ isTrack1 = true ; }
        if(rimpacVec[0]!=0.0 && rimpacVec[0]==rimpac2){ isTrack2 = true ; }
        if(rimpacVec[0]!=0.0 && rimpacVec[0]==rimpac3){ isTrack3 = true ; }
        if(rimpacVec[0]!=0.0 && rimpacVec[0]==rimpac4){ isTrack4 = true ; }       
	//	
	if(rimpacVec[1]!=0.0 && rimpacVec[1]==rimpac1){ isTrack1 = true ; }
        if(rimpacVec[1]!=0.0 && rimpacVec[1]==rimpac2){ isTrack2 = true ; }
        if(rimpacVec[1]!=0.0 && rimpacVec[1]==rimpac3){ isTrack3 = true ; }
        if(rimpacVec[1]!=0.0 && rimpacVec[1]==rimpac4){ isTrack4 = true ; }
     
  //----------------------------------------------------------------------
  // VERTEX
  
  int nvtx=0;
  int ntrkvtx = 0;
  int ntrkvtxU = 0;
  
  for(VertexCollection::const_iterator itVtx = vertices->begin();itVtx != vertices->end();++itVtx) {
    int vtxisfake = itVtx->isFake();
    if(vtxisfake==0) nvtx++;    
    else continue;

    //...Luiz 
    ntrkvtx = itVtx->nTracks();
    ntrkvtxU = itVtx->tracksSize();
    //itVtx->Print();
    //std::cout << " ntrkvtx = " << ntrkvtx << std::endl;
  }
  
  histosTH1F["hnvtx"]->Fill(nvtx);

  if(nvtx==1) histosTH1F["hntrkvtx"]->Fill(ntrkvtx);
  if(nvtx==1) histosTH1F["hntrkvtxU"]->Fill(ntrkvtxU);
       //...Luiz
       if(nvtx==0) histosTH1F["hntrkvtx0"]->Fill(ntrkvtx);
       if(nvtx==2) histosTH1F["hntrkvtx2"]->Fill(ntrkvtx);
       if(nvtx==3) histosTH1F["hntrkvtx3"]->Fill(ntrkvtx);
       if(nvtx==4) histosTH1F["hntrkvtx4"]->Fill(ntrkvtx);
       
       
  // considering all of nvtx ...Luiz
  // commented out
  //if(nvtx!=1) return;

  //...Luiz
  int isfake = vertices->begin()->isFake();
  
  double xvtx = vertices->begin()->x();
  double yvtx = vertices->begin()->y();
  double zvtx = vertices->begin()->z();
  double chi2vtx = vertices->begin()->normalizedChi2();
  double ndofvtx = vertices->begin()->ndof();
  
  histosTH1F["hvtxx"]->Fill(xvtx);
  histosTH1F["hvtxy"]->Fill(yvtx);
  histosTH1F["hvtxz"]->Fill(zvtx);
  histosTH1F["hvtxchi2"]->Fill(chi2vtx);

  /*
  // ...Luiz     ...this is interesting! primary vertex position & impact parameters
  math::XYZPoint pv2(xvtx,yvtx,zvtx);
  double pv2dxy = dxy(pv2);   ??
  double pv2dz  = dz(pv2);    ??
  histosTH1F["hpv2dxy"]->Fill(pv2dxy);
  histosTH1F["hpv2dz"]->Fill(pv2dz);
  */
  
  //----------------------------------------------------------------------
  // V0
  /*
  //...Kshort collection...Luiz
  bool isKshort = false;
  int nks=0;

  for(VertexCompositeCandidateCollection::const_iterator it_ks = kshorts->begin() ; it_ks != kshorts->end() ; ++it_ks){

	 nks++;
	 isKshort = nks;
	 double ksvertexx = it_ks->vx();
	 double ksvertexy = it_ks->vy();
	 double ksvertexz = it_ks->vz();
	 double kspt = it_ks->pt();
	 double kseta = it_ks->eta();
	 double ksphi = it_ks->phi();
	 double ksmass = it_ks->mass();
	 // ...do not mix with primary vertex !!!
	 //	 double ksradius = TMath::Sqrt((ksvertexx-xvtx)*(ksvertexx-xvtx)+(ksvertexy-yvtx)*(ksvertexy-yvtx));
	 //	 double ks3Dradius = TMath::Sqrt((ksvertexx-xvtx)*(ksvertexx-xvtx)+(ksvertexy-yvtx)*(ksvertexy-yvtx)+(ksvertexz-zvtx)*(ksvertexz-zvtx));
	 double ksradius = TMath::Sqrt((ksvertexx)*(ksvertexx)+(ksvertexy)*(ksvertexy));
	 double ks3Dradius = TMath::Sqrt((ksvertexx)*(ksvertexx)+(ksvertexy)*(ksvertexy)+(ksvertexz)*(ksvertexz));
	 double ksenergy = TMath::Sqrt(kspt*kspt+0.4976*0.4976);
	 double gammalorentzks = ksenergy/0.4976;
	 double kslifetime = ksradius/gammalorentzks;  
	 double ks3Dlifetime = ks3Dradius/gammalorentzks;
	 
	 histosTH1F["hkspt"]->Fill(kspt);
	 histosTH1F["hkseta"]->Fill(kseta);
	 histosTH1F["hksphi"]->Fill(ksphi);
	 histosTH1F["hksmass"]->Fill(ksmass);
	 //
	 if(nks == 1){histosTH1F["hksmassv1"]->Fill(ksmass);
                      histosTH1F["hksradiusv1"]->Fill(ksradius);
                      histosTH1F["hkslifetimev1"]->Fill(kslifetime);
                      histosTH1F["hks3Dradiusv1"]->Fill(ks3Dradius);
                      histosTH1F["hks3Dlifetimev1"]->Fill(ks3Dlifetime);
	 }
	 if(nks == 2){histosTH1F["hksmassv2"]->Fill(ksmass);
                      histosTH1F["hksradiusv2"]->Fill(ksradius);
                      histosTH1F["hkslifetimev2"]->Fill(kslifetime);
                      histosTH1F["hks3Dradiusv2"]->Fill(ks3Dradius);
                      histosTH1F["hks3Dlifetimev2"]->Fill(ks3Dlifetime);
	 }
	 //
	 histosTH1F["hksvertexx"]->Fill(ksvertexx);
	 histosTH1F["hksvertexy"]->Fill(ksvertexy);
	 histosTH1F["hksvertexz"]->Fill(ksvertexz);
	 histosTH1F["hksradius"]->Fill(ksradius);
	 histosTH1F["hks3Dradius"]->Fill(ks3Dradius);
	 histosTH1F["hkslifetime"]->Fill(kslifetime);
	 histosTH1F["hks3Dlifetime"]->Fill(ks3Dlifetime);
	 histosTH2F["h2dimksxy"]->Fill(ksvertexx,ksvertexy);
	 histosTH2F["h2dimksxz"]->Fill(ksvertexx,ksvertexz);
	 histosTH2F["h2dimksyz"]->Fill(ksvertexy,ksvertexz);
*/	 
	 /*
    	 std::cout << " * * * "<< std::endl;
	 std::cout << " nks = " << nks << std::endl;
	 std::cout << " ksvertexx = " << ksvertexx << std::endl;
	 std::cout << " ksvertexy = " << ksvertexy << std::endl;
	 std::cout << " ksvertexz = " << ksvertexz << std::endl;
	 std::cout << " xvtx = " << xvtx << std::endl;
	 std::cout << " yvtx = " << yvtx << std::endl;
	 std::cout << " zvtx = " << zvtx << std::endl;
	 std::cout << " ksmass = " << ksmass << std::endl;
	 std::cout << " kspt = " << kspt << std::endl;
	 std::cout << " ksradius = " << ksradius << std::endl;
	 */
/*
	 //it_ks->Print();
}
       //...end of Kshort
       histosTH1F["hnks"]->Fill(nks);
       histosTH2F["hntrknks"]->Fill(ntrk,nks);
       histosTH2F["hnvtxnks"]->Fill(nvtx,nks);
       histosTH2F["hntrknvtx"]->Fill(ntrk,nvtx);
       //std::cout << " --------------------------- " << std::endl;
       //std::cout << " nks  = " << nks << std::endl;
       //std::cout << " ntrk = " << ntrk << std::endl;
       //std::cout << " nvtx = " << nvtx << std::endl;
       //std::cout << " isKshort = " << isKshort << std::endl;
       //std::cout << " --------------------------- " << std::endl;
*/
/*       
   //...Lambda collection...Luiz
   ////bool isLambda = false;
   //std::cout << "isLambda ...boo! " << isLambda << std::endl;
   int nlam=0;

   for(VertexCompositeCandidateCollection::const_iterator it_lam = lambdas->begin() ; it_lam != lambdas->end() ; ++it_lam){

	 nlam++;
	 //...isLambda = nlam;  //compiler didn't like it for some mysterious reason! but liked the above isKshort
	 bool   isLambda = nlam;
	 double lamvertexx = it_lam->vx();
	 double lamvertexy = it_lam->vy();
	 double lamvertexz = it_lam->vz();
	 double lampt = it_lam->pt();
	 double lameta = it_lam->eta();
	 double lamphi = it_lam->phi();
	 double lammass = it_lam->mass();
	 // ...do not mix with primary vertex !!!
	 //	 double lamradius = TMath::Sqrt((lamvertexx-xvtx)*(lamvertexx-xvtx)+(lamvertexy-yvtx)*(lamvertexy-yvtx));
	 // 	 double lam3Dradius = TMath::Sqrt((lamvertexx-xvtx)*(lamvertexx-xvtx)+(lamvertexy-yvtx)*(lamvertexy-yvtx)+(lamvertexz-zvtx)*(lamvertexz-zvtx));
	 double lamradius = TMath::Sqrt((lamvertexx)*(lamvertexx)+(lamvertexy)*(lamvertexy));
 	 double lam3Dradius = TMath::Sqrt((lamvertexx)*(lamvertexx)+(lamvertexy)*(lamvertexy)+(lamvertexz)*(lamvertexz));
  	 double lamenergy = TMath::Sqrt(lampt*lampt+1.115683*1.115683);
	 double gammalorentzlam = lamenergy/1.115683;
	 double lamlifetime = lamradius/gammalorentzlam;  
	 double lam3Dlifetime = lam3Dradius/gammalorentzlam;  
	 histosTH1F["hlampt"]->Fill(lampt);
	 histosTH1F["hlameta"]->Fill(lameta);
	 histosTH1F["hlamphi"]->Fill(lamphi);
	 histosTH1F["hlammass"]->Fill(lammass);
	 histosTH1F["hlamvertexx"]->Fill(lamvertexx);
	 histosTH1F["hlamvertexy"]->Fill(lamvertexy);
	 histosTH1F["hlamvertexz"]->Fill(lamvertexz);
	 histosTH1F["hlamradius"]->Fill(lamradius);
	 histosTH1F["hlam3Dradius"]->Fill(lam3Dradius);
	 histosTH1F["hlamlifetime"]->Fill(lamlifetime);
	 histosTH1F["hlam3Dlifetime"]->Fill(lam3Dlifetime);
	 histosTH2F["h2dimlamxy"]->Fill(lamvertexx,lamvertexy);
	 histosTH2F["h2dimlamxz"]->Fill(lamvertexx,lamvertexz);
	 histosTH2F["h2dimlamyz"]->Fill(lamvertexy,lamvertexz);
	 //std::cout << " ksvertexx = " << ksvertexx << std::endl;
	 //std::cout << " ksvertexy = " << ksvertexy << std::endl;
	 //std::cout << " ksvertexz = " << ksvertexz << std::endl;
	 //std::cout << " ksmass = " << ksmass << std::endl;
	 //it_lam->Print();
       }
       //...end of Lambda
       histosTH1F["hnlam"]->Fill(nlam);
       histosTH2F["hntrknlam"]->Fill(ntrk,nlam);
       histosTH2F["hnvtxnlam"]->Fill(nvtx,nlam);
*/
       
  //--------------------------------------------------------------------------------
  //--------------------------------------------------------------------------------
  // PPS

  // 2018 setup
  //------------------------------------------------      
  // -z                    IP               +z
  //         sec45                  sec56
  //top:  24       4           104         124
  //ver:     23  3                 103 123      
  //bot:  25       5           105         125
  //
  //------------------------------------
  
  //  load track data;
  //  fwlite::Handle< vector<CTPPSLocalTrackLite> > RPtracks;
  //  RPtracks.getByLabel(iEvent, "ctppsLocalTrackLiteProducer");

  bool rp_valid_004 = false;
  bool rp_valid_005 = false;
  bool rp_valid_024 = false;
  bool rp_valid_025 = false;
  
  bool rp_valid_104 = false;
  bool rp_valid_105 = false;
  bool rp_valid_124 = false;
  bool rp_valid_125 = false;
  
  double xLN=100, xLF=100, yLN=100, yLF=100;
  double xRN=100, xRF=100, yRN=100, yRF=100;

  //here is a new set of alignment constants, this time extracted from the
  //physics runs:
  //alignment_corrections[24] = UnitHitData(false, -0.465, -0.689);
  //alignment_corrections[4] = UnitHitData(false, -0.210, -1.479);
  //alignment_corrections[104] = UnitHitData(false, +0.167, -0.916);
  //alignment_corrections[124] = UnitHitData(false, -0.450, +0.044);
  //alignment_corrections[25] = UnitHitData(false, -0.081, +0.009);
  //alignment_corrections[5] = UnitHitData(false, -0.112, +0.842);
  //alignment_corrections[105] = UnitHitData(false, +0.373, +1.312);
  //alignment_corrections[125] = UnitHitData(false, -0.574, +0.316);
  //The format is:
  //  UnitHitData(false, MEAN_X, MEAN_Y)

  //T, from L to R, as on Jan's slides
  double mean_x24  = -0.465;
  double mean_x4   = -0.210;
  double mean_x104 =  0.167;
  double mean_x124 = -0.450;
  
  //B, from L to R
  double mean_x25  = -0.081;
  double mean_x5   = -0.112;
  double mean_x105 =  0.373;
  double mean_x125 = -0.574;
  
  //T, from L to R
  double mean_y24  = -0.689;
  double mean_y4   = -1.479;
  double mean_y104 = -0.916;
  double mean_y124 =  0.044;
  
  //B, from L to R
  double mean_y25  = 0.009;
  double mean_y5   = 0.842;
  double mean_y105 = 1.312;
  double mean_y125 = 0.316;

  // process track data
  for (const auto &tr : *RPtracks){
    
    CTPPSDetId rpId(tr.getRPId());
    unsigned int rpDecId = 100*rpId.arm() + 10*rpId.station() + 1*rpId.rp();
    
    //    std::cout<<"rpDecId= "<<rpDecId<<std::endl;
    //std::cout<<" --- RP Id --- "<<std::endl;
    //std::cout<<"rpDecId= "<<rpDecId<<std::endl;
    //std::cout<<"rpId.arm = "<<rpId.arm()<<std::endl;
    //std::cout<<"rpId.station= "<<rpId.station()<<std::endl;
    //std::cout<<"rpId.rp= "<<rpId.rp()<<std::endl;
    //std::cout<<"                          "<<std::endl;
	
    if(rpDecId == 4){rp_valid_004 = true; xLN = tr.getX() + mean_x4; yLN = tr.getY() + mean_y4;}
    if(rpDecId == 5){rp_valid_005 = true; xLN = tr.getX() + mean_x5; yLN = tr.getY() + mean_y5;}

    if(rpDecId == 24){rp_valid_024 = true; xLF = tr.getX() + mean_x24; yLF = tr.getY() + mean_y24;}
    if(rpDecId == 25){rp_valid_025 = true; xLF = tr.getX() + mean_x25; yLF = tr.getY() + mean_y25;}
    
    if(rpDecId == 104){rp_valid_104 = true; xRN = tr.getX() + mean_x104; yRN = tr.getY() + mean_y104;}
    if(rpDecId == 105){rp_valid_105 = true; xRN = tr.getX() + mean_x105; yRN = tr.getY() + mean_y105;}

    if(rpDecId == 124){rp_valid_124 = true; xRF = tr.getX() + mean_x124; yRF = tr.getY() + mean_y124;}
    if(rpDecId == 125){rp_valid_125 = true; xRF = tr.getX() + mean_x125; yRF = tr.getY() + mean_y125;}

    //    if(rpDecId == 4){rp_valid_004 = true; xLN = tr.getX(); yLN = tr.getY();}
    //    if(rpDecId == 5){rp_valid_005 = true; xLN = tr.getX(); yLN = tr.getY();}

    //    if(rpDecId == 24){rp_valid_024 = true; xLF = tr.getX(); yLF = tr.getY();}
    //    if(rpDecId == 25){rp_valid_025 = true; xLF = tr.getX(); yLF = tr.getY();}
    
    //    if(rpDecId == 104){rp_valid_104 = true; xRN = tr.getX(); yRN = tr.getY();}
    //    if(rpDecId == 105){rp_valid_105 = true; xRN = tr.getX(); yRN = tr.getY();}

    //    if(rpDecId == 124){rp_valid_124 = true; xRF = tr.getX(); yRF = tr.getY();}
    //    if(rpDecId == 125){rp_valid_125 = true; xRF = tr.getX(); yRF = tr.getY();}
  }

  bool diag_top45_bot56 = rp_valid_024 && rp_valid_004 && rp_valid_105 && rp_valid_125;
  bool diag_bot45_top56 = rp_valid_025 && rp_valid_005 && rp_valid_104 && rp_valid_124;
  
  bool top45_top56      = rp_valid_024 && rp_valid_004 && rp_valid_104 && rp_valid_124;
  bool bot45_bot56      = rp_valid_025 && rp_valid_005 && rp_valid_105 && rp_valid_125;
  
  int nconf=0;
  if(diag_top45_bot56) nconf++;
  if(diag_bot45_top56) nconf++;
  if(top45_top56) nconf++;
  if(bot45_bot56) nconf++;
  
  histosTH1F["hnconf"]->Fill(nconf);
  if(nconf != 1) return;
  
  //  bool diag=false;
  //  if(diag_top45_bot56 || diag_bot45_top56) diag = true;
  
  //Topol
  //1 - TB, 2 - BT
  //3 - TT, 4 - BB       
  
  int tb=-1;
  
  if(diag_top45_bot56) tb=0;
  if(diag_bot45_top56) tb=1;
  if(top45_top56) tb=2;
  if(bot45_bot56) tb=3;
  
  histosTH1F["hconf"]->Fill(tb);
  

  // ----- single-arm kinematics reconstruction -----
  
  //  double ThxR, ThyR, ThxL, ThyL;//, xVtxL, xVtxR;
  double ThxR, ThyR, ThxL, ThyL, xVtxL, xVtxR;
  
  //  double D_x_L = - v_x_L_1_F * L_x_L_2_F + v_x_L_2_F * L_x_L_1_F;
  //sign convention
  double D_x_L = + v_x_L_1_F * L_x_L_2_F - v_x_L_2_F * L_x_L_1_F;
  ThxL = (v_x_L_1_F * xLF - v_x_L_2_F * xLN) / D_x_L;
  xVtxL = (- xLN * L_x_L_2_F + xLF * L_x_L_1_F) / D_x_L;
  
  double D_x_R = + v_x_R_1_F * L_x_R_2_F - v_x_R_2_F * L_x_R_1_F;
  ThxR = (v_x_R_1_F * xRF - v_x_R_2_F * xRN) / D_x_R;
  xVtxR = (+ xRN * L_x_R_2_F - xRF * L_x_R_1_F) / D_x_R;
  
  //  double th_y_L_1_F = - yLN / L_y_L_1_F;
  //  double th_y_L_2_F = - yLF / L_y_L_2_F;
  // sign convention
  double th_y_L_1_F = + yLN / L_y_L_1_F;
  double th_y_L_2_F = + yLF / L_y_L_2_F;
  ThyL = (th_y_L_1_F + th_y_L_2_F) / 2.;
  
  double th_y_R_1_F = + yRN / L_y_R_1_F;
  double th_y_R_2_F = + yRF / L_y_R_2_F;
  ThyR = (th_y_R_1_F + th_y_R_2_F) / 2.;

  // ----- theta reconstruction -----
  //  double th_sq = k.th_x*k.th_x + k.th_y*k.th_y;
  //  k.th = sqrt(th_sq);
  //  k.phi = atan2(k.th_y, k.th_x);
  
  // t reconstruction
  //  k.t_x = env.p*env.p * k.th_x * k.th_x;
  //  k.t_y = env.p*env.p * k.th_y * k.th_y;
  //  k.t = k.t_x + k.t_y;
  
  // Correct residual shifts in thx (only, not needed thy)  
  // 2015
  //	if(specialreco) ThxL=rec_proton_left->thx-5.04e-5;
  // 2018
  // was on during the run for Express
  // Gauss fit: shift xL+xR = -1.80371e-04
  // ThxR += 1.815e-04;
  
  //my calculations from shift in dpx/6500
  
  double a_off =  0.000002386 ; //TB
  double b_off = -0.000006593 ; //BT
  double c_off = -0.000007524 ; //TT
  double d_off =  0.000003268 ; //BB
  
  if(tb==0){ThxL += 0. ;           ThxR += a_off ;}//TB
  if(tb==1){ThxL += (b_off-c_off); ThxR += c_off ;}//BT
  if(tb==2){ThxL += 0. ;           ThxR += c_off ;}//TT
  if(tb==3){ThxL += (d_off-a_off); ThxR += a_off ;}//BB


  histosTH1F["hthxEla"]->Fill(ThxL+ThxR);
  histosTH1F["hthyEla"]->Fill(ThyL+ThyR);

  
  histosTH2F["hthx2DIM"]->Fill(ThxL,ThxR);
  histosTH2F["hthy2DIM"]->Fill(ThyL,ThyR);
  histosTH2F["hthythx2DIM"]->Fill(ThxL+ThxR,ThyL+ThyR);

  // smug alert!
  
  // |-t1|
  double tt1 = 6500*6500*(ThxR*ThxR+ThyR*ThyR);
  // |-t2|
  double tt2 = 6500*6500*(ThxL*ThxL+ThyL*ThyL);

  // Theta versus t1,t2         ...no cuts
    histosTH2F["hthyrtt1"]->Fill(tt1,ThyR);
    histosTH2F["hthyltt2"]->Fill(tt2,ThyL);
    histosTH2F["hthxrtt1"]->Fill(tt1,ThxR);
    histosTH2F["hthxltt2"]->Fill(tt2,ThxL);

  //---------------------------------------------Elastic off
  
  bool isElastic = false;
  
  // 5 sigma
  if(TMath::Abs(ThyL+ThyR)< 15e-6 && 
     TMath::Abs(ThxL+ThxR)<45e-6) isElastic=true;
  
  if(isElastic) return;
  //non-elastic
  histosTH1F["hthxEla2"]->Fill(ThxL+ThxR);
  histosTH1F["hthyEla2"]->Fill(ThyL+ThyR);

  //---------------------------------------------
  // after vtx cut, proton cut, elastic cut

  if(runnr == 319104) histosTH1F["hLS104"]->Fill(LS);
  if(runnr == 319124) histosTH1F["hLS124"]->Fill(LS);
  if(runnr == 319125) histosTH1F["hLS125"]->Fill(LS);
  if(runnr == 319159) histosTH1F["hLS159"]->Fill(LS);
  if(runnr == 319174) histosTH1F["hLS174"]->Fill(LS);
  if(runnr == 319175) histosTH1F["hLS175"]->Fill(LS);
  if(runnr == 319176) histosTH1F["hLS176"]->Fill(LS);
  if(runnr == 319177) histosTH1F["hLS177"]->Fill(LS);
  if(runnr == 319190) histosTH1F["hLS190"]->Fill(LS);
  
  if(runnr == 319222) histosTH1F["hLS222"]->Fill(LS);
  if(runnr == 319223) histosTH1F["hLS223"]->Fill(LS);
  if(runnr == 319254) histosTH1F["hLS254"]->Fill(LS);
  if(runnr == 319255) histosTH1F["hLS255"]->Fill(LS);
  if(runnr == 319256) histosTH1F["hLS256"]->Fill(LS);
  if(runnr == 319262) histosTH1F["hLS262"]->Fill(LS);
  if(runnr == 319263) histosTH1F["hLS263"]->Fill(LS);
  if(runnr == 319264) histosTH1F["hLS264"]->Fill(LS);
  if(runnr == 319265) histosTH1F["hLS265"]->Fill(LS);
  if(runnr == 319266) histosTH1F["hLS266"]->Fill(LS);
  if(runnr == 319267) histosTH1F["hLS267"]->Fill(LS);
  if(runnr == 319268) histosTH1F["hLS268"]->Fill(LS);
  if(runnr == 319270) histosTH1F["hLS270"]->Fill(LS);

  if(runnr == 319300) histosTH1F["hLS300"]->Fill(LS);
  if(runnr == 319311) histosTH1F["hLS311"]->Fill(LS);

  //-------------------------------
  //  if(!jsonLocal(runnr,LS)) return;

  //-------------------------------
  // Forward protons ...Luiz
  //for (const auto &rpProton : *ProtonsMultiRP ) {
  /*
  for (const auto &rpProton : *ProtonsSingleRP ) {
    
  CTPPSDetId rpId((*rpProton.contributingLocalTracks().begin())->getRPId());
    unsigned int decRPId = rpId.arm() * 100 + rpId.station() * 10 + rpId.rp();
    
  double proton_t = rpProton.t();
  double proton_xi = rpProton.xi();
  double proton_p = rpProton.p();
  double proton_pt = rpProton.pt();
  double proton_thx = rpProton.thetaX();
  double proton_thy = rpProton.thetaY();
  double proton_mass = rpProton.mass();

  std::cout<<" --- Forward Protons RP Id --- "<<std::endl;
  std::cout<<"decRPId= "<<decRPId<<std::endl;

  std::cout << " *** forward proton ***" <<std::endl;
  std::cout << " t    = " << proton_t <<std::endl;
  std::cout << " xi   = " << proton_xi <<std::endl;
  std::cout << " p    = " << proton_p <<std::endl;
  std::cout << " pt   = " << proton_pt <<std::endl;
  std::cout << " thx  = " << proton_thx <<std::endl;
  std::cout << " thy  = " << proton_thy <<std::endl;
  std::cout << " mass = " << proton_mass <<std::endl;
    }
  */
  //}
  
  //---------------------------------------------
  //CMS-TOTEM matching
  
  double TOTEMpy= 6500.*(ThyL+ThyR);
  double TOTEMpx=-6500.*(ThxL+ThxR);
  
  //  double TOTEMpt1= TMath::Sqrt(pow(TOTEMpx1,2)+pow(TOTEMpy1,2));
  //  double TOTEMpt2= TMath::Sqrt(pow(TOTEMpx2,2)+pow(TOTEMpy2,2));

  double TOTEMphiL = TMath::ATan2(ThyL,ThxL);
  double TOTEMphiR = TMath::ATan2(ThyR,ThxR);
  
  double TOTEMdphi = TOTEMphiL-TOTEMphiR;
  if(TOTEMdphi<0) TOTEMdphi = TOTEMdphi + 2*TMath::Pi(); // from (-2pi,2pi) to (0,2pi)
  if(TOTEMdphi>TMath::Pi()) TOTEMdphi = 2*TMath::Pi() - TOTEMdphi; // from (0,2pi) to (0,pi)
  
  //double CMSpx=pipiRec.Px();
  //double CMSpy=pipiRec.Py();

  //...Ferenc
  // (int & topology, pair<double,double> & pL, pair<double,double> & pR)
  // pL = pair<double,double>(-6500*ThxL, 6500*ThyL);
  // pR = pair<double,double>(-6500*ThxR, 6500*ThyR);
  double proton_left_px  = -6500*ThxL ;
  double proton_left_py  =  6500*ThyL ;
  double proton_right_px = -6500*ThxR ;
  double proton_right_py =  6500*ThyR ;

  //std::cout << " proton_left_py = " << proton_left_py << std::endl;
  //std::cout << " proton_right_py = " << proton_right_py << std::endl;
  
  //...proton transverse momentum
  double proton_left_pt  = TMath::Sqrt((-6500*ThxL)*(-6500*ThxL)+(6500*ThyL)*(6500*ThyL));
  double proton_right_pt = TMath::Sqrt((-6500*ThxR)*(-6500*ThxR)+(6500*ThyR)*(6500*ThyR));

    histosTH1F["hprotonlpx"]->Fill(proton_left_px);
    histosTH1F["hprotonlpy"]->Fill(proton_left_py);
    histosTH1F["hprotonrpx"]->Fill(proton_right_px);
    histosTH1F["hprotonrpy"]->Fill(proton_right_py);
  
    histosTH1F["hprotonlpt"]->Fill(proton_left_pt);
    histosTH1F["hprotonrpt"]->Fill(proton_right_pt);


    
  //...py efficiency

  //proton_left  = arm1 = 1
  //proton_right = arm2 = 2
	     
	     //std::cout<<" efvec[71][0] = "<< efvec[71][0] << " efvec[71][1] = "<< efvec[71][1] <<std::endl;

    double effic1 = -1.;
    double effic2 = -1.;
    double effic  = -1.;
    double effictotal = -1.;
	     
  for(int i = 0; i < 20800; i++)
  {
  //arm1...0
  if( efvec[i][4] == runnr && efvec[i][3] == 1 && proton_left_py >= (efvec[i][0]-0.0025) && proton_left_py < (efvec[i][0]+0.0025) )
    {
      effic1 = efvec[i][1];
    }
  //arm2...1
  if( efvec[i][4] == runnr && efvec[i][3] == 2 && proton_right_py >= (efvec[i][0]-0.0025) && proton_right_py < (efvec[i][0]+0.0025) )
    {
      effic2 = efvec[i][1];
    }    
  }
      effic = effic1*effic2;
  
  if( effic == 0. )
    {
      //effictotal = 1.;
      //...no contribution
      effictotal = 1000.;
    }
  else
    {
      effictotal = effic;
    }
  
  //std::cout<<" run# = "<<runnr<<" *************************************** "<<std::endl;	 
  //std::cout<<" proton_left_py = "<< proton_left_py <<" proton_right_py = "<< proton_right_py <<std::endl;
  //std::cout<<" effic1 = "<< effic1 <<" effic2 = "<< effic2 <<std::endl;


    
  //...4-momentum transfer squared
  // t2 = -p² (theta*_x_R² + theta*_y_R²)
  // t1 = -p² (theta*_x_L² + theta*_y_L²)
  // |-t2|
  double t2 = 6500*6500*(ThxR*ThxR+ThyR*ThyR);
  // |-t1|
  double t1 = 6500*6500*(ThxL*ThxL+ThyL*ThyL);

  // |-(t1+t2)|  ...is it OR or AND ?
  //double Thetax = ThxL+ThxR; //.....................WRONG!
  //double Thetay = ThyL+ThyR; //........................WRONG!
  //double t1t2 = 6500*6500*(Thetax*Thetax+Thetay*Thetay); // ...WRONG!
 
    histosTH1F["ht1"]->Fill(t1);
    histosTH1F["ht2"]->Fill(t2);
    histosTH2F["h2dimt1t2"]->Fill(t1,t2);
    //
    //histosTH1F["ht1t2"]->Fill(t1t2);

    double t12 = t1;
    histosTH1F["ht12"]->Fill(t12);
    t12 = t2;
    histosTH1F["ht12"]->Fill(t12);

    /*
    double t = t1 + t2 ; .............WRONG!
    histosTH1F["ht"]->Fill(t);
    */
    
    
    // t = -p^2 (Theta_x^2+Theta_y^2)        ...Luiz
    // t1 versus left pt    ...parabolic!
    histosTH2F["h2dimt1protonpt"]->Fill(proton_left_pt,t1);
    // t2 versus right pt
    histosTH2F["h2dimt2protonpt"]->Fill(proton_right_pt,t2);
    

    
    //...relative longitudinal momentum loss
    //
    //  xi = ( x - vx * x0 - Lx * Thx ) / Dx
    //
    //...choosing Lx = 0 ? ...why?
    //
    //  xi = ( x - vx * x0 ) / Dx
    //
      
    double xi1_1 = ( xLN - v_x_L_1_F * xVtxL - L_x_L_1_F * ThxL ) / D_x_L ;
    double xi1_2 = ( xLF - v_x_L_2_F * xVtxL - L_x_L_2_F * ThxL ) / D_x_L ;
    //
    double xi2_1 = ( xRN - v_x_R_1_F * xVtxR - L_x_R_1_F * ThxR ) / D_x_R ;
    double xi2_2 = ( xRF - v_x_R_2_F * xVtxR - L_x_R_2_F * ThxR ) / D_x_R ;
    //
    double xi1 = ( xi1_1 + xi1_2 )/2. ;
    double xi2 = ( xi2_1 + xi2_2 )/2. ;
    
    histosTH1F["hxi1"]->Fill(xi1);
    histosTH1F["hxi2"]->Fill(xi2);
    histosTH2F["h2dimxi1xi2"]->Fill(xi1,xi2);

    histosTH1F["hxi1_1"]->Fill(xi1_1);
    histosTH1F["hxi1_2"]->Fill(xi1_2);
    histosTH1F["hxi2_1"]->Fill(xi2_1);
    histosTH1F["hxi2_2"]->Fill(xi2_2);

    histosTH2F["h2dimxi1protonpt"]->Fill(proton_left_pt,xi2);
    histosTH2F["h2dimxi2protonpt"]->Fill(proton_right_pt,xi2);

    //double phi1 = TOTEMphiL;
    //double phi2 = TOTEMphiR;
    //double dphi = TOTEMdphi; // (0,pi)
    histosTH1F["hphi1"]->Fill(TOTEMphiL);
    histosTH1F["hphi2"]->Fill(TOTEMphiR);
    histosTH1F["hdphi"]->Fill(TOTEMdphi);
    histosTH1F["hdphig"]->Fill(TOTEMdphi*180/TMath::Pi());
    histosTH1F["hdphigcor"]->Fill(TOTEMdphi*180/TMath::Pi(),1/effictotal);

    // ...do not work!
    //double histomax = hdphig->GetMaximum();
    //histosTH1F["hdphif"]->Fill(TOTEMdphi*180/TMath::Pi()/histomax);

    // non-elastic
    //theta vs t
    histosTH2F["hthyrt1"]->Fill(t1,ThyR);
    histosTH2F["hthylt2"]->Fill(t2,ThyL);
    histosTH2F["hthxrt1"]->Fill(t1,ThxR);
    histosTH2F["hthxlt2"]->Fill(t2,ThxL);
    //theta vs xi
    histosTH2F["hthyrxi1"]->Fill(xi1,ThyR);
    histosTH2F["hthylxi2"]->Fill(xi2,ThyL);
    histosTH2F["hthxrxi1"]->Fill(xi1,ThxR);
    histosTH2F["hthxlxi2"]->Fill(xi2,ThxL);

  //---------------------------------------------------------------
  // V0
  // ...Kshort collection                           ...non-elastic!

  //...Kshort collection...Luiz
  bool isKshort = false;
  int nks=0;
  double ksmass = 0.0;

  //TMath::XYZPoint referencePOS(theBeamSpot->position());
  
  for(VertexCompositeCandidateCollection::const_iterator it_ks = kshorts->begin() ; it_ks != kshorts->end() ; ++it_ks){

	 nks++;
	 isKshort = nks;
	 double ksvertexx = it_ks->vx();
	 double ksvertexy = it_ks->vy();
	 double ksvertexz = it_ks->vz();
	 double kspt = it_ks->pt();
	 double kspz = it_ks->pz();      //...attention!
         //std::cout << " ++++++++++++++++++++++++ " <<std::endl;
         //std::cout << " kspz = "  << kspz <<std::endl;
	 double kseta = it_ks->eta();
	 double ksphi = it_ks->phi();
	 //double ksmass = it_ks->mass();
	 ksmass = it_ks->mass();
	 // ...consider the primary vertex !!!
	 // reference point is not the beam spot
 	 double ksradiusvtx = TMath::Sqrt((ksvertexx-xvtx)*(ksvertexx-xvtx)+(ksvertexy-yvtx)*(ksvertexy-yvtx));
 	 double ks3Dradiusvtx = TMath::Sqrt((ksvertexx-xvtx)*(ksvertexx-xvtx)+(ksvertexy-yvtx)*(ksvertexy-yvtx)+(ksvertexz-zvtx)*(ksvertexz-zvtx));
	 //...wrt BS
 	 double ksradius = TMath::Sqrt((ksvertexx-xBS)*(ksvertexx-xBS)+(ksvertexy-yBS)*(ksvertexy-yBS));
 	 double ks3Dradius = TMath::Sqrt((ksvertexx-xBS)*(ksvertexx-xBS)+(ksvertexy-yBS)*(ksvertexy-yBS)+(ksvertexz-zBS)*(ksvertexz-zBS));
	 //double ksradius = TMath::Sqrt((ksvertexx)*(ksvertexx)+(ksvertexy)*(ksvertexy));
	 //double ks3Dradius = TMath::Sqrt((ksvertexx)*(ksvertexx)+(ksvertexy)*(ksvertexy)+(ksvertexz)*(ksvertexz));
	 double ksenergy = TMath::Sqrt( kspt*kspt + 0.4976*0.4976 );
	 double gammalorentzks = ksenergy/0.4976;
	 double kslifetime = ksradius/gammalorentzks;  
	 double ks3Dlifetime = ks3Dradius/gammalorentzks;
	 double kslifetimevtx = ksradiusvtx/gammalorentzks;  
	 double ks3Dlifetimevtx = ks3Dradiusvtx/gammalorentzks;
         std::cout << " * * * * "<< std::endl;
         std::cout << " ksenergy = "<< ksenergy << std::endl;
         std::cout << " ks3Dlifetime = "<< ks3Dlifetime << std::endl;
         std::cout << " ks3Dlifetimevtx = "<< ks3Dlifetimevtx << std::endl;
	 // ...using pz of K0s
	 double ksenergyz = TMath::Sqrt( kspt*kspt + kspz*kspz + 0.4976*0.4976 );      //...attention!
	 double gammalorentzksz = ksenergyz/0.4976;
	 double kslifetimez = ksradius/gammalorentzksz;  
	 double ks3Dlifetimez = ks3Dradius/gammalorentzksz;
	 double kslifetimezvtx = ksradiusvtx/gammalorentzksz;  
	 double ks3Dlifetimezvtx = ks3Dradiusvtx/gammalorentzksz;
         std::cout << " ksenergyz = "<< ksenergyz << std::endl;
         std::cout << " ks3Dlifetimez = "<< ks3Dlifetimez << std::endl;
         std::cout << " ks3Dlifetimezvtx = "<< ks3Dlifetimezvtx << std::endl;
	 
         //std::cout << " ksenergy = " << ksenergy <<std::endl;
         //std::cout << " gammalorentzks = " << gammalorentzks <<std::endl;
         //std::cout << " ks3Dlifetime = " << ks3Dlifetime <<std::endl;
	 
	 histosTH1F["hkspt"]->Fill(kspt);
	 histosTH1F["hkseta"]->Fill(kseta);
	 histosTH1F["hksphi"]->Fill(ksphi);
	 histosTH1F["hksmass"]->Fill(ksmass);
	 histosTH1F["hksmasscor"]->Fill(ksmass,1/effictotal);
	 //
	 if(nks == 1){histosTH1F["hksmassv1"]->Fill(ksmass);
                      histosTH1F["hksradiusv1"]->Fill(ksradius);
                      histosTH1F["hkslifetimev1"]->Fill(kslifetime);
                      histosTH1F["hks3Dradiusv1"]->Fill(ks3Dradius);
                      histosTH1F["hks3Dlifetimev1"]->Fill(ks3Dlifetime);
	 }
	 /*
	 if(nks == 2){histosTH1F["hksmassv2"]->Fill(ksmass);
                      histosTH1F["hksradiusv2"]->Fill(ksradius);
                      histosTH1F["hkslifetimev2"]->Fill(kslifetime);
                      histosTH1F["hks3Dradiusv2"]->Fill(ks3Dradius);
                      histosTH1F["hks3Dlifetimev2"]->Fill(ks3Dlifetime);
	 }
	 */
	 if(ntrk==4 && TMath::Abs(pi1.Eta())<etaCut && TMath::Abs(pi2.Eta())<etaCut &&
       		       TMath::Abs(pi3.Eta())<etaCut && TMath::Abs(pi4.Eta())<etaCut &&
	               pi1.Pt()>ptCut && pi2.Pt()>ptCut && pi3.Pt()>ptCut && pi4.Pt()>ptCut){
	   if(totcharge == 0){
	   if(nks == 2){
	     if(nvtx==0){
	              histosTH1F["hksmassv2"]->Fill(ksmass);
	              histosTH1F["hksmassv2cor"]->Fill(ksmass,1/effictotal);
                      histosTH1F["hksradiusv2"]->Fill(ksradius);
                      histosTH1F["hks3Dradiusv2"]->Fill(ks3Dradius);
                      histosTH1F["hksradiusv2vtx"]->Fill(ksradiusvtx);
                      histosTH1F["hks3Dradiusv2vtx"]->Fill(ks3Dradiusvtx);
		      // cor
                      histosTH1F["hkslifetimev2"]->Fill(kslifetime);
                      histosTH1F["hks3Dlifetimev2"]->Fill(ks3Dlifetime);
                      histosTH1F["hks3Dlifetimev2cor"]->Fill(ks3Dlifetime,1/effictotal);
		      // vtx cor
                      histosTH1F["hkslifetimev2vtx"]->Fill(kslifetimevtx);
                      histosTH1F["hks3Dlifetimev2vtx"]->Fill(ks3Dlifetimevtx);
                      histosTH1F["hks3Dlifetimev2vtxcor"]->Fill(ks3Dlifetimevtx,1/effictotal);
		      // z cor
                      histosTH1F["hkslifetimev2z"]->Fill(kslifetimez);
                      histosTH1F["hks3Dlifetimev2z"]->Fill(ks3Dlifetimez);
                      histosTH1F["hks3Dlifetimev2zcor"]->Fill(ks3Dlifetimez,1/effictotal);
		      // z vtx cor
                      histosTH1F["hkslifetimev2zvtx"]->Fill(kslifetimezvtx);
                      histosTH1F["hks3Dlifetimev2zvtx"]->Fill(ks3Dlifetimezvtx);
                      histosTH1F["hks3Dlifetimev2zvtxcor"]->Fill(ks3Dlifetimezvtx,1/effictotal);
	     }
	     //std::cout << " * * * ntrk = 4 & Q = 0"<< std::endl;
	     //std::cout << " nks = " << nks << std::endl;
	     //std::cout << " nvtx = " << nvtx << std::endl;
	     }
	   }
	 }
	 //
	 histosTH1F["hksvertexx"]->Fill(ksvertexx);
	 histosTH1F["hksvertexy"]->Fill(ksvertexy);
	 histosTH1F["hksvertexz"]->Fill(ksvertexz);
	 histosTH1F["hksradius"]->Fill(ksradius);
	 histosTH1F["hks3Dradius"]->Fill(ks3Dradius);
	 histosTH1F["hksradiusvtx"]->Fill(ksradiusvtx);
	 histosTH1F["hks3Dradiusvtx"]->Fill(ks3Dradiusvtx);
	 // cor
	 histosTH1F["hkslifetime"]->Fill(kslifetime);
	 histosTH1F["hks3Dlifetime"]->Fill(ks3Dlifetime);
	 histosTH1F["hks3Dlifetimecor"]->Fill(ks3Dlifetime,1/effictotal);
	 // vtx cor
	 histosTH1F["hkslifetimevtx"]->Fill(kslifetimevtx);
	 histosTH1F["hks3Dlifetimevtx"]->Fill(ks3Dlifetimevtx);
	 histosTH1F["hks3Dlifetimevtxcor"]->Fill(ks3Dlifetimevtx,1/effictotal);
	 // z
	 histosTH1F["hkslifetimez"]->Fill(kslifetimez);
	 histosTH1F["hks3Dlifetimez"]->Fill(ks3Dlifetimez);
	 histosTH1F["hks3Dlifetimecorz"]->Fill(ks3Dlifetimez,1/effictotal);
	 // z vtx cor
	 histosTH1F["hkslifetimezvtx"]->Fill(kslifetimezvtx);
	 histosTH1F["hks3Dlifetimezvtx"]->Fill(ks3Dlifetimezvtx);
	 histosTH1F["hks3Dlifetimezvtxcor"]->Fill(ks3Dlifetimezvtx,1/effictotal);
	 // 2D
	 histosTH2F["h2dimksxy"]->Fill(ksvertexx,ksvertexy);
	 histosTH2F["h2dimksxz"]->Fill(ksvertexx,ksvertexz);
	 histosTH2F["h2dimksyz"]->Fill(ksvertexy,ksvertexz);
	 
	 
	 /*
    	 std::cout << " * * * "<< std::endl;
	 std::cout << " nks = " << nks << std::endl;
	 std::cout << " ksvertexx = " << ksvertexx << std::endl;
	 std::cout << " ksvertexy = " << ksvertexy << std::endl;
	 std::cout << " ksvertexz = " << ksvertexz << std::endl;
	 std::cout << " xvtx = " << xvtx << std::endl;
	 std::cout << " yvtx = " << yvtx << std::endl;
	 std::cout << " zvtx = " << zvtx << std::endl;
	 std::cout << " ksmass = " << ksmass << std::endl;
	 std::cout << " kspt = " << kspt << std::endl;
	 std::cout << " ksradius = " << ksradius << std::endl;
	 */
	 
	 //it_ks->Print();
       }
       //...end of Kshort
       histosTH1F["hnks"]->Fill(nks);
       histosTH2F["hntrknks"]->Fill(ntrk,nks);
       histosTH2F["hnvtxnks"]->Fill(nvtx,nks);
       histosTH2F["hntrknvtx"]->Fill(ntrk,nvtx);
       //std::cout << " --------------------------- " << std::endl;
       //std::cout << " nks  = " << nks << std::endl;
       //std::cout << " ntrk = " << ntrk << std::endl;
       //std::cout << " nvtx = " << nvtx << std::endl;
       //std::cout << " isKshort = " << isKshort << std::endl;
       //std::cout << " --------------------------- " << std::endl;

       //std::cout << " gammalorentzks 2 = " << gammalorentzks <<std::endl;
	  
   //...Lambda collection...Luiz
   ////bool isLambda = false;
   //std::cout << "isLambda ...boo! " << isLambda << std::endl;
   int nlam=0;

   for(VertexCompositeCandidateCollection::const_iterator it_lam = lambdas->begin() ; it_lam != lambdas->end() ; ++it_lam){

	 nlam++;
	 //...isLambda = nlam;  //compiler didn't like it for some mysterious reason! but liked the above isKshort
	 bool   isLambda = nlam;
	 double lamvertexx = it_lam->vx();
	 double lamvertexy = it_lam->vy();
	 double lamvertexz = it_lam->vz();
	 double lampt = it_lam->pt();
	 double lameta = it_lam->eta();
	 double lamphi = it_lam->phi();
	 double lammass = it_lam->mass();
	 // ...w.r.t primary vertex !!!
	 double lamradius = TMath::Sqrt((lamvertexx-xvtx)*(lamvertexx-xvtx)+(lamvertexy-yvtx)*(lamvertexy-yvtx));
	 double lam3Dradius = TMath::Sqrt((lamvertexx-xvtx)*(lamvertexx-xvtx)+(lamvertexy-yvtx)*(lamvertexy-yvtx)+(lamvertexz-zvtx)*(lamvertexz-zvtx));
	 //double lamradius = TMath::Sqrt((lamvertexx)*(lamvertexx)+(lamvertexy)*(lamvertexy));
 	 //double lam3Dradius = TMath::Sqrt((lamvertexx)*(lamvertexx)+(lamvertexy)*(lamvertexy)+(lamvertexz)*(lamvertexz));
  	 double lamenergy = TMath::Sqrt(lampt*lampt+1.115683*1.115683);
	 double gammalorentzlam = lamenergy/1.115683;
	 double lamlifetime = lamradius/gammalorentzlam;  
	 double lam3Dlifetime = lam3Dradius/gammalorentzlam;  
	 histosTH1F["hlampt"]->Fill(lampt);
	 histosTH1F["hlameta"]->Fill(lameta);
	 histosTH1F["hlamphi"]->Fill(lamphi);
	 histosTH1F["hlammass"]->Fill(lammass);
	 histosTH1F["hlamvertexx"]->Fill(lamvertexx);
	 histosTH1F["hlamvertexy"]->Fill(lamvertexy);
	 histosTH1F["hlamvertexz"]->Fill(lamvertexz);
	 histosTH1F["hlamradius"]->Fill(lamradius);
	 histosTH1F["hlam3Dradius"]->Fill(lam3Dradius);
	 histosTH1F["hlamlifetime"]->Fill(lamlifetime);
	 histosTH1F["hlam3Dlifetime"]->Fill(lam3Dlifetime);
	 histosTH2F["h2dimlamxy"]->Fill(lamvertexx,lamvertexy);
	 histosTH2F["h2dimlamxz"]->Fill(lamvertexx,lamvertexz);
	 histosTH2F["h2dimlamyz"]->Fill(lamvertexy,lamvertexz);
	 //std::cout << " lamvertexx = " << lamvertexx << std::endl;
	 //std::cout << " lamvertexy = " << lamvertexy << std::endl;
	 //std::cout << " lamvertexz = " << lamvertexz << std::endl;
	 //std::cout << " lammass = " << lammass << std::endl;
	 //it_lam->Print();
       }
       //...end of Lambda
       histosTH1F["hnlam"]->Fill(nlam);
       histosTH2F["hntrknlam"]->Fill(ntrk,nlam);
       histosTH2F["hnvtxnlam"]->Fill(nvtx,nlam);
    


    
  //--------------------------------------------------------------------------------
  // 2 tracks

  double CMSpx2=pipiRec.Px();
  double CMSpy2=pipiRec.Py();

    bool CTpxcut2 = TMath::Abs(CMSpx2+TOTEMpx)<0.15;    
    bool CTpycut2 = TMath::Abs(CMSpy2+TOTEMpy)<0.06;
    
    bool allCuts2 = (CTpxcut2 && CTpycut2);

  if(ntrk==2 && totcharge==0){

    //...Luiz
    //bool CTpxcut = TMath::Abs(CMSpx+TOTEMpx)<0.15;    
    //bool CTpycut = TMath::Abs(CMSpy+TOTEMpy)<0.06;
    
    //bool allCuts = CTpxcut && CTpycut;

    histosTH1F["hdpy2trk"]->Fill(CMSpy2+TOTEMpy);
    histosTH1F["hdpx2trk"]->Fill(CMSpx2+TOTEMpx);
    if(CTpxcut2) histosTH1F["hdpy2trkB"]->Fill(CMSpy2+TOTEMpy);
    if(CTpycut2) histosTH1F["hdpx2trkB"]->Fill(CMSpx2+TOTEMpx);
    
    histosTH2F["h2DIMdpy2trk"]->Fill(CMSpy2,TOTEMpy);
    histosTH2F["h2DIMdpx2trk"]->Fill(CMSpx2,TOTEMpx);
    
    //Mass 2 tracks
    //double mrec = pipiRec.M();
    double mrec2 = pipiRec.M();

    if(allCuts2){

      //histosTH1F["hm2"]->Fill(mrec);      
      histosTH1F["hm2"]->Fill(mrec2);      

      histosTH1F["hpt2"]->Fill(pi4pos1.Pt());      //????
      histosTH1F["hpt2"]->Fill(pi4neg1.Pt());
      
      histosTH1F["heta2"]->Fill(pi4pos1.Eta());
      histosTH1F["heta2"]->Fill(pi4neg1.Eta());
	
      histosTH1F["hdphi2"]->Fill(TOTEMdphi);
      
    }

  } 

  //--------------------------------------------------------------------------------
  // 4 tracks
  
  double CMSpx4=pi4Rec.Px();
  double CMSpy4=pi4Rec.Py();

    bool CTpxcut4 = TMath::Abs(CMSpx4+TOTEMpx)<0.15;    
    bool CTpycut4 = TMath::Abs(CMSpy4+TOTEMpy)<0.06;
    
    bool allCuts4 = (CTpxcut4 && CTpycut4);
    
  if(ntrk==4 && totcharge==0){

    //...Luiz
    //bool CTpxcut4 = TMath::Abs(CMSpx4+TOTEMpx)<0.15;
    //bool CTpycut4 = TMath::Abs(CMSpy4+TOTEMpy)<0.06;
    
    //bool allCuts4 = CTpxcut4 && CTpycut4;
    
    histosTH1F["hdpy4trk"]->Fill(CMSpy4+TOTEMpy);
    histosTH1F["hdpx4trk"]->Fill(CMSpx4+TOTEMpx);
    if(CTpxcut4) histosTH1F["hdpy4trkB"]->Fill(CMSpy4+TOTEMpy);
    if(CTpycut4) histosTH1F["hdpx4trkB"]->Fill(CMSpx4+TOTEMpx);
    
    histosTH2F["h2DIMdpy4trk"]->Fill(CMSpy4,TOTEMpy);
    histosTH2F["h2DIMdpx4trk"]->Fill(CMSpx4,TOTEMpx);

    //------------------
    //Mass 4 tracks

    double mrec4=pi4Rec.M();
    //std::cout << " mrec4 = " << mrec4 <<std::endl;

    if(allCuts4){

      //...no eta cut
      histosTH1F["hm4"]->Fill(mrec4);
      
      histosTH1F["hpt4"]->Fill(pi4pos1.Pt());
      histosTH1F["hpt4"]->Fill(pi4neg1.Pt());
      histosTH1F["hpt4"]->Fill(pi4pos2.Pt());
      histosTH1F["hpt4"]->Fill(pi4neg2.Pt());
      
      histosTH1F["heta4"]->Fill(pi4pos1.Eta());
      histosTH1F["heta4"]->Fill(pi4neg1.Eta());
      histosTH1F["heta4"]->Fill(pi4pos2.Eta());
      histosTH1F["heta4"]->Fill(pi4neg2.Eta());

      histosTH1F["hdphi4"]->Fill(TOTEMdphi);
     
    }

  }

    //-------------------------------------------
    //...Luiz
    //...pions
    fiducialRegion4 = (ntrk==4 && TMath::Abs(pi1.Eta())<etaCut && TMath::Abs(pi2.Eta())<etaCut &&
       		   TMath::Abs(pi3.Eta())<etaCut && TMath::Abs(pi4.Eta())<etaCut);  
    fiducialRegionPt4 = (ntrk==4 && pi1.Pt()>ptCut && pi2.Pt()>ptCut &&
       			   pi3.Pt()>ptCut && pi4.Pt()>ptCut);
    //...kaons
    ////fiducialRegionK4   = (ntrk==4 && TMath::Abs(k1.Eta())<etaCut && TMath::Abs(k2.Eta())<etaCut &&
    ////		   TMath::Abs(k3.Eta())<etaCut && TMath::Abs(k4.Eta())<etaCut);  
    ////fiducialRegionPtK4 = (ntrk==4 && k1.Pt()>ptCut && k2.Pt()>ptCut &&
    ////		   k3.Pt()>ptCut && k4.Pt()>ptCut);
    //-------------------------------------------

    //...Luiz
    double mrec=pipipipiRec.M();
    double mrecKKKK=kkkkRec.M();      

    //...Luiz
    histosTH2F["h2dimdpyAll"]->Fill(CMSpy4,TOTEMpy);
    histosTH1F["hdpyAll"]->Fill(CMSpy4+TOTEMpy);

    
    // pseudo-rapidity
    // non-elastic events
         histosTH1F["hetanela"]->Fill(pi1eta);
         histosTH1F["hetanela"]->Fill(pi2eta);
         histosTH1F["hetanela"]->Fill(pi3eta);
         histosTH1F["hetanela"]->Fill(pi4eta);
       //
       if(ntrk==4){
         histosTH1F["hetanela4"]->Fill(pi1eta);
         histosTH1F["hetanela4"]->Fill(pi2eta);
         histosTH1F["hetanela4"]->Fill(pi3eta);
         histosTH1F["hetanela4"]->Fill(pi4eta);
       }
       // 
       if(ntrk==4 && totcharge==0){
         histosTH1F["hetanela40"]->Fill(pi1eta);
         histosTH1F["hetanela40"]->Fill(pi2eta);
         histosTH1F["hetanela40"]->Fill(pi3eta);
         histosTH1F["hetanela40"]->Fill(pi4eta);
       }

       
         if(fiducialRegion4 && fiducialRegionPt4){
    	 histosTH2F["h2dimdpy"]->Fill(CMSpy4,TOTEMpy);
	 histosTH1F["hdpy"]->Fill(CMSpy4+TOTEMpy);
	 //
         histosTH1F["hetanela44"]->Fill(pi1eta);
         histosTH1F["hetanela44"]->Fill(pi2eta);
         histosTH1F["hetanela44"]->Fill(pi3eta);
         histosTH1F["hetanela44"]->Fill(pi4eta);
 	 }
         if(fiducialRegion4 && fiducialRegionPt4 && totcharge==0){
	   //histosTH2F["h2dimdpy"]->Fill(CMSpy4,TOTEMpy);
	 histosTH1F["hdpy0"]->Fill(CMSpy4+TOTEMpy);
	 //
         histosTH1F["hetanela440"]->Fill(pi1eta);
         histosTH1F["hetanela440"]->Fill(pi2eta);
         histosTH1F["hetanela440"]->Fill(pi3eta);
         histosTH1F["hetanela440"]->Fill(pi4eta);
	 }

    histosTH2F["h2dimdpxAll"]->Fill(CMSpx4,TOTEMpx);
    histosTH1F["hdpxAll"]->Fill(CMSpx4+TOTEMpx);
    
         if(fiducialRegion4 && fiducialRegionPt4){
	 histosTH2F["h2dimdpx"]->Fill(CMSpx4,TOTEMpx);
	 histosTH1F["hdpx"]->Fill(CMSpx4+TOTEMpx);
	 }
         if(fiducialRegion4 && fiducialRegionPt4 && totcharge==0){
	   //histosTH2F["h2dimdpx"]->Fill(CMSpx4,TOTEMpx);
	 histosTH1F["hdpx0"]->Fill(CMSpx4+TOTEMpx);
	 }

	 // ...hello compiler ...boo!
         ////if(fiducialRegionK4 && fiducialRegionPtK4){
	 ////  std::cout << "." ;
	 ////}

	 // ...checking 4pi-channel mass
	 //std::cout << " mrec  = " << mrec <<std::endl;
	 
    //...Luiz
    histosTH1F["hvtx"]->Fill(isfake);
       //...Luiz
       if(ntrk==4){
	 histosTH1F["hvtx2"]->Fill(isfake);
	 if(fiducialRegion4 && totcharge==0) histosTH1F["hvtx3"]->Fill(isfake);
       }    

       if(ntrk==4){
	 if(totcharge==0){
       histosTH1F["hrimpac1"]->Fill(rimpac1);
       histosTH1F["hrimpac2"]->Fill(rimpac2);
       histosTH1F["hrimpac3"]->Fill(rimpac3);
       histosTH1F["hrimpac4"]->Fill(rimpac4);
 	 }
       }

   //-------------------------------------------------------------------------------

       //...number of vertices accepted
       
       /* Ferenc's recommendation!! commented out!

       //...very important...needed for theVees
       //......not this--> if(nvtx!=0 || nvtx!=1) continue;
       // this:
       if(nvtx!=0){
	 if(nvtx!=1){
	   if(nvtx!=2) continue; //...Ok! it works!
	 }
       }
       */
       
       //if(nvtx!=1) continue;
  
  //--------------------------------------------------------------------------------
  // my stuff

       // M(1,2) M(3,4) M(1,3) M(2,4) M(1,4) M(2,3) 
       double mrecpi1pi2=pi1pi2Rec.M();
       double mrecpi3pi4=pi3pi4Rec.M();
       double mrecpi1pi3=pi1pi3Rec.M();
       double mrecpi2pi4=pi2pi4Rec.M();
       double mrecpi1pi4=pi1pi4Rec.M();
       double mrecpi2pi3=pi2pi3Rec.M();
       //
       double mrec1234=pi1234Rec.M();
       double mrec1324=pi1324Rec.M();
       double mrec1423=pi1423Rec.M();
       //
       double ptpi1pi2=pi1pi2Rec.Pt();
       double ptpi3pi4=pi3pi4Rec.Pt();
       double ptpi1pi3=pi1pi3Rec.Pt();
       double ptpi2pi4=pi2pi4Rec.Pt();
       double ptpi1pi4=pi1pi4Rec.Pt();
       double ptpi2pi3=pi2pi3Rec.Pt();
       //
       double pzpi1pi2=pi1pi2Rec.Pz();
       double pzpi3pi4=pi3pi4Rec.Pz();
       double pzpi1pi3=pi1pi3Rec.Pz();
       double pzpi2pi4=pi2pi4Rec.Pz();
       double pzpi1pi4=pi1pi4Rec.Pz();
       double pzpi2pi3=pi2pi3Rec.Pz();
       //
       double etapi1pi2=pi1pi2Rec.Eta();
       double etapi3pi4=pi3pi4Rec.Eta();
       double etapi1pi3=pi1pi3Rec.Eta();
       double etapi2pi4=pi2pi4Rec.Eta();
       double etapi1pi4=pi1pi4Rec.Eta();
       double etapi2pi3=pi2pi3Rec.Eta();

       /*
       // fixing the mass of the pion pair 
       double eneMK012 = TMath::Sqrt(ptpi1pi2*ptpi1pi2 + pzpi1pi2*pzpi1pi2 + m_k0*m_k0);
       double eneMK034 = TMath::Sqrt(ptpi3pi4*ptpi3pi4 + pzpi3pi4*pzpi3pi4 + m_k0*m_k0);

       histosTH1F["henemk012"]->Fill(eneMK012);
       histosTH1F["henemk034"]->Fill(eneMK034);
       */
       
     //...combining pions and kaons for the event selection type = 11 (one primary & one Vee)
	   /* 
     ...first combining, then select the Q_pairs=0
     pi1pi2 pi3k4
     pi1pi3 pi2k4
     pi2pi3 pi1k4

     pi1pi2 k3pi4
     pi1pi4 k3pi2
     pi2pi4 k3pi1

     pi1k2 pi3pi4
     pi3k2 pi1pi4
     pi4k2 pi1pi3

     k1pi2 pi3pi4
     k1pi3 pi2pi4
     k1pi4 pi2pi3
 	   */

       double mrecpi3k4=pi3k4Rec.M();
       double mrecpi2k4=pi2k4Rec.M(); 
       double mrecpi1k4=pi1k4Rec.M();
       double mreck3pi4=k3pi4Rec.M();
       double mreck3pi2=k3pi2Rec.M();
       double mreck3pi1=k3pi1Rec.M();
       double mrecpi1k2=pi1k2Rec.M();
       double mrecpi3k2=pi3k2Rec.M();
       double mrecpi4k2=pi4k2Rec.M();
       double mreck1pi2=k1pi2Rec.M();
       double mreck1pi3=k1pi3Rec.M();
       double mreck1pi4=k1pi4Rec.M();
       //
       double mrecKpi = 0.0;

       // KsK* or Kpi+pipi type:11 ??????????????????????????...WRONG!
       //k4
       double mrecpi1pi2_pi3k4= pi1pi2Rec.M() + pi3k4Rec.M() ;
       double mrecpi1pi3_pi2k4= pi1pi3Rec.M() + pi2k4Rec.M() ;
       double mrecpi2pi3_pi1k4= pi2pi3Rec.M() + pi1k4Rec.M() ;
       //k3
       double mrecpi1pi2_k3pi4= pi1pi2Rec.M() + k3pi4Rec.M() ;
       double mrecpi1pi4_k3pi2= pi1pi4Rec.M() + k3pi2Rec.M() ;
       double mrecpi2pi4_k3pi1= pi2pi4Rec.M() + k3pi1Rec.M() ;
       //k2
       double mrecpi1k2_pi3pi4= pi1k2Rec.M() + pi3pi4Rec.M() ;
       double mrecpi3k2_pi1pi4= pi3k2Rec.M() + pi1pi4Rec.M() ;
       double mrecpi4k2_pi1pi3= pi4k2Rec.M() + pi1pi3Rec.M() ;
       //k1
       double mreck1pi2_pi3pi4= k1pi2Rec.M() + pi3pi4Rec.M() ;
       double mreck1pi3_pi2pi4= k1pi3Rec.M() + pi2pi4Rec.M() ;
       double mreck1pi4_pi2pi3= k1pi4Rec.M() + pi2pi3Rec.M() ;
       //
       //...better using this one
       double mreck1pi2pi3pi4= k1pi2pi3pi4Rec.M() ;
       double mrecpi1k2pi3pi4= pi1k2pi3pi4Rec.M() ;
       double mrecpi1pi2k3pi4= pi1pi2k3pi4Rec.M() ;
       double mrecpi1pi2pi3k4= pi1pi2pi3k4Rec.M() ;
       
       
       // M(1,2) M(3,4) M(1,3) M(2,4) M(1,4) M(2,3) ...kaons only 
       double mreck1k2=k1k2Rec.M();
       double mreck3k4=k3k4Rec.M();
       double mreck1k3=k1k3Rec.M();
       double mreck2k4=k2k4Rec.M();
       double mreck1k4=k1k4Rec.M();
       double mreck2k3=k2k3Rec.M();
       //
       double ptk1k2=k1k2Rec.Pt();
       double ptk3k4=k3k4Rec.Pt();
       double ptk1k3=k1k3Rec.Pt();
       double ptk2k4=k2k4Rec.Pt();
       double ptk1k4=k1k4Rec.Pt();
       double ptk2k3=k2k3Rec.Pt();
       //
       double etak1k2=k1k2Rec.Eta();
       double etak3k4=k3k4Rec.Eta();
       double etak1k3=k1k3Rec.Eta();
       double etak2k4=k2k4Rec.Eta();
       double etak1k4=k1k4Rec.Eta();
       double etak2k3=k2k3Rec.Eta();


       double mreck2pi3=k2pi3Rec.M();
       double mreck2pi4=k2pi4Rec.M();
       //
       double mrecpi1k3=pi1k3Rec.M();
       //
       double mrecpi2k3=pi2k3Rec.M();


       double mreck1pi2k3pi4=k1pi2k3pi4Rec.M();
       double mreck1pi2pi3k4=k1pi2pi3k4Rec.M();
       double mrecpi1k2pi3k4=pi1k2pi3k4Rec.M();
       double mrecpi1k2k3pi4=pi1k2k3pi4Rec.M();

       double mreck1pi3k2pi4=k1pi3k2pi4Rec.M();
       double mreck1pi3pi2k4=k1pi3pi2k4Rec.M();
       double mrecpi1k3k2pi4=pi1k3k2pi4Rec.M();
       double mrecpi1k3pi2k4=pi1k3pi2k4Rec.M();

       double mreck1pi4k2pi3=k1pi4k2pi3Rec.M();
       double mreck1pi4pi2k3=k1pi4pi2k3Rec.M();
       double mrecpi1k4k2pi3=pi1k4k2pi3Rec.M();
       double mrecpi1k4pi2k3=pi1k4pi2k3Rec.M();

    // acceptance        ...Luiz
       double pt4pi = pipipipiRec.Pt();
    //
    // left t1 versus pt
    histosTH2F["h2dimt1pt"]->Fill(pt4pi,t1);
    // right t2 versus pt
    histosTH2F["h2dimt2pt"]->Fill(pt4pi,t2);
    // t12 versus right pt
    histosTH2F["h2dimt12pt"]->Fill(pt4pi,t12);

 
       //...dE/dx ...1 cut ...no elastic events!
	     int nnn=0;
             for(TrackCollection::const_iterator itTrack = tracks->begin();itTrack != tracks->end();++itTrack) {
	     //	       
	     reco::TrackRef aTrackRefno = reco::TrackRef(tracks, nnn);
             const reco::DeDxData& dedxPIXObjno = dedxPIXTrack[aTrackRefno];	   
             mydedxPIXno = dedxPIXObjno.dEdx();
	     pidno = myparticleType->guess( itTrack->p() , mydedxPIXno );
	     //
             histosTH2F["hdedxPIXno"]->Fill( itTrack->p() , mydedxPIXno );
	     //...only pions
             if(pidno==2){
             histosTH2F["hdedxPIXnopi"]->Fill( itTrack->p() , mydedxPIXno );
	     }
	     //...only kaons
	     if(pidno==3){
             histosTH2F["hdedxPIXnok"]->Fill( itTrack->p() , mydedxPIXno );
	     }
	     //...only protons
	     if(pidno==4){
             histosTH2F["hdedxPIXnop"]->Fill( itTrack->p() , mydedxPIXno );
	     }
	     //...only electrons
	     if(pidno==1){
             histosTH2F["hdedxPIXnoe"]->Fill( itTrack->p() , mydedxPIXno );
	     }
	     //...only unkown
	     if(pidno==0){
             histosTH2F["hdedxPIXnou"]->Fill( itTrack->p() , mydedxPIXno );
	     }
	     //...only K or pi
	     if(pidno==2 || pidno==3){
             histosTH2F["hdedxPIXnokpi"]->Fill( itTrack->p() , mydedxPIXno );
	     }
	     nnn++;
	     } //..end of dE/dx no cuts


  
  //--------------------------------------------------------------------------------
  // my cuts

  if(fiducialRegion4 && fiducialRegionPt4 && allCuts4){

       //...ntrk vs nks
       histosTH2F["hntrknksall4"]->Fill(ntrk,nks);
       //...nvtx vs nks
       histosTH2F["hnvtxnksall4"]->Fill(nvtx,nks);

    
         //...K0 mass window
         // win=20MeV
	 double masslow = 0.49;
         double masshigh = 0.51;
	 // win=40MeV
	 double masslow2 = 0.48;
         double masshigh2 = 0.52;
	 //
         //...rho mass window
	 // win=40MeV
	 double rholow = 0.75;
         double rhohigh = 0.79;
	 // win=70MeV
	 double rholow2 = 0.72;
         double rhohigh2 = 0.79;

	 // pipi
	 double sigpipik0   = 0.031;  // 31 MeV = 1*sigma, 62 MeV window = 2*sigma

	 // rhorho
	 double sigpipirho  = 0.118; // 118 MeV window ...based on Robert's AN sig=59 2sig=118 window=4*sig=236
	 double sigpipirho2 = 0.060; // 60 MeV = 1*sigma, 120 MeV = 2*sigma
	 double sigpipirho3 = 0.010; // 10 MeV = 1*sigma,  20 MeV = 2*sigma
	 double sigpipirho4 = 0.020; // 20 MeV = 1*sigma,  40 MeV = 2*sigma
	 double sigpipirho5 = 0.030; // 30 MeV = 1*sigma,  60 MeV = 2*sigma
	 double sigpipirho6 = 0.040; // 40 MeV = 1*sigma,  80 MeV = 2*sigma
	 
	 // limits for 4-momentum transfer squared
	 double limtA = 0.4  ;
	 double limtB = 0.3  ;
	 double limtC = 0.2  ;
	 double limtD = 0.17 ;
	 double limtE = 0.1  ;
	 double limtF = 0.08 ;
	 double limtG = 0.06 ;
	 double limtH = 0.05 ;

	 // from Breit-Wigner fits
	 //
	 // GAMMA(1370)=0.0 GeV  sig=0.035 GeV
	 // GAMMA(1549)=0.0825 GeV  sig=0.035 GeV
	 // GAMMA(1731)=0.0704 GeV  sig=0.030 GeV
	 //
	 //bool m1370 = TMath::Abs(mrec - 1.370) <= 0.1 ;   // sigma = GeV
	 //bool m1549 = TMath::Abs(mrec - 1.549) <= 0.116 ; // sigma = 0.116 GeV    ???????????????????
	 //bool m1731 = TMath::Abs(mrec - 1.731) <= 0.172 ; // sigma = 0.172 GeV
	 // ...using window=GAMMA
	 bool m1370 = TMath::Abs(mrec - 1.370) <= 0.030   ; // GAMMA/2 = 0.0030   ;...view inspection
	 bool m1549 = TMath::Abs(mrec - 1.549) <= 0.04125 ; // GAMMA/2 = 0.004125 ; sigma = 0.035 GeV    
	 bool m1731 = TMath::Abs(mrec - 1.731) <= 0.0352  ; // GAMMA/2 = 0.00352  ; sigma = 0.030 GeV

	 
	 //...cut 05       ...K0sK0s channel ...selection by mass   WHY NOT TO IMPOSE nvtx=0 for K0sK0s DISPLACED and
	 //                                                                           nvtx=1 for rhorho PROMPT         HERE? !!! DO IT !..........
	 if(totcharge==0){
	  if(nvtx==0){
	   //...only pions
	   if(pidarray[0]==2 && pidarray[1]==2 && pidarray[2]==2 && pidarray[3]==2){
	   //
	     ////if(pidv0array[0]==2 && pidv0array[1]==2 && pidv0array[2]==2 && pidv0array[3]==2){
	     //
	     //1234
	     if(charray[0]+charray[1] == 0)
	          {
		     histosTH1F["hm4rec2OS_pi1pi2t05"]->Fill(mrecpi1pi2);
	             if(mrecpi1pi2 < masshigh && mrecpi1pi2 > masslow){
		     histosTH1F["hm4rec2OS_pi1pi2m05"]->Fill(mrecpi1pi2);
                     histosTH1F["h2OSpt1205"]->Fill(ptpi1pi2);
                     histosTH1F["h2OSeta1205"]->Fill(etapi1pi2);                    // etaCut < 3.0 is not working!!!
                     histosTH2F["h2dim2OSpteta1205"]->Fill(ptpi1pi2,etapi1pi2);
		     }
                     histosTH1F["hm4rec2OS_pi3pi4t05"]->Fill(mrecpi3pi4);
		     if(mrecpi3pi4 < masshigh && mrecpi3pi4 > masslow){
                     histosTH1F["hm4rec2OS_pi3pi4m05"]->Fill(mrecpi3pi4);
                     histosTH1F["h2OSpt3405"]->Fill(ptpi3pi4);
                     histosTH1F["h2OSeta3405"]->Fill(etapi3pi4);
                     histosTH2F["h2dim2OSpteta3405"]->Fill(ptpi3pi4,etapi3pi4);	
		     }
		     /*
	   	     if(pidarray[0]==pidPion && pidarray[1]==pidPion &&
	   	        pidarray[2]==pidPion && pidarray[3]==pidPion){
		       //
		     histosTH1F["hm4rec2OSm1234t05pid"]->Fill(mrec);
		     }
		     */
		     histosTH1F["hm4rec2OSm1234t05"]->Fill(mrec);
               	     histosTH2F["h2dim2OSm12x34t05"]->Fill(mrecpi1pi2,mrecpi3pi4);
	             if(mrecpi1pi2 < masshigh && mrecpi1pi2 > masslow &&
			mrecpi3pi4 < masshigh && mrecpi3pi4 > masslow ){		     
		       	histosTH1F["hm4rec2OSm123405"]->Fill(mrec);
		       	histosTH1F["hm4rec2OSmrec123405"]->Fill(mrec1234);
			// testing mix-up channels
			if(!nks){
		       	histosTH1F["hm4rec2OSm123405nov"]->Fill(mrec);
			}
			if(nks==2){
		       	histosTH1F["hm4rec2OSm123405yesv"]->Fill(mrec);
			}
			if(nks==1){
		       	histosTH1F["hm4rec2OSm123405yes1"]->Fill(mrec);
			}
			//
			if(mrec > 1.50 && mrec < 1.58){			  
			  histosTH1F["hm4rec2OSm123405pi1pt"]->Fill(pi1pt);
			  histosTH1F["hm4rec2OSm123405pi2pt"]->Fill(pi2pt);
			  histosTH1F["hm4rec2OSm123405pi3pt"]->Fill(pi3pt);
			  histosTH1F["hm4rec2OSm123405pi4pt"]->Fill(pi4pt);
			}
               	        histosTH2F["h2dim2OSm12x3405"]->Fill(mrecpi1pi2,mrecpi3pi4);
		     }
	             if(mrecpi1pi2 < masshigh2 && mrecpi1pi2 > masslow2 &&
			mrecpi3pi4 < masshigh2 && mrecpi3pi4 > masslow2 ){
		       	histosTH1F["hm4rec2OSm1234052"]->Fill(mrec);
		     }
		     //...| M(pi+pi-) - M(K0) | < 31 MeV = *sigma, mass window = *sigma = 62 MeV ?????????????
	             if( ( TMath::Abs(mrecpi1pi2 - m_k0) < sigpipik0 ) &&
			 ( TMath::Abs(mrecpi3pi4 - m_k0) < sigpipik0 ) ){
		       	histosTH1F["hm4rec2OSm123405sig"]->Fill(mrec);
		     }
		     // K0
		     // Gamma = 12.02 MeV
		     // sigma = 5.1 MeV
		     //...new
		     // Gamma = 10.36 MeV
		     // sigma = 4.4 MeV
		     //
		     //...final study K0
  		     //...| M(pi+pi-) - M(K0) | < 5.1 MeV ==> mass window: 2*sigma = 10.2 MeV ...prior
  		     //...| M(pi+pi-) - M(K0) | < 4.4 MeV ==> mass window: 2*sigma = 8.8 MeV
	             if( ( TMath::Abs(mrecpi1pi2 - m_k0) < 0.0044 ) &&
			 ( TMath::Abs(mrecpi3pi4 - m_k0) < 0.0044 ) ){
		       	histosTH1F["hm4rec2OSm123405win10"]->Fill(mrec);
		     }
		     //...| M(pi+pi-) - M(K0) | < 10.2 MeV ==> mass window: 4*sigma = 20.4 MeV ...prior
		     //...| M(pi+pi-) - M(K0) | < 8.8 MeV ==> mass window: 4*sigma = 17.6 MeV
	             if( ( TMath::Abs(mrecpi1pi2 - m_k0) < 0.0088 ) &&
			 ( TMath::Abs(mrecpi3pi4 - m_k0) < 0.0088 ) ){
		       	histosTH1F["hm4rec2OSm123405win20"]->Fill(mrec);
		     }
		     //...| M(pi+pi-) - M(K0) | < 13.2 MeV ==> mass window: 6*sigma = 26.4 MeV
	             if( ( TMath::Abs(mrecpi1pi2 - m_k0) < 0.0132 ) &&
			 ( TMath::Abs(mrecpi3pi4 - m_k0) < 0.0132 ) ){
		       	histosTH1F["hm4rec2OSm123405win30"]->Fill(mrec);
		     }
		     
		     /*
		     //
		     // rho
		     // Gamma = 87.46 MeV
		     // sigma = 37.14 MeV
	             //
		     //...final study rho
    		     //...| M(pi+pi-) - M(rho) | < 37.14 MeV ==> mass window: 2*sigma = 74.3 MeV
	             if( ( TMath::Abs(mrecpi1pi2 - m_rho) < 0.03714 ) &&
			 ( TMath::Abs(mrecpi3pi4 - m_rho) < 0.03714 ) ){
		       	histosTH1F["hm4rec2OSr123405win74"]->Fill(mrec);
		     }
    		     //...| M(pi+pi-) - M(rho) | < 74.3 MeV ==> mass window: 4*sigma = 148.6 MeV
	             if( ( TMath::Abs(mrecpi1pi2 - m_rho) < 0.0743 ) &&
			 ( TMath::Abs(mrecpi3pi4 - m_rho) < 0.0743 ) ){
		       	histosTH1F["hm4rec2OSr123405win148"]->Fill(mrec);
		     }
    		     //...| M(pi+pi-) - M(rho) | < 18.6 MeV ==> mass window: sigma = 37.14 MeV
	             if( ( TMath::Abs(mrecpi1pi2 - m_rho) < 0.0186 ) &&
			 ( TMath::Abs(mrecpi3pi4 - m_rho) < 0.0186 ) ){
		       	histosTH1F["hm4rec2OSr123405winhalf"]->Fill(mrec);
		     }
		     // ...delta_phi < 90 degrees: parallel protons TTBB
		     // pT < 0.8 GeV/c
		     if( pipipipiRec.Pt() < 0.8 ){
	             if( ( TMath::Abs(mrecpi1pi2 - m_rho) < 0.0186 ) &&
			 ( TMath::Abs(mrecpi3pi4 - m_rho) < 0.0186 ) ){
		       	histosTH1F["hm4rec2OSr123405winpt"]->Fill(mrec);
		     }}
		     //
		     //...rho
		     if(mrecpi1pi2 < rhohigh && mrecpi1pi2 > rholow &&
			mrecpi3pi4 < rhohigh && mrecpi3pi4 > rholow ){		     
		       	histosTH1F["hm4rec2OSr123405"]->Fill(mrec);
		     }
		     //...| M(pi+pi-) - M(rho) | < 35 MeV = 1*sigma, mass window = 2*sigma = 70 MeV
		     if(mrecpi1pi2 < rhohigh2 && mrecpi1pi2 > rholow2 &&
			mrecpi3pi4 < rhohigh2 && mrecpi3pi4 > rholow2 ){		     
		       	histosTH1F["hm4rec2OSr1234052"]->Fill(mrec);
		     }
		     //...| M(pi+pi-) - M(rho) | < 118 MeV = 1*sigma, mass window = 2*sigma = 236 MeV
	             if( ( TMath::Abs(mrecpi1pi2 - m_rho) < sigpipirho ) &&
			 ( TMath::Abs(mrecpi3pi4 - m_rho) < sigpipirho ) ){
		       	histosTH1F["hm4rec2OSr123405sig"]->Fill(mrec);
		     }
		     //...| M(pi+pi-) - M(rho) | < 60 MeV = 1*sigma, mass window = 2*sigma = 120 MeV
	             if( ( TMath::Abs(mrecpi1pi2 - m_rho) < sigpipirho2 ) &&
			 ( TMath::Abs(mrecpi3pi4 - m_rho) < sigpipirho2 ) ){
		       	histosTH1F["hm4rec2OSr123405sig2"]->Fill(mrec);
		     }
		     //...| M(pi+pi-) - M(rho) | < 10 MeV = 1*sigma, mass window = 2*sigma = 20 MeV
	             if( ( TMath::Abs(mrecpi1pi2 - m_rho) < sigpipirho3 ) &&
			 ( TMath::Abs(mrecpi3pi4 - m_rho) < sigpipirho3 ) ){
		       	histosTH1F["hm4rec2OSr123405sig3"]->Fill(mrec);
		     }
		     //...| M(pi+pi-) - M(rho) | < 20 MeV = 1*sigma, mass window = 2*sigma = 40 MeV
	             if( ( TMath::Abs(mrecpi1pi2 - m_rho) < sigpipirho4 ) &&
			 ( TMath::Abs(mrecpi3pi4 - m_rho) < sigpipirho4 ) ){
		       	histosTH1F["hm4rec2OSr123405sig4"]->Fill(mrec);
		     }
		     //...| M(pi+pi-) - M(rho) | < 30 MeV = 1*sigma, mass window = 2*sigma = 60 MeV
	             if( ( TMath::Abs(mrecpi1pi2 - m_rho) < sigpipirho5 ) &&
			 ( TMath::Abs(mrecpi3pi4 - m_rho) < sigpipirho5 ) ){
		       	histosTH1F["hm4rec2OSr123405sig5"]->Fill(mrec);
		     }
		     //...| M(pi+pi-) - M(rho) | < 40 MeV = 1*sigma, mass window = 2*sigma = 80 MeV
	             if( ( TMath::Abs(mrecpi1pi2 - m_rho) < sigpipirho6 ) &&
			 ( TMath::Abs(mrecpi3pi4 - m_rho) < sigpipirho6 ) ){
		       	histosTH1F["hm4rec2OSr123405sig6"]->Fill(mrec);
		     }
		     */
    		  }//...end Q_pair=0 ...1234
	   //
	   //1324
	   if(charray[0]+charray[2] == 0)
		  {
		     histosTH1F["hm4rec2OS_pi1pi3t05"]->Fill(mrecpi1pi3);
	             if(mrecpi1pi3 < masshigh && mrecpi1pi3 > masslow){
		     histosTH1F["hm4rec2OS_pi1pi3m05"]->Fill(mrecpi1pi3);
                     histosTH1F["h2OSpt1305"]->Fill(ptpi1pi3);
                     histosTH1F["h2OSeta1305"]->Fill(etapi1pi3);
                     histosTH2F["h2dim2OSpteta1305"]->Fill(ptpi1pi3,etapi1pi3);
		     }
		     histosTH1F["hm4rec2OS_pi2pi4t05"]->Fill(mrecpi2pi4);
		     if(mrecpi2pi4 < masshigh && mrecpi2pi4 > masslow){
		     histosTH1F["hm4rec2OS_pi2pi4m05"]->Fill(mrecpi2pi4);
                     histosTH1F["h2OSpt2405"]->Fill(ptpi2pi4);
                     histosTH1F["h2OSeta2405"]->Fill(etapi2pi4);
                     histosTH2F["h2dim2OSpteta2405"]->Fill(ptpi2pi4,etapi2pi4);		     
 		     }
		     /*
	   	     if(pidarray[0]==pidPion && pidarray[1]==pidPion &&
	   	        pidarray[2]==pidPion && pidarray[3]==pidPion){
		       //
		     histosTH1F["hm4rec2OSm1324t05pid"]->Fill(mrec);
		     }
		     */
		     histosTH1F["hm4rec2OSm1324t05"]->Fill(mrec);
    		     histosTH2F["h2dim2OSm13x24t05"]->Fill(mrecpi1pi3,mrecpi2pi4);
	             if(mrecpi1pi3 < masshigh && mrecpi1pi3 > masslow &&
			mrecpi2pi4 < masshigh && mrecpi2pi4 > masslow ){		     
		        histosTH1F["hm4rec2OSm132405"]->Fill(mrec);
		       	histosTH1F["hm4rec2OSmrec132405"]->Fill(mrec1324);
			// testing mix-up channels
			if(!nks){
		       	histosTH1F["hm4rec2OSm132405nov"]->Fill(mrec);
			}
			if(nks==2){
		       	histosTH1F["hm4rec2OSm132405yesv"]->Fill(mrec);
			}
			if(nks==1){
		       	histosTH1F["hm4rec2OSm132405yes1"]->Fill(mrec);
			}
			//
 			if(mrec > 1.50 && mrec < 1.58){			  
			  histosTH1F["hm4rec2OSm132405pi1pt"]->Fill(pi1pt);
			  histosTH1F["hm4rec2OSm132405pi2pt"]->Fill(pi2pt);
			  histosTH1F["hm4rec2OSm132405pi3pt"]->Fill(pi3pt);
			  histosTH1F["hm4rec2OSm132405pi4pt"]->Fill(pi4pt);
			}
   		     histosTH2F["h2dim2OSm13x2405"]->Fill(mrecpi1pi3,mrecpi2pi4);
		     }
	             if(mrecpi1pi3 < masshigh2 && mrecpi1pi3 > masslow2 &&
			mrecpi2pi4 < masshigh2 && mrecpi2pi4 > masslow2 ){
		        histosTH1F["hm4rec2OSm1324052"]->Fill(mrec);
		     }
		     //...| M(pi+pi-) - M(K0) | < 31 MeV = 1*sigma, mass window = 2*sigma = 62 MeV
	             if( ( TMath::Abs(mrecpi1pi3 - m_k0) < sigpipik0 ) &&
			 ( TMath::Abs(mrecpi2pi4 - m_k0) < sigpipik0 ) ){
		       	histosTH1F["hm4rec2OSm132405sig"]->Fill(mrec);
		     }
		     // K0
		     // Gamma = 12.02 MeV
		     // sigma = 5.1 MeV
	             //
		     //...final study K0
 		     //...| M(pi+pi-) - M(K0) | < 5.1 MeV ==> mass window: 2*sigma = 10.2 MeV ...prior
  		     //...| M(pi+pi-) - M(K0) | < 4.4 MeV ==> mass window: 2*sigma = 8.8 MeV
	             if( ( TMath::Abs(mrecpi1pi3 - m_k0) < 0.0044 ) &&
			 ( TMath::Abs(mrecpi2pi4 - m_k0) < 0.0044 ) ){
		       	histosTH1F["hm4rec2OSm132405win10"]->Fill(mrec);
		     }
		     //...| M(pi+pi-) - M(K0) | < 10.2 MeV ==> mass window: 4*sigma = 20.4 MeV ...prior
		     //...| M(pi+pi-) - M(K0) | < 8.8 MeV ==> mass window: 4*sigma = 17.6 MeV
	             if( ( TMath::Abs(mrecpi1pi3 - m_k0) < 0.0088 ) &&
			 ( TMath::Abs(mrecpi2pi4 - m_k0) < 0.0088 ) ){
		       	histosTH1F["hm4rec2OSm132405win20"]->Fill(mrec);
		     }
		     //...| M(pi+pi-) - M(K0) | < 13.2 MeV ==> mass window: 6*sigma = 26.4 MeV
	             if( ( TMath::Abs(mrecpi1pi3 - m_k0) < 0.0132 ) &&
			 ( TMath::Abs(mrecpi2pi4 - m_k0) < 0.0132 ) ){
		       	histosTH1F["hm4rec2OSm132405win30"]->Fill(mrec);
		     }
		     
		     /*
		     //
		     // rho
		     // Gamma = 87.46 MeV
		     // sigma = 37.14 MeV
	             //
		     //...final study rho
    		     //...| M(pi+pi-) - M(rho) | < 37.14 MeV ==> mass window: 2*sigma = 74.3 MeV
	             if( ( TMath::Abs(mrecpi1pi3 - m_rho) < 0.03714 ) &&
			 ( TMath::Abs(mrecpi2pi4 - m_rho) < 0.03714 ) ){
		       	histosTH1F["hm4rec2OSr132405win74"]->Fill(mrec);
		     }
    		     //...| M(pi+pi-) - M(rho) | < 74.3 MeV ==> mass window: 4*sigma = 148.6 MeV
	             if( ( TMath::Abs(mrecpi1pi3 - m_rho) < 0.0743 ) &&
			 ( TMath::Abs(mrecpi2pi4 - m_rho) < 0.0743 ) ){
		       	histosTH1F["hm4rec2OSr132405win148"]->Fill(mrec);
		     }
    		     //...| M(pi+pi-) - M(rho) | < 18.6 MeV ==> mass window: sigma = 37.14 MeV
	             if( ( TMath::Abs(mrecpi1pi3 - m_rho) < 0.0186 ) &&
			 ( TMath::Abs(mrecpi2pi4 - m_rho) < 0.0186 ) ){
		       	histosTH1F["hm4rec2OSr132405winhalf"]->Fill(mrec);
		     }
		     // ...delta_phi < 90 degrees: parallel protons TTBB
		     // pT < 0.8 GeV/c
		     if( pipipipiRec.Pt() < 0.8 ){
	             if( ( TMath::Abs(mrecpi1pi3 - m_rho) < 0.0186 ) &&
			 ( TMath::Abs(mrecpi2pi4 - m_rho) < 0.0186 ) ){
		       	histosTH1F["hm4rec2OSr132405winpt"]->Fill(mrec);
		     }}
		     //
		     //...rho
	             if(mrecpi1pi3 < rhohigh && mrecpi1pi3 > rholow &&
			mrecpi2pi4 < rhohigh && mrecpi2pi4 > rholow ){		     
		        histosTH1F["hm4rec2OSr132405"]->Fill(mrec);
		     }		     
	             if(mrecpi1pi3 < rhohigh2 && mrecpi1pi3 > rholow2 &&
			mrecpi2pi4 < rhohigh2 && mrecpi2pi4 > rholow2 ){		     
		        histosTH1F["hm4rec2OSr1324052"]->Fill(mrec);
		     }
    		     //...| M(pi+pi-) - M(rho) | < 35 MeV = 1*sigma, mass window = 2*sigma = 70 MeV
    		     //...| M(pi+pi-) - M(rho) | < 138 MeV = 1*sigma, mass window = 2*sigma = 276 MeV
	             if( ( TMath::Abs(mrecpi1pi3 - m_rho) < sigpipirho ) &&
			 ( TMath::Abs(mrecpi2pi4 - m_rho) < sigpipirho ) ){
		       	histosTH1F["hm4rec2OSr132405sig"]->Fill(mrec);
		     }
    		     //...| M(pi+pi-) - M(rho) | < 60 MeV = 1*sigma, mass window = 2*sigma = 120 MeV
	             if( ( TMath::Abs(mrecpi1pi3 - m_rho) < sigpipirho2 ) &&
			 ( TMath::Abs(mrecpi2pi4 - m_rho) < sigpipirho2 ) ){
		       	histosTH1F["hm4rec2OSr132405sig2"]->Fill(mrec);
		     }
    		     //...| M(pi+pi-) - M(rho) | < 10 MeV = 1*sigma, mass window = 2*sigma = 20 MeV
	             if( ( TMath::Abs(mrecpi1pi3 - m_rho) < sigpipirho3 ) &&
			 ( TMath::Abs(mrecpi2pi4 - m_rho) < sigpipirho3 ) ){
		       	histosTH1F["hm4rec2OSr132405sig3"]->Fill(mrec);
		     }
    		     //...| M(pi+pi-) - M(rho) | < 20 MeV = 1*sigma, mass window = 2*sigma = 40 MeV
	             if( ( TMath::Abs(mrecpi1pi3 - m_rho) < sigpipirho4 ) &&
			 ( TMath::Abs(mrecpi2pi4 - m_rho) < sigpipirho4 ) ){
		       	histosTH1F["hm4rec2OSr132405sig4"]->Fill(mrec);
		     }
    		     //...| M(pi+pi-) - M(rho) | < 30 MeV = 1*sigma, mass window = 2*sigma = 60 MeV
	             if( ( TMath::Abs(mrecpi1pi3 - m_rho) < sigpipirho5 ) &&
			 ( TMath::Abs(mrecpi2pi4 - m_rho) < sigpipirho5 ) ){
		       	histosTH1F["hm4rec2OSr132405sig5"]->Fill(mrec);
		     }
    		     //...| M(pi+pi-) - M(rho) | < 40 MeV = 1*sigma, mass window = 2*sigma = 80 MeV
	             if( ( TMath::Abs(mrecpi1pi3 - m_rho) < sigpipirho6 ) &&
			 ( TMath::Abs(mrecpi2pi4 - m_rho) < sigpipirho6 ) ){
		       	histosTH1F["hm4rec2OSr132405sig6"]->Fill(mrec);
		     }
		     */
    		  }//...end Q_pair=0 ...1324
	   //
	   //1423
	   if(charray[0]+charray[3] == 0)
	          {	    
		     histosTH1F["hm4rec2OS_pi1pi4t05"]->Fill(mrecpi1pi4);
	             if(mrecpi1pi4 < masshigh && mrecpi1pi4 > masslow){
		     histosTH1F["hm4rec2OS_pi1pi4m05"]->Fill(mrecpi1pi4);
                     histosTH1F["h2OSpt1405"]->Fill(ptpi1pi4);
                     histosTH1F["h2OSeta1405"]->Fill(etapi1pi4);                // ... etaCut < 3.0 is not working!!! why ????????????
                     histosTH2F["h2dim2OSpteta1405"]->Fill(ptpi1pi4,etapi1pi4);
 		     }
		     histosTH1F["hm4rec2OS_pi2pi3t05"]->Fill(mrecpi2pi3);
		     if(mrecpi2pi3 < masshigh && mrecpi2pi3 > masslow){
		     histosTH1F["hm4rec2OS_pi2pi3m05"]->Fill(mrecpi2pi3);
                     histosTH1F["h2OSpt2305"]->Fill(ptpi2pi3);
                     histosTH1F["h2OSeta2305"]->Fill(etapi2pi3);                // .. Why hasn't it cut < 3.0 ?????????????????????????????
                     histosTH2F["h2dim2OSpteta2305"]->Fill(ptpi2pi3,etapi2pi3);		     
		     }
		     /*
	   	     if(pidarray[0]==pidPion && pidarray[1]==pidPion &&
	   	        pidarray[2]==pidPion && pidarray[3]==pidPion){
		       //
		     histosTH1F["hm4rec2OSm1423t05pid"]->Fill(mrec);
		     }
		     */
		     histosTH1F["hm4rec2OSm1423t05"]->Fill(mrec);
                     histosTH2F["h2dim2OSm14x23t05"]->Fill(mrecpi1pi4,mrecpi2pi3);
	             if(mrecpi1pi4 < masshigh && mrecpi1pi4 > masslow &&
			mrecpi2pi3 < masshigh && mrecpi2pi3 > masslow ){		     
		        histosTH1F["hm4rec2OSm142305"]->Fill(mrec);
		       	histosTH1F["hm4rec2OSmrec142305"]->Fill(mrec1423);
			// testing mix-up channels
			if(!nks){
		       	histosTH1F["hm4rec2OSm142305nov"]->Fill(mrec);
			}
			if(nks==2){
		       	histosTH1F["hm4rec2OSm142305yesv"]->Fill(mrec);
			}
			if(nks==1){
		       	histosTH1F["hm4rec2OSm142305yes1"]->Fill(mrec);
			}
			//
			if(mrec > 1.50 && mrec < 1.58){			  
			  histosTH1F["hm4rec2OSm142305pi1pt"]->Fill(pi1pt);
			  histosTH1F["hm4rec2OSm142305pi2pt"]->Fill(pi2pt);
			  histosTH1F["hm4rec2OSm142305pi3pt"]->Fill(pi3pt);
			  histosTH1F["hm4rec2OSm142305pi4pt"]->Fill(pi4pt);
			}		
                     histosTH2F["h2dim2OSm14x2305"]->Fill(mrecpi1pi4,mrecpi2pi3);
		     }
	             if(mrecpi1pi4 < masshigh2 && mrecpi1pi4 > masslow2 &&
			mrecpi2pi3 < masshigh2 && mrecpi2pi3 > masslow2 ){
		        histosTH1F["hm4rec2OSm1423052"]->Fill(mrec);
		     }
		     //...| M(pi+pi-) - M(K0) | < 31 MeV = 1*sigma, mass window = 2*sigma = 62 MeV
	             if( ( TMath::Abs(mrecpi1pi4 - m_k0) < sigpipik0 ) &&
			 ( TMath::Abs(mrecpi2pi3 - m_k0) < sigpipik0 ) ){
		       	histosTH1F["hm4rec2OSm142305sig"]->Fill(mrec);
		     }
		     
	   // ----------- K0 ------------study
	   // 1423 for M(K0) shape study ...based on inspection of M(14)xM(23).
	   // | M(pi2pi3) - M(K0) | < 50 MeV win=100MeV
	             if( TMath::Abs(mrecpi2pi3 - m_k0) < 0.050 ){
		     histosTH1F["hm4rec2OS14k0"]->Fill(mrecpi1pi4);
 		     }
	   // 1423 for M(K0) shape study ...based on inspection of M(14)xM(23). ...tasting!
	   // | M(pi2pi3) - M(K0) | < 50 MeV win=100MeV
	             if( (TMath::Abs(mrecpi2pi3 - m_k0) < 0.050) && (TMath::Abs(mrecpi1pi4 - m_k0) < 0.050) ){
		     histosTH1F["hm4rec2OS1423k0"]->Fill(mrec);
 		     }
     
		     // K0
		     // Gamma = 12.02 MeV
		     // sigma = 5.1 MeV
		     //
		     //...final study K0
 		     //...| M(pi+pi-) - M(K0) | < 5.1 MeV ==> mass window: 2*sigma = 10.2 MeV ...prior
  		     //...| M(pi+pi-) - M(K0) | < 4.4 MeV ==> mass window: 2*sigma = 8.8 MeV
	             if( ( TMath::Abs(mrecpi1pi4 - m_k0) < 0.0044 ) &&
			 ( TMath::Abs(mrecpi2pi3 - m_k0) < 0.0044 ) ){
		       	histosTH1F["hm4rec2OSm142305win10"]->Fill(mrec);
		     }
		     //...| M(pi+pi-) - M(K0) | < 10.2 MeV ==> mass window: 4*sigma = 20.4 MeV ...prior
		     //...| M(pi+pi-) - M(K0) | < 8.8 MeV ==> mass window: 4*sigma = 17.6 MeV
	             if( ( TMath::Abs(mrecpi1pi4 - m_k0) < 0.0088 ) &&
			 ( TMath::Abs(mrecpi2pi3 - m_k0) < 0.0088 ) ){
		       	histosTH1F["hm4rec2OSm142305win20"]->Fill(mrec);
		     }
		     //...| M(pi+pi-) - M(K0) | < 13.2 MeV ==> mass window: 6*sigma = 26.4 MeV
	             if( ( TMath::Abs(mrecpi1pi4 - m_k0) < 0.0132 ) &&
			 ( TMath::Abs(mrecpi2pi3 - m_k0) < 0.0132 ) ){
		       	histosTH1F["hm4rec2OSm142305win30"]->Fill(mrec);
		     }

		     
		     /*		          
	   // ----------- rho -----------study	          
	   // 1423 for M(rho) shape study ...beased on inspection of M(14)xM(23).
	   // | M(pi2pi3) - M(rho) | < 150 MeV win=300MeV
	             if( TMath::Abs(mrecpi2pi3 - m_rho) < 0.150 ){
		     histosTH1F["hm4rec2OS14rho"]->Fill(mrecpi1pi4);
 		     }

		     // rho
		     // Gamma = 87.46 MeV
		     // sigma = 37.14 MeV
	             //
		     //...final study rho
    		     //...| M(pi+pi-) - M(rho) | < 37.14 MeV ==> mass window: 2*sigma = 74.3 MeV
	             if( ( TMath::Abs(mrecpi1pi4 - m_rho) < 0.03714 ) &&
			 ( TMath::Abs(mrecpi2pi3 - m_rho) < 0.03714 ) ){
		       	histosTH1F["hm4rec2OSr142305win74"]->Fill(mrec);
		     }
    		     //...| M(pi+pi-) - M(rho) | < 74.3 MeV ==> mass window: 4*sigma = 148.6 MeV
	             if( ( TMath::Abs(mrecpi1pi4 - m_rho) < 0.0743 ) &&
			 ( TMath::Abs(mrecpi2pi3 - m_rho) < 0.0743 ) ){
		       	histosTH1F["hm4rec2OSr142305win148"]->Fill(mrec);
		     }
    		     //...| M(pi+pi-) - M(rho) | < 18.6 MeV ==> mass window: sigma = 37.14 MeV
	             if( ( TMath::Abs(mrecpi1pi4 - m_rho) < 0.0186 ) &&
			 ( TMath::Abs(mrecpi2pi3 - m_rho) < 0.0186 ) ){
		       	histosTH1F["hm4rec2OSr142305winhalf"]->Fill(mrec);
		     }
		     // ...delta_phi < 90 degrees: parallel protons TTBB
		     // pT < 0.8 GeV/c
		     if( pipipipiRec.Pt() < 0.8 ){
	             if( ( TMath::Abs(mrecpi1pi4 - m_rho) < 0.0186 ) &&
			 ( TMath::Abs(mrecpi2pi3 - m_rho) < 0.0186 ) ){
		       	histosTH1F["hm4rec2OSr142305winpt"]->Fill(mrec);
		     }}
		     //
		     //...rho
		     // win=40MeV
		     if(mrecpi1pi4 < rhohigh && mrecpi1pi4 > rholow &&
			mrecpi2pi3 < rhohigh && mrecpi2pi3 > rholow ){		     
		        histosTH1F["hm4rec2OSr142305"]->Fill(mrec);
		     }
		     // win=70MeV asymetric
		     if(mrecpi1pi4 < rhohigh2 && mrecpi1pi4 > rholow2 &&
			mrecpi2pi3 < rhohigh2 && mrecpi2pi3 > rholow2 ){		     
		        histosTH1F["hm4rec2OSr1423052"]->Fill(mrec);
		     }
    		     //...| M(pi+pi-) - M(rho) | < 118 MeV = 1*sigma, mass window = 2*sigma =  MeV
	             if( ( TMath::Abs(mrecpi1pi4 - m_rho) < sigpipirho ) &&
			 ( TMath::Abs(mrecpi2pi3 - m_rho) < sigpipirho ) ){
		       	histosTH1F["hm4rec2OSr142305sig"]->Fill(mrec);
		     }
		     //...| M(pi+pi-) - M(rho) | < 60 MeV = 1*sigma, mass window = 2*sigma = 120 MeV
	             if( ( TMath::Abs(mrecpi1pi4 - m_rho) < sigpipirho2 ) &&
			 ( TMath::Abs(mrecpi2pi3 - m_rho) < sigpipirho2 ) ){
		       	histosTH1F["hm4rec2OSr142305sig2"]->Fill(mrec);
		     }
		     //...| M(pi+pi-) - M(rho) | < 10 MeV = 1*sigma, mass window = 2*sigma = 20 MeV
	             if( ( TMath::Abs(mrecpi1pi4 - m_rho) < sigpipirho3 ) &&
			 ( TMath::Abs(mrecpi2pi3 - m_rho) < sigpipirho3 ) ){
		       	histosTH1F["hm4rec2OSr142305sig3"]->Fill(mrec);
		     }
		     //...| M(pi+pi-) - M(rho) | < 20 MeV = 1*sigma, mass window = 2*sigma = 40 MeV
	             if( ( TMath::Abs(mrecpi1pi4 - m_rho) < sigpipirho4 ) &&
			 ( TMath::Abs(mrecpi2pi3 - m_rho) < sigpipirho4 ) ){
		       	histosTH1F["hm4rec2OSr142305sig4"]->Fill(mrec);
		     }
		     //...| M(pi+pi-) - M(rho) | < 30 MeV = 1*sigma, mass window = 2*sigma = 60 MeV
	             if( ( TMath::Abs(mrecpi1pi4 - m_rho) < sigpipirho5 ) &&
			 ( TMath::Abs(mrecpi2pi3 - m_rho) < sigpipirho5 ) ){
		       	histosTH1F["hm4rec2OSr142305sig5"]->Fill(mrec);
		     }
		     //...| M(pi+pi-) - M(rho) | < 40 MeV = 1*sigma, mass window = 2*sigma = 80 MeV
	             if( ( TMath::Abs(mrecpi1pi4 - m_rho) < sigpipirho6 ) &&
			 ( TMath::Abs(mrecpi2pi3 - m_rho) < sigpipirho6 ) ){
		       	histosTH1F["hm4rec2OSr142305sig6"]->Fill(mrec);
		     }
		     */
		  }//...end Q_pair=0 ...1423
	   
	    }//...end of PID pions
	  }//...end of nvtx=0

	  
	  //                            ... rho rho ... 
	  if(nvtx==1){
	   //...only pions
	   if(pidarray[0]==2 && pidarray[1]==2 && pidarray[2]==2 && pidarray[3]==2){
	     //
	     //1234
	     if(charray[0]+charray[1] == 0)
	          {
		     histosTH1F["hm4rec2OSr1234t05"]->Fill(mrec);
               	     histosTH2F["h2dim2OSr12x34t05"]->Fill(mrecpi1pi2,mrecpi3pi4);
		     //
		     // rho
		     // Gamma = 87.46 MeV
		     // sigma = 37.14 MeV
	             //
		     //...final study rho
    		     //...| M(pi+pi-) - M(rho) | < 37.14 MeV ==> mass window: 2*sigma = 74.3 MeV
	             if( ( TMath::Abs(mrecpi1pi2 - m_rho) < 0.03714 ) &&
			 ( TMath::Abs(mrecpi3pi4 - m_rho) < 0.03714 ) ){
		       	histosTH1F["hm4rec2OSr123405win74"]->Fill(mrec);
		     }
    		     //...| M(pi+pi-) - M(rho) | < 74.3 MeV ==> mass window: 4*sigma = 148.6 MeV
	             if( ( TMath::Abs(mrecpi1pi2 - m_rho) < 0.0743 ) &&
			 ( TMath::Abs(mrecpi3pi4 - m_rho) < 0.0743 ) ){
		       	histosTH1F["hm4rec2OSr123405win148"]->Fill(mrec);
		     }
    		     //...| M(pi+pi-) - M(rho) | < 18.6 MeV ==> mass window: sigma = 37.14 MeV
	             if( ( TMath::Abs(mrecpi1pi2 - m_rho) < 0.0186 ) &&
			 ( TMath::Abs(mrecpi3pi4 - m_rho) < 0.0186 ) ){
		       	histosTH1F["hm4rec2OSr123405winhalf"]->Fill(mrec);
		     }
		     // ...delta_phi < 90 degrees: parallel protons TTBB
		     // pT < 0.8 GeV/c
		     if( pipipipiRec.Pt() < 0.8 ){
	             if( ( TMath::Abs(mrecpi1pi2 - m_rho) < 0.0186 ) &&
			 ( TMath::Abs(mrecpi3pi4 - m_rho) < 0.0186 ) ){
		       	histosTH1F["hm4rec2OSr123405winpt"]->Fill(mrec);
		     }}
		     //
		     //...rho
		     if(mrecpi1pi2 < rhohigh && mrecpi1pi2 > rholow &&
			mrecpi3pi4 < rhohigh && mrecpi3pi4 > rholow ){		     
		       	histosTH1F["hm4rec2OSr123405"]->Fill(mrec);
		     }
		     //...| M(pi+pi-) - M(rho) | < 35 MeV = 1*sigma, mass window = 2*sigma = 70 MeV
		     if(mrecpi1pi2 < rhohigh2 && mrecpi1pi2 > rholow2 &&
			mrecpi3pi4 < rhohigh2 && mrecpi3pi4 > rholow2 ){		     
		       	histosTH1F["hm4rec2OSr1234052"]->Fill(mrec);
		     }
		     //...| M(pi+pi-) - M(rho) | < 118 MeV = 1*sigma, mass window = 2*sigma = 236 MeV
	             if( ( TMath::Abs(mrecpi1pi2 - m_rho) < sigpipirho ) &&
			 ( TMath::Abs(mrecpi3pi4 - m_rho) < sigpipirho ) ){
		       	histosTH1F["hm4rec2OSr123405sig"]->Fill(mrec);
		     }
		     //...| M(pi+pi-) - M(rho) | < 60 MeV = 1*sigma, mass window = 2*sigma = 120 MeV
	             if( ( TMath::Abs(mrecpi1pi2 - m_rho) < sigpipirho2 ) &&
			 ( TMath::Abs(mrecpi3pi4 - m_rho) < sigpipirho2 ) ){
		       	histosTH1F["hm4rec2OSr123405sig2"]->Fill(mrec);
		     }
		     //...| M(pi+pi-) - M(rho) | < 10 MeV = 1*sigma, mass window = 2*sigma = 20 MeV
	             if( ( TMath::Abs(mrecpi1pi2 - m_rho) < sigpipirho3 ) &&
			 ( TMath::Abs(mrecpi3pi4 - m_rho) < sigpipirho3 ) ){
		       	histosTH1F["hm4rec2OSr123405sig3"]->Fill(mrec);
		     }
		     //...| M(pi+pi-) - M(rho) | < 20 MeV = 1*sigma, mass window = 2*sigma = 40 MeV
	             if( ( TMath::Abs(mrecpi1pi2 - m_rho) < sigpipirho4 ) &&
			 ( TMath::Abs(mrecpi3pi4 - m_rho) < sigpipirho4 ) ){
		       	histosTH1F["hm4rec2OSr123405sig4"]->Fill(mrec);
		     }
		     //...| M(pi+pi-) - M(rho) | < 30 MeV = 1*sigma, mass window = 2*sigma = 60 MeV
	             if( ( TMath::Abs(mrecpi1pi2 - m_rho) < sigpipirho5 ) &&
			 ( TMath::Abs(mrecpi3pi4 - m_rho) < sigpipirho5 ) ){
		       	histosTH1F["hm4rec2OSr123405sig5"]->Fill(mrec);
		     }
		     //...| M(pi+pi-) - M(rho) | < 40 MeV = 1*sigma, mass window = 2*sigma = 80 MeV
	             if( ( TMath::Abs(mrecpi1pi2 - m_rho) < sigpipirho6 ) &&
			 ( TMath::Abs(mrecpi3pi4 - m_rho) < sigpipirho6 ) ){
		       	histosTH1F["hm4rec2OSr123405sig6"]->Fill(mrec);
		     }
		   }//...end of Qpair=0...1234

	   //
	   //1324
	     if(charray[0]+charray[2] == 0)
		  {
		     histosTH1F["hm4rec2OSr1324t05"]->Fill(mrec);
               	     histosTH2F["h2dim2OSr13x24t05"]->Fill(mrecpi1pi3,mrecpi2pi4);
		     //
		     // rho
		     // Gamma = 87.46 MeV
		     // sigma = 37.14 MeV
	             //
		     //...final study rho
    		     //...| M(pi+pi-) - M(rho) | < 37.14 MeV ==> mass window: 2*sigma = 74.3 MeV
	             if( ( TMath::Abs(mrecpi1pi3 - m_rho) < 0.03714 ) &&
			 ( TMath::Abs(mrecpi2pi4 - m_rho) < 0.03714 ) ){
		       	histosTH1F["hm4rec2OSr132405win74"]->Fill(mrec);
		     }
    		     //...| M(pi+pi-) - M(rho) | < 74.3 MeV ==> mass window: 4*sigma = 148.6 MeV
	             if( ( TMath::Abs(mrecpi1pi3 - m_rho) < 0.0743 ) &&
			 ( TMath::Abs(mrecpi2pi4 - m_rho) < 0.0743 ) ){
		       	histosTH1F["hm4rec2OSr132405win148"]->Fill(mrec);
		     }
    		     //...| M(pi+pi-) - M(rho) | < 18.6 MeV ==> mass window: sigma = 37.14 MeV
	             if( ( TMath::Abs(mrecpi1pi3 - m_rho) < 0.0186 ) &&
			 ( TMath::Abs(mrecpi2pi4 - m_rho) < 0.0186 ) ){
		       	histosTH1F["hm4rec2OSr132405winhalf"]->Fill(mrec);
		     }
		     // ...delta_phi < 90 degrees: parallel protons TTBB
		     // pT < 0.8 GeV/c
		     if( pipipipiRec.Pt() < 0.8 ){
	             if( ( TMath::Abs(mrecpi1pi3 - m_rho) < 0.0186 ) &&
			 ( TMath::Abs(mrecpi2pi4 - m_rho) < 0.0186 ) ){
		       	histosTH1F["hm4rec2OSr132405winpt"]->Fill(mrec);
		     }}
		     //
		     //...rho
	             if(mrecpi1pi3 < rhohigh && mrecpi1pi3 > rholow &&
			mrecpi2pi4 < rhohigh && mrecpi2pi4 > rholow ){		     
		        histosTH1F["hm4rec2OSr132405"]->Fill(mrec);
		     }		     
	             if(mrecpi1pi3 < rhohigh2 && mrecpi1pi3 > rholow2 &&
			mrecpi2pi4 < rhohigh2 && mrecpi2pi4 > rholow2 ){		     
		        histosTH1F["hm4rec2OSr1324052"]->Fill(mrec);
		     }
    		     //...| M(pi+pi-) - M(rho) | < 35 MeV = 1*sigma, mass window = 2*sigma = 70 MeV
    		     //...| M(pi+pi-) - M(rho) | < 138 MeV = 1*sigma, mass window = 2*sigma = 276 MeV
	             if( ( TMath::Abs(mrecpi1pi3 - m_rho) < sigpipirho ) &&
			 ( TMath::Abs(mrecpi2pi4 - m_rho) < sigpipirho ) ){
		       	histosTH1F["hm4rec2OSr132405sig"]->Fill(mrec);
		     }
    		     //...| M(pi+pi-) - M(rho) | < 60 MeV = 1*sigma, mass window = 2*sigma = 120 MeV
	             if( ( TMath::Abs(mrecpi1pi3 - m_rho) < sigpipirho2 ) &&
			 ( TMath::Abs(mrecpi2pi4 - m_rho) < sigpipirho2 ) ){
		       	histosTH1F["hm4rec2OSr132405sig2"]->Fill(mrec);
		     }
    		     //...| M(pi+pi-) - M(rho) | < 10 MeV = 1*sigma, mass window = 2*sigma = 20 MeV
	             if( ( TMath::Abs(mrecpi1pi3 - m_rho) < sigpipirho3 ) &&
			 ( TMath::Abs(mrecpi2pi4 - m_rho) < sigpipirho3 ) ){
		       	histosTH1F["hm4rec2OSr132405sig3"]->Fill(mrec);
		     }
    		     //...| M(pi+pi-) - M(rho) | < 20 MeV = 1*sigma, mass window = 2*sigma = 40 MeV
	             if( ( TMath::Abs(mrecpi1pi3 - m_rho) < sigpipirho4 ) &&
			 ( TMath::Abs(mrecpi2pi4 - m_rho) < sigpipirho4 ) ){
		       	histosTH1F["hm4rec2OSr132405sig4"]->Fill(mrec);
		     }
    		     //...| M(pi+pi-) - M(rho) | < 30 MeV = 1*sigma, mass window = 2*sigma = 60 MeV
	             if( ( TMath::Abs(mrecpi1pi3 - m_rho) < sigpipirho5 ) &&
			 ( TMath::Abs(mrecpi2pi4 - m_rho) < sigpipirho5 ) ){
		       	histosTH1F["hm4rec2OSr132405sig5"]->Fill(mrec);
		     }
    		     //...| M(pi+pi-) - M(rho) | < 40 MeV = 1*sigma, mass window = 2*sigma = 80 MeV
	             if( ( TMath::Abs(mrecpi1pi3 - m_rho) < sigpipirho6 ) &&
			 ( TMath::Abs(mrecpi2pi4 - m_rho) < sigpipirho6 ) ){
		       	histosTH1F["hm4rec2OSr132405sig6"]->Fill(mrec);
		     }
		  }//...end of Qpair=0...1324
	   //
	   //1423
	   if(charray[0]+charray[3] == 0)
	          {	    
		     histosTH1F["hm4rec2OSr1423t05"]->Fill(mrec);
               	     histosTH2F["h2dim2OSr14x23t05"]->Fill(mrecpi1pi4,mrecpi2pi3);
		     //		     	          
	   // ----------- rho -----------study	          
	   // 1423 for M(rho) shape study ...beased on inspection of M(14)xM(23).
	   // | M(pi2pi3) - M(rho) | < 150 MeV win=300MeV
	             if( TMath::Abs(mrecpi2pi3 - m_rho) < 0.150 ){
		     histosTH1F["hm4rec2OS14rho"]->Fill(mrecpi1pi4);
 		     }

		     // rho
		     // Gamma = 87.46 MeV
		     // sigma = 37.14 MeV
	             //
		     //...final study rho
    		     //...| M(pi+pi-) - M(rho) | < 37.14 MeV ==> mass window: 2*sigma = 74.3 MeV
	             if( ( TMath::Abs(mrecpi1pi4 - m_rho) < 0.03714 ) &&
			 ( TMath::Abs(mrecpi2pi3 - m_rho) < 0.03714 ) ){
		       	histosTH1F["hm4rec2OSr142305win74"]->Fill(mrec);
		     }
    		     //...| M(pi+pi-) - M(rho) | < 74.3 MeV ==> mass window: 4*sigma = 148.6 MeV
	             if( ( TMath::Abs(mrecpi1pi4 - m_rho) < 0.0743 ) &&
			 ( TMath::Abs(mrecpi2pi3 - m_rho) < 0.0743 ) ){
		       	histosTH1F["hm4rec2OSr142305win148"]->Fill(mrec);
		     }
    		     //...| M(pi+pi-) - M(rho) | < 18.6 MeV ==> mass window: sigma = 37.14 MeV
	             if( ( TMath::Abs(mrecpi1pi4 - m_rho) < 0.0186 ) &&
			 ( TMath::Abs(mrecpi2pi3 - m_rho) < 0.0186 ) ){
		       	histosTH1F["hm4rec2OSr142305winhalf"]->Fill(mrec);
		     }
		     // ...delta_phi < 90 degrees: parallel protons TTBB
		     // pT < 0.8 GeV/c
		     if( pipipipiRec.Pt() < 0.8 ){
	             if( ( TMath::Abs(mrecpi1pi4 - m_rho) < 0.0186 ) &&
			 ( TMath::Abs(mrecpi2pi3 - m_rho) < 0.0186 ) ){
		       	histosTH1F["hm4rec2OSr142305winpt"]->Fill(mrec);
		     }}
		     //
		     //...rho
		     // win=40MeV
		     if(mrecpi1pi4 < rhohigh && mrecpi1pi4 > rholow &&
			mrecpi2pi3 < rhohigh && mrecpi2pi3 > rholow ){		     
		        histosTH1F["hm4rec2OSr142305"]->Fill(mrec);
		     }
		     // win=70MeV asymetric
		     if(mrecpi1pi4 < rhohigh2 && mrecpi1pi4 > rholow2 &&
			mrecpi2pi3 < rhohigh2 && mrecpi2pi3 > rholow2 ){		     
		        histosTH1F["hm4rec2OSr1423052"]->Fill(mrec);
		     }
    		     //...| M(pi+pi-) - M(rho) | < 118 MeV = 1*sigma, mass window = 2*sigma =  MeV
	             if( ( TMath::Abs(mrecpi1pi4 - m_rho) < sigpipirho ) &&
			 ( TMath::Abs(mrecpi2pi3 - m_rho) < sigpipirho ) ){
		       	histosTH1F["hm4rec2OSr142305sig"]->Fill(mrec);
		     }
		     //...| M(pi+pi-) - M(rho) | < 60 MeV = 1*sigma, mass window = 2*sigma = 120 MeV
	             if( ( TMath::Abs(mrecpi1pi4 - m_rho) < sigpipirho2 ) &&
			 ( TMath::Abs(mrecpi2pi3 - m_rho) < sigpipirho2 ) ){
		       	histosTH1F["hm4rec2OSr142305sig2"]->Fill(mrec);
		     }
		     //...| M(pi+pi-) - M(rho) | < 10 MeV = 1*sigma, mass window = 2*sigma = 20 MeV
	             if( ( TMath::Abs(mrecpi1pi4 - m_rho) < sigpipirho3 ) &&
			 ( TMath::Abs(mrecpi2pi3 - m_rho) < sigpipirho3 ) ){
		       	histosTH1F["hm4rec2OSr142305sig3"]->Fill(mrec);
		     }
		     //...| M(pi+pi-) - M(rho) | < 20 MeV = 1*sigma, mass window = 2*sigma = 40 MeV
	             if( ( TMath::Abs(mrecpi1pi4 - m_rho) < sigpipirho4 ) &&
			 ( TMath::Abs(mrecpi2pi3 - m_rho) < sigpipirho4 ) ){
		       	histosTH1F["hm4rec2OSr142305sig4"]->Fill(mrec);
		     }
		     //...| M(pi+pi-) - M(rho) | < 30 MeV = 1*sigma, mass window = 2*sigma = 60 MeV
	             if( ( TMath::Abs(mrecpi1pi4 - m_rho) < sigpipirho5 ) &&
			 ( TMath::Abs(mrecpi2pi3 - m_rho) < sigpipirho5 ) ){
		       	histosTH1F["hm4rec2OSr142305sig5"]->Fill(mrec);
		     }
		     //...| M(pi+pi-) - M(rho) | < 40 MeV = 1*sigma, mass window = 2*sigma = 80 MeV
	             if( ( TMath::Abs(mrecpi1pi4 - m_rho) < sigpipirho6 ) &&
			 ( TMath::Abs(mrecpi2pi3 - m_rho) < sigpipirho6 ) ){
		       	histosTH1F["hm4rec2OSr142305sig6"]->Fill(mrec);
		     }
		  }//...end Q_pair=0 ...1423

	    }//...end of PID pions
	     
	   }//...end of nvtx=1
	  
	 } //...end of cut05
	 


	//...cut 8..................theVees	   
	  
/* not now 1

	  //AA...using PID
      	  if(totcharge==0){
	       
	  //...using PID Pions & Kaons for selection=11

	     //...Luiz
	     //histosTH1F["hm4rec2OSvee"]->Fill(mrec);

	     if(isKshort){

	     //...one primary & one Vee && no Lambda //  K+pi- pi+pi-  or  K-pi+ pi+pi- 
	     if(nvtx==1 && nks==1 && nlam==0){

   not now 1 */

	  
      	   /* 
     ...first combining, then select the Q_pairs=0

     pi1pi2 pi3k4
     pi1pi3 pi2k4
     pi2pi3 pi1k4

     pi1pi2 k3pi4
     pi1pi4 k3pi2
     pi2pi4 k3pi1

     pi1k2 pi3pi4
     pi3k2 pi1pi4
     pi4k2 pi1pi3

     k1pi2 pi3pi4
     k1pi3 pi2pi4
     k1pi4 pi2pi3
 	   */
	       //double mrecKpi = 0.0 ;
	       //

	       //...d0 is the transverse impact parameter (dxy) w.r.t. IP
	       //...dz is the longitudinal impact parameter     w.r.t. IP
	       //...vtxdxy is the transverse impact parameter   w.r.t. primary vertex
	       //...vtxdz is the longitudinal impact parameter  w.r.t. primary vertex

/* not now 2  
               //
	       if(charray[0]+charray[1] == 0 && pidarray[2]==3 && pidarray[3]==2
		  && isTrack3 && isTrack4 )
		 {mrecKpi = mrecpi3k4 ; histosTH1F["hm4rec2OSvee11"]->Fill(mrecKpi);}
	       if(charray[0]+charray[2] == 0 && pidarray[1]==3 && pidarray[3]==2
		  && isTrack2 && isTrack4 )
		 {mrecKpi = mrecpi2k4 ; histosTH1F["hm4rec2OSvee11"]->Fill(mrecKpi);}
	       if(charray[1]+charray[2] == 0 && pidarray[0]==3 && pidarray[3]==2
		  && isTrack1 && isTrack4 )
		 {mrecKpi = mrecpi1k4 ; histosTH1F["hm4rec2OSvee11"]->Fill(mrecKpi);}
	       //
	       if(charray[0]+charray[1] == 0 && pidarray[2]==2 && pidarray[3]==3
		  && isTrack3 && isTrack4 )
		 {mrecKpi = mreck3pi4 ; histosTH1F["hm4rec2OSvee11"]->Fill(mrecKpi);}
	       if(charray[0]+charray[3] == 0 && pidarray[2]==2 && pidarray[1]==3
		  && isTrack3 && isTrack2 )
		 {mrecKpi = mreck3pi2 ; histosTH1F["hm4rec2OSvee11"]->Fill(mrecKpi);}
	       if(charray[1]+charray[3] == 0 && pidarray[2]==2 && pidarray[0]==3
		  && isTrack3 && isTrack1 )
		 {mrecKpi = mreck3pi1 ; histosTH1F["hm4rec2OSvee11"]->Fill(mrecKpi);}
	       //
	       if(charray[2]+charray[3] == 0 && pidarray[0]==3 && pidarray[1]==2
		  && isTrack1 && isTrack2 )
		 {mrecKpi = mrecpi1k2 ; histosTH1F["hm4rec2OSvee11"]->Fill(mrecKpi);}
	       if(charray[0]+charray[3] == 0 && pidarray[2]==3 && pidarray[1]==2
		  && isTrack3 && isTrack2 )
		 {mrecKpi = mrecpi3k2 ; histosTH1F["hm4rec2OSvee11"]->Fill(mrecKpi);}
	       if(charray[0]+charray[2] == 0 && pidarray[3]==3 && pidarray[1]==2
		  && isTrack4 && isTrack2 )
		 {mrecKpi = mrecpi4k2 ; histosTH1F["hm4rec2OSvee11"]->Fill(mrecKpi);}
	       //
	       if(charray[2]+charray[3] == 0 && pidarray[0]==2 && pidarray[1]==3
		  && isTrack1 && isTrack2 )
		 {mrecKpi = mreck1pi2 ; histosTH1F["hm4rec2OSvee11"]->Fill(mrecKpi);}
	       if(charray[1]+charray[3] == 0 && pidarray[0]==2 && pidarray[2]==3
		  && isTrack1 && isTrack3 )
		 {mrecKpi = mreck1pi3 ; histosTH1F["hm4rec2OSvee11"]->Fill(mrecKpi);}
	       if(charray[1]+charray[2] == 0 && pidarray[0]==2 && pidarray[3]==3
		  && isTrack1 && isTrack4 )
		 {mrecKpi = mreck1pi4 ; histosTH1F["hm4rec2OSvee11"]->Fill(mrecKpi);}


	       //A  
	       if(charray[0]+charray[1] == 0 && pidarray[0]==3 && pidarray[1]==3 && isTrack1 && isTrack2 ) histosTH1F["hm4rec2OS_pi1pi2vee11"]->Fill(mrecpi1pi2);
	       if(charray[2]+charray[3] == 0 && pidarray[2]==3 && pidarray[3]==2 && isTrack3 && isTrack4 ) histosTH1F["hm4rec2OS_pi3k4vee11"]->Fill(mrecpi3k4);    
	       ///histosTH2F["hm4dim2OS_pi1pi2_pi3k4vee11"]->Fill(mrecpi1pi2,mrecpi3k4);
	       //     
	       if(charray[0]+charray[2] == 0 && pidarray[0]==3 && pidarray[2]==3 && isTrack1 && isTrack3 ) histosTH1F["hm4rec2OS_pi1pi3vee11"]->Fill(mrecpi1pi3);
	       if(charray[1]+charray[3] == 0 && pidarray[1]==3 && pidarray[3]==2 && isTrack2 && isTrack4 ) histosTH1F["hm4rec2OS_pi2k4vee11"]->Fill(mrecpi2k4);
	       ///histosTH2F["hm4dim2OS_pi1pi3_pi2k4vee11"]->Fill(mrecpi1pi3,mrecpi2k4);
	       //
	       if(charray[1]+charray[2] == 0 && pidarray[1]==3 && pidarray[2]==3 && isTrack2 && isTrack3 ) histosTH1F["hm4rec2OS_pi2pi3vee11"]->Fill(mrecpi2pi3);
	       if(charray[0]+charray[3] == 0 && pidarray[0]==3 && pidarray[3]==2 && isTrack1 && isTrack4 ) histosTH1F["hm4rec2OS_pi1k4vee11"]->Fill(mrecpi1k4);
	       ///histosTH2F["hm4dim2OS_pi2pi3_pi1k4vee11"]->Fill(mrecpi2pi3,mrecpi1k4);

	       //B  
	       //if(charray[0]+charray[1] == 0) histosTH1F["hm4rec2OS_pi1pi2vee11"]->Fill(mrecpi1pi2);
	       if(charray[2]+charray[3] == 0 && pidarray[2]==2 && pidarray[3]==3 && isTrack3 && isTrack4 ) histosTH1F["hm4rec2OS_k3pi4vee11"]->Fill(mreck3pi4);
	       ///histosTH2F["hm4dim2OS_pi1pi2_k3pi4vee11"]->Fill(mrecpi1pi2,mreck3pi4);
	       //
	       if(charray[0]+charray[3] == 0 && pidarray[0]==3 && pidarray[3]==3 && isTrack1 && isTrack4 ) histosTH1F["hm4rec2OS_pi1pi4vee11"]->Fill(mrecpi1pi4);
	       if(charray[2]+charray[1] == 0 && pidarray[2]==2 && pidarray[1]==3 && isTrack3 && isTrack2 ) histosTH1F["hm4rec2OS_k3pi2vee11"]->Fill(mreck3pi2);
	       ///histosTH2F["hm4dim2OS_pi1pi4_k3pi2vee11"]->Fill(mrecpi1pi4,mreck3pi2);
	       //
	       if(charray[1]+charray[3] == 0 && pidarray[1]==3 && pidarray[3]==3 && isTrack2 && isTrack4 ) histosTH1F["hm4rec2OS_pi2pi4vee11"]->Fill(mrecpi2pi4);
	       if(charray[2]+charray[0] == 0 && pidarray[2]==2 && pidarray[0]==3 && isTrack3 && isTrack1 ) histosTH1F["hm4rec2OS_k3pi1vee11"]->Fill(mreck3pi1);
	       ///histosTH2F["hm4dim2OS_pi2pi4_k3pi1vee11"]->Fill(mrecpi2pi4,mreck3pi1);

	       //C
	       if(charray[0]+charray[1] == 0 && pidarray[0]==3 && pidarray[1]==2 && isTrack1 && isTrack2 ) histosTH1F["hm4rec2OS_pi1k2vee11"]->Fill(mrecpi1k2);
	       if(charray[2]+charray[3] == 0 && pidarray[2]==3 && pidarray[3]==3 && isTrack3 && isTrack4 ) histosTH1F["hm4rec2OS_pi3pi4vee11"]->Fill(mrecpi3pi4);
      	       ///histosTH2F["hm4dim2OS_pi1k2_pi3pi4vee11"]->Fill(mrecpi1k2,mrecpi3pi4);
	       //
	       if(charray[2]+charray[1] == 0 && pidarray[2]==3 && pidarray[1]==2 && isTrack3 && isTrack2 ) histosTH1F["hm4rec2OS_pi3k2vee11"]->Fill(mrecpi3k2);
	       //if(charray[0]+charray[3] == 0) histosTH1F["hm4rec2OS_pi1pi4vee11"]->Fill(mrecpi1pi4);
	       ///histosTH2F["hm4dim2OS_pi3k2_pi1pi4vee11"]->Fill(mrecpi3k2,mrecpi1pi4);
	       //
	       if(charray[3]+charray[1] == 0 && pidarray[3]==3 && pidarray[1]==2 && isTrack4 && isTrack2 ) histosTH1F["hm4rec2OS_pi4k2vee11"]->Fill(mrecpi4k2);
	       //if(charray[0]+charray[2] == 0) histosTH1F["hm4rec2OS_pi1pi3vee11"]->Fill(mrecpi1pi3);
	       ///histosTH2F["hm4dim2OS_pi4k2_pi1pi3vee11"]->Fill(mrecpi4k2,mrecpi1pi3);

	       //D
	       if(charray[0]+charray[1] == 0 && pidarray[0]==2 && pidarray[1]==3 && isTrack1 && isTrack2 ) histosTH1F["hm4rec2OS_k1pi2vee11"]->Fill(mreck1pi2);
	       //if(charray[2]+charray[3] == 0) histosTH1F["hm4rec2OS_pi3pi4vee11"]->Fill(mrecpi3pi4);
	       ///histosTH2F["hm4dim2OS_k1pi2_pi3pi4vee11"]->Fill(mreck1pi2,mrecpi3pi4);
	       //
	       if(charray[0]+charray[2] == 0 && pidarray[0]==2 && pidarray[2]==3 && isTrack1 && isTrack3 ) histosTH1F["hm4rec2OS_k1pi3vee11"]->Fill(mreck1pi3);
	       //if(charray[1]+charray[3] == 0) histosTH1F["hm4rec2OS_pi2pi4vee11"]->Fill(mrecpi2pi4);
	       ///histosTH2F["hm4dim2OS_k1pi3_pi2pi4vee11"]->Fill(mreck1pi3,mrecpi2pi4);
	       //
	       if(charray[0]+charray[3] == 0 && pidarray[0]==2 && pidarray[3]==3 && isTrack1 && isTrack4 ) histosTH1F["hm4rec2OS_k1pi4vee11"]->Fill(mreck1pi4);
	       //if(charray[1]+charray[2] == 0) histosTH1F["hm4rec2OS_pi2pi3vee11"]->Fill(mrecpi2pi3);
	       ///histosTH2F["hm4dim2OS_k1pi3_pi2pi4vee11"]->Fill(mreck1pi3,mrecpi2pi4);

	     } //end of nvtx=1 nks=1 nlam=0
	     //} //...end of PID Pions & Kaons

	  //...using PID Pions for selection=02 or 01
	     if(pidarray[0]==3 && pidarray[1]==3 && pidarray[2]==3 && pidarray[3]==3)
	     {
	     //...no primary & two Vees
	     if(nvtx==0 && nks==2){
	     histosTH1F["hm4rec2OSvee02"]->Fill(mrec);
		 if(charray[0]+charray[1] == 0)
	         {
	       histosTH1F["hm4rec2OS_pi1pi2vee02"]->Fill(mrecpi1pi2);
	       histosTH1F["hm4rec2OS_pi3pi4vee02"]->Fill(mrecpi3pi4);
	       histosTH2F["hm4dim2OS_pi1pi2_pi3pi4vee02"]->Fill(mrecpi1pi2,mrecpi3pi4);
	         }else if(charray[0]+charray[2] == 0){
	       histosTH1F["hm4rec2OS_pi1pi3vee02"]->Fill(mrecpi1pi3);
	       histosTH1F["hm4rec2OS_pi2pi4vee02"]->Fill(mrecpi2pi4);
	       histosTH2F["hm4dim2OS_pi1pi3_pi2pi4vee02"]->Fill(mrecpi1pi3,mrecpi2pi4);
		 }
	       } //end of nvtx=0 nks=2

	     //...no primary & 1 Vee
             if(nvtx==0 && nks==1){
	     histosTH1F["hm4rec2OSvee01"]->Fill(mrec); 
		 if(charray[0]+charray[1] == 0)
	         {
	       histosTH1F["hm4rec2OS_pi1pi2vee01"]->Fill(mrecpi1pi2);
	       histosTH1F["hm4rec2OS_pi3pi4vee01"]->Fill(mrecpi3pi4);
	       histosTH2F["hm4dim2OS_pi1pi2_pi3pi4vee01"]->Fill(mrecpi1pi2,mrecpi3pi4);
	         }else if(charray[0]+charray[2] == 0){
	       histosTH1F["hm4rec2OS_pi1pi3vee01"]->Fill(mrecpi1pi3);
	       histosTH1F["hm4rec2OS_pi2pi4vee01"]->Fill(mrecpi2pi4);
	       histosTH2F["hm4dim2OS_pi1pi3_pi2pi4vee01"]->Fill(mrecpi1pi3,mrecpi2pi4);
		 }
	       } //end of nvtx=0 nks=1

	       } //...end of PID Pions	     
	     } //...end of isKshort
	   } //...end of totalcharge=0	    
	   //AA...end of PID

not now 2 */

	   //BB...PID Pions or Kaons
     	   if(totcharge==0){

	     //...ntrk vs nks
	     histosTH2F["hntrknksq0"]->Fill(ntrk,nks);
	     //...nvtx vs nks
	     histosTH2F["hnvtxnksq0"]->Fill(nvtx,nks);

	     if(isKshort){

	       //...veeno11
	       //...type:11 is double counting ???
	     if(nvtx==1 && nks==1){

               //double rlimit = 0.5;
	       // k4
	       if(pidarray[0]==2 && pidarray[1]==2 && pidarray[2]==2 && pidarray[3]==3){
	       if(charray[0]+charray[1] == 0 
		  && isTrack3 && isTrack4 )
		 {mrecKpi = mrecpi3k4 ; histosTH1F["hm4rec2OSveeno11Kpi"]->Fill(mrecKpi);
		                        histosTH1F["hm4rec2OSveeno11x"]->Fill(mrecpi1pi2_pi3k4);
		                        histosTH1F["hm4rec2OSveeno11"]->Fill(mrecpi1pi2pi3k4);}
	       if(charray[0]+charray[2] == 0 
		  && isTrack2 && isTrack4  )
		 {mrecKpi = mrecpi2k4 ; histosTH1F["hm4rec2OSveeno11Kpi"]->Fill(mrecKpi);
	                                histosTH1F["hm4rec2OSveeno11x"]->Fill(mrecpi1pi3_pi2k4);
	                                histosTH1F["hm4rec2OSveeno11"]->Fill(mrecpi1pi2pi3k4);}
	       if(charray[1]+charray[2] == 0 
		  && isTrack1 && isTrack4 )
		 {mrecKpi = mrecpi1k4 ; histosTH1F["hm4rec2OSveeno11Kpi"]->Fill(mrecKpi);
		                        histosTH1F["hm4rec2OSveeno11x"]->Fill(mrecpi2pi3_pi1k4);
		                        histosTH1F["hm4rec2OSveeno11"]->Fill(mrecpi1pi2pi3k4);}
	       }
	       // k3
       	       if(pidarray[0]==2 && pidarray[1]==2 && pidarray[2]==3 && pidarray[3]==2){
	       if(charray[0]+charray[1] == 0 
		  && isTrack3 && isTrack4 )
		 {mrecKpi = mreck3pi4 ; histosTH1F["hm4rec2OSveeno11Kpi"]->Fill(mrecKpi);
		                        histosTH1F["hm4rec2OSveeno11x"]->Fill(mrecpi1pi2_k3pi4);
		                        histosTH1F["hm4rec2OSveeno11"]->Fill(mrecpi1pi2k3pi4);}
	       if(charray[0]+charray[3] == 0 
		  && isTrack3 && isTrack2 )
		 {mrecKpi = mreck3pi2 ; histosTH1F["hm4rec2OSveeno11Kpi"]->Fill(mrecKpi);
		                        histosTH1F["hm4rec2OSveeno11x"]->Fill(mrecpi1pi4_k3pi2);
		                        histosTH1F["hm4rec2OSveeno11"]->Fill(mrecpi1pi2k3pi4);}
	       if(charray[1]+charray[3] == 0 
		  && isTrack3 && isTrack1 )
		 {mrecKpi = mreck3pi1 ; histosTH1F["hm4rec2OSveeno11Kpi"]->Fill(mrecKpi);
		                        histosTH1F["hm4rec2OSveeno11x"]->Fill(mrecpi2pi4_k3pi1);
		                        histosTH1F["hm4rec2OSveeno11"]->Fill(mrecpi1pi2k3pi4);}
	       }
	       // k2
       	       if(pidarray[0]==2 && pidarray[1]==3 && pidarray[2]==2 && pidarray[3]==2){
	       if(charray[2]+charray[3] == 0 
		  && isTrack1 && isTrack2 )
		 {mrecKpi = mrecpi1k2 ; histosTH1F["hm4rec2OSveeno11Kpi"]->Fill(mrecKpi);
		                        histosTH1F["hm4rec2OSveeno11x"]->Fill(mrecpi1k2_pi3pi4);
		                        histosTH1F["hm4rec2OSveeno11"]->Fill(mrecpi1k2pi3pi4);}
	       if(charray[0]+charray[3] == 0 
		  && isTrack3 && isTrack2 )
		 {mrecKpi = mrecpi3k2 ; histosTH1F["hm4rec2OSveeno11Kpi"]->Fill(mrecKpi);
		                        histosTH1F["hm4rec2OSveeno11x"]->Fill(mrecpi3k2_pi1pi4);
		                        histosTH1F["hm4rec2OSveeno11"]->Fill(mrecpi1k2pi3pi4);}
	       if(charray[0]+charray[2] == 0 
		  && isTrack4 && isTrack2 )
		 {mrecKpi = mrecpi4k2 ; histosTH1F["hm4rec2OSveeno11Kpi"]->Fill(mrecKpi);
		                        histosTH1F["hm4rec2OSveeno11x"]->Fill(mrecpi4k2_pi1pi3);
		                        histosTH1F["hm4rec2OSveeno11"]->Fill(mrecpi1k2pi3pi4);}
	       }
	       // k1
      	       if(pidarray[0]==3 && pidarray[1]==2 && pidarray[2]==2 && pidarray[3]==2){
	       if(charray[2]+charray[3] == 0 
		  && isTrack1 && isTrack2 )
		 {mrecKpi = mreck1pi2 ; histosTH1F["hm4rec2OSveeno11Kpi"]->Fill(mrecKpi);
		                        histosTH1F["hm4rec2OSveeno11x"]->Fill(mreck1pi2_pi3pi4);
		                        histosTH1F["hm4rec2OSveeno11"]->Fill(mreck1pi2pi3pi4);}
	       if(charray[1]+charray[3] == 0
		  && isTrack1 && isTrack3 )
		 {mrecKpi = mreck1pi3 ; histosTH1F["hm4rec2OSveeno11Kpi"]->Fill(mrecKpi);
		                        histosTH1F["hm4rec2OSveeno11x"]->Fill(mreck1pi3_pi2pi4);
		                        histosTH1F["hm4rec2OSveeno11"]->Fill(mreck1pi2pi3pi4);}
	       if(charray[1]+charray[2] == 0 
		  && isTrack1 && isTrack4 )
		 {mrecKpi = mreck1pi4 ; histosTH1F["hm4rec2OSveeno11Kpi"]->Fill(mrecKpi);
		                        histosTH1F["hm4rec2OSveeno11x"]->Fill(mreck1pi4_pi2pi3);
		                        histosTH1F["hm4rec2OSveeno11"]->Fill(mreck1pi2pi3pi4);}
	       }

	       //...attention here!
	       
	       //A  k4
	       if(pidarray[0]==2 && pidarray[1]==2 && pidarray[2]==2 && pidarray[3]==3){
	       if(charray[0]+charray[1] == 0 && isTrack3 && isTrack4 ) {
		 histosTH1F["hm4rec2OS_pi1pi2k4veeno11"]->Fill(mrecpi1pi2); //a
	         histosTH1F["hm4rec2OS_pi3k4veeno11"]->Fill(mrecpi3k4);}    
	       ///histosTH2F["hm4dim2OS_pi1pi2_pi3k4veeno11"]->Fill(mrecpi1pi2,mrecpi3k4);
	       //     
	       if(charray[0]+charray[2] == 0 && isTrack2 && isTrack4 ) {
		 histosTH1F["hm4rec2OS_pi1pi3k4veeno11"]->Fill(mrecpi1pi3); //b
	         histosTH1F["hm4rec2OS_pi2k4veeno11"]->Fill(mrecpi2k4);}
	       ///histosTH2F["hm4dim2OS_pi1pi3_pi2k4veeno11"]->Fill(mrecpi1pi3,mrecpi2k4);
	       //
	       if(charray[1]+charray[2] == 0 && isTrack1 && isTrack4 ) {
		 histosTH1F["hm4rec2OS_pi2pi3k4veeno11"]->Fill(mrecpi2pi3); //c
	         histosTH1F["hm4rec2OS_pi1k4veeno11"]->Fill(mrecpi1k4);}
	       ///histosTH2F["hm4dim2OS_pi2pi3_pi1k4veeno11"]->Fill(mrecpi2pi3,mrecpi1k4);
	       }
	       //B  k3
       	       if(pidarray[0]==2 && pidarray[1]==2 && pidarray[2]==3 && pidarray[3]==2){
	       if(charray[0]+charray[1] == 0 && isTrack3 && isTrack4 ) {
		 histosTH1F["hm4rec2OS_pi1pi2k3veeno11"]->Fill(mrecpi1pi2); //ax2
	         histosTH1F["hm4rec2OS_k3pi4veeno11"]->Fill(mreck3pi4);}
	       ///histosTH2F["hm4dim2OS_pi1pi2_k3pi4veeno11"]->Fill(mrecpi1pi2,mreck3pi4);
	       //
	       if(charray[0]+charray[3] == 0 && isTrack3 && isTrack2 ) {
		 histosTH1F["hm4rec2OS_pi1pi4k3veeno11"]->Fill(mrecpi1pi4); //d
	         histosTH1F["hm4rec2OS_k3pi2veeno11"]->Fill(mreck3pi2);}
	       ///histosTH2F["hm4dim2OS_pi1pi4_k3pi2veeno11"]->Fill(mrecpi1pi4,mreck3pi2);
	       //
	       if(charray[1]+charray[3] == 0 && isTrack3 && isTrack1 ) {
		 histosTH1F["hm4rec2OS_pi2pi4k3veeno11"]->Fill(mrecpi2pi4); //e
	         histosTH1F["hm4rec2OS_k3pi1veeno11"]->Fill(mreck3pi1);}
	       ///histosTH2F["hm4dim2OS_pi2pi4_k3pi1veeno11"]->Fill(mrecpi2pi4,mreck3pi1);
	       }
	       //C  k2
      	       if(pidarray[0]==2 && pidarray[1]==3 && pidarray[2]==2 && pidarray[3]==2){
	       if(charray[2]+charray[3] == 0 && isTrack1 && isTrack2 ) {
		 histosTH1F["hm4rec2OS_pi1k2veeno11"]->Fill(mrecpi1k2);
                 histosTH1F["hm4rec2OS_pi3pi4k2veeno11"]->Fill(mrecpi3pi4);} //f
      	       ///histosTH2F["hm4dim2OS_pi1k2_pi3pi4veeno11"]->Fill(mrecpi1k2,mrecpi3pi4);
	       //
	       if(charray[0]+charray[3] == 0 && isTrack3 && isTrack2 ) {
		 histosTH1F["hm4rec2OS_pi3k2veeno11"]->Fill(mrecpi3k2);
	         histosTH1F["hm4rec2OS_pi1pi4k2veeno11"]->Fill(mrecpi1pi4);} //dx2
	       ///histosTH2F["hm4dim2OS_pi3k2_pi1pi4veeno11"]->Fill(mrecpi3k2,mrecpi1pi4);
	       //
	       if(charray[0]+charray[2] == 0 && isTrack4 && isTrack2 ) {
		 histosTH1F["hm4rec2OS_pi4k2veeno11"]->Fill(mrecpi4k2);
	         histosTH1F["hm4rec2OS_pi1pi3k2veeno11"]->Fill(mrecpi1pi3);} //bx2
	       ///histosTH2F["hm4dim2OS_pi4k2_pi1pi3veeno11"]->Fill(mrecpi4k2,mrecpi1pi3);
	       }
	       //D  k1
      	       if(pidarray[0]==3 && pidarray[1]==2 && pidarray[2]==2 && pidarray[3]==2){
	       if(charray[2]+charray[3] == 0 && isTrack1 && isTrack2 ) {
		 histosTH1F["hm4rec2OS_k1pi2veeno11"]->Fill(mreck1pi2);
	         histosTH1F["hm4rec2OS_pi3pi4k1veeno11"]->Fill(mrecpi3pi4);} //fx2
	       ///histosTH2F["hm4dim2OS_k1pi2_pi3pi4veeno11"]->Fill(mreck1pi2,mrecpi3pi4);
	       //
	       if(charray[1]+charray[3] == 0 && isTrack1 && isTrack3 ) {
		 histosTH1F["hm4rec2OS_k1pi3veeno11"]->Fill(mreck1pi3);
	         histosTH1F["hm4rec2OS_pi2pi4k1veeno11"]->Fill(mrecpi2pi4);} //ex2
	       ///histosTH2F["hm4dim2OS_k1pi3_pi2pi4veeno11"]->Fill(mreck1pi3,mrecpi2pi4);
	       //
	       if(charray[1]+charray[2] == 0 && isTrack1 && isTrack4 ) {
		 histosTH1F["hm4rec2OS_k1pi4veeno11"]->Fill(mreck1pi4);
	         histosTH1F["hm4rec2OS_pi2pi3k1veeno11"]->Fill(mrecpi2pi3);} //cx2
	       ///histosTH2F["hm4dim2OS_k1pi3_pi2pi4veeno11"]->Fill(mreck1pi3,mrecpi2pi4);
	       }
	     } //end of nvtx=1 nks=1


	     //...veeno02
             if(nvtx==0 && nks==2){

	       //...only pions ...BETTER USING PID
             if(pidarray[0]==2 && pidarray[1]==2 && pidarray[2]==2 && pidarray[3]==2){
	       ////if(pidv0array[0]==2 && pidv0array[1]==2 && pidv0array[2]==2 && pidv0array[3]==2){
	       
	     histosTH1F["hm4rec2OSveeno02"]->Fill(mrec);      
             histosTH1F["hm4rec2OSveeno02cor"]->Fill(mrec,1/effictotal);
	     if(ksmass >= 0.466 && ksmass <= 0.530){
             histosTH1F["hm4rec2OSveeno02col"]->Fill(mrec);
             histosTH1F["hm4rec2OSveeno02colcor"]->Fill(mrec,1/effictotal);
	     }
		 if(charray[0]+charray[1] == 0)
	         {
	       histosTH1F["hm4rec2OS_pi1pi2veeno02"]->Fill(mrecpi1pi2);
	       histosTH1F["hm4rec2OS_pi3pi4veeno02"]->Fill(mrecpi3pi4);
	       histosTH2F["hm4dim2OS_pi1pi2_pi3pi4veeno02"]->Fill(mrecpi1pi2,mrecpi3pi4);
	         }else if(charray[0]+charray[2] == 0){
	       histosTH1F["hm4rec2OS_pi1pi3veeno02"]->Fill(mrecpi1pi3);
	       histosTH1F["hm4rec2OS_pi2pi4veeno02"]->Fill(mrecpi2pi4);
	       histosTH2F["hm4dim2OS_pi1pi3_pi2pi4veeno02"]->Fill(mrecpi1pi3,mrecpi2pi4);
		 }
	 // |-(t1+t2)|
         // ...OR
	 if(t1 <= limtA || t2 <= limtA ){histosTH1F["hm4rec2OSveeno02t12a"]->Fill(mrec);}//0.4
	 if(t1 <= limtB || t2 <= limtB ){histosTH1F["hm4rec2OSveeno02t12b"]->Fill(mrec);}//0.3
	 if(t1 <= limtC || t2 <= limtC ){histosTH1F["hm4rec2OSveeno02t12c"]->Fill(mrec);}//0.2
	 if(t1 <= limtD || t2 <= limtD ){histosTH1F["hm4rec2OSveeno02t12d"]->Fill(mrec);}//0.17
	 if(t1 <= limtE || t2 <= limtE ){histosTH1F["hm4rec2OSveeno02t12e"]->Fill(mrec);}//0.1
	 if(t1 <= limtF || t2 <= limtF ){histosTH1F["hm4rec2OSveeno02t12f"]->Fill(mrec);}//0.08
	 if(t1 <= limtG || t2 <= limtG ){histosTH1F["hm4rec2OSveeno02t12g"]->Fill(mrec);}//0.06
	 if(t1 <= limtH || t2 <= limtH ){histosTH1F["hm4rec2OSveeno02t12h"]->Fill(mrec);}//0.05
         // ...OR complementary
         if(t1 > limtA || t2 > limtA ){histosTH1F["hm4rec2OSveeno02t12aC"]->Fill(mrec);}
         if(t1 > limtB || t2 > limtB ){histosTH1F["hm4rec2OSveeno02t12bC"]->Fill(mrec);}
         if(t1 > limtC || t2 > limtC ){histosTH1F["hm4rec2OSveeno02t12cC"]->Fill(mrec);}
         if(t1 > limtD || t2 > limtD ){histosTH1F["hm4rec2OSveeno02t12dC"]->Fill(mrec);}
         if(t1 > limtE || t2 > limtE ){histosTH1F["hm4rec2OSveeno02t12eC"]->Fill(mrec);}
         if(t1 > limtF || t2 > limtF ){histosTH1F["hm4rec2OSveeno02t12fC"]->Fill(mrec);}
         if(t1 > limtG || t2 > limtG ){histosTH1F["hm4rec2OSveeno02t12gC"]->Fill(mrec);}
         if(t1 > limtH || t2 > limtH ){histosTH1F["hm4rec2OSveeno02t12hC"]->Fill(mrec);}
         // ...AND
	 if(t1 <= limtA && t2 <= limtA ){histosTH1F["hm4rec2OSveeno02t12aa"]->Fill(mrec);}//0.4
	 if(t1 <= limtB && t2 <= limtB ){histosTH1F["hm4rec2OSveeno02t12bb"]->Fill(mrec);}//0.3
	 if(t1 <= limtC && t2 <= limtC ){histosTH1F["hm4rec2OSveeno02t12cc"]->Fill(mrec);}//0.2
	 if(t1 <= limtD && t2 <= limtD ){histosTH1F["hm4rec2OSveeno02t12dd"]->Fill(mrec);}//0.17
	 if(t1 <= limtE && t2 <= limtE ){histosTH1F["hm4rec2OSveeno02t12ee"]->Fill(mrec);}//0.1
	 if(t1 <= limtF && t2 <= limtF ){histosTH1F["hm4rec2OSveeno02t12ff"]->Fill(mrec);}//0.08
	 if(t1 <= limtG && t2 <= limtG ){histosTH1F["hm4rec2OSveeno02t12gg"]->Fill(mrec);}//0.06
	 if(t1 <= limtH && t2 <= limtH ){histosTH1F["hm4rec2OSveeno02t12hh"]->Fill(mrec);}//0.05
         // ...AND complementary
         if(t1 > limtA && t2 > limtA ){histosTH1F["hm4rec2OSveeno02t12aaC"]->Fill(mrec);}
         if(t1 > limtB && t2 > limtB ){histosTH1F["hm4rec2OSveeno02t12bbC"]->Fill(mrec);}
         if(t1 > limtC && t2 > limtC ){histosTH1F["hm4rec2OSveeno02t12ccC"]->Fill(mrec);}
         if(t1 > limtD && t2 > limtD ){histosTH1F["hm4rec2OSveeno02t12ddC"]->Fill(mrec);}
         if(t1 > limtE && t2 > limtE ){histosTH1F["hm4rec2OSveeno02t12eeC"]->Fill(mrec);}
         if(t1 > limtF && t2 > limtF ){histosTH1F["hm4rec2OSveeno02t12ffC"]->Fill(mrec);}
         if(t1 > limtG && t2 > limtG ){histosTH1F["hm4rec2OSveeno02t12ggC"]->Fill(mrec);}
         if(t1 > limtH && t2 > limtH ){histosTH1F["hm4rec2OSveeno02t12hhC"]->Fill(mrec);}
	      
	 // dphi
	 if(m1370){ histosTH1F["hdphig1370"]->Fill(TOTEMdphi*180/TMath::Pi());}
	 if(m1549){ histosTH1F["hdphig1549"]->Fill(TOTEMdphi*180/TMath::Pi());}
	 if(m1731){ histosTH1F["hdphig1731"]->Fill(TOTEMdphi*180/TMath::Pi());}

	 // (dphi,mrec) map
	 histosTH2F["h2dimdphigm"]->Fill(TOTEMdphi*180/TMath::Pi(),mrec);
	 // (t12,mrec) map                  ................
	 histosTH2F["h2dimt12m"]->Fill(t12,mrec);
		
	     }//...end of PID=pions
	     } //end of nvtx=0 nks=2

	     
	     //...veeno01
             if(nvtx==0 && nks==1){

	       //...only pions ...BETTER USING PID
	     if(pidarray[0]==2 && pidarray[1]==2 && pidarray[2]==2 && pidarray[3]==2){
	       ////if(pidv0array[0]==2 && pidv0array[1]==2 && pidv0array[2]==2 && pidv0array[3]==2){
	       
	     histosTH1F["hm4rec2OSveeno01"]->Fill(mrec);      
		 if(charray[0]+charray[1] == 0)
	         {
	       histosTH1F["hm4rec2OS_pi1pi2veeno01"]->Fill(mrecpi1pi2);
	       histosTH1F["hm4rec2OS_pi3pi4veeno01"]->Fill(mrecpi3pi4);
	       histosTH2F["hm4dim2OS_pi1pi2_pi3pi4veeno01"]->Fill(mrecpi1pi2,mrecpi3pi4);
	         }else if(charray[0]+charray[2] == 0){
	       histosTH1F["hm4rec2OS_pi1pi3veeno01"]->Fill(mrecpi1pi3);
	       histosTH1F["hm4rec2OS_pi2pi4veeno01"]->Fill(mrecpi2pi4);
	       histosTH2F["hm4dim2OS_pi1pi3_pi2pi4veeno01"]->Fill(mrecpi1pi3,mrecpi2pi4);
		 }
	 // |-(t1+t2)|
         // ...OR
	 if(t1 <= limtA || t2 <= limtA ){histosTH1F["hm4rec2OSveeno01t12a"]->Fill(mrec);}//0.4
	 if(t1 <= limtB || t2 <= limtB ){histosTH1F["hm4rec2OSveeno01t12b"]->Fill(mrec);}//0.3
	 if(t1 <= limtC || t2 <= limtC ){histosTH1F["hm4rec2OSveeno01t12c"]->Fill(mrec);}//0.2
	 if(t1 <= limtD || t2 <= limtD ){histosTH1F["hm4rec2OSveeno01t12d"]->Fill(mrec);}//0.17
	 if(t1 <= limtE || t2 <= limtE ){histosTH1F["hm4rec2OSveeno01t12e"]->Fill(mrec);}//0.1
	 if(t1 <= limtF || t2 <= limtF ){histosTH1F["hm4rec2OSveeno01t12f"]->Fill(mrec);}//0.08
	 if(t1 <= limtG || t2 <= limtG ){histosTH1F["hm4rec2OSveeno01t12g"]->Fill(mrec);}//0.06
	 if(t1 <= limtH || t2 <= limtH ){histosTH1F["hm4rec2OSveeno01t12h"]->Fill(mrec);}//0.05
         // ...OR complementary
         if(t1 > limtA || t2 > limtA ){histosTH1F["hm4rec2OSveeno01t12aC"]->Fill(mrec);}
         if(t1 > limtB || t2 > limtB ){histosTH1F["hm4rec2OSveeno01t12bC"]->Fill(mrec);}
         if(t1 > limtC || t2 > limtC ){histosTH1F["hm4rec2OSveeno01t12cC"]->Fill(mrec);}
         if(t1 > limtD || t2 > limtD ){histosTH1F["hm4rec2OSveeno01t12dC"]->Fill(mrec);}
         if(t1 > limtE || t2 > limtE ){histosTH1F["hm4rec2OSveeno01t12eC"]->Fill(mrec);}
         if(t1 > limtF || t2 > limtF ){histosTH1F["hm4rec2OSveeno01t12fC"]->Fill(mrec);}
         if(t1 > limtG || t2 > limtG ){histosTH1F["hm4rec2OSveeno01t12gC"]->Fill(mrec);}
         if(t1 > limtH || t2 > limtH ){histosTH1F["hm4rec2OSveeno01t12hC"]->Fill(mrec);}
         // ...AND
	 if(t1 <= limtA && t2 <= limtA ){histosTH1F["hm4rec2OSveeno01t12aa"]->Fill(mrec);}//0.4
	 if(t1 <= limtB && t2 <= limtB ){histosTH1F["hm4rec2OSveeno01t12bb"]->Fill(mrec);}//0.3
	 if(t1 <= limtC && t2 <= limtC ){histosTH1F["hm4rec2OSveeno01t12cc"]->Fill(mrec);}//0.2
	 if(t1 <= limtD && t2 <= limtD ){histosTH1F["hm4rec2OSveeno01t12dd"]->Fill(mrec);}//0.17
	 if(t1 <= limtE && t2 <= limtE ){histosTH1F["hm4rec2OSveeno01t12ee"]->Fill(mrec);}//0.1
	 if(t1 <= limtF && t2 <= limtF ){histosTH1F["hm4rec2OSveeno01t12ff"]->Fill(mrec);}//0.08
	 if(t1 <= limtG && t2 <= limtG ){histosTH1F["hm4rec2OSveeno01t12gg"]->Fill(mrec);}//0.06
	 if(t1 <= limtH && t2 <= limtH ){histosTH1F["hm4rec2OSveeno01t12hh"]->Fill(mrec);}//0.05
         // ...AND complementary
         if(t1 > limtA && t2 > limtA ){histosTH1F["hm4rec2OSveeno01t12aaC"]->Fill(mrec);}
         if(t1 > limtB && t2 > limtB ){histosTH1F["hm4rec2OSveeno01t12bbC"]->Fill(mrec);}
         if(t1 > limtC && t2 > limtC ){histosTH1F["hm4rec2OSveeno01t12ccC"]->Fill(mrec);}
         if(t1 > limtD && t2 > limtD ){histosTH1F["hm4rec2OSveeno01t12ddC"]->Fill(mrec);}
         if(t1 > limtE && t2 > limtE ){histosTH1F["hm4rec2OSveeno01t12eeC"]->Fill(mrec);}
         if(t1 > limtF && t2 > limtF ){histosTH1F["hm4rec2OSveeno01t12ffC"]->Fill(mrec);}
         if(t1 > limtG && t2 > limtG ){histosTH1F["hm4rec2OSveeno01t12ggC"]->Fill(mrec);}
         if(t1 > limtH && t2 > limtH ){histosTH1F["hm4rec2OSveeno01t12hhC"]->Fill(mrec);}

	       }//...end of PID=pions
	      } //end of nvtx=0 nks=1
	     
	     } //...end of isKshort
	   } //BB...end of totalcharge=0

	   

	   //CC...no PID
     	   if(totcharge==0){

	     if(isKshort){

	     //...veeno02
             if(nvtx==0 && nks==2){

	       //...only pions ...BETTER USING PID
	       //             if(pidarray[0]==2 && pidarray[1]==2 && pidarray[2]==2 && pidarray[3]==2){
	       ////if(pidv0array[0]==2 && pidv0array[1]==2 && pidv0array[2]==2 && pidv0array[3]==2){
	       
	     histosTH1F["hm4rec2OSveenopid02"]->Fill(mrec);      
             histosTH1F["hm4rec2OSveenopid02cor"]->Fill(mrec,1/effictotal);
	     if(ksmass >= 0.466 && ksmass <= 0.530){
             histosTH1F["hm4rec2OSveenopid02col"]->Fill(mrec);
             histosTH1F["hm4rec2OSveenopid02colcor"]->Fill(mrec,1/effictotal);
	     }
	     if(charray[0]+charray[1] == 0)
	         {
	       histosTH1F["hm4rec2OS_pi1pi2veenopid02"]->Fill(mrecpi1pi2);
	       histosTH1F["hm4rec2OS_pi3pi4veenopid02"]->Fill(mrecpi3pi4);
	       histosTH2F["hm4dim2OS_pi1pi2_pi3pi4veenopid02"]->Fill(mrecpi1pi2,mrecpi3pi4);
	         }else if(charray[0]+charray[2] == 0){
	       histosTH1F["hm4rec2OS_pi1pi3veenopid02"]->Fill(mrecpi1pi3);
	       histosTH1F["hm4rec2OS_pi2pi4veenopid02"]->Fill(mrecpi2pi4);
	       histosTH2F["hm4dim2OS_pi1pi3_pi2pi4veenopid02"]->Fill(mrecpi1pi3,mrecpi2pi4);
		 }		 
	 // |-(t1+t2)|
         // ...OR
	 if(t1 <= limtA || t2 <= limtA ){histosTH1F["hm4rec2OSveenopid02t12a"]->Fill(mrec);}//0.4
	 if(t1 <= limtB || t2 <= limtB ){histosTH1F["hm4rec2OSveenopid02t12b"]->Fill(mrec);}//0.3
	 if(t1 <= limtC || t2 <= limtC ){histosTH1F["hm4rec2OSveenopid02t12c"]->Fill(mrec);}//0.2
	 if(t1 <= limtD || t2 <= limtD ){histosTH1F["hm4rec2OSveenopid02t12d"]->Fill(mrec);}//0.17
	 if(t1 <= limtE || t2 <= limtE ){histosTH1F["hm4rec2OSveenopid02t12e"]->Fill(mrec);}//0.1
	 if(t1 <= limtF || t2 <= limtF ){histosTH1F["hm4rec2OSveenopid02t12f"]->Fill(mrec);}//0.08
	 if(t1 <= limtG || t2 <= limtG ){histosTH1F["hm4rec2OSveenopid02t12g"]->Fill(mrec);}//0.06
	 if(t1 <= limtH || t2 <= limtH ){histosTH1F["hm4rec2OSveenopid02t12h"]->Fill(mrec);}//0.05
         // ...OR complementary
         if(t1 > limtA || t2 > limtA ){histosTH1F["hm4rec2OSveenopid02t12aC"]->Fill(mrec);}
         if(t1 > limtB || t2 > limtB ){histosTH1F["hm4rec2OSveenopid02t12bC"]->Fill(mrec);}
         if(t1 > limtC || t2 > limtC ){histosTH1F["hm4rec2OSveenopid02t12cC"]->Fill(mrec);}
         if(t1 > limtD || t2 > limtD ){histosTH1F["hm4rec2OSveenopid02t12dC"]->Fill(mrec);}
         if(t1 > limtE || t2 > limtE ){histosTH1F["hm4rec2OSveenopid02t12eC"]->Fill(mrec);}
         if(t1 > limtF || t2 > limtF ){histosTH1F["hm4rec2OSveenopid02t12fC"]->Fill(mrec);}
         if(t1 > limtG || t2 > limtG ){histosTH1F["hm4rec2OSveenopid02t12gC"]->Fill(mrec);}
         if(t1 > limtH || t2 > limtH ){histosTH1F["hm4rec2OSveenopid02t12hC"]->Fill(mrec);}
         // ...AND
	 if(t1 <= limtA && t2 <= limtA ){histosTH1F["hm4rec2OSveenopid02t12aa"]->Fill(mrec);}//0.4
	 if(t1 <= limtB && t2 <= limtB ){histosTH1F["hm4rec2OSveenopid02t12bb"]->Fill(mrec);}//0.3
	 if(t1 <= limtC && t2 <= limtC ){histosTH1F["hm4rec2OSveenopid02t12cc"]->Fill(mrec);}//0.2
	 if(t1 <= limtD && t2 <= limtD ){histosTH1F["hm4rec2OSveenopid02t12dd"]->Fill(mrec);}//0.17
	 if(t1 <= limtE && t2 <= limtE ){histosTH1F["hm4rec2OSveenopid02t12ee"]->Fill(mrec);}//0.1
	 if(t1 <= limtF && t2 <= limtF ){histosTH1F["hm4rec2OSveenopid02t12ff"]->Fill(mrec);}//0.08
	 if(t1 <= limtG && t2 <= limtG ){histosTH1F["hm4rec2OSveenopid02t12gg"]->Fill(mrec);}//0.06
	 if(t1 <= limtH && t2 <= limtH ){histosTH1F["hm4rec2OSveenopid02t12hh"]->Fill(mrec);}//0.05
         // ...AND complementary
         if(t1 > limtA && t2 > limtA ){histosTH1F["hm4rec2OSveenopid02t12aaC"]->Fill(mrec);}
         if(t1 > limtB && t2 > limtB ){histosTH1F["hm4rec2OSveenopid02t12bbC"]->Fill(mrec);}
         if(t1 > limtC && t2 > limtC ){histosTH1F["hm4rec2OSveenopid02t12ccC"]->Fill(mrec);}
         if(t1 > limtD && t2 > limtD ){histosTH1F["hm4rec2OSveenopid02t12ddC"]->Fill(mrec);}
         if(t1 > limtE && t2 > limtE ){histosTH1F["hm4rec2OSveenopid02t12eeC"]->Fill(mrec);}
         if(t1 > limtF && t2 > limtF ){histosTH1F["hm4rec2OSveenopid02t12ffC"]->Fill(mrec);}
         if(t1 > limtG && t2 > limtG ){histosTH1F["hm4rec2OSveenopid02t12ggC"]->Fill(mrec);}
         if(t1 > limtH && t2 > limtH ){histosTH1F["hm4rec2OSveenopid02t12hhC"]->Fill(mrec);}

	 // dphi
	 if(m1370){ histosTH1F["hdphig1370nopid"]->Fill(TOTEMdphi*180/TMath::Pi());}
	 if(m1549){ histosTH1F["hdphig1549nopid"]->Fill(TOTEMdphi*180/TMath::Pi());}
	 if(m1731){ histosTH1F["hdphig1731nopid"]->Fill(TOTEMdphi*180/TMath::Pi());}
	 // dphi ...effic correction
	 if(m1370){ histosTH1F["hdphig1370nopidcor"]->Fill(TOTEMdphi*180/TMath::Pi(),1/effictotal);}
	 if(m1549){ histosTH1F["hdphig1549nopidcor"]->Fill(TOTEMdphi*180/TMath::Pi(),1/effictotal);}
	 if(m1731){ histosTH1F["hdphig1731nopidcor"]->Fill(TOTEMdphi*180/TMath::Pi(),1/effictotal);}

	 // (dphi,mrec) map
	 histosTH2F["h2dimdphigmnopid"]->Fill(TOTEMdphi*180/TMath::Pi(),mrec);
	 // (t12,mrec) map   .....................................
	 histosTH2F["h2dimt12mnopid"]->Fill(t12,mrec);
	 
	 //}//...end of PID=pions
	     } //end of nvtx=0 nks=2

	     	     
	     //...veeno01
             if(nvtx==0 && nks==1){

	       //...only pions ...BETTER USING PID
	       //if(pidarray[0]==2 && pidarray[1]==2 && pidarray[2]==2 && pidarray[3]==2){
	       ////if(pidv0array[0]==2 && pidv0array[1]==2 && pidv0array[2]==2 && pidv0array[3]==2){
	       
	     histosTH1F["hm4rec2OSveenopid01"]->Fill(mrec);      
		 if(charray[0]+charray[1] == 0)
	         {
	       histosTH1F["hm4rec2OS_pi1pi2veenopid01"]->Fill(mrecpi1pi2);
	       histosTH1F["hm4rec2OS_pi3pi4veenopid01"]->Fill(mrecpi3pi4);
	       histosTH2F["hm4dim2OS_pi1pi2_pi3pi4veenopid01"]->Fill(mrecpi1pi2,mrecpi3pi4);
	         }else if(charray[0]+charray[2] == 0){
	       histosTH1F["hm4rec2OS_pi1pi3veenopid01"]->Fill(mrecpi1pi3);
	       histosTH1F["hm4rec2OS_pi2pi4veenopid01"]->Fill(mrecpi2pi4);
	       histosTH2F["hm4dim2OS_pi1pi3_pi2pi4veenopid01"]->Fill(mrecpi1pi3,mrecpi2pi4);
		 }
	 // |-(t1+t2)|
         // ...OR
	 if(t1 <= limtA || t2 <= limtA ){histosTH1F["hm4rec2OSveenopid01t12a"]->Fill(mrec);}//0.4
	 if(t1 <= limtB || t2 <= limtB ){histosTH1F["hm4rec2OSveenopid01t12b"]->Fill(mrec);}//0.3
	 if(t1 <= limtC || t2 <= limtC ){histosTH1F["hm4rec2OSveenopid01t12c"]->Fill(mrec);}//0.2
	 if(t1 <= limtD || t2 <= limtD ){histosTH1F["hm4rec2OSveenopid01t12d"]->Fill(mrec);}//0.17
	 if(t1 <= limtE || t2 <= limtE ){histosTH1F["hm4rec2OSveenopid01t12e"]->Fill(mrec);}//0.1
	 if(t1 <= limtF || t2 <= limtF ){histosTH1F["hm4rec2OSveenopid01t12f"]->Fill(mrec);}//0.08
	 if(t1 <= limtG || t2 <= limtG ){histosTH1F["hm4rec2OSveenopid01t12g"]->Fill(mrec);}//0.06
	 if(t1 <= limtH || t2 <= limtH ){histosTH1F["hm4rec2OSveenopid01t12h"]->Fill(mrec);}//0.05
         // ...OR complementary
         if(t1 > limtA || t2 > limtA ){histosTH1F["hm4rec2OSveenopid01t12aC"]->Fill(mrec);}
         if(t1 > limtB || t2 > limtB ){histosTH1F["hm4rec2OSveenopid01t12bC"]->Fill(mrec);}
         if(t1 > limtC || t2 > limtC ){histosTH1F["hm4rec2OSveenopid01t12cC"]->Fill(mrec);}
         if(t1 > limtD || t2 > limtD ){histosTH1F["hm4rec2OSveenopid01t12dC"]->Fill(mrec);}
         if(t1 > limtE || t2 > limtE ){histosTH1F["hm4rec2OSveenopid01t12eC"]->Fill(mrec);}
         if(t1 > limtF || t2 > limtF ){histosTH1F["hm4rec2OSveenopid01t12fC"]->Fill(mrec);}
         if(t1 > limtG || t2 > limtG ){histosTH1F["hm4rec2OSveenopid01t12gC"]->Fill(mrec);}
         if(t1 > limtH || t2 > limtH ){histosTH1F["hm4rec2OSveenopid01t12hC"]->Fill(mrec);}
         // ...AND
	 if(t1 <= limtA && t2 <= limtA ){histosTH1F["hm4rec2OSveenopid01t12aa"]->Fill(mrec);}//0.4
	 if(t1 <= limtB && t2 <= limtB ){histosTH1F["hm4rec2OSveenopid01t12bb"]->Fill(mrec);}//0.3
	 if(t1 <= limtC && t2 <= limtC ){histosTH1F["hm4rec2OSveenopid01t12cc"]->Fill(mrec);}//0.2
	 if(t1 <= limtD && t2 <= limtD ){histosTH1F["hm4rec2OSveenopid01t12dd"]->Fill(mrec);}//0.17
	 if(t1 <= limtE && t2 <= limtE ){histosTH1F["hm4rec2OSveenopid01t12ee"]->Fill(mrec);}//0.1
	 if(t1 <= limtF && t2 <= limtF ){histosTH1F["hm4rec2OSveenopid01t12ff"]->Fill(mrec);}//0.08
	 if(t1 <= limtG && t2 <= limtG ){histosTH1F["hm4rec2OSveenopid01t12gg"]->Fill(mrec);}//0.06
	 if(t1 <= limtH && t2 <= limtH ){histosTH1F["hm4rec2OSveenopid01t12hh"]->Fill(mrec);}//0.05
         // ...AND complementary
         if(t1 > limtA && t2 > limtA ){histosTH1F["hm4rec2OSveenopid01t12aaC"]->Fill(mrec);}
         if(t1 > limtB && t2 > limtB ){histosTH1F["hm4rec2OSveenopid01t12bbC"]->Fill(mrec);}
         if(t1 > limtC && t2 > limtC ){histosTH1F["hm4rec2OSveenopid01t12ccC"]->Fill(mrec);}
         if(t1 > limtD && t2 > limtD ){histosTH1F["hm4rec2OSveenopid01t12ddC"]->Fill(mrec);}
         if(t1 > limtE && t2 > limtE ){histosTH1F["hm4rec2OSveenopid01t12eeC"]->Fill(mrec);}
         if(t1 > limtF && t2 > limtF ){histosTH1F["hm4rec2OSveenopid01t12ffC"]->Fill(mrec);}
         if(t1 > limtG && t2 > limtG ){histosTH1F["hm4rec2OSveenopid01t12ggC"]->Fill(mrec);}
         if(t1 > limtH && t2 > limtH ){histosTH1F["hm4rec2OSveenopid01t12hhC"]->Fill(mrec);}
		 
	 //}//...end of no PID=pions
	      } //end of nvtx=0 nks=1
	     
	     } //...end of isKshort
	   } //CC...end of totalcharge=0
 

	   //-----------end of cut 8
	   

  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

	 //...cut k1  ...K+K-K+K-

	 double sigkkphi  = 0.010;
	 double sigkkphi2 = 0.020;
	 double sigkkphi3 = 0.030;


	 if(totcharge==0){

	   //...kaons
	   
	   histosTH1F["hm4rec2OSk"]->Fill(mrecKKKK);

	   //..........................???????????????????  1 or 2 primaries ?
	   if(nvtx==1){
	     //
	   histosTH1F["hm4rec2OS2k"]->Fill(mrecKKKK);
	   
           if(pidarray[0]==3 && pidarray[1]==3 && pidarray[2]==3 && pidarray[3]==3){
	   histosTH1F["hm4rec2OS2kpid"]->Fill(mrecKKKK);
	                } //...end of PID=kaons

	   //1234
	   if(charray[0]+charray[1] == 0)
	          {
	             histosTH1F["hm4rec2OS_k1k2m"]->Fill(mreck1k2);
		     //
                     ////histosTH1F["h2OSpt12k"]->Fill(ptk1k2);
                     ////histosTH1F["h2OSeta12k"]->Fill(etak1k2);
                     ////histosTH2F["h2dim2OSpteta12k"]->Fill(ptk1k2,etak1k2);
		     //
		     histosTH1F["hm4rec2OS_k3k4m"]->Fill(mreck3k4);
		     //
                     ////histosTH1F["h2OSpt34k"]->Fill(ptk3k4);
                     ////histosTH1F["h2OSeta34k"]->Fill(etak3k4);
                     ////histosTH2F["h2dim2OSpteta34k"]->Fill(ptk3k4,etak3k4);		     
		     // 
	               histosTH1F["hm4rec2OSm1234k"]->Fill(mrecKKKK);
		       histosTH2F["h2dim2OSm12x34k"]->Fill(mreck1k2,mreck3k4);

	   // ----------- phi(1020) -----------study  
	   // 1234 for M(phi) shape study ...beased on inspection of M(12)xM(34).
	   // | M(k1k2) - M(phi) | < 25 MeV win=50MeV
	             if( TMath::Abs(mreck1k2 - m_phi) < 0.025 ){
		     histosTH1F["hm4rec2OS34phi"]->Fill(mreck3k4);
 		     }
		       
		     //
		     // phi(1020)
		     // Gamma = 10.95 MeV
		     // sigma = 4.65 MeV
	             //
		     //...final study phi(1020)
    		     //...| M(K+K-) - M(phi) | < 4.65 MeV ==> mass window: 2*sigma = 9.3 MeV
	             if( ( TMath::Abs(mreck1k2 - m_phi) < 0.00465 ) &&
			 ( TMath::Abs(mreck3k4 - m_phi) < 0.00465 ) ){
		       	histosTH1F["hm4rec2OSk123405win9"]->Fill(mrecKKKK);
		     }
    		     //...| M(K+K-) - M(phi) | < 9.3 MeV ==> mass window: 4*sigma = 18.6 MeV
	             if( ( TMath::Abs(mreck1k2 - m_phi) < 0.0093 ) &&
			 ( TMath::Abs(mreck3k4 - m_phi) < 0.0093 ) ){
		       	histosTH1F["hm4rec2OSk123405win18"]->Fill(mrecKKKK);
		     }
		     //
		     //
		     // X(1070)
		     // Gamma = 10.02 MeV
		     // sigma = 4.25 MeV
	             //
		     //...final study X(1070)
    		     //...| M(K+K-) - M(1070) | < 4.25 MeV ==> mass window: 2*sigma = 8.5 MeV
	             if( ( TMath::Abs(mreck1k2 - m_x1070) < 0.00425 ) &&
			 ( TMath::Abs(mreck3k4 - m_x1070) < 0.00425 ) ){
		       	histosTH1F["hm4rec2OSk123470win8"]->Fill(mrecKKKK);
		     }
    		     //...| M(K+K-) - M(1070) | < 8.5 MeV ==> mass window: 4*sigma = 17.0 MeV
	             if( ( TMath::Abs(mreck1k2 - m_x1070) < 0.0085 ) &&
			 ( TMath::Abs(mreck3k4 - m_x1070) < 0.0085 ) ){
		       	histosTH1F["hm4rec2OSk123470win17"]->Fill(mrecKKKK);
		     }
		     //
		     //
		     // X(1235)
		     // Gamma = 87.0 MeV
		     // sigma = 36.9 MeV
	             //
		     //...final study X(1235)
    		     //...| M(K+K-) - M(1235) | < 36.9 MeV ==> mass window: 2*sigma = 73.9 MeV
	             if( ( TMath::Abs(mreck1k2 - m_x1235) < 0.0369 ) &&
			 ( TMath::Abs(mreck3k4 - m_x1235) < 0.0369 ) ){
		       	histosTH1F["hm4rec2OSk1234235win74"]->Fill(mrecKKKK);
		     }
    		     //...| M(K+K-) - M(1235) | < 73.9 MeV ==> mass window: 4*sigma = 147.8 MeV
	             if( ( TMath::Abs(mreck1k2 - m_x1235) < 0.0739 ) &&
			 ( TMath::Abs(mreck3k4 - m_x1235) < 0.0739 ) ){
		       	histosTH1F["hm4rec2OSk1234235win148"]->Fill(mrecKKKK);
		     }
		     //

		     
		     //... phi
		     //...| M(K+K-) - M(phi) | < 10 MeV = 1*sigma, mass window = 2*sigma = 20 MeV
	             if( ( TMath::Abs(mreck1k2 - m_phi) < sigkkphi ) &&
			 ( TMath::Abs(mreck3k4 - m_phi) < sigkkphi ) ){
		       	histosTH1F["hm4rec2OSk1234sig"]->Fill(mrecKKKK);
                        if(pidarray[0]==3 && pidarray[1]==3 && pidarray[2]==3 && pidarray[3]==3){
			histosTH1F["hm4rec2OSk1234sigpid"]->Fill(mrecKKKK);}
		     }
		     //...| M(K+K-) - M(phi) | < 20 MeV = 1*sigma, mass window = 2*sigma = 40 MeV
	             if( ( TMath::Abs(mreck1k2 - m_phi) < sigkkphi2 ) &&
			 ( TMath::Abs(mreck3k4 - m_phi) < sigkkphi2 ) ){
		       	histosTH1F["hm4rec2OSk1234sig2"]->Fill(mrecKKKK);
                        if(pidarray[0]==3 && pidarray[1]==3 && pidarray[2]==3 && pidarray[3]==3){
			histosTH1F["hm4rec2OSk1234sig2pid"]->Fill(mrecKKKK);}
		     }
		     //...| M(K+K-) - M(phi) | < 30 MeV = 1*sigma, mass window = 2*sigma = 60 MeV
	             if( ( TMath::Abs(mreck1k2 - m_phi) < sigkkphi3 ) &&
			 ( TMath::Abs(mreck3k4 - m_phi) < sigkkphi3 ) ){
		       	histosTH1F["hm4rec2OSk1234sig3"]->Fill(mrecKKKK);
                        if(pidarray[0]==3 && pidarray[1]==3 && pidarray[2]==3 && pidarray[3]==3){
			histosTH1F["hm4rec2OSk1234sig3pid"]->Fill(mrecKKKK);}
		     }
		  }
	   //1324
	   if(charray[0]+charray[2] == 0)
		  {
	             histosTH1F["hm4rec2OS_k1k3m"]->Fill(mreck1k3);
                     //
		     ////histosTH1F["h2OSpt13k"]->Fill(ptk1k3);
                     ////histosTH1F["h2OSeta13k"]->Fill(etak1k3);
                     ////histosTH2F["h2dim2OSpteta13k"]->Fill(ptk1k3,etak1k3);
 		     //
		     histosTH1F["hm4rec2OS_k2k4m"]->Fill(mreck2k4);
		     //
		     ////histosTH1F["h2OSpt24k"]->Fill(ptk2k4);
                     ////histosTH1F["h2OSeta24k"]->Fill(etak2k4);
                     ////histosTH2F["h2dim2OSpteta24k"]->Fill(ptk2k4,etak2k4);		     
 		     //
	               histosTH1F["hm4rec2OSm1324k"]->Fill(mrecKKKK);
		       histosTH2F["h2dim2OSm13x24k"]->Fill(mreck1k3,mreck2k4);
		     //

		     //
		     // phi(1020)
		     // Gamma = 10.95 MeV
		     // sigma = 4.65 MeV
	             //
		     //...final study phi(1020)
    		     //...| M(K+K-) - M(phi) | < 4.65 MeV ==> mass window: 2*sigma = 9.3 MeV
	             if( ( TMath::Abs(mreck1k3 - m_phi) < 0.00465 ) &&
			 ( TMath::Abs(mreck2k4 - m_phi) < 0.00465 ) ){
		       	histosTH1F["hm4rec2OSk132405win9"]->Fill(mrecKKKK);
		     }
    		     //...| M(K+K-) - M(phi) | < 9.3 MeV ==> mass window: 4*sigma = 18.6 MeV
	             if( ( TMath::Abs(mreck1k3 - m_phi) < 0.0093 ) &&
			 ( TMath::Abs(mreck2k4 - m_phi) < 0.0093 ) ){
		       	histosTH1F["hm4rec2OSk132405win18"]->Fill(mrecKKKK);
		     }
		     //
		     //
		     // X(1070)
		     // Gamma = 10.02 MeV
		     // sigma = 4.25 MeV
	             //
		     //...final study X(1070)
    		     //...| M(K+K-) - M(1070) | < 4.25 MeV ==> mass window: 2*sigma = 8.5 MeV
	             if( ( TMath::Abs(mreck1k3 - m_x1070) < 0.00425 ) &&
			 ( TMath::Abs(mreck2k4 - m_x1070) < 0.00425 ) ){
		       	histosTH1F["hm4rec2OSk132470win8"]->Fill(mrecKKKK);
		     }
    		     //...| M(K+K-) - M(1070) | < 8.5 MeV ==> mass window: 4*sigma = 17.0 MeV
	             if( ( TMath::Abs(mreck1k3 - m_x1070) < 0.0085 ) &&
			 ( TMath::Abs(mreck2k4 - m_x1070) < 0.0085 ) ){
		       	histosTH1F["hm4rec2OSk132470win17"]->Fill(mrecKKKK);
		     }
		     //
		     //
		     // X(1235)
		     // Gamma = 87.0 MeV
		     // sigma = 36.9 MeV
	             //
		     //...final study X(1235)
    		     //...| M(K+K-) - M(1235) | < 36.9 MeV ==> mass window: 2*sigma = 73.9 MeV
	             if( ( TMath::Abs(mreck1k3 - m_x1235) < 0.0369 ) &&
			 ( TMath::Abs(mreck2k4 - m_x1235) < 0.0369 ) ){
		       	histosTH1F["hm4rec2OSk1324235win74"]->Fill(mrecKKKK);
		     }
    		     //...| M(K+K-) - M(1235) | < 73.9 MeV ==> mass window: 4*sigma = 147.8 MeV
	             if( ( TMath::Abs(mreck1k3 - m_x1235) < 0.0739 ) &&
			 ( TMath::Abs(mreck2k4 - m_x1235) < 0.0739 ) ){
		       	histosTH1F["hm4rec2OSk1324235win148"]->Fill(mrecKKKK);
		     }
		     //



		     //...phi
		     //...| M(K+K-) - M(phi) | < 10 MeV = 1*sigma, mass window = 2*sigma = 20 MeV
	             if( ( TMath::Abs(mreck1k3 - m_phi) < sigkkphi ) &&
			 ( TMath::Abs(mreck2k4 - m_phi) < sigkkphi ) ){
		       	histosTH1F["hm4rec2OSk1324sig"]->Fill(mrecKKKK);
                        if(pidarray[0]==3 && pidarray[1]==3 && pidarray[2]==3 && pidarray[3]==3){
			histosTH1F["hm4rec2OSk1324sigpid"]->Fill(mrecKKKK);}
		     }
		     //...| M(K+K-) - M(phi) | < 20 MeV = 1*sigma, mass window = 2*sigma = 40 MeV
	             if( ( TMath::Abs(mreck1k3 - m_phi) < sigkkphi2 ) &&
			 ( TMath::Abs(mreck2k4 - m_phi) < sigkkphi2 ) ){
		       	histosTH1F["hm4rec2OSk1324sig2"]->Fill(mrecKKKK);
                        if(pidarray[0]==3 && pidarray[1]==3 && pidarray[2]==3 && pidarray[3]==3){
			histosTH1F["hm4rec2OSk1324sig2pid"]->Fill(mrecKKKK);}
		     }
		     //...| M(K+K-) - M(phi) | < 30 MeV = 1*sigma, mass window = 2*sigma = 60 MeV
	             if( ( TMath::Abs(mreck1k3 - m_phi) < sigkkphi3 ) &&
			 ( TMath::Abs(mreck2k4 - m_phi) < sigkkphi3 ) ){
		       	histosTH1F["hm4rec2OSk1324sig3"]->Fill(mrecKKKK);
                        if(pidarray[0]==3 && pidarray[1]==3 && pidarray[2]==3 && pidarray[3]==3){
			histosTH1F["hm4rec2OSk1324sig3pid"]->Fill(mrecKKKK);}
		     }
		  }
	   //1423
	   if(charray[0]+charray[3] == 0)
	          {	    
	             histosTH1F["hm4rec2OS_k1k4m"]->Fill(mreck1k4);
                     //
		     ////histosTH1F["h2OSpt14k"]->Fill(ptk1k4);
                     ////histosTH1F["h2OSeta14k"]->Fill(etak1k4);
                     ////histosTH2F["h2dim2OSpteta14k"]->Fill(ptk1k4,etak1k4);
		     //
		     histosTH1F["hm4rec2OS_k2k3m"]->Fill(mreck2k3);
                     //
		     ////histosTH1F["h2OSpt23k"]->Fill(ptk2k3);
                     ////histosTH1F["h2OSeta23k"]->Fill(etak2k3);
                     ////histosTH2F["h2dim2OSpteta23k"]->Fill(ptk2k3,etak2k3);		     
		     // 
	               histosTH1F["hm4rec2OSm1423k"]->Fill(mrecKKKK);
		       histosTH2F["h2dim2OSm14x23k"]->Fill(mreck1k4,mreck2k3);

	   // ----------- phi(1020) -----------study  
	   // 1423 for M(phi) shape study ...based on inspection of M(14)xM(23).   ...tasting!!
	   // | M(k2k3) - M(phi) | < 20 MeV win=40MeV
	             if( TMath::Abs(mreck2k3 - m_phi) < 0.020 ){
		     histosTH1F["hm4rec2OS14phi"]->Fill(mreck1k4);
 		     }
		     
	   // ----------- X(1070) -----------study A : M(1090) ...center 1090      ...see the map!
	   // 1423 for M(1090) shape study ...based on inspection of M(14)xM(23).
	   // | M(k1k4) - M(1090) | < 30 MeV win=60MeV
	             if( TMath::Abs(mreck1k4 - 1.090 ) < 0.030 ){
		     histosTH1F["hm4rec2OS23exo70"]->Fill(mreck2k3);
 		     }
	   // ----------- X(1070) -----------study B : M(1070) ...center 1070
	   // 1423 for M(1070) shape study ...based on inspection of M(14)xM(23).
	   // | M(k2k3) - M(1070) | < 25 MeV win=50MeV
	             if( TMath::Abs(mreck2k3 - 1.070 ) < 0.025 ){
		     histosTH1F["hm4rec2OS14exo70"]->Fill(mreck1k4);
 		     }
		       
	   // ----------- exotic at 1220 unkown ----------study
	   // 1423 for M(1220) shape study ...based on inspection of M(14)xM(23).
	   // | M(k2k3) - M(1220) | < 90 MeV win=180MeV
	             if( TMath::Abs(mreck2k3 - 1.220 ) < 0.090 ){
		     histosTH1F["hm4rec2OS14exo1220"]->Fill(mreck1k4);
 		     }


		     //
		     // phi(1020)
		     // Gamma = 10.95 MeV
		     // sigma = 4.65 MeV
	             //
		     //...final study phi(1020)
    		     //...| M(K+K-) - M(phi) | < 4.65 MeV ==> mass window: 2*sigma = 9.3 MeV
	             if( ( TMath::Abs(mreck1k4 - m_phi) < 0.00465 ) &&
			 ( TMath::Abs(mreck2k3 - m_phi) < 0.00465 ) ){
		       	histosTH1F["hm4rec2OSk142305win9"]->Fill(mrecKKKK);
		     }
    		     //...| M(K+K-) - M(phi) | < 9.3 MeV ==> mass window: 4*sigma = 18.6 MeV
	             if( ( TMath::Abs(mreck1k4 - m_phi) < 0.0093 ) &&
			 ( TMath::Abs(mreck2k3 - m_phi) < 0.0093 ) ){
		       	histosTH1F["hm4rec2OSk142305win18"]->Fill(mrecKKKK);
		     }
		     //
		     //
		     // X(1070)
		     // Gamma = 10.02 MeV
		     // sigma = 4.25 MeV
	             //
		     //...final study X(1070)
    		     //...| M(K+K-) - M(1070) | < 4.25 MeV ==> mass window: 2*sigma = 8.5 MeV
	             if( ( TMath::Abs(mreck1k4 - m_x1070) < 0.00425 ) &&
			 ( TMath::Abs(mreck2k3 - m_x1070) < 0.00425 ) ){
		       	histosTH1F["hm4rec2OSk142370win8"]->Fill(mrecKKKK);
		     }
    		     //...| M(K+K-) - M(1070) | < 8.5 MeV ==> mass window: 4*sigma = 17.0 MeV
	             if( ( TMath::Abs(mreck1k4 - m_x1070) < 0.0085 ) &&
			 ( TMath::Abs(mreck2k3 - m_x1070) < 0.0085 ) ){
		       	histosTH1F["hm4rec2OSk142370win17"]->Fill(mrecKKKK);
		     }
		     //
		     //
		     // X(1235)
		     // Gamma = 87.0 MeV
		     // sigma = 36.9 MeV
	             //
		     //...final study X(1235)
    		     //...| M(K+K-) - M(1235) | < 36.9 MeV ==> mass window: 2*sigma = 73.9 MeV
	             if( ( TMath::Abs(mreck1k4 - m_x1235) < 0.0369 ) &&
			 ( TMath::Abs(mreck2k3 - m_x1235) < 0.0369 ) ){
		       	histosTH1F["hm4rec2OSk1423235win74"]->Fill(mrecKKKK);
		     }
    		     //...| M(K+K-) - M(1235) | < 73.9 MeV ==> mass window: 4*sigma = 147.8 MeV
	             if( ( TMath::Abs(mreck1k4 - m_x1235) < 0.0739 ) &&
			 ( TMath::Abs(mreck2k3 - m_x1235) < 0.0739 ) ){
		       	histosTH1F["hm4rec2OSk1423235win148"]->Fill(mrecKKKK);
		     }
		     //


		     
		     //...phi
		     //...| M(K+K-) - M(phi) | < 10 MeV = 1*sigma, mass window = 2*sigma = 20 MeV
	             if( ( TMath::Abs(mreck1k4 - m_phi) < sigkkphi ) &&
			 ( TMath::Abs(mreck2k3 - m_phi) < sigkkphi ) ){
		       	histosTH1F["hm4rec2OSk1423sig"]->Fill(mrecKKKK);
                        if(pidarray[0]==3 && pidarray[1]==3 && pidarray[2]==3 && pidarray[3]==3){
			histosTH1F["hm4rec2OSk1423sigpid"]->Fill(mrecKKKK);}
		     }
		     //...| M(K+K-) - M(phi) | < 20 MeV = 1*sigma, mass window = 2*sigma = 40 MeV
	             if( ( TMath::Abs(mreck1k4 - m_phi) < sigkkphi2 ) &&
			 ( TMath::Abs(mreck2k3 - m_phi) < sigkkphi2 ) ){
		       	histosTH1F["hm4rec2OSk1423sig2"]->Fill(mrecKKKK);
                        if(pidarray[0]==3 && pidarray[1]==3 && pidarray[2]==3 && pidarray[3]==3){
			histosTH1F["hm4rec2OSk1423sig2pid"]->Fill(mrecKKKK);}
		     }
		     //...| M(K+K-) - M(phi) | < 30 MeV = 1*sigma, mass window = 2*sigma = 60 MeV
	             if( ( TMath::Abs(mreck1k4 - m_phi) < sigkkphi3 ) &&
			 ( TMath::Abs(mreck2k3 - m_phi) < sigkkphi3 ) ){
		       	histosTH1F["hm4rec2OSk1423sig3"]->Fill(mrecKKKK);
                        if(pidarray[0]==3 && pidarray[1]==3 && pidarray[2]==3 && pidarray[3]==3){
			histosTH1F["hm4rec2OSk1423sig3pid"]->Fill(mrecKKKK);}
		     }

		     //exotics
		     //...exotic at 1230MeV unkown ...TEST
		     //...| M(K+K-) - M(phi) | < 57.5 MeV = 1*sigma, mass window  2*sigma = 115 MeV
	             if( ( TMath::Abs(mreck1k4 - 1.230) < 0.0575 ) &&
			 ( TMath::Abs(mreck2k3 - 1.230) < 0.0575 ) ){
		       	histosTH1F["hm4rec2OS1423exo"]->Fill(mrecKKKK);
                        if(pidarray[0]==3 && pidarray[1]==3 && pidarray[2]==3 && pidarray[3]==3){
			histosTH1F["hm4rec2OS1423exopid"]->Fill(mrecKKKK);}
		     }
		     //...exotic at 1085MeV unkown ...TEST
		     //...| M(K+K-) - M(phi) | < 40 MeV = 1*sigma, mass window  2*sigma = 80 MeV
	             if( ( TMath::Abs(mreck1k4 - 1.085) < 0.040 ) &&
			 ( TMath::Abs(mreck2k3 - 1.085) < 0.040 ) ){
		       	histosTH1F["hm4rec2OS1423exo85"]->Fill(mrecKKKK);
                        if(pidarray[0]==3 && pidarray[1]==3 && pidarray[2]==3 && pidarray[3]==3){
			histosTH1F["hm4rec2OS1423exo85pid"]->Fill(mrecKKKK);}
		     }
		  }

	   } //...end of nvtx=1
      	 } //...end of cut k1 totalcharge



	 /*

	   //...K*Kbar* channel

	   //1234
           k1pi2k3pi4Rec  = k1  + pi2 + k3 + pi4;
           k1pi2pi3k4Rec  = k1  + pi2 + pi3 + k4;
           pi1k2k3pi4Rec  = pi1  + k2 + k3 + pi4;
           pi1k2pi3k4Rec  = pi1  + k2 + pi3 + k4;
	   //1324
           k1pi3k2pi4Rec  = k1  + pi3 + k2 + pi4;
           k1pi3pi2k4Rec  = k1  + pi3 + pi2 + k4;
           pi1k3k2pi4Rec  = pi1  + k3 + k2 + pi4;
           pi1k3pi2k4Rec  = pi1  + k3 + pi2 + k4;
	   //1423
           k1pi4k2pi3Rec  = k1  + pi4 + k2 + pi3;
           k1pi4pi2k3Rec  = k1  + pi4 + pi2 + k3;
           pi1k4k2pi3Rec  = pi1  + k4 + k2 + pi3;
           pi1k4pi2k3Rec  = pi1  + k4 + pi2 + k3;

	 */


	 
//...cut K*Kbar* or K+pi-K-pi+    ...nvtx==1 nks=0 nlam=0
       
	 if(totcharge==0){

	   if(nvtx==1 && nks==0 && nlam==0){

//...1234
	     //k1pi2k3pi4
		 
 	     if( charray[0]== +1 && charray[1]== -1 && charray[2]== -1 && charray[3]== +1 &&
		 pidarray[0]== 3 && pidarray[1]== 2 && pidarray[2]== 3 && pidarray[3]== 2 ){

	       mrecKstar    = mreck1pi2;
	       mrecKstarbar = mreck3pi4;
	       
	         histosTH1F["hm4rec2OSkstar"]->Fill(mrecKstar);
	         histosTH1F["hm4rec2OSkstarbar"]->Fill(mrecKstarbar);

                 //histosTH1F["hm4rec2OSk1pi2"]->Fill(mreck1pi2);
		 //histosTH1F["hm4rec2OSk3pi4b"]->Fill(mreck3pi4);
	         //histosTH1F["hm4rec2OSk1p2k3p4b"]->Fill(mreck1pi2k3pi4);
		 
		 //histosTH2F["h2dim2OSk1p2k3p4b"]->Fill(mreck1pi2,mreck3pi4);
		 histosTH2F["h2dim2OSKK"]->Fill(mrecKstar,mrecKstarbar);

		 //...K*Kbar* study 0.85 - 1.05 GeV at Kbar*
		 //
		 //...| M(Kbar*) - M(895) | < 50.0 MeV ==> 100.0 MeV window ...study
		 if( ( TMath::Abs(mrecKstarbar - m_kstar) < 0.050 ) ){
		   histosTH1F["hm4rec2OSKK100"]->Fill(mrecKstar);
		 }
		 //

		 // Gamma = 34.3 MeV
		 // sigma = 14.6 MeV
		 // 2*sigma = 29.1 MeV
		 // 4*sigma = 58.3 MeV
		 
		 //...| M(K*Kbar*) - M(895) | < 14.6 MeV ==> mass window: 2*sigma = 29.1 MeV
		 if( ( TMath::Abs(mrecKstar - m_kstar) < 0.0146 ) &&
		     ( TMath::Abs(mrecKstarbar - m_kstar) < 0.0146 ) ){
		   histosTH1F["hm4rec2OSKKwin29"]->Fill(mreck1pi2k3pi4);
		 }
		 //...| M(K*Kbar*) - M(895) | < 29.1 MeV ==> mass window: 4*sigma = 58.3 MeV
		 if( ( TMath::Abs(mrecKstar - m_kstar) < 0.0291 ) &&
		     ( TMath::Abs(mrecKstarbar - m_kstar) < 0.0291 ) ){
		   histosTH1F["hm4rec2OSKKwin58"]->Fill(mreck1pi2k3pi4);
		 }
	     }

 	     if( charray[0]== -1 && charray[1]== +1 && charray[2]== +1 && charray[3]== -1 &&
		 pidarray[0]== 3 && pidarray[1]== 2 && pidarray[2]== 3 && pidarray[3]== 2 ){

	       mrecKstarbar = mreck1pi2;
	       mrecKstar    = mreck3pi4;
	       
	         histosTH1F["hm4rec2OSkstar"]->Fill(mrecKstar);
	         histosTH1F["hm4rec2OSkstarbar"]->Fill(mrecKstarbar);

                 //histosTH1F["hm4rec2OSk1pi2b"]->Fill(mreck1pi2);
		 //histosTH1F["hm4rec2OSk3pi4"]->Fill(mreck3pi4);
	         //histosTH1F["hm4rec2OSk1p2bk3p4"]->Fill(mreck1pi2k3pi4);
		 
		 //histosTH2F["h2dim2OSk1p2bk3p4"]->Fill(mreck1pi2,mreck3pi4);
		 histosTH2F["h2dim2OSKK"]->Fill(mrecKstar,mrecKstarbar);

		 //...K*Kbar* study 0.8 - 1.0 GeV at Kbar*
		 //
		 //...| M(Kbar*) - M(895) | < 50.0 MeV ==> 100.0 MeV window ...study
		 if( ( TMath::Abs(mrecKstarbar - m_kstar) < 0.050 ) ){
		   histosTH1F["hm4rec2OSKK100"]->Fill(mrecKstar);
		 }
		 //
		 //...| M(K*Kbar*) - M(895) | < 14.6 MeV ==> mass window: 2*sigma = 29.1 MeV
		 if( ( TMath::Abs(mrecKstar - m_kstar) < 0.0146 ) &&
		     ( TMath::Abs(mrecKstarbar - m_kstar) < 0.0146 ) ){
		   histosTH1F["hm4rec2OSKKwin29"]->Fill(mreck1pi2k3pi4);
		 }
		 //...| M(K*Kbar*) - M(895) | < 29.1 MeV ==> mass window: 4*sigma = 58.3 MeV
		 if( ( TMath::Abs(mrecKstar - m_kstar) < 0.0291 ) &&
		     ( TMath::Abs(mrecKstarbar - m_kstar) < 0.0291 ) ){
		   histosTH1F["hm4rec2OSKKwin58"]->Fill(mreck1pi2k3pi4);
		 }
	     }


	     //k1pi2pi3k4
	     
  	     if( charray[0]== +1 && charray[1]== -1 && charray[2]== +1 && charray[3]== -1 &&
		 pidarray[0]== 3 && pidarray[1]== 2 && pidarray[2]== 2 && pidarray[3]== 3 ){
	       
	         mrecKstar    = mreck1pi2;
	         mrecKstarbar = mrecpi3k4;

	         histosTH1F["hm4rec2OSkstar"]->Fill(mrecKstar);
	         histosTH1F["hm4rec2OSkstarbar"]->Fill(mrecKstarbar);

	         //histosTH1F["hm4rec2OSk1pi2"]->Fill(mreck1pi2);
		 //histosTH1F["hm4rec2OSpi3k4b"]->Fill(mrecpi3k4);
	         //histosTH1F["hm4rec2OSk1p2p3k4b"]->Fill(mreck1pi2pi3k4);
		 
		 //histosTH2F["h2dim2OSk1p2p3k4b"]->Fill(mreck1pi2,mrecpi3k4);
		 histosTH2F["h2dim2OSKK"]->Fill(mrecKstar,mrecKstarbar);
	         
		 //...K*Kbar* study 0.8 - 1.0 GeV at Kbar*
		 //
		 //...| M(Kbar*) - M(895) | < 50.0 MeV ==> 100.0 MeV window ...study
		 if( ( TMath::Abs(mrecKstarbar - m_kstar) < 0.050 ) ){
		   histosTH1F["hm4rec2OSKK100"]->Fill(mrecKstar);
		 }
		 //
		 //...| M(K*Kbar*) - M(895) | < 14.6 MeV ==> mass window: 2*sigma = 29.1 MeV 
		 if( ( TMath::Abs(mrecKstar - m_kstar) < 0.0146 ) &&
		     ( TMath::Abs(mrecKstarbar - m_kstar) < 0.0146 ) ){
		   histosTH1F["hm4rec2OSKKwin29"]->Fill(mreck1pi2pi3k4);
		 }
		 //...| M(K*Kbar*) - M(895) | < 29.1 MeV ==> mass window: 4*sigma = 58.3 MeV 
		 if( ( TMath::Abs(mrecKstar - m_kstar) < 0.0291 ) &&
		     ( TMath::Abs(mrecKstarbar - m_kstar) < 0.0291 ) ){
		   histosTH1F["hm4rec2OSKKwin58"]->Fill(mreck1pi2pi3k4);
		 }
	     }
    
  	     if( charray[0]== -1 && charray[1]== +1 && charray[2]== -1 && charray[3]== +1 &&
		 pidarray[0]== 3 && pidarray[1]== 2 && pidarray[2]== 2 && pidarray[3]== 3 ){
	       
	         mrecKstarbar = mreck1pi2;
	         mrecKstar    = mrecpi3k4;

	         histosTH1F["hm4rec2OSkstar"]->Fill(mrecKstar);
	         histosTH1F["hm4rec2OSkstarbar"]->Fill(mrecKstarbar);

	         //histosTH1F["hm4rec2OSk1pi2b"]->Fill(mreck1pi2);
		 //histosTH1F["hm4rec2OSpi3k4"]->Fill(mrecpi3k4);
	         //histosTH1F["hm4rec2OSk1p2bp3k4"]->Fill(mreck1pi2pi3k4);
		 
		 //histosTH2F["h2dim2OSk1p2bp3k4"]->Fill(mreck1pi2,mrecpi3k4);
		 histosTH2F["h2dim2OSKK"]->Fill(mrecKstar,mrecKstarbar);
	         
		 //...K*Kbar* study 0.8 - 1.0 GeV at Kbar*
		 //
		 //...| M(Kbar*) - M(895) | < 50.0 MeV ==> 100.0 MeV window ...study
		 if( ( TMath::Abs(mrecKstarbar - m_kstar) < 0.050 ) ){
		   histosTH1F["hm4rec2OSKK100"]->Fill(mrecKstar);
		 }
		 //
		 //...| M(K*Kbar*) - M(895) | < 14.6 MeV ==> mass window: 2*sigma = 29.1 MeV
		 if( ( TMath::Abs(mrecKstar - m_kstar) < 0.0146 ) &&
		     ( TMath::Abs(mrecKstarbar - m_kstar) < 0.0146 ) ){
		   histosTH1F["hm4rec2OSKKwin29"]->Fill(mreck1pi2pi3k4);
		 }
		 //...| M(K*Kbar*) - M(895) | < 29.1 MeV ==> mass window: 4*sigma = 58.3 MeV
		 if( ( TMath::Abs(mrecKstar - m_kstar) < 0.0291 ) &&
		     ( TMath::Abs(mrecKstarbar - m_kstar) < 0.0291 ) ){
		   histosTH1F["hm4rec2OSKKwin58"]->Fill(mreck1pi2pi3k4);
		 }
	     }
    

	     //pi1k2pi3k4
	     
  	     if( charray[0]== +1 && charray[1]== -1 && charray[2]== -1 && charray[3]== +1 &&
		 pidarray[0]== 2 && pidarray[1]== 3 && pidarray[2]== 2 && pidarray[3]== 3 ){
	       
	         mrecKstarbar = mrecpi1k2;
	         mrecKstar    = mrecpi3k4;

	         histosTH1F["hm4rec2OSkstar"]->Fill(mrecKstar);
	         histosTH1F["hm4rec2OSkstarbar"]->Fill(mrecKstarbar);

	         //histosTH1F["hm4rec2OSpi1k2b"]->Fill(mrecpi1k2);
		 //histosTH1F["hm4rec2OSpi3k4"]->Fill(mrecpi3k4);
	         //histosTH1F["hm4rec2OSp1k2bp3k4"]->Fill(mrecpi1k2pi3k4);
		 
		 //histosTH2F["h2dim2OSp1k2bp3k4"]->Fill(mrecpi1k2,mrecpi3k4);
		 histosTH2F["h2dim2OSKK"]->Fill(mrecKstar,mrecKstarbar);
	         
		 //...K*Kbar* study 0.8 - 1.0 GeV at Kbar*
		 //
		 //...| M(Kbar*) - M(895) | < 50.0 MeV ==> 100.0 MeV window ...study
		 if( ( TMath::Abs(mrecKstarbar - m_kstar) < 0.050 ) ){
		   histosTH1F["hm4rec2OSKK100"]->Fill(mrecKstar);
		 }
		 //
		 //...| M(K*Kbar*) - M(895) | < 14.6 MeV ==> mass window: 2*sigma = 29.1 MeV
		 if( ( TMath::Abs(mrecKstar - m_kstar) < 0.0146 ) &&
		     ( TMath::Abs(mrecKstarbar - m_kstar) < 0.0146 ) ){
		   histosTH1F["hm4rec2OSKKwin29"]->Fill(mrecpi1k2pi3k4);
		 }
		 //...| M(K*Kbar*) - M(895) | < 29.1 MeV ==> mass window: 4*sigma = 58.3 MeV
		 if( ( TMath::Abs(mrecKstar - m_kstar) < 0.0291 ) &&
		     ( TMath::Abs(mrecKstarbar - m_kstar) < 0.0291 ) ){
		   histosTH1F["hm4rec2OSKKwin58"]->Fill(mrecpi1k2pi3k4);
		 }
	     }
    
  	     if( charray[0]== -1 && charray[1]== +1 && charray[2]== +1 && charray[3]== -1 &&
		 pidarray[0]== 2 && pidarray[1]== 3 && pidarray[2]== 2 && pidarray[3]== 3 ){

	         mrecKstar    = mrecpi1k2;
	         mrecKstarbar = mrecpi3k4;

	         histosTH1F["hm4rec2OSkstar"]->Fill(mrecKstar);
	         histosTH1F["hm4rec2OSkstarbar"]->Fill(mrecKstarbar);

	         //histosTH1F["hm4rec2OSpi1k2"]->Fill(mrecpi1k2);
		 //histosTH1F["hm4rec2OSpi3k4b"]->Fill(mrecpi3k4);
	         //histosTH1F["hm4rec2OSp1k2p3k4b"]->Fill(mrecpi1k2pi3k4);
		 
		 //histosTH2F["h2dim2OSp1k2p3k4b"]->Fill(mrecpi1k2,mrecpi3k4);
		 histosTH2F["h2dim2OSKK"]->Fill(mrecKstar,mrecKstarbar);
	         
		 //...K*Kbar* study 0.8 - 1.0 GeV at Kbar*
		 //
		 //...| M(Kbar*) - M(895) | < 50.0 MeV ==> 100.0 MeV window ...study
		 if( ( TMath::Abs(mrecKstarbar - m_kstar) < 0.050 ) ){
		   histosTH1F["hm4rec2OSKK100"]->Fill(mrecKstar);
		 }
		 //
		 //...| M(K*Kbar*) - M(895) | < 14.6 MeV ==> mass window: 2*sigma = 29.1 MeV
		 if( ( TMath::Abs(mrecKstar - m_kstar) < 0.0146 ) &&
		     ( TMath::Abs(mrecKstarbar - m_kstar) < 0.0146 ) ){
		   histosTH1F["hm4rec2OSKKwin29"]->Fill(mrecpi1k2pi3k4);
		 }
		 //...| M(K*Kbar*) - M(895) | < 29.1 MeV ==> mass window: 4*sigma = 58.3 MeV
		 if( ( TMath::Abs(mrecKstar - m_kstar) < 0.0291 ) &&
		     ( TMath::Abs(mrecKstarbar - m_kstar) < 0.0291 ) ){
		   histosTH1F["hm4rec2OSKKwin58"]->Fill(mrecpi1k2pi3k4);
		 }
	     }
    
    
	     //pi1k2k3pi4
	     
  	     if( charray[0]== +1 && charray[1]== -1 && charray[2]== +1 && charray[3]== -1 &&
		 pidarray[0]== 2 && pidarray[1]== 3 && pidarray[2]== 3 && pidarray[3]== 2 ){
	       
	         mrecKstarbar = mrecpi1k2;
	         mrecKstar    = mreck3pi4;

	         histosTH1F["hm4rec2OSkstar"]->Fill(mrecKstar);
	         histosTH1F["hm4rec2OSkstarbar"]->Fill(mrecKstarbar);

	         //histosTH1F["hm4rec2OSpi1k2b"]->Fill(mrecpi1k2);
		 //histosTH1F["hm4rec2OSk3pi4"]->Fill(mreck3pi4);
	         //histosTH1F["hm4rec2OSp1k2bk3p4"]->Fill(mrecpi1k2k3pi4);
		 
		 //histosTH2F["h2dim2OSp1k2bk3p4"]->Fill(mrecpi1k2,mreck3pi4);
		 histosTH2F["h2dim2OSKK"]->Fill(mrecKstar,mrecKstarbar);
	         
		 //...K*Kbar* study 0.8 - 1.0 GeV at Kbar*
		 //
		 //...| M(Kbar*) - M(895) | < 50.0 MeV ==> 100.0 MeV window ...study
		 if( ( TMath::Abs(mrecKstarbar - m_kstar) < 0.050 ) ){
		   histosTH1F["hm4rec2OSKK100"]->Fill(mrecKstar);
		 }
		 //
 		 //...| M(K*Kbar*) - M(895) | < 14.6 MeV ==> mass window: 2*sigma = 29.1 MeV
		 if( ( TMath::Abs(mrecKstar - m_kstar) < 0.0146 ) &&
		     ( TMath::Abs(mrecKstarbar - m_kstar) < 0.0146 ) ){
		   histosTH1F["hm4rec2OSKKwin29"]->Fill(mrecpi1k2k3pi4);
		 }
 		 //...| M(K*Kbar*) - M(895) | < 29.1 MeV ==> mass window: 4*sigma = 58.3 MeV
		 if( ( TMath::Abs(mrecKstar - m_kstar) < 0.0291 ) &&
		     ( TMath::Abs(mrecKstarbar - m_kstar) < 0.0291 ) ){
		   histosTH1F["hm4rec2OSKKwin58"]->Fill(mrecpi1k2k3pi4);
		 }
	     }
    
 	     if( charray[0]== -1 && charray[1]== +1 && charray[2]== -1 && charray[3]== +1 &&
		 pidarray[0]== 2 && pidarray[1]== 3 && pidarray[2]== 3 && pidarray[3]== 2 ){

	         mrecKstar    = mrecpi1k2;
	         mrecKstarbar = mreck3pi4;

	         histosTH1F["hm4rec2OSkstar"]->Fill(mrecKstar);
	         histosTH1F["hm4rec2OSkstarbar"]->Fill(mrecKstarbar);

	         //histosTH1F["hm4rec2OSpi1k2"]->Fill(mrecpi1k2);
		 //histosTH1F["hm4rec2OSk3pi4b"]->Fill(mreck3pi4);
	         //histosTH1F["hm4rec2OSp1k2k3p4b"]->Fill(mrecpi1k2k3pi4);
		 
		 //histosTH2F["h2dim2OSp1k2k3p4b"]->Fill(mrecpi1k2,mreck3pi4);
		 histosTH2F["h2dim2OSKK"]->Fill(mrecKstar,mrecKstarbar);
	         
		 //...K*Kbar* study 0.8 - 1.0 GeV at Kbar*
		 //
		 //...| M(Kbar*) - M(895) | < 50.0 MeV ==> 100.0 MeV window ...study
		 if( ( TMath::Abs(mrecKstarbar - m_kstar) < 0.050 ) ){
		   histosTH1F["hm4rec2OSKK100"]->Fill(mrecKstar);
		 }
		 //
		 //...| M(K*Kbar*) - M(895) | < 14.6 MeV ==> mass window: 2*sigma = 29.1 MeV
		 if( ( TMath::Abs(mrecKstar - m_kstar) < 0.0146 ) &&
		     ( TMath::Abs(mrecKstarbar - m_kstar) < 0.0146 ) ){
		   histosTH1F["hm4rec2OSKKwin29"]->Fill(mrecpi1k2k3pi4);
		 }
		 //...| M(K*Kbar*) - M(895) | < 29.1 MeV ==> mass window: 4*sigma = 58.3 MeV
		 if( ( TMath::Abs(mrecKstar - m_kstar) < 0.0291 ) &&
		     ( TMath::Abs(mrecKstarbar - m_kstar) < 0.0291 ) ){
		   histosTH1F["hm4rec2OSKKwin58"]->Fill(mrecpi1k2k3pi4);
		 }
	     }
    

//...1324
	     //k1pi3k2pi4
	     
  	     if( charray[0]== +1 && charray[1]== -1 && charray[2]== -1 && charray[3]== +1 &&
		 pidarray[0]== 3 && pidarray[1]== 2 && pidarray[2]== 3 && pidarray[3]== 2 ){
	       
	         mrecKstar    = mreck1pi3;
	         mrecKstarbar = mreck2pi4;

	         histosTH1F["hm4rec2OSkstar"]->Fill(mrecKstar);
	         histosTH1F["hm4rec2OSkstarbar"]->Fill(mrecKstarbar);

	         //histosTH1F["hm4rec2OSk1pi3"]->Fill(mreck1pi3);
		 //histosTH1F["hm4rec2OSk2pi4b"]->Fill(mreck2pi4);
	         //histosTH1F["hm4rec2OSk1p3k2p4b"]->Fill(mreck1pi3k2pi4);
		 
		 //histosTH2F["h2dim2OSk1p3k2p4b"]->Fill(mreck1pi3,mreck2pi4);
		 histosTH2F["h2dim2OSKK"]->Fill(mrecKstar,mrecKstarbar);
	         
		 //...K*Kbar* study 0.8 - 1.0 GeV at Kbar*
		 //
		 //...| M(Kbar*) - M(895) | < 50.0 MeV ==> 100.0 MeV window ...study
		 if( ( TMath::Abs(mrecKstarbar - m_kstar) < 0.050 ) ){
		   histosTH1F["hm4rec2OSKK100"]->Fill(mrecKstar);
		 }
		 //
		 //...| M(K*Kbar*) - M(895) | < 14.6 MeV ==> mass window: 2*sigma = 29.1 MeV
		 if( ( TMath::Abs(mrecKstar - m_kstar) < 0.0146 ) &&
		     ( TMath::Abs(mrecKstarbar - m_kstar) < 0.0146 ) ){
		   histosTH1F["hm4rec2OSKKwin29"]->Fill(mreck1pi3k2pi4);
		 }
		 //...| M(K*Kbar*) - M(895) | < 29.1 MeV ==> mass window: 4*sigma = 58.3 MeV
		 if( ( TMath::Abs(mrecKstar - m_kstar) < 0.0291 ) &&
		     ( TMath::Abs(mrecKstarbar - m_kstar) < 0.0291 ) ){
		   histosTH1F["hm4rec2OSKKwin58"]->Fill(mreck1pi3k2pi4);
		 }
	     }
    
 	     if( charray[0]== -1 && charray[1]== +1 && charray[2]== +1 && charray[3]== -1 &&
		 pidarray[0]== 3 && pidarray[1]== 2 && pidarray[2]== 3 && pidarray[3]== 2 ){
	       
	         mrecKstarbar = mreck1pi3;
	         mrecKstar    = mreck2pi4;

	         histosTH1F["hm4rec2OSkstar"]->Fill(mrecKstar);
	         histosTH1F["hm4rec2OSkstarbar"]->Fill(mrecKstarbar);

	         //histosTH1F["hm4rec2OSk1pi3b"]->Fill(mreck1pi3);
		 //histosTH1F["hm4rec2OSk2pi4"]->Fill(mreck2pi4);
	         //histosTH1F["hm4rec2OSk1p3bk2p4"]->Fill(mreck1pi3k2pi4);
		 
		 //histosTH2F["h2dim2OSk1p3bk2p4"]->Fill(mreck1pi3,mreck2pi4);
		 histosTH2F["h2dim2OSKK"]->Fill(mrecKstar,mrecKstarbar);
	         
		 //...K*Kbar* study 0.8 - 1.0 GeV at Kbar*
		 //
		 //...| M(Kbar*) - M(895) | < 50.0 MeV ==> 100.0 MeV window ...study
		 if( ( TMath::Abs(mrecKstarbar - m_kstar) < 0.050 ) ){
		   histosTH1F["hm4rec2OSKK100"]->Fill(mrecKstar);
		 }
		 //
		 //...| M(K*Kbar*) - M(895) | < 14.6 MeV ==> mass window: 2*sigma = 29.1 MeV
		 if( ( TMath::Abs(mrecKstar - m_kstar) < 0.0146 ) &&
		     ( TMath::Abs(mrecKstarbar - m_kstar) < 0.0146 ) ){
		   histosTH1F["hm4rec2OSKKwin29"]->Fill(mreck1pi3k2pi4);
		 }
		 //...| M(K*Kbar*) - M(895) | < 29.1 MeV ==> mass window: 4*sigma = 58.3 MeV
		 if( ( TMath::Abs(mrecKstar - m_kstar) < 0.0291 ) &&
		     ( TMath::Abs(mrecKstarbar - m_kstar) < 0.0291 ) ){
		   histosTH1F["hm4rec2OSKKwin58"]->Fill(mreck1pi3k2pi4);
		 }
	     }
    
	     //k1pi3pi2k4
	     
  	     if( charray[0]== +1 && charray[1]== -1 && charray[2]== +1 && charray[3]== -1 &&
		 pidarray[0]== 3 && pidarray[1]== 2 && pidarray[2]== 2 && pidarray[3]== 3 ){
	       
	         mrecKstar    = mreck1pi3;
	         mrecKstarbar = mrecpi2k4;

	         histosTH1F["hm4rec2OSkstar"]->Fill(mrecKstar);
	         histosTH1F["hm4rec2OSkstarbar"]->Fill(mrecKstarbar);

	         //histosTH1F["hm4rec2OSk1pi3"]->Fill(mreck1pi3);
		 //histosTH1F["hm4rec2OSpi2k4b"]->Fill(mrecpi2k4);
	         //histosTH1F["hm4rec2OSk1p3p2k4b"]->Fill(mreck1pi3pi2k4);
		 
		 //histosTH2F["h2dim2OSk1p3p2k4b"]->Fill(mreck1pi3,mrecpi2k4);
		 histosTH2F["h2dim2OSKK"]->Fill(mrecKstar,mrecKstarbar);
	         
		 //...K*Kbar* study 0.8 - 1.0 GeV at Kbar*
		 //
		 //...| M(Kbar*) - M(895) | < 50.0 MeV ==> 100.0 MeV window ...study
		 if( ( TMath::Abs(mrecKstarbar - m_kstar) < 0.050 ) ){
		   histosTH1F["hm4rec2OSKK100"]->Fill(mrecKstar);
		 }
		 //
		 //...| M(K*Kbar*) - M(895) | < 14.6 MeV ==> mass window: 2*sigma = 29.1 MeV
		 if( ( TMath::Abs(mrecKstar - m_kstar) < 0.0146 ) &&
		     ( TMath::Abs(mrecKstarbar - m_kstar) < 0.0146 ) ){
		   histosTH1F["hm4rec2OSKKwin29"]->Fill(mreck1pi3pi2k4);
		 }
		 //...| M(K*Kbar*) - M(895) | < 29.1 MeV ==> mass window: 4*sigma = 58.3 MeV
		 if( ( TMath::Abs(mrecKstar - m_kstar) < 0.0291 ) &&
		     ( TMath::Abs(mrecKstarbar - m_kstar) < 0.0291 ) ){
		   histosTH1F["hm4rec2OSKKwin58"]->Fill(mreck1pi3pi2k4);
		 }
	     }
    
  	     if( charray[0]== -1 && charray[1]== +1 && charray[2]== -1 && charray[3]== +1 &&
		 pidarray[0]== 3 && pidarray[1]== 2 && pidarray[2]== 2 && pidarray[3]== 3 ){
	       
	         mrecKstarbar = mreck1pi3;
	         mrecKstar    = mrecpi2k4;

	         histosTH1F["hm4rec2OSkstar"]->Fill(mrecKstar);
	         histosTH1F["hm4rec2OSkstarbar"]->Fill(mrecKstarbar);

	         //histosTH1F["hm4rec2OSk1pi3b"]->Fill(mreck1pi3);
		 //histosTH1F["hm4rec2OSpi2k4"]->Fill(mrecpi2k4);
	         //histosTH1F["hm4rec2OSk1p3bpi2k4"]->Fill(mreck1pi3pi2k4);
		 
		 //histosTH2F["h2dim2OSk1p3bp2k4"]->Fill(mreck1pi3,mrecpi2k4);
		 histosTH2F["h2dim2OSKK"]->Fill(mrecKstar,mrecKstarbar);
	         
		 //...K*Kbar* study 0.8 - 1.0 GeV at Kbar*
		 //
		 //...| M(Kbar*) - M(895) | < 50.0 MeV ==> 100.0 MeV window ...study
		 if( ( TMath::Abs(mrecKstarbar - m_kstar) < 0.050 ) ){
		   histosTH1F["hm4rec2OSKK100"]->Fill(mrecKstar);
		 }
		 //
		 //...| M(K*Kbar*) - M(895) | < 14.6 MeV ==> mass window: 2*sigma = 29.1 MeV
		 if( ( TMath::Abs(mrecKstar - m_kstar) < 0.0146 ) &&
		     ( TMath::Abs(mrecKstarbar - m_kstar) < 0.0146 ) ){
		   histosTH1F["hm4rec2OSKKwin29"]->Fill(mreck1pi3pi2k4);
		 }
		 //...| M(K*Kbar*) - M(895) | < 29.1 MeV ==> mass window: 4*sigma = 58.3 MeV
		 if( ( TMath::Abs(mrecKstar - m_kstar) < 0.0291 ) &&
		     ( TMath::Abs(mrecKstarbar - m_kstar) < 0.0291 ) ){
		   histosTH1F["hm4rec2OSKKwin58"]->Fill(mreck1pi3pi2k4);
		 }
	     }
    
	     //pi1k3k2pi4
	     
  	     if( charray[0]== +1 && charray[1]== -1 && charray[2]== +1 && charray[3]== -1 &&
		 pidarray[0]== 2 && pidarray[1]== 3 && pidarray[2]== 3 && pidarray[3]== 2 ){
	       
	         mrecKstarbar = mrecpi1k3;
	         mrecKstar    = mreck2pi4;

	         histosTH1F["hm4rec2OSkstar"]->Fill(mrecKstar);
	         histosTH1F["hm4rec2OSkstarbar"]->Fill(mrecKstarbar);

	         //histosTH1F["hm4rec2OSpi1k3b"]->Fill(mrecpi1k3);
		 //histosTH1F["hm4rec2OSk2pi4"]->Fill(mreck2pi4);
	         //histosTH1F["hm4rec2OSp1k3bk2p4"]->Fill(mrecpi1k3k2pi4);
		 
		 //histosTH2F["h2dim2OSp1k3bk2p4"]->Fill(mrecpi1k3,mreck2pi4);
		 histosTH2F["h2dim2OSKK"]->Fill(mrecKstar,mrecKstarbar);
	         
		 //...K*Kbar* study 0.8 - 1.0 GeV at Kbar*
		 //
		 //...| M(Kbar*) - M(895) | < 50.0 MeV ==> 100.0 MeV window ...study
		 if( ( TMath::Abs(mrecKstarbar - m_kstar) < 0.050 ) ){
		   histosTH1F["hm4rec2OSKK100"]->Fill(mrecKstar);
		 }
		 //
		 //...| M(K*Kbar*) - M(895) | < 14.6 MeV ==> mass window: 2*sigma = 29.1 MeV
		 if( ( TMath::Abs(mrecKstar - m_kstar) < 0.0146 ) &&
		     ( TMath::Abs(mrecKstarbar - m_kstar) < 0.0146 ) ){
		   histosTH1F["hm4rec2OSKKwin29"]->Fill(mrecpi1k3k2pi4);
		 }
		 //...| M(K*Kbar*) - M(895) | < 29.1 MeV ==> mass window: 4*sigma = 58.3 MeV
		 if( ( TMath::Abs(mrecKstar - m_kstar) < 0.0291 ) &&
		     ( TMath::Abs(mrecKstarbar - m_kstar) < 0.0291 ) ){
		   histosTH1F["hm4rec2OSKKwin58"]->Fill(mrecpi1k3k2pi4);
		 }
	     }
    
  	     if( charray[0]== -1 && charray[1]== +1 && charray[2]== -1 && charray[3]== +1 &&
		 pidarray[0]== 2 && pidarray[1]== 3 && pidarray[2]== 3 && pidarray[3]== 2 ){
	       
	         mrecKstar    = mrecpi1k3;
	         mrecKstarbar = mreck2pi4;

	         histosTH1F["hm4rec2OSkstar"]->Fill(mrecKstar);
	         histosTH1F["hm4rec2OSkstarbar"]->Fill(mrecKstarbar);

	         //histosTH1F["hm4rec2OSpi1k3"]->Fill(mrecpi1k3);
		 //histosTH1F["hm4rec2OSk2pi4b"]->Fill(mreck2pi4);
	         //histosTH1F["hm4rec2OSp1k3k2p4b"]->Fill(mrecpi1k3k2pi4);
		 
		 //histosTH2F["h2dim2OSp1k3k2p4b"]->Fill(mrecpi1k3,mreck2pi4);
		 histosTH2F["h2dim2OSKK"]->Fill(mrecKstar,mrecKstarbar);
	         
		 //...K*Kbar* study 0.8 - 1.0 GeV at Kbar*
		 //
		 //...| M(Kbar*) - M(895) | < 50.0 MeV ==> 100.0 MeV window ...study
		 if( ( TMath::Abs(mrecKstarbar - m_kstar) < 0.050 ) ){
		   histosTH1F["hm4rec2OSKK100"]->Fill(mrecKstar);
		 }
		 //
		 //...| M(K*Kbar*) - M(895) | < 14.6 MeV ==> mass window: 2*sigma = 29.1 MeV
		 if( ( TMath::Abs(mrecKstar - m_kstar) < 0.0146 ) &&
		     ( TMath::Abs(mrecKstarbar - m_kstar) < 0.0146 ) ){
		   histosTH1F["hm4rec2OSKKwin29"]->Fill(mrecpi1k3k2pi4);
		 }
		 //...| M(K*Kbar*) - M(895) | < 29.1 MeV ==> mass window: 4*sigma = 58.3 MeV
		 if( ( TMath::Abs(mrecKstar - m_kstar) < 0.0291 ) &&
		     ( TMath::Abs(mrecKstarbar - m_kstar) < 0.0291 ) ){
		   histosTH1F["hm4rec2OSKKwin58"]->Fill(mrecpi1k3k2pi4);
		 }
	     }
    
	     //pi1k3pi2k4
	     
  	     if( charray[0]== +1 && charray[1]== -1 && charray[2]== -1 && charray[3]== +1 &&
		 pidarray[0]== 2 && pidarray[1]== 3 && pidarray[2]== 2 && pidarray[3]== 3 ){
	       
	         mrecKstarbar = mrecpi1k3;
	         mrecKstar    = mrecpi2k4;

	         histosTH1F["hm4rec2OSkstar"]->Fill(mrecKstar);
	         histosTH1F["hm4rec2OSkstarbar"]->Fill(mrecKstarbar);

	         //histosTH1F["hm4rec2OSpi1k3b"]->Fill(mrecpi1k3);
		 //histosTH1F["hm4rec2OSpi2k4"]->Fill(mrecpi2k4);
	         //histosTH1F["hm4rec2OSp1k3bp2k4"]->Fill(mrecpi1k3pi2k4);
		 
		 //histosTH2F["h2dim2OSp1k3bp2k4"]->Fill(mrecpi1k3,mrecpi2k4);
		 histosTH2F["h2dim2OSKK"]->Fill(mrecKstar,mrecKstarbar);
	         
		 //...K*Kbar* study 0.8 - 1.0 GeV at Kbar*
		 //
		 //...| M(Kbar*) - M(895) | < 50.0 MeV ==> 100.0 MeV window ...study
		 if( ( TMath::Abs(mrecKstarbar - m_kstar) < 0.050 ) ){
		   histosTH1F["hm4rec2OSKK100"]->Fill(mrecKstar);
		 }
		 //
		 //...| M(K*Kbar*) - M(895) | < 14.6 MeV ==> mass window: 2*sigma = 29.1 MeV
		 if( ( TMath::Abs(mrecKstar - m_kstar) < 0.0146 ) &&
		     ( TMath::Abs(mrecKstarbar - m_kstar) < 0.0146 ) ){
		   histosTH1F["hm4rec2OSKKwin29"]->Fill(mrecpi1k3pi2k4);
		 }
		 //...| M(K*Kbar*) - M(895) | < 29.1 MeV ==> mass window: 4*sigma = 58.3 MeV
		 if( ( TMath::Abs(mrecKstar - m_kstar) < 0.0291 ) &&
		     ( TMath::Abs(mrecKstarbar - m_kstar) < 0.0291 ) ){
		   histosTH1F["hm4rec2OSKKwin58"]->Fill(mrecpi1k3pi2k4);
		 }
	     }
    
  	     if( charray[0]== -1 && charray[1]== +1 && charray[2]== +1 && charray[3]== -1 &&
		 pidarray[0]== 2 && pidarray[1]== 3 && pidarray[2]== 2 && pidarray[3]== 3 ){
	       
	         mrecKstar    = mrecpi1k3;
	         mrecKstarbar = mrecpi2k4;

	         histosTH1F["hm4rec2OSkstar"]->Fill(mrecKstar);
	         histosTH1F["hm4rec2OSkstarbar"]->Fill(mrecKstarbar);

	         //histosTH1F["hm4rec2OSpi1k3"]->Fill(mrecpi1k3);
		 //histosTH1F["hm4rec2OSpi2k4b"]->Fill(mrecpi2k4);
	         //histosTH1F["hm4rec2OSp1k3p2k4b"]->Fill(mrecpi1k3pi2k4);
		 
		 //histosTH2F["h2dim2OSp1k3p2k4b"]->Fill(mrecpi1k3,mrecpi2k4);
		 histosTH2F["h2dim2OSKK"]->Fill(mrecKstar,mrecKstarbar);
	         
		 //...K*Kbar* study 0.8 - 1.0 GeV at Kbar*
		 //
		 //...| M(Kbar*) - M(895) | < 50.0 MeV ==> 100.0 MeV window ...study
		 if( ( TMath::Abs(mrecKstarbar - m_kstar) < 0.050 ) ){
		   histosTH1F["hm4rec2OSKK100"]->Fill(mrecKstar);
		 }
		 //
		 //...| M(K*Kbar*) - M(895) | < 14.6 MeV ==> mass window: 2*sigma = 29.1 MeV
		 if( ( TMath::Abs(mrecKstar - m_kstar) < 0.0146 ) &&
		     ( TMath::Abs(mrecKstarbar - m_kstar) < 0.0146 ) ){
		   histosTH1F["hm4rec2OSKKwin29"]->Fill(mrecpi1k3pi2k4);
		 }
		 //...| M(K*Kbar*) - M(895) | < 29.1 MeV ==> mass window: 4*sigma = 58.3 MeV
		 if( ( TMath::Abs(mrecKstar - m_kstar) < 0.0291 ) &&
		     ( TMath::Abs(mrecKstarbar - m_kstar) < 0.0291 ) ){
		   histosTH1F["hm4rec2OSKKwin58"]->Fill(mrecpi1k3pi2k4);
		 }
	     }

//...1423
	     //k1pi4k2pi3
	     
  	     if( charray[0]== +1 && charray[1]== -1 && charray[2]== -1 && charray[3]== +1 &&
		 pidarray[0]== 3 && pidarray[1]== 2 && pidarray[2]== 3 && pidarray[3]== 2 ){
	       
	         mrecKstar    = mreck1pi4;
	         mrecKstarbar = mreck2pi3;

	         histosTH1F["hm4rec2OSkstar"]->Fill(mrecKstar);
	         histosTH1F["hm4rec2OSkstarbar"]->Fill(mrecKstarbar);

	         //histosTH1F["hm4rec2OSk1pi4"]->Fill(mreck1pi4);
		 //histosTH1F["hm4rec2OSk2pi3b"]->Fill(mreck2pi3);
	         //histosTH1F["hm4rec2OSk1p4k2p3b"]->Fill(mreck1pi4k2pi3);
		 
		 //histosTH2F["h2dim2OSk1p4k2p3b"]->Fill(mreck1pi4,mreck2pi3);
		 histosTH2F["h2dim2OSKK"]->Fill(mrecKstar,mrecKstarbar);
	         
		 //...K*Kbar* study 0.8 - 1.0 GeV at Kbar*
		 //
		 //...| M(Kbar*) - M(895) | < 50.0 MeV ==> 100.0 MeV window ...study
		 if( ( TMath::Abs(mrecKstarbar - m_kstar) < 0.050 ) ){
		   histosTH1F["hm4rec2OSKK100"]->Fill(mrecKstar);
		 }
		 //
		 //...| M(K*Kbar*) - M(895) | < 14.6 MeV ==> mass window: 2*sigma = 29.1 MeV
		 if( ( TMath::Abs(mrecKstar - m_kstar) < 0.0146 ) &&
		     ( TMath::Abs(mrecKstarbar - m_kstar) < 0.0146 ) ){
		   histosTH1F["hm4rec2OSKKwin29"]->Fill(mreck1pi4k2pi3);
		 }
		 //...| M(K*Kbar*) - M(895) | < 29.1 MeV ==> mass window: 4*sigma = 58.3 MeV
		 if( ( TMath::Abs(mrecKstar - m_kstar) < 0.0291 ) &&
		     ( TMath::Abs(mrecKstarbar - m_kstar) < 0.0291 ) ){
		   histosTH1F["hm4rec2OSKKwin58"]->Fill(mreck1pi4k2pi3);
		 }
	     }
    
 	     if( charray[0]== -1 && charray[1]== +1 && charray[2]== +1 && charray[3]== -1 &&
		 pidarray[0]== 3 && pidarray[1]== 2 && pidarray[2]== 3 && pidarray[3]== 2 ){
	       
	         mrecKstarbar = mreck1pi4;
	         mrecKstar    = mreck2pi3;

	         histosTH1F["hm4rec2OSkstar"]->Fill(mrecKstar);
	         histosTH1F["hm4rec2OSkstarbar"]->Fill(mrecKstarbar);

	         //histosTH1F["hm4rec2OSk1pi4b"]->Fill(mreck1pi4);
		 //histosTH1F["hm4rec2OSk2pi3"]->Fill(mreck2pi3);
	         //histosTH1F["hm4rec2OSk1p4bk2p3"]->Fill(mreck1pi4k2pi3);
		 
		 //histosTH2F["h2dim2OSk1p4bk2p3"]->Fill(mreck1pi4,mreck2pi3);
		 histosTH2F["h2dim2OSKK"]->Fill(mrecKstar,mrecKstarbar);
	         
		 //...K*Kbar* study 0.8 - 1.0 GeV at Kbar*
		 //
		 //...| M(Kbar*) - M(895) | < 50.0 MeV ==> 100.0 MeV window ...study
		 if( ( TMath::Abs(mrecKstarbar - m_kstar) < 0.050 ) ){
		   histosTH1F["hm4rec2OSKK100"]->Fill(mrecKstar);
		 }
		 //
		 //...| M(K*Kbar*) - M(895) | < 14.6 MeV ==> mass window: 2*sigma = 29.1 MeV
		 if( ( TMath::Abs(mrecKstar - m_kstar) < 0.0146 ) &&
		     ( TMath::Abs(mrecKstarbar - m_kstar) < 0.0146 ) ){
		   histosTH1F["hm4rec2OSKKwin29"]->Fill(mreck1pi4k2pi3);
		 }
		 //...| M(K*Kbar*) - M(895) | < 29.1 MeV ==> mass window: 4*sigma = 58.3 MeV
		 if( ( TMath::Abs(mrecKstar - m_kstar) < 0.0291 ) &&
		     ( TMath::Abs(mrecKstarbar - m_kstar) < 0.0291 ) ){
		   histosTH1F["hm4rec2OSKKwin58"]->Fill(mreck1pi4k2pi3);
		 }
	     }
    
	     //k1pi4pi2k3
	     
  	     if( charray[0]== +1 && charray[1]== -1 && charray[2]== +1 && charray[3]== -1 &&
		 pidarray[0]== 3 && pidarray[1]== 2 && pidarray[2]== 2 && pidarray[3]== 3 ){
	       
	         mrecKstar    = mreck1pi4;
	         mrecKstarbar = mrecpi2k3;

	         histosTH1F["hm4rec2OSkstar"]->Fill(mrecKstar);
	         histosTH1F["hm4rec2OSkstarbar"]->Fill(mrecKstarbar);

	         //histosTH1F["hm4rec2OSk1pi4"]->Fill(mreck1pi4);
		 //histosTH1F["hm4rec2OSpi2k3b"]->Fill(mrecpi2k3);
	         //histosTH1F["hm4rec2OSk1p4p2k3b"]->Fill(mreck1pi4pi2k3);
		 
		 //histosTH2F["h2dim2OSk1p4p2k3b"]->Fill(mreck1pi4,mrecpi2k3);
		 histosTH2F["h2dim2OSKK"]->Fill(mrecKstar,mrecKstarbar);
	         
		 //...K*Kbar* study 0.8 - 1.0 GeV at Kbar*
		 //
		 //...| M(Kbar*) - M(895) | < 50.0 MeV ==> 100.0 MeV window ...study
		 if( ( TMath::Abs(mrecKstarbar - m_kstar) < 0.050 ) ){
		   histosTH1F["hm4rec2OSKK100"]->Fill(mrecKstar);
		 }
		 //
		 //...| M(K*Kbar*) - M(895) | < 14.6 MeV ==> mass window: 2*sigma = 29.1 MeV
		 if( ( TMath::Abs(mrecKstar - m_kstar) < 0.0146 ) &&
		     ( TMath::Abs(mrecKstarbar - m_kstar) < 0.0146 ) ){
		   histosTH1F["hm4rec2OSKKwin29"]->Fill(mreck1pi4pi2k3);
		 }
		 //...| M(K*Kbar*) - M(895) | < 29.1 MeV ==> mass window: 4*sigma = 58.3 MeV
		 if( ( TMath::Abs(mrecKstar - m_kstar) < 0.0291 ) &&
		     ( TMath::Abs(mrecKstarbar - m_kstar) < 0.0291 ) ){
		   histosTH1F["hm4rec2OSKKwin58"]->Fill(mreck1pi4pi2k3);
		 }
	     }
    
  	     if( charray[0]== -1 && charray[1]== +1 && charray[2]== -1 && charray[3]== +1 &&
		 pidarray[0]== 3 && pidarray[1]== 2 && pidarray[2]== 2 && pidarray[3]== 3 ){
	       
	         mrecKstarbar = mreck1pi4;
	         mrecKstar    = mrecpi2k3;

	         histosTH1F["hm4rec2OSkstar"]->Fill(mrecKstar);
	         histosTH1F["hm4rec2OSkstarbar"]->Fill(mrecKstarbar);

	         //histosTH1F["hm4rec2OSk1pi4b"]->Fill(mreck1pi4);
		 //histosTH1F["hm4rec2OSpi2k3"]->Fill(mrecpi2k3);
	         //histosTH1F["hm4rec2OSk1p4bp2k3"]->Fill(mreck1pi4pi2k3);
		 
		 //histosTH2F["h2dim2OSk1p4bp2k3"]->Fill(mreck1pi4,mrecpi2k3);
		 histosTH2F["h2dim2OSKK"]->Fill(mrecKstar,mrecKstarbar);
	         
		 //...K*Kbar* study 0.8 - 1.0 GeV at Kbar*
		 //
		 //...| M(Kbar*) - M(895) | < 50.0 MeV ==> 100.0 MeV window ...study
		 if( ( TMath::Abs(mrecKstarbar - m_kstar) < 0.050 ) ){
		   histosTH1F["hm4rec2OSKK100"]->Fill(mrecKstar);
		 }
		 //
		 //...| M(K*Kbar*) - M(895) | < 14.6 MeV ==> mass window: 2*sigma = 29.1 MeV
		 if( ( TMath::Abs(mrecKstar - m_kstar) < 0.0146 ) &&
		     ( TMath::Abs(mrecKstarbar - m_kstar) < 0.0146 ) ){
		   histosTH1F["hm4rec2OSKKwin29"]->Fill(mreck1pi4pi2k3);
		 }
		 //...| M(K*Kbar*) - M(895) | < 29.1 MeV ==> mass window: 4*sigma = 58.3 MeV
		 if( ( TMath::Abs(mrecKstar - m_kstar) < 0.0291 ) &&
		     ( TMath::Abs(mrecKstarbar - m_kstar) < 0.0291 ) ){
		   histosTH1F["hm4rec2OSKKwin58"]->Fill(mreck1pi4pi2k3);
		 }
	     }
    
	     //pi1k4k2pi3
	     
  	     if( charray[0]== +1 && charray[1]== -1 && charray[2]== +1 && charray[3]== -1 &&
		 pidarray[0]== 2 && pidarray[1]== 3 && pidarray[2]== 3 && pidarray[3]== 2 ){
	       
	         mrecKstarbar = mrecpi1k4;
	         mrecKstar    = mreck2pi3;

	         histosTH1F["hm4rec2OSkstar"]->Fill(mrecKstar);
	         histosTH1F["hm4rec2OSkstarbar"]->Fill(mrecKstarbar);

	         //histosTH1F["hm4rec2OSpi1k4b"]->Fill(mrecpi1k4);
		 //histosTH1F["hm4rec2OSk2pi3"]->Fill(mreck2pi3);
	         //histosTH1F["hm4rec2OSp1k4bk2p3"]->Fill(mrecpi1k4k2pi3);
		 
		 //histosTH2F["h2dim2OSp1k4bk2p3"]->Fill(mrecpi1k4,mreck2pi3);
		 histosTH2F["h2dim2OSKK"]->Fill(mrecKstar,mrecKstarbar);
	         
		 //...K*Kbar* study 0.8 - 1.0 GeV at Kbar*
		 //
		 //...| M(Kbar*) - M(895) | < 50.0 MeV ==> 100.0 MeV window ...study
		 if( ( TMath::Abs(mrecKstarbar - m_kstar) < 0.050 ) ){
		   histosTH1F["hm4rec2OSKK100"]->Fill(mrecKstar);
		 }
		 //
		 //...| M(K*Kbar*) - M(895) | < 14.6 MeV ==> mass window: 2*sigma = 29.1 MeV
		 if( ( TMath::Abs(mrecKstar - m_kstar) < 0.0146 ) &&
		     ( TMath::Abs(mrecKstarbar - m_kstar) < 0.0146 ) ){
		   histosTH1F["hm4rec2OSKKwin29"]->Fill(mrecpi1k4k2pi3);
		 }
		 //...| M(K*Kbar*) - M(895) | < 29.1 MeV ==> mass window: 4*sigma = 58.3 MeV
		 if( ( TMath::Abs(mrecKstar - m_kstar) < 0.0291 ) &&
		     ( TMath::Abs(mrecKstarbar - m_kstar) < 0.0291 ) ){
		   histosTH1F["hm4rec2OSKKwin58"]->Fill(mrecpi1k4k2pi3);
		 }
	     }
    
  	     if( charray[0]== -1 && charray[1]== +1 && charray[2]== -1 && charray[3]== +1 &&
		 pidarray[0]== 2 && pidarray[1]== 3 && pidarray[2]== 3 && pidarray[3]== 2 ){
	       
	         mrecKstar    = mrecpi1k4;
	         mrecKstarbar = mreck2pi3;

	         histosTH1F["hm4rec2OSkstar"]->Fill(mrecKstar);
	         histosTH1F["hm4rec2OSkstarbar"]->Fill(mrecKstarbar);

	         //histosTH1F["hm4rec2OSpi1k4"]->Fill(mrecpi1k4);
		 //histosTH1F["hm4rec2OSk2pi3b"]->Fill(mreck2pi3);
	         //histosTH1F["hm4rec2OSp1k4k2p3b"]->Fill(mrecpi1k4k2pi3);
		 
		 //histosTH2F["h2dim2OSp1k4k2p3b"]->Fill(mrecpi1k4,mreck2pi3);
		 histosTH2F["h2dim2OSKK"]->Fill(mrecKstar,mrecKstarbar);
	         
		 //...K*Kbar* study 0.8 - 1.0 GeV at Kbar*
		 //
		 //...| M(Kbar*) - M(895) | < 50.0 MeV ==> 100.0 MeV window ...study
		 if( ( TMath::Abs(mrecKstarbar - m_kstar) < 0.050 ) ){
		   histosTH1F["hm4rec2OSKK100"]->Fill(mrecKstar);
		 }
		 //
		 //...| M(K*Kbar*) - M(895) | < 14.6 MeV ==> mass window: 2*sigma = 29.1 MeV
		 if( ( TMath::Abs(mrecKstar - m_kstar) < 0.0146 ) &&
		     ( TMath::Abs(mrecKstarbar - m_kstar) < 0.0146 ) ){
		   histosTH1F["hm4rec2OSKKwin29"]->Fill(mrecpi1k4k2pi3);
		 }
		 //...| M(K*Kbar*) - M(895) | < 29.1 MeV ==> mass window: 4*sigma = 58.3 MeV
		 if( ( TMath::Abs(mrecKstar - m_kstar) < 0.0291 ) &&
		     ( TMath::Abs(mrecKstarbar - m_kstar) < 0.0291 ) ){
		   histosTH1F["hm4rec2OSKKwin58"]->Fill(mrecpi1k4k2pi3);
		 }
	     }
    
	     //pi1k4pi2k3
	     
  	     if( charray[0]== +1 && charray[1]== -1 && charray[2]== -1 && charray[3]== +1 &&
		 pidarray[0]== 2 && pidarray[1]== 3 && pidarray[2]== 2 && pidarray[3]== 3 ){
	       
	         mrecKstarbar = mrecpi1k4;
	         mrecKstar    = mrecpi2k3;

	         histosTH1F["hm4rec2OSkstar"]->Fill(mrecKstar);
	         histosTH1F["hm4rec2OSkstarbar"]->Fill(mrecKstarbar);

	         //histosTH1F["hm4rec2OSpi1k4b"]->Fill(mrecpi1k4);
		 //histosTH1F["hm4rec2OSpi2k3"]->Fill(mrecpi2k3);
	         //histosTH1F["hm4rec2OSp1k4bp2k3"]->Fill(mrecpi1k4pi2k3);
		 
		 //histosTH2F["h2dim2OSp1k4bp2k3"]->Fill(mrecpi1k4,mrecpi2k3);
		 histosTH2F["h2dim2OSKK"]->Fill(mrecKstar,mrecKstarbar);
	         
		 //...K*Kbar* study 0.8 - 1.0 GeV at Kbar*
		 //
		 //...| M(Kbar*) - M(895) | < 50.0 MeV ==> 100.0 MeV window ...study
		 if( ( TMath::Abs(mrecKstarbar - m_kstar) < 0.050 ) ){
		   histosTH1F["hm4rec2OSKK100"]->Fill(mrecKstar);
		 }
		 //
		 //...| M(K*Kbar*) - M(895) | < 14.6 MeV ==> mass window: 2*sigma = 29.1 MeV
		 if( ( TMath::Abs(mrecKstar - m_kstar) < 0.0146 ) &&
		     ( TMath::Abs(mrecKstarbar - m_kstar) < 0.0146 ) ){
		   histosTH1F["hm4rec2OSKKwin29"]->Fill(mrecpi1k4pi2k3);
		 }
		 //...| M(K*Kbar*) - M(895) | < 29.1 MeV ==> mass window: 4*sigma = 58.3 MeV
		 if( ( TMath::Abs(mrecKstar - m_kstar) < 0.0291 ) &&
		     ( TMath::Abs(mrecKstarbar - m_kstar) < 0.0291 ) ){
		   histosTH1F["hm4rec2OSKKwin58"]->Fill(mrecpi1k4pi2k3);
		 }
	     }
    
  	     if( charray[0]== -1 && charray[1]== +1 && charray[2]== +1 && charray[3]== -1 &&
		 pidarray[0]== 2 && pidarray[1]== 3 && pidarray[2]== 2 && pidarray[3]== 3 ){
	       
	         mrecKstarbar = mrecpi1k4;
	         mrecKstar    = mrecpi2k3;

	         histosTH1F["hm4rec2OSkstar"]->Fill(mrecKstar);
	         histosTH1F["hm4rec2OSkstarbar"]->Fill(mrecKstarbar);

	         //histosTH1F["hm4rec2OSpi1k4"]->Fill(mrecpi1k4);
		 //histosTH1F["hm4rec2OSpi2k3b"]->Fill(mrecpi2k3);
	         //histosTH1F["hm4rec2OSp1k4pi2k3b"]->Fill(mrecpi1k4pi2k3);
		 
		 //histosTH2F["h2dim2OSp1k4pi2k3b"]->Fill(mrecpi1k4,mrecpi2k3);
		 histosTH2F["h2dim2OSKK"]->Fill(mrecKstar,mrecKstarbar);
	         
		 //...K*Kbar* study 0.8 - 1.0 GeV at Kbar*
		 //
		 //...| M(Kbar*) - M(895) | < 50.0 MeV ==> 100.0 MeV window ...study
		 if( ( TMath::Abs(mrecKstarbar - m_kstar) < 0.050 ) ){
		   histosTH1F["hm4rec2OSKK100"]->Fill(mrecKstar);
		 }
		 //
		 //...| M(K*Kbar*) - M(895) | < 14.6 MeV ==> mass window: 2*sigma = 29.1 MeV
		 if( ( TMath::Abs(mrecKstar - m_kstar) < 0.0146 ) &&
		     ( TMath::Abs(mrecKstarbar - m_kstar) < 0.0146 ) ){
		   histosTH1F["hm4rec2OSKKwin29"]->Fill(mrecpi1k4pi2k3);
		 }
		 //...| M(K*Kbar*) - M(895) | < 29.1 MeV ==> mass window: 4*sigma = 58.3 MeV
		 if( ( TMath::Abs(mrecKstar - m_kstar) < 0.0291 ) &&
		     ( TMath::Abs(mrecKstarbar - m_kstar) < 0.0291 ) ){
		   histosTH1F["hm4rec2OSKKwin58"]->Fill(mrecpi1k4pi2k3);
		 }
	     }
    
	     /*
	     // K*Kbar*
       	     // Gamma =  MeV
	     // sigma =  MeV
             //
	     //...final study TEST
	     //...| M(K*Kbar*) - M(895) | < 14.6 MeV ==> mass window: 2*sigma = 29.1 MeV
             if( ( TMath::Abs(mrecKstar - m_kstar) < 0.0146 ) &&
		 ( TMath::Abs(mrecKstarbar - m_kstar) < 0.0146 ) ){
	       	histosTH1F["hm4rec2OSKKwin29"]->Fill(mrec);
	     }
	     */
	     
	    } //...end of nvtx=1 nks=0 nlam=0
          } //...end totalcharge ...end cut
	  


	 
  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


       //...cut 1    ...nvtx==1 or 2     ..............????????
       
	 if(totcharge==0){
	   //...Luiz
	   histosTH1F["hm4recOS"]->Fill(mrec);
	   //...extended to 10 GeV
	   histosTH1F["hm4recOS2"]->Fill(mrec);
         }
	  

       //...cut 2    ...Luiz

	 if(totcharge==0){                      // ...4-pions pi+pi-pi+pi-

	     histosTH1F["hm4rec2OS"]->Fill(mrec); 

	     //...dEdx ntrk=4 nvtx=1  ...fiducial Q=0
	     //
	     //////if(nvtx==1){
	     int nn=0;
             for(TrackCollection::const_iterator itTrack = tracks->begin();itTrack != tracks->end();++itTrack) {
	     //	       
	     reco::TrackRef aTrackRef4 = reco::TrackRef(tracks, nn);
             const reco::DeDxData& dedxPIXObj4 = dedxPIXTrack[aTrackRef4];	   
             mydedxPIX4 = dedxPIXObj4.dEdx();
	     //
             histosTH2F["hdedxPIX4"]->Fill( itTrack->p() , mydedxPIX4 );
	     //...only pions
             if(pidarray[0]==2 && pidarray[1]==2 && pidarray[2]==2 && pidarray[3]==2){
             histosTH2F["hdedxPIX4pi"]->Fill( itTrack->p() , mydedxPIX4 );
	     }
	     //...only kaons
	     if(pidarray[0]==3 && pidarray[1]==3 && pidarray[2]==3 && pidarray[3]==3){
             histosTH2F["hdedxPIX4k"]->Fill( itTrack->p() , mydedxPIX4 );
	     }
	     //...only protons
	     if(pidarray[0]==4 && pidarray[1]==4 && pidarray[2]==4 && pidarray[3]==4){
             histosTH2F["hdedxPIX4p"]->Fill( itTrack->p() , mydedxPIX4 );
	     }
	     //...only electrons
	     if(pidarray[0]==1 && pidarray[1]==1 && pidarray[2]==1 && pidarray[3]==1){
             histosTH2F["hdedxPIX4e"]->Fill( itTrack->p() , mydedxPIX4 );
	     }
	     //...only unkown
	     if(pidarray[0]==0 && pidarray[1]==0 && pidarray[2]==0 && pidarray[3]==0){
             histosTH2F["hdedxPIX4u"]->Fill( itTrack->p() , mydedxPIX4 );
	     }
	     //...only Kpi
	     if( (pidarray[0]==2 && pidarray[1]==2 && pidarray[2]==2 && pidarray[3]==3) ||
		 (pidarray[0]==2 && pidarray[1]==2 && pidarray[2]==3 && pidarray[3]==2) ||
		 (pidarray[0]==2 && pidarray[1]==3 && pidarray[2]==2 && pidarray[3]==2) ||
		 (pidarray[0]==3 && pidarray[1]==2 && pidarray[2]==2 && pidarray[3]==2) ){
             histosTH2F["hdedxPIX4kpi"]->Fill( itTrack->p() , mydedxPIX4 );
	     }
	     nn++;
	     } //...end of dE/dx cuts
	     //////} //...end of nvtx=1
	     
	     //...nvtx=1
	     if(nvtx==1){
	        histosTH1F["hm4rec2OS2"]->Fill(mrec);
		//...PID=2 pion
		if(pidarray[0]==2 && pidarray[1]==2 && pidarray[2]==2 && pidarray[3]==2){
	        histosTH1F["hm4rec2OS2pid"]->Fill(mrec);
	        }
		//...pidV0
	        if(pidv0array[0]==2 && pidv0array[1]==2 && pidv0array[2]==2 && pidv0array[3]==2){
	        histosTH1F["hm4rec2OS2pidv0"]->Fill(mrec);
	        }
		//...no V0    ...pure? 4-pion channel
       	        if(!nks){
		   histosTH1F["hm4rec2OS2nov0"]->Fill(mrec);
		 }
	     }
	     
	     //...nvtx=1 and nks=0  : type:10
	     //...does not improve the 4-pi channel, the difference is 22k events only
	     if(nvtx==1 && nks==0){
	        histosTH1F["hm4rec2OS2veeno10"]->Fill(mrec);
	     }

	     //...nvtx=2
	     if(nvtx==2){
	        histosTH1F["hm4rec2OS3"]->Fill(mrec);
	        //...no V0    ...pure 4-pion channel
       	        if(!nks){
		   histosTH1F["hm4rec2OS3nov0"]->Fill(mrec);
		 }
	     }
	     //...nvtx=2 and nks=0  : type:20
	     //...
	     if(nvtx==2 && nks==0){
	        histosTH1F["hm4rec2OS3veeno20"]->Fill(mrec);
	     }
	     
	 }//...end of cut 2

	 
	 
	 // ...Armenteros-Podolanski method                               ...Luiz
	 //

	 // ...first do a graphical cut and save it to ***apcut.root***
	 
         edm::FileInPath apcutInPath("PromptLUAna/PromptAnalyzer/python/apcut.root");
	 
         TFile* cutgfile = TFile::Open(apcutInPath.fullPath().c_str());
	 TCutG* apcut = (TCutG*)gROOT->FindObject("apcut");
	 TCutG* apcut2 = (TCutG*)gROOT->FindObject("apcut2");

	 if(totcharge==0){

	   //double pstar = 0.20597 ;

	   //1234
	   if(charray[0]+charray[1] == 0)
	        {
	   //................................................
	   // LAB 12 [0,1]
	   //
	   double v0px12 = trkpx[0]+trkpx[1];
	   double v0py12 = trkpy[0]+trkpy[1];
	   double v0pz12 = trkpz[0]+trkpz[1];
	   //
	   double v0P12 = TMath::Sqrt(v0px12*v0px12 + v0py12*v0py12 + v0pz12*v0pz12);
	   //
	   double trkpx1 = trkpx[0];
	   double trkpy1 = trkpy[0];
	   double trkpz1 = trkpz[0];
	   //
	   double trkpx2 = trkpx[1];
	   double trkpy2 = trkpy[1];
	   double trkpz2 = trkpz[1];
	   //
	   double Ptrk1 = TMath::Sqrt(trkpx1*trkpx1 + trkpy1*trkpy1 + trkpz1*trkpz1);
	   double Ptrk2 = TMath::Sqrt(trkpx2*trkpx2 + trkpy2*trkpy2 + trkpz2*trkpz2);
	   //
	   double costheta1A = (trkpx1*v0px12 + trkpy1*v0py12 + trkpz1*v0pz12)/(Ptrk1*v0P12);
	   double costheta2A = (trkpx2*v0px12 + trkpy2*v0py12 + trkpz2*v0pz12)/(Ptrk2*v0P12);
	   //
  	   double p1LA = TMath::Abs((trkpx1*v0px12 + trkpy1*v0py12 + trkpz1*v0pz12)/v0P12);
	   double p2LA = TMath::Abs((trkpx2*v0px12 + trkpy2*v0py12 + trkpz2*v0pz12)/v0P12);
	   //
	   double alpha12 = (p1LA-p2LA)/(p1LA+p2LA);
	   //
	   double qt1A = Ptrk1*TMath::Sqrt(1.0-costheta1A*costheta1A);
	   double qt2A = Ptrk2*TMath::Sqrt(1.0-costheta2A*costheta2A);
	   //
	   double L1A = rimpac1/TMath::Sqrt(1.0-costheta1A*costheta1A);
	   double L2A = rimpac2/TMath::Sqrt(1.0-costheta2A*costheta2A);
	   //
	   histosTH2F["hpod12"]->Fill(alpha12,qt1A);
	   histosTH2F["hpod12"]->Fill(alpha12,qt2A);
	   histosTH2F["hpod"]->Fill(alpha12,qt1A);
	   histosTH2F["hpod"]->Fill(alpha12,qt2A);
           //std::cout << " alpha12  = " << alpha12 << std::endl;
           //std::cout << " qt1A  = " << qt1A << std::endl;
           //std::cout << " qt2A  = " << qt2A << std::endl;
	   //
	   //	   double pt12 = (m_k0*m_k0)*(1-alpha12*alpha12)/4.0 - (m_pi*m_pi);
	   //if( -0.82784 < alpha12 && alpha12 < 0.82784 &&
	   //     0.9*pstar*pstar < qt1A*qt1A && qt1A*qt1A < 1.1*pstar*pstar )
	   /////if( -0.30 < alpha12 && alpha12 < 0.80 &&
	   /////     0.09 < qt1A && qt1A < 0.22 )
	   apcut->SetVarX("alpha12");
	   apcut->SetVarY("qt1A");
	   if(apcut->IsInside(alpha12,qt1A))
	     {
	   histosTH1F["hm4rec2OSpod12"]->Fill(mrecpi1pi2);
	   //-----------------------------
	   //-----------------------------
	   // LIFETIME
	   // impact parameter has to have the same track order based on momentum 
	   // L is decay distance
	   histosTH1F["hlife1A"]->Fill(L1A);
	   histosTH1F["hlife2A"]->Fill(L2A);
	   histosTH1F["hlife12"]->Fill(L1A);
	   histosTH1F["hlife12"]->Fill(L2A);
	   histosTH1F["hlife"]->Fill(L1A);
	   histosTH1F["hlife"]->Fill(L2A);
	   histosTH1F["hlifecor"]->Fill(L1A,1/effictotal);
	   histosTH1F["hlifecor"]->Fill(L2A,1/effictotal);
	   //------------------------------
	   //------------------------------
	   }
	   //
	   // rho/w
	   /////if(  0.50 < alpha12 && alpha12 < 1.0 &&
	   /////     0.22 < qt1A && qt1A < 0.40 )
	   apcut2->SetVarX("alpha12");
	   apcut2->SetVarY("qt1A");
	   if(apcut2->IsInside(alpha12,qt1A))
	     {
	   histosTH1F["hm4rec2OSpod12rw"]->Fill(mrecpi1pi2);
	   }
	   //	   
	   double costheta1AA = (pi1px*vee12px + pi1py*vee12py + pi1pz*vee12pz)/(Ppi1*Pvee12);
	   double costheta2AA = (pi2px*vee12px + pi2py*vee12py + pi2pz*vee12pz)/(Ppi2*Pvee12);
	   //std::cout << " .......................................... " << std::endl;
	   //std::cout << " costheta1AA  = " << costheta1AA << std::endl;
	   //std::cout << " costheta2AA  = " << costheta2AA << std::endl;
  	   double p1LAA = TMath::Abs((pi1px*vee12px + pi1py*vee12py + pi1pz*vee12pz)/Pvee12);
	   double p2LAA = TMath::Abs((pi2px*vee12px + pi2py*vee12py + pi2pz*vee12pz)/Pvee12);
	   //std::cout << " p1LAA  = " << p1LAA << std::endl;
	   //std::cout << " p2LAA  = " << p2LAA << std::endl;	   
	   double alpha12AA = (p1LAA-p2LAA)/(p1LAA+p2LAA);
	   double qt1AA = Ppi1*TMath::Sqrt(1.0-costheta1AA*costheta1AA);
	   double qt2AA = Ppi2*TMath::Sqrt(1.0-costheta2AA*costheta2AA);
	   //
	   histosTH2F["hpod12vee"]->Fill(alpha12AA,qt1AA);
	   histosTH2F["hpod12vee"]->Fill(alpha12AA,qt2AA);
	   histosTH2F["hpodvee"]->Fill(alpha12AA,qt1AA);
	   histosTH2F["hpodvee"]->Fill(alpha12AA,qt2AA);
           //std::cout << " alpha12AA  = " << alpha12AA << std::endl;
           //std::cout << " qt1AA  = " << qt1AA << std::endl;
           //std::cout << " qt2AA  = " << qt2AA << std::endl;
	   //
	   //	   double pt12 = (m_k0*m_k0)*(1-alpha12AA*alpha12AA)/4.0 - (m_pi*m_pi);
	   //	   if( -0.82784 < alpha12AA && alpha12AA < 0.82784 &&
	   //           0.9*pstar*pstar < qt1AA*qt1AA && qt1AA*qt1AA < 1.1*pstar*pstar )
	   /////if( -0.30 < alpha12AA && alpha12AA < 0.80 &&
	   /////     0.09 < qt1AA && qt1AA < 0.22 )
	   apcut->SetVarX("alpha12AA");
	   apcut->SetVarY("qt1AA");
	   if(apcut->IsInside(alpha12AA,qt1AA))
	     {
	   histosTH1F["hm4rec2OSpod12vee"]->Fill(mrecpi1pi2);
	   }
	   //
	   // rho/w
	   /////if(  0.50 < alpha12AA && alpha12AA < 1.0 &&
	   /////     0.22 < qt1AA && qt1AA < 0.40 )
	   apcut2->SetVarX("alpha12AA");
	   apcut2->SetVarY("qt1AA");
	   if(apcut2->IsInside(alpha12AA,qt1AA))
	     {
	   histosTH1F["hm4rec2OSpod12veerw"]->Fill(mrecpi1pi2);
	   }
	   // LAB 34 [2,3]
	   //
	   double v0px34 = trkpx[2]+trkpx[3];
	   double v0py34 = trkpy[2]+trkpy[3];
	   double v0pz34 = trkpz[2]+trkpz[3];
	   //
	   double v0P34 = TMath::Sqrt(v0px34*v0px34 + v0py34*v0py34 + v0pz34*v0pz34);
	   //
	   double trkpx3 = trkpx[2];
	   double trkpy3 = trkpy[2];
	   double trkpz3 = trkpz[2];
	   //
	   double trkpx4 = trkpx[3];
	   double trkpy4 = trkpy[3];
	   double trkpz4 = trkpz[3];
	   //
	   double Ptrk3 = TMath::Sqrt(trkpx3*trkpx3 + trkpy3*trkpy3 + trkpz3*trkpz3);
	   double Ptrk4 = TMath::Sqrt(trkpx4*trkpx4 + trkpy4*trkpy4 + trkpz4*trkpz4);
	   //
	   double costheta3A = (trkpx3*v0px34 + trkpy3*v0py34 + trkpz3*v0pz34)/(Ptrk3*v0P34);
	   double costheta4A = (trkpx4*v0px34 + trkpy4*v0py34 + trkpz4*v0pz34)/(Ptrk4*v0P34);
	   //
  	   double p3LA = TMath::Abs((trkpx3*v0px34 + trkpy3*v0py34 + trkpz3*v0pz34)/v0P34);
	   double p4LA = TMath::Abs((trkpx4*v0px34 + trkpy4*v0py34 + trkpz4*v0pz34)/v0P34);
	   //
	   double alpha34 = (p3LA-p4LA)/(p3LA+p4LA);
	   //
	   double qt3A = Ptrk3*TMath::Sqrt(1.0-costheta3A*costheta3A);
	   double qt4A = Ptrk4*TMath::Sqrt(1.0-costheta4A*costheta4A);
	   //
	   double L3A = rimpac3/TMath::Sqrt(1.0-costheta3A*costheta3A);
	   double L4A = rimpac4/TMath::Sqrt(1.0-costheta4A*costheta4A);
	   //
	   histosTH2F["hpod34"]->Fill(alpha34,qt3A);
	   histosTH2F["hpod34"]->Fill(alpha34,qt4A);
	   histosTH2F["hpod"]->Fill(alpha34,qt3A);
	   histosTH2F["hpod"]->Fill(alpha34,qt4A);
	   //
	   //	   double pt34 = (m_k0*m_k0)*(1-alpha34*alpha34)/4.0 - (m_pi*m_pi);
	   //if( -0.82784 < alpha34 && alpha34 < 0.82784 &&
	   //     0.9*pstar*pstar < qt3A*qt3A && qt3A*qt3A < 1.1*pstar*pstar )
	   /////if( -0.30 < alpha34 && alpha34 < 0.80 &&
	   /////     0.09 < qt3A && qt3A < 0.22 )
	   apcut->SetVarX("alpha34");
	   apcut->SetVarY("qt3A");
	   if(apcut->IsInside(alpha34,qt3A))
	     {
	   histosTH1F["hm4rec2OSpod34"]->Fill(mrecpi3pi4);
	   //-----------------------------
	   //-----------------------------
	   // LIFETIME
	   // impact parameter has to have the same track order based on momentum 
	   // L is decay distance
	   histosTH1F["hlife3A"]->Fill(L3A);
	   histosTH1F["hlife4A"]->Fill(L4A);
	   histosTH1F["hlife34"]->Fill(L3A);
	   histosTH1F["hlife34"]->Fill(L4A);
	   histosTH1F["hlife"]->Fill(L3A);
	   histosTH1F["hlife"]->Fill(L4A);
	   histosTH1F["hlifecor"]->Fill(L3A,1/effictotal);
	   histosTH1F["hlifecor"]->Fill(L4A,1/effictotal);
	   //------------------------------
	   //------------------------------
	   }
	   //
	   // rho/w
	   /////if(  0.50 < alpha34 && alpha34 < 1.0 &&
	   /////     0.22 < qt3A && qt3A < 0.40 )
	   apcut2->SetVarX("alpha34");
	   apcut2->SetVarY("qt3A");
	   if(apcut2->IsInside(alpha34,qt3A))
	     {
	   histosTH1F["hm4rec2OSpod34rw"]->Fill(mrecpi3pi4);
	   }
	   // pod1234
	   //if( -0.82784 < alpha12 && alpha12 < 0.82784 &&
	   //     0.9*pstar*pstar < qt1A*qt1A && qt1A*qt1A < 1.1*pstar*pstar &&
           //    -0.82784 < alpha34 && alpha34 < 0.82784 &&
	   //     0.9*pstar*pstar < qt3A*qt3A && qt3A*qt3A < 1.1*pstar*pstar )
	   /////if( -0.30 < alpha12 && alpha12 < 0.80 &&
	   /////     0.09 < qt1A && qt1A < 0.22 &&
	   /////    -0.30 < alpha34 && alpha34 < 0.80 &&
	   /////     0.09 < qt3A && qt3A < 0.22 )
	   if(apcut->IsInside(alpha12,qt1A) &&
	      apcut->IsInside(alpha34,qt3A))
	     {
	   histosTH1F["hm4rec2OSpod1234"]->Fill(mrec);
	   }
	   // rho/w
	   /////if(  0.50 < alpha12 && alpha12 < 1.0 &&
	   /////     0.22 < qt1A && qt1A < 0.40 &&
	   /////     0.50 < alpha34 && alpha34 < 1.0 &&
	   /////     0.22 < qt3A && qt3A < 0.40 )
	   if(apcut2->IsInside(alpha12,qt1A) &&
	      apcut2->IsInside(alpha34,qt3A))
	     {
	   histosTH1F["hm4rec2OSpod1234rw"]->Fill(mrec);
	   }
	   //
	   double costheta3AA = (pi3px*vee34px + pi3py*vee34py + pi3pz*vee34pz)/(Ppi3*Pvee34);
	   double costheta4AA = (pi4px*vee34px + pi4py*vee34py + pi4pz*vee34pz)/(Ppi4*Pvee34);
    	   double p3LAA = TMath::Abs((pi3px*vee34px + pi3py*vee34py + pi3pz*vee34pz)/Pvee34);
	   double p4LAA = TMath::Abs((pi4px*vee34px + pi4py*vee34py + pi4pz*vee34pz)/Pvee34);
	   double alpha34AA = (p3LAA-p4LAA)/(p3LAA+p4LAA);
	   double qt3AA =  Ppi3*TMath::Sqrt(1.0-costheta3AA*costheta3AA);
	   double qt4AA =  Ppi4*TMath::Sqrt(1.0-costheta4AA*costheta4AA);
	   //
	   histosTH2F["hpod34vee"]->Fill(alpha34AA,qt3AA);
	   histosTH2F["hpod34vee"]->Fill(alpha34AA,qt4AA);
	   histosTH2F["hpodvee"]->Fill(alpha34AA,qt3AA);
	   histosTH2F["hpodvee"]->Fill(alpha34AA,qt4AA);
	   //
	   //double pt34AA = (m_k0*m_k0)*(1-alpha34AA*alpha34AA)/4.0 - (m_pi*m_pi);
	   //if( -0.82784 < alpha34AA && alpha34AA < 0.82784 &&
	   //     0.9*pstar*pstar < qt3AA*qt3AA && qt3AA*qt3AA < 1.1*pstar*pstar )
	   /////if( -0.30 < alpha34AA && alpha34AA < 0.80 &&
	   /////     0.09 < qt3AA && qt3AA < 0.22 )
	   apcut->SetVarX("alpha34AA");
	   apcut->SetVarY("qt3AA");
	   if(apcut->IsInside(alpha34AA,qt3AA))
	     {
	   histosTH1F["hm4rec2OSpod34vee"]->Fill(mrecpi3pi4);
	   }
	   // rho/w
	   /////if(  0.50 < alpha34AA && alpha34AA < 1.0 &&
	   /////     0.22 < qt3AA && qt3AA < 0.40 )
	   apcut2->SetVarX("alpha34AA");
	   apcut2->SetVarY("qt3AA");
	   if(apcut2->IsInside(alpha34AA,qt3AA))
	     {
	   histosTH1F["hm4rec2OSpod34veerw"]->Fill(mrecpi3pi4);
	   }
	   // pod1234vee
	   //if( -0.82784 < alpha12AA && alpha12AA < 0.82784 &&
	   //     0.9*pstar*pstar < qt1AA*qt1AA && qt1AA*qt1AA < 1.1*pstar*pstar &&
           //    -0.82784 < alpha34AA && alpha34AA < 0.82784 &&
	   //     0.9*pstar*pstar < qt3AA*qt3AA && qt3AA*qt3AA < 1.1*pstar*pstar )
	   /////if( -0.30 < alpha12AA && alpha12AA < 0.80 &&
	   /////     0.09 < qt1AA && qt1AA < 0.22 &&
	   /////    -0.30 < alpha34AA && alpha34AA < 0.80 &&
	   /////     0.09 < qt3AA && qt3AA < 0.22 )
	   if(apcut->IsInside(alpha12AA,qt1AA) &&
	      apcut->IsInside(alpha34AA,qt3AA))
	     {
	   histosTH1F["hm4rec2OSpod1234vee"]->Fill(mrec);
	   }
	   // rho/w
	   /////if(  0.50 < alpha12AA && alpha12AA < 1.0 &&
	   /////     0.22 < qt1AA && qt1AA < 0.40 &&
	   /////     0.50 < alpha34AA && alpha34AA < 1.0 &&
	   /////     0.22 < qt3AA && qt3AA < 0.40 )
	   if(apcut2->IsInside(alpha12AA,qt1AA) &&
	      apcut2->IsInside(alpha34AA,qt3AA))
	     {
	   histosTH1F["hm4rec2OSpod1234veerw"]->Fill(mrec);
	   }
		}//...end Q_pair=0 ...1234

	   
	   //................................................
	   //1324
	   if(charray[0]+charray[2] == 0)
	        {
	   //................................................
	   // LAB 13 [0,2]
	   //
	   double v0px13 = trkpx[0]+trkpx[2];
	   double v0py13 = trkpy[0]+trkpy[2];
	   double v0pz13 = trkpz[0]+trkpz[2];
	   //
	   double v0P13 = TMath::Sqrt(v0px13*v0px13 + v0py13*v0py13 + v0pz13*v0pz13);
	   //
	   double trkpx1 = trkpx[0];
	   double trkpy1 = trkpy[0];
	   double trkpz1 = trkpz[0];
	   //
	   double trkpx3 = trkpx[2];
	   double trkpy3 = trkpy[2];
	   double trkpz3 = trkpz[2];
	   //
	   double Ptrk1 = TMath::Sqrt(trkpx1*trkpx1 + trkpy1*trkpy1 + trkpz1*trkpz1);
	   double Ptrk3 = TMath::Sqrt(trkpx3*trkpx3 + trkpy3*trkpy3 + trkpz3*trkpz3);
	   //
	   double costheta1B = (trkpx1*v0px13 + trkpy1*v0py13 + trkpz1*v0pz13)/(Ptrk1*v0P13);
	   double costheta3B = (trkpx3*v0px13 + trkpy3*v0py13 + trkpz3*v0pz13)/(Ptrk3*v0P13);
	   //
  	   double p1LB = TMath::Abs((trkpx1*v0px13 + trkpy1*v0py13 + trkpz1*v0pz13)/v0P13);
	   double p3LB = TMath::Abs((trkpx3*v0px13 + trkpy3*v0py13 + trkpz3*v0pz13)/v0P13);
	   //
	   double alpha13 = (p1LB-p3LB)/(p1LB+p3LB);
	   //
	   double qt1B = Ptrk1*TMath::Sqrt(1.0-costheta1B*costheta1B);
	   double qt3B = Ptrk3*TMath::Sqrt(1.0-costheta3B*costheta3B);
	   //
	   double L1B = rimpac1/TMath::Sqrt(1.0-costheta1B*costheta1B);
	   double L3B = rimpac3/TMath::Sqrt(1.0-costheta3B*costheta3B);
	   //
	   histosTH2F["hpod13"]->Fill(alpha13,qt1B);
	   histosTH2F["hpod13"]->Fill(alpha13,qt3B);
	   histosTH2F["hpod"]->Fill(alpha13,qt1B);
	   histosTH2F["hpod"]->Fill(alpha13,qt3B);
	   //
	   //	   double pt13 = (m_k0*m_k0)*(1-alpha13*alpha13)/4.0 - (m_pi*m_pi);
	   //if( -0.82784 < alpha13 && alpha13 < 0.82784 &&
	   //     0.9*pstar*pstar < qt1B*qt1B && qt1B*qt1B < 1.1*pstar*pstar )
	   /////if( -0.30 < alpha13 && alpha13 < 0.80 &&
	   /////     0.09 < qt1B && qt1B < 0.22 )
	   apcut->SetVarX("alpha13");
	   apcut->SetVarY("qt1B");
	   if(apcut->IsInside(alpha13,qt1B))
	     {
	   histosTH1F["hm4rec2OSpod13"]->Fill(mrecpi1pi3);
	   //-----------------------------
	   //-----------------------------
	   // LIFETIME
	   // impact parameter has to have the same track order based on momentum 
	   // L is decay distance
	   histosTH1F["hlife1B"]->Fill(L1B);
	   histosTH1F["hlife3B"]->Fill(L3B);
	   histosTH1F["hlife13"]->Fill(L1B);
	   histosTH1F["hlife13"]->Fill(L3B);
	   histosTH1F["hlife"]->Fill(L1B);
	   histosTH1F["hlife"]->Fill(L3B);
	   histosTH1F["hlifecor"]->Fill(L1B,1/effictotal);
	   histosTH1F["hlifecor"]->Fill(L3B,1/effictotal);
	   //------------------------------
	   //------------------------------
	   }
	   //
	   // rho/w
	   /////if(  0.50 < alpha13 && alpha13 < 1.0 &&
	   /////     0.22 < qt1B && qt1B < 0.40 )
	   apcut2->SetVarX("alpha13");
	   apcut2->SetVarY("qt1B");
	   if(apcut2->IsInside(alpha13,qt1B))
	     {
	   histosTH1F["hm4rec2OSpod13rw"]->Fill(mrecpi1pi3);
	   }
	   //   
	   double costheta1BB = (pi1px*vee13px + pi1py*vee13py + pi1pz*vee13pz)/(Ppi1*Pvee13);
	   double costheta3BB = (pi3px*vee13px + pi3py*vee13py + pi3pz*vee13pz)/(Ppi3*Pvee13);
  	   double p1LBB = TMath::Abs((pi1px*vee13px + pi1py*vee13py + pi1pz*vee13pz)/Pvee13);
	   double p3LBB = TMath::Abs((pi3px*vee13px + pi3py*vee13py + pi3pz*vee13pz)/Pvee13);
	   double alpha13BB = (p1LBB-p3LBB)/(p1LBB+p3LBB);
	   double qt1BB = Ppi1*TMath::Sqrt(1.0-costheta1BB*costheta1BB);
	   double qt3BB = Ppi3*TMath::Sqrt(1.0-costheta3BB*costheta3BB);
           //
	   histosTH2F["hpod13vee"]->Fill(alpha13BB,qt1BB);
	   histosTH2F["hpod13vee"]->Fill(alpha13BB,qt3BB);
	   histosTH2F["hpodvee"]->Fill(alpha13BB,qt1BB);
	   histosTH2F["hpodvee"]->Fill(alpha13BB,qt3BB);
	   //
	   //double pt13BB = (m_k0*m_k0)*(1-alpha13BB*alpha13BB)/4.0 - (m_pi*m_pi);
	   //if( -0.82784 < alpha13BB && alpha13BB < 0.82784 &&
	   //     0.9*pstar*pstar < qt1BB*qt1BB && qt1BB*qt1BB < 1.1*pstar*pstar )
	   /////if( -0.30 < alpha13BB && alpha13BB < 0.80 &&
	   /////     0.09 < qt1BB && qt1BB < 0.22 )
	   apcut->SetVarX("alpha13BB");
	   apcut->SetVarY("qt1BB");
	   if(apcut->IsInside(alpha13BB,qt1BB))
	     {
	   histosTH1F["hm4rec2OSpod13vee"]->Fill(mrecpi1pi3);
	   }
	   // rho/w
	   /////if(  0.50 < alpha13BB && alpha13BB < 1.0 &&
	   /////     0.22 < qt1BB && qt1BB < 0.40 )
	   apcut2->SetVarX("alpha13BB");
	   apcut2->SetVarY("qt1BB");
	   if(apcut2->IsInside(alpha13BB,qt1BB))
	     {
	   histosTH1F["hm4rec2OSpod13veerw"]->Fill(mrecpi1pi3);
	   }
	   // LAB 24 [1,3]
	   //
	   double v0px24 = trkpx[1]+trkpx[3];
	   double v0py24 = trkpy[1]+trkpy[3];
	   double v0pz24 = trkpz[1]+trkpz[3];
	   //
	   double v0P24 = TMath::Sqrt(v0px24*v0px24 + v0py24*v0py24 + v0pz24*v0pz24);
	   //
	   double trkpx2 = trkpx[1];
	   double trkpy2 = trkpy[1];
	   double trkpz2 = trkpz[1];
	   //
	   double trkpx4 = trkpx[3];
	   double trkpy4 = trkpy[3];
	   double trkpz4 = trkpz[3];
	   //
	   double Ptrk2 = TMath::Sqrt(trkpx2*trkpx2 + trkpy2*trkpy2 + trkpz2*trkpz2);
	   double Ptrk4 = TMath::Sqrt(trkpx4*trkpx4 + trkpy4*trkpy4 + trkpz4*trkpz4);
	   //
	   double costheta2B = (trkpx2*v0px24 + trkpy2*v0py24 + trkpz2*v0pz24)/(Ptrk2*v0P24);
	   double costheta4B = (trkpx4*v0px24 + trkpy4*v0py24 + trkpz4*v0pz24)/(Ptrk4*v0P24);
	   //
  	   double p2LB = TMath::Abs((trkpx2*v0px24 + trkpy2*v0py24 + trkpz2*v0pz24)/v0P24);
	   double p4LB = TMath::Abs((trkpx4*v0px24 + trkpy4*v0py24 + trkpz4*v0pz24)/v0P24);
	   //
	   double alpha24 = (p2LB-p4LB)/(p2LB+p4LB);
	   //
	   double qt2B = Ptrk2*TMath::Sqrt(1.0-costheta2B*costheta2B);
	   double qt4B = Ptrk4*TMath::Sqrt(1.0-costheta4B*costheta4B);
	   //
	   double L2B = rimpac2/TMath::Sqrt(1.0-costheta2B*costheta2B);
	   double L4B = rimpac4/TMath::Sqrt(1.0-costheta4B*costheta4B);
	   //
	   histosTH2F["hpod24"]->Fill(alpha24,qt2B);
	   histosTH2F["hpod24"]->Fill(alpha24,qt4B);
	   histosTH2F["hpod"]->Fill(alpha24,qt2B);
	   histosTH2F["hpod"]->Fill(alpha24,qt4B);
	   //
	   //double pt24 = (m_k0*m_k0)*(1-alpha24*alpha24)/4.0 - (m_pi*m_pi);
	   //if( -0.82784 < alpha24 && alpha24 < 0.82784 &&
	   //     0.9*pstar*pstar < qt3B*qt3B && qt3B*qt3B < 1.1*pstar*pstar )
	   /////if( -0.30 < alpha24 && alpha24 < 0.80 &&
	   /////     0.09 < qt2B && qt2B < 0.22 )
	   apcut->SetVarX("alpha24");
	   apcut->SetVarY("qt2B");
	   if(apcut->IsInside(alpha24,qt2B))
	     {
	   histosTH1F["hm4rec2OSpod24"]->Fill(mrecpi2pi4);
	   //-----------------------------
	   //-----------------------------
	   // LIFETIME
	   // impact parameter has to have the same track order based on momentum 
	   // L is decay distance
	   histosTH1F["hlife2B"]->Fill(L2B);
	   histosTH1F["hlife4B"]->Fill(L4B);
	   histosTH1F["hlife24"]->Fill(L2B);
	   histosTH1F["hlife24"]->Fill(L4B);
	   histosTH1F["hlife"]->Fill(L2B);
	   histosTH1F["hlife"]->Fill(L4B);
	   histosTH1F["hlifecor"]->Fill(L2B,1/effictotal);
	   histosTH1F["hlifecor"]->Fill(L4B,1/effictotal);
	   }
	   //------------------------------
	   //------------------------------
	   //
	   // rho/w
	   /////if(  0.50 < alpha24 && alpha24 < 1.0 &&
	   /////     0.22 < qt2B && qt2B < 0.40 )
	   apcut2->SetVarX("alpha24");
	   apcut2->SetVarY("qt2B");
	   if(apcut2->IsInside(alpha24,qt2B))
	     {
	   histosTH1F["hm4rec2OSpod24rw"]->Fill(mrecpi2pi4);
	   }
	   // pod1324
	   //if( -0.82784 < alpha13 && alpha13 < 0.82784 &&
	   //     0.9*pstar*pstar < qt1B*qt1B && qt1B*qt1B < 1.1*pstar*pstar &&
           //    -0.82784 < alpha24 && alpha24 < 0.82784 &&
	   //     0.9*pstar*pstar < qt2B*qt2B && qt2B*qt2B < 1.1*pstar*pstar )
	   /////if( -0.30 < alpha13 && alpha13 < 0.80 &&
	   /////     0.09 < qt1B && qt1B < 0.22 &&
           /////    -0.30 < alpha24 && alpha24 < 0.80 &&
	   /////     0.09 < qt2B && qt2B < 0.22 )
	   if(apcut->IsInside(alpha13,qt1B) &&
	      apcut->IsInside(alpha24,qt2B))
	     {
	   histosTH1F["hm4rec2OSpod1324"]->Fill(mrec);
	   }
	   // rho/w
	   /////if(  0.50 < alpha13 && alpha13 < 1.0 &&
	   /////     0.22 < qt1B && qt1B < 0.40 &&
           /////     0.50 < alpha24 && alpha24 < 1.0 &&
	   /////     0.22 < qt2B && qt2B < 0.40 )
	   if(apcut2->IsInside(alpha13,qt1B) &&
	      apcut2->IsInside(alpha24,qt2B))
	     {
	   histosTH1F["hm4rec2OSpod1324rw"]->Fill(mrec);
	   }
	   //
	   double costheta2BB = (pi2px*vee24px + pi2py*vee24py + pi2pz*vee24pz)/(Ppi2*Pvee24);
	   double costheta4BB = (pi4px*vee24px + pi4py*vee24py + pi4pz*vee24pz)/(Ppi4*Pvee24);
    	   double p2LBB = TMath::Abs((pi2px*vee24px + pi2py*vee24py + pi2pz*vee24pz)/Pvee24);
	   double p4LBB = TMath::Abs((pi4px*vee24px + pi4py*vee24py + pi4pz*vee24pz)/Pvee24);
	   double alpha24BB = (p2LBB-p4LBB)/(p2LBB+p4LBB);
	   double qt2BB =  Ppi2*TMath::Sqrt(1.0-costheta2BB*costheta2BB);
	   double qt4BB =  Ppi4*TMath::Sqrt(1.0-costheta4BB*costheta4BB);
	   //
	   histosTH2F["hpod24vee"]->Fill(alpha24BB,qt2BB);
	   histosTH2F["hpod24vee"]->Fill(alpha24BB,qt4BB);
	   histosTH2F["hpodvee"]->Fill(alpha24BB,qt2BB);
	   histosTH2F["hpodvee"]->Fill(alpha24BB,qt4BB);
	   //
	   //double pt24BB = (m_k0*m_k0)*(1-alpha24BB*alpha24BB)/4.0 - (m_pi*m_pi);
	   //if( -0.82784 < alpha24BB && alpha24BB < 0.82784 &&
	   //     0.9*pstar*pstar < qt2BB*qt2BB && qt2BB*qt2BB < 1.1*pstar*pstar )
	   /////if( -0.30 < alpha24BB && alpha24BB < 0.80 &&
	   /////     0.09 < qt2BB && qt2BB < 0.22 )
	   apcut->SetVarX("alpha24BB");
	   apcut->SetVarY("qt2BB");
	   if(apcut->IsInside(alpha24BB,qt2BB))
	     {
	   histosTH1F["hm4rec2OSpod24vee"]->Fill(mrecpi2pi4);
	   }
	   // rho/w
	   /////if(  0.50 < alpha24BB && alpha24BB < 1.0 &&
	   /////     0.22 < qt2BB && qt2BB < 0.40 )
	   apcut2->SetVarX("alpha24BB");
	   apcut2->SetVarY("qt2BB");
	   if(apcut2->IsInside(alpha24BB,qt2BB))
	     {
	   histosTH1F["hm4rec2OSpod24veerw"]->Fill(mrecpi2pi4);
	   }
	   // pod1324vee
	   //if( -0.82784 < alpha13BB && alpha13BB < 0.82784 &&
	   //     0.9*pstar*pstar < qt1BB*qt1BB && qt1BB*qt1BB < 1.1*pstar*pstar &&
           //    -0.82784 < alpha24BB && alpha24BB < 0.82784 &&
	   //     0.9*pstar*pstar < qt2BB*qt2BB && qt2BB*qt2BB < 1.1*pstar*pstar )
	   /////if( -0.30 < alpha13BB && alpha13BB < 0.80 &&
	   /////     0.09 < qt1BB && qt1BB < 0.22 &&
           /////    -0.30 < alpha24BB && alpha24BB < 0.80 &&
	   /////     0.09 < qt2BB && qt2BB < 0.22 )
	   if(apcut->IsInside(alpha13BB,qt1BB) &&
	      apcut->IsInside(alpha24BB,qt2BB))
	     {
	   histosTH1F["hm4rec2OSpod1324vee"]->Fill(mrec);
	   }
	   // rho/w
	   /////if(  0.50 < alpha13BB && alpha13BB < 1.0 &&
	   /////     0.22 < qt1BB && qt1BB < 0.40 &&
           /////     0.50 < alpha24BB && alpha24BB < 1.0 &&
	   /////     0.22 < qt2BB && qt2BB < 0.40 )
	   if(apcut2->IsInside(alpha13BB,qt1BB) &&
	      apcut2->IsInside(alpha24BB,qt2BB))
	     {
	   histosTH1F["hm4rec2OSpod1324veerw"]->Fill(mrec);
	   }
		}//...end Q_pair=0 ...1324

	   
	   //................................................
	   //1423
	   if(charray[0]+charray[3] == 0)
	        {
	   //................................................
	   // LAB 14 [0,3]
	   //
	   double v0px14 = trkpx[0]+trkpx[3];
	   double v0py14 = trkpy[0]+trkpy[3];
	   double v0pz14 = trkpz[0]+trkpz[3];
	   //
	   double v0P14 = TMath::Sqrt(v0px14*v0px14 + v0py14*v0py14 + v0pz14*v0pz14);
	   //
	   double trkpx1 = trkpx[0];
	   double trkpy1 = trkpy[0];
	   double trkpz1 = trkpz[0];
	   //
	   double trkpx4 = trkpx[3];
	   double trkpy4 = trkpy[3];
	   double trkpz4 = trkpz[3];
	   //
	   double Ptrk1 = TMath::Sqrt(trkpx1*trkpx1 + trkpy1*trkpy1 + trkpz1*trkpz1);
	   double Ptrk4 = TMath::Sqrt(trkpx4*trkpx4 + trkpy4*trkpy4 + trkpz4*trkpz4);
	   //
	   double costheta1C = (trkpx1*v0px14 + trkpy1*v0py14 + trkpz1*v0pz14)/(Ptrk1*v0P14);
	   double costheta4C = (trkpx4*v0px14 + trkpy4*v0py14 + trkpz4*v0pz14)/(Ptrk4*v0P14);
	   //
  	   double p1LC = TMath::Abs((trkpx1*v0px14 + trkpy1*v0py14 + trkpz1*v0pz14)/v0P14);
	   double p4LC = TMath::Abs((trkpx4*v0px14 + trkpy4*v0py14 + trkpz4*v0pz14)/v0P14);
	   //
	   double alpha14 = (p1LC-p4LC)/(p1LC+p4LC);
	   //
	   double qt1C = Ptrk1*TMath::Sqrt(1.0-costheta1C*costheta1C);
	   double qt4C = Ptrk4*TMath::Sqrt(1.0-costheta4C*costheta4C);
	   //
	   double L1C = rimpac1/TMath::Sqrt(1.0-costheta1C*costheta1C);
	   double L4C = rimpac4/TMath::Sqrt(1.0-costheta4C*costheta4C);
	   //
	   histosTH2F["hpod14"]->Fill(alpha14,qt1C);
	   histosTH2F["hpod14"]->Fill(alpha14,qt4C);
	   histosTH2F["hpod"]->Fill(alpha14,qt1C);
	   histosTH2F["hpod"]->Fill(alpha14,qt4C);
	   //
	   //double pt14 = (m_k0*m_k0)*(1-alpha14*alpha14)/4.0 - (m_pi*m_pi);
	   //if( -0.82784 < alpha14 && alpha14 < 0.82784 &&
	   //     0.9*pstar*pstar < qt1C*qt1C && qt1C*qt1C < 1.1*pstar*pstar )
	   /////if( -0.30 < alpha14 && alpha14 < 0.80 &&
	   /////     0.09 < qt1C && qt1C < 0.22 )
	   apcut->SetVarX("alpha14");
	   apcut->SetVarY("qt1C");
	   if(apcut->IsInside(alpha14,qt1C))
	     {
	   histosTH1F["hm4rec2OSpod14"]->Fill(mrecpi1pi4);
	   //-----------------------------
	   //-----------------------------
	   // LIFETIME
	   // impact parameter has to have the same track order based on momentum 
	   // L is decay distance
	   histosTH1F["hlife1C"]->Fill(L1C);
	   histosTH1F["hlife4C"]->Fill(L4C);
	   histosTH1F["hlife14"]->Fill(L1C);
	   histosTH1F["hlife14"]->Fill(L4C);
	   histosTH1F["hlife"]->Fill(L1C);
	   histosTH1F["hlife"]->Fill(L4C);
	   histosTH1F["hlifecor"]->Fill(L1C,1/effictotal);
	   histosTH1F["hlifecor"]->Fill(L4C,1/effictotal);
	   //------------------------------
	   //------------------------------
	   }
	   //
	   // rho/w
	   /////if(  0.50 < alpha14 && alpha14 < 1.0 &&
	   /////     0.22 < qt1C && qt1C < 0.40 )
	   apcut2->SetVarX("alpha14");
	   apcut2->SetVarY("qt1C");
	   if(apcut2->IsInside(alpha14,qt1C))
	     {
	   histosTH1F["hm4rec2OSpod14rw"]->Fill(mrecpi1pi4);
	   }
	   //
	   double costheta1CC = (pi1px*vee14px + pi1py*vee14py + pi1pz*vee14pz)/(Ppi1*Pvee14);
	   double costheta4CC = (pi4px*vee14px + pi4py*vee14py + pi4pz*vee14pz)/(Ppi4*Pvee14);
  	   double p1LCC = TMath::Abs((pi1px*vee14px + pi1py*vee14py + pi1pz*vee14pz)/Pvee14);
	   double p4LCC = TMath::Abs((pi4px*vee14px + pi4py*vee14py + pi4pz*vee14pz)/Pvee14);
	   double alpha14CC = (p1LCC-p4LCC)/(p1LCC+p4LCC);
	   double qt1CC = Ppi1*TMath::Sqrt(1.0-costheta1CC*costheta1CC);
	   double qt4CC = Ppi4*TMath::Sqrt(1.0-costheta4CC*costheta4CC);
	   //
	   histosTH2F["hpod14vee"]->Fill(alpha14CC,qt1CC);
	   histosTH2F["hpod14vee"]->Fill(alpha14CC,qt4CC);
	   histosTH2F["hpodvee"]->Fill(alpha14CC,qt1CC);
	   histosTH2F["hpodvee"]->Fill(alpha14CC,qt4CC);
	   //
	   //double pt14CC = (m_k0*m_k0)*(1-alpha14CC*alpha14CC)/4.0 - (m_pi*m_pi);
	   //if( -0.82784 < alpha14CC && alpha14CC < 0.82784 &&
	   //     0.9*pstar*pstar < qt1CC*qt1CC && qt1CC*qt1CC < 1.1*pstar*pstar )
	   /////if( -0.30 < alpha14CC && alpha14CC < 0.80 &&
	   /////     0.09 < qt1CC && qt1CC < 0.22 )
	   apcut->SetVarX("alpha14CC");
	   apcut->SetVarY("qt1CC");
	   if(apcut->IsInside(alpha14CC,qt1CC))
	     {
	   histosTH1F["hm4rec2OSpod14vee"]->Fill(mrecpi1pi4);
	   }
	   // rho/w
	   /////if(  0.50 < alpha14CC && alpha14CC < 1.0 &&
	   /////     0.22 < qt1CC && qt1CC < 0.40 )
	   apcut2->SetVarX("alpha14CC");
	   apcut2->SetVarY("qt1CC");
	   if(apcut2->IsInside(alpha14CC,qt1CC))
	     {
	   histosTH1F["hm4rec2OSpod14veerw"]->Fill(mrecpi1pi4);
	   }
	   // LAB 23 [1,2]
	   //
	   double v0px23 = trkpx[1]+trkpx[2];
	   double v0py23 = trkpy[1]+trkpy[2];
	   double v0pz23 = trkpz[1]+trkpz[2];
	   //
	   double v0P23 = TMath::Sqrt(v0px23*v0px23 + v0py23*v0py23 + v0pz23*v0pz23);
	   //
	   double trkpx2 = trkpx[1];
	   double trkpy2 = trkpy[1];
	   double trkpz2 = trkpz[1];
	   //
	   double trkpx3 = trkpx[2];
	   double trkpy3 = trkpy[2];
	   double trkpz3 = trkpz[2];
	   //
	   double Ptrk2 = TMath::Sqrt(trkpx2*trkpx2 + trkpy2*trkpy2 + trkpz2*trkpz2);
	   double Ptrk3 = TMath::Sqrt(trkpx3*trkpx3 + trkpy3*trkpy3 + trkpz3*trkpz3);
	   //
	   double costheta2C = (trkpx2*v0px23 + trkpy2*v0py23 + trkpz2*v0pz23)/(Ptrk2*v0P23);
	   double costheta3C = (trkpx3*v0px23 + trkpy3*v0py23 + trkpz3*v0pz23)/(Ptrk3*v0P23);
	   //
  	   double p2LC = TMath::Abs((trkpx2*v0px23 + trkpy2*v0py23 + trkpz2*v0pz23)/v0P23);
	   double p3LC = TMath::Abs((trkpx3*v0px23 + trkpy3*v0py23 + trkpz3*v0pz23)/v0P23);
	   //
	   double alpha23 = (p2LC-p3LC)/(p2LC+p3LC);
	   //
	   double qt2C = Ptrk2*TMath::Sqrt(1.0-costheta2C*costheta2C);
	   double qt3C = Ptrk3*TMath::Sqrt(1.0-costheta3C*costheta3C);
	   //
	   double L2C = rimpac2/TMath::Sqrt(1.0-costheta2C*costheta2C);
	   double L3C = rimpac3/TMath::Sqrt(1.0-costheta3C*costheta3C);
	   //
	   histosTH2F["hpod23"]->Fill(alpha23,qt2C);
	   histosTH2F["hpod23"]->Fill(alpha23,qt3C);
	   histosTH2F["hpod"]->Fill(alpha23,qt2C);
	   histosTH2F["hpod"]->Fill(alpha23,qt3C);
	   //
	   //double pt23 = (m_k0*m_k0)*(1-alpha23*alpha23)/4.0 - (m_pi*m_pi);
	   //if( -0.82784 < alpha23 && alpha23 < 0.82784 &&
	   //     0.9*pstar*pstar < qt2C*qt2C && qt2C*qt2C < 1.1*pstar*pstar )
	   /////if( -0.30 < alpha23 && alpha23 < 0.80 &&
	   /////     0.09 < qt2C && qt2C < 0.22 )
	   apcut->SetVarX("alpha23");
	   apcut->SetVarY("qt2C");
	   if(apcut->IsInside(alpha23,qt2C))
	     {
	   histosTH1F["hm4rec2OSpod23"]->Fill(mrecpi2pi3);
	   //-----------------------------
	   //-----------------------------
	   // LIFETIME
	   // impact parameter has to have the same track order based on momentum 
	   // L is decay distance
	   histosTH1F["hlife2C"]->Fill(L2C);
	   histosTH1F["hlife3C"]->Fill(L3C);
	   histosTH1F["hlife23"]->Fill(L2C);
	   histosTH1F["hlife23"]->Fill(L3C);
	   histosTH1F["hlife"]->Fill(L2C);
	   histosTH1F["hlife"]->Fill(L3C);
	   histosTH1F["hlifecor"]->Fill(L2C,1/effictotal);
	   histosTH1F["hlifecor"]->Fill(L3C,1/effictotal);
	   //------------------------------
	   //------------------------------
	   }
	   //
	   // rho/w
	   /////if(  0.50 < alpha23 && alpha23 < 1.0 &&
	   /////     0.22 < qt2C && qt2C < 0.40 )
	   apcut2->SetVarX("alpha23");
	   apcut2->SetVarY("qt2C");
	   if(apcut2->IsInside(alpha23,qt2C))
	     {
	   histosTH1F["hm4rec2OSpod23rw"]->Fill(mrecpi2pi3);
	   }
	   // pod1423
	   //if( -0.82784 < alpha14 && alpha14 < 0.82784 &&
	   //     0.9*pstar*pstar < qt1C*qt1C && qt1C*qt1C < 1.1*pstar*pstar &&
           //    -0.82784 < alpha23  && alpha23 < 0.82784 &&
	   //     0.9*pstar*pstar < qt2C*qt2C && qt2C*qt2C < 1.1*pstar*pstar )
	   /////if( -0.30 < alpha14 && alpha14 < 0.80 &&
	   /////     0.09 < qt1C && qt1C < 0.22 &&
           /////    -0.30 < alpha23 && alpha23 < 0.80 &&
	   /////     0.09 < qt2C && qt2C < 0.22 )
	   if(apcut->IsInside(alpha14,qt1C) &&
	      apcut->IsInside(alpha23,qt2C))
	     {
	   histosTH1F["hm4rec2OSpod1423"]->Fill(mrec);
	   }
	   // rho/w
	   /////if(  0.50 < alpha14 && alpha14 < 1.0 &&
	   /////     0.22 < qt1C && qt1C < 0.40 &&
           /////     0.50 < alpha23 && alpha23 < 1.0 &&
	   /////     0.22 < qt2C && qt2C < 0.40 )
           if(apcut2->IsInside(alpha14,qt1C) &&
	      apcut2->IsInside(alpha23,qt2C))
	     {
	   histosTH1F["hm4rec2OSpod1423rw"]->Fill(mrec);
	   }
	   //
	   double costheta2CC = (pi2px*vee23px + pi2py*vee23py + pi2pz*vee23pz)/(Ppi2*Pvee23);
	   double costheta3CC = (pi3px*vee23px + pi3py*vee23py + pi3pz*vee23pz)/(Ppi3*Pvee23);
    	   double p2LCC = TMath::Abs((pi2px*vee23px + pi2py*vee23py + pi2pz*vee23pz)/Pvee23);
	   double p3LCC = TMath::Abs((pi3px*vee23px + pi3py*vee23py + pi3pz*vee23pz)/Pvee23);
	   double alpha23CC = (p2LCC-p3LCC)/(p2LCC+p3LCC);
	   double qt2CC =  Ppi2*TMath::Sqrt(1.0-costheta2CC*costheta2CC);
	   double qt3CC =  Ppi3*TMath::Sqrt(1.0-costheta3CC*costheta3CC);
	   //
	   histosTH2F["hpod23vee"]->Fill(alpha23CC,qt2CC);
	   histosTH2F["hpod23vee"]->Fill(alpha23CC,qt3CC);
	   histosTH2F["hpodvee"]->Fill(alpha23CC,qt2CC);
	   histosTH2F["hpodvee"]->Fill(alpha23CC,qt3CC);
	   //
	   //double pt23CC = (m_k0*m_k0)*(1-alpha23CC*alpha23CC)/4.0 - (m_pi*m_pi);
	   //if( -0.82784 < alpha23CC && alpha23CC < 0.82784 &&
	   //     0.9*pstar*pstar < qt2CC*qt2CC && qt2CC*qt2CC < 1.1*pstar*pstar )
	   /////if( -0.30 < alpha23CC && alpha23CC < 0.80 &&
	   /////     0.09 < qt2CC && qt2CC < 0.22 )
	   apcut->SetVarX("alpha23CC");
	   apcut->SetVarY("qt2CC");
	   if(apcut->IsInside(alpha23CC,qt2CC))
	     {
	   histosTH1F["hm4rec2OSpod23vee"]->Fill(mrecpi2pi3);
	   }
	   // rho/w
	   /////if(  0.50 < alpha23CC && alpha23CC < 1.0 &&
	   /////     0.22 < qt2CC && qt2CC < 0.40 )
	   apcut2->SetVarX("alpha23CC");
	   apcut2->SetVarY("qt2CC");
	   if(apcut2->IsInside(alpha23CC,qt2CC))
	     {
	   histosTH1F["hm4rec2OSpod23veerw"]->Fill(mrecpi2pi3);
	   }
	   // pod1423vee
	   //if( -0.82784 < alpha14CC && alpha14CC < 0.82784 &&
	   //     0.9*pstar*pstar < qt1CC*qt1CC && qt1CC*qt1CC < 1.1*pstar*pstar &&
           //    -0.82784 < alpha23CC  && alpha23CC < 0.82784 &&
	   //     0.9*pstar*pstar < qt2CC*qt2CC && qt2CC*qt2CC < 1.1*pstar*pstar )
	   /////if( -0.30 < alpha14CC && alpha14CC < 0.80 &&
	   /////     0.09 < qt1CC && qt1CC < 0.22 &&
           /////    -0.30 < alpha23CC && alpha23CC < 0.80 &&
	   /////     0.09 < qt2CC && qt2CC < 0.22 )
	   if(apcut->IsInside(alpha14CC,qt1CC) &&
	      apcut->IsInside(alpha23CC,qt2CC))
	     {
	   histosTH1F["hm4rec2OSpod1423vee"]->Fill(mrec);
	   }
	   // rho/w
	   /////if(  0.50 < alpha14CC && alpha14CC < 1.0 &&
	   /////     0.22 < qt1CC && qt1CC < 0.40 &&
           /////     0.50 < alpha23CC && alpha23CC < 1.0 &&
	   /////     0.22 < qt2CC && qt2CC < 0.40 )
	   if(apcut2->IsInside(alpha14CC,qt1CC) &&
	      apcut2->IsInside(alpha23CC,qt2CC))
	     {
	   histosTH1F["hm4rec2OSpod1423veerw"]->Fill(mrec);
	   }
		}//...end Q_pair=0 ...1423
	   
	 } // Q=0

	 cutgfile->Close();
	 // ...end of AP
	 
 
  } //...end of fiducialRegion4 && allCuts4

    
  //------------------------------------

#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   //   iEvent.getByLabel("example",pIn);
   //...Luiz
   iEvent.getByToken(exampletoken, pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif

}


// ------------ method called once each job just before starting event loop  ------------
void
PromptAnalyzer::beginJob()
{
    //...py efficiency  ...Luiz
    //............................reading py versus efficiency
       
    edm::FileInPath fileInPath("PromptLUAna/PromptAnalyzer/plugins/pyEffic.his");
    ifstream hisfile(fileInPath.fullPath().c_str());

    std::string line;
    double word;

    //...defined at the beginning
    //std::vector<std::vector<double>> efvec;

    if(hisfile)
    {
        while(getline(hisfile, line, '\n'))        
        {
            //create a temporary vector that will contain all the columns
            std::vector<double> tempVec;
            
            std::istringstream ss(line);
            
            //read word by word(or double by double) 
            while(ss >> word)
            {
                //std::cout<<"word:"<<word<<std::endl;
                //add the word to the temporary vector 
                tempVec.push_back(word);
            
            }      
            
            //now all the words from the current line has been added to the temporary vector 
            efvec.emplace_back(tempVec);
        }    
    }
    else {
           std::cout << "Couldn't open file\n";
    }
    hisfile.close();

    /*
    std::cout<<" ############### "<<std::endl;
 
    for (int i = 0; i < 20800; i++)
      {
      for (int j = 0; j < 5; j++)
	{
	  std::cout<< efvec[i][j] << " ";
        }
      std::cout<<std::endl;
      }
    */



   // Armenteros-Podolanski analysis                                   ...Luiz

   ////edm::FileInPath apcutInPath("PromptLUAna/PromptAnalyzer/python/apcut.root");
   //std::cout << " apfile  = " << apcutInPath << std::endl;
   ////TFile g(apcutInPath.fullPath().c_str());
    
    ////TFile* ff = TFile::Open("/afs/cern.ch/user/l/lregisem/CMSSW_10_6_18/src/PromptLUAna/PromptAnalyzer/python/apcut.root");

   //TCutG* apcut = (TCutG*)gROOT->FindObject("apcut");
   //TCutG* apcut2 = (TCutG*)gROOT->FindObject("apcut2");

   //std::cout << " apfilefull  = " << apcutInPath.fullPath().c_str() << std::endl;


   
  //...Luiz
  edm::Service<TFileService> fs;
  
  int nbins_eta = 80;
  int nbins_pt = 100;
  int nbins_phi = 64;


  // TTree memmory 1 TB ...Luiz  
  //TTree::SetMaxTreeSize( 1000000000000LL );

  
  histosTH1F["hrun"] = fs->make<TH1F>("hrun","run#",215,319100,319315);


  histosTH1F["hpt"] = fs->make<TH1F>("hpt","p_{T}",nbins_pt,0,5);
  histosTH1F["heta"] = fs->make<TH1F>("heta","#eta",nbins_eta,-4,4);
  histosTH1F["hphi"] = fs->make<TH1F>("hphi","#varphi",nbins_phi,-3.2,3.2);
  histosTH1F["halgo"] = fs->make<TH1F>("halgo","Algo",15,0,15.);
  histosTH1F["hnhits"] = fs->make<TH1F>("hnhits","nhits pix+strip",40,0,40.);

  histosTH1F["hlooper"] = fs->make<TH1F>("hlooper","isLooper",5,0,5);  
  histosTH1F["hchi2"] = fs->make<TH1F>("hchi2","normalized #chi^{2}",1050,-50,1000.);   
  histosTH1F["hd0"] = fs->make<TH1F>("hd0","d0",    2000,-10,10.);
  histosTH1F["hdz"] = fs->make<TH1F>("hdz","dz",    500,-100,100.);

  histosTH1F["hd0BS"] = fs->make<TH1F>("hd0BS","d0",2000,-10,10.);
  histosTH1F["hdzBS"] = fs->make<TH1F>("hdzBS","dz",500,-100,100.);
  
  histosTH1F["hntrk0"] = fs->make<TH1F>("hntrk0","Ntrk",150,0,150);
  histosTH1F["hntrk"] = fs->make<TH1F>("hntrk","Ntrk for nPixelHits>0",150,0,150);
  
  histosTH1F["hnvtx"] = fs->make<TH1F>("hnvtx","Nvtx",10,0,10);
  histosTH1F["hvtxx"] = fs->make<TH1F>("hvtxx","X vtx",1000,-1.,1.);
  histosTH1F["hvtxy"] = fs->make<TH1F>("hvtxy","Y vtx",1000,-1.,1.);
  histosTH1F["hvtxz"] = fs->make<TH1F>("hvtxz","Z vtx",300,-30.,30.);
  histosTH1F["hvtxchi2"] = fs->make<TH1F>("hvtxchi2","chi2 vtx",1100,-100.,1000.);

  histosTH1F["hetanela"] = fs->make<TH1F>("hetanela","#eta",nbins_eta,-4,4);
  histosTH1F["hetanela4"] = fs->make<TH1F>("hetanela4","#eta",nbins_eta,-4,4);
  histosTH1F["hetanela40"] = fs->make<TH1F>("hetanela40","#eta",nbins_eta,-4,4);
  histosTH1F["hetanela44"] = fs->make<TH1F>("hetanela44","#eta",nbins_eta,-4,4);
  histosTH1F["hetanela440"] = fs->make<TH1F>("hetanela440","#eta",nbins_eta,-4,4);

  //---------------------------
  // CTPPS
  
  histosTH1F["hnconf"] = fs->make<TH1F>("hnconf", "Number of configurations (TB or BT or TT or BB)" , 5, 0., 5.);
  
  vector<string> labRP;
  labRP.push_back("TB"); labRP.push_back("BT"); labRP.push_back("TT"); labRP.push_back("BB");
  
  histosTH1F["hconf"]  = fs->make<TH1F>("hconf"," ",labRP.size(),0,labRP.size());
  for(size_t k = 0; k < labRP.size(); ++k){histosTH1F["hconf"]->GetXaxis()->SetBinLabel((k+1),labRP[k].c_str()); }  
  
  //--------------------------------
  // PPS
  
  histosTH1F["hthyEla"] = fs->make<TH1F>("hthyEla"  ,"#theta_{Y}^{L}+#theta_{Y}^{R}", 2000 , -0.0004 , 0.0004);
  histosTH1F["hthxEla"] = fs->make<TH1F>("hthxEla"  ,"#theta_{X}^{L}+#theta_{X}^{R}", 2000 , -0.0004 , 0.0004);

  histosTH1F["hthyEla2"] = fs->make<TH1F>("hthyEla2"  ,"#theta_{Y}^{L}+#theta_{Y}^{R}", 2000 , -0.0004 , 0.0004);
  histosTH1F["hthxEla2"] = fs->make<TH1F>("hthxEla2"  ,"#theta_{X}^{L}+#theta_{X}^{R}", 2000 , -0.0004 , 0.0004);
  
  histosTH2F["hthx2DIM"] = fs->make<TH2F>("hthx2DIM"  , "#theta_{X}^{R} vs #theta_{X}^{L}" ,400,-0.0004,0.0004,400,-0.0004,0.0004);  
  histosTH2F["hthy2DIM"] = fs->make<TH2F>("hthy2DIM"  , "#theta_{Y}^{R} vs #theta_{Y}^{L}" ,400,-0.0004,0.0004,400,-0.0004,0.0004);  
  histosTH2F["hthythx2DIM"] = fs->make<TH2F>("hthythx2DIM"  , "#theta_{Y} vs #theta_{X}" ,800,-0.0004,0.0004,800,-0.0004,0.0004);
 
  // CT matching 

  histosTH1F["hdpy2trk"] = fs->make<TH1F>("hdpy2trk","2trk, p^{CMS}_{Y}+p^{TOTEM}_{Y}",500,-0.5,0.5);  
  histosTH1F["hdpx2trk"] = fs->make<TH1F>("hdpx2trk","2trk, p^{CMS}_{X}+p^{TOTEM}_{X}",500,-0.5,0.5);

  histosTH1F["hdpy2trkB"] = fs->make<TH1F>("hdpy2trkB","2trk, p^{CMS}_{Y}+p^{TOTEM}_{Y}",500,-0.5,0.5);
  histosTH1F["hdpx2trkB"] = fs->make<TH1F>("hdpx2trkB","2trk, p^{CMS}_{X}+p^{TOTEM}_{X}",500,-0.5,0.5);

  histosTH2F["h2DIMdpy2trk"] = fs->make<TH2F>("h2DIMdpy2trk","2trk, p^{TOTEM}_{Y} vs p^{CMS}_{Y}",200,-2.,2.,200,-2.,2.);  
  histosTH2F["h2DIMdpx2trk"] = fs->make<TH2F>("h2DIMdpx2trk","p^{TOTEM}_{X} vs p^{CMS}_{X}",200,-2.,2.,200,-2.,2.);  

  histosTH1F["hdpy4trk"] = fs->make<TH1F>("hdpy4trk","p^{CMS}_{Y}+p^{TOTEM}_{Y}",500,-0.5,0.5);
  histosTH1F["hdpx4trk"] = fs->make<TH1F>("hdpx4trk","p^{CMS}_{X}+p^{TOTEM}_{X}",500,-0.5,0.5);

  histosTH1F["hdpy4trkB"] = fs->make<TH1F>("hdpy4trkB","p^{CMS}_{Y}+p^{TOTEM}_{Y}",500,-0.5,0.5);
  histosTH1F["hdpx4trkB"] = fs->make<TH1F>("hdpx4trkB","p^{CMS}_{X}+p^{TOTEM}_{X}",500,-0.5,0.5);
  
  histosTH2F["h2DIMdpy4trk"] = fs->make<TH2F>("h2DIMdpy4trk","4trk, p^{TOTEM}_{Y} vs p^{CMS}_{Y}",200,-2.,2.,200,-2.,2.);  
  histosTH2F["h2DIMdpx4trk"] = fs->make<TH2F>("h2DIMdpx4trk","4trk, p^{TOTEM}_{X} vs p^{CMS}_{X}",200,-2.,2.,200,-2.,2.);  

  //-------------------------------------
  // LS validation

  histosTH1F["hLS104"] = fs->make<TH1F>("hLS104","LS 319104",1850,0.,1850.);
  histosTH1F["hLS124"] = fs->make<TH1F>("hLS124","LS 319124",1850,0.,1850.);
  histosTH1F["hLS125"] = fs->make<TH1F>("hLS125","LS 319125",1850,0.,1850.);
  histosTH1F["hLS159"] = fs->make<TH1F>("hLS159","LS 319159",1850,0.,1850.);
  histosTH1F["hLS174"] = fs->make<TH1F>("hLS174","LS 319174",1850,0.,1850.);
  histosTH1F["hLS175"] = fs->make<TH1F>("hLS175","LS 319175",1850,0.,1850.);
  histosTH1F["hLS176"] = fs->make<TH1F>("hLS176","LS 319176",1850,0.,1850.);
  histosTH1F["hLS177"] = fs->make<TH1F>("hLS177","LS 319177",1850,0.,1850.);
  histosTH1F["hLS190"] = fs->make<TH1F>("hLS190","LS 319190",1850,0.,1850.);

  histosTH1F["hLS222"] = fs->make<TH1F>("hLS222","LS 319222",1850,0.,1850.);
  histosTH1F["hLS223"] = fs->make<TH1F>("hLS223","LS 319223",1850,0.,1850.);
  histosTH1F["hLS254"] = fs->make<TH1F>("hLS254","LS 319254",1850,0.,1850.);
  histosTH1F["hLS255"] = fs->make<TH1F>("hLS255","LS 319255",1850,0.,1850.);
  histosTH1F["hLS256"] = fs->make<TH1F>("hLS256","LS 319256",1850,0.,1850.);
  histosTH1F["hLS262"] = fs->make<TH1F>("hLS262","LS 319262",1850,0.,1850.);
  histosTH1F["hLS263"] = fs->make<TH1F>("hLS263","LS 319263",1850,0.,1850.);
  histosTH1F["hLS264"] = fs->make<TH1F>("hLS264","LS 319264",1850,0.,1850.);
  histosTH1F["hLS265"] = fs->make<TH1F>("hLS265","LS 319265",1850,0.,1850.);
  histosTH1F["hLS266"] = fs->make<TH1F>("hLS266","LS 319266",1850,0.,1850.);
  histosTH1F["hLS267"] = fs->make<TH1F>("hLS267","LS 319267",1850,0.,1850.);
  histosTH1F["hLS268"] = fs->make<TH1F>("hLS268","LS 319268",1850,0.,1850.);
  histosTH1F["hLS270"] = fs->make<TH1F>("hLS270","LS 319270",1850,0.,1850.);

  histosTH1F["hLS300"] = fs->make<TH1F>("hLS300","LS 319300",1850,0.,1850.);
  histosTH1F["hLS311"] = fs->make<TH1F>("hLS311","LS 319311",1850,0.,1850.);

  //-------------------------------------
  // Mass spectra
  
  //int massbins=1000;
  //
  //histosTH1F["hm2"]    = fs->make<TH1F>("hm2", "M_{#pi^{+}#pi^{-}} (GeV)", massbins,0.,10.);
  //histosTH1F["hm4"]    = fs->make<TH1F>("hm4", "M_{#pi^{+}#pi^{+}#pi^{-}#pi^{-}} (GeV)", massbins,0.,10.);
  //...Luiz
  histosTH1F["hm2"]    = fs->make<TH1F>("hm2", "M_{#pi^{+}#pi^{-}} (GeV)",1000,0.,10.);
  histosTH1F["hm4"]    = fs->make<TH1F>("hm4", "M_{#pi^{+}#pi^{+}#pi^{-}#pi^{-}} (GeV)",1000,0.,10.);
 
  histosTH1F["hpt2"] = fs->make<TH1F>("hpt2","p_{T}",40,0.,2.);
  histosTH1F["heta2"]= fs->make<TH1F>("heta2","#eta",50,-5.,5.);
  
  histosTH1F["hpt4"] = fs->make<TH1F>("hpt4","p_{T}",40,0.,2.);
  histosTH1F["heta4"]= fs->make<TH1F>("heta4","#eta",50,-5.,5.);
  
  //-----
  
  histosTH1F["hdphi2"] = fs->make<TH1F>("hdphi2","#Delta#varphi_{LR}",320,0,TMath::Pi());
  histosTH1F["hdphi4"] = fs->make<TH1F>("hdphi4","#Delta#varphi_{LR}",320,0,TMath::Pi());
  
  std::cout<<"booked all."<<std::endl;

  //-------------------------------------------------------------------------------
  //...Luiz
  //...my histograms

  histosTH1F["heta4pi"]= fs->make<TH1F>("heta4pi","#eta",2000,-5.,5.);

  //...10-micron resolution ...Luiz          4000 --> 2000 is there a limit here ?
  //**histosTH1F["hpv2dxy"] = fs->make<TH1F>("hpv2dxy","pv2dxy transverse impact parameter w.r.t. pv",2000,-1.,1.);
  //**histosTH1F["hpv2dz"] = fs->make<TH1F>("hpv2dz","pv2dz longitudinal impact parameter w.r.t. pv",2000,-1.,1.);

  histosTH1F["hntrk4q0"] = fs->make<TH1F>("hntrk4q0","Ntrk for nPixelHits>0 Q=0",150,0,150);

  histosTH1F["hntrkvtx"] = fs->make<TH1F>("hntrkvtx","Ntrkvtx",150,0,150);
  histosTH1F["hntrkvtxU"] = fs->make<TH1F>("hntrkvtxU","NtrkvtxU",150,0,150);
  histosTH1F["hntrkvtx0"] = fs->make<TH1F>("hntrkvtx0","Ntrkvtx0",150,0,150);
  histosTH1F["hntrkvtx2"] = fs->make<TH1F>("hntrkvtx2","Ntrkvtx2",150,0,150);
  histosTH1F["hntrkvtx3"] = fs->make<TH1F>("hntrkvtx3","Ntrkvtx3",150,0,150);
  histosTH1F["hntrkvtx4"] = fs->make<TH1F>("hntrkvtx4","Ntrkvtx4",150,0,150);
  //histosTH1F["hntrkntrkvtx2"] = fs->make<TH1F>("hntrkntrkvtx2","Ntrk for Ntrkvtx=2",150,0,150);
  //histosTH1F["hntrk2ntrkvtx"] = fs->make<TH1F>("hntrk2ntrkvtx","Ntrkvtx for Ntrk=2",150,0,150);
  
  //histosTH2F["hntrkntrkvtx"] = new TH2F("hntrkntrkvtx","Ntrk vs Ntrkvtx",150,0,150,150,0,150);
  
  histosTH1F["hvtx"] = fs->make<TH1F>("hvtx","vtx.isFake()",2,0,2);
  //...Luiz
  histosTH1F["hvtx2"] = fs->make<TH1F>("hvtx2","vtx.isFake() 4 tracks",2,0,2);
  histosTH1F["hvtx3"] = fs->make<TH1F>("hvtx3","vtx.isFake() 4 tracks both |#eta|<2.5 and OS",2,0,2);

  //...Kshorts
  histosTH1F["hnks"] = fs->make<TH1F>("hnks","N Kshorts",10,0,10);
  histosTH2F["hntrknks"] = fs->make<TH2F>("hntrknks","# of K0s Vees vs # of Tracks",150,0,150,10,0,10);
  histosTH2F["hnvtxnks"] = fs->make<TH2F>("hnvtxnks","# of K0s Vees vs # of Vertices",150,0,150,10,0,10);
  //...all4
  histosTH2F["hntrknksall4"] = fs->make<TH2F>("hntrknksall4","# of K0s Vees vs # of Tracks all4",150,0,150,10,0,10);
  histosTH2F["hnvtxnksall4"] = fs->make<TH2F>("hnvtxnksall4","# of K0s Vees vs # of Vertices all4",150,0,150,10,0,10);
  //...Q=0
  histosTH2F["hntrknksq0"] = fs->make<TH2F>("hntrknksq0","# of K0s Vees vs # of Tracks Q=0",150,0,150,10,0,10);
  histosTH2F["hnvtxnksq0"] = fs->make<TH2F>("hnvtxnksq0","# of K0s Vees vs # of Vertices Q=0",150,0,150,10,0,10);
  //
  histosTH2F["hntrknvtx"] = fs->make<TH2F>("hntrknvtx","# of Vertices vs # of Tracks",150,0,150,150,0,150);
  histosTH1F["hksvertexx"] = fs->make<TH1F>("hksvertexx","K0s X vertex",1200,-30.,30.);
  histosTH1F["hksvertexy"] = fs->make<TH1F>("hksvertexy","K0s Y vertex",1200,-30.,30.);
  histosTH1F["hksvertexz"] = fs->make<TH1F>("hksvertexz","K0s Z vertex",1200,-50.,50.);
  //
  histosTH1F["hksradius"] = fs->make<TH1F>("hksradius","K0s vertex radius",20000,0.,100.);
  histosTH1F["hksradiusvtx"] = fs->make<TH1F>("hksradiusvtx","K0s vertex radius vtx",20000,0.,100.);
  histosTH1F["hksradiusv1"] = fs->make<TH1F>("hksradiusv1","K0s vertex radius v1",20000,0.,100.);
  histosTH1F["hksradiusv2"] = fs->make<TH1F>("hksradiusv2","K0s vertex radius v2",20000,0.,100.);
  histosTH1F["hksradiusv2vtx"] = fs->make<TH1F>("hksradiusv2vtx","K0s vertex radius v2 vtx",20000,0.,100.);
  histosTH1F["hks3Dradius"] = fs->make<TH1F>("hks3Dradius","K0s vertex radius 3D",20000,0.,100.);
  histosTH1F["hks3Dradiusvtx"] = fs->make<TH1F>("hks3Dradiusvtx","K0s vertex radius 3D vtx",20000,0.,100.);
  histosTH1F["hks3Dradiusv1"] = fs->make<TH1F>("hks3Dradiusv1","K0s vertex radius 3D v1",20000,0.,100.);
  histosTH1F["hks3Dradiusv2"] = fs->make<TH1F>("hks3Dradiusv2","K0s vertex radius 3D v2",20000,0.,100.);
  histosTH1F["hks3Dradiusv2vtx"] = fs->make<TH1F>("hks3Dradiusv2vtx","K0s vertex radius 3D v2 vtx",20000,0.,100.);
  //
  // lifetime
  histosTH1F["hkslifetime"] = fs->make<TH1F>("hkslifetime","K0s lifetime",20000,0.,100.);
  histosTH1F["hkslifetimev1"] = fs->make<TH1F>("hkslifetimev1","K0s lifetime v1",20000,0.,100.);
  histosTH1F["hkslifetimev2"] = fs->make<TH1F>("hkslifetimev2","K0s lifetime v2",20000,0.,100.);
  // vtx
  histosTH1F["hkslifetimevtx"] = fs->make<TH1F>("hkslifetimevtx","K0s lifetime vtx",20000,0.,100.);
  //histosTH1F["hkslifetimev1"] = fs->make<TH1F>("hkslifetimev1","K0s lifetime v1",20000,0.,100.);
  histosTH1F["hkslifetimev2vtx"] = fs->make<TH1F>("hkslifetimev2vtx","K0s lifetime v2 vtx",20000,0.,100.);
  // z
  histosTH1F["hkslifetimez"] = fs->make<TH1F>("hkslifetimez","K0s lifetime z",20000,0.,100.);
  histosTH1F["hkslifetimev2z"] = fs->make<TH1F>("hkslifetimev2z","K0s lifetime v2 z",20000,0.,100.);
  // z vtx
  histosTH1F["hkslifetimezvtx"] = fs->make<TH1F>("hkslifetimezvtx","K0s lifetime z vtx",20000,0.,100.);
  histosTH1F["hkslifetimev2zvtx"] = fs->make<TH1F>("hkslifetimev2zvtx","K0s lifetime v2 z vtx",20000,0.,100.);
  // 3D cor
  histosTH1F["hks3Dlifetime"] = fs->make<TH1F>("hks3Dlifetime","K0s lifetime 3D",20000,0.,100.);
  histosTH1F["hks3Dlifetimecor"] = fs->make<TH1F>("hks3Dlifetimecor","K0s lifetime 3D cor",20000,0.,100.);
  // 3D vtx cor
  histosTH1F["hks3Dlifetimevtx"] = fs->make<TH1F>("hks3Dlifetimevtx","K0s lifetime 3D vtx",20000,0.,100.);
  histosTH1F["hks3Dlifetimevtxcor"] = fs->make<TH1F>("hks3Dlifetimevtxcor","K0s lifetime 3D vtx cor",20000,0.,100.);
  // 3D z
  histosTH1F["hks3Dlifetimez"] = fs->make<TH1F>("hks3Dlifetimez","K0s lifetime 3D z",20000,0.,100.);
  histosTH1F["hks3Dlifetimecorz"] = fs->make<TH1F>("hks3Dlifetimecorz","K0s lifetime 3D z",20000,0.,100.);
  // 3D z vtx cor
  histosTH1F["hks3Dlifetimezvtx"] = fs->make<TH1F>("hks3Dlifetimezvtx","K0s lifetime 3D z vtx",20000,0.,100.);
  histosTH1F["hks3Dlifetimezvtxcor"] = fs->make<TH1F>("hks3Dlifetimezvtxcor","K0s lifetime 3D z vtx cor",20000,0.,100.);
  // 3D v2
  histosTH1F["hks3Dlifetimev1"] = fs->make<TH1F>("hks3Dlifetimev1","K0s lifetime 3D v1",20000,0.,100.);
  histosTH1F["hks3Dlifetimev2"] = fs->make<TH1F>("hks3Dlifetimev2","K0s lifetime 3D v2",20000,0.,100.);
  histosTH1F["hks3Dlifetimev2cor"] = fs->make<TH1F>("hks3Dlifetimev2cor","K0s lifetime 3D v2",20000,0.,100.);
  // 3D v2 vtx cor
  //histosTH1F["hks3Dlifetimev1"] = fs->make<TH1F>("hks3Dlifetimev1","K0s lifetime 3D v1",20000,0.,100.);
  histosTH1F["hks3Dlifetimev2vtx"] = fs->make<TH1F>("hks3Dlifetimev2vtx","K0s lifetime 3D v2 vtx",20000,0.,100.);
  histosTH1F["hks3Dlifetimev2vtxcor"] = fs->make<TH1F>("hks3Dlifetimev2vtxcor","K0s lifetime 3D v2 vtx cor",20000,0.,100.);
  // 3D v2 z
  histosTH1F["hks3Dlifetimev2z"] = fs->make<TH1F>("hks3Dlifetimev2z","K0s lifetime 3D v2 z",20000,0.,100.);
  histosTH1F["hks3Dlifetimev2zcor"] = fs->make<TH1F>("hks3Dlifetimev2zcor","K0s lifetime 3D v2 z cor",20000,0.,100.);
  // 3D v2 z vtx cor
  histosTH1F["hks3Dlifetimev2zvtx"] = fs->make<TH1F>("hks3Dlifetimev2zvtx","K0s lifetime 3D v2 z vtx",20000,0.,100.);
  histosTH1F["hks3Dlifetimev2zvtxcor"] = fs->make<TH1F>("hks3Dlifetimev2zvtxcor","K0s lifetime 3D v2 z vtx cor",20000,0.,100.);
  //
  //...2D
  histosTH2F["h2dimksxy"] = fs->make<TH2F>("h2dimksxy","K0s X vs Y vtx",300,-30.,30.,300,-30.,30.);
  histosTH2F["h2dimksxz"] = fs->make<TH2F>("h2dimksxz","K0s X vs Z vtx",300,-30.,30.,300,-30.,30.);
  histosTH2F["h2dimksyz"] = fs->make<TH2F>("h2dimksyz","K0s Y vs Z vtx",300,-30.,30.,300,-30.,30.);
  //
  histosTH1F["hkspt"] = fs->make<TH1F>("hkspt","K0s pt",100,0.,5.);
  histosTH1F["hkseta"] = fs->make<TH1F>("hkseta","K0s #eta",80,-4.,4.);
  histosTH1F["hksphi"] = fs->make<TH1F>("hksphi","K0s #varphi",64,-3.2,3.2);
  histosTH1F["hksmass"] = fs->make<TH1F>("hksmass","K0s mass",2500,0.,5.);
  histosTH1F["hksmasscor"] = fs->make<TH1F>("hksmasscor","K0s mass",2500,0.,5.);
  //
  histosTH1F["hksmassv1"] = fs->make<TH1F>("hksmassv1","K0s mass 1 vertex",2500,0.,5.);
  histosTH1F["hksmassv2"] = fs->make<TH1F>("hksmassv2","K0sK0s mass 2 vertices",2500,0.,5.);
  histosTH1F["hksmassv2cor"] = fs->make<TH1F>("hksmassv2cor","K0sK0s mass 2 vertices efficiency correction",2500,0.,5.);
  //**histosTH1F["hksmassv3"] = fs->make<TH1F>("hksmassv3","K0s mass 3 vertices",250,0.,5.);
  //
  //**histosTH1F["hksrad"] = fs->make<TH1F>("hksrad","K0s radius",800,0,40.);
  //**histosTH1F["hksrad4"] = fs->make<TH1F>("hksrad4","K0s radius4",800,0,40.);
  //**histosTH1F["hksrad0"] = fs->make<TH1F>("hksrad0","K0s radius0",800,0,40.);
  //
  histosTH1F["hrimpac1"] = fs->make<TH1F>("hrimpac1","3D impact parameter 1",10000,0,10.);
  histosTH1F["hrimpac2"] = fs->make<TH1F>("hrimpac2","3D impact parameter 2",10000,0,10.);
  histosTH1F["hrimpac3"] = fs->make<TH1F>("hrimpac3","3D impact parameter 3",10000,0,10.);
  histosTH1F["hrimpac4"] = fs->make<TH1F>("hrimpac4","3D impact parameter 4",10000,0,10.);
  //
  //**histosTH1F["hntag1"] = fs->make<TH1F>("hntag1","test ntag1",10,0,10);
  //**histosTH1F["hntag2"] = fs->make<TH1F>("hntag2","test ntag2",10,0,10);
  //**histosTH1F["hntag3"] = fs->make<TH1F>("hntag3","test ntag3",10,0,10);
  
  //...Lambdas
  histosTH1F["hnlam"] = fs->make<TH1F>("hnlam","N Lambdas",10,0,10);
  histosTH2F["hntrknlam"] = fs->make<TH2F>("hntrknlam","# of #Lambda Vees vs # of Tracks",150,0,150,10,0,10);
  histosTH2F["hnvtxnlam"] = fs->make<TH2F>("hnvtxnlam","# of #Lambda Vees vs # of Vertices",150,0,150,10,0,10);
  histosTH1F["hlamvertexx"] = fs->make<TH1F>("hlamvertexx","#Lambda X vertex",120,-30.,30.);
  histosTH1F["hlamvertexy"] = fs->make<TH1F>("hlamvertexy","#Lambda Y vertex",120,-30.,30.);
  histosTH1F["hlamvertexz"] = fs->make<TH1F>("hlamvertexz","#Lambda Z vertex",120,-30.,30.);
  histosTH1F["hlamradius"] = fs->make<TH1F>("hlamradius","#Lambda vertex radius",20000,0.,100.);
  histosTH1F["hlam3Dradius"] = fs->make<TH1F>("hlam3Dradius","#Lambda vertex 3D radius",20000,0.,100.);
  histosTH1F["hlamlifetime"] = fs->make<TH1F>("hlamlifetime","#Lambda vertex lifetime",20000,0.,100.);
  histosTH1F["hlam3Dlifetime"] = fs->make<TH1F>("hlam3Dlifetime","#Lambda vertex 3D lifetime",20000,0.,100.);
  //...2D
  histosTH2F["h2dimlamxy"] = fs->make<TH2F>("h2dimlamxy","#Lambda X vs Y vtx",300,-30.,30.,300,-30.,30.);
  histosTH2F["h2dimlamxz"] = fs->make<TH2F>("h2dimlamxz","#Lambda X vs Z vtx",300,-30.,30.,300,-30.,30.);
  histosTH2F["h2dimlamyz"] = fs->make<TH2F>("h2dimlamyz","#Lambda Y vs Z vtx",300,-30.,30.,300,-30.,30.);
  //
  histosTH1F["hlampt"] = fs->make<TH1F>("hlampt","#Lambda pt",100,0.,5.);
  histosTH1F["hlameta"] = fs->make<TH1F>("hlameta","#Lambda #eta",80,-4.,4.);
  histosTH1F["hlamphi"] = fs->make<TH1F>("hlamphi","#Lambda #varphi",64,-3.2,3.2);
  histosTH1F["hlammass"] = fs->make<TH1F>("hlammass","#Lambda mass",250,0.,5.);

  //...already defined
  //histosTH1F["hvtxx"] = fs->make<TH1F>("hvtxx","X vtx",1000,-1.,1.);
  //histosTH1F["hvtxy"] = fs->make<TH1F>("hvtxy","Y vtx",1000,-1.,1.);
  //**histosTH1F["hvtxx4"] = fs->make<TH1F>("hvtxx4","X vtx",1000,-1.,1.);
  //**histosTH1F["hvtxy4"] = fs->make<TH1F>("hvtxy4","Y vtx",1000,-1.,1.);
  //histosTH1F["hvtxz"] = fs->make<TH1F>("hvtxz","Z vtx",300,-30.,30.);
  //**histosTH1F["hvtxz4"] = fs->make<TH1F>("hvtxz4","Z vtx",300,-30.,30.);

  //...pair position
  //**histosTH1F["hxpi1pi2"] = fs->make<TH1F>("hxpi1pi2","X pi1pi2",10000,-100.,100.);
  //**histosTH1F["hypi1pi2"] = fs->make<TH1F>("hypi1pi2","Y pi1pi2",10000,-100.,100.);
  //**histosTH1F["hr12"] = fs->make<TH1F>("hr12","R pi1pi2",10000,0.,200.);

  //...Luiz
  //***histosTH2F["hvtx2dimxy"] = fs->make<TH2F>("hvtx2dimxy","X vs Y vtx",1000,-1.,1.,1000,-1.,1.);
  //***histosTH2F["hvtx2dimxz"] = fs->make<TH2F>("hvtx2dimxz","X vs Z vtx",1000,-1.,1.,300,-3.,3.);
  //***histosTH2F["hvtx2dimyz"] = fs->make<TH2F>("hvtx2dimyz","Y vs Z vtx",1000,-1.,1.,300,-3.,3.);

  int massbins=250;
  
  //**histosTH1F["hm"] = fs->make<TH1F>("hm","M_{4#pi} ",massbins,0,5.);
  //...Luiz
  //**histosTH1F["hmxicut"] = fs->make<TH1F>("hmxicut","M_{4#pi} ",massbins,0,5.);
  
  //**histosTH1F["hm4rec"] = fs->make<TH1F>("hm4rec","M_{4#pi} ",massbins,0,5.);
  //**histosTH1F["hm4recbis"] = fs->make<TH1F>("hm4recbis","M_{4#pi}",2*massbins,0,5.);

  //...Luiz
  //**histosTH1F["hm4recPPPP"] = fs->make<TH1F>("hm4recPPPP","M_{4#pi} ",massbins,0,5.);
  //**histosTH1F["hm4recKKKK"] = fs->make<TH1F>("hm4recKKKK","M_{4K} ",massbins,0,5.);
 
  //...transverse impact parameter histograms d0 (dxy)
  histosTH1F["hd01"] = fs->make<TH1F>("hd01","d01 ntrk=4 OS",5000,-5.,5.);
  histosTH1F["hd02"] = fs->make<TH1F>("hd02","d02 ntrk=4 OS",5000,-5.,5.);
  histosTH1F["hd03"] = fs->make<TH1F>("hd03","d03 ntrk=4 OS",5000,-5.,5.);
  histosTH1F["hd04"] = fs->make<TH1F>("hd04","d04 ntrk=4 OS",5000,-5.,5.);
  //...longitudinal impact parameter histograms dz
  histosTH1F["hdz1"] = fs->make<TH1F>("hdz1","dz1 ntrk=4 OS",4000,-20.,20.);
  histosTH1F["hdz2"] = fs->make<TH1F>("hdz2","dz2 ntrk=4 OS",4000,-20.,20.);
  histosTH1F["hdz3"] = fs->make<TH1F>("hdz3","dz3 ntrk=4 OS",4000,-20.,20.);
  histosTH1F["hdz4"] = fs->make<TH1F>("hdz4","dz4 ntrk=4 OS",4000,-20.,20.);
  //..transverse impact parameter histograms vtxdxy
  histosTH1F["hvtxdxy1"] = fs->make<TH1F>("hvtxdxy1","vtxdxy1 ntrk=4 OS",5000,-5.,5.);
  histosTH1F["hvtxdxy2"] = fs->make<TH1F>("hvtxdxy2","vtxdxy2 ntrk=4 OS",5000,-5.,5.);
  histosTH1F["hvtxdxy3"] = fs->make<TH1F>("hvtxdxy3","vtxdxy3 ntrk=4 OS",5000,-5.,5.);
  histosTH1F["hvtxdxy4"] = fs->make<TH1F>("hvtxdxy4","vtxdxy4 ntrk=4 OS",5000,-5.,5.);
  //...longitudinal impact parameter histograms vtxdz
  histosTH1F["hvtxdz1"] = fs->make<TH1F>("hvtxdz1","vtxdz1 ntrk=4 OS",5000,-5.,5.);
  histosTH1F["hvtxdz2"] = fs->make<TH1F>("hvtxdz2","vtxdz2 ntrk=4 OS",5000,-5.,5.);
  histosTH1F["hvtxdz3"] = fs->make<TH1F>("hvtxdz3","vtxdz3 ntrk=4 OS",5000,-5.,5.);
  histosTH1F["hvtxdz4"] = fs->make<TH1F>("hvtxdz4","vtxdz4 ntrk=4 OS",5000,-5.,5.);

  //...OS-SS
  histosTH1F["hm4recOS"] = fs->make<TH1F>("hm4recOS","M_{4#pi} OS",massbins,0,5.);
  histosTH1F["hm4recOS2"] = fs->make<TH1F>("hm4recOS2","M_{4#pi} OS",2.0*massbins,0,10.);
  //
  //**histosTH1F["hm4recSS"] = fs->make<TH1F>("hm4recSS","M_{4#pi} SS",massbins,0,5.);
  //
  //**histosTH1F["hm4recOS_diag"] = fs->make<TH1F>("hm4recOS_diag","M_{4#pi} TB/BT OS",massbins,0,5.);
  //**histosTH1F["hm4recSS_diag"] = fs->make<TH1F>("hm4recSS_diag","M_{4#pi} TB/BT SS",massbins,0,5.);
  //**histosTH1F["hm4recOS_ttbb"] = fs->make<TH1F>("hm4recOS_ttbb","M_{4#pi} TT/BB OS",massbins,0,5.);
  //**histosTH1F["hm4recSS_ttbb"] = fs->make<TH1F>("hm4recSS_ttbb","M_{4#pi} TT/BB SS",massbins,0,5.);

  //...2OS-2SS
  histosTH1F["hm4rec2OS"] = fs->make<TH1F>("hm4rec2OS","M_{4#pi} OS",massbins,0,5.);
  histosTH1F["hm4rec2OScor"] = fs->make<TH1F>("hm4rec2OScor","M_{4#pi} OS",massbins,0,5.);
  histosTH1F["hm4rec2OSk"] = fs->make<TH1F>("hm4rec2OSk","M_{4K} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OSrejKp"] = fs->make<TH1F>("hm4rec2OSrejKp","M_{4#pi} OS rejecting at least one K or p",massbins,0,5.);
  //**histosTH1F["hm4rec2OS2rejKp"] = fs->make<TH1F>("hm4rec2OS2rejKp","M_{4#pi} OS nvtx=1 rejecting at least one K or p",massbins,0,5.);
  //**histosTH1F["hm4rec2OSrejKpu"] = fs->make<TH1F>("hm4rec2OSrejKpu","M_{4#pi} OS rejecting at least one K or p or unknown",massbins,0,5.);
  //**histosTH1F["hm4rec2OS2rejKpu"] = fs->make<TH1F>("hm4rec2OS2rejKpu","M_{4#pi} OS nvtx=1 rejecting at least one K or p or unknown",massbins,0,5.);
  //
  //...mass selection
  //
  //**histosTH1F["hm4rec2OSm1234"] = fs->make<TH1F>("hm4rec2OSm1234","M_{4#pi} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OSm1324"] = fs->make<TH1F>("hm4rec2OSm1324","M_{4#pi} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OSm1423"] = fs->make<TH1F>("hm4rec2OSm1423","M_{4#pi} OS",massbins,0,5.);
  //
  //**histosTH1F["hm4rec2OSm123400"] = fs->make<TH1F>("hm4rec2OSm123400","M_{4#pi} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OSm132400"] = fs->make<TH1F>("hm4rec2OSm132400","M_{4#pi} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OSm142300"] = fs->make<TH1F>("hm4rec2OSm142300","M_{4#pi} OS",massbins,0,5.);
  //
  //**histosTH1F["hm4rec2OSm1234000"] = fs->make<TH1F>("hm4rec2OSm1234000","M_{4#pi} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OSm1324000"] = fs->make<TH1F>("hm4rec2OSm1324000","M_{4#pi} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OSm1423000"] = fs->make<TH1F>("hm4rec2OSm1423000","M_{4#pi} OS",massbins,0,5.);
  //
  //**histosTH1F["hm4rec2OSm123404"] = fs->make<TH1F>("hm4rec2OSm123404","M_{4#pi} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OSm132404"] = fs->make<TH1F>("hm4rec2OSm132404","M_{4#pi} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OSm142304"] = fs->make<TH1F>("hm4rec2OSm142304","M_{4#pi} OS",massbins,0,5.);
  //
  histosTH1F["hm4rec2OSm123405"] = fs->make<TH1F>("hm4rec2OSm123405","M_{4#pi} OS",massbins,0,5.);
  histosTH1F["hm4rec2OSm132405"] = fs->make<TH1F>("hm4rec2OSm132405","M_{4#pi} OS",massbins,0,5.);
  histosTH1F["hm4rec2OSm142305"] = fs->make<TH1F>("hm4rec2OSm142305","M_{4#pi} OS",massbins,0,5.);
  //
  histosTH1F["hm4rec2OSm1234052"] = fs->make<TH1F>("hm4rec2OSm1234052","M_{4#pi} OS range",massbins,0,5.);
  histosTH1F["hm4rec2OSm1324052"] = fs->make<TH1F>("hm4rec2OSm1324052","M_{4#pi} OS range",massbins,0,5.);
  histosTH1F["hm4rec2OSm1423052"] = fs->make<TH1F>("hm4rec2OSm1423052","M_{4#pi} OS range",massbins,0,5.);
  //
  histosTH1F["hm4rec2OSm123405sig"] = fs->make<TH1F>("hm4rec2OSm123405sig","M_{4#pi} OS sigma",massbins,0,5.);
  histosTH1F["hm4rec2OSm132405sig"] = fs->make<TH1F>("hm4rec2OSm132405sig","M_{4#pi} OS sigma",massbins,0,5.);
  histosTH1F["hm4rec2OSm142305sig"] = fs->make<TH1F>("hm4rec2OSm142305sig","M_{4#pi} OS sigma",massbins,0,5.);


  // final studies : mass windows
  //K0
  histosTH1F["hm4rec2OSm123405win10"] = fs->make<TH1F>("hm4rec2OSm123405win10","M_{4#pi} OS study 9 MeV",massbins,0,5.);
  histosTH1F["hm4rec2OSm132405win10"] = fs->make<TH1F>("hm4rec2OSm132405win10","M_{4#pi} OS study 9 MeV",massbins,0,5.);
  histosTH1F["hm4rec2OSm142305win10"] = fs->make<TH1F>("hm4rec2OSm142305win10","M_{4#pi} OS study 9 MeV",massbins,0,5.);
  //
  histosTH1F["hm4rec2OSm123405win20"] = fs->make<TH1F>("hm4rec2OSm123405win20","M_{4#pi} OS study 18 MeV",massbins,0,5.);
  histosTH1F["hm4rec2OSm132405win20"] = fs->make<TH1F>("hm4rec2OSm132405win20","M_{4#pi} OS study 18 MeV",massbins,0,5.);
  histosTH1F["hm4rec2OSm142305win20"] = fs->make<TH1F>("hm4rec2OSm142305win20","M_{4#pi} OS study 18 MeV",massbins,0,5.);
  //
  histosTH1F["hm4rec2OSm123405win30"] = fs->make<TH1F>("hm4rec2OSm123405win30","M_{4#pi} OS study 27 MeV",massbins,0,5.);
  histosTH1F["hm4rec2OSm132405win30"] = fs->make<TH1F>("hm4rec2OSm132405win30","M_{4#pi} OS study 27 MeV",massbins,0,5.);
  histosTH1F["hm4rec2OSm142305win30"] = fs->make<TH1F>("hm4rec2OSm142305win30","M_{4#pi} OS study 27 MeV",massbins,0,5.);
  //
  //rho(770)
  histosTH1F["hm4rec2OSr123405win74"] = fs->make<TH1F>("hm4rec2OSr123405win74","M_{4#pi} OS study 74 MeV",massbins,0,5.);  
  histosTH1F["hm4rec2OSr132405win74"] = fs->make<TH1F>("hm4rec2OSr132405win74","M_{4#pi} OS study 74 MeV",massbins,0,5.);  
  histosTH1F["hm4rec2OSr142305win74"] = fs->make<TH1F>("hm4rec2OSr142305win74","M_{4#pi} OS study 74 MeV",massbins,0,5.);  
  //
  histosTH1F["hm4rec2OSr123405win148"] = fs->make<TH1F>("hm4rec2OSr123405win148","M_{4#pi} OS study 148 MeV",massbins,0,5.);  
  histosTH1F["hm4rec2OSr132405win148"] = fs->make<TH1F>("hm4rec2OSr132405win148","M_{4#pi} OS study 148 MeV",massbins,0,5.);  
  histosTH1F["hm4rec2OSr142305win148"] = fs->make<TH1F>("hm4rec2OSr142305win148","M_{4#pi} OS study 148 MeV",massbins,0,5.);  
  //
  histosTH1F["hm4rec2OSr123405winhalf"] = fs->make<TH1F>("hm4rec2OSr123405winhalf","M_{4#pi} OS study 37 MeV",massbins,0,5.);  
  histosTH1F["hm4rec2OSr132405winhalf"] = fs->make<TH1F>("hm4rec2OSr132405winhalf","M_{4#pi} OS study 37 MeV",massbins,0,5.);  
  histosTH1F["hm4rec2OSr142305winhalf"] = fs->make<TH1F>("hm4rec2OSr142305winhalf","M_{4#pi} OS study 37 MeV",massbins,0,5.);  
  //
  histosTH1F["hm4rec2OSr123405winpt"] = fs->make<TH1F>("hm4rec2OSr123405winpt","M_{4#pi} OS study 37 MeV pT<0.8",massbins,0,5.);  
  histosTH1F["hm4rec2OSr132405winpt"] = fs->make<TH1F>("hm4rec2OSr132405winpt","M_{4#pi} OS study 37 MeV pT<0.8",massbins,0,5.);  
  histosTH1F["hm4rec2OSr142305winpt"] = fs->make<TH1F>("hm4rec2OSr142305winpt","M_{4#pi} OS study 37 MeV pT<0.8",massbins,0,5.);  
  //
  //phi(1020)
  histosTH1F["hm4rec2OSk123405win9"] = fs->make<TH1F>("hm4rec2OSk123405win9","M_{4#pi} OS study 9 MeV",massbins,0,5.);  
  histosTH1F["hm4rec2OSk132405win9"] = fs->make<TH1F>("hm4rec2OSk132405win9","M_{4#pi} OS study 9 MeV",massbins,0,5.);  
  histosTH1F["hm4rec2OSk142305win9"] = fs->make<TH1F>("hm4rec2OSk142305win9","M_{4#pi} OS study 9 MeV",massbins,0,5.);  
  //
  histosTH1F["hm4rec2OSk123405win18"] = fs->make<TH1F>("hm4rec2OSk123405win18","M_{4#pi} OS study 18 MeV",massbins,0,5.);  
  histosTH1F["hm4rec2OSk132405win18"] = fs->make<TH1F>("hm4rec2OSk132405win18","M_{4#pi} OS study 18 MeV",massbins,0,5.);  
  histosTH1F["hm4rec2OSk142305win18"] = fs->make<TH1F>("hm4rec2OSk142305win18","M_{4#pi} OS study 18 MeV",massbins,0,5.);  
  //
  //X(1070)
  histosTH1F["hm4rec2OSk123470win8"] = fs->make<TH1F>("hm4rec2OSk123470win8","M_{4#pi} OS study 8 MeV",massbins,0,5.);  
  histosTH1F["hm4rec2OSk132470win8"] = fs->make<TH1F>("hm4rec2OSk132470win8","M_{4#pi} OS study 8 MeV",massbins,0,5.);  
  histosTH1F["hm4rec2OSk142370win8"] = fs->make<TH1F>("hm4rec2OSk142370win8","M_{4#pi} OS study 8 MeV",massbins,0,5.);  
  //
  histosTH1F["hm4rec2OSk123470win17"] = fs->make<TH1F>("hm4rec2OSk123470win17","M_{4#pi} OS study 17 MeV",massbins,0,5.);  
  histosTH1F["hm4rec2OSk132470win17"] = fs->make<TH1F>("hm4rec2OSk132470win17","M_{4#pi} OS study 17 MeV",massbins,0,5.);  
  histosTH1F["hm4rec2OSk142370win17"] = fs->make<TH1F>("hm4rec2OSk142370win17","M_{4#pi} OS study 17 MeV",massbins,0,5.);  
  //
  //X(1235)
  histosTH1F["hm4rec2OSk1234235win74"] = fs->make<TH1F>("hm4rec2OSk1234235win74","M_{4#pi} OS study 74 MeV",massbins,0,5.); 
  histosTH1F["hm4rec2OSk1324235win74"] = fs->make<TH1F>("hm4rec2OSk1324235win74","M_{4#pi} OS study 74 MeV",massbins,0,5.); 
  histosTH1F["hm4rec2OSk1423235win74"] = fs->make<TH1F>("hm4rec2OSk1423235win74","M_{4#pi} OS study 74 MeV",massbins,0,5.); 
  //
  histosTH1F["hm4rec2OSk1234235win148"] = fs->make<TH1F>("hm4rec2OSk1234235win148","M_{4#pi} OS study 148 MeV",massbins,0,5.); 
  histosTH1F["hm4rec2OSk1324235win148"] = fs->make<TH1F>("hm4rec2OSk1324235win148","M_{4#pi} OS study 148 MeV",massbins,0,5.); 
  histosTH1F["hm4rec2OSk1423235win148"] = fs->make<TH1F>("hm4rec2OSk1423235win148","M_{4#pi} OS study 148 MeV",massbins,0,5.); 
  //
  
   //
   // studies
   //
   //...rho study ...M(14)xM(23)
  histosTH1F["hm4rec2OS14rho"] = fs->make<TH1F>("hm4rec2OS14rho","M_{4#pi} OS #rho#rho study 300 MeV",massbins,0,5.);
   //...K0 study ...M(14)xM(23)
  histosTH1F["hm4rec2OS14k0"] = fs->make<TH1F>("hm4rec2OS14k0","M_{4#pi} OS K0K0 study 100 MeV",massbins,0,5.);
  histosTH1F["hm4rec2OS1423k0"] = fs->make<TH1F>("hm4rec2OS1423k0","M_{4#pi} OS 1423 study 100 MeV",massbins,0,5.);
   //...phi study ...M(12)xM(34) !!!
  histosTH1F["hm4rec2OS34phi"] = fs->make<TH1F>("hm4rec2OS34phi","M_{4K} OS #phi#phi study 50 MeV",massbins,0,5.);
   //...phi study ...M(14)xM(23) !!!
  histosTH1F["hm4rec2OS14phi"] = fs->make<TH1F>("hm4rec2OS14phi","M_{4K} OS #phi#phi study 40 MeV",massbins,0,5.);
   //...exo study 1070 ...M(14)xM(23)
  histosTH1F["hm4rec2OS14exo70"] = fs->make<TH1F>("hm4rec2OS14exo70","M_{4K} OS 1070 study B 50 MeV",massbins,0,5.);
  histosTH1F["hm4rec2OS23exo70"] = fs->make<TH1F>("hm4rec2OS23exo70","M_{4K} OS 1070 study A 60 MeV",massbins,0,5.);
   //...exo study 1220 ...M(14)xM(23)
  histosTH1F["hm4rec2OS14exo1220"] = fs->make<TH1F>("hm4rec2OS14exo1220","M_{4K} OS 1220 study 180 MeV",massbins,0,5.);
   //
   //
  histosTH1F["hm4rec2OS1423exo"] = fs->make<TH1F>("hm4rec2OS1423exo","M_{4K} OS at 1230 230 MeV ",massbins,0,5.);
  histosTH1F["hm4rec2OS1423exopid"] = fs->make<TH1F>("hm4rec2OS1423exopid","M_{4K} OS at 1230 230 MeV PID",massbins,0,5.);
   //...exo study 1085 ...M(14)xM(23)
  histosTH1F["hm4rec2OS14exo85"] = fs->make<TH1F>("hm4rec2OS14exo85","M_{4K} OS 1085 study",massbins,0,5.);
  histosTH1F["hm4rec2OS1423exo85"] = fs->make<TH1F>("hm4rec2OS1423exo85","M_{4K} OS at 1085 80 MeV ",massbins,0,5.);
  histosTH1F["hm4rec2OS1423exo85pid"] = fs->make<TH1F>("hm4rec2OS1423exo85pid","M_{4K} OS at 1085 80 MeV MeV PID",massbins,0,5.);
  //
  histosTH1F["hm4rec2OSr123405"] = fs->make<TH1F>("hm4rec2OSr123405","M_{4#pi} OS #rho#rho",massbins,0,5.);
  histosTH1F["hm4rec2OSr132405"] = fs->make<TH1F>("hm4rec2OSr132405","M_{4#pi} OS #rho#rho",massbins,0,5.);
  histosTH1F["hm4rec2OSr142305"] = fs->make<TH1F>("hm4rec2OSr142305","M_{4#pi} OS #rho#rho",massbins,0,5.);
   //
  histosTH1F["hm4rec2OSr1234052"] = fs->make<TH1F>("hm4rec2OSr1234052","M_{4#pi} OS #rho#rho range",massbins,0,5.);
  histosTH1F["hm4rec2OSr1324052"] = fs->make<TH1F>("hm4rec2OSr1324052","M_{4#pi} OS #rho#rho range",massbins,0,5.);
  histosTH1F["hm4rec2OSr1423052"] = fs->make<TH1F>("hm4rec2OSr1423052","M_{4#pi} OS #rho#rho range",massbins,0,5.);
   //
  histosTH1F["hm4rec2OSr123405sig"] = fs->make<TH1F>("hm4rec2OSr123405sig","M_{4#pi} OS #rho#rho sigma",massbins,0,5.);
  histosTH1F["hm4rec2OSr132405sig"] = fs->make<TH1F>("hm4rec2OSr132405sig","M_{4#pi} OS #rho#rho sigma",massbins,0,5.);
  histosTH1F["hm4rec2OSr142305sig"] = fs->make<TH1F>("hm4rec2OSr142305sig","M_{4#pi} OS #rho#rho sigma",massbins,0,5.);
  //
  histosTH1F["hm4rec2OSr123405sig2"] = fs->make<TH1F>("hm4rec2OSr123405sig2","M_{4#pi} OS #rho#rho 60 MeV",massbins,0,5.);
  histosTH1F["hm4rec2OSr132405sig2"] = fs->make<TH1F>("hm4rec2OSr132405sig2","M_{4#pi} OS #rho#rho 60 MeV",massbins,0,5.);
  histosTH1F["hm4rec2OSr142305sig2"] = fs->make<TH1F>("hm4rec2OSr142305sig2","M_{4#pi} OS #rho#rho 60 MeV",massbins,0,5.);
  //
  histosTH1F["hm4rec2OSr123405sig3"] = fs->make<TH1F>("hm4rec2OSr123405sig3","M_{4#pi} OS #rho#rho 10 MeV",massbins,0,5.);
  histosTH1F["hm4rec2OSr132405sig3"] = fs->make<TH1F>("hm4rec2OSr132405sig3","M_{4#pi} OS #rho#rho 10 MeV",massbins,0,5.);
  histosTH1F["hm4rec2OSr142305sig3"] = fs->make<TH1F>("hm4rec2OSr142305sig3","M_{4#pi} OS #rho#rho 10 MeV",massbins,0,5.);
  //
  histosTH1F["hm4rec2OSr123405sig4"] = fs->make<TH1F>("hm4rec2OSr123405sig4","M_{4#pi} OS #rho#rho 20 MeV",massbins,0,5.);
  histosTH1F["hm4rec2OSr132405sig4"] = fs->make<TH1F>("hm4rec2OSr132405sig4","M_{4#pi} OS #rho#rho 20 MeV",massbins,0,5.);
  histosTH1F["hm4rec2OSr142305sig4"] = fs->make<TH1F>("hm4rec2OSr142305sig4","M_{4#pi} OS #rho#rho 20 MeV",massbins,0,5.);
  //
  histosTH1F["hm4rec2OSr123405sig5"] = fs->make<TH1F>("hm4rec2OSr123405sig5","M_{4#pi} OS #rho#rho 30 MeV",massbins,0,5.);
  histosTH1F["hm4rec2OSr132405sig5"] = fs->make<TH1F>("hm4rec2OSr132405sig5","M_{4#pi} OS #rho#rho 30 MeV",massbins,0,5.);
  histosTH1F["hm4rec2OSr142305sig5"] = fs->make<TH1F>("hm4rec2OSr142305sig5","M_{4#pi} OS #rho#rho 30 MeV",massbins,0,5.);
  //
  histosTH1F["hm4rec2OSr123405sig6"] = fs->make<TH1F>("hm4rec2OSr123405sig6","M_{4#pi} OS #rho#rho 40 MeV",massbins,0,5.);
  histosTH1F["hm4rec2OSr132405sig6"] = fs->make<TH1F>("hm4rec2OSr132405sig6","M_{4#pi} OS #rho#rho 40 MeV",massbins,0,5.);
  histosTH1F["hm4rec2OSr142305sig6"] = fs->make<TH1F>("hm4rec2OSr142305sig6","M_{4#pi} OS #rho#rho 40 MeV",massbins,0,5.);
  //
  histosTH1F["hm4rec2OSmrec123405"] = fs->make<TH1F>("hm4rec2OSmrec123405","M_{4#pi} OS mrec1234",massbins,0,5.);
  histosTH1F["hm4rec2OSmrec132405"] = fs->make<TH1F>("hm4rec2OSmrec132405","M_{4#pi} OS mrec1324",massbins,0,5.);
  histosTH1F["hm4rec2OSmrec142305"] = fs->make<TH1F>("hm4rec2OSmrec142305","M_{4#pi} OS mrec1423",massbins,0,5.);

  //
  // testing mix-up channels ...no 2V0
  histosTH1F["hm4rec2OSm123405nov"] = fs->make<TH1F>("hm4rec2OSm123405nov","M_{4#pi} OS no 2V0",massbins,0,5.);
  histosTH1F["hm4rec2OSm132405nov"] = fs->make<TH1F>("hm4rec2OSm132405nov","M_{4#pi} OS no 2V0",massbins,0,5.);
  histosTH1F["hm4rec2OSm142305nov"] = fs->make<TH1F>("hm4rec2OSm142305nov","M_{4#pi} OS no 2V0",massbins,0,5.);
  // testing mix-up channels ...yes 2V0
  histosTH1F["hm4rec2OSm123405yesv"] = fs->make<TH1F>("hm4rec2OSm123405yesv","M_{4#pi} OS yes 2V0",massbins,0,5.);
  histosTH1F["hm4rec2OSm132405yesv"] = fs->make<TH1F>("hm4rec2OSm132405yesv","M_{4#pi} OS yes 2V0",massbins,0,5.);
  histosTH1F["hm4rec2OSm142305yesv"] = fs->make<TH1F>("hm4rec2OSm142305yesv","M_{4#pi} OS yes 2V0",massbins,0,5.);
  // testing mix-up channels ...yes 1V0
  histosTH1F["hm4rec2OSm123405yes1"] = fs->make<TH1F>("hm4rec2OSm123405yes1","M_{4#pi} OS yes 2V0",massbins,0,5.);
  histosTH1F["hm4rec2OSm132405yes1"] = fs->make<TH1F>("hm4rec2OSm132405yes1","M_{4#pi} OS yes 2V0",massbins,0,5.);
  histosTH1F["hm4rec2OSm142305yes1"] = fs->make<TH1F>("hm4rec2OSm142305yes1","M_{4#pi} OS yes 2V0",massbins,0,5.);
  // fixed
  histosTH1F["hm4rec2OSm123405fix12"] = fs->make<TH1F>("hm4rec2OSm123405fix12","M_{4#pi} OS",massbins,0,5.);
  histosTH1F["hm4rec2OSm123405fix34"] = fs->make<TH1F>("hm4rec2OSm123405fix34","M_{4#pi} OS",massbins,0,5.);
  // fixed
  histosTH1F["hm4rec2OSm132405fix13"] = fs->make<TH1F>("hm4rec2OSm132405fix13","M_{4#pi} OS",massbins,0,5.);
  histosTH1F["hm4rec2OSm132405fix24"] = fs->make<TH1F>("hm4rec2OSm132405fix24","M_{4#pi} OS",massbins,0,5.);
  // fixed
  histosTH1F["hm4rec2OSm142305fix14"] = fs->make<TH1F>("hm4rec2OSm142305fix14","M_{4#pi} OS",massbins,0,5.);
  histosTH1F["hm4rec2OSm142305fix23"] = fs->make<TH1F>("hm4rec2OSm142305fix23","M_{4#pi} OS",massbins,0,5.);
  //
  histosTH1F["hm4rec2OSm1234t05"] = fs->make<TH1F>("hm4rec2OSm1234t05","M_{4#pi} OS",massbins,0,5.);
  histosTH1F["hm4rec2OSm1324t05"] = fs->make<TH1F>("hm4rec2OSm1324t05","M_{4#pi} OS",massbins,0,5.);
  histosTH1F["hm4rec2OSm1423t05"] = fs->make<TH1F>("hm4rec2OSm1423t05","M_{4#pi} OS",massbins,0,5.);
  //
  histosTH1F["hm4rec2OSr1234t05"] = fs->make<TH1F>("hm4rec2OSr1234t05","M_{4#pi} OS",massbins,0,5.);
  histosTH1F["hm4rec2OSr1324t05"] = fs->make<TH1F>("hm4rec2OSr1324t05","M_{4#pi} OS",massbins,0,5.);
  histosTH1F["hm4rec2OSr1423t05"] = fs->make<TH1F>("hm4rec2OSr1423t05","M_{4#pi} OS",massbins,0,5.);
  //
  //**histosTH1F["hm4rec2OSm1234t05pid"] = fs->make<TH1F>("hm4rec2OSm1234t05pid","M_{4#pi} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OSm1324t05pid"] = fs->make<TH1F>("hm4rec2OSm1324t05pid","M_{4#pi} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OSm1423t05pid"] = fs->make<TH1F>("hm4rec2OSm1423t05pid","M_{4#pi} OS",massbins,0,5.);
  //
  //...pT tracks
  histosTH1F["hm4rec2OSm123405pi1pt"] = fs->make<TH1F>("hm4rec2OSm123405pi1pt","M_{4#pi} OS pi1 pT",massbins,0,5.);
  histosTH1F["hm4rec2OSm123405pi2pt"] = fs->make<TH1F>("hm4rec2OSm123405pi2pt","M_{4#pi} OS pi2 pT",massbins,0,5.);
  histosTH1F["hm4rec2OSm123405pi3pt"] = fs->make<TH1F>("hm4rec2OSm123405pi3pt","M_{4#pi} OS pi3 pT",massbins,0,5.);
  histosTH1F["hm4rec2OSm123405pi4pt"] = fs->make<TH1F>("hm4rec2OSm123405pi4pt","M_{4#pi} OS pi4 pT",massbins,0,5.);
  //
  histosTH1F["hm4rec2OSm132405pi1pt"] = fs->make<TH1F>("hm4rec2OSm132405pi1pt","M_{4#pi} OS pi1 pT",massbins,0,5.);
  histosTH1F["hm4rec2OSm132405pi2pt"] = fs->make<TH1F>("hm4rec2OSm132405pi2pt","M_{4#pi} OS pi2 pT",massbins,0,5.);
  histosTH1F["hm4rec2OSm132405pi3pt"] = fs->make<TH1F>("hm4rec2OSm132405pi3pt","M_{4#pi} OS pi3 pT",massbins,0,5.);
  histosTH1F["hm4rec2OSm132405pi4pt"] = fs->make<TH1F>("hm4rec2OSm132405pi4pt","M_{4#pi} OS pi4 pT",massbins,0,5.);
  //
  histosTH1F["hm4rec2OSm142305pi1pt"] = fs->make<TH1F>("hm4rec2OSm142305pi1pt","M_{4#pi} OS pi1 pT",massbins,0,5.);
  histosTH1F["hm4rec2OSm142305pi2pt"] = fs->make<TH1F>("hm4rec2OSm142305pi2pt","M_{4#pi} OS pi2 pT",massbins,0,5.);
  histosTH1F["hm4rec2OSm142305pi3pt"] = fs->make<TH1F>("hm4rec2OSm142305pi3pt","M_{4#pi} OS pi3 pT",massbins,0,5.);
  histosTH1F["hm4rec2OSm142305pi4pt"] = fs->make<TH1F>("hm4rec2OSm142305pi4pt","M_{4#pi} OS pi4 pT",massbins,0,5.);
  //...end of pT tracks
  //
  //**histosTH1F["hm4rec2OSm123406"] = fs->make<TH1F>("hm4rec2OSm123406","M_{4#pi} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OSm132406"] = fs->make<TH1F>("hm4rec2OSm132406","M_{4#pi} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OSm142306"] = fs->make<TH1F>("hm4rec2OSm142306","M_{4#pi} OS",massbins,0,5.);
  //
  //...cut k1
  histosTH1F["hm4rec2OSm1234k"] = fs->make<TH1F>("hm4rec2OSm1234k","M_{4K} OS",massbins,0,5.);
  histosTH1F["hm4rec2OSm1324k"] = fs->make<TH1F>("hm4rec2OSm1324k","M_{4K} OS",massbins,0,5.);
  histosTH1F["hm4rec2OSm1423k"] = fs->make<TH1F>("hm4rec2OSm1423k","M_{4K} OS",massbins,0,5.);
  //...sig 10 MeV
  histosTH1F["hm4rec2OSk1234sig"] = fs->make<TH1F>("hm4rec2OSk1234sig","M_{4K} OS sig 10 MeV",massbins,0,5.);
  histosTH1F["hm4rec2OSk1324sig"] = fs->make<TH1F>("hm4rec2OSk1324sig","M_{4K} OS sig 10 MeV",massbins,0,5.);
  histosTH1F["hm4rec2OSk1423sig"] = fs->make<TH1F>("hm4rec2OSk1423sig","M_{4K} OS sig 10 MeV",massbins,0,5.);
  //
  histosTH1F["hm4rec2OSk1423exo"] = fs->make<TH1F>("hm4rec2OSk1423exo","M_{4K} OS sig 58 MeV",massbins,0,5.);
  //
  histosTH1F["hm4rec2OSk1234sigpid"] = fs->make<TH1F>("hm4rec2OSk1234sigpid","M_{4K} OS sig 10 MeV pid",massbins,0,5.);
  histosTH1F["hm4rec2OSk1324sigpid"] = fs->make<TH1F>("hm4rec2OSk1324sigpid","M_{4K} OS sig 10 MeV pid",massbins,0,5.);
  histosTH1F["hm4rec2OSk1423sigpid"] = fs->make<TH1F>("hm4rec2OSk1423sigpid","M_{4K} OS sig 10 MeV pid",massbins,0,5.);
  //
  histosTH1F["hm4rec2OSk1423exopid"] = fs->make<TH1F>("hm4rec2OSk1423exopid","M_{4K} OS sig 58 MeV pid",massbins,0,5.);
  //..sig 20 MeV
  histosTH1F["hm4rec2OSk1234sig2"] = fs->make<TH1F>("hm4rec2OSk1234sig2","M_{4K} OS sig 20 MeV",massbins,0,5.);
  histosTH1F["hm4rec2OSk1324sig2"] = fs->make<TH1F>("hm4rec2OSk1324sig2","M_{4K} OS sig 20 MeV",massbins,0,5.);
  histosTH1F["hm4rec2OSk1423sig2"] = fs->make<TH1F>("hm4rec2OSk1423sig2","M_{4K} OS sig 20 MeV",massbins,0,5.);
  histosTH1F["hm4rec2OSk1234sig2pid"] = fs->make<TH1F>("hm4rec2OSk1234sig2pid","M_{4K} OS sig 20 MeV pid",massbins,0,5.);
  histosTH1F["hm4rec2OSk1324sig2pid"] = fs->make<TH1F>("hm4rec2OSk1324sig2pid","M_{4K} OS sig 20 MeV pid",massbins,0,5.);
  histosTH1F["hm4rec2OSk1423sig2pid"] = fs->make<TH1F>("hm4rec2OSk1423sig2pid","M_{4K} OS sig 20 MeV pid",massbins,0,5.);
  //..sig 30 MeV
  histosTH1F["hm4rec2OSk1234sig3"] = fs->make<TH1F>("hm4rec2OSk1234sig3","M_{4K} OS sig 30 MeV",massbins,0,5.);
  histosTH1F["hm4rec2OSk1324sig3"] = fs->make<TH1F>("hm4rec2OSk1324sig3","M_{4K} OS sig 30 MeV",massbins,0,5.);
  histosTH1F["hm4rec2OSk1423sig3"] = fs->make<TH1F>("hm4rec2OSk1423sig3","M_{4K} OS sig 30 MeV",massbins,0,5.);
  histosTH1F["hm4rec2OSk1234sig3pid"] = fs->make<TH1F>("hm4rec2OSk1234sig3pid","M_{4K} OS sig 30 MeV pid",massbins,0,5.);
  histosTH1F["hm4rec2OSk1324sig3pid"] = fs->make<TH1F>("hm4rec2OSk1324sig3pid","M_{4K} OS sig 30 MeV pid",massbins,0,5.);
  histosTH1F["hm4rec2OSk1423sig3pid"] = fs->make<TH1F>("hm4rec2OSk1423sig3pid","M_{4K} OS sig 30 MeV pid",massbins,0,5.);
  //
  //...mass selection ...pions ...cut 0
  //***histosTH2F["h2dim2OSm12x34"]  = fs->make<TH2F>("h2dim2OSm12x34","M_{#pi_{1}#pi_{2}} vs M_{#pi_{3}#pi_{4}} K0s mass windows",massbins,0,5.,massbins,0,5.);
  //***histosTH2F["h2dim2OSm12x34t"]  = fs->make<TH2F>("h2dim2OSm12x34t","M_{#pi_{1}#pi_{2}} vs M_{#pi_{3}#pi_{4}} K0s",massbins,0,5.,massbins,0,5.);
  //***histosTH2F["h2dim2OSpteta12"] = fs->make<TH2F>("h2dim2OSpteta12","pt_{#pi_{1}#pi_{2}} vs #eta_{#pi_{1}#pi_{2}} K0s mass windows",100,0,5.,100,-5.,5.);
  //***histosTH2F["h2dim2OSpteta34"] = fs->make<TH2F>("h2dim2OSpteta34","pt_{#pi_{3}#pi_{4}} vs #eta_{#pi_{3}#pi_{4}} K0s mass windows",100,0,5.,100,-5.,5.);
  //**histosTH1F["h2OSpt12"]  = fs->make<TH1F>("h2OSpt12","pt_{#pi_{1}#pi_{2}} K0s mass windows",100,0,5.);
  //**histosTH1F["h2OSpt34"]  = fs->make<TH1F>("h2OSpt34","pt_{#pi_{3}#pi_{4}} K0s mass windows",100,0,5.);
  //**histosTH1F["h2OSeta12"] = fs->make<TH1F>("h2OSeta12","#eta_{#pi_{1}#pi_{2}} K0s mass windows",100,-5.,5.);
  //**histosTH1F["h2OSeta34"] = fs->make<TH1F>("h2OSeta34","#eta_{#pi_{3}#pi_{4}} K0s mass windows",100,-5.,5.);
  //
  //***histosTH2F["h2dim2OSm13x24"]  = fs->make<TH2F>("h2dim2OSm13x24","M_{#pi_{1}#pi_{3}} vs M_{#pi_{2}#pi_{4}} K0s mass windows",massbins,0,5.,massbins,0,5.);
  //***histosTH2F["h2dim2OSm13x24t"]  = fs->make<TH2F>("h2dim2OSm13x24t","M_{#pi_{1}#pi_{3}} vs M_{#pi_{2}#pi_{4}} K0s",massbins,0,5.,massbins,0,5.);
  //***histosTH2F["h2dim2OSpteta13"] = fs->make<TH2F>("h2dim2OSpteta13","pt_{#pi_{1}#pi_{3}} vs #eta_{#pi_{1}#pi_{3}} K0s mass windows",100,0,5.,100,-5.,5.);
  //***histosTH2F["h2dim2OSpteta24"] = fs->make<TH2F>("h2dim2OSpteta24","pt_{#pi_{2}#pi_{4}} vs #eta_{#pi_{2}#pi_{4}} K0s mass windows",100,0,5.,100,-5.,5.);
  //**histosTH1F["h2OSpt13"]  = fs->make<TH1F>("h2OSpt13","pt_{#pi_{1}#pi_{3}} K0s mass windows",100,0,5.);
  //**histosTH1F["h2OSpt24"]  = fs->make<TH1F>("h2OSpt24","pt_{#pi_{2}#pi_{4}} K0s mass windows",100,0,5.);
  //**histosTH1F["h2OSeta13"] = fs->make<TH1F>("h2OSeta13","#eta_{#pi_{1}#pi_{3}} K0s mass windows",100,-5.,5.);
  //**histosTH1F["h2OSeta24"] = fs->make<TH1F>("h2OSeta24","#eta_{#pi_{2}#pi_{4}} K0s mass windows",100,-5.,5.);
  //  
  //***histosTH2F["h2dim2OSm14x23"]  = fs->make<TH2F>("h2dim2OSm14x23","M_{#pi_{1}#pi_{4}} vs M_{#pi_{2}#pi_{3}} K0s mass windows",massbins,0,5.,massbins,0,5.);
  //***histosTH2F["h2dim2OSm14x23t"]  = fs->make<TH2F>("h2dim2OSm14x23t","M_{#pi_{1}#pi_{4}} vs M_{#pi_{2}#pi_{3}} K0s",massbins,0,5.,massbins,0,5.);
  //***histosTH2F["h2dim2OSpteta14"] = fs->make<TH2F>("h2dim2OSpteta14","pt_{#pi_{1}#pi_{4}} vs #eta_{#pi_{1}#pi_{4}} K0s mass windows",100,0,5.,100,-5.,5.);
  //***histosTH2F["h2dim2OSpteta23"] = fs->make<TH2F>("h2dim2OSpteta23","pt_{#pi_{2}#pi_{3}} vs #eta_{#pi_{2}#pi_{3}} K0s mass windows",100,0,5.,100,-5.,5.);
  //**histosTH1F["h2OSpt14"]  = fs->make<TH1F>("h2OSpt14","pt_{#pi_{1}#pi_{4}} K0s mass windows",100,0,5.);
  //**histosTH1F["h2OSpt23"]  = fs->make<TH1F>("h2OSpt23","pt_{#pi_{2}#pi_{3}} K0s mass windows",100,0,5.);
  //**histosTH1F["h2OSeta14"] = fs->make<TH1F>("h2OSeta14","#eta_{#pi_{1}#pi_{4}} K0s mass windows",100,-5.,5.);
  //**histosTH1F["h2OSeta23"] = fs->make<TH1F>("h2OSeta23","#eta_{#pi_{2}#pi_{3}} K0s mass windows",100,-5.,5.);
  //
  //...mass selection ...pions ...cut 00
  //...mass selection ...pions ...cut 000
  //
  //...mass selection ...pions ...cut 04
  //***histosTH2F["h2dim2OSm12x3404"]  = fs->make<TH2F>("h2dim2OSm12x3404","M_{#pi_{1}#pi_{2}} vs M_{#pi_{3}#pi_{4}} K0s mass windows",massbins,0,5.,massbins,0,5.);
  //***histosTH2F["h2dim2OSm12x34t04"]  = fs->make<TH2F>("h2dim2OSm12x34t04","M_{#pi_{1}#pi_{2}} vs M_{#pi_{3}#pi_{4}} K0s",massbins,0,5.,massbins,0,5.);
  //***histosTH2F["h2dim2OSpteta1204"] = fs->make<TH2F>("h2dim2OSpteta1204","pt_{#pi_{1}#pi_{2}} vs #eta_{#pi_{1}#pi_{2}} K0s mass windows",100,0,5.,100,-5.,5.);
  //***histosTH2F["h2dim2OSpteta3404"] = fs->make<TH2F>("h2dim2OSpteta3404","pt_{#pi_{3}#pi_{4}} vs #eta_{#pi_{3}#pi_{4}} K0s mass windows",100,0,5.,100,-5.,5.);
  //**histosTH1F["h2OSpt1204"]  = fs->make<TH1F>("h2OSpt1204","pt_{#pi_{1}#pi_{2}} K0s mass windows",100,0,5.);
  //**histosTH1F["h2OSpt3404"]  = fs->make<TH1F>("h2OSpt3404","pt_{#pi_{3}#pi_{4}} K0s mass windows",100,0,5.);
  //**histosTH1F["h2OSeta1204"] = fs->make<TH1F>("h2OSeta1204","#eta_{#pi_{1}#pi_{2}} K0s mass windows",100,-5.,5.);
  //**histosTH1F["h2OSeta3404"] = fs->make<TH1F>("h2OSeta3404","#eta_{#pi_{3}#pi_{4}} K0s mass windows",100,-5.,5.);
  //
  //***histosTH2F["h2dim2OSm13x2404"]  = fs->make<TH2F>("h2dim2OSm13x2404","M_{#pi_{1}#pi_{3}} vs M_{#pi_{2}#pi_{4}} K0s mass windows",massbins,0,5.,massbins,0,5.);
  //***histosTH2F["h2dim2OSm13x24t04"]  = fs->make<TH2F>("h2dim2OSm13x24t04","M_{#pi_{1}#pi_{3}} vs M_{#pi_{2}#pi_{4}} K0s",massbins,0,5.,massbins,0,5.);
  //***histosTH2F["h2dim2OSpteta1304"] = fs->make<TH2F>("h2dim2OSpteta1304","pt_{#pi_{1}#pi_{3}} vs #eta_{#pi_{1}#pi_{3}} K0s mass windows",100,0,5.,100,-5.,5.);
  //***histosTH2F["h2dim2OSpteta2404"] = fs->make<TH2F>("h2dim2OSpteta2404","pt_{#pi_{2}#pi_{4}} vs #eta_{#pi_{2}#pi_{4}} K0s mass windows",100,0,5.,100,-5.,5.);
  //**histosTH1F["h2OSpt1304"]  = fs->make<TH1F>("h2OSpt1304","pt_{#pi_{1}#pi_{3}} K0s mass windows",100,0,5.);
  //**histosTH1F["h2OSpt2404"]  = fs->make<TH1F>("h2OSpt2404","pt_{#pi_{2}#pi_{4}} K0s mass windows",100,0,5.);
  //**histosTH1F["h2OSeta1304"] = fs->make<TH1F>("h2OSeta1304","#eta_{#pi_{1}#pi_{3}} K0s mass windows",100,-5.,5.);
  //**histosTH1F["h2OSeta2404"] = fs->make<TH1F>("h2OSeta2404","#eta_{#pi_{2}#pi_{4}} K0s mass windows",100,-5.,5.);
  //  
  //***histosTH2F["h2dim2OSm14x2304"]  = fs->make<TH2F>("h2dim2OSm14x2304","M_{#pi_{1}#pi_{4}} vs M_{#pi_{2}#pi_{3}} K0s mass windows",massbins,0,5.,massbins,0,5.);
  //***histosTH2F["h2dim2OSm14x23t04"]  = fs->make<TH2F>("h2dim2OSm14x23t04","M_{#pi_{1}#pi_{4}} vs M_{#pi_{2}#pi_{3}} K0s",massbins,0,5.,massbins,0,5.);
  //***histosTH2F["h2dim2OSpteta1404"] = fs->make<TH2F>("h2dim2OSpteta1404","pt_{#pi_{1}#pi_{4}} vs #eta_{#pi_{1}#pi_{4}} K0s mass windows",100,0,5.,100,-5.,5.);
  //***histosTH2F["h2dim2OSpteta2304"] = fs->make<TH2F>("h2dim2OSpteta2304","pt_{#pi_{2}#pi_{3}} vs #eta_{#pi_{2}#pi_{3}} K0s mass windows",100,0,5.,100,-5.,5.);
  //**histosTH1F["h2OSpt1404"]  = fs->make<TH1F>("h2OSpt1404","pt_{#pi_{1}#pi_{4}} K0s mass windows",100,0,5.);
  //**histosTH1F["h2OSpt2304"]  = fs->make<TH1F>("h2OSpt2304","pt_{#pi_{2}#pi_{3}} K0s mass windows",100,0,5.);
  //**histosTH1F["h2OSeta1404"] = fs->make<TH1F>("h2OSeta1404","#eta_{#pi_{1}#pi_{4}} K0s mass windows",100,-5.,5.);
  //**histosTH1F["h2OSeta2304"] = fs->make<TH1F>("h2OSeta2304","#eta_{#pi_{2}#pi_{3}} K0s mass windows",100,-5.,5.);
  //
  //...mass selection ...pions ...cut 05
  histosTH2F["h2dim2OSm12x3405"]  = fs->make<TH2F>("h2dim2OSm12x3405","M_{#pi_{1}#pi_{2}} vs M_{#pi_{3}#pi_{4}} K0s mass windows",massbins,0,5.,massbins,0,5.);
  histosTH2F["h2dim2OSm12x34t05"]  = fs->make<TH2F>("h2dim2OSm12x34t05","M_{#pi_{1}#pi_{2}} vs M_{#pi_{3}#pi_{4}} K0s",massbins,0,5.,massbins,0,5.);
  histosTH2F["h2dim2OSpteta1205"] = fs->make<TH2F>("h2dim2OSpteta1205","pt_{#pi_{1}#pi_{2}} vs #eta_{#pi_{1}#pi_{2}} K0s mass windows",100,0,5.,100,-5.,5.);
  histosTH2F["h2dim2OSpteta3405"] = fs->make<TH2F>("h2dim2OSpteta3405","pt_{#pi_{3}#pi_{4}} vs #eta_{#pi_{3}#pi_{4}} K0s mass windows",100,0,5.,100,-5.,5.);
  histosTH1F["h2OSpt1205"]  = fs->make<TH1F>("h2OSpt1205","pt_{#pi_{1}#pi_{2}} K0s mass windows",100,0,5.);
  histosTH1F["h2OSpt3405"]  = fs->make<TH1F>("h2OSpt3405","pt_{#pi_{3}#pi_{4}} K0s mass windows",100,0,5.);
  histosTH1F["h2OSeta1205"] = fs->make<TH1F>("h2OSeta1205","#eta_{#pi_{1}#pi_{2}} K0s mass windows",100,-5.,5.);
  histosTH1F["h2OSeta3405"] = fs->make<TH1F>("h2OSeta3405","#eta_{#pi_{3}#pi_{4}} K0s mass windows",100,-5.,5.);
  //
  histosTH2F["h2dim2OSm13x2405"]  = fs->make<TH2F>("h2dim2OSm13x2405","M_{#pi_{1}#pi_{3}} vs M_{#pi_{2}#pi_{4}} K0s mass windows",massbins,0,5.,massbins,0,5.);
  histosTH2F["h2dim2OSm13x24t05"]  = fs->make<TH2F>("h2dim2OSm13x24t05","M_{#pi_{1}#pi_{3}} vs M_{#pi_{2}#pi_{4}} K0s",massbins,0,5.,massbins,0,5.);
  histosTH2F["h2dim2OSpteta1305"] = fs->make<TH2F>("h2dim2OSpteta1305","pt_{#pi_{1}#pi_{3}} vs #eta_{#pi_{1}#pi_{3}} K0s mass windows",100,0,5.,100,-5.,5.);
  histosTH2F["h2dim2OSpteta2405"] = fs->make<TH2F>("h2dim2OSpteta2405","pt_{#pi_{2}#pi_{4}} vs #eta_{#pi_{2}#pi_{4}} K0s mass windows",100,0,5.,100,-5.,5.);
  histosTH1F["h2OSpt1305"]  = fs->make<TH1F>("h2OSpt1305","pt_{#pi_{1}#pi_{3}} K0s mass windows",100,0,5.);
  histosTH1F["h2OSpt2405"]  = fs->make<TH1F>("h2OSpt2405","pt_{#pi_{2}#pi_{4}} K0s mass windows",100,0,5.);
  histosTH1F["h2OSeta1305"] = fs->make<TH1F>("h2OSeta1305","#eta_{#pi_{1}#pi_{3}} K0s mass windows",100,-5.,5.);
  histosTH1F["h2OSeta2405"] = fs->make<TH1F>("h2OSeta2405","#eta_{#pi_{2}#pi_{4}} K0s mass windows",100,-5.,5.);
  //  
  histosTH2F["h2dim2OSm14x2305"]  = fs->make<TH2F>("h2dim2OSm14x2305","M_{#pi_{1}#pi_{4}} vs M_{#pi_{2}#pi_{3}} K0s mass windows",massbins,0,5.,massbins,0,5.);
  histosTH2F["h2dim2OSm14x23t05"]  = fs->make<TH2F>("h2dim2OSm14x23t05","M_{#pi_{1}#pi_{4}} vs M_{#pi_{2}#pi_{3}} K0s",massbins,0,5.,massbins,0,5.);
  histosTH2F["h2dim2OSpteta1405"] = fs->make<TH2F>("h2dim2OSpteta1405","pt_{#pi_{1}#pi_{4}} vs #eta_{#pi_{1}#pi_{4}} K0s mass windows",100,0,5.,100,-5.,5.);
  histosTH2F["h2dim2OSpteta2305"] = fs->make<TH2F>("h2dim2OSpteta2305","pt_{#pi_{2}#pi_{3}} vs #eta_{#pi_{2}#pi_{3}} K0s mass windows",100,0,5.,100,-5.,5.);
  histosTH1F["h2OSpt1405"]  = fs->make<TH1F>("h2OSpt1405","pt_{#pi_{1}#pi_{4}} K0s mass windows",100,0,5.);
  histosTH1F["h2OSpt2305"]  = fs->make<TH1F>("h2OSpt2305","pt_{#pi_{2}#pi_{3}} K0s mass windows",100,0,5.);
  histosTH1F["h2OSeta1405"] = fs->make<TH1F>("h2OSeta1405","#eta_{#pi_{1}#pi_{4}} K0s mass windows",100,-5.,5.);
  histosTH1F["h2OSeta2305"] = fs->make<TH1F>("h2OSeta2305","#eta_{#pi_{2}#pi_{3}} K0s mass windows",100,-5.,5.);
  //
  histosTH2F["h2dim2OSr12x34t05"]  = fs->make<TH2F>("h2dim2OSr12x34t05","M_{#pi_{1}#pi_{2}} vs M_{#pi_{3}#pi_{4}} K0s",massbins,0,5.,massbins,0,5.);
  histosTH2F["h2dim2OSr13x24t05"]  = fs->make<TH2F>("h2dim2OSr13x24t05","M_{#pi_{1}#pi_{3}} vs M_{#pi_{2}#pi_{4}} K0s",massbins,0,5.,massbins,0,5.);
  histosTH2F["h2dim2OSr14x23t05"]  = fs->make<TH2F>("h2dim2OSr14x23t05","M_{#pi_{1}#pi_{4}} vs M_{#pi_{2}#pi_{3}} K0s",massbins,0,5.,massbins,0,5.);
  //
  //...mass selection ...pions ...cut 06
  //***histosTH2F["h2dim2OSm12x3406"]  = fs->make<TH2F>("h2dim2OSm12x3406","M_{#pi_{1}#pi_{2}} vs M_{#pi_{3}#pi_{4}} K0s mass windows",massbins,0,5.,massbins,0,5.);
  //***histosTH2F["h2dim2OSpteta1206"] = fs->make<TH2F>("h2dim2OSpteta1206","pt_{#pi_{1}#pi_{2}} vs #eta_{#pi_{1}#pi_{2}} K0s mass windows",100,0,5.,100,-5.,5.);
  //***histosTH2F["h2dim2OSpteta3406"] = fs->make<TH2F>("h2dim2OSpteta3406","pt_{#pi_{3}#pi_{4}} vs #eta_{#pi_{3}#pi_{4}} K0s mass windows",100,0,5.,100,-5.,5.);
  //**histosTH1F["h2OSpt1206"]  = fs->make<TH1F>("h2OSpt1206","pt_{#pi_{1}#pi_{2}} K0s mass windows",100,0,5.);
  //**histosTH1F["h2OSpt3406"]  = fs->make<TH1F>("h2OSpt3406","pt_{#pi_{3}#pi_{4}} K0s mass windows",100,0,5.);
  //**histosTH1F["h2OSeta1206"] = fs->make<TH1F>("h2OSeta1206","#eta_{#pi_{1}#pi_{2}} K0s mass windows",100,-5.,5.);
  //**histosTH1F["h2OSeta3406"] = fs->make<TH1F>("h2OSeta3406","#eta_{#pi_{3}#pi_{4}} K0s mass windows",100,-5.,5.);
  //
  //***histosTH2F["h2dim2OSm13x2406"]  = fs->make<TH2F>("h2dim2OSm13x2406","M_{#pi_{1}#pi_{3}} vs M_{#pi_{2}#pi_{4}} K0s mass windows",massbins,0,5.,massbins,0,5.);
  //***histosTH2F["h2dim2OSpteta1306"] = fs->make<TH2F>("h2dim2OSpteta1306","pt_{#pi_{1}#pi_{3}} vs #eta_{#pi_{1}#pi_{3}} K0s mass windows",100,0,5.,100,-5.,5.);
  //***histosTH2F["h2dim2OSpteta2406"] = fs->make<TH2F>("h2dim2OSpteta2406","pt_{#pi_{2}#pi_{4}} vs #eta_{#pi_{2}#pi_{4}} K0s mass windows",100,0,5.,100,-5.,5.);
  //**histosTH1F["h2OSpt1306"]  = fs->make<TH1F>("h2OSpt1306","pt_{#pi_{1}#pi_{3}} K0s mass windows",100,0,5.);
  //**histosTH1F["h2OSpt2406"]  = fs->make<TH1F>("h2OSpt2406","pt_{#pi_{2}#pi_{4}} K0s mass windows",100,0,5.);
  //**histosTH1F["h2OSeta1306"] = fs->make<TH1F>("h2OSeta1306","#eta_{#pi_{1}#pi_{3}} K0s mass windows",100,-5.,5.);
  //**histosTH1F["h2OSeta2406"] = fs->make<TH1F>("h2OSeta2406","#eta_{#pi_{2}#pi_{4}} K0s mass windows",100,-5.,5.);
  //  
  //***histosTH2F["h2dim2OSm14x2306"]  = fs->make<TH2F>("h2dim2OSm14x2306","M_{#pi_{1}#pi_{4}} vs M_{#pi_{2}#pi_{3}} K0s mass windows",massbins,0,5.,massbins,0,5.);
  //***histosTH2F["h2dim2OSpteta1406"] = fs->make<TH2F>("h2dim2OSpteta1406","pt_{#pi_{1}#pi_{4}} vs #eta_{#pi_{1}#pi_{4}} K0s mass windows",100,0,5.,100,-5.,5.);
  //***histosTH2F["h2dim2OSpteta2306"] = fs->make<TH2F>("h2dim2OSpteta2306","pt_{#pi_{2}#pi_{3}} vs #eta_{#pi_{2}#pi_{3}} K0s mass windows",100,0,5.,100,-5.,5.);
  //**histosTH1F["h2OSpt1406"]  = fs->make<TH1F>("h2OSpt1406","pt_{#pi_{1}#pi_{4}} K0s mass windows",100,0,5.);
  //**histosTH1F["h2OSpt2306"]  = fs->make<TH1F>("h2OSpt2306","pt_{#pi_{2}#pi_{3}} K0s mass windows",100,0,5.);
  //**histosTH1F["h2OSeta1406"] = fs->make<TH1F>("h2OSeta1406","#eta_{#pi_{1}#pi_{4}} K0s mass windows",100,-5.,5.);
  //**histosTH1F["h2OSeta2306"] = fs->make<TH1F>("h2OSeta2306","#eta_{#pi_{2}#pi_{3}} K0s mass windows",100,-5.,5.);
  //
  //...cut k1 
  histosTH2F["h2dim2OSm12x34k"]  = fs->make<TH2F>("h2dim2OSm12x34k","M_{K_{1}K_{2}} vs M_{K_{3}K_{4}} K^{+}K^{-}",massbins,0,5.,massbins,0,5.);
  //***histosTH2F["h2dim2OSpteta12k"] = fs->make<TH2F>("h2dim2OSpteta12k","pt_{K_{1}K_{2}} vs #eta_{K_{1}K_{2}} K^{+}K^{-}",100,0,5.,100,-5.,5.);
  //***histosTH2F["h2dim2OSpteta34k"] = fs->make<TH2F>("h2dim2OSpteta34k","pt_{K_{3}K_{4}} vs #eta_{K_{3}K_{4}} K^{+}K^{-}",100,0,5.,100,-5.,5.);
  //**histosTH1F["h2OSpt12k"]  = fs->make<TH1F>("h2OSpt12k","pt_{K_{1}K_{2}} K^{+}K^{-}",100,0,5.);
  //**histosTH1F["h2OSpt34k"]  = fs->make<TH1F>("h2OSpt34k","pt_{K_{3}K_{4}} K^{+}K^{-}",100,0,5.);
  //**histosTH1F["h2OSeta12k"] = fs->make<TH1F>("h2OSeta12k","#eta_{K_{1}K_{2}} K^{+}K^{-}",100,-5.,5.);
  //**histosTH1F["h2OSeta34k"] = fs->make<TH1F>("h2OSeta34k","#eta_{K_{3}K_{4}} K^{+}K^{-}",100,-5.,5.);
  //
  histosTH2F["h2dim2OSm13x24k"]  = fs->make<TH2F>("h2dim2OSm13x24k","M_{K_{1}K_{3}} vs M_{K_{2}K_{4}} K^{+}K^{-}",massbins,0,5.,massbins,0,5.);
  //***histosTH2F["h2dim2OSpteta13k"] = fs->make<TH2F>("h2dim2OSpteta13k","pt_{K_{1}K_{3}} vs #eta_{K_{1}K_{3}} K^{+}K^{-}",100,0,5.,100,-5.,5.);
  //***histosTH2F["h2dim2OSpteta24k"] = fs->make<TH2F>("h2dim2OSpteta24k","pt_{K_{2}K_{4}} vs #eta_{K_{2}K_{4}} K^{+}K^{-}",100,0,5.,100,-5.,5.);
  //**histosTH1F["h2OSpt13k"]  = fs->make<TH1F>("h2OSpt13k","pt_{K_{1}K_{3}} K^{+}K^{-}",100,0,5.);
  //**histosTH1F["h2OSpt24k"]  = fs->make<TH1F>("h2OSpt24k","pt_{K_{2}K_{4}} K^{+}K^{-}",100,0,5.);
  //**histosTH1F["h2OSeta13k"] = fs->make<TH1F>("h2OSeta13k","#eta_{K_{1}K_{3}} K^{+}K^{-}",100,-5.,5.);
  //**histosTH1F["h2OSeta24k"] = fs->make<TH1F>("h2OSeta24k","#eta_{K_{2}K_{4}} K^{+}K^{-}",100,-5.,5.);
  //  
  histosTH2F["h2dim2OSm14x23k"]  = fs->make<TH2F>("h2dim2OSm14x23k","M_{K_{1}K_{4}} vs M_{K_{2}K_{3}} K^{+}K^{-}",massbins,0,5.,massbins,0,5.);
  //***histosTH2F["h2dim2OSpteta14k"] = fs->make<TH2F>("h2dim2OSpteta14k","pt_{K_{1}K_{4}} vs #eta_{K_{1}K_{4}} K^{+}K^{-}",100,0,5.,100,-5.,5.);
  //***histosTH2F["h2dim2OSpteta23k"] = fs->make<TH2F>("h2dim2OSpteta23k","pt_{K_{2}K_{3}} vs #eta_{K_{2}K_{3}} K^{+}K^{-}",100,0,5.,100,-5.,5.);
  //**histosTH1F["h2OSpt14k"]  = fs->make<TH1F>("h2OSpt14k","pt_{K_{1}K_{4}} K^{+}K^{-}",100,0,5.);
  //**histosTH1F["h2OSpt23k"]  = fs->make<TH1F>("h2OSpt23k","pt_{K_{2}K_{3}} K^{+}K^{-}",100,0,5.);
  //**histosTH1F["h2OSeta14k"] = fs->make<TH1F>("h2OSeta14k","#eta_{K_{1}K_{4}} K^{+}K^{-}",100,-5.,5.);
  //**histosTH1F["h2OSeta23k"] = fs->make<TH1F>("h2OSeta23k","#eta_{K_{2}K_{3}} K^{+}K^{-}",100,-5.,5.);
    
  //**histosTH1F["hm4rec2OSvee"] = fs->make<TH1F>("hm4rec2OSvee","M_{4#pi} OS",massbins,0,5.);
  //
  //**histosTH1F["hm4rec2OSvee11"] = fs->make<TH1F>("hm4rec2OSvee11","M_{K#pi} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OSvee02"] = fs->make<TH1F>("hm4rec2OSvee02","M_{4#pi} OS",2*massbins,0,10.);
  //**histosTH1F["hm4rec2OSvee01"] = fs->make<TH1F>("hm4rec2OSvee01","M_{2#pi} OS",massbins,0,5.);
  //
  histosTH1F["hm4rec2OSveeno11Kpi"] = fs->make<TH1F>("hm4rec2OSveeno11Kpi","M_{Kpi} OS",massbins,0,5.);
  histosTH1F["hm4rec2OSveeno11"] = fs->make<TH1F>("hm4rec2OSveeno11","M_{K0sK0s} OS",massbins,0,5.);
  histosTH1F["hm4rec2OSveeno11x"] = fs->make<TH1F>("hm4rec2OSveeno11x","M_{K0sK0s} OS",massbins,0,5.);
  histosTH1F["hm4rec2OSveeno02"] = fs->make<TH1F>("hm4rec2OSveeno02","M_{K0sK0s} OS",massbins,0,5.);
  histosTH1F["hm4rec2OSveeno02cor"] = fs->make<TH1F>("hm4rec2OSveeno02cor","M_{K0sK0s} OS efficiency correction",massbins,0,5.);
  histosTH1F["hm4rec2OSveeno02col"] = fs->make<TH1F>("hm4rec2OSveeno02col","M_{K0sK0s} OS kscoll",massbins,0,5.);
  histosTH1F["hm4rec2OSveeno02colcor"] = fs->make<TH1F>("hm4rec2OSveeno02colcor","M_{K0sK0s} OS efficiency correction kscoll",massbins,0,5.);
  histosTH1F["hm4rec2OSveeno01"] = fs->make<TH1F>("hm4rec2OSveeno01","M_{K0sK0s} OS",massbins,0,5.);
  //
  histosTH1F["hm4rec2OSveenopid02"] = fs->make<TH1F>("hm4rec2OSveenopid02","M_{K0sK0s} OS",massbins,0,5.);
  histosTH1F["hm4rec2OSveenopid02cor"] = fs->make<TH1F>("hm4rec2OSveenopid02cor","M_{K0sK0s} OS efficiency correction",massbins,0,5.);
  histosTH1F["hm4rec2OSveenopid02col"] = fs->make<TH1F>("hm4rec2OSveenopid02col","M_{K0sK0s} OS kscoll",massbins,0,5.);
  histosTH1F["hm4rec2OSveenopid02colcor"] = fs->make<TH1F>("hm4rec2OSveenopid02colcor","M_{K0sK0s} OS efficiency correction kscoll",massbins,0,5.);
  histosTH1F["hm4rec2OSveenopid01"] = fs->make<TH1F>("hm4rec2OSveenopid01","M_{K0sK0s} OS",massbins,0,5.);
  //
  // ...OR ...02
  histosTH1F["hm4rec2OSveeno02t12a"] = fs->make<TH1F>("hm4rec2OSveeno02t12a","M_{K0sK0s} t1+t2 < 0.4",massbins,0,5.);
  histosTH1F["hm4rec2OSveeno02t12b"] = fs->make<TH1F>("hm4rec2OSveeno02t12b","M_{K0sK0s} t1+t2 < 0.3",massbins,0,5.);
  histosTH1F["hm4rec2OSveeno02t12c"] = fs->make<TH1F>("hm4rec2OSveeno02t12c","M_{K0sK0s} t1+t2 < 0.2",massbins,0,5.);
  histosTH1F["hm4rec2OSveeno02t12d"] = fs->make<TH1F>("hm4rec2OSveeno02t12d","M_{K0sK0s} t1+t2 < 0.17",massbins,0,5.);
  histosTH1F["hm4rec2OSveeno02t12e"] = fs->make<TH1F>("hm4rec2OSveeno02t12e","M_{K0sK0s} t1+t2 < 0.1",massbins,0,5.);
  histosTH1F["hm4rec2OSveeno02t12f"] = fs->make<TH1F>("hm4rec2OSveeno02t12f","M_{K0sK0s} t1+t2 < 0.08",massbins,0,5.);
  histosTH1F["hm4rec2OSveeno02t12g"] = fs->make<TH1F>("hm4rec2OSveeno02t12g","M_{K0sK0s} t1+t2 < 0.06",massbins,0,5.);
  histosTH1F["hm4rec2OSveeno02t12h"] = fs->make<TH1F>("hm4rec2OSveeno02t12h","M_{K0sK0s} t1+t2 < 0.05",massbins,0,5.);
  // ...OR complementary
  histosTH1F["hm4rec2OSveeno02t12aC"] = fs->make<TH1F>("hm4rec2OSveeno02t12aC","M_{K0sK0s} t1+t2 > 0.4",massbins,0,5.);
  histosTH1F["hm4rec2OSveeno02t12bC"] = fs->make<TH1F>("hm4rec2OSveeno02t12bC","M_{K0sK0s} t1+t2 > 0.3",massbins,0,5.);
  histosTH1F["hm4rec2OSveeno02t12cC"] = fs->make<TH1F>("hm4rec2OSveeno02t12cC","M_{K0sK0s} t1+t2 > 0.2",massbins,0,5.);
  histosTH1F["hm4rec2OSveeno02t12dC"] = fs->make<TH1F>("hm4rec2OSveeno02t12dC","M_{K0sK0s} t1+t2 > 0.17",massbins,0,5.);
  histosTH1F["hm4rec2OSveeno02t12eC"] = fs->make<TH1F>("hm4rec2OSveeno02t12eC","M_{K0sK0s} t1+t2 > 0.1",massbins,0,5.);
  histosTH1F["hm4rec2OSveeno02t12fC"] = fs->make<TH1F>("hm4rec2OSveeno02t12fC","M_{K0sK0s} t1+t2 > 0.08",massbins,0,5.);
  histosTH1F["hm4rec2OSveeno02t12gC"] = fs->make<TH1F>("hm4rec2OSveeno02t12gC","M_{K0sK0s} t1+t2 > 0.06",massbins,0,5.);
  histosTH1F["hm4rec2OSveeno02t12hC"] = fs->make<TH1F>("hm4rec2OSveeno02t12hC","M_{K0sK0s} t1+t2 > 0.05",massbins,0,5.);
  // ...AND
  histosTH1F["hm4rec2OSveeno02t12aa"] = fs->make<TH1F>("hm4rec2OSveeno02t12aa","M_{K0sK0s} t1+t2 < 0.4",massbins,0,5.);
  histosTH1F["hm4rec2OSveeno02t12bb"] = fs->make<TH1F>("hm4rec2OSveeno02t12bb","M_{K0sK0s} t1+t2 < 0.3",massbins,0,5.);
  histosTH1F["hm4rec2OSveeno02t12cc"] = fs->make<TH1F>("hm4rec2OSveeno02t12cc","M_{K0sK0s} t1+t2 < 0.2",massbins,0,5.);
  histosTH1F["hm4rec2OSveeno02t12dd"] = fs->make<TH1F>("hm4rec2OSveeno02t12dd","M_{K0sK0s} t1+t2 < 0.17",massbins,0,5.);
  histosTH1F["hm4rec2OSveeno02t12ee"] = fs->make<TH1F>("hm4rec2OSveeno02t12ee","M_{K0sK0s} t1+t2 < 0.1",massbins,0,5.);
  histosTH1F["hm4rec2OSveeno02t12ff"] = fs->make<TH1F>("hm4rec2OSveeno02t12ff","M_{K0sK0s} t1+t2 < 0.08",massbins,0,5.);
  histosTH1F["hm4rec2OSveeno02t12gg"] = fs->make<TH1F>("hm4rec2OSveeno02t12gg","M_{K0sK0s} t1+t2 < 0.06",massbins,0,5.);
  histosTH1F["hm4rec2OSveeno02t12hh"] = fs->make<TH1F>("hm4rec2OSveeno02t12hh","M_{K0sK0s} t1+t2 < 0.05",massbins,0,5.);
  // ...AND complementary
  histosTH1F["hm4rec2OSveeno02t12aaC"] = fs->make<TH1F>("hm4rec2OSveeno02t12aaC","M_{K0sK0s} t1+t2 > 0.4",massbins,0,5.);
  histosTH1F["hm4rec2OSveeno02t12bbC"] = fs->make<TH1F>("hm4rec2OSveeno02t12bbC","M_{K0sK0s} t1+t2 > 0.3",massbins,0,5.);
  histosTH1F["hm4rec2OSveeno02t12ccC"] = fs->make<TH1F>("hm4rec2OSveeno02t12ccC","M_{K0sK0s} t1+t2 > 0.2",massbins,0,5.);
  histosTH1F["hm4rec2OSveeno02t12ddC"] = fs->make<TH1F>("hm4rec2OSveeno02t12ddC","M_{K0sK0s} t1+t2 > 0.17",massbins,0,5.);
  histosTH1F["hm4rec2OSveeno02t12eeC"] = fs->make<TH1F>("hm4rec2OSveeno02t12eeC","M_{K0sK0s} t1+t2 > 0.1",massbins,0,5.);
  histosTH1F["hm4rec2OSveeno02t12ffC"] = fs->make<TH1F>("hm4rec2OSveeno02t12ffC","M_{K0sK0s} t1+t2 > 0.08",massbins,0,5.);
  histosTH1F["hm4rec2OSveeno02t12ggC"] = fs->make<TH1F>("hm4rec2OSveeno02t12ggC","M_{K0sK0s} t1+t2 > 0.06",massbins,0,5.);
  histosTH1F["hm4rec2OSveeno02t12hhC"] = fs->make<TH1F>("hm4rec2OSveeno02t12hhC","M_{K0sK0s} t1+t2 > 0.05",massbins,0,5.);
  //
  // ...OR ...02 ...NO PID
  histosTH1F["hm4rec2OSveenopid02t12a"] = fs->make<TH1F>("hm4rec2OSveenopid02t12a","M_{K0sK0s} t1+t2 < 0.4",massbins,0,5.);
  histosTH1F["hm4rec2OSveenopid02t12b"] = fs->make<TH1F>("hm4rec2OSveenopid02t12b","M_{K0sK0s} t1+t2 < 0.3",massbins,0,5.);
  histosTH1F["hm4rec2OSveenopid02t12c"] = fs->make<TH1F>("hm4rec2OSveenopid02t12c","M_{K0sK0s} t1+t2 < 0.2",massbins,0,5.);
  histosTH1F["hm4rec2OSveenopid02t12d"] = fs->make<TH1F>("hm4rec2OSveenopid02t12d","M_{K0sK0s} t1+t2 < 0.17",massbins,0,5.);
  histosTH1F["hm4rec2OSveenopid02t12e"] = fs->make<TH1F>("hm4rec2OSveenopid02t12e","M_{K0sK0s} t1+t2 < 0.1",massbins,0,5.);
  histosTH1F["hm4rec2OSveenopid02t12f"] = fs->make<TH1F>("hm4rec2OSveenopid02t12f","M_{K0sK0s} t1+t2 < 0.08",massbins,0,5.);
  histosTH1F["hm4rec2OSveenopid02t12g"] = fs->make<TH1F>("hm4rec2OSveenopid02t12g","M_{K0sK0s} t1+t2 < 0.06",massbins,0,5.);
  histosTH1F["hm4rec2OSveenopid02t12h"] = fs->make<TH1F>("hm4rec2OSveenopid02t12h","M_{K0sK0s} t1+t2 < 0.05",massbins,0,5.);
  // ...OR complementary
  histosTH1F["hm4rec2OSveenopid02t12aC"] = fs->make<TH1F>("hm4rec2OSveenopid02t12aC","M_{K0sK0s} t1+t2 > 0.4",massbins,0,5.);
  histosTH1F["hm4rec2OSveenopid02t12bC"] = fs->make<TH1F>("hm4rec2OSveenopid02t12bC","M_{K0sK0s} t1+t2 > 0.3",massbins,0,5.);
  histosTH1F["hm4rec2OSveenopid02t12cC"] = fs->make<TH1F>("hm4rec2OSveenopid02t12cC","M_{K0sK0s} t1+t2 > 0.2",massbins,0,5.);
  histosTH1F["hm4rec2OSveenopid02t12dC"] = fs->make<TH1F>("hm4rec2OSveenopid02t12dC","M_{K0sK0s} t1+t2 > 0.17",massbins,0,5.);
  histosTH1F["hm4rec2OSveenopid02t12eC"] = fs->make<TH1F>("hm4rec2OSveenopid02t12eC","M_{K0sK0s} t1+t2 > 0.1",massbins,0,5.);
  histosTH1F["hm4rec2OSveenopid02t12fC"] = fs->make<TH1F>("hm4rec2OSveenopid02t12fC","M_{K0sK0s} t1+t2 > 0.08",massbins,0,5.);
  histosTH1F["hm4rec2OSveenopid02t12gC"] = fs->make<TH1F>("hm4rec2OSveenopid02t12gC","M_{K0sK0s} t1+t2 > 0.06",massbins,0,5.);
  histosTH1F["hm4rec2OSveenopid02t12hC"] = fs->make<TH1F>("hm4rec2OSveenopid02t12hC","M_{K0sK0s} t1+t2 > 0.05",massbins,0,5.);
  // ...AND
  histosTH1F["hm4rec2OSveenopid02t12aa"] = fs->make<TH1F>("hm4rec2OSveenopid02t12aa","M_{K0sK0s} t1+t2 < 0.4",massbins,0,5.);
  histosTH1F["hm4rec2OSveenopid02t12bb"] = fs->make<TH1F>("hm4rec2OSveenopid02t12bb","M_{K0sK0s} t1+t2 < 0.3",massbins,0,5.);
  histosTH1F["hm4rec2OSveenopid02t12cc"] = fs->make<TH1F>("hm4rec2OSveenopid02t12cc","M_{K0sK0s} t1+t2 < 0.2",massbins,0,5.);
  histosTH1F["hm4rec2OSveenopid02t12dd"] = fs->make<TH1F>("hm4rec2OSveenopid02t12dd","M_{K0sK0s} t1+t2 < 0.17",massbins,0,5.);
  histosTH1F["hm4rec2OSveenopid02t12ee"] = fs->make<TH1F>("hm4rec2OSveenopid02t12ee","M_{K0sK0s} t1+t2 < 0.1",massbins,0,5.);
  histosTH1F["hm4rec2OSveenopid02t12ff"] = fs->make<TH1F>("hm4rec2OSveenopid02t12ff","M_{K0sK0s} t1+t2 < 0.08",massbins,0,5.);
  histosTH1F["hm4rec2OSveenopid02t12gg"] = fs->make<TH1F>("hm4rec2OSveenopid02t12gg","M_{K0sK0s} t1+t2 < 0.06",massbins,0,5.);
  histosTH1F["hm4rec2OSveenopid02t12hh"] = fs->make<TH1F>("hm4rec2OSveenopid02t12hh","M_{K0sK0s} t1+t2 < 0.05",massbins,0,5.);
  // ...AND complementary
  histosTH1F["hm4rec2OSveenopid02t12aaC"] = fs->make<TH1F>("hm4rec2OSveenopid02t12aaC","M_{K0sK0s} t1+t2 > 0.4",massbins,0,5.);
  histosTH1F["hm4rec2OSveenopid02t12bbC"] = fs->make<TH1F>("hm4rec2OSveenopid02t12bbC","M_{K0sK0s} t1+t2 > 0.3",massbins,0,5.);
  histosTH1F["hm4rec2OSveenopid02t12ccC"] = fs->make<TH1F>("hm4rec2OSveenopid02t12ccC","M_{K0sK0s} t1+t2 > 0.2",massbins,0,5.);
  histosTH1F["hm4rec2OSveenopid02t12ddC"] = fs->make<TH1F>("hm4rec2OSveenopid02t12ddC","M_{K0sK0s} t1+t2 > 0.17",massbins,0,5.);
  histosTH1F["hm4rec2OSveenopid02t12eeC"] = fs->make<TH1F>("hm4rec2OSveenopid02t12eeC","M_{K0sK0s} t1+t2 > 0.1",massbins,0,5.);
  histosTH1F["hm4rec2OSveenopid02t12ffC"] = fs->make<TH1F>("hm4rec2OSveenopid02t12ffC","M_{K0sK0s} t1+t2 > 0.08",massbins,0,5.);
  histosTH1F["hm4rec2OSveenopid02t12ggC"] = fs->make<TH1F>("hm4rec2OSveenopid02t12ggC","M_{K0sK0s} t1+t2 > 0.06",massbins,0,5.);
  histosTH1F["hm4rec2OSveenopid02t12hhC"] = fs->make<TH1F>("hm4rec2OSveenopid02t12hhC","M_{K0sK0s} t1+t2 > 0.05",massbins,0,5.);
  //
  // ...OR ...01
  histosTH1F["hm4rec2OSveeno01t12a"] = fs->make<TH1F>("hm4rec2OSveeno01t12a","M_{K0sK0s} t1+t2 < 0.4",massbins,0,5.);
  histosTH1F["hm4rec2OSveeno01t12b"] = fs->make<TH1F>("hm4rec2OSveeno01t12b","M_{K0sK0s} t1+t2 < 0.3",massbins,0,5.);
  histosTH1F["hm4rec2OSveeno01t12c"] = fs->make<TH1F>("hm4rec2OSveeno01t12c","M_{K0sK0s} t1+t2 < 0.2",massbins,0,5.);
  histosTH1F["hm4rec2OSveeno01t12d"] = fs->make<TH1F>("hm4rec2OSveeno01t12d","M_{K0sK0s} t1+t2 < 0.17",massbins,0,5.);
  histosTH1F["hm4rec2OSveeno01t12e"] = fs->make<TH1F>("hm4rec2OSveeno01t12e","M_{K0sK0s} t1+t2 < 0.1",massbins,0,5.);
  histosTH1F["hm4rec2OSveeno01t12f"] = fs->make<TH1F>("hm4rec2OSveeno01t12f","M_{K0sK0s} t1+t2 < 0.08",massbins,0,5.);
  histosTH1F["hm4rec2OSveeno01t12g"] = fs->make<TH1F>("hm4rec2OSveeno01t12g","M_{K0sK0s} t1+t2 < 0.06",massbins,0,5.);
  histosTH1F["hm4rec2OSveeno01t12h"] = fs->make<TH1F>("hm4rec2OSveeno01t12h","M_{K0sK0s} t1+t2 < 0.05",massbins,0,5.);
  // ...OR complementary
  histosTH1F["hm4rec2OSveeno01t12aC"] = fs->make<TH1F>("hm4rec2OSveeno01t12aC","M_{K0sK0s} t1+t2 > 0.4",massbins,0,5.);
  histosTH1F["hm4rec2OSveeno01t12bC"] = fs->make<TH1F>("hm4rec2OSveeno01t12bC","M_{K0sK0s} t1+t2 > 0.3",massbins,0,5.);
  histosTH1F["hm4rec2OSveeno01t12cC"] = fs->make<TH1F>("hm4rec2OSveeno01t12cC","M_{K0sK0s} t1+t2 > 0.2",massbins,0,5.);
  histosTH1F["hm4rec2OSveeno01t12dC"] = fs->make<TH1F>("hm4rec2OSveeno01t12dC","M_{K0sK0s} t1+t2 > 0.17",massbins,0,5.);
  histosTH1F["hm4rec2OSveeno01t12eC"] = fs->make<TH1F>("hm4rec2OSveeno01t12eC","M_{K0sK0s} t1+t2 > 0.1",massbins,0,5.);
  histosTH1F["hm4rec2OSveeno01t12fC"] = fs->make<TH1F>("hm4rec2OSveeno01t12fC","M_{K0sK0s} t1+t2 > 0.08",massbins,0,5.);
  histosTH1F["hm4rec2OSveeno01t12gC"] = fs->make<TH1F>("hm4rec2OSveeno01t12gC","M_{K0sK0s} t1+t2 > 0.06",massbins,0,5.);
  histosTH1F["hm4rec2OSveeno01t12hC"] = fs->make<TH1F>("hm4rec2OSveeno01t12hC","M_{K0sK0s} t1+t2 > 0.05",massbins,0,5.);
  // ...AND
  histosTH1F["hm4rec2OSveeno01t12aa"] = fs->make<TH1F>("hm4rec2OSveeno01t12aa","M_{K0sK0s} t1+t2 < 0.4",massbins,0,5.);
  histosTH1F["hm4rec2OSveeno01t12bb"] = fs->make<TH1F>("hm4rec2OSveeno01t12bb","M_{K0sK0s} t1+t2 < 0.3",massbins,0,5.);
  histosTH1F["hm4rec2OSveeno01t12cc"] = fs->make<TH1F>("hm4rec2OSveeno01t12cc","M_{K0sK0s} t1+t2 < 0.2",massbins,0,5.);
  histosTH1F["hm4rec2OSveeno01t12dd"] = fs->make<TH1F>("hm4rec2OSveeno01t12dd","M_{K0sK0s} t1+t2 < 0.17",massbins,0,5.);
  histosTH1F["hm4rec2OSveeno01t12ee"] = fs->make<TH1F>("hm4rec2OSveeno01t12ee","M_{K0sK0s} t1+t2 < 0.1",massbins,0,5.);
  histosTH1F["hm4rec2OSveeno01t12ff"] = fs->make<TH1F>("hm4rec2OSveeno01t12ff","M_{K0sK0s} t1+t2 < 0.08",massbins,0,5.);
  histosTH1F["hm4rec2OSveeno01t12gg"] = fs->make<TH1F>("hm4rec2OSveeno01t12gg","M_{K0sK0s} t1+t2 < 0.06",massbins,0,5.);
  histosTH1F["hm4rec2OSveeno01t12hh"] = fs->make<TH1F>("hm4rec2OSveeno01t12hh","M_{K0sK0s} t1+t2 < 0.05",massbins,0,5.);
  // ...AND complementary
  histosTH1F["hm4rec2OSveeno01t12aaC"] = fs->make<TH1F>("hm4rec2OSveeno01t12aaC","M_{K0sK0s} t1+t2 > 0.4",massbins,0,5.);
  histosTH1F["hm4rec2OSveeno01t12bbC"] = fs->make<TH1F>("hm4rec2OSveeno01t12bbC","M_{K0sK0s} t1+t2 > 0.3",massbins,0,5.);
  histosTH1F["hm4rec2OSveeno01t12ccC"] = fs->make<TH1F>("hm4rec2OSveeno01t12ccC","M_{K0sK0s} t1+t2 > 0.2",massbins,0,5.);
  histosTH1F["hm4rec2OSveeno01t12ddC"] = fs->make<TH1F>("hm4rec2OSveeno01t12ddC","M_{K0sK0s} t1+t2 > 0.17",massbins,0,5.);
  histosTH1F["hm4rec2OSveeno01t12eeC"] = fs->make<TH1F>("hm4rec2OSveeno01t12eeC","M_{K0sK0s} t1+t2 > 0.1",massbins,0,5.);
  histosTH1F["hm4rec2OSveeno01t12ffC"] = fs->make<TH1F>("hm4rec2OSveeno01t12ffC","M_{K0sK0s} t1+t2 > 0.08",massbins,0,5.);
  histosTH1F["hm4rec2OSveeno01t12ggC"] = fs->make<TH1F>("hm4rec2OSveeno01t12ggC","M_{K0sK0s} t1+t2 > 0.06",massbins,0,5.);
  histosTH1F["hm4rec2OSveeno01t12hhC"] = fs->make<TH1F>("hm4rec2OSveeno01t12hhC","M_{K0sK0s} t1+t2 > 0.05",massbins,0,5.);
  //
  // ...OR ...01 ...NO PID
  histosTH1F["hm4rec2OSveenopid01t12a"] = fs->make<TH1F>("hm4rec2OSveenopid01t12a","M_{K0sK0s} t1+t2 < 0.4",massbins,0,5.);
  histosTH1F["hm4rec2OSveenopid01t12b"] = fs->make<TH1F>("hm4rec2OSveenopid01t12b","M_{K0sK0s} t1+t2 < 0.3",massbins,0,5.);
  histosTH1F["hm4rec2OSveenopid01t12c"] = fs->make<TH1F>("hm4rec2OSveenopid01t12c","M_{K0sK0s} t1+t2 < 0.2",massbins,0,5.);
  histosTH1F["hm4rec2OSveenopid01t12d"] = fs->make<TH1F>("hm4rec2OSveenopid01t12d","M_{K0sK0s} t1+t2 < 0.17",massbins,0,5.);
  histosTH1F["hm4rec2OSveenopid01t12e"] = fs->make<TH1F>("hm4rec2OSveenopid01t12e","M_{K0sK0s} t1+t2 < 0.1",massbins,0,5.);
  histosTH1F["hm4rec2OSveenopid01t12f"] = fs->make<TH1F>("hm4rec2OSveenopid01t12f","M_{K0sK0s} t1+t2 < 0.08",massbins,0,5.);
  histosTH1F["hm4rec2OSveenopid01t12g"] = fs->make<TH1F>("hm4rec2OSveenopid01t12g","M_{K0sK0s} t1+t2 < 0.06",massbins,0,5.);
  histosTH1F["hm4rec2OSveenopid01t12h"] = fs->make<TH1F>("hm4rec2OSveenopid01t12h","M_{K0sK0s} t1+t2 < 0.05",massbins,0,5.);
  // ...OR complementary
  histosTH1F["hm4rec2OSveenopid01t12aC"] = fs->make<TH1F>("hm4rec2OSveenopid01t12aC","M_{K0sK0s} t1+t2 > 0.4",massbins,0,5.);
  histosTH1F["hm4rec2OSveenopid01t12bC"] = fs->make<TH1F>("hm4rec2OSveenopid01t12bC","M_{K0sK0s} t1+t2 > 0.3",massbins,0,5.);
  histosTH1F["hm4rec2OSveenopid01t12cC"] = fs->make<TH1F>("hm4rec2OSveenopid01t12cC","M_{K0sK0s} t1+t2 > 0.2",massbins,0,5.);
  histosTH1F["hm4rec2OSveenopid01t12dC"] = fs->make<TH1F>("hm4rec2OSveenopid01t12dC","M_{K0sK0s} t1+t2 > 0.17",massbins,0,5.);
  histosTH1F["hm4rec2OSveenopid01t12eC"] = fs->make<TH1F>("hm4rec2OSveenopid01t12eC","M_{K0sK0s} t1+t2 > 0.1",massbins,0,5.);
  histosTH1F["hm4rec2OSveenopid01t12fC"] = fs->make<TH1F>("hm4rec2OSveenopid01t12fC","M_{K0sK0s} t1+t2 > 0.08",massbins,0,5.);
  histosTH1F["hm4rec2OSveenopid01t12gC"] = fs->make<TH1F>("hm4rec2OSveenopid01t12gC","M_{K0sK0s} t1+t2 > 0.06",massbins,0,5.);
  histosTH1F["hm4rec2OSveenopid01t12hC"] = fs->make<TH1F>("hm4rec2OSveenopid01t12hC","M_{K0sK0s} t1+t2 > 0.05",massbins,0,5.);
  // ...AND
  histosTH1F["hm4rec2OSveenopid01t12aa"] = fs->make<TH1F>("hm4rec2OSveenopid01t12aa","M_{K0sK0s} t1+t2 < 0.4",massbins,0,5.);
  histosTH1F["hm4rec2OSveenopid01t12bb"] = fs->make<TH1F>("hm4rec2OSveenopid01t12bb","M_{K0sK0s} t1+t2 < 0.3",massbins,0,5.);
  histosTH1F["hm4rec2OSveenopid01t12cc"] = fs->make<TH1F>("hm4rec2OSveenopid01t12cc","M_{K0sK0s} t1+t2 < 0.2",massbins,0,5.);
  histosTH1F["hm4rec2OSveenopid01t12dd"] = fs->make<TH1F>("hm4rec2OSveenopid01t12dd","M_{K0sK0s} t1+t2 < 0.17",massbins,0,5.);
  histosTH1F["hm4rec2OSveenopid01t12ee"] = fs->make<TH1F>("hm4rec2OSveenopid01t12ee","M_{K0sK0s} t1+t2 < 0.1",massbins,0,5.);
  histosTH1F["hm4rec2OSveenopid01t12ff"] = fs->make<TH1F>("hm4rec2OSveenopid01t12ff","M_{K0sK0s} t1+t2 < 0.08",massbins,0,5.);
  histosTH1F["hm4rec2OSveenopid01t12gg"] = fs->make<TH1F>("hm4rec2OSveenopid01t12gg","M_{K0sK0s} t1+t2 < 0.06",massbins,0,5.);
  histosTH1F["hm4rec2OSveenopid01t12hh"] = fs->make<TH1F>("hm4rec2OSveenopid01t12hh","M_{K0sK0s} t1+t2 < 0.05",massbins,0,5.);
  // ...AND complementary
  histosTH1F["hm4rec2OSveenopid01t12aaC"] = fs->make<TH1F>("hm4rec2OSveenopid01t12aaC","M_{K0sK0s} t1+t2 > 0.4",massbins,0,5.);
  histosTH1F["hm4rec2OSveenopid01t12bbC"] = fs->make<TH1F>("hm4rec2OSveenopid01t12bbC","M_{K0sK0s} t1+t2 > 0.3",massbins,0,5.);
  histosTH1F["hm4rec2OSveenopid01t12ccC"] = fs->make<TH1F>("hm4rec2OSveenopid01t12ccC","M_{K0sK0s} t1+t2 > 0.2",massbins,0,5.);
  histosTH1F["hm4rec2OSveenopid01t12ddC"] = fs->make<TH1F>("hm4rec2OSveenopid01t12ddC","M_{K0sK0s} t1+t2 > 0.17",massbins,0,5.);
  histosTH1F["hm4rec2OSveenopid01t12eeC"] = fs->make<TH1F>("hm4rec2OSveenopid01t12eeC","M_{K0sK0s} t1+t2 > 0.1",massbins,0,5.);
  histosTH1F["hm4rec2OSveenopid01t12ffC"] = fs->make<TH1F>("hm4rec2OSveenopid01t12ffC","M_{K0sK0s} t1+t2 > 0.08",massbins,0,5.);
  histosTH1F["hm4rec2OSveenopid01t12ggC"] = fs->make<TH1F>("hm4rec2OSveenopid01t12ggC","M_{K0sK0s} t1+t2 > 0.06",massbins,0,5.);
  histosTH1F["hm4rec2OSveenopid01t12hhC"] = fs->make<TH1F>("hm4rec2OSveenopid01t12hhC","M_{K0sK0s} t1+t2 > 0.05",massbins,0,5.);
  //
  
  //
  //**histosTH1F["hm4rec2OSvee11a"] = fs->make<TH1F>("hm4rec2OSvee11a","M_{K#pi} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OSvee11b"] = fs->make<TH1F>("hm4rec2OSvee11b","M_{K#pi} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OSvee11c"] = fs->make<TH1F>("hm4rec2OSvee11c","M_{K#pi} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OSvee11d"] = fs->make<TH1F>("hm4rec2OSvee11d","M_{K#pi} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OSvee11e"] = fs->make<TH1F>("hm4rec2OSvee11e","M_{K#pi} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OSvee11f"] = fs->make<TH1F>("hm4rec2OSvee11f","M_{K#pi} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OSvee11g"] = fs->make<TH1F>("hm4rec2OSvee11g","M_{K#pi} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OSvee11h"] = fs->make<TH1F>("hm4rec2OSvee11h","M_{K#pi} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OSvee11i"] = fs->make<TH1F>("hm4rec2OSvee11i","M_{K#pi} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OSvee11j"] = fs->make<TH1F>("hm4rec2OSvee11j","M_{K#pi} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OSvee11k"] = fs->make<TH1F>("hm4rec2OSvee11k","M_{K#pi} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OSvee11m"] = fs->make<TH1F>("hm4rec2OSvee11m","M_{K#pi} OS",massbins,0,5.);
  //
  //**histosTH1F["hm4rec2OSveeno11a"] = fs->make<TH1F>("hm4rec2OSveeno11a","M_{K#pi} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OSveeno11b"] = fs->make<TH1F>("hm4rec2OSveeno11b","M_{K#pi} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OSveeno11c"] = fs->make<TH1F>("hm4rec2OSveeno11c","M_{K#pi} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OSveeno11d"] = fs->make<TH1F>("hm4rec2OSveeno11d","M_{K#pi} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OSveeno11e"] = fs->make<TH1F>("hm4rec2OSveeno11e","M_{K#pi} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OSveeno11f"] = fs->make<TH1F>("hm4rec2OSveeno11f","M_{K#pi} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OSveeno11g"] = fs->make<TH1F>("hm4rec2OSveeno11g","M_{K#pi} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OSveeno11h"] = fs->make<TH1F>("hm4rec2OSveeno11h","M_{K#pi} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OSveeno11i"] = fs->make<TH1F>("hm4rec2OSveeno11i","M_{K#pi} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OSveeno11j"] = fs->make<TH1F>("hm4rec2OSveeno11j","M_{K#pi} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OSveeno11k"] = fs->make<TH1F>("hm4rec2OSveeno11k","M_{K#pi} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OSveeno11m"] = fs->make<TH1F>("hm4rec2OSveeno11m","M_{K#pi} OS",massbins,0,5.);

  //**histosTH1F["hm4rec2OSvee9"] = fs->make<TH1F>("hm4rec2OSvee9","M_{4#pi} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OSvee90"] = fs->make<TH1F>("hm4rec2OSvee90","M_{4#pi} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OSvee91"] = fs->make<TH1F>("hm4rec2OSvee91","M_{4#pi} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OSvee92"] = fs->make<TH1F>("hm4rec2OSvee92","M_{4#pi} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OSvtx0"] = fs->make<TH1F>("hm4rec2OSvtx0","M_{4#pi} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OSvtx01"] = fs->make<TH1F>("hm4rec2OSvtx01","M_{4#pi} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OSvtx02"] = fs->make<TH1F>("hm4rec2OSvtx02","M_{4#pi} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OSvtx11"] = fs->make<TH1F>("hm4rec2OSvtx11","M_{4#pi} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OSvtx1"] = fs->make<TH1F>("hm4rec2OSvtx1","M_{4#pi} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OSvtx2"] = fs->make<TH1F>("hm4rec2OSvtx2","M_{4#pi} OS",massbins,0,5.);
  
  histosTH1F["hm4rec2OS2"] = fs->make<TH1F>("hm4rec2OS2","M_{4#pi} OS",massbins,0,5.);
  histosTH1F["hm4rec2OS2pid"] = fs->make<TH1F>("hm4rec2OS2pid","M_{4#pi} OS pid=pion",massbins,0,5.);
  histosTH1F["hm4rec2OS2pidv0"] = fs->make<TH1F>("hm4rec2OS2pidv0","M_{4#pi} OS pidV0=pion",massbins,0,5.);
  histosTH1F["hm4rec2OS2nov0"] = fs->make<TH1F>("hm4rec2OS2nov0","M_{4#pi} OS no V0",massbins,0,5.);
  histosTH1F["hm4rec2OS2veeno10"] = fs->make<TH1F>("hm4rec2OS2veeno10","M_{4#pi} OS type:10",massbins,0,5.);
  histosTH1F["hm4rec2OS2k"] = fs->make<TH1F>("hm4rec2OS2k","M_{4K} OS",massbins,0,5.);
  histosTH1F["hm4rec2OS2kpid"] = fs->make<TH1F>("hm4rec2OS2kpid","M_{4K} OS pid=kaons",massbins,0,5.);
  // 12 34 13 24  ...for now
  //**histosTH1F["hm4rec2OS_pipi"] = fs->make<TH1F>("hm4rec2OS_pipi","M_{#pi#pi} OS",massbins,0,5.);

  //**histosTH1F["hm4rec2OS_pi1pi2"] = fs->make<TH1F>("hm4rec2OS_pi1pi2","M_{#pi_{1}#pi_{2}} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi3pi4"] = fs->make<TH1F>("hm4rec2OS_pi3pi4","M_{#pi_{3}#pi_{4}} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi1pi3"] = fs->make<TH1F>("hm4rec2OS_pi1pi3","M_{#pi_{1}#pi_{3}} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi2pi4"] = fs->make<TH1F>("hm4rec2OS_pi2pi4","M_{#pi_{2}#pi_{4}} OS",massbins,0,5.);
  //

  // 2 primary vertices
  histosTH1F["hm4rec2OS3"] = fs->make<TH1F>("hm4rec2OS3","M_{4#pi} OS 2 primary vertices",massbins,0,5.);
  histosTH1F["hm4rec2OS3nov0"] = fs->make<TH1F>("hm4rec2OS3nov0","M_{4#pi} OS no V0 2 primary vertices",massbins,0,5.);
  histosTH1F["hm4rec2OS3veeno20"] = fs->make<TH1F>("hm4rec2OS3veeno20","M_{4#pi} OS type:20 2 primary vertices",massbins,0,5.);

  //...mass selection
  //**histosTH1F["hm4rec2OS_pi1pi2m"] = fs->make<TH1F>("hm4rec2OS_pi1pi2m","M_{#pi_{1}#pi_{2}} K0s mass selection",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi3pi4m"] = fs->make<TH1F>("hm4rec2OS_pi3pi4m","M_{#pi_{3}#pi_{4}} K0s mass selection",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi1pi3m"] = fs->make<TH1F>("hm4rec2OS_pi1pi3m","M_{#pi_{1}#pi_{3}} K0s mass selection",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi2pi4m"] = fs->make<TH1F>("hm4rec2OS_pi2pi4m","M_{#pi_{2}#pi_{4}} K0s mass selection",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi1pi4m"] = fs->make<TH1F>("hm4rec2OS_pi1pi4m","M_{#pi_{1}#pi_{4}} K0s mass selection",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi2pi3m"] = fs->make<TH1F>("hm4rec2OS_pi2pi3m","M_{#pi_{2}#pi_{3}} K0s mass selection",massbins,0,5.);
  //
  //**histosTH1F["hm4rec2OS_pi1pi2t"] = fs->make<TH1F>("hm4rec2OS_pi1pi2t","M_{#pi_{1}#pi_{2}} K0s",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi3pi4t"] = fs->make<TH1F>("hm4rec2OS_pi3pi4t","M_{#pi_{3}#pi_{4}} K0s",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi1pi3t"] = fs->make<TH1F>("hm4rec2OS_pi1pi3t","M_{#pi_{1}#pi_{3}} K0s",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi2pi4t"] = fs->make<TH1F>("hm4rec2OS_pi2pi4t","M_{#pi_{2}#pi_{4}} K0s",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi1pi4t"] = fs->make<TH1F>("hm4rec2OS_pi1pi4t","M_{#pi_{1}#pi_{4}} K0s",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi2pi3t"] = fs->make<TH1F>("hm4rec2OS_pi2pi3t","M_{#pi_{2}#pi_{3}} K0s",massbins,0,5.);
  //
  //**histosTH1F["hm4rec2OS_pi1pi2m00"] = fs->make<TH1F>("hm4rec2OS_pi1pi2m00","M_{#pi_{1}#pi_{2}} K0s mass selection",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi3pi4m00"] = fs->make<TH1F>("hm4rec2OS_pi3pi4m00","M_{#pi_{3}#pi_{4}} K0s mass selection",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi1pi3m00"] = fs->make<TH1F>("hm4rec2OS_pi1pi3m00","M_{#pi_{1}#pi_{3}} K0s mass selection",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi2pi4m00"] = fs->make<TH1F>("hm4rec2OS_pi2pi4m00","M_{#pi_{2}#pi_{4}} K0s mass selection",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi1pi4m00"] = fs->make<TH1F>("hm4rec2OS_pi1pi4m00","M_{#pi_{1}#pi_{4}} K0s mass selection",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi2pi3m00"] = fs->make<TH1F>("hm4rec2OS_pi2pi3m00","M_{#pi_{2}#pi_{3}} K0s mass selection",massbins,0,5.);
  //
  //**histosTH1F["hm4rec2OS_pi1pi2m000"] = fs->make<TH1F>("hm4rec2OS_pi1pi2m000","M_{#pi_{1}#pi_{2}} K0s mass selection",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi3pi4m000"] = fs->make<TH1F>("hm4rec2OS_pi3pi4m000","M_{#pi_{3}#pi_{4}} K0s mass selection",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi1pi3m000"] = fs->make<TH1F>("hm4rec2OS_pi1pi3m000","M_{#pi_{1}#pi_{3}} K0s mass selection",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi2pi4m000"] = fs->make<TH1F>("hm4rec2OS_pi2pi4m000","M_{#pi_{2}#pi_{4}} K0s mass selection",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi1pi4m000"] = fs->make<TH1F>("hm4rec2OS_pi1pi4m000","M_{#pi_{1}#pi_{4}} K0s mass selection",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi2pi3m000"] = fs->make<TH1F>("hm4rec2OS_pi2pi3m000","M_{#pi_{2}#pi_{3}} K0s mass selection",massbins,0,5.);
  //
  //**histosTH1F["hm4rec2OS_pi1pi2m04"] = fs->make<TH1F>("hm4rec2OS_pi1pi2m04","M_{#pi_{1}#pi_{2}} K0s mass selection",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi3pi4m04"] = fs->make<TH1F>("hm4rec2OS_pi3pi4m04","M_{#pi_{3}#pi_{4}} K0s mass selection",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi1pi3m04"] = fs->make<TH1F>("hm4rec2OS_pi1pi3m04","M_{#pi_{1}#pi_{3}} K0s mass selection",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi2pi4m04"] = fs->make<TH1F>("hm4rec2OS_pi2pi4m04","M_{#pi_{2}#pi_{4}} K0s mass selection",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi1pi4m04"] = fs->make<TH1F>("hm4rec2OS_pi1pi4m04","M_{#pi_{1}#pi_{4}} K0s mass selection",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi2pi3m04"] = fs->make<TH1F>("hm4rec2OS_pi2pi3m04","M_{#pi_{2}#pi_{3}} K0s mass selection",massbins,0,5.);
  //
  //**histosTH1F["hm4rec2OS_pi1pi2t04"] = fs->make<TH1F>("hm4rec2OS_pi1pi2t04","M_{#pi_{1}#pi_{2}} K0s",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi3pi4t04"] = fs->make<TH1F>("hm4rec2OS_pi3pi4t04","M_{#pi_{3}#pi_{4}} K0s",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi1pi3t04"] = fs->make<TH1F>("hm4rec2OS_pi1pi3t04","M_{#pi_{1}#pi_{3}} K0s",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi2pi4t04"] = fs->make<TH1F>("hm4rec2OS_pi2pi4t04","M_{#pi_{2}#pi_{4}} K0s",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi1pi4t04"] = fs->make<TH1F>("hm4rec2OS_pi1pi4t04","M_{#pi_{1}#pi_{4}} K0s",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi2pi3t04"] = fs->make<TH1F>("hm4rec2OS_pi2pi3t04","M_{#pi_{2}#pi_{3}} K0s",massbins,0,5.);
  //
  histosTH1F["hm4rec2OS_pi1pi2m05"] = fs->make<TH1F>("hm4rec2OS_pi1pi2m05","M_{#pi_{1}#pi_{2}} K0s mass selection",massbins,0,5.);
  histosTH1F["hm4rec2OS_pi3pi4m05"] = fs->make<TH1F>("hm4rec2OS_pi3pi4m05","M_{#pi_{3}#pi_{4}} K0s mass selection",massbins,0,5.);
  histosTH1F["hm4rec2OS_pi1pi3m05"] = fs->make<TH1F>("hm4rec2OS_pi1pi3m05","M_{#pi_{1}#pi_{3}} K0s mass selection",massbins,0,5.);
  histosTH1F["hm4rec2OS_pi2pi4m05"] = fs->make<TH1F>("hm4rec2OS_pi2pi4m05","M_{#pi_{2}#pi_{4}} K0s mass selection",massbins,0,5.);
  histosTH1F["hm4rec2OS_pi1pi4m05"] = fs->make<TH1F>("hm4rec2OS_pi1pi4m05","M_{#pi_{1}#pi_{4}} K0s mass selection",massbins,0,5.);
  histosTH1F["hm4rec2OS_pi2pi3m05"] = fs->make<TH1F>("hm4rec2OS_pi2pi3m05","M_{#pi_{2}#pi_{3}} K0s mass selection",massbins,0,5.);
  //
  histosTH1F["hm4rec2OS_pi1pi2t05"] = fs->make<TH1F>("hm4rec2OS_pi1pi2t05","M_{#pi_{1}#pi_{2}} K0s",massbins,0,5.);
  histosTH1F["hm4rec2OS_pi3pi4t05"] = fs->make<TH1F>("hm4rec2OS_pi3pi4t05","M_{#pi_{3}#pi_{4}} K0s",massbins,0,5.);
  histosTH1F["hm4rec2OS_pi1pi3t05"] = fs->make<TH1F>("hm4rec2OS_pi1pi3t05","M_{#pi_{1}#pi_{3}} K0s",massbins,0,5.);
  histosTH1F["hm4rec2OS_pi2pi4t05"] = fs->make<TH1F>("hm4rec2OS_pi2pi4t05","M_{#pi_{2}#pi_{4}} K0s",massbins,0,5.);
  histosTH1F["hm4rec2OS_pi1pi4t05"] = fs->make<TH1F>("hm4rec2OS_pi1pi4t05","M_{#pi_{1}#pi_{4}} K0s",massbins,0,5.);
  histosTH1F["hm4rec2OS_pi2pi3t05"] = fs->make<TH1F>("hm4rec2OS_pi2pi3t05","M_{#pi_{2}#pi_{3}} K0s",massbins,0,5.);
  //
  //**histosTH1F["hm4rec2OS_pi1pi2m06"] = fs->make<TH1F>("hm4rec2OS_pi1pi2m06","M_{#pi_{1}#pi_{2}} K0s mass selection",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi3pi4m06"] = fs->make<TH1F>("hm4rec2OS_pi3pi4m06","M_{#pi_{3}#pi_{4}} K0s mass selection",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi1pi3m06"] = fs->make<TH1F>("hm4rec2OS_pi1pi3m06","M_{#pi_{1}#pi_{3}} K0s mass selection",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi2pi4m06"] = fs->make<TH1F>("hm4rec2OS_pi2pi4m06","M_{#pi_{2}#pi_{4}} K0s mass selection",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi1pi4m06"] = fs->make<TH1F>("hm4rec2OS_pi1pi4m06","M_{#pi_{1}#pi_{4}} K0s mass selection",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi2pi3m06"] = fs->make<TH1F>("hm4rec2OS_pi2pi3m06","M_{#pi_{2}#pi_{3}} K0s mass selection",massbins,0,5.);
  //
  //..cut k1
  histosTH1F["hm4rec2OS_k1k2m"] = fs->make<TH1F>("hm4rec2OS_k1k2m","M_{K_{1}K_{2}} K^{+}K^{-}",massbins,0,5.);
  histosTH1F["hm4rec2OS_k3k4m"] = fs->make<TH1F>("hm4rec2OS_k3k4m","M_{K_{3}K_{4}} K^{+}K^{-}",massbins,0,5.);
  histosTH1F["hm4rec2OS_k1k3m"] = fs->make<TH1F>("hm4rec2OS_k1k3m","M_{K_{1}K_{3}} K^{+}K^{-}",massbins,0,5.);
  histosTH1F["hm4rec2OS_k2k4m"] = fs->make<TH1F>("hm4rec2OS_k2k4m","M_{K_{2}K_{4}} K^{+}K^{-}",massbins,0,5.);
  histosTH1F["hm4rec2OS_k1k4m"] = fs->make<TH1F>("hm4rec2OS_k1k4m","M_{K_{1}K_{4}} K^{+}K^{-}",massbins,0,5.);
  histosTH1F["hm4rec2OS_k2k3m"] = fs->make<TH1F>("hm4rec2OS_k2k3m","M_{K_{2}K_{3}} K^{+}K^{-}",massbins,0,5.);
  
  //...v2
  //**histosTH1F["hm4rec2OS_pi1pi2v2"] = fs->make<TH1F>("hm4rec2OS_pi1pi2v2","M_{#pi_{1}#pi_{2}} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi3pi4v2"] = fs->make<TH1F>("hm4rec2OS_pi3pi4v2","M_{#pi_{3}#pi_{4}} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi1pi3v2"] = fs->make<TH1F>("hm4rec2OS_pi1pi3v2","M_{#pi_{1}#pi_{3}} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi2pi4v2"] = fs->make<TH1F>("hm4rec2OS_pi2pi4v2","M_{#pi_{2}#pi_{4}} OS",massbins,0,5.);
  //...2dim
  //***histosTH2F["hm4dim2OS_pi1pi2_pi3pi4"] = fs->make<TH2F>("hm4dim2OS_pi1pi2_pi3pi4","M_{#pi_{1}#pi_{2}} vs M_{#pi_{3}#pi_{4}} OS",massbins,0,5.,massbins,0,5.);
  //***histosTH2F["hm4dim2OS_pi1pi3_pi2pi4"] = fs->make<TH2F>("hm4dim2OS_pi1pi3_pi2pi4","M_{#pi_{1}#pi_{3}} vs M_{#pi_{2}#pi_{4}} OS",massbins,0,5.,massbins,0,5.);
  //...v2
  //***histosTH2F["hm4dim2OS_pi1pi2_pi3pi4v2"] = fs->make<TH2F>("hm4dim2OS_pi1pi2_pi3pi4v2","M_{#pi_{1}#pi_{2}} vs M_{#pi_{3}#pi_{4}} OS",massbins,0,5.,massbins,0,5.);
  //***histosTH2F["hm4dim2OS_pi1pi3_pi2pi4v2"] = fs->make<TH2F>("hm4dim2OS_pi1pi3_pi2pi4v2","M_{#pi_{1}#pi_{3}} vs M_{#pi_{2}#pi_{4}} OS",massbins,0,5.,massbins,0,5.);
  
  //...testing vee9  ........cut 9
  //**histosTH1F["hm4rec2OS_pi1pi2vee9"] = fs->make<TH1F>("hm4rec2OS_pi1pi2vee9","M_{#pi_{1}#pi_{2}} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi3pi4vee9"] = fs->make<TH1F>("hm4rec2OS_pi3pi4vee9","M_{#pi_{3}#pi_{4}} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi1pi3vee9"] = fs->make<TH1F>("hm4rec2OS_pi1pi3vee9","M_{#pi_{1}#pi_{3}} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi2pi4vee9"] = fs->make<TH1F>("hm4rec2OS_pi2pi4vee9","M_{#pi_{2}#pi_{4}} OS",massbins,0,5.);
  //
  //**histosTH1F["hm4rec2OS_pi1pi2vee90"] = fs->make<TH1F>("hm4rec2OS_pi1pi2vee90","M_{#pi_{1}#pi_{2}} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi3pi4vee90"] = fs->make<TH1F>("hm4rec2OS_pi3pi4vee90","M_{#pi_{3}#pi_{4}} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi1pi3vee90"] = fs->make<TH1F>("hm4rec2OS_pi1pi3vee90","M_{#pi_{1}#pi_{3}} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi2pi4vee90"] = fs->make<TH1F>("hm4rec2OS_pi2pi4vee90","M_{#pi_{2}#pi_{4}} OS",massbins,0,5.);
  //
  //**histosTH1F["hm4rec2OS_pi1pi2vee91"] = fs->make<TH1F>("hm4rec2OS_pi1pi2vee91","M_{#pi_{1}#pi_{2}} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi3pi4vee91"] = fs->make<TH1F>("hm4rec2OS_pi3pi4vee91","M_{#pi_{3}#pi_{4}} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi1pi3vee91"] = fs->make<TH1F>("hm4rec2OS_pi1pi3vee91","M_{#pi_{1}#pi_{3}} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi2pi4vee91"] = fs->make<TH1F>("hm4rec2OS_pi2pi4vee91","M_{#pi_{2}#pi_{4}} OS",massbins,0,5.);
  //
  //**histosTH1F["hm4rec2OS_pi1pi2vee92"] = fs->make<TH1F>("hm4rec2OS_pi1pi2vee92","M_{#pi_{1}#pi_{2}} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi3pi4vee92"] = fs->make<TH1F>("hm4rec2OS_pi3pi4vee92","M_{#pi_{3}#pi_{4}} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi1pi3vee92"] = fs->make<TH1F>("hm4rec2OS_pi1pi3vee92","M_{#pi_{1}#pi_{3}} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi2pi4vee92"] = fs->make<TH1F>("hm4rec2OS_pi2pi4vee92","M_{#pi_{2}#pi_{4}} OS",massbins,0,5.);
  //
  //**histosTH1F["hm4rec2OS_pi1pi2vtx0"] = fs->make<TH1F>("hm4rec2OS_pi1pi2vtx0","M_{#pi_{1}#pi_{2}} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi3pi4vtx0"] = fs->make<TH1F>("hm4rec2OS_pi3pi4vtx0","M_{#pi_{3}#pi_{4}} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi1pi3vtx0"] = fs->make<TH1F>("hm4rec2OS_pi1pi3vtx0","M_{#pi_{1}#pi_{3}} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi2pi4vtx0"] = fs->make<TH1F>("hm4rec2OS_pi2pi4vtx0","M_{#pi_{2}#pi_{4}} OS",massbins,0,5.);
   //
  //**histosTH1F["hm4rec2OS_pi1pi2vtx01"] = fs->make<TH1F>("hm4rec2OS_pi1pi2vtx01","M_{#pi_{1}#pi_{2}} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi3pi4vtx01"] = fs->make<TH1F>("hm4rec2OS_pi3pi4vtx01","M_{#pi_{3}#pi_{4}} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi1pi3vtx01"] = fs->make<TH1F>("hm4rec2OS_pi1pi3vtx01","M_{#pi_{1}#pi_{3}} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi2pi4vtx01"] = fs->make<TH1F>("hm4rec2OS_pi2pi4vtx01","M_{#pi_{2}#pi_{4}} OS",massbins,0,5.);
   //
  //**histosTH1F["hm4rec2OS_pi1pi2vtx02"] = fs->make<TH1F>("hm4rec2OS_pi1pi2vtx02","M_{#pi_{1}#pi_{2}} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi3pi4vtx02"] = fs->make<TH1F>("hm4rec2OS_pi3pi4vtx02","M_{#pi_{3}#pi_{4}} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi1pi3vtx02"] = fs->make<TH1F>("hm4rec2OS_pi1pi3vtx02","M_{#pi_{1}#pi_{3}} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi2pi4vtx02"] = fs->make<TH1F>("hm4rec2OS_pi2pi4vtx02","M_{#pi_{2}#pi_{4}} OS",massbins,0,5.);
   //
  //**histosTH1F["hm4rec2OS_pi1pi2vtx11"] = fs->make<TH1F>("hm4rec2OS_pi1pi2vtx11","M_{#pi_{1}#pi_{2}} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi3pi4vtx11"] = fs->make<TH1F>("hm4rec2OS_pi3pi4vtx11","M_{#pi_{3}#pi_{4}} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi1pi3vtx11"] = fs->make<TH1F>("hm4rec2OS_pi1pi3vtx11","M_{#pi_{1}#pi_{3}} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi2pi4vtx11"] = fs->make<TH1F>("hm4rec2OS_pi2pi4vtx11","M_{#pi_{2}#pi_{4}} OS",massbins,0,5.);
  //
  //**histosTH1F["hm4rec2OS_pi1pi2vtx1"] = fs->make<TH1F>("hm4rec2OS_pi1pi2vtx1","M_{#pi_{1}#pi_{2}} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi3pi4vtx1"] = fs->make<TH1F>("hm4rec2OS_pi3pi4vtx1","M_{#pi_{3}#pi_{4}} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi1pi3vtx1"] = fs->make<TH1F>("hm4rec2OS_pi1pi3vtx1","M_{#pi_{1}#pi_{3}} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi2pi4vtx1"] = fs->make<TH1F>("hm4rec2OS_pi2pi4vtx1","M_{#pi_{2}#pi_{4}} OS",massbins,0,5.);
  //
  //**histosTH1F["hm4rec2OS_pi1pi2vtx2"] = fs->make<TH1F>("hm4rec2OS_pi1pi2vtx2","M_{#pi_{1}#pi_{2}} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi3pi4vtx2"] = fs->make<TH1F>("hm4rec2OS_pi3pi4vtx2","M_{#pi_{3}#pi_{4}} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi1pi3vtx2"] = fs->make<TH1F>("hm4rec2OS_pi1pi3vtx2","M_{#pi_{1}#pi_{3}} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi2pi4vtx2"] = fs->make<TH1F>("hm4rec2OS_pi2pi4vtx2","M_{#pi_{2}#pi_{4}} OS",massbins,0,5.);
  //
  
   //...A
   //**histosTH1F["hm4rec2OS_pi1pi2vee11"] = fs->make<TH1F>("hm4rec2OS_pi1pi2vee11","M_{#pi_{1}#pi_{2}} OS",massbins,0,5.);
   //**histosTH1F["hm4rec2OS_pi3k4vee11"] = fs->make<TH1F>("hm4rec2OS_pi3k4vee11","M_{#pi_{3}K_{4}} OS",massbins,0,5.);
   //  
   //**histosTH1F["hm4rec2OS_pi1pi3vee11"] = fs->make<TH1F>("hm4rec2OS_pi1pi3vee11","M_{#pi_{1}#pi_{3}} OS",massbins,0,5.);
   //**histosTH1F["hm4rec2OS_pi2k4vee11"] = fs->make<TH1F>("hm4rec2OS_pi2k4vee11","M_{#pi_{2}K_{4}} OS",massbins,0,5.);
   //
   //**histosTH1F["hm4rec2OS_pi2pi3vee11"] = fs->make<TH1F>("hm4rec2OS_pi2pi3vee11","M_{#pi_{2}#pi_{3}} OS",massbins,0,5.);
   //**histosTH1F["hm4rec2OS_pi1k4vee11"] = fs->make<TH1F>("hm4rec2OS_pi1k4vee11","M_{#pi_{1}K_{4}} OS",massbins,0,5.);
   //...B  
   //**histosTH1F["hm4rec2OS_k3pi4vee11"] = fs->make<TH1F>("hm4rec2OS_k3pi4vee11","M_{K_{3}#pi_{4}} OS",massbins,0,5.);
   //
   //**histosTH1F["hm4rec2OS_pi1pi4vee11"] = fs->make<TH1F>("hm4rec2OS_pi1pi4vee11","M_{#pi_{1}#pi_{4}} OS",massbins,0,5.);
   //**histosTH1F["hm4rec2OS_k3pi2vee11"] = fs->make<TH1F>("hm4rec2OS_k3pi2vee11","M_{K_{3}#pi_{2}} OS",massbins,0,5.);
   //
   //**histosTH1F["hm4rec2OS_pi2pi4vee11"] = fs->make<TH1F>("hm4rec2OS_pi2pi4vee11","M_{#pi_{2}#pi_{4}} OS",massbins,0,5.);
   //**histosTH1F["hm4rec2OS_k3pi1vee11"] = fs->make<TH1F>("hm4rec2OS_k3pi1vee11","M_{K_{3}#pi_{1}} OS",massbins,0,5.);
   //...C
   //**histosTH1F["hm4rec2OS_pi1k2vee11"] = fs->make<TH1F>("hm4rec2OS_pi1k2vee11","M_{#pi_{1}K_{2}} OS",massbins,0,5.);
   //**histosTH1F["hm4rec2OS_pi3pi4vee11"] = fs->make<TH1F>("hm4rec2OS_pi3pi4vee11","M_{#pi_{3}#pi_{4}} OS",massbins,0,5.);
   //
   //**histosTH1F["hm4rec2OS_pi3k2vee11"] = fs->make<TH1F>("hm4rec2OS_pi3k2vee11","M_{#pi_{3}K_{2}} OS",massbins,0,5.);
   //
   //**histosTH1F["hm4rec2OS_pi4k2vee11"] = fs->make<TH1F>("hm4rec2OS_pi4k2vee11","M_{#pi_{4}K_{2}} OS",massbins,0,5.);
   //...D
   //**histosTH1F["hm4rec2OS_k1pi2vee11"] = fs->make<TH1F>("hm4rec2OS_k1pi2vee11","M_{K_{1}#pi_{2}} OS",massbins,0,5.);
   //
   //**histosTH1F["hm4rec2OS_k1pi3vee11"] = fs->make<TH1F>("hm4rec2OS_k1pi3vee11","M_{K_{1}#pi_{3}} OS",massbins,0,5.);
   //
   //**histosTH1F["hm4rec2OS_k1pi4vee11"] = fs->make<TH1F>("hm4rec2OS_k1pi4vee11","M_{K_{1}#pi_{4}} OS",massbins,0,5.);
   //
  
   //...A...no PID  k4
   histosTH1F["hm4rec2OS_pi1pi2k4veeno11"] = fs->make<TH1F>("hm4rec2OS_pi1pi2k4veeno11","M_{#pi_{1}#pi_{2}} OS",massbins,0,5.);
   histosTH1F["hm4rec2OS_pi3k4veeno11"] = fs->make<TH1F>("hm4rec2OS_pi3k4veeno11","M_{#pi_{3}K_{4}} OS",massbins,0,5.);
   //  
   histosTH1F["hm4rec2OS_pi1pi3k4veeno11"] = fs->make<TH1F>("hm4rec2OS_pi1pi3k4veeno11","M_{#pi_{1}#pi_{3}} OS",massbins,0,5.);
   histosTH1F["hm4rec2OS_pi2k4veeno11"] = fs->make<TH1F>("hm4rec2OS_pi2k4veeno11","M_{#pi_{2}K_{4}} OS",massbins,0,5.);
   //
   histosTH1F["hm4rec2OS_pi2pi3k4veeno11"] = fs->make<TH1F>("hm4rec2OS_pi2pi3k4veeno11","M_{#pi_{2}#pi_{3}} OS",massbins,0,5.);
   histosTH1F["hm4rec2OS_pi1k4veeno11"] = fs->make<TH1F>("hm4rec2OS_pi1k4veeno11","M_{#pi_{1}K_{4}} OS",massbins,0,5.);
   //
   //...B...no PID  k3
   histosTH1F["hm4rec2OS_k3pi4veeno11"] = fs->make<TH1F>("hm4rec2OS_k3pi4veeno11","M_{K_{3}#pi_{4}} OS",massbins,0,5.);
   histosTH1F["hm4rec2OS_pi1pi2k3veeno11"] = fs->make<TH1F>("hm4rec2OS_pi1pi2k3veeno11","M_{#pi_{1}#pi_{2}} OS",massbins,0,5.);
   //
   histosTH1F["hm4rec2OS_pi1pi4k3veeno11"] = fs->make<TH1F>("hm4rec2OS_pi1pi4k3veeno11","M_{#pi_{1}#pi_{4}} OS",massbins,0,5.);
   histosTH1F["hm4rec2OS_k3pi2veeno11"] = fs->make<TH1F>("hm4rec2OS_k3pi2veeno11","M_{K_{3}#pi_{2}} OS",massbins,0,5.);
   //
   histosTH1F["hm4rec2OS_pi2pi4k3veeno11"] = fs->make<TH1F>("hm4rec2OS_pi2pi4k3veeno11","M_{#pi_{2}#pi_{4}} OS",massbins,0,5.);
   histosTH1F["hm4rec2OS_k3pi1veeno11"] = fs->make<TH1F>("hm4rec2OS_k3pi1veeno11","M_{K_{3}#pi_{1}} OS",massbins,0,5.);
   //
   //...C...no PID  k2
   histosTH1F["hm4rec2OS_pi1k2veeno11"] = fs->make<TH1F>("hm4rec2OS_pi1k2veeno11","M_{#pi_{1}K_{2}} OS",massbins,0,5.);
   histosTH1F["hm4rec2OS_pi3pi4k2veeno11"] = fs->make<TH1F>("hm4rec2OS_pi3pi4k2veeno11","M_{#pi_{3}#pi_{4}} OS",massbins,0,5.);
   //
   histosTH1F["hm4rec2OS_pi3k2veeno11"] = fs->make<TH1F>("hm4rec2OS_pi3k2veeno11","M_{#pi_{3}K_{2}} OS",massbins,0,5.);
   histosTH1F["hm4rec2OS_pi1pi4k2veeno11"] = fs->make<TH1F>("hm4rec2OS_pi1pi4k2veeno11","M_{#pi_{1}#pi_{4}} OS",massbins,0,5.);
   //
   histosTH1F["hm4rec2OS_pi4k2veeno11"] = fs->make<TH1F>("hm4rec2OS_pi4k2veeno11","M_{#pi_{4}K_{2}} OS",massbins,0,5.);
   histosTH1F["hm4rec2OS_pi1pi3k2veeno11"] = fs->make<TH1F>("hm4rec2OS_pi1pi3k2veeno11","M_{#pi_{1}#pi_{3}} OS",massbins,0,5.);
   //
   //...D...no PID  k1
   histosTH1F["hm4rec2OS_k1pi2veeno11"] = fs->make<TH1F>("hm4rec2OS_k1pi2veeno11","M_{K_{1}#pi_{2}} OS",massbins,0,5.);
   histosTH1F["hm4rec2OS_pi3pi4k1veeno11"] = fs->make<TH1F>("hm4rec2OS_pi3pi4k1veeno11","M_{#pi_{3}#pi_{4}} OS",massbins,0,5.);
   //
   histosTH1F["hm4rec2OS_k1pi3veeno11"] = fs->make<TH1F>("hm4rec2OS_k1pi3veeno11","M_{K_{1}#pi_{3}} OS",massbins,0,5.);
   histosTH1F["hm4rec2OS_pi2pi4k1veeno11"] = fs->make<TH1F>("hm4rec2OS_pi2pi4k1veeno11","M_{#pi_{2}#pi_{4}} OS",massbins,0,5.);
   //
   histosTH1F["hm4rec2OS_k1pi4veeno11"] = fs->make<TH1F>("hm4rec2OS_k1pi4veeno11","M_{K_{1}#pi_{4}} OS",massbins,0,5.);
   histosTH1F["hm4rec2OS_pi2pi3k1veeno11"] = fs->make<TH1F>("hm4rec2OS_pi2pi3k1veeno11","M_{#pi_{2}#pi_{3}} OS",massbins,0,5.);
   
  //...vee02
  //**histosTH1F["hm4rec2OS_pi1pi2vee02"] = fs->make<TH1F>("hm4rec2OS_pi1pi2vee02","M_{#pi_{1}#pi_{2}} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi3pi4vee02"] = fs->make<TH1F>("hm4rec2OS_pi3pi4vee02","M_{#pi_{3}#pi_{4}} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi1pi3vee02"] = fs->make<TH1F>("hm4rec2OS_pi1pi3vee02","M_{#pi_{1}#pi_{3}} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi2pi4vee02"] = fs->make<TH1F>("hm4rec2OS_pi2pi4vee02","M_{#pi_{2}#pi_{4}} OS",massbins,0,5.);
   //...vee01
  //**histosTH1F["hm4rec2OS_pi1pi2vee01"] = fs->make<TH1F>("hm4rec2OS_pi1pi2vee01","M_{#pi_{1}#pi_{2}} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi3pi4vee01"] = fs->make<TH1F>("hm4rec2OS_pi3pi4vee01","M_{#pi_{3}#pi_{4}} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi1pi3vee01"] = fs->make<TH1F>("hm4rec2OS_pi1pi3vee01","M_{#pi_{1}#pi_{3}} OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_pi2pi4vee01"] = fs->make<TH1F>("hm4rec2OS_pi2pi4vee01","M_{#pi_{2}#pi_{4}} OS",massbins,0,5.);
  
  //...2dim vee11

  //...2dim vee02
  //***histosTH2F["hm4dim2OS_pi1pi2_pi3pi4vee02"] = fs->make<TH2F>("hm4dim2OS_pi1pi2_pi3pi4vee02","M_{#pi_{1}#pi_{2}} vs M_{#pi_{3}#pi_{4}} OS",massbins,0,5.,massbins,0,5.);
  //***histosTH2F["hm4dim2OS_pi1pi3_pi2pi4vee02"] = fs->make<TH2F>("hm4dim2OS_pi1pi3_pi2pi4vee02","M_{#pi_{1}#pi_{3}} vs M_{#pi_{2}#pi_{4}} OS",massbins,0,5.,massbins,0,5.);
  //...2dim vee01
  //***histosTH2F["hm4dim2OS_pi1pi2_pi3pi4vee01"] = fs->make<TH2F>("hm4dim2OS_pi1pi2_pi3pi4vee01","M_{#pi_{1}#pi_{2}} vs M_{#pi_{3}#pi_{4}} OS",massbins,0,5.,massbins,0,5.);
  //***histosTH2F["hm4dim2OS_pi1pi3_pi2pi4vee01"] = fs->make<TH2F>("hm4dim2OS_pi1pi3_pi2pi4vee01","M_{#pi_{1}#pi_{3}} vs M_{#pi_{2}#pi_{4}} OS",massbins,0,5.,massbins,0,5.);
  
  //...veeno02
  histosTH1F["hm4rec2OS_pi1pi2veeno02"] = fs->make<TH1F>("hm4rec2OS_pi1pi2veeno02","M_{#pi_{1}#pi_{2}} OS",massbins,0,5.);
  histosTH1F["hm4rec2OS_pi3pi4veeno02"] = fs->make<TH1F>("hm4rec2OS_pi3pi4veeno02","M_{#pi_{3}#pi_{4}} OS",massbins,0,5.);
  histosTH1F["hm4rec2OS_pi1pi3veeno02"] = fs->make<TH1F>("hm4rec2OS_pi1pi3veeno02","M_{#pi_{1}#pi_{3}} OS",massbins,0,5.);
  histosTH1F["hm4rec2OS_pi2pi4veeno02"] = fs->make<TH1F>("hm4rec2OS_pi2pi4veeno02","M_{#pi_{2}#pi_{4}} OS",massbins,0,5.);
  //...veeno01
  histosTH1F["hm4rec2OS_pi1pi2veeno01"] = fs->make<TH1F>("hm4rec2OS_pi1pi2veeno01","M_{#pi_{1}#pi_{2}} OS",massbins,0,5.);
  histosTH1F["hm4rec2OS_pi3pi4veeno01"] = fs->make<TH1F>("hm4rec2OS_pi3pi4veeno01","M_{#pi_{3}#pi_{4}} OS",massbins,0,5.);
  histosTH1F["hm4rec2OS_pi1pi3veeno01"] = fs->make<TH1F>("hm4rec2OS_pi1pi3veeno01","M_{#pi_{1}#pi_{3}} OS",massbins,0,5.);
  histosTH1F["hm4rec2OS_pi2pi4veeno01"] = fs->make<TH1F>("hm4rec2OS_pi2pi4veeno01","M_{#pi_{2}#pi_{4}} OS",massbins,0,5.);
  
  //...veenopid02
  histosTH1F["hm4rec2OS_pi1pi2veenopid02"] = fs->make<TH1F>("hm4rec2OS_pi1pi2veenopid02","M_{#pi_{1}#pi_{2}} OS",massbins,0,5.);
  histosTH1F["hm4rec2OS_pi3pi4veenopid02"] = fs->make<TH1F>("hm4rec2OS_pi3pi4veenopid02","M_{#pi_{3}#pi_{4}} OS",massbins,0,5.);
  histosTH1F["hm4rec2OS_pi1pi3veenopid02"] = fs->make<TH1F>("hm4rec2OS_pi1pi3veenopid02","M_{#pi_{1}#pi_{3}} OS",massbins,0,5.);
  histosTH1F["hm4rec2OS_pi2pi4veenopid02"] = fs->make<TH1F>("hm4rec2OS_pi2pi4veenopid02","M_{#pi_{2}#pi_{4}} OS",massbins,0,5.);
  //...veenopid01
  histosTH1F["hm4rec2OS_pi1pi2veenopid01"] = fs->make<TH1F>("hm4rec2OS_pi1pi2veenopid01","M_{#pi_{1}#pi_{2}} OS",massbins,0,5.);
  histosTH1F["hm4rec2OS_pi3pi4veenopid01"] = fs->make<TH1F>("hm4rec2OS_pi3pi4veenopid01","M_{#pi_{3}#pi_{4}} OS",massbins,0,5.);
  histosTH1F["hm4rec2OS_pi1pi3veenopid01"] = fs->make<TH1F>("hm4rec2OS_pi1pi3veenopid01","M_{#pi_{1}#pi_{3}} OS",massbins,0,5.);
  histosTH1F["hm4rec2OS_pi2pi4veenopid01"] = fs->make<TH1F>("hm4rec2OS_pi2pi4veenopid01","M_{#pi_{2}#pi_{4}} OS",massbins,0,5.);
  
  //...2dim veeno11

  //...2dim veeno02
  histosTH2F["hm4dim2OS_pi1pi2_pi3pi4veeno02"] = fs->make<TH2F>("hm4dim2OS_pi1pi2_pi3pi4veeno02","M_{#pi_{1}#pi_{2}} vs M_{#pi_{3}#pi_{4}} OS",massbins,0,5.,massbins,0,5.);
  histosTH2F["hm4dim2OS_pi1pi3_pi2pi4veeno02"] = fs->make<TH2F>("hm4dim2OS_pi1pi3_pi2pi4veeno02","M_{#pi_{1}#pi_{3}} vs M_{#pi_{2}#pi_{4}} OS",massbins,0,5.,massbins,0,5.);
  //...2dim veeno01
  histosTH2F["hm4dim2OS_pi1pi2_pi3pi4veeno01"] = fs->make<TH2F>("hm4dim2OS_pi1pi2_pi3pi4veeno01","M_{#pi_{1}#pi_{2}} vs M_{#pi_{3}#pi_{4}} OS",massbins,0,5.,massbins,0,5.);
  histosTH2F["hm4dim2OS_pi1pi3_pi2pi4veeno01"] = fs->make<TH2F>("hm4dim2OS_pi1pi3_pi2pi4veeno01","M_{#pi_{1}#pi_{3}} vs M_{#pi_{2}#pi_{4}} OS",massbins,0,5.,massbins,0,5.);
  //...2dim veeno11

  //...2dim veenopid02
  histosTH2F["hm4dim2OS_pi1pi2_pi3pi4veenopid02"] = fs->make<TH2F>("hm4dim2OS_pi1pi2_pi3pi4veenopid02","M_{#pi_{1}#pi_{2}} vs M_{#pi_{3}#pi_{4}} OS",massbins,0,5.,massbins,0,5.);
  histosTH2F["hm4dim2OS_pi1pi3_pi2pi4veenopid02"] = fs->make<TH2F>("hm4dim2OS_pi1pi3_pi2pi4veenopid02","M_{#pi_{1}#pi_{3}} vs M_{#pi_{2}#pi_{4}} OS",massbins,0,5.,massbins,0,5.);
  //...2dim veeno01
  histosTH2F["hm4dim2OS_pi1pi2_pi3pi4veenopid01"] = fs->make<TH2F>("hm4dim2OS_pi1pi2_pi3pi4veenopid01","M_{#pi_{1}#pi_{2}} vs M_{#pi_{3}#pi_{4}} OS",massbins,0,5.,massbins,0,5.);
  histosTH2F["hm4dim2OS_pi1pi3_pi2pi4veenopid01"] = fs->make<TH2F>("hm4dim2OS_pi1pi3_pi2pi4veenopid01","M_{#pi_{1}#pi_{3}} vs M_{#pi_{2}#pi_{4}} OS",massbins,0,5.,massbins,0,5.);
    
  //...........
  //
  //...Kaons
  //**histosTH1F["hm4rec2OS_k1k2"] = fs->make<TH1F>("hm4rec2OS_k1k2","M_{k_{1}k_{2}} OS",2*massbins,0,5.);
  //**histosTH1F["hm4rec2OS_k3k4"] = fs->make<TH1F>("hm4rec2OS_k3k4","M_{k_{3}k_{4}} OS",2*massbins,0,5.);
  //**histosTH1F["hm4rec2OS_k1k3"] = fs->make<TH1F>("hm4rec2OS_k1k3","M_{k_{1}k_{3}} OS",2*massbins,0,5.);
  //**histosTH1F["hm4rec2OS_k2k4"] = fs->make<TH1F>("hm4rec2OS_k2k4","M_{k_{2}k_{4}} OS",2*massbins,0,5.);
  //...v2
  //**histosTH1F["hm4rec2OS_k1k2v2"] = fs->make<TH1F>("hm4rec2OS_k1k2v2","M_{k_{1}k_{2}} OS",2*massbins,0,5.);
  //**histosTH1F["hm4rec2OS_k3k4v2"] = fs->make<TH1F>("hm4rec2OS_k3k4v2","M_{k_{3}k_{4}} OS",2*massbins,0,5.);
  //**histosTH1F["hm4rec2OS_k1k3v2"] = fs->make<TH1F>("hm4rec2OS_k1k3v2","M_{k_{1}k_{3}} OS",2*massbins,0,5.);
  //**histosTH1F["hm4rec2OS_k2k4v2"] = fs->make<TH1F>("hm4rec2OS_k2k4v2","M_{k_{2}k_{4}} OS",2*massbins,0,5.);
  //...2dim
  //***histosTH2F["hm4dim2OS_k1k2_k3k4"] = fs->make<TH2F>("hm4dim2OS_k1k2_k3k4","M_{k_{1}k_{2}} vs M_{k_{3}k_{4}} OS",2*massbins,0,5.,2*massbins,0,5.);
  //***histosTH2F["hm4dim2OS_k1k3_k2k4"] = fs->make<TH2F>("hm4dim2OS_k1k3_k2k4","M_{k_{1}k_{3}} vs M_{k_{2}k_{4}} OS",2*massbins,0,5.,2*massbins,0,5.);
  //...v2
  //***histosTH2F["hm4dim2OS_k1k2_k3k4v2"] = fs->make<TH2F>("hm4dim2OS_k1k2_k3k4v2","M_{k_{1}k_{2}} vs M_{k_{3}k_{4}} OS",2*massbins,0,5.,2*massbins,0,5.);
  //***histosTH2F["hm4dim2OS_k1k3_k2k4v2"] = fs->make<TH2F>("hm4dim2OS_k1k3_k2k4v2","M_{k_{1}k_{3}} vs M_{k_{2}k_{4}} OS",2*massbins,0,5.,2*massbins,0,5.);
  //...Kaons
  //
  //**histosTH1F["hm4rec2SS"] = fs->make<TH1F>("hm4rec2SS","M_{4#pi} SS",massbins,0,5.);
  //...2OSdiag
  //**histosTH1F["hm4rec2OS_diag"] = fs->make<TH1F>("hm4rec2OS_diag","M_{4#pi} TB/BT OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_diag2"] = fs->make<TH1F>("hm4rec2OS_diag2","M_{4#pi} TB/BT OS",1.60*massbins,0.0,8.0);
  //**histosTH1F["hm4rec2OS_diag3"] = fs->make<TH1F>("hm4rec2OS_diag3","M_{4#pi} TB/BT OS",0.50*massbins,0.0,2.5);
  //**histosTH1F["hm4rec2OS_diag4"] = fs->make<TH1F>("hm4rec2OS_diag4","M_{4#pi} TB/BT OS",0.24*massbins,2.5,4.0);
  //**histosTH1F["hm4rec2OS_diag5"] = fs->make<TH1F>("hm4rec2OS_diag5","M_{4#pi} TB/BT OS",0.32*massbins,4.0,8.0);

 //0-2.5 (125bins), 2.5-4(60bins), 4-8(80bins)
 double xmin1 = 0.;
 double xmax1 = 2.5;
 const int nbins1 = 125;

 double xmin2 = 2.5;
 double xmax2 = 4.;
 const int nbins2 = 60;

 double xmin3 = 4.;
 double xmax3 = 8.;
 const int nbins3 = 80;

 //**double bwidth1 = (xmax1 - xmin1)/nbins1;
 //**double bwidth2 = (xmax2 - xmin2)/nbins2;
 //**double bwidth3 = (xmax3 - xmin3)/nbins3;

 const int nbinstot = nbins1 + nbins2 + nbins3;
 //...Luiz
 //**double edges[nbinstot+1] ;
 
 //nbinstot++;

 int nbins=0;

 //**for( int i=0; i<nbins1; i++){ edges[nbins] = xmin1 + bwidth1 * i; nbins++;}
 //**for( int i=0; i<nbins2; i++){ edges[nbins] = xmin2 + bwidth2 * i; nbins++;}
 //...Luiz
 //**for( int i=0; i<=nbins3; i++){ edges[nbins] = xmin3 + bwidth3 * i; nbins++;}

 //**histosTH1F["hm4rec2OS_ttbb2varbin"] = fs->make<TH1F>("hm4rec2OS_ttbb2varbin","TTBB variable bins",nbinstot,edges);
 //**histosTH1F["hm4rec2OS_diag2varbin"] = fs->make<TH1F>("hm4rec2OS_diag2varbin","DIAG variable bins",nbinstot,edges);

  //...Pions
  //**histosTH1F["hm4rec2OS_diag_pi1pi2"] = fs->make<TH1F>("hm4rec2OS_diag_pi1pi2","M_{#pi_{1}#pi_{2}} OS",2.0*massbins,0,10.);
  //**histosTH1F["hm4rec2OS_diag_pi3pi4"] = fs->make<TH1F>("hm4rec2OS_diag_pi3pi4","M_{#pi_{3}#pi_{4}} OS",2.0*massbins,0,10.);
  //**histosTH1F["hm4rec2OS_diag_pi1pi3"] = fs->make<TH1F>("hm4rec2OS_diag_pi1pi3","M_{#pi_{1}#pi_{3}} OS",2.0*massbins,0,10.);
  //**histosTH1F["hm4rec2OS_diag_pi2pi4"] = fs->make<TH1F>("hm4rec2OS_diag_pi2pi4","M_{#pi_{2}#pi_{4}} OS",2.0*massbins,0,10.);
  //...Kaons
  //**histosTH1F["hm4rec2OS_diag_k1k2"] = fs->make<TH1F>("hm4rec2OS_diag_k1k2","M_{k_{1}k_{2}} OS",2.0*massbins,0,10.);
  //**histosTH1F["hm4rec2OS_diag_k3k4"] = fs->make<TH1F>("hm4rec2OS_diag_k3k4","M_{k_{3}k_{4}} OS",2.0*massbins,0,10.);
  //**histosTH1F["hm4rec2OS_diag_k1k3"] = fs->make<TH1F>("hm4rec2OS_diag_k1k3","M_{k_{1}k_{3}} OS",2.0*massbins,0,10.);
  //**histosTH1F["hm4rec2OS_diag_k2k4"] = fs->make<TH1F>("hm4rec2OS_diag_k2k4","M_{k_{2}k_{4}} OS",2.0*massbins,0,10.);
  
  //...2SSdiag
  //**histosTH1F["hm4rec2SS_diag"] = fs->make<TH1F>("hm4rec2SS_diag","M_{4#pi} TB/BT SS",massbins,0,5.);
  //
  //...2OSttbb
  //**histosTH1F["hm4rec2OS_ttbb"] = fs->make<TH1F>("hm4rec2OS_ttbb","M_{4#pi} TT/BB OS",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_ttbb2"] = fs->make<TH1F>("hm4rec2OS_ttbb2","M_{4#pi} TT/BB OS",1.60*massbins,0.0,8.0);
  //**histosTH1F["hm4rec2OS_ttbb3"] = fs->make<TH1F>("hm4rec2OS_ttbb3","M_{4#pi} TT/BB OS",0.50*massbins,0.0,2.5);
  //**histosTH1F["hm4rec2OS_ttbb4"] = fs->make<TH1F>("hm4rec2OS_ttbb4","M_{4#pi} TT/BB OS",0.24*massbins,2.5,4.0);
  //**histosTH1F["hm4rec2OS_ttbb5"] = fs->make<TH1F>("hm4rec2OS_ttbb5","M_{4#pi} TT/BB OS",0.32*massbins,4.0,8.0);
  
  //...Pions
  //**histosTH1F["hm4rec2OS_ttbb_pi1pi2"] = fs->make<TH1F>("hm4rec2OS_ttbb_pi1pi2","M_{#pi_{1}#pi_{2}} OS",2.0*massbins,0,10.);
  //**histosTH1F["hm4rec2OS_ttbb_pi3pi4"] = fs->make<TH1F>("hm4rec2OS_ttbb_pi3pi4","M_{#pi_{3}#pi_{4}} OS",2.0*massbins,0,10.);
  //**histosTH1F["hm4rec2OS_ttbb_pi1pi3"] = fs->make<TH1F>("hm4rec2OS_ttbb_pi1pi3","M_{#pi_{1}#pi_{3}} OS",2.0*massbins,0,10.);
  //**histosTH1F["hm4rec2OS_ttbb_pi2pi4"] = fs->make<TH1F>("hm4rec2OS_ttbb_pi2pi4","M_{#pi_{2}#pi_{4}} OS",2.0*massbins,0,10.);
  //...Kaons
  //**histosTH1F["hm4rec2OS_ttbb_k1k2"] = fs->make<TH1F>("hm4rec2OS_ttbb_k1k2","M_{k_{1}k_{2}} OS",2.0*massbins,0,10.);
  //**histosTH1F["hm4rec2OS_ttbb_k3k4"] = fs->make<TH1F>("hm4rec2OS_ttbb_k3k4","M_{k_{3}k_{4}} OS",2.0*massbins,0,10.);
  //**histosTH1F["hm4rec2OS_ttbb_k1k3"] = fs->make<TH1F>("hm4rec2OS_ttbb_k1k3","M_{k_{1}k_{3}} OS",2.0*massbins,0,10.);
  //**histosTH1F["hm4rec2OS_ttbb_k2k4"] = fs->make<TH1F>("hm4rec2OS_ttbb_k2k4","M_{k_{2}k_{4}} OS",2.0*massbins,0,10.);
  //...2SSttbb
  //**histosTH1F["hm4rec2SS_ttbb"] = fs->make<TH1F>("hm4rec2SS_ttbb","M_{4#pi} TT/BB SS",massbins,0,5.);
  //
  
  //...Luiz
  //**histosTH1F["hm4rec2OS_diag_trkP"] = fs->make<TH1F>("hm4rec2OS_diag_trkP","M_{4#pi} TB/BT OS, py_{#pi_{1}}py_{#pi_{2}}>0",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_diag_trkM"] = fs->make<TH1F>("hm4rec2OS_diag_trkM","M_{4#pi} TB/BT OS, py_{#pi_{1}}py_{#pi_{2}}<0",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_ttbb_trkP"] = fs->make<TH1F>("hm4rec2OS_ttbb_trkP","M_{4#pi} TT/BB OS, py_{#pi_{1}}py_{#pi_{2}}>0",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_ttbb_trkM"] = fs->make<TH1F>("hm4rec2OS_ttbb_trkM","M_{4#pi} TT/BB OS, py_{#pi_{1}}py_{#pi_{2}}<0",massbins,0,5.);
  
  //**histosTH1F["hm4rec2OS_diag_pypxP"] = fs->make<TH1F>("hm4rec2OS_diag_pypxP","M_{4#pi} TB/BT OS, |py/px|_{4#pi} > 1",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_diag_pypxM"] = fs->make<TH1F>("hm4rec2OS_diag_pypxM","M_{4#pi} TB/BT OS, |py/px|_{4#pi} < 1",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_ttbb_pypxP"] = fs->make<TH1F>("hm4rec2OS_ttbb_pypxP","M_{4#pi} TT/BB OS, |py/px|_{4#pi} > 1",massbins,0,5.);
  //**histosTH1F["hm4rec2OS_ttbb_pypxM"] = fs->make<TH1F>("hm4rec2OS_ttbb_pypxM","M_{4#pi} TT/BB OS, |py/px|_{4#pi} < 1",massbins,0,5.);

  //**histosTH1F["hm4rec3OS"] = fs->make<TH1F>("hm4rec3OS","M_{4#pi} OS",massbins,0,5.);
  //**histosTH1F["hm4rec3SS"] = fs->make<TH1F>("hm4rec3SS","M_{4#pi} SS",massbins,0,5.);
  //**histosTH1F["hm4rec3OS_diag"] = fs->make<TH1F>("hm4rec3OS_diag","M_{4#pi} TB/BT OS",massbins,0,5.);
  //**histosTH1F["hm4rec3SS_diag"] = fs->make<TH1F>("hm4rec3SS_diag","M_{4#pi} TB/BT SS",massbins,0,5.);
  //**histosTH1F["hm4rec3OS_ttbb"] = fs->make<TH1F>("hm4rec3OS_ttbb","M_{4#pi} TT/BB OS",massbins,0,5.);
  //**histosTH1F["hm4rec3SS_ttbb"] = fs->make<TH1F>("hm4rec3SS_ttbb","M_{4#pi} TT/BB SS",massbins,0,5.);

  //...Luiz
  //**histosTH1F["hm4rec3OS_diag_trkP"] = fs->make<TH1F>("hm4rec3OS_diag_trkP","M_{4#pi} TB/BT OS, py_{#pi_{1}}py_{#pi_{2}}>0",massbins,0,5.);
  //**histosTH1F["hm4rec3OS_diag_trkM"] = fs->make<TH1F>("hm4rec3OS_diag_trkM","M_{4#pi} TB/BT OS, py_{#pi_{1}}py_{#pi_{2}}<0",massbins,0,5.);
  //**histosTH1F["hm4rec3OS_ttbb_trkP"] = fs->make<TH1F>("hm4rec3OS_ttbb_trkP","M_{4#pi} TT/BB OS, py_{#pi_{1}}py_{#pi_{2}}>0",massbins,0,5.);
  //**histosTH1F["hm4rec3OS_ttbb_trkM"] = fs->make<TH1F>("hm4rec3OS_ttbb_trkM","M_{4#pi} TT/BB OS, py_{#pi_{1}}py_{#pi_{2}}<0",massbins,0,5.);
  
  //**histosTH1F["hm4rec3OS_diag_pypxP"] = fs->make<TH1F>("hm4rec3OS_diag_pypxP","M_{4#pi} TB/BT OS, |py/px|_{4#pi} > 1",massbins,0,5.);
  //**histosTH1F["hm4rec3OS_diag_pypxM"] = fs->make<TH1F>("hm4rec3OS_diag_pypxM","M_{4#pi} TB/BT OS, |py/px|_{4#pi} < 1",massbins,0,5.);
  //**histosTH1F["hm4rec3OS_ttbb_pypxP"] = fs->make<TH1F>("hm4rec3OS_ttbb_pypxP","M_{4#pi} TT/BB OS, |py/px|_{4#pi} > 1",massbins,0,5.);
  //**histosTH1F["hm4rec3OS_ttbb_pypxM"] = fs->make<TH1F>("hm4rec3OS_ttbb_pypxM","M_{4#pi} TT/BB OS, |py/px|_{4#pi} < 1",massbins,0,5.);
  
  //**histosTH1F["hm4rec4OS"] = fs->make<TH1F>("hm4rec4OS","M_{4#pi} OS",massbins,0,5.);
  //**histosTH1F["hm4rec4SS"] = fs->make<TH1F>("hm4rec4SS","M_{4#pi} SS",massbins,0,5.);
  //**histosTH1F["hm4rec4OS_diag"] = fs->make<TH1F>("hm4rec4OS_diag","M_{4#pi} TB/BT OS",massbins,0,5.);
  //**histosTH1F["hm4rec4SS_diag"] = fs->make<TH1F>("hm4rec4SS_diag","M_{4#pi} TB/BT SS",massbins,0,5.);
  //**histosTH1F["hm4rec4OS_ttbb"] = fs->make<TH1F>("hm4rec4OS_ttbb","M_{4#pi} TT/BB OS",massbins,0,5.);
  //**histosTH1F["hm4rec4SS_ttbb"] = fs->make<TH1F>("hm4rec4SS_ttbb","M_{4#pi} TT/BB SS",massbins,0,5.);

  //...Luiz
  //**histosTH1F["hm4rec4OS_diag_trkP"] = fs->make<TH1F>("hm4rec4OS_diag_trkP","M_{4#pi}TB/BT OS, py_{#pi_{1}}py_{#pi_{2}}>0",massbins,0,5.);
  //**histosTH1F["hm4rec4OS_diag_trkM"] = fs->make<TH1F>("hm4rec4OS_diag_trkM","M_{4#pi}TB/BT OS, py_{#pi_{1}}py_{#pi_{2}}<0",massbins,0,5.);  
  //**histosTH1F["hm4rec4OS_ttbb_trkP"] = fs->make<TH1F>("hm4rec4OS_ttbb_trkP","M_{4#pi}TT/BB OS, py_{#pi_{1}}py_{#pi_{2}}>0",massbins,0,5.);
  //**histosTH1F["hm4rec4OS_ttbb_trkM"] = fs->make<TH1F>("hm4rec4OS_ttbb_trkM","M_{4#pi}TT/BB OS, py_{#pi_{1}}py_{#pi_{2}}<0",massbins,0,5.);

  //**histosTH1F["hm4rec4OS_diag_pypxP"] = fs->make<TH1F>("hm4rec4OS_diag_pypxP","M_{4#pi} TB/BT OS, |py/px|_{4#pi} > 1",massbins,0,5.);
  //**histosTH1F["hm4rec4OS_diag_pypxM"] = fs->make<TH1F>("hm4rec4OS_diag_pypxM","M_{4#pi} TB/BT OS, |py/px|_{4#pi} < 1",massbins,0,5.);
  //**histosTH1F["hm4rec4OS_ttbb_pypxP"] = fs->make<TH1F>("hm4rec4OS_ttbb_pypxP","M_{4#pi} TT/BB OS, |py/px|_{4#pi} > 1",massbins,0,5.);
  //**histosTH1F["hm4rec4OS_ttbb_pypxM"] = fs->make<TH1F>("hm4rec4OS_ttbb_pypxM","M_{4#pi} TT/BB OS, |py/px|_{4#pi} < 1",massbins,0,5.);

  //**histosTH1F["hm4rec5OS"] = fs->make<TH1F>("hm4rec5OS","M_{4#pi} OS",massbins,0,5.);
  //**histosTH1F["hm4rec5SS"] = fs->make<TH1F>("hm4rec5SS","M_{4#pi} SS",massbins,0,5.);
  //**histosTH1F["hm4rec5OS_diag"] = fs->make<TH1F>("hm4rec5OS_diag","M_{4#pi} TB/BT OS",massbins,0,5.);
  //**histosTH1F["hm4rec5SS_diag"] = fs->make<TH1F>("hm4rec5SS_diag","M_{4#pi} TB/BT SS",massbins,0,5.);
  //**histosTH1F["hm4rec5OS_ttbb"] = fs->make<TH1F>("hm4rec5OS_ttbb","M_{4#pi} TT/BB OS",massbins,0,5.);
  //**histosTH1F["hm4rec5SS_ttbb"] = fs->make<TH1F>("hm4rec5SS_ttbb","M_{4#pi} TT/BB SS",massbins,0,5.);

  //...Luiz
  //**histosTH1F["hm4rec5OS_diag_trkP"] = fs->make<TH1F>("hm4rec5OS_diag_trkP","M_{4#pi} TB/BT OS, py_{#pi_{1}}py_{#pi_{2}}>0",massbins,0,5.);
  //**histosTH1F["hm4rec5OS_diag_trkM"] = fs->make<TH1F>("hm4rec5OS_diag_trkM","M_{4#pi} TB/BT OS, py_{#pi_{1}}py_{#pi_{2}}<0",massbins,0,5.);  
  //**histosTH1F["hm4rec5OS_ttbb_trkP"] = fs->make<TH1F>("hm4rec5OS_ttbb_trkP","M_{4#pi} TB/BT OS, py_{#pi_{1}}py_{#pi_{2}}>0",massbins,0,5.);
  //**histosTH1F["hm4rec5OS_ttbb_trkM"] = fs->make<TH1F>("hm4rec5OS_ttbb_trkM","M_{4#pi} TB/BT OS, py_{#pi_{1}}py_{#pi_{2}}<0",massbins,0,5.);

  //**histosTH1F["hm4rec5OS_diag_pypxP"] = fs->make<TH1F>("hm4rec5OS_diag_pypxP","M_{4#pi} TB/BT OS, |py/px|_{4#pi} > 1",massbins,0,5.);
  //**histosTH1F["hm4rec5OS_diag_pypxM"] = fs->make<TH1F>("hm4rec5OS_diag_pypxM","M_{4#pi} TB/BT OS, |py/px|_{4#pi} < 1",massbins,0,5.);
  //**histosTH1F["hm4rec5OS_ttbb_pypxP"] = fs->make<TH1F>("hm4rec5OS_ttbb_pypxP","M_{4#pi} TT/BB OS, |py/px|_{4#pi} > 1",massbins,0,5.);
  //**histosTH1F["hm4rec5OS_ttbb_pypxM"] = fs->make<TH1F>("hm4rec5OS_ttbb_pypxM","M_{4#pi} TT/BB OS, |py/px|_{4#pi} < 1",massbins,0,5.);

  //**histosTH1F["hm4rec6OS"] = fs->make<TH1F>("hm4rec6OS","M_{4#pi} OS",massbins,0,5.);
  //**histosTH1F["hm4rec6SS"] = fs->make<TH1F>("hm4rec6SS","M_{4#pi} SS",massbins,0,5.);
  //**histosTH1F["hm4rec6OS_diag"] = fs->make<TH1F>("hm4rec6OS_diag","M_{4#pi} TB/BT OS",massbins,0,5.);
  //**histosTH1F["hm4rec6SS_diag"] = fs->make<TH1F>("hm4rec6SS_diag","M_{4#pi} TB/BT SS",massbins,0,5.);
  //**histosTH1F["hm4rec6OS_ttbb"] = fs->make<TH1F>("hm4rec6OS_ttbb","M_{4#pi} TT/BB OS",massbins,0,5.);
  //**histosTH1F["hm4rec6SS_ttbb"] = fs->make<TH1F>("hm4rec6SS_ttbb","M_{4#pi} TT/BB SS",massbins,0,5.);

  //...Luiz
  //**histosTH1F["hm4rec6OS_diag_trkP"] = fs->make<TH1F>("hm4rec6OS_diag_trkP","M_{4#pi} TB/BT OS, py_{#pi_{1}}py_{#pi_{2}}>0",massbins,0,5.);
  //**histosTH1F["hm4rec6OS_diag_trkM"] = fs->make<TH1F>("hm4rec6OS_diag_trkM","M_{4#pi} TB/BT OS, py_{#pi_{1}}py_{#pi_{2}}<0",massbins,0,5.);  
  //**histosTH1F["hm4rec6OS_ttbb_trkP"] = fs->make<TH1F>("hm4rec6OS_ttbb_trkP","M_{4#pi} TT/BB OS, py_{#pi_{1}}py_{#pi_{2}}>0",massbins,0,5.);
  //**histosTH1F["hm4rec6OS_ttbb_trkM"] = fs->make<TH1F>("hm4rec6OS_ttbb_trkM","M_{4#pi} TT/BB OS, py_{#pi_{1}}py_{#pi_{2}}<0",massbins,0,5.);

  //**histosTH1F["hm4recHFvetoOS"] = fs->make<TH1F>("hm4recHFvetoOS","M_{4#pi} HFv OS",massbins,0,5.);
  //**histosTH1F["hm4recHFvetoSS"] = fs->make<TH1F>("hm4recHFvetoSS","M_{4#pi} HFv SS",massbins,0,5.);
  
  //**histosTH1F["hm4rec45OS"] = fs->make<TH1F>("hm4rec45OS","M_{4#pi} OS",massbins,0,5.);
  //**histosTH1F["hm4rec45SS"] = fs->make<TH1F>("hm4rec45SS","M_{4#pi} SS",massbins,0,5.);
  //**histosTH1F["hm4rec4515OS"] = fs->make<TH1F>("hm4rec4515OS","M_{4#pi} OS",massbins,0,5.);
  //**histosTH1F["hm4rec4515SS"] = fs->make<TH1F>("hm4rec4515SS","M_{4#pi} SS",massbins,0,5.);

  //**histosTH1F["hm4rec9919"] = fs->make<TH1F>("hm4rec9919","M_{4#pi} 9919",massbins,0,5.);
  //**histosTH1F["hm4rec9922"] = fs->make<TH1F>("hm4rec9922","M_{4#pi} 9919,9922",massbins,0,5.);
  //**histosTH1F["hm4rec9971"] = fs->make<TH1F>("hm4rec9971","M_{4#pi} 9971",massbins,0,5.);
  //**histosTH1F["hm4rec9978"] = fs->make<TH1F>("hm4rec9978","M_{4#pi} 9978",massbins,0,5.);

  //**histosTH1F["hnclusters"] = fs->make<TH1F>("hnclusters","nPixelClusters",500,0,500.);
  //**histosTH1F["hnclusters2"] = fs->make<TH1F>("hnclusters2","nStripClusters",500,0,500.);
  //**histosTH1F["hnclustersOSdiag"] = fs->make<TH1F>("hnclustersOSdiag","nPixelClusters",500,0,500.);
  //**histosTH1F["hnclusters2OSdiag"] = fs->make<TH1F>("hnclusters2OSdiag","nStripClusters",500,0,500.);
  
  //histosTH1F["halgo"] = fs->make<TH1F>("halgo","Algo",15,0,15.);
  //histosTH1F["hnhits"] = fs->make<TH1F>("hnhits","nhits pix+strip",40,0,40.);
  //histosTH1F["hchi2"] = fs->make<TH1F>("hchi2","normalized #chi^{2}",1050,-50,1000.);   

  //...Luiz
  //histosTH1F["hdz"] = fs->make<TH1F>("hdz","dz",2000,-100,100.);
  //histosTH1F["hd0"] = fs->make<TH1F>("hd0","d0",2000,-100,100.);
  
  //**histosTH1F["halgov"] = fs->make<TH1F>("halgov","Algo",15,0,15.);
  //**histosTH1F["hnhitsv"] = fs->make<TH1F>("hnhitsv","nhits pixel",40,0,40.);
  //**histosTH1F["hchi2v"] = fs->make<TH1F>("hchi2v","normalized #chi^{2} vtx-fitted",550,-50,500.);   

  //...Luiz
  //**histosTH1F["hdzv"] = fs->make<TH1F>("hdzv","dz vtx-fitted",1000,-100,100.);
  //**histosTH1F["hd0v"] = fs->make<TH1F>("hd0v","d0 vtx-fitted",2000,-20,20.);

  //**histosTH1F["hchi2fin"] = fs->make<TH1F>("hchi2fin","normalized #chi^{2} vtx-fitted",550,-50,500.);   

  //...Luiz
  //**histosTH1F["hdzfin"] = fs->make<TH1F>("hdzfin","dz vtx-fitted",1000,-100,100.);
  //**histosTH1F["hd0fin"] = fs->make<TH1F>("hd0fin","d0 vtx-fitted",2000,-20,20.);
  
  //...Luiz
  //**histosTH1F["hdeltaR"] = fs->make<TH1F>("hdeltaR","#DeltaR trk-trk",200,0,10.);
  //**histosTH1F["hdeltaR2"] = fs->make<TH1F>("hdeltaR2","#DeltaR trk-trk",200,0,10.);

  //-----------------
  histosTH2F["h2dimdpyAll"] = fs->make<TH2F>("h2dimdpyAll", "p_{y}^{TOTEM} vs p_{y}^{CMS}",200,-2.,2.,200,-2., 2.);
  histosTH2F["h2dimdpy"] = fs->make<TH2F>("h2dimdpy","p_{y}^{TOTEM} vs p_{y}^{CMS}",200,-2.,2.,200,-2.,2.);
  //***histosTH2F["h2dimdpy_diag"] = fs->make<TH2F>("h2dimdpy_diag","p_{y}^{TOTEM} vs p_{y}^{CMS} diag",100,-2.,2.,100,-2.,2.);
  //***histosTH2F["h2dimdpy_ttbb"] = fs->make<TH2F>("h2dimdpy_ttbb","p_{y}^{TOTEM} vs p_{y}^{CMS} TT/BB",100,-2.,2.,100,-2.,2.);
  
  //...Luiz
  histosTH1F["hdpyAll"] = fs->make<TH1F>("hdpyAll"  ,"#Deltap_{Y} CMS-TOTEM",500,-0.5,0.5);
  histosTH1F["hdpy"] = fs->make<TH1F>("hdpy"     ,"#Deltap_{Y} CMS-TOTEM",500,-0.5,0.5);
  histosTH1F["hdpy0"] = fs->make<TH1F>("hdpy0"     ,"#Deltap_{Y} CMS-TOTEM Q=0",500,-0.5,0.5);
  //**histosTH1F["hdpy_diag"] = fs->make<TH1F>("hdpy_diag","#Deltap_{Y} CMS-TOTEM TB/BT",500,-0.5,0.5);
  //**histosTH1F["hdpy_ttbb"] = fs->make<TH1F>("hdpy_ttbb","#Deltap_{Y} CMS-TOTEM TT/BB",500,-0.5,0.5);
  
  histosTH2F["h2dimdpxAll"] = fs->make<TH2F>("h2dimdpxAll", "p_{x}^{TOTEM} vs p_{x}^{CMS}",200,-2.,2.,200,-2., 2.);
  histosTH2F["h2dimdpx"] = fs->make<TH2F>("h2dimdpx","p_{x}^{TOTEM} vs p_{x}^{CMS}",200,-2.,2.,200,-2.,2.);
  //***histosTH2F["h2dimdpx_diag"] = fs->make<TH2F>("h2dimdpx_diag","p_{x}^{TOTEM} vs p_{x}^{CMS} diag",100,-2.,2.,100,-2.,2.);
  //***histosTH2F["h2dimdpx_ttbb"] = fs->make<TH2F>("h2dimdpx_ttbb","p_{x}^{TOTEM} vs p_{x}^{CMS} TT/BB",100,-2.,2.,100,-2.,2.);
  
  //...Luiz
  histosTH1F["hdpxAll"] = fs->make<TH1F>("hdpxAll", "#Deltap_{X} CMS-TOTEM",500,-0.5,0.5);
  histosTH1F["hdpx"] = fs->make<TH1F>("hdpx", "#Deltap_{X} CMS-TOTEM",500,-0.5,0.5);
  histosTH1F["hdpx0"] = fs->make<TH1F>("hdpx0", "#Deltap_{X} CMS-TOTEM Q=0",500,-0.5,0.5);
  //**histosTH1F["hdpx_diag"] = fs->make<TH1F>("hdpx_diag", "#Deltap_{X} CMS-TOTEM TB/BT",500,-0.5,0.5);
  //**histosTH1F["hdpx_ttbb"] = fs->make<TH1F>("hdpx_ttbb", "#Deltap_{X} CMS-TOTEM TT/BB",500,-0.5,0.5);

  //...checking!
  //**histosTH1F["hcmspx"] = fs->make<TH1F>("hcmspx","CMSpx",500,-0.5,0.5);
  //**histosTH1F["hcmspy"] = fs->make<TH1F>("hcmspy","CMSpy",500,-0.5,0.5);
  //**histosTH1F["hcmspxk"] = fs->make<TH1F>("hcmspxk","CMSpxK",500,-0.5,0.5);
  //**histosTH1F["hcmspyk"] = fs->make<TH1F>("hcmspyk","CMSpyK",500,-0.5,0.5);

  //------------------
  //***histosTH2F["h2dimxVtxRL"] = fs->make<TH2F>("h2dimxVtxRL","xVtxL vs xVtxR (m)",1000,-0.004,0.001,1000,-0.004,0.001);
  //***histosTH2F["h2dimxVtxcmsR"] = fs->make<TH2F>("h2dimxVtxcmsR","xVtxCMS vs xVtxR (cm)",300,-0.3,0.3,400,-0.3,0.5);
  //***histosTH2F["h2dimxVtxcmsL"] = fs->make<TH2F>("h2dimxVtxcmsL","xVtxCMS vs xVtxL (cm)",300,-0.3,0.3,400,-0.3,0.5);
  //***histosTH2F["h2dimxVtxcmsRL"] = fs->make<TH2F>("h2dimxVtxcmsRL","xVtxCMS vs xVtxRL (cm)",300,-0.3,0.3,400,-0.3,0.5);

  //***histosTH2F["h2dimxVtxcmsR2"] = fs->make<TH2F>("h2dimxVtxcmsR2","xVtxCMS vs xVtxR (cm) (|xVtxL-xVtxR|<3e-5)",300,-0.3,0.3,400,-0.3,0.5);
  //***histosTH2F["h2dimxVtxcmsL2"] = fs->make<TH2F>("h2dimxVtxcmsL2","xVtxCMS vs xVtxL (cm) (|xVtxL-xVtxR|<3e-5)",300,-0.3,0.3,400,-0.3,0.5);
  //***histosTH2F["h2dimxVtxcmsRL2"] = fs->make<TH2F>("h2dimxVtxcmsRL2","xVtxCMS vs xVtxRL (cm)",300,-0.3,0.3,400,-0.3,0.5);

  //***histosTH2F["h2dimxVtx_zVtx_CT"] = fs->make<TH2F>("h2dimxVtx_zVtx_CT","xVtxCMS-xVtxTOTEM vs zVtx (cm)",300,-20.,20.,400,-0.3,0.5);
  //***histosTH2F["h2dimxVtx_zVtx_C"] = fs->make<TH2F>("h2dimxVtx_zVtx_C","xVtxCMS vs zVtx (cm)",300,-20.,20.,400,-0.3,0.5);
  //***histosTH2F["h2dimxVtx_zVtx_T"] = fs->make<TH2F>("h2dimxVtx_zVtx_T","xVtxTOTEM vs zVtx (cm)",300,-20.,20.,400,-0.3,0.5);

  //**histosTH1F["hxVtxRL"] = fs->make<TH1F>("hxVtxRL","xVtxR-xVtxL (m)",300,-0.0003,0.0003);

  //...Luiz
  //**histosTH1F["hxVtxcmsR"] = fs->make<TH1F>("hxVtxcmsR","xVtxCMS-xVtxR (cm)",500,-0.5,0.5);
  //**histosTH1F["hxVtxcmsL"] = fs->make<TH1F>("hxVtxcmsL","xVtxCMS-xVtxL (cm)",500,-0.5,0.5);
  //**histosTH1F["hxVtxcmsRL"] = fs->make<TH1F>("hxVtxcmsRL","xVtxCMS-xVtxTOTEM (cm)",500,-0.5,0.5);

  //**histosTH1F["hxVtxRL_diag"] = fs->make<TH1F>("hxVtxRL_diag","xVtxR-xVtxL (m)",300,-0.0003,0.0003);

  //...Luiz
  //**histosTH1F["hxVtxcmsR_diag"] = fs->make<TH1F>("hxVtxcmsR_diag","xVtxCMS-xVtxR (cm)",500,-0.5,0.5);
  //**histosTH1F["hxVtxcmsL_diag"] = fs->make<TH1F>("hxVtxcmsL_diag","xVtxCMS-xVtxL (cm)",500,-0.5,0.5);
  //**histosTH1F["hxVtxcmsRL_diag"] = fs->make<TH1F>("hxVtxcmsRL_diag","xVtxCMS-xVtxTOTEM (cm)",500,-0.5,0.5);

  //**histosTH1F["hxVtxRL_ttbb"] = fs->make<TH1F>("hxVtxRL_ttbb","xVtxR-xVtxL (m)",300,-0.0003,0.0003);

  //...Luiz
  //**histosTH1F["hxVtxcmsR_ttbb"] = fs->make<TH1F>("hxVtxcmsR_ttbb","xVtxCMS-xVtxR (cm)",500,-0.5,0.5);
  //**histosTH1F["hxVtxcmsL_ttbb"] = fs->make<TH1F>("hxVtxcmsL_ttbb","xVtxCMS-xVtxL (cm)",500,-0.5,0.5);
  //**histosTH1F["hxVtxcmsRL_ttbb"] = fs->make<TH1F>("hxVtxcmsRL_ttbb","xVtxCMS-xVtxTOTEM (cm)",500,-0.5,0.5);
  
  //...Luiz

  //  histosTH2F["hdedx"] = fs->make<TH2F>("hdedx","dE/dx vs p", 1000, 0.,20.,1000, 0.,200.);
  //...no cuts
  histosTH2F["hdedxPIX"] = fs->make<TH2F>("hdedxPIX","dE/dx vs p PIX", 1000, 0.,10.,1000, 0.,50.);
  histosTH2F["hdedxPIXno"] = fs->make<TH2F>("hdedxPIXno","dE/dx vs p PIX", 1000, 0.,10.,1000, 0.,50.);
  histosTH2F["hdedxPIXnopi"] = fs->make<TH2F>("hdedxPIXnopi","dE/dx vs p PIX pions", 1000, 0.,10.,1000, 0.,50.);
  histosTH2F["hdedxPIXnok"] = fs->make<TH2F>("hdedxPIXnok","dE/dx vs p PIX kaons", 1000, 0.,10.,1000, 0.,50.);
  histosTH2F["hdedxPIXnop"] = fs->make<TH2F>("hdedxPIXnop","dE/dx vs p PIX protons", 1000, 0.,10.,1000, 0.,50.);
  histosTH2F["hdedxPIXnou"] = fs->make<TH2F>("hdedxPIXnou","dE/dx vs p PIX unkown", 1000, 0.,10.,1000, 0.,50.);
  histosTH2F["hdedxPIXnoe"] = fs->make<TH2F>("hdedxPIXnoe","dE/dx vs p PIX electrons", 1000, 0.,10.,1000, 0.,50.);
  histosTH2F["hdedxPIXnokpi"] = fs->make<TH2F>("hdedxPIXnokpi","dE/dx vs p PIX K or #pi", 1000, 0.,10.,1000, 0.,50.);
  //...fiducial allcuts Q=0 nvtx=1
  histosTH2F["hdedxPIX4"] = fs->make<TH2F>("hdedxPIX4","dE/dx vs p PIX 4-track", 1000, 0.,10.,1000, 0.,50.);
  histosTH2F["hdedxPIX4pi"] = fs->make<TH2F>("hdedxPIX4pi","dE/dx vs p PIX 4-track pions", 1000, 0.,10.,1000, 0.,50.);
  histosTH2F["hdedxPIX4k"] = fs->make<TH2F>("hdedxPIX4k","dE/dx vs p PIX 4-track kaons", 1000, 0.,10.,1000, 0.,50.);
  histosTH2F["hdedxPIX4p"] = fs->make<TH2F>("hdedxPIX4p","dE/dx vs p PIX 4-track protons", 1000, 0.,10.,1000, 0.,50.);
  histosTH2F["hdedxPIX4u"] = fs->make<TH2F>("hdedxPIX4u","dE/dx vs p PIX 4-track unkown", 1000, 0.,10.,1000, 0.,50.);
  histosTH2F["hdedxPIX4e"] = fs->make<TH2F>("hdedxPIX4e","dE/dx vs p PIX 4-track electrons", 1000, 0.,10.,1000, 0.,50.);
  histosTH2F["hdedxPIX4kpi"] = fs->make<TH2F>("hdedxPIX4kpi","dE/dx vs p PIX 4-track K^{#pm}#pi^{#mp}#pi^{+}#pi^{-}", 1000, 0.,10.,1000, 0.,50.);
  //  histosTH2F["hdedxpixel"] = fs->make<TH2F>("hdedx","dE/dx vs p", 1000, 0.,20.,1000, 0.,200.);
  //***histosTH2F["hdedxvee11"] = fs->make<TH2F>("hdedxvee11","dE/dx vs p type:11", 1000, 0.,20.,1000, 0.,200.);
  //***histosTH2F["hdedxvee02"] = fs->make<TH2F>("hdedxvee02","dE/dx vs p type:02", 1000, 0.,20.,1000, 0.,200.);
  //***histosTH2F["hdedxvee01"] = fs->make<TH2F>("hdedxvee01","dE/dx vs p type:01", 1000, 0.,20.,1000, 0.,200.);
  //***histosTH2F["hdedxrejKp"] = fs->make<TH2F>("hdedxrejKp","dE/dx vs p rejecting K or p", 1000, 0.,20.,1000, 0.,200.);
  //***histosTH2F["hdedx2rejKp"] = fs->make<TH2F>("hdedx2rejKp","dE/dx vs p nvtx=1 rejecting K or p", 1000, 0.,20.,1000, 0.,200.);
  //***histosTH2F["hdedxrejKpu"] = fs->make<TH2F>("hdedxrejKpu","dE/dx vs p rejecting K or p or unknown", 1000, 0.,20.,1000, 0.,200.);
  //***histosTH2F["hdedx2rejKpu"] = fs->make<TH2F>("hdedx2rejKpu","dE/dx vs p nvtx=1 rejecting K or p or unknown", 1000, 0.,20.,1000, 0.,200.);

  //...Luiz
  //***histosTH2F["hlndedx"]  = fs->make<TH2F>("hlndedx","ln dE/dx vs p", 500, 0.,5.,1000, 0.,5.);
  //***histosTH2F["hl10dedx"] = fs->make<TH2F>("hl10dedx","log10 dE/dx vs p", 500, 0.,5.,1000, 0.,5.);

  /*
  // eneMK0
  histosTH1F["henemk012"] = fs->make<TH1F>("henemk012","eneMK0 #pi1#pi2", 5000, 0.,50.);
  histosTH1F["henemk034"] = fs->make<TH1F>("henemk034","eneMK0 #pi3#pi4", 5000, 0.,50.);
  */
  
  //...proton transverse momentum
  histosTH1F["hprotonlpt"] = fs->make<TH1F>("hprotonlpt","proton left p_{T}",2000,-20.,20.);
  histosTH1F["hprotonrpt"] = fs->make<TH1F>("hprotonrpt","proton right p_(T)",2000,-20.,20.);
  //
  histosTH1F["hprotonlpx"] = fs->make<TH1F>("hprotonlpx","proton left p_{x}",2000,-20.,20.);
  histosTH1F["hprotonlpy"] = fs->make<TH1F>("hprotonlpy","proton left p_{y}",2000,-20.,20.);
  histosTH1F["hprotonrpx"] = fs->make<TH1F>("hprotonrpx","proton right p_{x}",2000,-20.,20.);
  histosTH1F["hprotonrpy"] = fs->make<TH1F>("hprotonrpy","proton right p_{y}",2000,-20.,20.);

  //...4-momentum transfer squared
  histosTH1F["ht1"] = fs->make<TH1F>("ht1","|-t1|",1000,0.,5.);
  histosTH1F["ht2"] = fs->make<TH1F>("ht2","|-t2|",1000,0.,5.);
  //histosTH1F["ht1t2"] = fs->make<TH1F>("ht1t2","|-(t1+t2)|",1000,0.,5.);
  histosTH1F["ht12"] = fs->make<TH1F>("ht12","|-(t1+t2)|",1000,0.,5.);
  //histosTH1F["ht"] = fs->make<TH1F>("ht","|-t|",1000,0.,5.); //.............WRONG!
  //
  histosTH2F["h2dimt1t2"]  = fs->make<TH2F>("h2dimt1t2","t1 vs t2", 1000,0.,5.,1000,0.,5.);
  // t
  histosTH2F["h2dimt1protonpt"]  = fs->make<TH2F>("h2dimt1protonpt","t1 vs proton left pt", 1000,0.,5.,1000,0.,5.);
  histosTH2F["h2dimt2protonpt"]  = fs->make<TH2F>("h2dimt2protonpt","t2 vs proton right pt", 1000,0.,5.,1000,0.,5.);
  // xi
  histosTH2F["h2dimxi1protonpt"]  = fs->make<TH2F>("h2dimxi1protonpt","xi1 vs proton left pt", 1000,0.,5.,1000,0.,5.);
  histosTH2F["h2dimxi2protonpt"]  = fs->make<TH2F>("h2dimxi2protonpt","xi2 vs proton right pt", 1000,0.,5.,1000,0.,5.);
  //
  // acceptance
  //histosTH2F["h2dimt1pt"]  = fs->make<TH2F>("h2dimt1pt","t1 vs pt", 1000,0.,5.,1000,0.,5.);
  //histosTH2F["h2dimt2pt"]  = fs->make<TH2F>("h2dimt2pt","t2 vs pt", 1000,0.,5.,1000,0.,5.);
  //histosTH2F["h2dimt12pt"] = fs->make<TH2F>("h2dimt12pt","t12 vs pt", 1000,0.,5.,1000,0.,5.);
  histosTH2F["h2dimt1pt"]  = fs->make<TH2F>("h2dimt1pt","t1 vs pt", 1000,0.,50.,1000,0.,50.);
  histosTH2F["h2dimt2pt"]  = fs->make<TH2F>("h2dimt2pt","t2 vs pt", 1000,0.,50.,1000,0.,50.);
  histosTH2F["h2dimt12pt"] = fs->make<TH2F>("h2dimt12pt","t12 vs pt", 1000,0.,50.,1000,0.,50.);
 
  //
  // xi 200,-0.5,0.5
  //
  //...relative momentum loss
  histosTH1F["hxi1"] = fs->make<TH1F>("hxi1","xi1",1000,-0.05,0.05);
  histosTH1F["hxi2"] = fs->make<TH1F>("hxi2","xi2",1000,-0.05,0.05);
  //
  histosTH2F["h2dimxi1xi2"]  = fs->make<TH2F>("h2dimxi1xi2","xi1 vs xi2",1000,-0.05,0.05,1000,-0.05,0.05);
  //
  histosTH1F["hxi1_1"] = fs->make<TH1F>("hxi1_1","xi1_1",100,-0.05,0.05);
  histosTH1F["hxi1_2"] = fs->make<TH1F>("hxi1_2","xi1_2",100,-0.05,0.05);
  histosTH1F["hxi2_1"] = fs->make<TH1F>("hxi2_1","xi2_1",100,-0.05,0.05);
  histosTH1F["hxi2_2"] = fs->make<TH1F>("hxi2_2","xi2_2",100,-0.05,0.05);

  // phi1, phi2, dphi
  histosTH1F["hphi1"] = fs->make<TH1F>("hphi1","phi1",320,-3.2,3.2);
  histosTH1F["hphi2"] = fs->make<TH1F>("hphi2","phi2",320,-3.2,3.2);
  histosTH1F["hdphi"] = fs->make<TH1F>("hdphi","dphi in rads",320,0.,3.2);
  histosTH1F["hdphig"] = fs->make<TH1F>("hdphig","dphi in degrees",180,0.,180.);
  histosTH1F["hdphigcor"] = fs->make<TH1F>("hdphigcor","dphi in degrees efficiency",180,0.,180.);
  //histosTH1F["hdphif"] = fs->make<TH1F>("hdphif","dphi fraction in degrees",1800,0.,180.);

  // dphi with mass cut ...veeno02
  histosTH1F["hdphig1370"] = fs->make<TH1F>("hdphig1370","dphi f0(1370)",180,0.,180.);
  histosTH1F["hdphig1549"] = fs->make<TH1F>("hdphig1549","dphi f0(1500)",180,0.,180.);
  histosTH1F["hdphig1731"] = fs->make<TH1F>("hdphig1731","dphi f0(1710)",180,0.,180.);

  histosTH2F["h2dimdphigm"]  = fs->make<TH2F>("h2dimdphigm","dphi vs m",180,0.,180.,250,0.,5.);
  histosTH2F["h2dimt12m"]  = fs->make<TH2F>("h2dimt12m","|-t12| vs m",1000,0.,5.,250,0.,5.);
  // dphi with mass cut ...veeno02 ...NO PID
  histosTH1F["hdphig1370nopid"] = fs->make<TH1F>("hdphig1370nopid","dphi f0(1370)",180,0.,180.);
  histosTH1F["hdphig1549nopid"] = fs->make<TH1F>("hdphig1549nopid","dphi f0(1500)",180,0.,180.);
  histosTH1F["hdphig1731nopid"] = fs->make<TH1F>("hdphig1731nopid","dphi f0(1710)",180,0.,180.);
  //...effic correction
  histosTH1F["hdphig1370nopidcor"] = fs->make<TH1F>("hdphig1370nopidcor","dphi f0(1370) effic",180,0.,180.);
  histosTH1F["hdphig1549nopidcor"] = fs->make<TH1F>("hdphig1549nopidcor","dphi f0(1500) 1549 effic",180,0.,180.);
  histosTH1F["hdphig1731nopidcor"] = fs->make<TH1F>("hdphig1731nopidcor","dphi f0(1710) 1731 effic",180,0.,180.);

  histosTH2F["h2dimdphigmnopid"]  = fs->make<TH2F>("h2dimdphigmnopid","dphi vs m",180,0.,180.,250,0.,5.);
  histosTH2F["h2dimt12mnopid"]  = fs->make<TH2F>("h2dimt12mnopid","|-t12| vs m",1000,0.,5.,250,0.,5.);


  //...K*Kbar*
  //
  histosTH1F["hm4rec2OSkstar"] = fs->make<TH1F>("hm4rec2OSkstar","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.);
  histosTH1F["hm4rec2OSkstarbar"] = fs->make<TH1F>("hm4rec2OSkstarbar","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.);
  histosTH1F["hm4rec2OSKK100"] = fs->make<TH1F>("hm4rec2OSKK100","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}* study 100 MeV",massbins,0,5.);
  histosTH1F["hm4rec2OSKKwin29"] = fs->make<TH1F>("hm4rec2OSKKwin29","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}* 2#sigma = 29 MeV",massbins,0,5.);
  histosTH1F["hm4rec2OSKKwin58"] = fs->make<TH1F>("hm4rec2OSKKwin58","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}* 4#sigma = 58 MeV",massbins,0,5.);
  histosTH2F["h2dim2OSKK"] = fs->make<TH2F>("h2dim2OSKK","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.,massbins,0,5.);

  
  //...Armenteros-Podolanski
  // LAB system - tracks
  histosTH2F["hpod12"] = fs->make<TH2F>("hpod12","Armenteros-Podolanski vee12",2000,-1.,1.,2000,0.,1.);
  histosTH2F["hpod34"] = fs->make<TH2F>("hpod34","Armenteros-Podolanski vee34",2000,-1.,1.,2000,0.,1.);
  histosTH2F["hpod13"] = fs->make<TH2F>("hpod13","Armenteros-Podolanski vee13",2000,-1.,1.,2000,0.,1.);
  histosTH2F["hpod24"] = fs->make<TH2F>("hpod24","Armenteros-Podolanski vee24",2000,-1.,1.,2000,0.,1.);
  histosTH2F["hpod14"] = fs->make<TH2F>("hpod14","Armenteros-Podolanski vee14",2000,-1.,1.,2000,0.,1.);
  histosTH2F["hpod23"] = fs->make<TH2F>("hpod23","Armenteros-Podolanski vee23",2000,-1.,1.,2000,0.,1.);
  histosTH2F["hpod"] = fs->make<TH2F>("hpod","Armenteros-Podolanski analysis",2000,-1.,1.,2000,0.,1.);
  // LAB system - vees
  histosTH2F["hpod12vee"] = fs->make<TH2F>("hpod12vee","Armenteros-Podolanski vee12",2000,-1.,1.,2000,0.,1.);
  histosTH2F["hpod34vee"] = fs->make<TH2F>("hpod34vee","Armenteros-Podolanski vee34",2000,-1.,1.,2000,0.,1.);
  histosTH2F["hpod13vee"] = fs->make<TH2F>("hpod13vee","Armenteros-Podolanski vee13",2000,-1.,1.,2000,0.,1.);
  histosTH2F["hpod24vee"] = fs->make<TH2F>("hpod24vee","Armenteros-Podolanski vee24",2000,-1.,1.,2000,0.,1.);
  histosTH2F["hpod14vee"] = fs->make<TH2F>("hpod14vee","Armenteros-Podolanski vee14",2000,-1.,1.,2000,0.,1.);
  histosTH2F["hpod23vee"] = fs->make<TH2F>("hpod23vee","Armenteros-Podolanski vee23",2000,-1.,1.,2000,0.,1.);
  histosTH2F["hpodvee"] = fs->make<TH2F>("hpodvee","Armenteros-Podolanski analysis",2000,-1.,1.,2000,0.,1.);
  // ...daughter mass - tracks
  histosTH1F["hm4rec2OSpod12"] = fs->make<TH1F>("hm4rec2OSpod12","Armenteros-Podolanski trk_{1}trk_{2}",10*massbins,0.,5.);
  histosTH1F["hm4rec2OSpod34"] = fs->make<TH1F>("hm4rec2OSpod34","Armenteros-Podolanski trk_{3}trk_{4}",10*massbins,0.,5.);
  histosTH1F["hm4rec2OSpod13"] = fs->make<TH1F>("hm4rec2OSpod13","Armenteros-Podolanski trk_{1}trk_{3}",10*massbins,0.,5.);
  histosTH1F["hm4rec2OSpod24"] = fs->make<TH1F>("hm4rec2OSpod24","Armenteros-Podolanski trk_{2}trk_{4}",10*massbins,0.,5.);
  histosTH1F["hm4rec2OSpod14"] = fs->make<TH1F>("hm4rec2OSpod14","Armenteros-Podolanski trk_{1}trk_{4}",10*massbins,0.,5.);
  histosTH1F["hm4rec2OSpod23"] = fs->make<TH1F>("hm4rec2OSpod23","Armenteros-Podolanski trk_{2}trk_{3}",10*massbins,0.,5.);
  // ...daughter mass - vees
  histosTH1F["hm4rec2OSpod12vee"] = fs->make<TH1F>("hm4rec2OSpod12vee","Armenteros-Podolanski #pi_{1}#pi_{2}",10*massbins,0.,5.);
  histosTH1F["hm4rec2OSpod34vee"] = fs->make<TH1F>("hm4rec2OSpod34vee","Armenteros-Podolanski #pi_{3}#pi_{4}",10*massbins,0.,5.);
  histosTH1F["hm4rec2OSpod13vee"] = fs->make<TH1F>("hm4rec2OSpod13vee","Armenteros-Podolanski #pi_{1}#pi_{3}",10*massbins,0.,5.);
  histosTH1F["hm4rec2OSpod24vee"] = fs->make<TH1F>("hm4rec2OSpod24vee","Armenteros-Podolanski #pi_{2}#pi_{4}",10*massbins,0.,5.);
  histosTH1F["hm4rec2OSpod14vee"] = fs->make<TH1F>("hm4rec2OSpod14vee","Armenteros-Podolanski #pi_{1}#pi_{4}",10*massbins,0.,5.);
  histosTH1F["hm4rec2OSpod23vee"] = fs->make<TH1F>("hm4rec2OSpod23vee","Armenteros-Podolanski #pi_{2}#pi_{3}",10*massbins,0.,5.);
  // ...parent mass - tracks
  histosTH1F["hm4rec2OSpod1234"] = fs->make<TH1F>("hm4rec2OSpod1234","Armenteros-Podolanski trk_{1}trk_{2}#otimestrk_{3}trk_{4}",massbins,0.,5.); 
  histosTH1F["hm4rec2OSpod1324"] = fs->make<TH1F>("hm4rec2OSpod1324","Armenteros-Podolanski trk_{1}trk_{3}#otimestrk_{2}trk_{4}",massbins,0.,5.); 
  histosTH1F["hm4rec2OSpod1423"] = fs->make<TH1F>("hm4rec2OSpod1423","Armenteros-Podolanski trk_{1}trk_{4}#otimestrk_{2}trk_{3}",massbins,0.,5.); 
  // ...parent mass - vees
  histosTH1F["hm4rec2OSpod1234vee"] = fs->make<TH1F>("hm4rec2OSpod1234vee","Armenteros-Podolanski #pi_{1}#pi_{2}#otimes#pi_{3}#pi_{4}",massbins,0.,5.); 
  histosTH1F["hm4rec2OSpod1324vee"] = fs->make<TH1F>("hm4rec2OSpod1324vee","Armenteros-Podolanski #pi_{1}#pi_{3}#otimes#pi_{2}#pi_{4}",massbins,0.,5.); 
  histosTH1F["hm4rec2OSpod1423vee"] = fs->make<TH1F>("hm4rec2OSpod1423vee","Armenteros-Podolanski #pi_{1}#pi_{4}#otimes#pi_{2}#pi_{3}",massbins,0.,5.); 
  //
  // rho/w
  // LAB system - tracks
  histosTH2F["hpod12rw"] = fs->make<TH2F>("hpod12rw","Armenteros-Podolanski vee12",2000,-1.,1.,2000,0.,1.);
  histosTH2F["hpod34rw"] = fs->make<TH2F>("hpod34rw","Armenteros-Podolanski vee34",2000,-1.,1.,2000,0.,1.);
  histosTH2F["hpod13rw"] = fs->make<TH2F>("hpod13rw","Armenteros-Podolanski vee13",2000,-1.,1.,2000,0.,1.);
  histosTH2F["hpod24rw"] = fs->make<TH2F>("hpod24rw","Armenteros-Podolanski vee24",2000,-1.,1.,2000,0.,1.);
  histosTH2F["hpod14rw"] = fs->make<TH2F>("hpod14rw","Armenteros-Podolanski vee14",2000,-1.,1.,2000,0.,1.);
  histosTH2F["hpod23rw"] = fs->make<TH2F>("hpod23rw","Armenteros-Podolanski vee23",2000,-1.,1.,2000,0.,1.);
  histosTH2F["hpodrw"] = fs->make<TH2F>("hpodrw","Armenteros-Podolanski analysis",2000,-1.,1.,2000,0.,1.);
  // LAB system - vees
  histosTH2F["hpod12veerw"] = fs->make<TH2F>("hpod12veerw","Armenteros-Podolanski vee12",2000,-1.,1.,2000,0.,1.);
  histosTH2F["hpod34veerw"] = fs->make<TH2F>("hpod34veerw","Armenteros-Podolanski vee34",2000,-1.,1.,2000,0.,1.);
  histosTH2F["hpod13veerw"] = fs->make<TH2F>("hpod13veerw","Armenteros-Podolanski vee13",2000,-1.,1.,2000,0.,1.);
  histosTH2F["hpod24veerw"] = fs->make<TH2F>("hpod24veerw","Armenteros-Podolanski vee24",2000,-1.,1.,2000,0.,1.);
  histosTH2F["hpod14veerw"] = fs->make<TH2F>("hpod14veerw","Armenteros-Podolanski vee14",2000,-1.,1.,2000,0.,1.);
  histosTH2F["hpod23veerw"] = fs->make<TH2F>("hpod23veerw","Armenteros-Podolanski vee23",2000,-1.,1.,2000,0.,1.);
  histosTH2F["hpodveerw"] = fs->make<TH2F>("hpodveewr","Armenteros-Podolanski analysis",2000,-1.,1.,2000,0.,1.);
  // ...daughter mass - tracks
  histosTH1F["hm4rec2OSpod12rw"] = fs->make<TH1F>("hm4rec2OSpod12rw","Armenteros-Podolanski trk_{1}trk_{2}",10*massbins,0.,5.);
  histosTH1F["hm4rec2OSpod34rw"] = fs->make<TH1F>("hm4rec2OSpod34rw","Armenteros-Podolanski trk_{3}trk_{4}",10*massbins,0.,5.);
  histosTH1F["hm4rec2OSpod13rw"] = fs->make<TH1F>("hm4rec2OSpod13rw","Armenteros-Podolanski trk_{1}trk_{3}",10*massbins,0.,5.);
  histosTH1F["hm4rec2OSpod24rw"] = fs->make<TH1F>("hm4rec2OSpod24rw","Armenteros-Podolanski trk_{2}trk_{4}",10*massbins,0.,5.);
  histosTH1F["hm4rec2OSpod14rw"] = fs->make<TH1F>("hm4rec2OSpod14rw","Armenteros-Podolanski trk_{1}trk_{4}",10*massbins,0.,5.);
  histosTH1F["hm4rec2OSpod23rw"] = fs->make<TH1F>("hm4rec2OSpod23rw","Armenteros-Podolanski trk_{2}trk_{3}",10*massbins,0.,5.);
  // ...daughter mass - vees
  histosTH1F["hm4rec2OSpod12veerw"] = fs->make<TH1F>("hm4rec2OSpod12veerw","Armenteros-Podolanski #pi_{1}#pi_{2}",10*massbins,0.,5.);
  histosTH1F["hm4rec2OSpod34veerw"] = fs->make<TH1F>("hm4rec2OSpod34veerw","Armenteros-Podolanski #pi_{3}#pi_{4}",10*massbins,0.,5.);
  histosTH1F["hm4rec2OSpod13veerw"] = fs->make<TH1F>("hm4rec2OSpod13veerw","Armenteros-Podolanski #pi_{1}#pi_{3}",10*massbins,0.,5.);
  histosTH1F["hm4rec2OSpod24veerw"] = fs->make<TH1F>("hm4rec2OSpod24veerw","Armenteros-Podolanski #pi_{2}#pi_{4}",10*massbins,0.,5.);
  histosTH1F["hm4rec2OSpod14veerw"] = fs->make<TH1F>("hm4rec2OSpod14veerw","Armenteros-Podolanski #pi_{1}#pi_{4}",10*massbins,0.,5.);
  histosTH1F["hm4rec2OSpod23veerw"] = fs->make<TH1F>("hm4rec2OSpod23veerw","Armenteros-Podolanski #pi_{2}#pi_{3}",10*massbins,0.,5.);
  // ...parent mass - tracks
  histosTH1F["hm4rec2OSpod1234rw"] = fs->make<TH1F>("hm4rec2OSpod1234rw","Armenteros-Podolanski trk_{1}trk_{2}#otimestrk_{3}trk_{4}",massbins,0.,5.); 
  histosTH1F["hm4rec2OSpod1324rw"] = fs->make<TH1F>("hm4rec2OSpod1324rw","Armenteros-Podolanski trk_{1}trk_{3}#otimestrk_{2}trk_{4}",massbins,0.,5.); 
  histosTH1F["hm4rec2OSpod1423rw"] = fs->make<TH1F>("hm4rec2OSpod1423rw","Armenteros-Podolanski trk_{1}trk_{4}#otimestrk_{2}trk_{3}",massbins,0.,5.); 
  // ...parent mass - vees
  histosTH1F["hm4rec2OSpod1234veerw"] = fs->make<TH1F>("hm4rec2OSpod1234veerw","Armenteros-Podolanski #pi_{1}#pi_{2}#otimes#pi_{3}#pi_{4}",massbins,0.,5.); 
  histosTH1F["hm4rec2OSpod1324veerw"] = fs->make<TH1F>("hm4rec2OSpod1324veerw","Armenteros-Podolanski #pi_{1}#pi_{3}#otimes#pi_{2}#pi_{4}",massbins,0.,5.); 
  histosTH1F["hm4rec2OSpod1423veerw"] = fs->make<TH1F>("hm4rec2OSpod1423veerw","Armenteros-Podolanski #pi_{1}#pi_{4}#otimes#pi_{2}#pi_{3}",massbins,0.,5.); 
  // LIFETIME
  histosTH1F["hlife1A"] = fs->make<TH1F>("hlife1A","Armenteros-Podolanski lifetime 12",2000,0.,20.);
  histosTH1F["hlife2A"] = fs->make<TH1F>("hlife2A","Armenteros-Podolanski lifetime 12",2000,0.,20.);
  histosTH1F["hlife12"] = fs->make<TH1F>("hlife12","Armenteros-Podolanski lifetime 12",2000,0.,20.);
  //
  histosTH1F["hlife3A"] = fs->make<TH1F>("hlife3A","Armenteros-Podolanski lifetime 34",2000,0.,20.);
  histosTH1F["hlife4A"] = fs->make<TH1F>("hlife4A","Armenteros-Podolanski lifetime 34",2000,0.,20.);
  histosTH1F["hlife34"] = fs->make<TH1F>("hlife34","Armenteros-Podolanski lifetime 34",2000,0.,20.);
  //
  histosTH1F["hlife1B"] = fs->make<TH1F>("hlife1B","Armenteros-Podolanski lifetime 13",2000,0.,20.);
  histosTH1F["hlife3B"] = fs->make<TH1F>("hlife3B","Armenteros-Podolanski lifetime 13",2000,0.,20.);
  histosTH1F["hlife13"] = fs->make<TH1F>("hlife13","Armenteros-Podolanski lifetime 13",2000,0.,20.);
  //
  histosTH1F["hlife2B"] = fs->make<TH1F>("hlife2B","Armenteros-Podolanski lifetime 24",2000,0.,20.);
  histosTH1F["hlife4B"] = fs->make<TH1F>("hlife4B","Armenteros-Podolanski lifetime 24",2000,0.,20.);
  histosTH1F["hlife24"] = fs->make<TH1F>("hlife24","Armenteros-Podolanski lifetime 24",2000,0.,20.);
  //
  histosTH1F["hlife1C"] = fs->make<TH1F>("hlife1C","Armenteros-Podolanski lifetime 14",2000,0.,20.);
  histosTH1F["hlife4C"] = fs->make<TH1F>("hlife4C","Armenteros-Podolanski lifetime 14",2000,0.,20.);
  histosTH1F["hlife14"] = fs->make<TH1F>("hlife14","Armenteros-Podolanski lifetime 14",2000,0.,20.);
  //
  histosTH1F["hlife2C"] = fs->make<TH1F>("hlife2C","Armenteros-Podolanski lifetime 23",2000,0.,20.);
  histosTH1F["hlife3C"] = fs->make<TH1F>("hlife3C","Armenteros-Podolanski lifetime 23",2000,0.,20.);
  histosTH1F["hlife23"] = fs->make<TH1F>("hlife23","Armenteros-Podolanski lifetime 23",2000,0.,20.);
  //
  histosTH1F["hlife"] = fs->make<TH1F>("hlife","Armenteros-Podolanski lifetime 12",2000,0.,20.);
  histosTH1F["hlifecor"] = fs->make<TH1F>("hlifecor","Armenteros-Podolanski lifetime 12",2000,0.,20.);
  
/*
//1234

histosTH1F["hm4rec2OSk1pi2"] = fs->make<TH1F>("hm4rec2OSk1pi2","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.);
histosTH1F["hm4rec2OSk3pi4b"] = fs->make<TH1F>("hm4rec2OSk3pi4b","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.);
histosTH1F["hm4rec2OSk1p2k3p4b"] = fs->make<TH1F>("hm4rec2OSk1p2k3p4b","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.);
histosTH2F["h2dim2OSk1p2k3p4b"] = fs->make<TH2F>("h2dim2OSk1p2k3p4b","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.,massbins,0,5.);

histosTH1F["hm4rec2OSk1pi2b"] = fs->make<TH1F>("hm4rec2OSk1pi2b","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.);
histosTH1F["hm4rec2OSk3pi4"] = fs->make<TH1F>("hm4rec2OSk3pi4","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.);
histosTH1F["hm4rec2OSk1p2bk3p4"] = fs->make<TH1F>("hm4rec2OSk1p2bk3p4","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.);
histosTH2F["h2dim2OSk1p2bk3p4"] = fs->make<TH2F>("h2dim2OSk1p2bk3p4","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.,massbins,0,5.);

	   //hm4rec2OSk1pi2
histosTH1F["hm4rec2OSpi3k4b"] = fs->make<TH1F>("hm4rec2OSpi3k4b","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.);
histosTH1F["hm4rec2OSk1p2p3k4b"] = fs->make<TH1F>("hm4rec2OSk1p2p3k4b","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.);
histosTH2F["h2dim2OSk1p2p3k4b"] = fs->make<TH2F>("h2dim2OSk1p2p3k4b","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.,massbins,0,5.);

           //hm4rec2OSk1pi2b
histosTH1F["hm4rec2OSpi3k4"] = fs->make<TH1F>("hm4rec2OSpi3k4","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.);
histosTH1F["hm4rec2OSk1p2bp3k4"] = fs->make<TH1F>("hm4rec2OSk1p2bp3k4","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.);
histosTH2F["h2dim2OSk1p2bp3k4"] = fs->make<TH2F>("h2dim2OSk1p2bp3k4","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.,massbins,0,5.);

histosTH1F["hm4rec2OSpi1k2b"] = fs->make<TH1F>("hm4rec2OSpi1k2b","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.);
	   //hm4rec2OSpi3k4
histosTH1F["hm4rec2OSp1k2bp3k4"] = fs->make<TH1F>("hm4rec2OSp1k2bp3k4","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.);
histosTH2F["h2dim2OSp1k2bp3k4"] = fs->make<TH2F>("h2dim2OSp1k2bp3k4","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.,massbins,0,5.);

histosTH1F["hm4rec2OSpi1k2"] = fs->make<TH1F>("hm4rec2OSpi1k2","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.);
histosTH1F["hm4rec2OSpi3k4b"] = fs->make<TH1F>("hm4rec2OSpi3k4b","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.);
histosTH1F["hm4rec2OSp1k2p3k4b"] = fs->make<TH1F>("hm4rec2OSp1k2p3k4b","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.);
histosTH2F["h2dim2OSp1k2p3k4b"] = fs->make<TH2F>("h2dim2OSp1k2p3k4b","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.,massbins,0,5.);

           //hm4rec2OSpi1k2
           //hm4rec2OSk3pi4b
histosTH1F["hm4rec2OSp1k2k3p4b"] = fs->make<TH1F>("hm4rec2OSp1k2k3p4b","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.);
histosTH2F["h2dim2OSp1k2k3p4b"] = fs->make<TH2F>("h2dim2OSp1k2k3p4b","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.,massbins,0,5.);

           //hm4rec2OSpi1k2b
           //hm4rec2OSk3pi4
histosTH1F["hm4rec2OSp1k2bk3p4"] = fs->make<TH1F>("hm4rec2OSp1k2bk3p4","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.);
histosTH2F["h2dim2OSp1k2bk3p4"] = fs->make<TH2F>("h2dim2OSp1k2bk3p4","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.,massbins,0,5.);


//1324

histosTH1F["hm4rec2OSk1pi3"] = fs->make<TH1F>("hm4rec2OSk1pi3","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.);
histosTH1F["hm4rec2OSk2pi4b"] = fs->make<TH1F>("hm4rec2OSk2pi4b","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.);
histosTH1F["hm4rec2OSk1p3k2p4b"] = fs->make<TH1F>("hm4rec2OSk1p3k2p4b","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.);
histosTH2F["h2dim2OSk1p3k2p4b"] = fs->make<TH2F>("h2dim2OSk1p3k2p4b","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.,massbins,0,5.);
						  
histosTH1F["hm4rec2OSk1pi3b"] = fs->make<TH1F>("hm4rec2OSk1pi3b","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.);
histosTH1F["hm4rec2OSk2pi4"] = fs->make<TH1F>("hm4rec2OSk2pi4","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.);
histosTH1F["hm4rec2OSk1p3bk2p4"] = fs->make<TH1F>("hm4rec2OSk1p3bk2p4","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.);
						  
histosTH2F["h2dim2OSk1p3bk2p4"] = fs->make<TH2F>("h2dim2OSk1p3bk2p4","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.,massbins,0,5.);
						  
           //hm4rec2OSk1pi3
histosTH1F["hm4rec2OSpi2k4b"] = fs->make<TH1F>("hm4rec2OSpi2k4b","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.);
histosTH1F["hm4rec2OSk1p3pi2k4b"] = fs->make<TH1F>("hm4rec2OSk1p3pi2k4b","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.);
histosTH2F["h2dim2OSk1p3p2k4b"] = fs->make<TH2F>("h2dim2OSk1p3p2k4b","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.,massbins,0,5.);
						  
           //hm4rec2OSk1pi3b
histosTH1F["hm4rec2OSpi2k4"] = fs->make<TH1F>("hm4rec2OSpi2k4","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.);
histosTH1F["hm4rec2OSk1p3bpi2k4"] = fs->make<TH1F>("hm4rec2OSk1p3bpi2k4","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.);
histosTH2F["h2dim2OSk1p3bp2k4"] = fs->make<TH2F>("h2dim2OSk1p3bp2k4","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.,massbins,0,5.);
						  
histosTH1F["hm4rec2OSpi1k3b"] = fs->make<TH1F>("hm4rec2OSpi1k3b","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.);
           //hm4rec2OSpi2k4
histosTH1F["hm4rec2OSp1k3bp2k4"] = fs->make<TH1F>("hm4rec2OSp1k3bp2k4","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.);
histosTH2F["h2dim2OSp1k3bp2k4"] = fs->make<TH2F>("h2dim2OSp1k3bp2k4","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.,massbins,0,5.);
						  
histosTH1F["hm4rec2OSpi1k3"] = fs->make<TH1F>("hm4rec2OSpi1k3","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.);
           //hm4rec2OSpi2k4b
histosTH1F["hm4rec2OSp1k3p2k4b"] = fs->make<TH1F>("hm4rec2OSp1k3p2k4b","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.);
histosTH2F["h2dim2OSp1k3p2k4b"] = fs->make<TH2F>("h2dim2OSp1k3p2k4b","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.,massbins,0,5.);


//1423

histosTH1F["hm4rec2OSk1pi4"] = fs->make<TH1F>("hm4rec2OSk1pi4","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.);
histosTH1F["hm4rec2OSk2pi3b"] = fs->make<TH1F>("hm4rec2OSk2pi3b","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.);
histosTH1F["hm4rec2OSk1p4k2p3b"] = fs->make<TH1F>("hm4rec2OSk1p4k2p3b","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.);
histosTH2F["h2dim2OSk1p4k2p3b"] = fs->make<TH2F>("h2dim2OSk1p4k2p3b","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.,massbins,0,5.);

histosTH1F["hm4rec2OSk1pi4b"] = fs->make<TH1F>("hm4rec2OSk1pi4b","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.);
histosTH1F["hm4rec2OSk2pi3"] = fs->make<TH1F>("hm4rec2OSk2pi3","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.);
histosTH1F["hm4rec2OSk1p4bk2p3"] = fs->make<TH1F>("hm4rec2OSk1p4bk2p3","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.);
histosTH2F["h2dim2OSk1p4bk2p3"] = fs->make<TH2F>("h2dim2OSk1p4bk2p3","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.,massbins,0,5.);
						  
           //hm4rec2OSk1pi4
histosTH1F["hm4rec2OSpi2k3b"] = fs->make<TH1F>("hm4rec2OSpi2k3b","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.);
histosTH1F["hm4rec2OSk1p4p2k3b"] = fs->make<TH1F>("hm4rec2OSk1p4p2k3b","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.);
histosTH2F["h2dim2OSk1p4p2k3b"] = fs->make<TH2F>("h2dim2OSk1p4p2k3b","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.,massbins,0,5.);
						  
histosTH1F["hm4rec2OSk1pi4b"] = fs->make<TH1F>("hm4rec2OSk1pi4b","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.);
histosTH1F["hm4rec2OSpi2k3"] = fs->make<TH1F>("hm4rec2OSpi2k3","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.);
histosTH1F["hm4rec2OSk1p4bp2k3"] = fs->make<TH1F>("hm4rec2OSk1p4bp2k3","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.);
histosTH2F["h2dim2OSk1p4bp2k3"] = fs->make<TH2F>("h2dim2OSk1p4bp2k3","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.,massbins,0,5.);
						  
histosTH1F["hm4rec2OSpi1k4b"] = fs->make<TH1F>("hm4rec2OSpi1k4b","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.);
           //hm4rec2OSk2pi3
histosTH1F["hm4rec2OSp1k4bk2p3"] = fs->make<TH1F>("hm4rec2OSp1k4bk2p3","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.);
histosTH2F["h2dim2OSp1k4bk2p3"] = fs->make<TH2F>("h2dim2OSp1k4bk2p3","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.,massbins,0,5.);
						  
histosTH1F["hm4rec2OSpi1k4"] = fs->make<TH1F>("hm4rec2OSpi1k4","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.);
           //hm4rec2OSk2pi3b
histosTH1F["hm4rec2OSp1k4k2p3b"] = fs->make<TH1F>("hm4rec2OSp1k4k2p3b","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.);
histosTH2F["h2dim2OSp1k4k2p3b"] = fs->make<TH2F>("h2dim2OSp1k4k2p3b","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.,massbins,0,5.);
						  
histosTH1F["hm4rec2OSpi1k4b"] = fs->make<TH1F>("hm4rec2OSpi1k4b","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.);
           //hm4rec2OSpi2k3
histosTH1F["hm4rec2OSp1k4bp2k3"] = fs->make<TH1F>("hm4rec2OSp1k4bp2k3","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.);
histosTH2F["h2dim2OSp1k4bp2k3"] = fs->make<TH2F>("h2dim2OSp1k4bp2k3","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.,massbins,0,5.);
						  
           //hm4rec2OSpi1k4
           //hm4rec2OSpi2k3b
histosTH1F["hm4rec2OSp1k4p2k3b"] = fs->make<TH1F>("hm4rec2OSp1k4p2k3b","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.);
histosTH2F["h2dim2OSp1k4pi2k3b"] = fs->make<TH2F>("h2dim2OSp1k4pi2k3b","M_{K^{#pm}#pi^{#mp}} OS K*#bar{K}*",massbins,0,5.,massbins,0,5.);
*/

/*
    histosTH2F["hthyrtt1"]->Fill(tt1,ThyR);
    histosTH2F["hthyltt2"]->Fill(tt2,ThyL);
    histosTH2F["hthxrtt1"]->Fill(tt1,ThxR);
    histosTH2F["hthxltt2"]->Fill(tt2,ThxL);
*/

// no cuts
 histosTH2F["hthyrtt1"] = fs->make<TH2F>("hthyrtt1","#Theta_{y}^{R} vs |-t1| no cuts",2000,0,5.,2000,-0.0004,0.0004);
 histosTH2F["hthyltt2"] = fs->make<TH2F>("hthyltt2","#Theta_{y}^{L} vs |-t2| no cuts",2000,0,5.,2000,-0.0004,0.0004);
 histosTH2F["hthxrtt1"] = fs->make<TH2F>("hthxrtt1","#Theta_{x}^{R} vs |-t1| no cuts",2000,0,5.,2000,-0.0004,0.0004);
 histosTH2F["hthxltt2"] = fs->make<TH2F>("hthxltt2","#Theta_{x}^{L} vs |-t2| no cuts",2000,0,5.,2000,-0.0004,0.0004);
 //
 // inelastic
 histosTH2F["hthyrt1"]  = fs->make<TH2F>("hthyrt1","#Theta_{y}^{R} vs |-t1| inelastic",2000,0,5.,2000,-0.0004,0.0004);
 histosTH2F["hthylt2"]  = fs->make<TH2F>("hthylt2","#Theta_{y}^{L} vs |-t2| inelastic",2000,0,5.,2000,-0.0004,0.0004);
 histosTH2F["hthxrt1"]  = fs->make<TH2F>("hthxrt1","#Theta_{x}^{R} vs |-t1| inelastic",2000,0,5.,2000,-0.0004,0.0004);
 histosTH2F["hthxlt2"]  = fs->make<TH2F>("hthxlt2","#Theta_{x}^{L} vs |-t2| inelastic",2000,0,5.,2000,-0.0004,0.0004);
 // inelastic
 histosTH2F["hthyrxi1"]  = fs->make<TH2F>("hthyrxi1","#Theta_{y}^{R} vs #xi_{1} inelastic",2000,0,5.,2000,-0.0004,0.0004);
 histosTH2F["hthylxi2"]  = fs->make<TH2F>("hthylxi2","#Theta_{y}^{L} vs #xi_{2} inelastic",2000,0,5.,2000,-0.0004,0.0004);
 histosTH2F["hthxrxi1"]  = fs->make<TH2F>("hthxrxi1","#Theta_{x}^{R} vs #xi_{1} inelastic",2000,0,5.,2000,-0.0004,0.0004);
 histosTH2F["hthxlxi2"]  = fs->make<TH2F>("hthxlxi2","#Theta_{x}^{L} vs #xi_{2} inelastic",2000,0,5.,2000,-0.0004,0.0004);
    
  // <<<<<<<<<<

  std::cout<<"booked all of Luiz' histograms."<<std::endl;
  
  //--------------end of my histograms

}

// ------------ method called once each job just after ending the event loop  ------------
void
PromptAnalyzer::endJob()
{
  // this does not work ...Luiz
  // Output file
  // TFile* output = new TFile(outputFileName.c_str(),"RECREATE");
  ////TFile* output = new TFile("output.root","RECREATE");
  ////output->cd();

  //// don't include it with TFileService ...Write() and Close() are done automatically!
  ////for(map<string,TH1F*>::iterator it_histo = histosTH1F.begin();it_histo != histosTH1F.end(); ++it_histo)(*it_histo).second->Write();
  ////for(map<string,TH2F*>::iterator it_histo = histosTH2F.begin();it_histo != histosTH2F.end(); ++it_histo)(*it_histo).second->Write();

  //  for(map<string,TProfile*>::iterator it_histo = histosTProf.begin();it_histo != histosTProf.end(); ++it_histo)(*it_histo).second->Write();
  //  for(map<string,TH3F*>::iterator it_histo = histosTH3F.begin();it_histo != histosTH3F.end(); ++it_histo)(*it_histo).second->Write();
  
  ////output->Close();
  std::cout<<"ciao ciao..."<<std::endl;
}
//--------------------------------------------------
bool PromptAnalyzer::jsonLocal(int runnr, int ls){

  int accept = false;
  
  if(runnr == 319104 && ls >= 22 && ls <=176) accept = true; 
  if(runnr == 319124){
    if(ls >= 151 && ls <= 186) accept = true;
    if(ls >= 192 && ls <= 277) accept = true; //changed from 149 276
  }
  if(runnr == 319125 && ls >= 25 && ls <= 191) accept = true; 
  if(runnr == 319159 && ls >= 202 && ls <= 617) accept = true; 
  if(runnr == 319174 && ls >= 24 && ls <= 70) accept = true; 
  if(runnr == 319175 && ls >= 1 && ls <= 139) accept = true; 
  if(runnr == 319176 && ls >= 1 && ls <= 1799) accept = true; 
  if(runnr == 319177){
    if(ls >= 11 && ls <= 190) accept = true;
    if(ls >= 215 && ls <= 223) accept = true;
  }
  if(runnr == 319190){
    if(ls >= 39 && ls <= 125) accept = true;
    if(ls >= 147 && ls <= 309) accept = true;
  }
  
  if(runnr == 319222){
    if(ls >= 192 && ls <= 230) accept = true;
    if(ls >= 233 && ls <= 294) accept = true;
  }
  if(runnr == 319223 && ls >= 5 && ls <= 131) accept = true; 
  if(runnr == 319254 && ls >= 199 && ls <= 262) accept = true; 
  if(runnr == 319255 && ls >= 1 && ls <= 164) accept = true; 
  if(runnr == 319256){
    if(ls >= 1 && ls <= 38) accept = true;
    if(ls >= 41 && ls <= 726) accept = true;
  }
  if(runnr == 319262){
    if(ls == 10) accept = true;
    if(ls >= 15 && ls <= 16) accept = true;
    if(ls >= 20 && ls <= 23) accept = true;
    if(ls >= 29 && ls <= 34) accept = true;
    if(ls >= 39 && ls <= 40) accept = true;
    if(ls >= 46 && ls <= 58) accept = true;
    if(ls >= 61 && ls <= 78) accept = true;
    if(ls >= 82 && ls <= 88) accept = true;
    if(ls >= 90 && ls <= 123) accept = true;
    if(ls >= 129 && ls <= 358) accept = true;
  }
  if(runnr == 319263 && ls >= 1 && ls <= 364) accept = true; 
  if(runnr == 319264 && ls >= 1 && ls <= 57) accept = true;  
  if(runnr == 319265 && ls >= 1 && ls <= 396) accept = true;  
  if(runnr == 319266){
    if(ls >= 1 && ls <= 18) accept = true;
    if(ls >= 20 && ls <= 26) accept = true;
  }
  if(runnr == 319267 && ls >= 1 && ls <= 204) accept = true; 
  if(runnr == 319268){
    if(ls >= 1 && ls <= 185) accept = true;
    if(ls >= 187 && ls <= 462) accept = true;
    }
  if(runnr == 319270 && ls >= 1 && ls <= 205) accept = true; 

  if(runnr == 319300){
    if(ls >= 57 && ls <= 194) accept = true;
    if(ls >= 203 && ls <= 604) accept = true;
    if(ls >= 606 && ls <= 871) accept = true;
    if(ls >= 874 && ls <= 987) accept = true;
    if(ls >= 990 && ls <= 1127) accept = true;
  }
  if(runnr == 319311){
    if(ls >= 60 && ls <= 76) accept = true;
    if(ls >= 78 && ls <= 275) accept = true;
    if(ls >= 282 && ls <= 300) accept = true;
    if(ls >= 302 && ls <= 526) accept = true;
    if(ls >= 530 && ls <= 829) accept = true;
    if(ls >= 839 && ls <= 1236) accept = true;
    if(ls >= 1238 && ls <= 1489) accept = true;
    if(ls >= 1491 && ls <= 1713) accept = true;
  }

    return accept;

}

// ------------ method called when starting to processes a run  ------------
void 
PromptAnalyzer::beginRun(edm::Run const& run, edm::EventSetup const& es)
{
  
  bool changed(true);
  if (hltConfig_.init(run, es, "HLT",changed)) {
    hltConfig_.dump("Triggers");
    hltConfig_.dump("PrescaleTable"); 
    
  }

  // TTree memmory 1 TB ...Luiz  
  //TTree::SetMaxTreeSize( 1000000000000LL );

}

// ------------ method called when ending the processing of a run  ------------
void 
PromptAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
  //...Luiz
  //cutfile.Close();
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PromptAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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

//define this as a plug-in
DEFINE_FWK_MODULE(PromptAnalyzer);

// -*- C++ -*-
//
// Package:    Bspi
// Class:      Bspi
// 
//=================================================
// original author:  Jhovanny Andres Mejia        |
//         created:  October 2021                 |
//         <jhovanny.andres.mejia.guisao@cern.ch> | 
// Implemented by :  Mohammad Mobassir Ameen      |
//         <mohammad.mobassir.ameen@cern.ch>      |
//=================================================

// system include files
#include <memory>

// user include files
#include "myAnalyzers/BspiPAT/src/Bspi.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/FWLite/interface/EventBase.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "DataFormats/MuonReco/interface/MuonChamberMatch.h"

#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Candidate/interface/ShallowCloneCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Math/interface/Error.h"
#include "RecoVertex/VertexTools/interface/VertexDistance.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/CLHEP/interface/Migration.h"

#include "RecoVertex/VertexPrimitives/interface/BasicSingleVertexState.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"

#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

//
// constants, enums and typedefs
//

  typedef math::Error<3>::type CovarianceMatrix;

//
// static data member definitions
//
  const double PI = 3.141592653589793;

//
// constructors and destructor
//
Bspi::Bspi(const edm::ParameterSet& iConfig)
  :
  dimuon_Label(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("dimuons"))),
  trakCollection_label(consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("Trak"))),
  trakCollection_label_lostTracks(consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("Trak_lowpt"))),            
  genCands_(consumes<reco::GenParticleCollection>(iConfig.getParameter < edm::InputTag > ("GenParticles"))), 
  //genCands_(consumes<reco::GenParticleCollection>(iConfig.getParameter < edm::InputTag > ("prunedGenParticles"))), 
  packedGenToken_(consumes<pat::PackedGenParticleCollection>(iConfig.getParameter <edm::InputTag> ("packedGenParticles"))), 
  primaryVertices_Label(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
  triggerResults_Label(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
  BSLabel_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("bslabel"))),

  OnlyBest_(iConfig.getParameter<bool>("OnlyBest")),
  isMC_(iConfig.getParameter<bool>("isMC")),
  OnlyGen_(iConfig.getParameter<bool>("OnlyGen")),
  Trkmass_(iConfig.getParameter<double>("Trkmass")),
  TrkTrkMasscut_(iConfig.getParameter<std::vector<double> >("TrkTrkMasscut")),
  BarebMasscut_(iConfig.getParameter<std::vector<double> >("BarebMasscut")),
  bMasscut_(iConfig.getParameter<std::vector<double> >("bMasscut")),

  tree_(0), 

  mumC2(0), mumAngT(0), mumNHits(0), mumNPHits(0),
  mupC2(0), mupAngT(0), mupNHits(0), mupNPHits(0),
  mumdxy(0), mupdxy(0), mumdz(0), mupdz(0),
  muon_dca(0),

  tri_Dim25(0), tri_JpsiTk(0), tri_JpsiTkTk(0),

  mu1soft(0), mu2soft(0), mu1tight(0), mu2tight(0), 
  mu1PF(0), mu2PF(0), mu1loose(0), mu2loose(0),

  // *******************************************************
  
  nB(0), nMu(0),

  //Bspi_mass_vertex(0),
  //Bspion_mass(0),
  Bc_mass(0), Bc_pt(0), Bc_px(0), Bc_py(0), Bc_pz(0), Bc_eta(0), Bc_phi(0),Bc_ct(0), Bc_charge(0),
  Bcvtxcl(0), Bc_prob(0),

  deltaMass(0),
  T_Bspion_mass(0), Notmatch_Bspion_mass(0),

  Bc_DecayVtxX(0),    Bc_DecayVtxY(0),   Bc_DecayVtxZ(0),
  Bc_DecayVtxXE(0),   Bc_DecayVtxYE(0),  Bc_DecayVtxZE(0),

  //Bc_DecayVtx_vtxfit_X(0),   Bc_DecayVtx_vtxfit_Y(0),  Bc_DecayVtx_vtxfit_Z(0),
  //Bc_DecayVtx_vtxfit_XE(0),   Bc_DecayVtx_vtxfit_YE(0),  Bc_DecayVtx_vtxfit_ZE(0),
  //Bc_DecayVtx_vtxfit_XYE(0),   Bc_DecayVtx_vtxfit_XZE(0),  Bc_DecayVtx_vtxfit_YZE(0),


  B_mass(0), B_px(0), B_py(0), B_pz(0),
  B_pt(0), B_eta(0), B_phi(0),

  //B_pion_px(0), B_pion_py(0), B_pion_pz(0),
  pion_track_normchi2(0), pion_NumHits(0), pion_NumPixelHits(0),
  pion_dxy(0), pion_dz(0), 
  pion_dxy_(0), pion_dz_(0), 
  pion_NTrackerLayers(0), pion_NPixelLayers(0),
  B_pion_px_track(0), B_pion_py_track(0), B_pion_pz_track(0),
  B_pion_pt_track(0),
  B_pion_charge(0),
 
  B_phi_mass(0), 
  B_phi_px(0), B_phi_py(0), B_phi_pz(0),
  B_phi_pt(0), B_phi_eta(0), B_phi_phi(0),
  

  B_phi_px1(0), B_phi_py1(0), B_phi_pz1(0), B_phi_pt1(0), 
  B_phi_eta1(0), B_phi_phi1(0), 
  B_phi_px2(0), B_phi_py2(0), B_phi_pz2(0), B_phi_pt2(0),
  B_phi_eta2(0), B_phi_phi2(0), 

  B_phi_px1_track(0), B_phi_py1_track(0), B_phi_pz1_track(0), 
  B_phi_px2_track(0), B_phi_py2_track(0), B_phi_pz2_track(0),

  B_phi_charge1(0), B_phi_charge2(0),
  k1dxy(0), k2dxy(0), k1dz(0), k2dz(0),
  k1dxy_e(0), k2dxy_e(0), k1dz_e(0), k2dz_e(0),
  k1InnerHits(0), k2InnerHits(0),

  B_J_mass(0), B_J_px(0), B_J_py(0), B_J_pz(0),
  B_J_pt(0), B_J_eta(0), B_J_phi(0),
  
  B_J_pt1(0),
  B_J_eta1(0), B_J_phi1(0),
  B_J_px1(0), B_J_py1(0), B_J_pz1(0), 
  B_J_pt2(0), 
  B_J_eta2(0), B_J_phi2(0),
  B_J_px2(0), B_J_py2(0), B_J_pz2(0), 
  B_J_charge1(0), B_J_charge2(0),

  nVtx(0),
  priVtxX(0), priVtxY(0), priVtxZ(0), priVtxXE(0), priVtxYE(0), priVtxZE(0), priVtxCL(0),
  priVtxXYE(0), priVtxXZE(0), priVtxYZE(0),
  
  // ************************ ****************************************************
  
  pVtxIPX(0),  pVtxIPY(0),  pVtxIPZ(0),
  pVtxIPXE(0),  pVtxIPYE(0),  pVtxIPZE(0),  pVtxIPCL(0),
  pVtxIPXYE(0),  pVtxIPXZE(0),  pVtxIPYZE(0),
  
  B_l3d(0),  B_l3dE(0),  /*B_lxy(0),*/Bc_lxy(0), B_lxyE(0),
  B_cosalpha(0),   /*B_cosalphaxy(0),*/Bc_cosalphaxy(0), alpha(0),  //B_treco(0),   B_trecoe(0),  B_trecoxy(0), B_trecoxye(0),
  B_pvip(0), B_pviperr(0), B_pvips(0), B_pvlzip(0), B_pvlziperr(0), B_pvlzips(0),
  B_pv2ip(0), B_pv2iperr(0), B_pv2ips(0), B_pv2lzip(0), B_pv2lziperr(0), B_pv2lzips(0),

  B_l3d_pv2(0),  B_l3dE_pv2(0),
  
  // ************************ ****************************************************

  Bc_chi2(0), B_chi2(0), B_J_chi2(0), 
  B_Prob(0), B_J_Prob(0), 

  B_DecayVtxX(0),     B_DecayVtxY(0),     B_DecayVtxZ(0),
  B_DecayVtxXE(0),    B_DecayVtxYE(0),    B_DecayVtxZE(0),
  B_DecayVtxXYE(0),   B_DecayVtxXZE(0),   B_DecayVtxYZE(0),

  deltaRmum(0), deltaRmup(0), deltaRkp(0), deltaRkm(0), deltaRpion(0),
  istruemum(0), istruemup(0), istruekp(0), istruekm(0), istruebs(0), istruepion(0), istruebc(0),

  run(0), event(0),
  lumiblock(0)

{
   //now do what ever initialization is needed
}


Bspi::~Bspi()
{

}


// ------------ method called to for each event  ------------
void Bspi::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{  
  using std::vector;
  using namespace edm;
  using namespace reco;
  using namespace std;
  
  //*********************************
  // Get event content information
  //*********************************  

 // Kinematic fit
  edm::ESHandle<TransientTrackBuilder> theB; 
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB); 

  edm::Handle< View<pat::PackedCandidate> > thePATTrackHandle;
  iEvent.getByToken(trakCollection_label,thePATTrackHandle);

  edm::Handle< View<pat::PackedCandidate> > thePATTrackHandle_lostTracks;
  iEvent.getByToken(trakCollection_label_lostTracks,thePATTrackHandle_lostTracks); 

  edm::Handle< View<pat::Muon> > thePATMuonHandle; 
  iEvent.getByToken(dimuon_Label,thePATMuonHandle);

  edm::Handle<reco::GenParticleCollection> pruned;
  //edm::Handle<pat::PackedGenParticle> pruned; 
  iEvent.getByToken(genCands_, pruned);
  
  edm::Handle<pat::PackedGenParticleCollection> packed;
  iEvent.getByToken(packedGenToken_,packed);

  // ************* VERY IMPORTANT if you need to used Low pt tracks*****************************************************
  // add "lostTracks" to "packedPFCandidates"
  std::vector<const pat::PackedCandidate *> thePATTrackHandleT ;
  for (const pat::PackedCandidate &p1 : *thePATTrackHandle ) thePATTrackHandleT.push_back(&p1);
  for (const pat::PackedCandidate &p2 : *thePATTrackHandle_lostTracks ) thePATTrackHandleT.push_back(&p2);
  //cout << thePATTrackHandleT.size() << " = " << thePATTrackHandle->size() << " + " << thePATTrackHandle_lostTracks->size() << endl;
  // ***********************************************************************************************   

  
  //*********************************
  // Get gen level information
  //*********************************  
  /* 
  gen_bc_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_bc_vtx.SetXYZ(0.,0.,0.);
  gen_bs_p4.SetPtEtaPhiM(0.,0.,0.,0.); //new var
  gen_bs_vtx.SetXYZ(0.,0.,0.); // new var
  //gen_phi_p4.SetPtEtaPhiM(0.,0.,0.,0.); 
  gen_jpsi_vtx.SetXYZ(0.,0.,0.);
  gen_muon1_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_muon2_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_jpsi_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_kaon1_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_kaon2_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_pion_p4.SetPtEtaPhiM(0.,0.,0.,0.); //new var
  
  gen_b_ct = -9999.;
  */
  if ((isMC_ || OnlyGen_) && pruned.isValid() ) 
    {
      int foundit = 0;
      for (size_t i=0; i<pruned->size(); i++) 
	{
	  foundit = 0;
	  const reco::Candidate *dau = &(*pruned)[i];
	  
	  if ( (abs(dau->pdgId()) == 541) ) // Bc = 541 
	    {                              
	      foundit++;
	      //std::cout << "Found Bc" << std::endl;
	      //gen_bc_pt->push_back(dau->mass());
	      gen_bc_p4.SetPtEtaPhiM(dau->pt(),dau->eta(),dau->phi(),dau->mass());
	      gen_bc_vtx.SetXYZ(dau->vx(),dau->vy(),dau->vz());
	      //std::cout << "saving information" << std::endl;

	     
	      for ( size_t bs=0; bs<dau->numberOfDaughters(); bs++ ) 
		{
		  const reco::Candidate *gdau = dau->daughter(bs);
		  if ( gdau->pdgId() == 531 ) //&& (dau->status() == 2) ) 		Bs = 531 , came from Bc
		    { 
		      foundit++;
		      //std::cout << "Found Bs" << std::endl;
		      gen_bs_p4.SetPtEtaPhiM(gdau->pt(),gdau->eta(),gdau->phi(),gdau->mass());
		      gen_bs_vtx.SetXYZ(gdau->vx(),gdau->vy(),gdau->vz());
		      gen_b_ct = GetLifetime(gen_bs_p4,gen_bc_vtx,gen_bs_vtx);
		      // ------- start JPsi ------
		      
		      for (size_t k=0; k<gdau->numberOfDaughters(); k++) 
			{
			  const reco::Candidate *ggdau = gdau->daughter(k);

			  if (ggdau->pdgId() == 443 ) //&& gdau->status()==2)				JPsi = 443, came from Bs not from Bc
			    {

			      foundit++;

			      int nm=0;
			      for (size_t l=0; l<ggdau->numberOfDaughters(); l++) 
				{
				  const reco::Candidate *mm = ggdau->daughter(l);
				  // mu-
				  if (mm->pdgId()==13) 
				    { 
				      foundit++;     
				      //std::cout << "Found Muon" << std::endl;
				      if (mm->status()!=1) // What does "(mm->status()!=1)" mean?
					{				
					  for (size_t m=0; m<mm->numberOfDaughters(); m++) 
					    {
					      const reco::Candidate *mu = mm->daughter(m);
					      if (mu->pdgId()==13 ) //&& mu->status()==1) 
						{ 
						  nm++;
						  //std::cout<< "muon1 OK" << std::endl;
						  gen_muon1_p4.SetPtEtaPhiM(mu->pt(),mu->eta(),mu->phi(),mu->mass());
						  break;
						}
					    }
					}
				      else 
					{
					  //std::cout <<"muon1 OK else"  << std::endl;
					  gen_muon1_p4.SetPtEtaPhiM(mm->pt(),mm->eta(),mm->phi(),mm->mass());
					  nm++;
					}
				    }
				  // mu+
				  if (mm->pdgId()==-13) 
				    { 
				      foundit++;
				      if (mm->status()!=1) // ??????????????
					{			
					  for (size_t m=0; m<mm->numberOfDaughters(); m++) 
					    {
					      const reco::Candidate *mu = mm->daughter(m);
					      if (mu->pdgId()==-13 ) //&& mu->status()==1) 
						{ 
						  nm++;
						  //std::cout << "muon2 OK" << std::endl;
						  gen_muon2_p4.SetPtEtaPhiM(mu->pt(),mu->eta(),mu->phi(),mu->mass());
						  break;
						}
					    }
					} 
				      else 
					{
					  //std::cout << "muon2 OK else" << std::endl;
					  gen_muon2_p4.SetPtEtaPhiM(mm->pt(),mm->eta(),mm->phi(),mm->mass());
					  nm++;
					}
				    }
				}
			      
			      if (nm==2) 
				{
				  gen_jpsi_p4.SetPtEtaPhiM(ggdau->pt(),ggdau->eta(),ggdau->phi(),ggdau->mass());
				  gen_jpsi_vtx.SetXYZ(ggdau->vx(),ggdau->vy(),ggdau->vz());
				}
			      else 
				{
				  foundit-=nm;
				}
			      
			      //gen_jpsi_p4.SetPtEtaPhiM(ggdau->pt(),ggdau->eta(),ggdau->phi(),ggdau->mass());
			      // note from Jhovanny: Consider to add this line. 
			      //if (foundit == 5 ) break;  // Chic->Jpsi decay found!
			    } //end of: if (ggdau->pdgId()== 443 ) 
			  
			  // end of k -> checking
			  //--------  end of JPsi --------
			  
			  //-------- start of Phi -> K K -------
			  for (size_t lk=0; lk<packed->size(); lk++) 
			    {
			      const reco::Candidate * dauInPrunedColl = (*packed)[lk].mother(0);
			      int stable_id = (*packed)[lk].pdgId();
			      if (dauInPrunedColl != nullptr && isAncestor(gdau,dauInPrunedColl)) 
				{
				  if(stable_id == 321) 
				    {
				      foundit++;
				      //std::cout << "Found Kaon1" << std::endl;
				      gen_kaon1_p4.SetPtEtaPhiM((*packed)[lk].pt(),(*packed)[lk].eta(),(*packed)[lk].phi(),(*packed)[lk].mass());
				    }
				  if(stable_id == -321)
				    { 
				      foundit++;
				      //std::cout << "Found Kaon2" << std::endl;
				      gen_kaon2_p4.SetPtEtaPhiM((*packed)[lk].pt(),(*packed)[lk].eta(),(*packed)[lk].phi(),(*packed)[lk].mass());
				    }
				}
			    }
			  //-------- end of Phi -> K K --------
			} // for (size_t k
		    } // if (abs(dau->pdgId())==531 )
			// note from Jhovanny: Consider to add this line. 
			//if (foundit == 7 ) break;  // just one decay of this kind is expected
		  // for bs -> X
		  // This is working --->
		  // for pion
		  // ****************************************************************************
		  // Note from Jhovanny:
		  // I think the problem is here.
		  // For the pion you need to use the "packed" container, in the same way as with kaons 
		  // (use my "isAncestor" function and check if it comes from Bc)
		  // *****************************************************************************
		  for (size_t pi=0; pi<packed->size(); pi++) {
		  	const reco::Candidate * dauInPrunedCollp = (*packed)[pi].mother(0);	
		  	int stable_id_p = (*packed)[pi].pdgId();
		  		if (dauInPrunedCollp != nullptr && isAncestor(dau,dauInPrunedCollp)) {
					if(abs(stable_id_p == 211)) {foundit++;
					//std::cout << "Found Pion" << std::endl;
					gen_pion_p4.SetPtEtaPhiM((*packed)[pi].pt(),(*packed)[pi].eta(),(*packed)[pi].phi(),(*packed)[pi].mass());
					//break;// Jhovanny comment to add break here.
					}
				}
		  } // for pi
		//Done by me//-----------------------------------------------------------------------------------------------------------
		  /*for ( size_t pi=0; pi<dau->numberOfDaughters(); pi++ ) {
		  const reco::Candidate *gdau = dau->daughter(pi);
		  if(TMath::Abs(gdau->pdgId()) == 211 ) 
		    {
		      foundit++;
		      //std::cout << "Found pion" << std::endl;
		      gen_pion_p4.SetPtEtaPhiM(gdau->pt(),gdau->eta(), gdau->phi(), gdau->mass());	
		      //std::cout << "saving pion information" << std::endl;
		    } // for if(abs(dua->pdgId()) == 211 ) 
		  } // for pi */
		//------------------------------------------------------------------------------------------------------------
		} // for bs -> checking: --> checked and working now.
	    }// if (abs(dau->pdgId())==541)
	  if (foundit>=8) break;	 //1-Bc, 2-Bs, 3-JPsi, 4-mu1, 5-mu2, 6-kaon1, 7-kaon2, 8-pion
	} // for i
        /*      
	if (foundit!=8) {
	gen_bc_p4.SetPtEtaPhiM(0.,0.,0.,0.);
	gen_bc_vtx.SetXYZ(0.,0.,0.);
	gen_bs_p4.SetPtEtaPhiM(0.,0.,0.,0.);
	gen_bs_vtx.SetXYZ(0.,0.,0.);
	//gen_phi_p4.SetPtEtaPhiM(0.,0.,0.,0.);
	gen_jpsi_vtx.SetXYZ(0.,0.,0.);
	gen_muon1_p4.SetPtEtaPhiM(0.,0.,0.,0.);
	gen_muon2_p4.SetPtEtaPhiM(0.,0.,0.,0.);
	gen_jpsi_p4.SetPtEtaPhiM(0.,0.,0.,0.);
	gen_kaon1_p4.SetPtEtaPhiM(0.,0.,0.,0.);
	gen_kaon2_p4.SetPtEtaPhiM(0.,0.,0.,0.);
	gen_pion_p4.SetPtEtaPhiM(0.,0.,0.,0.);
	gen_b_ct = -9999.;
	//std::cout << "Does not found the given decay " << run << "," << event << " foundit=" << foundit << std::endl; // sanity check
	}
        */     
    }// end of isMC_ || OnlyGen_
  /*
  nB = 0; nMu = 0;
  if( isMC_) {
	tree_->Fill();
	std::cout << "GenLevel tree is filling" << std::endl;
	return;
  }
  */

  //*********************************
  //Now we get the primary vertex 
  //*********************************

  reco::Vertex bestVtx;
  reco::Vertex BcVtx;

  edm::Handle<std::vector<reco::Vertex> > recVtxs;
  iEvent.getByToken(primaryVertices_Label, recVtxs);
  edm::Handle<reco::VertexCollection> pvHandle_;
  iEvent.getByToken(primaryVertices_Label,pvHandle_);

  // Getting the first primary vertex of the container
  bestVtx = *(recVtxs->begin());
  
  priVtxX = bestVtx.x();
  priVtxY = bestVtx.y();
  priVtxZ = bestVtx.z();
  priVtxXE = bestVtx.covariance(0, 0);
  priVtxYE = bestVtx.covariance(1, 1);
  priVtxZE = bestVtx.covariance(2, 2);
  priVtxXYE = bestVtx.covariance(0, 1);
  priVtxXZE = bestVtx.covariance(0, 2);
  priVtxYZE = bestVtx.covariance(1, 2);
  priVtxCL = ChiSquaredProbability((double)(bestVtx.chi2()),(double)(bestVtx.ndof())); 
  
  nVtx = recVtxs->size();  
 
  lumiblock = iEvent.id().luminosityBlock();
  run = iEvent.id().run();
  event = iEvent.id().event();

//***************************************************************************************************
// Let's begin by looking for Jpsi->mu+mu-
//***************************************************************************************************

  unsigned int nMu_tmp = thePATMuonHandle->size();
  
  for(View<pat::Muon>::const_iterator iMuon1 = thePATMuonHandle->begin(); iMuon1 != thePATMuonHandle->end(); ++iMuon1) 
    {
      for(View<pat::Muon>::const_iterator iMuon2 = iMuon1+1; iMuon2 != thePATMuonHandle->end(); ++iMuon2) 
	{
	  if(iMuon1==iMuon2) continue;
	  
	  //opposite charge 
	  if( (iMuon1->charge())*(iMuon2->charge()) == 1) continue;

	  const pat::Muon *patMuonP = 0;
	  const pat::Muon *patMuonM = 0;
	  TrackRef glbTrackP;	  
	  TrackRef glbTrackM;	  
	  
	  if(iMuon1->charge() == 1){ patMuonP = &(*iMuon1); glbTrackP = iMuon1->track();}
	  if(iMuon1->charge() == -1){patMuonM = &(*iMuon1); glbTrackM = iMuon1->track();}
	  
	  if(iMuon2->charge() == 1) {patMuonP = &(*iMuon2); glbTrackP = iMuon2->track();}
	  if(iMuon2->charge() == -1){patMuonM = &(*iMuon2); glbTrackM = iMuon2->track();}
	  
	  if( glbTrackP.isNull() || glbTrackM.isNull() ) 
	    {
	      //std::cout << "continue due to no track ref" << endl;
	      continue;
	    }

	  if(iMuon1->track()->pt()<4.0) continue;
	  if(iMuon2->track()->pt()<4.0) continue;
	  //if(fabs(iMuon1->eta())>2.2 || fabs(iMuon2->eta())>2.2) continue;

	  if(!(glbTrackM->quality(reco::TrackBase::highPurity))) continue;
	  if(!(glbTrackP->quality(reco::TrackBase::highPurity))) continue;
        
	  //Let's check the vertex and mass
	  reco::TransientTrack muon1TT((*theB).build(glbTrackP));
	  reco::TransientTrack muon2TT((*theB).build(glbTrackM));

	  // *****  Trajectory states to calculate DCA for the 2 muons *********************
	  FreeTrajectoryState mu1State = muon1TT.impactPointTSCP().theState();
	  FreeTrajectoryState mu2State = muon2TT.impactPointTSCP().theState();

	  if( !muon1TT.impactPointTSCP().isValid() || !muon2TT.impactPointTSCP().isValid() ) continue;

	  // Measure distance between tracks at their closest approach
	  ClosestApproachInRPhi cApp;
	  cApp.calculate(mu1State, mu2State);
	  if( !cApp.status() ) continue;
	  float dca = fabs( cApp.distance() );	  
	  //if (dca < 0. || dca > 0.5) continue;
	  //cout<<" closest approach  "<<dca<<endl;

	  
	  //The mass of a muon and the insignificant mass sigma 
	  //to avoid singularities in the covariance matrix.
	  ParticleMass muon_mass = 0.10565837; //pdg mass
	  ParticleMass psi_mass = 3.096916;
	  float muon_sigma = muon_mass*1.e-6;
	  //float psi_sigma = psi_mass*1.e-6;

	  //Creating a KinematicParticleFactory
	  KinematicParticleFactoryFromTransientTrack pFactory;
	  VirtualKinematicParticleFactory vFactory;
	  
	  //initial chi2 and ndf before kinematic fits.
	  float chi = 0.;
	  float ndf = 0.;
	  vector<RefCountedKinematicParticle> muonParticles;
	  try {
	    muonParticles.push_back(pFactory.particle(muon1TT,muon_mass,chi,ndf,muon_sigma));
	    muonParticles.push_back(pFactory.particle(muon2TT,muon_mass,chi,ndf,muon_sigma));
	  }
	  catch(...) { 
	    std::cout<<" Exception caught ... continuing 1 "<<std::endl; 
	    continue;
	  }

	  KinematicParticleVertexFitter fitter;   

	  RefCountedKinematicTree psiVertexFitTree;
	  try {
	    psiVertexFitTree = fitter.fit(muonParticles); 
	  }
	  catch (...) { 
	    std::cout<<" Exception caught ... continuing 2 "<<std::endl; 
	    continue;
	  }

	  if (!psiVertexFitTree->isValid()) 
	    {
	      //std::cout << "caught an exception in the psi vertex fit" << std::endl;
	      continue; 
	    }

	  psiVertexFitTree->movePointerToTheTop();
	  
	  RefCountedKinematicParticle psi_vFit_noMC = psiVertexFitTree->currentParticle();
	  RefCountedKinematicVertex psi_vFit_vertex_noMC = psiVertexFitTree->currentDecayVertex();
	  
	  if( psi_vFit_vertex_noMC->chiSquared() < 0 )
	    {
	      //std::cout << "negative chisq from psi fit" << endl;
	      continue;
	    }

	   //some loose cuts go here

	   if(psi_vFit_vertex_noMC->chiSquared()>26.) continue;
	   if(psi_vFit_noMC->currentState().mass()<2.9 || psi_vFit_noMC->currentState().mass()>3.3) continue;

	   double J_Prob_tmp   = TMath::Prob(psi_vFit_vertex_noMC->chiSquared(),(int)psi_vFit_vertex_noMC->degreesOfFreedom());
	   if(J_Prob_tmp<0.01)
	     {
	       continue;
	     }

	   //Now that we have a J/psi candidate, we look for phi candidates

	   for(View<pat::PackedCandidate>::const_iterator iTrack1 = thePATTrackHandle->begin();
	   iTrack1 != thePATTrackHandle->end(); ++iTrack1 )
	     {
	       //quality cuts track1
               if(iTrack1->charge()==0) continue;
	       //if(fabs(iTrack1->pdgId())!=211) continue;
	       if(iTrack1->pt()<0.95) continue;
	       //if(iTrack1->pt()<0.5) continue;
	       if(!(iTrack1->trackHighPurity())) continue;
	       if(iTrack1->numberOfPixelHits()<1)continue;
	       if(iTrack1->numberOfHits()<5)continue;

	       for(View<pat::PackedCandidate>::const_iterator iTrack2 = iTrack1+1;
	       iTrack2 != thePATTrackHandle->end(); ++iTrack2 ) 
		 {
		   //quality cuts track2
		   if(iTrack1==iTrack2) continue;
		   if(iTrack2->charge()==0) continue;
		   //if(fabs(iTrack2->pdgId())!=211) continue;
		   if(iTrack2->pt()<0.95) continue;
		   //if(iTrack2->pt()<0.5) continue;
		   if(!(iTrack2->trackHighPurity())) continue;
		   if(iTrack2->numberOfPixelHits()<1)continue;
		   if(iTrack2->numberOfHits()<5)continue;

		   if(iTrack1->charge() == iTrack2->charge()) continue;

		   //Now let's checks if our muons do not use the same tracks as we are using now
		   if ( IsTheSame(*iTrack1,*iMuon1) || IsTheSame(*iTrack1,*iMuon2) ) continue;
		   if ( IsTheSame(*iTrack2,*iMuon1) || IsTheSame(*iTrack2,*iMuon2) ) continue;

		   //Now let's see if these two tracks make a vertex
		   reco::TransientTrack pion1TT((*theB).build(iTrack1->pseudoTrack()));
		   reco::TransientTrack pion2TT((*theB).build(iTrack2->pseudoTrack()));

		   //ParticleMass kaon_mass = 0.493677;
		   //float kaon_sigma = kaon_mass*1.e-6;
		   //ParticleMass pion_mass = 0.13957018;
		   ParticleMass pion_mass = Trkmass_;
		   float pion_sigma = pion_mass*1.e-6;

		   // ***************************
		   // pipi invariant mass (before kinematic vertex fit)
		   // ***************************
		   TLorentzVector pion14V,pion24V,pipi4V, Jpsi4V; 
		   pion14V.SetXYZM(iTrack1->px(),iTrack1->py(),iTrack1->pz(),pion_mass);
		   pion24V.SetXYZM(iTrack2->px(),iTrack2->py(),iTrack2->pz(),pion_mass);

		   pipi4V=pion14V+pion24V;
		   // in the case of the BstoJpsiphi It could be:
		   //if(pipi4V.M()<0.970 || pipi4V.M()>1.070) continue;
		   // in the case of the Jpsipipi It could be (0.450)
		   //if(pipi4V.M()<0.350) continue;
		   if(pipi4V.M()<TrkTrkMasscut_[0] || pipi4V.M()>TrkTrkMasscut_[1]) continue;

		   //initial chi2 and ndf before kinematic fits.
		   float chi = 0.;
		   float ndf = 0.;
	
		   // ************************************************
		   // Bs invariant mass (before kinematic vertex fit)
		   // ************************************************

		   Jpsi4V.SetXYZM(psi_vFit_noMC->currentState().globalMomentum().x(),psi_vFit_noMC->currentState().globalMomentum().y(),psi_vFit_noMC->currentState().globalMomentum().z(),psi_vFit_noMC->currentState().mass());

		   // in the case of the BstoJpsiphi It could be:
		   //if ( (pipi4V + Jpsi4V).M()<4.4 || (pipi4V + Jpsi4V).M()>6.3 ) continue;
		   // in the case of the Jpsipipi It could be
		   //if ( (pipi4V + Jpsi4V).M()<3.2 || (pipi4V + Jpsi4V).M()>4.4 ) continue;
		   if ( (pipi4V + Jpsi4V).M()<BarebMasscut_[0] || (pipi4V + Jpsi4V).M()>BarebMasscut_[1] ) continue;

		   vector<RefCountedKinematicParticle> vFitMCParticles;
		   vFitMCParticles.push_back(pFactory.particle(muon1TT,muon_mass,chi,ndf,muon_sigma));
		   vFitMCParticles.push_back(pFactory.particle(muon2TT,muon_mass,chi,ndf,muon_sigma));
		   vFitMCParticles.push_back(pFactory.particle(pion1TT,pion_mass,chi,ndf,pion_sigma));
		   vFitMCParticles.push_back(pFactory.particle(pion2TT,pion_mass,chi,ndf,pion_sigma));

                   // JPsi mass constraint is applied in the final Bs fit,                                                               
                   MultiTrackKinematicConstraint *  j_psi_c = new  TwoTrackMassKinematicConstraint(psi_mass);
                   KinematicConstrainedVertexFitter kcvFitter;
                   RefCountedKinematicTree vertexFitTree = kcvFitter.fit(vFitMCParticles, j_psi_c);
                   if (!vertexFitTree->isValid()) {
                     //std::cout << "caught an exception in the B vertex fit with MC" << std::endl;                                        
                     continue;
                   }

		   vertexFitTree->movePointerToTheTop();
		   RefCountedKinematicParticle bCandMC = vertexFitTree->currentParticle();
		   RefCountedKinematicVertex bDecayVertexMC = vertexFitTree->currentDecayVertex();

		   reco::TransientTrack BsTrack = bCandMC->refittedTransientTrack();

		   if (!bDecayVertexMC->vertexIsValid()){
		     //std::cout << "B MC fit vertex is not valid" << endl;
		     continue;
		   }

		   // in the case of the BstoJpsiphi It could be:
		   //if(bCandMC->currentState().mass()<5.0 || bCandMC->currentState().mass()>6.0) continue;
		   // in the case of the Jpsipipi It could be
		   //if(bCandMC->currentState().mass()<3.3 || bCandMC->currentState().mass()>4.0) continue;
		   if(bCandMC->currentState().mass()<bMasscut_[0] || bCandMC->currentState().mass()>bMasscut_[1]) continue;
 
		   if(bDecayVertexMC->chiSquared()<0 || bDecayVertexMC->chiSquared()>50 ) 
		     {
		       //std::cout << " continue from negative chi2 = " << bDecayVertexMC->chiSquared() << endl;
		       continue;
		     }		   

		   double B_Prob_tmp  = TMath::Prob(bDecayVertexMC->chiSquared(),(int)bDecayVertexMC->degreesOfFreedom());
		   if(B_Prob_tmp<0.01)
		     {
		       continue;
		     }
		    /*-----------------------------------------------------------------------
		    // Chosing the closest PV in Z direction to the Bs trajectory projection
		    double dzMin = 1000000.;
		    const reco::VertexCollection* vertices = recVtxs.product();
		    for(reco::VertexCollection::const_iterator  primVertex = vertices->begin(); primVertex!= vertices->end(); primVertex++) 
		    {
			// std::cout "prim vertex z: " << primVertex->z() << std::endl;
			if (abs(dzMin) > abs(BsTrack.track().dz(primVertex->position())))
        		{
          			bestVtx = *(primVertex);
          			dzMin = BsTrack.track().dz(primVertex->position());
        		}
      		    } 

		    reco::Vertex bestVertexBSC = getPVConstrainedToBS(iEvent, iSetup, bestVtx);
		    TVector3 primaryVertex(bestVtx.x(),bestVtx.y(),0.);
		    TVector3 bDecayVertexPosition(bDecayVertexMC->position().x(),
                          bDecayVertexMC->position().y(),
                          0.);

		    double Lxy_tmp = bDecayVertexPosition.Mag() - primaryVertex.Dot(bDecayVertexPosition) / bDecayVertexPosition.Mag();
		    ------------------------------------------------------------------------------*/
		   // get children from final B fit

		   vertexFitTree->movePointerToTheFirstChild();
                   RefCountedKinematicParticle mu1CandMC = vertexFitTree->currentParticle();

                   vertexFitTree->movePointerToTheNextChild();
                   RefCountedKinematicParticle mu2CandMC = vertexFitTree->currentParticle();

                   vertexFitTree->movePointerToTheNextChild();
                   RefCountedKinematicParticle T1CandMC = vertexFitTree->currentParticle();

                   vertexFitTree->movePointerToTheNextChild();
                   RefCountedKinematicParticle T2CandMC = vertexFitTree->currentParticle();

		   KinematicParameters psiMu1KP = mu1CandMC->currentState().kinematicParameters();
		   KinematicParameters psiMu2KP = mu2CandMC->currentState().kinematicParameters();		   

		   KinematicParameters phiPi1KP = T1CandMC->currentState().kinematicParameters();
		   KinematicParameters phiPi2KP = T2CandMC->currentState().kinematicParameters();

		   const reco::Muon *recoMuonM = patMuonM;
		   const reco::Muon *recoMuonP = patMuonP;

		   TVector3 primaryVertex(bestVtx.x(),bestVtx.y(),0.);
		
		   // Loking for Pion ---------------------------------------------------------

		   for (unsigned int i = 0, n = thePATTrackHandleT.size(); i<n ; ++i )
		     {
		       const pat::PackedCandidate *iTrack3 =  (thePATTrackHandleT)[i];
		       //quality cuts track3
		       if(iTrack3->charge()==0) continue;
		       if(fabs(iTrack3->pdgId())!=211) continue;
		       if(iTrack3->pt()<0.5) continue;
		       //if(iTrack3->pt()<0.2) continue;
		       if(!(iTrack3->trackHighPurity())) continue;
		       if(iTrack3->numberOfPixelHits()<1)continue;
		       if(iTrack3->numberOfHits()<5)continue;
		       		       
		       //Now let's checks if our muons do not use the same tracks as we are using now
		       if ( IsTheSame(*iTrack3,*iMuon1) || IsTheSame(*iTrack3,*iMuon2) ) continue;

		       //***** Let's recontruct the pion track
		       reco::Track globalTrack3;
                       globalTrack3 = iTrack3->pseudoTrack();
		       
		       reco::TransientTrack pionTT((*theB).build(iTrack3->pseudoTrack()));
		       if(!pionTT.isValid()) continue;

		       // ***************************
		       // Bc invariant mass ()
		       // ***************************
		       
		       ParticleMass pion_mass = 0.13957018;
		       float pion_sigma = pion_mass*1.e-6;
		       ParticleMass bs_mass = 5.36679;
		       
		       TLorentzVector Bs4V,pion4V,Bc4V; 
		       
		       pion4V.SetXYZM(iTrack3->px(),iTrack3->py(),iTrack3->pz(),pion_mass);
		       Bs4V.SetXYZM(bCandMC->currentState().globalMomentum().x(),bCandMC->currentState().globalMomentum().y(),bCandMC->currentState().globalMomentum().z(),bCandMC->currentState().mass());
		       
		       ParticleMass PDG_BC_MASS = 6.2751; 
		       Bc4V=pion4V+Bs4V;
		       //cout<<"mass Bc: "<<Bc4V.M()<<endl;
		       //if(fabs(Bc4V.M() - PDG_BC_MASS > 0.6)) continue;
		       //if ((Bc4V.M() - Bs4V.M() + bs_mass < 5.0) || (Bc4V.M() - Bs4V.M() + bs_mass > 7.5)) continue;
		       //if ((Bc4V.M() - Bs4V.M() + bs_mass < 6.0) || (Bc4V.M() - Bs4V.M() + bs_mass > 6.6)) continue;
		       //if ((Bc4V.M() - Bs4V.M() + bs_mass < 5.0) || (Bc4V.M() - Bs4V.M() + bs_mass > 6.6)) continue;
		       //if ((Bc4V.M() - Bs4V.M() + bs_mass < 5.8) || (Bc4V.M() - Bs4V.M() + bs_mass > 6.6)) continue;
		       if ((Bc4V.M() < 5.8) || (Bc4V.M() > 6.6)) continue;
		       //if ((Bc4V.M() < 6.0) || (Bc4V.M() > 6.6)) continue;
		       //TLorentzVector bc4V;


		       // using VirtualKinematicParticleFactory vFactory for Bs
		       float Bs_dof  = bDecayVertexMC->degreesOfFreedom();
		       float Bs_chi2 = bDecayVertexMC->chiSquared();

		       vector<RefCountedKinematicParticle> vFitMCParticlesBspi;
		       vFitMCParticlesBspi.push_back(vFactory.particle(bCandMC->currentState(),Bs_chi2,Bs_dof,bCandMC));
		       vFitMCParticlesBspi.push_back(pFactory.particle(pionTT,pion_mass,chi,ndf,pion_sigma));

		       KinematicParticleVertexFitter kcvFitterBspi;
		       RefCountedKinematicTree vertexFitTreeBspi = kcvFitterBspi.fit(vFitMCParticlesBspi);
		       if (!vertexFitTreeBspi->isValid()) {
			 //std::cout << "caught an exception in the Bc vertex" << std::endl;
			 continue;
		       }
			
		       vertexFitTreeBspi->movePointerToTheTop();
		       RefCountedKinematicParticle bspiCandMC = vertexFitTreeBspi->currentParticle();
		       RefCountedKinematicVertex bspiDecayVertexMC = vertexFitTreeBspi->currentDecayVertex();
                      
		       TVector3 bcDecayVertexPosition(0.,0.,0.);
		       bcDecayVertexPosition.SetXYZ(bspiDecayVertexMC->position().x(),bspiDecayVertexMC->position().y(),bspiDecayVertexMC->position().z());
			/*
		       bc4V.SetPtEtaPhiM(bspiCandMC->currentState().globalMomentum().perp(),
              		   bspiCandMC->currentState().globalMomentum().eta(),
              		   bspiCandMC->currentState().globalMomentum().phi(),
              		   bspiCandMC->currentState().mass());
			*/
		       if (!bspiDecayVertexMC->vertexIsValid()){
			 //std::cout << "B MC fit vertex is not valid" << endl;
			 continue;
		       }

		       double Bc_mass_cjp_tmp = bspiCandMC->currentState().mass();
		       //if(fabs(Bc_mass_cjp_tmp - PDG_BC_MASS) > 0.2) continue;
		       //if(bspiCandMC->currentState().mass()<6.1 || bspiCandMC->currentState().mass()>6.7) continue;
		       //if(bspiCandMC->currentState().mass()<6.1 || bspiCandMC->currentState().mass()>6.5) continue;
		       //if(bspiCandMC->currentState().mass()<6.1 || bspiCandMC->currentState().mass()>6.3) continue;
		       //if(bspiCandMC->currentState().mass()<6.13 || bspiCandMC->currentState().mass()>6.25) continue;
		       //if(bspiCandMC->currentState().mass()<6.15 || bspiCandMC->currentState().mass()>6.22) continue;
		       //if(bspiCandMC->currentState().mass()<5.8 || bspiCandMC->currentState().mass()>6.6) continue;
		       //if(bspiCandMC->currentState().mass()<6.19 || bspiCandMC->currentState().mass()>6.22) continue; // Good Bc mass peak but problem with pion matching
		       if(bspiCandMC->currentState().mass()<6.26|| bspiCandMC->currentState().mass()>6.28) continue; // Good Bc and pion is also matching- This mass window applied for MC and data- using this I also working with BDT.
		       //if(bspiCandMC->currentState().mass()<5.4 || bspiCandMC->currentState().mass()>7.2) continue;   // seems Good Bc mass but pion is also working-in this case also
		       //if(bspiCandMC->currentState().mass()<5.0 || bspiCandMC->currentState().mass()>7.0) continue;   // seems Good Bc mass but pion is also working-in this case also
		       //if(bspiCandMC->currentState().mass()<6.15 || bspiCandMC->currentState().mass()>6.4) continue; 
		       //if(bspiDecayVertexMC->chiSquared()<0) continue;

		       double Bc_Prob_tmp  = TMath::Prob(bspiDecayVertexMC->chiSquared(),(int)bspiDecayVertexMC->degreesOfFreedom());
		       //if(Bc_Prob_tmp<0.01)
		       //{
		       	  //continue;
		       //}

		       ROOT::Math::XYZPoint BcDecayPoint( (*bspiDecayVertexMC).position().x(), (*bspiDecayVertexMC).position().y(), (*bspiDecayVertexMC).position().z() );	

		       //std::cout << "bspi chi2 " << bspiDecayVertexMC->chiSquared() << std::endl;
		       //std::cout << "bspi mass " << bspiCandMC->currentState().mass() << std::endl;
			/* =========
		       //get children from final Bspi fit (I am not absolutly sure about if it is necesary)		   
		       vertexFitTreeBspi->movePointerToTheFirstChild();
		       RefCountedKinematicParticle BsCandfrombspi = vertexFitTreeBspi->currentParticle();
		       
		       vertexFitTreeBspi->movePointerToTheNextChild();
		       RefCountedKinematicParticle pionCand = vertexFitTreeBspi->currentParticle();

		       KinematicParameters PionTrk = pionCand->currentState().kinematicParameters();
		
		       TVector3 reco_trk_p3;
        	       reco_trk_p3.SetXYZ(globalTrack3.px(),globalTrack3.py(),globalTrack3.pz());
			============= */
		       double Bc_ct_tmp = GetLifetime(Bc4V, primaryVertex, bcDecayVertexPosition);


		       
		       // fill candidate variables now		       
		       if(nB==0){
			 nMu  = nMu_tmp;
			 // cout<< "*Number of Muons : " << nMu_tmp << endl;
		       } // end nB==0

		       //Bspi_mass_vertex->push_back( bspiCandMC->currentState().mass() );
		       //=========== Information of Bc meson
		       deltaMass->push_back( Bc4V.M() - Bs4V.M() ); //+ bs_mass );
		       //Bspion_mass->push_back( Bc4V.M() );
		       Bc_mass->push_back(Bc4V.M());
		       //Bc_mass->push_back(Bc_mass_cjp_tmp);
		       //Bc_mass->push_back(bspiCandMC->currentState().mass());
                       Bc_pt->push_back(bspiCandMC->currentState().globalMomentum().perp());
                       //Bc_px->push_back(Bc4V.Px());
                       //Bc_py->push_back(Bc4V.Py());
                       //Bc_pz->push_back(Bc4V.Pz());
                       Bc_px->push_back(bspiCandMC->currentState().globalMomentum().x());
                       Bc_py->push_back(bspiCandMC->currentState().globalMomentum().y());
                       Bc_pz->push_back(bspiCandMC->currentState().globalMomentum().z());
                       Bc_eta->push_back(bspiCandMC->currentState().globalMomentum().eta());
                       Bc_phi->push_back(bspiCandMC->currentState().globalMomentum().phi());
                       Bc_ct->push_back(Bc_ct_tmp);
                       Bc_charge->push_back(bspiCandMC->currentState().particleCharge());
		       //Bcvtxcl->push_back(ChiSquaredProbability((double)(BcVtx.chi2()),(double)(BcVtx.ndof())));
		       Bcvtxcl->push_back(ChiSquaredProbability((double)(bspiDecayVertexMC->chiSquared()),int(rint(bspiDecayVertexMC->degreesOfFreedom()))));
		       Bcvtxcl->push_back(priVtxCL);
		       Bc_prob->push_back(Bc_Prob_tmp);
		       Bc_chi2->push_back(bspiDecayVertexMC->chiSquared());

		       Bc_DecayVtxX->push_back(bspiDecayVertexMC->position().x());
                       Bc_DecayVtxY->push_back(bspiDecayVertexMC->position().y());
                       Bc_DecayVtxZ->push_back(bspiDecayVertexMC->position().z());
                       Bc_DecayVtxXE->push_back(bspiDecayVertexMC->error().cxx());
                       Bc_DecayVtxYE->push_back(bspiDecayVertexMC->error().cyy());
                       Bc_DecayVtxZE->push_back(bspiDecayVertexMC->error().czz());

                       //Bc_DecayVtx_vtxfit_X ->push_back(    BcVtx.x() );
                       //Bc_DecayVtx_vtxfit_Y ->push_back(    BcVtx.y() );
                       //Bc_DecayVtx_vtxfit_Z ->push_back(    BcVtx.z() );
                       //Bc_DecayVtx_vtxfit_XE->push_back(    BcVtx.covariance(0, 0) );
                       //Bc_DecayVtx_vtxfit_YE->push_back(    BcVtx.covariance(1, 1) );
                       //Bc_DecayVtx_vtxfit_ZE->push_back(    BcVtx.covariance(2, 2) );
                       //Bc_DecayVtx_vtxfit_XYE->push_back(   BcVtx.covariance(0, 1) );
                       //Bc_DecayVtx_vtxfit_XZE->push_back(   BcVtx.covariance(0, 2) );
                       //Bc_DecayVtx_vtxfit_YZE->push_back(   BcVtx.covariance(1, 2) );
		       
		     
    		       //========= Information of Bs meson 
		       B_mass->push_back(bCandMC->currentState().mass());
		       B_px->push_back(bCandMC->currentState().globalMomentum().x());
		       B_py->push_back(bCandMC->currentState().globalMomentum().y());
		       B_pz->push_back(bCandMC->currentState().globalMomentum().z());
		       B_pt->push_back(bCandMC->currentState().globalMomentum().perp());
		       B_eta->push_back(bCandMC->currentState().globalMomentum().eta());
		       B_phi->push_back(bCandMC->currentState().globalMomentum().phi());
		       
		       //========= Information of Phi meson  
		       B_phi_mass->push_back( pipi4V.M() );
		       B_phi_px->push_back( pipi4V.Px() );
		       B_phi_py->push_back( pipi4V.Py() );
		       B_phi_pz->push_back( pipi4V.Pz() );
		       B_phi_pt->push_back( pipi4V.Pt() );
		       B_phi_eta->push_back( pipi4V.Eta() );
		       B_phi_phi->push_back( pipi4V.Phi() );
		       
		       // =========== Information of JPsi meson ===========
		       B_J_mass->push_back( psi_vFit_noMC->currentState().mass() );
		       B_J_px->push_back( psi_vFit_noMC->currentState().globalMomentum().x() );
		       B_J_py->push_back( psi_vFit_noMC->currentState().globalMomentum().y() );
		       B_J_pz->push_back( psi_vFit_noMC->currentState().globalMomentum().z() );
		       B_J_pt->push_back( psi_vFit_noMC->currentState().globalMomentum().perp() );
		       B_J_eta->push_back( psi_vFit_noMC->currentState().globalMomentum().eta() );
		       B_J_phi->push_back( psi_vFit_noMC->currentState().globalMomentum().phi() );

		       //========= Getting Pion Information ================

		       //B_pion_px->push_back(PionTrk.momentum().x());
		       //B_pion_py->push_back(PionTrk.momentum().y());
		       //B_pion_pz->push_back(PionTrk.momentum().z());
		       //B_pion_px->push_back(pion4V.Px());
		       //B_pion_py->push_back(pion4V.Py());
		       //B_pion_pz->push_back(pion4V.Px());

		       B_pion_px_track->push_back(iTrack3->px());
		       B_pion_py_track->push_back(iTrack3->py());
		       B_pion_pz_track->push_back(iTrack3->pz());
		       B_pion_pt_track->push_back(iTrack3->pt());
		       //B_pion_pt_track->push_back(pion4V.Pt());
		       //B_pion_charge->push_back(pionCand->currentState().particleCharge());
		       B_pion_charge->push_back(iTrack3->charge());
		       pion_track_normchi2  ->push_back(globalTrack3.normalizedChi2());
		       pion_NumHits->push_back(globalTrack3.numberOfValidHits());
                       pion_NumPixelHits->push_back(globalTrack3.hitPattern().numberOfValidPixelHits()); 
                       pion_dxy->push_back(globalTrack3.dxy(BcDecayPoint)); 
                       pion_dxy_->push_back(iTrack3->dxy()); 
                       pion_dz->push_back(globalTrack3.dz(BcDecayPoint));
                       pion_dz_->push_back(iTrack3->dz());
		       pion_NTrackerLayers->push_back ( globalTrack3.hitPattern().trackerLayersWithMeasurement() );
                       pion_NPixelLayers->push_back ( globalTrack3.hitPattern().pixelLayersWithMeasurement() );


		       // You can get the momentum components (for muons and kaon) from the final B childrens or of the original Tracks. Here, a example for the kaons:
		       B_phi_px1->push_back(phiPi1KP.momentum().x());
		       B_phi_py1->push_back(phiPi1KP.momentum().y());
		       B_phi_pz1->push_back(phiPi1KP.momentum().z());
		       B_phi_pz1->push_back(phiPi1KP.momentum().z());
		       B_phi_pt1->push_back(phiPi1KP.momentum().perp());
		       B_phi_eta1->push_back(phiPi1KP.momentum().eta());
		       B_phi_phi1->push_back(phiPi1KP.momentum().phi());
		       B_phi_px1_track->push_back(iTrack1->px());
		       B_phi_py1_track->push_back(iTrack1->py());
		       B_phi_pz1_track->push_back(iTrack1->pz());
		       B_phi_charge1->push_back(T1CandMC->currentState().particleCharge());
		       
		       B_phi_px2->push_back(phiPi2KP.momentum().x());
		       B_phi_py2->push_back(phiPi2KP.momentum().y());
		       B_phi_pz2->push_back(phiPi2KP.momentum().z());
		       B_phi_pt2->push_back(phiPi2KP.momentum().perp());
		       B_phi_eta2->push_back(phiPi2KP.momentum().eta());
		       B_phi_phi2->push_back(phiPi2KP.momentum().phi());
		       B_phi_px2_track->push_back(iTrack2->px());
		       B_phi_py2_track->push_back(iTrack2->py());
		       B_phi_pz2_track->push_back(iTrack2->pz());
		       B_phi_charge2->push_back(T2CandMC->currentState().particleCharge());
		       
		       B_J_px1->push_back(psiMu1KP.momentum().x());
		       B_J_py1->push_back(psiMu1KP.momentum().y());
		       B_J_pz1->push_back(psiMu1KP.momentum().z());
		       B_J_pt1->push_back(psiMu1KP.momentum().perp());
		       B_J_eta1->push_back(psiMu1KP.momentum().eta());
		       B_J_phi1->push_back(psiMu1KP.momentum().phi());
		       B_J_charge1->push_back(mu1CandMC->currentState().particleCharge());
		       
		       B_J_px2->push_back(psiMu2KP.momentum().x());
		       B_J_py2->push_back(psiMu2KP.momentum().y());
		       B_J_pz2->push_back(psiMu2KP.momentum().z());
		       B_J_pt2->push_back(psiMu2KP.momentum().perp());
		       B_J_eta2->push_back(psiMu1KP.momentum().eta());
		       B_J_phi2->push_back(psiMu1KP.momentum().phi());
		       B_J_charge2->push_back(mu2CandMC->currentState().particleCharge());
		       
		       B_J_chi2->push_back(psi_vFit_vertex_noMC->chiSquared());
		       B_chi2->push_back(bDecayVertexMC->chiSquared());
		       
		       //double B_Prob_tmp       = TMath::Prob(bDecayVertexMC->chiSquared(),(int)bDecayVertexMC->degreesOfFreedom());
		       //double J_Prob_tmp   = TMath::Prob(psi_vFit_vertex_noMC->chiSquared(),(int)psi_vFit_vertex_noMC->degreesOfFreedom());
		       B_Prob    ->push_back(B_Prob_tmp);
		       B_J_Prob  ->push_back(J_Prob_tmp);
		       
		       B_DecayVtxX ->push_back((*bDecayVertexMC).position().x());    
		       B_DecayVtxY ->push_back((*bDecayVertexMC).position().y());
		       B_DecayVtxZ ->push_back((*bDecayVertexMC).position().z());
		       B_DecayVtxXE ->push_back(bDecayVertexMC->error().cxx());   
		       B_DecayVtxYE ->push_back(bDecayVertexMC->error().cyy());   
		       B_DecayVtxZE ->push_back(bDecayVertexMC->error().czz());
		       B_DecayVtxXYE ->push_back(bDecayVertexMC->error().cyx());
		       B_DecayVtxXZE ->push_back(bDecayVertexMC->error().czx());
		       B_DecayVtxYZE ->push_back(bDecayVertexMC->error().czy());
		       
		       // ********************* muon-trigger-machint ****************
		       
		       const pat::Muon* muon1 = &(*iMuon1);
		       const pat::Muon* muon2 = &(*iMuon2);
		       
		       int tri_Dim25_tmp = 0, tri_JpsiTk_tmp = 0,  tri_JpsiTkTk_tmp = 0;
		       
		       if (muon1->triggerObjectMatchByPath("HLT_Dimuon25_Jpsi_v*")!=nullptr && muon2->triggerObjectMatchByPath("HLT_Dimuon25_Jpsi_v*")!=nullptr) tri_Dim25_tmp = 1;
		       if (muon1->triggerObjectMatchByPath("HLT_DoubleMu4_JpsiTrk_Displaced_v*")!=nullptr && muon2->triggerObjectMatchByPath("HLT_DoubleMu4_JpsiTrk_Displaced_v*")!=nullptr) tri_JpsiTk_tmp = 1;
		       if (muon1->triggerObjectMatchByPath("HLT_DoubleMu4_JpsiTrkTrk_Displaced_v*")!=nullptr && muon2->triggerObjectMatchByPath("HLT_DoubleMu4_JpsiTrkTrk_Displaced_v*")!=nullptr) tri_JpsiTkTk_tmp = 1;
		       
		       tri_Dim25->push_back( tri_Dim25_tmp );	       
		       tri_JpsiTk->push_back( tri_JpsiTk_tmp );
		       tri_JpsiTkTk->push_back( tri_JpsiTkTk_tmp );
		       
		       // ************ Different muons Id, and other properties  ****************
		       
		       mu1soft->push_back(iMuon1->isSoftMuon(bestVtx) );
		       mu2soft->push_back(iMuon2->isSoftMuon(bestVtx) );
		       mu1tight->push_back(iMuon1->isTightMuon(bestVtx) );
		       mu2tight->push_back(iMuon2->isTightMuon(bestVtx) );
		       mu1PF->push_back(iMuon1->isPFMuon());
		       mu2PF->push_back(iMuon2->isPFMuon());
		       mu1loose->push_back(muon::isLooseMuon(*iMuon1));
		       mu2loose->push_back(muon::isLooseMuon(*iMuon2));
		       
		       mumC2->push_back( glbTrackM->normalizedChi2() );
		       mumAngT->push_back( muon::isGoodMuon(*recoMuonM,muon::TMOneStationTight) ); 
		       mumNHits->push_back( glbTrackM->numberOfValidHits() );
		       mumNPHits->push_back( glbTrackM->hitPattern().numberOfValidPixelHits() );	       
		       mupC2->push_back( glbTrackP->normalizedChi2() );
		       mupAngT->push_back( muon::isGoodMuon(*recoMuonP,muon::TMOneStationTight) ); 
		       mupNHits->push_back( glbTrackP->numberOfValidHits() );
		       mupNPHits->push_back( glbTrackP->hitPattern().numberOfValidPixelHits() );
		       mumdxy->push_back(glbTrackM->dxy(bestVtx.position()) );
		       mupdxy->push_back(glbTrackP->dxy(bestVtx.position()) );
		       mumdz->push_back(glbTrackM->dz(bestVtx.position()) );
		       mupdz->push_back(glbTrackP->dz(bestVtx.position()) );
		       muon_dca->push_back(dca);
		       
		       k1dxy->push_back(iTrack1->dxy());
		       k2dxy->push_back(iTrack2->dxy());
		       k1dz->push_back(iTrack1->dz());
		       k2dz->push_back(iTrack2->dz());
		       
		       k1dxy_e->push_back(iTrack1->dxyError());
		       k2dxy_e->push_back(iTrack2->dxyError());
		       k1dz_e->push_back(iTrack1->dzError());
		       k2dz_e->push_back(iTrack2->dzError());
		       
		       k1InnerHits->push_back(iTrack1->lostInnerHits());
		       k2InnerHits->push_back(iTrack2->lostInnerHits());		   
		       
		       nB++;	       

		       int pvIndex = -1;
		       saveIP(vertexFitTree,*pvHandle_.product(), pvIndex);
		       //std::cout<<" pvindex "<<pvIndex<<std::endl;

		       if(isMC_)saveTruthMatch(iEvent);	       
		       
		       /////////////////////////////////////////////////
		       
		       //pionParticles.clear();
		       muonParticles.clear();
		       vFitMCParticles.clear();
		       vFitMCParticlesBspi.clear();
		       
		     }
		 }
	     }
	}
    }
  
   //fill the tree and clear the vectors
   if (nB > 0 ) 
     {
       //std::cout << "filling tree" << endl;
       tree_->Fill();
     }

   // *********

   nB = 0; nMu = 0;

   //triggersL = 0; 
   //Bspi_mass_vertex->clear();
   //Bspion_mass->clear();
   Bc_mass->clear();  Bc_pt->clear(); Bc_px->clear(); Bc_py->clear(); Bc_pz->clear(); Bc_eta->clear(); Bc_phi->clear(); Bc_ct->clear(); Bc_charge->clear();
   Bcvtxcl->clear(); Bc_prob->clear();

   deltaMass->clear();
   T_Bspion_mass->clear(); Notmatch_Bspion_mass->clear();

   Bc_DecayVtxX->clear();    Bc_DecayVtxY->clear();   Bc_DecayVtxZ->clear();
   Bc_DecayVtxXE->clear();   Bc_DecayVtxYE->clear();  Bc_DecayVtxZE->clear();

   //Bc_DecayVtx_vtxfit_X->clear();   Bc_DecayVtx_vtxfit_Y->clear();  Bc_DecayVtx_vtxfit_Z->clear();
   //Bc_DecayVtx_vtxfit_XE->clear();   Bc_DecayVtx_vtxfit_YE->clear();  Bc_DecayVtx_vtxfit_ZE->clear();
   //Bc_DecayVtx_vtxfit_XYE->clear();   Bc_DecayVtx_vtxfit_XZE->clear();  Bc_DecayVtx_vtxfit_YZE->clear();

   B_mass->clear();    B_px->clear();    B_py->clear();    B_pz->clear();
   B_pt->clear(); B_eta->clear(); B_phi->clear();	

   B_phi_mass->clear(); 
   B_phi_px->clear(); B_phi_py->clear(); B_phi_pz->clear();
   B_phi_pt->clear(); B_phi_eta->clear(); B_phi_phi->clear();

   B_J_mass->clear();  B_J_px->clear();  B_J_py->clear();  B_J_pz->clear();
   B_J_pt->clear();  B_J_eta->clear();  B_J_phi->clear();
   
   //B_pion_px->clear(); B_pion_py->clear(); B_pion_pz->clear();
   pion_track_normchi2->clear(); pion_NumHits->clear(); pion_NumPixelHits->clear();
   pion_dxy->clear(); pion_dz->clear();
   pion_dxy_->clear(); pion_dz_->clear();
   pion_NTrackerLayers->clear(); pion_NPixelLayers->clear();
   B_pion_px_track->clear(); B_pion_py_track->clear(); B_pion_pz_track->clear();
   B_pion_pt_track->clear(); 
   B_pion_charge->clear();

   B_phi_px1->clear(); B_phi_py1->clear(); B_phi_pz1->clear(); B_phi_pt1->clear(); B_phi_charge1->clear(); 
   B_phi_eta1->clear(); B_phi_phi1->clear(); 
   B_phi_px2->clear(); B_phi_py2->clear(); B_phi_pz2->clear(); B_phi_pt2->clear(); B_phi_charge2->clear(); 
   B_phi_eta2->clear(); B_phi_phi2->clear(); 

   B_phi_px1_track->clear(); B_phi_py1_track->clear(); B_phi_pz1_track->clear();
   B_phi_px2_track->clear(); B_phi_py2_track->clear(); B_phi_pz2_track->clear();

   B_J_px1->clear();  B_J_py1->clear();  B_J_pz1->clear(), B_J_charge1->clear();
   B_J_px2->clear();  B_J_py2->clear();  B_J_pz2->clear(), B_J_charge2->clear();
   B_J_pt1->clear();
   B_J_eta1->clear(); B_J_phi1->clear();
   B_J_pt2->clear();
   B_J_eta2->clear(); B_J_phi2->clear();

   Bc_chi2->clear(); B_chi2->clear(); B_J_chi2->clear();
   B_Prob->clear(); B_J_Prob->clear(); 

   B_DecayVtxX->clear();     B_DecayVtxY->clear();     B_DecayVtxZ->clear();
   B_DecayVtxXE->clear();    B_DecayVtxYE->clear();    B_DecayVtxZE->clear();
   B_DecayVtxXYE->clear();   B_DecayVtxXZE->clear();   B_DecayVtxYZE->clear();

   nVtx = 0; 
   priVtxX = 0;     priVtxY = 0;     priVtxZ = 0; 
   priVtxXE = 0;    priVtxYE = 0;    priVtxZE = 0; priVtxCL = 0;
   priVtxXYE = 0;   priVtxXZE = 0;   priVtxYZE = 0; 

   //******************************************************************
   
   pVtxIPX->clear();  pVtxIPY->clear();  pVtxIPZ->clear();
   pVtxIPXE->clear();  pVtxIPYE->clear();  pVtxIPZE->clear();  pVtxIPCL->clear();
   pVtxIPXYE->clear();  pVtxIPXZE->clear();  pVtxIPYZE->clear();
   
   B_l3d->clear();  B_l3dE->clear();  /*B_lxy->clear();*/ Bc_lxy->clear(); B_lxyE->clear();
   B_cosalpha->clear();   /*B_cosalphaxy->clear();*/Bc_cosalphaxy->clear(); alpha->clear();  //B_treco->clear();   B_trecoe->clear();  B_trecoxy->clear(); B_trecoxye->clear();
   B_pvip->clear(); B_pviperr->clear(); B_pvips->clear(); B_pvlzip->clear(); B_pvlziperr->clear(); B_pvlzips->clear();
   B_pv2ip->clear(); B_pv2iperr->clear(); B_pv2ips->clear(); B_pv2lzip->clear(); B_pv2lziperr->clear(); B_pv2lzips->clear();
   
   B_l3d_pv2->clear();  B_l3dE_pv2->clear();
   
   //******************************************************************

   k1dxy->clear(); k2dxy->clear(); k1dz->clear(); k2dz->clear();
   k1dxy_e->clear(); k2dxy_e->clear(); k1dz_e->clear(); k2dz_e->clear();
   k1InnerHits->clear(); k2InnerHits->clear(); 

   mumC2->clear();
   mumAngT->clear(); mumNHits->clear(); mumNPHits->clear();
   mupC2->clear();
   mupAngT->clear(); mupNHits->clear(); mupNPHits->clear();
   mumdxy->clear(); mupdxy->clear(); mumdz->clear(); mupdz->clear(); muon_dca->clear();

   tri_Dim25->clear(); tri_JpsiTk->clear(); tri_JpsiTkTk->clear(); 

   mu1soft->clear(); mu2soft->clear(); mu1tight->clear(); mu2tight->clear();
   mu1PF->clear(); mu2PF->clear(); mu1loose->clear(); mu2loose->clear(); 
 
   deltaRmum->clear(); deltaRmup->clear(); deltaRkp->clear(); deltaRkm->clear(); deltaRpion->clear();
   istruemum->clear(); istruemup->clear(); istruekp->clear(); istruekm->clear(); istruebs->clear(); istruepion->clear(); istruebc->clear();
}

bool Bspi::IsTheSame(const pat::GenericParticle& tk, const pat::Muon& mu){
  double DeltaEta = fabs(mu.eta()-tk.eta());
  double DeltaP   = fabs(mu.p()-tk.p());
  if (DeltaEta < 0.02 && DeltaP < 0.02) return true;
  return false;
}

bool Bspi::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle) {
    if (ancestor == particle ) return true;
    for (size_t i=0; i< particle->numberOfMothers(); i++) {
        if (isAncestor(ancestor,particle->mother(i))) return true;
    }
    return false;
}

double Bspi::GetLifetime(TLorentzVector b_p4, TVector3 production_vtx, TVector3 decay_vtx) {
   TVector3 pv_dv = decay_vtx - production_vtx;
   TVector3 b_p3  = b_p4.Vect();
   pv_dv.SetZ(0.);
   b_p3.SetZ(0.);
   Double_t lxy   = pv_dv.Dot(b_p3)/b_p3.Mag();
   return lxy*b_p4.M()/b_p3.Mag();
}

double Bspi::calEta (double Px, double Py, double Pz) {
	double P = sqrt(Px*Px + Py*Py + Pz*Pz);
	return 0.5*log((P + Pz) / (P - Pz));
}

double Bspi::calPhi (double Px, double Py, double Pz) {
	double phi = atan(Py / Px);
	if (Px < 0 && Py < 0) phi = phi - PI;
	if (Px < 0 && Py > 0) phi = phi + PI;
	return phi;
}

double Bspi::calEtaPhiDistance (double Px1, double Py1, double Pz1,
				double Px2, double Py2, double Pz2) {
	double phi1 = calPhi (Px1,Py1,Pz1);
	double eta1 = calEta (Px1,Py1,Pz1);
	double phi2 = calPhi (Px2,Py2,Pz2);
	double eta2 = calEta (Px2,Py2,Pz2);
	return sqrt((eta1-eta2) * (eta1-eta2) + (phi1-phi2) * (phi1-phi2));
}

void Bspi::saveTruthMatch(const edm::Event& iEvent) {
	double deltaEtaPhimum;	
	double deltaEtaPhimup;	
	double deltaEtaPhikp;	
	double deltaEtaPhikm;	
	//double deltaEtaPhibs;	
	double deltaEtaPhipion;	
for (std::vector<int>::size_type j = 0; j < Bc_mass->size(); j++) {
//		double TruthMatchPionMaxR_ = 0.004;
//	for (std::vector<int>::size_type j = 0; j < B_mass->size(); j++) { 
	
		double TruthMatchMuonMaxR_ = 0.004;
		double TruthMatchKaonMaxR_ = 0.003;
		double TruthMatchPionMaxR_ = 0.004;

		//--------------------------------
		// truth match with mu-
		//--------------------------------
		deltaEtaPhimum = calEtaPhiDistance(gen_muon1_p4.Px(), gen_muon1_p4.Py(), gen_muon1_p4.Pz(),
						   B_J_px2->at(j), B_J_py2->at(j), B_J_pz2->at(j));

		deltaRmum->push_back( deltaEtaPhimum );

		if (deltaEtaPhimum < TruthMatchMuonMaxR_) {
                        istruemum->push_back(true);
			//std::cout<<"Getting correct -Ve Muon "<<std::endl;
                } else {
                        istruemum->push_back(false);
			//std::cout<<"Getting wrong -Ve Muon "<<std::endl;
                }

		//--------------------------------
		// truth matching with mu+
		//--------------------------------
		deltaEtaPhimup = calEtaPhiDistance(gen_muon2_p4.Px(), gen_muon2_p4.Py(), gen_muon2_p4.Pz(),
                                                   B_J_px1->at(j), B_J_py1->at(j), B_J_pz1->at(j));

                deltaRmup->push_back( deltaEtaPhimup );

                if (deltaEtaPhimup < TruthMatchMuonMaxR_) {
                        istruemup->push_back(true);
			//std::cout<<"Getting correct +Ve Muon "<<std::endl;
                } else {
                        istruemup->push_back(false);
			//std::cout<<"Getting wrong +Ve Muon "<<std::endl;
                }

		//--------------------------------
		// truth match with Kaon+
		//--------------------------------
		deltaEtaPhikp = calEtaPhiDistance(gen_kaon1_p4.Px(), gen_kaon1_p4.Py(), gen_kaon1_p4.Pz(),
						B_phi_px1->at(j), B_phi_py1->at(j), B_phi_pz1->at(j));

		//h_deltaR->Fill(deltaEtaPhi);	
		deltaRkp->push_back( deltaEtaPhikp );

		if (deltaEtaPhikp < TruthMatchKaonMaxR_) {
			istruekp->push_back(true);
			//std::cout<<"Getting correct +Ve Kaon "<<std::endl;
		} else {
			istruekp->push_back(false);
			//std::cout<<"Getting wrong +Ve Kaon "<<std::endl;
		}

		//---------------------------------
		// truth match with Kaoun-
		//---------------------------------
		deltaEtaPhikm = calEtaPhiDistance(gen_kaon2_p4.Px(), gen_kaon2_p4.Py(), gen_kaon2_p4.Pz(),
                                                  B_phi_px2->at(j), B_phi_py2->at(j), B_phi_pz2->at(j));

		deltaRkm->push_back( deltaEtaPhikm );

                if (deltaEtaPhikm < TruthMatchKaonMaxR_) {
                        istruekm->push_back(true);
			//std::cout<<"Getting correct -Ve Kaon "<<std::endl;
                } else {
                        istruekm->push_back(false);
			//std::cout<<"Getting wrong -Ve Kaon "<<std::endl;
		}
		//---------------------------------
		// truth match with Bs or Bs bar
		//--------------------------------
		if ( istruemum->back() && istruemup->back() && istruekm->back() && istruekp->back() ) {
			istruebs->push_back(true);
			//std::cout<<"Getting correct Bs "<<std::endl;
		} else {
			istruebs->push_back(false);
			//std::cout<<"Getting wrong Bs "<<std::endl;
		}
	//}
		//----------------------------------
		// truth match with pion track
		//----------------------------------
		deltaEtaPhipion = calEtaPhiDistance(gen_pion_p4.Px(), gen_pion_p4.Py(), gen_pion_p4.Pz(),
                                                    B_pion_px_track->at(j), B_pion_py_track->at(j), B_pion_pz_track->at(j));
		//deltaEtaPhipion = calEtaPhiDistance(gen_pion_p4.Px(), gen_pion_p4.Py(), gen_pion_p4.Pz(),
                  //                                  B_pion_px->at(j), B_pion_py->at(j), B_pion_pz->at(j));

		//std::cout<< "delatR of pion================================================>>>>>>>>>>>>>>>>>>>>>>>>>>> " << deltaEtaPhipion << std::endl;
		deltaRpion->push_back( deltaEtaPhipion );

		if (deltaEtaPhipion < TruthMatchPionMaxR_) {
		//if (deltaEtaPhipion < 0.1)  // Entery = 3
		//if (deltaEtaPhipion < 0.2)  // Entery = 7
		//if (deltaEtaPhipion < 0.3)  // Entery = 13
		//if (deltaEtaPhipion < 0.4)  // Entery = 23
		//if (deltaEtaPhipion < 0.5)  // Entery = 27
		//if (deltaEtaPhipion < 0.6)  // Entery = 32
		//if (deltaEtaPhipion < 0.7)  // Entery = 34
		//if (deltaEtaPhipion < 0.8)  // Entery = 39 // <--
		//if (deltaEtaPhipion < 0.9)  // Entery = 39
		//if (deltaEtaPhipion < 1.0)  // Entery = 39
		//if (deltaEtaPhipion < 0.09) // Entery = 1
                        istruepion->push_back(true);
			//std::cout<<"Getting correct pion "<<std::endl;
                } else {
                        istruepion->push_back(false);
			//std::cout<<"Getting wrong pion "<<std::endl;
                }

		//----------------------------------
		// truth match with Bc
		//----------------------------------
		//if ( istruemum->back() && istruemup->back() && istruekm->back() && istruekp->back() && istruepion->back() ) 
		if ( istruebs->back() && istruepion->back() ) {
                        istruebc->push_back(true);
			std::cout<<"========XXXXXXXXXXXXXXXXXXXXXXXXXXX>>>>>>>>>>>>>>>>> GETTING CORRECT Bc MESON <<<<<<<<<<<<<<<<<<<<XXXXXXXXXXXXXXXXXXXXX========= "<<std::endl;
			//T_Bspion_mass->push_back(Bspion_mass->at(l));
			T_Bspion_mass->push_back(Bc_mass->at(j));
                } else {
                        istruebc->push_back(false);
			Notmatch_Bspion_mass->push_back(Bc_mass->at(j));
                        //std::cout<<"Getting wrong Bc "<<std::endl;
                }

	}

}

//*****************************************************************************************************************************

void Bspi::saveIP(const RefCountedKinematicTree& vertexFitTreeBspi,
		     const reco::VertexCollection& vertices, int & pvIndex){

	vertexFitTreeBspi->movePointerToTheTop();
	RefCountedKinematicParticle bspiCandMC = vertexFitTreeBspi->currentParticle();
        RefCountedKinematicVertex bspiDecayVertexMC = vertexFitTreeBspi->currentDecayVertex();
	
	//std::cout<<" inside the save IP vertices size "<<vertices.size()<<std::endl;	
	auto candTransientTrack = bspiCandMC->refittedTransientTrack();

	// find the first primary vertex
	const reco::Vertex* bestVertex_t(0);
	int bestVertexIndex(-1);
	double minDistance(999.);
	for ( unsigned int i = 0; i<vertices.size(); ++i ){
		const auto & vertex = vertices.at(i);
		auto impactParameter3D = IPTools::absoluteImpactParameter3D(candTransientTrack, vertex);
		//std::cout<<" inside the vertex loop "<<i<<"\t"<<impactParameter3D.second.value()<<endl;
		if (impactParameter3D.first and impactParameter3D.second.value() < minDistance){
			minDistance = impactParameter3D.second.value();
			bestVertex_t = &vertex;
			bestVertexIndex = i;
			//std::cout<<" index list "<<i<< "\t"<<minDistance<<std::endl;
		}
	}
	pvIndex = bestVertexIndex;
	//std::cout<<" inside the save IP index of vtx "<<bestVertexIndex<<std::endl;
	pVtxIPX->push_back( bestVertex_t->x());
	pVtxIPY->push_back( bestVertex_t->y());
  	pVtxIPZ->push_back( bestVertex_t->z());
  	pVtxIPXE->push_back( bestVertex_t->covariance(0, 0) );
  	pVtxIPYE->push_back( bestVertex_t->covariance(1, 1) );
  	pVtxIPZE->push_back( bestVertex_t->covariance(2, 2) );
  	pVtxIPXYE->push_back(bestVertex_t->covariance(0, 1) );
  	pVtxIPXZE->push_back(bestVertex_t->covariance(0, 2) );  
  	pVtxIPYZE->push_back( bestVertex_t->covariance(1, 2) );
  	pVtxIPCL->push_back(  ChiSquaredProbability((double)(bestVertex_t->chi2()),(double)(bestVertex_t->ndof())) );

	//find second best vertex
	const reco::Vertex* bestVertex2(0);
	//int bestVertexIndex2(-1);
	double minDistance2(999.);
	for ( unsigned int i = 0; i<vertices.size(); ++i ){
		const auto & vertex = vertices.at(i);
		auto impactParameter3D = IPTools::absoluteImpactParameter3D(candTransientTrack, vertex);
		if (impactParameter3D.first and impactParameter3D.second.value() < minDistance2 and impactParameter3D.second.value() > minDistance){
			minDistance2 = impactParameter3D.second.value();
			bestVertex2 = &vertex;
			//bestVertexIndex2 = i;
		}
	}	
	//std::cout<<" inside the save IP second vertex "<<std::endl;
	//  if (! bestVertex) continue;
	auto impactParameter3D = IPTools::absoluteImpactParameter3D(candTransientTrack, *bestVertex_t);
	auto impactParameterZ  = IPTools::signedDecayLength3D(candTransientTrack, GlobalVector(0,0,1), *bestVertex_t);
	//pv = bestVertex;
	//pvIndex = bestVertexIndex;
	double longitudinalImpactParameter(0.0), longitudinalImpactParameterErr(0.0);
	double distaceOfClosestApproach(0.0), distaceOfClosestApproachErr(0.0) ;
	if (impactParameterZ.first) {
		longitudinalImpactParameter    = impactParameterZ.second.value();
		longitudinalImpactParameterErr = impactParameterZ.second.error();
	}
	if(impactParameter3D.first) {
		distaceOfClosestApproach       = impactParameter3D.second.value();
		distaceOfClosestApproachErr    = impactParameter3D.second.error();
	}
	B_pvip->push_back(distaceOfClosestApproach);
	B_pviperr->push_back(distaceOfClosestApproachErr);
  	B_pvips->push_back(distaceOfClosestApproachErr/distaceOfClosestApproach);
  	B_pvlzip->push_back(longitudinalImpactParameter);
  	B_pvlziperr->push_back(longitudinalImpactParameterErr);
  	B_pvlzips->push_back(longitudinalImpactParameterErr/longitudinalImpactParameter);

	VertexDistance3D distance3D;
  	auto dist = distance3D.distance(*bestVertex_t, bspiDecayVertexMC->vertexState() );
  	double decayLength(-1.), decayLengthErr(0);
  	decayLength    = dist.value();
  	decayLengthErr = dist.error();
  
  	VertexDistanceXY distanceXY;
  	auto distXY = distanceXY.distance(*bestVertex_t, bspiDecayVertexMC->vertexState() );

  	B_l3d ->push_back(decayLength);
  	B_l3dE ->push_back(decayLengthErr);
  	//B_lxy ->push_back(distXY.value());
  	Bc_lxy ->push_back(distXY.value());
  	B_lxyE ->push_back(distXY.error());

	if (bestVertex2){
		double longitudinalImpactParameter2(0.0), longitudinalImpactParameter2Err(0.0);
    		double distaceOfClosestApproach2(0.0), distaceOfClosestApproach2Err(0.0) ;
    		auto impactParameter3D2 = IPTools::absoluteImpactParameter3D(candTransientTrack, *bestVertex2);
    		auto impactParameterZ2  = IPTools::signedDecayLength3D(candTransientTrack, GlobalVector(0,0,1), *bestVertex2);
		if (impactParameterZ2.first) {
      			longitudinalImpactParameter2    = impactParameterZ2.second.value();
      			longitudinalImpactParameter2Err = impactParameterZ2.second.error();
   		 }
		if (impactParameter3D2.first) {
      			distaceOfClosestApproach2       = impactParameter3D2.second.value();
      			distaceOfClosestApproach2Err    = impactParameter3D2.second.error();
    		}
		B_pv2ip->push_back(distaceOfClosestApproach2);
    		B_pv2iperr->push_back(distaceOfClosestApproach2Err);
    		B_pv2ips->push_back(distaceOfClosestApproach2Err/distaceOfClosestApproach2);
    		B_pv2lzip->push_back(longitudinalImpactParameter2);
    		B_pv2lziperr->push_back(longitudinalImpactParameter2Err);
    		B_pv2lzips->push_back(longitudinalImpactParameter2Err/longitudinalImpactParameter2);

		// compute decay length
		VertexDistance3D distance3D;
    		auto dist = distance3D.distance(*bestVertex2, bspiDecayVertexMC->vertexState() );
    		B_l3d_pv2 ->push_back(dist.value());
    		B_l3dE_pv2 ->push_back(dist.error());

	}

	TVector3 plab(bspiCandMC->currentState().globalMomentum().x(),
		      bspiCandMC->currentState().globalMomentum().y(),
		      bspiCandMC->currentState().globalMomentum().z());
	TVector3 p1(bestVertex_t->x(), bestVertex_t->y(), bestVertex_t->z());
	TVector3 p2(bspiDecayVertexMC->vertexState().position().x(), 
	            bspiDecayVertexMC->vertexState().position().y(), 
	            bspiDecayVertexMC->vertexState().position().z());
	TVector3 pDiff = p2-p1;
	TVector3 pDiffXY = TVector3(pDiff.X(), pDiff.Y(), 0.);
	TVector3 ptrans  = TVector3(plab.X(), plab.Y(), 0.);
	double cosAlpha(-999.),cosAlphaXY(-999.); //,decayTime(-999.),decayTimeError(-999.), decayTimeXY(-999.),decayTimeXYError(-999.);
	cosAlpha  = plab.Dot(pDiff) / (plab.Mag() * pDiff.Mag());
	cosAlphaXY  = ptrans.Dot(pDiffXY) / (ptrans.Mag() * pDiffXY.Mag());
	B_cosalpha -> push_back(cosAlpha);
	//B_cosalphaxy -> push_back(cosAlphaXY);
	Bc_cosalphaxy -> push_back(cosAlphaXY);
  	alpha->push_back(TMath::ACos(cosAlpha));
//
}

//*****************************************************************************************************************************

// ------------ method called once each job just before starting event loop  ------------

void 
Bspi::beginJob()
{
  std::cout << "Beginning analyzer job with value of isMC_ = " << isMC_ << std::endl;

  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("ntuple","Bc->BsPi ntuple");
  //tree_ = fs->make<TTree>("genTree","Bc->BsPi genLevel");
  //h_deltaR = fs->make<TH1D>("deltaR", "DR", 1000, 0.0, 0.1),

  tree_->Branch("nB",&nB,"nB/i");
  tree_->Branch("nMu",&nMu,"nMu/i");
	  
  //tree_->Branch("Bspi_massi_vertex", &Bspi_mass_vertex);
  //tree_->Branch("Bspion_mass", &Bspion_mass);
  tree_->Branch("Bc_mass",&Bc_mass);
  tree_->Branch("Bc_pt",&Bc_pt);
  tree_->Branch("Bc_px",&Bc_px);
  tree_->Branch("Bc_py",&Bc_py);
  tree_->Branch("Bc_pz",&Bc_pz);
  tree_->Branch("Bc_eta",&Bc_eta);
  tree_->Branch("Bc_phi",&Bc_phi);
  tree_->Branch("Bc_ct", &Bc_ct);
  tree_->Branch("Bc_charge",&Bc_charge);
  tree_->Branch("Bcvtxcl",&Bcvtxcl);
  tree_->Branch("Bc_prob",&Bc_prob);

  tree_->Branch("deltaMass", &deltaMass);
  tree_->Branch("T_Bspion_mass", &T_Bspion_mass);
  tree_->Branch("Notmatch_Bspion_mass", &Notmatch_Bspion_mass);

  tree_->Branch("Bc_DecayVtxX"       , &Bc_DecayVtxX          );
  tree_->Branch("Bc_DecayVtxY"       , &Bc_DecayVtxY          );
  tree_->Branch("Bc_DecayVtxZ"       , &Bc_DecayVtxZ          );
  tree_->Branch("Bc_DecayVtxXE"      , &Bc_DecayVtxXE         );
  tree_->Branch("Bc_DecayVtxYE"      , &Bc_DecayVtxYE         );
  tree_->Branch("Bc_DecayVtxZE"      , &Bc_DecayVtxZE         );

  //tree_->Branch("Bc_DecayVtx_vtxfit_X"       , &Bc_DecayVtx_vtxfit_X          );
  //tree_->Branch("Bc_DecayVtx_vtxfit_Y"       , &Bc_DecayVtx_vtxfit_Y          );
  //tree_->Branch("Bc_DecayVtx_vtxfit_Z"       , &Bc_DecayVtx_vtxfit_Z          );
  //tree_->Branch("Bc_DecayVtx_vtxfit_XE"       , &Bc_DecayVtx_vtxfit_XE          );
  //tree_->Branch("Bc_DecayVtx_vtxfit_YE"       , &Bc_DecayVtx_vtxfit_YE          );
  //tree_->Branch("Bc_DecayVtx_vtxfit_ZE"       , &Bc_DecayVtx_vtxfit_ZE         );
  //tree_->Branch("Bc_DecayVtx_vtxfit_XYE"       , &Bc_DecayVtx_vtxfit_XYE          );
  //tree_->Branch("Bc_DecayVtx_vtxfit_XZE"       , &Bc_DecayVtx_vtxfit_XZE          );
  //tree_->Branch("Bc_DecayVtx_vtxfit_YZE"       , &Bc_DecayVtx_vtxfit_YZE          );

  tree_->Branch("B_mass", &B_mass);
  tree_->Branch("B_px", &B_px);
  tree_->Branch("B_py", &B_py);
  tree_->Branch("B_pz", &B_pz);
  tree_->Branch("B_pt", &B_pt);
  tree_->Branch("B_eta", &B_eta);
  tree_->Branch("B_phi", &B_phi);

  tree_->Branch("B_phi_mass", &B_phi_mass);
  tree_->Branch("B_phi_px", &B_phi_px);
  tree_->Branch("B_phi_py", &B_phi_py);
  tree_->Branch("B_phi_pz", &B_phi_pz);
  tree_->Branch("B_phi_pt", &B_phi_pt);
  tree_->Branch("B_phi_eta", &B_phi_eta);
  tree_->Branch("B_phi_phi", &B_phi_phi);

  tree_->Branch("B_J_mass", &B_J_mass);
  tree_->Branch("B_J_px", &B_J_px);
  tree_->Branch("B_J_py", &B_J_py);
  tree_->Branch("B_J_pz", &B_J_pz);
  tree_->Branch("B_J_pt", &B_J_pt);
  tree_->Branch("B_J_eta", &B_J_eta);
  tree_->Branch("B_J_phi", &B_J_phi);
  
  //tree_->Branch("B_pion_px", &B_pion_px);
  //tree_->Branch("B_pion_py", &B_pion_py);
  //tree_->Branch("B_pion_pz", &B_pion_pz);
  tree_->Branch("pion_track_normchi2"   , &pion_track_normchi2      );
  tree_->Branch("pion_NumHits"        , &pion_NumHits           );
  tree_->Branch("pion_NumPixelHits"        , &pion_NumPixelHits           );

  tree_->Branch("B_pion_px_track", &B_pion_px_track);
  tree_->Branch("B_pion_py_track", &B_pion_py_track);
  tree_->Branch("B_pion_pz_track", &B_pion_pz_track);
  tree_->Branch("B_pion_pt_track", &B_pion_pt_track);
  tree_->Branch("B_pion_charge", &B_pion_charge); 
  tree_->Branch("pion_dxy", &pion_dxy);
  tree_->Branch("pion_dz", &pion_dz);
  tree_->Branch("pion_dxy_", &pion_dxy_);
  tree_->Branch("pion_dz_", &pion_dz_);
  tree_->Branch("pion_NTrackerLayers", &pion_NTrackerLayers);
  tree_->Branch("pion_NPixelLayers", &pion_NPixelLayers);

  tree_->Branch("B_phi_px1", &B_phi_px1);
  tree_->Branch("B_phi_py1", &B_phi_py1);
  tree_->Branch("B_phi_pz1", &B_phi_pz1);
  tree_->Branch("B_phi_pt1", &B_phi_pt1);
  tree_->Branch("B_phi_eta1", &B_phi_eta1);
  tree_->Branch("B_phi_phi1", &B_phi_phi1);
  tree_->Branch("B_phi_px1_track", &B_phi_px1_track);
  tree_->Branch("B_phi_py1_track", &B_phi_py1_track);
  tree_->Branch("B_phi_pz1_track", &B_phi_pz1_track);
  tree_->Branch("B_phi_charge1", &B_phi_charge1); 
 
  tree_->Branch("B_phi_px2", &B_phi_px2);
  tree_->Branch("B_phi_py2", &B_phi_py2);
  tree_->Branch("B_phi_pz2", &B_phi_pz2);
  tree_->Branch("B_phi_pt2", &B_phi_pt2);
  tree_->Branch("B_phi_eta2", &B_phi_eta2);
  tree_->Branch("B_phi_phi2", &B_phi_phi2);
  tree_->Branch("B_phi_px2_track", &B_phi_px2_track);
  tree_->Branch("B_phi_py2_track", &B_phi_py2_track);
  tree_->Branch("B_phi_pz2_track", &B_phi_pz2_track);
  tree_->Branch("B_phi_charge2", &B_phi_charge2);

  tree_->Branch("B_J_px1", &B_J_px1);
  tree_->Branch("B_J_py1", &B_J_py1);
  tree_->Branch("B_J_pz1", &B_J_pz1);
  tree_->Branch("B_J_pt1", &B_J_pt1);
  tree_->Branch("B_J_eta1", &B_J_eta1);
  tree_->Branch("B_J_phi1", &B_J_phi1);
  tree_->Branch("B_J_charge1", &B_J_charge1);

  tree_->Branch("B_J_px2", &B_J_px2);
  tree_->Branch("B_J_py2", &B_J_py2);
  tree_->Branch("B_J_pz2", &B_J_pz2);
  tree_->Branch("B_J_pt2", &B_J_pt2);
  tree_->Branch("B_J_eta2", &B_J_eta2);
  tree_->Branch("B_J_phi2", &B_J_phi2);
  tree_->Branch("B_J_charge2", &B_J_charge2);

  tree_->Branch("Bc_chi2",    &Bc_chi2);
  tree_->Branch("B_chi2",    &B_chi2);
  tree_->Branch("B_J_chi2",  &B_J_chi2);

  tree_->Branch("B_Prob",    &B_Prob);
  tree_->Branch("B_J_Prob",  &B_J_Prob);
     
  tree_->Branch("B_DecayVtxX",     &B_DecayVtxX);
  tree_->Branch("B_DecayVtxY",     &B_DecayVtxY);
  tree_->Branch("B_DecayVtxZ",     &B_DecayVtxZ);
  tree_->Branch("B_DecayVtxXE",    &B_DecayVtxXE);
  tree_->Branch("B_DecayVtxYE",    &B_DecayVtxYE);
  tree_->Branch("B_DecayVtxZE",    &B_DecayVtxZE);
  tree_->Branch("B_DecayVtxXYE",    &B_DecayVtxXYE);
  tree_->Branch("B_DecayVtxXZE",    &B_DecayVtxXZE);
  tree_->Branch("B_DecayVtxYZE",    &B_DecayVtxYZE);

  tree_->Branch("priVtxX",&priVtxX, "priVtxX/f");
  tree_->Branch("priVtxY",&priVtxY, "priVtxY/f");
  tree_->Branch("priVtxZ",&priVtxZ, "priVtxZ/f");
  tree_->Branch("priVtxXE",&priVtxXE, "priVtxXE/f");
  tree_->Branch("priVtxYE",&priVtxYE, "priVtxYE/f");
  tree_->Branch("priVtxZE",&priVtxZE, "priVtxZE/f");
  tree_->Branch("priVtxXYE",&priVtxXYE, "priVtxXYE/f");
  tree_->Branch("priVtxXZE",&priVtxXZE, "priVtxXZE/f");
  tree_->Branch("priVtxYZE",&priVtxYZE, "priVtxYZE/f");
  tree_->Branch("priVtxCL",&priVtxCL, "priVtxCL/f");
  
  tree_->Branch("nVtx",       &nVtx);
  tree_->Branch("run",        &run,       "run/I");
  tree_->Branch("event",        &event,     "event/I");
  tree_->Branch("lumiblock",&lumiblock,"lumiblock/I");
    
  // *************************
  tree_->Branch("pVtxIPX",     &pVtxIPX);
  tree_->Branch("pVtxIPY",     &pVtxIPY);
  tree_->Branch("pVtxIPZ",     &pVtxIPZ);
  tree_->Branch("pVtxIPXE",     &pVtxIPXE);
  tree_->Branch("pVtxIPYE",     &pVtxIPYE);
  tree_->Branch("pVtxIPZE",     &pVtxIPZE);
  tree_->Branch("pVtxIPXYE",     &pVtxIPXYE);
  tree_->Branch("pVtxIPXZE",     &pVtxIPXZE);
  tree_->Branch("pVtxIPYZE",     &pVtxIPYZE);
  tree_->Branch("pVtxIPCL",     &pVtxIPCL);  
  
  tree_->Branch("B_pvip",&B_pvip);
  tree_->Branch("B_pviperr",&B_pviperr);
  tree_->Branch("B_pvips",&B_pvips);
  tree_->Branch("B_pvlzip",&B_pvlzip);
  tree_->Branch("B_pvlziperr",&B_pvlziperr);
  tree_->Branch("B_pvlzips",&B_pvlzips);
  tree_->Branch("B_pv2ip",&B_pv2ip);
  tree_->Branch("B_pv2iperr",&B_pv2iperr);
  tree_->Branch("B_pv2ips",&B_pv2ips);
  tree_->Branch("B_pv2lzip",&B_pv2lzip);
  tree_->Branch("B_pv2lziperr",&B_pv2lziperr);
  tree_->Branch("B_pv2lzips",&B_pv2lzips);
  tree_->Branch("B_l3d_pv2",&B_l3d_pv2);
  tree_->Branch("B_l3dE_pv2",&B_l3dE_pv2);
  tree_->Branch("B_l3d",&B_l3d);
  tree_->Branch("B_l3dE",&B_l3dE);
  //tree_->Branch("B_lxy", &B_lxy);
  tree_->Branch("Bc_lxy", &Bc_lxy);
  tree_->Branch("B_lxyE",&B_lxyE);
  tree_->Branch("B_cosalpha",&B_cosalpha);  
  //tree_->Branch("B_cosalphaxy",&B_cosalphaxy);
  tree_->Branch("Bc_cosalphaxy",&Bc_cosalphaxy);
  tree_->Branch("alpha",&alpha);
  // *************************
  
  tree_->Branch("k1dxy",&k1dxy);
  tree_->Branch("k2dxy",&k2dxy);
  tree_->Branch("k1dz",&k1dz);
  tree_->Branch("k2dz",&k2dz);

  tree_->Branch("k1dxy_e",&k1dxy_e);
  tree_->Branch("k2dxy_e",&k2dxy_e);
  tree_->Branch("k1dz_e",&k1dz_e);
  tree_->Branch("k2dz_e",&k2dz_e);

  tree_->Branch("k1InnerHits",&k1InnerHits);
  tree_->Branch("k2InnerHits",&k2InnerHits);

  tree_->Branch("mumC2",&mumC2);  
  tree_->Branch("mumAngT",&mumAngT);
  tree_->Branch("mumNHits",&mumNHits);
  tree_->Branch("mumNPHits",&mumNPHits);
  tree_->Branch("mupC2",&mupC2);  
  tree_->Branch("mupAngT",&mupAngT);
  tree_->Branch("mupNHits",&mupNHits);
  tree_->Branch("mupNPHits",&mupNPHits);
  tree_->Branch("mumdxy",&mumdxy);
  tree_->Branch("mupdxy",&mupdxy);
  tree_->Branch("mumdz",&mumdz);
  tree_->Branch("mupdz",&mupdz);
  tree_->Branch("muon_dca",&muon_dca);

  tree_->Branch("tri_Dim25",&tri_Dim25);
  tree_->Branch("tri_JpsiTk",&tri_JpsiTk);
  tree_->Branch("tri_JpsiTkTk",&tri_JpsiTkTk); 

  tree_->Branch("mu1soft",&mu1soft);
  tree_->Branch("mu2soft",&mu2soft);
  tree_->Branch("mu1tight",&mu1tight);
  tree_->Branch("mu2tight",&mu2tight);
  tree_->Branch("mu1PF",&mu1PF);
  tree_->Branch("mu2PF",&mu2PF);
  tree_->Branch("mu1loose",&mu1loose);
  tree_->Branch("mu2loose",&mu2loose);

  // gen
  if (isMC_) {
     tree_->Branch("gen_bc_p4",     "TLorentzVector",  &gen_bc_p4);
     tree_->Branch("gen_bc_vtx",    "TVector3",        &gen_bc_vtx);
     tree_->Branch("gen_bs_p4",     "TLorentzVector",  &gen_bs_p4);
     tree_->Branch("gen_bs_vtx",    "TVector3",        &gen_bs_vtx);
     //tree_->Branch("gen_phi_p4",   "TLorentzVector",  &gen_phi_p4);
     tree_->Branch("gen_jpsi_vtx",  "TVector3",        &gen_jpsi_vtx);
     tree_->Branch("gen_muon1_p4",  "TLorentzVector",  &gen_muon1_p4);
     tree_->Branch("gen_muon2_p4",  "TLorentzVector",  &gen_muon2_p4);
     tree_->Branch("gen_jpsi_p4",   "TLorentzVector",  &gen_jpsi_p4);
     tree_->Branch("gen_kaon1_p4",  "TLorentzVector",  &gen_kaon1_p4);
     tree_->Branch("gen_kaon2_p4",  "TLorentzVector",  &gen_kaon2_p4);
     tree_->Branch("gen_pion_p4",     "TLorentzVector",  &gen_pion_p4);
     tree_->Branch("gen_b_ct",     &gen_b_ct,        "gen_b_ct/F");
  }

     tree_->Branch("deltaRmum",   &deltaRmum  );
     tree_->Branch("deltaRmup",   &deltaRmup  );
     tree_->Branch("deltaRkp",   &deltaRkp  );
     tree_->Branch("deltaRkm",   &deltaRkm  );
     tree_->Branch("deltaRpion",   &deltaRpion  );

     tree_->Branch("istruemum",   &istruemum  );
     tree_->Branch("istruemup",   &istruemup  );
     tree_->Branch("istruekp",   &istruekp  );
     tree_->Branch("istruekm",   &istruekm  );
     tree_->Branch("istruebs",   &istruebs  );
     tree_->Branch("istruepion",   &istruepion  );
     tree_->Branch("istruebc",   &istruebc  );

}


// ------------ method called once each job just after ending the event loop  ------------
void Bspi::endJob() {
  tree_->GetDirectory()->cd();
  tree_->Write();
}

//define this as a plug-in
DEFINE_FWK_MODULE(Bspi);

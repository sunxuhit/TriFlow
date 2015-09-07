#include "StTriFlowV0.h"
#include "StTriFlowConstants.h"
#include "StTriFlowCut.h"
#include "StRoot/StPicoDstMaker/StPicoDst.h"
#include "StRoot/StPicoDstMaker/StPicoEvent.h"
#include "StRoot/StPicoDstMaker/StPicoTrack.h"
#include "StRoot/StAlexPhiMesonEvent/StAlexPhiMesonEvent.h"
#include "StRoot/StV0Event/StV0Event.h"
#include "StRoot/StV0TofCorrection/StV0TofCorrection.h"
#include <vector>
#include "TLorentzVector.h"
#include "StThreeVectorF.hh"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "StarClassLibrary/StThreeVector.hh"
#include "StarClassLibrary/SystemOfUnits.h"
#include "TMath.h"
#include "TObject.h"
#include "TVector3.h"
#include "StTriFlow2ndVertexFinder.h"
#include "StLorentzVectorD.hh"

ClassImp(StTriFlowV0)

//------------------------------------------------------------------------------------------------------------------
StTriFlowV0::StTriFlowV0(Int_t energy)
{
  mEnergy = energy;
}

StTriFlowV0::~StTriFlowV0()
{
}

//------------------------------------------------------------------------------------------------------------------

void StTriFlowV0::InitPhi()
{
  mTriFlowCut = new StTriFlowCut(mEnergy);
  TString HistName = "Mass2_pt";
  h_Mass2 = new TH2F(HistName.Data(),HistName.Data(),20,0.2,5.0,200,0.98,1.08);

  for(Int_t cent = 0; cent < TriFlow::Bin_Centrality; cent++)
  {
    for(Int_t vz = 0; vz < TriFlow::Bin_VertexZ; vz++)
    {
      for(Int_t phi_psi = 0; phi_psi < TriFlow::Bin_Phi_Psi; phi_psi++)
      {
        mEventCounter2[cent][vz][phi_psi] = 0;
	clear_phi(cent,vz,phi_psi);
      }
    }
  }

  mXuPhiMesonEvent = new StAlexPhiMesonEvent();
  mTree_Phi = new TTree("XuPhiMesonEvent","XuPhiMesonEvent");
  mTree_Phi->Branch("phi_flow_branch","StAlexPhiMesonEvent",&mXuPhiMesonEvent);
  mTree_Phi->SetAutoSave(5000000);
}

void StTriFlowV0::InitLambda()
{
  mTriFlowCut = new StTriFlowCut(mEnergy);
  mTofCorr = new StV0TofCorrection();
  TString HistName = "Mass2_pt";
  h_Mass2 = new TH2F(HistName.Data(),HistName.Data(),20,0.2,5.0,200,1.06,1.20);
  HistName = "Mass2_pt_sub";
  h_Mass2_sub = new TH2F(HistName.Data(),HistName.Data(),20,0.2,5.0,200,1.06,1.20);
  HistName = "Mass2_pt_K0s";
  h_Mass2_K0s = new TH2F(HistName.Data(),HistName.Data(),20,0.2,5.0,200,0.40,0.60);
  HistName = "Mass2_p";
  h_Mass2_p = new TH2F(HistName.Data(),HistName.Data(),200,-5.0,5.0,200,-0.2,1.2);

  for(Int_t cent = 0; cent < TriFlow::Bin_Centrality; cent++)
  {
    for(Int_t vz = 0; vz < TriFlow::Bin_VertexZ; vz++)
    {
      for(Int_t phi_psi = 0; phi_psi < TriFlow::Bin_Phi_Psi; phi_psi++)
      {
        mEventCounter2[cent][vz][phi_psi] = 0;
	clear_Lambda(cent,vz,phi_psi);
      }
    }
  }

  mV0Event = new StV0Event();
  mTree_Lambda = new TTree("LambdaEvent","LambdaEvent");
  mTree_Lambda->Branch("Lambda_flow_branch","StV0Event",&mV0Event);
  mTree_Lambda->SetAutoSave(5000000);
}

void StTriFlowV0::InitAntiLambda()
{
  mTriFlowCut = new StTriFlowCut(mEnergy);
  mTofCorr = new StV0TofCorrection();
  TString HistName = "Mass2_pt";
  h_Mass2 = new TH2F(HistName.Data(),HistName.Data(),20,0.2,5.0,200,1.06,1.20);
  HistName = "Mass2_pt_sub";
  h_Mass2_sub = new TH2F(HistName.Data(),HistName.Data(),20,0.2,5.0,200,1.06,1.20);
  HistName = "Mass2_pt_K0s";
  h_Mass2_K0s = new TH2F(HistName.Data(),HistName.Data(),20,0.2,5.0,200,0.40,0.60);
  HistName = "Mass2_p";
  h_Mass2_p = new TH2F(HistName.Data(),HistName.Data(),200,-5.0,5.0,200,-0.2,1.2);

  for(Int_t cent = 0; cent < TriFlow::Bin_Centrality; cent++)
  {
    for(Int_t vz = 0; vz < TriFlow::Bin_VertexZ; vz++)
    {
      for(Int_t phi_psi = 0; phi_psi < TriFlow::Bin_Phi_Psi; phi_psi++)
      {
        mEventCounter2[cent][vz][phi_psi] = 0;
	clear_AntiLambda(cent,vz,phi_psi);
      }
    }
  }

  mV0Event = new StV0Event();
  mTree_AntiLambda = new TTree("antiLambdaEvent","antiLambdaEvent");
  mTree_AntiLambda->Branch("antiLambda_flow_branch","StV0Event",&mV0Event);
  mTree_AntiLambda->SetAutoSave(5000000);
}

void StTriFlowV0::InitK0S()
{
  mTriFlowCut = new StTriFlowCut(mEnergy);
  mTofCorr = new StV0TofCorrection();
  h_Mass2 = new TH2F("h_Mass2","h_Mass2",20,0.2,5.0,200,0.4,0.6);
  h_Mass2_p = new TH2F("h_Mass2_p","h_Mass2_p",200,-5.0,5.0,200,-0.2,1.2);
  h_Mass2_Lambda = new TH2F("h_Mass2_Lambda","h_Mass2_Lambda",20,0.2,5.0,200,1.06,1.2);
  h_Mass2_AntiLambda = new TH2F("h_Mass2_AntiLambda","h_Mass2_AntiLambda",20,0.2,5.0,200,1.06,1.2);
  h_Mass2_sub = new TH2F("h_Mass2_sub","h_Mass2_sub",20,0.2,5.0,200,0.4,0.6);

  for(Int_t cent = 0; cent < TriFlow::Bin_Centrality; cent++)
  {
    for(Int_t vz = 0; vz < TriFlow::Bin_VertexZ; vz++)
    {
      for(Int_t phi_psi = 0; phi_psi < TriFlow::Bin_Phi_Psi; phi_psi++)
      {
        mEventCounter2[cent][vz][phi_psi] = 0;
	clear_K0S(cent,vz,phi_psi);
      }
    }
  }

  mV0Event = new StV0Event();
  mTree_K0S = new TTree("K0SEvent","K0SEvent");
  mTree_K0S->Branch("K0S_flow_branch","StV0Event",&mV0Event);
  mTree_K0S->SetAutoSave(5000000);
}
//------------------------------------------------------------------------------------------------------------------

void StTriFlowV0::WritePhiMass2()
{
  h_Mass2->Write();
  mTree_Phi->Write("",TObject::kOverwrite);
}

void StTriFlowV0::WriteLambdaMass2()
{
  h_Mass2->Write();
  h_Mass2_sub->Write();
  h_Mass2_K0s->Write();
  h_Mass2_p->Write();
  mTree_Lambda->Write("",TObject::kOverwrite);
}

void StTriFlowV0::WriteAntiLambdaMass2()
{
  h_Mass2->Write();
  h_Mass2_sub->Write();
  h_Mass2_K0s->Write();
  h_Mass2_p->Write();
  mTree_AntiLambda->Write("",TObject::kOverwrite);
}

void StTriFlowV0::WriteK0SMass2()
{
  h_Mass2->Write();
  h_Mass2_p->Write();
  h_Mass2_Lambda->Write();
  h_Mass2_AntiLambda->Write();
  h_Mass2_sub->Write();
  mTree_K0S->Write("",TObject::kOverwrite);
}
//------------------------------------------------------------------------------------------------------------------

void StTriFlowV0::clear_phi(Int_t cent9, Int_t Bin_vz, Int_t Bin_Psi2)
{
  mPrimaryvertex[cent9][Bin_vz][Bin_Psi2].clear();
  mRefMult[cent9][Bin_vz][Bin_Psi2].clear();
  mCentrality[cent9][Bin_vz][Bin_Psi2].clear();
  mRunId[cent9][Bin_vz][Bin_Psi2].clear();
  mN_prim[cent9][Bin_vz][Bin_Psi2].clear();
  mN_non_prim[cent9][Bin_vz][Bin_Psi2].clear();
  mN_Tof_match[cent9][Bin_vz][Bin_Psi2].clear();
  mZDCx[cent9][Bin_vz][Bin_Psi2].clear();
  mBBCx[cent9][Bin_vz][Bin_Psi2].clear();
  mVzVpd[cent9][Bin_vz][Bin_Psi2].clear();

  for(Int_t j = 0; j < TriFlow::EtaGap_total; j++)
  {
    mQ2East[cent9][Bin_vz][Bin_Psi2][j].clear();
    mQ2West[cent9][Bin_vz][Bin_Psi2][j].clear();
    mQ3East[cent9][Bin_vz][Bin_Psi2][j].clear();
    mQ3West[cent9][Bin_vz][Bin_Psi2][j].clear();
  }

  for(Int_t Bin_Event = 0; Bin_Event < TriFlow::Buffer_depth; Bin_Event++)
  {
    for(Int_t charge = 0; charge < 2; charge++)
    {
      MEKey key = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,charge);
      mHelix_Kaon[key].clear();
      mMomentum[key].clear();
      mMass2[key].clear();
      mDca[key].clear();
      mNHitsFit[key].clear();
      mNSigmaKaon[key].clear();
    }
  }
  mEventCounter2[cent9][Bin_vz][Bin_Psi2] = 0;
}

void StTriFlowV0::size_phi(Int_t cent9, Int_t Bin_vz, Int_t Bin_Psi2)
{
  LOG_INFO << "Event Buffer: Centrality = " << cent9 << ", VertexZ = " << Bin_vz << ", Psi2 = " << Bin_Psi2 << endm;
  LOG_INFO << "Buffer_depth = " << mEventCounter2[cent9][Bin_vz][Bin_Psi2] << endm;

  LOG_INFO << "Size of primaryVertex = " << mPrimaryvertex[cent9][Bin_vz][Bin_Psi2].size() << endm;;
  LOG_INFO << "Size of refMult       = " << mRefMult[cent9][Bin_vz][Bin_Psi2].size() << endm;;
  LOG_INFO << "---------------------------------------------------------------------------" << endm;

  for(Int_t Bin_Event = 0; Bin_Event < mEventCounter2[cent9][Bin_vz][Bin_Psi2]; Bin_Event++)
  {
    MEKey key = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,0);
    LOG_INFO << "Event Number " << Bin_Event << ":" << endm; 
    LOG_INFO << "Positive Particle:" << endm;
    LOG_INFO << "  Size of Helix_Kplus  = " << mHelix_Kaon[key].size() << endm;;
    LOG_INFO << "  Size of Momentum     = " << mMomentum[key].size() << endm;
    LOG_INFO << "  Size of Mass2        = " << mMass2[key].size() << endm;
    LOG_INFO << "  Size of dca          = " << mDca[key].size() << endm;
    LOG_INFO << "  Size of nHitsFit     = " << mNHitsFit[key].size() << endm;
    LOG_INFO << "  Size of nSigmaKaon   = " << mNSigmaKaon[key].size() << endm;

    key = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,1);
    LOG_INFO << "Negative Particle:" << endm;
    LOG_INFO << "  Size of Helix_Kminus = " << mHelix_Kaon[key].size() << endm;
    LOG_INFO << "  Size of Momentum     = " << mMomentum[key].size() << endm;
    LOG_INFO << "  Size of Mass2        = " << mMass2[key].size() << endm;
    LOG_INFO << "  Size of dca          = " << mDca[key].size() << endm;
    LOG_INFO << "  Size of nHitsFit     = " << mNHitsFit[key].size() << endm;
    LOG_INFO << "  Size of nSigmaKaon   = " << mNSigmaKaon[key].size() << endm;
    LOG_INFO << "---------------------------------------------------------------------------" << endm;
  }
}

void StTriFlowV0::clear_Lambda(Int_t cent9, Int_t Bin_vz, Int_t Bin_Psi2)
{
  mPrimaryvertex[cent9][Bin_vz][Bin_Psi2].clear();
  mRefMult[cent9][Bin_vz][Bin_Psi2].clear();
  mCentrality[cent9][Bin_vz][Bin_Psi2].clear();
  mRunId[cent9][Bin_vz][Bin_Psi2].clear();
  mN_prim[cent9][Bin_vz][Bin_Psi2].clear();
  mN_non_prim[cent9][Bin_vz][Bin_Psi2].clear();
  mN_Tof_match[cent9][Bin_vz][Bin_Psi2].clear();
  mZDCx[cent9][Bin_vz][Bin_Psi2].clear();
  mBBCx[cent9][Bin_vz][Bin_Psi2].clear();
  mVzVpd[cent9][Bin_vz][Bin_Psi2].clear();
  mField[cent9][Bin_vz][Bin_Psi2].clear();

  for(Int_t j = 0; j < TriFlow::EtaGap_total; j++)
  {
    mQ2East[cent9][Bin_vz][Bin_Psi2][j].clear();
    mQ2West[cent9][Bin_vz][Bin_Psi2][j].clear();
    mQ3East[cent9][Bin_vz][Bin_Psi2][j].clear();
    mQ3West[cent9][Bin_vz][Bin_Psi2][j].clear();
  }

  for(Int_t Bin_Event = 0; Bin_Event < TriFlow::Buffer_depth; Bin_Event++)
  {
    for(Int_t charge = 0; charge < 2; charge++)
    {
      if(charge == 0) // proton
      {
	MEKey key = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,charge);
	mHelix_Proton[key].clear();
	mMomentum[key].clear();
	mMass2[key].clear();
	mDca[key].clear();
	mNHitsFit[key].clear();
	mNSigmaProton[key].clear();
	mLPTrack[key].clear();
	mTofFlag[key].clear();
	mTofTime[key].clear();
	mTofBeta[key].clear();
	mTofHit[key].clear();
      }
      if(charge == 1) // pi-
      {
	MEKey key = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,charge);
	mHelix_Pion[key].clear();
	mMomentum[key].clear();
	mMass2[key].clear();
	mDca[key].clear();
	mNHitsFit[key].clear();
	mNSigmaPion[key].clear();
	mLPTrack[key].clear();
	mTofFlag[key].clear();
	mTofTime[key].clear();
	mTofBeta[key].clear();
	mTofHit[key].clear();
      }
    }
  }
  mEventCounter2[cent9][Bin_vz][Bin_Psi2] = 0;
}

void StTriFlowV0::clear_AntiLambda(Int_t cent9, Int_t Bin_vz, Int_t Bin_Psi2)
{
  mPrimaryvertex[cent9][Bin_vz][Bin_Psi2].clear();
  mRefMult[cent9][Bin_vz][Bin_Psi2].clear();
  mCentrality[cent9][Bin_vz][Bin_Psi2].clear();
  mRunId[cent9][Bin_vz][Bin_Psi2].clear();
  mN_prim[cent9][Bin_vz][Bin_Psi2].clear();
  mN_non_prim[cent9][Bin_vz][Bin_Psi2].clear();
  mN_Tof_match[cent9][Bin_vz][Bin_Psi2].clear();
  mZDCx[cent9][Bin_vz][Bin_Psi2].clear();
  mBBCx[cent9][Bin_vz][Bin_Psi2].clear();
  mVzVpd[cent9][Bin_vz][Bin_Psi2].clear();
  mField[cent9][Bin_vz][Bin_Psi2].clear();

  for(Int_t j = 0; j < TriFlow::EtaGap_total; j++)
  {
    mQ2East[cent9][Bin_vz][Bin_Psi2][j].clear();
    mQ2West[cent9][Bin_vz][Bin_Psi2][j].clear();
    mQ3East[cent9][Bin_vz][Bin_Psi2][j].clear();
    mQ3West[cent9][Bin_vz][Bin_Psi2][j].clear();
  }

  for(Int_t Bin_Event = 0; Bin_Event < TriFlow::Buffer_depth; Bin_Event++)
  {
    for(Int_t charge = 0; charge < 2; charge++)
    {
      if(charge == 1) // anti-proton
      {
	MEKey key = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,charge);
	mHelix_Proton[key].clear();
	mMomentum[key].clear();
	mMass2[key].clear();
	mDca[key].clear();
	mNHitsFit[key].clear();
	mNSigmaProton[key].clear();
	mLPTrack[key].clear();
	mTofFlag[key].clear();
	mTofTime[key].clear();
	mTofBeta[key].clear();
	mTofHit[key].clear();
      }
      if(charge == 0) // pi+
      {
	MEKey key = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,charge);
	mHelix_Pion[key].clear();
	mMomentum[key].clear();
	mMass2[key].clear();
	mDca[key].clear();
	mNHitsFit[key].clear();
	mNSigmaPion[key].clear();
	mLPTrack[key].clear();
	mTofFlag[key].clear();
	mTofTime[key].clear();
	mTofBeta[key].clear();
	mTofHit[key].clear();
      }
    }
  }
  mEventCounter2[cent9][Bin_vz][Bin_Psi2] = 0;
}

void StTriFlowV0::clear_K0S(Int_t cent9, Int_t Bin_vz, Int_t Bin_Psi2)
{
  mPrimaryvertex[cent9][Bin_vz][Bin_Psi2].clear();
  mRefMult[cent9][Bin_vz][Bin_Psi2].clear();
  mCentrality[cent9][Bin_vz][Bin_Psi2].clear();
  mRunId[cent9][Bin_vz][Bin_Psi2].clear();
  mN_prim[cent9][Bin_vz][Bin_Psi2].clear();
  mN_non_prim[cent9][Bin_vz][Bin_Psi2].clear();
  mN_Tof_match[cent9][Bin_vz][Bin_Psi2].clear();
  mZDCx[cent9][Bin_vz][Bin_Psi2].clear();
  mBBCx[cent9][Bin_vz][Bin_Psi2].clear();
  mVzVpd[cent9][Bin_vz][Bin_Psi2].clear();
  mField[cent9][Bin_vz][Bin_Psi2].clear();

  for(Int_t j = 0; j < TriFlow::EtaGap_total; j++)
  {
    mQ2East[cent9][Bin_vz][Bin_Psi2][j].clear();
    mQ2West[cent9][Bin_vz][Bin_Psi2][j].clear();
    mQ3East[cent9][Bin_vz][Bin_Psi2][j].clear();
    mQ3West[cent9][Bin_vz][Bin_Psi2][j].clear();
  }

  for(Int_t Bin_Event = 0; Bin_Event < TriFlow::Buffer_depth; Bin_Event++)
  {
    for(Int_t charge = 0; charge < 2; charge++)
    {
      MEKey key = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,charge);
      mHelix_Pion[key].clear();
      mMomentum[key].clear();
      mMass2[key].clear();
      mDca[key].clear();
      mNHitsFit[key].clear();
      mNSigmaProton[key].clear();
      mLPTrack[key].clear();
      mTofFlag[key].clear();
      mTofTime[key].clear();
      mTofBeta[key].clear();
      mTofHit[key].clear();
    }
  }
  mEventCounter2[cent9][Bin_vz][Bin_Psi2] = 0;
}
//------------------------------------------------------------------------------------------------------------------

void StTriFlowV0::doPhi(Int_t Flag_ME, Int_t cent9, Int_t Bin_vz, Int_t Bin_Psi2) // 0: Same Event, 1: Mix Event
{
  if(Flag_ME == 0) // same event
  {
    for(Int_t Bin_Event = 0; Bin_Event < mEventCounter2[cent9][Bin_vz][Bin_Psi2]; Bin_Event++)
    {
      // event header
      mXuPhiMesonEvent->clearTrackList();
      mXuPhiMesonEvent->setPrimaryVertex(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mXuPhiMesonEvent->setRunId(mRunId[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mXuPhiMesonEvent->setRefMult(mRefMult[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mXuPhiMesonEvent->setCentrality(mCentrality[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mXuPhiMesonEvent->setN_prim(mN_prim[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mXuPhiMesonEvent->setN_non_prim(mN_non_prim[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mXuPhiMesonEvent->setN_Tof_match(mN_Tof_match[cent9][Bin_vz][Bin_Psi2][Bin_Event]);

      for(Int_t j = 0; j < TriFlow::EtaGap_total; j++)
      {
	// QVector
	mXuPhiMesonEvent->setQ2East(mQ2East[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
	mXuPhiMesonEvent->setQ2West(mQ2West[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
	mXuPhiMesonEvent->setQ3East(mQ3East[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
	mXuPhiMesonEvent->setQ3West(mQ3West[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
	// Number of Tracks
	mXuPhiMesonEvent->setNumTrackEast(mNumTrackEast[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
	mXuPhiMesonEvent->setNumTrackWest(mNumTrackWest[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
      }

      mXuPhiMesonEvent->setZDCx(mZDCx[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mXuPhiMesonEvent->setBBCx(mBBCx[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mXuPhiMesonEvent->setVzVpd(mVzVpd[cent9][Bin_vz][Bin_Psi2][Bin_Event]);

      // start to select phi candidate in a event
      MEKey key_plus  = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,0);
      MEKey key_minus = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,1);

      TLorentzVector ltrackA, ltrackB;
      for(Int_t n_kplus = 0; n_kplus < mHelix_Kaon[key_plus].size(); n_kplus++) // first track loop over K+ candidates
      {
	StThreeVectorF p_vecA = mHelix_Kaon[key_plus][n_kplus].cat(mHelix_Kaon[key_plus][n_kplus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event]));  // primary momentum
	p_vecA *= mMomentum[key_plus][n_kplus];
	ltrackA.SetXYZM(p_vecA.x(),p_vecA.y(),p_vecA.z(),TriFlow::mMassKaon);

	for(Int_t n_kminus = 0; n_kminus < mHelix_Kaon[key_minus].size(); n_kminus++) // second track loop over K- candidates
	{
	  StThreeVectorF p_vecB = mHelix_Kaon[key_minus][n_kminus].cat(mHelix_Kaon[key_minus][n_kminus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event]));  // primary momentum
	  p_vecB *= mMomentum[key_minus][n_kminus];
	  ltrackB.SetXYZM(p_vecB.x(),p_vecB.y(),p_vecB.z(),TriFlow::mMassKaon);

	  TLorentzVector trackAB      = ltrackA+ltrackB;
	  Double_t InvMassAB          = trackAB.M();
	  Double_t pt = trackAB.Perp();

	  // fill phi candidate into mTree_Phi
	  if(InvMassAB > TriFlow::mMassKaon*2 && InvMassAB < 1.05) 
	  {
	    mXuPhiMesonTrack = mXuPhiMesonEvent->createTrack();
	    mXuPhiMesonTrack->setMass2A(mMass2[key_plus][n_kplus]); // K+
	    mXuPhiMesonTrack->setMass2B(mMass2[key_minus][n_kminus]); // K-
	    mXuPhiMesonTrack->setNSigKaonA(mNSigmaKaon[key_plus][n_kplus]); // K+
	    mXuPhiMesonTrack->setNSigKaonB(mNSigmaKaon[key_minus][n_kminus]); // K-
	    mXuPhiMesonTrack->setDcaA(mDca[key_plus][n_kplus]); // K+
	    mXuPhiMesonTrack->setDcaB(mDca[key_minus][n_kminus]); // K-
	    mXuPhiMesonTrack->setTrackA(ltrackA); // K+
	    mXuPhiMesonTrack->setTrackB(ltrackB); // K-
	    mXuPhiMesonTrack->setFlagA(Bin_Event); // K+
	    mXuPhiMesonTrack->setFlagB(Bin_Event); // K-
	  }

	  // Fill histogram with InvMassAB information
	  h_Mass2->Fill(pt,InvMassAB);
	}
      }
    }
    mTree_Phi->Fill();
  }

  if(Flag_ME == 1) // mixed event
  {
    for(Int_t Bin_Event_A = 0; Bin_Event_A < mEventCounter2[cent9][Bin_vz][Bin_Psi2]-1; Bin_Event_A++)
    {
      MEKey key_A_plus  = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event_A,0);
      MEKey key_A_minus = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event_A,1);
      for(Int_t Bin_Event_B = Bin_Event_A+1; Bin_Event_B < mEventCounter2[cent9][Bin_vz][Bin_Psi2]; Bin_Event_B++)
      {
	MEKey key_B_plus  = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event_B,0);
	MEKey key_B_minus = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event_B,1);

	if(Bin_Event_A == 0 && Bin_Event_B == 1)
	{
	  Int_t Bin_Event = Bin_Event_A;
	  // event header
	  mXuPhiMesonEvent->clearTrackList();
	  mXuPhiMesonEvent->setPrimaryVertex(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mXuPhiMesonEvent->setRunId(mRunId[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mXuPhiMesonEvent->setRefMult(mRefMult[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mXuPhiMesonEvent->setCentrality(mCentrality[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mXuPhiMesonEvent->setN_prim(mN_prim[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mXuPhiMesonEvent->setN_non_prim(mN_non_prim[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mXuPhiMesonEvent->setN_Tof_match(mN_Tof_match[cent9][Bin_vz][Bin_Psi2][Bin_Event]);

	  for(Int_t j = 0; j < TriFlow::EtaGap_total; j++)
	  {
	    // QVector
	    mXuPhiMesonEvent->setQ2East(mQ2East[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
	    mXuPhiMesonEvent->setQ2West(mQ2West[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
	    mXuPhiMesonEvent->setQ3East(mQ3East[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
	    mXuPhiMesonEvent->setQ3West(mQ3West[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);

	    // Number of Tracks
	    mXuPhiMesonEvent->setNumTrackEast(mNumTrackEast[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
	    mXuPhiMesonEvent->setNumTrackWest(mNumTrackWest[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
	  }

	  mXuPhiMesonEvent->setZDCx(mZDCx[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mXuPhiMesonEvent->setBBCx(mBBCx[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mXuPhiMesonEvent->setVzVpd(mVzVpd[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	}

	TLorentzVector ltrackA, ltrackB;

	// start to mix events
	// mix K+ candidates from A event with K- candidates from B event
	for(Int_t n_kplus = 0; n_kplus < mHelix_Kaon[key_A_plus].size(); n_kplus++) // first track loop over K+ candidates from event A
	{
	  StThreeVectorF p_vecA(mHelix_Kaon[key_A_plus][n_kplus].cat(mHelix_Kaon[key_A_plus][n_kplus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event_A])));  // primary momentum
	  p_vecA *= mMomentum[key_A_plus][n_kplus];
	  ltrackA.SetXYZM(p_vecA.x(),p_vecA.y(),p_vecA.z(),TriFlow::mMassKaon); // K+

	  for(Int_t n_kminus = 0; n_kminus < mHelix_Kaon[key_B_minus].size(); n_kminus++) // second track loop over K- candidates from event B
	  {
	    StThreeVectorF p_vecB(mHelix_Kaon[key_B_minus][n_kminus].cat(mHelix_Kaon[key_B_minus][n_kminus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event_B])));  // primary momentum
	    p_vecB *= mMomentum[key_B_minus][n_kminus];
	    ltrackB.SetXYZM(p_vecB.x(),p_vecB.y(),p_vecB.z(),TriFlow::mMassKaon); // K-

	    TLorentzVector trackAB      = ltrackA+ltrackB;
	    Double_t InvMassAB          = trackAB.M();
	    Double_t pt = trackAB.Perp();

	    // fill phi candidate background into mTree_Phi
	    if(InvMassAB > TriFlow::mMassKaon*2 && InvMassAB < 1.05) 
	    {
	      mXuPhiMesonTrack = mXuPhiMesonEvent->createTrack();
	      mXuPhiMesonTrack->setMass2A(mMass2[key_A_plus][n_kplus]); // K+
	      mXuPhiMesonTrack->setMass2B(mMass2[key_B_minus][n_kminus]); // K-
	      mXuPhiMesonTrack->setNSigKaonA(mNSigmaKaon[key_A_plus][n_kplus]); // K+
	      mXuPhiMesonTrack->setNSigKaonB(mNSigmaKaon[key_B_minus][n_kminus]); // K-
	      mXuPhiMesonTrack->setDcaA(mDca[key_A_plus][n_kplus]); // K+
	      mXuPhiMesonTrack->setDcaB(mDca[key_B_minus][n_kminus]); // K-
	      mXuPhiMesonTrack->setTrackA(ltrackA); // K+
	      mXuPhiMesonTrack->setTrackB(ltrackB); // K-
	      mXuPhiMesonTrack->setFlagA(Bin_Event_A); // K+
	      mXuPhiMesonTrack->setFlagB(Bin_Event_B); // K-
	    }

	    // Fill histogram with InvMassAB information
	    h_Mass2->Fill(pt,InvMassAB);
	  }
	}

	// mix K- candidates from A event with K+ candidates from B event
	for(Int_t n_kminus = 0; n_kminus < mHelix_Kaon[key_A_minus].size(); n_kminus++) // first track loop over K- candidates from event A
	{
	  StThreeVectorF p_vecA(mHelix_Kaon[key_A_minus][n_kminus].cat(mHelix_Kaon[key_A_minus][n_kminus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event_A])));  // primary momentum
	  p_vecA *= mMomentum[key_A_minus][n_kminus];
	  ltrackA.SetXYZM(p_vecA.x(),p_vecA.y(),p_vecA.z(),TriFlow::mMassKaon); // K-

	  for(Int_t n_kplus = 0; n_kplus < mHelix_Kaon[key_B_plus].size(); n_kplus++) // second track loop over K+ candidates from event B
	  {
	    StThreeVectorF p_vecB(mHelix_Kaon[key_B_plus][n_kplus].cat(mHelix_Kaon[key_B_plus][n_kplus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event_B])));  // primary momentum
	    p_vecB *= mMomentum[key_B_plus][n_kplus];
	    ltrackB.SetXYZM(p_vecB.x(),p_vecB.y(),p_vecB.z(),TriFlow::mMassKaon); // K+

	    TLorentzVector trackAB      = ltrackA+ltrackB;
	    Double_t InvMassAB          = trackAB.M();
	    Double_t pt = trackAB.Perp();

	    // fill phi candidate background into mTree_Phi
	    if(InvMassAB > TriFlow::mMassKaon*2 && InvMassAB < 1.05) 
	    {
	      mXuPhiMesonTrack = mXuPhiMesonEvent->createTrack();
	      mXuPhiMesonTrack->setMass2A(mMass2[key_B_plus][n_kplus]); // K+
	      mXuPhiMesonTrack->setMass2B(mMass2[key_A_minus][n_kminus]); // K-
	      mXuPhiMesonTrack->setNSigKaonA(mNSigmaKaon[key_B_plus][n_kplus]); // K+
	      mXuPhiMesonTrack->setNSigKaonB(mNSigmaKaon[key_A_minus][n_kminus]); // K-
	      mXuPhiMesonTrack->setDcaA(mDca[key_B_plus][n_kplus]); // K+
	      mXuPhiMesonTrack->setDcaB(mDca[key_A_minus][n_kminus]); // K-
	      mXuPhiMesonTrack->setTrackA(ltrackB); // K+
	      mXuPhiMesonTrack->setTrackB(ltrackA); // K-
	      mXuPhiMesonTrack->setFlagA(Bin_Event_B); // K+
	      mXuPhiMesonTrack->setFlagB(Bin_Event_A); // K-
	    }

	    // Fill histogram with InvMassAB information
	    h_Mass2->Fill(pt,InvMassAB);
	  }
	}
      }
    }
    mTree_Phi->Fill();
  }
}

//------------------------------------------------------------------------------------------------------------------

void StTriFlowV0::doLambda(Int_t Flag_ME, Int_t cent9, Int_t Bin_vz, Int_t Bin_Psi2) // 0: Same Event, 1: Mix Event
{
  if(Flag_ME == 0) // same event
  {
    StTriFlow2ndVertexFinder *VertexFinder = new StTriFlow2ndVertexFinder();

    StPhysicalHelixD helixA, helixB;
    TLorentzVector ltrackA, ltrackB, ltrackC;
    Float_t EventVertexXA,EventVertexYA,EventVertexZA;
    StThreeVectorF vectorprim, vectorAB;
    Float_t MomentumA, MomentumB;
    Float_t dcaA, dcaB, dcaAB;
    Float_t VerdistX, VerdistY;

    for(Int_t Bin_Event = 0; Bin_Event < mEventCounter2[cent9][Bin_vz][Bin_Psi2]; Bin_Event++)
    {
      MEKey key_proton   = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,0); // p
      MEKey key_pi_minus = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,1); // pi-

      EventVertexXA = mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event].x();
      EventVertexYA = mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event].y();
      EventVertexZA = mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event].z();
      vectorprim.set(EventVertexXA,EventVertexYA,EventVertexZA);

      // event header
      mV0Event->clearTrackList();
      mV0Event->setPrimaryVertex(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mV0Event->setRunId(mRunId[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mV0Event->setRefMult(mRefMult[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mV0Event->setCentrality(mCentrality[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mV0Event->setN_prim(mN_prim[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mV0Event->setN_non_prim(mN_non_prim[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mV0Event->setN_Tof_match(mN_Tof_match[cent9][Bin_vz][Bin_Psi2][Bin_Event]);

      for(Int_t j = 0; j < TriFlow::EtaGap_total; j++)
      {
	// QVector
	mV0Event->setQ2East(mQ2East[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
	mV0Event->setQ2West(mQ2West[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
	mV0Event->setQ3East(mQ3East[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
	mV0Event->setQ3West(mQ3West[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
	// Number of Tracks
	mV0Event->setNumTrackEast(mNumTrackEast[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
	mV0Event->setNumTrackWest(mNumTrackWest[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
      }

      mV0Event->setZDCx(mZDCx[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mV0Event->setBBCx(mBBCx[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mV0Event->setVzVpd(mVzVpd[cent9][Bin_vz][Bin_Psi2][Bin_Event]);

      // start to select Lambda candidate in a event
      for(Int_t n_proton = 0; n_proton < mHelix_Proton[key_proton].size(); n_proton++) // first track loop over proton candidates
      {
        helixA = mHelix_Proton[key_proton][n_proton]; // Proton global track
	MomentumA = mMomentum[key_proton][n_proton];
	dcaA = mDca[key_proton][n_proton];

	for(Int_t n_pi_minus = 0; n_pi_minus < mHelix_Pion[key_pi_minus].size(); n_pi_minus++) // second track loop over pion minus candidates
	{
	  helixB = mHelix_Pion[key_pi_minus][n_pi_minus]; // Pion Minus global track
	  MomentumB = mMomentum[key_pi_minus][n_pi_minus];
	  dcaB = mDca[key_pi_minus][n_pi_minus];

	  VertexFinder->Find2ndVertex(helixA,helixB,vectorprim,MomentumA,MomentumB,ltrackA,ltrackB,VerdistX,VerdistY,vectorAB,dcaAB,0); // Lambda mode

	  // Invariant mass calculations
	  TLorentzVector trackAB = ltrackA+ltrackB; // mother particle
	  Double_t InvMassAB     = trackAB.M(); // invariant mass of mother particle

	  //-----------------------------------------------------------------------------
	  // get K0s by misidentification
	  ltrackC.SetXYZM(ltrackA.X(),ltrackA.Y(),ltrackA.Z(),TriFlow::mMassPion); // set Lorentz vector for pion from Event A

	  // Invariant mass calculations K0s
	  TLorentzVector trackCB      = ltrackC+ltrackB; // mother particle
	  Double_t InvMassCB          = trackCB.M(); // invariant mass of mother particle
	  //-----------------------------------------------------------------------------

	  if(fabs(dcaA) > mDca_proton && fabs(dcaB) > mDca_pion && fabs(dcaAB) < mDcaAB && VerdistX > mDecayLength && VerdistY < mDcaV0 && InvMassAB > mInvLambda_low && InvMassAB < mInvLambda_high)
	  {
	    Float_t Mass2_proton = mMass2[key_proton][n_proton];
	    Float_t Mass2_pion   = mMass2[key_pi_minus][n_pi_minus];

	    // final mass2 cut
	    Float_t Mass2_low_proton = 0.7;
	    Float_t Mass2_up_proton  = 1.1;
	    Float_t Mass2_low_pi_minus;
	    Float_t Mass2_up_pi_minus;
	    if(mMomentum[key_pi_minus][n_pi_minus] < 0.75)
	    {
	      Mass2_low_pi_minus = -0.013;
	      Mass2_up_pi_minus  =  0.047;
	    }
	    if(mMomentum[key_pi_minus][n_pi_minus] >= 0.75)
	    {
	      Mass2_low_pi_minus =  0.053 - 0.088*mMomentum[key_pi_minus][n_pi_minus];
	      Mass2_up_pi_minus  = -0.004 + 0.068*mMomentum[key_pi_minus][n_pi_minus];
	    }

	    StLorentzVectorD SttrackAB(trackAB.X(),trackAB.Y(),trackAB.Z(),trackAB.E());

	    // StV0TofCorrection for proton
	    if(mTofFlag[key_proton][n_proton] > 0 && mTofTime[key_proton][n_proton] != 0)
	    {
	      mTofCorr->setVectors3D(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event])(vectorAB)(mTofHit[key_proton][n_proton]);
	      mTofCorr->setMotherTracks(SttrackAB);
	      Float_t Time_before = mTofTime[key_proton][n_proton];
	      Float_t Beta_before = mTofBeta[key_proton][n_proton];
	      Float_t Time_after = Time_before;
	      Float_t Beta_after = Beta_before;
	      mTofCorr->correctBeta(helixA,Time_after,Beta_after);
	      if(Beta_after != 0) Mass2_proton = mMomentum[key_proton][n_proton]*mMomentum[key_proton][n_proton]*(1.0/(Beta_after*Beta_after) - 1.0);
	      mTofCorr->clearContainers();
	    }
	    // StV0TofCorrection for pion_minus
	    if(mTofFlag[key_pi_minus][n_pi_minus] > 0 && mTofTime[key_pi_minus][n_pi_minus] != 0)
	    {
	      mTofCorr->setVectors3D(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event])(vectorAB)(mTofHit[key_pi_minus][n_pi_minus]);
	      mTofCorr->setMotherTracks(SttrackAB);
	      Float_t Time_before = mTofTime[key_pi_minus][n_pi_minus];
	      Float_t Beta_before = mTofBeta[key_pi_minus][n_pi_minus];
	      Float_t Time_after = Time_before;
	      Float_t Beta_after = Beta_before;
	      mTofCorr->correctBeta(helixB,Time_after,Beta_after);
	      if(Beta_after != 0) Mass2_pion = mMomentum[key_pi_minus][n_pi_minus]*mMomentum[key_pi_minus][n_pi_minus]*(1.0/(Beta_after*Beta_after) - 1.0);
	      mTofCorr->clearContainers();
	    }

	    if(
		 ((Mass2_proton > Mass2_low_proton && Mass2_proton < Mass2_up_proton) || Mass2_proton < -10.0) // proton mass2 cut
	      && ((Mass2_pion > Mass2_low_pi_minus && Mass2_pion < Mass2_up_pi_minus) || Mass2_pion < -10.0) // pion_minus mass2 cut
	      )
	    {
	      h_Mass2_K0s->Fill(trackCB.Perp(),InvMassCB);
	      h_Mass2->Fill(trackAB.Perp(),InvMassAB);

	      // fill Lambda candidate into mTree_Lambda | the InvMass cut already done
	      mV0Track = mV0Event->createTrack();
	      mV0Track->setMass2A(Mass2_proton); // proton
	      mV0Track->setMass2B(Mass2_pion); // pi_minus
	      mV0Track->setNSigA(mNSigmaProton[key_proton][n_proton]); // proton
	      mV0Track->setNSigB(mNSigmaPion[key_pi_minus][n_pi_minus]); // pi_minus
	      mV0Track->setDcaA(mDca[key_proton][n_proton]); // proton
	      mV0Track->setDcaB(mDca[key_pi_minus][n_pi_minus]); // pi_minus
	      mV0Track->setGTrackA(ltrackA); // proton
	      mV0Track->setGTrackB(ltrackB); // pi_minus
	      mV0Track->setPTrackA(mLPTrack[key_proton][n_proton]); // proton
	      mV0Track->setPTrackB(mLPTrack[key_pi_minus][n_pi_minus]); // pi_minus
	      mV0Track->setFlagA(Bin_Event); // proton
	      mV0Track->setFlagB(Bin_Event); // pi_minus
	      mV0Track->setDcaAB(dcaAB);
	      mV0Track->setDecayLength(VerdistX);
	      mV0Track->setDcaV0(VerdistY);

	      if(!(InvMassCB > 0.498-3*0.005 && InvMassCB < 0.498+3*0.005))
	      {
		h_Mass2_sub->Fill(trackAB.Perp(),InvMassAB);
	      }
	    }
	  }
	}
      }
    }
    mTree_Lambda->Fill();
  }

  if(Flag_ME == 1) // mixed event
  {
    StTriFlow2ndVertexFinder *VertexFinder = new StTriFlow2ndVertexFinder();

    StPhysicalHelixD helixA, helixB;
    TLorentzVector ltrackA, ltrackB, ltrackC;
    StThreeVectorF vectorprim, vectorprimB, vectordiff, vectorAB;
    Float_t EventVertexXA, EventVertexYA, EventVertexZA, EventVertexXB, EventVertexYB, EventVertexZB, vertexAB_dist;
    Float_t MomentumA, MomentumB;
    Float_t dcaA, dcaB, dcaAB;
    Float_t VerdistX, VerdistY;

    for(Int_t Bin_Event_A = 0; Bin_Event_A < mEventCounter2[cent9][Bin_vz][Bin_Psi2]-1; Bin_Event_A++)
    {
      MEKey key_A_proton   = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event_A,0);
      MEKey key_A_pi_minus = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event_A,1);

      EventVertexXA = mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event_A].x();
      EventVertexYA = mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event_A].y();
      EventVertexZA = mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event_A].z();
      vectorprim.set(EventVertexXA,EventVertexYA,EventVertexZA);
      for(Int_t Bin_Event_B = Bin_Event_A+1; Bin_Event_B < mEventCounter2[cent9][Bin_vz][Bin_Psi2]; Bin_Event_B++)
      {
	MEKey key_B_proton   = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event_B,0);
	MEKey key_B_pi_minus = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event_B,1);

	EventVertexXB = mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event_B].x();
	EventVertexYB = mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event_B].y();
	EventVertexZB = mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event_B].z();
	vectorprimB.set(EventVertexXB,EventVertexYB,EventVertexZB);

	vectordiff = (vectorprim - vectorprimB); // difference between VertexA and VertexB
        vertexAB_dist  = vectordiff.mag(); // distance between eventA and eventB vertex

	if(Bin_Event_A == 0 && Bin_Event_B == 1)
	{
	  Int_t Bin_Event = Bin_Event_A;
	  // event header
	  mV0Event->clearTrackList();
	  mV0Event->setPrimaryVertex(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mV0Event->setRunId(mRunId[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mV0Event->setRefMult(mRefMult[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mV0Event->setCentrality(mCentrality[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mV0Event->setN_prim(mN_prim[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mV0Event->setN_non_prim(mN_non_prim[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mV0Event->setN_Tof_match(mN_Tof_match[cent9][Bin_vz][Bin_Psi2][Bin_Event]);

	  for(Int_t j = 0; j < TriFlow::EtaGap_total; j++)
	  {
	    // QVector
	    mV0Event->setQ2East(mQ2East[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
	    mV0Event->setQ2West(mQ2West[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
	    mV0Event->setQ3East(mQ3East[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
	    mV0Event->setQ3West(mQ3West[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);

	    // Number of Tracks
	    mV0Event->setNumTrackEast(mNumTrackEast[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
	    mV0Event->setNumTrackWest(mNumTrackWest[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
	  }

	  mV0Event->setZDCx(mZDCx[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mV0Event->setBBCx(mBBCx[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mV0Event->setVzVpd(mVzVpd[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	}

	// start to mix events
	// mix proton candidates from A event with pi_minus candidates from B event
	for(Int_t n_proton = 0; n_proton < mHelix_Proton[key_A_proton].size(); n_proton++) // first track loop over proton candidates from event A
	{
	  helixA = mHelix_Proton[key_A_proton][n_proton]; // proton
	  MomentumA = mMomentum[key_A_proton][n_proton];
	  dcaA = mDca[key_A_proton][n_proton];

	  for(Int_t n_pi_minus = 0; n_pi_minus < mHelix_Pion[key_B_pi_minus].size(); n_pi_minus++) // second track loop over pi_minus candidates from event B
	  {
	    Float_t BField = mField[cent9][Bin_vz][Bin_Psi2][Bin_Event_B];
	    Float_t Charge = mHelix_Pion[key_B_pi_minus][n_pi_minus].charge(BField);
	    StThreeVector<double> Momentum = mHelix_Pion[key_B_pi_minus][n_pi_minus].momentum(BField);
	    StThreeVector<double> Origin   = mHelix_Pion[key_B_pi_minus][n_pi_minus].origin();
	    helixB = StPhysicalHelixD(Momentum,Origin+vectordiff,BField,Charge); // pi_minus | redefine the helix of pi_minus after move the primary vertex from eventB to eventA
	    MomentumB = mMomentum[key_B_pi_minus][n_pi_minus];
	    dcaB = mDca[key_B_pi_minus][n_pi_minus];

	    VertexFinder->Find2ndVertex(helixA,helixB,vectorprim,MomentumA,MomentumB,ltrackA,ltrackB,VerdistX,VerdistY,vectorAB,dcaAB,0); // Lambda mode

	    // Invariant mass calculations
	    TLorentzVector trackAB = ltrackA+ltrackB; // mother particle
	    Double_t InvMassAB     = trackAB.M(); // invariant mass of mother particle

	    //-----------------------------------------------------------------------------
	    // get K0s by misidentification
	    TLorentzVector  ltrackC;
	    ltrackC.SetXYZM(ltrackA.X(),ltrackA.Y(),ltrackA.Z(),TriFlow::mMassPion); // set Lorentz vector for pion from Event A

	    // Invariant mass calculations K0s
	    TLorentzVector trackCB      = ltrackC+ltrackB; // mother particle
	    Double_t InvMassCB          = trackCB.M(); // invariant mass of mother particle
	    //-----------------------------------------------------------------------------

	    if(fabs(dcaA) > mDca_proton && fabs(dcaB) > mDca_pion && fabs(dcaAB) < mDcaAB && VerdistX > mDecayLength && VerdistY < mDcaV0 && InvMassAB > mInvLambda_low && InvMassAB < mInvLambda_high)
	    {
	      Float_t Mass2_proton = mMass2[key_A_proton][n_proton];
	      Float_t Mass2_pion   = mMass2[key_B_pi_minus][n_pi_minus];

	      // final mass2 cut
	      Float_t Mass2_low_proton = 0.7;
	      Float_t Mass2_up_proton  = 1.1;
	      Float_t Mass2_low_pi_minus;
	      Float_t Mass2_up_pi_minus;
	      if(mMomentum[key_B_pi_minus][n_pi_minus] < 0.75)
	      {
		Mass2_low_pi_minus = -0.013;
		Mass2_up_pi_minus  =  0.047;
	      }
	      if(mMomentum[key_B_pi_minus][n_pi_minus] >= 0.75)
	      {
		Mass2_low_pi_minus =  0.053 - 0.088*mMomentum[key_B_pi_minus][n_pi_minus];
		Mass2_up_pi_minus  = -0.004 + 0.068*mMomentum[key_B_pi_minus][n_pi_minus];
	      }

	      StLorentzVectorD SttrackAB(trackAB.X(),trackAB.Y(),trackAB.Z(),trackAB.E());

	      // StV0TofCorrection for proton from Event A
	      if(mTofFlag[key_A_proton][n_proton] > 0 && mTofTime[key_A_proton][n_proton] != 0)
	      {
		mTofCorr->setVectors3D(vectorprim)(vectorAB)(mTofHit[key_A_proton][n_proton]);
		mTofCorr->setMotherTracks(SttrackAB);
		Float_t Time_before = mTofTime[key_A_proton][n_proton];
		Float_t Beta_before = mTofBeta[key_A_proton][n_proton];
		Float_t Time_after = Time_before;
		Float_t Beta_after = Beta_before;
		mTofCorr->correctBeta(helixA,Time_after,Beta_after);
		if(Beta_after != 0) Mass2_proton = mMomentum[key_A_proton][n_proton]*mMomentum[key_A_proton][n_proton]*(1.0/(Beta_after*Beta_after) - 1.0);
		mTofCorr->clearContainers();
	      }
	      // StV0TofCorrection for pion_minus from Event B
	      if(mTofFlag[key_B_pi_minus][n_pi_minus] > 0 && mTofTime[key_B_pi_minus][n_pi_minus] != 0)
	      {
		mTofCorr->setVectors3D(vectorprim)(vectorAB)(mTofHit[key_B_pi_minus][n_pi_minus]+vectordiff); 
		mTofCorr->setMotherTracks(SttrackAB);
		Float_t Time_before = mTofTime[key_B_pi_minus][n_pi_minus];
		Float_t Beta_before = mTofBeta[key_B_pi_minus][n_pi_minus];
		Float_t Time_after = Time_before;
		Float_t Beta_after = Beta_before;
		mTofCorr->correctBeta(helixB,Time_after,Beta_after);
		if(Beta_after != 0) Mass2_pion = mMomentum[key_B_pi_minus][n_pi_minus]*mMomentum[key_B_pi_minus][n_pi_minus]*(1.0/(Beta_after*Beta_after) - 1.0);
		mTofCorr->clearContainers();
	      }

	      if(
		   ((Mass2_proton > Mass2_low_proton && Mass2_proton < Mass2_up_proton) || Mass2_proton < -10.0) // proton mass2 cut
		&& ((Mass2_pion > Mass2_low_pi_minus && Mass2_pion < Mass2_up_pi_minus) || Mass2_pion < -10.0) // pion_minus mass2 cut
		)
	      {
		h_Mass2_K0s->Fill(trackCB.Perp(),InvMassCB);
		h_Mass2->Fill(trackAB.Perp(),InvMassAB);

		// fill Lambda candidate into mTree_Lambda
		mV0Track = mV0Event->createTrack();
		mV0Track->setMass2A(mMass2[key_A_proton][n_proton]); // proton
		mV0Track->setMass2B(mMass2[key_B_pi_minus][n_pi_minus]); // pi_minus
		mV0Track->setNSigA(mNSigmaProton[key_A_proton][n_proton]); // proton
		mV0Track->setNSigB(mNSigmaPion[key_B_pi_minus][n_pi_minus]); // pi_minus
		mV0Track->setDcaA(mDca[key_A_proton][n_proton]); // proton
		mV0Track->setDcaB(mDca[key_B_pi_minus][n_pi_minus]); // pi_minus
		mV0Track->setGTrackA(ltrackA); // proton
		mV0Track->setGTrackB(ltrackB); // pi_minus
		mV0Track->setPTrackA(mLPTrack[key_A_proton][n_proton]); // proton
		mV0Track->setPTrackB(mLPTrack[key_B_pi_minus][n_pi_minus]); // pi_minus
		mV0Track->setFlagA(Bin_Event_A); // proton
		mV0Track->setFlagB(Bin_Event_B); // pi_minus
		mV0Track->setDcaAB(dcaAB);
		mV0Track->setDecayLength(VerdistX);
		mV0Track->setDcaV0(VerdistY);

		if(!(InvMassCB > 0.498-3*0.005 && InvMassCB < 0.498+3*0.005))
		{
		  h_Mass2_sub->Fill(trackAB.Perp(),InvMassAB);
		}
	      }
	    }
	  }
	}

	// mix pi_minus candidates from A event with proton candidates from B event
	for(Int_t n_pi_minus = 0; n_pi_minus < mHelix_Pion[key_A_pi_minus].size(); n_pi_minus++) // first track loop over pi_minus candidates from event A
	{
	  helixA = mHelix_Pion[key_A_pi_minus][n_pi_minus]; // pi_minus
	  MomentumA = mMomentum[key_A_pi_minus][n_pi_minus];
	  dcaA = mDca[key_A_pi_minus][n_pi_minus];

	  for(Int_t n_proton = 0; n_proton < mHelix_Proton[key_B_proton].size(); n_proton++) // second track loop over proton candidates from event B
	  {
	    Float_t BField = mField[cent9][Bin_vz][Bin_Psi2][Bin_Event_B];
	    Float_t Charge = mHelix_Proton[key_B_proton][n_proton].charge(BField);
	    StThreeVector<double> Momentum = mHelix_Proton[key_B_proton][n_proton].momentum(BField);
	    StThreeVector<double> Origin   = mHelix_Proton[key_B_proton][n_proton].origin();
	    helixB = StPhysicalHelixD(Momentum,Origin+vectordiff,BField,Charge); // proton | redefine the helix of pi_minus after move the primary vertex from eventB to eventA
	    MomentumB = mMomentum[key_B_proton][n_proton];
	    dcaB = mDca[key_B_proton][n_proton];

	    VertexFinder->Find2ndVertex(helixB,helixA,vectorprim,MomentumB,MomentumA,ltrackB,ltrackA,VerdistX,VerdistY,vectorAB,dcaAB,0); // Lambda mode

	    // Invariant mass calculations
	    TLorentzVector trackAB = ltrackA+ltrackB; // mother particle
	    Double_t InvMassAB     = trackAB.M(); // invariant mass of mother particle

	    //-----------------------------------------------------------------------------
	    // get K0s by misidentification
	    TLorentzVector  ltrackC;
	    ltrackC.SetXYZM(ltrackB.X(),ltrackB.Y(),ltrackB.Z(),TriFlow::mMassPion); // set Lorentz vector for pion from Event B

	    // Invariant mass calculations K0s
	    TLorentzVector trackCA      = ltrackC+ltrackA; // mother particle
	    Double_t InvMassCA          = trackCA.M(); // invariant mass of mother particle
	    //-----------------------------------------------------------------------------

	    if(fabs(dcaA) > mDca_pion && fabs(dcaB) > mDca_proton && fabs(dcaAB) < mDcaAB && VerdistX > mDecayLength && VerdistY < mDcaV0 && InvMassAB > mInvLambda_low && InvMassAB < mInvLambda_high)
	    {
	      Float_t Mass2_proton = mMass2[key_B_proton][n_proton];
	      Float_t Mass2_pion   = mMass2[key_A_pi_minus][n_pi_minus];

	      // final mass2 cut
	      Float_t Mass2_low_proton = 0.7;
	      Float_t Mass2_up_proton  = 1.1;
	      Float_t Mass2_low_pi_minus;
	      Float_t Mass2_up_pi_minus;
	      if(mMomentum[key_A_pi_minus][n_pi_minus] < 0.75)
	      {
		Mass2_low_pi_minus = -0.013;
		Mass2_up_pi_minus  =  0.047;
	      }
	      if(mMomentum[key_A_pi_minus][n_pi_minus] >= 0.75)
	      {
		Mass2_low_pi_minus =  0.053 - 0.088*mMomentum[key_A_pi_minus][n_pi_minus];
		Mass2_up_pi_minus  = -0.004 + 0.068*mMomentum[key_A_pi_minus][n_pi_minus];
	      }

	      StLorentzVectorD SttrackAB(trackAB.X(),trackAB.Y(),trackAB.Z(),trackAB.E());

	      // StV0TofCorrection for proton from Event B
	      if(mTofFlag[key_B_proton][n_proton] > 0 && mTofTime[key_B_proton][n_proton] != 0)
	      {
		mTofCorr->setVectors3D(vectorprim)(vectorAB)(mTofHit[key_B_proton][n_proton]+vectordiff);
		mTofCorr->setMotherTracks(SttrackAB);
		Float_t Time_before = mTofTime[key_B_proton][n_proton];
		Float_t Beta_before = mTofBeta[key_B_proton][n_proton];
		Float_t Time_after = Time_before;
		Float_t Beta_after = Beta_before;
		mTofCorr->correctBeta(helixB,Time_after,Beta_after);
		if(Beta_after != 0) Mass2_proton = mMomentum[key_B_proton][n_proton]*mMomentum[key_B_proton][n_proton]*(1.0/(Beta_after*Beta_after) - 1.0);
		mTofCorr->clearContainers();
	      }
	      // StV0TofCorrection for pion_minus from Event A
	      if(mTofFlag[key_A_pi_minus][n_pi_minus] > 0 && mTofTime[key_A_pi_minus][n_pi_minus] != 0)
	      {
		mTofCorr->setVectors3D(vectorprim)(vectorAB)(mTofHit[key_A_pi_minus][n_pi_minus]); 
		mTofCorr->setMotherTracks(SttrackAB);
		Float_t Time_before = mTofTime[key_A_pi_minus][n_pi_minus];
		Float_t Beta_before = mTofBeta[key_A_pi_minus][n_pi_minus];
		Float_t Time_after = Time_before;
		Float_t Beta_after = Beta_before;
		mTofCorr->correctBeta(helixA,Time_after,Beta_after);
		if(Beta_after != 0) Mass2_pion = mMomentum[key_A_pi_minus][n_pi_minus]*mMomentum[key_A_pi_minus][n_pi_minus]*(1.0/(Beta_after*Beta_after) - 1.0);
		mTofCorr->clearContainers();
	      }

	      if(
		   ((Mass2_proton > Mass2_low_proton && Mass2_proton < Mass2_up_proton) || Mass2_proton < -10.0) // proton mass2 cut
		&& ((Mass2_pion > Mass2_low_pi_minus && Mass2_pion < Mass2_up_pi_minus) || Mass2_pion < -10.0) // pion_minus mass2 cut
		)
	      {
		h_Mass2_K0s->Fill(trackCA.Perp(),InvMassCA);
		h_Mass2->Fill(trackAB.Perp(),InvMassAB);

		// fill Lambda candidate into mTree_Phi
		mV0Track = mV0Event->createTrack();
		mV0Track->setMass2A(mMass2[key_B_proton][n_proton]); // proton
		mV0Track->setMass2B(mMass2[key_A_pi_minus][n_pi_minus]); // pi_minus
		mV0Track->setNSigA(mNSigmaProton[key_B_proton][n_proton]); // proton
		mV0Track->setNSigB(mNSigmaPion[key_A_pi_minus][n_pi_minus]); // pi_minus
		mV0Track->setDcaA(mDca[key_B_proton][n_proton]); // proton
		mV0Track->setDcaB(mDca[key_A_pi_minus][n_pi_minus]); // pi_minus
		mV0Track->setGTrackA(ltrackB); // proton
		mV0Track->setGTrackB(ltrackA); // pi_minus
		mV0Track->setPTrackA(mLPTrack[key_B_proton][n_proton]); // proton
		mV0Track->setPTrackB(mLPTrack[key_A_pi_minus][n_pi_minus]); // pi_minus
		mV0Track->setFlagA(Bin_Event_B); // proton
		mV0Track->setFlagB(Bin_Event_A); // pi_minus
		mV0Track->setDcaAB(dcaAB);
		mV0Track->setDecayLength(VerdistX);
		mV0Track->setDcaV0(VerdistY);

		if(!(InvMassCA > 0.498-3*0.005 && InvMassCA < 0.498+3*0.005))
		{
		  h_Mass2_sub->Fill(trackAB.Perp(),InvMassAB);
		}
	      }
	    }
	  }
	}
      }
    }
    mTree_Lambda->Fill();
  }
}

//------------------------------------------------------------------------------------------------------------------

void StTriFlowV0::doAntiLambda(Int_t Flag_ME, Int_t cent9, Int_t Bin_vz, Int_t Bin_Psi2) // 0: Same Event, 1: Mix Event
{
  if(Flag_ME == 0) // same event
  {
    StTriFlow2ndVertexFinder *VertexFinder = new StTriFlow2ndVertexFinder();

    StPhysicalHelixD helixA, helixB;
    TLorentzVector ltrackA, ltrackB, ltrackC;
    Float_t EventVertexXA,EventVertexYA,EventVertexZA;
    StThreeVectorF vectorprim, vectorAB;
    Float_t MomentumA, MomentumB;
    Float_t dcaA, dcaB, dcaAB;
    Float_t VerdistX, VerdistY;

    for(Int_t Bin_Event = 0; Bin_Event < mEventCounter2[cent9][Bin_vz][Bin_Psi2]; Bin_Event++)
    {
      MEKey key_antiproton   = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,1); // anti-proton 
      MEKey key_pi_plus = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,0); // pi+

      EventVertexXA = mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event].x();
      EventVertexYA = mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event].y();
      EventVertexZA = mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event].z();
      vectorprim.set(EventVertexXA,EventVertexYA,EventVertexZA);

      // event header
      mV0Event->clearTrackList();
      mV0Event->setPrimaryVertex(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mV0Event->setRunId(mRunId[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mV0Event->setRefMult(mRefMult[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mV0Event->setCentrality(mCentrality[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mV0Event->setN_prim(mN_prim[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mV0Event->setN_non_prim(mN_non_prim[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mV0Event->setN_Tof_match(mN_Tof_match[cent9][Bin_vz][Bin_Psi2][Bin_Event]);

      for(Int_t j = 0; j < TriFlow::EtaGap_total; j++)
      {
	// QVector
	mV0Event->setQ2East(mQ2East[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
	mV0Event->setQ2West(mQ2West[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
	mV0Event->setQ3East(mQ3East[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
	mV0Event->setQ3West(mQ3West[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
	// Number of Tracks
	mV0Event->setNumTrackEast(mNumTrackEast[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
	mV0Event->setNumTrackWest(mNumTrackWest[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
      }

      mV0Event->setZDCx(mZDCx[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mV0Event->setBBCx(mBBCx[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mV0Event->setVzVpd(mVzVpd[cent9][Bin_vz][Bin_Psi2][Bin_Event]);

      // start to select antiLambda candidate in a event
      for(Int_t n_antiproton = 0; n_antiproton < mHelix_Proton[key_antiproton].size(); n_antiproton++) // first track loop over antiproton candidates
      {
        helixA = mHelix_Proton[key_antiproton][n_antiproton]; // anti-Proton global track
	MomentumA = mMomentum[key_antiproton][n_antiproton];
	dcaA = mDca[key_antiproton][n_antiproton];

	for(Int_t n_pi_plus = 0; n_pi_plus < mHelix_Pion[key_pi_plus].size(); n_pi_plus++) // second track loop over pion plus candidates
	{
	  helixB = mHelix_Pion[key_pi_plus][n_pi_plus]; // Pion Plus global track
	  MomentumB = mMomentum[key_pi_plus][n_pi_plus];
	  dcaB = mDca[key_pi_plus][n_pi_plus];

	  VertexFinder->Find2ndVertex(helixA,helixB,vectorprim,MomentumA,MomentumB,ltrackA,ltrackB,VerdistX,VerdistY,vectorAB,dcaAB,1); // antiLambda mode

	  // Invariant mass calculations
	  TLorentzVector trackAB = ltrackA+ltrackB; // mother particle
	  Double_t InvMassAB     = trackAB.M(); // invariant mass of mother particle

	  //-----------------------------------------------------------------------------
	  // get K0s by misidentification
	  ltrackC.SetXYZM(ltrackA.X(),ltrackA.Y(),ltrackA.Z(),TriFlow::mMassPion); // set Lorentz vector for pion from Event A

	  // Invariant mass calculations K0s
	  TLorentzVector trackCB      = ltrackC+ltrackB; // mother particle
	  Double_t InvMassCB          = trackCB.M(); // invariant mass of mother particle
	  //-----------------------------------------------------------------------------

	  if(fabs(dcaA) > mDca_proton && fabs(dcaB) > mDca_pion && fabs(dcaAB) < mDcaAB && VerdistX > mDecayLength && VerdistY < mDcaV0 && InvMassAB > mInvLambda_low && InvMassAB < mInvLambda_high)
	  {
	    Float_t Mass2_proton = mMass2[key_antiproton][n_antiproton];
	    Float_t Mass2_pion   = mMass2[key_pi_plus][n_pi_plus];

	    // final mass2 cut
	    Float_t Mass2_low_proton = 0.7;
	    Float_t Mass2_up_proton  = 1.1;
	    Float_t Mass2_low_pi_minus;
	    Float_t Mass2_up_pi_minus;
	    if(mMomentum[key_pi_plus][n_pi_plus] < 0.75)
	    {
	      Mass2_low_pi_minus = -0.013;
	      Mass2_up_pi_minus  =  0.047;
	    }
	    if(mMomentum[key_pi_plus][n_pi_plus] >= 0.75)
	    {
	      Mass2_low_pi_minus =  0.053 - 0.088*mMomentum[key_pi_plus][n_pi_plus];
	      Mass2_up_pi_minus  = -0.004 + 0.068*mMomentum[key_pi_plus][n_pi_plus];
	    }

	    StLorentzVectorD SttrackAB(trackAB.X(),trackAB.Y(),trackAB.Z(),trackAB.E());

	    // StV0TofCorrection for anti-proton
	    if(mTofFlag[key_antiproton][n_antiproton] > 0 && mTofTime[key_antiproton][n_antiproton] != 0)
	    {
	      mTofCorr->setVectors3D(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event])(vectorAB)(mTofHit[key_antiproton][n_antiproton]);
	      mTofCorr->setMotherTracks(SttrackAB);
	      Float_t Time_before = mTofTime[key_antiproton][n_antiproton];
	      Float_t Beta_before = mTofBeta[key_antiproton][n_antiproton];
	      Float_t Time_after = Time_before;
	      Float_t Beta_after = Beta_before;
	      mTofCorr->correctBeta(helixA,Time_after,Beta_after);
	      if(Beta_after != 0) Mass2_proton = mMomentum[key_antiproton][n_antiproton]*mMomentum[key_antiproton][n_antiproton]*(1.0/(Beta_after*Beta_after) - 1.0);
	      mTofCorr->clearContainers();
	    }
	    // StV0TofCorrection for pion_plus
	    if(mTofFlag[key_pi_plus][n_pi_plus] > 0 && mTofTime[key_pi_plus][n_pi_plus] != 0)
	    {
	      mTofCorr->setVectors3D(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event])(vectorAB)(mTofHit[key_pi_plus][n_pi_plus]);
	      mTofCorr->setMotherTracks(SttrackAB);
	      Float_t Time_before = mTofTime[key_pi_plus][n_pi_plus];
	      Float_t Beta_before = mTofBeta[key_pi_plus][n_pi_plus];
	      Float_t Time_after = Time_before;
	      Float_t Beta_after = Beta_before;
	      mTofCorr->correctBeta(helixB,Time_after,Beta_after);
	      if(Beta_after != 0) Mass2_pion = mMomentum[key_pi_plus][n_pi_plus]*mMomentum[key_pi_plus][n_pi_plus]*(1.0/(Beta_after*Beta_after) - 1.0);
	      mTofCorr->clearContainers();
	    }

	    if(
		 ((Mass2_proton > Mass2_low_proton && Mass2_proton < Mass2_up_proton) || Mass2_proton < -10.0) // proton mass2 cut
	      && ((Mass2_pion > Mass2_low_pi_minus && Mass2_pion < Mass2_up_pi_minus) || Mass2_pion < -10.0) // pion_minus mass2 cut
	      )
	    {
	      h_Mass2_K0s->Fill(trackCB.Perp(),InvMassCB);
	      h_Mass2->Fill(trackAB.Perp(),InvMassAB);

	      // fill Lambda candidate into mTree_AntiLambda | the InvMass cut already done
	      mV0Track = mV0Event->createTrack();
	      mV0Track->setMass2A(Mass2_proton); // antiproton
	      mV0Track->setMass2B(Mass2_pion); // pi_plus
	      mV0Track->setNSigA(mNSigmaProton[key_antiproton][n_antiproton]); // antiproton
	      mV0Track->setNSigB(mNSigmaPion[key_pi_plus][n_pi_plus]); // pi_plus
	      mV0Track->setDcaA(mDca[key_antiproton][n_antiproton]); // antiproton
	      mV0Track->setDcaB(mDca[key_pi_plus][n_pi_plus]); // pi_plus
	      mV0Track->setGTrackA(ltrackA); // antiproton
	      mV0Track->setGTrackB(ltrackB); // pi_plus
	      mV0Track->setPTrackA(mLPTrack[key_antiproton][n_antiproton]); // antiproton
	      mV0Track->setPTrackB(mLPTrack[key_pi_plus][n_pi_plus]); // pi_plus
	      mV0Track->setFlagA(Bin_Event); // antiproton
	      mV0Track->setFlagB(Bin_Event); // pi_plus
	      mV0Track->setDcaAB(dcaAB);
	      mV0Track->setDecayLength(VerdistX);
	      mV0Track->setDcaV0(VerdistY);

	      if(!(InvMassCB > 0.498-3*0.005 && InvMassCB < 0.498+3*0.005))
	      {
		h_Mass2_sub->Fill(trackAB.Perp(),InvMassAB);
	      }
	    }
	  }
	}
      }
    }
    mTree_AntiLambda->Fill();
  }

  if(Flag_ME == 1) // mixed event
  {
    StTriFlow2ndVertexFinder *VertexFinder = new StTriFlow2ndVertexFinder();

    StPhysicalHelixD helixA, helixB;
    TLorentzVector ltrackA, ltrackB, ltrackC;
    StThreeVectorF vectorprim, vectorprimB, vectordiff, vectorAB;
    Float_t EventVertexXA, EventVertexYA, EventVertexZA, EventVertexXB, EventVertexYB, EventVertexZB, vertexAB_dist;
    Float_t MomentumA, MomentumB;
    Float_t dcaA, dcaB, dcaAB;
    Float_t VerdistX, VerdistY;

    for(Int_t Bin_Event_A = 0; Bin_Event_A < mEventCounter2[cent9][Bin_vz][Bin_Psi2]-1; Bin_Event_A++)
    {
      MEKey key_A_antiproton   = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event_A,1); // anti-proton
      MEKey key_A_pi_plus = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event_A,0); // pi+

      EventVertexXA = mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event_A].x();
      EventVertexYA = mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event_A].y();
      EventVertexZA = mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event_A].z();
      vectorprim.set(EventVertexXA,EventVertexYA,EventVertexZA);
      for(Int_t Bin_Event_B = Bin_Event_A+1; Bin_Event_B < mEventCounter2[cent9][Bin_vz][Bin_Psi2]; Bin_Event_B++)
      {
	MEKey key_B_antiproton   = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event_B,1); // anti-proton
	MEKey key_B_pi_plus = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event_B,0); // pi+

	EventVertexXB = mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event_B].x();
	EventVertexYB = mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event_B].y();
	EventVertexZB = mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event_B].z();
	vectorprimB.set(EventVertexXB,EventVertexYB,EventVertexZB);

	vectordiff = (vectorprim - vectorprimB); // difference between VertexA and VertexB
        vertexAB_dist  = vectordiff.mag(); // distance between eventA and eventB vertex

	if(Bin_Event_A == 0 && Bin_Event_B == 1)
	{
	  Int_t Bin_Event = Bin_Event_A;
	  // event header
	  mV0Event->clearTrackList();
	  mV0Event->setPrimaryVertex(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mV0Event->setRunId(mRunId[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mV0Event->setRefMult(mRefMult[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mV0Event->setCentrality(mCentrality[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mV0Event->setN_prim(mN_prim[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mV0Event->setN_non_prim(mN_non_prim[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mV0Event->setN_Tof_match(mN_Tof_match[cent9][Bin_vz][Bin_Psi2][Bin_Event]);

	  for(Int_t j = 0; j < TriFlow::EtaGap_total; j++)
	  {
	    // QVector
	    mV0Event->setQ2East(mQ2East[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
	    mV0Event->setQ2West(mQ2West[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
	    mV0Event->setQ3East(mQ3East[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
	    mV0Event->setQ3West(mQ3West[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);

	    // Number of Tracks
	    mV0Event->setNumTrackEast(mNumTrackEast[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
	    mV0Event->setNumTrackWest(mNumTrackWest[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
	  }

	  mV0Event->setZDCx(mZDCx[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mV0Event->setBBCx(mBBCx[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mV0Event->setVzVpd(mVzVpd[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	}

	// start to mix events
	// mix anti-proton candidates from A event with pi_plus candidates from B event
	for(Int_t n_antiproton = 0; n_antiproton < mHelix_Proton[key_A_antiproton].size(); n_antiproton++) // first track loop over anti-proton candidates from event A
	{
	  helixA = mHelix_Proton[key_A_antiproton][n_antiproton]; // anti-proton
	  MomentumA = mMomentum[key_A_antiproton][n_antiproton];
	  dcaA = mDca[key_A_antiproton][n_antiproton];

	  for(Int_t n_pi_plus = 0; n_pi_plus < mHelix_Pion[key_B_pi_plus].size(); n_pi_plus++) // second track loop over pi_plus candidates from event B
	  {
	    Float_t BField = mField[cent9][Bin_vz][Bin_Psi2][Bin_Event_B];
	    Float_t Charge = mHelix_Pion[key_B_pi_plus][n_pi_plus].charge(BField);
	    StThreeVector<double> Momentum = mHelix_Pion[key_B_pi_plus][n_pi_plus].momentum(BField);
	    StThreeVector<double> Origin   = mHelix_Pion[key_B_pi_plus][n_pi_plus].origin();
	    helixB = StPhysicalHelixD(Momentum,Origin+vectordiff,BField,Charge); // pi_plus
	    MomentumB = mMomentum[key_B_pi_plus][n_pi_plus];
	    dcaB = mDca[key_B_pi_plus][n_pi_plus];

	    VertexFinder->Find2ndVertex(helixA,helixB,vectorprim,MomentumA,MomentumB,ltrackA,ltrackB,VerdistX,VerdistY,vectorAB,dcaAB,1); // antiLambda mode

	    // Invariant mass calculations
	    TLorentzVector trackAB = ltrackA+ltrackB; // mother particle
	    Double_t InvMassAB     = trackAB.M(); // invariant mass of mother particle

	    //-----------------------------------------------------------------------------
	    // get K0s by misidentification
	    TLorentzVector  ltrackC;
	    ltrackC.SetXYZM(ltrackA.X(),ltrackA.Y(),ltrackA.Z(),TriFlow::mMassPion); // set Lorentz vector for pion from Event A

	    // Invariant mass calculations K0s
	    TLorentzVector trackCB      = ltrackC+ltrackB; // mother particle
	    Double_t InvMassCB          = trackCB.M(); // invariant mass of mother particle
	    //-----------------------------------------------------------------------------

	    if(fabs(dcaA) > mDca_proton && fabs(dcaB) > mDca_pion && fabs(dcaAB) < mDcaAB && VerdistX > mDecayLength && VerdistY < mDcaV0 && InvMassAB > mInvLambda_low && InvMassAB < mInvLambda_high)
	    {
	      Float_t Mass2_proton = mMass2[key_A_antiproton][n_antiproton];
	      Float_t Mass2_pion   = mMass2[key_B_pi_plus][n_pi_plus];

	      // final mass2 cut
	      Float_t Mass2_low_proton = 0.7;
	      Float_t Mass2_up_proton  = 1.1;
	      Float_t Mass2_low_pi_minus;
	      Float_t Mass2_up_pi_minus;
	      if(mMomentum[key_B_pi_plus][n_pi_plus] < 0.75)
	      {
		Mass2_low_pi_minus = -0.013;
		Mass2_up_pi_minus  =  0.047;
	      }
	      if(mMomentum[key_B_pi_plus][n_pi_plus] >= 0.75)
	      {
		Mass2_low_pi_minus =  0.053 - 0.088*mMomentum[key_B_pi_plus][n_pi_plus];
		Mass2_up_pi_minus  = -0.004 + 0.068*mMomentum[key_B_pi_plus][n_pi_plus];
	      }

	      StLorentzVectorD SttrackAB(trackAB.X(),trackAB.Y(),trackAB.Z(),trackAB.E());

	      // StV0TofCorrection for proton from Event A
	      if(mTofFlag[key_A_antiproton][n_antiproton] > 0 && mTofTime[key_A_antiproton][n_antiproton] != 0)
	      {
		mTofCorr->setVectors3D(vectorprim)(vectorAB)(mTofHit[key_A_antiproton][n_antiproton]);
		mTofCorr->setMotherTracks(SttrackAB);
		Float_t Time_before = mTofTime[key_A_antiproton][n_antiproton];
		Float_t Beta_before = mTofBeta[key_A_antiproton][n_antiproton];
		Float_t Time_after = Time_before;
		Float_t Beta_after = Beta_before;
		mTofCorr->correctBeta(helixA,Time_after,Beta_after);
		if(Beta_after != 0) Mass2_proton = mMomentum[key_A_antiproton][n_antiproton]*mMomentum[key_A_antiproton][n_antiproton]*(1.0/(Beta_after*Beta_after) - 1.0);
		mTofCorr->clearContainers();
	      }
	      // StV0TofCorrection for pion_minus from Event B
	      if(mTofFlag[key_B_pi_plus][n_pi_plus] > 0 && mTofTime[key_B_pi_plus][n_pi_plus] != 0)
	      {
		mTofCorr->setVectors3D(vectorprim)(vectorAB)(mTofHit[key_B_pi_plus][n_pi_plus]+vectordiff); 
		mTofCorr->setMotherTracks(SttrackAB);
		Float_t Time_before = mTofTime[key_B_pi_plus][n_pi_plus];
		Float_t Beta_before = mTofBeta[key_B_pi_plus][n_pi_plus];
		Float_t Time_after = Time_before;
		Float_t Beta_after = Beta_before;
		mTofCorr->correctBeta(helixB,Time_after,Beta_after);
		if(Beta_after != 0) Mass2_pion = mMomentum[key_B_pi_plus][n_pi_plus]*mMomentum[key_B_pi_plus][n_pi_plus]*(1.0/(Beta_after*Beta_after) - 1.0);
		mTofCorr->clearContainers();
	      }

	      if(
		   ((Mass2_proton > Mass2_low_proton && Mass2_proton < Mass2_up_proton) || Mass2_proton < -10.0) // proton mass2 cut
		&& ((Mass2_pion > Mass2_low_pi_minus && Mass2_pion < Mass2_up_pi_minus) || Mass2_pion < -10.0) // pion_minus mass2 cut
		)
	      {
		h_Mass2_K0s->Fill(trackCB.Perp(),InvMassCB);
		h_Mass2->Fill(trackAB.Perp(),InvMassAB);

		// fill Lambda candidate into mTree_Phi
		mV0Track = mV0Event->createTrack();
		mV0Track->setMass2A(mMass2[key_A_antiproton][n_antiproton]); // antiproton
		mV0Track->setMass2B(mMass2[key_B_pi_plus][n_pi_plus]); // pi_plus
		mV0Track->setNSigA(mNSigmaProton[key_A_antiproton][n_antiproton]); // antiproton
		mV0Track->setNSigB(mNSigmaPion[key_B_pi_plus][n_pi_plus]); // pi_plus
		mV0Track->setDcaA(mDca[key_A_antiproton][n_antiproton]); // antiproton
		mV0Track->setDcaB(mDca[key_B_pi_plus][n_pi_plus]); // pi_plus
		mV0Track->setGTrackA(ltrackA); // antiproton
		mV0Track->setGTrackB(ltrackB); // pi_plus
		mV0Track->setPTrackA(mLPTrack[key_A_antiproton][n_antiproton]); // antiproton
		mV0Track->setPTrackB(mLPTrack[key_B_pi_plus][n_pi_plus]); // pi_plus
		mV0Track->setFlagA(Bin_Event_A); // antiproton
		mV0Track->setFlagB(Bin_Event_B); // pi_plus
		mV0Track->setDcaAB(dcaAB);
		mV0Track->setDecayLength(VerdistX);
		mV0Track->setDcaV0(VerdistY);

		if(!(InvMassCB > 0.498-3*0.005 && InvMassCB < 0.498+3*0.005))
		{
		  h_Mass2_sub->Fill(trackAB.Perp(),InvMassAB);
		}
	      }
	    }
	  }
	}

	// mix pi_plus candidates from A event with anti-proton candidates from B event
	for(Int_t n_pi_plus = 0; n_pi_plus < mHelix_Pion[key_A_pi_plus].size(); n_pi_plus++) // first track loop over pi_plus candidates from event A
	{
	  helixA = mHelix_Pion[key_A_pi_plus][n_pi_plus]; // pi_plus
	  MomentumA = mMomentum[key_A_pi_plus][n_pi_plus];
	  dcaA = mDca[key_A_pi_plus][n_pi_plus];

	  for(Int_t n_antiproton = 0; n_antiproton < mHelix_Proton[key_B_antiproton].size(); n_antiproton++) // second track loop over anti-proton candidates from event B
	  {
	    Float_t BField = mField[cent9][Bin_vz][Bin_Psi2][Bin_Event_B];
	    Float_t Charge = mHelix_Proton[key_B_antiproton][n_antiproton].charge(BField);
	    StThreeVector<double> Momentum = mHelix_Proton[key_B_antiproton][n_antiproton].momentum(BField);
	    StThreeVector<double> Origin   = mHelix_Proton[key_B_antiproton][n_antiproton].origin();
	    helixB = StPhysicalHelixD(Momentum,Origin+vectordiff,BField,Charge); // anti-proton
	    MomentumB = mMomentum[key_B_antiproton][n_antiproton];
	    dcaB = mDca[key_B_antiproton][n_antiproton];

	    VertexFinder->Find2ndVertex(helixB,helixA,vectorprim,MomentumB,MomentumA,ltrackB,ltrackA,VerdistX,VerdistY,vectorAB,dcaAB,1); // Lambda mode

	    // Invariant mass calculations
	    TLorentzVector trackAB = ltrackA+ltrackB; // mother particle
	    Double_t InvMassAB     = trackAB.M(); // invariant mass of mother particle

	    //-----------------------------------------------------------------------------
	    // get K0s by misidentification
	    TLorentzVector  ltrackC;
	    ltrackC.SetXYZM(ltrackB.X(),ltrackB.Y(),ltrackB.Z(),TriFlow::mMassPion); // set Lorentz vector for pion from Event B

	    // Invariant mass calculations K0s
	    TLorentzVector trackCA      = ltrackC+ltrackA; // mother particle
	    Double_t InvMassCA          = trackCA.M(); // invariant mass of mother particle
	    //-----------------------------------------------------------------------------

	    if(fabs(dcaA) > mDca_pion && fabs(dcaB) > mDca_proton && fabs(dcaAB) < mDcaAB && VerdistX > mDecayLength && VerdistY < mDcaV0 && InvMassAB > mInvLambda_low && InvMassAB < mInvLambda_high)
	    {
	      Float_t Mass2_proton = mMass2[key_B_antiproton][n_antiproton];
	      Float_t Mass2_pion   = mMass2[key_A_pi_plus][n_pi_plus];

	      // final mass2 cut
	      Float_t Mass2_low_proton = 0.7;
	      Float_t Mass2_up_proton  = 1.1;
	      Float_t Mass2_low_pi_minus;
	      Float_t Mass2_up_pi_minus;
	      if(mMomentum[key_A_pi_plus][n_pi_plus] < 0.75)
	      {
		Mass2_low_pi_minus = -0.013;
		Mass2_up_pi_minus  =  0.047;
	      }
	      if(mMomentum[key_A_pi_plus][n_pi_plus] >= 0.75)
	      {
		Mass2_low_pi_minus =  0.053 - 0.088*mMomentum[key_A_pi_plus][n_pi_plus];
		Mass2_up_pi_minus  = -0.004 + 0.068*mMomentum[key_A_pi_plus][n_pi_plus];
	      }

	      StLorentzVectorD SttrackAB(trackAB.X(),trackAB.Y(),trackAB.Z(),trackAB.E());

	      // StV0TofCorrection for anti-proton from Event B
	      if(mTofFlag[key_B_antiproton][n_antiproton] > 0 && mTofTime[key_B_antiproton][n_antiproton] != 0)
	      {
		mTofCorr->setVectors3D(vectorprim)(vectorAB)(mTofHit[key_B_antiproton][n_antiproton]+vectordiff);
		mTofCorr->setMotherTracks(SttrackAB);
		Float_t Time_before = mTofTime[key_B_antiproton][n_antiproton];
		Float_t Beta_before = mTofBeta[key_B_antiproton][n_antiproton];
		Float_t Time_after = Time_before;
		Float_t Beta_after = Beta_before;
		mTofCorr->correctBeta(helixB,Time_after,Beta_after);
		if(Beta_after != 0) Mass2_proton = mMomentum[key_B_antiproton][n_antiproton]*mMomentum[key_B_antiproton][n_antiproton]*(1.0/(Beta_after*Beta_after) - 1.0);
		mTofCorr->clearContainers();
	      }
	      // StV0TofCorrection for pion_plus from Event A
	      if(mTofFlag[key_A_pi_plus][n_pi_plus] > 0 && mTofTime[key_A_pi_plus][n_pi_plus] != 0)
	      {
		mTofCorr->setVectors3D(vectorprim)(vectorAB)(mTofHit[key_A_pi_plus][n_pi_plus]); 
		mTofCorr->setMotherTracks(SttrackAB);
		Float_t Time_before = mTofTime[key_A_pi_plus][n_pi_plus];
		Float_t Beta_before = mTofBeta[key_A_pi_plus][n_pi_plus];
		Float_t Time_after = Time_before;
		Float_t Beta_after = Beta_before;
		mTofCorr->correctBeta(helixA,Time_after,Beta_after);
		if(Beta_after != 0) Mass2_pion = mMomentum[key_A_pi_plus][n_pi_plus]*mMomentum[key_A_pi_plus][n_pi_plus]*(1.0/(Beta_after*Beta_after) - 1.0);
		mTofCorr->clearContainers();
	      }

	      if(
		   ((Mass2_proton > Mass2_low_proton && Mass2_proton < Mass2_up_proton) || Mass2_proton < -10.0) // proton mass2 cut
		&& ((Mass2_pion > Mass2_low_pi_minus && Mass2_pion < Mass2_up_pi_minus) || Mass2_pion < -10.0) // pion_minus mass2 cut
		)
	      {
		h_Mass2_K0s->Fill(trackCA.Perp(),InvMassCA);
		h_Mass2->Fill(trackAB.Perp(),InvMassAB);

		// fill Lambda candidate into mTree_Phi
		mV0Track = mV0Event->createTrack();
		mV0Track->setMass2A(mMass2[key_B_antiproton][n_antiproton]); // anti-proton
		mV0Track->setMass2B(mMass2[key_A_pi_plus][n_pi_plus]); // pi_plus
		mV0Track->setNSigA(mNSigmaProton[key_B_antiproton][n_antiproton]); // anti-proton
		mV0Track->setNSigB(mNSigmaPion[key_A_pi_plus][n_pi_plus]); // pi_plus
		mV0Track->setDcaA(mDca[key_B_antiproton][n_antiproton]); // anti-proton
		mV0Track->setDcaB(mDca[key_A_pi_plus][n_pi_plus]); // pi_plus
		mV0Track->setGTrackA(ltrackB); // anti-proton
		mV0Track->setGTrackB(ltrackA); // pi_plus
		mV0Track->setPTrackA(mLPTrack[key_B_antiproton][n_antiproton]); // anti-proton
		mV0Track->setPTrackB(mLPTrack[key_A_pi_plus][n_pi_plus]); // pi_plus
		mV0Track->setFlagA(Bin_Event_B); // anti-proton
		mV0Track->setFlagB(Bin_Event_A); // pi_plus
		mV0Track->setDcaAB(dcaAB);
		mV0Track->setDecayLength(VerdistX);
		mV0Track->setDcaV0(VerdistY);

		if(!(InvMassCA > 0.498-3*0.005 && InvMassCA < 0.498+3*0.005))
		{
		  h_Mass2_sub->Fill(trackAB.Perp(),InvMassAB);
		}
	      }
	    }
	  }
	}
      }
    }
    mTree_AntiLambda->Fill();
  }
}

//------------------------------------------------------------------------------------------------------------------

void StTriFlowV0::doK0S(Int_t Flag_ME, Int_t cent9, Int_t Bin_vz, Int_t Bin_Psi2) // 0: Same Event, 1: Mix Event
{
  if(Flag_ME == 0) // same event
  {
    StTriFlow2ndVertexFinder *VertexFinder = new StTriFlow2ndVertexFinder();

    StPhysicalHelixD helixA, helixB;
    TLorentzVector ltrackA, ltrackB, ltrackC;
    Float_t EventVertexXA,EventVertexYA,EventVertexZA;
    StThreeVectorF vectorprim, vectorAB;
    Float_t MomentumA, MomentumB;
    Float_t dcaA, dcaB, dcaAB;
    Float_t VerdistX, VerdistY;

    for(Int_t Bin_Event = 0; Bin_Event < mEventCounter2[cent9][Bin_vz][Bin_Psi2]; Bin_Event++)
    {
      MEKey key_pi_plus  = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,0); // pi+
      MEKey key_pi_minus = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,1); // pi-

      EventVertexXA = mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event].x();
      EventVertexYA = mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event].y();
      EventVertexZA = mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event].z();
      vectorprim.set(EventVertexXA,EventVertexYA,EventVertexZA);

      // event header
      mV0Event->clearTrackList();
      mV0Event->setPrimaryVertex(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mV0Event->setRunId(mRunId[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mV0Event->setRefMult(mRefMult[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mV0Event->setCentrality(mCentrality[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mV0Event->setN_prim(mN_prim[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mV0Event->setN_non_prim(mN_non_prim[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mV0Event->setN_Tof_match(mN_Tof_match[cent9][Bin_vz][Bin_Psi2][Bin_Event]);

      for(Int_t j = 0; j < TriFlow::EtaGap_total; j++)
      {
	// QVector
	mV0Event->setQ2East(mQ2East[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
	mV0Event->setQ2West(mQ2West[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
	mV0Event->setQ3East(mQ3East[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
	mV0Event->setQ3West(mQ3West[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
	// Number of Tracks
	mV0Event->setNumTrackEast(mNumTrackEast[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
	mV0Event->setNumTrackWest(mNumTrackWest[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
      }

      mV0Event->setZDCx(mZDCx[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mV0Event->setBBCx(mBBCx[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
      mV0Event->setVzVpd(mVzVpd[cent9][Bin_vz][Bin_Psi2][Bin_Event]);

      // start to select K0S candidate in a event
      for(Int_t n_pi_plus = 0; n_pi_plus < mHelix_Pion[key_pi_plus].size(); n_pi_plus++) // first track loop over pion plus candidates
      {
        helixA = mHelix_Pion[key_pi_plus][n_pi_plus]; // Pion Plus global track
	MomentumA = mMomentum[key_pi_plus][n_pi_plus];
	dcaA = mDca[key_pi_plus][n_pi_plus];

	for(Int_t n_pi_minus = 0; n_pi_minus < mHelix_Pion[key_pi_minus].size(); n_pi_minus++) // second track loop over pion minus candidates
	{
	  helixB = mHelix_Pion[key_pi_minus][n_pi_minus]; // Pion Minus global track
	  MomentumB = mMomentum[key_pi_minus][n_pi_minus];
	  dcaB = mDca[key_pi_minus][n_pi_minus];

	  VertexFinder->Find2ndVertex(helixA,helixB,vectorprim,MomentumA,MomentumB,ltrackA,ltrackB,VerdistX,VerdistY,vectorAB,dcaAB,2); // K0S mode

	  // Invariant mass calculations
	  TLorentzVector trackAB = ltrackA+ltrackB; // mother particle
	  Double_t InvMassAB     = trackAB.M(); // invariant mass of mother particle

	  //-----------------------------------------------------------------------------
	  // get Lambda by misidentification
	  TLorentzVector  ltrackP, ltrackPbar;
	  ltrackP.SetXYZM(ltrackA.X(),ltrackA.Y(),ltrackA.Z(),TriFlow::mMassProton); // set Lorentz vector of pi+ to p
	  ltrackPbar.SetXYZM(ltrackB.X(),ltrackB.Y(),ltrackB.Z(),TriFlow::mMassProton); // set Lorentz vector of pi- to p

	  // Invariant mass calculations K0s
	  TLorentzVector trackLambda      = ltrackP+ltrackB; // p + pi-
	  TLorentzVector trackAntiLambda      = ltrackPbar+ltrackA; // pbar + pi+
	  Double_t InvMassLambda = trackLambda.M(); // invariant mass of Lambda
	  Double_t InvMassAntiLambda = trackAntiLambda.M(); // invariant mass of Lambda
	  //-----------------------------------------------------------------------------

	  if(fabs(dcaA) > mDca_pion && fabs(dcaB) > mDca_pion && fabs(dcaAB) < mDcaAB && VerdistX > mDecayLength && VerdistY < mDcaV0 && InvMassAB > mInvK0S_low && InvMassAB < mInvK0S_high)
	  {
	    Float_t Mass2_pi_plus  = mMass2[key_pi_plus][n_pi_plus];
	    Float_t Mass2_pi_minus = mMass2[key_pi_minus][n_pi_minus];

	    // final mass2 cut
	    Float_t Mass2_low_pi_plus;
	    Float_t Mass2_up_pi_plus;
	    if(mMomentum[key_pi_plus][n_pi_plus] < 0.75)
	    {
	      Mass2_low_pi_plus = -0.013;
	      Mass2_up_pi_plus  =  0.047;
	    }
	    if(mMomentum[key_pi_plus][n_pi_plus] >= 0.75)
	    {
	      Mass2_low_pi_plus =  0.053 - 0.088*mMomentum[key_pi_plus][n_pi_plus];
	      Mass2_up_pi_plus  = -0.004 + 0.068*mMomentum[key_pi_plus][n_pi_plus];
	    }

	    Float_t Mass2_low_pi_minus;
	    Float_t Mass2_up_pi_minus;
	    if(mMomentum[key_pi_minus][n_pi_minus] < 0.75)
	    {
	      Mass2_low_pi_minus = -0.013;
	      Mass2_up_pi_minus  =  0.047;
	    }
	    if(mMomentum[key_pi_minus][n_pi_minus] >= 0.75)
	    {
	      Mass2_low_pi_minus =  0.053 - 0.088*mMomentum[key_pi_minus][n_pi_minus];
	      Mass2_up_pi_minus  = -0.004 + 0.068*mMomentum[key_pi_minus][n_pi_minus];
	    }

	    StLorentzVectorD SttrackAB(trackAB.X(),trackAB.Y(),trackAB.Z(),trackAB.E());

	    // StV0TofCorrection for pion_plus
	    if(mTofFlag[key_pi_plus][n_pi_plus] > 0 && mTofTime[key_pi_plus][n_pi_plus] != 0)
	    {
	      mTofCorr->setVectors3D(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event])(vectorAB)(mTofHit[key_pi_plus][n_pi_plus]);
	      mTofCorr->setMotherTracks(SttrackAB);
	      Float_t Time_before = mTofTime[key_pi_plus][n_pi_plus];
	      Float_t Beta_before = mTofBeta[key_pi_plus][n_pi_plus];
	      Float_t Time_after = Time_before;
	      Float_t Beta_after = Beta_before;
	      mTofCorr->correctBeta(helixA,Time_after,Beta_after);
	      if(Beta_after != 0) Mass2_pi_plus = mMomentum[key_pi_plus][n_pi_plus]*mMomentum[key_pi_plus][n_pi_plus]*(1.0/(Beta_after*Beta_after) - 1.0);
	      mTofCorr->clearContainers();
	    }
	    // StV0TofCorrection for pion_minus
	    if(mTofFlag[key_pi_minus][n_pi_minus] > 0 && mTofTime[key_pi_minus][n_pi_minus] != 0)
	    {
	      mTofCorr->setVectors3D(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event])(vectorAB)(mTofHit[key_pi_minus][n_pi_minus]);
	      mTofCorr->setMotherTracks(SttrackAB);
	      Float_t Time_before = mTofTime[key_pi_minus][n_pi_minus];
	      Float_t Beta_before = mTofBeta[key_pi_minus][n_pi_minus];
	      Float_t Time_after = Time_before;
	      Float_t Beta_after = Beta_before;
	      mTofCorr->correctBeta(helixB,Time_after,Beta_after);
	      if(Beta_after != 0) Mass2_pi_minus = mMomentum[key_pi_minus][n_pi_minus]*mMomentum[key_pi_minus][n_pi_minus]*(1.0/(Beta_after*Beta_after) - 1.0);
	      mTofCorr->clearContainers();
	    }

	    if(
	         ((Mass2_pi_plus  > Mass2_low_pi_plus  && Mass2_pi_plus  < Mass2_up_pi_plus ) || Mass2_pi_plus  < -10.0) // pion_plus  mass2 cut
	      && ((Mass2_pi_minus > Mass2_low_pi_minus && Mass2_pi_minus < Mass2_up_pi_minus) || Mass2_pi_minus < -10.0) // pion_minus mass2 cut
	      )
	    {
	      h_Mass2->Fill(trackAB.Perp(),InvMassAB);
	      h_Mass2_Lambda->Fill(trackLambda.Perp(),InvMassLambda);
	      h_Mass2_AntiLambda->Fill(trackAntiLambda.Perp(),InvMassAntiLambda);
	      if(!((InvMassLambda > 1.1157-0.006*3 && InvMassLambda < 1.1157+0.006*3) || (InvMassAntiLambda > 1.1157-0.006*3 && InvMassAntiLambda < 1.1157+0.006*3)))
		h_Mass2_sub->Fill(trackAB.Perp(),InvMassAB);

	      // fill K0S candidate into mTree_K0S | the InvMass cut already done
	      mV0Track = mV0Event->createTrack();
	      mV0Track->setMass2A(Mass2_pi_plus); // pi_plus
	      mV0Track->setMass2B(Mass2_pi_minus); // pi_minus
	      mV0Track->setNSigA(mNSigmaPion[key_pi_plus][n_pi_plus]); // pi_plus
	      mV0Track->setNSigB(mNSigmaPion[key_pi_minus][n_pi_minus]); // pi_minus
	      mV0Track->setDcaA(mDca[key_pi_plus][n_pi_plus]); // pi_plus
	      mV0Track->setDcaB(mDca[key_pi_minus][n_pi_minus]); // pi_minus
	      mV0Track->setGTrackA(ltrackA); // pi_plus
	      mV0Track->setGTrackB(ltrackB); // pi_minus
	      mV0Track->setPTrackA(mLPTrack[key_pi_plus][n_pi_plus]); // pi_plus
	      mV0Track->setPTrackB(mLPTrack[key_pi_minus][n_pi_minus]); // pi_minus
	      mV0Track->setFlagA(Bin_Event); // pi_plus
	      mV0Track->setFlagB(Bin_Event); // pi_minus
	      mV0Track->setDcaAB(dcaAB);
	      mV0Track->setDecayLength(VerdistX);
	      mV0Track->setDcaV0(VerdistY);
	    }
	  }
	}
      }
    }
    mTree_K0S->Fill();
  }

  if(Flag_ME == 1) // mixed event
  {
    StTriFlow2ndVertexFinder *VertexFinder = new StTriFlow2ndVertexFinder();

    StPhysicalHelixD helixA, helixB;
    TLorentzVector ltrackA, ltrackB, ltrackC;
    StThreeVectorF vectorprim, vectorprimB, vectordiff, vectorAB;
    Float_t EventVertexXA, EventVertexYA, EventVertexZA, EventVertexXB, EventVertexYB, EventVertexZB, vertexAB_dist;
    Float_t MomentumA, MomentumB;
    Float_t dcaA, dcaB, dcaAB;
    Float_t VerdistX, VerdistY;

    for(Int_t Bin_Event_A = 0; Bin_Event_A < mEventCounter2[cent9][Bin_vz][Bin_Psi2]-1; Bin_Event_A++)
    {
      MEKey key_A_pi_plus  = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event_A,0);
      MEKey key_A_pi_minus = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event_A,1);

      EventVertexXA = mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event_A].x();
      EventVertexYA = mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event_A].y();
      EventVertexZA = mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event_A].z();
      vectorprim.set(EventVertexXA,EventVertexYA,EventVertexZA);
      for(Int_t Bin_Event_B = Bin_Event_A+1; Bin_Event_B < mEventCounter2[cent9][Bin_vz][Bin_Psi2]; Bin_Event_B++)
      {
	MEKey key_B_pi_plus  = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event_B,0);
	MEKey key_B_pi_minus = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event_B,1);

	EventVertexXB = mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event_B].x();
	EventVertexYB = mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event_B].y();
	EventVertexZB = mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event_B].z();
	vectorprimB.set(EventVertexXB,EventVertexYB,EventVertexZB);

	vectordiff = (vectorprim - vectorprimB); // difference between VertexA and VertexB
        vertexAB_dist  = vectordiff.mag(); // distance between eventA and eventB vertex

	if(Bin_Event_A == 0 && Bin_Event_B == 1)
	{
	  Int_t Bin_Event = Bin_Event_A;
	  // event header
	  mV0Event->clearTrackList();
	  mV0Event->setPrimaryVertex(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mV0Event->setRunId(mRunId[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mV0Event->setRefMult(mRefMult[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mV0Event->setCentrality(mCentrality[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mV0Event->setN_prim(mN_prim[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mV0Event->setN_non_prim(mN_non_prim[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mV0Event->setN_Tof_match(mN_Tof_match[cent9][Bin_vz][Bin_Psi2][Bin_Event]);

	  for(Int_t j = 0; j < TriFlow::EtaGap_total; j++)
	  {
	    // QVector
	    mV0Event->setQ2East(mQ2East[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
	    mV0Event->setQ2West(mQ2West[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
	    mV0Event->setQ3East(mQ3East[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
	    mV0Event->setQ3West(mQ3West[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);

	    // Number of Tracks
	    mV0Event->setNumTrackEast(mNumTrackEast[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
	    mV0Event->setNumTrackWest(mNumTrackWest[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
	  }

	  mV0Event->setZDCx(mZDCx[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mV0Event->setBBCx(mBBCx[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	  mV0Event->setVzVpd(mVzVpd[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
	}

	// start to mix events
	// mix pi_plus candidates from A event with pi_minus candidates from B event
	for(Int_t n_pi_plus = 0; n_pi_plus < mHelix_Pion[key_A_pi_plus].size(); n_pi_plus++) // first track loop over pi_plus candidates from event A
	{
	  helixA = mHelix_Pion[key_A_pi_plus][n_pi_plus]; // pi_plus
	  MomentumA = mMomentum[key_A_pi_plus][n_pi_plus];
	  dcaA = mDca[key_A_pi_plus][n_pi_plus];

	  for(Int_t n_pi_minus = 0; n_pi_minus < mHelix_Pion[key_B_pi_minus].size(); n_pi_minus++) // second track loop over pi_minus candidates from event B
	  {
	    Float_t BField = mField[cent9][Bin_vz][Bin_Psi2][Bin_Event_B];
	    Float_t Charge = mHelix_Pion[key_B_pi_minus][n_pi_minus].charge(BField);
	    StThreeVector<double> Momentum = mHelix_Pion[key_B_pi_minus][n_pi_minus].momentum(BField);
	    StThreeVector<double> Origin   = mHelix_Pion[key_B_pi_minus][n_pi_minus].origin();
	    helixB = StPhysicalHelixD(Momentum,Origin+vectordiff,BField,Charge); // pi_minus | redefine the helix of pi_minus after move the primary vertex from eventB to eventA
	    MomentumB = mMomentum[key_B_pi_minus][n_pi_minus];
	    dcaB = mDca[key_B_pi_minus][n_pi_minus];

	    VertexFinder->Find2ndVertex(helixA,helixB,vectorprim,MomentumA,MomentumB,ltrackA,ltrackB,VerdistX,VerdistY,vectorAB,dcaAB,2); // K0S mode

	    // Invariant mass calculations
	    TLorentzVector trackAB = ltrackA+ltrackB; // mother particle
	    Double_t InvMassAB     = trackAB.M(); // invariant mass of mother particle

	    //-----------------------------------------------------------------------------
	    // get Lambda and antiLambda by misidentification
	    TLorentzVector  ltrackP, ltrackPbar;
	    ltrackP.SetXYZM(ltrackA.X(),ltrackA.Y(),ltrackA.Z(),TriFlow::mMassProton); // set Lorentz vector of pi+ to p
	    ltrackPbar.SetXYZM(ltrackB.X(),ltrackB.Y(),ltrackB.Z(),TriFlow::mMassProton); // set Lorentz vector of pi- to p

	    // Invariant mass calculations Lambda and antiLambda
	    TLorentzVector trackLambda      = ltrackP+ltrackB; // p + pi-
	    TLorentzVector trackAntiLambda      = ltrackPbar+ltrackA; // pbar + pi+
	    Double_t InvMassLambda = trackLambda.M(); // invariant mass of Lambda
	    Double_t InvMassAntiLambda = trackAntiLambda.M(); // invariant mass of Lambda
	    //-----------------------------------------------------------------------------

	    if(fabs(dcaA) > mDca_pion && fabs(dcaB) > mDca_pion && fabs(dcaAB) < mDcaAB && VerdistX > mDecayLength && VerdistY < mDcaV0 && InvMassAB > mInvK0S_low && InvMassAB < mInvK0S_high)
	    {
	      Float_t Mass2_pi_plus  = mMass2[key_A_pi_plus][n_pi_plus];
	      Float_t Mass2_pi_minus = mMass2[key_B_pi_minus][n_pi_minus];

	      // final mass2 cut
	      Float_t Mass2_low_pi_plus;
	      Float_t Mass2_up_pi_plus;
	      if(mMomentum[key_A_pi_plus][n_pi_plus] < 0.75)
	      {
		Mass2_low_pi_plus = -0.013;
		Mass2_up_pi_plus  =  0.047;
	      }
	      if(mMomentum[key_A_pi_plus][n_pi_plus] >= 0.75)
	      {
		Mass2_low_pi_plus =  0.053 - 0.088*mMomentum[key_A_pi_plus][n_pi_plus];
		Mass2_up_pi_plus  = -0.004 + 0.068*mMomentum[key_A_pi_plus][n_pi_plus];
	      }

	      Float_t Mass2_low_pi_minus;
	      Float_t Mass2_up_pi_minus;
	      if(mMomentum[key_B_pi_minus][n_pi_minus] < 0.75)
	      {
		Mass2_low_pi_minus = -0.013;
		Mass2_up_pi_minus  =  0.047;
	      }
	      if(mMomentum[key_B_pi_minus][n_pi_minus] >= 0.75)
	      {
		Mass2_low_pi_minus =  0.053 - 0.088*mMomentum[key_B_pi_minus][n_pi_minus];
		Mass2_up_pi_minus  = -0.004 + 0.068*mMomentum[key_B_pi_minus][n_pi_minus];
	      }

	      StLorentzVectorD SttrackAB(trackAB.X(),trackAB.Y(),trackAB.Z(),trackAB.E());

	      // StV0TofCorrection for proton from Event A
	      if(mTofFlag[key_A_pi_plus][n_pi_plus] > 0 && mTofTime[key_A_pi_plus][n_pi_plus] != 0)
	      {
		mTofCorr->setVectors3D(vectorprim)(vectorAB)(mTofHit[key_A_pi_plus][n_pi_plus]);
		mTofCorr->setMotherTracks(SttrackAB);
		Float_t Time_before = mTofTime[key_A_pi_plus][n_pi_plus];
		Float_t Beta_before = mTofBeta[key_A_pi_plus][n_pi_plus];
		Float_t Time_after = Time_before;
		Float_t Beta_after = Beta_before;
		mTofCorr->correctBeta(helixA,Time_after,Beta_after);
		if(Beta_after != 0) Mass2_pi_plus = mMomentum[key_A_pi_plus][n_pi_plus]*mMomentum[key_A_pi_plus][n_pi_plus]*(1.0/(Beta_after*Beta_after) - 1.0);
		mTofCorr->clearContainers();
	      }
	      // StV0TofCorrection for pion_minus from Event B
	      if(mTofFlag[key_B_pi_minus][n_pi_minus] > 0 && mTofTime[key_B_pi_minus][n_pi_minus] != 0)
	      {
		mTofCorr->setVectors3D(vectorprim)(vectorAB)(mTofHit[key_B_pi_minus][n_pi_minus]+vectordiff); 
		mTofCorr->setMotherTracks(SttrackAB);
		Float_t Time_before = mTofTime[key_B_pi_minus][n_pi_minus];
		Float_t Beta_before = mTofBeta[key_B_pi_minus][n_pi_minus];
		Float_t Time_after = Time_before;
		Float_t Beta_after = Beta_before;
		mTofCorr->correctBeta(helixB,Time_after,Beta_after);
		if(Beta_after != 0) Mass2_pi_minus = mMomentum[key_B_pi_minus][n_pi_minus]*mMomentum[key_B_pi_minus][n_pi_minus]*(1.0/(Beta_after*Beta_after) - 1.0);
		mTofCorr->clearContainers();
	      }

	      if(
		   ((Mass2_pi_plus  > Mass2_low_pi_plus  && Mass2_pi_plus  < Mass2_up_pi_plus)  || Mass2_pi_plus  < -10.0) // pion_plus  mass2 cut
		&& ((Mass2_pi_minus > Mass2_low_pi_minus && Mass2_pi_minus < Mass2_up_pi_minus) || Mass2_pi_minus < -10.0) // pion_minus mass2 cut
		)
	      {
		h_Mass2->Fill(trackAB.Perp(),InvMassAB);
		h_Mass2_Lambda->Fill(trackLambda.Perp(),InvMassLambda);
		h_Mass2_AntiLambda->Fill(trackAntiLambda.Perp(),InvMassAntiLambda);
		if(!((InvMassLambda > 1.1157-0.006*3 && InvMassLambda < 1.1157+0.006*3) || (InvMassAntiLambda > 1.1157-0.006*3 && InvMassAntiLambda < 1.1157+0.006*3)))
		  h_Mass2_sub->Fill(trackAB.Perp(),InvMassAB);

		// fill Lambda candidate into mTree_Lambda
		mV0Track = mV0Event->createTrack();
		mV0Track->setMass2A(mMass2[key_A_pi_plus][n_pi_plus]); // pi_plus
		mV0Track->setMass2B(mMass2[key_B_pi_minus][n_pi_minus]); // pi_minus
		mV0Track->setNSigA(mNSigmaPion[key_A_pi_plus][n_pi_plus]); // pi_plus
		mV0Track->setNSigB(mNSigmaPion[key_B_pi_minus][n_pi_minus]); // pi_minus
		mV0Track->setDcaA(mDca[key_A_pi_plus][n_pi_plus]); // pi_plus
		mV0Track->setDcaB(mDca[key_B_pi_minus][n_pi_minus]); // pi_minus
		mV0Track->setGTrackA(ltrackA); // pi_plus
		mV0Track->setGTrackB(ltrackB); // pi_minus
		mV0Track->setPTrackA(mLPTrack[key_A_pi_plus][n_pi_plus]); // pi_plus
		mV0Track->setPTrackB(mLPTrack[key_B_pi_minus][n_pi_minus]); // pi_minus
		mV0Track->setFlagA(Bin_Event_A); // pi_plus
		mV0Track->setFlagB(Bin_Event_B); // pi_minus
		mV0Track->setDcaAB(dcaAB);
		mV0Track->setDecayLength(VerdistX);
		mV0Track->setDcaV0(VerdistY);
	      }
	    }
	  }
	}

	// mix pi_minus candidates from A event with pi_plus candidates from B event
	for(Int_t n_pi_minus = 0; n_pi_minus < mHelix_Pion[key_A_pi_minus].size(); n_pi_minus++) // first track loop over pi_minus candidates from event A
	{
	  helixA = mHelix_Pion[key_A_pi_minus][n_pi_minus]; // pi_minus
	  MomentumA = mMomentum[key_A_pi_minus][n_pi_minus];
	  dcaA = mDca[key_A_pi_minus][n_pi_minus];

	  for(Int_t n_pi_plus = 0; n_pi_plus < mHelix_Pion[key_B_pi_plus].size(); n_pi_plus++) // second track loop over proton candidates from event B
	  {
	    Float_t BField = mField[cent9][Bin_vz][Bin_Psi2][Bin_Event_B];
	    Float_t Charge = mHelix_Pion[key_B_pi_plus][n_pi_plus].charge(BField);
	    StThreeVector<double> Momentum = mHelix_Pion[key_B_pi_plus][n_pi_plus].momentum(BField);
	    StThreeVector<double> Origin   = mHelix_Pion[key_B_pi_plus][n_pi_plus].origin();
	    helixB = StPhysicalHelixD(Momentum,Origin+vectordiff,BField,Charge); // proton | redefine the helix of pi_minus after move the primary vertex from eventB to eventA
	    MomentumB = mMomentum[key_B_pi_plus][n_pi_plus];
	    dcaB = mDca[key_B_pi_plus][n_pi_plus];

	    VertexFinder->Find2ndVertex(helixB,helixA,vectorprim,MomentumB,MomentumA,ltrackB,ltrackA,VerdistX,VerdistY,vectorAB,dcaAB,2); // K0S mode

	    // Invariant mass calculations
	    TLorentzVector trackAB = ltrackA+ltrackB; // mother particle
	    Double_t InvMassAB     = trackAB.M(); // invariant mass of mother particle

	    //-----------------------------------------------------------------------------
	    // get Lambda and antiLambda by misidentification
	    TLorentzVector  ltrackP, ltrackPbar;
	    ltrackP.SetXYZM(ltrackB.X(),ltrackB.Y(),ltrackB.Z(),TriFlow::mMassProton); // set Lorentz vector of pi+ to p
	    ltrackPbar.SetXYZM(ltrackA.X(),ltrackA.Y(),ltrackA.Z(),TriFlow::mMassProton); // set Lorentz vector of pi- to pbar

	    // Invariant mass calculations Lambda and antiLambda
	    TLorentzVector trackLambda      = ltrackP+ltrackA; // p + pi-
	    TLorentzVector trackAntiLambda      = ltrackPbar+ltrackB; // pbar + pi+
	    Double_t InvMassLambda = trackLambda.M(); // invariant mass of Lambda
	    Double_t InvMassAntiLambda = trackAntiLambda.M(); // invariant mass of Lambda
	    //-----------------------------------------------------------------------------

	    if(fabs(dcaA) > mDca_pion && fabs(dcaB) > mDca_pion && fabs(dcaAB) < mDcaAB && VerdistX > mDecayLength && VerdistY < mDcaV0 && InvMassAB > mInvK0S_low && InvMassAB < mInvK0S_high)
	    {
	      Float_t Mass2_pi_plus  = mMass2[key_B_pi_plus][n_pi_plus];
	      Float_t Mass2_pi_minus = mMass2[key_A_pi_minus][n_pi_minus];

	      // final mass2 cut
	      Float_t Mass2_low_pi_plus;
	      Float_t Mass2_up_pi_plus;
	      if(mMomentum[key_B_pi_plus][n_pi_plus] < 0.75)
	      {
		Mass2_low_pi_plus = -0.013;
		Mass2_up_pi_plus  =  0.047;
	      }
	      if(mMomentum[key_B_pi_plus][n_pi_plus] >= 0.75)
	      {
		Mass2_low_pi_plus =  0.053 - 0.088*mMomentum[key_B_pi_plus][n_pi_plus];
		Mass2_up_pi_plus  = -0.004 + 0.068*mMomentum[key_B_pi_plus][n_pi_plus];
	      }

	      Float_t Mass2_low_pi_minus;
	      Float_t Mass2_up_pi_minus;
	      if(mMomentum[key_A_pi_minus][n_pi_minus] < 0.75)
	      {
		Mass2_low_pi_minus = -0.013;
		Mass2_up_pi_minus  =  0.047;
	      }
	      if(mMomentum[key_A_pi_minus][n_pi_minus] >= 0.75)
	      {
		Mass2_low_pi_minus =  0.053 - 0.088*mMomentum[key_A_pi_minus][n_pi_minus];
		Mass2_up_pi_minus  = -0.004 + 0.068*mMomentum[key_A_pi_minus][n_pi_minus];
	      }

	      StLorentzVectorD SttrackAB(trackAB.X(),trackAB.Y(),trackAB.Z(),trackAB.E());

	      // StV0TofCorrection for pi_plus from Event B
	      if(mTofFlag[key_B_pi_plus][n_pi_plus] > 0 && mTofTime[key_B_pi_plus][n_pi_plus] != 0)
	      {
		mTofCorr->setVectors3D(vectorprim)(vectorAB)(mTofHit[key_B_pi_plus][n_pi_plus]+vectordiff);
		mTofCorr->setMotherTracks(SttrackAB);
		Float_t Time_before = mTofTime[key_B_pi_plus][n_pi_plus];
		Float_t Beta_before = mTofBeta[key_B_pi_plus][n_pi_plus];
		Float_t Time_after = Time_before;
		Float_t Beta_after = Beta_before;
		mTofCorr->correctBeta(helixB,Time_after,Beta_after);
		if(Beta_after != 0) Mass2_pi_plus = mMomentum[key_B_pi_plus][n_pi_plus]*mMomentum[key_B_pi_plus][n_pi_plus]*(1.0/(Beta_after*Beta_after) - 1.0);
		mTofCorr->clearContainers();
	      }
	      // StV0TofCorrection for pion_minus from Event A
	      if(mTofFlag[key_A_pi_minus][n_pi_minus] > 0 && mTofTime[key_A_pi_minus][n_pi_minus] != 0)
	      {
		mTofCorr->setVectors3D(vectorprim)(vectorAB)(mTofHit[key_A_pi_minus][n_pi_minus]); 
		mTofCorr->setMotherTracks(SttrackAB);
		Float_t Time_before = mTofTime[key_A_pi_minus][n_pi_minus];
		Float_t Beta_before = mTofBeta[key_A_pi_minus][n_pi_minus];
		Float_t Time_after = Time_before;
		Float_t Beta_after = Beta_before;
		mTofCorr->correctBeta(helixA,Time_after,Beta_after);
		if(Beta_after != 0) Mass2_pi_minus = mMomentum[key_A_pi_minus][n_pi_minus]*mMomentum[key_A_pi_minus][n_pi_minus]*(1.0/(Beta_after*Beta_after) - 1.0);
		mTofCorr->clearContainers();
	      }

	      if(
		   ((Mass2_pi_plus  > Mass2_low_pi_plus  && Mass2_pi_plus  < Mass2_up_pi_plus)  || Mass2_pi_plus  < -10.0) // pion_plus  mass2 cut
		&& ((Mass2_pi_minus > Mass2_low_pi_minus && Mass2_pi_minus < Mass2_up_pi_minus) || Mass2_pi_minus < -10.0) // pion_minus mass2 cut
		)
	      {
		h_Mass2->Fill(trackAB.Perp(),InvMassAB);
		h_Mass2_Lambda->Fill(trackLambda.Perp(),InvMassLambda);
		h_Mass2_AntiLambda->Fill(trackAntiLambda.Perp(),InvMassAntiLambda);
		if(!((InvMassLambda > 1.1157-0.006*3 && InvMassLambda < 1.1157+0.006*3) || (InvMassAntiLambda > 1.1157-0.006*3 && InvMassAntiLambda < 1.1157+0.006*3)))
		  h_Mass2_sub->Fill(trackAB.Perp(),InvMassAB);

		// fill Lambda candidate into mTree_Phi
		mV0Track = mV0Event->createTrack();
		mV0Track->setMass2A(mMass2[key_B_pi_plus][n_pi_plus]); // pi_plus
		mV0Track->setMass2B(mMass2[key_A_pi_minus][n_pi_minus]); // pi_minus
		mV0Track->setNSigA(mNSigmaPion[key_B_pi_plus][n_pi_plus]); // pi_plus
		mV0Track->setNSigB(mNSigmaPion[key_A_pi_minus][n_pi_minus]); // pi_minus
		mV0Track->setDcaA(mDca[key_B_pi_plus][n_pi_plus]); // pi_plus
		mV0Track->setDcaB(mDca[key_A_pi_minus][n_pi_minus]); // pi_minus
		mV0Track->setGTrackA(ltrackB); // pi_plus
		mV0Track->setGTrackB(ltrackA); // pi_minus
		mV0Track->setPTrackA(mLPTrack[key_B_pi_plus][n_pi_plus]); // pi_plus
		mV0Track->setPTrackB(mLPTrack[key_A_pi_minus][n_pi_minus]); // pi_minus
		mV0Track->setFlagA(Bin_Event_B); // pi_plus
		mV0Track->setFlagB(Bin_Event_A); // pi_minus
		mV0Track->setDcaAB(dcaAB);
		mV0Track->setDecayLength(VerdistX);
		mV0Track->setDcaV0(VerdistY);
	      }
	    }
	  }
	}
      }
    }
    mTree_K0S->Fill();
  }
}

//------------------------------------------------------------------------------------------------------------------

void StTriFlowV0::MixEvent_Phi(Int_t Flag_ME, StPicoDst *pico, Int_t cent9, Float_t vz, Float_t Psi2)
{
  StPicoEvent *event = (StPicoEvent*)pico->event();

  Int_t Bin_vz, Bin_Psi2;

  Float_t vz_start = TriFlow::mVzMaxMap[event->energy()];
  Float_t vz_bin = 2*vz_start/TriFlow::Bin_VertexZ;

  Float_t psi2_start = TMath::Pi()/2.0;
  Float_t psi2_bin = 2*psi2_start/TriFlow::Bin_Phi_Psi;

  for(Int_t i = 0; i < TriFlow::Bin_VertexZ; i++)
  {
    if((vz > -1.0*vz_start+i*vz_bin) && (vz <= -1.0*vz_start+(i+1)*vz_bin))
    {
      Bin_vz = i;
    }
  }
  for(Int_t i = 0; i < TriFlow::Bin_Phi_Psi; i++)
  {
    if((Psi2 > -1.0*psi2_start+i*psi2_bin) && (Psi2 <= -1.0*psi2_start+(i+1)*psi2_bin))
    {
      Bin_Psi2 = i;
    }
  }

  Int_t Bin_Event = mEventCounter2[cent9][Bin_vz][Bin_Psi2];

  const Double_t MAGFIELDFACTOR = kilogauss;
  const Int_t nTracks = pico->numberOfTracks();

  // store Enent Information
  mPrimaryvertex[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<StThreeVectorF>(event->primaryVertex()));
  mRefMult[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(event->refMult()));
  mCentrality[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(cent9));
  mRunId[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(event->runId()));
  mN_prim[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(mNumber_prim));
  mN_non_prim[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(mNumber_non_prim));
  mN_Tof_match[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(mNumber_Tof_match));
  mZDCx[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Float_t>(event->ZDCx()));
  mBBCx[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Float_t>(event->BBCx()));
  mVzVpd[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Float_t>(event->vzVpd()));
  mNumTracks[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<UShort_t>(pico->numberOfTracks()));
  for(Int_t j = 0; j < TriFlow::EtaGap_total; j++)
  {
    mQ2East[cent9][Bin_vz][Bin_Psi2][j].push_back(static_cast<TVector2>(mQVector2East[j]));
    mQ2West[cent9][Bin_vz][Bin_Psi2][j].push_back(static_cast<TVector2>(mQVector2West[j]));
    mQ3East[cent9][Bin_vz][Bin_Psi2][j].push_back(static_cast<TVector2>(mQVector3East[j]));
    mQ3West[cent9][Bin_vz][Bin_Psi2][j].push_back(static_cast<TVector2>(mQVector3West[j]));
    mNumTrackEast[cent9][Bin_vz][Bin_Psi2][j].push_back(static_cast<Int_t>(mTrackEtaEast[j]));
    mNumTrackWest[cent9][Bin_vz][Bin_Psi2][j].push_back(static_cast<Int_t>(mTrackEtaWest[j]));
  }

  // store Track Information
  for(Int_t i = 0; i < nTracks; i ++) // loop over all particles in event
  {
    StPicoTrack *track = pico->track(i);

    if(mTriFlowCut->passTrackPhi(track))
    {
      Float_t Mass2 = mTriFlowCut->getMass2(track);
      Float_t scale_nSigma_factor = TriFlow::mSigScaleMap[event->energy()];
      Float_t Polarity = static_cast<Float_t>(track->charge());
      Float_t momentum = track->pMom().mag();
      Float_t Mass2_low;
      Float_t Mass2_up;
      if(momentum < 0.5)
      {
        Mass2_low = 0.4*0.4;
	Mass2_up = 0.6*0.6;
      }
      if(momentum >= 0.5)
      {
	Mass2_low = 0.277205 - 0.0812931*momentum;
	Mass2_up = 0.215517 + 0.076801*momentum;
      }

      Int_t charge = 0; // k+
      if(Polarity < 0) charge = 1; // k-


      if(mTriFlowCut->passSigKaonCut(track,scale_nSigma_factor))
      {
	if(
	    (momentum < 0.65 && ((Mass2 > Mass2_low && Mass2 < Mass2_up) || Mass2 < -10.0)) // dE/dx + ToF
	    || (momentum >= 0.65 && (Mass2 > Mass2_low && Mass2 < Mass2_up)) // dE/dx + ToF(always)
	  )
	{
	  MEKey key = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,charge);
	  mMass2[key].push_back(static_cast<Float_t>(mTriFlowCut->getMass2(track))); // mass2
	  mDca[key].push_back(static_cast<Float_t>(track->dca()*track->charge())); // dca*charge 
	  mNHitsFit[key].push_back(static_cast<Float_t>(track->nHitsFit())); // nHitsFit
	  mNSigmaKaon[key].push_back(static_cast<Float_t>((track->nSigmaKaon())*scale_nSigma_factor)); // nSigmaKaon
	  mHelix_Kaon[key].push_back(static_cast<StPhysicalHelixD>(StPhysicalHelixD(track->pMom(),event->primaryVertex(),event->bField()*MAGFIELDFACTOR,track->charge())));// get helix from the pMom 
	  mMomentum[key].push_back(static_cast<Float_t>(track->pMom().mag()));// get helix from the pMom 
	}
      }
    }
  }

  mEventCounter2[cent9][Bin_vz][Bin_Psi2]++;

  if(Flag_ME == 0) // same event
  {
    doPhi(Flag_ME,cent9,Bin_vz,Bin_Psi2);
    clear_phi(cent9,Bin_vz,Bin_Psi2);
  }

  if(Flag_ME == 1) // mix event
  {
    if(mEventCounter2[cent9][Bin_vz][Bin_Psi2] == TriFlow::Buffer_depth)
    {
      doPhi(Flag_ME,cent9,Bin_vz,Bin_Psi2);
      clear_phi(cent9,Bin_vz,Bin_Psi2);
    }
  }
}

//------------------------------------------------------------------------------------------------------------------

void StTriFlowV0::MixEvent_Lambda(Int_t Flag_ME, StPicoDst *pico, Int_t cent9, Float_t vz, Float_t Psi2) // proton + pi-
{
  StPicoEvent *event = (StPicoEvent*)pico->event();

  Int_t Bin_vz = 0;
  Int_t Bin_Psi2 = 0;

  Float_t vz_start = TriFlow::mVzMaxMap[event->energy()];
  Float_t vz_bin = 2*vz_start/TriFlow::Bin_VertexZ;

  Float_t psi2_start = TMath::Pi()/2.0;
  Float_t psi2_bin = 2*psi2_start/TriFlow::Bin_Phi_Psi;

  for(Int_t i = 0; i < TriFlow::Bin_VertexZ; i++)
  {
    if((vz > -1.0*vz_start+i*vz_bin) && (vz <= -1.0*vz_start+(i+1)*vz_bin))
    {
      Bin_vz = i;
    }
  }
  for(Int_t i = 0; i < TriFlow::Bin_Phi_Psi; i++)
  {
    if((Psi2 > -1.0*psi2_start+i*psi2_bin) && (Psi2 <= -1.0*psi2_start+(i+1)*psi2_bin))
    {
      Bin_Psi2 = i;
    }
  }

  Int_t Bin_Event = mEventCounter2[cent9][Bin_vz][Bin_Psi2];

  const Double_t MAGFIELDFACTOR = kilogauss;
  const Int_t nTracks = pico->numberOfTracks();

  // store Event Information
  mPrimaryvertex[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<StThreeVectorF>(event->primaryVertex()));
  mRefMult[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(event->refMult()));
  mCentrality[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(cent9));
  mRunId[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(event->runId()));
  mN_prim[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(mNumber_prim));
  mN_non_prim[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(mNumber_non_prim));
  mN_Tof_match[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(mNumber_Tof_match));
  mZDCx[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Float_t>(event->ZDCx()));
  mBBCx[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Float_t>(event->BBCx()));
  mVzVpd[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Float_t>(event->vzVpd()));
  mField[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Float_t>(event->bField()*MAGFIELDFACTOR));
  mNumTracks[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<UShort_t>(pico->numberOfTracks()));
  for(Int_t j = 0; j < TriFlow::EtaGap_total; j++)
  {
    mQ2East[cent9][Bin_vz][Bin_Psi2][j].push_back(static_cast<TVector2>(mQVector2East[j]));
    mQ2West[cent9][Bin_vz][Bin_Psi2][j].push_back(static_cast<TVector2>(mQVector2West[j]));
    mQ3East[cent9][Bin_vz][Bin_Psi2][j].push_back(static_cast<TVector2>(mQVector3East[j]));
    mQ3West[cent9][Bin_vz][Bin_Psi2][j].push_back(static_cast<TVector2>(mQVector3West[j]));
    mNumTrackEast[cent9][Bin_vz][Bin_Psi2][j].push_back(static_cast<Int_t>(mTrackEtaEast[j]));
    mNumTrackWest[cent9][Bin_vz][Bin_Psi2][j].push_back(static_cast<Int_t>(mTrackEtaWest[j]));
  }

  // store Track Information
  for(Int_t i = 0; i < nTracks; i ++) // loop over all particles in event
  {
    StPicoTrack *track = pico->track(i);

    if(mTriFlowCut->passTrackV0(track)) // global pt > 0.1
    {
      Float_t Mass2 = mTriFlowCut->getV0Mass2(track);
      Float_t scale_nSigma_factor = TriFlow::mSigScaleMap[event->energy()];
      Float_t Polarity = static_cast<Float_t>(track->charge());
      Float_t momentum = track->gMom().mag();

      // initial mass2 cut
      Float_t Mass2_low_proton = 0.45;
      Float_t Mass2_up_proton  = 1.5;
      Float_t Mass2_low_pi_minus = -0.5;
      Float_t Mass2_up_pi_minus  = 0.1;

      h_Mass2_p->Fill(momentum*Polarity,Mass2);

      Int_t charge = 0; // p 
      if(Polarity < 0) charge = 1; // pi-

      if(charge == 0 && mTriFlowCut->passSigProntonCut(track,scale_nSigma_factor)) // p
      {
	if(
	    (Mass2 > Mass2_low_proton && Mass2 < Mass2_up_proton) || Mass2 < -10.0 // dE/dx + ToF
	  )
	{
	  MEKey key = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,charge);
	  mMass2[key].push_back(static_cast<Float_t>(mTriFlowCut->getV0Mass2(track))); // mass2
	  mDca[key].push_back(static_cast<Float_t>(track->dca()*track->charge())); // dca*charge 
	  mNHitsFit[key].push_back(static_cast<Float_t>(track->nHitsFit())); // nHitsFit
	  mNSigmaProton[key].push_back(static_cast<Float_t>((track->nSigmaProton())*scale_nSigma_factor)); // nSigmaProton
	  mHelix_Proton[key].push_back(static_cast<StPhysicalHelixD>(StPhysicalHelixD(track->gMom(),track->origin(),event->bField()*MAGFIELDFACTOR,track->charge())));// get helix from the gMom 
	  mMomentum[key].push_back(static_cast<Float_t>(track->gMom().mag()));// get helix from the gMom 
	  TLorentzVector lptrack(track->pMom().x(),track->pMom().y(),track->pMom().z(),TriFlow::mMassProton);
//	  cout << "gMom.mag() = " << track->gMom().mag() << endl;
//	  cout << "pMom.mag() = " << track->pMom().mag() << endl;
//	  cout << endl;
	  mLPTrack[key].push_back(static_cast<TLorentzVector>(lptrack));

	  // StV0TofCorr
	  mTofFlag[key].push_back(static_cast<Int_t>(track->btofMatchFlag()));
	  mTofTime[key].push_back(static_cast<Float_t>(track->btof()));
	  if(track->btofMatchFlag() > 0 && track->btof() != 0)
	  {
	    mTofHit[key].push_back(static_cast<StThreeVectorF>(track->btofHisPos()));
	    mTofBeta[key].push_back(static_cast<Float_t>(track->btofBeta()));
	  }
	  else
	  {
	    StThreeVectorF FakeHit(-999.9,-999.9,-999.9);
	    Float_t FakeBeta = -999.9;
	    mTofHit[key].push_back(static_cast<StThreeVectorF>(FakeHit));
	    mTofBeta[key].push_back(static_cast<Float_t>(FakeBeta));
	  }
	}
      }
      if(charge == 1 && mTriFlowCut->passSigPionCut(track,scale_nSigma_factor)) // pi-
      {
	if(
	    ((Mass2 > Mass2_low_pi_minus && Mass2 < Mass2_up_pi_minus) || Mass2 < -10.0) // pre-mass cut: dE/dx + ToF
	  && (track->dca() > mDca_pion_Pre) // pre-dca cut
	  )
	{
	  MEKey key = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,charge);
	  mMass2[key].push_back(static_cast<Float_t>(mTriFlowCut->getV0Mass2(track))); // mass2
	  mDca[key].push_back(static_cast<Float_t>(track->dca()*track->charge())); // dca*charge 
	  mNHitsFit[key].push_back(static_cast<Float_t>(track->nHitsFit())); // nHitsFit
	  mNSigmaPion[key].push_back(static_cast<Float_t>((track->nSigmaPion())*scale_nSigma_factor)); // nSigmaPion
	  mHelix_Pion[key].push_back(static_cast<StPhysicalHelixD>(StPhysicalHelixD(track->gMom(),track->origin(),event->bField()*MAGFIELDFACTOR,track->charge())));// get helix from the gMom 
	  mMomentum[key].push_back(static_cast<Float_t>(track->gMom().mag()));// get helix from the gMom 
	  TLorentzVector lptrack(track->pMom().x(),track->pMom().y(),track->pMom().z(),TriFlow::mMassPion);
	  mLPTrack[key].push_back(static_cast<TLorentzVector>(lptrack));


	  // StV0TofCorr
	  mTofFlag[key].push_back(static_cast<Int_t>(track->btofMatchFlag()));
	  mTofTime[key].push_back(static_cast<Float_t>(track->btof()));
	  if(track->btofMatchFlag() > 0 && track->btof() != 0)
	  {
	    mTofHit[key].push_back(static_cast<StThreeVectorF>(track->btofHisPos()));
	    mTofBeta[key].push_back(static_cast<Float_t>(track->btofBeta()));
	  }
	  else
	  {
	    StThreeVectorF FakeHit(-999.9,-999.9,-999.9);
	    Float_t FakeBeta = -999.9;
	    mTofHit[key].push_back(static_cast<StThreeVectorF>(FakeHit));
	    mTofBeta[key].push_back(static_cast<Float_t>(FakeBeta));
	  }
	}
      }
    }
  }

  mEventCounter2[cent9][Bin_vz][Bin_Psi2]++;

  if(Flag_ME == 0) // same event
  {
    doLambda(Flag_ME,cent9,Bin_vz,Bin_Psi2);
    clear_Lambda(cent9,Bin_vz,Bin_Psi2);
  }

  if(Flag_ME == 1) // mix event
  {
    if(mEventCounter2[cent9][Bin_vz][Bin_Psi2] == TriFlow::Buffer_depth)
    {
      doLambda(Flag_ME,cent9,Bin_vz,Bin_Psi2);
      clear_Lambda(cent9,Bin_vz,Bin_Psi2);
    }
  }
}

void StTriFlowV0::MixEvent_AntiLambda(Int_t Flag_ME, StPicoDst *pico, Int_t cent9, Float_t vz, Float_t Psi2) // anti-proton + pi+
{
  StPicoEvent *event = (StPicoEvent*)pico->event();

  Int_t Bin_vz = 0;
  Int_t Bin_Psi2 = 0;

  Float_t vz_start = TriFlow::mVzMaxMap[event->energy()];
  Float_t vz_bin = 2*vz_start/TriFlow::Bin_VertexZ;

  Float_t psi2_start = TMath::Pi()/2.0;
  Float_t psi2_bin = 2*psi2_start/TriFlow::Bin_Phi_Psi;

  for(Int_t i = 0; i < TriFlow::Bin_VertexZ; i++)
  {
    if((vz > -1.0*vz_start+i*vz_bin) && (vz <= -1.0*vz_start+(i+1)*vz_bin))
    {
      Bin_vz = i;
    }
  }
  for(Int_t i = 0; i < TriFlow::Bin_Phi_Psi; i++)
  {
    if((Psi2 > -1.0*psi2_start+i*psi2_bin) && (Psi2 <= -1.0*psi2_start+(i+1)*psi2_bin))
    {
      Bin_Psi2 = i;
    }
  }

  Int_t Bin_Event = mEventCounter2[cent9][Bin_vz][Bin_Psi2];

  const Double_t MAGFIELDFACTOR = kilogauss;
  const Int_t nTracks = pico->numberOfTracks();

  // store Event Information
  mPrimaryvertex[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<StThreeVectorF>(event->primaryVertex()));
  mRefMult[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(event->refMult()));
  mCentrality[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(cent9));
  mRunId[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(event->runId()));
  mN_prim[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(mNumber_prim));
  mN_non_prim[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(mNumber_non_prim));
  mN_Tof_match[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(mNumber_Tof_match));
  mZDCx[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Float_t>(event->ZDCx()));
  mBBCx[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Float_t>(event->BBCx()));
  mVzVpd[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Float_t>(event->vzVpd()));
  mField[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Float_t>(event->bField()*MAGFIELDFACTOR));
  mNumTracks[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<UShort_t>(pico->numberOfTracks()));
  for(Int_t j = 0; j < TriFlow::EtaGap_total; j++)
  {
    mQ2East[cent9][Bin_vz][Bin_Psi2][j].push_back(static_cast<TVector2>(mQVector2East[j]));
    mQ2West[cent9][Bin_vz][Bin_Psi2][j].push_back(static_cast<TVector2>(mQVector2West[j]));
    mQ3East[cent9][Bin_vz][Bin_Psi2][j].push_back(static_cast<TVector2>(mQVector3East[j]));
    mQ3West[cent9][Bin_vz][Bin_Psi2][j].push_back(static_cast<TVector2>(mQVector3West[j]));
    mNumTrackEast[cent9][Bin_vz][Bin_Psi2][j].push_back(static_cast<Int_t>(mTrackEtaEast[j]));
    mNumTrackWest[cent9][Bin_vz][Bin_Psi2][j].push_back(static_cast<Int_t>(mTrackEtaWest[j]));
  }

  // store Track Information
  for(Int_t i = 0; i < nTracks; i ++) // loop over all particles in event
  {
    StPicoTrack *track = pico->track(i);

    if(mTriFlowCut->passTrackV0(track)) // global pt > 0.1
    {
      Float_t Mass2 = mTriFlowCut->getV0Mass2(track);
      Float_t scale_nSigma_factor = TriFlow::mSigScaleMap[event->energy()];
      Float_t Polarity = static_cast<Float_t>(track->charge());
      Float_t momentum = track->gMom().mag();

      // initial mass2 cut
      Float_t Mass2_low_proton = 0.45;
      Float_t Mass2_up_proton  = 1.5;
      Float_t Mass2_low_pi_minus = -0.5;
      Float_t Mass2_up_pi_minus  = 0.1;

      h_Mass2_p->Fill(momentum*Polarity,Mass2);

      Int_t charge = 0; // pi+ 
      if(Polarity < 0) charge = 1; // anti-proton

      if(charge == 1 && mTriFlowCut->passSigProntonCut(track,scale_nSigma_factor)) // anti-proton 
      {
	if(
	    (Mass2 > Mass2_low_proton && Mass2 < Mass2_up_proton) || Mass2 < -10.0 // dE/dx + ToF
	  )
	{
	  MEKey key = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,charge);
	  mMass2[key].push_back(static_cast<Float_t>(mTriFlowCut->getV0Mass2(track))); // mass2
	  mDca[key].push_back(static_cast<Float_t>(track->dca()*track->charge())); // dca*charge 
	  mNHitsFit[key].push_back(static_cast<Float_t>(track->nHitsFit())); // nHitsFit
	  mNSigmaProton[key].push_back(static_cast<Float_t>((track->nSigmaProton())*scale_nSigma_factor)); // nSigmaProton
	  mHelix_Proton[key].push_back(static_cast<StPhysicalHelixD>(StPhysicalHelixD(track->gMom(),track->origin(),event->bField()*MAGFIELDFACTOR,track->charge())));// get helix from the gMom 
	  mMomentum[key].push_back(static_cast<Float_t>(track->gMom().mag()));// get helix from the gMom 
	  TLorentzVector lptrack(track->pMom().x(),track->pMom().y(),track->pMom().z(),TriFlow::mMassProton);
	  mLPTrack[key].push_back(static_cast<TLorentzVector>(lptrack));

	  // StV0TofCorr
	  mTofFlag[key].push_back(static_cast<Int_t>(track->btofMatchFlag()));
	  mTofTime[key].push_back(static_cast<Float_t>(track->btof()));
	  if(track->btofMatchFlag() > 0 && track->btof() != 0)
	  {
	    mTofHit[key].push_back(static_cast<StThreeVectorF>(track->btofHisPos()));
	    mTofBeta[key].push_back(static_cast<Float_t>(track->btofBeta()));
	  }
	  else
	  {
	    StThreeVectorF FakeHit(-999.9,-999.9,-999.9);
	    Float_t FakeBeta = -999.9;
	    mTofHit[key].push_back(static_cast<StThreeVectorF>(FakeHit));
	    mTofBeta[key].push_back(static_cast<Float_t>(FakeBeta));
	  }
	}
      }
      if(charge == 0 && mTriFlowCut->passSigPionCut(track,scale_nSigma_factor)) // pi+
      {
	if(
	    ((Mass2 > Mass2_low_pi_minus && Mass2 < Mass2_up_pi_minus) || Mass2 < -10.0) // pre-mass cut: dE/dx + ToF
	  && (track->dca() > mDca_pion_Pre) // pre-dca cut
	  )
	{
	  MEKey key = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,charge);
	  mMass2[key].push_back(static_cast<Float_t>(mTriFlowCut->getV0Mass2(track))); // mass2
	  mDca[key].push_back(static_cast<Float_t>(track->dca()*track->charge())); // dca*charge 
	  mNHitsFit[key].push_back(static_cast<Float_t>(track->nHitsFit())); // nHitsFit
	  mNSigmaPion[key].push_back(static_cast<Float_t>((track->nSigmaPion())*scale_nSigma_factor)); // nSigmaPion
	  mHelix_Pion[key].push_back(static_cast<StPhysicalHelixD>(StPhysicalHelixD(track->gMom(),track->origin(),event->bField()*MAGFIELDFACTOR,track->charge())));// get helix from the gMom 
	  mMomentum[key].push_back(static_cast<Float_t>(track->gMom().mag()));// get helix from the gMom 
	  TLorentzVector lptrack(track->pMom().x(),track->pMom().y(),track->pMom().z(),TriFlow::mMassPion);
	  mLPTrack[key].push_back(static_cast<TLorentzVector>(lptrack));


	  // StV0TofCorr
	  mTofFlag[key].push_back(static_cast<Int_t>(track->btofMatchFlag()));
	  mTofTime[key].push_back(static_cast<Float_t>(track->btof()));
	  if(track->btofMatchFlag() > 0 && track->btof() != 0)
	  {
	    mTofHit[key].push_back(static_cast<StThreeVectorF>(track->btofHisPos()));
	    mTofBeta[key].push_back(static_cast<Float_t>(track->btofBeta()));
	  }
	  else
	  {
	    StThreeVectorF FakeHit(-999.9,-999.9,-999.9);
	    Float_t FakeBeta = -999.9;
	    mTofHit[key].push_back(static_cast<StThreeVectorF>(FakeHit));
	    mTofBeta[key].push_back(static_cast<Float_t>(FakeBeta));
	  }
	}
      }
    }
  }

  mEventCounter2[cent9][Bin_vz][Bin_Psi2]++;

  if(Flag_ME == 0) // same event
  {
    doAntiLambda(Flag_ME,cent9,Bin_vz,Bin_Psi2);
    clear_AntiLambda(cent9,Bin_vz,Bin_Psi2);
  }

  if(Flag_ME == 1) // mix event
  {
    if(mEventCounter2[cent9][Bin_vz][Bin_Psi2] == TriFlow::Buffer_depth)
    {
      doAntiLambda(Flag_ME,cent9,Bin_vz,Bin_Psi2);
      clear_AntiLambda(cent9,Bin_vz,Bin_Psi2);
    }
  }
}

//------------------------------------------------------------------------------------------------------------------

void StTriFlowV0::MixEvent_K0S(Int_t Flag_ME, StPicoDst *pico, Int_t cent9, Float_t vz, Float_t Psi2) // pi+ + pi-
{
  StPicoEvent *event = (StPicoEvent*)pico->event();

  Int_t Bin_vz = 0;
  Int_t Bin_Psi2 = 0;

  Float_t vz_start = TriFlow::mVzMaxMap[event->energy()];
  Float_t vz_bin = 2*vz_start/TriFlow::Bin_VertexZ;

  Float_t psi2_start = TMath::Pi()/2.0;
  Float_t psi2_bin = 2*psi2_start/TriFlow::Bin_Phi_Psi;

  for(Int_t i = 0; i < TriFlow::Bin_VertexZ; i++)
  {
    if((vz > -1.0*vz_start+i*vz_bin) && (vz <= -1.0*vz_start+(i+1)*vz_bin))
    {
      Bin_vz = i;
    }
  }
  for(Int_t i = 0; i < TriFlow::Bin_Phi_Psi; i++)
  {
    if((Psi2 > -1.0*psi2_start+i*psi2_bin) && (Psi2 <= -1.0*psi2_start+(i+1)*psi2_bin))
    {
      Bin_Psi2 = i;
    }
  }

  Int_t Bin_Event = mEventCounter2[cent9][Bin_vz][Bin_Psi2];

  const Double_t MAGFIELDFACTOR = kilogauss;
  const Int_t nTracks = pico->numberOfTracks();

  // store Event Information
  mPrimaryvertex[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<StThreeVectorF>(event->primaryVertex()));
  mRefMult[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(event->refMult()));
  mCentrality[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(cent9));
  mRunId[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(event->runId()));
  mN_prim[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(mNumber_prim));
  mN_non_prim[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(mNumber_non_prim));
  mN_Tof_match[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(mNumber_Tof_match));
  mZDCx[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Float_t>(event->ZDCx()));
  mBBCx[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Float_t>(event->BBCx()));
  mVzVpd[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Float_t>(event->vzVpd()));
  mField[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Float_t>(event->bField()*MAGFIELDFACTOR));
  mNumTracks[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<UShort_t>(pico->numberOfTracks()));
  for(Int_t j = 0; j < TriFlow::EtaGap_total; j++)
  {
    mQ2East[cent9][Bin_vz][Bin_Psi2][j].push_back(static_cast<TVector2>(mQVector2East[j]));
    mQ2West[cent9][Bin_vz][Bin_Psi2][j].push_back(static_cast<TVector2>(mQVector2West[j]));
    mQ3East[cent9][Bin_vz][Bin_Psi2][j].push_back(static_cast<TVector2>(mQVector3East[j]));
    mQ3West[cent9][Bin_vz][Bin_Psi2][j].push_back(static_cast<TVector2>(mQVector3West[j]));
    mNumTrackEast[cent9][Bin_vz][Bin_Psi2][j].push_back(static_cast<Int_t>(mTrackEtaEast[j]));
    mNumTrackWest[cent9][Bin_vz][Bin_Psi2][j].push_back(static_cast<Int_t>(mTrackEtaWest[j]));
  }

  // store Track Information
  for(Int_t i = 0; i < nTracks; i ++) // loop over all particles in event
  {
    StPicoTrack *track = pico->track(i);

    if(mTriFlowCut->passTrackV0(track)) // global pt > 0.1
    {
      Float_t Mass2 = mTriFlowCut->getV0Mass2(track);
      Float_t scale_nSigma_factor = TriFlow::mSigScaleMap[event->energy()];
      Float_t Polarity = static_cast<Float_t>(track->charge());
      Float_t momentum = track->gMom().mag();

      // initial mass2 cut
      Float_t Mass2_low_pi_minus = -0.5;
      Float_t Mass2_up_pi_minus  = 0.1;

      h_Mass2_p->Fill(momentum*Polarity,Mass2);

      Int_t charge = 0; // pi+ 
      if(Polarity < 0) charge = 1; // pi-

      if(mTriFlowCut->passSigPionCut(track,scale_nSigma_factor)) // pi+/pi-
      {
	if(
	    ((Mass2 > Mass2_low_pi_minus && Mass2 < Mass2_up_pi_minus) || Mass2 < -10.0) // pre-mass cut: dE/dx + ToF
	  && (track->dca() > mDca_pion_Pre) // pre-dca cut
	  )
	{
	  MEKey key = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,charge);
	  mMass2[key].push_back(static_cast<Float_t>(mTriFlowCut->getV0Mass2(track))); // mass2
	  mDca[key].push_back(static_cast<Float_t>(track->dca()*track->charge())); // dca*charge 
	  mNHitsFit[key].push_back(static_cast<Float_t>(track->nHitsFit())); // nHitsFit
	  mNSigmaPion[key].push_back(static_cast<Float_t>((track->nSigmaPion())*scale_nSigma_factor)); // nSigmaPion
	  mHelix_Pion[key].push_back(static_cast<StPhysicalHelixD>(StPhysicalHelixD(track->gMom(),track->origin(),event->bField()*MAGFIELDFACTOR,track->charge())));// get helix from the gMom 
	  mMomentum[key].push_back(static_cast<Float_t>(track->gMom().mag()));// get helix from the gMom 
	  TLorentzVector lptrack(track->pMom().x(),track->pMom().y(),track->pMom().z(),TriFlow::mMassPion);
	  mLPTrack[key].push_back(static_cast<TLorentzVector>(lptrack));


	  // StV0TofCorr
	  mTofFlag[key].push_back(static_cast<Int_t>(track->btofMatchFlag()));
	  mTofTime[key].push_back(static_cast<Float_t>(track->btof()));
	  if(track->btofMatchFlag() > 0 && track->btof() != 0)
	  {
	    mTofHit[key].push_back(static_cast<StThreeVectorF>(track->btofHisPos()));
	    mTofBeta[key].push_back(static_cast<Float_t>(track->btofBeta()));
	  }
	  else
	  {
	    StThreeVectorF FakeHit(-999.9,-999.9,-999.9);
	    Float_t FakeBeta = -999.9;
	    mTofHit[key].push_back(static_cast<StThreeVectorF>(FakeHit));
	    mTofBeta[key].push_back(static_cast<Float_t>(FakeBeta));
	  }
	}
      }
    }
  }

  mEventCounter2[cent9][Bin_vz][Bin_Psi2]++;

  if(Flag_ME == 0) // same event
  {
    doK0S(Flag_ME,cent9,Bin_vz,Bin_Psi2);
    clear_K0S(cent9,Bin_vz,Bin_Psi2);
  }

  if(Flag_ME == 1) // mix event
  {
    if(mEventCounter2[cent9][Bin_vz][Bin_Psi2] == TriFlow::Buffer_depth)
    {
      doK0S(Flag_ME,cent9,Bin_vz,Bin_Psi2);
      clear_K0S(cent9,Bin_vz,Bin_Psi2);
    }
  }
}

//------------------------------------------------------------------------------------------------------------------
void StTriFlowV0::SetTopoCut(Float_t dca_proton, Float_t dca_pion, Float_t dcaAB, Float_t decaylength, Float_t dcaV0, Float_t dca_pion_pre, Float_t InvLambda_high)
{
  mDca_proton     = dca_proton;
  mDca_pion       = dca_pion;
  mDcaAB          = dcaAB;
  mDecayLength     = decaylength;
  mDcaV0          = dcaV0;
  mDca_pion_Pre   = dca_pion_pre;
  mInvLambda_low  = TriFlow::mMassProton+TriFlow::mMassPion;
  mInvLambda_high = InvLambda_high;
}

void StTriFlowV0::SetTopoCutK0S(Float_t dca_pion, Float_t dcaAB, Float_t decaylength, Float_t dcaV0, Float_t dca_pion_pre, Float_t InvK0S_high)
{
  mDca_pion       = dca_pion;
  mDcaAB          = dcaAB;
  mDecayLength     = decaylength;
  mDcaV0          = dcaV0;
  mDca_pion_Pre   = dca_pion_pre;
  mInvK0S_low  = TriFlow::mMassPion+TriFlow::mMassPion;
  mInvK0S_high = InvK0S_high;
}

void StTriFlowV0::PrintTopoCut()
{
  cout << "dca_proton = " << mDca_proton << endl;
  cout << "dca_pion   = " << mDca_pion   << endl;
  cout << "decaylength = " << mDecayLength << endl;
  cout << "dcaV0      = " << mDcaV0      << endl;
}
//------------------------------------------------------------------------------------------------------------------
// pass event information from Maker
void StTriFlowV0::clearEvent()
{
  mNumber_prim = 0;
  mNumber_non_prim = 0;
  mNumber_Tof_match = 0;

  for(Int_t j = 0; j < TriFlow::EtaGap_total; j++)
  {
    mQVector2East[j].Set(-999.9,-999.9);
    mQVector2West[j].Set(-999.9,-999.9);
    mQVector3East[j].Set(-999.9,-999.9);
    mQVector3West[j].Set(-999.9,-999.9);
  }
}

void StTriFlowV0::passEvent(Int_t N_prim, Int_t N_non_prim, Int_t N_Tof_match)
{
  mNumber_prim = N_prim;
  mNumber_non_prim = N_non_prim;
  mNumber_Tof_match = N_Tof_match;
}

void StTriFlowV0::passEventPlane2East(TVector2 Q2East_0, TVector2 Q2East_1, TVector2 Q2East_2, TVector2 Q2East_3)
{
  mQVector2East[0] = Q2East_0;
  mQVector2East[1] = Q2East_1;
  mQVector2East[2] = Q2East_2;
  mQVector2East[3] = Q2East_3;
}

void StTriFlowV0::passEventPlane2West(TVector2 Q2West_0, TVector2 Q2West_1, TVector2 Q2West_2, TVector2 Q2West_3)
{
  mQVector2West[0] = Q2West_0;
  mQVector2West[1] = Q2West_1;
  mQVector2West[2] = Q2West_2;
  mQVector2West[3] = Q2West_3;
}

void StTriFlowV0::passEventPlane3East(TVector2 Q3East_0, TVector2 Q3East_1, TVector2 Q3East_2, TVector2 Q3East_3)
{
  mQVector3East[0] = Q3East_0;
  mQVector3East[1] = Q3East_1;
  mQVector3East[2] = Q3East_2;
  mQVector3East[3] = Q3East_3;
}

void StTriFlowV0::passEventPlane3West(TVector2 Q3West_0, TVector2 Q3West_1, TVector2 Q3West_2, TVector2 Q3West_3)
{
  mQVector3West[0] = Q3West_0;
  mQVector3West[1] = Q3West_1;
  mQVector3West[2] = Q3West_2;
  mQVector3West[3] = Q3West_3;
}

void StTriFlowV0::passNumTrackEast(Int_t NumTrackEast_0, Int_t NumTrackEast_1, Int_t NumTrackEast_2, Int_t NumTrackEast_3)
{
  mTrackEtaEast[0] = NumTrackEast_0;
  mTrackEtaEast[1] = NumTrackEast_1;
  mTrackEtaEast[2] = NumTrackEast_2;
  mTrackEtaEast[3] = NumTrackEast_3;
}

void StTriFlowV0::passNumTrackWest(Int_t NumTrackWest_0, Int_t NumTrackWest_1, Int_t NumTrackWest_2, Int_t NumTrackWest_3)
{
  mTrackEtaWest[0] = NumTrackWest_0;
  mTrackEtaWest[1] = NumTrackWest_1;
  mTrackEtaWest[2] = NumTrackWest_2;
  mTrackEtaWest[3] = NumTrackWest_3;
}
//------------------------------------------------------------------------------------------------------------------

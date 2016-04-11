#include "StTriFlowMaker.h"
#include "StTriFlowConstants.h"
#include "StTriFlowCut.h"
#include "StTriFlowProManger.h"
#include "StTriFlowCorrection.h"
#include "StTriFlowHistoManger.h"
#include "StTriFlowV0.h"
#include "StRoot/StPicoDstMaker/StPicoDst.h"
#include "StRoot/StPicoDstMaker/StPicoEvent.h"
#include "StRoot/StPicoDstMaker/StPicoTrack.h"
#include "StRoot/StPicoDstMaker/StPicoV0.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"
#include "StRoot/StRunIdEventsDb/StRunIdEventsDb.h"
#include "StRoot/StCombPID/StCombPID.h"
#include "StThreeVectorF.hh"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TVector3.h"
#include "TMath.h"
#include "StMessMgr.h"
#include <algorithm>

ClassImp(StTriFlowMaker)

StRefMultCorr* StTriFlowMaker::mRefMultCorr = NULL;
//-----------------------------------------------------------------------------
StTriFlowMaker::StTriFlowMaker(const char* name, StPicoDstMaker *picoMaker, const Int_t jobCounter, const Int_t Mode, const Int_t energy, const Int_t flag_ME)
  : StMaker(name)
{
  mPicoDstMaker = picoMaker;
  mPicoDst = 0;
  mMode = Mode;
  mEnergy = energy;
  mFlag_ME = flag_ME;


  if(mMode == 0)
  {
    mOutPut_ReCenterPar = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/RecenterParameter/file_%s_ReCenterPar_",TriFlow::Energy[energy].Data(),TriFlow::Energy[energy].Data());
    mOutPut_ReCenterPar += jobCounter;
    mOutPut_ReCenterPar += ".root";
  }
  if(mMode == 1)
  {
    mOutPut_Corr_Shift = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/Correction/Shift/file_%s_Corr_Shift_",TriFlow::Energy[energy].Data(),TriFlow::Energy[energy].Data()); 
    mOutPut_Corr_Shift += jobCounter;
    mOutPut_Corr_Shift += ".root";

    mOutPut_Corr_ReCenter = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/Correction/ReCenter/file_%s_Corr_ReCenter_",TriFlow::Energy[energy].Data(),TriFlow::Energy[energy].Data()); 
    mOutPut_Corr_ReCenter += jobCounter;
    mOutPut_Corr_ReCenter += ".root";
  }
  if(mMode == 2)
  {
    mOutPut_ChargedFLow = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/ChargedFlow/file_%s_ChargedFlow_",TriFlow::Energy[energy].Data(),TriFlow::Energy[energy].Data()); 
    mOutPut_ChargedFLow += jobCounter;
    mOutPut_ChargedFLow += ".root";
  }
  if(mMode == 3)
  {
    mOutPut_M2_nSigPion = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/Mass2_nSigmaPion/file_%s_M2_nSigPion_%s_etagap_",TriFlow::Energy[energy].Data(),TriFlow::Energy[energy].Data(),TriFlow::Centrality_01[TriFlow::Centrality_start].Data()); 
    mOutPut_M2_nSigPion += TriFlow::EtaGap_start;
    mOutPut_M2_nSigPion += TriFlow::EtaGap_stop-1;
    mOutPut_M2_nSigPion += "_";
    mOutPut_M2_nSigPion += jobCounter;
    mOutPut_M2_nSigPion += ".root";

    mOutPut_Yields_nSigPion = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/Mass2_nSigmaPion/Yields/file_%s_Yields_nSigPion_%s_etagap_",TriFlow::Energy[energy].Data(),TriFlow::Energy[energy].Data(),TriFlow::Centrality_01[TriFlow::Centrality_start].Data()); 
    mOutPut_Yields_nSigPion += TriFlow::EtaGap_start;
    mOutPut_Yields_nSigPion += TriFlow::EtaGap_stop-1;
    mOutPut_Yields_nSigPion += "_";
    mOutPut_Yields_nSigPion += jobCounter;
    mOutPut_Yields_nSigPion += ".root";
  }
  if(mMode == 4)
  {
    mOutPut_M2_Proton = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/Mass2_Proton/file_%s_M2_Proton_%s_etagap_",TriFlow::Energy[energy].Data(),TriFlow::Energy[energy].Data(),TriFlow::Centrality_01[TriFlow::Centrality_start].Data()); 
    mOutPut_M2_Proton += TriFlow::EtaGap_start;
    mOutPut_M2_Proton += TriFlow::EtaGap_stop-1;
    mOutPut_M2_Proton += "_";
    mOutPut_M2_Proton += jobCounter;
    mOutPut_M2_Proton += ".root";

    mOutPut_Yields_Proton = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/Mass2_Proton/Yields/file_%s_Yields_Proton_%s_etagap_",TriFlow::Energy[energy].Data(),TriFlow::Energy[energy].Data(),TriFlow::Centrality_01[TriFlow::Centrality_start].Data()); 
    mOutPut_Yields_Proton += TriFlow::EtaGap_start;
    mOutPut_Yields_Proton += TriFlow::EtaGap_stop-1;
    mOutPut_Yields_Proton += "_";
    mOutPut_Yields_Proton += jobCounter;
    mOutPut_Yields_Proton += ".root";
  }
  if(mMode == 5)
  {
    mOutPut_Phi = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/Phi/file_%s_Phi_%s_",TriFlow::Energy[energy].Data(),TriFlow::Energy[energy].Data(),TriFlow::MixEvent[mFlag_ME].Data()); 
    mOutPut_Phi += jobCounter;
    mOutPut_Phi += ".root";
  }
  if(mMode == 6)
  {
    mOutPut_Lambda = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/Lambda/file_%s_Lambda_%s_",TriFlow::Energy[energy].Data(),TriFlow::Energy[energy].Data(),TriFlow::MixEvent[mFlag_ME].Data()); 
    mOutPut_Lambda += jobCounter;
    mOutPut_Lambda += ".root";
  }
  if(mMode == 7)
  {
    mOutPut_AntiLambda = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/AntiLambda/file_%s_AntiLambda_%s_",TriFlow::Energy[energy].Data(),TriFlow::Energy[energy].Data(),TriFlow::MixEvent[mFlag_ME].Data()); 
    mOutPut_AntiLambda += jobCounter;
    mOutPut_AntiLambda += ".root";
  }
  if(mMode == 8)
  {
    mOutPut_K0S = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/K0S/file_%s_K0S_%s_",TriFlow::Energy[energy].Data(),TriFlow::Energy[energy].Data(),TriFlow::MixEvent[mFlag_ME].Data()); 
    mOutPut_K0S += jobCounter;
    mOutPut_K0S += ".root";
  }
}

//----------------------------------------------------------------------------- 
StTriFlowMaker::~StTriFlowMaker()
{ /*  */ }

//----------------------------------------------------------------------------- 
Int_t StTriFlowMaker::Init() 
{
  if(!mRefMultCorr)
  {
    mRefMultCorr = CentralityMaker::instance()->getRefMultCorr();
  }
  mTriFlowCut = new StTriFlowCut(mEnergy);
  mTriFlowProManger = new StTriFlowProManger();
  mTriFlowCorrection = new StTriFlowCorrection(mEnergy);
  mTriFlowHistoManger = new StTriFlowHistoManger(mEnergy);
  mCombPID = new StCombPID();

  if(mMode == 0)
  {
    mFile_ReCenterPar = new TFile(mOutPut_ReCenterPar.Data(),"RECREATE");
    mFile_ReCenterPar->cd();
    mTriFlowProManger->InitReCenter();
    mTriFlowHistoManger->InitQA_Detector();
  }

  if(mMode == 1)
  {
    mUsedTrackCounter = 0;
    mTriFlowCorrection->InitReCenterCorrection(mEnergy);
    mFile_Corr_Shift = new TFile(mOutPut_Corr_Shift.Data(),"RECREATE");
    mTriFlowProManger->InitShift();
    mFile_Corr_ReCenter = new TFile(mOutPut_Corr_ReCenter.Data(),"RECREATE");
    mFile_Corr_ReCenter->cd();
    mTriFlowCorrection->InitNtuple();
  }
  if(mMode == 2)
  {
    mFile_ChargedFlow = new TFile(mOutPut_ChargedFLow.Data(),"RECREATE");
    mTriFlowCorrection->InitReCenterCorrection(mEnergy);
    mTriFlowCorrection->InitShiftCorrection(mEnergy);
    mTriFlowProManger->InitChargedFlow();
  }
  if(mMode == 3)
  {
    mFile_Yields_nSigPion = new TFile(mOutPut_Yields_nSigPion.Data(),"RECREATE");
    mTriFlowHistoManger->InitYields_nSigPion();
    mFile_M2_nSigPion = new TFile(mOutPut_M2_nSigPion.Data(),"RECREATE");
    mFile_M2_nSigPion->cd();
    mTriFlowHistoManger->InitHist();
    mTriFlowCorrection->InitReCenterCorrection(mEnergy);
    mTriFlowCorrection->InitShiftCorrection(mEnergy);
  }
  if(mMode == 4)
  {
    mFile_Yields_Proton = new TFile(mOutPut_Yields_Proton.Data(),"RECREATE");
    mTriFlowHistoManger->InitYields_Proton();
    mFile_M2_Proton = new TFile(mOutPut_M2_Proton.Data(),"RECREATE");
    mFile_M2_Proton->cd();
    mTriFlowHistoManger->InitProton();
    mTriFlowCorrection->InitReCenterCorrection(mEnergy);
    mTriFlowCorrection->InitShiftCorrection(mEnergy);
  }
  if(mMode == 5)
  {
    mTriFlowV0 = new StTriFlowV0(mEnergy);
    mFile_Phi = new TFile(mOutPut_Phi.Data(),"RECREATE");
    mFile_Phi->cd();
    mTriFlowV0->InitPhi();
    mTriFlowCorrection->InitReCenterCorrection(mEnergy);
    mTriFlowCorrection->InitShiftCorrection(mEnergy);
  }
  if(mMode == 6)
  {
    mTriFlowV0 = new StTriFlowV0(mEnergy);
    mFile_Lambda = new TFile(mOutPut_Lambda.Data(),"RECREATE");
    mFile_Lambda->cd();
    mTriFlowV0->InitLambda();
    mTriFlowV0->SetTopoCut(0.1,0.7,1.0,2.0,1.3,0.7,1.20); // dca_proton, dca_pion, dcaAB, decaylength, dcaV0, dca_pion_pre, InvLambda_low = TriFlow::mMassProton+TriFlow::mMassPion, InvLambda_high
    mTriFlowCorrection->InitReCenterCorrection(mEnergy);
    mTriFlowCorrection->InitShiftCorrection(mEnergy);
  }
  if(mMode == 7)
  {
    mTriFlowV0 = new StTriFlowV0(mEnergy);
    mFile_AntiLambda = new TFile(mOutPut_AntiLambda.Data(),"RECREATE");
    mFile_AntiLambda->cd();
    mTriFlowV0->InitAntiLambda();
    mTriFlowV0->SetTopoCut(0.1,0.7,1.0,2.0,1.3,0.7,1.20); // dca_proton, dca_pion, dcaAB, decaylength, dcaV0, dca_pion_pre, InvLambda_low = TriFlow::mMassProton+TriFlow::mMassPion, InvLambda_high
    mTriFlowCorrection->InitReCenterCorrection(mEnergy);
    mTriFlowCorrection->InitShiftCorrection(mEnergy);
  }
  if(mMode == 8)
  {
    mTriFlowV0 = new StTriFlowV0(mEnergy);
    mFile_K0S = new TFile(mOutPut_K0S.Data(),"RECREATE");
    mFile_K0S->cd();
    mTriFlowV0->InitK0S();
    mTriFlowV0->SetTopoCutK0S(0.7,1.0,2.0,1.3,0.7,0.6); // dca_pion, dcaAB, decaylength, dcaV0, dca_pion_pre, InvK0S_low = TriFlow::mMassPion+TriFlow::mMassPion, InvK0S_high
    mTriFlowCorrection->InitReCenterCorrection(mEnergy);
    mTriFlowCorrection->InitShiftCorrection(mEnergy);
  }

  return kStOK;
}

//----------------------------------------------------------------------------- 
Int_t StTriFlowMaker::Finish() 
{
  if(mMode == 0)
  {
    if(mOutPut_ReCenterPar != "")
    {
      mFile_ReCenterPar->cd();
      mTriFlowProManger->WriteReCenter();
      mTriFlowHistoManger->WriteQA_Detector();
      mFile_ReCenterPar->Close();
    }
  }
  if(mMode == 1)
  {
    if(mOutPut_Corr_ReCenter != "")
    {
      mFile_Corr_ReCenter->cd();
      mTriFlowCorrection->writeNtuple();
      mFile_Corr_ReCenter->Close();
    }
    if(mOutPut_Corr_Shift != "")
    {
      mFile_Corr_Shift->cd();
      mTriFlowProManger->WriteShift();
      mFile_Corr_Shift->Close();
    }
  }
  if(mMode == 2)
  {
    if(mOutPut_ChargedFLow != "")
    {
      mFile_ChargedFlow->cd();
      mTriFlowProManger->WriteChargedFlow();
      mFile_ChargedFlow->Close();
    }
  }
  if(mMode == 3)
  {
    if(mOutPut_M2_nSigPion != "")
    {
      mFile_M2_nSigPion->cd();
      mTriFlowHistoManger->WriteHist();
      mFile_M2_nSigPion->Close();
    }
    if(mOutPut_Yields_nSigPion != "")
    {
      mFile_Yields_nSigPion->cd();
      mTriFlowHistoManger->WriteYileds_nSigPion();
      mFile_Yields_nSigPion->Close();
    }
  }
  if(mMode == 4)
  {
    if(mOutPut_M2_Proton != "")
    {
      mFile_M2_Proton->cd();
      mTriFlowHistoManger->WriteProton();
      mTriFlowHistoManger->WriteQA();
      mFile_M2_Proton->Close();
    }
    if(mOutPut_Yields_Proton != "")
    {
      mFile_Yields_Proton->cd();
      mTriFlowHistoManger->WriteYileds_Proton();
      mFile_Yields_Proton->Close();
    }
  }
  if(mMode == 5)
  {
    if(mOutPut_Phi != "")
    {
      mFile_Phi->cd();
      mTriFlowV0->WritePhiMass2();
      mFile_Phi->Close();
    }
  }
  if(mMode == 6)
  {
    if(mOutPut_Lambda != "")
    {
      mFile_Lambda->cd();
      mTriFlowV0->WriteLambdaMass2();
      mFile_Lambda->Close();
    }
  }
  if(mMode == 7)
  {
    if(mOutPut_AntiLambda != "")
    {
      mFile_AntiLambda->cd();
      mTriFlowV0->WriteAntiLambdaMass2();
      mFile_AntiLambda->Close();
    }
  }
  if(mMode == 8)
  {
    if(mOutPut_K0S != "")
    {
      mFile_K0S->cd();
      mTriFlowV0->WriteK0SMass2();
      mFile_K0S->Close();
    }
  }

  return kStOK;
}

//----------------------------------------------------------------------------- 
void StTriFlowMaker::Clear(Option_t *opt) {
}

//----------------------------------------------------------------------------- 
Int_t StTriFlowMaker::Make() 
{
  if(!mPicoDstMaker) 
  {
    LOG_WARN << " No PicoDstMaker! Skip! " << endm;
    return kStWarn;
  }

  mPicoDst = mPicoDstMaker->picoDst();
  if(!mPicoDst) 
  {
    LOG_WARN << " No PicoDst! Skip! " << endm;
    return kStWarn;
  }

  mPicoEvent = (StPicoEvent*)mPicoDst->event();
  if(!mPicoEvent)
  {
    LOG_WARN << " No PicoEvent! Skip! " << endm;
    return kStWarn;
  }

  // RefMult
  Int_t runId = mPicoEvent->runId();
  Int_t refMult = mPicoEvent->refMult();
  Float_t vz = mPicoEvent->primaryVertex().z();
  Float_t zdcX = mPicoEvent->ZDCx();
  mRefMultCorr->init(runId);
  if(mEnergy == 0) mRefMultCorr->initEvent(refMult,vz,zdcX); // for 200 GeV
  if(mEnergy != 0) mRefMultCorr->initEvent(refMult,vz); // for BES Energy
//  cout << "StTriFLowMaker:" << endl;
//  cout << mRefMultCorr->getWeight() << endl;

  // vz sign
  Int_t vz_sign;
  if(vz > 0.0)
  {
    vz_sign = 0;
  }
  else
  {
    vz_sign = 1;
  }

  // runIndex
  mRunIdEventsDb = StRunIdEventsDb::Instance((Float_t)mPicoEvent->energy(),(Float_t)mPicoEvent->year());
  const Int_t runIndex = mRunIdEventsDb->getRunIdIndex(runId); // expensive
//  cout << runIndex << endl;
//  cout << mRunIdEventsDb->getTotalNrRunIds() << endl;

  // Event Cut
  if(mTriFlowCut->passEventCut(mPicoDst)) // event cut
  {
    const Int_t nTracks = mPicoDst->numberOfTracks();
    const Int_t cent9 = mRefMultCorr->getCentralityBin9();
//    if(cent9 < 0) cout << cent9 << endl;
    const Double_t reweight = mRefMultCorr->getWeight();
    const Int_t nToFMatched = mTriFlowCut->getMatchedToF();

//    cout << "nTracks = " << nTracks << endl;
    for(Int_t i = 0; i < nTracks; i++) // track loop
    {
      StPicoTrack *track = (StPicoTrack*)mPicoDst->track(i);

      if(mMode == 0)
      {
	Float_t eta = track->pMom().pseudoRapidity();
	if(fabs(eta) < TriFlow::mEtaMax && track->dca() < 3.0)
	{
	  Float_t Mass2 = mTriFlowCut->getMass2(track);
	  Float_t dEdx = track->dEdx();
	  Float_t p = track->pMom().mag();
	  mTriFlowHistoManger->FillQA_Detector(dEdx,Mass2,p);
	}
      }
      if(mTriFlowCut->passTrackEP(track)) // track cut
      {
	if(mMode == 0) // fill re-center parameter
	{
	  Float_t pt = track->pMom().perp();

	  if(mTriFlowCorrection->passTrackFull(track))
	  {
	    TVector2 q2Vector_Full = mTriFlowCorrection->calq2Vector(track);
	    TVector2 q3Vector_Full = mTriFlowCorrection->calq3Vector(track);
	    mTriFlowProManger->FillTrackFull(q2Vector_Full,q3Vector_Full,cent9,runIndex,vz_sign,pt);
	  }

	  for(Int_t j = 0; j < 4; j++) // eta_gap loop
	  {
	    if(mTriFlowCorrection->passTrackEtaEast(track,j,0)) // neg eta sub
	    {
	      TVector2 q2Vector_East = mTriFlowCorrection->calq2Vector(track);
	      TVector2 q3Vector_East = mTriFlowCorrection->calq3Vector(track);
	      mTriFlowProManger->FillTrackEast(q2Vector_East,q3Vector_East,cent9,runIndex,vz_sign,j,pt);
	    }
	    if(mTriFlowCorrection->passTrackEtaWest(track,j,0)) // pos eta sub
	    {
	      TVector2 q2Vector_West = mTriFlowCorrection->calq2Vector(track);
	      TVector2 q3Vector_West = mTriFlowCorrection->calq3Vector(track);
	      mTriFlowProManger->FillTrackWest(q2Vector_West,q3Vector_West,cent9,runIndex,vz_sign,j,pt);
	    }
	  }
	}

	if(mMode == 1 || mMode == 2 || mMode == 3 || mMode == 4 || mMode == 5 || mMode == 6 || mMode == 7 || mMode == 8) // calculate Q Vector after recentering for full event and eta sub event
	{
	  if(mTriFlowCorrection->passTrackFull(track))
	  {
	    mTriFlowCorrection->addTrack_Full(track,cent9,runIndex,vz_sign);
	    mUsedTrackCounter++;
	  }
	  for(Int_t j = 0; j < 4; j++) // eta_gap loop
	  {
	    if(mTriFlowCorrection->passTrackEtaEast(track,j,0)) // neg eta sub
	    {
	      mTriFlowCorrection->addTrack_East(track,cent9,runIndex,vz_sign,j);
	    }
	    if(mTriFlowCorrection->passTrackEtaWest(track,j,0)) // pos eta sub
	    {
	      mTriFlowCorrection->addTrack_West(track,cent9,runIndex,vz_sign,j);
	    }
	  }
	}
      }
    }
    if(mMode == 1) // calculate Q vector after recentering for Random Sub Event
    {
      Int_t iTrack[mUsedTrackCounter];
      Float_t ranCounter = (Float_t)mUsedTrackCounter/2.0 - 1;
      for(Int_t i = 0; i < mUsedTrackCounter; i++)
      {
        iTrack[i] = i;
      }
      random_shuffle(iTrack,iTrack+mUsedTrackCounter);
      mUsedTrackCounter = 0;
      for(Int_t i = 0; i < nTracks; i++) // track loop
      {
	StPicoTrack *track = (StPicoTrack*)mPicoDst->track(i);
	if(mTriFlowCut->passTrackEP(track)) // track cut
	{
	  if(mTriFlowCorrection->passTrackFull(track))
	  {
	    if((Float_t)iTrack[mUsedTrackCounter] > ranCounter) // Sub Event A
	    {
	      mTriFlowCorrection->addTrack_A(track,cent9,runIndex,vz_sign);
	    }
	    else // Sub Event B
	    {
	      mTriFlowCorrection->addTrack_B(track,cent9,runIndex,vz_sign);
	    }
	    mUsedTrackCounter++;
	  }
	}
      }
    }
    if(mMode == 1) // re-center and calculate shift parameter for EP and SP
    {
      mTriFlowCorrection->fillNtuple(mPicoDst,cent9,nToFMatched,runIndex);

      // full event shift parameter
      if(mTriFlowCorrection->passTrackFullNumCut())
      {
	for(Int_t k = 0; k < 5; k++) // ShiftOrder loop
	{
	  // Event Plane method
	  TVector2 Psi2Vector_Full_EP = mTriFlowCorrection->calPsi2_Full_EP(k);
	  TVector2 Psi3Vector_Full_EP = mTriFlowCorrection->calPsi3_Full_EP(k);
	  mTriFlowProManger->FillEventFull_EP(Psi2Vector_Full_EP,Psi3Vector_Full_EP,cent9,runIndex,vz_sign,k);

	  // Scalor Product method
	  TVector2 Psi2Vector_Full_SP = mTriFlowCorrection->calPsi2_Full_SP(k);
	  TVector2 Psi3Vector_Full_SP = mTriFlowCorrection->calPsi3_Full_SP(k);
	  mTriFlowProManger->FillEventFull_SP(Psi2Vector_Full_SP,Psi3Vector_Full_SP,cent9,runIndex,vz_sign,k);
	}
      }

      // eta sub event shift parameter
      for(Int_t j = 0; j < 4; j ++)
      {
        if(mTriFlowCorrection->passTrackEtaNumCut(j))
	{
	  for(Int_t k = 0; k < 5; k++)
	  {
	    // Event Plane method
	    TVector2 Psi2Vector_East_EP = mTriFlowCorrection->calPsi2_East_EP(j,k);
	    TVector2 Psi3Vector_East_EP = mTriFlowCorrection->calPsi3_East_EP(j,k);
	    mTriFlowProManger->FillEventEast_EP(Psi2Vector_East_EP,Psi3Vector_East_EP,cent9,runIndex,vz_sign,j,k);

	    TVector2 Psi2Vector_West_EP = mTriFlowCorrection->calPsi2_West_EP(j,k);
	    TVector2 Psi3Vector_West_EP = mTriFlowCorrection->calPsi3_West_EP(j,k);
	    mTriFlowProManger->FillEventWest_EP(Psi2Vector_West_EP,Psi3Vector_West_EP,cent9,runIndex,vz_sign,j,k);

	    // Scalor Product method
	    TVector2 Psi2Vector_East_SP = mTriFlowCorrection->calPsi2_East_SP(j,k);
	    TVector2 Psi3Vector_East_SP = mTriFlowCorrection->calPsi3_East_SP(j,k);
	    mTriFlowProManger->FillEventEast_SP(Psi2Vector_East_SP,Psi3Vector_East_SP,cent9,runIndex,vz_sign,j,k);

	    TVector2 Psi2Vector_West_SP = mTriFlowCorrection->calPsi2_West_SP(j,k);
	    TVector2 Psi3Vector_West_SP = mTriFlowCorrection->calPsi3_West_SP(j,k);
	    mTriFlowProManger->FillEventWest_SP(Psi2Vector_West_SP,Psi3Vector_West_SP,cent9,runIndex,vz_sign,j,k);
	  }
	}
      }
      mTriFlowCorrection->clear();
      mUsedTrackCounter = 0;
    }

    if(mMode == 2) // Charged Flow
    {
      // calculate Event Plane for EP method or QVector for SP method
      // eta sub method
      for(Int_t j = 0; j < 4; j ++)
      {
	if(mTriFlowCorrection->passTrackEtaNumCut(j))
	{
	  // Event Plane method
	  Float_t Psi2_East = mTriFlowCorrection->calShiftAngle2East_EP(runIndex,cent9,vz_sign,j);
	  Float_t Psi2_West = mTriFlowCorrection->calShiftAngle2West_EP(runIndex,cent9,vz_sign,j);
	  Float_t Res2_EP = mTriFlowCorrection->getResolution2_EP(cent9,j);
	  Float_t Psi3_East = mTriFlowCorrection->calShiftAngle3East_EP(runIndex,cent9,vz_sign,j);
	  Float_t Psi3_West = mTriFlowCorrection->calShiftAngle3West_EP(runIndex,cent9,vz_sign,j);
	  Float_t Res3_EP = mTriFlowCorrection->getResolution3_EP(cent9,j);

	  // Scalor Product method
	  TVector2 Q2_East = mTriFlowCorrection->calQVector2East_SP(runIndex,cent9,vz_sign,j);
	  TVector2 Q2_West = mTriFlowCorrection->calQVector2West_SP(runIndex,cent9,vz_sign,j);
	  Float_t Res2_SP = mTriFlowCorrection->getResolution2_SP(cent9,j);
	  TVector2 Q3_East = mTriFlowCorrection->calQVector3East_SP(runIndex,cent9,vz_sign,j);
	  TVector2 Q3_West = mTriFlowCorrection->calQVector3West_SP(runIndex,cent9,vz_sign,j);
	  Float_t Res3_SP = mTriFlowCorrection->getResolution3_SP(cent9,j);

	  for(Int_t i = 0; i < nTracks; i++) // track loop
	  {
	    StPicoTrack *track = (StPicoTrack*)mPicoDst->track(i);
	    if(mTriFlowCut->passTrackCut(track)) // track cut
	    {
	      if(mTriFlowCorrection->passTrackEtaEast(track,j,1))
	      {
		Float_t pt = track->pMom().perp();
	        // Event Plane method
		Float_t phi_East = track->pMom().phi();
		Float_t flow_2_EP = TMath::Cos(2.0*(phi_East-Psi2_West))/Res2_EP;
		Float_t flow_3_EP = TMath::Cos(3.0*(phi_East-Psi3_West))/Res3_EP;
//		Float_t flow_2_EP = TMath::Cos(2.0*(phi_East-Psi2_West));
//		Float_t flow_3_EP = TMath::Cos(3.0*(phi_East-Psi3_West));

	        // Scalar Product method
		TVector2 q2Vector_East_SP = mTriFlowCorrection->calq2Vector(track);
		Float_t flow_2_SP = q2Vector_East_SP*Q2_West/Res2_SP;
//		Float_t flow_2_SP = q2Vector_East_SP*Q2_West/Q2_West.Mod();
		TVector2 q3Vector_East_SP = mTriFlowCorrection->calq3Vector(track);
		Float_t flow_3_SP = q3Vector_East_SP*Q3_West/Res3_SP;
//		Float_t flow_3_SP = q3Vector_East_SP*Q3_West/Q3_West.Mod();

		mTriFlowProManger->FillEtaCharged2Flow(pt,flow_2_EP,Res2_EP,flow_2_SP,Res2_SP,cent9,j,reweight);
		mTriFlowProManger->FillEtaCharged3Flow(pt,flow_3_EP,Res3_EP,flow_3_SP,Res3_SP,cent9,j,reweight);
	      }
	      if(mTriFlowCorrection->passTrackEtaWest(track,j,1))
	      {
		Float_t pt = track->pMom().perp();
	        // Event Plane method
		Float_t phi_West = track->pMom().phi();
		Float_t flow_2_EP = TMath::Cos(2.0*(phi_West-Psi2_East))/Res2_EP;
		Float_t flow_3_EP = TMath::Cos(3.0*(phi_West-Psi3_East))/Res3_EP;
//		Float_t flow_2_EP = TMath::Cos(2.0*(phi_West-Psi2_East));
//		Float_t flow_3_EP = TMath::Cos(3.0*(phi_West-Psi3_East));

	        // Scalar Product method
		TVector2 q2Vector_West_SP = mTriFlowCorrection->calq2Vector(track);
		Float_t flow_2_SP = q2Vector_West_SP*Q2_East/Res2_SP;
//		Float_t flow_2_SP = q2Vector_West_SP*Q2_East/Q2_East.Mod();
		TVector2 q3Vector_West_SP = mTriFlowCorrection->calq3Vector(track);
		Float_t flow_3_SP = q3Vector_West_SP*Q3_East/Res3_SP;
//		Float_t flow_3_SP = q3Vector_West_SP*Q3_East/Q3_East.Mod();

		mTriFlowProManger->FillEtaCharged2Flow(pt,flow_2_EP,Res2_EP,flow_2_SP,Res2_SP,cent9,j,reweight);
		mTriFlowProManger->FillEtaCharged3Flow(pt,flow_3_EP,Res3_EP,flow_3_SP,Res3_SP,cent9,j,reweight);
	      }
	    }
	  }
	}
      }
      // random sub event
      if(mTriFlowCorrection->passTrackFullNumCut())
      {
	Float_t Res2_Full_EP = mTriFlowCorrection->getResolution2_Full_EP(cent9);
	Float_t Res3_Full_EP = mTriFlowCorrection->getResolution3_Full_EP(cent9);
	Float_t Res2_Full_SP = mTriFlowCorrection->getResolution2_Full_SP(cent9);
	Float_t Res3_Full_SP = mTriFlowCorrection->getResolution3_Full_SP(cent9);
	for(Int_t i = 0; i < nTracks; i++)
	{
	  StPicoTrack *track = (StPicoTrack*)mPicoDst->track(i);
	  if(mTriFlowCut->passTrackCut(track)) // track cut
	  {
	    Float_t pt = track->pMom().perp();
	    // Event Plane method
	    Float_t phi = track->pMom().phi();

	    Float_t Psi2 = mTriFlowCorrection->calShiftAngle2Full_EP(runIndex,cent9,vz_sign,track);
	    Float_t flow_2_full_EP = TMath::Cos(2.0*(phi-Psi2))/Res2_Full_EP;
//	    Float_t flow_2_full_EP = TMath::Cos(2.0*(phi-Psi2));

	    Float_t Psi3 = mTriFlowCorrection->calShiftAngle3Full_EP(runIndex,cent9,vz_sign,track);
	    Float_t flow_3_full_EP = TMath::Cos(3.0*(phi-Psi3))/Res3_Full_EP;
//	    Float_t flow_3_full_EP = TMath::Cos(3.0*(phi-Psi3));

	    // Scalar Product method
	    TVector2 q2Vector_Full_SP = mTriFlowCorrection->calq2Vector(track);
	    TVector2 Q2Vector_Full_SP = mTriFlowCorrection->calQVector2Full_SP(runIndex,cent9,vz_sign,track);
	    Float_t flow_2_full_SP = q2Vector_Full_SP*Q2Vector_Full_SP/Res2_Full_SP;
//	    Float_t flow_2_full_SP = q2Vector_Full_SP*Q2Vector_Full_SP/Q2Vector_Full_SP.Mod();

	    TVector2 q3Vector_Full_SP = mTriFlowCorrection->calq3Vector(track);
	    TVector2 Q3Vector_Full_SP = mTriFlowCorrection->calQVector3Full_SP(runIndex,cent9,vz_sign,track);
	    Float_t flow_3_full_SP = q3Vector_Full_SP*Q3Vector_Full_SP/Res3_Full_SP;
//	    Float_t flow_3_full_SP = q3Vector_Full_SP*Q3Vector_Full_SP/Q3Vector_Full_SP.Mod();

	    mTriFlowProManger->FillRanCharged2Flow(pt,flow_2_full_EP,Res2_Full_EP,flow_2_full_SP,Res2_Full_SP,cent9,reweight);
	    mTriFlowProManger->FillRanCharged3Flow(pt,flow_3_full_EP,Res3_Full_EP,flow_3_full_SP,Res3_Full_SP,cent9,reweight);
	  }
	}
      }
      mTriFlowCorrection->clear();
    }

    if(mMode == 3)
    {
      // particle identification for pion and kaon
      // eta sub method
      for(Int_t j = TriFlow::EtaGap_start; j < TriFlow::EtaGap_stop; j ++)
      {
	if(mTriFlowCorrection->passTrackEtaNumCut(j))
	{
	  // Event Plane method
	  Float_t Psi2_East = mTriFlowCorrection->calShiftAngle2East_EP(runIndex,cent9,vz_sign,j);
	  Float_t Psi2_West = mTriFlowCorrection->calShiftAngle2West_EP(runIndex,cent9,vz_sign,j);
	  Float_t Res2_EP = mTriFlowCorrection->getResolution2_EP(cent9,j);
	  Float_t Psi3_East = mTriFlowCorrection->calShiftAngle3East_EP(runIndex,cent9,vz_sign,j);
	  Float_t Psi3_West = mTriFlowCorrection->calShiftAngle3West_EP(runIndex,cent9,vz_sign,j);
	  Float_t Res3_EP = mTriFlowCorrection->getResolution3_EP(cent9,j);

	  Int_t i_cut = 0;

	  for(Int_t i_dca = 0; i_dca < 3; i_dca++)
	  {
	    for(Int_t i_nHitsFit = 0; i_nHitsFit < 2; i_nHitsFit++)
	    {
	      for(Int_t i = 0; i < nTracks; i++) // track loop
	      {
		StPicoTrack *track = (StPicoTrack*)mPicoDst->track(i);

		if(mTriFlowCut->passTrackCutSys(track,i_dca,i_nHitsFit)) // track cut
		{
		  if(mTriFlowCut->passPIDCut(track))
		  {
		    Float_t pt = track->pMom().perp();
		    Float_t Mass2 = mTriFlowCut->getMass2(track);
		    Float_t nSigmaPion = track->nSigmaPion();
		    Float_t scale_nSigma_factor = TriFlow::mSigScaleMap[mPicoEvent->energy()];
		    Int_t nCharge = track->charge();
		    Int_t charge_bin;
		    if(nCharge > 0)
		    {
		      charge_bin = 0;
		    }
		    if(nCharge < 0)
		    {
		      charge_bin = 1;
		    }

		    Float_t ex_scale_factor = TriFlow::mExScaleMap[mPicoEvent->energy()];
		    mCombPID->setInitValues(pt,nSigmaPion*scale_nSigma_factor,Mass2,ex_scale_factor);
		    Float_t New_X = mCombPID->getNewX();
		    Float_t New_Y = mCombPID->getNewY();

		    // Fill Histogram
		    if(mTriFlowCorrection->passTrackEtaEast(track,j,1))
		    {
		      Float_t phi_East = track->pMom().phi();
		      Float_t phi_Psi2 = phi_East - Psi2_West;
		      Float_t phi_Psi3 = phi_East - Psi3_West;
		      //		  cout << "phi = " << phi_East << endl;
		      //		  cout << "psi = " << Psi3_West << endl;
		      mTriFlowHistoManger->FillHist(pt,cent9,charge_bin,j,phi_Psi2,Res2_EP,phi_Psi3,Res3_EP,New_X,New_Y,reweight,i_cut);
		    }
		    if(mTriFlowCorrection->passTrackEtaWest(track,j,1))
		    {
		      Float_t phi_West = track->pMom().phi();
		      Float_t phi_Psi2 = phi_West - Psi2_East;
		      Float_t phi_Psi3 = phi_West - Psi3_East;
		      //		  cout << "phi = " << phi_West << endl;
		      //		  cout << "psi = " << Psi3_East << endl;
		      mTriFlowHistoManger->FillHist(pt,cent9,charge_bin,j,phi_Psi2,Res2_EP,phi_Psi3,Res3_EP,New_X,New_Y,reweight,i_cut);
		    }
		    mTriFlowHistoManger->FillYields_PiK(cent9,charge_bin,j,New_X,New_Y,reweight,i_cut);
		  }
		  mTriFlowHistoManger->FillToFLocal(track);
		}
	      }
	      i_cut++;
	    }
	  }
	}
      }
      mTriFlowCorrection->clear();
    }

    if(mMode == 4)
    {
      // particle identification for proton
      // eta sub method
      for(Int_t j = TriFlow::EtaGap_start; j < TriFlow::EtaGap_stop; j++)
      {
	if(mTriFlowCorrection->passTrackEtaNumCut(j))
	{
	  // Event Plane method
	  Float_t Psi2_East = mTriFlowCorrection->calShiftAngle2East_EP(runIndex,cent9,vz_sign,j);
	  Float_t Psi2_West = mTriFlowCorrection->calShiftAngle2West_EP(runIndex,cent9,vz_sign,j);
	  Float_t Res2_EP = mTriFlowCorrection->getResolution2_EP(cent9,j);
	  Float_t Psi3_East = mTriFlowCorrection->calShiftAngle3East_EP(runIndex,cent9,vz_sign,j);
	  Float_t Psi3_West = mTriFlowCorrection->calShiftAngle3West_EP(runIndex,cent9,vz_sign,j);
	  Float_t Res3_EP = mTriFlowCorrection->getResolution3_EP(cent9,j);

	  Int_t i_cut = 0;

	  for(Int_t i_dca = 0; i_dca < 3; i_dca++)
	  {
	    for(Int_t i_nHitsFit = 0; i_nHitsFit < 2; i_nHitsFit++)
	    {
	      for(Int_t i_proton = 0; i_proton < 3; i_proton++)
	      {
		for(Int_t i = 0; i < nTracks; i++) // track loop
		{
		  StPicoTrack *track = (StPicoTrack*)mPicoDst->track(i);

		  if(mTriFlowCut->passTrackCutSys(track,i_dca,i_nHitsFit)) // track cut
		  {
		    if(mTriFlowCut->passPIDCut(track))
		    {
		      Float_t pt = track->pMom().perp();
		      Float_t p = track->pMom().mag();
		      Float_t Mass2 = mTriFlowCut->getMass2(track);
		      Float_t dEdx = track->dEdx();
		      Int_t nCharge = track->charge();
		      Int_t charge_bin;
		      if(nCharge > 0)
		      {
			charge_bin = 0;
		      }
		      if(nCharge < 0)
		      {
			charge_bin = 1;
		      }

		      mTriFlowHistoManger->FillQA_before(j,Mass2,dEdx,p*nCharge);

		      Float_t scale_nSigma_factor = TriFlow::mSigScaleMap[mPicoEvent->energy()];
		      if(mTriFlowCut->passSigProntonCutSys(track,scale_nSigma_factor,i_proton)) // nSigmaProton cut
		      {
			if(mTriFlowCorrection->passTrackEtaEast(track,j,1))
			{
			  Float_t phi_East = track->pMom().phi();
			  Float_t phi_Psi2 = phi_East - Psi2_West;
			  Float_t phi_Psi3 = phi_East - Psi3_West;
			  //		  cout << "phi = " << phi_East << endl;
			  //		  cout << "psi = " << Psi3_West << endl;
			  mTriFlowHistoManger->FillProton(pt,cent9,charge_bin,j,phi_Psi2,Res2_EP,phi_Psi3,Res3_EP,Mass2,reweight,i_cut);
			}
			if(mTriFlowCorrection->passTrackEtaWest(track,j,1))
			{
			  Float_t phi_West = track->pMom().phi();
			  Float_t phi_Psi2 = phi_West - Psi2_East;
			  Float_t phi_Psi3 = phi_West - Psi3_East;
			  //		  cout << "phi = " << phi_West << endl;
			  //		  cout << "psi = " << Psi3_East << endl;
			  mTriFlowHistoManger->FillProton(pt,cent9,charge_bin,j,phi_Psi2,Res2_EP,phi_Psi3,Res3_EP,Mass2,reweight,i_cut);
			}
			mTriFlowHistoManger->FillQA_after(j,Mass2,dEdx,p*nCharge);
			mTriFlowHistoManger->FillYields_Proton(cent9,charge_bin,j,Mass2,reweight,i_cut);
		      }
		    }
		  }
		}
		i_cut++;
	      }
	    }
	  }
	}
      }
      mTriFlowCorrection->clear();
      mTriFlowHistoManger->FillQA_Event(refMult,vz);
    }

    if(mMode == 5)
    { // phi meson
      Float_t Psi2_East;
      TVector2 Q2East[TriFlow::EtaGap_total], Q2West[TriFlow::EtaGap_total];
      TVector2 Q3East[TriFlow::EtaGap_total], Q3West[TriFlow::EtaGap_total]; 
      Int_t NumTrackEast[TriFlow::EtaGap_total], NumTrackWest[TriFlow::EtaGap_total]; 
      for(Int_t j = 0; j < TriFlow::EtaGap_total; j++)
      {
	Q2East[j].Set(-999.9,-999.9); // initialize Q Vector to unreasonable value
	Q2West[j].Set(-999.9,-999.9);
	Q3East[j].Set(-999.9,-999.9);
	Q3West[j].Set(-999.9,-999.9);
	NumTrackEast[j] = 0;
	NumTrackWest[j] = 0;
	if(mTriFlowCorrection->passTrackEtaNumCut(j))
	{
	  // get QVector of sub event
	  Q2East[j] = mTriFlowCorrection->getQVector(j,0,0); // 0 = eta_gap, 1 = flow type, 2 = east/west
	  Q2West[j] = mTriFlowCorrection->getQVector(j,0,1); // 0 = eta_gap, 1 = flow type, 2 = east/west
	  Q3East[j] = mTriFlowCorrection->getQVector(j,1,0); // 0 = eta_gap, 1 = flow type, 2 = east/west
	  Q3West[j] = mTriFlowCorrection->getQVector(j,1,1); // 0 = eta_gap, 1 = flow type, 2 = east/west
	  NumTrackEast[j] = mTriFlowCorrection->getNumTrack(j,0); // 0 = eta_gap, 1 = east/west
	  NumTrackWest[j] = mTriFlowCorrection->getNumTrack(j,1); // 0 = eta_gap, 1 = east/west
	}
      }

      if(mTriFlowCorrection->passTrackEtaNumCut(0))
      {
	// Event Plane method
	Psi2_East = mTriFlowCorrection->calShiftAngle2East_EP(runIndex,cent9,vz_sign,0);

	// get N_prim, N_non_prim, N_Tof_match
	Int_t N_prim = mTriFlowCut->getNpirm();
	Int_t N_non_prim = mTriFlowCut->getNnonprim();
	Int_t N_Tof_match = mTriFlowCut->getMatchedToF();

	// pass the event information to StTriFlowV0
	mTriFlowV0->clearEvent();
	mTriFlowV0->passEvent(N_prim, N_non_prim, N_Tof_match);

	// 2nd sub event plane
	mTriFlowV0->passEventPlane2East(Q2East[0],Q2East[1],Q2East[2],Q2East[3]);
	mTriFlowV0->passEventPlane2West(Q2West[0],Q2West[1],Q2West[2],Q2West[3]);

	// 3rd sub event plane
	mTriFlowV0->passEventPlane3East(Q3East[0],Q3East[1],Q3East[2],Q3East[3]);
	mTriFlowV0->passEventPlane3West(Q3West[0],Q3West[1],Q3West[2],Q3West[3]);

	// Number of Track in East and West part of TPC
	mTriFlowV0->passNumTrackEast(NumTrackEast[0],NumTrackEast[1],NumTrackEast[2],NumTrackEast[3]);
	mTriFlowV0->passNumTrackWest(NumTrackWest[0],NumTrackWest[1],NumTrackWest[2],NumTrackWest[3]);

	mTriFlowV0->MixEvent_Phi(mFlag_ME,mPicoDst,cent9,vz,Psi2_East);
      }
      mTriFlowCorrection->clear();
    }

    if(mMode == 6 || mMode == 7)
    { // Lambda and anti-Lambda
      Float_t Psi2_East;
      TVector2 Q2East[TriFlow::EtaGap_total], Q2West[TriFlow::EtaGap_total];
      TVector2 Q3East[TriFlow::EtaGap_total], Q3West[TriFlow::EtaGap_total]; 
      Int_t NumTrackEast[TriFlow::EtaGap_total], NumTrackWest[TriFlow::EtaGap_total]; 
      for(Int_t j = 0; j < TriFlow::EtaGap_total; j++)
      {
	Q2East[j].Set(-999.9,-999.9); // initialize Q Vector to unreasonable value
	Q2West[j].Set(-999.9,-999.9);
	Q3East[j].Set(-999.9,-999.9);
	Q3West[j].Set(-999.9,-999.9);
	NumTrackEast[j] = 0;
	NumTrackWest[j] = 0;
	if(mTriFlowCorrection->passTrackEtaNumCut(j))
	{
	  // get QVector of sub event
	  Q2East[j] = mTriFlowCorrection->getQVector(j,0,0); // 0 = eta_gap, 1 = flow type, 2 = east/west
	  Q2West[j] = mTriFlowCorrection->getQVector(j,0,1); // 0 = eta_gap, 1 = flow type, 2 = east/west
	  Q3East[j] = mTriFlowCorrection->getQVector(j,1,0); // 0 = eta_gap, 1 = flow type, 2 = east/west
	  Q3West[j] = mTriFlowCorrection->getQVector(j,1,1); // 0 = eta_gap, 1 = flow type, 2 = east/west
	  NumTrackEast[j] = mTriFlowCorrection->getNumTrack(j,0); // 0 = eta_gap, 1 = east/west
	  NumTrackWest[j] = mTriFlowCorrection->getNumTrack(j,1); // 0 = eta_gap, 1 = east/west
	}
      }

      if(mTriFlowCorrection->passTrackEtaNumCut(0))
      {
	// Event Plane method
	Psi2_East = mTriFlowCorrection->calShiftAngle2East_EP(runIndex,cent9,vz_sign,0);

	// get N_prim, N_non_prim, N_Tof_match
	Int_t N_prim = mTriFlowCut->getNpirm();
	Int_t N_non_prim = mTriFlowCut->getNnonprim();
	Int_t N_Tof_match = mTriFlowCut->getMatchedToF();

	// pass the event information to StTriFlowV0
	mTriFlowV0->clearEvent();
	mTriFlowV0->passEvent(N_prim, N_non_prim, N_Tof_match);

	// 2nd sub event plane
	mTriFlowV0->passEventPlane2East(Q2East[0],Q2East[1],Q2East[2],Q2East[3]);
	mTriFlowV0->passEventPlane2West(Q2West[0],Q2West[1],Q2West[2],Q2West[3]);

	// 3rd sub event plane
	mTriFlowV0->passEventPlane3East(Q3East[0],Q3East[1],Q3East[2],Q3East[3]);
	mTriFlowV0->passEventPlane3West(Q3West[0],Q3West[1],Q3West[2],Q3West[3]);

	// Number of Track in East and West part of TPC
	mTriFlowV0->passNumTrackEast(NumTrackEast[0],NumTrackEast[1],NumTrackEast[2],NumTrackEast[3]);
	mTriFlowV0->passNumTrackWest(NumTrackWest[0],NumTrackWest[1],NumTrackWest[2],NumTrackWest[3]);

//	cout << "cent9 = " << cent9 << endl;
	if(mMode == 6) mTriFlowV0->MixEvent_Lambda(mFlag_ME,mPicoDst,cent9,vz,Psi2_East);
	if(mMode == 7) mTriFlowV0->MixEvent_AntiLambda(mFlag_ME,mPicoDst,cent9,vz,Psi2_East);
      }
      mTriFlowCorrection->clear();
    }
    if(mMode == 8)
    { // K0S
      Float_t Psi2_East;
      TVector2 Q2East[TriFlow::EtaGap_total], Q2West[TriFlow::EtaGap_total];
      TVector2 Q3East[TriFlow::EtaGap_total], Q3West[TriFlow::EtaGap_total]; 
      Int_t NumTrackEast[TriFlow::EtaGap_total], NumTrackWest[TriFlow::EtaGap_total]; 
      for(Int_t j = 0; j < TriFlow::EtaGap_total; j++)
      {
	Q2East[j].Set(-999.9,-999.9); // initialize Q Vector to unreasonable value
	Q2West[j].Set(-999.9,-999.9);
	Q3East[j].Set(-999.9,-999.9);
	Q3West[j].Set(-999.9,-999.9);
	NumTrackEast[j] = 0;
	NumTrackWest[j] = 0;
	if(mTriFlowCorrection->passTrackEtaNumCut(j))
	{
	  // get QVector of sub event
	  Q2East[j] = mTriFlowCorrection->getQVector(j,0,0); // 0 = eta_gap, 1 = flow type, 2 = east/west
	  Q2West[j] = mTriFlowCorrection->getQVector(j,0,1); // 0 = eta_gap, 1 = flow type, 2 = east/west
	  Q3East[j] = mTriFlowCorrection->getQVector(j,1,0); // 0 = eta_gap, 1 = flow type, 2 = east/west
	  Q3West[j] = mTriFlowCorrection->getQVector(j,1,1); // 0 = eta_gap, 1 = flow type, 2 = east/west
	  NumTrackEast[j] = mTriFlowCorrection->getNumTrack(j,0); // 0 = eta_gap, 1 = east/west
	  NumTrackWest[j] = mTriFlowCorrection->getNumTrack(j,1); // 0 = eta_gap, 1 = east/west
	}
      }

      if(mTriFlowCorrection->passTrackEtaNumCut(0))
      {
	// Event Plane method
	Psi2_East = mTriFlowCorrection->calShiftAngle2East_EP(runIndex,cent9,vz_sign,0);

	// get N_prim, N_non_prim, N_Tof_match
	Int_t N_prim = mTriFlowCut->getNpirm();
	Int_t N_non_prim = mTriFlowCut->getNnonprim();
	Int_t N_Tof_match = mTriFlowCut->getMatchedToF();

	// pass the event information to StTriFlowV0
	mTriFlowV0->clearEvent();
	mTriFlowV0->passEvent(N_prim, N_non_prim, N_Tof_match);

	// 2nd sub event plane
	mTriFlowV0->passEventPlane2East(Q2East[0],Q2East[1],Q2East[2],Q2East[3]);
	mTriFlowV0->passEventPlane2West(Q2West[0],Q2West[1],Q2West[2],Q2West[3]);

	// 3rd sub event plane
	mTriFlowV0->passEventPlane3East(Q3East[0],Q3East[1],Q3East[2],Q3East[3]);
	mTriFlowV0->passEventPlane3West(Q3West[0],Q3West[1],Q3West[2],Q3West[3]);

	// Number of Track in East and West part of TPC
	mTriFlowV0->passNumTrackEast(NumTrackEast[0],NumTrackEast[1],NumTrackEast[2],NumTrackEast[3]);
	mTriFlowV0->passNumTrackWest(NumTrackWest[0],NumTrackWest[1],NumTrackWest[2],NumTrackWest[3]);

//	cout << "cent9 = " << cent9 << endl;
	mTriFlowV0->MixEvent_K0S(mFlag_ME,mPicoDst,cent9,vz,Psi2_East);
      }
      mTriFlowCorrection->clear();
    }
  }

  return kStOK;
}


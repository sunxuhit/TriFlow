#include "StTriFlowResolution.h"
#include "StTriFlowConstants.h"
#include "StTriFlowHistManger.h"
#include "StMessMgr.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TMath.h"
#include "TString.h"
#include "TProfile2D.h"
#include "TProfile.h"

ClassImp(StTriFlowResolution)

TString StTriFlowResolution::mVStr[2] = {"pos","neg"};

//----------------------------------------------------------------------------------------------------------

StTriFlowResolution::StTriFlowResolution(const Int_t jobCounter, const Int_t energy)
{
  mEnergy = energy;
  mInPut_Corr_ReCenter = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/Correction/ReCenter/file_%s_Corr_ReCenter_%d.root",TriFlow::Energy[energy].Data(),TriFlow::Energy[energy].Data(),jobCounter);

  mOutPut_Resolution = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/Resolution/file_%s_Resolution_%d.root",TriFlow::Energy[energy].Data(),TriFlow::Energy[energy].Data(),jobCounter);

  mOutPut_EventPlane = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/EventPlane/file_%s_EventPlane_%d.root",TriFlow::Energy[energy].Data(),TriFlow::Energy[energy].Data(),jobCounter);
}

StTriFlowResolution::~StTriFlowResolution()
{
}

//----------------------------------------------------------------------------------------------------------

Int_t StTriFlowResolution::Init()
{
  mFile_Resolution = new TFile(mOutPut_Resolution.Data(),"RECREATE");
  InitReCenterNutple(mInPut_Corr_ReCenter);
  InitShiftCorrection(mEnergy);
  InitResolution();

  mFile_EventPlane = new TFile(mOutPut_EventPlane.Data(),"RECREATE");
  mTriFlowHistManger = new StTriFlowHistManger();
  mTriFlowHistManger->InitEventPlane();

  return kTRUE;
}

//----------------------------------------------------------------------------------------------------------

Int_t StTriFlowResolution::Finish()
{
  if(mOutPut_Resolution != "")
  {
    mFile_Resolution->cd();
    WriteResolution();
    mFile_Resolution->Close();
    mInPutFile_Shift->Close();
  }
  if(mOutPut_EventPlane != "")
  {
    mFile_EventPlane->cd();
    mTriFlowHistManger->WriteEventPlane();
    mFile_EventPlane->Close();
  }

  return kFALSE;
}

//----------------------------------------------------------------------------------------------------------

void StTriFlowResolution::InitReCenterNutple(TString InPutName)
{
  TFile *InPutFile = TFile::Open(InPutName.Data());

  mNtuple = (TNtuple*)InPutFile->Get("Ntuple");

  // initialize Ntuple
  mNtuple->SetBranchAddress("runId",&mRunId);
  mNtuple->SetBranchAddress("eventId",&mEventId);
  mNtuple->SetBranchAddress("RefMult",&mRefMult);
  mNtuple->SetBranchAddress("ZDCx",&mZDCx);
  mNtuple->SetBranchAddress("BBCx",&mBBCx);
  mNtuple->SetBranchAddress("vzVpd",&mVzVpd);
  mNtuple->SetBranchAddress("Centrality9",&mCentrality9);
  mNtuple->SetBranchAddress("vx",&mVx);
  mNtuple->SetBranchAddress("vy",&mVy);
  mNtuple->SetBranchAddress("vz",&mVz);
  mNtuple->SetBranchAddress("nToFMatched",&mNToFMatched);

  // Event Plane method
  // 2nd East
  mNtuple->SetBranchAddress("mQ2X_East_EP_0",&mQ2X_East_EP[0]);
  mNtuple->SetBranchAddress("mQ2X_East_EP_1",&mQ2X_East_EP[1]);
  mNtuple->SetBranchAddress("mQ2X_East_EP_2",&mQ2X_East_EP[2]);
  mNtuple->SetBranchAddress("mQ2X_East_EP_3",&mQ2X_East_EP[3]);

  mNtuple->SetBranchAddress("mQ2Y_East_EP_0",&mQ2Y_East_EP[0]);
  mNtuple->SetBranchAddress("mQ2Y_East_EP_1",&mQ2Y_East_EP[1]);
  mNtuple->SetBranchAddress("mQ2Y_East_EP_2",&mQ2Y_East_EP[2]);
  mNtuple->SetBranchAddress("mQ2Y_East_EP_3",&mQ2Y_East_EP[3]);
  // 3rd East
  mNtuple->SetBranchAddress("mQ3X_East_EP_0",&mQ3X_East_EP[0]);
  mNtuple->SetBranchAddress("mQ3X_East_EP_1",&mQ3X_East_EP[1]);
  mNtuple->SetBranchAddress("mQ3X_East_EP_2",&mQ3X_East_EP[2]);
  mNtuple->SetBranchAddress("mQ3X_East_EP_3",&mQ3X_East_EP[3]);

  mNtuple->SetBranchAddress("mQ3Y_East_EP_0",&mQ3Y_East_EP[0]);
  mNtuple->SetBranchAddress("mQ3Y_East_EP_1",&mQ3Y_East_EP[1]);
  mNtuple->SetBranchAddress("mQ3Y_East_EP_2",&mQ3Y_East_EP[2]);
  mNtuple->SetBranchAddress("mQ3Y_East_EP_3",&mQ3Y_East_EP[3]);
  // 2nd West
  mNtuple->SetBranchAddress("mQ2X_West_EP_0",&mQ2X_West_EP[0]);
  mNtuple->SetBranchAddress("mQ2X_West_EP_1",&mQ2X_West_EP[1]);
  mNtuple->SetBranchAddress("mQ2X_West_EP_2",&mQ2X_West_EP[2]);
  mNtuple->SetBranchAddress("mQ2X_West_EP_3",&mQ2X_West_EP[3]);

  mNtuple->SetBranchAddress("mQ2Y_West_EP_0",&mQ2Y_West_EP[0]);
  mNtuple->SetBranchAddress("mQ2Y_West_EP_1",&mQ2Y_West_EP[1]);
  mNtuple->SetBranchAddress("mQ2Y_West_EP_2",&mQ2Y_West_EP[2]);
  mNtuple->SetBranchAddress("mQ2Y_West_EP_3",&mQ2Y_West_EP[3]);
  // 3rd West
  mNtuple->SetBranchAddress("mQ3X_West_EP_0",&mQ3X_West_EP[0]);
  mNtuple->SetBranchAddress("mQ3X_West_EP_1",&mQ3X_West_EP[1]);
  mNtuple->SetBranchAddress("mQ3X_West_EP_2",&mQ3X_West_EP[2]);
  mNtuple->SetBranchAddress("mQ3X_West_EP_3",&mQ3X_West_EP[3]);

  mNtuple->SetBranchAddress("mQ3Y_West_EP_0",&mQ3Y_West_EP[0]);
  mNtuple->SetBranchAddress("mQ3Y_West_EP_1",&mQ3Y_West_EP[1]);
  mNtuple->SetBranchAddress("mQ3Y_West_EP_2",&mQ3Y_West_EP[2]);
  mNtuple->SetBranchAddress("mQ3Y_West_EP_3",&mQ3Y_West_EP[3]);

  // 2nd full event plane
  mNtuple->SetBranchAddress("mQ2X_Full_EP",&mQ2X_Full_EP);
  mNtuple->SetBranchAddress("mQ2Y_Full_EP",&mQ2Y_Full_EP);
  // 3rd full event plane
  mNtuple->SetBranchAddress("mQ3X_Full_EP",&mQ3X_Full_EP);
  mNtuple->SetBranchAddress("mQ3Y_Full_EP",&mQ3Y_Full_EP);

  // Scalor Product method
  // 2nd East
  mNtuple->SetBranchAddress("mQ2X_East_SP_0",&mQ2X_East_SP[0]);
  mNtuple->SetBranchAddress("mQ2X_East_SP_1",&mQ2X_East_SP[1]);
  mNtuple->SetBranchAddress("mQ2X_East_SP_2",&mQ2X_East_SP[2]);
  mNtuple->SetBranchAddress("mQ2X_East_SP_3",&mQ2X_East_SP[3]);

  mNtuple->SetBranchAddress("mQ2Y_East_SP_0",&mQ2Y_East_SP[0]);
  mNtuple->SetBranchAddress("mQ2Y_East_SP_1",&mQ2Y_East_SP[1]);
  mNtuple->SetBranchAddress("mQ2Y_East_SP_2",&mQ2Y_East_SP[2]);
  mNtuple->SetBranchAddress("mQ2Y_East_SP_3",&mQ2Y_East_SP[3]);
  // 3rd East
  mNtuple->SetBranchAddress("mQ3X_East_SP_0",&mQ3X_East_SP[0]);
  mNtuple->SetBranchAddress("mQ3X_East_SP_1",&mQ3X_East_SP[1]);
  mNtuple->SetBranchAddress("mQ3X_East_SP_2",&mQ3X_East_SP[2]);
  mNtuple->SetBranchAddress("mQ3X_East_SP_3",&mQ3X_East_SP[3]);

  mNtuple->SetBranchAddress("mQ3Y_East_SP_0",&mQ3Y_East_SP[0]);
  mNtuple->SetBranchAddress("mQ3Y_East_SP_1",&mQ3Y_East_SP[1]);
  mNtuple->SetBranchAddress("mQ3Y_East_SP_2",&mQ3Y_East_SP[2]);
  mNtuple->SetBranchAddress("mQ3Y_East_SP_3",&mQ3Y_East_SP[3]);
  // 2nd West
  mNtuple->SetBranchAddress("mQ2X_West_SP_0",&mQ2X_West_SP[0]);
  mNtuple->SetBranchAddress("mQ2X_West_SP_1",&mQ2X_West_SP[1]);
  mNtuple->SetBranchAddress("mQ2X_West_SP_2",&mQ2X_West_SP[2]);
  mNtuple->SetBranchAddress("mQ2X_West_SP_3",&mQ2X_West_SP[3]);

  mNtuple->SetBranchAddress("mQ2Y_West_SP_0",&mQ2Y_West_SP[0]);
  mNtuple->SetBranchAddress("mQ2Y_West_SP_1",&mQ2Y_West_SP[1]);
  mNtuple->SetBranchAddress("mQ2Y_West_SP_2",&mQ2Y_West_SP[2]);
  mNtuple->SetBranchAddress("mQ2Y_West_SP_3",&mQ2Y_West_SP[3]);
  // 3rd West
  mNtuple->SetBranchAddress("mQ3X_West_SP_0",&mQ3X_West_SP[0]);
  mNtuple->SetBranchAddress("mQ3X_West_SP_1",&mQ3X_West_SP[1]);
  mNtuple->SetBranchAddress("mQ3X_West_SP_2",&mQ3X_West_SP[2]);
  mNtuple->SetBranchAddress("mQ3X_West_SP_3",&mQ3X_West_SP[3]);

  mNtuple->SetBranchAddress("mQ3Y_West_SP_0",&mQ3Y_West_SP[0]);
  mNtuple->SetBranchAddress("mQ3Y_West_SP_1",&mQ3Y_West_SP[1]);
  mNtuple->SetBranchAddress("mQ3Y_West_SP_2",&mQ3Y_West_SP[2]);
  mNtuple->SetBranchAddress("mQ3Y_West_SP_3",&mQ3Y_West_SP[3]);

  // 2nd full event plane
  mNtuple->SetBranchAddress("mQ2X_Full_SP",&mQ2X_Full_SP);
  mNtuple->SetBranchAddress("mQ2Y_Full_SP",&mQ2Y_Full_SP);
  // 3rd full event plane
  mNtuple->SetBranchAddress("mQ3X_Full_SP",&mQ3X_Full_SP);
  mNtuple->SetBranchAddress("mQ3Y_Full_SP",&mQ3Y_Full_SP);

  // Counter
  // East
  mNtuple->SetBranchAddress("mQCounter_East_0",&mQCounter_East[0]);
  mNtuple->SetBranchAddress("mQCounter_East_1",&mQCounter_East[1]);
  mNtuple->SetBranchAddress("mQCounter_East_2",&mQCounter_East[2]);
  mNtuple->SetBranchAddress("mQCounter_East_3",&mQCounter_East[3]);
  // West
  mNtuple->SetBranchAddress("mQCounter_West_0",&mQCounter_West[0]);
  mNtuple->SetBranchAddress("mQCounter_West_1",&mQCounter_West[1]);
  mNtuple->SetBranchAddress("mQCounter_West_2",&mQCounter_West[2]);
  mNtuple->SetBranchAddress("mQCounter_West_3",&mQCounter_West[3]);
  // full
  mNtuple->SetBranchAddress("mQCounter_Full",&mQCounter_Full);
  mNtuple->SetBranchAddress("mQCounter_Full_East",&mQCounter_Full_East);
  mNtuple->SetBranchAddress("mQCounter_Full_West",&mQCounter_Full_West);

  mNtuple->SetBranchAddress("runIndex",&mRunIndex);

  // random sub event plane
  mNtuple->SetBranchAddress("mQ2X_A_EP",&mQ2X_A_EP);
  mNtuple->SetBranchAddress("mQ2Y_A_EP",&mQ2Y_A_EP);
  mNtuple->SetBranchAddress("mQ3X_A_EP",&mQ3X_A_EP);
  mNtuple->SetBranchAddress("mQ3Y_A_EP",&mQ3Y_A_EP);
  mNtuple->SetBranchAddress("mQ2X_A_SP",&mQ2X_A_SP);
  mNtuple->SetBranchAddress("mQ2Y_A_SP",&mQ2Y_A_SP);
  mNtuple->SetBranchAddress("mQ3X_A_SP",&mQ3X_A_SP);
  mNtuple->SetBranchAddress("mQ3Y_A_SP",&mQ3Y_A_SP);
  mNtuple->SetBranchAddress("mQCounter_A",&mQCounter_A);

  mNtuple->SetBranchAddress("mQ2X_B_EP",&mQ2X_B_EP);
  mNtuple->SetBranchAddress("mQ2Y_B_EP",&mQ2Y_B_EP);
  mNtuple->SetBranchAddress("mQ3X_B_EP",&mQ3X_B_EP);
  mNtuple->SetBranchAddress("mQ3Y_B_EP",&mQ3Y_B_EP);
  mNtuple->SetBranchAddress("mQ2X_B_SP",&mQ2X_B_SP);
  mNtuple->SetBranchAddress("mQ2Y_B_SP",&mQ2Y_B_SP);
  mNtuple->SetBranchAddress("mQ3X_B_SP",&mQ3X_B_SP);
  mNtuple->SetBranchAddress("mQ3Y_B_SP",&mQ3Y_B_SP);
  mNtuple->SetBranchAddress("mQCounter_B",&mQCounter_B);

  mN_Entires = (Int_t)mNtuple->GetEntries();
}

//----------------------------------------------------------------------------------------------------------

void StTriFlowResolution::InitShiftCorrection(Int_t mEnergy)
{
  TString InPutFile_Shift = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/Correction/Shift/file_%s_Corr_Shift.root",TriFlow::Energy[mEnergy].Data(),TriFlow::Energy[mEnergy].Data());
  mInPutFile_Shift = TFile::Open(InPutFile_Shift.Data());
}

//----------------------------------------------------------------------------------------------------------

void StTriFlowResolution::InitResolution()
{
  for(Int_t j = 0; j < 4; j++)
  {
    TString Res2_EP = Form("Res2_EtaGap_%d_EP",j);
    p_mRes2_EP[j] = new TProfile(Res2_EP.Data(),Res2_EP.Data(),9,-0.5,8.5);
    TString Res3_EP = Form("Res3_EtaGap_%d_EP",j);
    p_mRes3_EP[j] = new TProfile(Res3_EP.Data(),Res3_EP.Data(),9,-0.5,8.5);

    TString Res2_SP = Form("Res2_EtaGap_%d_SP",j);
    p_mRes2_SP[j] = new TProfile(Res2_SP.Data(),Res2_SP.Data(),9,-0.5,8.5);
    TString Res3_SP = Form("Res3_EtaGap_%d_SP",j);
    p_mRes3_SP[j] = new TProfile(Res3_SP.Data(),Res3_SP.Data(),9,-0.5,8.5);
  }

  TString Res2_ran_EP = "Res2_Ran_EP";
  p_mRes2_ran_EP = new TProfile(Res2_ran_EP.Data(),Res2_ran_EP.Data(),9,-0.5,8.5);
  TString Res3_ran_EP = "Res3_Ran_EP";
  p_mRes3_ran_EP = new TProfile(Res3_ran_EP.Data(),Res3_ran_EP.Data(),9,-0.5,8.5);

  TString Res2_ran_SP = "Res2_Ran_SP";
  p_mRes2_ran_SP = new TProfile(Res2_ran_SP.Data(),Res2_ran_SP.Data(),9,-0.5,8.5);
  TString Res3_ran_SP = "Res3_Ran_SP";
  p_mRes3_ran_SP = new TProfile(Res3_ran_SP.Data(),Res3_ran_SP.Data(),9,-0.5,8.5);
}

//----------------------------------------------------------------------------------------------------------
// eta sub start
// Event Plane method
// 2nd
TVector2 StTriFlowResolution::calShiftAngle2East_EP(Int_t vz_sign, Int_t eta_gap)
{
  Float_t Psi_ReCenter = TMath::ATan2(mQ2Y_East_EP[eta_gap],mQ2X_East_EP[eta_gap])/2.0;
  Float_t mean_sin[5], mean_cos[5];
  Float_t delta_Psi = 0.0;
  TVector2 Psi_Shift(0.0,0.0);

  for(Int_t k = 0; k < 5; k++) // Shift Order loop
  {
    TString ProName_cos, ProName_sin;
    TProfile2D *p_cos, *p_sin;

    ProName_cos = Form("CosPsi2_Vertex_%s_EtaGap_%d_Order_%d_East_EP",mVStr[vz_sign].Data(),eta_gap,k);
    p_cos = (TProfile2D*)mInPutFile_Shift->Get(ProName_cos.Data());
    mean_cos[k] = p_cos->GetBinContent(p_cos->FindBin((Double_t)mRunIndex,(Double_t)mCentrality9));

    ProName_sin = Form("SinPsi2_Vertex_%s_EtaGap_%d_Order_%d_East_EP",mVStr[vz_sign].Data(),eta_gap,k);
    p_sin = (TProfile2D*)mInPutFile_Shift->Get(ProName_sin.Data());
    mean_sin[k] = p_sin->GetBinContent(p_sin->FindBin((Double_t)mRunIndex,(Double_t)mCentrality9));

    delta_Psi += (1.0/2.0)*(2.0/(Float_t)(k+1))*(-1.0*mean_sin[k]*TMath::Cos(TriFlow::mShiftOrder2[k]*Psi_ReCenter)+mean_cos[k]*TMath::Sin(TriFlow::mShiftOrder2[k]*Psi_ReCenter));
  }

  Psi_Shift.Set(Psi_ReCenter,delta_Psi);

  return Psi_Shift;
}

TVector2 StTriFlowResolution::calShiftAngle2West_EP(Int_t vz_sign, Int_t eta_gap)
{
  Float_t Psi_ReCenter = TMath::ATan2(mQ2Y_West_EP[eta_gap],mQ2X_West_EP[eta_gap])/2.0;
  Float_t mean_sin[5], mean_cos[5];
  Float_t delta_Psi = 0.0;
  TVector2 Psi_Shift(0.0,0.0);

  for(Int_t k = 0; k < 5; k++) // Shift Order loop
  {
    TString ProName_cos, ProName_sin;
    TProfile2D *p_cos, *p_sin;

    ProName_cos = Form("CosPsi2_Vertex_%s_EtaGap_%d_Order_%d_West_EP",mVStr[vz_sign].Data(),eta_gap,k);
    p_cos = (TProfile2D*)mInPutFile_Shift->Get(ProName_cos.Data());
    mean_cos[k] = p_cos->GetBinContent(p_cos->FindBin((Double_t)mRunIndex,(Double_t)mCentrality9));

    ProName_sin = Form("SinPsi2_Vertex_%s_EtaGap_%d_Order_%d_West_EP",mVStr[vz_sign].Data(),eta_gap,k);
    p_sin = (TProfile2D*)mInPutFile_Shift->Get(ProName_sin.Data());
    mean_sin[k] = p_sin->GetBinContent(p_sin->FindBin((Double_t)mRunIndex,(Double_t)mCentrality9));

    delta_Psi += (1.0/2.0)*(2.0/(Float_t)(k+1))*(-1.0*mean_sin[k]*TMath::Cos(TriFlow::mShiftOrder2[k]*Psi_ReCenter)+mean_cos[k]*TMath::Sin(TriFlow::mShiftOrder2[k]*Psi_ReCenter));
  }

  Psi_Shift.Set(Psi_ReCenter,delta_Psi);

  return Psi_Shift;
}
// 3rd
TVector2 StTriFlowResolution::calShiftAngle3East_EP(Int_t vz_sign, Int_t eta_gap)
{
  Float_t Psi_ReCenter = TMath::ATan2(mQ3Y_East_EP[eta_gap],mQ3X_East_EP[eta_gap])/3.0;
  Float_t mean_sin[5], mean_cos[5];
  Float_t delta_Psi = 0.0;
  TVector2 Psi_Shift(0.0,0.0);

  for(Int_t k = 0; k < 5; k++) // Shift Order loop
  {
    TString ProName_cos, ProName_sin;
    TProfile2D *p_cos, *p_sin;

    ProName_cos = Form("CosPsi3_Vertex_%s_EtaGap_%d_Order_%d_East_EP",mVStr[vz_sign].Data(),eta_gap,k);
    p_cos = (TProfile2D*)mInPutFile_Shift->Get(ProName_cos.Data());
    mean_cos[k] = p_cos->GetBinContent(p_cos->FindBin((Double_t)mRunIndex,(Double_t)mCentrality9));

    ProName_sin = Form("SinPsi3_Vertex_%s_EtaGap_%d_Order_%d_East_EP",mVStr[vz_sign].Data(),eta_gap,k);
    p_sin = (TProfile2D*)mInPutFile_Shift->Get(ProName_sin.Data());
    mean_sin[k] = p_sin->GetBinContent(p_sin->FindBin((Double_t)mRunIndex,(Double_t)mCentrality9));

    delta_Psi += (1.0/3.0)*(2.0/(Float_t)(k+1))*(-1.0*mean_sin[k]*TMath::Cos(TriFlow::mShiftOrder3[k]*Psi_ReCenter)+mean_cos[k]*TMath::Sin(TriFlow::mShiftOrder3[k]*Psi_ReCenter));
  }

  Psi_Shift.Set(Psi_ReCenter,delta_Psi);

  return Psi_Shift;
}

TVector2 StTriFlowResolution::calShiftAngle3West_EP(Int_t vz_sign, Int_t eta_gap)
{
  Float_t Psi_ReCenter = TMath::ATan2(mQ3Y_West_EP[eta_gap],mQ3X_West_EP[eta_gap])/3.0;
  Float_t mean_sin[5], mean_cos[5];
  Float_t delta_Psi = 0.0;
  TVector2 Psi_Shift(0.0,0.0);

  for(Int_t k = 0; k < 5; k++) // Shift Order loop
  {
    TString ProName_cos, ProName_sin;
    TProfile2D *p_cos, *p_sin;

    ProName_cos = Form("CosPsi3_Vertex_%s_EtaGap_%d_Order_%d_West_EP",mVStr[vz_sign].Data(),eta_gap,k);
    p_cos = (TProfile2D*)mInPutFile_Shift->Get(ProName_cos.Data());
    mean_cos[k] = p_cos->GetBinContent(p_cos->FindBin((Double_t)mRunIndex,(Double_t)mCentrality9));

    ProName_sin = Form("SinPsi3_Vertex_%s_EtaGap_%d_Order_%d_West_EP",mVStr[vz_sign].Data(),eta_gap,k);
    p_sin = (TProfile2D*)mInPutFile_Shift->Get(ProName_sin.Data());
    mean_sin[k] = p_sin->GetBinContent(p_sin->FindBin((Double_t)mRunIndex,(Double_t)mCentrality9));

    delta_Psi += (1.0/3.0)*(2.0/(Float_t)(k+1))*(-1.0*mean_sin[k]*TMath::Cos(TriFlow::mShiftOrder3[k]*Psi_ReCenter)+mean_cos[k]*TMath::Sin(TriFlow::mShiftOrder3[k]*Psi_ReCenter));
  }

  Psi_Shift.Set(Psi_ReCenter,delta_Psi);

  return Psi_Shift;
}

// Scalor Product method
// 2nd
TVector2 StTriFlowResolution::calShiftAngle2East_SP(Int_t vz_sign, Int_t eta_gap)
{
  Float_t Psi_ReCenter = TMath::ATan2(mQ2Y_East_SP[eta_gap],mQ2X_East_SP[eta_gap])/2.0;
  Float_t mean_sin[5], mean_cos[5];
  Float_t delta_Psi = 0.0;
  TVector2 Psi_Shift(0.0,0.0);

  for(Int_t k = 0; k < 5; k++) // Shift Order loop
  {
    TString ProName_cos, ProName_sin;
    TProfile2D *p_cos, *p_sin;

    ProName_cos = Form("CosPsi2_Vertex_%s_EtaGap_%d_Order_%d_East_SP",mVStr[vz_sign].Data(),eta_gap,k);
    p_cos = (TProfile2D*)mInPutFile_Shift->Get(ProName_cos.Data());
    mean_cos[k] = p_cos->GetBinContent(p_cos->FindBin((Double_t)mRunIndex,(Double_t)mCentrality9));

    ProName_sin = Form("SinPsi2_Vertex_%s_EtaGap_%d_Order_%d_East_SP",mVStr[vz_sign].Data(),eta_gap,k);
    p_sin = (TProfile2D*)mInPutFile_Shift->Get(ProName_sin.Data());
    mean_sin[k] = p_sin->GetBinContent(p_sin->FindBin((Double_t)mRunIndex,(Double_t)mCentrality9));

    delta_Psi += (1.0/2.0)*(2.0/(Float_t)(k+1))*(-1.0*mean_sin[k]*TMath::Cos(TriFlow::mShiftOrder2[k]*Psi_ReCenter)+mean_cos[k]*TMath::Sin(TriFlow::mShiftOrder2[k]*Psi_ReCenter));
  }

  Psi_Shift.Set(Psi_ReCenter,delta_Psi);

  return Psi_Shift;
}

TVector2 StTriFlowResolution::calShiftAngle2West_SP(Int_t vz_sign, Int_t eta_gap)
{
  Float_t Psi_ReCenter = TMath::ATan2(mQ2Y_West_SP[eta_gap],mQ2X_West_SP[eta_gap])/2.0;
  Float_t mean_sin[5], mean_cos[5];
  Float_t delta_Psi = 0.0;
  TVector2 Psi_Shift(0.0,0.0);

  for(Int_t k = 0; k < 5; k++) // Shift Order loop
  {
    TString ProName_cos, ProName_sin;
    TProfile2D *p_cos, *p_sin;

    ProName_cos = Form("CosPsi2_Vertex_%s_EtaGap_%d_Order_%d_West_SP",mVStr[vz_sign].Data(),eta_gap,k);
    p_cos = (TProfile2D*)mInPutFile_Shift->Get(ProName_cos.Data());
    mean_cos[k] = p_cos->GetBinContent(p_cos->FindBin((Double_t)mRunIndex,(Double_t)mCentrality9));

    ProName_sin = Form("SinPsi2_Vertex_%s_EtaGap_%d_Order_%d_West_SP",mVStr[vz_sign].Data(),eta_gap,k);
    p_sin = (TProfile2D*)mInPutFile_Shift->Get(ProName_sin.Data());
    mean_sin[k] = p_sin->GetBinContent(p_sin->FindBin((Double_t)mRunIndex,(Double_t)mCentrality9));

    delta_Psi += (1.0/2.0)*(2.0/(Float_t)(k+1))*(-1.0*mean_sin[k]*TMath::Cos(TriFlow::mShiftOrder2[k]*Psi_ReCenter)+mean_cos[k]*TMath::Sin(TriFlow::mShiftOrder2[k]*Psi_ReCenter));
  }

  Psi_Shift.Set(Psi_ReCenter,delta_Psi);

  return Psi_Shift;
}
// 3rd
TVector2 StTriFlowResolution::calShiftAngle3East_SP(Int_t vz_sign, Int_t eta_gap)
{
  Float_t Psi_ReCenter = TMath::ATan2(mQ3Y_East_SP[eta_gap],mQ3X_East_SP[eta_gap])/3.0;
  Float_t mean_sin[5], mean_cos[5];
  Float_t delta_Psi = 0.0;
  TVector2 Psi_Shift(0.0,0.0);

  for(Int_t k = 0; k < 5; k++) // Shift Order loop
  {
    TString ProName_cos, ProName_sin;
    TProfile2D *p_cos, *p_sin;

    ProName_cos = Form("CosPsi3_Vertex_%s_EtaGap_%d_Order_%d_East_SP",mVStr[vz_sign].Data(),eta_gap,k);
    p_cos = (TProfile2D*)mInPutFile_Shift->Get(ProName_cos.Data());
    mean_cos[k] = p_cos->GetBinContent(p_cos->FindBin((Double_t)mRunIndex,(Double_t)mCentrality9));

    ProName_sin = Form("SinPsi3_Vertex_%s_EtaGap_%d_Order_%d_East_SP",mVStr[vz_sign].Data(),eta_gap,k);
    p_sin = (TProfile2D*)mInPutFile_Shift->Get(ProName_sin.Data());
    mean_sin[k] = p_sin->GetBinContent(p_sin->FindBin((Double_t)mRunIndex,(Double_t)mCentrality9));

    delta_Psi += (1.0/3.0)*(2.0/(Float_t)(k+1))*(-1.0*mean_sin[k]*TMath::Cos(TriFlow::mShiftOrder3[k]*Psi_ReCenter)+mean_cos[k]*TMath::Sin(TriFlow::mShiftOrder3[k]*Psi_ReCenter));
  }

  Psi_Shift.Set(Psi_ReCenter,delta_Psi);

  return Psi_Shift;
}

TVector2 StTriFlowResolution::calShiftAngle3West_SP(Int_t vz_sign, Int_t eta_gap)
{
  Float_t Psi_ReCenter = TMath::ATan2(mQ3Y_West_SP[eta_gap],mQ3X_West_SP[eta_gap])/3.0;
  Float_t mean_sin[5], mean_cos[5];
  Float_t delta_Psi = 0.0;
  TVector2 Psi_Shift(0.0,0.0);

  for(Int_t k = 0; k < 5; k++) // Shift Order loop
  {
    TString ProName_cos, ProName_sin;
    TProfile2D *p_cos, *p_sin;

    ProName_cos = Form("CosPsi3_Vertex_%s_EtaGap_%d_Order_%d_West_SP",mVStr[vz_sign].Data(),eta_gap,k);
    p_cos = (TProfile2D*)mInPutFile_Shift->Get(ProName_cos.Data());
    mean_cos[k] = p_cos->GetBinContent(p_cos->FindBin((Double_t)mRunIndex,(Double_t)mCentrality9));

    ProName_sin = Form("SinPsi3_Vertex_%s_EtaGap_%d_Order_%d_West_SP",mVStr[vz_sign].Data(),eta_gap,k);
    p_sin = (TProfile2D*)mInPutFile_Shift->Get(ProName_sin.Data());
    mean_sin[k] = p_sin->GetBinContent(p_sin->FindBin((Double_t)mRunIndex,(Double_t)mCentrality9));

    delta_Psi += (1.0/3.0)*(2.0/(Float_t)(k+1))*(-1.0*mean_sin[k]*TMath::Cos(TriFlow::mShiftOrder3[k]*Psi_ReCenter)+mean_cos[k]*TMath::Sin(TriFlow::mShiftOrder3[k]*Psi_ReCenter));
  }

  Psi_Shift.Set(Psi_ReCenter,delta_Psi);

  return Psi_Shift;
}
// eta sub end
//----------------------------------------------------------------------------------------------------------
// random sub start
// Event Plane method
// 2nd
TVector2 StTriFlowResolution::calShiftAngle2A_EP(Int_t vz_sign)
{
  Float_t Psi_ReCenter = TMath::ATan2(mQ2Y_A_EP,mQ2X_A_EP)/2.0;
  Float_t mean_sin[5], mean_cos[5];
  Float_t delta_Psi = 0.0;
  TVector2 Psi_Shift(0.0,0.0);

  for(Int_t k = 0; k < 5; k++) // Shift Order loop
  {
    TString ProName_cos, ProName_sin;
    TProfile2D *p_cos, *p_sin;

    ProName_cos = Form("CosPsi2_Vertex_%s_Order_%d_Full_EP",mVStr[vz_sign].Data(),k);
    p_cos = (TProfile2D*)mInPutFile_Shift->Get(ProName_cos.Data());
    mean_cos[k] = p_cos->GetBinContent(p_cos->FindBin((Double_t)mRunIndex,(Double_t)mCentrality9));

    ProName_sin = Form("SinPsi2_Vertex_%s_Order_%d_Full_EP",mVStr[vz_sign].Data(),k);
    p_sin = (TProfile2D*)mInPutFile_Shift->Get(ProName_sin.Data());
    mean_sin[k] = p_sin->GetBinContent(p_sin->FindBin((Double_t)mRunIndex,(Double_t)mCentrality9));

    delta_Psi += (1.0/2.0)*(2.0/(Float_t)(k+1))*(-1.0*mean_sin[k]*TMath::Cos(TriFlow::mShiftOrder2[k]*Psi_ReCenter)+mean_cos[k]*TMath::Sin(TriFlow::mShiftOrder2[k]*Psi_ReCenter));
  }

  Psi_Shift.Set(Psi_ReCenter,delta_Psi);

  return Psi_Shift;
}

TVector2 StTriFlowResolution::calShiftAngle2B_EP(Int_t vz_sign)
{
  Float_t Psi_ReCenter = TMath::ATan2(mQ2Y_B_EP,mQ2X_B_EP)/2.0;
  Float_t mean_sin[5], mean_cos[5];
  Float_t delta_Psi = 0.0;
  TVector2 Psi_Shift(0.0,0.0);

  for(Int_t k = 0; k < 5; k++) // Shift Order loop
  {
    TString ProName_cos, ProName_sin;
    TProfile2D *p_cos, *p_sin;

    ProName_cos = Form("CosPsi2_Vertex_%s_Order_%d_Full_EP",mVStr[vz_sign].Data(),k);
    p_cos = (TProfile2D*)mInPutFile_Shift->Get(ProName_cos.Data());
    mean_cos[k] = p_cos->GetBinContent(p_cos->FindBin((Double_t)mRunIndex,(Double_t)mCentrality9));

    ProName_sin = Form("SinPsi2_Vertex_%s_Order_%d_Full_EP",mVStr[vz_sign].Data(),k);
    p_sin = (TProfile2D*)mInPutFile_Shift->Get(ProName_sin.Data());
    mean_sin[k] = p_sin->GetBinContent(p_sin->FindBin((Double_t)mRunIndex,(Double_t)mCentrality9));

    delta_Psi += (1.0/2.0)*(2.0/(Float_t)(k+1))*(-1.0*mean_sin[k]*TMath::Cos(TriFlow::mShiftOrder2[k]*Psi_ReCenter)+mean_cos[k]*TMath::Sin(TriFlow::mShiftOrder2[k]*Psi_ReCenter));
  }

  Psi_Shift.Set(Psi_ReCenter,delta_Psi);

  return Psi_Shift;
}

TVector2 StTriFlowResolution::calShiftAngle2Full_EP(Int_t vz_sign)
{
  Float_t Psi_ReCenter = TMath::ATan2(mQ2Y_Full_EP,mQ2X_Full_EP)/2.0;
  Float_t mean_sin[5], mean_cos[5];
  Float_t delta_Psi = 0.0;
  TVector2 Psi_Shift(0.0,0.0);

  for(Int_t k = 0; k < 5; k++) // Shift Order loop
  {
    TString ProName_cos, ProName_sin;
    TProfile2D *p_cos, *p_sin;

    ProName_cos = Form("CosPsi2_Vertex_%s_Order_%d_Full_EP",mVStr[vz_sign].Data(),k);
    p_cos = (TProfile2D*)mInPutFile_Shift->Get(ProName_cos.Data());
    mean_cos[k] = p_cos->GetBinContent(p_cos->FindBin((Double_t)mRunIndex,(Double_t)mCentrality9));

    ProName_sin = Form("SinPsi2_Vertex_%s_Order_%d_Full_EP",mVStr[vz_sign].Data(),k);
    p_sin = (TProfile2D*)mInPutFile_Shift->Get(ProName_sin.Data());
    mean_sin[k] = p_sin->GetBinContent(p_sin->FindBin((Double_t)mRunIndex,(Double_t)mCentrality9));

    delta_Psi += (1.0/2.0)*(2.0/(Float_t)(k+1))*(-1.0*mean_sin[k]*TMath::Cos(TriFlow::mShiftOrder2[k]*Psi_ReCenter)+mean_cos[k]*TMath::Sin(TriFlow::mShiftOrder2[k]*Psi_ReCenter));
  }

  Psi_Shift.Set(Psi_ReCenter,delta_Psi);

  return Psi_Shift;
}
// 3rd
TVector2 StTriFlowResolution::calShiftAngle3A_EP(Int_t vz_sign)
{
  Float_t Psi_ReCenter = TMath::ATan2(mQ3Y_A_EP,mQ3X_A_EP)/3.0;
  Float_t mean_sin[5], mean_cos[5];
  Float_t delta_Psi = 0.0;
  TVector2 Psi_Shift(0.0,0.0);

  for(Int_t k = 0; k < 5; k++) // Shift Order loop
  {
    TString ProName_cos, ProName_sin;
    TProfile2D *p_cos, *p_sin;

    ProName_cos = Form("CosPsi3_Vertex_%s_Order_%d_Full_EP",mVStr[vz_sign].Data(),k);
    p_cos = (TProfile2D*)mInPutFile_Shift->Get(ProName_cos.Data());
    mean_cos[k] = p_cos->GetBinContent(p_cos->FindBin((Double_t)mRunIndex,(Double_t)mCentrality9));

    ProName_sin = Form("SinPsi3_Vertex_%s_Order_%d_Full_EP",mVStr[vz_sign].Data(),k);
    p_sin = (TProfile2D*)mInPutFile_Shift->Get(ProName_sin.Data());
    mean_sin[k] = p_sin->GetBinContent(p_sin->FindBin((Double_t)mRunIndex,(Double_t)mCentrality9));

    delta_Psi += (1.0/3.0)*(2.0/(Float_t)(k+1))*(-1.0*mean_sin[k]*TMath::Cos(TriFlow::mShiftOrder3[k]*Psi_ReCenter)+mean_cos[k]*TMath::Sin(TriFlow::mShiftOrder3[k]*Psi_ReCenter));
  }

  Psi_Shift.Set(Psi_ReCenter,delta_Psi);

  return Psi_Shift;
}

TVector2 StTriFlowResolution::calShiftAngle3B_EP(Int_t vz_sign)
{
  Float_t Psi_ReCenter = TMath::ATan2(mQ3Y_B_EP,mQ3X_B_EP)/3.0;
  Float_t mean_sin[5], mean_cos[5];
  Float_t delta_Psi = 0.0;
  TVector2 Psi_Shift(0.0,0.0);

  for(Int_t k = 0; k < 5; k++) // Shift Order loop
  {
    TString ProName_cos, ProName_sin;
    TProfile2D *p_cos, *p_sin;

    ProName_cos = Form("CosPsi3_Vertex_%s_Order_%d_Full_EP",mVStr[vz_sign].Data(),k);
    p_cos = (TProfile2D*)mInPutFile_Shift->Get(ProName_cos.Data());
    mean_cos[k] = p_cos->GetBinContent(p_cos->FindBin((Double_t)mRunIndex,(Double_t)mCentrality9));

    ProName_sin = Form("SinPsi3_Vertex_%s_Order_%d_Full_EP",mVStr[vz_sign].Data(),k);
    p_sin = (TProfile2D*)mInPutFile_Shift->Get(ProName_sin.Data());
    mean_sin[k] = p_sin->GetBinContent(p_sin->FindBin((Double_t)mRunIndex,(Double_t)mCentrality9));

    delta_Psi += (1.0/3.0)*(2.0/(Float_t)(k+1))*(-1.0*mean_sin[k]*TMath::Cos(TriFlow::mShiftOrder3[k]*Psi_ReCenter)+mean_cos[k]*TMath::Sin(TriFlow::mShiftOrder3[k]*Psi_ReCenter));
  }

  Psi_Shift.Set(Psi_ReCenter,delta_Psi);

  return Psi_Shift;
}

TVector2 StTriFlowResolution::calShiftAngle3Full_EP(Int_t vz_sign)
{
  Float_t Psi_ReCenter = TMath::ATan2(mQ3Y_Full_EP,mQ3X_Full_EP)/3.0;
  Float_t mean_sin[5], mean_cos[5];
  Float_t delta_Psi = 0.0;
  TVector2 Psi_Shift(0.0,0.0);

  for(Int_t k = 0; k < 5; k++) // Shift Order loop
  {
    TString ProName_cos, ProName_sin;
    TProfile2D *p_cos, *p_sin;

    ProName_cos = Form("CosPsi3_Vertex_%s_Order_%d_Full_EP",mVStr[vz_sign].Data(),k);
    p_cos = (TProfile2D*)mInPutFile_Shift->Get(ProName_cos.Data());
    mean_cos[k] = p_cos->GetBinContent(p_cos->FindBin((Double_t)mRunIndex,(Double_t)mCentrality9));

    ProName_sin = Form("SinPsi3_Vertex_%s_Order_%d_Full_EP",mVStr[vz_sign].Data(),k);
    p_sin = (TProfile2D*)mInPutFile_Shift->Get(ProName_sin.Data());
    mean_sin[k] = p_sin->GetBinContent(p_sin->FindBin((Double_t)mRunIndex,(Double_t)mCentrality9));

    delta_Psi += (1.0/3.0)*(2.0/(Float_t)(k+1))*(-1.0*mean_sin[k]*TMath::Cos(TriFlow::mShiftOrder3[k]*Psi_ReCenter)+mean_cos[k]*TMath::Sin(TriFlow::mShiftOrder3[k]*Psi_ReCenter));
  }

  Psi_Shift.Set(Psi_ReCenter,delta_Psi);

  return Psi_Shift;
}

// Scalar Product method
// 2nd
TVector2 StTriFlowResolution::calShiftAngle2A_SP(Int_t vz_sign)
{
  Float_t Psi_ReCenter = TMath::ATan2(mQ2Y_A_SP,mQ2X_A_SP)/2.0;
  Float_t mean_sin[5], mean_cos[5];
  Float_t delta_Psi = 0.0;
  TVector2 Psi_Shift(0.0,0.0);

  for(Int_t k = 0; k < 5; k++) // Shift Order loop
  {
    TString ProName_cos, ProName_sin;
    TProfile2D *p_cos, *p_sin;

    ProName_cos = Form("CosPsi2_Vertex_%s_Order_%d_Full_SP",mVStr[vz_sign].Data(),k);
    p_cos = (TProfile2D*)mInPutFile_Shift->Get(ProName_cos.Data());
    mean_cos[k] = p_cos->GetBinContent(p_cos->FindBin((Double_t)mRunIndex,(Double_t)mCentrality9));

    ProName_sin = Form("SinPsi2_Vertex_%s_Order_%d_Full_SP",mVStr[vz_sign].Data(),k);
    p_sin = (TProfile2D*)mInPutFile_Shift->Get(ProName_sin.Data());
    mean_sin[k] = p_sin->GetBinContent(p_sin->FindBin((Double_t)mRunIndex,(Double_t)mCentrality9));

    delta_Psi += (1.0/2.0)*(2.0/(Float_t)(k+1))*(-1.0*mean_sin[k]*TMath::Cos(TriFlow::mShiftOrder2[k]*Psi_ReCenter)+mean_cos[k]*TMath::Sin(TriFlow::mShiftOrder2[k]*Psi_ReCenter));
  }

  Psi_Shift.Set(Psi_ReCenter,delta_Psi);

  return Psi_Shift;
}

TVector2 StTriFlowResolution::calShiftAngle2B_SP(Int_t vz_sign)
{
  Float_t Psi_ReCenter = TMath::ATan2(mQ2Y_B_SP,mQ2X_B_SP)/2.0;
  Float_t mean_sin[5], mean_cos[5];
  Float_t delta_Psi = 0.0;
  TVector2 Psi_Shift(0.0,0.0);

  for(Int_t k = 0; k < 5; k++) // Shift Order loop
  {
    TString ProName_cos, ProName_sin;
    TProfile2D *p_cos, *p_sin;

    ProName_cos = Form("CosPsi2_Vertex_%s_Order_%d_Full_SP",mVStr[vz_sign].Data(),k);
    p_cos = (TProfile2D*)mInPutFile_Shift->Get(ProName_cos.Data());
    mean_cos[k] = p_cos->GetBinContent(p_cos->FindBin((Double_t)mRunIndex,(Double_t)mCentrality9));

    ProName_sin = Form("SinPsi2_Vertex_%s_Order_%d_Full_SP",mVStr[vz_sign].Data(),k);
    p_sin = (TProfile2D*)mInPutFile_Shift->Get(ProName_sin.Data());
    mean_sin[k] = p_sin->GetBinContent(p_sin->FindBin((Double_t)mRunIndex,(Double_t)mCentrality9));

    delta_Psi += (1.0/2.0)*(2.0/(Float_t)(k+1))*(-1.0*mean_sin[k]*TMath::Cos(TriFlow::mShiftOrder2[k]*Psi_ReCenter)+mean_cos[k]*TMath::Sin(TriFlow::mShiftOrder2[k]*Psi_ReCenter));
  }

  Psi_Shift.Set(Psi_ReCenter,delta_Psi);

  return Psi_Shift;
}

TVector2 StTriFlowResolution::calShiftAngle2Full_SP(Int_t vz_sign)
{
  Float_t Psi_ReCenter = TMath::ATan2(mQ2Y_Full_SP,mQ2X_Full_SP)/2.0;
  Float_t mean_sin[5], mean_cos[5];
  Float_t delta_Psi = 0.0;
  TVector2 Psi_Shift(0.0,0.0);

  for(Int_t k = 0; k < 5; k++) // Shift Order loop
  {
    TString ProName_cos, ProName_sin;
    TProfile2D *p_cos, *p_sin;

    ProName_cos = Form("CosPsi2_Vertex_%s_Order_%d_Full_SP",mVStr[vz_sign].Data(),k);
    p_cos = (TProfile2D*)mInPutFile_Shift->Get(ProName_cos.Data());
    mean_cos[k] = p_cos->GetBinContent(p_cos->FindBin((Double_t)mRunIndex,(Double_t)mCentrality9));

    ProName_sin = Form("SinPsi2_Vertex_%s_Order_%d_Full_SP",mVStr[vz_sign].Data(),k);
    p_sin = (TProfile2D*)mInPutFile_Shift->Get(ProName_sin.Data());
    mean_sin[k] = p_sin->GetBinContent(p_sin->FindBin((Double_t)mRunIndex,(Double_t)mCentrality9));

    delta_Psi += (1.0/2.0)*(2.0/(Float_t)(k+1))*(-1.0*mean_sin[k]*TMath::Cos(TriFlow::mShiftOrder2[k]*Psi_ReCenter)+mean_cos[k]*TMath::Sin(TriFlow::mShiftOrder2[k]*Psi_ReCenter));
  }

  Psi_Shift.Set(Psi_ReCenter,delta_Psi);

  return Psi_Shift;
}
// 3rd
TVector2 StTriFlowResolution::calShiftAngle3A_SP(Int_t vz_sign)
{
  Float_t Psi_ReCenter = TMath::ATan2(mQ3Y_A_SP,mQ3X_A_SP)/3.0;
  Float_t mean_sin[5], mean_cos[5];
  Float_t delta_Psi = 0.0;
  TVector2 Psi_Shift(0.0,0.0);

  for(Int_t k = 0; k < 5; k++) // Shift Order loop
  {
    TString ProName_cos, ProName_sin;
    TProfile2D *p_cos, *p_sin;

    ProName_cos = Form("CosPsi3_Vertex_%s_Order_%d_Full_SP",mVStr[vz_sign].Data(),k);
    p_cos = (TProfile2D*)mInPutFile_Shift->Get(ProName_cos.Data());
    mean_cos[k] = p_cos->GetBinContent(p_cos->FindBin((Double_t)mRunIndex,(Double_t)mCentrality9));

    ProName_sin = Form("SinPsi3_Vertex_%s_Order_%d_Full_SP",mVStr[vz_sign].Data(),k);
    p_sin = (TProfile2D*)mInPutFile_Shift->Get(ProName_sin.Data());
    mean_sin[k] = p_sin->GetBinContent(p_sin->FindBin((Double_t)mRunIndex,(Double_t)mCentrality9));

    delta_Psi += (1.0/3.0)*(2.0/(Float_t)(k+1))*(-1.0*mean_sin[k]*TMath::Cos(TriFlow::mShiftOrder3[k]*Psi_ReCenter)+mean_cos[k]*TMath::Sin(TriFlow::mShiftOrder3[k]*Psi_ReCenter));
  }

  Psi_Shift.Set(Psi_ReCenter,delta_Psi);

  return Psi_Shift;
}

TVector2 StTriFlowResolution::calShiftAngle3B_SP(Int_t vz_sign)
{
  Float_t Psi_ReCenter = TMath::ATan2(mQ3Y_B_SP,mQ3X_B_SP)/3.0;
  Float_t mean_sin[5], mean_cos[5];
  Float_t delta_Psi = 0.0;
  TVector2 Psi_Shift(0.0,0.0);

  for(Int_t k = 0; k < 5; k++) // Shift Order loop
  {
    TString ProName_cos, ProName_sin;
    TProfile2D *p_cos, *p_sin;

    ProName_cos = Form("CosPsi3_Vertex_%s_Order_%d_Full_SP",mVStr[vz_sign].Data(),k);
    p_cos = (TProfile2D*)mInPutFile_Shift->Get(ProName_cos.Data());
    mean_cos[k] = p_cos->GetBinContent(p_cos->FindBin((Double_t)mRunIndex,(Double_t)mCentrality9));

    ProName_sin = Form("SinPsi3_Vertex_%s_Order_%d_Full_SP",mVStr[vz_sign].Data(),k);
    p_sin = (TProfile2D*)mInPutFile_Shift->Get(ProName_sin.Data());
    mean_sin[k] = p_sin->GetBinContent(p_sin->FindBin((Double_t)mRunIndex,(Double_t)mCentrality9));

    delta_Psi += (1.0/3.0)*(2.0/(Float_t)(k+1))*(-1.0*mean_sin[k]*TMath::Cos(TriFlow::mShiftOrder3[k]*Psi_ReCenter)+mean_cos[k]*TMath::Sin(TriFlow::mShiftOrder3[k]*Psi_ReCenter));
  }

  Psi_Shift.Set(Psi_ReCenter,delta_Psi);

  return Psi_Shift;
}

TVector2 StTriFlowResolution::calShiftAngle3Full_SP(Int_t vz_sign)
{
  Float_t Psi_ReCenter = TMath::ATan2(mQ3Y_Full_SP,mQ3X_Full_SP)/3.0;
  Float_t mean_sin[5], mean_cos[5];
  Float_t delta_Psi = 0.0;
  TVector2 Psi_Shift(0.0,0.0);

  for(Int_t k = 0; k < 5; k++) // Shift Order loop
  {
    TString ProName_cos, ProName_sin;
    TProfile2D *p_cos, *p_sin;

    ProName_cos = Form("CosPsi3_Vertex_%s_Order_%d_Full_SP",mVStr[vz_sign].Data(),k);
    p_cos = (TProfile2D*)mInPutFile_Shift->Get(ProName_cos.Data());
    mean_cos[k] = p_cos->GetBinContent(p_cos->FindBin((Double_t)mRunIndex,(Double_t)mCentrality9));

    ProName_sin = Form("SinPsi3_Vertex_%s_Order_%d_Full_SP",mVStr[vz_sign].Data(),k);
    p_sin = (TProfile2D*)mInPutFile_Shift->Get(ProName_sin.Data());
    mean_sin[k] = p_sin->GetBinContent(p_sin->FindBin((Double_t)mRunIndex,(Double_t)mCentrality9));

    delta_Psi += (1.0/3.0)*(2.0/(Float_t)(k+1))*(-1.0*mean_sin[k]*TMath::Cos(TriFlow::mShiftOrder3[k]*Psi_ReCenter)+mean_cos[k]*TMath::Sin(TriFlow::mShiftOrder3[k]*Psi_ReCenter));
  }

  Psi_Shift.Set(Psi_ReCenter,delta_Psi);

  return Psi_Shift;
}
// random sub stop
//----------------------------------------------------------------------------------------------------------

Float_t StTriFlowResolution::Psi_Shift(Float_t Psi, Float_t order)
{
  Float_t Psi_Corr = Psi;
  if(Psi > TMath::Pi()/order)
  {
    Psi_Corr = Psi - 2.0*TMath::Pi()/order;
  }
  if(Psi < -1.0*TMath::Pi()/order)
  {
    Psi_Corr = Psi + 2.0*TMath::Pi()/order;
  }

  return Psi_Corr;
}

//----------------------------------------------------------------------------------------------------------

void StTriFlowResolution::getResolution()
{
  for(Int_t n_entires = 0; n_entires < mN_Entires; n_entires++)
  {
    Double_t res2_EP[4], res3_EP[4];
    Double_t res2_SP[4], res3_SP[4];
    Double_t res2_ran_EP, res3_ran_EP;
    Double_t res2_ran_SP, res3_ran_SP;
    mNtuple->GetEntry(n_entires);

    // vz sign
    Int_t vz_sign;
    if(mVz > 0.0)
    {
      vz_sign = 0;
    }
    else
    {
      vz_sign = 1;
    }

    // eta sub
    for(Int_t j = 0; j < 4; j++)
    {
      if(mQCounter_East[j] > TriFlow::mTrackMin && mQCounter_West[j] > TriFlow::mTrackMin)
      {
	// Event Plane method
	// 2nd resolution
        TVector2 Psi2_East_EP = calShiftAngle2East_EP(vz_sign,j);
	Float_t  Psi2_Shift_East_EP_raw = Psi2_East_EP.X() + Psi2_East_EP.Y();
	Float_t  Psi2_Shift_East_EP = Psi_Shift(Psi2_Shift_East_EP_raw,2.0);

        TVector2 Psi2_West_EP = calShiftAngle2West_EP(vz_sign,j);
	Float_t  Psi2_Shift_West_EP_raw = Psi2_West_EP.X() + Psi2_West_EP.Y();
	Float_t  Psi2_Shift_West_EP = Psi_Shift(Psi2_Shift_West_EP_raw,2.0);
	res2_EP[j] = TMath::Cos(2.0*((Double_t)Psi2_Shift_East_EP-(Double_t)Psi2_Shift_West_EP));

	// 3rd resolution
        TVector2 Psi3_East_EP = calShiftAngle3East_EP(vz_sign,j);
	Float_t  Psi3_Shift_East_EP_raw = Psi3_East_EP.X() + Psi3_East_EP.Y();
	Float_t  Psi3_Shift_East_EP = Psi_Shift(Psi3_Shift_East_EP_raw,3.0);

        TVector2 Psi3_West_EP = calShiftAngle3West_EP(vz_sign,j);
	Float_t  Psi3_Shift_West_EP_raw = Psi3_West_EP.X() + Psi3_West_EP.Y();
	Float_t  Psi3_Shift_West_EP = Psi_Shift(Psi3_Shift_West_EP_raw,3.0);
	res3_EP[j] = TMath::Cos(3.0*((Double_t)Psi3_Shift_East_EP-(Double_t)Psi3_Shift_West_EP));

	p_mRes2_EP[j]->Fill((Double_t)mCentrality9,res2_EP[j]);
	p_mRes3_EP[j]->Fill((Double_t)mCentrality9,res3_EP[j]);
	mTriFlowHistManger->FillEventPlane_East_EP((Float_t)Psi2_East_EP.X(),(Float_t)Psi2_Shift_East_EP,(Float_t)Psi3_East_EP.X(),(Float_t)Psi3_Shift_East_EP,(Int_t)mCentrality9,j);
	mTriFlowHistManger->FillEventPlane_West_EP((Float_t)Psi2_West_EP.X(),(Float_t)Psi2_Shift_West_EP,(Float_t)Psi3_West_EP.X(),(Float_t)Psi3_Shift_West_EP,(Int_t)mCentrality9,j);

	// Scalor Product method
	// 2nd resolution
        TVector2 Psi2_East_SP = calShiftAngle2East_SP(vz_sign,j);
	Float_t  Psi2_Shift_East_SP_raw = Psi2_East_SP.X() + Psi2_East_SP.Y();
	Float_t  Psi2_Shift_East_SP = Psi_Shift(Psi2_Shift_East_SP_raw,2.0);
	Float_t  Q2_East = TMath::Sqrt(mQ2X_East_SP[j]*mQ2X_East_SP[j] + mQ2Y_East_SP[j]*mQ2Y_East_SP[j]);
	Float_t  Q2X_East = Q2_East*TMath::Cos(2.0*Psi2_Shift_East_SP);
	Float_t  Q2Y_East = Q2_East*TMath::Sin(2.0*Psi2_Shift_East_SP);

        TVector2 Psi2_West_SP = calShiftAngle2West_SP(vz_sign,j);
	Float_t  Psi2_Shift_West_SP_raw = Psi2_West_SP.X() + Psi2_West_SP.Y();
	Float_t  Psi2_Shift_West_SP = Psi_Shift(Psi2_Shift_West_SP_raw,2.0);
	Float_t  Q2_West = TMath::Sqrt(mQ2X_West_SP[j]*mQ2X_West_SP[j] + mQ2Y_West_SP[j]*mQ2Y_West_SP[j]);
	Float_t  Q2X_West = Q2_West*TMath::Cos(2.0*Psi2_Shift_West_SP);
	Float_t  Q2Y_West = Q2_West*TMath::Sin(2.0*Psi2_Shift_West_SP);
	res2_SP[j] = Q2X_East*Q2X_West+Q2Y_East*Q2Y_West;

	// 3rd resolution
        TVector2 Psi3_East_SP = calShiftAngle3East_SP(vz_sign,j);
	Float_t  Psi3_Shift_East_SP_raw = Psi3_East_SP.X() + Psi3_East_SP.Y();
	Float_t  Psi3_Shift_East_SP = Psi_Shift(Psi3_Shift_East_SP_raw,3.0);
	Float_t  Q3_East = TMath::Sqrt(mQ3X_East_SP[j]*mQ3X_East_SP[j] + mQ3Y_East_SP[j]*mQ3Y_East_SP[j]);
	Float_t  Q3X_East = Q3_East*TMath::Cos(3.0*Psi3_Shift_East_SP);
	Float_t  Q3Y_East = Q3_East*TMath::Sin(3.0*Psi3_Shift_East_SP);

        TVector2 Psi3_West_SP = calShiftAngle3West_SP(vz_sign,j);
	Float_t  Psi3_Shift_West_SP_raw = Psi3_West_SP.X() + Psi3_West_SP.Y();
	Float_t  Psi3_Shift_West_SP = Psi_Shift(Psi3_Shift_West_SP_raw,3.0);
	Float_t  Q3_West = TMath::Sqrt(mQ3X_West_SP[j]*mQ3X_West_SP[j] + mQ3Y_West_SP[j]*mQ3Y_West_SP[j]);
	Float_t  Q3X_West = Q3_West*TMath::Cos(3.0*Psi3_Shift_West_SP);
	Float_t  Q3Y_West = Q3_West*TMath::Sin(3.0*Psi3_Shift_West_SP);
	res3_SP[j] = Q3X_East*Q3X_West+Q3Y_East*Q3Y_West;

	p_mRes2_SP[j]->Fill((Double_t)mCentrality9,res2_SP[j]);
	p_mRes3_SP[j]->Fill((Double_t)mCentrality9,res3_SP[j]);
	mTriFlowHistManger->FillEventPlane_East_SP((Float_t)Psi2_East_SP.X(),(Float_t)Psi2_Shift_East_SP,(Float_t)Psi3_East_SP.X(),(Float_t)Psi3_Shift_East_SP,(Int_t)mCentrality9,j);
	mTriFlowHistManger->FillEventPlane_West_SP((Float_t)Psi2_West_SP.X(),(Float_t)Psi2_Shift_West_SP,(Float_t)Psi3_West_SP.X(),(Float_t)Psi3_Shift_West_SP,(Int_t)mCentrality9,j);
      }
    }
    // random sub
    if(mQCounter_Full > TriFlow::mTrackMin_Full && mQCounter_Full_East > 0 && mQCounter_Full_West > 0)
    {
      // Event Plane method
      // 2nd resolution
      TVector2 Psi2_A_EP = calShiftAngle2A_EP(vz_sign);
      Float_t  Psi2_Shift_A_EP_raw = Psi2_A_EP.X() + Psi2_A_EP.Y();
      Float_t  Psi2_Shift_A_EP = Psi_Shift(Psi2_Shift_A_EP_raw,2.0);

      TVector2 Psi2_B_EP = calShiftAngle2B_EP(vz_sign);
      Float_t  Psi2_Shift_B_EP_raw = Psi2_B_EP.X() + Psi2_B_EP.Y();
      Float_t  Psi2_Shift_B_EP = Psi_Shift(Psi2_Shift_B_EP_raw,2.0);
      res2_ran_EP = TMath::Cos(2.0*((Double_t)Psi2_Shift_A_EP-(Double_t)Psi2_Shift_B_EP));

      TVector2 Psi2_Full_EP = calShiftAngle2Full_EP(vz_sign);
      Float_t  Psi2_Shift_Full_EP_raw = Psi2_Full_EP.X() + Psi2_Full_EP.Y();
      Float_t  Psi2_Shift_Full_EP = Psi_Shift(Psi2_Shift_Full_EP_raw,2.0);

      // 3rd resolution
      TVector2 Psi3_A_EP = calShiftAngle3A_EP(vz_sign);
      Float_t  Psi3_Shift_A_EP_raw = Psi3_A_EP.X() + Psi3_A_EP.Y();
      Float_t  Psi3_Shift_A_EP = Psi_Shift(Psi3_Shift_A_EP_raw,3.0);

      TVector2 Psi3_B_EP = calShiftAngle3B_EP(vz_sign);
      Float_t  Psi3_Shift_B_EP_raw = Psi3_B_EP.X() + Psi3_B_EP.Y();
      Float_t  Psi3_Shift_B_EP = Psi_Shift(Psi3_Shift_B_EP_raw,3.0);
      res3_ran_EP = TMath::Cos(3.0*((Double_t)Psi3_Shift_A_EP-(Double_t)Psi3_Shift_B_EP));

      TVector2 Psi3_Full_EP = calShiftAngle3Full_EP(vz_sign);
      Float_t  Psi3_Shift_Full_EP_raw = Psi3_Full_EP.X() + Psi3_Full_EP.Y();
      Float_t  Psi3_Shift_Full_EP = Psi_Shift(Psi3_Shift_Full_EP_raw,3.0);

      p_mRes2_ran_EP->Fill((Double_t)mCentrality9,res2_ran_EP);
      p_mRes3_ran_EP->Fill((Double_t)mCentrality9,res3_ran_EP);
      mTriFlowHistManger->FillEventPlane_A_EP((Float_t)Psi2_A_EP.X(),(Float_t)Psi2_Shift_A_EP,(Float_t)Psi3_A_EP.X(),(Float_t)Psi3_Shift_A_EP,(Int_t)mCentrality9);
      mTriFlowHistManger->FillEventPlane_B_EP((Float_t)Psi2_B_EP.X(),(Float_t)Psi2_Shift_B_EP,(Float_t)Psi3_B_EP.X(),(Float_t)Psi3_Shift_B_EP,(Int_t)mCentrality9);
      mTriFlowHistManger->FillEventPlane_Full_EP((Float_t)Psi2_Full_EP.X(),(Float_t)Psi2_Shift_Full_EP,(Float_t)Psi3_Full_EP.X(),(Float_t)Psi3_Shift_Full_EP,(Int_t)mCentrality9);

      // Scalor Product method
      // 2nd resolution
      TVector2 Psi2_A_SP = calShiftAngle2A_SP(vz_sign);
      Float_t  Psi2_Shift_A_SP_raw = Psi2_A_SP.X() + Psi2_A_SP.Y();
      Float_t  Psi2_Shift_A_SP = Psi_Shift(Psi2_Shift_A_SP_raw,2.0);
      Float_t  Q2_A = TMath::Sqrt(mQ2X_A_SP*mQ2X_A_SP + mQ2Y_A_SP*mQ2Y_A_SP);
      Float_t  Q2X_A = Q2_A*TMath::Cos(2.0*Psi2_Shift_A_SP);
      Float_t  Q2Y_A = Q2_A*TMath::Sin(2.0*Psi2_Shift_A_SP);

      TVector2 Psi2_B_SP = calShiftAngle2B_SP(vz_sign);
      Float_t  Psi2_Shift_B_SP_raw = Psi2_B_SP.X() + Psi2_B_SP.Y();
      Float_t  Psi2_Shift_B_SP = Psi_Shift(Psi2_Shift_B_SP_raw,2.0);
      Float_t  Q2_B = TMath::Sqrt(mQ2X_B_SP*mQ2X_B_SP + mQ2Y_B_SP*mQ2Y_B_SP);
      Float_t  Q2X_B = Q2_B*TMath::Cos(2.0*Psi2_Shift_B_SP);
      Float_t  Q2Y_B = Q2_B*TMath::Sin(2.0*Psi2_Shift_B_SP);
      res2_ran_SP = Q2X_A*Q2X_B+Q2Y_A*Q2Y_B;

      TVector2 Psi2_Full_SP = calShiftAngle2Full_SP(vz_sign);
      Float_t  Psi2_Shift_Full_SP_raw = Psi2_Full_SP.X() + Psi2_Full_SP.Y();
      Float_t  Psi2_Shift_Full_SP = Psi_Shift(Psi2_Shift_Full_SP_raw,2.0);

      // 3rd resolution
      TVector2 Psi3_A_SP = calShiftAngle3A_SP(vz_sign);
      Float_t  Psi3_Shift_A_SP_raw = Psi3_A_SP.X() + Psi3_A_SP.Y();
      Float_t  Psi3_Shift_A_SP = Psi_Shift(Psi3_Shift_A_SP_raw,3.0);
      Float_t  Q3_A = TMath::Sqrt(mQ3X_A_SP*mQ3X_A_SP + mQ3Y_A_SP*mQ3Y_A_SP);
      Float_t  Q3X_A = Q3_A*TMath::Cos(3.0*Psi3_Shift_A_SP);
      Float_t  Q3Y_A = Q3_A*TMath::Sin(3.0*Psi3_Shift_A_SP);

      TVector2 Psi3_B_SP = calShiftAngle3B_SP(vz_sign);
      Float_t  Psi3_Shift_B_SP_raw = Psi3_B_SP.X() + Psi3_B_SP.Y();
      Float_t  Psi3_Shift_B_SP = Psi_Shift(Psi3_Shift_B_SP_raw,3.0);
      Float_t  Q3_B = TMath::Sqrt(mQ3X_B_SP*mQ3X_B_SP + mQ3Y_B_SP*mQ3Y_B_SP);
      Float_t  Q3X_B = Q3_B*TMath::Cos(3.0*Psi3_Shift_B_SP);
      Float_t  Q3Y_B = Q3_B*TMath::Sin(3.0*Psi3_Shift_B_SP);
      res3_ran_SP = Q3X_A*Q3X_B+Q3Y_A*Q3Y_B;

      TVector2 Psi3_Full_SP = calShiftAngle3Full_SP(vz_sign);
      Float_t  Psi3_Shift_Full_SP_raw = Psi3_Full_SP.X() + Psi3_Full_SP.Y();
      Float_t  Psi3_Shift_Full_SP = Psi_Shift(Psi3_Shift_Full_SP_raw,3.0);

      p_mRes2_ran_SP->Fill((Double_t)mCentrality9,res2_ran_SP);
      p_mRes3_ran_SP->Fill((Double_t)mCentrality9,res3_ran_SP);
      mTriFlowHistManger->FillEventPlane_A_SP((Float_t)Psi2_A_SP.X(),(Float_t)Psi2_Shift_A_SP,(Float_t)Psi3_A_SP.X(),(Float_t)Psi3_Shift_A_SP,(Int_t)mCentrality9);
      mTriFlowHistManger->FillEventPlane_B_SP((Float_t)Psi2_B_SP.X(),(Float_t)Psi2_Shift_B_SP,(Float_t)Psi3_B_SP.X(),(Float_t)Psi3_Shift_B_SP,(Int_t)mCentrality9);
      mTriFlowHistManger->FillEventPlane_Full_SP((Float_t)Psi2_Full_SP.X(),(Float_t)Psi2_Shift_Full_SP,(Float_t)Psi3_Full_SP.X(),(Float_t)Psi3_Shift_Full_SP,(Int_t)mCentrality9);
    }
  }
}

//----------------------------------------------------------------------------------------------------------

void StTriFlowResolution::WriteResolution()
{
  for(Int_t j = 0; j < 4; j++)
  {
    p_mRes2_EP[j]->Write();
    p_mRes3_EP[j]->Write();
    p_mRes2_SP[j]->Write();
    p_mRes3_SP[j]->Write();
  }
  p_mRes2_ran_EP->Write();
  p_mRes3_ran_EP->Write();
  p_mRes2_ran_SP->Write();
  p_mRes3_ran_SP->Write();
}

//----------------------------------------------------------------------------------------------------------

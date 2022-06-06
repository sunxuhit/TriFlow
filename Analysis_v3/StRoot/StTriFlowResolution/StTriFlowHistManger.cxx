#include "StTriFlowHistManger.h"
#include "TMath.h"
#include "TString.h"
#include "TH1F.h"

ClassImp(StTriFlowHistManger)

//--------------------------------------------------------------------------------------

StTriFlowHistManger::StTriFlowHistManger()
{
}

StTriFlowHistManger::~StTriFlowHistManger()
{
}

//--------------------------------------------------------------------------------------

void StTriFlowHistManger::InitEventPlane()
{
  for(Int_t i = 0; i < 9; i++) // Centrality
  {
    for(Int_t j = 0; j < 4; j++) // eta_gap
    {
      TString HistName;
      // Event Plane method
      // ReCenter Correction
      HistName = Form("Psi2_Centrality_%d_EtaGap_%d_ReCenter_East_EP",i,j);
      h_mPsi2_East_ReCenter_EP[i][j] = new TH1F(HistName.Data(),HistName.Data(),360,-1.0*TMath::Pi()/2.0,TMath::Pi()/2.0);
      HistName = Form("Psi2_Centrality_%d_EtaGap_%d_ReCenter_West_EP",i,j);
      h_mPsi2_West_ReCenter_EP[i][j] = new TH1F(HistName.Data(),HistName.Data(),360,-1.0*TMath::Pi()/2.0,TMath::Pi()/2.0);
      HistName = Form("Psi3_Centrality_%d_EtaGap_%d_ReCenter_East_EP",i,j);
      h_mPsi3_East_ReCenter_EP[i][j] = new TH1F(HistName.Data(),HistName.Data(),360,-1.0*TMath::Pi()/3.0,TMath::Pi()/3.0);
      HistName = Form("Psi3_Centrality_%d_EtaGap_%d_ReCenter_West_EP",i,j);
      h_mPsi3_West_ReCenter_EP[i][j] = new TH1F(HistName.Data(),HistName.Data(),360,-1.0*TMath::Pi()/3.0,TMath::Pi()/3.0);

      // Shift Correction
      HistName = Form("Psi2_Centrality_%d_EtaGap_%d_Shift_East_EP",i,j);
      h_mPsi2_East_Shift_EP[i][j] = new TH1F(HistName.Data(),HistName.Data(),360,-1.0*TMath::Pi()/2.0,TMath::Pi()/2.0);
      HistName = Form("Psi2_Centrality_%d_EtaGap_%d_Shift_West_EP",i,j);
      h_mPsi2_West_Shift_EP[i][j] = new TH1F(HistName.Data(),HistName.Data(),360,-1.0*TMath::Pi()/2.0,TMath::Pi()/2.0);
      HistName = Form("Psi3_Centrality_%d_EtaGap_%d_Shift_East_EP",i,j);
      h_mPsi3_East_Shift_EP[i][j] = new TH1F(HistName.Data(),HistName.Data(),360,-1.0*TMath::Pi()/3.0,TMath::Pi()/3.0);
      HistName = Form("Psi3_Centrality_%d_EtaGap_%d_Shift_West_EP",i,j);
      h_mPsi3_West_Shift_EP[i][j] = new TH1F(HistName.Data(),HistName.Data(),360,-1.0*TMath::Pi()/3.0,TMath::Pi()/3.0);

      // Scalor Product method
      // ReCenter Correction
      HistName = Form("Psi2_Centrality_%d_EtaGap_%d_ReCenter_East_SP",i,j);
      h_mPsi2_East_ReCenter_SP[i][j] = new TH1F(HistName.Data(),HistName.Data(),360,-1.0*TMath::Pi()/2.0,TMath::Pi()/2.0);
      HistName = Form("Psi2_Centrality_%d_EtaGap_%d_ReCenter_West_SP",i,j);
      h_mPsi2_West_ReCenter_SP[i][j] = new TH1F(HistName.Data(),HistName.Data(),360,-1.0*TMath::Pi()/2.0,TMath::Pi()/2.0);
      HistName = Form("Psi3_Centrality_%d_EtaGap_%d_ReCenter_East_SP",i,j);
      h_mPsi3_East_ReCenter_SP[i][j] = new TH1F(HistName.Data(),HistName.Data(),360,-1.0*TMath::Pi()/3.0,TMath::Pi()/3.0);
      HistName = Form("Psi3_Centrality_%d_EtaGap_%d_ReCenter_West_SP",i,j);
      h_mPsi3_West_ReCenter_SP[i][j] = new TH1F(HistName.Data(),HistName.Data(),360,-1.0*TMath::Pi()/3.0,TMath::Pi()/3.0);

      // Shift Correction
      HistName = Form("Psi2_Centrality_%d_EtaGap_%d_Shift_East_SP",i,j);
      h_mPsi2_East_Shift_SP[i][j] = new TH1F(HistName.Data(),HistName.Data(),360,-1.0*TMath::Pi()/2.0,TMath::Pi()/2.0);
      HistName = Form("Psi2_Centrality_%d_EtaGap_%d_Shift_West_SP",i,j);
      h_mPsi2_West_Shift_SP[i][j] = new TH1F(HistName.Data(),HistName.Data(),360,-1.0*TMath::Pi()/2.0,TMath::Pi()/2.0);
      HistName = Form("Psi3_Centrality_%d_EtaGap_%d_Shift_East_SP",i,j);
      h_mPsi3_East_Shift_SP[i][j] = new TH1F(HistName.Data(),HistName.Data(),360,-1.0*TMath::Pi()/3.0,TMath::Pi()/3.0);
      HistName = Form("Psi3_Centrality_%d_EtaGap_%d_Shift_West_SP",i,j);
      h_mPsi3_West_Shift_SP[i][j] = new TH1F(HistName.Data(),HistName.Data(),360,-1.0*TMath::Pi()/3.0,TMath::Pi()/3.0);
    }
    TString HistName_Full;
    // Event Plane method
    // ReCenter Correction
    HistName_Full = Form("Psi2_Centrality_%d_ReCenter_A_EP",i);
    h_mPsi2_A_ReCenter_EP[i] = new TH1F(HistName_Full.Data(),HistName_Full.Data(),360,-1.0*TMath::Pi()/2.0,TMath::Pi()/2.0);
    HistName_Full = Form("Psi2_Centrality_%d_ReCenter_B_EP",i);
    h_mPsi2_B_ReCenter_EP[i] = new TH1F(HistName_Full.Data(),HistName_Full.Data(),360,-1.0*TMath::Pi()/2.0,TMath::Pi()/2.0);
    HistName_Full = Form("Psi2_Centrality_%d_ReCenter_Full_EP",i);
    h_mPsi2_Full_ReCenter_EP[i] = new TH1F(HistName_Full.Data(),HistName_Full.Data(),360,-1.0*TMath::Pi()/2.0,TMath::Pi()/2.0);
    HistName_Full = Form("Psi3_Centrality_%d_ReCenter_A_EP",i);
    h_mPsi3_A_ReCenter_EP[i] = new TH1F(HistName_Full.Data(),HistName_Full.Data(),360,-1.0*TMath::Pi()/3.0,TMath::Pi()/3.0);
    HistName_Full = Form("Psi3_Centrality_%d_ReCenter_B_EP",i);
    h_mPsi3_B_ReCenter_EP[i] = new TH1F(HistName_Full.Data(),HistName_Full.Data(),360,-1.0*TMath::Pi()/3.0,TMath::Pi()/3.0);
    HistName_Full = Form("Psi3_Centrality_%d_ReCenter_Full_EP",i);
    h_mPsi3_Full_ReCenter_EP[i] = new TH1F(HistName_Full.Data(),HistName_Full.Data(),360,-1.0*TMath::Pi()/3.0,TMath::Pi()/3.0);

    // Shift Correction
    HistName_Full = Form("Psi2_Centrality_%d_Shift_A_EP",i);
    h_mPsi2_A_Shift_EP[i] = new TH1F(HistName_Full.Data(),HistName_Full.Data(),360,-1.0*TMath::Pi()/2.0,TMath::Pi()/2.0);
    HistName_Full = Form("Psi2_Centrality_%d_Shift_B_EP",i);
    h_mPsi2_B_Shift_EP[i] = new TH1F(HistName_Full.Data(),HistName_Full.Data(),360,-1.0*TMath::Pi()/2.0,TMath::Pi()/2.0);
    HistName_Full = Form("Psi2_Centrality_%d_Shift_Full_EP",i);
    h_mPsi2_Full_Shift_EP[i] = new TH1F(HistName_Full.Data(),HistName_Full.Data(),360,-1.0*TMath::Pi()/2.0,TMath::Pi()/2.0);
    HistName_Full = Form("Psi3_Centrality_%d_Shift_A_EP",i);
    h_mPsi3_A_Shift_EP[i] = new TH1F(HistName_Full.Data(),HistName_Full.Data(),360,-1.0*TMath::Pi()/3.0,TMath::Pi()/3.0);
    HistName_Full = Form("Psi3_Centrality_%d_Shift_B_EP",i);
    h_mPsi3_B_Shift_EP[i] = new TH1F(HistName_Full.Data(),HistName_Full.Data(),360,-1.0*TMath::Pi()/3.0,TMath::Pi()/3.0);
    HistName_Full = Form("Psi3_Centrality_%d_Shift_Full_EP",i);
    h_mPsi3_Full_Shift_EP[i] = new TH1F(HistName_Full.Data(),HistName_Full.Data(),360,-1.0*TMath::Pi()/3.0,TMath::Pi()/3.0);
    
    // Scalar Product method
    // ReCenter Correction
    HistName_Full = Form("Psi2_Centrality_%d_ReCenter_A_SP",i);
    h_mPsi2_A_ReCenter_SP[i] = new TH1F(HistName_Full.Data(),HistName_Full.Data(),360,-1.0*TMath::Pi()/2.0,TMath::Pi()/2.0);
    HistName_Full = Form("Psi2_Centrality_%d_ReCenter_B_SP",i);
    h_mPsi2_B_ReCenter_SP[i] = new TH1F(HistName_Full.Data(),HistName_Full.Data(),360,-1.0*TMath::Pi()/2.0,TMath::Pi()/2.0);
    HistName_Full = Form("Psi2_Centrality_%d_ReCenter_Full_SP",i);
    h_mPsi2_Full_ReCenter_SP[i] = new TH1F(HistName_Full.Data(),HistName_Full.Data(),360,-1.0*TMath::Pi()/2.0,TMath::Pi()/2.0);
    HistName_Full = Form("Psi3_Centrality_%d_ReCenter_A_SP",i);
    h_mPsi3_A_ReCenter_SP[i] = new TH1F(HistName_Full.Data(),HistName_Full.Data(),360,-1.0*TMath::Pi()/3.0,TMath::Pi()/3.0);
    HistName_Full = Form("Psi3_Centrality_%d_ReCenter_B_SP",i);
    h_mPsi3_B_ReCenter_SP[i] = new TH1F(HistName_Full.Data(),HistName_Full.Data(),360,-1.0*TMath::Pi()/3.0,TMath::Pi()/3.0);
    HistName_Full = Form("Psi3_Centrality_%d_ReCenter_Full_SP",i);
    h_mPsi3_Full_ReCenter_SP[i] = new TH1F(HistName_Full.Data(),HistName_Full.Data(),360,-1.0*TMath::Pi()/3.0,TMath::Pi()/3.0);

    // Shift Correction
    HistName_Full = Form("Psi2_Centrality_%d_Shift_A_SP",i);
    h_mPsi2_A_Shift_SP[i] = new TH1F(HistName_Full.Data(),HistName_Full.Data(),360,-1.0*TMath::Pi()/2.0,TMath::Pi()/2.0);
    HistName_Full = Form("Psi2_Centrality_%d_Shift_B_SP",i);
    h_mPsi2_B_Shift_SP[i] = new TH1F(HistName_Full.Data(),HistName_Full.Data(),360,-1.0*TMath::Pi()/2.0,TMath::Pi()/2.0);
    HistName_Full = Form("Psi2_Centrality_%d_Shift_Full_SP",i);
    h_mPsi2_Full_Shift_SP[i] = new TH1F(HistName_Full.Data(),HistName_Full.Data(),360,-1.0*TMath::Pi()/2.0,TMath::Pi()/2.0);
    HistName_Full = Form("Psi3_Centrality_%d_Shift_A_SP",i);
    h_mPsi3_A_Shift_SP[i] = new TH1F(HistName_Full.Data(),HistName_Full.Data(),360,-1.0*TMath::Pi()/3.0,TMath::Pi()/3.0);
    HistName_Full = Form("Psi3_Centrality_%d_Shift_B_SP",i);
    h_mPsi3_B_Shift_SP[i] = new TH1F(HistName_Full.Data(),HistName_Full.Data(),360,-1.0*TMath::Pi()/3.0,TMath::Pi()/3.0);
    HistName_Full = Form("Psi3_Centrality_%d_Shift_Full_SP",i);
    h_mPsi3_Full_Shift_SP[i] = new TH1F(HistName_Full.Data(),HistName_Full.Data(),360,-1.0*TMath::Pi()/3.0,TMath::Pi()/3.0);
  }
  for(Int_t j = 0; j < 4; j++)
  {
    TString HistName_minBias;
    // Event Plane method
    // ReCenter Correction
    HistName_minBias = Form("Psi2_minBias_EtaGap_%d_ReCenter_East_EP",j);
    h_mPsi2_East_ReCenter_EP[9][j] = new TH1F(HistName_minBias.Data(),HistName_minBias.Data(),360,-1.0*TMath::Pi()/2.0,TMath::Pi()/2.0);
    HistName_minBias = Form("Psi2_minBias_EtaGap_%d_ReCenter_West_EP",j);
    h_mPsi2_West_ReCenter_EP[9][j] = new TH1F(HistName_minBias.Data(),HistName_minBias.Data(),360,-1.0*TMath::Pi()/2.0,TMath::Pi()/2.0);
    HistName_minBias = Form("Psi3_minBias_EtaGap_%d_ReCenter_East_EP",j);
    h_mPsi3_East_ReCenter_EP[9][j] = new TH1F(HistName_minBias.Data(),HistName_minBias.Data(),360,-1.0*TMath::Pi()/3.0,TMath::Pi()/3.0);
    HistName_minBias = Form("Psi3_minBias_EtaGap_%d_ReCenter_West_EP",j);
    h_mPsi3_West_ReCenter_EP[9][j] = new TH1F(HistName_minBias.Data(),HistName_minBias.Data(),360,-1.0*TMath::Pi()/3.0,TMath::Pi()/3.0);

    // Shift Correction
    HistName_minBias = Form("Psi2_minBias_EtaGap_%d_Shift_East_EP",j);
    h_mPsi2_East_Shift_EP[9][j] = new TH1F(HistName_minBias.Data(),HistName_minBias.Data(),360,-1.0*TMath::Pi()/2.0,TMath::Pi()/2.0);
    HistName_minBias = Form("Psi2_minBias_EtaGap_%d_Shift_West_EP",j);
    h_mPsi2_West_Shift_EP[9][j] = new TH1F(HistName_minBias.Data(),HistName_minBias.Data(),360,-1.0*TMath::Pi()/2.0,TMath::Pi()/2.0);
    HistName_minBias = Form("Psi3_minBias_EtaGap_%d_Shift_East_EP",j);
    h_mPsi3_East_Shift_EP[9][j] = new TH1F(HistName_minBias.Data(),HistName_minBias.Data(),360,-1.0*TMath::Pi()/3.0,TMath::Pi()/3.0);
    HistName_minBias = Form("Psi3_minBias_EtaGap_%d_Shift_West_EP",j);
    h_mPsi3_West_Shift_EP[9][j] = new TH1F(HistName_minBias.Data(),HistName_minBias.Data(),360,-1.0*TMath::Pi()/3.0,TMath::Pi()/3.0);

    // Scalor Product method
    // ReCenter Correction
    HistName_minBias = Form("Psi2_minBias_EtaGap_%d_ReCenter_East_SP",j);
    h_mPsi2_East_ReCenter_SP[9][j] = new TH1F(HistName_minBias.Data(),HistName_minBias.Data(),360,-1.0*TMath::Pi()/2.0,TMath::Pi()/2.0);
    HistName_minBias = Form("Psi2_minBias_EtaGap_%d_ReCenter_West_SP",j);
    h_mPsi2_West_ReCenter_SP[9][j] = new TH1F(HistName_minBias.Data(),HistName_minBias.Data(),360,-1.0*TMath::Pi()/2.0,TMath::Pi()/2.0);
    HistName_minBias = Form("Psi3_minBias_EtaGap_%d_ReCenter_East_SP",j);
    h_mPsi3_East_ReCenter_SP[9][j] = new TH1F(HistName_minBias.Data(),HistName_minBias.Data(),360,-1.0*TMath::Pi()/3.0,TMath::Pi()/3.0);
    HistName_minBias = Form("Psi3_minBias_EtaGap_%d_ReCenter_West_SP",j);
    h_mPsi3_West_ReCenter_SP[9][j] = new TH1F(HistName_minBias.Data(),HistName_minBias.Data(),360,-1.0*TMath::Pi()/3.0,TMath::Pi()/3.0);

    // Shift Correction
    HistName_minBias = Form("Psi2_minBias_EtaGap_%d_Shift_East_SP",j);
    h_mPsi2_East_Shift_SP[9][j] = new TH1F(HistName_minBias.Data(),HistName_minBias.Data(),360,-1.0*TMath::Pi()/2.0,TMath::Pi()/2.0);
    HistName_minBias = Form("Psi2_minBias_EtaGap_%d_Shift_West_SP",j);
    h_mPsi2_West_Shift_SP[9][j] = new TH1F(HistName_minBias.Data(),HistName_minBias.Data(),360,-1.0*TMath::Pi()/2.0,TMath::Pi()/2.0);
    HistName_minBias = Form("Psi3_minBias_EtaGap_%d_Shift_East_SP",j);
    h_mPsi3_East_Shift_SP[9][j] = new TH1F(HistName_minBias.Data(),HistName_minBias.Data(),360,-1.0*TMath::Pi()/3.0,TMath::Pi()/3.0);
    HistName_minBias = Form("Psi3_minBias_EtaGap_%d_Shift_West_SP",j);
    h_mPsi3_West_Shift_SP[9][j] = new TH1F(HistName_minBias.Data(),HistName_minBias.Data(),360,-1.0*TMath::Pi()/3.0,TMath::Pi()/3.0);
  }
  TString HistName_minBias_Full;
  // Event Plane method
  // ReCenter Correction
  HistName_minBias_Full = "Psi2_minBias_ReCenter_A_EP";
  h_mPsi2_A_ReCenter_EP[9] = new TH1F(HistName_minBias_Full.Data(),HistName_minBias_Full.Data(),360,-1.0*TMath::Pi()/2.0,TMath::Pi()/2.0);
  HistName_minBias_Full = "Psi2_minBias_ReCenter_B_EP";
  h_mPsi2_B_ReCenter_EP[9] = new TH1F(HistName_minBias_Full.Data(),HistName_minBias_Full.Data(),360,-1.0*TMath::Pi()/2.0,TMath::Pi()/2.0);
  HistName_minBias_Full = "Psi2_minBias_ReCenter_Full_EP";
  h_mPsi2_Full_ReCenter_EP[9] = new TH1F(HistName_minBias_Full.Data(),HistName_minBias_Full.Data(),360,-1.0*TMath::Pi()/2.0,TMath::Pi()/2.0);
  HistName_minBias_Full = "Psi3_minBias_ReCenter_A_EP";
  h_mPsi3_A_ReCenter_EP[9] = new TH1F(HistName_minBias_Full.Data(),HistName_minBias_Full.Data(),360,-1.0*TMath::Pi()/3.0,TMath::Pi()/3.0);
  HistName_minBias_Full = "Psi3_minBias_ReCenter_B_EP";
  h_mPsi3_B_ReCenter_EP[9] = new TH1F(HistName_minBias_Full.Data(),HistName_minBias_Full.Data(),360,-1.0*TMath::Pi()/3.0,TMath::Pi()/3.0);
  HistName_minBias_Full = "Psi3_minBias_ReCenter_Full_EP";
  h_mPsi3_Full_ReCenter_EP[9] = new TH1F(HistName_minBias_Full.Data(),HistName_minBias_Full.Data(),360,-1.0*TMath::Pi()/3.0,TMath::Pi()/3.0);

  // Shift Correction
  HistName_minBias_Full = "Psi2_minBias_Shift_A_EP";
  h_mPsi2_A_Shift_EP[9] = new TH1F(HistName_minBias_Full.Data(),HistName_minBias_Full.Data(),360,-1.0*TMath::Pi()/2.0,TMath::Pi()/2.0);
  HistName_minBias_Full = "Psi2_minBias_Shift_B_EP";
  h_mPsi2_B_Shift_EP[9] = new TH1F(HistName_minBias_Full.Data(),HistName_minBias_Full.Data(),360,-1.0*TMath::Pi()/2.0,TMath::Pi()/2.0);
  HistName_minBias_Full = "Psi2_minBias_Shift_Full_EP";
  h_mPsi2_Full_Shift_EP[9] = new TH1F(HistName_minBias_Full.Data(),HistName_minBias_Full.Data(),360,-1.0*TMath::Pi()/2.0,TMath::Pi()/2.0);
  HistName_minBias_Full = "Psi3_minBias_Shift_A_EP";
  h_mPsi3_A_Shift_EP[9] = new TH1F(HistName_minBias_Full.Data(),HistName_minBias_Full.Data(),360,-1.0*TMath::Pi()/3.0,TMath::Pi()/3.0);
  HistName_minBias_Full = "Psi3_minBias_Shift_B_EP";
  h_mPsi3_B_Shift_EP[9] = new TH1F(HistName_minBias_Full.Data(),HistName_minBias_Full.Data(),360,-1.0*TMath::Pi()/3.0,TMath::Pi()/3.0);
  HistName_minBias_Full = "Psi3_minBias_Shift_Full_EP";
  h_mPsi3_Full_Shift_EP[9] = new TH1F(HistName_minBias_Full.Data(),HistName_minBias_Full.Data(),360,-1.0*TMath::Pi()/3.0,TMath::Pi()/3.0);

  // Scalar Product method
  // ReCenter Correction
  HistName_minBias_Full = "Psi2_minBias_ReCenter_A_SP";
  h_mPsi2_A_ReCenter_SP[9] = new TH1F(HistName_minBias_Full.Data(),HistName_minBias_Full.Data(),360,-1.0*TMath::Pi()/2.0,TMath::Pi()/2.0);
  HistName_minBias_Full = "Psi2_minBias_ReCenter_B_SP";
  h_mPsi2_B_ReCenter_SP[9] = new TH1F(HistName_minBias_Full.Data(),HistName_minBias_Full.Data(),360,-1.0*TMath::Pi()/2.0,TMath::Pi()/2.0);
  HistName_minBias_Full = "Psi2_minBias_ReCenter_Full_SP";
  h_mPsi2_Full_ReCenter_SP[9] = new TH1F(HistName_minBias_Full.Data(),HistName_minBias_Full.Data(),360,-1.0*TMath::Pi()/2.0,TMath::Pi()/2.0);
  HistName_minBias_Full = "Psi3_minBias_ReCenter_A_SP";
  h_mPsi3_A_ReCenter_SP[9] = new TH1F(HistName_minBias_Full.Data(),HistName_minBias_Full.Data(),360,-1.0*TMath::Pi()/3.0,TMath::Pi()/3.0);
  HistName_minBias_Full = "Psi3_minBias_ReCenter_B_SP";
  h_mPsi3_B_ReCenter_SP[9] = new TH1F(HistName_minBias_Full.Data(),HistName_minBias_Full.Data(),360,-1.0*TMath::Pi()/3.0,TMath::Pi()/3.0);
  HistName_minBias_Full = "Psi3_minBias_ReCenter_Full_SP";
  h_mPsi3_Full_ReCenter_SP[9] = new TH1F(HistName_minBias_Full.Data(),HistName_minBias_Full.Data(),360,-1.0*TMath::Pi()/3.0,TMath::Pi()/3.0);

  // Shift Correction
  HistName_minBias_Full = "Psi2_minBias_Shift_A_SP";
  h_mPsi2_A_Shift_SP[9] = new TH1F(HistName_minBias_Full.Data(),HistName_minBias_Full.Data(),360,-1.0*TMath::Pi()/2.0,TMath::Pi()/2.0);
  HistName_minBias_Full = "Psi2_minBias_Shift_B_SP";
  h_mPsi2_B_Shift_SP[9] = new TH1F(HistName_minBias_Full.Data(),HistName_minBias_Full.Data(),360,-1.0*TMath::Pi()/2.0,TMath::Pi()/2.0);
  HistName_minBias_Full = "Psi2_minBias_Shift_Full_SP";
  h_mPsi2_Full_Shift_SP[9] = new TH1F(HistName_minBias_Full.Data(),HistName_minBias_Full.Data(),360,-1.0*TMath::Pi()/2.0,TMath::Pi()/2.0);
  HistName_minBias_Full = "Psi3_minBias_Shift_A_SP";
  h_mPsi3_A_Shift_SP[9] = new TH1F(HistName_minBias_Full.Data(),HistName_minBias_Full.Data(),360,-1.0*TMath::Pi()/3.0,TMath::Pi()/3.0);
  HistName_minBias_Full = "Psi3_minBias_Shift_B_SP";
  h_mPsi3_B_Shift_SP[9] = new TH1F(HistName_minBias_Full.Data(),HistName_minBias_Full.Data(),360,-1.0*TMath::Pi()/3.0,TMath::Pi()/3.0);
  HistName_minBias_Full = "Psi3_minBias_Shift_Full_SP";
  h_mPsi3_Full_Shift_SP[9] = new TH1F(HistName_minBias_Full.Data(),HistName_minBias_Full.Data(),360,-1.0*TMath::Pi()/3.0,TMath::Pi()/3.0);
}

//--------------------------------------------------------------------------------------
// Event Plane method
void StTriFlowHistManger::FillEventPlane_East_EP(Float_t Psi2_R, Float_t Psi2_S, Float_t Psi3_R, Float_t Psi3_S, Int_t cent9, Int_t j)
{
  // Centrality
  h_mPsi2_East_ReCenter_EP[cent9][j]->Fill(Psi2_R);
  h_mPsi2_East_Shift_EP[cent9][j]->Fill(Psi2_S);
  h_mPsi3_East_ReCenter_EP[cent9][j]->Fill(Psi3_R);
  h_mPsi3_East_Shift_EP[cent9][j]->Fill(Psi3_S);

  // minBias
  h_mPsi2_East_ReCenter_EP[9][j]->Fill(Psi2_R);
  h_mPsi2_East_Shift_EP[9][j]->Fill(Psi2_S);
  h_mPsi3_East_ReCenter_EP[9][j]->Fill(Psi3_R);
  h_mPsi3_East_Shift_EP[9][j]->Fill(Psi3_S);
}

void StTriFlowHistManger::FillEventPlane_West_EP(Float_t Psi2_R, Float_t Psi2_S, Float_t Psi3_R, Float_t Psi3_S, Int_t cent9, Int_t j)
{
  // Centrality
  h_mPsi2_West_ReCenter_EP[cent9][j]->Fill(Psi2_R);
  h_mPsi2_West_Shift_EP[cent9][j]->Fill(Psi2_S);
  h_mPsi3_West_ReCenter_EP[cent9][j]->Fill(Psi3_R);
  h_mPsi3_West_Shift_EP[cent9][j]->Fill(Psi3_S);

  // minBias
  h_mPsi2_West_ReCenter_EP[9][j]->Fill(Psi2_R);
  h_mPsi2_West_Shift_EP[9][j]->Fill(Psi2_S);
  h_mPsi3_West_ReCenter_EP[9][j]->Fill(Psi3_R);
  h_mPsi3_West_Shift_EP[9][j]->Fill(Psi3_S);
}

void StTriFlowHistManger::FillEventPlane_A_EP(Float_t Psi2_R, Float_t Psi2_S, Float_t Psi3_R, Float_t Psi3_S, Int_t cent9)
{
  // Centrality
  h_mPsi2_A_ReCenter_EP[cent9]->Fill(Psi2_R);
  h_mPsi2_A_Shift_EP[cent9]->Fill(Psi2_S);
  h_mPsi3_A_ReCenter_EP[cent9]->Fill(Psi3_R);
  h_mPsi3_A_Shift_EP[cent9]->Fill(Psi3_S);

  // minBias
  h_mPsi2_A_ReCenter_EP[9]->Fill(Psi2_R);
  h_mPsi2_A_Shift_EP[9]->Fill(Psi2_S);
  h_mPsi3_A_ReCenter_EP[9]->Fill(Psi3_R);
  h_mPsi3_A_Shift_EP[9]->Fill(Psi3_S);
}

void StTriFlowHistManger::FillEventPlane_B_EP(Float_t Psi2_R, Float_t Psi2_S, Float_t Psi3_R, Float_t Psi3_S, Int_t cent9)
{
  // Centrality
  h_mPsi2_B_ReCenter_EP[cent9]->Fill(Psi2_R);
  h_mPsi2_B_Shift_EP[cent9]->Fill(Psi2_S);
  h_mPsi3_B_ReCenter_EP[cent9]->Fill(Psi3_R);
  h_mPsi3_B_Shift_EP[cent9]->Fill(Psi3_S);

  // minBias
  h_mPsi2_B_ReCenter_EP[9]->Fill(Psi2_R);
  h_mPsi2_B_Shift_EP[9]->Fill(Psi2_S);
  h_mPsi3_B_ReCenter_EP[9]->Fill(Psi3_R);
  h_mPsi3_B_Shift_EP[9]->Fill(Psi3_S);
}

void StTriFlowHistManger::FillEventPlane_Full_EP(Float_t Psi2_R, Float_t Psi2_S, Float_t Psi3_R, Float_t Psi3_S, Int_t cent9)
{
  // Centrality
  h_mPsi2_Full_ReCenter_EP[cent9]->Fill(Psi2_R);
  h_mPsi2_Full_Shift_EP[cent9]->Fill(Psi2_S);
  h_mPsi3_Full_ReCenter_EP[cent9]->Fill(Psi3_R);
  h_mPsi3_Full_Shift_EP[cent9]->Fill(Psi3_S);

  // minBias
  h_mPsi2_Full_ReCenter_EP[9]->Fill(Psi2_R);
  h_mPsi2_Full_Shift_EP[9]->Fill(Psi2_S);
  h_mPsi3_Full_ReCenter_EP[9]->Fill(Psi3_R);
  h_mPsi3_Full_Shift_EP[9]->Fill(Psi3_S);
}

//--------------------------------------------------------------------------------------
// Scalor Product method
void StTriFlowHistManger::FillEventPlane_East_SP(Float_t Psi2_R, Float_t Psi2_S, Float_t Psi3_R, Float_t Psi3_S, Int_t cent9, Int_t j)
{
  // Centrality
  h_mPsi2_East_ReCenter_SP[cent9][j]->Fill(Psi2_R);
  h_mPsi2_East_Shift_SP[cent9][j]->Fill(Psi2_S);
  h_mPsi3_East_ReCenter_SP[cent9][j]->Fill(Psi3_R);
  h_mPsi3_East_Shift_SP[cent9][j]->Fill(Psi3_S);

  // minBias
  h_mPsi2_East_ReCenter_SP[9][j]->Fill(Psi2_R);
  h_mPsi2_East_Shift_SP[9][j]->Fill(Psi2_S);
  h_mPsi3_East_ReCenter_SP[9][j]->Fill(Psi3_R);
  h_mPsi3_East_Shift_SP[9][j]->Fill(Psi3_S);
}

void StTriFlowHistManger::FillEventPlane_West_SP(Float_t Psi2_R, Float_t Psi2_S, Float_t Psi3_R, Float_t Psi3_S, Int_t cent9, Int_t j)
{
  // Centrality
  h_mPsi2_West_ReCenter_SP[cent9][j]->Fill(Psi2_R);
  h_mPsi2_West_Shift_SP[cent9][j]->Fill(Psi2_S);
  h_mPsi3_West_ReCenter_SP[cent9][j]->Fill(Psi3_R);
  h_mPsi3_West_Shift_SP[cent9][j]->Fill(Psi3_S);

  // minBias
  h_mPsi2_West_ReCenter_SP[9][j]->Fill(Psi2_R);
  h_mPsi2_West_Shift_SP[9][j]->Fill(Psi2_S);
  h_mPsi3_West_ReCenter_SP[9][j]->Fill(Psi3_R);
  h_mPsi3_West_Shift_SP[9][j]->Fill(Psi3_S);
}

void StTriFlowHistManger::FillEventPlane_A_SP(Float_t Psi2_R, Float_t Psi2_S, Float_t Psi3_R, Float_t Psi3_S, Int_t cent9)
{
  // Centrality
  h_mPsi2_A_ReCenter_SP[cent9]->Fill(Psi2_R);
  h_mPsi2_A_Shift_SP[cent9]->Fill(Psi2_S);
  h_mPsi3_A_ReCenter_SP[cent9]->Fill(Psi3_R);
  h_mPsi3_A_Shift_SP[cent9]->Fill(Psi3_S);

  // minBias
  h_mPsi2_A_ReCenter_SP[9]->Fill(Psi2_R);
  h_mPsi2_A_Shift_SP[9]->Fill(Psi2_S);
  h_mPsi3_A_ReCenter_SP[9]->Fill(Psi3_R);
  h_mPsi3_A_Shift_SP[9]->Fill(Psi3_S);
}

void StTriFlowHistManger::FillEventPlane_B_SP(Float_t Psi2_R, Float_t Psi2_S, Float_t Psi3_R, Float_t Psi3_S, Int_t cent9)
{
  // Centrality
  h_mPsi2_B_ReCenter_SP[cent9]->Fill(Psi2_R);
  h_mPsi2_B_Shift_SP[cent9]->Fill(Psi2_S);
  h_mPsi3_B_ReCenter_SP[cent9]->Fill(Psi3_R);
  h_mPsi3_B_Shift_SP[cent9]->Fill(Psi3_S);

  // minBias
  h_mPsi2_B_ReCenter_SP[9]->Fill(Psi2_R);
  h_mPsi2_B_Shift_SP[9]->Fill(Psi2_S);
  h_mPsi3_B_ReCenter_SP[9]->Fill(Psi3_R);
  h_mPsi3_B_Shift_SP[9]->Fill(Psi3_S);
}

void StTriFlowHistManger::FillEventPlane_Full_SP(Float_t Psi2_R, Float_t Psi2_S, Float_t Psi3_R, Float_t Psi3_S, Int_t cent9)
{
  // Centrality
  h_mPsi2_Full_ReCenter_SP[cent9]->Fill(Psi2_R);
  h_mPsi2_Full_Shift_SP[cent9]->Fill(Psi2_S);
  h_mPsi3_Full_ReCenter_SP[cent9]->Fill(Psi3_R);
  h_mPsi3_Full_Shift_SP[cent9]->Fill(Psi3_S);

  // minBias
  h_mPsi2_Full_ReCenter_SP[9]->Fill(Psi2_R);
  h_mPsi2_Full_Shift_SP[9]->Fill(Psi2_S);
  h_mPsi3_Full_ReCenter_SP[9]->Fill(Psi3_R);
  h_mPsi3_Full_Shift_SP[9]->Fill(Psi3_S);
}

//--------------------------------------------------------------------------------------

void StTriFlowHistManger::WriteEventPlane()
{
  for(Int_t i = 0; i < 10; i++)
  {
    for(Int_t j = 0; j < 4; j++)
    {
      h_mPsi2_East_ReCenter_EP[i][j]->Write();
      h_mPsi2_West_ReCenter_EP[i][j]->Write();
      h_mPsi3_East_ReCenter_EP[i][j]->Write();
      h_mPsi3_West_ReCenter_EP[i][j]->Write();

      h_mPsi2_East_Shift_EP[i][j]->Write();
      h_mPsi2_West_Shift_EP[i][j]->Write();
      h_mPsi3_East_Shift_EP[i][j]->Write();
      h_mPsi3_West_Shift_EP[i][j]->Write();

      h_mPsi2_East_ReCenter_SP[i][j]->Write();
      h_mPsi2_West_ReCenter_SP[i][j]->Write();
      h_mPsi3_East_ReCenter_SP[i][j]->Write();
      h_mPsi3_West_ReCenter_SP[i][j]->Write();

      h_mPsi2_East_Shift_SP[i][j]->Write();
      h_mPsi2_West_Shift_SP[i][j]->Write();
      h_mPsi3_East_Shift_SP[i][j]->Write();
      h_mPsi3_West_Shift_SP[i][j]->Write();
    }
    h_mPsi2_A_ReCenter_EP[i]->Write();
    h_mPsi2_B_ReCenter_EP[i]->Write();
    h_mPsi2_Full_ReCenter_EP[i]->Write();
    h_mPsi3_A_ReCenter_EP[i]->Write();
    h_mPsi3_B_ReCenter_EP[i]->Write();
    h_mPsi3_Full_ReCenter_EP[i]->Write();

    h_mPsi2_A_Shift_EP[i]->Write();
    h_mPsi2_B_Shift_EP[i]->Write();
    h_mPsi2_Full_Shift_EP[i]->Write();
    h_mPsi3_A_Shift_EP[i]->Write();
    h_mPsi3_B_Shift_EP[i]->Write();
    h_mPsi3_Full_Shift_EP[i]->Write();

    h_mPsi2_A_ReCenter_SP[i]->Write();
    h_mPsi2_B_ReCenter_SP[i]->Write();
    h_mPsi2_Full_ReCenter_SP[i]->Write();
    h_mPsi3_A_ReCenter_SP[i]->Write();
    h_mPsi3_B_ReCenter_SP[i]->Write();
    h_mPsi3_Full_ReCenter_SP[i]->Write();

    h_mPsi2_A_Shift_SP[i]->Write();
    h_mPsi2_B_Shift_SP[i]->Write();
    h_mPsi2_Full_Shift_SP[i]->Write();
    h_mPsi3_A_Shift_SP[i]->Write();
    h_mPsi3_B_Shift_SP[i]->Write();
    h_mPsi3_Full_Shift_SP[i]->Write();
  }
}

//--------------------------------------------------------------------------------------


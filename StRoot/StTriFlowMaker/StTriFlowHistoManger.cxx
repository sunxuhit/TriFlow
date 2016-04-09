#include "StTriFlowHistoManger.h"
#include "StTriFlowConstants.h"
#include "StTriFlowCut.h"
#include "StRoot/StPicoDstMaker/StPicoTrack.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TString.h"
#include "TMath.h"

ClassImp(StTriFlowHistoManger)

//-------------------------------------------------------------------------------------------

StTriFlowHistoManger::StTriFlowHistoManger(Int_t energy)
{
  mEnergy = energy;
}

//-------------------------------------------------------------------------------------------

StTriFlowHistoManger::~StTriFlowHistoManger()
{
}

//-------------------------------------------------------------------------------------------

void StTriFlowHistoManger::InitHist()
{
  // flow
  for(Int_t i = 0; i < TriFlow::pt_total; i++) // pt bin
  {
    for(Int_t j = TriFlow::Centrality_start; j < TriFlow::Centrality_stop; j++) // centrality bin
    {
      for(Int_t k = TriFlow::Charge_start; k < TriFlow::Charge_stop; k++) // charge bin
      {
	for(Int_t l = TriFlow::EtaGap_start; l < TriFlow::EtaGap_stop; l++) // eta gap bin
	{
	  for(Int_t m = 0; m < 7; m++) // phi-psi bin
	  {
	    TString HistName;
	    HistName = Form("pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_phi_Psi_%d_2nd_EP",i,j,k,l,m);
	    h_mMass2_nSigmaPion2_EP[i][j][k][l][m] = new TH2F(HistName.Data(),HistName.Data(),400,TriFlow::x_low[i],TriFlow::x_up[i],400,TriFlow::y_low[i],TriFlow::y_up[i]);
	    h_mMass2_nSigmaPion2_EP[i][j][k][l][m]->Sumw2();
//	    cout << HistName.Data() << endl;
	    HistName = Form("pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_phi_Psi_%d_3rd_EP",i,j,k,l,m);
	    h_mMass2_nSigmaPion3_EP[i][j][k][l][m] = new TH2F(HistName.Data(),HistName.Data(),400,TriFlow::x_low[i],TriFlow::x_up[i],400,TriFlow::y_low[i],TriFlow::y_up[i]);
	    h_mMass2_nSigmaPion3_EP[i][j][k][l][m]->Sumw2();
//	    cout << HistName.Data() << endl;
	  }
	}
      }
    }
  }

  // raw pt spectra
  for(Int_t i = 0; i < TriFlow::pt_total; i++) // pt bin
  {
    for(Int_t j = TriFlow::Centrality_start; j < TriFlow::Centrality_stop; j++) // centrality bin
    {
      for(Int_t k = TriFlow::Charge_start; k < TriFlow::Charge_stop; k++) // charge bin
      {
	for(Int_t l = TriFlow::EtaGap_start; l < TriFlow::EtaGap_stop; l++) // eta gap bin
	{
	  TString SpecName;
	  SpecName = Form("Spectra_pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_pik_2nd",i,j,k,l);
	  h_pt_spectra2_pik[i][j][k][l] = new TH2F(SpecName.Data(),SpecName.Data(),400,TriFlow::x_low[i],TriFlow::x_up[i],400,TriFlow::y_low[i],TriFlow::y_up[i]);
	  h_pt_spectra2_pik[i][j][k][l]->Sumw2();

	  SpecName = Form("Spectra_pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_pik_3rd",i,j,k,l);
	  h_pt_spectra3_pik[i][j][k][l] = new TH2F(SpecName.Data(),SpecName.Data(),400,TriFlow::x_low[i],TriFlow::x_up[i],400,TriFlow::y_low[i],TriFlow::y_up[i]);
	  h_pt_spectra3_pik[i][j][k][l]->Sumw2();
	}
      }
    }
  }
  TString HistName_ToF;
  HistName_ToF = "ToFYLocal_Mass2";
  h_mToFYLocal_Mass2 = new TH2F(HistName_ToF.Data(),HistName_ToF.Data(),400,-0.1,1.5,400,-6.0,6.0);
  HistName_ToF = "ToFZLocal_Mass2";
  h_mToFZLocal_Mass2 = new TH2F(HistName_ToF.Data(),HistName_ToF.Data(),400,-0.1,1.5,400,-6.0,6.0);
  mTriFlowCut = new StTriFlowCut(mEnergy);

  TString HistName_Counter = "Mass2_pt";
  h_mMass2_pt = new TH2F(HistName_Counter.Data(),HistName_Counter.Data(),100,0.15,5.0,100,-0.2,1.2);

//  h_phi_psi2 = new TH1F("phi_psi2","phi_psi2",200,-3.0*TMath::Pi()/2.0,3.0*TMath::Pi()/2.0);
//  h_phi_psi3 = new TH1F("phi_psi3","phi_psi3",200,-4.0*TMath::Pi()/3.0,4.0*TMath::Pi()/3.0);
//  h_yield_phi_psi2 = new TH1F("yield_phi_psi2","yield_phi_psi2",7,0.0,TMath::Pi()/2.0);
//  h_yield_phi_psi3 = new TH1F("yield_phi_psi3","yield_phi_psi3",7,0.0,TMath::Pi()/3.0);
}

//-------------------------------------------------------------------------------------------

void StTriFlowHistoManger::InitProton()
{
  // flow
  for(Int_t i = 0; i < TriFlow::pt_total; i++) // pt bin
  {
    for(Int_t j = TriFlow::Centrality_start; j < TriFlow::Centrality_stop; j++) // centrality bin
    {
      for(Int_t k = TriFlow::Charge_start; k < TriFlow::Charge_stop; k++) // charge bin
      {
	for(Int_t l = TriFlow::EtaGap_start; l < TriFlow::EtaGap_stop; l++) // eta gap bin
	{
	  for(Int_t m = 0; m < 7; m++) // phi-psi bin
	  {
	    for(Int_t i_cut = 0; i_cut < 18; i_cut++)
	    {
	      TString KEY_Proton;

	      KEY_Proton = Form("pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_phi_Psi_%d_2nd_Proton_SysError_%d",i,j,k,l,m,i_cut);
	      h_mMass2_Proton2_EP[KEY_Proton] = new TH1F(KEY_Proton.Data(),KEY_Proton.Data(),400,-0.3,1.7);
	      h_mMass2_Proton2_EP[KEY_Proton]->Sumw2();

	      KEY_Proton = Form("pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_phi_Psi_%d_3rd_Proton_SysError_%d",i,j,k,l,m,i_cut);
	      h_mMass2_Proton3_EP[KEY_Proton] = new TH1F(KEY_Proton.Data(),KEY_Proton.Data(),400,-0.3,1.7);
	      h_mMass2_Proton3_EP[KEY_Proton]->Sumw2();
	    }
	  }
	}
      }
    }
  }

  // raw pt spectra
  for(Int_t i = 0; i < TriFlow::pt_total; i++) // pt bin
  {
    for(Int_t j = TriFlow::Centrality_start; j < TriFlow::Centrality_stop; j++) // centrality bin
    {
      for(Int_t k = TriFlow::Charge_start; k < TriFlow::Charge_stop; k++) // charge bin
      {
	for(Int_t l = TriFlow::EtaGap_start; l < TriFlow::EtaGap_stop; l++) // eta gap bin
	{
	  for(Int_t i_cut = 0; i_cut < 18; i_cut++)
	  {
	    TString KEY_Proton_Spec;
	    KEY_Proton_Spec = Form("Spectra_pt_%d_low_Centrality_%d_Charge_%d_EtaGap_%d_Proton_2nd_SysError_%d",i,j,k,l,i_cut);
	    h_pt_spectra2_proton[KEY_Proton_Spec] = new TH1F(KEY_Proton_Spec.Data(),KEY_Proton_Spec.Data(),400,-0.3,1.7);
	    h_pt_spectra2_proton[KEY_Proton_Spec]->Sumw2();

	    KEY_Proton_Spec = Form("Spectra_pt_%d_high_Centrality_%d_Charge_%d_EtaGap_%d_Proton_2nd_SysError_%d",i,j,k,l,i_cut);
	    h_pt_spectra2_proton[KEY_Proton_Spec] = new TH1F(KEY_Proton_Spec.Data(),KEY_Proton_Spec.Data(),400,-0.3,1.7);
	    h_pt_spectra2_proton[KEY_Proton_Spec]->Sumw2();

	    KEY_Proton_Spec = Form("Spectra_pt_%d_low_Centrality_%d_Charge_%d_EtaGap_%d_Proton_3rd_SysError_%d",i,j,k,l,i_cut);
	    h_pt_spectra3_proton[KEY_Proton_Spec] = new TH1F(KEY_Proton_Spec.Data(),KEY_Proton_Spec.Data(),400,-0.3,1.7);
	    h_pt_spectra3_proton[KEY_Proton_Spec]->Sumw2();

	    KEY_Proton_Spec = Form("Spectra_pt_%d_high_Centrality_%d_Charge_%d_EtaGap_%d_Proton_3rd_SysError_%d",i,j,k,l,i_cut);
	    h_pt_spectra3_proton[KEY_Proton_Spec] = new TH1F(KEY_Proton_Spec.Data(),KEY_Proton_Spec.Data(),400,-0.3,1.7);
	    h_pt_spectra3_proton[KEY_Proton_Spec]->Sumw2();
	  }
	}
      }
    }
  }
  for(Int_t l = TriFlow::EtaGap_start; l < TriFlow::EtaGap_stop; l++) // eta gap bin
  {
    TString HistName;
    HistName = Form("dEdx_pq_EtaGap_%d_before",l);
    h_mDEdx_pq_before[l] = new TH2F(HistName.Data(),HistName.Data(),100,-4.0,4.0,100,0,10);
    HistName = Form("dEdx_pq_EtaGap_%d_after",l);
    h_mDEdx_pq_after[l] = new TH2F(HistName.Data(),HistName.Data(),100,-4.0,4.0,100,0,10);
    HistName = Form("Mass2_pq_EtaGap_%d_before",l);
    h_mMass2_pq_before[l] = new TH2F(HistName.Data(),HistName.Data(),100,-4.0,4.0,100,-0.3,1.7);
    HistName = Form("Mass2_pq_EtaGap_%d_after",l);
    h_mMass2_pq_after[l] = new TH2F(HistName.Data(),HistName.Data(),100,-4.0,4.0,100,-0.3,1.7);
  }

  h_mRefMult = new TH1F("RefMult","RefMult",1000,-0.5,999.5);
  h_mVz = new TH1F("VertexZ","VertexZ",500,-80.0,80.0);
}

//-------------------------------------------------------------------------------------------
void StTriFlowHistoManger::InitYields_nSigPion()
{
  for(Int_t j = 0; j < 9; j++) // centrality bin
  {
    for(Int_t k = TriFlow::Charge_start; k < TriFlow::Charge_stop; k++) // charge bin
    {
      for(Int_t l = TriFlow::EtaGap_start; l < TriFlow::EtaGap_stop; l++) // eta gap bin
      {
	TString HistName = Form("Centrality_%d_Charge_%d_EtaGap_%d_Yields_PiK_EP",j,k,l);
	h_mMass2_nSigmaPion_Yields_PiK_EP[j][k][l] = new TH2F(HistName.Data(),HistName.Data(),400,TriFlow::x_low[12],TriFlow::x_up[12],400,TriFlow::y_low[12],TriFlow::y_up[12]);
	h_mMass2_nSigmaPion_Yields_PiK_EP[j][k][l]->Sumw2();
      }
    }
  }
}

void StTriFlowHistoManger::InitYields_Proton()
{
  for(Int_t j = 0; j < 9; j++) // centrality bin
  {
    for(Int_t k = TriFlow::Charge_start; k < TriFlow::Charge_stop; k++) // charge bin
    {
      for(Int_t l = TriFlow::EtaGap_start; l < TriFlow::EtaGap_stop; l++) // eta gap bin
      {
	for(Int_t i_cut = 0; i_cut < 18; i_cut++)
	{
	  TString KEY_Proton_Yield = Form("Centrality_%d_Charge_%d_EtaGap_%d_Yields_Proton_SysError_%d",j,k,l,i_cut);
	  h_mMass2_Yields_Proton_EP[KEY_Proton_Yield] = new TH1F(KEY_Proton_Yield.Data(),KEY_Proton_Yield.Data(),400,-0.3,1.7);
	  h_mMass2_Yields_Proton_EP[KEY_Proton_Yield]->Sumw2();
	}
      }
    }
  }
}
//-------------------------------------------------------------------------------------------

//-------------------------------------------------------------------------------------------
void StTriFlowHistoManger::InitQA_Detector()
{
  h_mDEdx = new TH2F("h_mDEdx","h_mDEdx",1000,0,4.0,1000,0,40);
  h_mMass2 = new TH2F("h_mMass2","h_mMass2",1000,0,4.0,1000,-0.3,1.7);
}
//-------------------------------------------------------------------------------------------

void StTriFlowHistoManger::FillHist(Float_t pt, Int_t Cent9, Int_t charge_bin, Int_t eta_gap, Float_t phi_psi2, Float_t Res2, Float_t phi_psi3, Float_t Res3, Float_t New_X, Float_t New_Y, Double_t reweight)
{
  if(Res2 > 0.0)
  {
    for(Int_t i = 0; i < TriFlow::pt_total; i++) // pt_bin
    {
      if(pt > TriFlow::pt_low[i] && pt < TriFlow::pt_up[i])
      {
	for(Int_t j = TriFlow::Centrality_start; j < TriFlow::Centrality_stop; j++) // centrality bin
	{
	  if(Cent9 >= TriFlow::cent_low[j] && Cent9 <= TriFlow::cent_up[j])
	  {
	    for(Int_t psi_bin = 0; psi_bin < 3; psi_bin++)
	    {
	      if(phi_psi2 >= TriFlow::Psi2_low[psi_bin] && phi_psi2 < TriFlow::Psi2_up[psi_bin])
	      {
//		cout << "phi_psi2 = " << phi_psi2 << endl;
		Double_t phi_psi2_final = phi_psi2 - (psi_bin-1)*2.0*TMath::Pi()/2.0;
//		cout << "phi_psi2_final = " << phi_psi2_final << endl;
		for(Int_t m = 0; m < 7; m++) // phi-psi2 bin
		{
		  if(TMath::Abs(phi_psi2_final) >= TriFlow::phi_Psi2_low[m] && TMath::Abs(phi_psi2_final) < TriFlow::phi_Psi2_up[m])
		  {
		    // flow
		    h_mMass2_nSigmaPion2_EP[i][j][charge_bin][eta_gap][m]->Fill(New_X,New_Y,(reweight/Res2));
		    // raw pt spectra
		    h_pt_spectra2_pik[i][j][charge_bin][eta_gap]->Fill(New_X,New_Y,reweight);
//		    cout << "m = " << m << endl;
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }

  if(Res3 > 0.0)
  {
    for(Int_t i = 0; i < TriFlow::pt_total; i++) // pt_bin
    {
      if(pt > TriFlow::pt_low[i] && pt < TriFlow::pt_up[i])
      {
	for(Int_t j = TriFlow::Centrality_start; j < TriFlow::Centrality_stop; j++) // centrality bin
	{
	  if(Cent9 >= TriFlow::cent_low[j] && Cent9 <= TriFlow::cent_up[j])
	  {
	    for(Int_t psi_bin = 0; psi_bin < 5; psi_bin++)
	    {
	      if(phi_psi3 >= TriFlow::Psi3_low[psi_bin] && phi_psi3 < TriFlow::Psi3_up[psi_bin])
	      {
//		cout << "phi_psi3 = " << phi_psi3 << endl;
		Double_t phi_psi3_final = phi_psi3 - (psi_bin-2)*2.0*TMath::Pi()/3.0;
//		cout << "phi_psi3_final = " << phi_psi3_final << endl;
		for(Int_t m = 0; m < 7; m++) // phi-psi3 bin
		{
		  if(TMath::Abs(phi_psi3_final) >= TriFlow::phi_Psi3_low[m] && TMath::Abs(phi_psi3_final) < TriFlow::phi_Psi3_up[m])
		  {
		    // flow
		    h_mMass2_nSigmaPion3_EP[i][j][charge_bin][eta_gap][m]->Fill(New_X,New_Y,(reweight/Res3));
		    // raw pt spectra
		    h_pt_spectra3_pik[i][j][charge_bin][eta_gap]->Fill(New_X,New_Y,reweight);
//		    cout << "phi_psi3_bin = " << m << endl;
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
}

//-------------------------------------------------------------------------------------------

void StTriFlowHistoManger::FillProton(Float_t pt, Int_t Cent9, Int_t charge_bin, Int_t eta_gap, Float_t phi_psi2, Float_t Res2, Float_t phi_psi3, Float_t Res3, Float_t Mass2, Double_t reweight, Int_t i_cut)
{
  if(Res2 > 0.0)
  {
    for(Int_t i = 0; i < TriFlow::pt_total; i++) // pt_bin
    {
      if(pt > TriFlow::pt_low[i] && pt < TriFlow::pt_up[i])
      {
	for(Int_t j = TriFlow::Centrality_start; j < TriFlow::Centrality_stop; j++) // centrality bin
	{
	  if(Cent9 >= TriFlow::cent_low[j] && Cent9 <= TriFlow::cent_up[j])
	  {
	    for(Int_t psi_bin = 0; psi_bin < 3; psi_bin++)
	    {
	      if(phi_psi2 >= TriFlow::Psi2_low[psi_bin] && phi_psi2 < TriFlow::Psi2_up[psi_bin])
	      {
//		cout << "phi_psi2 = " << phi_psi2 << endl;
		Double_t phi_psi2_final = phi_psi2 - (psi_bin-1)*2.0*TMath::Pi()/2.0;
//		cout << "phi_psi2_final = " << phi_psi2_final << endl;
		for(Int_t m = 0; m < 7; m++) // phi-psi2 bin
		{
		  if(TMath::Abs(phi_psi2_final) >= TriFlow::phi_Psi2_low[m] && TMath::Abs(phi_psi2_final) < TriFlow::phi_Psi2_up[m])
		  {
		    // flow
		    TString KEY_Proton = Form("pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_phi_Psi_%d_2nd_Proton_SysError_%d",i,j,charge_bin,eta_gap,m,i_cut);
		    h_mMass2_Proton2_EP[KEY_Proton]->Fill(Mass2,(reweight/Res2));
		    // raw pt spectra
		    if(pt < 0.5*(TriFlow::pt_low[i]+TriFlow::pt_up[i]))
		    {
		      TString KEY_Proton_Spec = Form("Spectra_pt_%d_low_Centrality_%d_Charge_%d_EtaGap_%d_Proton_2nd_SysError_%d",i,j,charge_bin,eta_gap,i_cut);
		      h_pt_spectra2_proton[KEY_Proton_Spec]->Fill(Mass2,reweight);
		    }
		    else
		    {
		      TString KEY_Proton_Spec = Form("Spectra_pt_%d_high_Centrality_%d_Charge_%d_EtaGap_%d_Proton_2nd_SysError_%d",i,j,charge_bin,eta_gap,i_cut);
		      h_pt_spectra2_proton[KEY_Proton_Spec]->Fill(Mass2,reweight);
		    }
//		    cout << "m = " << m << endl;
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }

  if(Res3 > 0.0)
  {
    for(Int_t i = 0; i < TriFlow::pt_total; i++) // pt_bin
    {
      if(pt > TriFlow::pt_low[i] && pt < TriFlow::pt_up[i])
      {
	for(Int_t j = TriFlow::Centrality_start; j < TriFlow::Centrality_stop; j++) // centrality bin
	{
	  if(Cent9 >= TriFlow::cent_low[j] && Cent9 <= TriFlow::cent_up[j])
	  {
	    for(Int_t psi_bin = 0; psi_bin < 5; psi_bin++)
	    {
	      if(phi_psi3 >= TriFlow::Psi3_low[psi_bin] && phi_psi3 < TriFlow::Psi3_up[psi_bin])
	      {
//		cout << "phi_psi3 = " << phi_psi3 << endl;
		Double_t phi_psi3_final = phi_psi3 - (psi_bin-2)*2.0*TMath::Pi()/3.0;
//		cout << "phi_psi3_final = " << phi_psi3_final << endl;
		for(Int_t m = 0; m < 7; m++) // phi-psi3 bin
		{
		  if(TMath::Abs(phi_psi3_final) >= TriFlow::phi_Psi3_low[m] && TMath::Abs(phi_psi3_final) < TriFlow::phi_Psi3_up[m])
		  {
		    // flow
		    TString KEY_Proton = Form("pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_phi_Psi_%d_3rd_Proton_SysError_%d",i,j,charge_bin,eta_gap,m,i_cut);
		    h_mMass2_Proton3_EP[KEY_Proton]->Fill(Mass2,(reweight/Res3));
		    // raw pt spectra
		    if(pt < 0.5*(TriFlow::pt_low[i]+TriFlow::pt_up[i]))
		    {
		      TString KEY_Proton_Spec = Form("Spectra_pt_%d_low_Centrality_%d_Charge_%d_EtaGap_%d_Proton_3rd_SysError_%d",i,j,charge_bin,eta_gap,i_cut);
		      h_pt_spectra3_proton[KEY_Proton_Spec]->Fill(Mass2,reweight);
		    }
		    else
		    {
		      TString KEY_Proton_Spec = Form("Spectra_pt_%d_high_Centrality_%d_Charge_%d_EtaGap_%d_Proton_3rd_SysError_%d",i,j,charge_bin,eta_gap,i_cut);
		      h_pt_spectra3_proton[KEY_Proton_Spec]->Fill(Mass2,reweight);
		    }
//		    cout << "phi_psi3_bin = " << m << endl;
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
}

//-------------------------------------------------------------------------------------------

void StTriFlowHistoManger::FillYields_PiK(Int_t Cent9, Int_t charge_bin, Int_t eta_gap, Float_t New_X, Float_t New_Y, Double_t reweight)
{
  h_mMass2_nSigmaPion_Yields_PiK_EP[Cent9][charge_bin][eta_gap]->Fill(New_X,New_Y,reweight);
}

void StTriFlowHistoManger::FillYields_Proton(Int_t Cent9, Int_t charge_bin, Int_t eta_gap, Float_t Mass2, Double_t reweight, Int_t i_cut)
{
  TString KEY_Proton_Yield = Form("Centrality_%d_Charge_%d_EtaGap_%d_Yields_Proton_SysError_%d",Cent9,charge_bin,eta_gap,i_cut);
  h_mMass2_Yields_Proton_EP[KEY_Proton_Yield]->Fill(Mass2,reweight);
}
//-------------------------------------------------------------------------------------------

void StTriFlowHistoManger::FillQA_before(Int_t eta_gap, Float_t Mass2, Float_t dEdx, Float_t pq)
{
  h_mDEdx_pq_before[eta_gap]->Fill(pq,dEdx);
  h_mMass2_pq_before[eta_gap]->Fill(pq,Mass2);
}

void StTriFlowHistoManger::FillQA_after(Int_t eta_gap, Float_t Mass2, Float_t dEdx, Float_t pq)
{
  h_mDEdx_pq_after[eta_gap]->Fill(pq,dEdx);
  h_mMass2_pq_after[eta_gap]->Fill(pq,Mass2);
}

void StTriFlowHistoManger::FillQA_Detector(Float_t dEdx, Float_t Mass2, Float_t p)
{
  h_mDEdx->Fill(p,dEdx);
  h_mMass2->Fill(p,Mass2);
}

//-------------------------------------------------------------------------------------------

void StTriFlowHistoManger::FillToFLocal(StPicoTrack *track)
{
  Float_t ToFYLocal = track->btofYLocal();
  Float_t ToFZLocal = track->btofZLocal();
  Float_t Mass2 = mTriFlowCut->getMass2(track);
  h_mToFYLocal_Mass2->Fill(Mass2,ToFYLocal);
  h_mToFZLocal_Mass2->Fill(Mass2,ToFZLocal);

  Float_t pt = track->pMom().perp();
  h_mMass2_pt->Fill(pt,Mass2);
}

//-------------------------------------------------------------------------------------------

void StTriFlowHistoManger::FillQA_Event(Int_t RefMult, Float_t Vz)
{
  h_mRefMult->Fill(RefMult);
  h_mVz->Fill(Vz);
}
//-------------------------------------------------------------------------------------------

void StTriFlowHistoManger::WriteHist()
{
  // flow
  for(Int_t i = 0; i < TriFlow::pt_total; i++) // pt bin
  {
    for(Int_t j = TriFlow::Centrality_start; j < TriFlow::Centrality_stop; j++) // centrality bin
    {
      for(Int_t k = TriFlow::Charge_start; k < TriFlow::Charge_stop; k++) // charge bin
      {
	for(Int_t l = TriFlow::EtaGap_start; l < TriFlow::EtaGap_stop; l++) // eta gap bin
	{
	  for(Int_t m = 0; m < 7; m ++) // phi-psi bin
	  {
	    h_mMass2_nSigmaPion2_EP[i][j][k][l][m]->Write();
	    h_mMass2_nSigmaPion3_EP[i][j][k][l][m]->Write();
	  }
	}
      }
    }
  }
  // pt spectra
  for(Int_t i = 0; i < TriFlow::pt_total; i++) // pt_bin
  {
    for(Int_t j = TriFlow::Centrality_start; j < TriFlow::Centrality_stop; j++) // centrality bin
    {
      for(Int_t k = TriFlow::Charge_start; k < TriFlow::Charge_stop; k++) // charge bin
      {
	for(Int_t l = TriFlow::EtaGap_start; l < TriFlow::EtaGap_stop; l++) // eta gap bin
	{
	  h_pt_spectra2_pik[i][j][k][l]->Write();
	  h_pt_spectra3_pik[i][j][k][l]->Write();
	}
      }
    }
  }
  h_mToFYLocal_Mass2->Write();
  h_mToFZLocal_Mass2->Write();
  h_mMass2_pt->Write();
//  h_phi_psi2->Write();
//  h_phi_psi3->Write();
//  h_yield_phi_psi2->Write();
//  h_yield_phi_psi3->Write();
}

//-------------------------------------------------------------------------------------------

void StTriFlowHistoManger::WriteProton()
{
  // flow
  for(Int_t i = 0; i < TriFlow::pt_total; i++) // pt bin
  {
    for(Int_t j = TriFlow::Centrality_start; j < TriFlow::Centrality_stop; j++) // centrality bin
    {
      for(Int_t k = TriFlow::Charge_start; k < TriFlow::Charge_stop; k++) // charge bin
      {
	for(Int_t l = TriFlow::EtaGap_start; l < TriFlow::EtaGap_stop; l++) // eta gap bin
	{
	  for(Int_t m = 0; m < 7; m ++) // phi-psi bin
	  {
	    for(Int_t i_cut = 0; i_cut < 18; i_cut++)
	    {
	      TString KEY_Proton;

	      KEY_Proton = Form("pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_phi_Psi_%d_2nd_Proton_SysError_%d",i,j,k,l,m,i_cut);
	      h_mMass2_Proton2_EP[KEY_Proton]->Write();

	      KEY_Proton = Form("pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_phi_Psi_%d_3rd_Proton_SysError_%d",i,j,k,l,m,i_cut);
	      h_mMass2_Proton3_EP[KEY_Proton]->Write();
	    }
	  }
	}
      }
    }
  }

  // raw pt spectra
  for(Int_t i = 0; i < TriFlow::pt_total; i++) // pt bin
  {
    for(Int_t j = TriFlow::Centrality_start; j < TriFlow::Centrality_stop; j++) // centrality bin
    {
      for(Int_t k = TriFlow::Charge_start; k < TriFlow::Charge_stop; k++) // charge bin
      {
	for(Int_t l = TriFlow::EtaGap_start; l < TriFlow::EtaGap_stop; l++) // eta gap bin
	{
	  for(Int_t i_cut = 0; i_cut < 18; i_cut++)
	  {
	    TString KEY_Proton_Spec;
	    KEY_Proton_Spec = Form("Spectra_pt_%d_low_Centrality_%d_Charge_%d_EtaGap_%d_Proton_2nd_SysError_%d",i,j,k,l,i_cut);
	    h_pt_spectra2_proton[KEY_Proton_Spec]->Write();

	    KEY_Proton_Spec = Form("Spectra_pt_%d_high_Centrality_%d_Charge_%d_EtaGap_%d_Proton_2nd_SysError_%d",i,j,k,l,i_cut);
	    h_pt_spectra2_proton[KEY_Proton_Spec]->Write();

	    KEY_Proton_Spec = Form("Spectra_pt_%d_low_Centrality_%d_Charge_%d_EtaGap_%d_Proton_3rd_SysError_%d",i,j,k,l,i_cut);
	    h_pt_spectra3_proton[KEY_Proton_Spec]->Write();

	    KEY_Proton_Spec = Form("Spectra_pt_%d_high_Centrality_%d_Charge_%d_EtaGap_%d_Proton_3rd_SysError_%d",i,j,k,l,i_cut);
	    h_pt_spectra3_proton[KEY_Proton_Spec]->Write();
	  }
	}
      }
    }
  }
}

//-------------------------------------------------------------------------------------------

void StTriFlowHistoManger::WriteYileds_nSigPion()
{
  for(Int_t cent = 0; cent < 9; cent++)
  {
    for(Int_t k = TriFlow::Charge_start; k < TriFlow::Charge_stop; k++) // charge bin
    {
      for(Int_t l = TriFlow::EtaGap_start; l < TriFlow::EtaGap_stop; l++) // eta gap bin
      {
	h_mMass2_nSigmaPion_Yields_PiK_EP[cent][k][l]->Write();
      }
    }
  }
}

void StTriFlowHistoManger::WriteYileds_Proton()
{
  for(Int_t cent = 0; cent < 9; cent++)
  {
    for(Int_t k = TriFlow::Charge_start; k < TriFlow::Charge_stop; k++) // charge bin
    {
      for(Int_t l = TriFlow::EtaGap_start; l < TriFlow::EtaGap_stop; l++) // eta gap bin
      {
	for(Int_t i_cut = 0; i_cut < 18; i_cut++)
	{
	  TString KEY_Proton_Yield = Form("Centrality_%d_Charge_%d_EtaGap_%d_Yields_Proton_SysError_%d",j,k,l,i_cut);
	  h_mMass2_Yields_Proton_EP[KEY_Proton_Yield]->Write();
	}
      }
    }
  }
}
//-------------------------------------------------------------------------------------------

void StTriFlowHistoManger::WriteQA()
{
  for(Int_t l = TriFlow::EtaGap_start; l < TriFlow::EtaGap_stop; l++) // eta gap bin
  {
    h_mDEdx_pq_before[l]->Write();
    h_mMass2_pq_before[l]->Write();
    h_mDEdx_pq_after[l]->Write();
    h_mMass2_pq_after[l]->Write();
  }
  h_mRefMult->Write();
  h_mVz->Write();
}

void StTriFlowHistoManger::WriteQA_Detector()
{
  h_mDEdx->Write();
  h_mMass2->Write();
}
//-------------------------------------------------------------------------------------------

#include "StStrangenessHistoManger.h"
#include "StStrangenessCons.h"
#include "TString.h"
#include "TMath.h"
#include "TH1F.h"

ClassImp(StStrangenessHistoManger)

//-------------------------------------------------------------
StStrangenessHistoManger::StStrangenessHistoManger()
{
}

StStrangenessHistoManger::~StStrangenessHistoManger()
{
}
//-------------------------------------------------------------

void StStrangenessHistoManger::Init(Int_t X_flag, Int_t mode) // 0 for Same Event, 1 for Mixed Event
{
  // flow analysis
  for(Int_t i = 0; i < Strangeness::pt_total_phi; i++) // pt bin | TODO: increase pt_bin to 8 GeV/c
  {
    for(Int_t j = Strangeness::Centrality_start; j < Strangeness::Centrality_stop; j++) // centrality bin
    {
      for(Int_t l = Strangeness::EtaGap_start; l < Strangeness::EtaGap_stop; l++) // eta gap bin
      {
	for(Int_t m = 0; m < Strangeness::Phi_Psi_total; m++) // phi-psi bin
	{
	  TString Mode[2] = {"SE","ME"};
	  TString HistName;
	  HistName = Form("pt_%d_Centrality_%d_EtaGap_%d_phi_Psi_%d_2nd_%s_%s",i,j,l,m,Strangeness::Partype[mode].Data(),Mode[X_flag].Data());
	  h_mMass2_EP[i][j][l][m] = new TH1F(HistName.Data(),HistName.Data(),200,Strangeness::InvMass_low[mode],Strangeness::InvMass_high[mode]);
	  h_mMass2_EP[i][j][l][m]->Sumw2();
	  HistName = Form("pt_%d_Centrality_%d_EtaGap_%d_phi_Psi_%d_3rd_%s_%s",i,j,l,m,Strangeness::Partype[mode].Data(),Mode[X_flag].Data());
	  h_mMass3_EP[i][j][l][m] = new TH1F(HistName.Data(),HistName.Data(),200,Strangeness::InvMass_low[mode],Strangeness::InvMass_high[mode]);
	  h_mMass3_EP[i][j][l][m]->Sumw2();

	  // subtract K0s
	  HistName = Form("pt_%d_Centrality_%d_EtaGap_%d_phi_Psi_%d_2nd_%s_%s_sub",i,j,l,m,Strangeness::Partype[mode].Data(),Mode[X_flag].Data());
	  h_mMass2_EP_sub[i][j][l][m] = new TH1F(HistName.Data(),HistName.Data(),200,Strangeness::InvMass_low[mode],Strangeness::InvMass_high[mode]);
	  h_mMass2_EP_sub[i][j][l][m]->Sumw2();
	  HistName = Form("pt_%d_Centrality_%d_EtaGap_%d_phi_Psi_%d_3rd_%s_%s_sub",i,j,l,m,Strangeness::Partype[mode].Data(),Mode[X_flag].Data());
	  h_mMass3_EP_sub[i][j][l][m] = new TH1F(HistName.Data(),HistName.Data(),200,Strangeness::InvMass_low[mode],Strangeness::InvMass_high[mode]);
	  h_mMass3_EP_sub[i][j][l][m]->Sumw2();
	}
      }
    }
  }

  // raw pt spectra | TODO: use finer pt_bin
  for(Int_t i = 0; i < Strangeness::pt_total_phi; i++) // pt bin
  {
    for(Int_t j = Strangeness::Centrality_start; j < Strangeness::Centrality_stop; j++) // centrality bin
    {
      for(Int_t l = Strangeness::EtaGap_start; l < Strangeness::EtaGap_stop; l++) // eta gap bin
      {
	TString Mode[2] = {"SE","ME"};
	TString HistName; 
	HistName = Form("Spec_pt_%d_Centrality_%d_EtaGap_%d_%s_%s",i,j,l,Strangeness::Partype[mode].Data(),Mode[X_flag].Data());
	h_mMass_Spec[i][j][l] = new TH1F(HistName.Data(),HistName.Data(),200,Strangeness::InvMass_low[mode],Strangeness::InvMass_high[mode]);
	h_mMass_Spec[i][j][l]->Sumw2();

	// subtract K0s
	HistName = Form("Spec_pt_%d_Centrality_%d_EtaGap_%d_%s_%s_sub",i,j,l,Strangeness::Partype[mode].Data(),Mode[X_flag].Data());
	h_mMass_Spec_sub[i][j][l] = new TH1F(HistName.Data(),HistName.Data(),200,Strangeness::InvMass_low[mode],Strangeness::InvMass_high[mode]);
	h_mMass_Spec_sub[i][j][l]->Sumw2();
      }
    }
  }
  
  // Yields
  for(Int_t j = 0; j < 9; j++) // centrality bin
  {
    for(Int_t l = Strangeness::EtaGap_start; l < Strangeness::EtaGap_stop; l++) // eta gap bin
    {
      TString Mode[2] = {"SE","ME"};
      TString HistName;
      HistName = Form("Yields_Centrality_%d_EtaGap_%d_%s_%s",j,l,Strangeness::Partype[mode].Data(),Mode[X_flag].Data());
      h_mMass_Yields[j][l] = new TH1F(HistName.Data(),HistName.Data(),200,Strangeness::InvMass_low[mode],Strangeness::InvMass_high[mode]);
      h_mMass_Yields[j][l]->Sumw2();

      // subtract K0s
      HistName = Form("Yields_Centrality_%d_EtaGap_%d_%s_%s_sub",j,l,Strangeness::Partype[mode].Data(),Mode[X_flag].Data());
      h_mMass_Yields_sub[j][l] = new TH1F(HistName.Data(),HistName.Data(),200,Strangeness::InvMass_low[mode],Strangeness::InvMass_high[mode]);
      h_mMass_Yields_sub[j][l]->Sumw2();
    }
  }
}
//-------------------------------------------------------------
void StStrangenessHistoManger::Fill(Float_t pt, Int_t Cent9, Int_t eta_gap, Float_t phi_psi2, Float_t Res2, Float_t phi_psi3, Float_t Res3, Float_t InvMass, Double_t reweight)
{
  if(Res2 > 0.0)
  {
    for(Int_t i = 0; i < Strangeness::pt_total_phi; i++) // pt_bin
    {
      if(pt > Strangeness::pt_low_phi[i] && pt < Strangeness::pt_up_phi[i])
      {
	for(Int_t j = Strangeness::Centrality_start; j < Strangeness::Centrality_stop; j++) // centrality bin
	{
	  if(Cent9 >= Strangeness::cent_low[j] && Cent9 <= Strangeness::cent_up[j])
	  {
	    for(Int_t psi_bin = 0; psi_bin < 3; psi_bin++)
	    {
	      if(phi_psi2 >= Strangeness::Psi2_low[psi_bin] && phi_psi2 < Strangeness::Psi2_up[psi_bin])
	      {
		Double_t phi_psi2_final = phi_psi2 - (psi_bin-1)*2.0*TMath::Pi()/2.0;
		for(Int_t m = 0; m < Strangeness::Phi_Psi_total; m++) // phi-psi2 bin
		{
		  if(TMath::Abs(phi_psi2_final) >= Strangeness::phi_Psi2_low[m] && TMath::Abs(phi_psi2_final) < Strangeness::phi_Psi2_up[m])
		  {
		    // flow
		    h_mMass2_EP[i][j][eta_gap][m]->Fill(InvMass,(reweight/Res2));
		    // raw pt spectra
		    h_mMass_Spec[i][j][eta_gap]->Fill(InvMass,reweight);
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
    for(Int_t i = 0; i < Strangeness::pt_total_phi; i++) // pt_bin
    {
      if(pt > Strangeness::pt_low_phi[i] && pt < Strangeness::pt_up_phi[i])
      {
	for(Int_t j = Strangeness::Centrality_start; j < Strangeness::Centrality_stop; j++) // centrality bin
	{
	  if(Cent9 >= Strangeness::cent_low[j] && Cent9 <= Strangeness::cent_up[j])
	  {
	    for(Int_t psi_bin = 0; psi_bin < 5; psi_bin++)
	    {
	      if(phi_psi3 >= Strangeness::Psi3_low[psi_bin] && phi_psi3 < Strangeness::Psi3_up[psi_bin])
	      {
		Double_t phi_psi3_final = phi_psi3 - (psi_bin-2)*2.0*TMath::Pi()/3.0;
		for(Int_t m = 0; m < Strangeness::Phi_Psi_total; m++) // phi-psi3 bin
		{
		  if(TMath::Abs(phi_psi3_final) >= Strangeness::phi_Psi3_low[m] && TMath::Abs(phi_psi3_final) < Strangeness::phi_Psi3_up[m])
		  {
		    // flow
		    h_mMass3_EP[i][j][eta_gap][m]->Fill(InvMass,(reweight/Res3));
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
  h_mMass_Yields[Cent9][eta_gap]->Fill(InvMass,reweight);
}
void StStrangenessHistoManger::Fill_sub(Float_t pt, Int_t Cent9, Int_t eta_gap, Float_t phi_psi2, Float_t Res2, Float_t phi_psi3, Float_t Res3, Float_t InvMass, Double_t reweight)
{
  if(Res2 > 0.0)
  {
    for(Int_t i = 0; i < Strangeness::pt_total_phi; i++) // pt_bin
    {
      if(pt > Strangeness::pt_low_phi[i] && pt < Strangeness::pt_up_phi[i])
      {
	for(Int_t j = Strangeness::Centrality_start; j < Strangeness::Centrality_stop; j++) // centrality bin
	{
	  if(Cent9 >= Strangeness::cent_low[j] && Cent9 <= Strangeness::cent_up[j])
	  {
	    for(Int_t psi_bin = 0; psi_bin < 3; psi_bin++)
	    {
	      if(phi_psi2 >= Strangeness::Psi2_low[psi_bin] && phi_psi2 < Strangeness::Psi2_up[psi_bin])
	      {
		Double_t phi_psi2_final = phi_psi2 - (psi_bin-1)*2.0*TMath::Pi()/2.0;
		for(Int_t m = 0; m < Strangeness::Phi_Psi_total; m++) // phi-psi2 bin
		{
		  if(TMath::Abs(phi_psi2_final) >= Strangeness::phi_Psi2_low[m] && TMath::Abs(phi_psi2_final) < Strangeness::phi_Psi2_up[m])
		  {
		    // flow
		    h_mMass2_EP_sub[i][j][eta_gap][m]->Fill(InvMass,(reweight/Res2));
		    // raw pt spectra
		    h_mMass_Spec_sub[i][j][eta_gap]->Fill(InvMass,reweight);
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
    for(Int_t i = 0; i < Strangeness::pt_total_phi; i++) // pt_bin
    {
      if(pt > Strangeness::pt_low_phi[i] && pt < Strangeness::pt_up_phi[i])
      {
	for(Int_t j = Strangeness::Centrality_start; j < Strangeness::Centrality_stop; j++) // centrality bin
	{
	  if(Cent9 >= Strangeness::cent_low[j] && Cent9 <= Strangeness::cent_up[j])
	  {
	    for(Int_t psi_bin = 0; psi_bin < 5; psi_bin++)
	    {
	      if(phi_psi3 >= Strangeness::Psi3_low[psi_bin] && phi_psi3 < Strangeness::Psi3_up[psi_bin])
	      {
		Double_t phi_psi3_final = phi_psi3 - (psi_bin-2)*2.0*TMath::Pi()/3.0;
		for(Int_t m = 0; m < Strangeness::Phi_Psi_total; m++) // phi-psi3 bin
		{
		  if(TMath::Abs(phi_psi3_final) >= Strangeness::phi_Psi3_low[m] && TMath::Abs(phi_psi3_final) < Strangeness::phi_Psi3_up[m])
		  {
		    // flow
		    h_mMass3_EP_sub[i][j][eta_gap][m]->Fill(InvMass,(reweight/Res3));
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
  h_mMass_Yields_sub[Cent9][eta_gap]->Fill(InvMass,reweight);
}
//-------------------------------------------------------------
void StStrangenessHistoManger::Write()
{
  // flow
  for(Int_t i = 0; i < Strangeness::pt_total_phi; i++) // pt bin
  {
    for(Int_t j = Strangeness::Centrality_start; j < Strangeness::Centrality_stop; j++) // centrality bin
    {
      for(Int_t l = Strangeness::EtaGap_start; l < Strangeness::EtaGap_stop; l++) // eta gap bin
      {
	for(Int_t m = 0; m < Strangeness::Phi_Psi_total; m ++) // phi-psi bin
	{
	  h_mMass2_EP[i][j][l][m]->Write();
	  h_mMass3_EP[i][j][l][m]->Write();
	  h_mMass2_EP_sub[i][j][l][m]->Write();
	  h_mMass3_EP_sub[i][j][l][m]->Write();
	}
      }
    }
  }

  // Yields
  for(Int_t j = 0; j < 9; j++) // centrality bin
  {
    for(Int_t l = Strangeness::EtaGap_start; l < Strangeness::EtaGap_stop; l++) // eta gap bin
    {
      h_mMass_Yields[j][l]->Write();
      h_mMass_Yields_sub[j][l]->Write();
    }
  }

  // raw pt spectra
  for(Int_t i = 0; i < Strangeness::pt_total_phi; i++) // pt bin
  {
    for(Int_t j = Strangeness::Centrality_start; j < Strangeness::Centrality_stop; j++) // centrality bin
    {
      for(Int_t l = Strangeness::EtaGap_start; l < Strangeness::EtaGap_stop; l++) // eta gap bin
      {
	h_mMass_Spec[i][j][l]->Write();
	h_mMass_Spec_sub[i][j][l]->Write();
      }
    }
  }
}

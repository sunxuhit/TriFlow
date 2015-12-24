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
  for(Int_t i = 0; i < Strangeness::pt_total_phi; i++) // pt bin 
  {
    for(Int_t j = Strangeness::Centrality_start; j < Strangeness::Centrality_stop; j++) // centrality bin
    {
      for(Int_t l = Strangeness::EtaGap_start; l < Strangeness::EtaGap_stop; l++) // eta gap bin
      {
	for(Int_t m = 0; m < Strangeness::Phi_Psi_total; m++) // phi-psi bin
	{
	  for(Int_t i_cut = 0; i_cut < Strangeness::N_cuts; i_cut++) // SysErrors
	  {
	    TString Mode[2] = {"SE","ME"};
	    TString HistName;
	    HistName = Form("pt_%d_Centrality_%d_EtaGap_%d_phi_Psi_%d_2nd_%s_%s_SysErrors_%d",i,j,l,m,Strangeness::Partype[mode].Data(),Mode[X_flag].Data(),i_cut);
	    h_mMass2_EP[i][j][l][m][i_cut] = new TH1F(HistName.Data(),HistName.Data(),200,Strangeness::InvMass_low[mode],Strangeness::InvMass_high[mode]);
	    h_mMass2_EP[i][j][l][m][i_cut]->Sumw2();
	    HistName = Form("pt_%d_Centrality_%d_EtaGap_%d_phi_Psi_%d_3rd_%s_%s_SysErrors_%d",i,j,l,m,Strangeness::Partype[mode].Data(),Mode[X_flag].Data(),i_cut);
	    h_mMass3_EP[i][j][l][m][i_cut] = new TH1F(HistName.Data(),HistName.Data(),200,Strangeness::InvMass_low[mode],Strangeness::InvMass_high[mode]);
	    h_mMass3_EP[i][j][l][m][i_cut]->Sumw2();

	    // subtract K0s
	    HistName = Form("pt_%d_Centrality_%d_EtaGap_%d_phi_Psi_%d_2nd_%s_%s_SysErrors_%d_sub",i,j,l,m,Strangeness::Partype[mode].Data(),Mode[X_flag].Data(),i_cut);
	    h_mMass2_EP_sub[i][j][l][m][i_cut] = new TH1F(HistName.Data(),HistName.Data(),200,Strangeness::InvMass_low[mode],Strangeness::InvMass_high[mode]);
	    h_mMass2_EP_sub[i][j][l][m][i_cut]->Sumw2();
	    HistName = Form("pt_%d_Centrality_%d_EtaGap_%d_phi_Psi_%d_3rd_%s_%s_SysErrors_%d_sub",i,j,l,m,Strangeness::Partype[mode].Data(),Mode[X_flag].Data(),i_cut);
	    h_mMass3_EP_sub[i][j][l][m][i_cut] = new TH1F(HistName.Data(),HistName.Data(),200,Strangeness::InvMass_low[mode],Strangeness::InvMass_high[mode]);
	    h_mMass3_EP_sub[i][j][l][m][i_cut]->Sumw2();
	  }
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
	for(Int_t i_cut = 0; i_cut < Strangeness::N_cuts; i_cut++) // SysErrors
	{
	  for(Int_t i_pT = 0; i_pT < 2; i_pT++)
	  {
	    TString Mode[2] = {"SE","ME"};
	    TString pT[2] = {"low","high"};
	    TString HistName; 
	    HistName = Form("Spec_pt_%d_%s_Centrality_%d_EtaGap_%d_%s_%s_SysErrors_%d",i,pT[i_pT].Data(),j,l,Strangeness::Partype[mode].Data(),Mode[X_flag].Data(),i_cut);
	    h_mMass_Spec[i][i_pT][j][l][i_cut] = new TH1F(HistName.Data(),HistName.Data(),200,Strangeness::InvMass_low[mode],Strangeness::InvMass_high[mode]);
	    h_mMass_Spec[i][i_pT][j][l][i_cut]->Sumw2();

	    // subtract K0s
	    HistName = Form("Spec_pt_%d_%s_Centrality_%d_EtaGap_%d_%s_%s_SysErrors_%d_sub",i,pT[i_pT].Data(),j,l,Strangeness::Partype[mode].Data(),Mode[X_flag].Data(),i_cut);
	    h_mMass_Spec_sub[i][i_pT][j][l][i_cut] = new TH1F(HistName.Data(),HistName.Data(),200,Strangeness::InvMass_low[mode],Strangeness::InvMass_high[mode]);
	    h_mMass_Spec_sub[i][i_pT][j][l][i_cut]->Sumw2();
	  }
	}
      }
    }
  }
  
  // Yields
  for(Int_t j = 0; j < 9; j++) // centrality bin
  {
    for(Int_t l = Strangeness::EtaGap_start; l < Strangeness::EtaGap_stop; l++) // eta gap bin
    {
      for(Int_t i_cut = 0; i_cut < Strangeness::N_cuts; i_cut++)
      {
	TString Mode[2] = {"SE","ME"};
	TString HistName;
	HistName = Form("Yields_Centrality_%d_EtaGap_%d_%s_%s_SysErrors_%d",j,l,Strangeness::Partype[mode].Data(),Mode[X_flag].Data(),i_cut);
	h_mMass_Yields[j][l][i_cut] = new TH1F(HistName.Data(),HistName.Data(),200,Strangeness::InvMass_low[mode],Strangeness::InvMass_high[mode]);
	h_mMass_Yields[j][l][i_cut]->Sumw2();

	// subtract K0s
	HistName = Form("Yields_Centrality_%d_EtaGap_%d_%s_%s_SysErrors_%d_sub",j,l,Strangeness::Partype[mode].Data(),Mode[X_flag].Data(),i_cut);
	h_mMass_Yields_sub[j][l][i_cut] = new TH1F(HistName.Data(),HistName.Data(),200,Strangeness::InvMass_low[mode],Strangeness::InvMass_high[mode]);
	h_mMass_Yields_sub[j][l][i_cut]->Sumw2();
      }
    }
  }
}
//-------------------------------------------------------------
void StStrangenessHistoManger::Fill(Float_t pt, Int_t Cent9, Int_t eta_gap, Float_t phi_psi2, Float_t Res2, Float_t phi_psi3, Float_t Res3, Float_t InvMass, Double_t reweight, Int_t i_cut)
{
  if(Res2 > 0.0)
  {
    for(Int_t i = 0; i < Strangeness::pt_total_phi; i++) // pt_bin
    {
      if(pt >= Strangeness::pt_low_phi[i] && pt < Strangeness::pt_up_phi[i])
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
		    h_mMass2_EP[i][j][eta_gap][m][i_cut]->Fill(InvMass,(reweight/Res2));
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
      if(pt >= Strangeness::pt_low_phi[i] && pt < Strangeness::pt_up_phi[i])
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
		    h_mMass3_EP[i][j][eta_gap][m][i_cut]->Fill(InvMass,(reweight/Res3));
		    cout << "i_cut = " << i_cut << ", pT = " << pt << ", InvMass_lTrack = " << InvMass << endl;
		    // raw pt spectra
		    if(pt < 0.5*(Strangeness::pt_low_phi[i]+Strangeness::pt_up_phi[i])) 
		    {
		      h_mMass_Spec[i][0][j][eta_gap][i_cut]->Fill(InvMass,reweight);
		    }
		    else
		    {
		      h_mMass_Spec[i][1][j][eta_gap][i_cut]->Fill(InvMass,reweight);
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
  h_mMass_Yields[Cent9][eta_gap][i_cut]->Fill(InvMass,reweight);
}
void StStrangenessHistoManger::Fill_sub(Float_t pt, Int_t Cent9, Int_t eta_gap, Float_t phi_psi2, Float_t Res2, Float_t phi_psi3, Float_t Res3, Float_t InvMass, Double_t reweight, Int_t i_cut)
{
  if(Res2 > 0.0)
  {
    for(Int_t i = 0; i < Strangeness::pt_total_phi; i++) // pt_bin
    {
      if(pt >= Strangeness::pt_low_phi[i] && pt < Strangeness::pt_up_phi[i])
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
		    h_mMass2_EP_sub[i][j][eta_gap][m][i_cut]->Fill(InvMass,(reweight/Res2));
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
      if(pt >= Strangeness::pt_low_phi[i] && pt < Strangeness::pt_up_phi[i])
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
		    h_mMass3_EP_sub[i][j][eta_gap][m][i_cut]->Fill(InvMass,(reweight/Res3));
		    // raw pt spectra
		    if(pt < 0.5*(Strangeness::pt_low_phi[i]+Strangeness::pt_up_phi[i])) 
		    {
		      h_mMass_Spec_sub[i][0][j][eta_gap][i_cut]->Fill(InvMass,reweight);
		    }
		    else
		    {
		      h_mMass_Spec_sub[i][1][j][eta_gap][i_cut]->Fill(InvMass,reweight);
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
  h_mMass_Yields_sub[Cent9][eta_gap][i_cut]->Fill(InvMass,reweight);
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
	  for(Int_t i_cut = 0; i_cut < Strangeness::N_cuts; i_cut++) // SysErrors
	  {
	    h_mMass2_EP[i][j][l][m][i_cut]->Write();
	    h_mMass3_EP[i][j][l][m][i_cut]->Write();
	    h_mMass2_EP_sub[i][j][l][m][i_cut]->Write();
	    h_mMass3_EP_sub[i][j][l][m][i_cut]->Write();
	  }
	}
      }
    }
  }

  // Yields
  for(Int_t j = 0; j < 9; j++) // centrality bin
  {
    for(Int_t l = Strangeness::EtaGap_start; l < Strangeness::EtaGap_stop; l++) // eta gap bin
    {
      for(Int_t i_cut = 0; i_cut < Strangeness::N_cuts; i_cut++) // SysErrors
      {
	h_mMass_Yields[j][l][i_cut]->Write();
	h_mMass_Yields_sub[j][l][i_cut]->Write();
      }
    }
  }

  // raw pt spectra
  for(Int_t i = 0; i < Strangeness::pt_total_phi; i++) // pt bin
  {
    for(Int_t j = Strangeness::Centrality_start; j < Strangeness::Centrality_stop; j++) // centrality bin
    {
      for(Int_t l = Strangeness::EtaGap_start; l < Strangeness::EtaGap_stop; l++) // eta gap bin
      {
	for(Int_t i_cut = 0; i_cut < Strangeness::N_cuts; i_cut++) // SysErrors
	{
	  for(Int_t i_pT = 0; i_pT < 2; i_pT++)
	  {
	    h_mMass_Spec[i][i_pT][j][l][i_cut]->Write();
	    h_mMass_Spec_sub[i][i_pT][j][l][i_cut]->Write();
	  }
	}
      }
    }
  }
}

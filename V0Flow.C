#include "TString.h"
#include "TFile.h"
#include "TH1F.h"
#include <map>
#include "TCanvas.h"
#include <iostream>
#include "draw.h"
#include "TF1.h"
#include "functions.h"
#include <vector>
#include "TMath.h"
#include "TLegend.h"

// mEnergy: 0, 200 GeV | 1, 39 GeV
// mPID: 0, phi | 1, Lambda | 2, anti-Lambda | 3, K0S
static const TString Energy[2] = {"200GeV","39GeV"};
static const TString PID[4] = {"Phi","Lambda","antiLambda","K0S"};
static const TString Order[2] = {"2nd","3rd"};
static const Float_t Norm_Start[4] = {1.04,1.14,1.14,0.41};
static const Float_t Norm_Stop[4]  = {1.05,1.19,1.19,0.46};
static const Float_t BW_Start[4] = {0.994,1.0,1.0,1.0};
static const Float_t BW_Stop[4]  = {1.050,1.0,1.0,1.0};
static const Float_t InvMass[4] = {1.019,1.116,1.116,0.498};
static const Float_t Width[4]   = {0.00426,0.0016,0.0016,0.0016};
static const Float_t nSigGaus = 1.0;
static const Double_t PI_max[2] = {TMath::Pi()/2.0,TMath::Pi()/3.0};
static const Float_t nSigV0 = 2.0;

static const Int_t pt_total = 25; // pT loop
static const Int_t pt_start = 0;
static const Int_t pt_stop  = 25;

// pt bin
//                                       0 , 1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 ,10 ,11 ,12 ,13 ,14 ,15 ,16 ,17 ,18 ,19 ,20 ,21 ,22, 23, 24
static const Float_t pt_low_raw[25] = {0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.4,5.8,6.2,6.6,7.2};
static const Float_t pt_up_raw[25]  = {0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.4,5.8,6.2,6.6,7.2,8.0};

static const Int_t pt_rebin = 11;
static const Float_t pt_low[pt_rebin] = {0.2,0.8,1.0,1.2,1.6,2.0,2.4,2.8,3.4,4.2,5.8};
static const Float_t pt_up[pt_rebin]  = {0.8,1.0,1.2,1.6,2.0,2.4,2.8,3.4,4.2,5.8,8.0};
static const Int_t pt_rebin_start[pt_rebin] = {0,3,4,5,7, 9,11,13,15,17,21};
static const Int_t pt_rebin_stop[pt_rebin]  = {2,3,4,6,8,10,12,14,16,20,24};
static const Int_t pt_rebin_first = 0;
static const Int_t pt_rebin_last  = 11;

static const Int_t Cent_total = 4; // Centrality loop
static const Int_t Cent_start = 0;
static const Int_t Cent_stop  = 1;

static const Int_t Eta_total = 4; // Eta loop
static const Int_t Eta_start = 0;
static const Int_t Eta_stop  = 1;

static const Int_t phi_total = 7; // phi loop
static const Int_t phi_start = 0;
static const Int_t phi_stop  = 7;

static const Int_t Sys_total = 20; // Systematic loop
static const Int_t Sys_start = 10;
static const Int_t Sys_stop  = 15;

typedef std::map<TString,TH1F*> TH1FMap;
typedef std::map<TString,TF1*> TF1Map;
typedef std::map<TString,std::vector<Float_t>> vecFMap;

void V0Flow(Int_t mEnergy = 0, Int_t mPID = 0, Int_t mOrder = 1)
{
  TString InPutFile_SE = Form("./Data/AuAu%s/%s/Yields_SE_%s.root",Energy[mEnergy].Data(),PID[mPID].Data(),Energy[mEnergy].Data());
  TString InPutFile_ME = Form("./Data/AuAu%s/%s/Yields_ME_%s.root",Energy[mEnergy].Data(),PID[mPID].Data(),Energy[mEnergy].Data());
  TFile *File_SE = TFile::Open(InPutFile_SE.Data());
  TFile *File_ME = TFile::Open(InPutFile_ME.Data());

  // read in histogram for same event and mixed event
  // calculate SE - ME
  TH1FMap h_mMass_SE, h_mMass_ME, h_mMass_SM;
  for(Int_t i_pt = pt_start; i_pt < pt_stop; i_pt++) // pt loop
  {
    for(Int_t i_cent = Cent_start; i_cent < Cent_stop; i_cent++) // Centrality loop
    {
      for(Int_t i_eta = Eta_start; i_eta < Eta_stop; i_eta++) // EtaGap loop
      {
	for(Int_t i_phi = phi_start; i_phi < phi_stop; i_phi++) // phi loop
	{
	  for(Int_t i_sys = Sys_start; i_sys < Sys_stop; i_sys++) // Systematic loop
	  {
	    TString KEY_SE = Form("pt_%d_Centrality_%d_EtaGap_%d_phi_Psi_%d_%s_%s_SE_SysErrors_%d",i_pt,i_cent,i_eta,i_phi,Order[mOrder].Data(),PID[mPID].Data(),i_sys);
	    h_mMass_SE[KEY_SE.Data()] = (TH1F*)File_SE->Get(KEY_SE.Data())->Clone(); 
	    Int_t Norm_bin_start = h_mMass_SE[KEY_SE.Data()]->FindBin(Norm_Start[mPID]);
	    Int_t Norm_bin_stop  = h_mMass_SE[KEY_SE.Data()]->FindBin(Norm_Stop[mPID]);
	    Float_t Inte_SE = h_mMass_SE[KEY_SE.Data()]->Integral(Norm_bin_start,Norm_bin_stop);

	    TString KEY_ME = Form("pt_%d_Centrality_%d_EtaGap_%d_phi_Psi_%d_%s_%s_ME_SysErrors_%d",i_pt,i_cent,i_eta,i_phi,Order[mOrder].Data(),PID[mPID].Data(),i_sys);
	    h_mMass_ME[KEY_ME.Data()] = (TH1F*)File_ME->Get(KEY_ME.Data())->Clone(); 
	    Float_t Inte_ME = h_mMass_ME[KEY_ME.Data()]->Integral(Norm_bin_start,Norm_bin_stop);
	    h_mMass_ME[KEY_ME.Data()]->Scale(Inte_SE/Inte_ME);

	    TString KEY_SM = Form("pt_%d_Centrality_%d_EtaGap_%d_phi_Psi_%d_%s_%s_SM_SysErrors_%d",i_pt,i_cent,i_eta,i_phi,Order[mOrder].Data(),PID[mPID].Data(),i_sys);
	    h_mMass_SM[KEY_SM.Data()] = (TH1F*)h_mMass_SE[KEY_SE]->Clone();
	    h_mMass_SM[KEY_SM.Data()]->Add(h_mMass_ME[KEY_ME.Data()],-1.0);
	  }
	}
      }
    }
  }

  /*
  // QA Plots for SE vs. ME
  TString KEY = "pt_4_Centrality_0_EtaGap_0_phi_Psi_3_3rd_Phi_SE_SysErrors_10";
  h_mMass_SE[KEY.Data()]->Draw("PE");

  KEY = "pt_4_Centrality_0_EtaGap_0_phi_Psi_3_3rd_Phi_ME_SysErrors_10";
  h_mMass_ME[KEY.Data()]->SetLineColor(2);
  h_mMass_ME[KEY.Data()]->SetFillColor(2);
  h_mMass_ME[KEY.Data()]->SetFillStyle(3002);
  h_mMass_ME[KEY.Data()]->Draw("h same");

  KEY = "pt_4_Centrality_0_EtaGap_0_phi_Psi_3_3rd_Phi_SM_SysErrors_10";
  h_mMass_SM[KEY.Data()]->SetLineColor(4);
  h_mMass_SM[KEY.Data()]->SetFillColor(4);
  h_mMass_SM[KEY.Data()]->SetFillStyle(3004);
  h_mMass_SM[KEY.Data()]->Draw("h same");
  */


  /*
  // QA Plots for pT bins
  TCanvas *c_pT = new TCanvas("c_pT","c_pT",10,10,1400,1400);
  c_pT->Divide(5,5);
  for(Int_t i_pt = pt_start; i_pt < pt_stop; i_pt++) // pt loop
  {
    c_pT->cd(i_pt+1);
    c_pT->cd(i_pt+1)->SetLeftMargin(0.15);
    c_pT->cd(i_pt+1)->SetBottomMargin(0.15);
    c_pT->cd(i_pt+1)->SetTicks(1,1);
    c_pT->cd(i_pt+1)->SetGrid(0,0);
    TString KEY_SE = Form("pt_%d_Centrality_0_EtaGap_0_phi_Psi_2_%s_%s_SE_SysErrors_14",i_pt,Order[mOrder].Data(),PID[mPID].Data());
    h_mMass_SE[KEY_SE.Data()]->Draw();

    TString KEY_ME = Form("pt_%d_Centrality_0_EtaGap_0_phi_Psi_2_%s_%s_ME_SysErrors_14",i_pt,Order[mOrder].Data(),PID[mPID].Data());
    h_mMass_ME[KEY_ME.Data()]->SetLineColor(2);
    h_mMass_ME[KEY_ME.Data()]->SetFillColor(2);
    h_mMass_ME[KEY_ME.Data()]->SetFillStyle(3002);
    h_mMass_ME[KEY_ME.Data()]->Draw("h same");

    TString KEY_SM = Form("pt_%d_Centrality_0_EtaGap_0_phi_Psi_2_%s_%s_SM_SysErrors_14",i_pt,Order[mOrder].Data(),PID[mPID].Data());
    h_mMass_SM[KEY_SM.Data()]->SetLineColor(4);
    h_mMass_SM[KEY_SM.Data()]->SetFillColor(4);
    h_mMass_SM[KEY_SM.Data()]->SetFillStyle(3004);
    h_mMass_SM[KEY_SM.Data()]->Draw("h same");

    TString pT_range = Form("[%.2f,%.2f]",pt_low_raw[i_pt],pt_up_raw[i_pt]);
    plotTopLegend((char*)pT_range.Data(),0.2,0.7,0.08,1,0.0,42,1);
  }
  */

  // pT rebin
  TH1FMap h_mMass; // rebinned InvMass distribution, SE-ME

  for(Int_t i_cent = Cent_start; i_cent < Cent_stop; i_cent++) // Centrality loop
  {
    for(Int_t i_eta = Eta_start; i_eta < Eta_stop; i_eta++) // EtaGap loop
    {
      for(Int_t i_phi = phi_start; i_phi < phi_stop; i_phi++) // phi loop
      {
	for(Int_t i_sys = Sys_start; i_sys < Sys_stop; i_sys++) // Systematic loop
	{
	  for(Int_t pt_bin = pt_rebin_first; pt_bin < pt_rebin_last; pt_bin++) // pt loop
	  {
	    TString KEY = Form("pt_%d_Centrality_%d_EtaGap_%d_phi_Psi_%d_%s_%s_SM_SysErrors_%d",pt_bin,i_cent,i_eta,i_phi,Order[mOrder].Data(),PID[mPID].Data(),i_sys);
	    for(Int_t i_pt = pt_rebin_start[pt_bin]; i_pt <= pt_rebin_stop[pt_bin]; i_pt++)
	    {
	      TString KEY_SM = Form("pt_%d_Centrality_%d_EtaGap_%d_phi_Psi_%d_%s_%s_SM_SysErrors_%d",i_pt,i_cent,i_eta,i_phi,Order[mOrder].Data(),PID[mPID].Data(),i_sys);
//	      cout << "KEY= " << KEY.Data() << ", KEY_SM = " << KEY_SM.Data() << endl;
	      if(i_pt == pt_rebin_start[pt_bin])
	      {
		h_mMass[KEY] = (TH1F*)h_mMass_SM[KEY_SM]->Clone();
	      }
	      else
	      {
		h_mMass[KEY]->Add(h_mMass_SM[KEY_SM],1.0);
	      }
	    }
	  }
	}
      }
    }
  }
  
  /*
  // QA Plots for pT rebins
  TCanvas *c_pT_rebin = new TCanvas("c_pT_rebin","c_pT_rebin",10,10,1400,1400);
  c_pT_rebin->Divide(5,5);
  for(Int_t i_pt = pt_rebin_first; i_pt < pt_rebin_last; i_pt++) // pt loop
  {
    c_pT_rebin->cd(pt_rebin_start[i_pt]+1);
    c_pT_rebin->cd(pt_rebin_start[i_pt]+1)->SetLeftMargin(0.15);
    c_pT_rebin->cd(pt_rebin_start[i_pt]+1)->SetBottomMargin(0.15);
    c_pT_rebin->cd(pt_rebin_start[i_pt]+1)->SetTicks(1,1);
    c_pT_rebin->cd(pt_rebin_start[i_pt]+1)->SetGrid(0,0);
    for(Int_t i_sys = Sys_start; i_sys < Sys_stop; i_sys++) // Systematic loop
    {
      TString KEY = Form("pt_%d_Centrality_0_EtaGap_0_phi_Psi_2_%s_%s_SM_SysErrors_%d",i_pt,Order[mOrder].Data(),PID[mPID].Data(),i_sys);
      h_mMass[KEY]->SetMarkerColor(1+i_sys-Sys_start);
      h_mMass[KEY]->SetMarkerStyle(20);
      h_mMass[KEY]->SetMarkerSize(0.4);
      if(i_sys == Sys_start) h_mMass[KEY]->DrawCopy("pE");
      else h_mMass[KEY]->DrawCopy("pE same");
    }

    TString pT_range = Form("[%.2f,%.2f]",pt_low[i_pt],pt_up[i_pt]);
    plotTopLegend((char*)pT_range.Data(),0.2,0.7,0.08,1,0.0,42,1);
  }
  */

  if(mPID == 0) // Polynomial fit subtraction is only needed for phi meson
  {
    // Poly + Breit Wignar fit to phi integrated InvMass
    TH1FMap h_mMass_phi;
    vecFMap ParFit_phi;

    for(Int_t i_pt = pt_rebin_first; i_pt < pt_rebin_last; i_pt++) // pt loop
    {
      for(Int_t i_cent = Cent_start; i_cent < Cent_stop; i_cent++) // Centrality loop
      {
	for(Int_t i_eta = Eta_start; i_eta < Eta_stop; i_eta++) // EtaGap loop
	{
	  for(Int_t i_sys = Sys_start; i_sys < Sys_stop; i_sys++) // Systematic loop
	  {
	    TString KEY_phi = Form("pt_%d_Centrality_%d_EtaGap_%d_phi_Psi_%d_%s_%s_SM_SysErrors_%d",i_pt,i_cent,i_eta,phi_start,Order[mOrder].Data(),PID[mPID].Data(),i_sys);
	    for(Int_t i_phi = phi_start; i_phi < phi_stop; i_phi++) // phi loop
	    {
	      TString KEY = Form("pt_%d_Centrality_%d_EtaGap_%d_phi_Psi_%d_%s_%s_SM_SysErrors_%d",i_pt,i_cent,i_eta,i_phi,Order[mOrder].Data(),PID[mPID].Data(),i_sys);
	      if(i_phi == phi_start) h_mMass_phi[KEY_phi] = (TH1F*)h_mMass[KEY]->Clone();
	      else h_mMass_phi[KEY_phi]->Add(h_mMass[KEY],1.0);
	    }
	    TF1 *f_bw = new TF1("f_bw",PolyBreitWigner,BW_Start[mPID],BW_Stop[mPID],5);
	    for(Int_t i_par = 0; i_par < 5; i_par++)
	    {
	      f_bw->ReleaseParameter(i_par);
	    }
	    f_bw->SetParameter(0,1.019);
	    f_bw->SetParLimits(0,1.014,1.024);
	    f_bw->SetParameter(1,0.0055);
	    f_bw->SetParameter(2,10000);
	    f_bw->SetParameter(3,-6000);
	    f_bw->SetParameter(4,0.5);
	    f_bw->SetRange(BW_Start[mPID],BW_Stop[mPID]);
	    ParFit_phi[KEY_phi].clear();
	    h_mMass_phi[KEY_phi]->Fit(f_bw,"NQR");
	    for(Int_t n_par = 0; n_par < 5; n_par++)
	    {
	      ParFit_phi[KEY_phi].push_back(static_cast<Float_t>(f_bw->GetParameter(n_par)));
	    }
	  }
	}
      }
    }

    /*
    // QA plots for Poly+Breit_Wignar fits for phi integrated InvMass
    TCanvas *c_mMass_phi = new TCanvas("c_mMass_phi","c_mMass_phi",10,10,1400,1400);
    c_mMass_phi->Divide(5,5);
    for(Int_t i_pt = pt_rebin_first; i_pt < pt_rebin_last; i_pt++) // pt loop
    {
      c_mMass_phi->cd(pt_rebin_start[i_pt]+1);
      c_mMass_phi->cd(pt_rebin_start[i_pt]+1)->SetLeftMargin(0.15);
      c_mMass_phi->cd(pt_rebin_start[i_pt]+1)->SetBottomMargin(0.15);
      c_mMass_phi->cd(pt_rebin_start[i_pt]+1)->SetTicks(1,1);
      c_mMass_phi->cd(pt_rebin_start[i_pt]+1)->SetGrid(0,0);
      TString KEY_phi = Form("pt_%d_Centrality_0_EtaGap_0_phi_Psi_%d_%s_%s_SM_SysErrors_14",i_pt,phi_start,Order[mOrder].Data(),PID[mPID].Data());
      h_mMass_phi[KEY_phi]->SetMarkerColor(1);
      h_mMass_phi[KEY_phi]->SetMarkerStyle(24);
      h_mMass_phi[KEY_phi]->SetMarkerSize(0.8);
      h_mMass_phi[KEY_phi]->Draw("PE");

      TF1 *f_bw = new TF1("f_bw",PolyBreitWigner,BW_Start[mPID],BW_Stop[mPID],5);
      for(Int_t i_par = 0; i_par < 5; i_par++)
      {
	f_bw->SetParameter(i_par,ParFit_phi[KEY_phi][i_par]);
      }
      f_bw->SetLineColor(2);
      f_bw->SetLineStyle(1);
      f_bw->SetLineWidth(2);
      f_bw->DrawCopy("l same");

      TString FuncName_poly = "ploy" + KEY_phi;
      TF1 *f_poly = new TF1(FuncName_poly.Data(),Poly,BW_Start[mPID],BW_Stop[mPID],2);
      f_poly->SetParameter(0,ParFit_phi[KEY_phi][3]);
      f_poly->SetParameter(1,ParFit_phi[KEY_phi][4]);
      f_poly->SetLineColor(4);
      f_poly->SetLineStyle(2);
      f_poly->SetLineWidth(4);
      f_poly->DrawCopy("l same");

      TString pT_range = Form("[%.2f,%.2f]",pt_low[i_pt],pt_up[i_pt]);
      plotTopLegend((char*)pT_range.Data(),0.2,0.7,0.08,1,0.0,42,1);
      PlotLine(0.98,1.05,0.0,0.0,1,2,2);
    }
    */

    // Poly+bw fits for phi differential InvMass
    for(Int_t i_pt = pt_rebin_first; i_pt < pt_rebin_last; i_pt++) // pt loop
    {
      for(Int_t i_cent = Cent_start; i_cent < Cent_stop; i_cent++) // Centrality loop
      {
	for(Int_t i_eta = Eta_start; i_eta < Eta_stop; i_eta++) // EtaGap loop
	{
	  for(Int_t i_sys = Sys_start; i_sys < Sys_stop; i_sys++) // Systematic loop
	  {
	    TString KEY_phi = Form("pt_%d_Centrality_%d_EtaGap_%d_phi_Psi_%d_%s_%s_SM_SysErrors_%d",i_pt,i_cent,i_eta,phi_start,Order[mOrder].Data(),PID[mPID].Data(),i_sys);
	    for(Int_t i_phi = phi_start; i_phi < phi_stop; i_phi++) // phi loop
	    {
	      TString KEY = Form("pt_%d_Centrality_%d_EtaGap_%d_phi_Psi_%d_%s_%s_SM_SysErrors_%d",i_pt,i_cent,i_eta,i_phi,Order[mOrder].Data(),PID[mPID].Data(),i_sys);
	      TF1 *f_bw = new TF1("f_bw",PolyBreitWigner,BW_Start[mPID],BW_Stop[mPID],5); //Poly+bw fits
	      for(Int_t i_par = 0; i_par < 5; i_par++)
	      {
		f_bw->ReleaseParameter(i_par);
	      }
	      f_bw->FixParameter(0,ParFit_phi[KEY_phi][0]);
	      f_bw->FixParameter(1,ParFit_phi[KEY_phi][1]);
	      f_bw->SetParameter(2,ParFit_phi[KEY_phi][2]/7.0);
	      f_bw->SetParameter(3,ParFit_phi[KEY_phi][3]);
	      f_bw->SetParameter(4,ParFit_phi[KEY_phi][4]);
	      f_bw->SetRange(BW_Start[mPID],BW_Stop[mPID]);
	      h_mMass[KEY]->Fit(f_bw,"NQR");

	      TString FuncName_poly = "ploy" + KEY_phi;
	      TF1 *f_poly = new TF1(FuncName_poly.Data(),Poly,BW_Start[mPID],BW_Stop[mPID],2);
	      f_poly->FixParameter(0,f_bw->GetParameter(3));
	      f_poly->FixParameter(1,f_bw->GetParameter(4));

	      h_mMass[KEY]->Add(f_poly,-1.0); // subtract linear background for phi differential InvMass
	    }
	  }
	}
      }
    }

    /*
    // QA plots for phi differential InvMass after linear background subtraction
    TCanvas *c_mMass_phi_diff = new TCanvas("c_mMass_phi_diff","c_mMass_phi_diff",10,10,1400,1400);
    c_mMass_phi_diff->Divide(5,5);
    for(Int_t i_pt = pt_rebin_first; i_pt < pt_rebin_last; i_pt++) // pt loop
    {
      c_mMass_phi_diff->cd(pt_rebin_start[i_pt]+1);
      c_mMass_phi_diff->cd(pt_rebin_start[i_pt]+1)->SetLeftMargin(0.15);
      c_mMass_phi_diff->cd(pt_rebin_start[i_pt]+1)->SetBottomMargin(0.15);
      c_mMass_phi_diff->cd(pt_rebin_start[i_pt]+1)->SetTicks(1,1);
      c_mMass_phi_diff->cd(pt_rebin_start[i_pt]+1)->SetGrid(0,0);
      TString KEY_phi = Form("pt_%d_Centrality_0_EtaGap_0_phi_Psi_3_%s_%s_SM_SysErrors_14",i_pt,Order[mOrder].Data(),PID[mPID].Data());
      h_mMass[KEY_phi]->SetMarkerColor(1);
      h_mMass[KEY_phi]->SetMarkerStyle(24);
      h_mMass[KEY_phi]->SetMarkerSize(0.8);
      h_mMass[KEY_phi]->Draw("PE");

      TString pT_range = Form("[%.2f,%.2f]",pt_low[i_pt],pt_up[i_pt]);
      plotTopLegend((char*)pT_range.Data(),0.2,0.7,0.08,1,0.0,42,1);
      PlotLine(0.98,1.05,0.0,0.0,1,2,2);
    }
    */
  }

  TH1FMap h_mMass_gauss; // phi integrated InvMass after linear background subtraction for guassian and bw fits to extract yields 
  vecFMap ParGaus, ParBW;

  for(Int_t i_pt = pt_rebin_first; i_pt < pt_rebin_last; i_pt++) // pt loop
  {
    for(Int_t i_cent = Cent_start; i_cent < Cent_stop; i_cent++) // Centrality loop
    {
      for(Int_t i_eta = Eta_start; i_eta < Eta_stop; i_eta++) // EtaGap loop
      {
	for(Int_t i_sys = Sys_start; i_sys < Sys_stop; i_sys++) // Systematic loop
	{
	  TString KEY_phi = Form("pt_%d_Centrality_%d_EtaGap_%d_phi_Psi_%d_%s_%s_SM_SysErrors_%d",i_pt,i_cent,i_eta,phi_start,Order[mOrder].Data(),PID[mPID].Data(),i_sys);
	  for(Int_t i_phi = phi_start; i_phi < phi_stop; i_phi++) // phi loop
	  {
	    TString KEY = Form("pt_%d_Centrality_%d_EtaGap_%d_phi_Psi_%d_%s_%s_SM_SysErrors_%d",i_pt,i_cent,i_eta,i_phi,Order[mOrder].Data(),PID[mPID].Data(),i_sys);
	    if(i_phi == phi_start) h_mMass_gauss[KEY_phi] = (TH1F*)h_mMass[KEY]->Clone();
	    else h_mMass_gauss[KEY_phi]->Add(h_mMass[KEY],1.0);
	  }
	  TF1 *f_gauss = new TF1("f_gauss",Gaussion,BW_Start[mPID],BW_Stop[mPID],3);
	  f_gauss->SetParameter(0,InvMass[mPID]);
	  f_gauss->SetParameter(1,Width[mPID]);
	  f_gauss->SetParameter(2,1000);
	  f_gauss->SetRange(InvMass[mPID]-nSigGaus*Width[mPID],InvMass[mPID]+nSigGaus*Width[mPID]);
	  h_mMass_gauss[KEY_phi]->Fit(f_gauss,"MQNR");
	  ParGaus[KEY_phi].push_back(static_cast<Float_t>(f_gauss->GetParameter(0)));
	  ParGaus[KEY_phi].push_back(static_cast<Float_t>(f_gauss->GetParameter(1)));
	  ParGaus[KEY_phi].push_back(static_cast<Float_t>(f_gauss->GetParameter(2)));

	  TF1 *f_bw = new TF1("f_bw",BreitWigner,BW_Start[mPID],BW_Stop[mPID],3);
	  f_bw->SetParameter(0,InvMass[mPID]);
	  f_bw->SetParameter(1,Width[mPID]);
	  f_bw->SetParameter(2,1000);
	  f_bw->SetRange(BW_Start[mPID],BW_Stop[mPID]);
	  h_mMass_gauss[KEY_phi]->Fit(f_bw,"MQNR");
	  ParBW[KEY_phi].push_back(static_cast<Float_t>(f_bw->GetParameter(0)));
	  ParBW[KEY_phi].push_back(static_cast<Float_t>(f_bw->GetParameter(1)));
	  ParBW[KEY_phi].push_back(static_cast<Float_t>(f_bw->GetParameter(2)));
	}
      }
    }
  }

  /*
  // QA guassian and bw fits to phi integrated InvMass
  TCanvas *c_mMass_guass = new TCanvas("c_mMass_guass","c_mMass_guass",10,10,1400,1400);
  c_mMass_guass->Divide(5,5);
  for(Int_t i_pt = pt_rebin_first; i_pt < pt_rebin_last; i_pt++) // pt loop
  {
    c_mMass_guass->cd(pt_rebin_start[i_pt]+1);
    c_mMass_guass->cd(pt_rebin_start[i_pt]+1)->SetLeftMargin(0.15);
    c_mMass_guass->cd(pt_rebin_start[i_pt]+1)->SetBottomMargin(0.15);
    c_mMass_guass->cd(pt_rebin_start[i_pt]+1)->SetTicks(1,1);
    c_mMass_guass->cd(pt_rebin_start[i_pt]+1)->SetGrid(0,0);
    TString KEY_phi = Form("pt_%d_Centrality_0_EtaGap_0_phi_Psi_%d_%s_%s_SM_SysErrors_14",i_pt,phi_start,Order[mOrder].Data(),PID[mPID].Data());
    h_mMass_gauss[KEY_phi]->SetMarkerColor(1);
    h_mMass_gauss[KEY_phi]->SetMarkerStyle(24);
    h_mMass_gauss[KEY_phi]->SetMarkerSize(0.8);
    h_mMass_gauss[KEY_phi]->Draw("PE");
    TF1 *f_gauss = new TF1("f_gauss",Gaussion,BW_Start[mPID],BW_Stop[mPID],3);
    f_gauss->FixParameter(0,ParGaus[KEY_phi][0]);
    f_gauss->FixParameter(1,ParGaus[KEY_phi][1]);
    f_gauss->FixParameter(2,ParGaus[KEY_phi][2]);
    f_gauss->SetRange(InvMass[mPID]-nSigGaus*Width[mPID],InvMass[mPID]+nSigGaus*Width[mPID]);
    f_gauss->SetLineColor(2);
    f_gauss->Draw("l same");

    TString pT_range = Form("[%.2f,%.2f]",pt_low[i_pt],pt_up[i_pt]);
    plotTopLegend((char*)pT_range.Data(),0.2,0.7,0.08,1,0.0,42,1);
    PlotLine(0.98,1.05,0.0,0.0,1,2,2);
  }

  TCanvas *c_mMass_bw = new TCanvas("c_mMass_bw","c_mMass_bw",10,10,1400,1400);
  c_mMass_bw->Divide(5,5);
  for(Int_t i_pt = pt_rebin_first; i_pt < pt_rebin_last; i_pt++) // pt loop
  {
    c_mMass_bw->cd(pt_rebin_start[i_pt]+1);
    c_mMass_bw->cd(pt_rebin_start[i_pt]+1)->SetLeftMargin(0.15);
    c_mMass_bw->cd(pt_rebin_start[i_pt]+1)->SetBottomMargin(0.15);
    c_mMass_bw->cd(pt_rebin_start[i_pt]+1)->SetTicks(1,1);
    c_mMass_bw->cd(pt_rebin_start[i_pt]+1)->SetGrid(0,0);
    TString KEY_phi = Form("pt_%d_Centrality_0_EtaGap_0_phi_Psi_%d_%s_%s_SM_SysErrors_14",i_pt,phi_start,Order[mOrder].Data(),PID[mPID].Data());
    h_mMass_gauss[KEY_phi]->SetMarkerColor(1);
    h_mMass_gauss[KEY_phi]->SetMarkerStyle(24);
    h_mMass_gauss[KEY_phi]->SetMarkerSize(0.8);
    h_mMass_gauss[KEY_phi]->Draw("PE");
    TF1 *f_bw = new TF1("f_bw",BreitWigner,BW_Start[mPID],BW_Stop[mPID],3);
    f_bw->SetParameter(0,ParBW[KEY_phi][0]);
    f_bw->SetParameter(1,ParBW[KEY_phi][1]);
    f_bw->SetParameter(2,ParBW[KEY_phi][2]);
    f_bw->SetRange(BW_Start[mPID],BW_Stop[mPID]);
    f_bw->SetLineColor(2);
    f_bw->Draw("l same");

    TString pT_range = Form("[%.2f,%.2f]",pt_low[i_pt],pt_up[i_pt]);
    plotTopLegend((char*)pT_range.Data(),0.2,0.7,0.08,1,0.0,42,1);
    PlotLine(0.98,1.05,0.0,0.0,1,2,2);
  }
  */

  // calculate counts and errors for phi bin with gaussian fits and breit wigner fits
  TH1FMap h_mCounts;
  for(Int_t i_pt = pt_rebin_first; i_pt < pt_rebin_last; i_pt++) // pt loop
  {
    for(Int_t i_cent = Cent_start; i_cent < Cent_stop; i_cent++) // Centrality loop
    {
      for(Int_t i_eta = Eta_start; i_eta < Eta_stop; i_eta++) // EtaGap loop
      {
	for(Int_t i_sys = Sys_start; i_sys < Sys_stop; i_sys++) // Systematic loop
	{
	  TString KEY_Gaus = Form("Gaus_pt_%d_Centrality_%d_EtaGap_%d_%s_%s_SM_SysErrors_%d",i_pt,i_cent,i_eta,Order[mOrder].Data(),PID[mPID].Data(),i_sys); // gaussian fits
	  h_mCounts[KEY_Gaus] = new TH1F(KEY_Gaus.Data(),KEY_Gaus.Data(),7,0.0,PI_max[mOrder]);
	  TString KEY_BW = Form("BW_pt_%d_Centrality_%d_EtaGap_%d_%s_%s_SM_SysErrors_%d",i_pt,i_cent,i_eta,Order[mOrder].Data(),PID[mPID].Data(),i_sys); // breit wigner fits
	  h_mCounts[KEY_BW] = new TH1F(KEY_BW.Data(),KEY_BW.Data(),7,0.0,PI_max[mOrder]);
	  TString KEY_phi = Form("pt_%d_Centrality_%d_EtaGap_%d_phi_Psi_%d_%s_%s_SM_SysErrors_%d",i_pt,i_cent,i_eta,phi_start,Order[mOrder].Data(),PID[mPID].Data(),i_sys);
	  for(Int_t i_phi = phi_start; i_phi < phi_stop; i_phi++) // phi loop
	  {
	    // gaussian counting
	    TString KEY = Form("pt_%d_Centrality_%d_EtaGap_%d_phi_Psi_%d_%s_%s_SM_SysErrors_%d",i_pt,i_cent,i_eta,i_phi,Order[mOrder].Data(),PID[mPID].Data(),i_sys);
	    Float_t counts = 0.0;
	    Float_t errors = 0.0;
	    Float_t bin_center = PI_max[mOrder]/14.0+i_phi*PI_max[mOrder]/7.0;
	    Int_t bin_start = h_mMass[KEY]->FindBin(ParGaus[KEY_phi][0]-nSigV0*ParGaus[KEY_phi][1]);
	    Int_t bin_stop  = h_mMass[KEY]->FindBin(ParGaus[KEY_phi][0]+nSigV0*ParGaus[KEY_phi][1]);
	    for(Int_t i_bin = bin_start; i_bin <= bin_stop; i_bin++)
	    {
	      counts += h_mMass[KEY]->GetBinContent(i_bin);
	      errors += h_mMass[KEY]->GetBinError(i_bin)*h_mMass[KEY]->GetBinError(i_bin);
	    }
	    h_mCounts[KEY_Gaus]->SetBinContent(h_mCounts[KEY_Gaus]->FindBin(bin_center),counts);
	    h_mCounts[KEY_Gaus]->SetBinError(h_mCounts[KEY_Gaus]->FindBin(bin_center),TMath::Sqrt(errors));

	    // breit wigner integrating
	    TF1 *f_bw = new TF1("f_bw",BreitWigner,BW_Start[mPID],BW_Stop[mPID],3);;
	    f_bw->FixParameter(0,ParBW[KEY_phi][0]);
	    f_bw->FixParameter(1,ParBW[KEY_phi][1]);
	    f_bw->SetParameter(2,1000);
	    f_bw->SetRange(BW_Start[mPID],BW_Stop[mPID]);
	    h_mMass[KEY]->Fit(f_bw,"NMQR");
	    Float_t bin_width = h_mMass[KEY]->GetBinWidth(1);
	    Float_t Inte_start = ParGaus[KEY_phi][0]-nSigV0*ParGaus[KEY_phi][1]-0.5*bin_width;
	    Float_t Inte_stop  = ParGaus[KEY_phi][0]+nSigV0*ParGaus[KEY_phi][1]+0.5*bin_width;
	    Float_t counts_bw = f_bw->Integral(Inte_start,Inte_stop)/bin_width;
	    Float_t errors_bw = f_bw->IntegralError(Inte_start,Inte_stop)/bin_width;
	    h_mCounts[KEY_BW]->SetBinContent(h_mCounts[KEY_BW]->FindBin(bin_center),counts_bw);
	    h_mCounts[KEY_BW]->SetBinError(h_mCounts[KEY_BW]->FindBin(bin_center),errors_bw);
	  }
	  //TODO: add flow fit to extract raw flow
	}
      }
    }
  }

  /*
  // QA InvMass vs. phi for gaussian and breit wigner fits
  TString KEY_phi = Form("pt_7_Centrality_0_EtaGap_0_phi_Psi_%d_%s_%s_SM_SysErrors_14",phi_start,Order[mOrder].Data(),PID[mPID].Data());
  TCanvas *c_mMass_psi = new TCanvas("c_mMass_psi","c_mMass_psi",10,10,900,900);
  c_mMass_psi->Divide(3,3);
  for(Int_t i_phi = phi_start; i_phi < phi_stop; i_phi++)
  {
    c_mMass_psi->cd(i_phi+1)->SetLeftMargin(0.15);
    c_mMass_psi->cd(i_phi+1)->SetBottomMargin(0.15);
    c_mMass_psi->cd(i_phi+1)->SetTicks(1,1);
    c_mMass_psi->cd(i_phi+1)->SetGrid(0,0);
  }
  for(Int_t i_phi = phi_start; i_phi < phi_stop; i_phi++)
  {
    c_mMass_psi->cd(i_phi+1);
    TString KEY = Form("pt_7_Centrality_0_EtaGap_0_phi_Psi_%d_%s_%s_SM_SysErrors_14",i_phi,Order[mOrder].Data(),PID[mPID].Data());
    h_mMass[KEY]->SetMarkerColor(1);
    h_mMass[KEY]->SetMarkerStyle(24);
    h_mMass[KEY]->SetMarkerSize(0.8);
    h_mMass[KEY]->Draw("hE");

    TF1 *f_bw = new TF1("f_bw",BreitWigner,BW_Start[mPID],BW_Stop[mPID],3);;
    f_bw->FixParameter(0,ParBW[KEY_phi][0]);
    f_bw->FixParameter(1,ParBW[KEY_phi][1]);
    f_bw->SetParameter(2,1000);
    f_bw->SetRange(BW_Start[mPID],BW_Stop[mPID]);
    h_mMass[KEY]->Fit(f_bw,"NMQR");
    f_bw->SetLineColor(2);
    f_bw->Draw("l same");

    Float_t x1 = ParGaus[KEY_phi][0] - nSigV0*ParGaus[KEY_phi][1];
    Float_t x2 = ParGaus[KEY_phi][0] + nSigV0*ParGaus[KEY_phi][1];
    Float_t y = h_mMass[KEY]->GetBinContent(h_mMass[KEY]->FindBin(ParGaus[KEY_phi][0]));
    PlotLine(x1,x1,0,y,4,2,2);
    PlotLine(x2,x2,0,y,4,2,2);
    PlotLine(0.98,1.05,0.0,0.0,1,2,2);
  }

  c_mMass_psi->cd(8);
  h_mMass_gauss[KEY_phi]->SetMarkerColor(1);
  h_mMass_gauss[KEY_phi]->SetMarkerStyle(24);
  h_mMass_gauss[KEY_phi]->SetMarkerSize(0.8);
  h_mMass_gauss[KEY_phi]->Draw("pE");
  TF1 *f_gauss_temp = new TF1("f_gauss_temp",Gaussion,BW_Start[mPID],BW_Stop[mPID],3);
  f_gauss_temp->FixParameter(0,ParGaus[KEY_phi][0]);
  f_gauss_temp->FixParameter(1,ParGaus[KEY_phi][1]);
  f_gauss_temp->FixParameter(2,ParGaus[KEY_phi][2]);
  f_gauss_temp->SetRange(InvMass[mPID]-nSigGaus*Width[mPID],InvMass[mPID]+nSigGaus*Width[mPID]);
  f_gauss_temp->SetLineColor(4);
  f_gauss_temp->SetLineWidth(2);
  f_gauss_temp->Draw("l same");
  PlotLine(0.98,1.05,0.0,0.0,1,2,2);

  c_mMass_psi->cd(9);
  TString KEY_Gaus = Form("Gaus_pt_7_Centrality_0_EtaGap_0_%s_%s_SM_SysErrors_14",Order[mOrder].Data(),PID[mPID].Data());
  h_mCounts[KEY_Gaus]->SetLineColor(4);
  h_mCounts[KEY_Gaus]->SetMarkerColor(4);
  h_mCounts[KEY_Gaus]->SetMarkerStyle(24);
  h_mCounts[KEY_Gaus]->SetMarkerSize(0.8);
  Float_t Inte_gaus = h_mCounts[KEY_Gaus]->Integral();
  h_mCounts[KEY_Gaus]->Scale(1.0/Inte_gaus);
  h_mCounts[KEY_Gaus]->Draw("pE");

  TString KEY_BW = Form("BW_pt_7_Centrality_0_EtaGap_0_%s_%s_SM_SysErrors_14",Order[mOrder].Data(),PID[mPID].Data());
  h_mCounts[KEY_BW]->SetLineColor(2);
  h_mCounts[KEY_BW]->SetMarkerColor(2);
  h_mCounts[KEY_BW]->SetMarkerStyle(24);
  h_mCounts[KEY_BW]->SetMarkerSize(0.8);
  Float_t Inte_bw = h_mCounts[KEY_BW]->Integral();
  h_mCounts[KEY_BW]->Scale(1.0/Inte_bw);
  h_mCounts[KEY_BW]->Draw("pE same");

  TLegend *leg_temp = new TLegend(0.5,0.6,0.8,0.8);
  leg_temp->SetFillColor(10);
  leg_temp->SetBorderSize(0.0);
  leg_temp->AddEntry(h_mCounts[KEY_Gaus],"bin counting","p");
  leg_temp->AddEntry(h_mCounts[KEY_BW],"breit wigner","p");
  leg_temp->Draw("same");
  */
}

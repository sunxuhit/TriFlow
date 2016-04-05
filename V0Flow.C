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
#include "TGaxis.h"
#include "TProfile.h"
#include "TGraphAsymmErrors.h"

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
static const Double_t PI_max[2] = {TMath::Pi()/2.0,TMath::Pi()/3.0};
static const Float_t nSigV0 = 2.0;
static const Float_t Flow_Order[2] = {2.0,3.0};

static const Int_t pt_total = 25; // pT loop
static const Int_t pt_start = 0;
static const Int_t pt_stop  = 25;
static const Int_t pt_QA    = 4;

// pt bin
//                                       0 , 1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 ,10 ,11 ,12 ,13 ,14 ,15 ,16 ,17 ,18 ,19 ,20 ,21 ,22, 23, 24
static const Float_t pt_low_raw[25] = {0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.4,5.8,6.2,6.6,7.2};
static const Float_t pt_up_raw[25]  = {0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.4,5.8,6.2,6.6,7.2,8.0};

// phi pt rebin
static const Int_t pt_rebin = 11;
static const Float_t pt_low[pt_rebin] = {0.2,0.8,1.0,1.2,1.6,2.0,2.4,2.8,3.4,4.2,5.8};
static const Float_t pt_up[pt_rebin]  = {0.8,1.0,1.2,1.6,2.0,2.4,2.8,3.4,4.2,5.8,8.0};
static const Int_t pt_rebin_start[pt_rebin] = {0,3,4,5,7, 9,11,13,15,17,21};
static const Int_t pt_rebin_stop[pt_rebin]  = {2,3,4,6,8,10,12,14,16,20,24};
static const Int_t pt_rebin_first = 0;
static const Int_t pt_rebin_last  = 11;

/*
// phi pt rebin consistency check
static const Int_t pt_rebin = 8;
static Float_t pt_low[pt_rebin] = {0.2,0.8,1.0,1.2,1.6,2.0,2.4,2.8};
static Float_t pt_up[pt_rebin]  = {0.8,1.0,1.2,1.6,2.0,2.4,2.8,4.2};
static Int_t pt_rebin_start[pt_rebin] = {0,3,4,5,7, 9,11,13};
static Int_t pt_rebin_stop[pt_rebin]  = {2,3,4,6,8,10,12,15};
static const Int_t pt_rebin_first = 0;
static const Int_t pt_rebin_last  = 8;
*/

static const Int_t Cent_total = 4; // Centrality loop
static const Int_t Cent_start = 0;
static const Int_t Cent_stop  = 1;
static const Int_t cent_low[4] = {0,7,4,0}; // 0 = 0-80%, 1 = 0-10%, 2 = 10-40%, 3 = 40-80%
static const Int_t cent_up[4]  = {8,8,6,3}; // 0 = 0-80%, 1 = 0-10%, 2 = 10-40%, 3 = 40-80%

static const Int_t Eta_total = 4; // Eta loop
static const Int_t Eta_start = 0;
static const Int_t Eta_stop  = 1;

static const Int_t phi_total = 7; // phi loop
static const Int_t phi_start = 0;
static const Int_t phi_stop  = 7;

static const Int_t Sys_total = 20; // Systematic loop
static const Int_t Sys_start = 0;
static const Int_t Sys_stop  = 20;

typedef std::map<TString,TH1F*> TH1FMap;
typedef std::map<TString,TProfile*> TProMap;
typedef std::map<TString,TGraphAsymmErrors*> TGraMap;
typedef std::map<TString,std::vector<Float_t>> vecFMap;

void V0Flow(Int_t mEnergy = 0, Int_t mPID = 0, Int_t mOrder = 1)
{
  TGaxis::SetMaxDigits(4);

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
	    h_mMass_SE[KEY_SE] = (TH1F*)File_SE->Get(KEY_SE.Data())->Clone(); 
	    Int_t Norm_bin_start = h_mMass_SE[KEY_SE]->FindBin(Norm_Start[mPID]);
	    Int_t Norm_bin_stop  = h_mMass_SE[KEY_SE]->FindBin(Norm_Stop[mPID]);
	    Float_t Inte_SE = h_mMass_SE[KEY_SE]->Integral(Norm_bin_start,Norm_bin_stop);

	    TString KEY_ME = Form("pt_%d_Centrality_%d_EtaGap_%d_phi_Psi_%d_%s_%s_ME_SysErrors_%d",i_pt,i_cent,i_eta,i_phi,Order[mOrder].Data(),PID[mPID].Data(),i_sys);
	    h_mMass_ME[KEY_ME] = (TH1F*)File_ME->Get(KEY_ME.Data())->Clone(); 
	    Float_t Inte_ME = h_mMass_ME[KEY_ME]->Integral(Norm_bin_start,Norm_bin_stop);
	    h_mMass_ME[KEY_ME]->Scale(Inte_SE/Inte_ME);

	    TString KEY_SM = Form("pt_%d_Centrality_%d_EtaGap_%d_phi_Psi_%d_%s_%s_SM_SysErrors_%d",i_pt,i_cent,i_eta,i_phi,Order[mOrder].Data(),PID[mPID].Data(),i_sys);
	    h_mMass_SM[KEY_SM] = (TH1F*)h_mMass_SE[KEY_SE]->Clone();
	    h_mMass_SM[KEY_SM]->Add(h_mMass_ME[KEY_ME],-1.0);
	  }
	}
      }
    }
  }

  /*
  // QA Plots for SE vs. ME
  TString KEY = Form("pt_%d_Centrality_%d_EtaGap_%d_phi_Psi_3_3rd_%s_SE_SysErrors_%d",pt_QA,Cent_start,Eta_start,PID[mPID].Data(),Sys_start);
  h_mMass_SE[KEY]->DrawCopy("PE");

  KEY = Form("pt_%d_Centrality_%d_EtaGap_%d_phi_Psi_3_3rd_%s_ME_SysErrors_%d",pt_QA,Cent_start,Eta_start,PID[mPID].Data(),Sys_start);
  h_mMass_ME[KEY]->SetLineColor(2);
  h_mMass_ME[KEY]->SetFillColor(2);
  h_mMass_ME[KEY]->SetFillStyle(3002);
  h_mMass_ME[KEY]->DrawCopy("h same");

  KEY = Form("pt_%d_Centrality_%d_EtaGap_%d_phi_Psi_3_3rd_%s_SM_SysErrors_%d",pt_QA,Cent_start,Eta_start,PID[mPID].Data(),Sys_start);
  h_mMass_SM[KEY]->SetLineColor(4);
  h_mMass_SM[KEY]->SetFillColor(4);
  h_mMass_SM[KEY]->SetFillStyle(3004);
  h_mMass_SM[KEY]->DrawCopy("h same");
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
    TString KEY_SE = Form("pt_%d_Centrality_%d_EtaGap_%d_phi_Psi_%d_%s_%s_SE_SysErrors_%d",i_pt,Cent_start,Eta_start,phi_start,Order[mOrder].Data(),PID[mPID].Data(),Sys_start);
    h_mMass_SE[KEY_SE]->DrawCopy();

    TString KEY_ME = Form("pt_%d_Centrality_%d_EtaGap_%d_phi_Psi_%d_%s_%s_ME_SysErrors_%d",i_pt,Cent_start,Eta_start,phi_start,Order[mOrder].Data(),PID[mPID].Data(),Sys_start);
    h_mMass_ME[KEY_ME]->SetLineColor(2);
    h_mMass_ME[KEY_ME]->SetFillColor(2);
    h_mMass_ME[KEY_ME]->SetFillStyle(3002);
    h_mMass_ME[KEY_ME]->DrawCopy("h same");

    TString KEY_SM = Form("pt_%d_Centrality_%d_EtaGap_%d_phi_Psi_%d_%s_%s_SM_SysErrors_%d",i_pt,Cent_start,Eta_start,phi_start,Order[mOrder].Data(),PID[mPID].Data(),Sys_start);
    h_mMass_SM[KEY_SM]->SetLineColor(4);
    h_mMass_SM[KEY_SM]->SetFillColor(4);
    h_mMass_SM[KEY_SM]->SetFillStyle(3004);
    h_mMass_SM[KEY_SM]->DrawCopy("h same");

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
      TString KEY = Form("pt_%d_Centrality_%d_EtaGap_%d_phi_Psi_%d_%s_%s_SM_SysErrors_%d",i_pt,Cent_start,Eta_start,phi_start,Order[mOrder].Data(),PID[mPID].Data(),i_sys);
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
      TString KEY_phi = Form("pt_%d_Centrality_%d_EtaGap_%d_phi_Psi_%d_%s_%s_SM_SysErrors_%d",i_pt,Cent_start,Eta_start,phi_start,Order[mOrder].Data(),PID[mPID].Data(),Sys_start);
      h_mMass_phi[KEY_phi]->SetMarkerColor(1);
      h_mMass_phi[KEY_phi]->SetMarkerStyle(24);
      h_mMass_phi[KEY_phi]->SetMarkerSize(0.8);
      h_mMass_phi[KEY_phi]->DrawCopy("PE");

      TF1 *f_bw = new TF1("f_bw",PolyBreitWigner,BW_Start[mPID],BW_Stop[mPID],5);
      for(Int_t i_par = 0; i_par < 5; i_par++)
      {
	f_bw->SetParameter(i_par,ParFit_phi[KEY_phi][i_par]);
      }
      f_bw->SetLineColor(2);
      f_bw->SetLineStyle(1);
      f_bw->SetLineWidth(2);
      f_bw->DrawCopy("l same");

      TF1 *f_poly = new TF1("f_poly",Poly,BW_Start[mPID],BW_Stop[mPID],2);
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

	      TF1 *f_poly = new TF1("f_poly",Poly,BW_Start[mPID],BW_Stop[mPID],2);
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
      TString KEY_phi = Form("pt_%d_Centrality_%d_EtaGap_%d_phi_Psi_3_%s_%s_SM_SysErrors_%d",i_pt,Cent_start,Eta_start,Order[mOrder].Data(),PID[mPID].Data(),Sys_start);
      h_mMass[KEY_phi]->SetMarkerColor(1);
      h_mMass[KEY_phi]->SetMarkerStyle(24);
      h_mMass[KEY_phi]->SetMarkerSize(0.8);
      h_mMass[KEY_phi]->DrawCopy("PE");

      TString pT_range = Form("[%.2f,%.2f]",pt_low[i_pt],pt_up[i_pt]);
      plotTopLegend((char*)pT_range.Data(),0.2,0.7,0.08,1,0.0,42,1);
      PlotLine(0.98,1.05,0.0,0.0,1,2,2);
    }
    */
  }

  TH1FMap h_mMass_total; // phi integrated InvMass after linear background subtraction for bw fits to extract yields 
  vecFMap ParBW;

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
	    if(i_phi == phi_start) h_mMass_total[KEY_phi] = (TH1F*)h_mMass[KEY]->Clone();
	    else h_mMass_total[KEY_phi]->Add(h_mMass[KEY],1.0);
	  }
	  TF1 *f_bw = new TF1("f_bw",BreitWigner,BW_Start[mPID],BW_Stop[mPID],3);
	  f_bw->SetParameter(0,InvMass[mPID]);
	  f_bw->SetParLimits(0,InvMass[mPID]-0.001,InvMass[mPID]+0.001);
	  f_bw->SetParameter(1,Width[mPID]);
	  f_bw->SetParameter(2,1000);
	  f_bw->SetRange(BW_Start[mPID],BW_Stop[mPID]);
	  h_mMass_total[KEY_phi]->Fit(f_bw,"MQNR");
	  ParBW[KEY_phi].clear();
	  ParBW[KEY_phi].push_back(static_cast<Float_t>(f_bw->GetParameter(0)));
	  ParBW[KEY_phi].push_back(static_cast<Float_t>(f_bw->GetParameter(1)));
	  ParBW[KEY_phi].push_back(static_cast<Float_t>(f_bw->GetParameter(2)));
	}
      }
    }
  }

  /*
  // QA: bw fits to phi integrated InvMass
  TCanvas *c_mMass_bw = new TCanvas("c_mMass_bw","c_mMass_bw",10,10,1400,1400);
  c_mMass_bw->Divide(5,5);
  for(Int_t i_pt = pt_rebin_first; i_pt < pt_rebin_last; i_pt++) // pt loop
  {
    c_mMass_bw->cd(pt_rebin_start[i_pt]+1);
    c_mMass_bw->cd(pt_rebin_start[i_pt]+1)->SetLeftMargin(0.15);
    c_mMass_bw->cd(pt_rebin_start[i_pt]+1)->SetBottomMargin(0.15);
    c_mMass_bw->cd(pt_rebin_start[i_pt]+1)->SetTicks(1,1);
    c_mMass_bw->cd(pt_rebin_start[i_pt]+1)->SetGrid(0,0);
    TString KEY_phi = Form("pt_%d_Centrality_%d_EtaGap_%d_phi_Psi_%d_%s_%s_SM_SysErrors_%d",i_pt,Cent_start,Eta_start,phi_start,Order[mOrder].Data(),PID[mPID].Data(),Sys_start);
    h_mMass_total[KEY_phi]->SetMarkerColor(1);
    h_mMass_total[KEY_phi]->SetMarkerStyle(24);
    h_mMass_total[KEY_phi]->SetMarkerSize(0.8);
    h_mMass_total[KEY_phi]->DrawCopy("PE");
    TF1 *f_bw = new TF1("f_bw",BreitWigner,BW_Start[mPID],BW_Stop[mPID],3);
    f_bw->SetParameter(0,ParBW[KEY_phi][0]);
    f_bw->SetParameter(1,ParBW[KEY_phi][1]);
    f_bw->SetParameter(2,ParBW[KEY_phi][2]);
    f_bw->SetRange(BW_Start[mPID],BW_Stop[mPID]);
    f_bw->SetLineColor(2);
    f_bw->DrawCopy("l same");

    TString pT_range = Form("[%.2f,%.2f]",pt_low[i_pt],pt_up[i_pt]);
    plotTopLegend((char*)pT_range.Data(),0.2,0.7,0.08,1,0.0,42,1);
    PlotLine(0.98,1.05,0.0,0.0,1,2,2);
  }
  */

  // calculate counts and errors for phi bin with breit wigner fits with bin counting and integrating
  TH1FMap h_mCounts, h_mRawFlow;
  vecFMap ParFlow_Gaus, ParFlow_BW;
  for(Int_t i_cent = Cent_start; i_cent < Cent_stop; i_cent++) // Centrality loop
  {
    for(Int_t i_eta = Eta_start; i_eta < Eta_stop; i_eta++) // EtaGap loop
    {
      for(Int_t i_sys = Sys_start; i_sys < Sys_stop; i_sys++) // Systematic loop
      {
	TString KEY_RawFlow_Gaus = Form("RawFlow_Gaus_Centrality_%d_EtaGap_%d_%s_%s_SM_SysErrors_%d",i_cent,i_eta,Order[mOrder].Data(),PID[mPID].Data(),i_sys); // bin counting
	h_mRawFlow[KEY_RawFlow_Gaus] = new TH1F(KEY_RawFlow_Gaus.Data(),KEY_RawFlow_Gaus.Data(),100,-0.05,9.95);
	TString KEY_RawFlow_BW = Form("RawFlow_BW_Centrality_%d_EtaGap_%d_%s_%s_SM_SysErrors_%d",i_cent,i_eta,Order[mOrder].Data(),PID[mPID].Data(),i_sys); // breit wigner fits
	h_mRawFlow[KEY_RawFlow_BW] = new TH1F(KEY_RawFlow_BW.Data(),KEY_RawFlow_BW.Data(),100,-0.05,9.95);
	for(Int_t i_pt = pt_rebin_first; i_pt < pt_rebin_last; i_pt++) // pt loop
	{
	  TString KEY_Gaus = Form("Gaus_pt_%d_Centrality_%d_EtaGap_%d_%s_%s_SM_SysErrors_%d",i_pt,i_cent,i_eta,Order[mOrder].Data(),PID[mPID].Data(),i_sys); // bin counting
	  h_mCounts[KEY_Gaus] = new TH1F(KEY_Gaus.Data(),KEY_Gaus.Data(),7,0.0,PI_max[mOrder]);
	  TString KEY_BW = Form("BW_pt_%d_Centrality_%d_EtaGap_%d_%s_%s_SM_SysErrors_%d",i_pt,i_cent,i_eta,Order[mOrder].Data(),PID[mPID].Data(),i_sys); // breit wigner fits
	  h_mCounts[KEY_BW] = new TH1F(KEY_BW.Data(),KEY_BW.Data(),7,0.0,PI_max[mOrder]);
	  TString KEY_phi = Form("pt_%d_Centrality_%d_EtaGap_%d_phi_Psi_%d_%s_%s_SM_SysErrors_%d",i_pt,i_cent,i_eta,phi_start,Order[mOrder].Data(),PID[mPID].Data(),i_sys);
	  for(Int_t i_phi = phi_start; i_phi < phi_stop; i_phi++) // phi loop
	  {
	    // bin counting
	    TString KEY = Form("pt_%d_Centrality_%d_EtaGap_%d_phi_Psi_%d_%s_%s_SM_SysErrors_%d",i_pt,i_cent,i_eta,i_phi,Order[mOrder].Data(),PID[mPID].Data(),i_sys);
	    Float_t counts = 0.0;
	    Float_t errors = 0.0;
	    Float_t bin_center = PI_max[mOrder]/14.0+i_phi*PI_max[mOrder]/7.0;
	    Int_t bin_start = h_mMass[KEY]->FindBin(ParBW[KEY_phi][0]-nSigV0*ParBW[KEY_phi][1]);
	    Int_t bin_stop  = h_mMass[KEY]->FindBin(ParBW[KEY_phi][0]+nSigV0*ParBW[KEY_phi][1]);
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
	    f_bw->SetParameter(2,ParBW[KEY_phi][2]/7.0);
	    f_bw->SetRange(BW_Start[mPID],BW_Stop[mPID]);
	    h_mMass[KEY]->Fit(f_bw,"NMQR");
	    Float_t bin_width = h_mMass[KEY]->GetBinWidth(1);
	    Float_t Inte_start = ParBW[KEY_phi][0]-nSigV0*ParBW[KEY_phi][1]-0.5*bin_width;
	    Float_t Inte_stop  = ParBW[KEY_phi][0]+nSigV0*ParBW[KEY_phi][1]+0.5*bin_width;
	    Float_t counts_bw = f_bw->Integral(Inte_start,Inte_stop)/bin_width;
	    Float_t errors_bw = f_bw->IntegralError(Inte_start,Inte_stop)/bin_width;
	    h_mCounts[KEY_BW]->SetBinContent(h_mCounts[KEY_BW]->FindBin(bin_center),counts_bw);
	    h_mCounts[KEY_BW]->SetBinError(h_mCounts[KEY_BW]->FindBin(bin_center),errors_bw);
	  }
	  Float_t pt_mean = (pt_low[i_pt]+pt_up[i_pt])/2.0;

	  TF1 *f_phi_gaus = new TF1("f_phi_gaus",flow,0.0,PI_max[mOrder],3);
	  f_phi_gaus->SetParameter(0,0.2);
	  f_phi_gaus->SetParameter(1,0.2);
	  f_phi_gaus->FixParameter(2,Flow_Order[mOrder]);
	  h_mCounts[KEY_Gaus]->Fit(f_phi_gaus,"NQM");
	  ParFlow_Gaus[KEY_Gaus].clear();
	  ParFlow_Gaus[KEY_Gaus].push_back(static_cast<Float_t>(f_phi_gaus->GetParameter(0)));
	  ParFlow_Gaus[KEY_Gaus].push_back(static_cast<Float_t>(f_phi_gaus->GetParameter(1)));
	  ParFlow_Gaus[KEY_Gaus].push_back(static_cast<Float_t>(f_phi_gaus->GetParameter(2)));
	  h_mRawFlow[KEY_RawFlow_Gaus]->SetBinContent(h_mRawFlow[KEY_RawFlow_Gaus]->FindBin(pt_mean),f_phi_gaus->GetParameter(1));
	  h_mRawFlow[KEY_RawFlow_Gaus]->SetBinError(h_mRawFlow[KEY_RawFlow_Gaus]->FindBin(pt_mean),f_phi_gaus->GetParError(1));

	  TF1 *f_phi_bw = new TF1("f_phi_bw",flow,0.0,PI_max[mOrder],3);
	  f_phi_bw->SetParameter(0,2.0);
	  f_phi_bw->SetParameter(1,1.0);
	  f_phi_bw->FixParameter(2,Flow_Order[mOrder]);
	  h_mCounts[KEY_BW]->Fit(f_phi_bw,"NQM");
	  ParFlow_BW[KEY_BW].clear();
	  ParFlow_BW[KEY_BW].push_back(static_cast<Float_t>(f_phi_bw->GetParameter(0)));
	  ParFlow_BW[KEY_BW].push_back(static_cast<Float_t>(f_phi_bw->GetParameter(1)));
	  ParFlow_BW[KEY_BW].push_back(static_cast<Float_t>(f_phi_bw->GetParameter(2)));
	  h_mRawFlow[KEY_RawFlow_BW]->SetBinContent(h_mRawFlow[KEY_RawFlow_BW]->FindBin(pt_mean),f_phi_bw->GetParameter(1));
	  h_mRawFlow[KEY_RawFlow_BW]->SetBinError(h_mRawFlow[KEY_RawFlow_BW]->FindBin(pt_mean),f_phi_bw->GetParError(1));
	}
      }
    }
  }

  /*
  // QA InvMass vs. phi for gaussian and breit wigner fits
  TString KEY_phi = Form("pt_%d_Centrality_%d_EtaGap_%d_phi_Psi_%d_%s_%s_SM_SysErrors_%d",pt_QA,Cent_start,Eta_start,phi_start,Order[mOrder].Data(),PID[mPID].Data(),Sys_start);
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
    TString KEY = Form("pt_%d_Centrality_%d_EtaGap_%d_phi_Psi_%d_%s_%s_SM_SysErrors_%d",pt_QA,Cent_start,Eta_start,i_phi,Order[mOrder].Data(),PID[mPID].Data(),Sys_start);
    h_mMass[KEY]->SetMarkerColor(1);
    h_mMass[KEY]->SetMarkerStyle(24);
    h_mMass[KEY]->SetMarkerSize(0.8);
    h_mMass[KEY]->DrawCopy("hE");

    TF1 *f_bw = new TF1("f_bw",BreitWigner,BW_Start[mPID],BW_Stop[mPID],3);;
    f_bw->FixParameter(0,ParBW[KEY_phi][0]);
    f_bw->FixParameter(1,ParBW[KEY_phi][1]);
    f_bw->SetParameter(2,ParBW[KEY_phi][2]/7.0);
    f_bw->SetRange(BW_Start[mPID],BW_Stop[mPID]);
    h_mMass[KEY]->Fit(f_bw,"NMQR");
    f_bw->SetLineColor(2);
    f_bw->DrawCopy("l same");

    Float_t x1 = ParBW[KEY_phi][0] - nSigV0*ParBW[KEY_phi][1];
    Float_t x2 = ParBW[KEY_phi][0] + nSigV0*ParBW[KEY_phi][1];
    Float_t y = h_mMass[KEY]->GetBinContent(h_mMass[KEY]->FindBin(ParBW[KEY_phi][0]));
    PlotLine(x1,x1,0,y,4,2,2);
    PlotLine(x2,x2,0,y,4,2,2);
    PlotLine(0.98,1.05,0.0,0.0,1,2,2);
  }

  c_mMass_psi->cd(8);
  h_mMass_total[KEY_phi]->SetMarkerColor(1);
  h_mMass_total[KEY_phi]->SetMarkerStyle(24);
  h_mMass_total[KEY_phi]->SetMarkerSize(0.8);
  h_mMass_total[KEY_phi]->DrawCopy("pE");
  TF1 *f_bw_QA = new TF1("f_bw_QA",BreitWigner,BW_Start[mPID],BW_Stop[mPID],3);;
  f_bw_QA->FixParameter(0,ParBW[KEY_phi][0]);
  f_bw_QA->FixParameter(1,ParBW[KEY_phi][1]);
  f_bw_QA->FixParameter(2,ParBW[KEY_phi][2]);
  f_bw_QA->SetRange(BW_Start[mPID],BW_Stop[mPID]);
  f_bw_QA->SetLineColor(2);
  f_bw_QA->Draw("l same");
  PlotLine(0.98,1.05,0.0,0.0,1,2,2);

  c_mMass_psi->cd(9);
  TString KEY_Gaus = Form("Gaus_pt_%d_Centrality_%d_EtaGap_%d_%s_%s_SM_SysErrors_%d",pt_QA,Cent_start,Eta_start,Order[mOrder].Data(),PID[mPID].Data(),Sys_start);
  h_mCounts[KEY_Gaus]->SetLineColor(4);
  h_mCounts[KEY_Gaus]->SetMarkerColor(4);
  h_mCounts[KEY_Gaus]->SetMarkerStyle(24);
  h_mCounts[KEY_Gaus]->SetMarkerSize(0.8);
  Float_t Inte_gaus = h_mCounts[KEY_Gaus]->Integral();
  h_mCounts[KEY_Gaus]->Scale(1.0/Inte_gaus);
  h_mCounts[KEY_Gaus]->DrawCopy("pE");
  TF1 *f_phi_gaus = new TF1("f_phi_gaus",flow,0.0,PI_max[mOrder],3);
  f_phi_gaus->FixParameter(0,ParFlow_Gaus[KEY_Gaus][0]/Inte_gaus);
//  f_phi_gaus->FixParameter(0,ParFlow_Gaus[KEY_Gaus][0]);
  f_phi_gaus->FixParameter(1,ParFlow_Gaus[KEY_Gaus][1]);
  f_phi_gaus->FixParameter(2,ParFlow_Gaus[KEY_Gaus][2]);
  f_phi_gaus->SetLineColor(4);
  f_phi_gaus->DrawCopy("l same");

  TString KEY_BW = Form("BW_pt_%d_Centrality_%d_EtaGap_%d_%s_%s_SM_SysErrors_%d",pt_QA,Cent_start,Eta_start,Order[mOrder].Data(),PID[mPID].Data(),Sys_start);
  h_mCounts[KEY_BW]->SetLineColor(2);
  h_mCounts[KEY_BW]->SetMarkerColor(2);
  h_mCounts[KEY_BW]->SetMarkerStyle(24);
  h_mCounts[KEY_BW]->SetMarkerSize(0.8);
  Float_t Inte_bw = h_mCounts[KEY_BW]->Integral();
  h_mCounts[KEY_BW]->Scale(1.0/Inte_bw);
  h_mCounts[KEY_BW]->DrawCopy("pE same");
  TF1 *f_phi_bw = new TF1("f_phi_bw",flow,0.0,PI_max[mOrder],3);
  f_phi_bw->FixParameter(0,ParFlow_BW[KEY_BW][0]/Inte_bw);
//  f_phi_bw->FixParameter(0,ParFlow_BW[KEY_BW][0]);
  f_phi_bw->FixParameter(1,ParFlow_BW[KEY_BW][1]);
  f_phi_bw->FixParameter(2,ParFlow_BW[KEY_BW][2]);
  f_phi_bw->SetLineColor(2);
  f_phi_bw->DrawCopy("l same");

  TLegend *leg_temp = new TLegend(0.5,0.6,0.8,0.8);
  leg_temp->SetFillColor(10);
  leg_temp->SetBorderSize(0.0);
  leg_temp->AddEntry(h_mCounts[KEY_Gaus],"bin counting","p");
  leg_temp->AddEntry(h_mCounts[KEY_BW],"breit wigner","p");
  leg_temp->Draw("same");
  */

  /*
  // QA raw flow vs. pt for gaussian fits and breit wigner fits
  TCanvas *c_raw_flow = new TCanvas("c_raw_flow","c_raw_flow",10,10,800,800);
  c_raw_flow->SetLeftMargin(0.15);
  c_raw_flow->SetBottomMargin(0.15);
  c_raw_flow->SetTicks(1,1);
  c_raw_flow->SetGrid(0,0);
  TH1F *h_play_raw = new TH1F("h_play_raw","h_play_raw",100,0.0,10.0);
  for(Int_t i_bin = 0; i_bin < 100; i_bin++)
  {
    h_play_raw->SetBinContent(i_bin+1,-10.0);
    h_play_raw->SetBinError(i_bin+1,1.0);
  }
  h_play_raw->SetTitle("");
  h_play_raw->SetStats(0);
  h_play_raw->GetXaxis()->SetRangeUser(0.0,8.0);
  h_play_raw->GetXaxis()->SetNdivisions(505,'N');
  h_play_raw->GetXaxis()->SetLabelSize(0.03);
  h_play_raw->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_play_raw->GetXaxis()->SetTitleSize(0.05);
  h_play_raw->GetXaxis()->SetTitleOffset(1.2);
  h_play_raw->GetXaxis()->CenterTitle();

  h_play_raw->GetYaxis()->SetRangeUser(-0.01,0.04);
  h_play_raw->GetYaxis()->SetNdivisions(505,'N');
  h_play_raw->GetYaxis()->SetTitle("v_{3}^{raw}");
  h_play_raw->GetYaxis()->SetTitleSize(0.05);
  h_play_raw->GetYaxis()->SetLabelSize(0.03);
  h_play_raw->GetYaxis()->CenterTitle();
  h_play_raw->DrawCopy("pE");
  PlotLine(0.0,8.0,0.0,0.0,1,2,2);
  for(Int_t i_cent = Cent_start; i_cent < Cent_stop; i_cent++) // Centrality loop
  {
    for(Int_t i_eta = Eta_start; i_eta < Eta_stop; i_eta++) // EtaGap loop
    {
      for(Int_t i_sys = Sys_start; i_sys < Sys_stop; i_sys++) // Systematic loop
      {
	TString KEY_RawFlow_Gaus = Form("RawFlow_Gaus_Centrality_%d_EtaGap_%d_%s_%s_SM_SysErrors_%d",i_cent,i_eta,Order[mOrder].Data(),PID[mPID].Data(),i_sys); // gaussian fits
	h_mRawFlow[KEY_RawFlow_Gaus]->SetMarkerStyle(24);
	h_mRawFlow[KEY_RawFlow_Gaus]->SetMarkerColor(4);
	h_mRawFlow[KEY_RawFlow_Gaus]->DrawCopy("PE same");
	TString KEY_RawFlow_BW = Form("RawFlow_BW_Centrality_%d_EtaGap_%d_%s_%s_SM_SysErrors_%d",i_cent,i_eta,Order[mOrder].Data(),PID[mPID].Data(),i_sys); // breit wigner fits
	h_mRawFlow[KEY_RawFlow_BW]->SetMarkerStyle(24);
	h_mRawFlow[KEY_RawFlow_BW]->SetMarkerColor(2);
	h_mRawFlow[KEY_RawFlow_BW]->DrawCopy("PE same");
      }
    }
  }
  TString KEY_RawFlow_Gaus_QA = Form("RawFlow_Gaus_Centrality_0_EtaGap_0_%s_%s_SM_SysErrors_%d",Order[mOrder].Data(),PID[mPID].Data(),Sys_start); // gaussian fits
  TString KEY_RawFlow_BW_QA = Form("RawFlow_BW_Centrality_0_EtaGap_0_%s_%s_SM_SysErrors_%d",Order[mOrder].Data(),PID[mPID].Data(),Sys_start); // gaussian fits
  TLegend *leg_raw_flow = new TLegend(0.5,0.6,0.8,0.8);
  leg_raw_flow->SetFillColor(10);
  leg_raw_flow->SetBorderSize(0.0);
  leg_raw_flow->AddEntry(h_mRawFlow[KEY_RawFlow_Gaus_QA],"bin counting","p");
  leg_raw_flow->AddEntry(h_mRawFlow[KEY_RawFlow_BW_QA],"breit wigner","p");
  leg_raw_flow->Draw("same");
//  c_raw_flow->SaveAs("./raw_flow.eps");
  */

  // resolution correction
  // extract yields
  TH1FMap h_mYield_SE, h_mYield_ME, h_mYield;
  for(Int_t i_cent = 0; i_cent < 9; i_cent++) // Centrality loop
  {
    for(Int_t i_eta = Eta_start; i_eta < Eta_stop; i_eta++) // EtaGap loop
    {
      for(Int_t i_sys = Sys_start; i_sys < Sys_stop; i_sys++) // Systematic loop
      {
	TString KEY_Yield_SE = Form("Yields_Centrality_%d_EtaGap_%d_%s_SE_SysErrors_%d",i_cent,i_eta,PID[mPID].Data(),i_sys);
	h_mYield_SE[KEY_Yield_SE] = (TH1F*)File_SE->Get(KEY_Yield_SE.Data())->Clone(); 
	Int_t Norm_bin_start = h_mYield_SE[KEY_Yield_SE]->FindBin(Norm_Start[mPID]);
	Int_t Norm_bin_stop  = h_mYield_SE[KEY_Yield_SE]->FindBin(Norm_Stop[mPID]);
	Float_t Inte_SE = h_mYield_SE[KEY_Yield_SE]->Integral(Norm_bin_start,Norm_bin_stop);

	TString KEY_Yield_ME = Form("Yields_Centrality_%d_EtaGap_%d_%s_ME_SysErrors_%d",i_cent,i_eta,PID[mPID].Data(),i_sys);
	h_mYield_ME[KEY_Yield_ME] = (TH1F*)File_ME->Get(KEY_Yield_ME.Data())->Clone(); 
	Float_t Inte_ME = h_mYield_ME[KEY_Yield_ME]->Integral(Norm_bin_start,Norm_bin_stop);
	h_mYield_ME[KEY_Yield_ME]->Scale(Inte_SE/Inte_ME);

	TString KEY_Yield = Form("Yields_Centrality_%d_EtaGap_%d_%s_SysErrors_%d",i_cent,i_eta,PID[mPID].Data(),i_sys);
	h_mYield[KEY_Yield] = (TH1F*)h_mYield_SE[KEY_Yield_SE]->Clone();
	h_mYield[KEY_Yield]->Add(h_mYield_ME[KEY_Yield_ME],-1.0);
      }
    }
  }

  /*
  //QA Yields vs Centrality
  TCanvas *c_Yields = new TCanvas("c_Yields","c_Yields",10,10,900,900);
  c_Yields->Divide(3,3);
  for(Int_t i_cent = 0; i_cent < 9; i_cent++)
  {
    c_Yields->cd(i_cent+1);
    c_Yields->cd(i_cent+1)->SetLeftMargin(0.20);
    c_Yields->cd(i_cent+1)->SetBottomMargin(0.20);
    c_Yields->cd(i_cent+1)->SetTicks(1,1);
    c_Yields->cd(i_cent+1)->SetGrid(0,0);
    TString KEY_Yield_SE = Form("Yields_Centrality_%d_EtaGap_%d_%s_SE_SysErrors_%d",i_cent,Eta_start,PID[mPID].Data(),Sys_start);
    h_mYield_SE[KEY_Yield_SE]->SetTitle("");
    h_mYield_SE[KEY_Yield_SE]->SetStats(0);
    h_mYield_SE[KEY_Yield_SE]->SetMarkerStyle(24);
    h_mYield_SE[KEY_Yield_SE]->SetMarkerSize(1.0);
    h_mYield_SE[KEY_Yield_SE]->GetXaxis()->SetNdivisions(505,'N');
    h_mYield_SE[KEY_Yield_SE]->GetXaxis()->SetTitle("InvMass (GeV/c^{2})");
    h_mYield_SE[KEY_Yield_SE]->GetXaxis()->CenterTitle();
    h_mYield_SE[KEY_Yield_SE]->GetYaxis()->SetTitle("Counts");
    h_mYield_SE[KEY_Yield_SE]->GetYaxis()->CenterTitle();
    h_mYield_SE[KEY_Yield_SE]->GetYaxis()->SetTitleOffset(1.2);
    h_mYield_SE[KEY_Yield_SE]->DrawCopy("pE");

    TString KEY_Yield_ME = Form("Yields_Centrality_%d_EtaGap_%d_%s_ME_SysErrors_%d",i_cent,Eta_start,PID[mPID].Data(),Sys_start);
    h_mYield_ME[KEY_Yield_ME]->SetFillColor(2);
    h_mYield_ME[KEY_Yield_ME]->SetFillStyle(3002);
    h_mYield_ME[KEY_Yield_ME]->DrawCopy("h Same");

    TString KEY_Yield = Form("Yields_Centrality_%d_EtaGap_%d_%s_SysErrors_%d",i_cent,Eta_start,PID[mPID].Data(),Sys_start);
    h_mYield[KEY_Yield]->SetFillColor(4);
    h_mYield[KEY_Yield]->SetFillStyle(3001);
    h_mYield[KEY_Yield]->DrawCopy("h Same");
  }
  */

  if(mPID == 0) // Polynomial fit subtraction is only needed for phi meson
  {
    // Poly + Breit Wignar fit to Yields
    TH1FMap h_mYield_SM; // for QA plot only
    vecFMap ParYield_SM;

    for(Int_t i_cent = 0; i_cent < 9; i_cent++) // Centrality loop
    {
      for(Int_t i_eta = Eta_start; i_eta < Eta_stop; i_eta++) // EtaGap loop
      {
	for(Int_t i_sys = Sys_start; i_sys < Sys_stop; i_sys++) // Systematic loop
	{
	  TString KEY_Yield = Form("Yields_Centrality_%d_EtaGap_%d_%s_SysErrors_%d",i_cent,i_eta,PID[mPID].Data(),i_sys);
	  h_mYield_SM[KEY_Yield] = (TH1F*)h_mYield[KEY_Yield]->Clone();
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
	  ParYield_SM[KEY_Yield].clear();
	  h_mYield[KEY_Yield]->Fit(f_bw,"NQR");
	  for(Int_t n_par = 0; n_par < 5; n_par++)
	  {
	    ParYield_SM[KEY_Yield].push_back(static_cast<Float_t>(f_bw->GetParameter(n_par)));
	  }

	  TF1 *f_poly = new TF1("f_poly",Poly,BW_Start[mPID],BW_Stop[mPID],2);
	  f_poly->FixParameter(0,f_bw->GetParameter(3));
	  f_poly->FixParameter(1,f_bw->GetParameter(4));
	  h_mYield[KEY_Yield]->Add(f_poly,-1.0); // subtract linear background for phi differential InvMass
	}
      }
    }

    /*
    //QA: Yields subtract linear background
    TCanvas *c_Yields = new TCanvas("c_Yields","c_Yields",10,10,900,900);
    c_Yields->Divide(3,3);
    for(Int_t i_cent = 0; i_cent < 9; i_cent++)
    {
      c_Yields->cd(i_cent+1);
      c_Yields->cd(i_cent+1)->SetLeftMargin(0.20);
      c_Yields->cd(i_cent+1)->SetBottomMargin(0.20);
      c_Yields->cd(i_cent+1)->SetTicks(1,1);
      c_Yields->cd(i_cent+1)->SetGrid(0,0);

      TString KEY_Yield = Form("Yields_Centrality_%d_EtaGap_0_%s_SysErrors_%d",i_cent,PID[mPID].Data(),Sys_start);
      h_mYield[KEY_Yield]->SetTitle("");
      h_mYield[KEY_Yield]->SetStats(0);
      h_mYield[KEY_Yield]->SetMarkerStyle(24);
      h_mYield[KEY_Yield]->SetMarkerColor(kGray+3);
      h_mYield[KEY_Yield]->SetMarkerSize(0.8);
      h_mYield[KEY_Yield]->GetXaxis()->SetNdivisions(505,'N');
      h_mYield[KEY_Yield]->GetXaxis()->SetTitle("InvMass (GeV/c^{2})");
      h_mYield[KEY_Yield]->GetXaxis()->CenterTitle();
      h_mYield[KEY_Yield]->GetYaxis()->SetTitle("Counts");
      h_mYield[KEY_Yield]->GetYaxis()->CenterTitle();
      h_mYield[KEY_Yield]->GetYaxis()->SetTitleOffset(1.2);
      h_mYield[KEY_Yield]->DrawCopy("pE");
      PlotLine(0.98,1.05,0.0,0.0,1,2,2);

      h_mYield_SM[KEY_Yield]->SetMarkerStyle(24);
      h_mYield_SM[KEY_Yield]->SetMarkerColor(4);
      h_mYield_SM[KEY_Yield]->SetMarkerSize(0.8);
      h_mYield_SM[KEY_Yield]->DrawCopy("pE Same");

      TF1 *f_bw = new TF1("f_bw",PolyBreitWigner,BW_Start[mPID],BW_Stop[mPID],5);
      for(Int_t i_par = 0; i_par < 5; i_par++)
      {
	f_bw->ReleaseParameter(i_par);
      }
      f_bw->FixParameter(0,ParYield_SM[KEY_Yield][0]);
      f_bw->FixParameter(1,ParYield_SM[KEY_Yield][1]);
      f_bw->FixParameter(2,ParYield_SM[KEY_Yield][2]);
      f_bw->FixParameter(3,ParYield_SM[KEY_Yield][3]);
      f_bw->FixParameter(4,ParYield_SM[KEY_Yield][4]);
      f_bw->SetRange(BW_Start[mPID],BW_Stop[mPID]);
      f_bw->SetLineColor(2);
      f_bw->SetLineWidth(2);
      f_bw->SetLineStyle(1);
      f_bw->DrawCopy("l same");

      TF1 *f_poly = new TF1("f_poly",Poly,BW_Start[mPID],BW_Stop[mPID],2);
      f_poly->FixParameter(0,ParYield_SM[KEY_Yield][3]);
      f_poly->FixParameter(1,ParYield_SM[KEY_Yield][4]);
      f_poly->SetLineColor(2);
      f_poly->SetLineWidth(2);
      f_poly->SetLineStyle(2);
      f_poly->DrawCopy("l same");
    }
    */
  }

  // calculate total yields for each centrality bin via gaussian and breit wigner fits
  vecFMap ParYield_BW;
  vecFMap yields_Gaus, yields_BW;
  for(Int_t i_cent = 0; i_cent < 9; i_cent++) // Centrality loop
  {
    for(Int_t i_eta = Eta_start; i_eta < Eta_stop; i_eta++) // EtaGap loop
    {
      for(Int_t i_sys = Sys_start; i_sys < Sys_stop; i_sys++) // Systematic loop
      {
	TString KEY_Yield = Form("Yields_Centrality_%d_EtaGap_%d_%s_SysErrors_%d",i_cent,i_eta,PID[mPID].Data(),i_sys);
	TF1 *f_yields_bw = new TF1("f_yields_bw",BreitWigner,BW_Start[mPID],BW_Stop[mPID],3);
	f_yields_bw->SetParameter(0,InvMass[mPID]);
	f_yields_bw->SetParLimits(0,InvMass[mPID]-0.001,InvMass[mPID]+0.001);
	f_yields_bw->SetParameter(1,Width[mPID]);
	f_yields_bw->SetParameter(2,1000);
	f_yields_bw->SetRange(BW_Start[mPID],BW_Stop[mPID]);
	h_mYield[KEY_Yield]->Fit(f_yields_bw,"MQNR");
	ParYield_BW[KEY_Yield].clear();
	ParYield_BW[KEY_Yield].push_back(static_cast<Float_t>(f_yields_bw->GetParameter(0)));
	ParYield_BW[KEY_Yield].push_back(static_cast<Float_t>(f_yields_bw->GetParameter(1)));
	ParYield_BW[KEY_Yield].push_back(static_cast<Float_t>(f_yields_bw->GetParameter(2)));

	// counting for guassian
	Float_t counts_gaus = 0.0;
	Float_t errors_gaus = 0.0;
	Int_t bin_start = h_mYield[KEY_Yield]->FindBin(ParYield_BW[KEY_Yield][0]-nSigV0*ParYield_BW[KEY_Yield][1]);
	Int_t bin_stop  = h_mYield[KEY_Yield]->FindBin(ParYield_BW[KEY_Yield][0]+nSigV0*ParYield_BW[KEY_Yield][1]);
	for(Int_t i_bin = bin_start; i_bin <= bin_stop; i_bin++)
	{
	  counts_gaus += h_mYield[KEY_Yield]->GetBinContent(i_bin);
	  errors_gaus += h_mYield[KEY_Yield]->GetBinError(i_bin)*h_mYield[KEY_Yield]->GetBinError(i_bin);
	}
	yields_Gaus[KEY_Yield].clear();
	yields_Gaus[KEY_Yield].push_back(static_cast<Float_t>(counts_gaus));
	yields_Gaus[KEY_Yield].push_back(static_cast<Float_t>(TMath::Sqrt(errors_gaus)));

	// integrating for breit wigner
	Float_t bin_width = h_mYield[KEY_Yield]->GetBinWidth(1);
	Float_t Inte_start = ParYield_BW[KEY_Yield][0]-nSigV0*ParYield_BW[KEY_Yield][1]-0.5*bin_width;
	Float_t Inte_stop  = ParYield_BW[KEY_Yield][0]+nSigV0*ParYield_BW[KEY_Yield][1]+0.5*bin_width;
	Float_t counts_bw = f_yields_bw->Integral(Inte_start,Inte_stop)/bin_width;
	Float_t errors_bw = f_yields_bw->IntegralError(Inte_start,Inte_stop)/bin_width;
	yields_BW[KEY_Yield].clear();
	yields_BW[KEY_Yield].push_back(static_cast<Float_t>(counts_bw));
	yields_BW[KEY_Yield].push_back(static_cast<Float_t>(errors_bw));
      }
    }
  }

  /*
  // QA: different counting method: bin counting vs breit wigner integrating
  TCanvas *c_Yields_counts = new TCanvas("c_Yields_counts","c_Yields_counts",10,10,900,900);
  c_Yields_counts->Divide(3,3);
  for(Int_t i_cent = 0; i_cent < 9; i_cent++)
  {
    c_Yields_counts->cd(i_cent+1);
    c_Yields_counts->cd(i_cent+1)->SetLeftMargin(0.20);
    c_Yields_counts->cd(i_cent+1)->SetBottomMargin(0.20);
    c_Yields_counts->cd(i_cent+1)->SetTicks(1,1);
    c_Yields_counts->cd(i_cent+1)->SetGrid(0,0);

    TString KEY_Yield = Form("Yields_Centrality_%d_EtaGap_0_%s_SysErrors_%d",i_cent,PID[mPID].Data(),Sys_start);
    h_mYield[KEY_Yield]->SetTitle("");
    h_mYield[KEY_Yield]->SetStats(0);
    h_mYield[KEY_Yield]->SetMarkerStyle(24);
    h_mYield[KEY_Yield]->SetMarkerColor(kGray+3);
    h_mYield[KEY_Yield]->SetMarkerSize(0.8);
    h_mYield[KEY_Yield]->GetXaxis()->SetNdivisions(505,'N');
    h_mYield[KEY_Yield]->GetXaxis()->SetTitle("InvMass (GeV/c^{2})");
    h_mYield[KEY_Yield]->GetXaxis()->CenterTitle();
    h_mYield[KEY_Yield]->GetYaxis()->SetTitle("Counts");
    h_mYield[KEY_Yield]->GetYaxis()->CenterTitle();
    h_mYield[KEY_Yield]->GetYaxis()->SetTitleOffset(1.2);
    h_mYield[KEY_Yield]->DrawCopy("pE");
    PlotLine(0.98,1.05,0.0,0.0,1,2,2);

    TF1 *f_yields_bw = new TF1("f_yields_bw",BreitWigner,BW_Start[mPID],BW_Stop[mPID],3);
    f_yields_bw->SetParameter(0,ParYield_BW[KEY_Yield][0]);
    f_yields_bw->SetParameter(1,ParYield_BW[KEY_Yield][1]);
    f_yields_bw->SetParameter(2,ParYield_BW[KEY_Yield][2]);
    f_yields_bw->SetRange(BW_Start[mPID],BW_Stop[mPID]);
    f_yields_bw->SetLineColor(2);
    f_yields_bw->SetLineStyle(1);
    f_yields_bw->SetLineWidth(2);
    f_yields_bw->DrawCopy("l same");

    Float_t x1 = ParYield_BW[KEY_Yield][0] - nSigV0*ParYield_BW[KEY_Yield][1];
    Float_t x2 = ParYield_BW[KEY_Yield][0] + nSigV0*ParYield_BW[KEY_Yield][1];
    Float_t y = h_mYield[KEY_Yield]->GetBinContent(h_mYield[KEY_Yield]->FindBin(ParYield_BW[KEY_Yield][0]));
    PlotLine(x1,x1,0,y,4,2,2);
    PlotLine(x2,x2,0,y,4,2,2);
  }
  */

  // calculate final resolution correction factors and correct flow
  TString InPutFile_Res = Form("./Data/AuAu%s/file_%s_Resolution.root",Energy[mEnergy].Data(),Energy[mEnergy].Data());
  TFile *File_Res = TFile::Open(InPutFile_Res.Data());
  TString Res_Order[2] = {"Res2","Res3"};
  TProMap p_mRes;
  vecFMap ResValue;
  for(Int_t i_eta = Eta_start; i_eta < Eta_stop; i_eta++) // eta_gap loop
  {
    TString KEY_eta = Form("%s_EtaGap_%d_EP",Res_Order[mOrder].Data(),i_eta);
    p_mRes[KEY_eta] = (TProfile*)File_Res->Get(KEY_eta.Data()); // read in resolution file
    for(Int_t i_sys = Sys_start; i_sys < Sys_stop; i_sys++) // Systematic errors loop
    {
      for(Int_t i_cent = Cent_start; i_cent < Cent_stop; i_cent++) // centrality bin loop
      {
	Float_t yields_total_gaus = 0.0;
	Float_t yields_total_bw   = 0.0;
	for(Int_t cent = cent_low[i_cent]; cent <= cent_up[i_cent]; cent++) // calculate resolution and total yields in selected centrality bin
	{
	  if(p_mRes[KEY_eta]->GetBinContent(cent+1) > 0) 
	  {
	    TString KEY_Yield = Form("Yields_Centrality_%d_EtaGap_%d_%s_SysErrors_%d",cent,i_eta,PID[mPID].Data(),i_sys);
	    ResValue[KEY_Yield].push_back(static_cast<Float_t>(TMath::Sqrt(p_mRes[KEY_eta]->GetBinContent(cent+1))));
	    yields_total_gaus += yields_Gaus[KEY_Yield][0];
	    yields_total_bw += yields_BW[KEY_Yield][0];
	  }
	}

	TString KEY_ResCorr = Form("Res_%s_Centrality_%d_EtaGap_%d_SysErrors_%d",Order[mOrder].Data(),i_cent,i_eta,i_sys); // KEY for final resolution correction factor
	Float_t mean_res_gaus = 0.0;
	Float_t mean_res_bw = 0.0;
	for(Int_t cent = cent_low[i_cent]; cent <= cent_up[i_cent]; cent++) // calculate final resolution correction factor <1/R(centrality)>
	{
	  TString KEY_Yield = Form("Yields_Centrality_%d_EtaGap_%d_%s_SysErrors_%d",cent,i_eta,PID[mPID].Data(),i_sys);
	  mean_res_gaus += yields_Gaus[KEY_Yield][0]/(ResValue[KEY_Yield][0]*yields_total_gaus);
	  mean_res_bw += yields_BW[KEY_Yield][0]/(ResValue[KEY_Yield][0]*yields_total_bw);
	}
	cout << "i_eta = " << i_eta << ", i_sys = " << i_sys << ", centrality_bin = " << i_cent << ", mean_res_gaus = " << mean_res_gaus << endl;
	cout << "i_eta = " << i_eta << ", i_sys = " << i_sys << ", centrality_bin = " << i_cent << ", mean_res_bw = " << mean_res_bw << endl;

	TString KEY_RawFlow_Gaus = Form("RawFlow_Gaus_Centrality_%d_EtaGap_%d_%s_%s_SM_SysErrors_%d",i_cent,i_eta,Order[mOrder].Data(),PID[mPID].Data(),i_sys); // gaussian fits
	h_mRawFlow[KEY_RawFlow_Gaus]->Scale(mean_res_gaus);
	TString KEY_RawFlow_BW = Form("RawFlow_BW_Centrality_%d_EtaGap_%d_%s_%s_SM_SysErrors_%d",i_cent,i_eta,Order[mOrder].Data(),PID[mPID].Data(),i_sys); // breit wigner fits
	h_mRawFlow[KEY_RawFlow_BW]->Scale(mean_res_bw);
      }
    }
  }

  /*
  // QA flow vs. pt for gaussian fits and breit wigner fits
  TCanvas *c_flow = new TCanvas("c_flow","c_flow",10,10,800,800);
  c_flow->SetLeftMargin(0.15);
  c_flow->SetBottomMargin(0.15);
  c_flow->SetTicks(1,1);
  c_flow->SetGrid(0,0);
  TH1F *h_play = new TH1F("h_play","h_play",100,0.0,10.0);
  for(Int_t i_bin = 0; i_bin < 100; i_bin++)
  {
    h_play->SetBinContent(i_bin+1,-10.0);
    h_play->SetBinError(i_bin+1,1.0);
  }
  h_play->SetTitle("");
  h_play->SetStats(0);
  h_play->GetXaxis()->SetRangeUser(0.0,8.0);
  h_play->GetXaxis()->SetNdivisions(505,'N');
  h_play->GetXaxis()->SetLabelSize(0.03);
  h_play->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_play->GetXaxis()->SetTitleSize(0.05);
  h_play->GetXaxis()->SetTitleOffset(1.2);
  h_play->GetXaxis()->CenterTitle();

  h_play->GetYaxis()->SetRangeUser(-0.01,0.20);
  h_play->GetYaxis()->SetNdivisions(505,'N');
  h_play->GetYaxis()->SetTitle("v_{3}");
  h_play->GetYaxis()->SetTitleSize(0.05);
  h_play->GetYaxis()->SetLabelSize(0.03);
  h_play->GetYaxis()->CenterTitle();
  h_play->DrawCopy("pE");
  PlotLine(0.0,8.0,0.0,0.0,1,2,2);
  for(Int_t i_cent = Cent_start; i_cent < Cent_stop; i_cent++) // Centrality loop
  {
    for(Int_t i_eta = Eta_start; i_eta < Eta_stop; i_eta++) // EtaGap loop
    {
      for(Int_t i_sys = Sys_start; i_sys < Sys_stop; i_sys++) // Systematic loop
      {
	TString KEY_RawFlow_Gaus = Form("RawFlow_Gaus_Centrality_%d_EtaGap_%d_%s_%s_SM_SysErrors_%d",i_cent,i_eta,Order[mOrder].Data(),PID[mPID].Data(),i_sys); // gaussian fits
	h_mRawFlow[KEY_RawFlow_Gaus]->SetMarkerStyle(24);
	h_mRawFlow[KEY_RawFlow_Gaus]->SetMarkerColor(4);
	h_mRawFlow[KEY_RawFlow_Gaus]->DrawCopy("PE same");
        TString KEY_RawFlow_BW = Form("RawFlow_BW_Centrality_%d_EtaGap_%d_%s_%s_SM_SysErrors_%d",i_cent,i_eta,Order[mOrder].Data(),PID[mPID].Data(),i_sys); // breit wigner fits
	h_mRawFlow[KEY_RawFlow_BW]->SetMarkerStyle(24);
	h_mRawFlow[KEY_RawFlow_BW]->SetMarkerColor(2);
	h_mRawFlow[KEY_RawFlow_BW]->DrawCopy("PE same");
      }
    }
  }
  TString KEY_RawFlow_Gaus = Form("RawFlow_Gaus_Centrality_0_EtaGap_0_%s_%s_SM_SysErrors_%d",Order[mOrder].Data(),PID[mPID].Data(),Sys_start); // gaussian fits
  TString KEY_RawFlow_BW = Form("RawFlow_BW_Centrality_0_EtaGap_0_%s_%s_SM_SysErrors_%d",Order[mOrder].Data(),PID[mPID].Data(),Sys_start); // breit wigner fits
  TLegend *leg_flow = new TLegend(0.5,0.6,0.8,0.8);
  leg_flow->SetFillColor(10);
  leg_flow->SetBorderSize(0.0);
  leg_flow->AddEntry(h_mRawFlow[KEY_RawFlow_Gaus],"bin counting","p");
  leg_flow->AddEntry(h_mRawFlow[KEY_RawFlow_BW],"breit wigner","p");
  leg_flow->Draw("same");
//  c_flow->SaveAs("./flow.eps");
  */

  // read in pt spectra
  TString InPutFile_Pt = Form("./OutPut/AuAu%s/%s/h_pt.root",Energy[mEnergy].Data(),PID[mPID].Data());
  TFile *File_Spec = TFile::Open(InPutFile_Pt.Data());
  TH1FMap h_mPt;
  for(Int_t i_cent = Cent_start; i_cent < Cent_stop; i_cent++) // Centrality loop
  {
    for(Int_t i_eta = Eta_start; i_eta < Eta_stop; i_eta++) // EtaGap loop
    {
      for(Int_t i_sys = Sys_start; i_sys < Sys_stop; i_sys++) // Systematic loop
      {
	TString KEY_pT_counts = Form("Count_Spec_Centrality_%d_EtaGap_%d_%s_SysErrors_%d",i_cent,i_eta,PID[mPID].Data(),i_sys);
	h_mPt[KEY_pT_counts] = (TH1F*)File_Spec->Get(KEY_pT_counts.Data());
	TString KEY_pT_inte = Form("Inte_Spec_Centrality_%d_EtaGap_%d_%s_SysErrors_%d",i_cent,i_eta,PID[mPID].Data(),i_sys);
	h_mPt[KEY_pT_inte] = (TH1F*)File_Spec->Get(KEY_pT_inte.Data());
      }
    }
  }

  /*
  // QA: pt spectra
  TCanvas *c_Pt = new TCanvas("c_Pt","c_Pt",10,10,800,800);
  c_Pt->cd();
  c_Pt->cd()->SetLeftMargin(0.20);
  c_Pt->cd()->SetBottomMargin(0.20);
  c_Pt->cd()->SetTicks(1,1);
  c_Pt->cd()->SetGrid(0,0);
  c_Pt->cd()->SetLogy();
  TString KEY_pT_counts_QA = Form("Count_Spec_Centrality_%d_EtaGap_%d_%s_SysErrors_%d",Cent_start,Eta_start,PID[mPID].Data(),Sys_start);
  TH1F *h_play = new TH1F("h_play","h_play",110,-1.0,10.0);
  for(Int_t i_bin = 0; i_bin < 110; i_bin++)
  {
    h_play->SetBinContent(i_bin+1,-1000.0);
    h_play->SetBinError(i_bin+1,1.0);
  }
  h_play->SetTitle("");
  h_play->SetStats(0);
  h_play->GetXaxis()->SetNdivisions(505,'N');
  h_play->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_play->GetXaxis()->CenterTitle();
  h_play->GetYaxis()->SetTitle("dN/dp_{T}");
  h_play->GetYaxis()->CenterTitle();
  h_play->GetYaxis()->SetTitleOffset(1.2);
  h_play->GetYaxis()->SetRangeUser(0.1,2.0*h_mPt[KEY_pT_counts_QA]->GetMaximum());
  h_play->DrawCopy("pE");
  h_mPt[KEY_pT_counts_QA]->SetMarkerStyle(24);
  h_mPt[KEY_pT_counts_QA]->SetMarkerColor(4);
  h_mPt[KEY_pT_counts_QA]->SetMarkerSize(1.0);
  h_mPt[KEY_pT_counts_QA]->DrawCopy("pE same");

  TString KEY_pT_inte_QA = Form("Inte_Spec_Centrality_%d_EtaGap_%d_%s_SysErrors_%d",Cent_start,Eta_start,PID[mPID].Data(),Sys_start);
  h_mPt[KEY_pT_inte_QA]->SetMarkerColor(2);
  h_mPt[KEY_pT_inte_QA]->SetMarkerStyle(24);
  h_mPt[KEY_pT_inte_QA]->SetMarkerSize(1.0);
  h_mPt[KEY_pT_inte_QA]->DrawCopy("pE same");
  */

  // pt spectra binning
  Float_t pt_low_spec[50], pt_up_spec[50], pt_width[50], pt_center[50];
  for(Int_t i_pt = pt_start; i_pt < pt_stop; i_pt++)
  {
    pt_low_spec[2*i_pt] = pt_low_raw[i_pt];
    pt_up_spec[2*i_pt]  = 0.5*(pt_low_raw[i_pt]+pt_up_raw[i_pt]);
    pt_width[2*i_pt] = pt_up_spec[2*i_pt]-pt_low_spec[2*i_pt];
    pt_center[2*i_pt] = 0.5*(pt_up_spec[2*i_pt]+pt_low_spec[2*i_pt]);
    
    pt_low_spec[2*i_pt+1] = 0.5*(pt_low_raw[i_pt]+pt_up_raw[i_pt]);
    pt_up_spec[2*i_pt+1]  = pt_up_raw[i_pt];
    pt_width[2*i_pt+1] = pt_up_spec[2*i_pt+1]-pt_low_spec[2*i_pt+1];
    pt_center[2*i_pt+1] = 0.5*(pt_up_spec[2*i_pt+1]+pt_low_spec[2*i_pt+1]);
  }

  vecFMap mean_pt; // mean pt with systematic errors
  for(Int_t i_cent = Cent_start; i_cent < Cent_stop; i_cent++) // Centrality loop
  {
    for(Int_t i_eta = Eta_start; i_eta < Eta_stop; i_eta++) // EtaGap loop
    {
      for(Int_t i_sys = Sys_start; i_sys < Sys_stop; i_sys++) // Systematic loop
      {
	TString KEY_RawFlow_Gaus_Save = Form("RawFlow_Gaus_Centrality_%d_EtaGap_%d_%s_%s_SM_SysErrors_%d",i_cent,i_eta,Order[mOrder].Data(),PID[mPID].Data(),i_sys); // gaussian fits
	mean_pt[KEY_RawFlow_Gaus_Save].clear();
	TString KEY_RawFlow_BW_Save = Form("RawFlow_BW_Centrality_%d_EtaGap_%d_%s_%s_SM_SysErrors_%d",i_cent,i_eta,Order[mOrder].Data(),PID[mPID].Data(),i_sys); // breit wigner fits
	mean_pt[KEY_RawFlow_BW_Save].clear();
	for(Int_t i_pt = pt_rebin_first; i_pt < pt_rebin_last; i_pt++) // loop over rebinned pT
	{
//	  cout << "i_pt = " << i_pt << endl;
	  Float_t mean_pt_counts = 0.0;
	  Float_t spec_counts    = 0.0;
	  Float_t mean_pt_inte   = 0.0;
	  Float_t spec_inte      = 0.0;
	  for(Int_t i_pt_raw = pt_rebin_start[i_pt]; i_pt_raw <= pt_rebin_stop[i_pt]; i_pt_raw++) // loop over raw pT bin
	  {
//	    cout << "i_pt_raw = " << i_pt_raw << endl;
	    for(Int_t i_pt_spec = 2*i_pt_raw; i_pt_spec <= 2*i_pt_raw+1; i_pt_spec++)
	    {
//	      cout << "i_pt_spec = " << i_pt_spec << endl;
	      TString KEY_pT_counts = Form("Count_Spec_Centrality_%d_EtaGap_%d_%s_SysErrors_%d",i_cent,i_eta,PID[mPID].Data(),i_sys);
	      mean_pt_counts += pt_center[i_pt_spec]*h_mPt[KEY_pT_counts]->GetBinContent(i_pt_spec+1)*pt_width[i_pt_spec];
	      spec_counts    += h_mPt[KEY_pT_counts]->GetBinContent(i_pt_spec+1)*pt_width[i_pt_spec];

	      TString KEY_pT_inte = Form("Inte_Spec_Centrality_%d_EtaGap_%d_%s_SysErrors_%d",i_cent,i_eta,PID[mPID].Data(),i_sys);
	      mean_pt_inte += pt_center[i_pt_spec]*h_mPt[KEY_pT_inte]->GetBinContent(i_pt_spec+1)*pt_width[i_pt_spec];
	      spec_inte    += h_mPt[KEY_pT_inte]->GetBinContent(i_pt_spec+1)*pt_width[i_pt_spec];
	    }
	  }
	  mean_pt[KEY_RawFlow_Gaus_Save].push_back(static_cast<Float_t>(mean_pt_counts/spec_counts));
	  mean_pt[KEY_RawFlow_BW_Save].push_back(static_cast<Float_t>(mean_pt_inte/spec_inte));
	}
      }
    }
  }

  // set final pt and flow to one TGraphAsymmErrors
  TGraMap g_mFlow;
  for(Int_t i_cent = Cent_start; i_cent < Cent_stop; i_cent++) // Centrality loop
  {
    for(Int_t i_eta = Eta_start; i_eta < Eta_stop; i_eta++) // EtaGap loop
    {
      for(Int_t i_sys = Sys_start; i_sys < Sys_stop; i_sys++) // Systematic loop
      {
	TString KEY_RawFlow_Gaus_Save = Form("RawFlow_Gaus_Centrality_%d_EtaGap_%d_%s_%s_SM_SysErrors_%d",i_cent,i_eta,Order[mOrder].Data(),PID[mPID].Data(),i_sys); // gaussian fits
	g_mFlow[KEY_RawFlow_Gaus_Save] = new TGraphAsymmErrors();
	TString KEY_RawFlow_BW_Save = Form("RawFlow_BW_Centrality_%d_EtaGap_%d_%s_%s_SM_SysErrors_%d",i_cent,i_eta,Order[mOrder].Data(),PID[mPID].Data(),i_sys); // breit wigner fits
	g_mFlow[KEY_RawFlow_BW_Save] = new TGraphAsymmErrors();
	for(Int_t i_pt = pt_rebin_first; i_pt < pt_rebin_last; i_pt++)
	{
	  Float_t pt_mean = (pt_low[i_pt]+pt_up[i_pt])/2.0;

	  Float_t gauss_content = h_mRawFlow[KEY_RawFlow_Gaus_Save]->GetBinContent(h_mRawFlow[KEY_RawFlow_Gaus_Save]->FindBin(pt_mean)); // bin counting
	  Float_t gauss_error   = h_mRawFlow[KEY_RawFlow_Gaus_Save]->GetBinError(h_mRawFlow[KEY_RawFlow_Gaus_Save]->FindBin(pt_mean));
	  g_mFlow[KEY_RawFlow_Gaus_Save]->SetPoint(i_pt,mean_pt[KEY_RawFlow_Gaus_Save][i_pt],gauss_content);
	  g_mFlow[KEY_RawFlow_Gaus_Save]->SetPointError(i_pt,0.0,0.0,gauss_error,gauss_error);
	  TString Name_gaus = Form("Flow_%s_%s_Centrality_%d_EtaGap_%d_SysErrors_%d_Gaus",Order[mOrder].Data(),PID[mPID].Data(),i_cent,i_eta,i_sys); // gaussian fits
	  g_mFlow[KEY_RawFlow_Gaus_Save]->SetName(Name_gaus.Data());

	  Float_t bw_content = h_mRawFlow[KEY_RawFlow_BW_Save]->GetBinContent(h_mRawFlow[KEY_RawFlow_BW_Save]->FindBin(pt_mean)); // breit wigner fits
	  Float_t bw_error   = h_mRawFlow[KEY_RawFlow_BW_Save]->GetBinError(h_mRawFlow[KEY_RawFlow_BW_Save]->FindBin(pt_mean));
	  g_mFlow[KEY_RawFlow_BW_Save]->SetPoint(i_pt,pt_mean,bw_content);
	  g_mFlow[KEY_RawFlow_BW_Save]->SetPointError(i_pt,0.0,0.0,bw_error,bw_error);
	  TString Name_bw = Form("Flow_%s_%s_Centrality_%d_EtaGap_%d_SysErrors_%d_Gaus",Order[mOrder].Data(),PID[mPID].Data(),i_cent,i_eta,i_sys); // gaussian fits
	  g_mFlow[KEY_RawFlow_BW_Save]->SetName(Name_bw.Data());
	}
      }
    }
  }

  // QA flow vs. pt for gaussian fits and breit wigner fits with TGraphAsymmErrors
  TCanvas *c_flow = new TCanvas("c_flow","c_flow",10,10,800,800);
  c_flow->SetLeftMargin(0.15);
  c_flow->SetBottomMargin(0.15);
  c_flow->SetTicks(1,1);
  c_flow->SetGrid(0,0);
  TH1F *h_play = new TH1F("h_play","h_play",100,0.0,10.0);
  for(Int_t i_bin = 0; i_bin < 100; i_bin++)
  {
    h_play->SetBinContent(i_bin+1,-10.0);
    h_play->SetBinError(i_bin+1,1.0);
  }
  h_play->SetTitle("");
  h_play->SetStats(0);
  h_play->GetXaxis()->SetRangeUser(0.0,8.0);
  h_play->GetXaxis()->SetNdivisions(505,'N');
  h_play->GetXaxis()->SetLabelSize(0.03);
  h_play->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_play->GetXaxis()->SetTitleSize(0.05);
  h_play->GetXaxis()->SetTitleOffset(1.2);
  h_play->GetXaxis()->CenterTitle();

  h_play->GetYaxis()->SetRangeUser(-0.01,0.20);
  h_play->GetYaxis()->SetNdivisions(505,'N');
  h_play->GetYaxis()->SetTitle("v_{3}");
  h_play->GetYaxis()->SetTitleSize(0.05);
  h_play->GetYaxis()->SetLabelSize(0.03);
  h_play->GetYaxis()->CenterTitle();
  h_play->DrawCopy("pE");
  PlotLine(0.0,8.0,0.0,0.0,1,2,2);
  for(Int_t i_cent = Cent_start; i_cent < Cent_stop; i_cent++) // Centrality loop
  {
    for(Int_t i_eta = Eta_start; i_eta < Eta_stop; i_eta++) // EtaGap loop
    {
      for(Int_t i_sys = Sys_start; i_sys < Sys_stop; i_sys++) // Systematic loop
      {
	TString KEY_RawFlow_Gaus = Form("RawFlow_Gaus_Centrality_%d_EtaGap_%d_%s_%s_SM_SysErrors_%d",i_cent,i_eta,Order[mOrder].Data(),PID[mPID].Data(),i_sys); // gaussian fits
	g_mFlow[KEY_RawFlow_Gaus]->SetMarkerStyle(24);
	g_mFlow[KEY_RawFlow_Gaus]->SetMarkerColor(4);
	g_mFlow[KEY_RawFlow_Gaus]->Draw("PE same");
        TString KEY_RawFlow_BW = Form("RawFlow_BW_Centrality_%d_EtaGap_%d_%s_%s_SM_SysErrors_%d",i_cent,i_eta,Order[mOrder].Data(),PID[mPID].Data(),i_sys); // breit wigner fits
	g_mFlow[KEY_RawFlow_BW]->SetMarkerStyle(24);
	g_mFlow[KEY_RawFlow_BW]->SetMarkerColor(2);
	g_mFlow[KEY_RawFlow_BW]->Draw("PE same");
      }
    }
  }
  TString KEY_RawFlow_Gaus = Form("RawFlow_Gaus_Centrality_0_EtaGap_0_%s_%s_SM_SysErrors_%d",Order[mOrder].Data(),PID[mPID].Data(),Sys_start); // gaussian fits
  TString KEY_RawFlow_BW = Form("RawFlow_BW_Centrality_0_EtaGap_0_%s_%s_SM_SysErrors_%d",Order[mOrder].Data(),PID[mPID].Data(),Sys_start); // breit wigner fits
  TLegend *leg_flow = new TLegend(0.5,0.6,0.7,0.7);
  leg_flow->SetFillColor(10);
  leg_flow->SetBorderSize(0.0);
  leg_flow->AddEntry(g_mFlow[KEY_RawFlow_Gaus],"bin counting","p");
  leg_flow->AddEntry(g_mFlow[KEY_RawFlow_BW],"breit wigner","p");
  leg_flow->Draw("same");
//  c_flow->SaveAs("./flow.eps");

  TString OutPutFile = Form("./OutPut/AuAu%s/%s/flow_%s.root",Energy[mEnergy].Data(),PID[mPID].Data(),Order[mOrder].Data());
  TFile *File_OutPut = new TFile(OutPutFile.Data(),"RECREATE");
  File_OutPut->cd();
  for(Int_t i_cent = Cent_start; i_cent < Cent_stop; i_cent++) // Centrality loop
  {
    for(Int_t i_eta = Eta_start; i_eta < Eta_stop; i_eta++) // EtaGap loop
    {
      for(Int_t i_sys = Sys_start; i_sys < Sys_stop; i_sys++) // Systematic loop
      {
	TString KEY_RawFlow_Gaus_Save = Form("RawFlow_Gaus_Centrality_%d_EtaGap_%d_%s_%s_SM_SysErrors_%d",i_cent,i_eta,Order[mOrder].Data(),PID[mPID].Data(),i_sys); // gaussian fits
	g_mFlow[KEY_RawFlow_Gaus_Save]->SetMarkerStyle(24);
	g_mFlow[KEY_RawFlow_Gaus_Save]->SetMarkerColor(4);
	g_mFlow[KEY_RawFlow_Gaus_Save]->Write();
	TString KEY_RawFlow_BW_Save = Form("RawFlow_BW_Centrality_%d_EtaGap_%d_%s_%s_SM_SysErrors_%d",i_cent,i_eta,Order[mOrder].Data(),PID[mPID].Data(),i_sys); // breit wigner fits
	g_mFlow[KEY_RawFlow_BW_Save]->SetMarkerStyle(24);
	g_mFlow[KEY_RawFlow_BW_Save]->SetMarkerColor(2);
	g_mFlow[KEY_RawFlow_BW_Save]->Write();
      }
    }
  }
  h_play->Write();
  File_OutPut->Close();
}

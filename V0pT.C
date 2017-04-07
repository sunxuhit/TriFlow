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
static const TString PID[4] = {"Phi","Lambda","AntiLambda","K0S"};
static const TString Order[2] = {"2nd","3rd"};
static const Float_t Norm_Start[4] = {1.04,1.14,1.14,0.56};
static const Float_t Norm_Stop[4]  = {1.05,1.19,1.19,0.60};
static const Float_t BW_Start[4] = {0.994,1.106,1.106,0.472};
static const Float_t BW_Stop[4]  = {1.050,1.122,1.122,0.524};
static const Float_t InvMass[4] = {1.019,1.116,1.116,0.498};
static const Float_t Width[4]   = {0.00426,0.0016,0.0016,0.0016};
static const Double_t PI_max[2] = {TMath::Pi()/2.0,TMath::Pi()/3.0};
static const Float_t nSigV0 = 2.0;
static const Float_t Flow_Order[2] = {2.0,3.0};

static const TString pt_range[2] = {"low","high"};
static const Int_t pt_total = 25; // pT loop
static const Int_t pt_start = 0;
static const Int_t pt_stop  = 25;
static const Int_t pt_QA    = 4;

// pt bin
//                                       0 , 1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 ,10 ,11 ,12 ,13 ,14 ,15 ,16 ,17 ,18 ,19 ,20 ,21 ,22, 23, 24
static const Float_t pt_low_raw[25] = {0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.4,5.8,6.2,6.6,7.2};
static const Float_t pt_up_raw[25]  = {0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.4,5.8,6.2,6.6,7.2,8.0};

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
typedef std::map<TString,std::vector<Float_t>> vecFMap;

#ifndef _PlotQA_
#define _PlotQA_ 1
#endif

void V0pT(Int_t mEnergy = 0, Int_t mPID = 0)
{
  TGaxis::SetMaxDigits(4);

  TString InPutFile_SE = Form("./Data/AuAu%s/%s/Yields_SE_%s.root",Energy[mEnergy].Data(),PID[mPID].Data(),Energy[mEnergy].Data());
  TString InPutFile_ME = Form("./Data/AuAu%s/%s/Yields_ME_%s.root",Energy[mEnergy].Data(),PID[mPID].Data(),Energy[mEnergy].Data());
  TFile *File_SE = TFile::Open(InPutFile_SE.Data());
  TFile *File_ME = TFile::Open(InPutFile_ME.Data());

  // read in histogram for same event and mixed event
  // calculate SE - ME
  TH1FMap h_mSpec_SE, h_mSpec_ME, h_mSpec;
  for(Int_t i_pt = pt_start; i_pt < pt_stop; i_pt++) // pt loop
  {
    for(Int_t i_range = 0; i_range < 2; i_range++)
    {
      for(Int_t i_cent = Cent_start; i_cent < Cent_stop; i_cent++) // Centrality loop
      {
	for(Int_t i_eta = Eta_start; i_eta < Eta_stop; i_eta++) // EtaGap loop
	{
	  for(Int_t i_sys = Sys_start; i_sys < Sys_stop; i_sys++) // Systematic loop
	  {
	    TString KEY_SE = Form("Spec_pt_%d_%s_Centrality_%d_EtaGap_%d_%s_SE_SysErrors_%d",i_pt,pt_range[i_range].Data(),i_cent,i_eta,PID[mPID].Data(),i_sys);
	    TString Hist_SE = KEY_SE;
	    if(mPID == 3) Hist_SE = KEY_SE+"_sub";
	    h_mSpec_SE[KEY_SE] = (TH1F*)File_SE->Get(Hist_SE.Data())->Clone(); 
	    Int_t Norm_bin_start = h_mSpec_SE[KEY_SE]->FindBin(Norm_Start[mPID]);
	    Int_t Norm_bin_stop  = h_mSpec_SE[KEY_SE]->FindBin(Norm_Stop[mPID]);
	    Float_t Inte_SE = h_mSpec_SE[KEY_SE]->Integral(Norm_bin_start,Norm_bin_stop);

	    TString KEY_ME = Form("Spec_pt_%d_%s_Centrality_%d_EtaGap_%d_%s_ME_SysErrors_%d",i_pt,pt_range[i_range].Data(),i_cent,i_eta,PID[mPID].Data(),i_sys);
	    TString Hist_ME = KEY_ME;
	    if(mPID == 3) Hist_ME = KEY_ME+"_sub";
	    h_mSpec_ME[KEY_ME] = (TH1F*)File_ME->Get(Hist_ME.Data())->Clone(); 
	    Float_t Inte_ME = h_mSpec_ME[KEY_ME]->Integral(Norm_bin_start,Norm_bin_stop);
	    h_mSpec_ME[KEY_ME]->Scale(Inte_SE/Inte_ME);

	    TString KEY = Form("Spec_pt_%d_%s_Centrality_%d_EtaGap_%d_%s_SysErrors_%d",i_pt,pt_range[i_range].Data(),i_cent,i_eta,PID[mPID].Data(),i_sys);
	    h_mSpec[KEY] = (TH1F*)h_mSpec_SE[KEY_SE]->Clone();
	    h_mSpec[KEY]->Add(h_mSpec_ME[KEY_ME],-1.0);
	  }
	}
      }
    }
  }

#if _PlotQA_
  // QA Plots for SE vs. ME
  TString KEY_QA = Form("Spec_pt_%d_low_Centrality_0_EtaGap_0_%s_SE_SysErrors_%d",pt_QA,PID[mPID].Data(),Sys_start);
  h_mSpec_SE[KEY_QA]->DrawCopy("PE");

  KEY_QA = Form("Spec_pt_%d_low_Centrality_0_EtaGap_0_%s_ME_SysErrors_%d",pt_QA,PID[mPID].Data(),Sys_start);
  h_mSpec_ME[KEY_QA]->SetLineColor(2);
  h_mSpec_ME[KEY_QA]->SetFillColor(2);
  h_mSpec_ME[KEY_QA]->SetFillStyle(3002);
  h_mSpec_ME[KEY_QA]->DrawCopy("h same");

  KEY_QA = Form("Spec_pt_%d_low_Centrality_0_EtaGap_0_%s_SysErrors_%d",pt_QA,PID[mPID].Data(),Sys_start);
  h_mSpec[KEY_QA]->SetLineColor(4);
  h_mSpec[KEY_QA]->SetFillColor(4);
  h_mSpec[KEY_QA]->SetFillStyle(3004);
  h_mSpec[KEY_QA]->DrawCopy("h same");
#endif

#if _PlotQA_
  // QA Plots for pT bins
  TCanvas *c_pT_low = new TCanvas("c_pT_low","c_pT_low",10,10,1400,1400);
  c_pT_low->Divide(5,5);
  for(Int_t i_pt = pt_start; i_pt < pt_stop; i_pt++) // pt loop
  {
    c_pT_low->cd(i_pt+1);
    c_pT_low->cd(i_pt+1)->SetLeftMargin(0.15);
    c_pT_low->cd(i_pt+1)->SetBottomMargin(0.15);
    c_pT_low->cd(i_pt+1)->SetTicks(1,1);
    c_pT_low->cd(i_pt+1)->SetGrid(0,0);

    TString KEY_pT_QA = Form("Spec_pt_%d_low_Centrality_0_EtaGap_0_%s_SE_SysErrors_%d",i_pt,PID[mPID].Data(),Sys_start);
    h_mSpec_SE[KEY_pT_QA]->DrawCopy("PE");

    KEY_pT_QA = Form("Spec_pt_%d_low_Centrality_0_EtaGap_0_%s_ME_SysErrors_%d",i_pt,PID[mPID].Data(),Sys_start);
    h_mSpec_ME[KEY_pT_QA]->SetLineColor(2);
    h_mSpec_ME[KEY_pT_QA]->SetFillColor(2);
    h_mSpec_ME[KEY_pT_QA]->SetFillStyle(3002);
    h_mSpec_ME[KEY_pT_QA]->DrawCopy("h same");

    KEY_pT_QA = Form("Spec_pt_%d_low_Centrality_0_EtaGap_0_%s_SysErrors_%d",i_pt,PID[mPID].Data(),Sys_start);
    h_mSpec[KEY_pT_QA]->SetLineColor(4);
    h_mSpec[KEY_pT_QA]->SetFillColor(4);
    h_mSpec[KEY_pT_QA]->SetFillStyle(3004);
    h_mSpec[KEY_pT_QA]->DrawCopy("h same");

    TString pT_range = Form("[%.2f,%.2f]",pt_low_raw[i_pt],0.5*(pt_low_raw[i_pt]+pt_up_raw[i_pt]));
    plotTopLegend((char*)pT_range.Data(),0.2,0.7,0.08,1,0.0,42,1);
  }

  TCanvas *c_pT_high = new TCanvas("c_pT_high","c_pT_high",10,10,1400,1400);
  c_pT_high->Divide(5,5);
  for(Int_t i_pt = pt_start; i_pt < pt_stop; i_pt++) // pt loop
  {
    c_pT_high->cd(i_pt+1);
    c_pT_high->cd(i_pt+1)->SetLeftMargin(0.15);
    c_pT_high->cd(i_pt+1)->SetBottomMargin(0.15);
    c_pT_high->cd(i_pt+1)->SetTicks(1,1);
    c_pT_high->cd(i_pt+1)->SetGrid(0,0);

    TString KEY_pT_QA = Form("Spec_pt_%d_high_Centrality_0_EtaGap_0_%s_SE_SysErrors_%d",i_pt,PID[mPID].Data(),Sys_start);
    h_mSpec_SE[KEY_pT_QA]->DrawCopy("PE");

    KEY_pT_QA = Form("Spec_pt_%d_high_Centrality_0_EtaGap_0_%s_ME_SysErrors_%d",i_pt,PID[mPID].Data(),Sys_start);
    h_mSpec_ME[KEY_pT_QA]->SetLineColor(2);
    h_mSpec_ME[KEY_pT_QA]->SetFillColor(2);
    h_mSpec_ME[KEY_pT_QA]->SetFillStyle(3002);
    h_mSpec_ME[KEY_pT_QA]->DrawCopy("h same");

    KEY_pT_QA = Form("Spec_pt_%d_high_Centrality_0_EtaGap_0_%s_SysErrors_%d",i_pt,PID[mPID].Data(),Sys_start);
    h_mSpec[KEY_pT_QA]->SetLineColor(4);
    h_mSpec[KEY_pT_QA]->SetFillColor(4);
    h_mSpec[KEY_pT_QA]->SetFillStyle(3004);
    h_mSpec[KEY_pT_QA]->DrawCopy("h same");

    TString pT_range = Form("[%.2f,%.2f]",0.5*(pt_low_raw[i_pt]+pt_up_raw[i_pt]),pt_up_raw[i_pt]);
    plotTopLegend((char*)pT_range.Data(),0.2,0.7,0.08,1,0.0,42,1);
  }
#endif

  if(mPID == 0) // subtract linear background only needed for phi meson
  {
    // Poly + Breit Wignar fit to InvMass
    vecFMap ParFit;
    TH1FMap h_mSpec_QA;

    for(Int_t i_pt = pt_start; i_pt < pt_stop; i_pt++) // pt loop
    {
      for(Int_t i_range = 0; i_range < 2; i_range++)
      {
	for(Int_t i_cent = Cent_start; i_cent < Cent_stop; i_cent++) // Centrality loop
	{
	  for(Int_t i_eta = Eta_start; i_eta < Eta_stop; i_eta++) // EtaGap loop
	  {
	    for(Int_t i_sys = Sys_start; i_sys < Sys_stop; i_sys++) // Systematic loop
	    {
	      TString KEY = Form("Spec_pt_%d_%s_Centrality_%d_EtaGap_%d_%s_SysErrors_%d",i_pt,pt_range[i_range].Data(),i_cent,i_eta,PID[mPID].Data(),i_sys);
	      h_mSpec_QA[KEY] = (TH1F*)h_mSpec[KEY]->Clone();
	      TF1 *f_bw = new TF1("f_bw",PolyBreitWigner,BW_Start[mPID],BW_Stop[mPID],5);
	      for(Int_t i_par = 0; i_par < 5; i_par++)
	      {
		f_bw->ReleaseParameter(i_par);
	      }
	      f_bw->SetParameter(0,InvMass[mPID]);
	      f_bw->SetParLimits(0,1.014,1.024);
	      f_bw->SetParameter(1,Width[mPID]);
	      f_bw->SetParameter(2,10000);
	      f_bw->SetParameter(3,-6000);
	      f_bw->SetParameter(4,0.5);
	      f_bw->SetRange(BW_Start[mPID],BW_Stop[mPID]);
	      ParFit[KEY].clear();
	      h_mSpec[KEY]->Fit(f_bw,"NQR");
	      for(Int_t n_par = 0; n_par < 5; n_par++)
	      {
		ParFit[KEY].push_back(static_cast<Float_t>(f_bw->GetParameter(n_par)));
	      }

	      TF1 *f_poly = new TF1("f_poly",Poly,BW_Start[mPID],BW_Stop[mPID],2); // poly fit for linear background
	      f_poly->SetParameter(0,ParFit[KEY][3]);
	      f_poly->SetParameter(1,ParFit[KEY][4]);

	      h_mSpec[KEY]->Add(f_poly,-1.0); // subtract linear background for phi differential InvMass
	    }
	  }
	}
      }
    }

#if _PlotQA_
    // QA plots for Poly+Breit_Wignar fits for phi integrated InvMass
    TCanvas *c_mSpec_sub_low = new TCanvas("c_mSpec_sub_low","c_mSpec_sub_low",10,10,1400,1400);
    c_mSpec_sub_low->Divide(5,5);
    for(Int_t i_pt = pt_start; i_pt < pt_stop; i_pt++) // pt loop
    {
      c_mSpec_sub_low->cd(i_pt+1);
      c_mSpec_sub_low->cd(i_pt+1)->SetLeftMargin(0.15);
      c_mSpec_sub_low->cd(i_pt+1)->SetBottomMargin(0.15);
      c_mSpec_sub_low->cd(i_pt+1)->SetTicks(1,1);
      c_mSpec_sub_low->cd(i_pt+1)->SetGrid(0,0);
      TString KEY_QA = Form("Spec_pt_%d_low_Centrality_%d_EtaGap_%d_%s_SysErrors_%d",i_pt,Cent_start,Eta_start,PID[mPID].Data(),Sys_start);
      h_mSpec[KEY_QA]->SetMarkerColor(1);
      h_mSpec[KEY_QA]->SetMarkerStyle(24);
      h_mSpec[KEY_QA]->SetMarkerSize(0.8);
      h_mSpec[KEY_QA]->DrawCopy("PE");

      h_mSpec_QA[KEY_QA]->SetMarkerColor(4);
      h_mSpec_QA[KEY_QA]->SetMarkerStyle(24);
      h_mSpec_QA[KEY_QA]->SetMarkerSize(0.8);
      h_mSpec_QA[KEY_QA]->DrawCopy("PE same");

      TF1 *f_bw = new TF1("f_bw",PolyBreitWigner,BW_Start[mPID],BW_Stop[mPID],5);
      for(Int_t i_par = 0; i_par < 5; i_par++)
      {
	f_bw->SetParameter(i_par,ParFit[KEY_QA][i_par]);
      }
      f_bw->SetLineColor(2);
      f_bw->SetLineStyle(1);
      f_bw->SetLineWidth(2);
      f_bw->DrawCopy("l same");

      TF1 *f_poly = new TF1("f_poly",Poly,BW_Start[mPID],BW_Stop[mPID],2);
      f_poly->SetParameter(0,ParFit[KEY_QA][3]);
      f_poly->SetParameter(1,ParFit[KEY_QA][4]);
      f_poly->SetLineColor(2);
      f_poly->SetLineStyle(2);
      f_poly->SetLineWidth(4);
      f_poly->DrawCopy("l same");

      TString pT_range = Form("[%.2f,%.2f]",pt_low_raw[i_pt],0.5*(pt_low_raw[i_pt]+pt_up_raw[i_pt]));
      plotTopLegend((char*)pT_range.Data(),0.2,0.7,0.08,1,0.0,42,1);
      PlotLine(0.98,1.05,0.0,0.0,1,2,2);
    }

    TCanvas *c_mSpec_sub_high = new TCanvas("c_mSpec_sub_high","c_mSpec_sub_high",10,10,1400,1400);
    c_mSpec_sub_high->Divide(5,5);
    for(Int_t i_pt = pt_start; i_pt < pt_stop; i_pt++) // pt loop
    {
      c_mSpec_sub_high->cd(i_pt+1);
      c_mSpec_sub_high->cd(i_pt+1)->SetLeftMargin(0.15);
      c_mSpec_sub_high->cd(i_pt+1)->SetBottomMargin(0.15);
      c_mSpec_sub_high->cd(i_pt+1)->SetTicks(1,1);
      c_mSpec_sub_high->cd(i_pt+1)->SetGrid(0,0);
      TString KEY_QA = Form("Spec_pt_%d_high_Centrality_%d_EtaGap_%d_%s_SysErrors_%d",i_pt,Cent_start,Eta_start,PID[mPID].Data(),Sys_start);
      h_mSpec[KEY_QA]->SetMarkerColor(1);
      h_mSpec[KEY_QA]->SetMarkerStyle(24);
      h_mSpec[KEY_QA]->SetMarkerSize(0.8);
      h_mSpec[KEY_QA]->DrawCopy("PE");

      h_mSpec_QA[KEY_QA]->SetMarkerColor(4);
      h_mSpec_QA[KEY_QA]->SetMarkerStyle(24);
      h_mSpec_QA[KEY_QA]->SetMarkerSize(0.8);
      h_mSpec_QA[KEY_QA]->DrawCopy("PE same");

      TF1 *f_bw = new TF1("f_bw",PolyBreitWigner,BW_Start[mPID],BW_Stop[mPID],5);
      for(Int_t i_par = 0; i_par < 5; i_par++)
      {
	f_bw->SetParameter(i_par,ParFit[KEY_QA][i_par]);
      }
      f_bw->SetLineColor(2);
      f_bw->SetLineStyle(1);
      f_bw->SetLineWidth(2);
      f_bw->DrawCopy("l same");

      TF1 *f_poly = new TF1("f_poly",Poly,BW_Start[mPID],BW_Stop[mPID],2);
      f_poly->SetParameter(0,ParFit[KEY_QA][3]);
      f_poly->SetParameter(1,ParFit[KEY_QA][4]);
      f_poly->SetLineColor(2);
      f_poly->SetLineStyle(2);
      f_poly->SetLineWidth(4);
      f_poly->DrawCopy("l same");

      TString pT_range = Form("[%.2f,%.2f]",0.5*(pt_low_raw[i_pt]+pt_up_raw[i_pt]),pt_up_raw[i_pt]);
      plotTopLegend((char*)pT_range.Data(),0.2,0.7,0.08,1,0.0,42,1);
      PlotLine(0.98,1.05,0.0,0.0,1,2,2);
    }
#endif
  }

  // calculate total yields for each pT bin via gaussian and breit wigner fits
  vecFMap ParYield;
  vecFMap yields_count, yields_inte;
  for(Int_t i_pt = pt_start; i_pt < pt_stop; i_pt++) // pt loop
  {
    for(Int_t i_range = 0; i_range < 2; i_range++)
    {
      for(Int_t i_cent = Cent_start; i_cent < Cent_stop; i_cent++) // Centrality loop
      {
	for(Int_t i_eta = Eta_start; i_eta < Eta_stop; i_eta++) // EtaGap loop
	{
	  for(Int_t i_sys = Sys_start; i_sys < Sys_stop; i_sys++) // Systematic loop
	  {
	    TString KEY = Form("Spec_pt_%d_%s_Centrality_%d_EtaGap_%d_%s_SysErrors_%d",i_pt,pt_range[i_range].Data(),i_cent,i_eta,PID[mPID].Data(),i_sys);
	    TF1 *f_yields_bw = new TF1("f_yields_bw",BreitWigner,BW_Start[mPID],BW_Stop[mPID],3);
	    f_yields_bw->SetParameter(0,InvMass[mPID]);
	    f_yields_bw->SetParLimits(0,InvMass[mPID]-0.005,InvMass[mPID]+0.005);
	    f_yields_bw->SetParameter(1,Width[mPID]);
	    f_yields_bw->SetParameter(2,h_mSpec[KEY]->GetMaximum());
	    f_yields_bw->SetRange(BW_Start[mPID],BW_Stop[mPID]);
	    h_mSpec[KEY]->Fit(f_yields_bw,"MQNR");
	    ParYield[KEY].clear();
	    ParYield[KEY].push_back(static_cast<Float_t>(f_yields_bw->GetParameter(0)));
	    ParYield[KEY].push_back(static_cast<Float_t>(f_yields_bw->GetParameter(1)));
	    ParYield[KEY].push_back(static_cast<Float_t>(f_yields_bw->GetParameter(2)));

	    // bin counting
	    Float_t counts_gaus = 0.0;
	    Float_t errors_gaus = 0.0;
	    Int_t bin_start = h_mSpec[KEY]->FindBin(ParYield[KEY][0]-nSigV0*ParYield[KEY][1]);
	    Int_t bin_stop  = h_mSpec[KEY]->FindBin(ParYield[KEY][0]+nSigV0*ParYield[KEY][1]);
	    for(Int_t i_bin = bin_start; i_bin <= bin_stop; i_bin++)
	    {
	      counts_gaus += h_mSpec[KEY]->GetBinContent(i_bin);
	      errors_gaus += h_mSpec[KEY]->GetBinError(i_bin)*h_mSpec[KEY]->GetBinError(i_bin);
	    }
	    yields_count[KEY].clear();
	    yields_count[KEY].push_back(static_cast<Float_t>(counts_gaus));
	    yields_count[KEY].push_back(static_cast<Float_t>(TMath::Sqrt(errors_gaus)));

	    // integrating for breit wigner
	    Float_t bin_width = h_mSpec[KEY]->GetBinWidth(1);
	    Float_t Inte_start = ParYield[KEY][0]-nSigV0*ParYield[KEY][1]-0.5*bin_width;
	    Float_t Inte_stop  = ParYield[KEY][0]+nSigV0*ParYield[KEY][1]+0.5*bin_width;
	    Float_t counts_bw = f_yields_bw->Integral(Inte_start,Inte_stop)/bin_width;
	    Float_t errors_bw = f_yields_bw->IntegralError(Inte_start,Inte_stop)/bin_width;
	    yields_inte[KEY].clear();
	    yields_inte[KEY].push_back(static_cast<Float_t>(counts_bw));
	    yields_inte[KEY].push_back(static_cast<Float_t>(errors_bw));
	  }
	}
      }
    }
  }

#if _PlotQA_
  // QA: different counting method: bin counting vs breit wigner integrating
  TCanvas *c_Yields_low = new TCanvas("c_Yields_low","c_Yields_low",10,10,900,900);
  c_Yields_low->Divide(5,5);
  for(Int_t i_pt = pt_start; i_pt < pt_stop; i_pt++)
  {
    c_Yields_low->cd(i_pt+1);
    c_Yields_low->cd(i_pt+1)->SetLeftMargin(0.20);
    c_Yields_low->cd(i_pt+1)->SetBottomMargin(0.20);
    c_Yields_low->cd(i_pt+1)->SetTicks(1,1);
    c_Yields_low->cd(i_pt+1)->SetGrid(0,0);

    TString KEY_QA = Form("Spec_pt_%d_low_Centrality_%d_EtaGap_%d_%s_SysErrors_%d",i_pt,Cent_start,Eta_start,PID[mPID].Data(),Sys_start);
    h_mSpec[KEY_QA]->SetTitle("");
    h_mSpec[KEY_QA]->SetStats(0);
    h_mSpec[KEY_QA]->SetMarkerStyle(24);
    h_mSpec[KEY_QA]->SetMarkerColor(kGray+3);
    h_mSpec[KEY_QA]->SetMarkerSize(0.8);
    h_mSpec[KEY_QA]->GetXaxis()->SetNdivisions(505,'N');
    h_mSpec[KEY_QA]->GetXaxis()->SetTitle("InvMass (GeV/c^{2})");
    h_mSpec[KEY_QA]->GetXaxis()->CenterTitle();
    h_mSpec[KEY_QA]->GetYaxis()->SetTitle("Counts");
    h_mSpec[KEY_QA]->GetYaxis()->CenterTitle();
    h_mSpec[KEY_QA]->GetYaxis()->SetTitleOffset(1.2);
    h_mSpec[KEY_QA]->DrawCopy("pE");
    PlotLine(0.98,1.05,0.0,0.0,1,2,2);

    TF1 *f_yields_bw = new TF1("f_yields_bw",BreitWigner,BW_Start[mPID],BW_Stop[mPID],3);
    f_yields_bw->SetParameter(0,ParYield[KEY_QA][0]);
    f_yields_bw->SetParameter(1,ParYield[KEY_QA][1]);
    f_yields_bw->SetParameter(2,ParYield[KEY_QA][2]);
    f_yields_bw->SetRange(BW_Start[mPID],BW_Stop[mPID]);
    f_yields_bw->SetLineColor(2);
    f_yields_bw->SetLineStyle(1);
    f_yields_bw->SetLineWidth(2);
    f_yields_bw->DrawCopy("l same");

    Float_t x1 = ParYield[KEY_QA][0] - nSigV0*ParYield[KEY_QA][1];
    Float_t x2 = ParYield[KEY_QA][0] + nSigV0*ParYield[KEY_QA][1];
    Float_t y = h_mSpec[KEY_QA]->GetBinContent(h_mSpec[KEY_QA]->FindBin(ParYield[KEY_QA][0]));
    PlotLine(x1,x1,0,y,4,2,2);
    PlotLine(x2,x2,0,y,4,2,2);
  }

  TCanvas *c_Yields_high = new TCanvas("c_Yields_high","c_Yields_high",10,10,900,900);
  c_Yields_high->Divide(5,5);
  for(Int_t i_pt = pt_start; i_pt < pt_stop; i_pt++)
  {
    c_Yields_high->cd(i_pt+1);
    c_Yields_high->cd(i_pt+1)->SetLeftMargin(0.20);
    c_Yields_high->cd(i_pt+1)->SetBottomMargin(0.20);
    c_Yields_high->cd(i_pt+1)->SetTicks(1,1);
    c_Yields_high->cd(i_pt+1)->SetGrid(0,0);

    TString KEY_QA = Form("Spec_pt_%d_high_Centrality_%d_EtaGap_%d_%s_SysErrors_%d",i_pt,Cent_start,Eta_start,PID[mPID].Data(),Sys_start);
    h_mSpec[KEY_QA]->SetTitle("");
    h_mSpec[KEY_QA]->SetStats(0);
    h_mSpec[KEY_QA]->SetMarkerStyle(24);
    h_mSpec[KEY_QA]->SetMarkerColor(kGray+3);
    h_mSpec[KEY_QA]->SetMarkerSize(0.8);
    h_mSpec[KEY_QA]->GetXaxis()->SetNdivisions(505,'N');
    h_mSpec[KEY_QA]->GetXaxis()->SetTitle("InvMass (GeV/c^{2})");
    h_mSpec[KEY_QA]->GetXaxis()->CenterTitle();
    h_mSpec[KEY_QA]->GetYaxis()->SetTitle("Counts");
    h_mSpec[KEY_QA]->GetYaxis()->CenterTitle();
    h_mSpec[KEY_QA]->GetYaxis()->SetTitleOffset(1.2);
    h_mSpec[KEY_QA]->DrawCopy("pE");
    PlotLine(0.98,1.05,0.0,0.0,1,2,2);

    TF1 *f_yields_bw = new TF1("f_yields_bw",BreitWigner,BW_Start[mPID],BW_Stop[mPID],3);
    f_yields_bw->SetParameter(0,ParYield[KEY_QA][0]);
    f_yields_bw->SetParameter(1,ParYield[KEY_QA][1]);
    f_yields_bw->SetParameter(2,ParYield[KEY_QA][2]);
    f_yields_bw->SetRange(BW_Start[mPID],BW_Stop[mPID]);
    f_yields_bw->SetLineColor(2);
    f_yields_bw->SetLineStyle(1);
    f_yields_bw->SetLineWidth(2);
    f_yields_bw->DrawCopy("l same");

    Float_t x1 = ParYield[KEY_QA][0] - nSigV0*ParYield[KEY_QA][1];
    Float_t x2 = ParYield[KEY_QA][0] + nSigV0*ParYield[KEY_QA][1];
    Float_t y = h_mSpec[KEY_QA]->GetBinContent(h_mSpec[KEY_QA]->FindBin(ParYield[KEY_QA][0]));
    PlotLine(x1,x1,0,y,4,2,2);
    PlotLine(x2,x2,0,y,4,2,2);
  }
#endif

  // declare histogram with different pT width
  Float_t pt_low[51], pt_up[51], pt_width[51];
  for(Int_t i_pt = pt_start; i_pt < pt_stop; i_pt++)
  {
    pt_low[2*i_pt] = pt_low_raw[i_pt];
    pt_up[2*i_pt]  = 0.5*(pt_low_raw[i_pt]+pt_up_raw[i_pt]);
    pt_width[2*i_pt] = pt_up[2*i_pt]-pt_low[2*i_pt];
    
    pt_low[2*i_pt+1] = 0.5*(pt_low_raw[i_pt]+pt_up_raw[i_pt]);
    pt_up[2*i_pt+1]  = pt_up_raw[i_pt];
    pt_width[2*i_pt+1] = pt_up[2*i_pt+1]-pt_low[2*i_pt+1];
  }
  pt_low[50] = 8.0; // make sure pT = 8.0 in the histogram
  pt_up[50] = 10.0;
  pt_width[50] = pt_up[50]-pt_low[50];

  TH1FMap h_mPt;
  for(Int_t i_cent = Cent_start; i_cent < Cent_stop; i_cent++) // Centrality loop
  {
    for(Int_t i_eta = Eta_start; i_eta < Eta_stop; i_eta++) // EtaGap loop
    {
      for(Int_t i_sys = Sys_start; i_sys < Sys_stop; i_sys++) // Systematic loop
      {
	TString KEY_pT_counts = Form("Count_Spec_Centrality_%d_EtaGap_%d_%s_SysErrors_%d",i_cent,i_eta,PID[mPID].Data(),i_sys);
	h_mPt[KEY_pT_counts] = new TH1F(KEY_pT_counts.Data(),KEY_pT_counts.Data(),2*pt_total,pt_low);
	TString KEY_pT_inte = Form("Inte_Spec_Centrality_%d_EtaGap_%d_%s_SysErrors_%d",i_cent,i_eta,PID[mPID].Data(),i_sys);
	h_mPt[KEY_pT_inte] = new TH1F(KEY_pT_inte.Data(),KEY_pT_inte.Data(),2*pt_total,pt_low);
	for(Int_t i_pt = pt_start; i_pt < pt_stop; i_pt++) // pt loop
	{
	  for(Int_t i_range = 0; i_range < 2; i_range++)
	  {
	    TString KEY = Form("Spec_pt_%d_%s_Centrality_%d_EtaGap_%d_%s_SysErrors_%d",i_pt,pt_range[i_range].Data(),i_cent,i_eta,PID[mPID].Data(),i_sys);
	    h_mPt[KEY_pT_counts]->SetBinContent(2*i_pt+1+i_range,yields_count[KEY][0]/pt_width[2*i_pt+i_range]);
	    h_mPt[KEY_pT_counts]->SetBinError(2*i_pt+1+i_range,yields_count[KEY][1]/TMath::Sqrt(pt_width[2*i_pt+i_range]));
	    h_mPt[KEY_pT_inte]->SetBinContent(2*i_pt+1+i_range,yields_inte[KEY][0]/pt_width[2*i_pt+i_range]);
	    h_mPt[KEY_pT_inte]->SetBinError(2*i_pt+1+i_range,yields_inte[KEY][1]/TMath::Sqrt(pt_width[2*i_pt+i_range]));
	  }
	}
      }
    }
  }

#if _PlotQA_
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
#endif

  // Save h_mPt
  TString OutPutFile = Form("./OutPut/AuAu%s/%s/h_pt.root",Energy[mEnergy].Data(),PID[mPID].Data());
  TFile *File_OutPut = new TFile(OutPutFile.Data(),"RECREATE");
  File_OutPut->cd();
  for(Int_t i_cent = Cent_start; i_cent < Cent_stop; i_cent++) // Centrality loop
  {
    for(Int_t i_eta = Eta_start; i_eta < Eta_stop; i_eta++) // EtaGap loop
    {
      for(Int_t i_sys = Sys_start; i_sys < Sys_stop; i_sys++) // Systematic loop
      {
	TString KEY_pT_counts = Form("Count_Spec_Centrality_%d_EtaGap_%d_%s_SysErrors_%d",i_cent,i_eta,PID[mPID].Data(),i_sys);
	h_mPt[KEY_pT_counts]->SetMarkerStyle(24);
	h_mPt[KEY_pT_counts]->SetMarkerColor(4);
	h_mPt[KEY_pT_counts]->Write();
	TString KEY_pT_inte = Form("Inte_Spec_Centrality_%d_EtaGap_%d_%s_SysErrors_%d",i_cent,i_eta,PID[mPID].Data(),i_sys);
	h_mPt[KEY_pT_inte]->SetMarkerStyle(24);
	h_mPt[KEY_pT_inte]->SetMarkerColor(2);
	h_mPt[KEY_pT_inte]->Write();
      }
    }
  }
  File_OutPut->Close();
}

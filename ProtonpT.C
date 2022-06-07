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
#include "./student_t_1d_single.h"

// mEnergy: 0, 200 GeV | 1, 39 GeV
static const TString Energy[2] = {"200GeV","39GeV"};
static const TString Order[2] = {"2nd","3rd"};
static const Double_t PI_max[2] = {TMath::Pi()/2.0,TMath::Pi()/3.0};
static const Float_t nSigProton = 3.0;
static const Float_t Flow_Order[2] = {2.0,3.0};

static const Int_t pt_total = 16; // pT loop
static const Int_t pt_start = 0;
static const Int_t pt_stop  = 16;
static const TString pt_range[2] = {"low","high"};

// pt bin
//                                      0 , 1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 ,10 ,11 ,12 ,13 ,14 ,15
static const Float_t pt_low_raw[16] = {0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4};
static const Float_t pt_up_raw[16]  = {0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,4.2};
//static const Float_t m2_cut_low[pt_total] = {0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.65,0.7,0.7,0.75,0.879,0.879,0.879,0.879,0.879};
//static const Float_t m2_cut_up[pt_total]  = {1.2,1.2,1.2,1.2,1.2,1.2,1.2,1.20,1.2,1.2,1.20,1.400,1.400,1.400,1.400,1.400};
static const Float_t m2_cut_low[pt_total] = {0.879,0.879,0.879,0.879,0.879,0.879,0.879,0.879,0.879,0.879,0.879,0.879,0.879,0.879,0.879,0.879};
static const Float_t m2_cut_up[pt_total]  = {1.400,1.400,1.400,1.400,1.400,1.400,1.400,1.400,1.400,1.400,1.400,1.400,1.400,1.400,1.400,1.400};

// pt rebin
static const Int_t pt_QA    = 4;
static const Float_t m2_yield_low  = 0.6;
static const Float_t m2_yield_high = 1.2;

static const Int_t Cent_total = 4; // Centrality loop
static const Int_t Cent_start = 0;
static const Int_t Cent_stop  = 1;
static const Int_t cent_low[4] = {0,7,4,0}; // 0 = 0-80%, 1 = 0-10%, 2 = 10-40%, 3 = 40-80%
static const Int_t cent_up[4]  = {8,8,6,3}; // 0 = 0-80%, 1 = 0-10%, 2 = 10-40%, 3 = 40-80%

static const Int_t Charge_total = 2;
static const Int_t Charge_start = 0;
static const Int_t Charge_stop  = 2;

static const Int_t Eta_total = 4; // Eta loop
static const Int_t Eta_start = 0;
static const Int_t Eta_stop  = 1;

static const Int_t phi_total = 7; // phi loop
static const Int_t phi_start = 0;
static const Int_t phi_stop  = 7;

static const Int_t Sys_total = 18; // Systematic loop
static const Int_t Sys_start = 0;
static const Int_t Sys_stop  = 18;

typedef std::map<TString,TH1F*> TH1FMap;
typedef std::map<TString,TProfile*> TProMap;
typedef std::map<TString,TGraphAsymmErrors*> TGraMap;
typedef std::map<TString,std::vector<Float_t>> vecFMap;

void ProtonpT(Int_t mEnergy = 0, Int_t mOrder = 1)
{
  TGaxis::SetMaxDigits(4);

  TString InPutFile_Spec = Form("./Data/AuAu%s/Proton/Flow_%s.root",Energy[mEnergy].Data(),Energy[mEnergy].Data());
  TFile *File_Spec = TFile::Open(InPutFile_Spec.Data());

  // read in  histogram for flow calculation
  TH1FMap h_mSpec;
  for(Int_t i_pt = 0; i_pt < pt_total; i_pt++) // pt bin
  {
    for(Int_t i_range = 0; i_range < 2; i_range++)
    {
      for(Int_t i_cent = Cent_start; i_cent < Cent_stop; i_cent++) // centrality bin
      {
	for(Int_t i_charge = Charge_start; i_charge < Charge_stop; i_charge++) // charge bin
	{
	  for(Int_t i_eta = Eta_start; i_eta < Eta_stop; i_eta++) // eta gap bin
	  {
	    for(Int_t i_sys = 0; i_sys < 18; i_sys++)
	    {
	      TString KEY_Proton_Spec = Form("Spectra_pt_%d_%s_Centrality_%d_Charge_%d_EtaGap_%d_Proton_%s_SysError_%d",i_pt,pt_range[i_range].Data(),i_cent,i_charge,i_eta,Order[mOrder].Data(),i_sys);
	      h_mSpec[KEY_Proton_Spec] = (TH1F*)File_Spec->Get(KEY_Proton_Spec.Data());
	    }
	  }
	}
      }
    }
  }

  /*
  // QA: spec vs pT
  TCanvas *c_pT_low = new TCanvas("c_pT_low","c_pT_low",10,10,1400,1400);
  c_pT_low->Divide(4,4);
  for(Int_t i_pt = pt_start; i_pt < pt_stop; i_pt++) // pt loop
  {
    c_pT_low->cd(i_pt+1);
    c_pT_low->cd(i_pt+1)->SetLeftMargin(0.15);
    c_pT_low->cd(i_pt+1)->SetBottomMargin(0.15);
    c_pT_low->cd(i_pt+1)->SetTicks(1,1);
    c_pT_low->cd(i_pt+1)->SetGrid(0,0);

    TString KEY_Proton_Spec_QA = Form("Spectra_pt_%d_%s_Centrality_%d_Charge_%d_EtaGap_%d_Proton_%s_SysError_%d",i_pt,pt_range[0].Data(),Cent_start,Charge_start,Eta_start,Order[mOrder].Data(),Sys_start);
    h_mSpec[KEY_Proton_Spec_QA]->SetTitle("");
    h_mSpec[KEY_Proton_Spec_QA]->SetStats(0);
    h_mSpec[KEY_Proton_Spec_QA]->GetXaxis()->SetTitle("m^{2} ((GeV/c^{2})^{2})");
    h_mSpec[KEY_Proton_Spec_QA]->GetXaxis()->CenterTitle();
    h_mSpec[KEY_Proton_Spec_QA]->GetXaxis()->SetTitleSize(0.06);
    h_mSpec[KEY_Proton_Spec_QA]->GetYaxis()->SetTitle("Counts");
    h_mSpec[KEY_Proton_Spec_QA]->GetYaxis()->CenterTitle();
    h_mSpec[KEY_Proton_Spec_QA]->GetYaxis()->SetTitleSize(0.06);
    h_mSpec[KEY_Proton_Spec_QA]->DrawCopy("PE");

    TString pT_range = Form("[%.2f,%.2f]",pt_low_raw[i_pt],0.5*(pt_low_raw[i_pt]+pt_up_raw[i_pt]));
    plotTopLegend((char*)pT_range.Data(),0.61,0.8,0.08,1,0.0,42,1);
    PlotLine(m2_cut_low[i_pt],m2_cut_low[i_pt],0.0,h_mSpec[KEY_Proton_Spec_QA]->GetMaximum()/2.0,4,2,2);
    PlotLine( m2_cut_up[i_pt], m2_cut_up[i_pt],0.0,h_mSpec[KEY_Proton_Spec_QA]->GetMaximum()/2.0,4,2,2);
  }

  TCanvas *c_pT_high = new TCanvas("c_pT_high","c_pT_high",10,10,1400,1400);
  c_pT_high->Divide(4,4);
  for(Int_t i_pt = pt_start; i_pt < pt_stop; i_pt++) // pt loop
  {
    c_pT_high->cd(i_pt+1);
    c_pT_high->cd(i_pt+1)->SetLeftMargin(0.15);
    c_pT_high->cd(i_pt+1)->SetBottomMargin(0.15);
    c_pT_high->cd(i_pt+1)->SetTicks(1,1);
    c_pT_high->cd(i_pt+1)->SetGrid(0,0);

    TString KEY_Proton_Spec_QA = Form("Spectra_pt_%d_%s_Centrality_%d_Charge_%d_EtaGap_%d_Proton_%s_SysError_%d",i_pt,pt_range[1].Data(),Cent_start,Charge_start,Eta_start,Order[mOrder].Data(),Sys_start);
    h_mSpec[KEY_Proton_Spec_QA]->SetTitle("");
    h_mSpec[KEY_Proton_Spec_QA]->SetStats(0);
    h_mSpec[KEY_Proton_Spec_QA]->GetXaxis()->SetTitle("m^{2} ((GeV/c^{2})^{2})");
    h_mSpec[KEY_Proton_Spec_QA]->GetXaxis()->CenterTitle();
    h_mSpec[KEY_Proton_Spec_QA]->GetXaxis()->SetTitleSize(0.06);
    h_mSpec[KEY_Proton_Spec_QA]->GetYaxis()->SetTitle("Counts");
    h_mSpec[KEY_Proton_Spec_QA]->GetYaxis()->CenterTitle();
    h_mSpec[KEY_Proton_Spec_QA]->GetYaxis()->SetTitleSize(0.06);
    h_mSpec[KEY_Proton_Spec_QA]->DrawCopy("PE");

    TString pT_range = Form("[%.2f,%.2f]",pt_low_raw[i_pt],0.5*(pt_low_raw[i_pt]+pt_up_raw[i_pt]));
    plotTopLegend((char*)pT_range.Data(),0.61,0.8,0.08,1,0.0,42,1);
    PlotLine(m2_cut_low[i_pt],m2_cut_low[i_pt],0.0,h_mSpec[KEY_Proton_Spec_QA]->GetMaximum()/2.0,4,2,2);
    PlotLine( m2_cut_up[i_pt], m2_cut_up[i_pt],0.0,h_mSpec[KEY_Proton_Spec_QA]->GetMaximum()/2.0,4,2,2);
  }
  */

  // get counts for all pT bin
  vecFMap yields;
  for(Int_t i_pt = 0; i_pt < pt_total; i_pt++) // pt bin
  {
    for(Int_t i_range = 0; i_range < 2; i_range++)
    {
      for(Int_t i_cent = Cent_start; i_cent < Cent_stop; i_cent++) // centrality bin
      {
	for(Int_t i_charge = Charge_start; i_charge < Charge_stop; i_charge++) // charge bin
	{
	  for(Int_t i_eta = Eta_start; i_eta < Eta_stop; i_eta++) // eta gap bin
	  {
	    for(Int_t i_sys = 0; i_sys < 18; i_sys++)
	    {
	      TString KEY_Proton_Spec = Form("Spectra_pt_%d_%s_Centrality_%d_Charge_%d_EtaGap_%d_Proton_%s_SysError_%d",i_pt,pt_range[i_range].Data(),i_cent,i_charge,i_eta,Order[mOrder].Data(),i_sys);
	      Int_t bin_start = h_mSpec[KEY_Proton_Spec]->FindBin(m2_cut_low[i_pt]);
	      Int_t bin_stop  = h_mSpec[KEY_Proton_Spec]->FindBin(m2_cut_up[i_pt]);
	      Float_t counts = 0.0;
	      Float_t errors = 0.0;
	      for(Int_t i_bin = bin_start; i_bin <= bin_stop; i_bin++)
	      {
		counts += h_mSpec[KEY_Proton_Spec]->GetBinContent(i_bin);
		errors += h_mSpec[KEY_Proton_Spec]->GetBinError(i_bin)*h_mSpec[KEY_Proton_Spec]->GetBinError(i_bin);
	      }
	      yields[KEY_Proton_Spec].clear();
	      yields[KEY_Proton_Spec].push_back(static_cast<Float_t>(counts));
	      yields[KEY_Proton_Spec].push_back(static_cast<Float_t>(TMath::Sqrt(errors)));
	    }
	  }
	}
      }
    }
  }

  // declare histogram with different pT width
  Float_t pt_low[33], pt_up[33], pt_width[33];
  for(Int_t i_pt = pt_start; i_pt < pt_stop; i_pt++)
  {
    pt_low[2*i_pt] = pt_low_raw[i_pt];
    pt_up[2*i_pt]  = 0.5*(pt_low_raw[i_pt]+pt_up_raw[i_pt]);
    pt_width[2*i_pt] = pt_up[2*i_pt]-pt_low[2*i_pt];
    
    pt_low[2*i_pt+1] = 0.5*(pt_low_raw[i_pt]+pt_up_raw[i_pt]);
    pt_up[2*i_pt+1]  = pt_up_raw[i_pt];
    pt_width[2*i_pt+1] = pt_up[2*i_pt+1]-pt_low[2*i_pt+1];
  }
  pt_low[32] = 4.2; // make sure pT = 4.2 in the histogram
  pt_up[32] = 5.0;
  pt_width[32] = pt_up[32]-pt_low[32];

  TH1FMap h_mPt;
  for(Int_t i_cent = Cent_start; i_cent < Cent_stop; i_cent++) // centrality bin
  {
    for(Int_t i_charge = Charge_start; i_charge < Charge_stop; i_charge++) // charge bin
    {
      for(Int_t i_eta = Eta_start; i_eta < Eta_stop; i_eta++) // eta gap bin
      {
	for(Int_t i_sys = 0; i_sys < 18; i_sys++)
	{
	  TString KEY_Proton_pT = Form("Spec_Centrality_%d_Charge_%d_EtaGap_%d_Proton_%s_SysErrors_%d",i_cent,i_charge,i_eta,Order[mOrder].Data(),i_sys);
	  h_mPt[KEY_Proton_pT] = new TH1F(KEY_Proton_pT.Data(),KEY_Proton_pT.Data(),2*pt_total,pt_low);
	  for(Int_t i_pt = 0; i_pt < pt_total; i_pt++) // pt bin
	  {
	    for(Int_t i_range = 0; i_range < 2; i_range++)
	    {
	      TString KEY_Proton_Spec = Form("Spectra_pt_%d_%s_Centrality_%d_Charge_%d_EtaGap_%d_Proton_%s_SysError_%d",i_pt,pt_range[i_range].Data(),i_cent,i_charge,i_eta,Order[mOrder].Data(),i_sys);
	      h_mPt[KEY_Proton_pT]->SetBinContent(2*i_pt+1+i_range,yields[KEY_Proton_Spec][0]/pt_width[2*i_pt+i_range]);
	      h_mPt[KEY_Proton_pT]->SetBinError(2*i_pt+1+i_range,yields[KEY_Proton_Spec][1]/TMath::Sqrt(pt_width[2*i_pt+i_range]));
	    }
	  }
	}
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
  TString KEY_Proton_pT_QA = Form("Spec_Centrality_%d_Charge_%d_EtaGap_%d_Proton_%s_SysErrors_%d",Cent_start,Charge_start,Eta_start,Order[mOrder].Data(),Sys_start);
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
  h_play->GetYaxis()->SetRangeUser(0.1,2.0*h_mPt[KEY_Proton_pT_QA]->GetMaximum());
  h_play->GetXaxis()->SetRangeUser(-0.1,4.2);
  h_play->DrawCopy("pE");
  h_mPt[KEY_Proton_pT_QA]->SetMarkerStyle(24);
  h_mPt[KEY_Proton_pT_QA]->SetMarkerColor(4);
  h_mPt[KEY_Proton_pT_QA]->SetMarkerSize(1.0);
  h_mPt[KEY_Proton_pT_QA]->DrawCopy("pE same");
  */

  // Save h_mPt
  TString OutPutFile = Form("./OutPut/AuAu%s/Proton/h_pt_%s.root",Energy[mEnergy].Data(),Order[mOrder].Data());
  TFile *File_OutPut = new TFile(OutPutFile.Data(),"RECREATE");
  File_OutPut->cd();
  for(Int_t i_cent = Cent_start; i_cent < Cent_stop; i_cent++) // centrality bin
  {
    for(Int_t i_charge = Charge_start; i_charge < Charge_stop; i_charge++) // charge bin
    {
      for(Int_t i_eta = Eta_start; i_eta < Eta_stop; i_eta++) // eta gap bin
      {
	for(Int_t i_sys = 0; i_sys < 18; i_sys++)
	{
	  TString KEY_Proton_pT = Form("Spec_Centrality_%d_Charge_%d_EtaGap_%d_Proton_%s_SysErrors_%d",i_cent,i_charge,i_eta,Order[mOrder].Data(),i_sys);
	  h_mPt[KEY_Proton_pT]->SetMarkerStyle(24);
	  h_mPt[KEY_Proton_pT]->SetMarkerColor(4);
	  h_mPt[KEY_Proton_pT]->SetMarkerSize(1.0);
	  h_mPt[KEY_Proton_pT]->Write();
	}
      }
    }
  }
  File_OutPut->Close();
  File_Spec->Close();
}

#include <iostream>
#include "TString.h"
#include "TFile.h"
#include "TH1F.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TLine.h"

Double_t PolyBreitWigner(Double_t *x_val, Double_t *par) 
{
  Double_t x = x_val[0];
  Double_t m0 = par[0];
  Double_t Gamma = par[1];
  Double_t Norm = par[2];

  Double_t denom = 2.0*TMath::Pi()*((x-m0)*(x-m0)+Gamma*Gamma/4.0);
  Double_t BW = Norm*Gamma/denom;

  Double_t Poly = par[3] + par[4]*x;

  Double_t y = BW + Poly;

  return y;
}

Double_t Poly(Double_t *x_val, Double_t *par) 
{
  Double_t x = x_val[0];
  Double_t y = par[0] + par[1]*x;

  return y;
}

Double_t BreitWigner(Double_t *x_val, Double_t *par) 
{
  Double_t x = x_val[0];
  Double_t m0 = par[0];
  Double_t Gamma = par[1];
  Double_t Norm = par[2];

  Double_t denom = 2.0*TMath::Pi()*((x-m0)*(x-m0)+Gamma*Gamma/4.0);
  Double_t BW = Norm*Gamma/denom;

  return BW;
}

Double_t flow_2(Double_t *x_val, Double_t *par)
{
  Double_t x, y;
  Double_t Ampl, v2;
  x = x_val[0];
  Ampl = par[0];
  v2 = par[1];

  y = Ampl*(1.0 + 2.0*v2*TMath::Cos(2.0*x));

  return y;
}

Double_t flow_3(Double_t *x_val, Double_t *par)
{
  Double_t x, y;
  Double_t Ampl, v3;
  x = x_val[0];
  Ampl = par[0];
  v3 = par[1];

  y = Ampl*(1.0 + 2.0*v3*TMath::Cos(3.0*x));

  return y;
}

Double_t Gaussion(Double_t *x_val, Double_t *par)
{
  Double_t x, y;
  x = x_val[0];

  Double_t mu = par[0];
  Double_t sigma = par[1];
  Double_t norm = par[2];

  y = norm*TMath::Exp(-1.0*(x-mu)*(x-mu)/(2.0*sigma*sigma))/(sigma*TMath::Sqrt(2*TMath::Pi()));

  return y;
}

void PlotLine(Double_t x1_val, Double_t x2_val, Double_t y1_val, Double_t y2_val, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle)
{
  TLine* Zero_line = new TLine();
  Zero_line -> SetX1(x1_val);
  Zero_line -> SetX2(x2_val);
  Zero_line -> SetY1(y1_val);
  Zero_line -> SetY2(y2_val);
  Zero_line -> SetLineWidth(LineWidth);
  Zero_line -> SetLineStyle(LineStyle);
  Zero_line -> SetLineColor(Line_Col);
  Zero_line -> Draw();
  //delete Zero_line;
}

static TString mBeamEnergy[7] = {"7GeV","11GeV","19GeV","27GeV","39GeV","62GeV","200GeV"};
static TString mMode_AMPT[2] = {"Default","StringMelting"};
static TString mMode_SM[2] = {"SE","ME"};

static Int_t mList_start[20] = {  1,101,201,301,401,501,601,701,801, 901,1001,1101,1201,1301,1401,1501,1601,1701,1801,1901};
static Int_t mList_stop[20]  = {100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000};

// pt bin
//                                 0 , 1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 ,10 ,11 ,12 ,13 ,14 ,15 ,16 ,17 ,18 ,19 ,10 ,21 ,22
static Double_t pt_low_phi[23] = {0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.4,5.8,6.2};
static Double_t pt_up_phi[23]  = {0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.4,5.8,6.2,6.6};

// phi-Psi bin
static Double_t phi_Psi2_low[7] = {0.0,TMath::Pi()/14.0,2.0*TMath::Pi()/14.0,3.0*TMath::Pi()/14.0,4.0*TMath::Pi()/14.0,5.0*TMath::Pi()/14.0,6.0*TMath::Pi()/14.0};
static Double_t phi_Psi2_up[7]  = {TMath::Pi()/14.0,2.0*TMath::Pi()/14.0,3.0*TMath::Pi()/14.0,4.0*TMath::Pi()/14.0,5.0*TMath::Pi()/14.0,6.0*TMath::Pi()/14.0,7.0*TMath::Pi()/14.0};
static Double_t phi_Psi3_low[7] = {0.0,TMath::Pi()/21.0,2.0*TMath::Pi()/21.0,3.0*TMath::Pi()/21.0,4.0*TMath::Pi()/21.0,5.0*TMath::Pi()/21.0,6.0*TMath::Pi()/21.0};
static Double_t phi_Psi3_up[7]  = {TMath::Pi()/21.0,2.0*TMath::Pi()/21.0,3.0*TMath::Pi()/21.0,4.0*TMath::Pi()/21.0,5.0*TMath::Pi()/21.0,6.0*TMath::Pi()/21.0,7.0*TMath::Pi()/21.0};


static TString mOrder[2] = {"2nd","3rd"};
static TString mCentrality[4] = {"0080","0010","1040","4080"};

static const Int_t pt_total_phi = 23;
static const Int_t Centrality_total = 4;
static const Int_t Centrality_start = 2;
static const Int_t Centrality_stop = 3;
static const Int_t Phi_Psi_total = 7;
static const Float_t nSigmaPhi = 2.0;
static const Double_t PI_max[2] = {TMath::Pi()/2.0,TMath::Pi()/3.0};
static const Int_t cent_low[4] = {0,7,4,0}; // 0 = 0-80%, 1 = 0-10%, 2 = 10-40%, 3 = 40-80%
static const Int_t cent_up[4]  = {8,8,6,3}; // 0 = 0-80%, 1 = 0-10%, 2 = 10-40%, 3 = 40-80%
static const Float_t BW_Start = 0.994;
static const Float_t BW_Stop  = 1.050;

// new pt bins => need to adjust for different energies
static const Int_t pt_total_New_phi = 23;
static Int_t pt_new_bin_start[pt_total_New_phi] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22};
static Int_t pt_new_bin_stop[pt_total_New_phi]  = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22};

// Energy = 0: 7GeV, 1: 11GeV, 2: 19GeV, 3: 27GeV, 4: 39GeV, 5: 62GeV | Mode = 0: Default, 1: String Melting
void PhiFlow(Int_t mEnergy = 4, Int_t mMode = 0) 
{
  TString inputfile_SE = Form("/home/xusun/Data/AMPT_%s/Flow/%s_%s/Phi/Flow_%s_SE.root",mMode_AMPT[mMode].Data(),mBeamEnergy[mEnergy].Data(),mMode_AMPT[mMode].Data(),mBeamEnergy[mEnergy].Data());
  TFile *File_input_SE = TFile::Open(inputfile_SE.Data());

  TString inputfile_ME = Form("/home/xusun/Data/AMPT_%s/Flow/%s_%s/Phi/Flow_%s_ME.root",mMode_AMPT[mMode].Data(),mBeamEnergy[mEnergy].Data(),mMode_AMPT[mMode].Data(),mBeamEnergy[mEnergy].Data());
  TFile *File_input_ME = TFile::Open(inputfile_ME.Data());

  // read in TH1F for flow calculation
  TH1F *h_flow_SE[2][Centrality_total][pt_total_phi][Phi_Psi_total];
  TH1F *h_flow_ME[2][Centrality_total][pt_total_phi][Phi_Psi_total];
  TH1F *h_flow_SM[2][Centrality_total][pt_total_phi][Phi_Psi_total];
  for(Int_t i_order = 0; i_order < 2; i_order++) // harmonic loop
  {
    for(Int_t i_cent = Centrality_start; i_cent < Centrality_stop; i_cent++) // centrality loop
    {
      for(Int_t i_pt = 0; i_pt < pt_total_phi; i_pt++) // pt loop
      {
	for(Int_t i_phi = 0; i_phi < Phi_Psi_total; i_phi++) // phi_Psi loop
	{
	  TString HistName = Form("Flow_phi_%s_%s_pt_%d_phi_Psi_%d",mOrder[i_order].Data(),mCentrality[i_cent].Data(),i_pt,i_phi);
	  h_flow_SE[i_order][i_cent][i_pt][i_phi] = (TH1F*)File_input_SE->Get(HistName.Data());
	  h_flow_ME[i_order][i_cent][i_pt][i_phi] = (TH1F*)File_input_ME->Get(HistName.Data());

	  // initialize Histogram for SE_ME
	  HistName = Form("Flow_phi_%s_%s_pt_%d_phi_Psi_%d_SM",mOrder[i_order].Data(),mCentrality[i_cent].Data(),i_pt,i_phi);
	  h_flow_SM[i_order][i_cent][i_pt][i_phi] = (TH1F*)h_flow_SE[i_order][i_cent][i_pt][i_phi]->Clone(HistName.Data());

	  // normalize the SE and ME
	  Int_t Bin_SE_Start = h_flow_SE[i_order][i_cent][i_pt][i_phi]->FindBin(1.04);
	  Int_t Bin_SE_Stop  = h_flow_SE[i_order][i_cent][i_pt][i_phi]->FindBin(1.05);
	  Float_t Inte_SE = h_flow_SE[i_order][i_cent][i_pt][i_phi]->Integral(Bin_SE_Start,Bin_SE_Stop-1);

	  Int_t Bin_ME_Start = h_flow_ME[i_order][i_cent][i_pt][i_phi]->FindBin(1.04);
	  Int_t Bin_ME_Stop  = h_flow_ME[i_order][i_cent][i_pt][i_phi]->FindBin(1.05);
	  Float_t Inte_ME = h_flow_ME[i_order][i_cent][i_pt][i_phi]->Integral(Bin_ME_Start,Bin_ME_Stop-1);
	  if(Inte_ME > 0)
	  {
	    h_flow_ME[i_order][i_cent][i_pt][i_phi]->Scale(Inte_SE/Inte_ME);
	    h_flow_SM[i_order][i_cent][i_pt][i_phi]->Add(h_flow_ME[i_order][i_cent][i_pt][i_phi],-1.0);
	  }
	}
      }
    }
  }

  /*
  // QA plots: pt bins
  TCanvas *c_pt[2][Centrality_total][Phi_Psi_total];
  for(Int_t i_order = 0; i_order < 2; i_order++)
  {
    for(Int_t i_cent = Centrality_start; i_cent < Centrality_stop; i_cent++)
    {
      for(Int_t i_phi = 0; i_phi < Phi_Psi_total; i_phi++)
      {
	TString CanName = Form("c_pt_%s_%s_phi_Psi_%d",mOrder[i_order].Data(),mCentrality[i_cent].Data(),i_phi);
	c_pt[i_order][i_cent][i_phi] = new TCanvas(CanName.Data(),CanName.Data(),10,10,1000,1000);
	c_pt[i_order][i_cent][i_phi]->Divide(5,5);
	for(Int_t i_pt = 0; i_pt < pt_total_New_phi; i_pt++)
	{
	  c_pt[i_order][i_cent][i_phi]->cd(pt_new_bin_start[i_pt]+1);
	  c_pt[i_order][i_cent][i_phi]->cd(pt_new_bin_start[i_pt]+1)->SetLeftMargin(0.15);
	  c_pt[i_order][i_cent][i_phi]->cd(pt_new_bin_start[i_pt]+1)->SetBottomMargin(0.15);
	  c_pt[i_order][i_cent][i_phi]->cd(pt_new_bin_start[i_pt]+1)->SetTicks(1,1);
	  c_pt[i_order][i_cent][i_phi]->cd(pt_new_bin_start[i_pt]+1)->SetGrid(0,0);
	  TString Title = Form("%1.1f < p_{T} < %1.1f",pt_low_phi[i_pt],pt_up_phi[i_pt]);
	  h_flow_SE[i_order][i_cent][i_pt][i_phi]->SetTitle(Title.Data());
	  h_flow_SE[i_order][i_cent][i_pt][i_phi]->SetStats(0);
	  h_flow_SE[i_order][i_cent][i_pt][i_phi]->GetXaxis()->SetNdivisions(505,'X');
	  h_flow_SE[i_order][i_cent][i_pt][i_phi]->GetYaxis()->SetNdivisions(505,'Y');
	  h_flow_SE[i_order][i_cent][i_pt][i_phi]->GetXaxis()->SetTitle("M(K^{+},K^{-}) (GeV/c^{2})");
	  h_flow_SE[i_order][i_cent][i_pt][i_phi]->GetYaxis()->SetTitle("Counts/Resolution");
	  h_flow_SE[i_order][i_cent][i_pt][i_phi]->GetXaxis()->SetTitleSize(0.04);
	  h_flow_SE[i_order][i_cent][i_pt][i_phi]->GetYaxis()->SetTitleSize(0.04);
	  h_flow_SE[i_order][i_cent][i_pt][i_phi]->GetXaxis()->CenterTitle();
	  h_flow_SE[i_order][i_cent][i_pt][i_phi]->GetYaxis()->CenterTitle();
	  h_flow_SE[i_order][i_cent][i_pt][i_phi]->SetLineColor(1);
	  h_flow_SE[i_order][i_cent][i_pt][i_phi]->SetMarkerStyle(24);
	  h_flow_SE[i_order][i_cent][i_pt][i_phi]->SetMarkerColor(1);
	  h_flow_SE[i_order][i_cent][i_pt][i_phi]->SetMarkerSize(0.8);
	  h_flow_SE[i_order][i_cent][i_pt][i_phi]->Draw("pE");
	  h_flow_ME[i_order][i_cent][i_pt][i_phi]->SetLineColor(2);
	  h_flow_ME[i_order][i_cent][i_pt][i_phi]->SetFillColor(2);
	  h_flow_ME[i_order][i_cent][i_pt][i_phi]->SetFillStyle(3002);
	  h_flow_ME[i_order][i_cent][i_pt][i_phi]->Draw("h same");
	  h_flow_SM[i_order][i_cent][i_pt][i_phi]->SetLineColor(4);
	  h_flow_SM[i_order][i_cent][i_pt][i_phi]->SetFillColor(4);
	  h_flow_SM[i_order][i_cent][i_pt][i_phi]->SetFillStyle(3003);
	  h_flow_SM[i_order][i_cent][i_pt][i_phi]->Draw("h same");
	}
	CanName = "./figures/" + CanName + ".eps";
	c_pt[i_order][i_cent][i_phi]->SaveAs(CanName.Data());
      }
    }
  }
  */

  // merge original pT bins to new pT bins
  TH1F *h_flow_SM_New[2][Centrality_total][pt_total_New_phi][Phi_Psi_total];
  for(Int_t i_order = 0; i_order < 2; i_order++)
  {
    for(Int_t i_cent = Centrality_start; i_cent < Centrality_stop; i_cent++) // centrality bin
    {
      for(Int_t i_phi = 0; i_phi < Phi_Psi_total; i_phi++) // phi-psi bin
      {
	for(Int_t i_pt = 0; i_pt < pt_total_New_phi; i_pt++)
	{
	  for(Int_t pt_bin = pt_new_bin_start[i_pt]; pt_bin <= pt_new_bin_stop[i_pt]; pt_bin++)
	  {
	    if(pt_bin == pt_new_bin_start[i_pt])
	    {
	      h_flow_SM_New[i_order][i_cent][i_pt][i_phi] = (TH1F*)h_flow_SM[i_order][i_cent][pt_bin][i_phi]->Clone();
	    }
	    else
	    {
	      h_flow_SM_New[i_order][i_cent][i_pt][i_phi]->Add(h_flow_SM[i_order][i_cent][pt_bin][i_phi],1.0);
	    }
	  }
	}
      }
    }
  }

  /*
  // QA plots: new pt bins
  TCanvas *c_pt_new[2][Centrality_total][Phi_Psi_total];
  for(Int_t i_order = 0; i_order < 2; i_order++)
  {
    for(Int_t i_cent = Centrality_start; i_cent < Centrality_stop; i_cent++)
    {
      for(Int_t i_phi = 0; i_phi < Phi_Psi_total; i_phi++)
      {
	TString CanName = Form("c_pt_new_%s_%s_phi_Psi_%d",mOrder[i_order].Data(),mCentrality[i_cent].Data(),i_phi);
	c_pt_new[i_order][i_cent][i_phi] = new TCanvas(CanName.Data(),CanName.Data(),10,10,1000,1000);
	c_pt_new[i_order][i_cent][i_phi]->Divide(5,5);
	for(Int_t i_pt = 0; i_pt < pt_total_New_phi; i_pt++)
	{
	  c_pt_new[i_order][i_cent][i_phi]->cd(pt_new_bin_start[i_pt]+1);
	  c_pt_new[i_order][i_cent][i_phi]->cd(pt_new_bin_start[i_pt]+1)->SetLeftMargin(0.15);
	  c_pt_new[i_order][i_cent][i_phi]->cd(pt_new_bin_start[i_pt]+1)->SetBottomMargin(0.15);
	  c_pt_new[i_order][i_cent][i_phi]->cd(pt_new_bin_start[i_pt]+1)->SetTicks(1,1);
	  c_pt_new[i_order][i_cent][i_phi]->cd(pt_new_bin_start[i_pt]+1)->SetGrid(0,0);
	  TString Title = Form("%1.1f < p_{T} < %1.1f",pt_low_phi[i_pt],pt_up_phi[i_pt]);
	  h_flow_SM_New[i_order][i_cent][i_pt][i_phi]->SetTitle(Title.Data());
	  h_flow_SM_New[i_order][i_cent][i_pt][i_phi]->SetStats(0);
	  h_flow_SM_New[i_order][i_cent][i_pt][i_phi]->GetXaxis()->SetNdivisions(505,'X');
	  h_flow_SM_New[i_order][i_cent][i_pt][i_phi]->GetYaxis()->SetNdivisions(505,'Y');
	  h_flow_SM_New[i_order][i_cent][i_pt][i_phi]->GetXaxis()->SetTitle("M(K^{+},K^{-}) (GeV/c^{2})");
	  h_flow_SM_New[i_order][i_cent][i_pt][i_phi]->GetYaxis()->SetTitle("Counts/Resolution");
	  h_flow_SM_New[i_order][i_cent][i_pt][i_phi]->GetXaxis()->SetTitleSize(0.04);
	  h_flow_SM_New[i_order][i_cent][i_pt][i_phi]->GetYaxis()->SetTitleSize(0.04);
	  h_flow_SM_New[i_order][i_cent][i_pt][i_phi]->GetXaxis()->CenterTitle();
	  h_flow_SM_New[i_order][i_cent][i_pt][i_phi]->GetYaxis()->CenterTitle();
	  h_flow_SM_New[i_order][i_cent][i_pt][i_phi]->SetMarkerStyle(24);
	  h_flow_SM_New[i_order][i_cent][i_pt][i_phi]->SetMarkerColor(2);
	  h_flow_SM_New[i_order][i_cent][i_pt][i_phi]->SetMarkerSize(0.8);
	  h_flow_SM_New[i_order][i_cent][i_pt][i_phi]->Draw("pE");
	}
	CanName = "./figures/" + CanName + ".eps";
	c_pt_new[i_order][i_cent][i_phi]->SaveAs(CanName.Data());
      }
    }
  }
  */

#if 0
  // Poly+BW fit for merged phi-Psi bin to extract fit parameters
  TH1F *h_flow_SM_total[2][Centrality_total][pt_total_New_phi];
  TF1  *f_PolyBW_total[2][Centrality_total][pt_total_New_phi];
  Float_t ParFit_total[2][Centrality_total][pt_total_New_phi][5];

  for(Int_t i_order = 0; i_order < 2; i_order++)
  {
    for(Int_t i_pt = 0; i_pt < pt_total_New_phi; i_pt++)
    {
      for(Int_t i_cent = Centrality_start; i_cent < Centrality_stop; i_cent++) // centrality bin
      {
	for(Int_t i_phi = 0; i_phi < Phi_Psi_total; i_phi++) // phi-psi bin => merge
	{
	  TString HistName;
	  if(i_phi == 0)
	  {
	    HistName = Form("pt_%d_%s_Total_%s_Phi_SM",i_pt,mCentrality[i_cent].Data(),mOrder[i_order].Data());
	    h_flow_SM_total[i_order][i_cent][i_pt] = (TH1F*)h_flow_SM_New[i_order][i_cent][i_pt][i_phi]->Clone(HistName.Data());
	    HistName = "f_"+HistName;
	    f_PolyBW_total[i_order][i_cent][i_pt] = new TF1(HistName.Data(),PolyBreitWigner,BW_Start,BW_Stop,5);
	    for(Int_t n_par = 0; n_par < 5; n_par++)
	    {
	      f_PolyBW_total[i_order][i_cent][i_pt]->ReleaseParameter(n_par);
	    }
	    f_PolyBW_total[i_order][i_cent][i_pt]->SetParameter(0,1.019);
	    f_PolyBW_total[i_order][i_cent][i_pt]->SetParLimits(0,1.014,1.024);
	    f_PolyBW_total[i_order][i_cent][i_pt]->SetParameter(1,0.0055);
	    f_PolyBW_total[i_order][i_cent][i_pt]->SetParameter(2,10000);
	    f_PolyBW_total[i_order][i_cent][i_pt]->SetParameter(3,-6000);
	    f_PolyBW_total[i_order][i_cent][i_pt]->SetParameter(4,0.5);
	    f_PolyBW_total[i_order][i_cent][i_pt]->SetRange(BW_Start,BW_Stop);
	  }
	  else
	  {
	    h_flow_SM_total[i_order][i_cent][i_pt]->Add(h_flow_SM_New[i_order][i_cent][i_pt][i_phi],1.0);
	  }
	}
	h_flow_SM_total[i_order][i_cent][i_pt]->Fit(f_PolyBW_total[i_order][i_cent][i_pt],"NQR");
	for(Int_t n_par = 0; n_par < 5; n_par++)
	{
	  ParFit_total[i_order][i_cent][i_pt][n_par] = f_PolyBW_total[i_order][i_cent][i_pt]->GetParameter(n_par);
	}
      }
    }
  }

  /*
  // QA plots: Poly+BW fit to extract fit parameters 
  TCanvas *c_PolyBW_total[2][Centrality_total];
  TF1 *f_Poly_total[2][Centrality_total][pt_total_New_phi];
  for(Int_t i_order = 0; i_order < 2; i_order++)
  {
    for(Int_t i_cent = Centrality_start; i_cent < Centrality_stop; i_cent++)
    {
      TString CanName = Form("c_PolyBW_total_%s_%s",mOrder[i_order].Data(),mCentrality[i_cent].Data());
      c_PolyBW_total[i_order][i_cent] = new TCanvas(CanName.Data(),CanName.Data(),10,10,1000,1000);
      c_PolyBW_total[i_order][i_cent]->Divide(5,5);
      for(Int_t i_pt = 0; i_pt < pt_total_New_phi; i_pt++)
      {
	c_PolyBW_total[i_order][i_cent]->cd(pt_new_bin_start[i_pt]+1);
	c_PolyBW_total[i_order][i_cent]->cd(pt_new_bin_start[i_pt]+1)->SetLeftMargin(0.15);
	c_PolyBW_total[i_order][i_cent]->cd(pt_new_bin_start[i_pt]+1)->SetBottomMargin(0.15);
	c_PolyBW_total[i_order][i_cent]->cd(pt_new_bin_start[i_pt]+1)->SetTicks(1,1);
	c_PolyBW_total[i_order][i_cent]->cd(pt_new_bin_start[i_pt]+1)->SetGrid(0,0);
	TString Title = Form("%1.1f < p_{T} < %1.1f",pt_low_phi[i_pt],pt_up_phi[i_pt]);
	h_flow_SM_total[i_order][i_cent][i_pt]->SetTitle(Title.Data());
	h_flow_SM_total[i_order][i_cent][i_pt]->SetStats(0);
	h_flow_SM_total[i_order][i_cent][i_pt]->GetXaxis()->SetNdivisions(505,'X');
	h_flow_SM_total[i_order][i_cent][i_pt]->GetYaxis()->SetNdivisions(505,'Y');
	h_flow_SM_total[i_order][i_cent][i_pt]->GetXaxis()->SetTitle("M(K^{+},K^{-}) (GeV/c^{2})");
	h_flow_SM_total[i_order][i_cent][i_pt]->GetYaxis()->SetTitle("Counts/Resolution");
	h_flow_SM_total[i_order][i_cent][i_pt]->GetXaxis()->SetTitleSize(0.04);
	h_flow_SM_total[i_order][i_cent][i_pt]->GetYaxis()->SetTitleSize(0.04);
	h_flow_SM_total[i_order][i_cent][i_pt]->GetXaxis()->CenterTitle();
	h_flow_SM_total[i_order][i_cent][i_pt]->GetYaxis()->CenterTitle();
	h_flow_SM_total[i_order][i_cent][i_pt]->SetMarkerStyle(24);
	h_flow_SM_total[i_order][i_cent][i_pt]->SetMarkerColor(1);
	h_flow_SM_total[i_order][i_cent][i_pt]->SetMarkerSize(0.8);
	h_flow_SM_total[i_order][i_cent][i_pt]->Draw("pE");
	f_PolyBW_total[i_order][i_cent][i_pt]->SetLineColor(2);
	f_PolyBW_total[i_order][i_cent][i_pt]->Draw("l same");
	TString FuncName = Form("f_Poly_Total_%s_%s_pt_%d",mOrder[i_order].Data(),mCentrality[i_cent].Data(),i_pt);
	f_Poly_total[i_order][i_cent][i_pt] = new TF1(FuncName.Data(),Poly,BW_Start,BW_Stop,2);
	f_Poly_total[i_order][i_cent][i_pt]->SetParameter(0,ParFit_total[i_order][i_cent][i_pt][3]);
	f_Poly_total[i_order][i_cent][i_pt]->SetParameter(1,ParFit_total[i_order][i_cent][i_pt][4]);
	f_Poly_total[i_order][i_cent][i_pt]->SetLineColor(4);
	f_Poly_total[i_order][i_cent][i_pt]->SetLineStyle(2);
	f_Poly_total[i_order][i_cent][i_pt]->SetLineWidth(2);
	f_Poly_total[i_order][i_cent][i_pt]->Draw("l same");
      }
      CanName = "./figures/" + CanName + ".eps";
      c_PolyBW_total[i_order][i_cent]->SaveAs(CanName.Data());
    }
  }
  */

  // Poly+BW Fit for phi-Psi bin and subtract linear background
  TF1 *f_PolyBW[2][Centrality_total][pt_total_New_phi][Phi_Psi_total];
  TF1 *f_Poly[2][Centrality_total][pt_total_New_phi][Phi_Psi_total];
  TH1F *h_flow_SM_Copy[2][Centrality_total][pt_total_New_phi][Phi_Psi_total];
  for(Int_t i_order = 0; i_order < 2; i_order++)
  {
    for(Int_t i_cent = Centrality_start; i_cent < Centrality_stop; i_cent++) // centrality bin
    {
      for(Int_t i_pt = 0; i_pt < pt_total_New_phi; i_pt++)
      {
	for(Int_t i_phi = 0; i_phi < Phi_Psi_total; i_phi++) // phi-psi bin
	{
	  TString HistName = Form("h_Copy_%s_%s_pt_%d_phi_Psi_%d",mOrder[i_order].Data(),mCentrality[i_cent].Data(),i_pt,i_phi);
	  h_flow_SM_Copy[i_order][i_cent][i_pt][i_phi] = (TH1F*)h_flow_SM_New[i_order][i_cent][i_pt][i_phi]->Clone(HistName.Data());

	  TString FuncName = Form("f_PolyBW_%s_%s_pt_%d_phi_Psi_%d",mOrder[i_order].Data(),mCentrality[i_cent].Data(),i_pt,i_phi);
	  f_PolyBW[i_order][i_cent][i_pt][i_phi] = new TF1(FuncName.Data(),PolyBreitWigner,BW_Start,BW_Stop,5);
	  for(Int_t n_par = 0; n_par < 5; n_par++)
	  {
	    f_PolyBW[i_order][i_cent][i_pt][i_phi]->ReleaseParameter(n_par);
	  }
	  f_PolyBW[i_order][i_cent][i_pt][i_phi]->FixParameter(0,ParFit_total[i_order][i_cent][i_pt][0]);
	  f_PolyBW[i_order][i_cent][i_pt][i_phi]->FixParameter(1,ParFit_total[i_order][i_cent][i_pt][1]);
	  f_PolyBW[i_order][i_cent][i_pt][i_phi]->SetParameter(2,ParFit_total[i_order][i_cent][i_pt][2]/7.0);
	  f_PolyBW[i_order][i_cent][i_pt][i_phi]->SetParameter(3,ParFit_total[i_order][i_cent][i_pt][3]);
	  f_PolyBW[i_order][i_cent][i_pt][i_phi]->SetParameter(4,ParFit_total[i_order][i_cent][i_pt][4]);
	  f_PolyBW[i_order][i_cent][i_pt][i_phi]->SetRange(BW_Start,BW_Stop);
	  
	  h_flow_SM_New[i_order][i_cent][i_pt][i_phi]->Fit(f_PolyBW[i_order][i_cent][i_pt][i_phi],"NQR");

	  FuncName = Form("f_Poly_%s_%s_phi_Psi_%d",mOrder[i_order].Data(),mCentrality[i_cent].Data(),i_phi);
	  f_Poly[i_order][i_cent][i_pt][i_phi] = new TF1(FuncName.Data(),Poly,BW_Start,BW_Stop,2);
	  f_Poly[i_order][i_cent][i_pt][i_phi]->FixParameter(0,ParFit_total[i_order][i_cent][i_pt][3]);
	  f_Poly[i_order][i_cent][i_pt][i_phi]->FixParameter(1,ParFit_total[i_order][i_cent][i_pt][4]);

	  h_flow_SM_New[i_order][i_cent][i_pt][i_phi]->Add(f_Poly[i_order][i_cent][i_pt][i_phi],-1.0);
	}
      }
    }
  }

  /*
  // QA plots: Poly+BW fits for phi-Psi bin
  TCanvas *c_PolyBW[2][Centrality_total][Phi_Psi_total];
  for(Int_t i_order = 0; i_order < 2; i_order++)
  {
    for(Int_t i_cent = Centrality_start; i_cent < Centrality_stop; i_cent++)
    {
      for(Int_t i_phi = 0; i_phi < Phi_Psi_total; i_phi++)
      {
	TString CanName = Form("c_PolyBW_%s_%s_phi_Psi_%d",mOrder[i_order].Data(),mCentrality[i_cent].Data(),i_phi);
	c_PolyBW[i_order][i_cent][i_phi] = new TCanvas(CanName.Data(),CanName.Data(),10,10,1000,1000);
	c_PolyBW[i_order][i_cent][i_phi]->Divide(5,5);
	for(Int_t i_pt = 0; i_pt < pt_total_New_phi; i_pt++)
	{
	  c_PolyBW[i_order][i_cent][i_phi]->cd(pt_new_bin_start[i_pt]+1);
	  c_PolyBW[i_order][i_cent][i_phi]->cd(pt_new_bin_start[i_pt]+1)->SetLeftMargin(0.15);
	  c_PolyBW[i_order][i_cent][i_phi]->cd(pt_new_bin_start[i_pt]+1)->SetBottomMargin(0.15);
	  c_PolyBW[i_order][i_cent][i_phi]->cd(pt_new_bin_start[i_pt]+1)->SetTicks(1,1);
	  c_PolyBW[i_order][i_cent][i_phi]->cd(pt_new_bin_start[i_pt]+1)->SetGrid(0,0);
	  TString Title = Form("%1.1f < p_{T} < %1.1f",pt_low_phi[i_pt],pt_up_phi[i_pt]);
	  h_flow_SM_Copy[i_order][i_cent][i_pt][i_phi]->SetTitle(Title.Data());
	  h_flow_SM_Copy[i_order][i_cent][i_pt][i_phi]->SetStats(0);
	  h_flow_SM_Copy[i_order][i_cent][i_pt][i_phi]->GetXaxis()->SetNdivisions(505,'X');
	  h_flow_SM_Copy[i_order][i_cent][i_pt][i_phi]->GetYaxis()->SetNdivisions(505,'Y');
	  h_flow_SM_Copy[i_order][i_cent][i_pt][i_phi]->GetXaxis()->SetTitle("M(K^{+},K^{-}) (GeV/c^{2})");
	  h_flow_SM_Copy[i_order][i_cent][i_pt][i_phi]->GetYaxis()->SetTitle("Counts/Resolution");
	  h_flow_SM_Copy[i_order][i_cent][i_pt][i_phi]->GetXaxis()->SetTitleSize(0.04);
	  h_flow_SM_Copy[i_order][i_cent][i_pt][i_phi]->GetYaxis()->SetTitleSize(0.04);
	  h_flow_SM_Copy[i_order][i_cent][i_pt][i_phi]->GetXaxis()->CenterTitle();
	  h_flow_SM_Copy[i_order][i_cent][i_pt][i_phi]->GetYaxis()->CenterTitle();
	  h_flow_SM_Copy[i_order][i_cent][i_pt][i_phi]->SetMarkerStyle(24);
	  h_flow_SM_Copy[i_order][i_cent][i_pt][i_phi]->SetMarkerColor(1);
	  h_flow_SM_Copy[i_order][i_cent][i_pt][i_phi]->SetMarkerSize(0.8);
	  h_flow_SM_Copy[i_order][i_cent][i_pt][i_phi]->Draw("pE");
	  f_PolyBW[i_order][i_cent][i_pt][i_phi]->SetLineColor(2);
	  f_PolyBW[i_order][i_cent][i_pt][i_phi]->Draw("l same");

	  f_Poly[i_order][i_cent][i_pt][i_phi]->SetLineColor(4);
	  f_Poly[i_order][i_cent][i_pt][i_phi]->SetLineStyle(2);
	  f_Poly[i_order][i_cent][i_pt][i_phi]->SetLineWidth(2);
	  f_Poly[i_order][i_cent][i_pt][i_phi]->Draw("l same");

//	  h_flow_SM_New[i_order][i_cent][i_pt][i_phi]->Draw("h same");
	}
	CanName = "./figures/" + CanName + ".eps";
	c_PolyBW[i_order][i_cent][i_phi]->SaveAs(CanName.Data());
      }
    }
  }
  */
#endif

  // Counting and BW fits to extract raw flow signal
  // merge phi-Psi distribution to extract fit range for signal
  TH1F *h_Sig_total[2][Centrality_total][pt_total_New_phi]; // merged distribution from phi-Psi distribution which subtracted linear background 
  TF1  *f_Sig_total[2][Centrality_total][pt_total_New_phi]; // Breit Wigner distribution
  Float_t ParFit_Sig_total[2][Centrality_total][pt_total_New_phi][3];
  Float_t Inte_start[2][Centrality_total][pt_total_New_phi];
  Float_t Inte_stop[2][Centrality_total][pt_total_New_phi];

  for(Int_t i_order = 0; i_order < 2; i_order++)
  {
    for(Int_t i_cent = Centrality_start; i_cent < Centrality_stop; i_cent++) // centrality bin
    {
      for(Int_t i_pt = 0; i_pt < pt_total_New_phi; i_pt++)
      {
	for(Int_t i_phi = 0; i_phi < Phi_Psi_total; i_phi++) // phi-psi bin => merge
	{
	  if(i_phi == 0)
	  {
	    TString HistName = Form("h_Sig_%s_%s_pt_%d",mOrder[i_order].Data(),mCentrality[i_cent].Data(),i_pt);
	    h_Sig_total[i_order][i_cent][i_pt] = (TH1F*)h_flow_SM_New[i_order][i_cent][i_pt][i_phi]->Clone(HistName.Data());
	    TString FuncName = Form("f_Sig_%s_%s_pt_%d",mOrder[i_order].Data(),mCentrality[i_cent].Data(),i_pt);
	    f_Sig_total[i_order][i_cent][i_pt] = new TF1(FuncName.Data(),BreitWigner,BW_Start,BW_Stop,3);
	    for(Int_t n_par = 0; n_par < 3; n_par++)
	    {
	      f_Sig_total[i_order][i_cent][i_pt]->ReleaseParameter(n_par);
	    }
//	    f_Sig_total[i_order][i_cent][i_pt]->SetParameter(0,ParFit_total[i_order][i_cent][i_pt][0]);
//	    f_Sig_total[i_order][i_cent][i_pt]->SetParameter(1,ParFit_total[i_order][i_cent][i_pt][1]);
//	    f_Sig_total[i_order][i_cent][i_pt]->SetParameter(2,ParFit_total[i_order][i_cent][i_pt][2]);
	    f_Sig_total[i_order][i_cent][i_pt]->SetParameter(0,1.019);
	    f_Sig_total[i_order][i_cent][i_pt]->SetParLimits(0,1.014,1.024);
	    f_Sig_total[i_order][i_cent][i_pt]->SetParameter(1,0.0055);
	    f_Sig_total[i_order][i_cent][i_pt]->SetParameter(2,10000);
	    f_Sig_total[i_order][i_cent][i_pt]->SetRange(BW_Start,BW_Stop);
	  }
	  else
	  {
	    h_Sig_total[i_order][i_cent][i_pt]->Add(h_flow_SM_New[i_order][i_cent][i_pt][i_phi],1.0);
	  }
	}
	h_Sig_total[i_order][i_cent][i_pt]->Fit(f_Sig_total[i_order][i_cent][i_pt],"NMQR");
	for(Int_t n_par = 0; n_par < 3; n_par++)
	{
	  ParFit_Sig_total[i_order][i_cent][i_pt][n_par] = f_Sig_total[i_order][i_cent][i_pt]->GetParameter(n_par);
	}
	// Integration range
	Inte_start[i_order][i_cent][i_pt] = ParFit_Sig_total[i_order][i_cent][i_pt][0] - nSigmaPhi*ParFit_Sig_total[i_order][i_cent][i_pt][1];
	Inte_stop[i_order][i_cent][i_pt]  = ParFit_Sig_total[i_order][i_cent][i_pt][0] + nSigmaPhi*ParFit_Sig_total[i_order][i_cent][i_pt][1];
      }
    }
  }

  // calculate total counts and errors for each phi-Psi bin by counting
  TH1F *h_Counts[2][Centrality_total][pt_total_New_phi]; // counts vs phi-Psi
  for(Int_t i_order = 0; i_order < 2; i_order++)
  {
    for(Int_t i_cent = Centrality_start; i_cent < Centrality_stop; i_cent++) // centrality bin
    {
      for(Int_t i_pt = 0; i_pt < pt_total_New_phi; i_pt++)
      {
	TString HistName = Form("h_Counts_%s_%s_pt_%d",mOrder[i_order].Data(),mCentrality[i_cent].Data(),i_pt);
	h_Counts[i_order][i_cent][i_pt] = new TH1F(HistName.Data(),HistName.Data(),7,0,PI_max[i_order]);
	for(Int_t i_phi = 0; i_phi < Phi_Psi_total; i_phi++) // phi-psi bin
	{
	  Float_t counts = 0.0;
	  Float_t errors = 0.0;
	  Float_t bin_center = PI_max[i_order]/14.0+i_phi*PI_max[i_order]/7.0;
	  Int_t bin_start = h_flow_SM_New[i_order][i_cent][i_pt][i_phi]->FindBin(Inte_start[i_order][i_cent][i_pt]);
	  Int_t bin_stop  = h_flow_SM_New[i_order][i_cent][i_pt][i_phi]->FindBin(Inte_stop[i_order][i_cent][i_pt]);
	  for(Int_t bin = bin_start; bin <= bin_stop; bin++)
	  {
	    counts += h_flow_SM_New[i_order][i_cent][i_pt][i_phi]->GetBinContent(bin);
	    errors += h_flow_SM_New[i_order][i_cent][i_pt][i_phi]->GetBinError(bin)*h_flow_SM_New[i_order][i_cent][i_pt][i_phi]->GetBinError(bin);
	  }
	  h_Counts[i_order][i_cent][i_pt]->SetBinContent(h_Counts[i_order][i_cent][i_pt]->FindBin(bin_center),counts);
	  h_Counts[i_order][i_cent][i_pt]->SetBinError(h_Counts[i_order][i_cent][i_pt]->FindBin(bin_center),TMath::Sqrt(errors));
	}
      }
    }
  }

  // QA plots: phi-Psi distribution for each pT bin
  TCanvas *c_phi_Psi[2][Centrality_total][pt_total_New_phi];
  for(Int_t i_order = 0; i_order < 2; i_order++)
  {
    for(Int_t i_cent = Centrality_start; i_cent < Centrality_stop; i_cent++)
    {
//      for(Int_t i_pt = 0; i_pt < pt_total_New_phi; i_pt++)
      for(Int_t i_pt = 1; i_pt < 6; i_pt++)
      {
	TString CanName = Form("c_phi_Psi_%s_%s_pt_%d",mOrder[i_order].Data(),mCentrality[i_cent].Data(),i_pt);
	c_phi_Psi[i_order][i_cent][i_pt] = new TCanvas(CanName.Data(),CanName.Data(),10,10,900,900);
	c_phi_Psi[i_order][i_cent][i_pt]->Divide(3,3);

	c_phi_Psi[i_order][i_cent][i_pt]->cd(8);
	c_phi_Psi[i_order][i_cent][i_pt]->cd(8)->SetLeftMargin(0.15);
	c_phi_Psi[i_order][i_cent][i_pt]->cd(8)->SetBottomMargin(0.15);
	c_phi_Psi[i_order][i_cent][i_pt]->cd(8)->SetTicks(1,1);
	c_phi_Psi[i_order][i_cent][i_pt]->cd(8)->SetGrid(0,0);
	h_Sig_total[i_order][i_cent][i_pt]->SetStats(0);
	h_Sig_total[i_order][i_cent][i_pt]->SetTitle("Merged Inv. Mass distribution");
	h_Sig_total[i_order][i_cent][i_pt]->GetXaxis()->SetNdivisions(505,'X');
	h_Sig_total[i_order][i_cent][i_pt]->GetYaxis()->SetNdivisions(505,'Y');
	h_Sig_total[i_order][i_cent][i_pt]->GetXaxis()->SetTitle("M(K^{+},K^{-}) (GeV/c^{2})");
	h_Sig_total[i_order][i_cent][i_pt]->GetYaxis()->SetTitle("Counts/Resolution");
	h_Sig_total[i_order][i_cent][i_pt]->GetXaxis()->SetTitleSize(0.04);
	h_Sig_total[i_order][i_cent][i_pt]->GetYaxis()->SetTitleSize(0.04);
	h_Sig_total[i_order][i_cent][i_pt]->GetXaxis()->CenterTitle();
	h_Sig_total[i_order][i_cent][i_pt]->GetYaxis()->CenterTitle();
//	h_Sig_total[i_order][i_cent][i_pt]->SetMarkerStyle(24);
//	h_Sig_total[i_order][i_cent][i_pt]->SetMarkerColor(1);
//	h_Sig_total[i_order][i_cent][i_pt]->SetMarkerSize(0.8);
	h_Sig_total[i_order][i_cent][i_pt]->Draw("pE");
	f_Sig_total[i_order][i_cent][i_pt]->SetLineColor(2);
	f_Sig_total[i_order][i_cent][i_pt]->Draw("l same");
	PlotLine(0.98,1.05,0.0,0.0,1,2,2);
	PlotLine(Inte_start[i_order][i_cent][i_pt],Inte_start[i_order][i_cent][i_pt],0.0,h_Sig_total[i_order][i_cent][i_pt]->GetMaximum(),4,2,2);
	PlotLine(Inte_stop[i_order][i_cent][i_pt],Inte_stop[i_order][i_cent][i_pt],0.0,h_Sig_total[i_order][i_cent][i_pt]->GetMaximum(),4,2,2);

	for(Int_t i_phi = 0; i_phi < Phi_Psi_total; i_phi++)
	{
	  c_phi_Psi[i_order][i_cent][i_pt]->cd(i_phi+1);
	  c_phi_Psi[i_order][i_cent][i_pt]->cd(i_phi+1)->SetLeftMargin(0.15);
	  c_phi_Psi[i_order][i_cent][i_pt]->cd(i_phi+1)->SetBottomMargin(0.15);
	  c_phi_Psi[i_order][i_cent][i_pt]->cd(i_phi+1)->SetTicks(1,1);
	  c_phi_Psi[i_order][i_cent][i_pt]->cd(i_phi+1)->SetGrid(0,0);
	  h_flow_SM_New[i_order][i_cent][i_pt][i_phi]->SetStats(0);
	  h_flow_SM_New[i_order][i_cent][i_pt][i_phi]->GetXaxis()->SetNdivisions(505,'X');
	  h_flow_SM_New[i_order][i_cent][i_pt][i_phi]->GetYaxis()->SetNdivisions(505,'Y');
	  h_flow_SM_New[i_order][i_cent][i_pt][i_phi]->GetXaxis()->SetTitle("M(K^{+},K^{-}) (GeV/c^{2})");
	  h_flow_SM_New[i_order][i_cent][i_pt][i_phi]->GetYaxis()->SetTitle("Counts/Resolution");
	  h_flow_SM_New[i_order][i_cent][i_pt][i_phi]->GetXaxis()->SetTitleSize(0.04);
	  h_flow_SM_New[i_order][i_cent][i_pt][i_phi]->GetYaxis()->SetTitleSize(0.04);
	  h_flow_SM_New[i_order][i_cent][i_pt][i_phi]->GetXaxis()->CenterTitle();
	  h_flow_SM_New[i_order][i_cent][i_pt][i_phi]->GetYaxis()->CenterTitle();
	  h_flow_SM_New[i_order][i_cent][i_pt][i_phi]->Draw("pE");
	  PlotLine(0.98,1.05,0.0,0.0,1,2,2);
	  PlotLine(Inte_start[i_order][i_cent][i_pt],Inte_start[i_order][i_cent][i_pt],0.0,h_flow_SM_New[i_order][i_cent][i_pt][i_phi]->GetMaximum(),4,2,2);
	  PlotLine(Inte_stop[i_order][i_cent][i_pt],Inte_stop[i_order][i_cent][i_pt],0.0,h_flow_SM_New[i_order][i_cent][i_pt][i_phi]->GetMaximum(),4,2,2);
	}
	CanName = "./figures/" + CanName + ".eps";
	c_phi_Psi[i_order][i_cent][i_pt]->SaveAs(CanName.Data());

	c_phi_Psi[i_order][i_cent][i_pt]->cd(9);
	c_phi_Psi[i_order][i_cent][i_pt]->cd(9)->SetLeftMargin(0.15);
	c_phi_Psi[i_order][i_cent][i_pt]->cd(9)->SetBottomMargin(0.15);
	c_phi_Psi[i_order][i_cent][i_pt]->cd(9)->SetTicks(1,1);
	c_phi_Psi[i_order][i_cent][i_pt]->cd(9)->SetGrid(0,0);
	h_Counts[i_order][i_cent][i_pt]->SetStats(0);
	h_Counts[i_order][i_cent][i_pt]->SetTitle("");
	h_Counts[i_order][i_cent][i_pt]->GetXaxis()->SetNdivisions(505,'X');
	h_Counts[i_order][i_cent][i_pt]->GetYaxis()->SetNdivisions(505,'Y');
	h_Counts[i_order][i_cent][i_pt]->GetXaxis()->SetTitle("#phi-#Psi");
	h_Counts[i_order][i_cent][i_pt]->GetYaxis()->SetTitle("Counts/Resolution");
	h_Counts[i_order][i_cent][i_pt]->GetXaxis()->SetTitleSize(0.04);
	h_Counts[i_order][i_cent][i_pt]->GetYaxis()->SetTitleSize(0.04);
	h_Counts[i_order][i_cent][i_pt]->GetXaxis()->CenterTitle();
	h_Counts[i_order][i_cent][i_pt]->GetYaxis()->CenterTitle();
	h_Counts[i_order][i_cent][i_pt]->SetMarkerStyle(24);
	h_Counts[i_order][i_cent][i_pt]->SetMarkerColor(1);
	h_Counts[i_order][i_cent][i_pt]->SetMarkerSize(0.8);
	h_Counts[i_order][i_cent][i_pt]->Draw("pE");
//	f_Sig_total[i_order][i_cent][i_pt]->SetLineColor(2);
//	f_Sig_total[i_order][i_cent][i_pt]->Draw("l same");
      }
    }
  }

  // read in TH1F for resolution correction
  TH1F *h_Centrality_SE[9];
  TH1F *h_Centrality_ME[9];
  TH1F *h_Centrality_SM[9];
  for(Int_t i_cent = 0; i_cent < 9; i_cent++) // harmonic loop
  {
    TString HistName = Form("h_mPhi_%d",i_cent);
    h_Centrality_SE[i_cent] = (TH1F*)File_input_SE->Get(HistName.Data());
    h_Centrality_ME[i_cent] = (TH1F*)File_input_ME->Get(HistName.Data());

    // initialize Histogram for SE_ME
    HistName = Form("h_mPhi_%d_SM",i_cent);
    h_Centrality_SM[i_cent] = (TH1F*)h_Centrality_SE[i_cent]->Clone(HistName.Data());

    // normalize the SE and ME
    Int_t Bin_SE_Start = h_Centrality_SE[i_cent]->FindBin(1.04);
    Int_t Bin_SE_Stop  = h_Centrality_SE[i_cent]->FindBin(1.05);
    Float_t Inte_SE = h_Centrality_SE[i_cent]->Integral(Bin_SE_Start,Bin_SE_Stop-1);

    Int_t Bin_ME_Start = h_Centrality_ME[i_cent]->FindBin(1.04);
    Int_t Bin_ME_Stop  = h_Centrality_ME[i_cent]->FindBin(1.05);
    Float_t Inte_ME = h_Centrality_ME[i_cent]->Integral(Bin_ME_Start,Bin_ME_Stop-1);
    if(Inte_ME > 0)
    {
      h_Centrality_ME[i_cent]->Scale(Inte_SE/Inte_ME);
      h_Centrality_SM[i_cent]->Add(h_Centrality_ME[i_cent],-1.0);
    }
  }
  /*
  h_Centrality_SE[7]->Draw("pE");
  h_Centrality_ME[7]->SetLineColor(2);
  h_Centrality_ME[7]->Draw("pE same");
  h_Centrality_SM[7]->SetLineColor(4);
  h_Centrality_SM[7]->Draw("PE same");
  */
}

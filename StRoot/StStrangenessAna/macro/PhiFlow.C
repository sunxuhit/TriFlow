#include "TFile.h"
#include "TH1F.h"
#include "TString.h"
#include "TMath.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TLatex.h"
#include "TLine.h"
#include "TPolyLine.h"

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

TLatex* plotTopLegend(char* label,Float_t x=-1,Float_t y=-1,Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1)
{
  // coordinates in NDC!
  // plots the string label in position x and y in NDC coordinates
  // size is the text size
  // color is the text color

  //    if(x<0||y<0)
  //    {   // defaults
  //      x=gPad->GetLeftMargin()*1.15;
  //      y=(1-gPad->GetTopMargin())*1.04;
  //    }
  TLatex* text=new TLatex(x,y,label);
  text->SetTextFont(font);
  text->SetTextSize(size);
  if(NDC == 1) text->SetNDC();
  text->SetTextColor(color);
  text->SetTextAngle(angle);
  text->Draw();
  return text;
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

static const TString Energy[2] = {"200GeV","39GeV"};
static const Int_t pt_total_phi = 16;
static const Int_t pt_total_New_phi = 8;
static const Int_t Centrality_total = 4;
static const Int_t Centrality_start = 0;
static const Int_t Centrality_stop = 1;
static const Int_t EtaGap_total = 4;
static const Int_t EtaGap_start = 0;
static const Int_t EtaGap_stop = 1;
static const Int_t Phi_Psi_total = 7;
static const Float_t nSigmaPhi = 2.5;
static const Double_t PI_max[2] = {TMath::Pi()/2.0,TMath::Pi()/3.0};

// pt bin
//                                 0 , 1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 ,10 ,11 ,12 ,13 ,14 ,15
/*
static Float_t pt_low_phi[16] = {0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4};
static Float_t pt_up_phi[16]  = {0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,4.2};
*/
static Float_t pt_low_phi[8] = {0.4,0.8,1.0,1.2,1.4,1.8,2.2,2.6};
static Float_t pt_up_phi[8]  = {0.8,1.0,1.2,1.4,1.8,2.2,2.6,3.2};

void PhiFlow(const Int_t energy = 1)
{
  // open input file for same event and mixed event
  TString inputfile_SE = Form("/project/projectdirs/star/xusun/OutPut/AuAu39GeV/Phi/flow_phi/merged_file/Yields_SE_%s.root",Energy[energy].Data());
  TFile *File_SE = TFile::Open(inputfile_SE.Data());
  TString inputfile_ME = Form("/project/projectdirs/star/xusun/OutPut/AuAu39GeV/Phi/flow_phi/merged_file/Yields_ME_%s.root",Energy[energy].Data());
  TFile *File_ME = TFile::Open(inputfile_ME.Data());

  // read histogram for same event and mixed event
  TH1F *h_mMass_Phi2_SE[pt_total_phi][Centrality_total][EtaGap_total][Phi_Psi_total];
  TH1F *h_mMass_Phi3_SE[pt_total_phi][Centrality_total][EtaGap_total][Phi_Psi_total];

  TH1F *h_mMass_Phi2_ME[pt_total_phi][Centrality_total][EtaGap_total][Phi_Psi_total];
  TH1F *h_mMass_Phi3_ME[pt_total_phi][Centrality_total][EtaGap_total][Phi_Psi_total];

  for(Int_t i = 0; i < pt_total_phi; i++) // pt bin
  {
    for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
    {
      for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
      {
	for(Int_t m = 0; m < Phi_Psi_total; m++) // phi-psi bin
	{
	  TString HistName;

	  HistName = Form("pt_%d_Centrality_%d_EtaGap_%d_phi_Psi_%d_2nd_Phi_SE",i,j,l,m);
	  h_mMass_Phi2_SE[i][j][l][m] = (TH1F*)File_SE->Get(HistName.Data())->Clone();
	  HistName = Form("pt_%d_Centrality_%d_EtaGap_%d_phi_Psi_%d_3rd_Phi_SE",i,j,l,m);
	  h_mMass_Phi3_SE[i][j][l][m] = (TH1F*)File_SE->Get(HistName.Data())->Clone();

	  HistName = Form("pt_%d_Centrality_%d_EtaGap_%d_phi_Psi_%d_2nd_Phi_ME",i,j,l,m);
	  h_mMass_Phi2_ME[i][j][l][m] = (TH1F*)File_ME->Get(HistName.Data())->Clone();
	  HistName = Form("pt_%d_Centrality_%d_EtaGap_%d_phi_Psi_%d_3rd_Phi_ME",i,j,l,m);
	  h_mMass_Phi3_ME[i][j][l][m] = (TH1F*)File_ME->Get(HistName.Data())->Clone();
	}
      }
    }
  }
  /*
  File_SE->Close();
  File_ME->Close();
  */

  // calculate scaling factor and subtract mixed event from same event
  TH1F *h_mMass_Phi2_SM[pt_total_phi][Centrality_total][EtaGap_total][Phi_Psi_total];
  TH1F *h_mMass_Phi3_SM[pt_total_phi][Centrality_total][EtaGap_total][Phi_Psi_total];
  for(Int_t i = 0; i < pt_total_phi; i++) // pt bin
  {
    for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
    {
      for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
      {
	for(Int_t m = 0; m < Phi_Psi_total; m++) // phi-psi bin
	{
	  Int_t bin2_SE_start = h_mMass_Phi2_SE[i][j][l][m]->FindBin(1.04);
	  Int_t bin2_SE_stop  = h_mMass_Phi2_SE[i][j][l][m]->FindBin(1.05);
	  Float_t Inte2_SE = h_mMass_Phi2_SE[i][j][l][m]->Integral(bin2_SE_start,bin2_SE_stop);
	  h_mMass_Phi2_SM[i][j][l][m] = (TH1F*)h_mMass_Phi2_SE[i][j][l][m]->Clone();

	  Int_t bin3_SE_start = h_mMass_Phi3_SE[i][j][l][m]->FindBin(1.04);
	  Int_t bin3_SE_stop  = h_mMass_Phi3_SE[i][j][l][m]->FindBin(1.05);
	  Float_t Inte3_SE = h_mMass_Phi3_SE[i][j][l][m]->Integral(bin3_SE_start,bin3_SE_stop);
	  h_mMass_Phi3_SM[i][j][l][m] = (TH1F*)h_mMass_Phi3_SE[i][j][l][m]->Clone();

	  Int_t bin2_ME_start = h_mMass_Phi2_ME[i][j][l][m]->FindBin(1.04);
	  Int_t bin2_ME_stop  = h_mMass_Phi2_ME[i][j][l][m]->FindBin(1.05);
	  Float_t Inte2_ME = h_mMass_Phi2_ME[i][j][l][m]->Integral(bin2_ME_start,bin2_ME_stop);
	  h_mMass_Phi2_ME[i][j][l][m]->Scale(Inte2_SE/Inte2_ME);
	  h_mMass_Phi2_SM[i][j][l][m]->Add(h_mMass_Phi2_ME[i][j][l][m],-1.0);

	  Int_t bin3_ME_start = h_mMass_Phi3_ME[i][j][l][m]->FindBin(1.04);
	  Int_t bin3_ME_stop  = h_mMass_Phi3_ME[i][j][l][m]->FindBin(1.05);
	  Float_t Inte3_ME = h_mMass_Phi3_ME[i][j][l][m]->Integral(bin3_ME_start,bin3_ME_stop);
	  h_mMass_Phi3_ME[i][j][l][m]->Scale(Inte3_SE/Inte3_ME);
	  h_mMass_Phi3_SM[i][j][l][m]->Add(h_mMass_Phi3_ME[i][j][l][m],-1.0);
	}
      }
    }
  }

  // merged original pT bins to get new pT bins
  TH1F *h_mMass_Phi2_SM_New_Bin[pt_total_New_phi][Centrality_total][EtaGap_total][Phi_Psi_total];
  TH1F *h_mMass_Phi3_SM_New_Bin[pt_total_New_phi][Centrality_total][EtaGap_total][Phi_Psi_total];

  Int_t pt_new_bin_start[pt_total_New_phi] = {1,3,4,5,6,8,10,12};
  Int_t pt_new_bin_stop[pt_total_New_phi]  = {2,3,4,5,7,9,11,14};
  for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
  {
    for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
    {
      for(Int_t m = 0; m < Phi_Psi_total; m++) // phi-psi bin
      {
	for(Int_t i = 0; i < pt_total_New_phi; i++)
	{
	  for(Int_t pt_bin = pt_new_bin_start[i]; pt_bin <= pt_new_bin_stop[i]; pt_bin++)
	  {
	    if(pt_bin == pt_new_bin_start[i])
	    {
	      h_mMass_Phi2_SM_New_Bin[i][j][l][m] = (TH1F*)h_mMass_Phi2_SM[pt_bin][j][l][m]->Clone();
	      h_mMass_Phi3_SM_New_Bin[i][j][l][m] = (TH1F*)h_mMass_Phi3_SM[pt_bin][j][l][m]->Clone();
	    }
	    else
	    {
	      h_mMass_Phi2_SM_New_Bin[i][j][l][m]->Add(h_mMass_Phi2_SM[pt_bin][j][l][m],1.0);
	      h_mMass_Phi3_SM_New_Bin[i][j][l][m]->Add(h_mMass_Phi3_SM[pt_bin][j][l][m],1.0);
	    }
	  }
	  if(i >= 6)
	  {
	    h_mMass_Phi2_SM_New_Bin[i][j][l][m]->Rebin(2);
	    h_mMass_Phi3_SM_New_Bin[i][j][l][m]->Rebin(2);
	  }
	}
      }
    }
  }


  /*
  // Draw new bins
  TCanvas *c_Phi2_New_Bin[Centrality_total][EtaGap_total][Phi_Psi_total];
  TCanvas *c_Phi3_New_Bin[Centrality_total][EtaGap_total][Phi_Psi_total];

  for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
  {
    for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
    {
      for(Int_t m = 0; m < Phi_Psi_total; m++) // phi-psi bin
      {
	TString CanName;

	CanName = Form("c_Phi2_SM_New_Bin_Centrality_%d_EtaGap_%d_phi_Psi_%d",j,l,m);
	c_Phi2_New_Bin[j][l][m] = new TCanvas(CanName.Data(),CanName.Data(),1400,10,900,900);
	c_Phi2_New_Bin[j][l][m]->Divide(3,3);
	for(Int_t i = 0; i < 8; i++)
	{
	  c_Phi2_New_Bin[j][l][m]->cd(i+1);
	  h_mMass_Phi2_SM_New_Bin[i][j][l][m]->SetTitle("");
	  h_mMass_Phi2_SM_New_Bin[i][j][l][m]->SetStats(0);
	  h_mMass_Phi2_SM_New_Bin[i][j][l][m]->DrawCopy("PE");
	}

	CanName = Form("c_Phi3_SM_New_Bin_Centrality_%d_EtaGap_%d_phi_Psi_%d",j,l,m);
	c_Phi3_New_Bin[j][l][m] = new TCanvas(CanName.Data(),CanName.Data(),1400,10,900,900);
	c_Phi3_New_Bin[j][l][m]->Divide(3,3);
	for(Int_t i = 0; i < 8; i++)
	{
	  c_Phi3_New_Bin[j][l][m]->cd(i+1);
	  h_mMass_Phi3_SM_New_Bin[i][j][l][m]->SetTitle("");
	  h_mMass_Phi3_SM_New_Bin[i][j][l][m]->SetStats(0);
	  h_mMass_Phi3_SM_New_Bin[i][j][l][m]->DrawCopy("PE");
	}
      }
    }
  }
  */

  // Poly+BW Fit for merged phi-psi bin
  TH1F *h_mMass_Phi2_SM_total[pt_total_New_phi][Centrality_total][EtaGap_total];
  TH1F *h_mMass_Phi3_SM_total[pt_total_New_phi][Centrality_total][EtaGap_total];

  TF1  *f_PolyBW_Phi2_total[pt_total_New_phi][Centrality_total][EtaGap_total];
  TF1  *f_PolyBW_Phi3_total[pt_total_New_phi][Centrality_total][EtaGap_total];

  Float_t ParFit2_total[pt_total_New_phi][Centrality_total][EtaGap_total][5];
  Float_t ParFit3_total[pt_total_New_phi][Centrality_total][EtaGap_total][5];

  for(Int_t i = 0; i < pt_total_New_phi; i++) // pt bin
  {
    for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
    {
      for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
      {
	for(Int_t m = 0; m < Phi_Psi_total; m++) // phi-psi bin
	{
	  TString HistName;
	  if(m == 0)
	  {
	    HistName = Form("pt_%d_Centrality_%d_EtaGap_%d_Total_2nd_Phi_SE",i,j,l);
	    h_mMass_Phi2_SM_total[i][j][l] = (TH1F*)h_mMass_Phi2_SM_New_Bin[i][j][l][m]->Clone(HistName.Data());
	    HistName = "f_"+HistName;
	    f_PolyBW_Phi2_total[i][j][l] = new TF1(HistName.Data(),PolyBreitWigner,0.994,1.05,5);
	    for(Int_t n_par = 0; n_par < 5; n_par++)
	    {
	      f_PolyBW_Phi2_total[i][j][l]->ReleaseParameter(n_par);
	    }
	    f_PolyBW_Phi2_total[i][j][l]->SetParameter(0,1.019);
//	    f_PolyBW_Phi2_total[i][j][l]->SetParLimits(0,1.014,1.024);
	    f_PolyBW_Phi2_total[i][j][l]->SetParameter(1,0.0055);
	    f_PolyBW_Phi2_total[i][j][l]->SetParameter(2,10000);
	    f_PolyBW_Phi2_total[i][j][l]->SetParameter(3,-6000);
	    f_PolyBW_Phi2_total[i][j][l]->SetParameter(4,0.5);
	    f_PolyBW_Phi2_total[i][j][l]->SetRange(0.994,1.05);

	    HistName = Form("pt_%d_Centrality_%d_EtaGap_%d_Total_3rd_Phi_SE",i,j,l);
	    h_mMass_Phi3_SM_total[i][j][l] = (TH1F*)h_mMass_Phi3_SM_New_Bin[i][j][l][m]->Clone(HistName.Data());
	    HistName = "f_"+HistName;
	    f_PolyBW_Phi3_total[i][j][l] = new TF1(HistName.Data(),PolyBreitWigner,0.994,1.05,5);
	    for(Int_t n_par = 0; n_par < 4; n_par++)
	    {
	      f_PolyBW_Phi3_total[i][j][l]->ReleaseParameter(n_par);
	    }
	    f_PolyBW_Phi3_total[i][j][l]->SetParameter(0,1.019);
//	    f_PolyBW_Phi3_total[i][j][l]->SetParLimits(0,1.014,1.024);
	    f_PolyBW_Phi3_total[i][j][l]->SetParameter(1,0.0055);
	    f_PolyBW_Phi3_total[i][j][l]->SetParameter(2,10000);
	    f_PolyBW_Phi3_total[i][j][l]->SetParameter(3,-6000);
	    f_PolyBW_Phi3_total[i][j][l]->SetParameter(4,0.5);
	    f_PolyBW_Phi3_total[i][j][l]->SetRange(0.994,1.05);
	  }
	  else
	  {
	    h_mMass_Phi2_SM_total[i][j][l]->Add(h_mMass_Phi2_SM_New_Bin[i][j][l][m],1.0);
	    h_mMass_Phi3_SM_total[i][j][l]->Add(h_mMass_Phi3_SM_New_Bin[i][j][l][m],1.0);
	  }
	}
	h_mMass_Phi2_SM_total[i][j][l]->Fit(f_PolyBW_Phi2_total[i][j][l],"NQR");
	h_mMass_Phi3_SM_total[i][j][l]->Fit(f_PolyBW_Phi3_total[i][j][l],"NQR");
	for(Int_t n_par = 0; n_par < 5; n_par++)
	{
	  ParFit2_total[i][j][l][n_par] = f_PolyBW_Phi2_total[i][j][l]->GetParameter(n_par);
	  ParFit3_total[i][j][l][n_par] = f_PolyBW_Phi3_total[i][j][l]->GetParameter(n_par);
	}
      }
    }
  }

  /*
  // Draw the merged histogram and fit line
  TCanvas *c_mMass_Phi2_SM_total[Centrality_total][EtaGap_total];
  TCanvas *c_mMass_Phi3_SM_total[Centrality_total][EtaGap_total];

  for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
  {
    for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
    {
      TString CanName;

      CanName = Form("c_Phi2_SM_Total_Centrality_%d_EtaGap_%d",j,l);
      c_mMass_Phi2_SM_total[j][l] = new TCanvas(CanName.Data(),CanName.Data(),1400,10,900,900);
      c_mMass_Phi2_SM_total[j][l]->Divide(3,3);
      for(Int_t i = 0; i < pt_total_New_phi; i++)
      {
	c_mMass_Phi2_SM_total[j][l]->cd(i+1);
	h_mMass_Phi2_SM_total[i][j][l]->Draw("PE");
	f_PolyBW_Phi2_total[i][j][l]->SetLineColor(2);
	f_PolyBW_Phi2_total[i][j][l]->Draw("l same");
      }

      CanName = Form("c_Phi3_SM_Total_Centrality_%d_EtaGap_%d",j,l);
      c_mMass_Phi3_SM_total[j][l] = new TCanvas(CanName.Data(),CanName.Data(),1400,10,900,900);
      c_mMass_Phi3_SM_total[j][l]->Divide(3,3);
      for(Int_t i = 0; i < pt_total_New_phi; i++)
      {
	c_mMass_Phi3_SM_total[j][l]->cd(i+1);
	h_mMass_Phi3_SM_total[i][j][l]->Draw("PE");
	f_PolyBW_Phi3_total[i][j][l]->SetLineColor(2);
	f_PolyBW_Phi3_total[i][j][l]->Draw("l same");
      }
    }
  }
  */

  // Poly+BW Fit for phi-psi bin
  TF1 *f_PolyBW_Phi2[pt_total_New_phi][Centrality_total][EtaGap_total][Phi_Psi_total];
  TF1 *f_PolyBW_Phi3[pt_total_New_phi][Centrality_total][EtaGap_total][Phi_Psi_total];
  for(Int_t i = 0; i < pt_total_New_phi; i++) // pt bin
  {
    for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
    {
      for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
      {
	for(Int_t m = 0; m < Phi_Psi_total; m++) // phi-psi bin
	{
	  TString HistName;

	  HistName = Form("f_pt_%d_Centrality_%d_EtaGap_%d_phi_Psi_%d_2nd_Phi_SE",i,j,l,m);
	  f_PolyBW_Phi2[i][j][l][m] = new TF1(HistName.Data(),PolyBreitWigner,0.994,1.05,5);
	  for(Int_t n_par = 0; n_par < 5; n_par++)
	  {
	    f_PolyBW_Phi2[i][j][l][m]->ReleaseParameter(n_par);
	  }
	  f_PolyBW_Phi2[i][j][l][m]->FixParameter(0,ParFit2_total[i][j][l][0]);
	  f_PolyBW_Phi2[i][j][l][m]->FixParameter(1,ParFit2_total[i][j][l][1]);
	  f_PolyBW_Phi2[i][j][l][m]->SetParameter(2,ParFit2_total[i][j][l][2]/7.0);
	  f_PolyBW_Phi2[i][j][l][m]->SetParameter(3,ParFit2_total[i][j][l][3]);
	  f_PolyBW_Phi2[i][j][l][m]->SetParameter(4,ParFit2_total[i][j][l][4]);
	  f_PolyBW_Phi2[i][j][l][m]->SetRange(0.994,1.05);

	  HistName = Form("f_pt_%d_Centrality_%d_EtaGap_%d_phi_Psi_%d_3rd_Phi_SE",i,j,l,m);
	  f_PolyBW_Phi3[i][j][l][m] = new TF1(HistName.Data(),PolyBreitWigner,0.994,1.05,5);
	  for(Int_t n_par = 0; n_par < 5; n_par++)
	  {
	    f_PolyBW_Phi3[i][j][l][m]->ReleaseParameter(n_par);
	  }
	  f_PolyBW_Phi3[i][j][l][m]->FixParameter(0,ParFit3_total[i][j][l][0]);
	  f_PolyBW_Phi3[i][j][l][m]->FixParameter(1,ParFit3_total[i][j][l][1]);
	  f_PolyBW_Phi3[i][j][l][m]->SetParameter(2,ParFit3_total[i][j][l][2]/7.0);
	  f_PolyBW_Phi3[i][j][l][m]->SetParameter(3,ParFit3_total[i][j][l][3]);
	  f_PolyBW_Phi3[i][j][l][m]->SetParameter(4,ParFit3_total[i][j][l][4]);
	  f_PolyBW_Phi3[i][j][l][m]->SetRange(0.994,1.05);

	  h_mMass_Phi2_SM_New_Bin[i][j][l][m]->Fit(f_PolyBW_Phi2[i][j][l][m],"NQR");
	  h_mMass_Phi3_SM_New_Bin[i][j][l][m]->Fit(f_PolyBW_Phi3[i][j][l][m],"NQR");
	}
      }
    }
  }

  /*
  TCanvas *c_mMass_Phi2_SM[Centrality_total][EtaGap_total][Phi_Psi_total];
  TCanvas *c_mMass_Phi3_SM[Centrality_total][EtaGap_total][Phi_Psi_total];

  for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
  {
    for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
    {
      for(Int_t m = 0; m < Phi_Psi_total; m++) // phi-psi bin
      {
	TString CanName;

	CanName = Form("c_Phi2_SM_Centrality_%d_EtaGap_%d_phi_Psi_%d",j,l,m);
	c_mMass_Phi2_SM[j][l][m] = new TCanvas(CanName.Data(),CanName.Data(),10,1400,900,900);
	c_mMass_Phi2_SM[j][l][m]->Divide(3,3);
	for(Int_t i = 0; i < pt_total_New_phi; i++)
	{
	  c_mMass_Phi2_SM[j][l][m]->cd(i+1)->SetLeftMargin(0.25);
	  c_mMass_Phi2_SM[j][l][m]->cd(i+1)->SetTicks(1,1);
	  c_mMass_Phi2_SM[j][l][m]->cd(i+1)->SetGrid(0,0);
	  c_mMass_Phi2_SM[j][l][m]->cd(i+1);
	  h_mMass_Phi2_SM_New_Bin[i][j][l][m]->SetTitle("");
	  h_mMass_Phi2_SM_New_Bin[i][j][l][m]->SetStats(0);
//	  h_mMass_Phi2_SM_New_Bin[i][j][l][m]->DrawCopy("PE");
	  f_PolyBW_Phi2[i][j][l][m]->SetLineColor(4);
	  f_PolyBW_Phi2[i][j][l][m]->Draw("l same");
	}

	CanName = Form("c_Phi3_SM_Centrality_%d_EtaGap_%d_phi_Psi_%d",j,l,m);
	c_mMass_Phi3_SM[j][l][m] = new TCanvas(CanName.Data(),CanName.Data(),10,1400,900,900);
	c_mMass_Phi3_SM[j][l][m]->Divide(3,3);
	for(Int_t i = 0; i < pt_total_New_phi; i++)
	{
	  c_mMass_Phi3_SM[j][l][m]->cd(i+1);
	  h_mMass_Phi3_SM_New_Bin[i][j][l][m]->SetTitle("");
	  h_mMass_Phi3_SM_New_Bin[i][j][l][m]->SetStats(0);
//	  h_mMass_Phi3_SM_New_Bin[i][j][l][m]->DrawCopy("PE");
	  f_PolyBW_Phi3[i][j][l][m]->SetLineColor(4);
	  f_PolyBW_Phi3[i][j][l][m]->Draw("l same");
	}
      }
    }
  }
  */

  TF1 *f_Poly_Phi2[pt_total_New_phi][Centrality_total][EtaGap_total][Phi_Psi_total];
  TF1 *f_Poly_Phi3[pt_total_New_phi][Centrality_total][EtaGap_total][Phi_Psi_total];

  for(Int_t i = 0; i < pt_total_New_phi; i++)
  {
    for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
    {
      for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
      {
	for(Int_t m = 0; m < Phi_Psi_total; m++) // phi-psi bin
	{
	  TString FuncName;

	  FuncName = Form("f_poly2_pt_%d_Centrality_%d_EtaGap_%d_phi_Psi_%d",i,j,l,m);
	  f_Poly_Phi2[i][j][l][m] = new TF1(FuncName.Data(),Poly,0.994,1.05,2);
	  f_Poly_Phi2[i][j][l][m]->FixParameter(0,f_PolyBW_Phi2[i][j][l][m]->GetParameter(3));
	  f_Poly_Phi2[i][j][l][m]->FixParameter(1,f_PolyBW_Phi2[i][j][l][m]->GetParameter(4));
	  h_mMass_Phi2_SM_New_Bin[i][j][l][m]->Add(f_Poly_Phi2[i][j][l][m],-1.0);

	  FuncName = Form("f_poly3_pt_%d_Centrality_%d_EtaGap_%d_phi_Psi_%d",i,j,l,m);
	  f_Poly_Phi3[i][j][l][m] = new TF1(FuncName.Data(),Poly,0.994,1.05,2);
	  f_Poly_Phi3[i][j][l][m]->FixParameter(0,f_PolyBW_Phi3[i][j][l][m]->GetParameter(3));
	  f_Poly_Phi3[i][j][l][m]->FixParameter(1,f_PolyBW_Phi3[i][j][l][m]->GetParameter(4));
	  h_mMass_Phi3_SM_New_Bin[i][j][l][m]->Add(f_Poly_Phi3[i][j][l][m],-1.0);
	}
      }
    }
  }

  /*
  // Draw the histogram after subtract poly background for phi-Psi bin
  for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
  {
    for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
    {
      for(Int_t m = 0; m < Phi_Psi_total; m++) // phi-psi bin
      {
	for(Int_t i = 0; i < pt_total_New_phi; i++)
	{
	  TString leg_pt = Form("p_{T} = %1.1f GeV/c - %1.1f GeV/c",pt_low_phi[i],pt_up_phi[i]);

	  c_mMass_Phi2_SM[j][l][m]->cd(i+1);
	  h_mMass_Phi2_SM_New_Bin[i][j][l][m]->SetLineColor(4);
	  h_mMass_Phi2_SM_New_Bin[i][j][l][m]->DrawCopy("PE");
	  // f_Poly_Phi2[i][j][l][m]->SetLineColor(2);
	  // f_Poly_Phi2[i][j][l][m]->SetLineStyle(2);
	  // f_Poly_Phi2[i][j][l][m]->Draw("l same");
	  plotTopLegend((char*)leg_pt.Data(),0.2,0.94,0.06,1,0.0,42,1);
	  Float_t x1 = ParFit2_total[i][j][l][0]-nSigmaPhi*ParFit2_total[i][j][l][1];
	  Float_t x2 = ParFit2_total[i][j][l][0]+nSigmaPhi*ParFit2_total[i][j][l][1];
	  Float_t y  = h_mMass_Phi2_SM_New_Bin[i][j][l][m]->GetBinContent(h_mMass_Phi2_SM_New_Bin[i][j][l][m]->FindBin(ParFit2_total[i][j][l][0]));
	  PlotLine(x1,x1,0,y,1,2,2);
	  PlotLine(x2,x2,0,y,1,2,2);

	  c_mMass_Phi3_SM[j][l][m]->cd(i+1);
	  h_mMass_Phi3_SM_New_Bin[i][j][l][m]->SetLineColor(4);
	  h_mMass_Phi3_SM_New_Bin[i][j][l][m]->DrawCopy("PE");
	  // f_Poly_Phi3[i][j][l][m]->SetLineColor(2);
	  // f_Poly_Phi3[i][j][l][m]->SetLineStyle(2);
	  // f_Poly_Phi3[i][j][l][m]->Draw("l same");
	  plotTopLegend((char*)leg_pt.Data(),0.2,0.94,0.06,1,0.0,42,1);
	  x1 = ParFit3_total[i][j][l][0]-nSigmaPhi*ParFit3_total[i][j][l][1];
	  x2 = ParFit3_total[i][j][l][0]+nSigmaPhi*ParFit3_total[i][j][l][1];
	  y  = h_mMass_Phi3_SM_New_Bin[i][j][l][m]->GetBinContent(h_mMass_Phi3_SM_New_Bin[i][j][l][m]->FindBin(ParFit3_total[i][j][l][0]));
	  PlotLine(x1,x1,0,y,1,2,2);
	  PlotLine(x2,x2,0,y,1,2,2);
	}
      }
    }
  }
  */

  // calculate total counts and errors for each phi-Psi bin
  TH1F *h_Counts2[pt_total_New_phi][Centrality_total][EtaGap_total];
  TH1F *h_Counts3[pt_total_New_phi][Centrality_total][EtaGap_total];
  for(Int_t i = 0; i < pt_total_New_phi; i++)
  {
    for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
    {
      for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
      {
	TString HistName;

	HistName = Form("Counts2_pt_%d_Centrality_%d_EtaGap_%d",i,j,l);
	h_Counts2[i][j][l] = new TH1F(HistName.Data(),HistName.Data(),100,0,PI_max[0]);
	for(Int_t m = 0; m < Phi_Psi_total; m++) // phi-psi bin
	{
	  Float_t counts = 0.0;
	  Float_t errors = 0.0;
	  Float_t bin_center = PI_max[0]/14.0+m*PI_max[0]/7.0;
	  Int_t bin_start = h_mMass_Phi2_SM_New_Bin[i][j][l][m]->FindBin(ParFit2_total[i][j][l][0]-nSigmaPhi*ParFit2_total[i][j][l][1]);
	  Int_t bin_stop  = h_mMass_Phi2_SM_New_Bin[i][j][l][m]->FindBin(ParFit2_total[i][j][l][0]+nSigmaPhi*ParFit2_total[i][j][l][1]);
	  for(Int_t bin = bin_start; bin <= bin_stop; bin++)
	  {
	    counts += h_mMass_Phi2_SM_New_Bin[i][j][l][m]->GetBinContent(bin);
	    errors += h_mMass_Phi2_SM_New_Bin[i][j][l][m]->GetBinError(bin)*h_mMass_Phi2_SM_New_Bin[i][j][l][m]->GetBinError(bin);
	  }
	  h_Counts2[i][j][l]->SetBinContent(h_Counts2[i][j][l]->FindBin(bin_center),counts);
	  h_Counts2[i][j][l]->SetBinError(h_Counts2[i][j][l]->FindBin(bin_center),TMath::Sqrt(errors));
	}

	HistName = Form("Counts3_pt_%d_Centrality_%d_EtaGap_%d",i,j,l);
	h_Counts3[i][j][l] = new TH1F(HistName.Data(),HistName.Data(),100,0,PI_max[1]);
	for(Int_t m = 0; m < Phi_Psi_total; m++) // phi-psi bin
	{
	  Float_t counts = 0.0;
	  Float_t errors = 0.0;
	  Float_t bin_center = PI_max[1]/14.0+m*PI_max[1]/7.0;
	  Int_t bin_start = h_mMass_Phi3_SM_New_Bin[i][j][l][m]->FindBin(ParFit3_total[i][j][l][0]-nSigmaPhi*ParFit3_total[i][j][l][1]);
	  Int_t bin_stop  = h_mMass_Phi3_SM_New_Bin[i][j][l][m]->FindBin(ParFit3_total[i][j][l][0]+nSigmaPhi*ParFit3_total[i][j][l][1]);
	  for(Int_t bin = bin_start; bin <= bin_stop; bin++)
	  {
	    counts += h_mMass_Phi3_SM_New_Bin[i][j][l][m]->GetBinContent(bin);
	    errors += h_mMass_Phi3_SM_New_Bin[i][j][l][m]->GetBinError(bin)*h_mMass_Phi3_SM_New_Bin[i][j][l][m]->GetBinError(bin);
	  }
	  h_Counts3[i][j][l]->SetBinContent(h_Counts3[i][j][l]->FindBin(bin_center),counts);
	  h_Counts3[i][j][l]->SetBinError(h_Counts3[i][j][l]->FindBin(bin_center),TMath::Sqrt(errors));
	}
      }
    }
  }

  // do cos fit to extract raw v2 and v3
  TF1 *f_phi2[pt_total_New_phi][Centrality_total][EtaGap_total];
  TF1 *f_phi3[pt_total_New_phi][Centrality_total][EtaGap_total];
  for(Int_t i = 0; i < pt_total_New_phi; i++)
  {
    for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
    {
      for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
      {
	TString Flow_phi;

	Flow_phi = Form("flow_pt_%d_Centrality_%d_etagap_%d_2nd_phi",i,j,l);
	f_phi2[i][j][l] = new TF1(Flow_phi.Data(),flow_2,0.0,PI_max[0],2);
	f_phi2[i][j][l]->SetParameter(0,2.0);
	f_phi2[i][j][l]->SetParameter(1,1.0);
	h_Counts2[i][j][l]->Fit(f_phi2[i][j][l],"NQI");

	Flow_phi = Form("flow_pt_%d_Centrality_%d_etagap_%d_3rd_phi",i,j,l);
	f_phi3[i][j][l] = new TF1(Flow_phi.Data(),flow_3,0.0,PI_max[1],2);
	f_phi3[i][j][l]->SetParameter(0,2.0);
	f_phi3[i][j][l]->SetParameter(1,1.0);
	h_Counts3[i][j][l]->Fit(f_phi3[i][j][l],"NQI");
      }
    }
  }

  /*
  // Draw phi-Psi distribution for all pt bin
  TCanvas *c2_phi_Psi[pt_total_New_phi][Centrality_total][EtaGap_total];
  TCanvas *c3_phi_Psi[pt_total_New_phi][Centrality_total][EtaGap_total];
  for(Int_t i = 0; i < pt_total_New_phi; i++)
  {
    for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
    {
      for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
      {
	TString CanName;

	CanName = Form("c2_phi_Psi_pt_%d_Centrality_%d_EtaGap_%d",i,j,l);
	c2_phi_Psi[i][j][l] = new TCanvas(CanName.Data(),CanName.Data(),1400,10,1200,600);
	c2_phi_Psi[i][j][l]->Divide(4,2);
	for(Int_t m = 0; m < Phi_Psi_total; m++) // phi-psi bin
	{
	  c2_phi_Psi[i][j][l]->cd(m+1);
	  c2_phi_Psi[i][j][l]->cd(m+1)->SetLeftMargin(0.2);
	  c2_phi_Psi[i][j][l]->cd(m+1)->SetTicks(1,1);
	  c2_phi_Psi[i][j][l]->cd(m+1)->SetGrid(0,0);
	  h_mMass_Phi2_SM_New_Bin[i][j][l][m]->SetStats(0);
	  h_mMass_Phi2_SM_New_Bin[i][j][l][m]->Draw("pE");
	}
	c2_phi_Psi[i][j][l]->cd(8);
	c2_phi_Psi[i][j][l]->cd(8)->SetLeftMargin(0.2);
	c2_phi_Psi[i][j][l]->cd(8)->SetTicks(1,1);
	c2_phi_Psi[i][j][l]->cd(8)->SetGrid(0,0);
	h_Counts2[i][j][l]->SetStats(0);
	h_Counts2[i][j][l]->Draw("pE");
	f_phi2[i][j][l]->SetLineColor(2);
	f_phi2[i][j][l]->Draw("l same");

	CanName = Form("c3_phi_Psi_pt_%d_Centrality_%d_EtaGap_%d",i,j,l);
	c3_phi_Psi[i][j][l] = new TCanvas(CanName.Data(),CanName.Data(),1400,10,1200,600);
	c3_phi_Psi[i][j][l]->Divide(4,2);
	for(Int_t m = 0; m < Phi_Psi_total; m++) // phi-psi bin
	{
	  c3_phi_Psi[i][j][l]->cd(m+1);
	  c3_phi_Psi[i][j][l]->cd(m+1)->SetLeftMargin(0.2);
	  c3_phi_Psi[i][j][l]->cd(m+1)->SetTicks(1,1);
	  c3_phi_Psi[i][j][l]->cd(m+1)->SetGrid(0,0);
	  h_mMass_Phi3_SM_New_Bin[i][j][l][m]->SetStats(0);
	  h_mMass_Phi3_SM_New_Bin[i][j][l][m]->Draw("pE");
	}
	c3_phi_Psi[i][j][l]->cd(8);
	c3_phi_Psi[i][j][l]->cd(8)->SetLeftMargin(0.2);
	c3_phi_Psi[i][j][l]->cd(8)->SetTicks(1,1);
	c3_phi_Psi[i][j][l]->cd(8)->SetGrid(0,0);
	h_Counts3[i][j][l]->SetStats(0);
	h_Counts3[i][j][l]->Draw("pE");
	f_phi3[i][j][l]->SetLineColor(2);
	f_phi3[i][j][l]->Draw("l same");
      }
    }
  }
  */

  // fill raw v2 and v3 into histogram
  TH1F *h_flow_2[Centrality_total][EtaGap_total];
  TH1F *h_flow_3[Centrality_total][EtaGap_total];
  for(Int_t j = Centrality_start; j < Centrality_stop; j++) // centrality bin
  {
    for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
    {
      TString HistName;

      HistName = Form("flow_2_Centrality_%d_EtaGap_%d_phi",j,l);
      h_flow_2[j][l] = new TH1F(HistName.Data(),HistName.Data(),100,0.0,3.4);
      for(Int_t i = 0; i < pt_total_New_phi; i++)
      {
	Int_t bin_center = h_flow_2[j][l]->FindBin((pt_low_phi[i]+pt_up_phi[i])/2.0);
	Float_t bin_content = f_phi2[i][j][l]->GetParameter(1);
	Float_t bin_error = f_phi2[i][j][l]->GetParError(1);
	h_flow_2[j][l]->SetBinContent(bin_center,bin_content);
	h_flow_2[j][l]->SetBinError(bin_center,bin_error);
      }

      HistName = Form("flow_3_Centrality_%d_EtaGap_%d_phi",j,l);
      h_flow_3[j][l] = new TH1F(HistName.Data(),HistName.Data(),100,0.0,3.4);
      for(Int_t i = 0; i < pt_total_New_phi; i++)
      {
	Int_t bin_center = h_flow_3[j][l]->FindBin((pt_low_phi[i]+pt_up_phi[i])/2.0);
	Float_t bin_content = f_phi3[i][j][l]->GetParameter(1);
	Float_t bin_error = f_phi3[i][j][l]->GetParError(1);
	h_flow_3[j][l]->SetBinContent(bin_center,bin_content);
	h_flow_3[j][l]->SetBinError(bin_center,bin_error);
      }
    }
  }

  //--------------------------------------------------------------------------------------------------
  // calculate the final resolution correction
  // read the file from Same Event and Mixed Event
  TH1F *h_Yields_Phi2_SE[9][EtaGap_total];
  TH1F *h_Yields_Phi3_SE[9][EtaGap_total];

  TH1F *h_Yields_Phi2_ME[9][EtaGap_total];
  TH1F *h_Yields_Phi3_ME[9][EtaGap_total];
  for(Int_t j = 0; j < 9; j++) // centrality bin
  {
    for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
    {
      TString HistName;

      HistName = Form("Centrality_%d_EtaGap_%d_2nd_Phi_SE",j,l);
      h_Yields_Phi2_SE[j][l] = (TH1F*)File_SE->Get(HistName.Data())->Clone();
      HistName = Form("Centrality_%d_EtaGap_%d_3rd_Phi_SE",j,l);
      h_Yields_Phi3_SE[j][l] = (TH1F*)File_SE->Get(HistName.Data())->Clone();

      HistName = Form("Centrality_%d_EtaGap_%d_2nd_Phi_ME",j,l);
      h_Yields_Phi2_ME[j][l] = (TH1F*)File_ME->Get(HistName.Data())->Clone();
      HistName = Form("Centrality_%d_EtaGap_%d_3rd_Phi_ME",j,l);
      h_Yields_Phi3_ME[j][l] = (TH1F*)File_ME->Get(HistName.Data())->Clone();
    }
  }

  TH1F *h_Yields_Phi2_SM[9][EtaGap_total];
  TH1F *h_Yields_Phi3_SM[9][EtaGap_total];
  for(Int_t j = 0; j < 9; j++) // centrality bin
  {
    for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
    {
      Int_t bin2_SE_start = h_Yields_Phi2_SE[j][l]->FindBin(1.04);
      Int_t bin2_SE_stop  = h_Yields_Phi2_SE[j][l]->FindBin(1.05);
      Float_t Inte2_SE = h_Yields_Phi2_SE[j][l]->Integral(bin2_SE_start,bin2_SE_stop);
      h_Yields_Phi2_SM[j][l] = (TH1F*)h_Yields_Phi2_SE[j][l]->Clone();

      Int_t bin3_SE_start = h_Yields_Phi3_SE[j][l]->FindBin(1.04);
      Int_t bin3_SE_stop  = h_Yields_Phi3_SE[j][l]->FindBin(1.05);
      Float_t Inte3_SE = h_Yields_Phi3_SE[j][l]->Integral(bin3_SE_start,bin3_SE_stop);
      h_Yields_Phi3_SM[j][l] = (TH1F*)h_Yields_Phi3_SE[j][l]->Clone();

      Int_t bin2_ME_start = h_Yields_Phi2_ME[j][l]->FindBin(1.04);
      Int_t bin2_ME_stop  = h_Yields_Phi2_ME[j][l]->FindBin(1.05);
      Float_t Inte2_ME = h_Yields_Phi2_ME[j][l]->Integral(bin2_ME_start,bin2_ME_stop);
      h_Yields_Phi2_ME[j][l]->Scale(Inte2_SE/Inte2_ME);
      h_Yields_Phi2_SM[j][l]->Add(h_Yields_Phi2_ME[j][l],-1.0);

      Int_t bin3_ME_start = h_Yields_Phi3_ME[j][l]->FindBin(1.04);
      Int_t bin3_ME_stop  = h_Yields_Phi3_ME[j][l]->FindBin(1.05);
      Float_t Inte3_ME = h_Yields_Phi3_ME[j][l]->Integral(bin3_ME_start,bin3_ME_stop);
      h_Yields_Phi3_ME[j][l]->Scale(Inte3_SE/Inte3_ME);
      h_Yields_Phi3_SM[j][l]->Add(h_Yields_Phi3_ME[j][l],-1.0);
    }
  }

  TF1  *f_PolyBW_Phi2_Yields[9][EtaGap_total];
  TF1  *f_PolyBW_Phi3_Yields[9][EtaGap_total];

  Float_t ParFit2_Yields[9][EtaGap_total][5];
  Float_t ParFit3_Yields[9][EtaGap_total][5];

  for(Int_t j = 0; j < 9; j++) // centrality bin
  {
    for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
    {
      TString HistName;

      HistName = Form("f_Centrality_%d_EtaGap_%d_Yields_2nd_Phi_SE",j,l);
      f_PolyBW_Phi2_Yields[j][l] = new TF1(HistName.Data(),PolyBreitWigner,0.994,1.05,5);
      for(Int_t n_par = 0; n_par < 5; n_par++)
      {
	f_PolyBW_Phi2_Yields[j][l]->ReleaseParameter(n_par);
      }
      f_PolyBW_Phi2_Yields[j][l]->SetParameter(0,1.019);
      f_PolyBW_Phi2_Yields[j][l]->SetParLimits(0,1.014,1.024);
      f_PolyBW_Phi2_Yields[j][l]->SetParameter(1,0.0055);
      f_PolyBW_Phi2_Yields[j][l]->SetParameter(2,10000);
      f_PolyBW_Phi2_Yields[j][l]->SetParameter(3,-6000);
      f_PolyBW_Phi2_Yields[j][l]->SetParameter(4,0.5);
      f_PolyBW_Phi2_Yields[j][l]->SetRange(0.994,1.05);

      HistName = Form("f_Centrality_%d_EtaGap_%d_Yields_3rd_Phi_SE",j,l);
      f_PolyBW_Phi3_Yields[j][l] = new TF1(HistName.Data(),PolyBreitWigner,0.994,1.05,5);
      for(Int_t n_par = 0; n_par < 4; n_par++)
      {
	f_PolyBW_Phi3_Yields[j][l]->ReleaseParameter(n_par);
      }
      f_PolyBW_Phi3_Yields[j][l]->SetParameter(0,1.019);
      f_PolyBW_Phi3_Yields[j][l]->SetParLimits(0,1.014,1.024);
      f_PolyBW_Phi3_Yields[j][l]->SetParameter(1,0.0055);
      f_PolyBW_Phi3_Yields[j][l]->SetParameter(2,10000);
      f_PolyBW_Phi3_Yields[j][l]->SetParameter(3,-6000);
      f_PolyBW_Phi3_Yields[j][l]->SetParameter(4,0.5);
      f_PolyBW_Phi3_Yields[j][l]->SetRange(0.994,1.05);

      h_Yields_Phi2_SM[j][l]->Fit(f_PolyBW_Phi2_Yields[j][l],"NQR");
      h_Yields_Phi3_SM[j][l]->Fit(f_PolyBW_Phi3_Yields[j][l],"NQR");
      for(Int_t n_par = 0; n_par < 5; n_par++)
      {
	ParFit2_Yields[j][l][n_par] = f_PolyBW_Phi2_Yields[j][l]->GetParameter(n_par);
	ParFit3_Yields[j][l][n_par] = f_PolyBW_Phi3_Yields[j][l]->GetParameter(n_par);
      }
    }
  }

  TF1 *f_Poly_Phi2_Yields[9][EtaGap_total];
  TF1 *f_Poly_Phi3_Yields[9][EtaGap_total];

  for(Int_t j = 0; j < 9; j++) // centrality bin
  {
    for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
    {
      TString FuncName;

      FuncName = Form("f_poly2_Centrality_%d_EtaGap_%d_Yields",j,l);
      f_Poly_Phi2_Yields[j][l] = new TF1(FuncName.Data(),Poly,0.994,1.05,2);
      f_Poly_Phi2_Yields[j][l]->FixParameter(0,ParFit2_Yields[j][l][3]);
      f_Poly_Phi2_Yields[j][l]->FixParameter(1,ParFit2_Yields[j][l][4]);
      h_Yields_Phi2_SM[j][l]->Add(f_Poly_Phi2_Yields[j][l],-1.0);

      FuncName = Form("f_poly3_Centrality_%d_EtaGap_%d_Yields",j,l);
      f_Poly_Phi3_Yields[j][l] = new TF1(FuncName.Data(),Poly,0.994,1.05,2);
      f_Poly_Phi3_Yields[j][l]->FixParameter(0,ParFit3_Yields[j][l][3]);
      f_Poly_Phi3_Yields[j][l]->FixParameter(1,ParFit3_Yields[j][l][4]);
      h_Yields_Phi3_SM[j][l]->Add(f_Poly_Phi3_Yields[j][l],-1.0);
    }
  }

  // calculate total counts and for each centrality bin and eta bin
  Float_t Counts2_Yields[9][EtaGap_total];
  Float_t Counts3_Yields[9][EtaGap_total];
  for(Int_t j = 0; j < 9; j++) // centrality bin
  {
    for(Int_t l = EtaGap_start; l < EtaGap_stop; l++) // eta gap bin
    {
      Counts2_Yields[j][l] = 0.0;
      Int_t bin_start = h_Yields_Phi2_SM[j][l]->FindBin(ParFit2_Yields[j][l][0]-nSigmaPhi*ParFit2_Yields[j][l][1]);
      Int_t bin_stop  = h_Yields_Phi2_SM[j][l]->FindBin(ParFit2_Yields[j][l][0]+nSigmaPhi*ParFit2_Yields[j][l][1]);
      for(Int_t bin = bin_start; bin <= bin_stop; bin++)
      {
	Counts2_Yields[j][l] += h_Yields_Phi2_SM[j][l]->GetBinContent(bin);
      }

      Counts3_Yields[j][l] = 0.0;
      bin_start = h_Yields_Phi3_SM[j][l]->FindBin(ParFit3_Yields[j][l][0]-nSigmaPhi*ParFit3_Yields[j][l][1]);
      bin_stop  = h_Yields_Phi3_SM[j][l]->FindBin(ParFit3_Yields[j][l][0]+nSigmaPhi*ParFit3_Yields[j][l][1]);
      for(Int_t bin = bin_start; bin <= bin_stop; bin++)
      {
	Counts3_Yields[j][l] += h_Yields_Phi3_SM[j][l]->GetBinContent(bin);
      }
    }
  }
}

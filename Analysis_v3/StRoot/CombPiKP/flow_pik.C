#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "draw.h"
#include "TString.h"
#include "TMath.h"
#include "TProfile.h"

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

// mEnergy: 0 for 200GeV, 1 for 39GeV
// mCharge: 0 for positive, 1 for negative
// mCentrality: 0 for 00-80, 1 for 00-10, 2 for 10-40, 3 for 40-80

void flow_pik(Int_t mEnergy = 0, Int_t mCharge = 0, Int_t mCentrality = 0)
{
  const Int_t pt_total = 16;
  const Int_t cent_total = 4; // 4
  const Int_t cent_start[4] = {0,1,2,3};
  const Int_t cent_stop[4]  = {1,2,3,4};
  const Int_t charge_total = 2; // 2
  const Int_t eta_total = 4; // 4
  const Int_t eta_start = 0;
  const Int_t eta_stop = 1;
  const Int_t phi_psi_total = 7; // 7
  const Int_t flow_total = 2; // 2
  TString Order[2] = {"2nd","3rd"};
  Double_t PI_max[2] = {TMath::Pi()/2.0,TMath::Pi()/3.0};
  TString Title[2] = {"#phi-#Psi_{2}","#phi-#Psi_{3}"};
  TString Title_Y[2] = {"v_{2}","v_{3}"};
  TString Title_pion[2] = {"#pi^{+}","#pi^{-}"};
  TString Title_kaon[2] = {"K^{+}","K^{-}"};
  TString Energy[2] = {"200GeV","39GeV"};
  TString Centrality[4] = {"0080","0010","1040","4080"};

  TString inputfile = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/Mass2_nSigmaPion/merged_file/flow_pik/Counts_%s_charge_%d_etagap_%d.root",Energy[mEnergy].Data(),Centrality[mCentrality].Data(),mCharge,eta_start);
  TFile *file_input = TFile::Open(inputfile.Data());

  // pt bin
  Float_t pt_low[pt_total] = {0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4};
  Float_t pt_up[pt_total]  = {0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,4.2};

  Float_t counts_pion[pt_total][cent_total][charge_total][eta_total][phi_psi_total][flow_total];
  Float_t errors_pion[pt_total][cent_total][charge_total][eta_total][phi_psi_total][flow_total];
  Float_t counts_kaon[pt_total][cent_total][charge_total][eta_total][phi_psi_total][flow_total];
  Float_t errors_kaon[pt_total][cent_total][charge_total][eta_total][phi_psi_total][flow_total];

  Float_t flow_pion[pt_total][cent_total][charge_total][eta_total][flow_total];
  Float_t Err_pion[pt_total][cent_total][charge_total][eta_total][flow_total];
  Float_t flow_kaon[pt_total][cent_total][charge_total][eta_total][flow_total];
  Float_t Err_kaon[pt_total][cent_total][charge_total][eta_total][flow_total];

  TH1F *h_counts_pion[pt_total][cent_total][charge_total][eta_total][flow_total]; // counts of pion
  TH1F *h_counts_kaon[pt_total][cent_total][charge_total][eta_total][flow_total]; // counts of kaon
  for(Int_t pt_bin = 0; pt_bin < pt_total; pt_bin++)
  {
    for(Int_t cent = cent_start[mCentrality]; cent < cent_stop[mCentrality]; cent++)
    {
      for(Int_t charge = mCharge; charge < mCharge+1; charge++)
      {
	for(Int_t eta_bin = eta_start; eta_bin < eta_stop; eta_bin++)
	{
	  for(Int_t flow_bin = 0; flow_bin < flow_total; flow_bin++)
	  {
	    TString HistName_pion = Form("pt_%d_Centrality_%d_charge_%d_etagap_%d_%s_counts_pion",pt_bin,cent,charge,eta_bin,Order[flow_bin].Data());
	    h_counts_pion[pt_bin][cent][charge][eta_bin][flow_bin] = (TH1F*)file_input->Get(HistName_pion.Data());
	    TString HistName_kaon = Form("pt_%d_Centrality_%d_charge_%d_etagap_%d_%s_counts_kaon",pt_bin,cent,charge,eta_bin,Order[flow_bin].Data());
	    h_counts_kaon[pt_bin][cent][charge][eta_bin][flow_bin] = (TH1F*)file_input->Get(HistName_kaon.Data());
	    for(Int_t phi_psi_bin = 0; phi_psi_bin < phi_psi_total; phi_psi_bin++)
	    {
	      counts_pion[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin] = h_counts_pion[pt_bin][cent][charge][eta_bin][flow_bin]->GetBinContent(phi_psi_bin+1);
	      errors_pion[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin] = h_counts_pion[pt_bin][cent][charge][eta_bin][flow_bin]->GetBinError(phi_psi_bin+1);
	      counts_kaon[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin] = h_counts_kaon[pt_bin][cent][charge][eta_bin][flow_bin]->GetBinContent(phi_psi_bin+1);
	      errors_kaon[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin] = h_counts_kaon[pt_bin][cent][charge][eta_bin][flow_bin]->GetBinError(phi_psi_bin+1);
	    }
	  }
	}
      }
    }
  }

  TH1F* pos_neg_dummy = new TH1F("pos_neg_dummy","pos_neg_dummy",500,0.,TMath::Pi());
  for(Int_t bin_x = 1; bin_x < pos_neg_dummy->GetNbinsX(); bin_x++)
  {
     pos_neg_dummy->SetBinContent(bin_x,-100.0);
  }
  pos_neg_dummy->SetStats(0);
  pos_neg_dummy->SetTitle("");
  pos_neg_dummy->GetXaxis()->SetTitleOffset(1.0);
  pos_neg_dummy->GetYaxis()->SetTitleOffset(1.2);
  pos_neg_dummy->GetXaxis()->SetLabelSize(0.05);
  pos_neg_dummy->GetYaxis()->SetLabelSize(0.05);
  pos_neg_dummy->GetXaxis()->SetTitleSize(0.06);
  pos_neg_dummy->GetYaxis()->SetTitleSize(0.06);
  pos_neg_dummy->GetXaxis()->SetNdivisions(505,'N');
  pos_neg_dummy->GetYaxis()->SetNdivisions(505,'N');
  pos_neg_dummy->GetXaxis()->CenterTitle();
  pos_neg_dummy->GetYaxis()->CenterTitle();

  // fit for pion and kaon
  TF1 *f_pion[pt_total][cent_total][charge_total][eta_total][flow_total];
  TF1 *f_kaon[pt_total][cent_total][charge_total][eta_total][flow_total];
  for(Int_t pt_bin = 0; pt_bin < pt_total; pt_bin++)
  {
    for(Int_t cent = cent_start[mCentrality]; cent < cent_stop[mCentrality]; cent++)
    {
      for(Int_t charge = mCharge; charge < mCharge+1; charge++)
      {
	for(Int_t eta_bin = eta_start; eta_bin < eta_stop; eta_bin++)
	{
	  for(Int_t flow_bin = 0; flow_bin < flow_total; flow_bin++)
	  {
	    if(flow_bin == 0)
	    {
	      TString Flow_pion = Form("flow_pt_%d_Centrality_%d_charge_%d_etagap_%d_%s_pion",pt_bin,cent,charge,eta_bin,Order[flow_bin].Data());
	      f_pion[pt_bin][cent][charge][eta_bin][flow_bin] = new TF1(Flow_pion.Data(),flow_2,0.0,PI_max[flow_bin],2);
	      f_pion[pt_bin][cent][charge][eta_bin][flow_bin]->SetParameter(0,10000);
	      f_pion[pt_bin][cent][charge][eta_bin][flow_bin]->SetParameter(1,0.2);
	      TString Flow_kaon = Form("flow_pt_%d_Centrality_%d_charge_%d_etagap_%d_%s_kaon",pt_bin,cent,charge,eta_bin,Order[flow_bin].Data());
	      f_kaon[pt_bin][cent][charge][eta_bin][flow_bin] = new TF1(Flow_kaon.Data(),flow_2,0.0,PI_max[flow_bin],2);
	      f_kaon[pt_bin][cent][charge][eta_bin][flow_bin]->SetParameter(0,1000);
	      f_kaon[pt_bin][cent][charge][eta_bin][flow_bin]->SetParameter(1,0.1);
	    }
	    if(flow_bin == 1)
	    {
	      TString Flow_pion = Form("flow_pt_%d_Centrality_%d_charge_%d_etagap_%d_%s_pion",pt_bin,cent,charge,eta_bin,Order[flow_bin].Data());
	      f_pion[pt_bin][cent][charge][eta_bin][flow_bin] = new TF1(Flow_pion.Data(),flow_3,0.0,PI_max[flow_bin],2);
	      f_pion[pt_bin][cent][charge][eta_bin][flow_bin]->SetParameter(0,1000);
	      f_pion[pt_bin][cent][charge][eta_bin][flow_bin]->SetParameter(1,0.1);
	      TString Flow_kaon = Form("flow_pt_%d_Centrality_%d_charge_%d_etagap_%d_%s_kaon",pt_bin,cent,charge,eta_bin,Order[flow_bin].Data());
	      f_kaon[pt_bin][cent][charge][eta_bin][flow_bin] = new TF1(Flow_kaon.Data(),flow_3,0.0,PI_max[flow_bin],2);
	      f_kaon[pt_bin][cent][charge][eta_bin][flow_bin]->SetParameter(0,1000);
	      f_kaon[pt_bin][cent][charge][eta_bin][flow_bin]->SetParameter(1,0.2);
	    }
	  }
	}
      }
    }
  }

  TCanvas *c_counts_pion[pt_total][cent_total][charge_total][eta_total][flow_total]; // counts of pion
  TCanvas *c_counts_kaon[pt_total][cent_total][charge_total][eta_total][flow_total]; // counts of kaon
  for(Int_t pt_bin = 0; pt_bin < pt_total; pt_bin++)
  {
    for(Int_t cent = cent_start[mCentrality]; cent < cent_stop[mCentrality]; cent++)
    {
      for(Int_t charge = mCharge; charge < mCharge+1; charge++)
      {
	for(Int_t eta_bin = eta_start; eta_bin < eta_stop; eta_bin++)
	{
	  for(Int_t flow_bin = 0; flow_bin < flow_total; flow_bin++)
	  {
	    TString CanName_pion = Form("c_pt_%d_Centrality_%d_charge_%d_etagap_%d_%s_counts_pion",pt_bin,cent,charge,eta_bin,Order[flow_bin].Data());
	    c_counts_pion[pt_bin][cent][charge][eta_bin][flow_bin] = new TCanvas(CanName_pion.Data(),CanName_pion.Data(),1400,10,800,800);
	    c_counts_pion[pt_bin][cent][charge][eta_bin][flow_bin]->cd();
	    c_counts_pion[pt_bin][cent][charge][eta_bin][flow_bin]->cd()->SetLeftMargin(0.15);
	    c_counts_pion[pt_bin][cent][charge][eta_bin][flow_bin]->cd()->SetBottomMargin(0.15);
	    c_counts_pion[pt_bin][cent][charge][eta_bin][flow_bin]->cd()->SetTicks(1,1);
	    c_counts_pion[pt_bin][cent][charge][eta_bin][flow_bin]->cd()->SetGrid(0,0);
	    pos_neg_dummy->GetXaxis()->SetTitle(Title[flow_bin].Data());
	    pos_neg_dummy->GetXaxis()->SetRangeUser(0.0,PI_max[flow_bin]);
	    pos_neg_dummy->GetYaxis()->SetTitle("Counts");
	    pos_neg_dummy->GetYaxis()->SetRangeUser(counts_pion[pt_bin][cent][charge][eta_bin][6][flow_bin]*0.99,counts_pion[pt_bin][cent][charge][eta_bin][0][flow_bin]*1.01);
	    pos_neg_dummy->DrawCopy("pE");
	    h_counts_pion[pt_bin][cent][charge][eta_bin][flow_bin]->Fit(f_pion[pt_bin][cent][charge][eta_bin][flow_bin],"QI");
	    flow_pion[pt_bin][cent][charge][eta_bin][flow_bin] = f_pion[pt_bin][cent][charge][eta_bin][flow_bin]->GetParameter(1);
	    Err_pion[pt_bin][cent][charge][eta_bin][flow_bin] = f_pion[pt_bin][cent][charge][eta_bin][flow_bin]->GetParError(1);
	    h_counts_pion[pt_bin][cent][charge][eta_bin][flow_bin]->Draw("pE same");

	    TString CanName_kaon = Form("c_pt_%d_Centrality_%d_charge_%d_etagap_%d_%s_counts_kaon",pt_bin,cent,charge,eta_bin,Order[flow_bin].Data());
	    c_counts_kaon[pt_bin][cent][charge][eta_bin][flow_bin] = new TCanvas(CanName_kaon.Data(),CanName_kaon.Data(),1400,10,800,800);
	    c_counts_kaon[pt_bin][cent][charge][eta_bin][flow_bin]->cd();
	    c_counts_kaon[pt_bin][cent][charge][eta_bin][flow_bin]->cd()->SetLeftMargin(0.15);
	    c_counts_kaon[pt_bin][cent][charge][eta_bin][flow_bin]->cd()->SetBottomMargin(0.15);
	    c_counts_kaon[pt_bin][cent][charge][eta_bin][flow_bin]->cd()->SetTicks(1,1);
	    c_counts_kaon[pt_bin][cent][charge][eta_bin][flow_bin]->cd()->SetGrid(0,0);
	    pos_neg_dummy->GetXaxis()->SetTitle(Title[flow_bin].Data());
	    pos_neg_dummy->GetXaxis()->SetRangeUser(0.0,PI_max[flow_bin]);
	    pos_neg_dummy->GetYaxis()->SetTitle("Counts");
	    pos_neg_dummy->GetYaxis()->SetRangeUser(counts_kaon[pt_bin][cent][charge][eta_bin][6][flow_bin]*0.99,counts_kaon[pt_bin][cent][charge][eta_bin][0][flow_bin]*1.01);
	    pos_neg_dummy->DrawCopy("pE");
	    h_counts_kaon[pt_bin][cent][charge][eta_bin][flow_bin]->Fit(f_kaon[pt_bin][cent][charge][eta_bin][flow_bin],"QI");
	    flow_kaon[pt_bin][cent][charge][eta_bin][flow_bin] = f_kaon[pt_bin][cent][charge][eta_bin][flow_bin]->GetParameter(1);
	    Err_kaon[pt_bin][cent][charge][eta_bin][flow_bin] = f_kaon[pt_bin][cent][charge][eta_bin][flow_bin]->GetParError(1);
	    h_counts_kaon[pt_bin][cent][charge][eta_bin][flow_bin]->Draw("pE same");
	  }
	}
      }
    }
  }

  // Resolution Correction
  Float_t mean_res_2_pion[charge_total][eta_total][2][cent_total]; // Energy: 0 = 200 GeV, 1 = 39 GeV
  Float_t mean_res_2_kaon[charge_total][eta_total][2][cent_total];
  Float_t mean_res_3_pion[charge_total][eta_total][2][cent_total];
  Float_t mean_res_3_kaon[charge_total][eta_total][2][cent_total];

  /*
  // 200 GeV
  // pos charge
  mean_res_2_pion[0][0][0] = 1.92072;
  mean_res_2_kaon[0][0][0] = 1.91034;
  mean_res_3_pion[0][0][0] = 4.16207;
  mean_res_3_kaon[0][0][0] = 4.07798;

  // neg charge
  mean_res_2_pion[1][0][0] = 1.92071;
  mean_res_2_kaon[1][0][0] = 1.91049;
  mean_res_3_pion[1][0][0] = 4.16049;
  mean_res_3_kaon[1][0][0] = 4.08106;
  */
  //
  // 200 GeV
  mean_res_2_pion[0][0][0][0] = 1.84219; // pos, 0.05, 200GeV, 00-80
  mean_res_2_kaon[0][0][0][0] = 1.82832; // pos, 0.05, 200GeV, 00-80
  mean_res_3_pion[0][0][0][0] = 4.03529; // pos, 0.05, 200GeV, 00-80
  mean_res_3_kaon[0][0][0][0] = 3.93780; // pos, 0.05, 200GeV, 00-80

  mean_res_2_pion[1][0][0][0] = 1.84208; // neg, 0.05, 200GeV, 00-80
  mean_res_2_kaon[1][0][0][0] = 1.82856; // neg, 0.05, 200GeV, 00-80
  mean_res_3_pion[1][0][0][0] = 4.03239; // neg, 0.05, 200GeV, 00-80
  mean_res_3_kaon[1][0][0][0] = 3.94038; // neg, 0.05, 200GeV, 00-80
  
  mean_res_2_pion[0][0][0][2] = 1.53783; // pos, 0.05, 200GeV, 10-40
  mean_res_2_kaon[0][0][0][2] = 1.53776; // pos, 0.05, 200GeV, 10-40
  mean_res_3_pion[0][0][0][2] = 3.49234; // pos, 0.05, 200GeV, 10-40
  mean_res_3_kaon[0][0][0][2] = 3.48602; // pos, 0.05, 200GeV, 10-40

  mean_res_2_pion[1][0][0][2] = 1.53783; // neg, 0.05, 200GeV, 10-40
  mean_res_2_kaon[1][0][0][2] = 1.53777; // neg, 0.05, 200GeV, 10-40
  mean_res_3_pion[1][0][0][2] = 3.49194; // neg, 0.05, 200GeV, 10-40
  mean_res_3_kaon[1][0][0][2] = 3.48631; // neg, 0.05, 200GeV, 10-40

  // 39 GeV
  mean_res_2_pion[0][0][1][0] = 2.52810; // pos, 0.05, 39GeV, 00-80
  mean_res_2_kaon[0][0][1][0] = 2.49843; // pos, 0.05, 39GeV, 00-80
  mean_res_3_pion[0][0][1][0] = 7.11682; // pos, 0.05, 39GeV, 00-80
  mean_res_3_kaon[0][0][1][0] = 6.75552; // pos, 0.05, 39GeV, 00-80
                                                                   
  mean_res_2_pion[1][0][1][0] = 2.52806; // neg, 0.05, 39GeV, 00-80
  mean_res_2_kaon[1][0][1][0] = 2.50093; // neg, 0.05, 39GeV, 00-80
  mean_res_3_pion[1][0][1][0] = 7.11947; // neg, 0.05, 39GeV, 00-80
  mean_res_3_kaon[1][0][1][0] = 6.80734; // neg, 0.05, 39GeV, 00-80

  mean_res_2_pion[0][0][1][2] = 2.02130; // pos, 0.05, 39GeV, 10-40
  mean_res_2_kaon[0][0][1][2] = 2.02092; // pos, 0.05, 39GeV, 10-40
  mean_res_3_pion[0][0][1][2] = 5.78088; // pos, 0.05, 39GeV, 10-40
  mean_res_3_kaon[0][0][1][2] = 5.76030; // pos, 0.05, 39GeV, 10-40
                                                                   
  mean_res_2_pion[1][0][1][2] = 2.02130; // neg, 0.05, 39GeV, 10-40
  mean_res_2_kaon[1][0][1][2] = 2.02099; // neg, 0.05, 39GeV, 10-40
  mean_res_3_pion[1][0][1][2] = 5.78104; // neg, 0.05, 39GeV, 10-40
  mean_res_3_kaon[1][0][1][2] = 5.76556; // neg, 0.05, 39GeV, 10-40

  //flow vs pt
  TGraphAsymmErrors *g_flow_pion[cent_total][charge_total][eta_total][flow_total]; // flow of pion vs pt
  TGraphAsymmErrors *g_flow_kaon[cent_total][charge_total][eta_total][flow_total]; // flow of kaon vs pt
  TH1F *h_flow_pion[cent_total][charge_total][eta_total][flow_total]; // flow of pion vs pt
  TH1F *h_flow_kaon[cent_total][charge_total][eta_total][flow_total]; // flow of kaon vs pt
  for(Int_t cent = cent_start[mCentrality]; cent < cent_stop[mCentrality]; cent++)
  {
    for(Int_t charge = mCharge; charge < mCharge+1; charge++)
    {
      for(Int_t eta_bin = eta_start; eta_bin < eta_stop; eta_bin++)
      {
	for(Int_t flow_bin = 0; flow_bin < flow_total; flow_bin++)
	{
	  TString g_Name_pion = Form("g_Centrality_%d_charge_%d_etagap_%d_%s_flow_pion",cent,charge,eta_bin,Order[flow_bin].Data());
	  TString g_Name_kaon = Form("g_Centrality_%d_charge_%d_etagap_%d_%s_flow_kaon",cent,charge,eta_bin,Order[flow_bin].Data());
	  g_flow_pion[cent][charge][eta_bin][flow_bin] = new TGraphAsymmErrors();
	  g_flow_pion[cent][charge][eta_bin][flow_bin]->SetName(g_Name_pion.Data());
	  g_flow_kaon[cent][charge][eta_bin][flow_bin] = new TGraphAsymmErrors();
	  g_flow_kaon[cent][charge][eta_bin][flow_bin]->SetName(g_Name_kaon.Data());
	  for(Int_t pt_bin = 0; pt_bin < pt_total; pt_bin++)
	  {
	    if(flow_bin == 0)
	    {
	      g_flow_pion[cent][charge][eta_bin][flow_bin]->SetPoint(pt_bin,(pt_low[pt_bin]+pt_up[pt_bin])/2.0,flow_pion[pt_bin][cent][charge][eta_bin][flow_bin]*mean_res_2_pion[charge][eta_bin][mEnergy][mCentrality]);
	      g_flow_pion[cent][charge][eta_bin][flow_bin]->SetPointError(pt_bin,0.0,0.0,Err_pion[pt_bin][cent][charge][eta_bin][flow_bin]*mean_res_2_pion[charge][eta_bin][mEnergy][mCentrality],Err_pion[pt_bin][cent][charge][eta_bin][flow_bin]*mean_res_2_pion[charge][eta_bin][mEnergy][mCentrality]);

	      g_flow_kaon[cent][charge][eta_bin][flow_bin]->SetPoint(pt_bin,(pt_low[pt_bin]+pt_up[pt_bin])/2.0,flow_kaon[pt_bin][cent][charge][eta_bin][flow_bin]*mean_res_2_kaon[charge][eta_bin][mEnergy][mCentrality]);
	      g_flow_kaon[cent][charge][eta_bin][flow_bin]->SetPointError(pt_bin,0.0,0.0,Err_kaon[pt_bin][cent][charge][eta_bin][flow_bin]*mean_res_2_kaon[charge][eta_bin][mEnergy][mCentrality],Err_kaon[pt_bin][cent][charge][eta_bin][flow_bin]*mean_res_2_kaon[charge][eta_bin][mEnergy][mCentrality]);
	    }
	    if(flow_bin == 1)
	    {
	      g_flow_pion[cent][charge][eta_bin][flow_bin]->SetPoint(pt_bin,(pt_low[pt_bin]+pt_up[pt_bin])/2.0,flow_pion[pt_bin][cent][charge][eta_bin][flow_bin]*mean_res_3_pion[charge][eta_bin][mEnergy][mCentrality]);
	      g_flow_pion[cent][charge][eta_bin][flow_bin]->SetPointError(pt_bin,0.0,0.0,Err_pion[pt_bin][cent][charge][eta_bin][flow_bin]*mean_res_3_pion[charge][eta_bin][mEnergy][mCentrality],Err_pion[pt_bin][cent][charge][eta_bin][flow_bin]*mean_res_3_pion[charge][eta_bin][mEnergy][mCentrality]);

	      g_flow_kaon[cent][charge][eta_bin][flow_bin]->SetPoint(pt_bin,(pt_low[pt_bin]+pt_up[pt_bin])/2.0,flow_kaon[pt_bin][cent][charge][eta_bin][flow_bin]*mean_res_3_kaon[charge][eta_bin][mEnergy][mCentrality]);
	      g_flow_kaon[cent][charge][eta_bin][flow_bin]->SetPointError(pt_bin,0.0,0.0,Err_kaon[pt_bin][cent][charge][eta_bin][flow_bin]*mean_res_3_kaon[charge][eta_bin][mEnergy][mCentrality],Err_kaon[pt_bin][cent][charge][eta_bin][flow_bin]*mean_res_3_kaon[charge][eta_bin][mEnergy][mCentrality]);
	    }
	  }
	}
      }
    }
  }
  // save flow
  TString output_pik = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/Mass2_nSigmaPion/merged_file/flow_pik/Flow_pik_%s_charge_%d_etagap_%d.root",Energy[mEnergy].Data(),Centrality[mCentrality].Data(),mCharge,eta_start);
  TFile *file_output = new TFile(output_pik.Data(),"RECREATE");
  file_output->cd();
  for(Int_t cent = cent_start[mCentrality]; cent < cent_stop[mCentrality]; cent++)
  {
    for(Int_t charge = mCharge; charge < mCharge+1; charge++)
    {
      for(Int_t eta_bin = eta_start; eta_bin < eta_stop; eta_bin++)
      {
	for(Int_t flow_bin = 0; flow_bin < flow_total; flow_bin++)
	{
	  g_flow_pion[cent][charge][eta_bin][flow_bin]->Write();
	  g_flow_kaon[cent][charge][eta_bin][flow_bin]->Write();
	}
      }
    }
  }
  file_output->Close();

  /*
  // figures for rehearsal
  TCanvas *c_flow2_pion[2], *c_flow2_kaon[2]; // 0 = particles, 1 = anti-particles
  for(Int_t charge = mCharge; charge < mCharge+1; charge++)
  {
    TString CanName_pion = Form("c_flow2_%s_charge_%d_etagap_%d_pion",Centrality[mCentrality].Data(),charge,eta_start);
    c_flow2_pion[charge] = new TCanvas(CanName_pion.Data(),CanName_pion.Data(),10,10,800,800);
    c_flow2_pion[charge]->cd();
    c_flow2_pion[charge]->cd()->SetLeftMargin(0.15);
    c_flow2_pion[charge]->cd()->SetBottomMargin(0.15);
    c_flow2_pion[charge]->cd()->SetTicks(1,1);
    c_flow2_pion[charge]->cd()->SetGrid(0,0);
    pos_neg_dummy->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    pos_neg_dummy->GetYaxis()->SetTitle(Title_Y[0]);
    pos_neg_dummy->GetXaxis()->SetRangeUser(0.0,3.5);
    pos_neg_dummy->GetYaxis()->SetRangeUser(0.0,0.2);
    pos_neg_dummy->DrawCopy("pE");
    Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_flow_pion[0][charge][eta_start][0],20,2,0.8);
    if(eta_start == 0)
    {
      Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_pion_alex[charge],20,1,0.8);
    }
    plotTopLegend((char*)Title_pion[charge].Data(),0.2,0.8,0.06,1,0.0,42,1);
    Draw_TGAE_Point_new_Symbol(0.2,0.15,0.0,0.0,0.0,0.0,20,2,0.8);
    plotTopLegend((char*)"My result",0.25,0.148,0.03,1,0.0,42,0);
    Draw_TGAE_Point_new_Symbol(0.2,0.14,0.0,0.0,0.0,0.0,20,1,0.8);
    plotTopLegend((char*)"Phys. Rev. C 88, 014902 (2013)",0.25,0.138,0.03,1,0.0,42,0);
    c_flow2_pion[charge]->SaveAs((Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/Mass2_nSigmaPion/merged_file/figures/%s.gif",Energy[mEnergy].Data(),CanName_pion.Data())).Data());

    TString CanName_kaon = Form("c_flow2_%s_charge_%d_etagap_%d_kaon",Centrality[mCentrality].Data(),charge,eta_start);
    c_flow2_kaon[charge] = new TCanvas(CanName_kaon.Data(),CanName_kaon.Data(),10,10,800,800);
    c_flow2_kaon[charge]->cd();
    c_flow2_kaon[charge]->cd()->SetLeftMargin(0.15);
    c_flow2_kaon[charge]->cd()->SetBottomMargin(0.15);
    c_flow2_kaon[charge]->cd()->SetTicks(1,1);
    c_flow2_kaon[charge]->cd()->SetGrid(0,0);
    pos_neg_dummy->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    pos_neg_dummy->GetYaxis()->SetTitle(Title_Y[0]);
    pos_neg_dummy->GetXaxis()->SetRangeUser(0.0,3.5);
    pos_neg_dummy->GetYaxis()->SetRangeUser(0.0,0.2);
    pos_neg_dummy->DrawCopy("pE");
    Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_flow_kaon[0][charge][eta_start][0],20,2,0.8);
    if(eta_start == 0)
    {
      Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_kaon_alex[charge],20,1,0.8);
    }
    plotTopLegend((char*)Title_kaon[charge].Data(),0.2,0.8,0.06,1,0.0,42,1);
    Draw_TGAE_Point_new_Symbol(0.2,0.15,0.0,0.0,0.0,0.0,20,2,0.8);
    plotTopLegend((char*)"My result",0.25,0.148,0.03,1,0.0,42,0);
    Draw_TGAE_Point_new_Symbol(0.2,0.14,0.0,0.0,0.0,0.0,20,1,0.8);
    plotTopLegend((char*)"Phys. Rev. C 88, 014902 (2013)",0.25,0.138,0.03,1,0.0,42,0);
    c_flow2_kaon[charge]->SaveAs(("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu200GeV/Mass2_nSigmaPion/merged_file/figures/"+CanName_kaon+".gif").Data());
  }
  
  TCanvas *c_flow3_pik[2], *c_flow3_pion, *c_flow3_kaon;
  for(Int_t charge = mCharge; charge < mCharge+1; charge++)
  {
    TString CanName_pik = Form("c_flow3_%s_charge_%d_etagap_%d_pik",Centrality[mCentrality].Data(),charge,eta_start);
    c_flow3_pik[charge] = new TCanvas(CanName_pik.Data(),CanName_pik.Data(),10,10,800,800);
    c_flow3_pik[charge]->cd();
    c_flow3_pik[charge]->cd()->SetLeftMargin(0.15);
    c_flow3_pik[charge]->cd()->SetBottomMargin(0.15);
    c_flow3_pik[charge]->cd()->SetTicks(1,1);
    c_flow3_pik[charge]->cd()->SetGrid(0,0);
    pos_neg_dummy->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    pos_neg_dummy->GetYaxis()->SetTitle(Title_Y[1]);
    pos_neg_dummy->GetXaxis()->SetRangeUser(0.0,3.5);
    pos_neg_dummy->GetYaxis()->SetRangeUser(0.0,0.1);
    pos_neg_dummy->DrawCopy("pE");
    Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_flow_pion[0][charge][eta_start][1],20,2,0.8);
    Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_flow_kaon[0][charge][eta_start][1],28,4,0.8);
    Draw_TGAE_Point_new_Symbol(0.2,0.09,0.0,0.0,0.0,0.0,20,2,0.8);
    plotTopLegend((char*)Title_pion[charge].Data(),0.25,0.088,0.03,1,0.0,42,0);
    Draw_TGAE_Point_new_Symbol(0.2,0.08,0.0,0.0,0.0,0.0,28,4,0.8);
    plotTopLegend((char*)Title_kaon[charge].Data(),0.25,0.078,0.03,1,0.0,42,0);
    c_flow3_pik[charge]->SaveAs(("./figures/"+CanName_pik+".gif").Data());
  }

  TString CanName_pion3 = Form("c_flow3_%s_etagap_%d_pion",Centrality[mCentrality].Data(),eta_start);
  c_flow3_pion = new TCanvas(CanName_pion3.Data(),CanName_pion3.Data(),10,10,800,800);
  c_flow3_pion->cd();
  c_flow3_pion->cd()->SetLeftMargin(0.15);
  c_flow3_pion->cd()->SetBottomMargin(0.15);
  c_flow3_pion->cd()->SetTicks(1,1);
  c_flow3_pion->cd()->SetGrid(0,0);
  pos_neg_dummy->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  pos_neg_dummy->GetYaxis()->SetTitle(Title_Y[1]);
  pos_neg_dummy->GetXaxis()->SetRangeUser(0.0,3.5);
  pos_neg_dummy->GetYaxis()->SetRangeUser(0.0,0.1);
  pos_neg_dummy->DrawCopy("pE");
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_flow_pion[0][0][eta_start][1],20,2,0.8);
  Draw_TGAE_Point_new_Symbol(0.2,0.09,0.0,0.0,0.0,0.0,20,2,0.8);
  plotTopLegend((char*)Title_pion[0].Data(),0.25,0.088,0.03,1,0.0,42,0);
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_flow_pion[0][1][eta_start][1],24,1,0.8);
  Draw_TGAE_Point_new_Symbol(0.2,0.08,0.0,0.0,0.0,0.0,24,1,0.8);
  plotTopLegend((char*)Title_pion[1].Data(),0.25,0.078,0.03,1,0.0,42,0);
  c_flow3_pion->SaveAs(("./figures/"+CanName_pion3+".gif"));

  TString CanName_kaon3 = Form("c_flow3_%s_etagap_%d_kaon",Centrality[mCentrality].Data(),eta_start);
  c_flow3_kaon = new TCanvas(CanName_kaon3.Data(),CanName_kaon3.Data(),10,10,800,800);
  c_flow3_kaon->cd();
  c_flow3_kaon->cd()->SetLeftMargin(0.15);
  c_flow3_kaon->cd()->SetBottomMargin(0.15);
  c_flow3_kaon->cd()->SetTicks(1,1);
  c_flow3_kaon->cd()->SetGrid(0,0);
  pos_neg_dummy->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  pos_neg_dummy->GetYaxis()->SetTitle(Title_Y[1]);
  pos_neg_dummy->GetXaxis()->SetRangeUser(0.0,3.5);
  pos_neg_dummy->GetYaxis()->SetRangeUser(0.0,0.1);
  pos_neg_dummy->DrawCopy("pE");
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_flow_kaon[0][0][eta_start][1],20,2,0.8);
  Draw_TGAE_Point_new_Symbol(0.2,0.09,0.0,0.0,0.0,0.0,20,2,0.8);
  plotTopLegend((char*)Title_kaon[0].Data(),0.25,0.088,0.03,1,0.0,42,0);
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_flow_kaon[0][1][eta_start][1],24,1,0.8);
  Draw_TGAE_Point_new_Symbol(0.2,0.08,0.0,0.0,0.0,0.0,24,1,0.8);
  plotTopLegend((char*)Title_kaon[1].Data(),0.25,0.078,0.03,1,0.0,42,0);
  c_flow3_kaon->SaveAs(("./figures/"+CanName_kaon3+".gif"));
  */
}

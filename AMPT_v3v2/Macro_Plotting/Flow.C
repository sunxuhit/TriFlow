#include "TString.h"
#include "TFile.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"

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

static TString Mode[2] = {"Default","StringMelting"};
static TString Energy[7] = {"7GeV","11GeV","19GeV","27GeV","39GeV","62GeV","200GeV"};
static TString Order[2] = {"2nd","3rd"};
static TString Centrality[4] = {"0080","0010","1040","4080"};
static TString ParType[10] = {"pi_plus","pi_minus","K_plus","K_minus","p","pbar","Lambda","Lambdabar","K0s","phi"};

void Flow(Int_t mEnergy = 4, Int_t mMode = 1) // 0: 7.7 GeV, 1: 11.5 GeV, 2: 19.6 GeV, 3: 27 GeV, 4: 39 GeV, 5: 62.4 GeV, 6: 200 GeV | 0: Default, 1: String Melting
{
  gStyle->SetTitleX(0.53);
  gStyle->SetTitleY(0.96);
  TString inputfile = Form("/home/xusun/Data/AMPT_%s/Flow/%s_%s/Flow_%s.root",Mode[mMode].Data(),Energy[mEnergy].Data(),Mode[mMode].Data(),Energy[mEnergy].Data()); // temperory file
  cout << "Input File: " << inputfile.Data() << endl;
  TFile *File_input = TFile::Open(inputfile.Data());

  // flow for pi, K, p, Lambda, K0s by using eta_sub event plane method
  // pi_plus,pi_minus,K_plus,K_minus,p,pbar,Lambda,Lambdabar,K0s,phi
  TProfile *p_mFlow[10][2][4]; // Particle types | 0 for 2nd, 1 for 3rd | 0 for 0-80%, 1 for 0-10%, 2 for 10-40%, 3 for 40-80%

  // pt spectra 
  TH1F *h_mPt[10][2][4]; // Particle types | 0 for 2nd, 1 for 3rd | 0 for 0-80%, 1 for 0-10%, 2 for 10-40%, 3 for 40-80%

  // Read in flow TProfile => transfer to TH1F | Read in pt Spectra
  for(Int_t i_par = 0; i_par < 10; i_par++)
  {
    for(Int_t i_order = 0; i_order < 2; i_order++)
    {
      for(Int_t i_cent = 0; i_cent < 4; i_cent++)
      {
	// v2 and v3 relative to event plane
	//--------------------------------------------------------------------------------------------------
	TString ProName = Form("Flow_%s_%s_%s",ParType[i_par].Data(),Order[i_order].Data(),Centrality[i_cent].Data()); // pi_plus
	p_mFlow[i_par][i_order][i_cent] = (TProfile*)File_input->Get(ProName.Data());
	TString HistName = Form("Pt_%s_%s_%s",ParType[i_par].Data(),Order[i_order].Data(),Centrality[i_cent].Data());
	h_mPt[i_par][i_order][i_cent] = (TH1F*)File_input->Get(HistName.Data());
	//--------------------------------------------------------------------------------------------------
      }
    }
  }

  // v2 of pi+ and p
  TCanvas *c_flow_2nd = new TCanvas("c_flow_2nd","c_flow_2nd",10,10,800,800);
  c_flow_2nd->cd();
  c_flow_2nd->cd()->SetLeftMargin(0.15);
  c_flow_2nd->cd()->SetBottomMargin(0.15);
  c_flow_2nd->cd()->SetTicks(1,1);
  TString Title = Form("Au+Au, %s, 0%%-70%%, AMPT_%s",Energy[mEnergy].Data(),Mode[mMode].Data());
  p_mFlow[4][0][0]->SetTitle(Title.Data());
  p_mFlow[4][0][0]->SetStats(0);
  p_mFlow[4][0][0]->SetLineColor(1);
  p_mFlow[4][0][0]->SetMarkerStyle(24);
  p_mFlow[4][0][0]->SetMarkerColor(2);
  p_mFlow[4][0][0]->SetMarkerSize(1.4);
  p_mFlow[4][0][0]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  p_mFlow[4][0][0]->GetYaxis()->SetTitle("v_{2}");
  p_mFlow[4][0][0]->GetXaxis()->CenterTitle();
  p_mFlow[4][0][0]->GetYaxis()->CenterTitle();
  p_mFlow[4][0][0]->GetXaxis()->SetTitleSize(0.06);
  p_mFlow[4][0][0]->GetYaxis()->SetTitleSize(0.06);
  p_mFlow[4][0][0]->GetXaxis()->SetRangeUser(0.0,3.0);
  p_mFlow[4][0][0]->GetYaxis()->SetRangeUser(-0.01,0.2);
  p_mFlow[4][0][0]->GetXaxis()->SetNdivisions(505,'X');
  p_mFlow[4][0][0]->GetYaxis()->SetNdivisions(505,'Y');
  p_mFlow[4][0][0]->Draw("pEX0");

  p_mFlow[0][0][0]->SetTitle("");
  p_mFlow[0][0][0]->SetStats(0);
  p_mFlow[0][0][0]->SetLineColor(1);
  p_mFlow[0][0][0]->SetMarkerStyle(24);
  p_mFlow[0][0][0]->SetMarkerColor(1);
  p_mFlow[0][0][0]->SetMarkerSize(1.4);
  p_mFlow[0][0][0]->Draw("pEX0 same");
  PlotLine(0.0,3.0,0.0,0.0,1,2,2);

  TLegend *leg = new TLegend(0.2,0.7,0.4,0.85);
  leg->SetFillColor(10);
  leg->SetBorderSize(0.0);
  leg->AddEntry(p_mFlow[4][0][0],"p","P");
  leg->AddEntry(p_mFlow[0][0][0],"#pi^{+}","P");
  leg->Draw("same");

  // v3 of pi+ and p
  TCanvas *c_flow_3rd = new TCanvas("c_flow_3rd","c_flow_3rd",10,10,800,800);
  c_flow_3rd->cd();
  c_flow_3rd->cd()->SetLeftMargin(0.15);
  c_flow_3rd->cd()->SetBottomMargin(0.15);
  c_flow_3rd->cd()->SetTicks(1,1);
  TString Title = Form("Au+Au, %s, 0%%-70%%, AMPT_%s",Energy[mEnergy].Data(),Mode[mMode].Data());
  p_mFlow[4][1][0]->SetTitle(Title.Data());
  p_mFlow[4][1][0]->SetStats(0);
  p_mFlow[4][1][0]->SetLineColor(1);
  p_mFlow[4][1][0]->SetMarkerStyle(24);
  p_mFlow[4][1][0]->SetMarkerColor(2);
  p_mFlow[4][1][0]->SetMarkerSize(1.4);
  p_mFlow[4][1][0]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  p_mFlow[4][1][0]->GetYaxis()->SetTitle("v_{3}");
  p_mFlow[4][1][0]->GetXaxis()->CenterTitle();
  p_mFlow[4][1][0]->GetYaxis()->CenterTitle();
  p_mFlow[4][1][0]->GetXaxis()->SetTitleSize(0.06);
  p_mFlow[4][1][0]->GetYaxis()->SetTitleSize(0.06);
  p_mFlow[4][1][0]->GetXaxis()->SetRangeUser(0.0,3.0);
  p_mFlow[4][1][0]->GetYaxis()->SetRangeUser(-0.01,0.2);
  p_mFlow[4][1][0]->GetXaxis()->SetNdivisions(505,'X');
  p_mFlow[4][1][0]->GetYaxis()->SetNdivisions(505,'Y');
  p_mFlow[4][1][0]->Draw("pEX0");

  p_mFlow[0][1][0]->SetTitle("");
  p_mFlow[0][1][0]->SetStats(0);
  p_mFlow[0][1][0]->SetLineColor(1);
  p_mFlow[0][1][0]->SetMarkerStyle(24);
  p_mFlow[0][1][0]->SetMarkerColor(1);
  p_mFlow[0][1][0]->SetMarkerSize(1.4);
  p_mFlow[0][1][0]->Draw("pEX0 same");
  PlotLine(0.0,3.0,0.0,0.0,1,2,2);

  TLegend *leg = new TLegend(0.2,0.7,0.4,0.85);
  leg->SetFillColor(10);
  leg->SetBorderSize(0.0);
  leg->AddEntry(p_mFlow[4][1][0],"p","P");
  leg->AddEntry(p_mFlow[0][1][0],"#pi^{+}","P");
  leg->Draw("same");

  c_flow_2nd->SaveAs("./figures/c_flow_2nd.eps");
  c_flow_3rd->SaveAs("./figures/c_flow_3rd.eps");
}

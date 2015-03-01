#include "TFile.h"
#include "TString.h"
#include "TH1F.h"
#include "TCanvas.h"
#include <iostream>
#include "TLine.h"

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

void Centrality(Int_t mEnergy = 4, Int_t mMode = 0, Int_t mScreen = 0) // mEnergy = 0: 7GeV, 1: 11GeV, 2: 19GeV, 3: 27GeV, 4: 39GeV, 5: 62GeV, 6: 200GeV | mMode = 0: Default, 1: StringMelting | mScreen = 0: 1mb, 1: 3mb, 2: 6mb
{
  TString Energy[7] = {"7GeV","11GeV","19GeV","27GeV","39GeV","62GeV","200GeV"};
  TString Mode_AMPT[2] ={"Default","StringMelting"};
  TString ScreenMass_AMPT[3] ={"1mb","3mb","6mb"};
  Float_t Centrality_start[9] = {0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.05, 0.0};
  Float_t Centrality_stop[9]  = {0.8,0.7,0.6,0.5,0.4,0.3,0.2, 0.1,0.05};
//  Float_t Centrality_start[9] = {0.0,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7};
//  Float_t Centrality_stop[9]  = {0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8};
  TString inputfile;
  if(mMode == 0)
  {
    inputfile = Form("/home/xusun/Data/AMPT_%s/Resolution/%s_Resolution/Resolution_%s.root",Mode_AMPT[mMode].Data(),Energy[mEnergy].Data(),Energy[mEnergy].Data());
  }
  if(mMode == 1)
  {
    inputfile = Form("/home/xusun/Data/AMPT_%s/Resolution/%s_Resolution/%s/Resolution_%s.root",Mode_AMPT[mMode].Data(),Energy[mEnergy].Data(),ScreenMass_AMPT[mScreen].Data(),Energy[mEnergy].Data());
  }
  TFile *file_input = TFile::Open(inputfile.Data());
  cout << inputfile.Data() << endl;
  TH1F *h_refMult = (TH1F*)file_input->Get("h_mRefMult");
  Float_t Inte = h_refMult->Integral();
  Int_t binx_max = h_refMult->GetNbinsX();
  Int_t bin_cent[10];
  for(Int_t i_binx = binx_max; i_binx > 0; i_binx--)
  {
    Float_t Inte_cent = h_refMult->Integral(i_binx,binx_max);
    Float_t ratio = (Float_t)Inte_cent/(Float_t)Inte;
    if(ratio == 0.0)
    {
      bin_cent[9] = i_binx;
    }
    for(Int_t i_cent = 0; i_cent < 9; i_cent++)
    {
      if(ratio > Centrality_start[i_cent] && ratio <= Centrality_stop[i_cent])
      {
	bin_cent[i_cent] = i_binx;
      }
    }
  }
  for(Int_t i_cent = 0; i_cent < 10; i_cent++)
  {
    cout << "bin number = " << bin_cent[i_cent] << ", bin center = " << h_refMult->GetBinCenter(bin_cent[i_cent]) << endl;
  }
  Float_t y_max = h_refMult->GetMaximum();

  TCanvas *c_play = new TCanvas("c_play","c_play",10,10,800,800);
  c_play->cd();
  c_play->SetLeftMargin(0.15);
  c_play->SetBottomMargin(0.15);
  c_play->SetTicks(1,1);
  c_play->SetGrid(0,0);
  h_refMult->GetXaxis()->SetRangeUser(0,1000);
  h_refMult->SetTitle("");
  h_refMult->SetStats(0);
  h_refMult->GetXaxis()->SetTitle("refMult");
  h_refMult->GetYaxis()->SetTitle("Counts");
  h_refMult->GetXaxis()->CenterTitle();
  h_refMult->GetYaxis()->CenterTitle();
  h_refMult->Draw("hE");
  for(Int_t i_cent = 0; i_cent < 10; i_cent++)
  {
    PlotLine(h_refMult->GetBinCenter(bin_cent[i_cent]),h_refMult->GetBinCenter(bin_cent[i_cent]),0.0,y_max/2.0,1,2,2);
  }
}

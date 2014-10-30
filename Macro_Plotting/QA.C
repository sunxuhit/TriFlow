#include "TString.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TProfile.h"
#include "TGraphAsymmErrors.h"
#include "TMath.h"
#include "TStyle.h"

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

void Draw_TGAE_new_Symbol(TGraphAsymmErrors* tgae, Int_t style, Int_t color, Float_t size)
{
    TString HistName;
    Float_t size_A = 1.35*size;
    Float_t size_B = size;
    Float_t size_C = size;
    Int_t alt_marker = style;
    if(style == 24)
    {
        alt_marker = 20;
        size_A = 1.35*size;
    }
    if(style == 25)
    {
        alt_marker = 21;
        size_A = 1.35*size;
    }
    if(style == 26)
    {
        alt_marker = 22;
        size_A = 1.4*size;
    }
    if(style == 23)
    {
        alt_marker = 23;
        size_A = 1.35*size;
    }
    if(style == 30 || style == 29)
    {
        alt_marker = 29;
        size_A = 1.55*size;
    }

    // black and filled outer marker
    HistName = "tgae_dummy_A";
    TGraphAsymmErrors* ge_clone_A = (TGraphAsymmErrors*)tgae->Clone(HistName.Data());
    ge_clone_A->SetMarkerSize(size_A);
    ge_clone_A->SetMarkerStyle(alt_marker);
    ge_clone_A->SetMarkerColor(1);
    ge_clone_A->SetLineColor(10);
    ge_clone_A->Draw("same PZ0");

    // white and filled inner marker
    HistName = "tgae_dummy_B";
    TGraphAsymmErrors* ge_clone_B = (TGraphAsymmErrors*)tgae->Clone(HistName.Data());
    ge_clone_B->SetMarkerSize(size_B);
    ge_clone_B->SetMarkerStyle(alt_marker);
    ge_clone_B->SetMarkerColor(10);
    ge_clone_B->SetLineColor(10);
    ge_clone_B->Draw("same PZ0");

    // color inner marker
    HistName = "tgae_dummy_C";
    TGraphAsymmErrors* ge_clone_C = (TGraphAsymmErrors*)tgae->Clone(HistName.Data());
    ge_clone_C->SetMarkerSize(size_C);
    ge_clone_C->SetMarkerStyle(style);
    ge_clone_C->SetMarkerColor(color);
    ge_clone_C->SetLineColor(1);
    ge_clone_C->Draw("same PZ0");
}

void Draw_TGAE_Point_new_Symbol(Double_t x_val, Double_t y_val, Double_t x_min_err, Double_t x_max_err, Double_t y_min_err, Double_t y_max_err, Int_t style, Int_t color, Float_t size)
{
    TGraphAsymmErrors* tgae = new TGraphAsymmErrors();
    tgae->SetPoint(0,x_val,y_val);
    tgae->SetPointError(0,x_min_err,x_max_err,y_min_err,y_max_err);

    TString HistName;
    Float_t size_A = 1.35*size;
    Float_t size_B = size;
    Float_t size_C = size;
    Int_t alt_marker = style;
    Int_t style_in = style;
    if(style == 24)
    {
        alt_marker = 20;
        size_A = 1.35*size;
    }
    if(style == 25)
    {
        alt_marker = 21;
        size_A = 1.35*size;
    }
    if(style == 26)
    {
        alt_marker = 22;
        size_A = 1.4*size;
    }
    if(style == 23)
    {
        alt_marker = 23;
        size_A = 1.35*size;
    }
    if(style == 30 || style == 29)
    {
        alt_marker = 29;
        size_A = 1.55*size;
    }
    if(style == 260)
    {
        alt_marker = 26;
        size_A = 1.15*size;
        style = 26;
    }
    if(style == 300)
    {
        alt_marker = 30;
        size_A = 1.3*size;
        style = 30;
    }

    // black and filled outer marker
    HistName = "tgae_dummy_A";
    TGraphAsymmErrors* ge_clone_A = (TGraphAsymmErrors*)tgae->Clone(HistName.Data());
    ge_clone_A->SetMarkerSize(size_A);
    ge_clone_A->SetMarkerStyle(alt_marker);
    ge_clone_A->SetMarkerColor(1);
    ge_clone_A->SetLineColor(10);
    if(style_in == 260 || style_in == 300) ge_clone_A->SetLineColor(1);
    ge_clone_A->Draw("same PZ0");

    if(!(style_in == 260 || style_in == 300))
    {
        // white and filled inner marker
        HistName = "tgae_dummy_B";
        TGraphAsymmErrors* ge_clone_B = (TGraphAsymmErrors*)tgae->Clone(HistName.Data());
        ge_clone_B->SetMarkerSize(size_B);
        ge_clone_B->SetMarkerStyle(alt_marker);
        ge_clone_B->SetMarkerColor(10);
        ge_clone_B->SetLineColor(10);
        ge_clone_B->Draw("same PZ0");
    }

    // color inner marker
    HistName = "tgae_dummy_C";
    TGraphAsymmErrors* ge_clone_C = (TGraphAsymmErrors*)tgae->Clone(HistName.Data());
    ge_clone_C->SetMarkerSize(size_C);
    ge_clone_C->SetMarkerStyle(style);
    ge_clone_C->SetMarkerColor(color);
    ge_clone_C->SetLineColor(1);
    ge_clone_C->Draw("same PZ0");
}

static TString Mode[2] = {"Default","StringMelting"};
static TString Energy[7] = {"7GeV","11GeV","19GeV","27GeV","39GeV","62GeV","200GeV"};
static TString Order[2] = {"2nd","3rd"};
static TString Centrality[4] = {"0080","0010","1040","4080"};
static TString ParType[10] = {"pi_plus","pi_minus","K_plus","K_minus","p","pbar","Lambda","Lambdabar","K0s","phi"};

void QA(Int_t mEnergy = 4, Int_t mMode = 0) // 0: 7.7 GeV, 1: 11.5 GeV, 2: 19.6 GeV, 3: 27 GeV, 4: 39 GeV, 5: 62.4 GeV, 6: 200 GeV | 0: Default, 1: String Melting
{
  gStyle->SetTitleY(0.97);
  TString inputfile = Form("/home/xusun/Data/AMPT_%s/Resolution/%s_Resolution/Resolution_%s.root",Mode[mMode].Data(),Energy[mEnergy].Data(),Energy[mEnergy].Data()); // temperory file
  cout << "Input File: " << inputfile.Data() << endl;
  TFile *File_input = TFile::Open(inputfile.Data());

  TH1F *h_mRefMult = (TH1F*)File_input->Get("h_mRefMult");
  TH1F *h_mPsi2_East = (TH1F*)File_input->Get("h_mPsi2_East");
  TH1F *h_mPsi3_West = (TH1F*)File_input->Get("h_mPsi3_West");
  TH2F *h_mPsi2 = (TH2F*)File_input->Get("h_mPsi2");
  TH2F *h_mPsi3 = (TH2F*)File_input->Get("h_mPsi3");

  // Resolution
  TProfile *p_mRes2 = (TProfile*)File_input->Get("p_mRes2"); 
  TProfile *p_mRes3 = (TProfile*)File_input->Get("p_mRes3"); 
  Float_t Centrality_start[9] = {0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.05, 0.0};
  Float_t Centrality_stop[9]  = {0.8,0.7,0.6,0.5,0.4,0.3,0.2, 0.1,0.05};
  TGraphAsymmErrors *g_mRes2 = new TGraphAsymmErrors();
  TGraphAsymmErrors *g_mRes3 = new TGraphAsymmErrors();
  for(Int_t i_bin = 1; i_bin < 10; i_bin++)
  {
    if(p_mRes2->GetBinContent(i_bin) > 0.0)
    {
      Float_t res2 = TMath::Sqrt(p_mRes2->GetBinContent(i_bin));
      Float_t err_res2 = p_mRes2->GetBinError(i_bin)/(0.5*TMath::Sqrt(p_mRes2->GetBinContent(i_bin)));
      g_mRes2->SetPoint(i_bin-1,50.0*(Centrality_start[i_bin-1]+Centrality_stop[i_bin-1]),res2);
      g_mRes2->SetPointError(i_bin-1,0.0,0.0,err_res2,err_res2);
    }
    if(p_mRes3->GetBinContent(i_bin) > 0.0)
    {
      Float_t res3 = TMath::Sqrt(p_mRes3->GetBinContent(i_bin));
      Float_t err_res3 = p_mRes3->GetBinError(i_bin)/(0.5*TMath::Sqrt(p_mRes3->GetBinContent(i_bin)));
      g_mRes3->SetPoint(i_bin-1,50.0*(Centrality_start[i_bin-1]+Centrality_stop[i_bin-1]),res3);
      g_mRes3->SetPointError(i_bin-1,0.0,0.0,err_res3,err_res3);
    }
  }

  // refMult
  TCanvas *c_refmult = new TCanvas("c_refmult","c_refmult",10,10,800,800);
  c_refmult->cd();
  c_refmult->cd()->SetLeftMargin(0.15);
  c_refmult->cd()->SetBottomMargin(0.15);
  c_refmult->cd()->SetTicks(1,1);
  h_mRefMult->SetTitle("");
  h_mRefMult->SetStats(0);
  h_mRefMult->GetXaxis()->SetRangeUser(0.0,500.0);
  h_mRefMult->GetXaxis()->SetTitle("refMult");
  h_mRefMult->GetYaxis()->SetTitle("N_{events}");
  h_mRefMult->GetXaxis()->SetTitleSize(0.06);
  h_mRefMult->GetYaxis()->SetTitleSize(0.06);
  h_mRefMult->GetXaxis()->CenterTitle();
  h_mRefMult->GetYaxis()->CenterTitle();
  h_mRefMult->GetXaxis()->SetNdivisions(505,'X');
  h_mRefMult->GetYaxis()->SetNdivisions(505,'Y');
  h_mRefMult->Draw("PE");
  TString leg = Form("Au+Au, %s",Energy[mEnergy].Data());
  plotTopLegend(leg.Data(),0.43,0.7,0.04,1,0.0,42,1);
  TString leg = Form("AMPT_%s",Mode[mMode].Data());
  plotTopLegend(leg.Data(),0.43,0.63,0.04,1,0.0,42,1);

  // Psi2_East
  TCanvas *c_Psi2_East = new TCanvas("c_Psi2_East","c_Psi2_East",10,10,800,800);
  c_Psi2_East->cd();
  c_Psi2_East->cd()->SetLeftMargin(0.15);
  c_Psi2_East->cd()->SetBottomMargin(0.15);
  c_Psi2_East->cd()->SetTicks(1,1);
  h_mPsi2_East->SetTitle("");
  h_mPsi2_East->SetStats(0);
  h_mPsi2_East->GetXaxis()->SetRangeUser(-2.0,2.0);
  h_mPsi2_East->GetXaxis()->SetTitle("#Psi_{2,East}");
  h_mPsi2_East->GetYaxis()->SetTitle("N_{events}");
  h_mPsi2_East->GetXaxis()->SetTitleSize(0.06);
  h_mPsi2_East->GetYaxis()->SetTitleSize(0.06);
  h_mPsi2_East->GetXaxis()->CenterTitle();
  h_mPsi2_East->GetYaxis()->CenterTitle();
  h_mPsi2_East->GetXaxis()->SetNdivisions(505,'X');
  h_mPsi2_East->GetYaxis()->SetNdivisions(505,'Y');
  h_mPsi2_East->Draw("hE");
  TString leg = Form("Au+Au, %s, 0%%-80%%",Energy[mEnergy].Data());
  plotTopLegend(leg.Data(),0.35,0.3,0.04,1,0.0,42,1);
  TString leg = Form("AMPT_%s",Mode[mMode].Data());
  plotTopLegend(leg.Data(),0.41,0.25,0.04,1,0.0,42,1);

  // Psi3_West
  TCanvas *c_Psi3_West = new TCanvas("c_Psi3_West","c_Psi3_West",10,10,800,800);
  c_Psi3_West->cd();
  c_Psi3_West->cd()->SetLeftMargin(0.15);
  c_Psi3_West->cd()->SetBottomMargin(0.15);
  c_Psi3_West->cd()->SetTicks(1,1);
  h_mPsi3_West->SetTitle("");
  h_mPsi3_West->SetStats(0);
  h_mPsi3_West->GetXaxis()->SetRangeUser(-2.0,2.0);
//  h_mPsi3_West->GetYaxis()->SetRangeUser(30000,100000);
  h_mPsi3_West->GetYaxis()->SetRangeUser(h_mPsi3_West->GetMaximum()/2.0,h_mPsi3_West->GetMaximum()*1.5);
  h_mPsi3_West->GetXaxis()->SetTitle("#Psi_{3,West}");
  h_mPsi3_West->GetYaxis()->SetTitle("N_{events}");
  h_mPsi3_West->GetXaxis()->SetTitleSize(0.06);
  h_mPsi3_West->GetYaxis()->SetTitleSize(0.06);
  h_mPsi3_West->GetXaxis()->CenterTitle();
  h_mPsi3_West->GetYaxis()->CenterTitle();
  h_mPsi3_West->GetXaxis()->SetNdivisions(505,'X');
  h_mPsi3_West->GetYaxis()->SetNdivisions(505,'Y');
  h_mPsi3_West->Draw("hE");
  TString leg = Form("Au+Au, %s, 0%%-80%%",Energy[mEnergy].Data());
  plotTopLegend(leg.Data(),0.35,0.7,0.04,1,0.0,42,1);
  TString leg = Form("AMPT_%s",Mode[mMode].Data());
  plotTopLegend(leg.Data(),0.41,0.65,0.04,1,0.0,42,1);

  TString Title;
  // Psi2_East vs Psi2_West
  TCanvas *c_Psi2 = new TCanvas("c_Psi2","c_Psi2",10,10,800,800);
  c_Psi2->cd();
  c_Psi2->cd()->SetLeftMargin(0.15);
  c_Psi2->cd()->SetBottomMargin(0.15);
  c_Psi2->cd()->SetRightMargin(0.15);
  c_Psi2->cd()->SetTicks(1,1);
//  c_Psi2->cd()->SetLogz();
  Title = Form("Au+Au,%s,0%%-80%%,AMPT_%s",Energy[mEnergy].Data(),Mode[mMode].Data());
  h_mPsi2->SetTitle(Title.Data());
  h_mPsi2->SetStats(0);
  h_mPsi2->GetXaxis()->SetRangeUser(-2.0,2.0);
  h_mPsi2->GetXaxis()->SetTitle("#Psi_{2,East}");
  h_mPsi2->GetYaxis()->SetTitle("#Psi_{2,West}");
  h_mPsi2->GetXaxis()->SetTitleSize(0.06);
  h_mPsi2->GetYaxis()->SetTitleSize(0.06);
  h_mPsi2->GetXaxis()->CenterTitle();
  h_mPsi2->GetYaxis()->CenterTitle();
  h_mPsi2->GetXaxis()->SetNdivisions(505,'X');
  h_mPsi2->GetYaxis()->SetNdivisions(505,'Y');
  h_mPsi2->Draw("colz");

  // Psi3_East vs Psi3_West
  TCanvas *c_Psi3 = new TCanvas("c_Psi3","c_Psi3",10,10,800,800);
  c_Psi3->cd();
  c_Psi3->cd()->SetLeftMargin(0.15);
  c_Psi3->cd()->SetBottomMargin(0.15);
  c_Psi3->cd()->SetRightMargin(0.15);
  c_Psi3->cd()->SetTicks(1,1);
//  c_Psi3->cd()->SetLogz();
  Title = Form("Au+Au,%s,0%%-80%%,AMPT_%s",Energy[mEnergy].Data(),Mode[mMode].Data());
  h_mPsi3->SetTitle(Title.Data());
  h_mPsi3->SetStats(0);
  h_mPsi3->GetXaxis()->SetRangeUser(-2.0,2.0);
  h_mPsi3->GetXaxis()->SetTitle("#Psi_{3,East}");
  h_mPsi3->GetYaxis()->SetTitle("#Psi_{3,West}");
  h_mPsi3->GetXaxis()->SetTitleSize(0.06);
  h_mPsi3->GetYaxis()->SetTitleSize(0.06);
  h_mPsi3->GetXaxis()->CenterTitle();
  h_mPsi3->GetYaxis()->CenterTitle();
  h_mPsi3->GetXaxis()->SetNdivisions(505,'X');
  h_mPsi3->GetYaxis()->SetNdivisions(505,'Y');
  h_mPsi3->Draw("colz");

  // Resolution
  TCanvas *c_Res = new TCanvas("c_Res","c_Res",10,10,800,800);
  c_Res->cd();
  c_Res->cd()->SetLeftMargin(0.15);
  c_Res->cd()->SetBottomMargin(0.15);
  c_Res->cd()->SetRightMargin(0.15);
  c_Res->cd()->SetTicks(1,1);
  TH1F *h_play = new TH1F("h_play","h_play",100,0,100);
  for(Int_t i_bin = 0; i_bin < 100; i_bin++)
  {
    h_play->SetBinContent(i_bin,-10.0);
  }
  Title = Form("Au+Au,%s,0%%-80%%,AMPT_%s",Energy[mEnergy].Data(),Mode[mMode].Data());
  h_play->SetTitle(Title.Data());
  h_play->SetStats(0);
  h_play->GetXaxis()->SetRangeUser(0,80);
  h_play->GetYaxis()->SetRangeUser(-0.05,0.7);
  h_play->GetXaxis()->SetTitle("Centrality(%)");
  h_play->GetYaxis()->SetTitle("Resolution");
  h_play->GetXaxis()->CenterTitle();
  h_play->GetYaxis()->CenterTitle();
  h_play->GetXaxis()->SetTitleSize(0.06);
  h_play->GetYaxis()->SetTitleSize(0.06);
  h_play->GetXaxis()->SetNdivisions(505,'X');
  h_play->GetYaxis()->SetNdivisions(505,'Y');
  h_play->Draw("pE");
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_mRes2,24,2,1.4);
  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_mRes3,24,kGray,1.4);
  PlotLine(0.0,80.0,0.0,0.0,1,2,2);
  Draw_TGAE_Point_new_Symbol(45,0.61,0.0,0.0,0,0.0,24,2,1.4);
  plotTopLegend("2^{nd} Event Plane",47,0.598,0.04,1,0.0,42,0);
  Draw_TGAE_Point_new_Symbol(45,0.56,0.0,0.0,0,0.0,24,kGray,1.4);
  plotTopLegend("3^{rd} Event Plane",47,0.548,0.04,1,0.0,42,0);

  // Save figures
  c_refmult->SaveAs("./figures/c_refmult.eps");
  c_Psi2_East->SaveAs("./figures/c_Psi2_East.eps");
  c_Psi3_West->SaveAs("./figures/c_Psi3_West.eps");
  c_Psi2->SaveAs("./figures/c_Psi2.eps");
  c_Psi3->SaveAs("./figures/c_Psi3.eps");
  c_Res->SaveAs("./figures/c_Res.eps");
}

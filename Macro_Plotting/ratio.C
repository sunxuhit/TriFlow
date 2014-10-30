#include "TString.h"
#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
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
static TString ParName[5] = {"#pi","K","p","#phi","#Lambda"};

// Calculate integrated v2 and v3
void ratio(Int_t mEnergy = 4, Int_t mMode = 0) // 0: 7.7 GeV, 1: 11.5 GeV, 2: 19.6 GeV, 3: 27 GeV, 4: 39 GeV, 5: 62.4 GeV, 6: 200 GeV | 0: Default, 1: String Melting
{
  gStyle->SetTitleX(0.55);
  gStyle->SetTitleY(0.98);
  TString inputfile = Form("/home/xusun/Data/AMPT_%s/InteFlow/%s_%s/InteFlow.root",Mode[mMode].Data(),Energy[mEnergy].Data(),Mode[mMode].Data());
  TFile *File_input = TFile::Open(inputfile.Data());

  // pi_plus,pi_minus,K_plus,K_minus,p,pbar,Lambda,Lambdabar,K0s,phi
  TH1F *h_ratio = (TH1F*)File_input->Get("h_ratio_0080");
  TGraphAsymmErrors *g_ratio = new TGraphAsymmErrors();
  g_ratio->SetPoint(0,0.5,h_ratio->GetBinContent(1)); // pi_plus
  g_ratio->SetPointError(0,0.0,0.0,h_ratio->GetBinError(1),h_ratio->GetBinError(1)); // pi_plus
  g_ratio->SetPoint(1,1.5,h_ratio->GetBinContent(3)); // K_plus
  g_ratio->SetPointError(1,0.0,0.0,h_ratio->GetBinError(3),h_ratio->GetBinError(5)); // K_plus
  g_ratio->SetPoint(2,2.5,h_ratio->GetBinContent(5)); // p
  g_ratio->SetPointError(2,0.0,0.0,h_ratio->GetBinError(5),h_ratio->GetBinError(5)); // p
  g_ratio->SetPoint(3,3.5,h_ratio->GetBinContent(10)); // phi
  g_ratio->SetPointError(3,0.0,0.0,h_ratio->GetBinError(10),h_ratio->GetBinError(10)); // phi 
  g_ratio->SetPoint(4,4.5,h_ratio->GetBinContent(7)); // Lambda
  g_ratio->SetPointError(4,0.0,0.0,h_ratio->GetBinError(7),h_ratio->GetBinError(7)); // Lambda

  TCanvas *c_play = new TCanvas("c_play","c_play",10,10,800,800);
  c_play->SetLeftMargin(0.2);
  c_play->SetBottomMargin(0.2);
  c_play->SetTicks(0,1);
  TH1F *h_play = new TH1F("h_play","h_play",100,0.0,10.0);
  for(Int_t i_bin = 1; i_bin < 101; i_bin++)
  {
    h_play->SetBinContent(i_bin,-10.0);
  }
  TString Title = Form("Au+Au, %s, 0%%-70%%, AMPT_%s",Energy[mEnergy].Data(),Mode[mMode].Data());
  h_play->SetTitle(Title.Data());
  h_play->SetStats(0);
  h_play->GetXaxis()->SetRangeUser(0.0,5.0);
  h_play->GetYaxis()->SetRangeUser(2.0,4.5);
  h_play->SetNdivisions(400,"X");
  h_play->SetNdivisions(505,"Y");
  h_play->GetXaxis()->SetLabelSize(0.0);
  h_play->GetYaxis()->SetTitle("v_{2}/v_{3}");
  h_play->GetYaxis()->SetTitleSize(0.06);
  h_play->GetYaxis()->CenterTitle();
  h_play->GetXaxis()->SetTickLength(0.0);
  h_play->Draw("pE");
  Draw_TGAE_new_Symbol(g_ratio,24,kGray,1.4);

  for(Int_t i_par = 0; i_par < 5; i_par++)
  {
    plotTopLegend((char*)ParName[i_par].Data(),0.26+0.14*i_par,0.117,0.06,1,0.0,42,1);
  }
  PlotLine(0.0,5.0,2.35,2.35,2,2,2);

  Draw_TGAE_Point_new_Symbol(0.5,4.3,0.0,0.0,0.0,0.0,24,kGray,1.4);
  plotTopLegend("AMPT",0.7,4.25,0.04,1,0.0,42,0);
  PlotLine(0.42,0.59,4.1,4.1,2,2,2);
  plotTopLegend("Data",0.7,4.07,0.04,1,0.0,42,0);

  c_play->SaveAs("./figures/c_ratio_AMPT.eps");
}

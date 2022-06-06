#include "TLine.h"
#include "TPolyLine.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TLatex.h"

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

void PlotHistLine2(TH1F* Histo, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle, Float_t x_start, Float_t x_stop, Float_t y_min, Float_t y_max)
{
    TPolyLine* Hist_line = new TPolyLine();
    Hist_line -> SetLineWidth(LineWidth);
    Hist_line -> SetLineStyle(LineStyle);
    Hist_line -> SetLineColor(Line_Col);
    for(Int_t i = 1; i < Histo->GetNbinsX(); i++)
    {
        Double_t x = Histo->GetBinCenter(i);
        Double_t y = Histo->GetBinContent(i);
        if(y >= y_min && y < y_max && y != 0 && x >= x_start && x <= x_stop)
        {
            //cout << "x = " << x << ", y = " << y << endl;
            Hist_line->SetNextPoint(x,y);
        }
    }
    Hist_line -> Draw();
    //delete Hist_line;
}


void PlotHistErrorBand(TH1F* Histo, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle, Int_t FillStyle, Int_t FillColor)
{
    Int_t N_points = 0;
    for(Int_t i = 1; i < Histo->GetNbinsX(); i++)
    {
        Double_t y = Histo->GetBinContent(i);
        if(y != 0)
        {
            N_points++;
        }
    }
    Int_t N_total_points = N_points*2+1;
    //cout << "N_total_points = " << N_total_points << endl;
    TGraph* Hist_line = new TGraph(N_total_points);
    Hist_line -> SetLineWidth(LineWidth);
    Hist_line -> SetLineStyle(LineStyle);
    Hist_line -> SetLineColor(Line_Col);
    Hist_line -> SetFillStyle(FillStyle);
    Hist_line -> SetFillColor(FillColor);
    Int_t N_point = 0;
    for(Int_t i = 1; i < Histo->GetNbinsX(); i++)
    {
        Double_t y     = Histo->GetBinContent(i);
        Double_t y_err = Histo->GetBinError(i);
        if(y != 0)
        {
            Double_t x = Histo->GetBinCenter(i);
            //cout << "N_point = " << N_point << ", x = " << x << ", y = " << y << ", y_err = " << y_err << endl;
            Hist_line->SetPoint(N_point,x,y-y_err);
            Hist_line->SetPoint(N_total_points-2-N_point,x,y+y_err);
            if(N_point == 0) Hist_line->SetPoint(N_total_points-1,x,y-y_err);
            N_point++;
        }
    }
    Hist_line -> Draw("f");
    //delete Hist_line;
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

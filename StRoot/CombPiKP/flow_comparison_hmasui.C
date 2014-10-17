#include "TFile.h"
#include "TString.h"
#include "TGraphAsymmErrors.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "draw.h"
#include "TMath.h"
#include "TF1.h"
#include <iostream>

Double_t ErrorAdd(Float_t x, Float_t y)
{
  return sqrt(x*x+y*y);
}

Double_t ErrTimes(Float_t x, Float_t y, Float_t dx, Float_t dy)
{
  return x*y*ErrorAdd(dx/x,dy/y);
}

Double_t ErrDiv(Float_t x, Float_t y, Float_t dx, Float_t dy)
{
  return x/y*ErrorAdd(dx/x,dy/y);
}

Double_t poly(Double_t *x_val, Double_t *par)
{
  Double_t x = x_val[0];
  Double_t y;

  y = par[0] + par[1]*x + par[2]*x*x + par[3]*x*x*x + par[4]*x*x*x*x + par[5]*x*x*x*x*x;

  return y;
}

// mEnergy: 0 for 200GeV, 1 for 39GeV
void flow_comparison_hmasui(Int_t mEnergy = 0)
{
  const Int_t charge_total = 2;
  const Int_t eta_total = 1;
  const Int_t flow_total = 2;
  const Int_t N_Species = 3;
  TString Order[2] = {"2nd","3rd"};
  TString Name[charge_total][N_Species] = {
                                            {"#pi^{+}","K^{+}","p"},
					    {"#pi^{-}","K^{-}","#bar{p}"}
                                          };
  TString TitleY[2] = {"v_{2}","v_{3}"};
//  TString EtaGap[eta_total] = {"#eta_{gap} = #pm 0.05", "#eta_{gap} = #pm 0.10"};
  TString EtaGap[eta_total] = {"#eta_{gap} = #pm 0.05"};
  Float_t rangeY[2] = {0.2,0.14};
  Double_t mass[N_Species] = {0.139,0.494,0.938};
  Double_t nq[N_Species] = {2.0,2.0,3.0};

  const Int_t Draw_Color[N_Species] = {2,kGray+2,kAzure-2};
  const Int_t Draw_Style[N_Species] = {29,21,20};
  TString Energy[2] = {"200GeV","39GeV"};

  TFile *file_pik_pos, *file_pik_neg, *file_proton;
  file_pik_pos = TFile::Open(Form("/project/projectdirs/star/xusun/OutPut/AuAu%s/Mass2_nSigmaPion/merged_file/flow_pik/Flow_pik_0080_charge_0_etagap_0.root",Energy[mEnergy].Data()));
  file_pik_neg = TFile::Open(Form("/project/projectdirs/star/xusun/OutPut/AuAu%s/Mass2_nSigmaPion/merged_file/flow_pik/Flow_pik_0080_charge_1_etagap_0.root",Energy[mEnergy].Data()));
  file_proton = TFile::Open(Form("/project/projectdirs/star/xusun/OutPut/AuAu%s/Mass2_Proton/merged_file/flow_proton/Flow_proton_0080_etagap_0.root",Energy[mEnergy].Data()));

  TFile *file_hmasui = TFile::Open("/project/projectdirs/star/xusun/OutPut/AuAu200GeV/Mass2_nSigmaPion/merged_file/flow_pik/v2pt_pikp_K0sLambda_preliminary_Results_Jul25_2012.root");
  TGraphAsymmErrors *g_flow_hmasui[N_Species+1];
  g_flow_hmasui[0] = (TGraphAsymmErrors*)file_hmasui->Get("gv2pt_5_0"); // pion
  g_flow_hmasui[1] = (TGraphAsymmErrors*)file_hmasui->Get("gv2pt_5_1"); // kaon
  g_flow_hmasui[2] = (TGraphAsymmErrors*)file_hmasui->Get("gv2pt_5_2"); // proton
  g_flow_hmasui[3] = (TGraphAsymmErrors*)file_hmasui->Get("gv2pt_5_3"); // k0s

  TF1 *f_flow_hmasui[N_Species];
  for(Int_t i = 0; i < N_Species; i++)
  {
    TString f_Name = Form("f_%d",i);
    f_flow_hmasui[i] = new TF1(f_Name.Data(),poly,0.4,3.0,6);
    f_flow_hmasui[i]->SetParameter(0,0.0);
    f_flow_hmasui[i]->SetParameter(1,0.00);
    f_flow_hmasui[i]->SetParameter(2,0.00);
    f_flow_hmasui[i]->SetParameter(3,0.00);
    f_flow_hmasui[i]->SetParameter(4,0.00);
    f_flow_hmasui[i]->SetParameter(5,0.00);
    f_flow_hmasui[i]->SetRange(0.3,2.8);
    g_flow_hmasui[i]->Fit(f_flow_hmasui[i],"R");
  }

  TGraphAsymmErrors *g_flow[charge_total][eta_total][flow_total][N_Species]; // charge: 0 = pos, 1 = neg | eta: 0 = +/- 0.05, 1 = +/- 0.1 | flow: 0 = 2nd, 1 = 3rd | species: 0 = pion, 1 = kaon, 2 = proton

  for(Int_t charge = 0; charge < charge_total; charge++)
  {
    for(Int_t eta_bin = 0; eta_bin < eta_total; eta_bin++)
    {
      for(Int_t flow_bin = 0; flow_bin < flow_total; flow_bin++)
      {
	for(Int_t n_species = 0; n_species < N_Species; n_species++)
	{
	  TString g_Name_pion = Form("g_Centrality_0_charge_%d_etagap_%d_%s_flow_pion",charge,eta_bin,Order[flow_bin].Data());
	  TString g_Name_kaon = Form("g_Centrality_0_charge_%d_etagap_%d_%s_flow_kaon",charge,eta_bin,Order[flow_bin].Data());
	  if(charge == 0)
	  {
	    g_flow[charge][eta_bin][flow_bin][0] = (TGraphAsymmErrors*)file_pik_pos->Get(g_Name_pion.Data());
	    g_flow[charge][eta_bin][flow_bin][1] = (TGraphAsymmErrors*)file_pik_pos->Get(g_Name_kaon.Data());
	  }
	  if(charge == 1)
	  {
	    g_flow[charge][eta_bin][flow_bin][0] = (TGraphAsymmErrors*)file_pik_neg->Get(g_Name_pion.Data());
	    g_flow[charge][eta_bin][flow_bin][1] = (TGraphAsymmErrors*)file_pik_neg->Get(g_Name_kaon.Data());
	  }
	  TString g_Name_proton = Form("g_Centrality_0_charge_%d_etagap_%d_%s_flow_proton",charge,eta_bin,Order[flow_bin].Data());
	  g_flow[charge][eta_bin][flow_bin][2] = (TGraphAsymmErrors*)file_proton->Get(g_Name_proton.Data());
	}
      }
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------
  Double_t label_size_factor = 1.2;
  // set canvas for v2

  TCanvas* c_v2_Ratio_v2_6 = new TCanvas("c_v2_Ratio_v2_6","c_v2_Ratio_v2_6",1400,10,1000,800);
  c_v2_Ratio_v2_6->SetFillColor(10);
  c_v2_Ratio_v2_6->SetTopMargin(0.0);
  c_v2_Ratio_v2_6->SetBottomMargin(0.3);
  c_v2_Ratio_v2_6->SetRightMargin(0.0);
  c_v2_Ratio_v2_6->SetLeftMargin(0.2);
  c_v2_Ratio_v2_6->Divide(3,4,0.0,0.0,10);
  TH1F* h_test_6p_2[12];

  for(Int_t j = 0; j < 12; j++)
  {
    TString HistName = "h_test_6p_2_";
    HistName += j;
    h_test_6p_2[j] = new TH1F(HistName.Data(),HistName.Data(),1000,-0.5,6);
    for(Int_t i = 0; i < h_test_6p_2[j]->GetNbinsX(); i++)
    {
      h_test_6p_2[j]->SetBinContent(i,-10);
    }
  }

  Int_t x_loop  = 0;
  Int_t y_loop2 = 0;
  Int_t y_loop[12] = {0,0,0,1,1,1,2,2,2,3,3,3};
  for(Int_t i = 0; i < 12; i++)
  {
    Float_t pad_factor = 1.0;
    if(i >= 6 && i <= 8) pad_factor = 1.1;
    if(i <= 2) {y_loop2 = 0;}
    else y_loop2 = 1;
    c_v2_Ratio_v2_6->cd(i+1);
    c_v2_Ratio_v2_6->cd(i+1);
    c_v2_Ratio_v2_6->cd(i+1)->SetTicks(1,1);
    c_v2_Ratio_v2_6->cd(i+1)->SetGrid(0,0);
    c_v2_Ratio_v2_6->cd(i+1)->SetFillColor(10);
    if(i >= 3 && i <= 5) c_v2_Ratio_v2_6->cd(i+1)->SetBottomMargin(0.01);
    if(i >= 9 && i <= 11) c_v2_Ratio_v2_6->cd(i+1)->SetBottomMargin(0.45);

    if(i < 3) c_v2_Ratio_v2_6->cd(i+1)->SetTopMargin(0.1);
    if(i == 2 || i == 5 || i == 8 || i == 11) c_v2_Ratio_v2_6->cd(i+1)->SetRightMargin(0.1);
    if(i == 3 || i == 6 || i == 9) x_loop = 0;
    Float_t x1_array[3] = {0,0.387,0.678};
    Float_t x2_array[3] = {x1_array[1],x1_array[2],1.0};
    Float_t y1_array[4] = {0.64,0.526+0.005,0.2,0.0};
    Float_t y2_array[4] = {1.0,y1_array[0],y1_array[1]-0.005,y1_array[2]};
    c_v2_Ratio_v2_6->cd(i+1)->SetPad(x1_array[x_loop],y1_array[y_loop[i]],x2_array[x_loop],y2_array[y_loop[i]]); // x1, y1, x2, y2
//    cout << "i = " << i << ", x1 = " << x1_array[x_loop] << ", x2 = " << x2_array[x_loop] << ", y1 = " << y1_array[y_loop[i]]
//      << ", y2 = " << y2_array[y_loop[i]] << endl;

    h_test_6p_2[i]->SetStats(0);
    h_test_6p_2[i]->SetTitle("");
    h_test_6p_2[i]->GetXaxis()->SetTitleOffset(1.2);
    h_test_6p_2[i]->GetYaxis()->SetTitleOffset(1.26);
    h_test_6p_2[i]->GetYaxis()->SetLabelOffset(0.01);
    h_test_6p_2[i]->GetXaxis()->SetLabelSize(0.09*label_size_factor);
    h_test_6p_2[i]->GetYaxis()->SetLabelSize(0.09*label_size_factor);
    h_test_6p_2[i]->GetXaxis()->SetTitleSize(0.09*label_size_factor);
    h_test_6p_2[i]->GetYaxis()->SetTitleSize(0.09*label_size_factor);
    h_test_6p_2[i]->GetXaxis()->SetNdivisions(505,'N');
    h_test_6p_2[i]->GetYaxis()->SetNdivisions(505,'N');
    h_test_6p_2[i]->GetXaxis()->CenterTitle();
    h_test_6p_2[i]->GetYaxis()->CenterTitle();
    h_test_6p_2[i]->GetXaxis()->SetRangeUser(-0.2,3.1);
    h_test_6p_2[i]->GetYaxis()->SetRangeUser(-0.02,0.22);
    h_test_6p_2[i]->GetYaxis()->SetTitle("");
    h_test_6p_2[i]->GetXaxis()->SetTitle("");

    if(i == 3 || i == 4 || i == 5) {h_test_6p_2[i]->SetTickLength(0.04,"X");}
    if(i == 2 || i == 8) {h_test_6p_2[i]->SetTickLength(0.04,"Y");}
    if(i == 1 || i == 7) {h_test_6p_2[i]->SetTickLength(0.042,"Y");}
    if(i == 4)  {h_test_6p_2[i]->SetTickLength(0.04,"Y");}
    if(i == 5)  {h_test_6p_2[i]->SetTickLength(0.04,"Y");}
    if(i == 9)  {h_test_6p_2[i]->SetTickLength(0.052,"Y");}
    if(i == 10) {h_test_6p_2[i]->SetTickLength(0.075,"Y");}
    if(i == 11) {h_test_6p_2[i]->SetTickLength(0.07,"Y");}

    if(i == 3 || i == 4 || i == 5 || i == 9 || i == 10 || i == 11)
    {
      h_test_6p_2[i]->GetYaxis()->SetRangeUser(0.88,1.07);
      h_test_6p_2[i]->GetYaxis()->SetNdivisions(503,'N');
    }
    if(i >= 9)
    {
      h_test_6p_2[i]->GetXaxis()->SetTitleOffset(1.05);
      h_test_6p_2[i]->GetXaxis()->SetLabelOffset(0.008);
      h_test_6p_2[i]->GetXaxis()->SetLabelSize(0.16*label_size_factor);
      h_test_6p_2[i]->GetXaxis()->SetTitleSize(0.16*label_size_factor);
    }
    if(i == 3)
    {
      h_test_6p_2[i]->GetYaxis()->SetTitleOffset(0.43);
      h_test_6p_2[i]->GetYaxis()->SetLabelOffset(0.01);
      h_test_6p_2[i]->GetYaxis()->SetLabelSize(0.27*label_size_factor);
      h_test_6p_2[i]->GetYaxis()->SetTitleSize(0.27*label_size_factor);
    }
    if(i == 6)
    {
      h_test_6p_2[i]->GetYaxis()->SetTitleOffset(1.1);
      h_test_6p_2[i]->GetYaxis()->SetLabelOffset(0.01);
      h_test_6p_2[i]->GetYaxis()->SetLabelSize(0.1*label_size_factor);
      h_test_6p_2[i]->GetYaxis()->SetTitleSize(0.1*label_size_factor);
    }
    if(i == 9)
    {
      h_test_6p_2[i]->GetYaxis()->SetTitleOffset(0.72);
      h_test_6p_2[i]->GetYaxis()->SetLabelOffset(0.01);
      h_test_6p_2[i]->GetYaxis()->SetLabelSize(0.16*label_size_factor);
      h_test_6p_2[i]->GetYaxis()->SetTitleSize(0.16*label_size_factor);
    }
    if(i == 0 || i == 6) h_test_6p_2[i]->GetYaxis()->SetTitle("v_{2}");
    if(i == 3 || i == 9) h_test_6p_2[i]->GetYaxis()->SetTitle("Ratio");
    if(i == 10) h_test_6p_2[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");

    h_test_6p_2[i]->DrawCopy("h");
    if(i == 0 || i == 1 || i == 2 || i == 6 || i == 7 || i == 8)
    {
      Float_t x_pos[3] = {0.3,0.07,0.068};
      Float_t y_pos[2] = {0.77,0.86};
//      plotTopLegend((char*)"blubb",x_pos[x_loop],y_pos[y_loop2],0.09*pad_factor,1,0.0,42,1);
    }
    PlotLine(-0.2,3.1,1.0,1.0,1,1,2); // x1_val, x2_val, y1_val, y2_val, Line_Col, LineWidth, LineStyle
    x_loop++;
  }

  //------------------------------------------------------------------------------------------------------------------------------------

  TGraphAsymmErrors *g_ratio_Xu[charge_total][N_Species];
  TGraphAsymmErrors *g_ratio_hmasui[N_Species];
  for(Int_t charge = 0; charge < charge_total; charge++)
  {
    for(Int_t n_species = 0; n_species < N_Species; n_species++)
    {
      g_ratio_Xu[charge][n_species] = new TGraphAsymmErrors();
      if(charge == 0)
      {
	g_ratio_hmasui[n_species] = new TGraphAsymmErrors();
	if(n_species != 2)
	{
	  for(Int_t i = 2; i < 18; i++)
	  {
	    Double_t x_hmasui, y_hmasui, err_hmasui;
	    g_flow_hmasui[n_species]->GetPoint(i,x_hmasui,y_hmasui);
	    err_hmasui = g_flow_hmasui[n_species]->GetErrorYhigh(i);

	    Double_t y_line, err_line;
	    y_line = f_flow_hmasui[n_species]->Eval(x_hmasui);
	    err_line = 0.0;

	    Double_t ratio_HL = y_hmasui/y_line;
	    Double_t err_HL = ErrDiv(y_hmasui,y_line,err_hmasui,err_line);

	    //	  if(n_species == 2) cout << "i = " << i << ", x = " << x_hmasui << endl;
	    //	  if(n_species == 2) cout << "i = " << i << ", ratio = " << ratio_HL << endl;

	    g_ratio_hmasui[n_species]->SetPoint(i-2,x_hmasui,ratio_HL);
	    g_ratio_hmasui[n_species]->SetPointError(i-2,0.0,0.0,err_HL,err_HL);
	  }
	}
	if(n_species == 2)
	{
	  for(Int_t i = 2; i < 23; i++)
	  {
	    if(i < 11 || i > 15)
	    {
	      Double_t x_hmasui, y_hmasui, err_hmasui;
	      g_flow_hmasui[n_species]->GetPoint(i,x_hmasui,y_hmasui);
	      err_hmasui = g_flow_hmasui[n_species]->GetErrorYhigh(i);

	      Double_t y_line, err_line;
	      y_line = f_flow_hmasui[n_species]->Eval(x_hmasui);
	      err_line = 0.0;

	      Double_t ratio_HL = y_hmasui/y_line;
	      Double_t err_HL = ErrDiv(y_hmasui,y_line,err_hmasui,err_line);

	      //	  if(n_species == 2) cout << "i = " << i << ", x = " << x_hmasui << endl;
	      //	  if(n_species == 2) cout << "i = " << i << ", ratio = " << ratio_HL << endl;

	      g_ratio_hmasui[n_species]->SetPoint(i-2,x_hmasui,ratio_HL);
	      g_ratio_hmasui[n_species]->SetPointError(i-2,0.0,0.0,err_HL,err_HL);
	    }
	  }
	}
      }

      for(Int_t i = 1; i < 13; i++)
      {
	Double_t x_Xu, y_Xu, err_Xu;
	g_flow[charge][0][0][n_species]->GetPoint(i,x_Xu,y_Xu);
	err_Xu = g_flow[charge][0][0][n_species]->GetErrorYhigh(i);

	Double_t y_line, err_line;
	y_line = f_flow_hmasui[n_species]->Eval(x_Xu);
	err_line = 0.0;

	Double_t ratio_XL = y_Xu/y_line;
	Double_t err_XL = ErrDiv(y_Xu,y_line,err_Xu,err_line);

//	cout << "pt = " << x_hmasui << ", ratio_XA = " << ratio_XA << ", err_XA = " << err_XA << endl;

	g_ratio_Xu[charge][n_species]->SetPoint(i-1,x_Xu,ratio_XL);
//	g_ratio[charge][n_species]->SetPointError(i,0.0,0.0,0.0,0.0);
	g_ratio_Xu[charge][n_species]->SetPointError(i-1,0.0,0.0,err_XL,err_XL);
      }
    }
  }

  //------------------------------------------------------------------------------------------------------------------------------------

  for(Int_t i = 0; i < 6; i++)
  {
    if(i < 3)
    {
      c_v2_Ratio_v2_6->cd(i+1);

      Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_flow[0][0][0][i],24,2,0.8);
      Draw_TGAE_Point_new_Symbol(0.1,0.18,0.0,0.0,0.0,0.0,24,2,1.0);
      plotTopLegend((char*)Name[0][i].Data(),0.2,0.177,0.06,1,0.0,42,0);

      Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_flow_hmasui[i],21,kGray+1,0.8);
      Draw_TGAE_Point_new_Symbol(0.1,0.16,0.0,0.0,0.0,0.0,21,kGray+1,1.0);
      plotTopLegend((char*)(Name[0][i]+"+"+Name[1][i]+" (Hiroshi)").Data(),0.2,0.157,0.06,1,0.0,42,0);

      if(i == 0) plotTopLegend((char*)"Au+Au, 200 GeV, 0-80%",0.6,0.177,0.06,1,0.0,42,0);

      if(i == 1) plotTopLegend((char*)EtaGap[0].Data(),0.6,0.18,0.06,1,0.0,42,0);

      c_v2_Ratio_v2_6->cd(i+4);
      Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_ratio_Xu[0][i],24,2,0.8);
      Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_ratio_hmasui[i],21,kGray+1,0.8);
    }
    else
    {
      c_v2_Ratio_v2_6->cd(i+4);

      Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_flow[1][0][0][i-3],24,1,0.8);
      Draw_TGAE_Point_new_Symbol(0.1,0.18,0.0,0.0,0.0,0.0,24,1,1.0);
      plotTopLegend((char*)Name[1][i-3].Data(),0.2,0.177,0.06,1,0.0,42,0);

      Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_flow_hmasui[i-3],21,kGray+1,0.8);
      Draw_TGAE_Point_new_Symbol(0.1,0.16,0.0,0.0,0.0,0.0,21,kGray+1,1.0);
      plotTopLegend((char*)(Name[0][i-3]+"+"+Name[1][i-3]+" (Hiroshi)").Data(),0.2,0.157,0.06,1,0.0,42,0);

      c_v2_Ratio_v2_6->cd(i+7);
      Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_ratio_Xu[1][i-3],24,1,0.8);
      Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_ratio_hmasui[i-3],21,kGray+1,0.8);
    }
  }
  c_v2_Ratio_v2_6->SaveAs("./figures/flow/c_v2_Ratio_v2_6.eps");
}

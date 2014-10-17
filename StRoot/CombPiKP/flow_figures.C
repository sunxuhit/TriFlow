#include "TFile.h"
#include "TString.h"
#include "TGraphAsymmErrors.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "draw.h"
#include "TMath.h"
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
// mEnergy: 0 for 200GeV, 1 for 39GeV
void flow_figures(Int_t mEnergy = 0)
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
  // set frame
  TH1F* pos_neg_dummy = new TH1F("pos_neg_dummy","pos_neg_dummy",500,-0.1,3.4);
  for(Int_t i = 0; i < pos_neg_dummy->GetNbinsX(); i++)
  {
    pos_neg_dummy->SetBinContent(i,-10);
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
  pos_neg_dummy->GetXaxis()->SetTitle("p_{T} GeV/c");

  //------------------------------------------------------------------------------------------------------------------------------------

  //------------------------------------------------------------------------------------------------------------------------------------
  Double_t label_size_factor = 1.2;
  /*
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
  for(Int_t i = 0; i < 15; i++)
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
      h_test_6p_2[i]->GetYaxis()->SetRangeUser(0.93,1.07);
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

  // v2 comparision
  TGraphAsymmErrors *g_Alex[charge_total][N_Species];

  // pion v2 from Alex
  TString input_pion_alex_0 = "/Users/mac/STAR/Data/BES/v2/Merge_PiP_v2_39GeV_histo_C_mult_0_ep_1_syst.root";
  TFile *file_pion_alex_0 = TFile::Open(input_pion_alex_0.Data());
  g_Alex[0][0] = (TGraphAsymmErrors*)file_pion_alex_0->Get("tgae_PiP_v2_39GeV_mult_0_ep_1_err_0");

  TString input_pion_alex_1 = "/Users/mac/STAR/Data/BES/v2/Merge_PiM_v2_39GeV_histo_C_mult_0_ep_1_syst.root";
  TFile *file_pion_alex_1 = TFile::Open(input_pion_alex_1.Data());
  g_Alex[1][0] = (TGraphAsymmErrors*)file_pion_alex_1->Get("tgae_PiM_v2_39GeV_mult_0_ep_1_err_0");

  // kaon v2 from Alex
  TString input_kaon_alex_0 = "/Users/mac/STAR/Data/BES/v2/Merge_KP_v2_39GeV_histo_C_mult_0_ep_1_syst.root";
  TFile *file_kaon_alex_0 = TFile::Open(input_kaon_alex_0.Data());
  g_Alex[0][1] = (TGraphAsymmErrors*)file_kaon_alex_0->Get("tgae_KP_v2_39GeV_mult_0_ep_1_err_0");

  TString input_kaon_alex_1 = "/Users/mac/STAR/Data/BES/v2/Merge_KM_v2_39GeV_histo_C_mult_0_ep_1_syst.root";
  TFile *file_kaon_alex_1 = TFile::Open(input_kaon_alex_1.Data());
  g_Alex[1][1] = (TGraphAsymmErrors*)file_kaon_alex_1->Get("tgae_KM_v2_39GeV_mult_0_ep_1_err_0");

  // proton v2 from Alex
  TString flow_proton_Alex_0 = "/Users/mac/STAR/Data/BES/v2/Merge_Proton_v2_39GeV_histo_C_mult_0_ep_1_syst.root";
  TFile *file_Alex_0 = TFile::Open(flow_proton_Alex_0.Data());
  g_Alex[0][2] = (TGraphAsymmErrors*)file_Alex_0->Get("tgae_Proton_v2_39GeV_mult_0_ep_1_err_0");

  TString flow_proton_Alex_1 = "/Users/mac/STAR/Data/BES/v2/Merge_antiProton_v2_39GeV_histo_C_mult_0_ep_1_syst.root";
  TFile *file_Alex_1 = TFile::Open(flow_proton_Alex_1.Data());
  g_Alex[1][2] = (TGraphAsymmErrors*)file_Alex_1->Get("tgae_antiProton_v2_39GeV_mult_0_ep_1_err_0");

  // Ratio
  TGraphAsymmErrors *g_ratio[charge_total][N_Species];
  for(Int_t charge = 0; charge < charge_total; charge++)
  {
    for(Int_t n_species = 0; n_species < N_Species; n_species++)
    {
      g_ratio[charge][n_species] = new TGraphAsymmErrors();
      for(Int_t i = 0; i < 11; i++)
      {
        Double_t x_Alex, y_Alex, err_Alex; 
	g_Alex[charge][n_species]->GetPoint(i,x_Alex,y_Alex);
	err_Alex = g_Alex[charge][n_species]->GetErrorYhigh(i);

	Double_t x_Xu, y_Xu, err_Xu;
	g_flow[charge][0][0][n_species]->GetPoint(i,x_Xu,y_Xu);
	err_Xu = g_flow[charge][0][0][n_species]->GetErrorYhigh(i);

	Double_t ratio_XA = y_Xu/y_Alex;
	Double_t err_XA = ErrDiv(y_Xu,y_Alex,err_Xu,err_Alex);

//	cout << "pt = " << x_Alex << ", ratio_XA = " << ratio_XA << ", err_XA = " << err_XA << endl;

	g_ratio[charge][n_species]->SetPoint(i,x_Alex,ratio_XA);
	g_ratio[charge][n_species]->SetPointError(i,0.0,0.0,0.0,0.0);
      }
    }
  }

  for(Int_t i = 0; i < 6; i++)
  {
    if(i < 3)
    {
      c_v2_Ratio_v2_6->cd(i+1);

      Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_flow[0][0][0][i],24,2,0.8);
      Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_Alex[0][i],24,1,0.8);

      if(i == 0)
      {
	Draw_TGAE_Point_new_Symbol(0.0,0.14,0.0,0.0,0.0,0.0,24,2,1.0);
	plotTopLegend((char*)"My Results",0.1,0.135,0.05,1,0.0,42,0);
	Draw_TGAE_Point_new_Symbol(0.0,0.16,0.0,0.0,0.0,0.0,24,1,1.0);
	plotTopLegend((char*)"Phys. Rev. C 88, 014902 (2013)",0.1,0.155,0.05,1,0.0,42,0);
      }
      if(i == 1)
      {
	plotTopLegend((char*)"200 GeV, 0-80%",0.1,0.155,0.05,1,0.0,42,0);
	plotTopLegend((char*)EtaGap[0].Data(),0.12,0.135,0.05,1,0.0,42,0);
      }

      plotTopLegend((char*)Name[0][i].Data(),0.2,0.18,0.08,1,0.0,42,0);

      c_v2_Ratio_v2_6->cd(i+4);
      Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_ratio[0][i],28,4,0.8);
    }
    else
    {
      c_v2_Ratio_v2_6->cd(i+4);
      Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_flow[1][0][0][i-3],24,2,0.8);
      Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_Alex[1][i-3],24,1,0.8);
      plotTopLegend((char*)Name[1][i-3].Data(),0.2,0.18,0.08,1,0.0,42,0);

      c_v2_Ratio_v2_6->cd(i+7);
      Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_ratio[1][i-3],28,4,0.8);
    }
  }
  c_v2_Ratio_v2_6->SaveAs("./figures/flow/c_v2_Ratio_v2_6.eps");

  //------------------------------------------------------------------------------------------------------------------------------------
  */
  // set canvas
  TCanvas* c_v3_Delta_v3_6[2];
  for(Int_t m = 0; m < 2; m++)
  {
    TString CanName = Form("c_v3_Delta_v3_6_%d",m);
    c_v3_Delta_v3_6[m] = new TCanvas(CanName.Data(),CanName.Data(),1400,10,1000,800);
    c_v3_Delta_v3_6[m]->SetFillColor(10);
    c_v3_Delta_v3_6[m]->SetTopMargin(0.0);
    c_v3_Delta_v3_6[m]->SetBottomMargin(0.3);
    c_v3_Delta_v3_6[m]->SetRightMargin(0.0);
    c_v3_Delta_v3_6[m]->SetLeftMargin(0.2);
    c_v3_Delta_v3_6[m]->Divide(3,4,0.0,0.0,10);
    TH1F* h_test_6p_3[12];

    for(Int_t j = 0; j < 12; j++)
    {
      TString HistName = "h_test_6p_3_";
      HistName += j;
      h_test_6p_3[j] = new TH1F(HistName.Data(),HistName.Data(),1000,-0.5,6);
      for(Int_t i = 0; i < h_test_6p_3[j]->GetNbinsX(); i++)
      {
	h_test_6p_3[j]->SetBinContent(i,-10);
      }
    }

    Int_t x_looop  = 0;
    Int_t y_looop2 = 0;
    Int_t y_looop[12] = {0,0,0,1,1,1,2,2,2,3,3,3};
    for(Int_t i = 0; i < 12; i++)
    {
      Float_t pad_factor = 1.0;
      if(i >= 6 && i <= 8) pad_factor = 1.1;
      if(i <= 2) {y_looop2 = 0;}
      else y_looop2 = 1;
      c_v3_Delta_v3_6[m]->cd(i+1);
      c_v3_Delta_v3_6[m]->cd(i+1);
      c_v3_Delta_v3_6[m]->cd(i+1)->SetTicks(1,1);
      c_v3_Delta_v3_6[m]->cd(i+1)->SetGrid(0,0);
      c_v3_Delta_v3_6[m]->cd(i+1)->SetFillColor(10);
      if(i >= 3 && i <= 5) c_v3_Delta_v3_6[m]->cd(i+1)->SetBottomMargin(0.01);
      if(i >= 9 && i <= 11) c_v3_Delta_v3_6[m]->cd(i+1)->SetBottomMargin(0.45);

      if(i < 3) c_v3_Delta_v3_6[m]->cd(i+1)->SetTopMargin(0.1);
      if(i == 2 || i == 5 || i == 8 || i == 11) c_v3_Delta_v3_6[m]->cd(i+1)->SetRightMargin(0.1);
      if(i == 3 || i == 6 || i == 9) x_looop = 0;
      Float_t x1_array[3] = {0,0.387,0.678};
      Float_t x2_array[3] = {x1_array[1],x1_array[2],1.0};
      Float_t y1_array[4] = {0.64,0.526+0.005,0.2,0.0};
      Float_t y2_array[4] = {1.0,y1_array[0],y1_array[1]-0.005,y1_array[2]};
      c_v3_Delta_v3_6[m]->cd(i+1)->SetPad(x1_array[x_looop],y1_array[y_looop[i]],x2_array[x_looop],y2_array[y_looop[i]]); // x1, y1, x2, y2
      //    cout << "i = " << i << ", x1 = " << x1_array[x_looop] << ", x2 = " << x2_array[x_looop] << ", y1 = " << y1_array[y_looop[i]]
      //      << ", y2 = " << y2_array[y_looop[i]] << endl;

      h_test_6p_3[i]->SetStats(0);
      h_test_6p_3[i]->SetTitle("");
      h_test_6p_3[i]->GetXaxis()->SetTitleOffset(1.2);
      h_test_6p_3[i]->GetYaxis()->SetTitleOffset(1.26);
      h_test_6p_3[i]->GetYaxis()->SetLabelOffset(0.01);
      h_test_6p_3[i]->GetXaxis()->SetLabelSize(0.09*label_size_factor);
      h_test_6p_3[i]->GetYaxis()->SetLabelSize(0.09*label_size_factor);
      h_test_6p_3[i]->GetXaxis()->SetTitleSize(0.09*label_size_factor);
      h_test_6p_3[i]->GetYaxis()->SetTitleSize(0.09*label_size_factor);
      h_test_6p_3[i]->GetXaxis()->SetNdivisions(505,'N');
      h_test_6p_3[i]->GetYaxis()->SetNdivisions(505,'N');
      h_test_6p_3[i]->GetXaxis()->CenterTitle();
      h_test_6p_3[i]->GetYaxis()->CenterTitle();
      h_test_6p_3[i]->GetXaxis()->SetRangeUser(-0.2,3.3);
      if(i < 3)
      {
	h_test_6p_3[i]->GetYaxis()->SetRangeUser(-0.02,0.21);
      }
      else
      {
	h_test_6p_3[i]->GetYaxis()->SetRangeUser(-0.02,0.14);
      }
      h_test_6p_3[i]->GetYaxis()->SetTitle("");
      h_test_6p_3[i]->GetXaxis()->SetTitle("");

      if(i == 3 || i == 4 || i == 5) {h_test_6p_3[i]->SetTickLength(0.04,"X");}
      if(i == 2 || i == 8) {h_test_6p_3[i]->SetTickLength(0.04,"Y");}
      if(i == 1 || i == 7) {h_test_6p_3[i]->SetTickLength(0.042,"Y");}
      if(i == 4)  {h_test_6p_3[i]->SetTickLength(0.04,"Y");}
      if(i == 5)  {h_test_6p_3[i]->SetTickLength(0.04,"Y");}
      if(i == 9)  {h_test_6p_3[i]->SetTickLength(0.052,"Y");}
      if(i == 10) {h_test_6p_3[i]->SetTickLength(0.075,"Y");}
      if(i == 11) {h_test_6p_3[i]->SetTickLength(0.07,"Y");}

      if(i == 3 || i == 4 || i == 5 || i == 9 || i == 10 || i == 11)
      {
	if(m == 0)
	{
	  h_test_6p_3[i]->GetYaxis()->SetRangeUser(-0.015,0.015);
	}
	if(m == 1)
	{
	  h_test_6p_3[i]->GetYaxis()->SetRangeUser(-0.003,0.007);
	}
	h_test_6p_3[i]->GetYaxis()->SetNdivisions(503,'N');
      }
      if(i >= 9)
      {
	h_test_6p_3[i]->GetXaxis()->SetTitleOffset(1.05);
	h_test_6p_3[i]->GetXaxis()->SetLabelOffset(0.008);
	h_test_6p_3[i]->GetXaxis()->SetLabelSize(0.14*label_size_factor);
	h_test_6p_3[i]->GetXaxis()->SetTitleSize(0.16*label_size_factor);
      }
      if(i == 3)
      {
	h_test_6p_3[i]->GetYaxis()->SetTitleOffset(0.43);
	h_test_6p_3[i]->GetYaxis()->SetLabelOffset(0.01);
	h_test_6p_3[i]->GetYaxis()->SetLabelSize(0.27*label_size_factor);
	h_test_6p_3[i]->GetYaxis()->SetTitleSize(0.27*label_size_factor);
      }
      if(i == 6)
      {
	h_test_6p_3[i]->GetYaxis()->SetTitleOffset(1.1);
	h_test_6p_3[i]->GetYaxis()->SetLabelOffset(0.01);
	h_test_6p_3[i]->GetYaxis()->SetLabelSize(0.1*label_size_factor);
	h_test_6p_3[i]->GetYaxis()->SetTitleSize(0.1*label_size_factor);
      }
      if(i == 9)
      {
	h_test_6p_3[i]->GetYaxis()->SetTitleOffset(0.72);
	h_test_6p_3[i]->GetYaxis()->SetLabelOffset(0.01);
	h_test_6p_3[i]->GetYaxis()->SetLabelSize(0.16*label_size_factor);
	h_test_6p_3[i]->GetYaxis()->SetTitleSize(0.16*label_size_factor);
      }
      if(i == 0) h_test_6p_3[i]->GetYaxis()->SetTitle("v_{2}");
      if(i == 6) h_test_6p_3[i]->GetYaxis()->SetTitle("v_{3}");
      if(i == 3) h_test_6p_3[i]->GetYaxis()->SetTitle("#Deltav_{2}");
      if(i == 9) h_test_6p_3[i]->GetYaxis()->SetTitle("#Deltav_{3}");
      if(i == 10) h_test_6p_3[i]->GetXaxis()->SetTitle("p_{T} (GeV/c)");

      h_test_6p_3[i]->DrawCopy("h");
      if(i == 0 || i == 1 || i == 2 || i == 6 || i == 7 || i == 8)
      {
	Float_t x_pos[3] = {0.3,0.07,0.068};
	Float_t y_pos[2] = {0.77,0.86};
	//      plotTopLegend((char*)"blubb",x_pos[x_looop],y_pos[y_looop2],0.09*pad_factor,1,0.0,42,1);
      }
      PlotLine(-0.2,3.1,0.0,0.0,1,1,2); // x1_val, x2_val, y1_val, y2_val, Line_Col, LineWidth, LineStyle
      x_looop++;
    }
  }
  //------------------------------------------------------------------------------------------------------------------------------------
  // v3 results

  // difference between particles and anti-particles
  TGraphAsymmErrors *g_diff_pm[eta_total][flow_total][N_Species];
  for(Int_t flow_bin = 0; flow_bin < flow_total; flow_bin++)
  {
    for(Int_t eta_bin = 0; eta_bin < eta_total; eta_bin++)
    {
      for(Int_t n_species = 0; n_species < N_Species; n_species++)
      {
	g_diff_pm[eta_bin][flow_bin][n_species] = new TGraphAsymmErrors();
	for(Int_t i = 0; i < 15; i++)
	{
	  Double_t x_pos, y_pos, err_pos;
	  g_flow[0][eta_bin][flow_bin][n_species]->GetPoint(i,x_pos,y_pos);
	  err_pos = g_flow[0][eta_bin][flow_bin][n_species]->GetErrorYhigh(i);

	  Double_t x_neg, y_neg, err_neg;
	  g_flow[1][eta_bin][flow_bin][n_species]->GetPoint(i,x_neg,y_neg);
	  err_neg = g_flow[1][eta_bin][flow_bin][n_species]->GetErrorYhigh(i);

	  Double_t y_diff = y_pos - y_neg;
	  Double_t err_diff = ErrorAdd(err_pos,err_neg);

	  g_diff_pm[eta_bin][flow_bin][n_species]->SetPoint(i,x_pos,y_diff);
	  g_diff_pm[eta_bin][flow_bin][n_species]->SetPointError(i,0.0,0.0,err_diff,err_diff);
	}
      }
    }
  }

  for(Int_t i = 0; i < 6; i++)
  {
    if(i < 3)
    {
      c_v3_Delta_v3_6[0]->cd(i+1);

      Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_flow[0][0][0][i],24,2,0.8);
      Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_flow[1][0][0][i],24,1,0.8);
//      Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_flow_hmasui[i],28,4,0.8);

      Draw_TGAE_Point_new_Symbol(0.1,0.18,0.0,0.0,0.0,0.0,24,2,1.0);
      plotTopLegend((char*)Name[0][i].Data(),0.2,0.177,0.06,1,0.0,42,0);

      Draw_TGAE_Point_new_Symbol(0.1,0.16,0.0,0.0,0.0,0.0,24,1,1.0);
      plotTopLegend((char*)Name[1][i].Data(),0.2,0.157,0.06,1,0.0,42,0);

//      Draw_TGAE_Point_new_Symbol(0.1,0.14,0.0,0.0,0.0,0.0,28,4,1.0);
//      plotTopLegend((char*)(Name[0][i]+"+"+Name[1][i]+" (Hiroshi)").Data(),0.2,0.137,0.06,1,0.0,42,0);
      if(i == 0)
      {
	plotTopLegend((char*)"Au+Au, 200 GeV, 0-80%",0.6,0.177,0.06,1,0.0,42,0);
	plotTopLegend((char*)EtaGap[0].Data(),0.62,0.157,0.06,1,0.0,42,0);
      }

      c_v3_Delta_v3_6[0]->cd(i+4);
      Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_diff_pm[0][0][i],28,4,0.8);
    }
    else
    {
      c_v3_Delta_v3_6[0]->cd(i+4);
      Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_flow[0][0][1][i-3],24,2,0.8);
      Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_flow[1][0][1][i-3],24,1,0.8);

      Draw_TGAE_Point_new_Symbol(0.1,0.12,0.0,0.0,0.0,0.0,24,2,1.0);
      plotTopLegend((char*)Name[0][i-3].Data(),0.2,0.117,0.07,1,0.0,42,0);
      Draw_TGAE_Point_new_Symbol(0.1,0.10,0.0,0.0,0.0,0.0,24,1,1.0);
      plotTopLegend((char*)Name[1][i-3].Data(),0.2,0.097,0.07,1,0.0,42,0);
      if(i == 3)
      {
	plotTopLegend((char*)"Au+Au, 200 GeV, 0-80%",0.6,0.117,0.07,1,0.0,42,0);
//	plotTopLegend((char*)EtaGap[1].Data(),0.62,0.097,0.07,1,0.0,42,0);
      }

      c_v3_Delta_v3_6[0]->cd(i+7);
      Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_diff_pm[0][1][i-3],28,4,0.8);
    }
  }
  c_v3_Delta_v3_6[0]->SaveAs("./figures/flow/c_v3_Delta_v3_6_pm.eps");

  /*
  // difference between different eta_gap
  TGraphAsymmErrors *g_diff_eta[charge_total][N_Species];
  for(Int_t charge = 0; charge < charge_total; charge++)
  {
    for(Int_t n_species = 0; n_species < N_Species; n_species++)
    {
      g_diff_eta[charge][n_species] = new TGraphAsymmErrors();
      for(Int_t i = 0; i < 15; i++)
      {
        Double_t x_0, y_0, err_0;
	g_flow[charge][0][1][n_species]->GetPoint(i,x_0,y_0);
	err_0 = g_flow[charge][0][1][n_species]->GetErrorYhigh(i);

        Double_t x_1, y_1, err_1;
	g_flow[charge][1][1][n_species]->GetPoint(i,x_1,y_1);
	err_1 = g_flow[charge][1][1][n_species]->GetErrorYhigh(i);

	Double_t y_diff = y_0 - y_1;
	Double_t err_diff = ErrorAdd(err_0,err_1);
	g_diff_eta[charge][n_species]->SetPoint(i,x_0,y_diff);
	g_diff_eta[charge][n_species]->SetPointError(i,0.0,0.0,0.0,0.0);
      }
    }
  }

  for(Int_t i = 0; i < 6; i++)
  {
    if(i < 3)
    {
      c_v3_Delta_v3_6[1]->cd(i+1);

      Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_flow[0][0][1][i],24,2,0.8);
      Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_flow[0][1][1][i],24,1,0.8);

      plotTopLegend((char*)Name[0][i].Data(),0.2,0.12,0.07,1,0.0,42,0);

      Draw_TGAE_Point_new_Symbol(0.1,0.10,0.0,0.0,0.0,0.0,24,2,1.0);
      plotTopLegend((char*)EtaGap[0].Data(),0.2,0.097,0.06,1,0.0,42,0);
      Draw_TGAE_Point_new_Symbol(0.1,0.08,0.0,0.0,0.0,0.0,24,1,1.0);
      plotTopLegend((char*)EtaGap[1].Data(),0.2,0.077,0.06,1,0.0,42,0);
      if(i == 0)
      {
	plotTopLegend((char*)"Au+Au, 200 GeV, 0-80%",0.7,0.12,0.06,1,0.0,42,0);
      }

      c_v3_Delta_v3_6[1]->cd(i+4);
      Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_diff_eta[0][i],28,4,0.8);
    }
    else
    {
      c_v3_Delta_v3_6[1]->cd(i+4);
      Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_flow[1][0][1][i-3],24,2,0.8);
      Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_flow[1][1][1][i-3],24,1,0.8);

      plotTopLegend((char*)Name[1][i-3].Data(),0.2,0.12,0.07,1,0.0,42,0);

      Draw_TGAE_Point_new_Symbol(0.1,0.10,0.0,0.0,0.0,0.0,24,2,1.0);
      plotTopLegend((char*)EtaGap[0].Data(),0.2,0.097,0.07,1,0.0,42,0);
      Draw_TGAE_Point_new_Symbol(0.1,0.08,0.0,0.0,0.0,0.0,24,1,1.0);
      plotTopLegend((char*)EtaGap[1].Data(),0.2,0.077,0.07,1,0.0,42,0);
      if(i == 3)
      {
	plotTopLegend((char*)"Au+Au, 200 GeV, 0-80%",0.7,0.12,0.07,1,0.0,42,0);
      }

      c_v3_Delta_v3_6[1]->cd(i+7);
      Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_diff_eta[1][i-3],28,4,0.8);
    }
  }
  c_v3_Delta_v3_6[1]->SaveAs("./figures/flow/c_v3_Delta_v3_6_etagap.eps");
  */

  // ncq scaling
  TGraphAsymmErrors *g_ncq_scale_pt[charge_total][eta_total][N_Species]; // pt/nq, v3/nq
  TGraphAsymmErrors *g_ncq_order_pt[charge_total][eta_total][N_Species]; // pt/nq, v3/nq^{3/2}
  TGraphAsymmErrors *g_flow_kt[charge_total][eta_total][N_Species]; // kt, v3
  TGraphAsymmErrors *g_ncq_scale_kt[charge_total][eta_total][N_Species]; // kt/nq, v3/nq
  TGraphAsymmErrors *g_ncq_order_kt[charge_total][eta_total][N_Species]; // kt/nq, v3/nq^{3/2}
  for(Int_t charge = 0; charge < charge_total; charge++)
  {
    for(Int_t eta_bin = 0; eta_bin < eta_total; eta_bin++)
    {
      for(Int_t n_species = 0; n_species < N_Species; n_species++)
      {
	g_ncq_scale_pt[charge][eta_bin][n_species] = new TGraphAsymmErrors();
	g_ncq_order_pt[charge][eta_bin][n_species] = new TGraphAsymmErrors();
	g_flow_kt[charge][eta_bin][n_species]      = new TGraphAsymmErrors();
	g_ncq_scale_kt[charge][eta_bin][n_species] = new TGraphAsymmErrors();
	g_ncq_order_kt[charge][eta_bin][n_species] = new TGraphAsymmErrors();
	for(Int_t i = 0; i < 15; i++)
	{
	  Double_t pt, v3, err;
	  g_flow[charge][eta_bin][1][n_species]->GetPoint(i,pt,v3);
	  err = g_flow[charge][eta_bin][1][n_species]->GetErrorYhigh(i);

	  g_ncq_scale_pt[charge][eta_bin][n_species]->SetPoint(i,pt/nq[n_species],v3/nq[n_species]);
	  g_ncq_scale_pt[charge][eta_bin][n_species]->SetPointError(i,0.0,0.0,err/nq[n_species],err/nq[n_species]);

	  Double_t Nq = TMath::Power(nq[n_species],3.0/2.0);
	  g_ncq_order_pt[charge][eta_bin][n_species]->SetPoint(i,pt/nq[n_species],v3/Nq);
	  g_ncq_order_pt[charge][eta_bin][n_species]->SetPointError(i,0.0,0.0,err/Nq,err/Nq);

	  Double_t kt = TMath::Sqrt(pt*pt+mass[n_species]*mass[n_species]) - mass[n_species];
	  g_flow_kt[charge][eta_bin][n_species]->SetPoint(i,kt,v3);
	  g_flow_kt[charge][eta_bin][n_species]->SetPointError(i,0.0,0.0,err,err);

	  g_ncq_scale_kt[charge][eta_bin][n_species]->SetPoint(i,kt/nq[n_species],v3/nq[n_species]);
	  g_ncq_scale_kt[charge][eta_bin][n_species]->SetPointError(i,0.0,0.0,err/nq[n_species],err/nq[n_species]);

	  g_ncq_order_kt[charge][eta_bin][n_species]->SetPoint(i,kt/nq[n_species],v3/Nq);
	  g_ncq_order_kt[charge][eta_bin][n_species]->SetPointError(i,0.0,0.0,err/Nq,err/Nq);
	}
      }
    }
  }

  // v3 vs pt
  TCanvas *c_v3_pt[eta_total];
  for(Int_t eta_bin = 0; eta_bin < eta_total; eta_bin++)
  {
    TString CanName = Form("c_v3_pt_%d",eta_bin);
    c_v3_pt[eta_bin] = new TCanvas(CanName.Data(),CanName.Data(),1400,10,1100,500);
    c_v3_pt[eta_bin]->Divide(2,1,0.0,0.0,10);;
    for(Int_t charge = 0; charge < charge_total; charge++)
    {
      c_v3_pt[eta_bin]->cd(charge+1);
      c_v3_pt[eta_bin]->cd(charge+1)->SetTicks(1,1);
      c_v3_pt[eta_bin]->cd(charge+1)->SetGrid(0,0);
      c_v3_pt[eta_bin]->cd(charge+1)->SetBottomMargin(0.15);
      if(charge == 0) 
      {
	c_v3_pt[eta_bin]->cd(charge+1)->SetLeftMargin(0.15);
	pos_neg_dummy->SetTickLength(0.034,"Y");
      }
      if(charge == 1) 
      {
	c_v3_pt[eta_bin]->cd(charge+1)->SetRightMargin(0.055);
	pos_neg_dummy->SetTickLength(0.04,"Y");
      }
      pos_neg_dummy->GetYaxis()->SetTitle(TitleY[1].Data());
      pos_neg_dummy->GetYaxis()->SetRangeUser(-0.02,0.14);
      pos_neg_dummy->GetXaxis()->SetRangeUser(-0.1,3.3);
      pos_neg_dummy->DrawCopy("h");
      PlotLine(-0.1,3.1,0.0,0.0,1,2,2);

      Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_flow[charge][eta_bin][1][0],Draw_Style[0],Draw_Color[0],0.8);
      Draw_TGAE_Point_new_Symbol(0.2,0.12,0.0,0.0,0.0,0.0,Draw_Style[0],Draw_Color[0],1.0);
      plotTopLegend((char*)Name[charge][0].Data(),0.3,0.117,0.06,1,0.0,42,0);

      Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_flow[charge][eta_bin][1][1],Draw_Style[1],Draw_Color[1],0.8);
      Draw_TGAE_Point_new_Symbol(0.2,0.10,0.0,0.0,0.0,0.0,Draw_Style[1],Draw_Color[1],1.0);
      plotTopLegend((char*)Name[charge][1].Data(),0.3,0.097,0.06,1,0.0,42,0);

      Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_flow[charge][eta_bin][1][2],Draw_Style[2],Draw_Color[2],0.8);
      Draw_TGAE_Point_new_Symbol(0.2,0.08,0.0,0.0,0.0,0.0,Draw_Style[2],Draw_Color[2],1.0);
      plotTopLegend((char*)Name[charge][2].Data(),0.3,0.077,0.06,1,0.0,42,0);

      if(charge == 0)
      {
	plotTopLegend((char*)"Au+Au, 200 GeV, 0-80%",1.4,0.12,0.05,1,0.0,42,0);
      }
      if(charge == 1)
      {
	plotTopLegend((char*)EtaGap[eta_bin].Data(),2.1,0.12,0.05,1,0.0,42,0);
      }
    }
    c_v3_pt[eta_bin]->SaveAs(("./figures/flow/"+CanName+".eps").Data());
  }

  // v3/nq vs pt/nq
  TCanvas *c_scale_pt[eta_total];
  for(Int_t eta_bin = 0; eta_bin < eta_total; eta_bin++)
  {
    TString CanName = Form("c_scale_pt_%d",eta_bin);
    c_scale_pt[eta_bin] = new TCanvas(CanName.Data(),CanName.Data(),1400,10,1100,500);
    c_scale_pt[eta_bin]->Divide(2,1,0.0,0.0,10);;
    for(Int_t charge = 0; charge < charge_total; charge++)
    {
      c_scale_pt[eta_bin]->cd(charge+1);
      c_scale_pt[eta_bin]->cd(charge+1)->SetTicks(1,1);
      c_scale_pt[eta_bin]->cd(charge+1)->SetGrid(0,0);
      c_scale_pt[eta_bin]->cd(charge+1)->SetBottomMargin(0.15);
      if(charge == 0) 
      {
	c_scale_pt[eta_bin]->cd(charge+1)->SetLeftMargin(0.15);
	pos_neg_dummy->SetTickLength(0.034,"Y");
      }
      if(charge == 1) 
      {
	c_scale_pt[eta_bin]->cd(charge+1)->SetRightMargin(0.055);
	pos_neg_dummy->SetTickLength(0.04,"Y");
      }
      pos_neg_dummy->GetYaxis()->SetTitle("v_{3}/n_{q}");
      pos_neg_dummy->GetXaxis()->SetTitle("p_{T}/n_{q} (GeV/c)");
      pos_neg_dummy->GetYaxis()->SetRangeUser(-0.02/2.0,0.14/2.0);
      pos_neg_dummy->GetXaxis()->SetRangeUser(-0.10/2.0,3.30/2.0);
      pos_neg_dummy->DrawCopy("h");
      PlotLine(-0.05,1.55,0.0,0.0,1,2,2);

      Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_ncq_scale_pt[charge][eta_bin][0],Draw_Style[0],Draw_Color[0],0.8);
      Draw_TGAE_Point_new_Symbol(0.2/2.0,0.12/2.0,0.0,0.0,0.0,0.0,Draw_Style[0],Draw_Color[0],1.0);
      plotTopLegend((char*)Name[charge][0].Data(),0.30/2.0,0.117/2.0,0.06,1,0.0,42,0);

      Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_ncq_scale_pt[charge][eta_bin][1],Draw_Style[1],Draw_Color[1],0.8);
      Draw_TGAE_Point_new_Symbol(0.2/2.0,0.10/2.0,0.0,0.0,0.0,0.0,Draw_Style[1],Draw_Color[1],1.0);
      plotTopLegend((char*)Name[charge][1].Data(),0.30/2.0,0.097/2.0,0.06,1,0.0,42,0);

      Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_ncq_scale_pt[charge][eta_bin][2],Draw_Style[2],Draw_Color[2],0.8);
      Draw_TGAE_Point_new_Symbol(0.2/2.0,0.08/2.0,0.0,0.0,0.0,0.0,Draw_Style[2],Draw_Color[2],1.0);
      plotTopLegend((char*)Name[charge][2].Data(),0.30/2.0,0.077/2.0,0.06,1,0.0,42,0);

      if(charge == 0)
      {
	plotTopLegend((char*)"Au+Au, 200 GeV, 0-80%",0.7,0.06,0.05,1,0.0,42,0);
      }
      if(charge == 1)
      {
	plotTopLegend((char*)EtaGap[eta_bin].Data(),1.05,0.06,0.05,1,0.0,42,0);
      }
    }
    c_scale_pt[eta_bin]->SaveAs(("./figures/flow/"+CanName+".eps").Data());
  }

  // v3/nq^{3/2} vs pt/nq
  TCanvas *c_order_pt[eta_total];
  for(Int_t eta_bin = 0; eta_bin < eta_total; eta_bin++)
  {
    TString CanName = Form("c_order_pt_%d",eta_bin);
    c_order_pt[eta_bin] = new TCanvas(CanName.Data(),CanName.Data(),1400,10,1100,500);
    c_order_pt[eta_bin]->Divide(2,1,0.0,0.0,10);;
    for(Int_t charge = 0; charge < charge_total; charge++)
    {
      c_order_pt[eta_bin]->cd(charge+1);
      c_order_pt[eta_bin]->cd(charge+1)->SetTicks(1,1);
      c_order_pt[eta_bin]->cd(charge+1)->SetGrid(0,0);
      c_order_pt[eta_bin]->cd(charge+1)->SetBottomMargin(0.15);
      if(charge == 0) 
      {
	c_order_pt[eta_bin]->cd(charge+1)->SetLeftMargin(0.15);
	pos_neg_dummy->SetTickLength(0.034,"Y");
      }
      if(charge == 1) 
      {
	c_order_pt[eta_bin]->cd(charge+1)->SetRightMargin(0.055);
	pos_neg_dummy->SetTickLength(0.04,"Y");
      }
      pos_neg_dummy->GetYaxis()->SetTitle("v_{3}/n_{q}^{3/2}");
      pos_neg_dummy->GetXaxis()->SetTitle("p_{T}/n_{q} (GeV/c)");
      pos_neg_dummy->GetYaxis()->SetRangeUser(-0.02/TMath::Power(2.0,3.0/2.0),0.14/TMath::Power(2.0,3.0/2.0));
      pos_neg_dummy->GetXaxis()->SetRangeUser(-0.10/2.0,3.30/2.0);
      pos_neg_dummy->DrawCopy("h");
      PlotLine(-0.05,1.55,0.0,0.0,1,2,2);

      Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_ncq_order_pt[charge][eta_bin][0],Draw_Style[0],Draw_Color[0],0.8);
      Draw_TGAE_Point_new_Symbol(0.2/2.0,0.12/TMath::Power(2.0,3.0/2.0),0.0,0.0,0.0,0.0,Draw_Style[0],Draw_Color[0],1.0);
      plotTopLegend((char*)Name[charge][0].Data(),0.30/2.0,0.117/TMath::Power(2.0,3.0/2.0),0.06,1,0.0,42,0);

      Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_ncq_order_pt[charge][eta_bin][1],Draw_Style[1],Draw_Color[1],0.8);
      Draw_TGAE_Point_new_Symbol(0.2/2.0,0.10/TMath::Power(2.0,3.0/2.0),0.0,0.0,0.0,0.0,Draw_Style[1],Draw_Color[1],1.0);
      plotTopLegend((char*)Name[charge][1].Data(),0.30/2.0,0.097/TMath::Power(2.0,3.0/2.0),0.06,1,0.0,42,0);

      Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_ncq_order_pt[charge][eta_bin][2],Draw_Style[2],Draw_Color[2],0.8);
      Draw_TGAE_Point_new_Symbol(0.2/2.0,0.08/TMath::Power(2.0,3.0/2.0),0.0,0.0,0.0,0.0,Draw_Style[2],Draw_Color[2],1.0);
      plotTopLegend((char*)Name[charge][2].Data(),0.30/2.0,0.077/TMath::Power(2.0,3.0/2.0),0.06,1,0.0,42,0);

      if(charge == 0)
      {
	plotTopLegend((char*)"Au+Au, 200 GeV, 0-80%",0.7,0.12/TMath::Power(2.0,3.0/2.0),0.05,1,0.0,42,0);
      }
      if(charge == 1)
      {
	plotTopLegend((char*)EtaGap[eta_bin].Data(),1.05,0.12/TMath::Power(2.0,3.0/2.0),0.05,1,0.0,42,0);
      }
    }
    c_order_pt[eta_bin]->SaveAs(("./figures/flow/"+CanName+".eps").Data());
  }

  // v3 vs kt
  TCanvas *c_v3_kt[eta_total];
  for(Int_t eta_bin = 0; eta_bin < eta_total; eta_bin++)
  {
    TString CanName = Form("c_v3_kt_%d",eta_bin);
    c_v3_kt[eta_bin] = new TCanvas(CanName.Data(),CanName.Data(),1400,10,1100,500);
    c_v3_kt[eta_bin]->Divide(2,1,0.0,0.0,10);;
    for(Int_t charge = 0; charge < charge_total; charge++)
    {
      c_v3_kt[eta_bin]->cd(charge+1);
      c_v3_kt[eta_bin]->cd(charge+1)->SetTicks(1,1);
      c_v3_kt[eta_bin]->cd(charge+1)->SetGrid(0,0);
      c_v3_kt[eta_bin]->cd(charge+1)->SetBottomMargin(0.15);
      if(charge == 0) 
      {
	c_v3_kt[eta_bin]->cd(charge+1)->SetLeftMargin(0.15);
	pos_neg_dummy->SetTickLength(0.034,"Y");
      }
      if(charge == 1) 
      {
	c_v3_kt[eta_bin]->cd(charge+1)->SetRightMargin(0.055);
	pos_neg_dummy->SetTickLength(0.04,"Y");
      }
      pos_neg_dummy->GetYaxis()->SetTitle(TitleY[1].Data());
      pos_neg_dummy->GetXaxis()->SetTitle("m_{T}-m_{0} (GeV/c^{2})");
      pos_neg_dummy->GetYaxis()->SetRangeUser(-0.02,0.14);
      pos_neg_dummy->GetXaxis()->SetRangeUser(-0.1,3.3);
      pos_neg_dummy->DrawCopy("h");
      PlotLine(-0.1,3.1,0.0,0.0,1,2,2);

      Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_flow_kt[charge][eta_bin][0],Draw_Style[0],Draw_Color[0],0.8);
      Draw_TGAE_Point_new_Symbol(0.2,0.12,0.0,0.0,0.0,0.0,Draw_Style[0],Draw_Color[0],1.0);
      plotTopLegend((char*)Name[charge][0].Data(),0.3,0.117,0.06,1,0.0,42,0);

      Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_flow_kt[charge][eta_bin][1],Draw_Style[1],Draw_Color[1],0.8);
      Draw_TGAE_Point_new_Symbol(0.2,0.10,0.0,0.0,0.0,0.0,Draw_Style[1],Draw_Color[1],1.0);
      plotTopLegend((char*)Name[charge][1].Data(),0.3,0.097,0.06,1,0.0,42,0);

      Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_flow_kt[charge][eta_bin][2],Draw_Style[2],Draw_Color[2],0.8);
      Draw_TGAE_Point_new_Symbol(0.2,0.08,0.0,0.0,0.0,0.0,Draw_Style[2],Draw_Color[2],1.0);
      plotTopLegend((char*)Name[charge][2].Data(),0.3,0.077,0.06,1,0.0,42,0);

      if(charge == 0)
      {
	plotTopLegend((char*)"Au+Au, 200 GeV, 0-80%",1.4,0.12,0.05,1,0.0,42,0);
      }
      if(charge == 1)
      {
	plotTopLegend((char*)EtaGap[eta_bin].Data(),2.1,0.12,0.05,1,0.0,42,0);
      }
    }
    c_v3_kt[eta_bin]->SaveAs(("./figures/flow/"+CanName+".eps").Data());
  }

  // v3/nq vs kt/nq
  TCanvas *c_scale_kt[eta_total];
  for(Int_t eta_bin = 0; eta_bin < eta_total; eta_bin++)
  {
    TString CanName = Form("c_scale_kt_%d",eta_bin);
    c_scale_kt[eta_bin] = new TCanvas(CanName.Data(),CanName.Data(),1400,10,1100,500);
    c_scale_kt[eta_bin]->Divide(2,1,0.0,0.0,10);;
    for(Int_t charge = 0; charge < charge_total; charge++)
    {
      c_scale_kt[eta_bin]->cd(charge+1);
      c_scale_kt[eta_bin]->cd(charge+1)->SetTicks(1,1);
      c_scale_kt[eta_bin]->cd(charge+1)->SetGrid(0,0);
      c_scale_kt[eta_bin]->cd(charge+1)->SetBottomMargin(0.15);
      if(charge == 0) 
      {
	c_scale_kt[eta_bin]->cd(charge+1)->SetLeftMargin(0.15);
	pos_neg_dummy->SetTickLength(0.034,"Y");
      }
      if(charge == 1) 
      {
	c_scale_kt[eta_bin]->cd(charge+1)->SetRightMargin(0.055);
	pos_neg_dummy->SetTickLength(0.04,"Y");
      }
      pos_neg_dummy->GetYaxis()->SetTitle("v_{3}/n_{q}");
      pos_neg_dummy->GetXaxis()->SetTitle("(m_{T}-m_{0})/n_{q} (GeV/c^{2})");
      pos_neg_dummy->GetYaxis()->SetRangeUser(-0.02/2.0,0.14/2.0);
      pos_neg_dummy->GetXaxis()->SetRangeUser(-0.10/2.0,3.30/2.0);
      pos_neg_dummy->DrawCopy("h");
      PlotLine(-0.05,1.55,0.0,0.0,1,2,2);

      Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_ncq_scale_kt[charge][eta_bin][0],Draw_Style[0],Draw_Color[0],0.8);
      Draw_TGAE_Point_new_Symbol(0.2/2.0,0.12/2.0,0.0,0.0,0.0,0.0,Draw_Style[0],Draw_Color[0],1.0);
      plotTopLegend((char*)Name[charge][0].Data(),0.30/2.0,0.117/2.0,0.06,1,0.0,42,0);

      Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_ncq_scale_kt[charge][eta_bin][1],Draw_Style[1],Draw_Color[1],0.8);
      Draw_TGAE_Point_new_Symbol(0.2/2.0,0.10/2.0,0.0,0.0,0.0,0.0,Draw_Style[1],Draw_Color[1],1.0);
      plotTopLegend((char*)Name[charge][1].Data(),0.30/2.0,0.097/2.0,0.06,1,0.0,42,0);

      Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_ncq_scale_kt[charge][eta_bin][2],Draw_Style[2],Draw_Color[2],0.8);
      Draw_TGAE_Point_new_Symbol(0.2/2.0,0.08/2.0,0.0,0.0,0.0,0.0,Draw_Style[2],Draw_Color[2],1.0);
      plotTopLegend((char*)Name[charge][2].Data(),0.30/2.0,0.077/2.0,0.06,1,0.0,42,0);

      if(charge == 0)
      {
	plotTopLegend((char*)"Au+Au, 200 GeV, 0-80%",0.7,0.06,0.05,1,0.0,42,0);
      }
      if(charge == 1)
      {
	plotTopLegend((char*)EtaGap[eta_bin].Data(),1.05,0.06,0.05,1,0.0,42,0);
      }
    }
    c_scale_kt[eta_bin]->SaveAs(("./figures/flow/"+CanName+".eps").Data());
  }

  // v3/nq^{3/2} vs kt/nq
  TCanvas *c_order_kt[eta_total];
  for(Int_t eta_bin = 0; eta_bin < eta_total; eta_bin++)
  {
    TString CanName = Form("c_order_kt_%d",eta_bin);
    c_order_kt[eta_bin] = new TCanvas(CanName.Data(),CanName.Data(),1400,10,1100,500);
    c_order_kt[eta_bin]->Divide(2,1,0.0,0.0,10);;
    for(Int_t charge = 0; charge < charge_total; charge++)
    {
      c_order_kt[eta_bin]->cd(charge+1);
      c_order_kt[eta_bin]->cd(charge+1)->SetTicks(1,1);
      c_order_kt[eta_bin]->cd(charge+1)->SetGrid(0,0);
      c_order_kt[eta_bin]->cd(charge+1)->SetBottomMargin(0.15);
      if(charge == 0) 
      {
	c_order_kt[eta_bin]->cd(charge+1)->SetLeftMargin(0.15);
	pos_neg_dummy->SetTickLength(0.034,"Y");
      }
      if(charge == 1) 
      {
	c_order_kt[eta_bin]->cd(charge+1)->SetRightMargin(0.055);
	pos_neg_dummy->SetTickLength(0.04,"Y");
      }
      pos_neg_dummy->GetYaxis()->SetTitle("v_{3}/n_{q}^{3/2}");
      pos_neg_dummy->GetXaxis()->SetTitle("(m_{T}-m_{0})/n_{q} (GeV/c^{2})");
      pos_neg_dummy->GetYaxis()->SetRangeUser(-0.02/TMath::Power(2.0,3.0/2.0),0.14/TMath::Power(2.0,3.0/2.0));
      pos_neg_dummy->GetXaxis()->SetRangeUser(-0.10/2.0,3.30/2.0);
      pos_neg_dummy->DrawCopy("h");
      PlotLine(-0.05,1.55,0.0,0.0,1,2,2);

      Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_ncq_order_kt[charge][eta_bin][0],Draw_Style[0],Draw_Color[0],0.8);
      Draw_TGAE_Point_new_Symbol(0.2/2.0,0.12/TMath::Power(2.0,3.0/2.0),0.0,0.0,0.0,0.0,Draw_Style[0],Draw_Color[0],1.0);
      plotTopLegend((char*)Name[charge][0].Data(),0.30/2.0,0.117/TMath::Power(2.0,3.0/2.0),0.06,1,0.0,42,0);

      Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_ncq_order_kt[charge][eta_bin][1],Draw_Style[1],Draw_Color[1],0.8);
      Draw_TGAE_Point_new_Symbol(0.2/2.0,0.10/TMath::Power(2.0,3.0/2.0),0.0,0.0,0.0,0.0,Draw_Style[1],Draw_Color[1],1.0);
      plotTopLegend((char*)Name[charge][1].Data(),0.30/2.0,0.097/TMath::Power(2.0,3.0/2.0),0.06,1,0.0,42,0);

      Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_ncq_order_kt[charge][eta_bin][2],Draw_Style[2],Draw_Color[2],0.8);
      Draw_TGAE_Point_new_Symbol(0.2/2.0,0.08/TMath::Power(2.0,3.0/2.0),0.0,0.0,0.0,0.0,Draw_Style[2],Draw_Color[2],1.0);
      plotTopLegend((char*)Name[charge][2].Data(),0.30/2.0,0.077/TMath::Power(2.0,3.0/2.0),0.06,1,0.0,42,0);

      if(charge == 0)
      {
	plotTopLegend((char*)"Au+Au, 200 GeV, 0-80%",0.7,0.12/TMath::Power(2.0,3.0/2.0),0.05,1,0.0,42,0);
      }
      if(charge == 1)
      {
	plotTopLegend((char*)EtaGap[eta_bin].Data(),1.05,0.12/TMath::Power(2.0,3.0/2.0),0.05,1,0.0,42,0);
      }
    }
    c_order_kt[eta_bin]->SaveAs(("./figures/flow/"+CanName+".eps").Data());
  }
}

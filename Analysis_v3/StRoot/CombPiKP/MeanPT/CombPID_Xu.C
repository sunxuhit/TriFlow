#include "TFile.h"
#include "TMath.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TString.h"
#include "TF2.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "Math/MinimizerOptions.h"
#include "TFitter.h"
#include "TGaxis.h"
#include "TColor.h"
#include <iostream>

#include "/u/xusun/STAR/Analysis_v3/StRoot/CombPiKP/student_t_2d_fit.h"
#include "/u/xusun/STAR/Analysis_v3/StRoot/CombPiKP/student_t_2d_single.h"
#include "/u/xusun/STAR/Analysis_v3/StRoot/CombPiKP/student_t_2d_double.h"
#include "/u/xusun/STAR/Analysis_v3/StRoot/CombPiKP/student_t_1d_single.h"
#include "/u/xusun/STAR/Analysis_v3/StRoot/CombPiKP/student_t_1d_double.h"
#include "/u/xusun/STAR/Analysis_v3/StRoot/CombPiKP/draw.h"
#include "/u/xusun/STAR/Analysis_v3/StRoot/CombPiKP/student_t_combPID.h"

// mEnergy: 0 for 200GeV, 1 for 39GeV
// mCharge: 0 for positive charge, 1 for negative charge

void CombPID_Xu(Int_t mode = 1,Int_t mEnergy = 0, Int_t mCharge = 0, Int_t mCentrality = 0)
{
  //**************************** Set graphic style ***************************************
  gStyle->SetPalette(1);
  gStyle->SetCanvasColor(10);
  gStyle->SetFrameFillColor(10);
  gStyle->SetGridColor(0);
  gStyle->SetGridStyle(0);
  gStyle->SetGridWidth(0);
  //gStyle->SetFillColor(4);
  TGaxis::SetMaxDigits(4);
  gStyle->SetPadTopMargin(0.15);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadRightMargin(0.14);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetLabelSize(0.07,"X");
  gStyle->SetLabelSize(0.07,"Y");
  gStyle->SetTitleSize(0.07,"X");
  gStyle->SetTitleSize(0.07,"Y");

  const Int_t NRGBs = 5;
  const Int_t NCont = 255;
  Double_t stops[NRGBs]  = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t reds[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t greens[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blues[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  Int_t  FI = TColor::CreateGradientColorTable(NRGBs, stops, reds, greens, blues, NCont);
  gStyle->SetNumberContours(NCont);
  //**************************************************************************************
  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(50000);

  //-----------------------------------------------------------------------------------------
  const Int_t pt_total = 16;
  const Int_t cent_total = 4; // 4
  const Int_t cent_start[4] = {0,1,2,3};
  const Int_t cent_stop[4]  = {1,2,3,4};
  const Int_t charge_total = 2; // 2
  const Int_t charge_start = (mCharge == 0) ? 0 : 1; // charge_start: 0 for pos, 1 for neg
  const Int_t charge_stop  = (mCharge == 0) ? 1 : 2; // charge_stop:  1 for pos, 2 for neg
  const Int_t eta_total = 4; // 4
  const Int_t eta_start = 0; // 4
  const Int_t eta_stop = 1; // 4
  const Int_t phi_psi_total = 7; // 7
  const Int_t flow_total = 2; // 2
  
  // Initial parameters 200 GeV
  Float_t order_pion[pt_total]   = {3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,2.5,3.00,2.50,3.0,2.50,3.0};
  Float_t order_proton[pt_total] = {3.0,3.0,3.0,3.0,3.0,2.5,2.5,3.0,3.0,3.0,2.5,2.55,2.55,3.0,2.55,3.0};
  // pion
  Float_t x_pi[pt_total] = {0.0025,-0.0025,-0.0025,-0.0025,-0.0025,-0.0075,-0.0125,-0.005,-0.005,-0.00625,0.00875,0.02375,0.0283,0.0344,0.0649,0.1135};
  Float_t y_pi[pt_total] = {0.0005,0.00125,0.00125,-0.00025,0.002,0.00825,0.01075,0.02175,0.02825,0.03137,0.02625,0.03375,0.0114,0.0236,0.0236,0.0825};
  Float_t width_pi[pt_total] = {0.0,0.0,0.01214,0.01992,0.02972,0.04194,0.05657,0.07499,0.09672,0.1205,0.137,0.173,0.2807,0.2807,0.2807,0.2807};
  Float_t nu_pi[pt_total] = {0.0,0.0,6.912,7.673,7.503,8.736,11.20,12.47,16.76,15.81,50,37.8,76.27,76.27,76.27,76.27}; // etagap = 0

  // kaon
  Float_t x_k[pt_total] = { 0.2275,  0.2275, 0.2225, 0.2225,0.2225, 0.2225, 0.2325, 0.2470, 0.2530, 0.2787, 0.3500, 0.3970,0.3519,0.4130,0.4924,0.8652};
  Float_t y_k[pt_total] = {-0.0045,-0.00325,0.00125,0.00325,0.0060,0.00375,0.00525,0.02175,0.00225,0.03137,0.02625,0.03375,0.0114,0.0236,0.0236,0.0825};
  Float_t width_k[pt_total] = {0.008031,0.01045,0.01612,0.02472,0.03588,0.04931,0.06609,0.08716,0.1064,0.1295,0.150,0.176,0.3187,0.3187,0.3187,0.3187};

  // proton
  Float_t x_p[pt_total] = {0.8625, 0.8625,  0.8675,  0.8625,0.8525, 0.8305, 0.8025, 0.7750, 0.7450, 0.7213, 0.7287, 0.7620, 0.7794, 0.8343, 0.9565, 0.8955};
  Float_t y_p[pt_total] = {0.0375,0.02075,-0.01125,-0.05975,-0.126,-0.2033,-0.2862,-0.3553,-0.4268,-0.4614,-0.5062,-0.5587,-0.5015,-0.5071,-0.4981,-0.5295};
//  Float_t width_p_x[pt_total] = {0.03359,0.03034,0.03311,0.04130,0.05238,0.06557,0.08030,0.09742,0.11610,0.1384,0.158,0.1984,0.3567,0.3567,0.3567,0.3567};
  Float_t width_p_x[pt_total] = {0.03359,0.03034,0.03311,0.04130,0.05238,0.07200,0.08030,0.09100,0.11610,0.1384,0.158,0.1984,0.3567,0.3567,0.36,0.3567};
  Float_t width_p_y[pt_total] = {0.01217,0.02191,0.03077,0.04085,0.05326,0.06540,0.07996,0.09776,0.11900,0.1419,0.169,0.2052,0.3111,0.3111,0.32,0.3111};
  Float_t nu_p[pt_total] = {8.975,9.936,7.798,6.460,5.970,5.905,6.109,6.991,8.212,10.58,20.,47.92,49.,50.2,55,50.2};

  /*
  // Initial parameters 39 GeV for 00-80
  Float_t order_pion[pt_total]   = {3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.00,2.80,2.55,3.00,3.0};
  Float_t order_proton[pt_total] = {3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,2.55,2.50,2.85,2.50,3.0};
  // pion
  Float_t x_pi[pt_total] = {0.0025,-0.0025,-0.0025,-0.0025,-0.0025,-0.0075,-0.0125,-0.005,-0.005,-0.00625,0.00875,0.02375,0.0305,0.0344,0.0649,0.1135};
  Float_t y_pi[pt_total] = {0.0005,0.00125,0.00125,-0.00025,0.002,0.00825,0.01075,0.02175,0.02825,0.03137,0.02625,0.03375,0.0114,0.0236,0.0236,0.0825};
  Float_t width_pi[pt_total] = {0.0,0.0,0.01214,0.01992,0.02972,0.04194,0.05657,0.07499,0.09672,0.1205,0.137,0.173,0.2807,0.2807,0.2807,0.2807};
  Float_t nu_pi[pt_total] = {0.0,0.0,6.912,7.673,7.503,8.736,11.20,12.47,16.76,15.81,50,37.8,76.27,76.27,76.27,76.27}; // etagap = 0

  // kaon
  Float_t x_k[pt_total] = { 0.2275,  0.2275, 0.2225, 0.2225,0.2225, 0.2225, 0.2325, 0.2470, 0.2530, 0.2787, 0.3560, 0.3970,0.3515,0.4130,0.4924,0.8652};
  Float_t y_k[pt_total] = {-0.0045,-0.00325,0.00125,0.00325,0.0060,0.00375,0.00525,0.02175,0.00225,0.03137,0.02625,0.03375,0.0114,0.0236,0.0236,0.0825};
  Float_t width_k[pt_total] = {0.008031,0.01045,0.01612,0.02472,0.03588,0.04931,0.06609,0.08716,0.1064,0.1295,0.150,0.176,0.3187,0.3187,0.3187,0.3187};

  // proton
  Float_t x_p[pt_total] = {0.8625, 0.8625,  0.8675,  0.8625,0.8525, 0.8305, 0.7825, 0.7750, 0.7450, 0.7213, 0.7287, 0.7680, 0.7794, 0.7794, 0.8285, 0.8955};
  Float_t y_p[pt_total] = {0.0375,0.02075,-0.01125,-0.05975,-0.126,-0.2033,-0.2862,-0.3553,-0.4268,-0.4614,-0.5062,-0.5587,-0.5015,-0.4832,-0.4441,-0.5295};
  Float_t width_p_x[pt_total] = {0.03359,0.03034,0.03311,0.04130,0.05238,0.06557,0.08030,0.09742,0.11610,0.1384,0.158,0.1984,0.3067,0.3567,0.3265,0.3567};
  Float_t width_p_y[pt_total] = {0.01217,0.02191,0.03077,0.04085,0.05326,0.06540,0.07996,0.09776,0.11900,0.1419,0.169,0.2052,0.3111,0.3111,0.3111,0.3111};
  Float_t nu_p[pt_total] = {8.975,9.936,7.798,6.460,5.970,5.905,6.509,6.991,8.212,10.58,20.,47.92,50.,55.5,55,50.2};
  */

  /*
  // Initial parameters 39 GeV for 10-40
  Float_t order_pion[pt_total]   = {3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,2.5,3.0,3.00,2.50,3.00,2.50,3.0};
  Float_t order_proton[pt_total] = {3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,2.5,3.0,2.55,2.55,2.55,2.55,3.0};
  // pion
  Float_t x_pi[pt_total] = {0.0025,-0.0025,-0.0025,-0.0025,-0.0025,-0.0075,-0.0125,-0.005,-0.005,-0.00625,0.00875,0.02375,0.0305,0.0344,0.0649,0.1135};
  Float_t y_pi[pt_total] = {0.0005,0.00125,0.00125,-0.00025,0.002,0.00825,0.01075,0.02175,0.02825,0.03137,0.02625,0.03375,0.0114,0.0236,0.0236,0.0825};
  Float_t width_pi[pt_total] = {0.0,0.0,0.01214,0.01992,0.02972,0.04194,0.05657,0.07499,0.09672,0.1205,0.137,0.173,0.2807,0.2807,0.2807,0.2807};
  Float_t nu_pi[pt_total] = {0.0,0.0,6.912,7.673,7.503,8.736,11.20,12.47,16.76,15.81,50,37.8,76.27,76.27,76.27,76.27}; // etagap = 0

  // kaon
  Float_t x_k[pt_total] = { 0.2275,  0.2275, 0.2225, 0.2225,0.2225, 0.2225, 0.2325, 0.2470, 0.2530, 0.2787, 0.3560, 0.3970,0.3515,0.4130,0.4924,0.8652};
  Float_t y_k[pt_total] = {-0.0045,-0.00325,0.00125,0.00325,0.0060,0.00375,0.00525,0.02175,0.00225,0.03137,0.02625,0.03375,0.0114,0.0236,0.0236,0.0825};
  Float_t width_k[pt_total] = {0.008031,0.01045,0.01612,0.02472,0.03588,0.04931,0.06609,0.08716,0.1064,0.1295,0.150,0.176,0.3187,0.3187,0.3187,0.3187};

  // proton
  Float_t x_p[pt_total] = {0.8625, 0.8625,  0.8675,  0.8625,0.8525, 0.8305, 0.7825, 0.7750, 0.7450, 0.7213, 0.7287, 0.7680, 0.7794, 0.7794, 0.8285, 0.8955};
  Float_t y_p[pt_total] = {0.0375,0.02075,-0.01125,-0.05975,-0.126,-0.2033,-0.2862,-0.3553,-0.4268,-0.4614,-0.5062,-0.5587,-0.5015,-0.4832,-0.4441,-0.5295};
  Float_t width_p_x[pt_total] = {0.03359,0.03034,0.03311,0.04130,0.05238,0.06557,0.08030,0.09742,0.11610,0.1384,0.158,0.1984,0.3067,0.3567,0.3265,0.3567};
  Float_t width_p_y[pt_total] = {0.01217,0.02191,0.03077,0.04085,0.05326,0.06540,0.07996,0.09776,0.11900,0.1419,0.169,0.2052,0.3111,0.3111,0.3111,0.3111};
  Float_t nu_p[pt_total] = {8.975,9.936,7.798,6.460,5.970,5.905,6.509,6.991,8.212,10.58,20.,47.92,50.,53.5,55,50.2};
  */

  // pt bin
  //                           0 , 1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 ,10 ,11 ,12 ,13 ,14 ,15
  Float_t pt_low[pt_total] = {0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4};
  Float_t pt_up[pt_total]  = {0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,4.2};

  // x and y range
  Float_t x_low[pt_total] = {-0.4, -0.6,-0.6, -0.6,-0.6,-0.6,-0.6,-0.8,-0.8,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0};
  Float_t x_up[pt_total]  = { 1.4,  1.4, 1.4,  1.4, 1.4, 1.4, 1.4, 1.6, 1.6, 2.0, 2.0, 2.0, 2.4, 2.4, 2.4, 2.4};

  Float_t y_low[pt_total] = {-0.15,-0.2,-0.3,-0.45,-0.6,-0.8,-1.0,-1.3,-1.5,-1.7,-1.8,-1.8,-2.0,-2.0,-2.0,-2.0};
  Float_t y_up[pt_total]  = { 0.15, 0.2, 0.2, 0.20, 0.3, 0.4, 0.6, 1.0, 1.0, 1.0, 1.2, 1.2, 1.4, 1.4, 1.4, 1.4};

  Double_t parfit[pt_total][cent_total][charge_total][eta_total][19];
  Int_t fit_bin = 13;

  TString Energy[2] = {"200GeV","39GeV"};
  TString Centrality1[4] = {"0-80%","0-10%","10-40%","40-80%"};
  TString Centrality2[4] = {"0-70%","0-10%","10-40%","40-70%"};
  TString Charge[2] = {"pos","neg"};
  TString EtaGap[4] = {"#eta_{gap} = 0.05","#eta_{gap} = 0.10","#eta_{gap} = 0.20","#eta_{gap} = 0.50"};
  TString Order[2] = {"2nd","3rd"};
  TString Centrality[4] = {"0080","0010","1040","4080"};
  TString Title[pt_total][cent_total][charge_total][eta_total];

  const Int_t fill_color_pion = kAzure-2;
  const Int_t fill_style_pion = 3002;
  const Int_t fill_color_kaon = kGray+2;
  const Int_t fill_style_kaon = 3003;
  const Int_t fill_color_proton = 2;
  const Int_t fill_style_proton = 3004;

  TString inputname = Form("/project/projectdirs/star/xusun/OutPut/AuAu%s/Mass2_nSigmaPion/merged_file/merged_file_%s_M2_nSigPion_%s_etagap_00.root",Energy[mEnergy].Data(),Energy[mEnergy].Data(),Centrality[mCentrality].Data()); 

  for(Int_t pt_bin = 0; pt_bin < pt_total; pt_bin ++)
  {
    for(Int_t cent = cent_start[mCentrality]; cent < cent_stop[mCentrality]; cent++)
    {
      for(Int_t charge = charge_start; charge < charge_stop; charge++)
      {
	for(Int_t eta_bin = eta_start; eta_bin < eta_stop; eta_bin++)
	{
	  Title[pt_bin][cent][charge][eta_bin] = Form("%s, %s, %1.1f-%1.1f GeV/c, %s, %s, total",Energy[mEnergy].Data(),Centrality1[cent].Data(), pt_low[pt_bin], pt_up[pt_bin],Charge[charge].Data(),EtaGap[eta_bin].Data());
	  for(Int_t i = 0; i < 19; i++)
	  {
	    parfit[pt_bin][cent][charge][eta_bin][i] = 0.0;
	  }
	}
      }
    }
  }
  //-----------------------------------------------------------------------------------------

  //-----------------------------------------------------------------------------------------
  // read Histogram
  TFile *input_student;
  input_student = TFile::Open(inputname.Data());
  TH2F *h_XY[pt_total][cent_total][charge_total][eta_total];
  for(Int_t pt_bin = 0; pt_bin < pt_total; pt_bin ++)
  {
    for(Int_t cent = cent_start[mCentrality]; cent < cent_stop[mCentrality]; cent++)
    {
      for(Int_t charge = charge_start; charge < charge_stop; charge++)
      {
	for(Int_t eta_bin = eta_start; eta_bin < eta_stop; eta_bin++)
	{
	  TString HistName = Form("Spectra_pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_pik_2nd",pt_bin,cent,charge,eta_bin);
	  h_XY[pt_bin][cent][charge][eta_bin] = (TH2F*)input_student->Get(HistName.Data());
          if(pt_bin > 8)
	  {
	    h_XY[pt_bin][cent][charge][eta_bin]->RebinX(2);
            h_XY[pt_bin][cent][charge][eta_bin]->RebinY(2);
            if(pt_bin > 12)
	    {
	      h_XY[pt_bin][cent][charge][eta_bin]->RebinX(2);
              h_XY[pt_bin][cent][charge][eta_bin]->RebinY(2);
	    }
	  }
	}
      }
    }
  }

  TCanvas *c_XY[cent_total][charge_total][eta_total];
  for(Int_t cent = cent_start[mCentrality]; cent < cent_stop[mCentrality]; cent++)
  {
    for(Int_t charge = charge_start; charge < charge_stop; charge++)
    {
      for(Int_t eta_bin = eta_start; eta_bin < eta_stop; eta_bin++)
      {
	TString CanName = Form("c_XY_cent_%d_charge_%d_etagap_%d",cent,charge,eta_bin);
	c_XY[cent][charge][eta_bin] = new TCanvas(CanName.Data(),CanName.Data(),1400,10,1200,1200);
	c_XY[cent][charge][eta_bin]->Divide(4,4);
	for(Int_t pt_bin = 0; pt_bin < pt_total; pt_bin++)
	{
	  c_XY[cent][charge][eta_bin]->cd(pt_bin+1);
	  c_XY[cent][charge][eta_bin]->cd(pt_bin+1)->SetLeftMargin(0.15);
	  c_XY[cent][charge][eta_bin]->cd(pt_bin+1)->SetRightMargin(0.15);
	  c_XY[cent][charge][eta_bin]->cd(pt_bin+1)->SetTopMargin(0.15);
	  c_XY[cent][charge][eta_bin]->cd(pt_bin+1)->SetBottomMargin(0.15);
	  c_XY[cent][charge][eta_bin]->cd(pt_bin+1)->SetLogz();
	  c_XY[cent][charge][eta_bin]->cd(pt_bin+1)->SetTicks(1,1);
	  h_XY[pt_bin][cent][charge][eta_bin]->SetStats(0);
	  h_XY[pt_bin][cent][charge][eta_bin]->SetTitle(Title[pt_bin][cent][charge][eta_bin].Data());
	  h_XY[pt_bin][cent][charge][eta_bin]->GetXaxis()->SetTitle("x(n#sigma_{#pi},m^{2})");
	  h_XY[pt_bin][cent][charge][eta_bin]->GetXaxis()->CenterTitle();
	  h_XY[pt_bin][cent][charge][eta_bin]->GetYaxis()->SetTitle("y(n#sigma_{#pi},m^{2})");
	  h_XY[pt_bin][cent][charge][eta_bin]->GetYaxis()->CenterTitle();
	  h_XY[pt_bin][cent][charge][eta_bin]->GetYaxis()->SetTitleOffset(1.6);
	  h_XY[pt_bin][cent][charge][eta_bin]->Draw("colz");
	}
      }
    }
  }
  //-----------------------------------------------------------------------------------------

  //-----------------------------------------------------------------------------------------
  // student_t function
  TF2 *f_student[pt_total][cent_total][charge_total][eta_total];
  for(Int_t pt_bin = 5; pt_bin < pt_total; pt_bin ++)
  {
    for(Int_t cent = cent_start[mCentrality]; cent < cent_stop[mCentrality]; cent++)
    {
      for(Int_t charge = charge_start; charge < charge_stop; charge++)
      {
	for(Int_t eta_bin = eta_start; eta_bin < eta_stop; eta_bin++)
	{
	    TString FuncName = Form("f_pt_%d_cent_%d_charge_%d_etagap_%d",pt_bin,cent,charge,eta_bin);
	    f_student[pt_bin][cent][charge][eta_bin] = new TF2(FuncName.Data(),student_t_2d_fit,x_low[pt_bin],x_up[pt_bin],y_low[pt_bin],y_up[pt_bin],19);
	    f_student[pt_bin][cent][charge][eta_bin]->SetNpx(400);
	    f_student[pt_bin][cent][charge][eta_bin]->SetNpy(400);

	    f_student[pt_bin][cent][charge][eta_bin]->SetParLimits(0,0,nu_pi[pt_bin]*10.0);
	    f_student[pt_bin][cent][charge][eta_bin]->SetParLimits(1,0,nu_pi[pt_bin]*10.0);
	    f_student[pt_bin][cent][charge][eta_bin]->SetParLimits(2,x_pi[pt_bin]-0.08,x_pi[pt_bin]+0.08);
	    f_student[pt_bin][cent][charge][eta_bin]->SetParLimits(3,y_pi[pt_bin]-0.08,y_pi[pt_bin]+0.08);
	    f_student[pt_bin][cent][charge][eta_bin]->SetParLimits(4,width_pi[pt_bin]*0.4,width_pi[pt_bin]*1.6);
	    f_student[pt_bin][cent][charge][eta_bin]->SetParLimits(5,width_pi[pt_bin]*0.4,width_pi[pt_bin]*1.6);
	    f_student[pt_bin][cent][charge][eta_bin]->SetParLimits(7,x_k[pt_bin]-0.18,x_k[pt_bin]+0.18);
	    f_student[pt_bin][cent][charge][eta_bin]->SetParLimits(8,y_k[pt_bin]-0.18,y_k[pt_bin]+0.18);
	    f_student[pt_bin][cent][charge][eta_bin]->SetParLimits(9,width_k[pt_bin]*0.4,width_k[pt_bin]*1.6);
	    f_student[pt_bin][cent][charge][eta_bin]->SetParLimits(10,width_k[pt_bin]*0.4,width_k[pt_bin]*1.6);
	    f_student[pt_bin][cent][charge][eta_bin]->SetParLimits(12,0,nu_p[pt_bin]*10.0);
	    //	f_student_total[pt_bin][cent][charge][eta_bin][flow_bin]->SetParLimits(13,0,nu_p[pt_bin]*10.0);
	    f_student[pt_bin][cent][charge][eta_bin]->SetParLimits(14,x_p[pt_bin]-0.08,x_p[pt_bin]+0.08);
	    f_student[pt_bin][cent][charge][eta_bin]->SetParLimits(15,y_p[pt_bin]-0.08,y_p[pt_bin]+0.08);
	    f_student[pt_bin][cent][charge][eta_bin]->SetParLimits(16,width_p_x[pt_bin]*0.4,width_p_x[pt_bin]*1.6);
	    f_student[pt_bin][cent][charge][eta_bin]->SetParLimits(17,width_p_y[pt_bin]*0.4,width_p_y[pt_bin]*1.6);

	    // parameter for pion start
	    f_student[pt_bin][cent][charge][eta_bin]->SetParameter(0,nu_pi[pt_bin]);
	    f_student[pt_bin][cent][charge][eta_bin]->SetParameter(1,nu_pi[pt_bin]);
	    f_student[pt_bin][cent][charge][eta_bin]->SetParameter(2,x_pi[pt_bin]);
	    f_student[pt_bin][cent][charge][eta_bin]->SetParameter(3,y_pi[pt_bin]);
	    f_student[pt_bin][cent][charge][eta_bin]->SetParameter(4,width_pi[pt_bin]);
	    f_student[pt_bin][cent][charge][eta_bin]->SetParameter(5,width_pi[pt_bin]);
	    f_student[pt_bin][cent][charge][eta_bin]->SetParameter(6,1);
	    // parameter for pion stop

	    // parameter for kaon start
	    f_student[pt_bin][cent][charge][eta_bin]->SetParameter(7,x_k[pt_bin]);
	    f_student[pt_bin][cent][charge][eta_bin]->SetParameter(8,y_k[pt_bin]);
	    f_student[pt_bin][cent][charge][eta_bin]->SetParameter(9,width_k[pt_bin]);
	    f_student[pt_bin][cent][charge][eta_bin]->SetParameter(10,width_k[pt_bin]);
	    f_student[pt_bin][cent][charge][eta_bin]->SetParameter(11,1);
	    // parameter for kaon stop

	    // parameter for proton start
	    f_student[pt_bin][cent][charge][eta_bin]->SetParameter(12,nu_p[pt_bin]);
	    f_student[pt_bin][cent][charge][eta_bin]->SetParameter(13,nu_p[pt_bin]);
	    f_student[pt_bin][cent][charge][eta_bin]->SetParameter(14,x_p[pt_bin]);
	    f_student[pt_bin][cent][charge][eta_bin]->SetParameter(15,y_p[pt_bin]);
	    f_student[pt_bin][cent][charge][eta_bin]->SetParameter(16,width_p_x[pt_bin]);
	    f_student[pt_bin][cent][charge][eta_bin]->SetParameter(17,width_p_y[pt_bin]);
	    f_student[pt_bin][cent][charge][eta_bin]->SetParameter(18,1);
	    // parameter for proton stop

	    Float_t pion_max = h_XY[pt_bin][cent][charge][eta_bin]->GetBinContent(h_XY[pt_bin][cent][charge][eta_bin]->GetXaxis()->FindBin(x_pi[pt_bin]),h_XY[pt_bin][cent][charge][eta_bin]->GetYaxis()->FindBin(y_pi[pt_bin]));
	    Float_t pion_fit = f_student[pt_bin][cent][charge][eta_bin]->Eval(x_pi[pt_bin],y_pi[pt_bin]);
	    Float_t pion_Nrom = pion_max/pion_fit;

	    Float_t kaon_max = h_XY[pt_bin][cent][charge][eta_bin]->GetBinContent(h_XY[pt_bin][cent][charge][eta_bin]->GetXaxis()->FindBin(x_k[pt_bin]),h_XY[pt_bin][cent][charge][eta_bin]->GetYaxis()->FindBin(y_k[pt_bin]));
	    Float_t kaon_fit = f_student[pt_bin][cent][charge][eta_bin]->Eval(x_k[pt_bin],y_k[pt_bin]);
	    Float_t kaon_Nrom = kaon_max/kaon_fit;

	    Float_t proton_max = h_XY[pt_bin][cent][charge][eta_bin]->GetBinContent(h_XY[pt_bin][cent][charge][eta_bin]->GetXaxis()->FindBin(x_p[pt_bin]),h_XY[pt_bin][cent][charge][eta_bin]->GetYaxis()->FindBin(y_p[pt_bin]));
	    Float_t proton_fit = f_student[pt_bin][cent][charge][eta_bin]->Eval(x_p[pt_bin],y_p[pt_bin]);
	    Float_t proton_Nrom = proton_max/proton_fit;
	    f_student[pt_bin][cent][charge][eta_bin]->SetParameter(6,pion_Nrom);
	    f_student[pt_bin][cent][charge][eta_bin]->SetParameter(11,kaon_Nrom);
	    f_student[pt_bin][cent][charge][eta_bin]->SetParameter(18,proton_Nrom);

	    f_student[pt_bin][cent][charge][eta_bin]->SetRange(x_pi[pt_bin]-(order_pion[pt_bin]+0.5)*width_pi[pt_bin],y_p[pt_bin],x_p[pt_bin],y_pi[pt_bin]+(order_pion[pt_bin]+0.5)*width_pi[pt_bin]);

	    f_student[pt_bin][cent][charge][eta_bin]->SetLineStyle(1);
	    f_student[pt_bin][cent][charge][eta_bin]->SetLineColor(2);
	}
      }
    }
  }

  // first fit
  for(Int_t cent = cent_start[mCentrality]; cent < cent_stop[mCentrality]; cent++)
  {
    for(Int_t charge = charge_start; charge < charge_stop; charge++)
    {
      for(Int_t eta_bin = eta_start; eta_bin < eta_stop; eta_bin++)
      {
	  for(Int_t pt_bin = 5; pt_bin < pt_total; pt_bin++)
	  {
	    if(mode == 0)
//	    if(mode == 0 && pt_bin == fit_bin)
	    {
	      h_XY[pt_bin][cent][charge][eta_bin]->Fit(f_student[pt_bin][cent][charge][eta_bin],"MQRN");
	    }
	    c_XY[cent][charge][eta_bin]->cd(pt_bin+1);
	    f_student[pt_bin][cent][charge][eta_bin]->SetContour(15);
	    f_student[pt_bin][cent][charge][eta_bin]->DrawCopy("cont3 same");
	    for(Int_t i = 0; i < 19; i++)
	    {
	      parfit[pt_bin][cent][charge][eta_bin][i] = f_student[pt_bin][cent][charge][eta_bin]->GetParameter(i);
	    }
	  }
      }
    }
  }
  TH1F *h_parfit_1st[pt_total][cent_total][charge_total][eta_total];
  for(Int_t pt_bin = 0; pt_bin < pt_total; pt_bin ++)
  {
    for(Int_t cent = cent_start[mCentrality]; cent < cent_stop[mCentrality]; cent++)
    {
      for(Int_t charge = charge_start; charge < charge_stop; charge++)
      {
	for(Int_t eta_bin = eta_start; eta_bin < eta_stop; eta_bin++)
	{
	  TString ParName = Form("parfit_pt_%d_Centrality_%d_charge_%d_etagap_%d_1st",pt_bin,cent,charge,eta_bin);
	  h_parfit_1st[pt_bin][cent][charge][eta_bin] = new TH1F(ParName.Data(),ParName.Data(),19,-0.5,18.5);
	  for(Int_t i = 0; i < 19; i++)
	  {
	    h_parfit_1st[pt_bin][cent][charge][eta_bin]->SetBinContent(i+1,parfit[pt_bin][cent][charge][eta_bin][i]);
	  }
	}
      }
    }
  }

  // second fit
  for(Int_t pt_bin = 5; pt_bin < pt_total; pt_bin++)
  {
    for(Int_t cent = cent_start[mCentrality]; cent < cent_stop[mCentrality]; cent++)
    {
      for(Int_t charge = charge_start; charge < charge_stop; charge++)
      {
	for(Int_t eta_bin = eta_start; eta_bin < eta_stop; eta_bin++)
	{
	  for(Int_t i = 0; i < 19; i++)
	  {
	    f_student[pt_bin][cent][charge][eta_bin]->ReleaseParameter(i);
	    f_student[pt_bin][cent][charge][eta_bin]->SetParameter(i,0.0);
	    f_student[pt_bin][cent][charge][eta_bin]->SetParError(i,0.0);
	    f_student[pt_bin][cent][charge][eta_bin]->SetParameter(i,parfit[pt_bin][cent][charge][eta_bin][i]);
	  }

	  f_student[pt_bin][cent][charge][eta_bin]->SetParLimits(0,0,342);
	  f_student[pt_bin][cent][charge][eta_bin]->SetParLimits(1,0,342);
	  f_student[pt_bin][cent][charge][eta_bin]->SetParLimits(2,x_pi[pt_bin]-0.08,x_pi[pt_bin]+0.08);
	  f_student[pt_bin][cent][charge][eta_bin]->SetParLimits(3,y_pi[pt_bin]-0.08,y_pi[pt_bin]+0.08);
	  f_student[pt_bin][cent][charge][eta_bin]->SetParLimits(7,x_k[pt_bin]-0.18,x_k[pt_bin]+0.18);
	  f_student[pt_bin][cent][charge][eta_bin]->SetParLimits(8,y_k[pt_bin]-0.18,y_k[pt_bin]+0.18);
	  if(pt_bin == 12) f_student[pt_bin][cent][charge][eta_bin]->SetParLimits(8,y_k[pt_bin]-0.36,y_k[pt_bin]+0.36);
	  f_student[pt_bin][cent][charge][eta_bin]->SetParLimits(12,0,342);
	  f_student[pt_bin][cent][charge][eta_bin]->SetParLimits(13,0,342);
	  f_student[pt_bin][cent][charge][eta_bin]->SetParLimits(14,x_p[pt_bin]-0.08,x_p[pt_bin]+0.08);
	  if(pt_bin >= 14) f_student[pt_bin][cent][charge][eta_bin]->SetParLimits(14,x_p[pt_bin]-0.16,x_p[pt_bin]+0.16);
	  f_student[pt_bin][cent][charge][eta_bin]->SetParLimits(15,y_p[pt_bin]-0.08,y_p[pt_bin]+0.08);
	  if(pt_bin >= 14) f_student[pt_bin][cent][charge][eta_bin]->SetParLimits(15,y_p[pt_bin]-0.16,y_p[pt_bin]+0.16);
	  f_student[pt_bin][cent][charge][eta_bin]->SetRange(x_pi[pt_bin]-order_pion[pt_bin]*parfit[pt_bin][cent][charge][eta_bin][4],y_p[pt_bin]-order_proton[pt_bin]*parfit[pt_bin][cent][charge][eta_bin][17],x_p[pt_bin]+order_proton[pt_bin]*parfit[pt_bin][cent][charge][eta_bin][16],y_pi[pt_bin]+order_pion[pt_bin]*parfit[pt_bin][cent][charge][eta_bin][5]);

	  f_student[pt_bin][cent][charge][eta_bin]->SetLineStyle(2);
	  f_student[pt_bin][cent][charge][eta_bin]->SetLineColor(4);
	}
      }
    }
  }
  for(Int_t cent = cent_start[mCentrality]; cent < cent_stop[mCentrality]; cent++)
  {
    for(Int_t charge = charge_start; charge < charge_stop; charge++)
    {
      for(Int_t eta_bin = eta_start; eta_bin < eta_stop; eta_bin++)
      {
	for(Int_t pt_bin = 5; pt_bin < pt_total; pt_bin++)
	{
	  if(mode == 0)
//	  if(mode == 0 && pt_bin == fit_bin)
	  {
	    cout << "pt_bin = " << pt_bin << ", charge = " << charge << ", eta_bin = " << eta_bin << endl;
	    h_XY[pt_bin][cent][charge][eta_bin]->Fit(f_student[pt_bin][cent][charge][eta_bin],"MRN");
	    Float_t chi2 = f_student[pt_bin][cent][charge][eta_bin]->GetChisquare();
	    Float_t Ndf  = f_student[pt_bin][cent][charge][eta_bin]->GetNDF();
	    cout << "chi2/Ndf = " << chi2 << "/" << Ndf << " = " << chi2/Ndf << endl;
	  }
	  c_XY[cent][charge][eta_bin]->cd(pt_bin+1);
	  f_student[pt_bin][cent][charge][eta_bin]->DrawCopy("cont3 same");
	  for(Int_t i = 0; i < 19; i++)
	  {
	    parfit[pt_bin][cent][charge][eta_bin][i] = f_student[pt_bin][cent][charge][eta_bin]->GetParameter(i);
	  }
	  f_student[pt_bin][cent][charge][eta_bin]->SetRange(x_low[pt_bin],y_low[pt_bin],x_up[pt_bin],y_up[pt_bin]);
	}
      }
    }
  }

  // save parfit to histogram
  TH1F *h_parfit_2nd[pt_total][cent_total][charge_total][eta_total];
  for(Int_t pt_bin = 0; pt_bin < pt_total; pt_bin ++)
  {
    for(Int_t cent = cent_start[mCentrality]; cent < cent_stop[mCentrality]; cent++)
    {
      for(Int_t charge = charge_start; charge < charge_stop; charge++)
      {
	for(Int_t eta_bin = eta_start; eta_bin < eta_stop; eta_bin++)
	{
	  TString ParName = Form("parfit_pt_%d_Centrality_%d_charge_%d_etagap_%d_2nd",pt_bin,cent,charge,eta_bin);
	  h_parfit_2nd[pt_bin][cent][charge][eta_bin] = new TH1F(ParName.Data(),ParName.Data(),19,-0.5,18.5);
	  for(Int_t i = 0; i < 19; i++)
	  {
	    h_parfit_2nd[pt_bin][cent][charge][eta_bin]->SetBinContent(i+1,parfit[pt_bin][cent][charge][eta_bin][i]);
	  }
	}
      }
    }
  }

  //-----------------------------------------------------------------------------------------

  //-----------------------------------------------------------------------------------------
  // Projection to X and Y

  TF2 *f_pion[pt_total][cent_total][charge_total][eta_total];
  TF2 *f_kaon[pt_total][cent_total][charge_total][eta_total];
  TF2 *f_proton[pt_total][cent_total][charge_total][eta_total];

  for(Int_t pt_bin = 5; pt_bin < pt_total; pt_bin++)
  {
    for(Int_t cent = cent_start[mCentrality]; cent < cent_stop[mCentrality]; cent++)
    {
      for(Int_t charge = charge_start; charge < charge_stop; charge++)
      {
	for(Int_t eta_bin = eta_start; eta_bin < eta_stop; eta_bin++)
	{
	  TString PionName = Form("f_pion_pt_%d_cent_%d_charge_%d_etagap_%d",pt_bin,cent,charge,eta_bin);
	  f_pion[pt_bin][cent][charge][eta_bin] = new TF2(PionName.Data(),student_t_2d_single,x_low[pt_bin],x_up[pt_bin],y_low[pt_bin],y_up[pt_bin],7);
	  f_pion[pt_bin][cent][charge][eta_bin]->SetNpx(400);
	  f_pion[pt_bin][cent][charge][eta_bin]->SetNpy(400);
	  f_pion[pt_bin][cent][charge][eta_bin]->SetContour(15);

	  TString KaonName = Form("f_kaon_pt_%d_cent_%d_charge_%d_etagap_%d",pt_bin,cent,charge,eta_bin);
	  f_kaon[pt_bin][cent][charge][eta_bin] = new TF2(KaonName.Data(),student_t_2d_single,x_low[pt_bin],x_up[pt_bin],y_low[pt_bin],y_up[pt_bin],7);
	  f_kaon[pt_bin][cent][charge][eta_bin]->SetNpx(400);
	  f_kaon[pt_bin][cent][charge][eta_bin]->SetNpy(400);
	  f_kaon[pt_bin][cent][charge][eta_bin]->SetContour(15);

	  TString ProtonName = Form("f_proton_pt_%d_cent_%d_charge_%d_etagap_%d",pt_bin,cent,charge,eta_bin);
	  f_proton[pt_bin][cent][charge][eta_bin] = new TF2(ProtonName.Data(),student_t_2d_single,x_low[pt_bin],x_up[pt_bin],y_low[pt_bin],y_up[pt_bin],7);
	  f_proton[pt_bin][cent][charge][eta_bin]->SetNpx(400);
	  f_proton[pt_bin][cent][charge][eta_bin]->SetNpy(400);
	  f_proton[pt_bin][cent][charge][eta_bin]->SetContour(15);

	  for(Int_t i = 0; i < 7; i++)
	  {
	    f_pion[pt_bin][cent][charge][eta_bin]->FixParameter(i,parfit[pt_bin][cent][charge][eta_bin][i]);
	    f_proton[pt_bin][cent][charge][eta_bin]->FixParameter(i,parfit[pt_bin][cent][charge][eta_bin][i+12]);
	  }
	  f_kaon[pt_bin][cent][charge][eta_bin]->FixParameter(0,(2.0*parfit[pt_bin][cent][charge][eta_bin][0]+parfit[pt_bin][cent][charge][eta_bin][12])/3.0);
	  f_kaon[pt_bin][cent][charge][eta_bin]->FixParameter(1,(2.0*parfit[pt_bin][cent][charge][eta_bin][1]+parfit[pt_bin][cent][charge][eta_bin][13])/3.0);
	  for(Int_t i = 2; i < 7; i++)
	  {
	    f_kaon[pt_bin][cent][charge][eta_bin]->FixParameter(i,parfit[pt_bin][cent][charge][eta_bin][i+5]);
	  }
	}
      }
    }
  }

  // Canvas
  TCanvas *c_X[cent_total][charge_total][eta_total];
  TCanvas *c_Y[cent_total][charge_total][eta_total];
  for(Int_t cent = cent_start[mCentrality]; cent < cent_stop[mCentrality]; cent++)
  {
    for(Int_t charge = charge_start; charge < charge_stop; charge++)
    {
      for(Int_t eta_bin = eta_start; eta_bin < eta_stop; eta_bin++)
      {
	TString CanName = Form("c_X_cent_%d_charge_%d_etagap_%d",cent,charge,eta_bin);
	c_X[cent][charge][eta_bin] = new TCanvas(CanName.Data(),CanName.Data(),1400,10,1200,1200);
	c_X[cent][charge][eta_bin]->Divide(4,4);
	CanName = Form("c_Y_cent_%d_charge_%d_etagap_%d",cent,charge,eta_bin);
	c_Y[cent][charge][eta_bin] = new TCanvas(CanName.Data(),CanName.Data(),1400,10,1200,1200);
	c_Y[cent][charge][eta_bin]->Divide(4,4);
	for(Int_t pt_bin = 0; pt_bin < pt_total; pt_bin++)
	{
	  c_X[cent][charge][eta_bin]->cd(pt_bin+1);
	  c_X[cent][charge][eta_bin]->cd(pt_bin+1)->SetLeftMargin(0.15);
	  c_X[cent][charge][eta_bin]->cd(pt_bin+1)->SetRightMargin(0.15);
	  c_X[cent][charge][eta_bin]->cd(pt_bin+1)->SetTopMargin(0.15);
	  c_X[cent][charge][eta_bin]->cd(pt_bin+1)->SetBottomMargin(0.15);
	  c_X[cent][charge][eta_bin]->cd(pt_bin+1)->SetLogz();
	  c_X[cent][charge][eta_bin]->cd(pt_bin+1)->SetTicks(1,1);

	  c_Y[cent][charge][eta_bin]->cd(pt_bin+1);
	  c_Y[cent][charge][eta_bin]->cd(pt_bin+1)->SetLeftMargin(0.15);
	  c_Y[cent][charge][eta_bin]->cd(pt_bin+1)->SetRightMargin(0.15);
	  c_Y[cent][charge][eta_bin]->cd(pt_bin+1)->SetTopMargin(0.15);
	  c_Y[cent][charge][eta_bin]->cd(pt_bin+1)->SetBottomMargin(0.15);
	  c_Y[cent][charge][eta_bin]->cd(pt_bin+1)->SetLogz();
	  c_Y[cent][charge][eta_bin]->cd(pt_bin+1)->SetTicks(1,1);
	}
      }
    }
  }

  // Projection X
  TH1F *h_X[pt_total][cent_total][charge_total][eta_total];
  TH1F *h_X_Fit[pt_total][cent_total][charge_total][eta_total];
  TH1F *h_X_Pion[pt_total][cent_total][charge_total][eta_total];
  TH1F *h_X_Kaon[pt_total][cent_total][charge_total][eta_total];
  TH1F *h_X_Proton[pt_total][cent_total][charge_total][eta_total];
  for(Int_t pt_bin = 0; pt_bin < pt_total; pt_bin++)
  {
    for(Int_t cent = cent_start[mCentrality]; cent < cent_stop[mCentrality]; cent++)
    {
      for(Int_t charge = charge_start; charge < charge_stop; charge++)
      {
	for(Int_t eta_bin = eta_start; eta_bin < eta_stop; eta_bin++)
	{
	  if(pt_bin < 5)
	  {
	    TString HistName = Form("pt_%d_Centrality_%d_charge_%d_etagap_%d_X",pt_bin,cent,charge,eta_bin);
	    h_X[pt_bin][cent][charge][eta_bin] = (TH1F*)h_XY[pt_bin][cent][charge][eta_bin]->ProjectionX(HistName.Data(),h_XY[pt_bin][cent][charge][eta_bin]->GetYaxis()->FindBin(y_low[pt_bin]),h_XY[pt_bin][cent][charge][eta_bin]->GetYaxis()->FindBin(y_up[pt_bin]));
            c_X[cent][charge][eta_bin]->cd(pt_bin+1);
	    h_X[pt_bin][cent][charge][eta_bin]->SetStats(0);
	    h_X[pt_bin][cent][charge][eta_bin]->DrawCopy("PE");
	  }
	  if(pt_bin >= 5)
	  {
	    TString HistName = Form("pt_%d_Centrality_%d_charge_%d_etagap_%d_X",pt_bin,cent,charge,eta_bin);
	    h_X[pt_bin][cent][charge][eta_bin] = (TH1F*)h_XY[pt_bin][cent][charge][eta_bin]->ProjectionX(HistName.Data(),h_XY[pt_bin][cent][charge][eta_bin]->GetYaxis()->FindBin(y_p[pt_bin]-order_proton[pt_bin]*parfit[pt_bin][cent][charge][eta_bin][17]),h_XY[pt_bin][cent][charge][eta_bin]->GetYaxis()->FindBin(y_pi[pt_bin]+order_pion[pt_bin]*parfit[pt_bin][cent][charge][eta_bin][5]));
            c_X[cent][charge][eta_bin]->cd(pt_bin+1);
	    h_X[pt_bin][cent][charge][eta_bin]->SetStats(0);
	    h_X[pt_bin][cent][charge][eta_bin]->DrawCopy("PE");

	    Int_t start_bin_X = ((TH2F*)f_student[pt_bin][cent][charge][eta_bin]->GetHistogram())->GetYaxis()->FindBin(y_p[pt_bin]-order_proton[pt_bin]*parfit[pt_bin][cent][charge][eta_bin][17]);
	    Int_t stop_bin_X  = ((TH2F*)f_student[pt_bin][cent][charge][eta_bin]->GetHistogram())->GetYaxis()->FindBin(y_pi[pt_bin]+order_pion[pt_bin]*parfit[pt_bin][cent][charge][eta_bin][5]);
	    Double_t Inte_X = h_X[pt_bin][cent][charge][eta_bin]->Integral(h_X[pt_bin][cent][charge][eta_bin]->FindBin(x_pi[pt_bin]-order_pion[pt_bin]*parfit[pt_bin][cent][charge][eta_bin][4]),h_X[pt_bin][cent][charge][eta_bin]->FindBin(x_p[pt_bin]+order_proton[pt_bin]*parfit[pt_bin][cent][charge][eta_bin][16]));

	    HistName = Form("pt_%d_Centrality_%d_charge_%d_etagap_%d_Fit_X",pt_bin,cent,charge,eta_bin);
            h_X_Fit[pt_bin][cent][charge][eta_bin] = (TH1F*)((TH2F*)f_student[pt_bin][cent][charge][eta_bin]->GetHistogram())->ProjectionX(HistName.Data(),start_bin_X,stop_bin_X);
         
            Double_t Inte_X_Fit = h_X_Fit[pt_bin][cent][charge][eta_bin]->Integral(h_X_Fit[pt_bin][cent][charge][eta_bin]->FindBin(x_pi[pt_bin]-order_pion[pt_bin]*parfit[pt_bin][cent][charge][eta_bin][4]),h_X_Fit[pt_bin][cent][charge][eta_bin]->FindBin(x_p[pt_bin]+order_proton[pt_bin]*parfit[pt_bin][cent][charge][eta_bin][16]));
	    Double_t scale_X = 1.0;
	    if(Inte_X_Fit > 0.0) scale_X = Inte_X/Inte_X_Fit;
	    h_X_Fit[pt_bin][cent][charge][eta_bin]->Scale(scale_X*h_X[pt_bin][cent][charge][eta_bin]->GetBinWidth(1)/h_X_Fit[pt_bin][cent][charge][eta_bin]->GetBinWidth(1));
            h_X_Fit[pt_bin][cent][charge][eta_bin]->SetLineColor(2);
	    h_X_Fit[pt_bin][cent][charge][eta_bin]->SetLineStyle(2);
	    h_X_Fit[pt_bin][cent][charge][eta_bin]->SetFillColor(0);

	    HistName = Form("pt_%d_Centrality_%d_charge_%d_etagap_%d_Pion_X",pt_bin,cent,charge,eta_bin);
	    h_X_Pion[pt_bin][cent][charge][eta_bin] = (TH1F*)((TH2F*)f_pion[pt_bin][cent][charge][eta_bin]->GetHistogram())->ProjectionX(HistName.Data(),start_bin_X,stop_bin_X);
	    h_X_Pion[pt_bin][cent][charge][eta_bin]->Scale(scale_X*h_X[pt_bin][cent][charge][eta_bin]->GetBinWidth(1)/h_X_Pion[pt_bin][cent][charge][eta_bin]->GetBinWidth(1));
	    h_X_Pion[pt_bin][cent][charge][eta_bin]->SetFillColor(fill_color_pion);
	    h_X_Pion[pt_bin][cent][charge][eta_bin]->SetLineColor(fill_color_pion);
	    h_X_Pion[pt_bin][cent][charge][eta_bin]->SetFillStyle(fill_style_pion);

	    HistName = Form("pt_%d_Centrality_%d_charge_%d_etagap_%d_Kaon_X",pt_bin,cent,charge,eta_bin);
	    h_X_Kaon[pt_bin][cent][charge][eta_bin] = (TH1F*)((TH2F*)f_kaon[pt_bin][cent][charge][eta_bin]->GetHistogram())->ProjectionX(HistName.Data(),start_bin_X,stop_bin_X);
	    h_X_Kaon[pt_bin][cent][charge][eta_bin]->Scale(scale_X*h_X[pt_bin][cent][charge][eta_bin]->GetBinWidth(1)/h_X_Kaon[pt_bin][cent][charge][eta_bin]->GetBinWidth(1));
	    h_X_Kaon[pt_bin][cent][charge][eta_bin]->SetFillColor(fill_color_kaon);
	    h_X_Kaon[pt_bin][cent][charge][eta_bin]->SetLineColor(fill_color_kaon);
	    h_X_Kaon[pt_bin][cent][charge][eta_bin]->SetFillStyle(fill_style_kaon);

	    HistName = Form("pt_%d_Centrality_%d_charge_%d_etagap_%d_Proton_X",pt_bin,cent,charge,eta_bin);
	    h_X_Proton[pt_bin][cent][charge][eta_bin] = (TH1F*)((TH2F*)f_proton[pt_bin][cent][charge][eta_bin]->GetHistogram())->ProjectionX(HistName.Data(),start_bin_X,stop_bin_X);
	    h_X_Proton[pt_bin][cent][charge][eta_bin]->Scale(scale_X*h_X[pt_bin][cent][charge][eta_bin]->GetBinWidth(1)/h_X_Proton[pt_bin][cent][charge][eta_bin]->GetBinWidth(1));
	    h_X_Proton[pt_bin][cent][charge][eta_bin]->SetFillColor(fill_color_proton);
	    h_X_Proton[pt_bin][cent][charge][eta_bin]->SetLineColor(fill_color_proton);
	    h_X_Proton[pt_bin][cent][charge][eta_bin]->SetFillStyle(fill_style_proton);

	    h_X_Pion[pt_bin][cent][charge][eta_bin]->DrawCopy("hE same");
	    h_X_Kaon[pt_bin][cent][charge][eta_bin]->DrawCopy("hE same");
	    h_X_Proton[pt_bin][cent][charge][eta_bin]->DrawCopy("hE same");
	    h_X_Fit[pt_bin][cent][charge][eta_bin]->DrawCopy("l same");

	    PlotLine(x_pi[pt_bin]-order_pion[pt_bin]*width_pi[pt_bin],x_pi[pt_bin]-order_pion[pt_bin]*width_pi[pt_bin],0,h_X[pt_bin][cent][charge][eta_bin]->GetMaximum()/1.5,1,2,2);
	    PlotLine(x_p[pt_bin]+order_proton[pt_bin]*width_p_x[pt_bin],x_p[pt_bin]+order_proton[pt_bin]*width_p_x[pt_bin],0,h_X[pt_bin][cent][charge][eta_bin]->GetMaximum()/1.5,1,2,2);
	  }
	}
      }
    }
  }

  // Projection Y
  TH1F *h_Y[pt_total][cent_total][charge_total][eta_total];
  TH1F *h_Y_Fit[pt_total][cent_total][charge_total][eta_total];
  TH1F *h_Y_Pion[pt_total][cent_total][charge_total][eta_total];
  TH1F *h_Y_Kaon[pt_total][cent_total][charge_total][eta_total];
  TH1F *h_Y_Proton[pt_total][cent_total][charge_total][eta_total];
  for(Int_t pt_bin = 0; pt_bin < pt_total; pt_bin++)
  {
    for(Int_t cent = cent_start[mCentrality]; cent < cent_stop[mCentrality]; cent++)
    {
      for(Int_t charge = charge_start; charge < charge_stop; charge++)
      {
	for(Int_t eta_bin = eta_start; eta_bin < eta_stop; eta_bin++)
	{
	  if(pt_bin < 5)
	  {
	    TString HistName = Form("pt_%d_Centrality_%d_charge_%d_etagap_%d_Y",pt_bin,cent,charge,eta_bin);
	    h_Y[pt_bin][cent][charge][eta_bin] = (TH1F*)h_XY[pt_bin][cent][charge][eta_bin]->ProjectionY(HistName.Data(),h_XY[pt_bin][cent][charge][eta_bin]->GetXaxis()->FindBin(x_low[pt_bin]),h_XY[pt_bin][cent][charge][eta_bin]->GetXaxis()->FindBin(x_up[pt_bin]));
	    c_Y[cent][charge][eta_bin]->cd(pt_bin+1);
	    h_Y[pt_bin][cent][charge][eta_bin]->SetStats(0);
	    h_Y[pt_bin][cent][charge][eta_bin]->DrawCopy("PE");
	  }
	  if(pt_bin >= 5)
	  {
	    TString HistName = Form("pt_%d_Centrality_%d_charge_%d_etagap_%d_Y",pt_bin,cent,charge,eta_bin);
	    h_Y[pt_bin][cent][charge][eta_bin] = (TH1F*)h_XY[pt_bin][cent][charge][eta_bin]->ProjectionY(HistName.Data(),h_XY[pt_bin][cent][charge][eta_bin]->GetXaxis()->FindBin(x_pi[pt_bin]-order_pion[pt_bin]*parfit[pt_bin][cent][charge][eta_bin][4]),h_XY[pt_bin][cent][charge][eta_bin]->GetXaxis()->FindBin(x_p[pt_bin]+order_proton[pt_bin]*parfit[pt_bin][cent][charge][eta_bin][16]));
	    c_Y[cent][charge][eta_bin]->cd(pt_bin+1);
	    h_Y[pt_bin][cent][charge][eta_bin]->SetStats(0);
	    h_Y[pt_bin][cent][charge][eta_bin]->DrawCopy("PE");

	    Int_t start_bin_Y = ((TH2F*)f_student[pt_bin][cent][charge][eta_bin]->GetHistogram())->GetXaxis()->FindBin(x_pi[pt_bin]-order_pion[pt_bin]*parfit[pt_bin][cent][charge][eta_bin][4]);
	    Int_t stop_bin_Y = ((TH2F*)f_student[pt_bin][cent][charge][eta_bin]->GetHistogram())->GetXaxis()->FindBin(x_p[pt_bin]+order_proton[pt_bin]*parfit[pt_bin][cent][charge][eta_bin][16]);
	    Double_t Inte_Y = h_Y[pt_bin][cent][charge][eta_bin]->Integral(h_Y[pt_bin][cent][charge][eta_bin]->FindBin(y_p[pt_bin]-order_proton[pt_bin]*parfit[pt_bin][cent][charge][eta_bin][17]),h_Y[pt_bin][cent][charge][eta_bin]->FindBin(y_pi[pt_bin]+order_pion[pt_bin]*parfit[pt_bin][cent][charge][eta_bin][5]));

	    HistName = Form("pt_%d_Centrality_%d_charge_%d_etagap_%d_Fit_Y",pt_bin,cent,charge,eta_bin);
	    h_Y_Fit[pt_bin][cent][charge][eta_bin] = (TH1F*)((TH2F*)f_student[pt_bin][cent][charge][eta_bin]->GetHistogram())->ProjectionY(HistName.Data(),start_bin_Y,stop_bin_Y);
	    Double_t Inte_Y_Fit = h_Y_Fit[pt_bin][cent][charge][eta_bin]->Integral(h_Y_Fit[pt_bin][cent][charge][eta_bin]->FindBin(y_p[pt_bin]-order_proton[pt_bin]*parfit[pt_bin][cent][charge][eta_bin][17]),h_Y_Fit[pt_bin][cent][charge][eta_bin]->FindBin(y_pi[pt_bin]+order_pion[pt_bin]*parfit[pt_bin][cent][charge][eta_bin][5]));
	    Double_t scale_Y = 1.0;
	    if(Inte_Y_Fit > 0.0) scale_Y = Inte_Y/Inte_Y_Fit;
	    h_Y_Fit[pt_bin][cent][charge][eta_bin]->Scale(scale_Y*h_Y[pt_bin][cent][charge][eta_bin]->GetBinWidth(1)/h_Y_Fit[pt_bin][cent][charge][eta_bin]->GetBinWidth(1));
	    h_Y_Fit[pt_bin][cent][charge][eta_bin]->GetXaxis()->SetTitle("y");
	    h_Y_Fit[pt_bin][cent][charge][eta_bin]->GetXaxis()->CenterTitle();
	    h_Y_Fit[pt_bin][cent][charge][eta_bin]->SetLineColor(2);
	    h_Y_Fit[pt_bin][cent][charge][eta_bin]->SetLineStyle(2);
	    h_Y_Fit[pt_bin][cent][charge][eta_bin]->SetFillColor(0);

	    HistName = Form("pt_%d_Centrality_%d_charge_%d_etagap_%d_Pion_Y",pt_bin,cent,charge,eta_bin);
	    h_Y_Pion[pt_bin][cent][charge][eta_bin] = (TH1F*)((TH2F*)f_pion[pt_bin][cent][charge][eta_bin]->GetHistogram())->ProjectionY(HistName.Data(),start_bin_Y,stop_bin_Y);
	    h_Y_Pion[pt_bin][cent][charge][eta_bin]->Scale(scale_Y*h_Y[pt_bin][cent][charge][eta_bin]->GetBinWidth(1)/h_Y_Pion[pt_bin][cent][charge][eta_bin]->GetBinWidth(1));
	    h_Y_Pion[pt_bin][cent][charge][eta_bin]->SetFillColor(fill_color_pion);
	    h_Y_Pion[pt_bin][cent][charge][eta_bin]->SetLineColor(fill_color_pion);
	    h_Y_Pion[pt_bin][cent][charge][eta_bin]->SetFillStyle(fill_style_pion);

	    HistName = Form("pt_%d_Centrality_%d_charge_%d_etagap_%d_Kaon_Y",pt_bin,cent,charge,eta_bin);
	    h_Y_Kaon[pt_bin][cent][charge][eta_bin] = (TH1F*)((TH2F*)f_kaon[pt_bin][cent][charge][eta_bin]->GetHistogram())->ProjectionY(HistName.Data(),start_bin_Y,stop_bin_Y);
	    h_Y_Kaon[pt_bin][cent][charge][eta_bin]->Scale(scale_Y*h_Y[pt_bin][cent][charge][eta_bin]->GetBinWidth(1)/h_Y_Kaon[pt_bin][cent][charge][eta_bin]->GetBinWidth(1));
	    h_Y_Kaon[pt_bin][cent][charge][eta_bin]->SetFillColor(fill_color_kaon);
	    h_Y_Kaon[pt_bin][cent][charge][eta_bin]->SetLineColor(fill_color_kaon);
	    h_Y_Kaon[pt_bin][cent][charge][eta_bin]->SetFillStyle(fill_style_kaon);

	    HistName = Form("pt_%d_Centrality_%d_charge_%d_etagap_%d_Proton_Y",pt_bin,cent,charge,eta_bin);
	    h_Y_Proton[pt_bin][cent][charge][eta_bin] = (TH1F*)((TH2F*)f_proton[pt_bin][cent][charge][eta_bin]->GetHistogram())->ProjectionY(HistName.Data(),start_bin_Y,stop_bin_Y);
	    h_Y_Proton[pt_bin][cent][charge][eta_bin]->Scale(scale_Y*h_Y[pt_bin][cent][charge][eta_bin]->GetBinWidth(1)/h_Y_Proton[pt_bin][cent][charge][eta_bin]->GetBinWidth(1));
	    h_Y_Proton[pt_bin][cent][charge][eta_bin]->SetFillColor(fill_color_proton);
	    h_Y_Proton[pt_bin][cent][charge][eta_bin]->SetLineColor(fill_color_proton);
	    h_Y_Proton[pt_bin][cent][charge][eta_bin]->SetFillStyle(fill_style_proton);

	    h_Y_Pion[pt_bin][cent][charge][eta_bin]->DrawCopy("hE same");
	    h_Y_Kaon[pt_bin][cent][charge][eta_bin]->DrawCopy("hE same");
	    h_Y_Proton[pt_bin][cent][charge][eta_bin]->DrawCopy("hE same");
	    h_Y_Fit[pt_bin][cent][charge][eta_bin]->DrawCopy("l same");

	    PlotLine(y_pi[pt_bin]+order_pion[pt_bin]*width_pi[pt_bin],y_pi[pt_bin]+order_pion[pt_bin]*width_pi[pt_bin],0,h_Y[pt_bin][cent][charge][eta_bin]->GetMaximum()/1.5,1,2,2);
	    PlotLine(y_p[pt_bin]-order_proton[pt_bin]*width_p_y[pt_bin],y_p[pt_bin]-order_proton[pt_bin]*width_p_y[pt_bin],0,h_Y[pt_bin][cent][charge][eta_bin]->GetMaximum()/1.5,1,2,2);
	  }
	}
      }
    }
  }
  //-----------------------------------------------------------------------------------------

  //-----------------------------------------------------------------------------------------
  // Subtraction proton
  TH2F *h_XY_sub[pt_total][cent_total][charge_total][eta_total];
  for(Int_t pt_bin = 0; pt_bin < pt_total; pt_bin ++)
  {
    for(Int_t cent = cent_start[mCentrality]; cent < cent_stop[mCentrality]; cent++)
    {
      for(Int_t charge = charge_start; charge < charge_stop; charge++)
      {
	for(Int_t eta_bin = eta_start; eta_bin < eta_stop; eta_bin++)
	{
	  TString HistName = Form("pt_%d_Centrality_%d_Charge_%d_etagap_%d_sub",pt_bin,cent,charge,eta_bin);
	  h_XY_sub[pt_bin][cent][charge][eta_bin] = (TH2F*)h_XY[pt_bin][cent][charge][eta_bin]->Clone(HistName.Data());
	}
      }
    }
  }
  TCanvas *c_XY_sub[cent_total][charge_total][eta_total];
  for(Int_t cent = cent_start[mCentrality]; cent < cent_stop[mCentrality]; cent++)
  {
    for(Int_t charge = charge_start; charge < charge_stop; charge++)
    {
      for(Int_t eta_bin = eta_start; eta_bin < eta_stop; eta_bin++)
      {
	TString CanName = Form("c_XY_cent_%d_charge_%d_etagap_%d_sub",cent,charge,eta_bin);
	c_XY_sub[cent][charge][eta_bin] = new TCanvas(CanName.Data(),CanName.Data(),1400,10,1200,1200);
	c_XY_sub[cent][charge][eta_bin]->Divide(4,4);
	for(Int_t pt_bin = 0; pt_bin < pt_total; pt_bin++)
	{
	  c_XY_sub[cent][charge][eta_bin]->cd(pt_bin+1);
	  c_XY_sub[cent][charge][eta_bin]->cd(pt_bin+1)->SetLeftMargin(0.15);
	  c_XY_sub[cent][charge][eta_bin]->cd(pt_bin+1)->SetRightMargin(0.15);
	  c_XY_sub[cent][charge][eta_bin]->cd(pt_bin+1)->SetTopMargin(0.15);
	  c_XY_sub[cent][charge][eta_bin]->cd(pt_bin+1)->SetBottomMargin(0.15);
	  c_XY_sub[cent][charge][eta_bin]->cd(pt_bin+1)->SetLogz();
	  c_XY_sub[cent][charge][eta_bin]->cd(pt_bin+1)->SetTicks(1,1);
	  h_XY_sub[pt_bin][cent][charge][eta_bin]->SetStats(0);
	  h_XY_sub[pt_bin][cent][charge][eta_bin]->SetTitle(Title[pt_bin][cent][charge][eta_bin].Data());
	  h_XY_sub[pt_bin][cent][charge][eta_bin]->GetXaxis()->SetTitle("x(n#sigma_{#pi},m^{2})");
	  h_XY_sub[pt_bin][cent][charge][eta_bin]->GetXaxis()->CenterTitle();
	  h_XY_sub[pt_bin][cent][charge][eta_bin]->GetYaxis()->SetTitle("y(n#sigma_{#pi},m^{2})");
	  h_XY_sub[pt_bin][cent][charge][eta_bin]->GetYaxis()->CenterTitle();
	  h_XY_sub[pt_bin][cent][charge][eta_bin]->GetYaxis()->SetTitleOffset(1.6);
	  if(pt_bin >=  5)
	  {
	    h_XY_sub[pt_bin][cent][charge][eta_bin]->Add(f_proton[pt_bin][cent][charge][eta_bin],-1.0);
	  }
	  h_XY_sub[pt_bin][cent][charge][eta_bin]->Draw("colz");
	  if(pt_bin >=  5)
	  {
	    // mass2 cut line start
	    Double_t x1, x2, y1, y2;
	    setInitValues((pt_low[pt_bin]+pt_up[pt_bin])/2.0,-3.0,0.60);
	    x1 = getNewX(); // x component of first point
	    y1 = getNewY(); // y component of first point
	    setInitValues((pt_low[pt_bin]+pt_up[pt_bin])/2.0,3.0,0.60);
	    x2 = getNewX(); // x component of second point
	    y2 = getNewY(); // y component of second point
	    TF1 *cutline = new TF1("cutline",getLine,x_low[pt_bin],x_up[pt_bin],4);
	    cutline->FixParameter(0,x1);
	    cutline->FixParameter(1,y1);
	    cutline->FixParameter(2,x2);
	    cutline->FixParameter(3,y2);
	    cutline->SetLineStyle(2);
	    cutline->SetLineColor(1);
	    cutline->SetLineWidth(2);
	    cutline->DrawCopy("l same");
	  }
	}
      }
    }
  }

  // mass2 cut
  TH2F *h_XY_cut[pt_total][cent_total][charge_total][eta_total];
  for(Int_t pt_bin = 0; pt_bin < pt_total; pt_bin ++)
  {
    for(Int_t cent = cent_start[mCentrality]; cent < cent_stop[mCentrality]; cent++)
    {
      for(Int_t charge = charge_start; charge < charge_stop; charge++)
      {
	for(Int_t eta_bin = eta_start; eta_bin < eta_stop; eta_bin++)
	{
	  TString HistName = Form("pt_%d_Centrality_%d_Charge_%d_etagap_%d_cut",pt_bin,cent,charge,eta_bin);
	  h_XY_cut[pt_bin][cent][charge][eta_bin] = (TH2F*)h_XY_sub[pt_bin][cent][charge][eta_bin]->Clone(HistName.Data());
	}
      }
    }
  }
  TCanvas *c_XY_cut[cent_total][charge_total][eta_total];
  for(Int_t cent = cent_start[mCentrality]; cent < cent_stop[mCentrality]; cent++)
  {
    for(Int_t charge = charge_start; charge < charge_stop; charge++)
    {
      for(Int_t eta_bin = eta_start; eta_bin < eta_stop; eta_bin++)
      {
	TString CanName = Form("c_XY_cent_%d_charge_%d_etagap_%d_cut",cent,charge,eta_bin);
	c_XY_cut[cent][charge][eta_bin] = new TCanvas(CanName.Data(),CanName.Data(),1400,10,1200,1200);
	c_XY_cut[cent][charge][eta_bin]->Divide(4,4);
	for(Int_t pt_bin = 0; pt_bin < pt_total; pt_bin++)
	{
	  c_XY_cut[cent][charge][eta_bin]->cd(pt_bin+1);
	  c_XY_cut[cent][charge][eta_bin]->cd(pt_bin+1)->SetLeftMargin(0.15);
	  c_XY_cut[cent][charge][eta_bin]->cd(pt_bin+1)->SetRightMargin(0.15);
	  c_XY_cut[cent][charge][eta_bin]->cd(pt_bin+1)->SetTopMargin(0.15);
	  c_XY_cut[cent][charge][eta_bin]->cd(pt_bin+1)->SetBottomMargin(0.15);
	  c_XY_cut[cent][charge][eta_bin]->cd(pt_bin+1)->SetLogz();
	  c_XY_cut[cent][charge][eta_bin]->cd(pt_bin+1)->SetTicks(1,1);
	  h_XY_cut[pt_bin][cent][charge][eta_bin]->SetStats(0);
	  h_XY_cut[pt_bin][cent][charge][eta_bin]->SetTitle(Title[pt_bin][cent][charge][eta_bin].Data());
	  h_XY_cut[pt_bin][cent][charge][eta_bin]->GetXaxis()->SetTitle("x(n#sigma_{#pi},m^{2})");
	  h_XY_cut[pt_bin][cent][charge][eta_bin]->GetXaxis()->CenterTitle();
	  h_XY_cut[pt_bin][cent][charge][eta_bin]->GetYaxis()->SetTitle("y(n#sigma_{#pi},m^{2})");
	  h_XY_cut[pt_bin][cent][charge][eta_bin]->GetYaxis()->CenterTitle();
	  h_XY_cut[pt_bin][cent][charge][eta_bin]->GetYaxis()->SetTitleOffset(1.6);
	  if(pt_bin >=  5)
	  {
	    // mass2 cut line start
	    Double_t x1, x2, y1, y2;
	    setInitValues((pt_low[pt_bin]+pt_up[pt_bin])/2.0,-3.0,0.60);
	    x1 = getNewX(); // x component of first point
	    y1 = getNewY(); // y component of first point
	    setInitValues((pt_low[pt_bin]+pt_up[pt_bin])/2.0,3.0,0.60);
	    x2 = getNewX(); // x component of second point
	    y2 = getNewY(); // y component of second point
	    TF1 *cutline = new TF1("cutline",getLine,x_low[pt_bin],x_up[pt_bin],4);
	    cutline->FixParameter(0,x1);
	    cutline->FixParameter(1,y1);
	    cutline->FixParameter(2,x2);
	    cutline->FixParameter(3,y2);

	    Int_t nbinx = h_XY_cut[pt_bin][cent][charge][eta_bin]->GetNbinsX();
	    Int_t nbiny = h_XY_cut[pt_bin][cent][charge][eta_bin]->GetNbinsY();

	    for(Int_t binx = 1; binx < nbinx; binx++)
	    {
	      Float_t binx_center = h_XY_cut[pt_bin][cent][charge][eta_bin]->GetXaxis()->GetBinCenter(binx);
	      for(Int_t biny = 1; biny < nbiny; biny++)
	      {
		Float_t biny_center = h_XY_cut[pt_bin][cent][charge][eta_bin]->GetYaxis()->GetBinCenter(biny);
		Float_t y_val = cutline->Eval(binx_center);
		if(biny_center < y_val)
		{
		  Int_t global_bin = h_XY_cut[pt_bin][cent][charge][eta_bin]->GetBin(binx,biny);
		  h_XY_cut[pt_bin][cent][charge][eta_bin]->SetBinContent(global_bin,0.0);
		}
	      }
	    }
	  }
	  h_XY_cut[pt_bin][cent][charge][eta_bin]->Draw("colz");
	}
      }
    }
  }
  //-----------------------------------------------------------------------------------------

  //-----------------------------------------------------------------------------------------
  // Pion and Kaon
  TF2 *f_double[pt_total][cent_total][charge_total][eta_total];
  for(Int_t pt_bin = 5; pt_bin < pt_total; pt_bin ++)
  {
    for(Int_t cent = cent_start[mCentrality]; cent < cent_stop[mCentrality]; cent++)
    {
      for(Int_t charge = charge_start; charge < charge_stop; charge++)
      {
	for(Int_t eta_bin = eta_start; eta_bin < eta_stop; eta_bin++)
	{
	  TString FuncName = Form("f_double_pt_%d_cent_%d_charge_%d_etagap_%d",pt_bin,cent,charge,eta_bin);
	  f_double[pt_bin][cent][charge][eta_bin] = new TF2(FuncName.Data(),student_t_2d_double,x_low[pt_bin],x_up[pt_bin],y_low[pt_bin],y_up[pt_bin],14);
	  f_double[pt_bin][cent][charge][eta_bin]->SetNpx(400);
	  f_double[pt_bin][cent][charge][eta_bin]->SetNpy(400);
	  f_double[pt_bin][cent][charge][eta_bin]->SetContour(15);
	  for(Int_t i = 0; i < 7; i++)
	  {
	    f_double[pt_bin][cent][charge][eta_bin]->SetParameter(i,parfit[pt_bin][cent][charge][eta_bin][i]);
	  }
	  f_double[pt_bin][cent][charge][eta_bin]->SetParameter(7,(2.0*parfit[pt_bin][cent][charge][eta_bin][0]+parfit[pt_bin][cent][charge][eta_bin][12])/3.0);
	  f_double[pt_bin][cent][charge][eta_bin]->SetParameter(8,(2.0*parfit[pt_bin][cent][charge][eta_bin][1]+parfit[pt_bin][cent][charge][eta_bin][13])/3.0);
	  for(Int_t i = 0; i < 5; i++)
	  {
	    f_double[pt_bin][cent][charge][eta_bin]->SetParameter(i+9,parfit[pt_bin][cent][charge][eta_bin][i+7]);
	  }
	  f_double[pt_bin][cent][charge][eta_bin]->SetLineStyle(1);
	  f_double[pt_bin][cent][charge][eta_bin]->SetLineColor(2);
	}
      }
    }
  }

  // Canvas
  TCanvas *c_X_cut[cent_total][charge_total][eta_total];
  TCanvas *c_Y_cut[cent_total][charge_total][eta_total];
  for(Int_t cent = cent_start[mCentrality]; cent < cent_stop[mCentrality]; cent++)
  {
    for(Int_t charge = charge_start; charge < charge_stop; charge++)
    {
      for(Int_t eta_bin = eta_start; eta_bin < eta_stop; eta_bin++)
      {
	TString CanName = Form("c_X_cent_%d_charge_%d_etagap_%d_cut",cent,charge,eta_bin);
	c_X_cut[cent][charge][eta_bin] = new TCanvas(CanName.Data(),CanName.Data(),1400,10,1200,1200);
	c_X_cut[cent][charge][eta_bin]->Divide(4,4);
	CanName = Form("c_Y_cent_%d_charge_%d_etagap_%d_cut",cent,charge,eta_bin);
	c_Y_cut[cent][charge][eta_bin] = new TCanvas(CanName.Data(),CanName.Data(),1400,10,1200,1200);
	c_Y_cut[cent][charge][eta_bin]->Divide(4,4);
	for(Int_t pt_bin = 0; pt_bin < pt_total; pt_bin++)
	{
	  c_X_cut[cent][charge][eta_bin]->cd(pt_bin+1);
	  c_X_cut[cent][charge][eta_bin]->cd(pt_bin+1)->SetLeftMargin(0.15);
	  c_X_cut[cent][charge][eta_bin]->cd(pt_bin+1)->SetRightMargin(0.15);
	  c_X_cut[cent][charge][eta_bin]->cd(pt_bin+1)->SetTopMargin(0.15);
	  c_X_cut[cent][charge][eta_bin]->cd(pt_bin+1)->SetBottomMargin(0.15);
	  c_X_cut[cent][charge][eta_bin]->cd(pt_bin+1)->SetLogz();
	  c_X_cut[cent][charge][eta_bin]->cd(pt_bin+1)->SetTicks(1,1);

	  c_Y_cut[cent][charge][eta_bin]->cd(pt_bin+1);
	  c_Y_cut[cent][charge][eta_bin]->cd(pt_bin+1)->SetLeftMargin(0.15);
	  c_Y_cut[cent][charge][eta_bin]->cd(pt_bin+1)->SetRightMargin(0.15);
	  c_Y_cut[cent][charge][eta_bin]->cd(pt_bin+1)->SetTopMargin(0.15);
	  c_Y_cut[cent][charge][eta_bin]->cd(pt_bin+1)->SetBottomMargin(0.15);
	  c_Y_cut[cent][charge][eta_bin]->cd(pt_bin+1)->SetLogz();
	  c_Y_cut[cent][charge][eta_bin]->cd(pt_bin+1)->SetTicks(1,1);
	}
      }
    }
  }

  // Projection X
  TH1F *h_X_cut[pt_total][cent_total][charge_total][eta_total];
  TH1F *h_X_cut_Fit[pt_total][cent_total][charge_total][eta_total];
  TH1F *h_X_cut_Pion[pt_total][cent_total][charge_total][eta_total];
  TH1F *h_X_cut_Kaon[pt_total][cent_total][charge_total][eta_total];
  for(Int_t pt_bin = 0; pt_bin < pt_total; pt_bin++)
  {
    for(Int_t cent = cent_start[mCentrality]; cent < cent_stop[mCentrality]; cent++)
    {
      for(Int_t charge = charge_start; charge < charge_stop; charge++)
      {
	for(Int_t eta_bin = eta_start; eta_bin < eta_stop; eta_bin++)
	{
	  if(pt_bin < 5)
	  {
	    TString HistName = Form("pt_%d_Centrality_%d_charge_%d_etagap_%d_X_cut",pt_bin,cent,charge,eta_bin);
	    h_X_cut[pt_bin][cent][charge][eta_bin] = (TH1F*)h_XY_cut[pt_bin][cent][charge][eta_bin]->ProjectionX(HistName.Data(),h_XY_cut[pt_bin][cent][charge][eta_bin]->GetYaxis()->FindBin(y_low[pt_bin]),h_XY_cut[pt_bin][cent][charge][eta_bin]->GetYaxis()->FindBin(y_up[pt_bin]));
	    c_X_cut[cent][charge][eta_bin]->cd(pt_bin+1);
	    h_X_cut[pt_bin][cent][charge][eta_bin]->SetStats(0);
	    h_X_cut[pt_bin][cent][charge][eta_bin]->DrawCopy("PE");
	  }
	  if(pt_bin >= 5)
	  {
	    TString HistName = Form("pt_%d_Centrality_%d_charge_%d_etagap_%d_X_cut",pt_bin,cent,charge,eta_bin);
	    h_X_cut[pt_bin][cent][charge][eta_bin] = (TH1F*)h_XY_cut[pt_bin][cent][charge][eta_bin]->ProjectionX(HistName.Data(),h_XY_cut[pt_bin][cent][charge][eta_bin]->GetYaxis()->FindBin(y_p[pt_bin]-order_proton[pt_bin]*parfit[pt_bin][cent][charge][eta_bin][17]),h_XY_cut[pt_bin][cent][charge][eta_bin]->GetYaxis()->FindBin(y_pi[pt_bin]+order_pion[pt_bin]*parfit[pt_bin][cent][charge][eta_bin][5]));
	    c_X_cut[cent][charge][eta_bin]->cd(pt_bin+1);
	    h_X_cut[pt_bin][cent][charge][eta_bin]->SetStats(0);
	    h_X_cut[pt_bin][cent][charge][eta_bin]->DrawCopy("PE");

	    Int_t start_bin_X = ((TH2F*)f_double[pt_bin][cent][charge][eta_bin]->GetHistogram())->GetYaxis()->FindBin(y_p[pt_bin]-order_proton[pt_bin]*parfit[pt_bin][cent][charge][eta_bin][17]);
	    Int_t stop_bin_X  = ((TH2F*)f_double[pt_bin][cent][charge][eta_bin]->GetHistogram())->GetYaxis()->FindBin(y_pi[pt_bin]+order_pion[pt_bin]*parfit[pt_bin][cent][charge][eta_bin][5]);
	    Double_t Inte_X = h_X_cut[pt_bin][cent][charge][eta_bin]->Integral(h_X_cut[pt_bin][cent][charge][eta_bin]->FindBin(x_pi[pt_bin]-order_pion[pt_bin]*parfit[pt_bin][cent][charge][eta_bin][4]),h_X_cut[pt_bin][cent][charge][eta_bin]->FindBin(x_p[pt_bin]+order_proton[pt_bin]*parfit[pt_bin][cent][charge][eta_bin][16]));


	    HistName = Form("pt_%d_Centrality_%d_charge_%d_etagap_%d_Fit_X",pt_bin,cent,charge,eta_bin);
	    h_X_cut_Fit[pt_bin][cent][charge][eta_bin] = (TH1F*)((TH2F*)f_double[pt_bin][cent][charge][eta_bin]->GetHistogram())->ProjectionX(HistName.Data(),start_bin_X,stop_bin_X);
	    Double_t Inte_X_Fit = h_X_cut_Fit[pt_bin][cent][charge][eta_bin]->Integral(h_X_cut_Fit[pt_bin][cent][charge][eta_bin]->FindBin(x_pi[pt_bin]-order_pion[pt_bin]*parfit[pt_bin][cent][charge][eta_bin][4]),h_X_cut_Fit[pt_bin][cent][charge][eta_bin]->FindBin(x_p[pt_bin]+order_proton[pt_bin]*parfit[pt_bin][cent][charge][eta_bin][16]));
	    Double_t scale_X = 1.0;
	    if(Inte_X_Fit > 0.0) scale_X = Inte_X/Inte_X_Fit;
	    h_X_cut_Fit[pt_bin][cent][charge][eta_bin]->Scale(scale_X*h_X_cut[pt_bin][cent][charge][eta_bin]->GetBinWidth(1)/h_X_cut_Fit[pt_bin][cent][charge][eta_bin]->GetBinWidth(1));
	    h_X_cut_Fit[pt_bin][cent][charge][eta_bin]->SetLineColor(2);
	    h_X_cut_Fit[pt_bin][cent][charge][eta_bin]->SetLineStyle(2);
	    h_X_cut_Fit[pt_bin][cent][charge][eta_bin]->SetFillColor(0);

	    HistName = Form("pt_%d_Centrality_%d_charge_%d_etagap_%d_Pion_X",pt_bin,cent,charge,eta_bin);
	    h_X_cut_Pion[pt_bin][cent][charge][eta_bin] = (TH1F*)((TH2F*)f_pion[pt_bin][cent][charge][eta_bin]->GetHistogram())->ProjectionX(HistName.Data(),start_bin_X,stop_bin_X);
	    h_X_cut_Pion[pt_bin][cent][charge][eta_bin]->Scale(scale_X*h_X_cut[pt_bin][cent][charge][eta_bin]->GetBinWidth(1)/h_X_cut_Pion[pt_bin][cent][charge][eta_bin]->GetBinWidth(1));
	    h_X_cut_Pion[pt_bin][cent][charge][eta_bin]->SetFillColor(fill_color_pion);
	    h_X_cut_Pion[pt_bin][cent][charge][eta_bin]->SetLineColor(fill_color_pion);
	    h_X_cut_Pion[pt_bin][cent][charge][eta_bin]->SetFillStyle(fill_style_pion);

	    HistName = Form("pt_%d_Centrality_%d_charge_%d_etagap_%d_Kaon_X",pt_bin,cent,charge,eta_bin);
	    h_X_cut_Kaon[pt_bin][cent][charge][eta_bin] = (TH1F*)((TH2F*)f_kaon[pt_bin][cent][charge][eta_bin]->GetHistogram())->ProjectionX(HistName.Data(),start_bin_X,stop_bin_X);
	    h_X_cut_Kaon[pt_bin][cent][charge][eta_bin]->Scale(scale_X*h_X_cut[pt_bin][cent][charge][eta_bin]->GetBinWidth(1)/h_X_cut_Kaon[pt_bin][cent][charge][eta_bin]->GetBinWidth(1));
	    h_X_cut_Kaon[pt_bin][cent][charge][eta_bin]->SetFillColor(fill_color_kaon);
	    h_X_cut_Kaon[pt_bin][cent][charge][eta_bin]->SetLineColor(fill_color_kaon);
	    h_X_cut_Kaon[pt_bin][cent][charge][eta_bin]->SetFillStyle(fill_style_kaon);

	    h_X_cut_Pion[pt_bin][cent][charge][eta_bin]->DrawCopy("hE same");
	    h_X_cut_Kaon[pt_bin][cent][charge][eta_bin]->DrawCopy("hE same");
	    h_X_cut_Fit[pt_bin][cent][charge][eta_bin]->DrawCopy("l same");
	  }
	}
      }
    }
  }

  // Projection Y
  TH1F *h_Y_cut[pt_total][cent_total][charge_total][eta_total];
  TH1F *h_Y_cut_Fit[pt_total][cent_total][charge_total][eta_total];
  TH1F *h_Y_cut_Pion[pt_total][cent_total][charge_total][eta_total];
  TH1F *h_Y_cut_Kaon[pt_total][cent_total][charge_total][eta_total];
  for(Int_t pt_bin = 0; pt_bin < pt_total; pt_bin++)
  {
    for(Int_t cent = cent_start[mCentrality]; cent < cent_stop[mCentrality]; cent++)
    {
      for(Int_t charge = charge_start; charge < charge_stop; charge++)
      {
	for(Int_t eta_bin = eta_start; eta_bin < eta_stop; eta_bin++)
	{
	  if(pt_bin < 5)
	  {
	    TString HistName = Form("pt_%d_Centrality_%d_charge_%d_etagap_%d_Y_cut",pt_bin,cent,charge,eta_bin);
	    h_Y_cut[pt_bin][cent][charge][eta_bin] = (TH1F*)h_XY_cut[pt_bin][cent][charge][eta_bin]->ProjectionY(HistName.Data(),h_XY_cut[pt_bin][cent][charge][eta_bin]->GetXaxis()->FindBin(x_low[pt_bin]),h_XY_cut[pt_bin][cent][charge][eta_bin]->GetXaxis()->FindBin(x_up[pt_bin]));
	    c_Y_cut[cent][charge][eta_bin]->cd(pt_bin+1);
	    h_Y_cut[pt_bin][cent][charge][eta_bin]->SetStats(0);
	    h_Y_cut[pt_bin][cent][charge][eta_bin]->DrawCopy("PE");
	  }
	  if(pt_bin >= 5)
	  {
	    TString HistName = Form("pt_%d_Centrality_%d_charge_%d_etagap_%d_Y_cut",pt_bin,cent,charge,eta_bin);
	    h_Y_cut[pt_bin][cent][charge][eta_bin] = (TH1F*)h_XY_cut[pt_bin][cent][charge][eta_bin]->ProjectionY(HistName.Data(),h_XY_cut[pt_bin][cent][charge][eta_bin]->GetXaxis()->FindBin(x_pi[pt_bin]-order_pion[pt_bin]*parfit[pt_bin][cent][charge][eta_bin][4]),h_XY_cut[pt_bin][cent][charge][eta_bin]->GetXaxis()->FindBin(x_p[pt_bin]+order_proton[pt_bin]*parfit[pt_bin][cent][charge][eta_bin][16]));
	    c_Y_cut[cent][charge][eta_bin]->cd(pt_bin+1);
	    h_Y_cut[pt_bin][cent][charge][eta_bin]->SetStats(0);
	    h_Y_cut[pt_bin][cent][charge][eta_bin]->DrawCopy("PE");

	    Int_t start_bin_Y = ((TH2F*)f_double[pt_bin][cent][charge][eta_bin]->GetHistogram())->GetXaxis()->FindBin(x_pi[pt_bin]-order_pion[pt_bin]*parfit[pt_bin][cent][charge][eta_bin][4]);
	    Int_t stop_bin_Y  = ((TH2F*)f_double[pt_bin][cent][charge][eta_bin]->GetHistogram())->GetYaxis()->FindBin(x_p[pt_bin]+order_proton[pt_bin]*parfit[pt_bin][cent][charge][eta_bin][16]);
	    Double_t Inte_Y = h_Y_cut[pt_bin][cent][charge][eta_bin]->Integral(h_Y_cut[pt_bin][cent][charge][eta_bin]->FindBin(y_p[pt_bin]-order_proton[pt_bin]*parfit[pt_bin][cent][charge][eta_bin][17]),h_Y_cut[pt_bin][cent][charge][eta_bin]->FindBin(y_pi[pt_bin]+order_pion[pt_bin]*parfit[pt_bin][cent][charge][eta_bin][5]));

	    HistName = Form("pt_%d_Centrality_%d_charge_%d_etagap_%d_Fit_Y",pt_bin,cent,charge,eta_bin);
	    h_Y_cut_Fit[pt_bin][cent][charge][eta_bin] = (TH1F*)((TH2F*)f_double[pt_bin][cent][charge][eta_bin]->GetHistogram())->ProjectionY(HistName.Data(),start_bin_Y,stop_bin_Y);
	    Double_t Inte_Y_Fit = h_Y_cut_Fit[pt_bin][cent][charge][eta_bin]->Integral(h_Y_cut_Fit[pt_bin][cent][charge][eta_bin]->FindBin(y_p[pt_bin]-order_proton[pt_bin]*parfit[pt_bin][cent][charge][eta_bin][17]),h_Y_cut_Fit[pt_bin][cent][charge][eta_bin]->FindBin(y_pi[pt_bin]+order_pion[pt_bin]*parfit[pt_bin][cent][charge][eta_bin][5]));
	    Double_t scale_Y = 1.0;
	    if(Inte_Y_Fit > 0.0) scale_Y = Inte_Y/Inte_Y_Fit;
	    h_Y_cut_Fit[pt_bin][cent][charge][eta_bin]->Scale(scale_Y*h_Y_cut[pt_bin][cent][charge][eta_bin]->GetBinWidth(1)/h_Y_cut_Fit[pt_bin][cent][charge][eta_bin]->GetBinWidth(1));
	    h_Y_cut_Fit[pt_bin][cent][charge][eta_bin]->SetLineColor(2);
	    h_Y_cut_Fit[pt_bin][cent][charge][eta_bin]->SetLineStyle(2);
	    h_Y_cut_Fit[pt_bin][cent][charge][eta_bin]->SetFillColor(0);

	    HistName = Form("pt_%d_Centrality_%d_charge_%d_etagap_%d_Pion_Y",pt_bin,cent,charge,eta_bin);
	    h_Y_cut_Pion[pt_bin][cent][charge][eta_bin] = (TH1F*)((TH2F*)f_pion[pt_bin][cent][charge][eta_bin]->GetHistogram())->ProjectionY(HistName.Data(),start_bin_Y,stop_bin_Y);
	    h_Y_cut_Pion[pt_bin][cent][charge][eta_bin]->Scale(scale_Y*h_Y_cut[pt_bin][cent][charge][eta_bin]->GetBinWidth(1)/h_Y_cut_Pion[pt_bin][cent][charge][eta_bin]->GetBinWidth(1));
	    h_Y_cut_Pion[pt_bin][cent][charge][eta_bin]->SetFillColor(fill_color_pion);
	    h_Y_cut_Pion[pt_bin][cent][charge][eta_bin]->SetLineColor(fill_color_pion);
	    h_Y_cut_Pion[pt_bin][cent][charge][eta_bin]->SetFillStyle(fill_style_pion);

	    HistName = Form("pt_%d_Centrality_%d_charge_%d_etagap_%d_Kaon_Y",pt_bin,cent,charge,eta_bin);
	    h_Y_cut_Kaon[pt_bin][cent][charge][eta_bin] = (TH1F*)((TH2F*)f_kaon[pt_bin][cent][charge][eta_bin]->GetHistogram())->ProjectionY(HistName.Data(),start_bin_Y,stop_bin_Y);
	    h_Y_cut_Kaon[pt_bin][cent][charge][eta_bin]->Scale(scale_Y*h_Y_cut[pt_bin][cent][charge][eta_bin]->GetBinWidth(1)/h_Y_cut_Kaon[pt_bin][cent][charge][eta_bin]->GetBinWidth(1));
	    h_Y_cut_Kaon[pt_bin][cent][charge][eta_bin]->SetFillColor(fill_color_kaon);
	    h_Y_cut_Kaon[pt_bin][cent][charge][eta_bin]->SetLineColor(fill_color_kaon);
	    h_Y_cut_Kaon[pt_bin][cent][charge][eta_bin]->SetFillStyle(fill_style_kaon);

	    h_Y_cut_Pion[pt_bin][cent][charge][eta_bin]->DrawCopy("hE same");
	    h_Y_cut_Kaon[pt_bin][cent][charge][eta_bin]->DrawCopy("hE same");
	    h_Y_cut_Fit[pt_bin][cent][charge][eta_bin]->DrawCopy("l same");
	  }
	}
      }
    }
  }
  //-----------------------------------------------------------------------------------------

  cout << "Save Histogram to the file!!" << endl;
  TString output_hist = Form("/project/projectdirs/star/xusun/OutPut/AuAu%s/Mass2_nSigmaPion/merged_file/pt_spectra/Mass2_nSigmaPion_cut_%s_charge_%d_etagap_%d.root",Energy[mEnergy].Data(),Centrality[mCentrality].Data(),charge_start,eta_start);
  TFile *File_hist = new TFile(output_hist.Data(),"RECREATE");
  File_hist->cd();
  for(Int_t pt_bin = 0; pt_bin < pt_total; pt_bin ++)
  {
    for(Int_t cent = cent_start[mCentrality]; cent < cent_stop[mCentrality]; cent++)
    {
      for(Int_t charge = charge_start; charge < charge_stop; charge++)
      {
	for(Int_t eta_bin = eta_start; eta_bin < eta_stop; eta_bin++)
	{
	  h_XY_cut[pt_bin][cent][charge][eta_bin]->Write();
	}
      }
    }
  }
  File_hist->Close();

  TString output_par = Form("/project/projectdirs/star/xusun/OutPut/AuAu%s/Mass2_nSigmaPion/merged_file/pt_spectra/parfit_%s_charge_%d_etagap_%d.root",Energy[mEnergy].Data(),Centrality[mCentrality].Data(),charge_start,eta_start);
  cout << output_par.Data() << endl;
  TFile *File_par = new TFile(output_par.Data(),"RECREATE");
  File_par->cd();
  for(Int_t pt_bin = 0; pt_bin < pt_total; pt_bin ++)
  {
    for(Int_t cent = cent_start[mCentrality]; cent < cent_stop[mCentrality]; cent++)
    {
      for(Int_t charge = charge_start; charge < charge_stop; charge++)
      {
	for(Int_t eta_bin = eta_start; eta_bin < eta_stop; eta_bin++)
	{
	  h_parfit_1st[pt_bin][cent][charge][eta_bin]->Write();
	  h_parfit_2nd[pt_bin][cent][charge][eta_bin]->Write();
	}
      }
    }
  }
  File_par->Close();

  for(Int_t cent = cent_start[mCentrality]; cent < cent_stop[mCentrality]; cent++)
  {
    for(Int_t charge = charge_start; charge < charge_stop; charge++)
    {
      for(Int_t eta_bin = eta_start; eta_bin < eta_stop; eta_bin++)
      {
	TString CanName = Form("/project/projectdirs/star/xusun/OutPut/AuAu%s/Mass2_nSigmaPion/merged_file/pt_spectra/figures/c_XY_cent_%d_charge_%d_etagap_%d_EP.gif",Energy[mEnergy].Data(),cent,charge,eta_bin);
	c_XY[cent][charge][eta_bin]->SaveAs(CanName.Data());

	TString CanName_X = Form("/project/projectdirs/star/xusun/OutPut/AuAu%s/Mass2_nSigmaPion/merged_file/pt_spectra/figures/c_X_cent_%d_charge_%d_etagap_%d_EP.gif",Energy[mEnergy].Data(),cent,charge,eta_bin);
	c_X[cent][charge][eta_bin]->SaveAs(CanName_X.Data());
	TString CanName_Y = Form("/project/projectdirs/star/xusun/OutPut/AuAu%s/Mass2_nSigmaPion/merged_file/pt_spectra/figures/c_Y_cent_%d_charge_%d_etagap_%d_EP.gif",Energy[mEnergy].Data(),cent,charge,eta_bin);
	c_Y[cent][charge][eta_bin]->SaveAs(CanName_Y.Data());

	TString CanName_sub = Form("/project/projectdirs/star/xusun/OutPut/AuAu%s/Mass2_nSigmaPion/merged_file/pt_spectra/figures/c_XY_cent_%d_charge_%d_etagap_%d_EP_sub.gif",Energy[mEnergy].Data(),cent,charge,eta_bin);
	c_XY_sub[cent][charge][eta_bin]->SaveAs(CanName_sub.Data());

	TString CanName_cut = Form("/project/projectdirs/star/xusun/OutPut/AuAu%s/Mass2_nSigmaPion/merged_file/pt_spectra/figures/c_XY_cent_%d_charge_%d_etagap_%d_EP_cut.gif",Energy[mEnergy].Data(),cent,charge,eta_bin);
	c_XY_cut[cent][charge][eta_bin]->SaveAs(CanName_cut.Data());
	TString CanName_cut_X = Form("/project/projectdirs/star/xusun/OutPut/AuAu%s/Mass2_nSigmaPion/merged_file/pt_spectra/figures/c_X_cent_%d_charge_%d_etagap_%d_EP_cut.gif",Energy[mEnergy].Data(),cent,charge,eta_bin);
	c_X_cut[cent][charge][eta_bin]->SaveAs(CanName_cut_X.Data());
	TString CanName_cut_Y = Form("/project/projectdirs/star/xusun/OutPut/AuAu%s/Mass2_nSigmaPion/merged_file/pt_spectra/figures/c_Y_cent_%d_charge_%d_etagap_%d_EP_cut.gif",Energy[mEnergy].Data(),cent,charge,eta_bin);
	c_Y_cut[cent][charge][eta_bin]->SaveAs(CanName_cut_Y.Data());

	//	    TString CanName_pik = Form("/u/xusun/STAR/Analysis_v3/StRoot/CombPiKP//figures/c_pik_cent_%d_charge_%d_etagap_%d_phi_psi_%d_%s_EP.gif",cent,charge,eta_bin,phi_psi_bin,Order[flow_bin].Data());
	//	    c_X_pik[cent][charge][eta_bin][phi_psi_bin][flow_bin]->SaveAs(CanName_pik.Data());
      }
    }
  }
  cout << "Work done!! Now it's time to close!!" << endl;
}

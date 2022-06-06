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
#include <map>
#include <vector>
#include "TF1.h"

#include "./student_t_2d_fit.h"
#include "./student_t_2d_single.h"
#include "./student_t_2d_double.h"
#include "./student_t_1d_single.h"
#include "./student_t_1d_double.h"
#include "./draw.h"
#include "./student_t_combPID.h"

static const TString Energy[2] = {"200GeV","39GeV"};
static const TString Order[2] = {"2nd","3rd"};
static const Double_t PI_max[2] = {TMath::Pi()/2.0,TMath::Pi()/3.0};
static const Float_t nSigProton = 3.0;
static const Float_t Flow_Order[2] = {2.0,3.0};

// pt bin
static const Int_t pt_total = 16;
//                                            0 , 1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 ,10 ,11 ,12 ,13 ,14 ,15
static const Float_t pt_low_raw[pt_total] = {0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4};
static const Float_t pt_up_raw[pt_total]  = {0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,4.2};

// pt rebin
static const Int_t pt_rebin_total = 16;
static const Float_t pt_low[pt_rebin_total] = {0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4};
static const Float_t pt_up[pt_rebin_total]  = {0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,4.2};
static const Int_t pt_rebin_start[pt_rebin_total] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
static const Int_t pt_rebin_stop[pt_rebin_total]  = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
static const Int_t pt_rebin_first = 0;
static const Int_t pt_rebin_last  = 16;
static const Int_t pt_QA    = 12;

// x and y range
static const Float_t x_low[pt_total] = {-0.4, -0.6,-0.6, -0.6,-0.6,-0.6,-0.6,-0.8,-0.8,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0};
static const Float_t x_up[pt_total]  = { 1.4,  1.4, 1.4,  1.4, 1.4, 1.4, 1.4, 1.6, 1.6, 2.0, 2.0, 2.0, 2.4, 2.4, 2.4, 2.4};
static const Float_t y_low[pt_total] = {-0.15,-0.2,-0.3,-0.45,-0.6,-0.8,-1.0,-1.3,-1.5,-1.7,-1.8,-1.8,-2.0,-2.0,-2.0,-2.0};
static const Float_t y_up[pt_total]  = { 0.15, 0.2, 0.2, 0.20, 0.3, 0.4, 0.6, 1.0, 1.0, 1.0, 1.2, 1.2, 1.4, 1.4, 1.4, 1.4};

static const Int_t Cent_total = 4; 
static const Int_t Cent_start = 0;
static const Int_t Cent_stop  = 1;

static const Int_t Eta_total = 4;
static const Int_t Eta_start = 0;
static const Int_t Eta_stop = 1;

static const Int_t phi_total = 7; // phi loop
static const Int_t phi_start = 0;
static const Int_t phi_stop  = 7;
static const Int_t phi_QA = 4;

static const Int_t Sys_total = 6;
static const Int_t Sys_start = 0;
static const Int_t Sys_stop  = 6;
static const Int_t Sys_QA = 3;

static const Int_t Proj_total = 3;
static const Int_t Proj_start = 0;
static const Int_t Proj_stop  = 3;
static const Int_t Proj_QA    = 2;
static const Float_t nSigProj[Proj_total] = {3.0,2.5,2.0};
static const Float_t nSigma = 2.0;

// Initial parameters 200 GeV
// pion
Float_t x_pi[pt_total] = {0.0025,-0.0025,-0.0025,-0.0025,-0.0025,-0.0075,-0.0125,-0.005,-0.005,-0.00625,0.00875,0.02375,0.0283,0.0344,0.0649,0.1135};
Float_t y_pi[pt_total] = {0.0005,0.00125,0.00125,-0.00025,0.002,0.00825,0.01075,0.02175,0.02825,0.03137,0.02625,0.03375,0.0114,0.0236,0.0236,0.0825};
// Float_t order_pion[pt_total]   = {3.0,3.0,3.0,3.0,3.0,2.55,2.55,3.0,3.0,3.0,2.5,3.00,3.00,3.00,3.00,3.0}; // positive
Float_t order_pion[pt_total]   = {3.0,3.0,3.0,3.0,3.0,2.50,2.55,3.0,3.0,3.0,2.5,3.00,2.55,2.50,2.55,3.0}; //negative 

// kaon
Float_t x_k[pt_total] = { 0.2275,  0.2275, 0.2225, 0.2225,0.2225, 0.2225, 0.2325, 0.2470, 0.2530, 0.2787, 0.3500, 0.3970,0.3519,0.4130,0.4924,0.8652};
Float_t y_k[pt_total] = {-0.0045,-0.00325,0.00125,0.00325,0.0060,0.00375,0.00525,0.02175,0.00225,0.03137,0.02625,0.03375,0.0114,0.0236,0.0236,0.0825};
// Float_t order_kaon[pt_total]   = {3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.05,3.05,3.00,3.0}; // positive
Float_t order_kaon[pt_total]   = {3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,2.55,3.00,3.00,3.0}; //negative

// proton
Float_t x_p[pt_total] = {0.8625, 0.8625,  0.8675,  0.8625,0.8525, 0.8305, 0.8025, 0.7750, 0.7450, 0.7213, 0.7287, 0.7620, 0.7794, 0.8343, 0.9565, 0.8955};
Float_t y_p[pt_total] = {0.0375,0.02075,-0.01125,-0.05975,-0.126,-0.2033,-0.2862,-0.3553,-0.4268,-0.4614,-0.5062,-0.5587,-0.5015,-0.5071,-0.4981,-0.5295};

typedef std::map<TString,TH1F*> TH1FMap;
typedef std::map<TString,TH2F*> TH2FMap;
typedef std::map<TString,TProfile*> TProMap;
typedef std::map<TString,TGraphAsymmErrors*> TGraMap;
typedef std::map<TString,std::vector<Float_t>> vecFMap;
typedef std::map<TString,TF1*> TF1Map;

#ifndef _PlotQA_
#define _PlotQA_ 0
#endif

// mMode:   0 for all pT bin, 1 for single pT bin
// mEnergy: 0 for 200GeV, 1 for 39GeV
// mCharge: 0 for positive, 1 for negative 
void CombPID_1d(Int_t mEnergy = 0, Int_t mCharge = 0, Int_t mOrder = 1)
{
  TString InPutFile = Form("./OutPut/AuAu%s/nSigmaPion/nSigmaPion_Sub_Charge_%d.root",Energy[mEnergy].Data(),mCharge);
  TFile *File_InPut = TFile::Open(InPutFile.Data());
  cout << "InPutFile: " << InPutFile.Data() << endl;

  TH2FMap h_mMass2_sub, h_mMass2;
  TH1FMap h_ParStudnet_total_1st, h_ParStudnet_total_2nd;
  vecFMap parfit_1st, parfit_2nd;
  parfit_1st.clear();
  parfit_2nd.clear();
  for(Int_t i_pt = pt_rebin_first; i_pt < pt_rebin_last; i_pt++) // pt bin
  {
    for(Int_t i_cent = Cent_start; i_cent < Cent_stop; i_cent++) // centrality bin
    {
      for(Int_t i_charge = mCharge; i_charge < mCharge+1; i_charge++) // charge bin
      {
	for(Int_t i_eta = Eta_start; i_eta < Eta_stop; i_eta++) // eta gap bin
	{
	  for(Int_t i_sys = Sys_start; i_sys < Sys_stop; i_sys++)
	  {
	    for(Int_t i_phi = phi_start; i_phi < phi_stop; i_phi++) // phi-psi bin
	    { // get histograms after proton subtraction
	      TString KEY_PiK = Form("pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_phi_Psi_%d_%s_PiK_SysError_%d",i_pt,i_cent,i_charge,i_eta,i_phi,Order[mOrder].Data(),i_sys);
	      h_mMass2_sub[KEY_PiK] = (TH2F*)File_InPut->Get(KEY_PiK.Data());
	      h_mMass2[KEY_PiK] = (TH2F*)h_mMass2_sub[KEY_PiK]->Clone();
	      if(i_pt > 4)
	      { // mass2 cut line start
		Double_t x1, x2, y1, y2;
		setInitValues((pt_low[i_pt]+pt_up[i_pt])/2.0,-3.0,0.60);
		x1 = getNewX(); // x component of first point
		y1 = getNewY(); // y component of first point
		setInitValues((pt_low[i_pt]+pt_up[i_pt])/2.0,3.0,0.60);
		x2 = getNewX(); // x component of second point
		y2 = getNewY(); // y component of second point
		TF1 *cutline = new TF1("cutline",getLine,x_low[i_pt],x_up[i_pt],4);
		cutline->FixParameter(0,x1);
		cutline->FixParameter(1,y1);
		cutline->FixParameter(2,x2);
		cutline->FixParameter(3,y2);

		Int_t nbinx = h_mMass2[KEY_PiK]->GetNbinsX();
		Int_t nbiny = h_mMass2[KEY_PiK]->GetNbinsY();

		for(Int_t binx = 1; binx < nbinx; binx++)
		{
		  Float_t binx_center = h_mMass2[KEY_PiK]->GetXaxis()->GetBinCenter(binx);
		  for(Int_t biny = 1; biny < nbiny; biny++)
		  {
		    Float_t biny_center = h_mMass2[KEY_PiK]->GetYaxis()->GetBinCenter(biny);
		    Float_t y_val = cutline->Eval(binx_center);
		    if(biny_center < y_val)
		    {
		      Int_t global_bin = h_mMass2[KEY_PiK]->GetBin(binx,biny);
		      h_mMass2[KEY_PiK]->SetBinContent(global_bin,0.0);
		    }
		  }
		}
	      }
	    }
	    if(i_pt > 4)
	    { // get fit parameters
	      TString KEY_PiK_Total = Form("pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_%s_PiK_SysError_%d",i_pt,i_cent,i_charge,i_eta,Order[mOrder].Data(),i_sys);
	      h_ParStudnet_total_1st[KEY_PiK_Total] = (TH1F*)File_InPut->Get((KEY_PiK_Total+"_1st").Data());
	      h_ParStudnet_total_2nd[KEY_PiK_Total] = (TH1F*)File_InPut->Get((KEY_PiK_Total+"_2nd").Data());
	      for(Int_t i_par = 0; i_par < 19; ++i_par)
	      {
		parfit_1st[KEY_PiK_Total].push_back(static_cast<Float_t>(h_ParStudnet_total_1st[KEY_PiK_Total]->GetBinContent(i_par+1)));
		parfit_2nd[KEY_PiK_Total].push_back(static_cast<Float_t>(h_ParStudnet_total_2nd[KEY_PiK_Total]->GetBinContent(i_par+1)));
	      }
	    }
	  }
	}
      }
    }
  }

#if _PlotQA_
  TCanvas *c_XY = new TCanvas("c_XY","c_XY",10,10,1600,800);
  c_XY->Divide(2,1);
  for(Int_t i_pad = 0; i_pad < 2; ++i_pad)
  {
    c_XY->cd(i_pad+1);
    c_XY->cd(i_pad+1)->SetLeftMargin(0.15);
    c_XY->cd(i_pad+1)->SetBottomMargin(0.15);
    c_XY->cd(i_pad+1)->SetTicks(1,1);
    c_XY->cd(i_pad+1)->SetGrid(0,0);
    c_XY->cd(i_pad+1)->SetLogz();
    TString KEY_PiK_QA = Form("pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_phi_Psi_%d_%s_PiK_SysError_%d",pt_QA,Cent_start,mCharge,Eta_start,phi_QA,Order[mOrder].Data(),Sys_QA); h_mMass2_sub[KEY_PiK_QA]->GetXaxis()->SetTitle("x(n#sigma_{#pi},m^{2})");
    if(i_pad == 0)
    {
      h_mMass2_sub[KEY_PiK_QA]->GetXaxis()->CenterTitle();
      h_mMass2_sub[KEY_PiK_QA]->GetYaxis()->SetTitle("y(n#sigma_{#pi},m^{2})");
      h_mMass2_sub[KEY_PiK_QA]->GetYaxis()->CenterTitle();
      h_mMass2_sub[KEY_PiK_QA]->Draw("colz");
    }
    if(i_pad == 1)
    {
      h_mMass2[KEY_PiK_QA]->GetXaxis()->SetTitle("x(n#sigma_{#pi},m^{2})");
      h_mMass2[KEY_PiK_QA]->GetXaxis()->CenterTitle();
      h_mMass2[KEY_PiK_QA]->GetYaxis()->SetTitle("y(n#sigma_{#pi},m^{2})");
      h_mMass2[KEY_PiK_QA]->GetYaxis()->CenterTitle();
      h_mMass2[KEY_PiK_QA]->Draw("colz");
    }
  }
#endif

  // integrated over phi and projected to x-axis
  TH2FMap h_mMass2_total;
  TH1FMap h_mXproj_total; 
  for(Int_t i_pt = pt_rebin_first; i_pt < pt_rebin_last; i_pt++) // pt bin
  {
    for(Int_t i_cent = Cent_start; i_cent < Cent_stop; i_cent++) // centrality bin
    {
      for(Int_t i_charge = mCharge; i_charge < mCharge+1; i_charge++) // charge bin
      {
	for(Int_t i_eta = Eta_start; i_eta < Eta_stop; i_eta++) // eta gap bin
	{
	  for(Int_t i_sys = Sys_start; i_sys < Sys_stop; i_sys++)
	  {
	    // integrated over phi
	    TString KEY_PiK_Total = Form("pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_%s_PiK_SysError_%d",i_pt,i_cent,i_charge,i_eta,Order[mOrder].Data(),i_sys);
	    for(Int_t i_phi = phi_start; i_phi < phi_stop; i_phi++) // phi-psi bin
	    {
	      TString KEY_PiK = Form("pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_phi_Psi_%d_%s_PiK_SysError_%d",i_pt,i_cent,i_charge,i_eta,i_phi,Order[mOrder].Data(),i_sys);
	      if(i_phi == phi_start) h_mMass2_total[KEY_PiK_Total] = (TH2F*)h_mMass2[KEY_PiK]->Clone();
	      else h_mMass2_total[KEY_PiK_Total]->Add(h_mMass2[KEY_PiK],1.0);
	    }

	    // projction to x-axis
	    for(Int_t i_proj = Proj_start; i_proj < Proj_stop; ++i_proj)
	    {
	      TString KEY_PiK_Tproj = Form("pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_%s_PiK_SysError_%d_Proj_%d",i_pt,i_cent,i_charge,i_eta,Order[mOrder].Data(),i_sys,i_proj);
	      if(i_pt > 4)
	      {
		Int_t bin_start = h_mMass2_total[KEY_PiK_Total]->GetYaxis()->FindBin(y_p[i_pt]-nSigProj[i_proj]*parfit_1st[KEY_PiK_Total][17]);
		Int_t bin_stop  = h_mMass2_total[KEY_PiK_Total]->GetYaxis()->FindBin(y_pi[i_pt]+nSigProj[i_proj]*parfit_1st[KEY_PiK_Total][5]);
		h_mXproj_total[KEY_PiK_Tproj] = (TH1F*)h_mMass2_total[KEY_PiK_Total]->ProjectionX(KEY_PiK_Tproj.Data(),bin_start,bin_stop);
	      }
	      else
	      {
		Int_t bin_start = h_mMass2_total[KEY_PiK_Total]->GetYaxis()->FindBin(y_low[i_pt]);
		Int_t bin_stop  = h_mMass2_total[KEY_PiK_Total]->GetYaxis()->FindBin(y_up[i_pt]);
		h_mXproj_total[KEY_PiK_Tproj] = (TH1F*)h_mMass2_total[KEY_PiK_Total]->ProjectionX(KEY_PiK_Tproj.Data(),bin_start,bin_stop);
	      }
	    }
	  }
	}
      }
    }
  }

#if _PlotQA_
  TCanvas *c_phi = new TCanvas("c_phi","c_phi",1200,1200);
  c_phi->Divide(3,3);
  for(Int_t i_pad = 0; i_pad < 9; ++i_pad)
  {
    c_phi->cd(i_pad+1);
    c_phi->cd(i_pad+1)->SetLeftMargin(0.15);
    c_phi->cd(i_pad+1)->SetBottomMargin(0.15);
    c_phi->cd(i_pad+1)->SetTicks(1,1);
    c_phi->cd(i_pad+1)->SetGrid(0,0);
    c_phi->cd(i_pad+1)->SetLogz();
    if(i_pad < phi_total)
    {
      TString KEY_PiK_QA = Form("pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_phi_Psi_%d_%s_PiK_SysError_%d",pt_QA,Cent_start,mCharge,Eta_start,i_pad,Order[mOrder].Data(),Sys_QA);
      h_mMass2[KEY_PiK_QA]->Draw("colz");
    }
    else if(i_pad == phi_total)
    {
      TString KEY_PiK_Total_QA = Form("pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_%s_PiK_SysError_%d",pt_QA,Cent_start,mCharge,Eta_start,Order[mOrder].Data(),Sys_QA);
      h_mMass2_total[KEY_PiK_Total_QA]->Draw("colz");
      if(pt_QA > 4)
      {
	Float_t y_start = y_p[pt_QA]-nSigProj[Proj_QA]*parfit_1st[KEY_PiK_Total_QA][17];
	PlotLine(x_low[pt_QA],x_up[pt_QA],y_start,y_start,2,2,2);
	Float_t y_stop  = y_pi[pt_QA]+nSigProj[Proj_QA]*parfit_1st[KEY_PiK_Total_QA][5];
	PlotLine(x_low[pt_QA],x_up[pt_QA],y_stop,y_stop,2,2,2);
      }
    }
    else
    {
      TString KEY_PiK_Tproj_QA = Form("pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_%s_PiK_SysError_%d_Proj_%d",pt_QA,Cent_start,mCharge,Eta_start,Order[mOrder].Data(),Sys_QA,Proj_QA);
      h_mXproj_total[KEY_PiK_Tproj_QA]->Draw("hE");
    }
  }
  // c_phi->SaveAs("c_phiQA.eps");
#endif

  // fit projction of pi, k integrated distribution on x-axis
  TF1Map f_1d_double_total, f_1d_pion_total, f_1d_kaon_total;
  vecFMap parfit;
  parfit.clear();
  for(Int_t i_pt = pt_rebin_first; i_pt < pt_rebin_last; i_pt++) // pt bin
  {
    for(Int_t i_cent = Cent_start; i_cent < Cent_stop; i_cent++) // centrality bin
    {
      for(Int_t i_charge = mCharge; i_charge < mCharge+1; i_charge++) // charge bin
      {
	for(Int_t i_eta = Eta_start; i_eta < Eta_stop; i_eta++) // eta gap bin
	{
	  for(Int_t i_sys = Sys_start; i_sys < Sys_stop; i_sys++)
	  {
	    TString KEY_PiK_Total = Form("pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_%s_PiK_SysError_%d",i_pt,i_cent,i_charge,i_eta,Order[mOrder].Data(),i_sys); // used for fit parameters
	    for(Int_t i_proj = Proj_start; i_proj < Proj_stop; ++i_proj)
	    {
	      if(i_pt <= 4) continue;

	      TString KEY_doule_total = Form("f_1d_double_pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_%s_PiK_SysError_%d_Proj_%d",i_pt,i_cent,i_charge,i_eta,Order[mOrder].Data(),i_sys,i_proj); // used for TF1s
	      f_1d_double_total[KEY_doule_total] = new TF1(KEY_doule_total.Data(),student_t_1d_double,x_low[i_pt],x_up[i_pt],8);

	      // pion
	      f_1d_double_total[KEY_doule_total]->SetParameter(0,parfit_2nd[KEY_PiK_Total][0]);
	      f_1d_double_total[KEY_doule_total]->SetParameter(1,parfit_2nd[KEY_PiK_Total][2]);
	      f_1d_double_total[KEY_doule_total]->SetParameter(2,parfit_2nd[KEY_PiK_Total][4]);
	      f_1d_double_total[KEY_doule_total]->SetParameter(3,1.0);

	      // kaon
	      f_1d_double_total[KEY_doule_total]->SetParameter(4,(2.0*parfit_2nd[KEY_PiK_Total][0]+parfit_2nd[KEY_PiK_Total][12])/3.0);
	      f_1d_double_total[KEY_doule_total]->SetParameter(5,parfit_2nd[KEY_PiK_Total][7]);
	      f_1d_double_total[KEY_doule_total]->SetParameter(6,parfit_2nd[KEY_PiK_Total][9]);
	      f_1d_double_total[KEY_doule_total]->SetParameter(7,1.0);

	      TString KEY_PiK_Tproj = Form("pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_%s_PiK_SysError_%d_Proj_%d",i_pt,i_cent,i_charge,i_eta,Order[mOrder].Data(),i_sys,i_proj); // used for TH1Fs and corresponding fit parameters
	      Float_t pi_max = h_mXproj_total[KEY_PiK_Tproj]->GetBinContent(h_mXproj_total[KEY_PiK_Tproj]->FindBin(parfit_2nd[KEY_PiK_Total][2]));
	      Float_t pi_func = f_1d_double_total[KEY_doule_total]->Eval(parfit_2nd[KEY_PiK_Total][2]);
	      Float_t Norm_pi = pi_max/pi_func;

	      Float_t k_max = h_mXproj_total[KEY_PiK_Tproj]->GetBinContent(h_mXproj_total[KEY_PiK_Tproj]->FindBin(parfit_2nd[KEY_PiK_Total][7]));
	      Float_t k_func = f_1d_double_total[KEY_doule_total]->Eval(parfit_2nd[KEY_PiK_Total][7]);
	      Float_t Norm_k = k_max/k_func;

	      f_1d_double_total[KEY_doule_total]->SetParameter(3,Norm_pi);
	      f_1d_double_total[KEY_doule_total]->SetParameter(7,Norm_k);
	      f_1d_double_total[KEY_doule_total]->SetRange(parfit_2nd[KEY_PiK_Total][2]-order_pion[i_pt]*parfit_2nd[KEY_PiK_Total][4],parfit_2nd[KEY_PiK_Total][7]+order_kaon[i_pt]*parfit_2nd[KEY_PiK_Total][9]);

	      // Fit
	      // cout << "i_pt = " << i_pt << ", charge = " << i_charge << ", eta_bin = " << i_eta <<  ", i_sys = " << i_sys << ", i_proj = " << i_proj << endl;
	      h_mXproj_total[KEY_PiK_Tproj]->Fit(f_1d_double_total[KEY_doule_total],"MQRN");

	      Float_t chi2 = f_1d_double_total[KEY_doule_total]->GetChisquare();
	      Float_t Ndf = f_1d_double_total[KEY_doule_total]->GetNDF();
	      // cout << "chi2/Ndf = " << chi2 << "/" << Ndf << " = " << chi2/Ndf << endl;
	      for(Int_t i_par = 0; i_par < 8; ++i_par)
	      {
		parfit[KEY_PiK_Tproj].push_back(static_cast<Float_t>(f_1d_double_total[KEY_doule_total]->GetParameter(i_par)));
	      }

	      f_1d_double_total[KEY_doule_total]->SetRange(x_low[i_pt],x_up[i_pt]);

	      TString KEY_pion_total = Form("f_1d_pion_pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_%s_PiK_SysError_%d_Proj_%d",i_pt,i_cent,i_charge,i_eta,Order[mOrder].Data(),i_sys,i_proj); // used for TF1s
	      f_1d_pion_total[KEY_pion_total] = new TF1(KEY_pion_total.Data(),student_t_1d_single,x_low[i_pt],x_up[i_pt],4);
	      f_1d_pion_total[KEY_pion_total]->FixParameter(0,f_1d_double_total[KEY_doule_total]->GetParameter(0));
	      f_1d_pion_total[KEY_pion_total]->FixParameter(1,f_1d_double_total[KEY_doule_total]->GetParameter(1));
	      f_1d_pion_total[KEY_pion_total]->FixParameter(2,f_1d_double_total[KEY_doule_total]->GetParameter(2));
	      f_1d_pion_total[KEY_pion_total]->FixParameter(3,f_1d_double_total[KEY_doule_total]->GetParameter(3));

	      TString KEY_kaon_total = Form("f_1d_kaon_pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_%s_PiK_SysError_%d_Proj_%d",i_pt,i_cent,i_charge,i_eta,Order[mOrder].Data(),i_sys,i_proj); // used for TF1s
	      f_1d_kaon_total[KEY_kaon_total] = new TF1(KEY_kaon_total.Data(),student_t_1d_single,x_low[i_pt],x_up[i_pt],4);
	      f_1d_kaon_total[KEY_kaon_total]->FixParameter(0,f_1d_double_total[KEY_doule_total]->GetParameter(4));
	      f_1d_kaon_total[KEY_kaon_total]->FixParameter(1,f_1d_double_total[KEY_doule_total]->GetParameter(5));
	      f_1d_kaon_total[KEY_kaon_total]->FixParameter(2,f_1d_double_total[KEY_doule_total]->GetParameter(6));
	      f_1d_kaon_total[KEY_kaon_total]->FixParameter(3,f_1d_double_total[KEY_doule_total]->GetParameter(7));
	    }
	  }
	}
      }
    }
  }

#if _PlotQA_
  TCanvas *c_pt_pik = new TCanvas("c_pt_pik","c_pt_pik",10,10,1200,1200);
  c_pt_pik->Divide(4,4);
  for(Int_t i_pt = pt_rebin_first; i_pt < pt_rebin_last; ++i_pt)
  {
    c_pt_pik->cd(i_pt+1);
    c_pt_pik->cd(i_pt+1)->SetLeftMargin(0.15);
    c_pt_pik->cd(i_pt+1)->SetBottomMargin(0.15);
    c_pt_pik->cd(i_pt+1)->SetTicks(1,1);
    c_pt_pik->cd(i_pt+1)->SetGrid(0,0);
    TString KEY_PiK_Tproj_QA = Form("pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_%s_PiK_SysError_%d_Proj_%d",i_pt,Cent_start,mCharge,Eta_start,Order[mOrder].Data(),Sys_QA,Proj_QA);
    h_mXproj_total[KEY_PiK_Tproj_QA]->SetStats(0);
    h_mXproj_total[KEY_PiK_Tproj_QA]->SetTitle(KEY_PiK_Tproj_QA.Data());

    h_mXproj_total[KEY_PiK_Tproj_QA]->GetXaxis()->SetNdivisions(505,'N');
    h_mXproj_total[KEY_PiK_Tproj_QA]->GetXaxis()->SetTitle("x (n#sigma_{#pi},m^{2})");
    h_mXproj_total[KEY_PiK_Tproj_QA]->GetXaxis()->SetTitleSize(0.05);
    h_mXproj_total[KEY_PiK_Tproj_QA]->GetXaxis()->CenterTitle();
    h_mXproj_total[KEY_PiK_Tproj_QA]->GetXaxis()->SetTitleOffset(0.9);

    h_mXproj_total[KEY_PiK_Tproj_QA]->GetYaxis()->SetNdivisions(505,'N');
    h_mXproj_total[KEY_PiK_Tproj_QA]->GetYaxis()->SetTitle("Counts/resolution");
    h_mXproj_total[KEY_PiK_Tproj_QA]->GetYaxis()->SetTitleSize(0.05);
    h_mXproj_total[KEY_PiK_Tproj_QA]->GetYaxis()->CenterTitle();
    h_mXproj_total[KEY_PiK_Tproj_QA]->GetYaxis()->SetTitleOffset(1.5);
    h_mXproj_total[KEY_PiK_Tproj_QA]->Draw("hE");
    TString pt_range = Form("[%1.1f,%1.1f]",pt_low[i_pt],pt_up[i_pt]);
    plotTopLegend((char*)pt_range.Data(),0.6,0.7,0.08,1,0.0,42,1);


    if(i_pt > 4)
    {
      TString KEY_doule_total_QA = Form("f_1d_double_pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_%s_PiK_SysError_%d_Proj_%d",i_pt,Cent_start,mCharge,Eta_start,Order[mOrder].Data(),Sys_QA,Proj_QA); // used for TF1s
      f_1d_double_total[KEY_doule_total_QA]->SetLineStyle(1);
      f_1d_double_total[KEY_doule_total_QA]->SetLineColor(2);
      f_1d_double_total[KEY_doule_total_QA]->SetFillStyle(0);
      f_1d_double_total[KEY_doule_total_QA]->Draw("l same");

      TString KEY_pion_total_QA = Form("f_1d_pion_pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_%s_PiK_SysError_%d_Proj_%d",i_pt,Cent_start,mCharge,Eta_start,Order[mOrder].Data(),Sys_QA,Proj_QA); // used for TF1s
      f_1d_pion_total[KEY_pion_total_QA]->SetFillColor(kAzure-2);
      f_1d_pion_total[KEY_pion_total_QA]->SetLineColor(kAzure-2);
      f_1d_pion_total[KEY_pion_total_QA]->SetLineStyle(1);
      f_1d_pion_total[KEY_pion_total_QA]->SetLineWidth(2);
      f_1d_pion_total[KEY_pion_total_QA]->Draw("h same");

      TString KEY_kaon_total_QA = Form("f_1d_kaon_pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_%s_PiK_SysError_%d_Proj_%d",i_pt,Cent_start,mCharge,Eta_start,Order[mOrder].Data(),Sys_QA,Proj_QA); // used for TF1s
      f_1d_kaon_total[KEY_kaon_total_QA]->SetFillColor(kGray+2);
      f_1d_kaon_total[KEY_kaon_total_QA]->SetLineColor(kGray+2);
      f_1d_kaon_total[KEY_kaon_total_QA]->SetLineStyle(1);
      f_1d_kaon_total[KEY_kaon_total_QA]->SetLineWidth(2);
      f_1d_kaon_total[KEY_kaon_total_QA]->SetFillStyle(3003);
      f_1d_kaon_total[KEY_kaon_total_QA]->Draw("h same");
    }
    else
    {
      Float_t Inte_pion_start = -0.05;
      Float_t Inte_pion_stop  = 0.05;
      Float_t Inte_kaon_start = 0.15;
      Float_t Inte_kaon_stop = 0.30;

      PlotLine(Inte_pion_start,Inte_pion_start,0,0.5*h_mXproj_total[KEY_PiK_Tproj_QA]->GetMaximum(),kAzure-2,2,2);
      PlotLine(Inte_pion_stop,Inte_pion_stop,0,0.5*h_mXproj_total[KEY_PiK_Tproj_QA]->GetMaximum(),kAzure-2,2,2);
      PlotLine(Inte_kaon_start,Inte_kaon_start,0,0.5*h_mXproj_total[KEY_PiK_Tproj_QA]->GetMaximum(),kGray+2,2,2);
      PlotLine(Inte_kaon_stop,Inte_kaon_stop,0,0.5*h_mXproj_total[KEY_PiK_Tproj_QA]->GetMaximum(),kGray+2,2,2);
    }
  }
#endif

  // projected phi-differential distritbuion to x-axis
  TH1FMap h_mXproj, h_mPion, h_mKaon; 
  TF1Map f_1d_double, f_1d_pion, f_1d_kaon;
  for(Int_t i_pt = pt_rebin_first; i_pt < pt_rebin_last; i_pt++) // pt bin
  {
    for(Int_t i_cent = Cent_start; i_cent < Cent_stop; i_cent++) // centrality bin
    {
      for(Int_t i_charge = mCharge; i_charge < mCharge+1; i_charge++) // charge bin
      {
	for(Int_t i_eta = Eta_start; i_eta < Eta_stop; i_eta++) // eta gap bin
	{
	  for(Int_t i_sys = Sys_start; i_sys < Sys_stop; i_sys++)
	  {
	    TString KEY_PiK_Total = Form("pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_%s_PiK_SysError_%d",i_pt,i_cent,i_charge,i_eta,Order[mOrder].Data(),i_sys);
	    for(Int_t i_phi = phi_start; i_phi < phi_stop; i_phi++) // phi-psi bin
	    {
	      TString KEY_PiK = Form("pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_phi_Psi_%d_%s_PiK_SysError_%d",i_pt,i_cent,i_charge,i_eta,i_phi,Order[mOrder].Data(),i_sys);
	      for(Int_t i_proj = Proj_start; i_proj < Proj_stop; ++i_proj)
	      { // projction to x-axis
		TString KEY_PiK_proj = Form("pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_phi_Psi_%d_%s_PiK_SysError_%d_Proj_%d",i_pt,i_cent,i_charge,i_eta,i_phi,Order[mOrder].Data(),i_sys,i_proj);
		TString KEY_YieldPion = Form("pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_phi_Psi_%d_%s_Pion_SysError_%d_Proj_%d",i_pt,i_cent,i_charge,i_eta,i_phi,Order[mOrder].Data(),i_sys,i_proj);
		TString KEY_YieldKaon = Form("pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_phi_Psi_%d_%s_Kaon_SysError_%d_Proj_%d",i_pt,i_cent,i_charge,i_eta,i_phi,Order[mOrder].Data(),i_sys,i_proj);
		if(i_pt > 4)
		{
		  Int_t bin_start = h_mMass2[KEY_PiK]->GetYaxis()->FindBin(y_p[i_pt]-nSigProj[i_proj]*parfit_1st[KEY_PiK_Total][17]);
		  Int_t bin_stop  = h_mMass2[KEY_PiK]->GetYaxis()->FindBin(y_pi[i_pt]+nSigProj[i_proj]*parfit_1st[KEY_PiK_Total][5]);
		  h_mXproj[KEY_PiK_proj] = (TH1F*)h_mMass2[KEY_PiK]->ProjectionX(KEY_PiK_proj.Data(),bin_start,bin_stop);
		  h_mPion[KEY_YieldPion] = (TH1F*)h_mMass2[KEY_PiK]->ProjectionX(KEY_YieldPion.Data(),bin_start,bin_stop);
		  h_mKaon[KEY_YieldKaon] = (TH1F*)h_mMass2[KEY_PiK]->ProjectionX(KEY_YieldKaon.Data(),bin_start,bin_stop);

		  TString KEY_doule = Form("f_1d_double_pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_phi_Psi_%d_%s_PiK_SysError_%d_Proj_%d",i_pt,i_cent,i_charge,i_eta,i_phi,Order[mOrder].Data(),i_sys,i_proj); // used for TF1s
		  f_1d_double[KEY_doule] = new TF1(KEY_doule.Data(),student_t_1d_double,x_low[i_pt],x_up[i_pt],8);

		  TString KEY_PiK_Tproj = Form("pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_%s_PiK_SysError_%d_Proj_%d",i_pt,i_cent,i_charge,i_eta,Order[mOrder].Data(),i_sys,i_proj); // used for TH1Fs and corresponding fit parameters
		  // pion
		  f_1d_double[KEY_doule]->FixParameter(0,parfit[KEY_PiK_Tproj][0]);
		  f_1d_double[KEY_doule]->FixParameter(1,parfit[KEY_PiK_Tproj][1]);
		  f_1d_double[KEY_doule]->FixParameter(2,parfit[KEY_PiK_Tproj][2]);
		  f_1d_double[KEY_doule]->SetParameter(3,parfit[KEY_PiK_Tproj][3]);

		  // kaon
		  f_1d_double[KEY_doule]->FixParameter(4,parfit[KEY_PiK_Tproj][4]);
		  f_1d_double[KEY_doule]->FixParameter(5,parfit[KEY_PiK_Tproj][5]);
		  f_1d_double[KEY_doule]->FixParameter(6,parfit[KEY_PiK_Tproj][6]);
		  f_1d_double[KEY_doule]->SetParameter(7,parfit[KEY_PiK_Tproj][7]);

		  Float_t pi_max = h_mXproj[KEY_PiK_proj]->GetBinContent(h_mXproj[KEY_PiK_proj]->FindBin(parfit[KEY_PiK_Tproj][1]));
		  Float_t pi_func = f_1d_double[KEY_doule]->Eval(parfit[KEY_PiK_Tproj][1]);
		  Float_t Norm_pi = pi_max/pi_func;

		  Float_t k_max = h_mXproj[KEY_PiK_proj]->GetBinContent(h_mXproj[KEY_PiK_proj]->FindBin(parfit[KEY_PiK_Tproj][5]));
		  Float_t k_func = f_1d_double[KEY_doule]->Eval(parfit[KEY_PiK_Tproj][5]);
		  Float_t Norm_k = k_max/k_func;

		  f_1d_double[KEY_doule]->SetParameter(3,Norm_pi);
		  f_1d_double[KEY_doule]->SetParameter(7,Norm_k);

		  f_1d_double[KEY_doule]->SetRange(parfit_2nd[KEY_PiK_Total][2]-order_pion[i_pt]*parfit_2nd[KEY_PiK_Total][4],parfit_2nd[KEY_PiK_Total][7]+order_kaon[i_pt]*parfit_2nd[KEY_PiK_Total][9]);

		  // Fit
		  // cout << "i_pt = " << i_pt << ", charge = " << i_charge << ", eta_bin = " << i_eta << ", i_phi = " << i_phi << ", i_sys = " << i_sys << ", i_proj = " << i_proj << endl;
		  h_mXproj[KEY_PiK_proj]->Fit(f_1d_double[KEY_doule],"MQRN");

		  Float_t chi2 = f_1d_double[KEY_doule]->GetChisquare();
		  Float_t Ndf = f_1d_double[KEY_doule]->GetNDF();
		  // cout << "chi2/Ndf = " << chi2 << "/" << Ndf << " = " << chi2/Ndf << endl;

		  f_1d_double[KEY_doule]->SetRange(x_low[i_pt],x_up[i_pt]);

		  TString KEY_pion = Form("f_1d_pion_pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_phi_Psi_%d_%s_PiK_SysError_%d_Proj_%d",i_pt,i_cent,i_charge,i_eta,i_phi,Order[mOrder].Data(),i_sys,i_proj); // used for TF1s
		  f_1d_pion[KEY_pion] = new TF1(KEY_pion.Data(),student_t_1d_single,x_low[i_pt],x_up[i_pt],4);
		  f_1d_pion[KEY_pion]->FixParameter(0,f_1d_double[KEY_doule]->GetParameter(0));
		  f_1d_pion[KEY_pion]->FixParameter(1,f_1d_double[KEY_doule]->GetParameter(1));
		  f_1d_pion[KEY_pion]->FixParameter(2,f_1d_double[KEY_doule]->GetParameter(2));
		  f_1d_pion[KEY_pion]->FixParameter(3,f_1d_double[KEY_doule]->GetParameter(3));

		  TString KEY_kaon = Form("f_1d_kaon_pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_phi_Psi_%d_%s_PiK_SysError_%d_Proj_%d",i_pt,i_cent,i_charge,i_eta,i_phi,Order[mOrder].Data(),i_sys,i_proj); // used for TF1s
		  f_1d_kaon[KEY_kaon] = new TF1(KEY_kaon.Data(),student_t_1d_single,x_low[i_pt],x_up[i_pt],4);
		  f_1d_kaon[KEY_kaon]->FixParameter(0,f_1d_double[KEY_doule]->GetParameter(4));
		  f_1d_kaon[KEY_kaon]->FixParameter(1,f_1d_double[KEY_doule]->GetParameter(5));
		  f_1d_kaon[KEY_kaon]->FixParameter(2,f_1d_double[KEY_doule]->GetParameter(6));
		  f_1d_kaon[KEY_kaon]->FixParameter(3,f_1d_double[KEY_doule]->GetParameter(7));
		}
		else
		{
		  Int_t bin_start = h_mMass2[KEY_PiK]->GetYaxis()->FindBin(y_low[i_pt]);
		  Int_t bin_stop  = h_mMass2[KEY_PiK]->GetYaxis()->FindBin(y_up[i_pt]);
		  h_mXproj[KEY_PiK_proj] = (TH1F*)h_mMass2[KEY_PiK]->ProjectionX(KEY_PiK_proj.Data(),bin_start,bin_stop);
		  h_mPion[KEY_YieldPion] = (TH1F*)h_mMass2[KEY_PiK]->ProjectionX(KEY_YieldPion.Data(),bin_start,bin_stop);
		  h_mKaon[KEY_YieldKaon] = (TH1F*)h_mMass2[KEY_PiK]->ProjectionX(KEY_YieldKaon.Data(),bin_start,bin_stop);
		}
	      }
	    }
	  }
	}
      }
    }
  }

#if _PlotQA_
  TCanvas *c_phi_proj = new TCanvas("c_phi_proj","c_phi_proj",1200,1200);
  c_phi_proj->Divide(3,3);
  for(Int_t i_phi = phi_start; i_phi < phi_stop; ++i_phi)
  {
    c_phi_proj->cd(i_phi+1);
    c_phi_proj->cd(i_phi+1)->SetLeftMargin(0.15);
    c_phi_proj->cd(i_phi+1)->SetBottomMargin(0.15);
    c_phi_proj->cd(i_phi+1)->SetTicks(1,1);
    c_phi_proj->cd(i_phi+1)->SetGrid(0,0);
    TString KEY_PiK_proj_QA = Form("pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_phi_Psi_%d_%s_PiK_SysError_%d_Proj_%d",pt_QA,Cent_start,mCharge,Eta_start,i_phi,Order[mOrder].Data(),Sys_QA,Proj_QA);
    h_mXproj[KEY_PiK_proj_QA]->SetStats(0);
    h_mXproj[KEY_PiK_proj_QA]->SetTitle(KEY_PiK_proj_QA.Data());

    h_mXproj[KEY_PiK_proj_QA]->GetXaxis()->SetNdivisions(505,'N');
    h_mXproj[KEY_PiK_proj_QA]->GetXaxis()->SetTitle("x (n#sigma_{#pi},m^{2})");
    h_mXproj[KEY_PiK_proj_QA]->GetXaxis()->SetTitleSize(0.05);
    h_mXproj[KEY_PiK_proj_QA]->GetXaxis()->CenterTitle();
    h_mXproj[KEY_PiK_proj_QA]->GetXaxis()->SetTitleOffset(0.9);

    h_mXproj[KEY_PiK_proj_QA]->GetYaxis()->SetNdivisions(505,'N');
    h_mXproj[KEY_PiK_proj_QA]->GetYaxis()->SetTitle("Counts/resolution");
    h_mXproj[KEY_PiK_proj_QA]->GetYaxis()->SetTitleSize(0.05);
    h_mXproj[KEY_PiK_proj_QA]->GetYaxis()->CenterTitle();
    h_mXproj[KEY_PiK_proj_QA]->GetYaxis()->SetTitleOffset(1.5);
    h_mXproj[KEY_PiK_proj_QA]->Draw("hE");
    if(pt_QA > 4)
    {
      TString KEY_doule_QA = Form("f_1d_double_pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_phi_Psi_%d_%s_PiK_SysError_%d_Proj_%d",pt_QA,Cent_start,mCharge,Eta_start,i_phi,Order[mOrder].Data(),Sys_start,Proj_QA); // used for TF1s
      f_1d_double[KEY_doule_QA]->SetLineStyle(1);
      f_1d_double[KEY_doule_QA]->SetLineColor(2);
      f_1d_double[KEY_doule_QA]->SetFillStyle(0);
      f_1d_double[KEY_doule_QA]->Draw("l same");

      TString KEY_pion_QA = Form("f_1d_pion_pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_phi_Psi_%d_%s_PiK_SysError_%d_Proj_%d",pt_QA,Cent_start,mCharge,Eta_start,i_phi,Order[mOrder].Data(),Sys_start,Proj_QA); // used for TF1s
      f_1d_pion[KEY_pion_QA]->SetFillColor(kAzure-2);
      f_1d_pion[KEY_pion_QA]->SetLineColor(kAzure-2);
      f_1d_pion[KEY_pion_QA]->SetLineStyle(1);
      f_1d_pion[KEY_pion_QA]->SetLineWidth(2);
      f_1d_pion[KEY_pion_QA]->Draw("h same");

      TString KEY_kaon_QA = Form("f_1d_kaon_pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_phi_Psi_%d_%s_PiK_SysError_%d_Proj_%d",pt_QA,Cent_start,mCharge,Eta_start,i_phi,Order[mOrder].Data(),Sys_start,Proj_QA); // used for TF1s
      f_1d_kaon[KEY_kaon_QA]->SetFillColor(kGray+2);
      f_1d_kaon[KEY_kaon_QA]->SetLineColor(kGray+2);
      f_1d_kaon[KEY_kaon_QA]->SetLineStyle(1);
      f_1d_kaon[KEY_kaon_QA]->SetLineWidth(2);
      f_1d_kaon[KEY_kaon_QA]->SetFillStyle(3003);
      f_1d_kaon[KEY_kaon_QA]->Draw("h same");
    }
    else
    {
      Float_t Inte_pion_start = -0.05;
      Float_t Inte_pion_stop  = 0.05;
      Float_t Inte_kaon_start = 0.15;
      Float_t Inte_kaon_stop = 0.30;

      PlotLine(Inte_pion_start,Inte_pion_start,0,0.5*h_mXproj[KEY_PiK_proj_QA]->GetMaximum(),kAzure-2,2,2);
      PlotLine(Inte_pion_stop,Inte_pion_stop,0,0.5*h_mXproj[KEY_PiK_proj_QA]->GetMaximum(),kAzure-2,2,2);
      PlotLine(Inte_kaon_start,Inte_kaon_start,0,0.5*h_mXproj[KEY_PiK_proj_QA]->GetMaximum(),kGray+2,2,2);
      PlotLine(Inte_kaon_stop,Inte_kaon_stop,0,0.5*h_mXproj[KEY_PiK_proj_QA]->GetMaximum(),kGray+2,2,2);
    }
  }
#endif

  // get yield of pions and kaons
  TH1FMap h_mCounts_pion, h_mCounts_kaon;
  for(Int_t i_pt = pt_rebin_first; i_pt < pt_rebin_last; i_pt++) // pt bin
  {
    for(Int_t i_cent = Cent_start; i_cent < Cent_stop; i_cent++) // centrality bin
    {
      for(Int_t i_charge = mCharge; i_charge < mCharge+1; i_charge++) // charge bin
      {
	for(Int_t i_eta = Eta_start; i_eta < Eta_stop; i_eta++) // eta gap bin
	{
	  for(Int_t i_sys = Sys_start; i_sys < Sys_stop; i_sys++)
	  {
	    for(Int_t i_proj = Proj_start; i_proj < Proj_stop; ++i_proj)
	    { 
	      TString KEY_Counts_Pion = Form("Pion_Counts_pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_%s_SysError_%d_Proj_%d",i_pt,i_cent,i_charge,i_eta,Order[mOrder].Data(),i_sys,i_proj); 
	      h_mCounts_pion[KEY_Counts_Pion] = new TH1F(KEY_Counts_Pion.Data(),KEY_Counts_Pion.Data(),7,0,PI_max[mOrder]);
	      TString KEY_Inte_Pion = Form("Pion_Inte_pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_%s_SysError_%d_Proj_%d",i_pt,i_cent,i_charge,i_eta,Order[mOrder].Data(),i_sys,i_proj); 
	      h_mCounts_pion[KEY_Inte_Pion] = new TH1F(KEY_Inte_Pion.Data(),KEY_Inte_Pion.Data(),7,0,PI_max[mOrder]);

	      TString KEY_Counts_Kaon = Form("Kaon_Counts_pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_%s_SysError_%d_Proj_%d",i_pt,i_cent,i_charge,i_eta,Order[mOrder].Data(),i_sys,i_proj); 
	      h_mCounts_kaon[KEY_Counts_Kaon] = new TH1F(KEY_Counts_Kaon.Data(),KEY_Counts_Kaon.Data(),7,0,PI_max[mOrder]);
	      TString KEY_Inte_Kaon = Form("Kaon_Inte_pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_%s_SysError_%d_Proj_%d",i_pt,i_cent,i_charge,i_eta,Order[mOrder].Data(),i_sys,i_proj); 
	      h_mCounts_kaon[KEY_Inte_Kaon] = new TH1F(KEY_Inte_Kaon.Data(),KEY_Inte_Kaon.Data(),7,0,PI_max[mOrder]);

	      TString KEY_PiK_Tproj = Form("pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_%s_PiK_SysError_%d_Proj_%d",i_pt,i_cent,i_charge,i_eta,Order[mOrder].Data(),i_sys,i_proj); // used for TH1Fs and corresponding fit parameters
	      for(Int_t i_phi = phi_start; i_phi < phi_stop; i_phi++) // phi-psi bin
	      {
		TString KEY_YieldPion = Form("pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_phi_Psi_%d_%s_Pion_SysError_%d_Proj_%d",i_pt,i_cent,i_charge,i_eta,i_phi,Order[mOrder].Data(),i_sys,i_proj);
		TString KEY_pion = Form("f_1d_pion_pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_phi_Psi_%d_%s_PiK_SysError_%d_Proj_%d",i_pt,i_cent,i_charge,i_eta,i_phi,Order[mOrder].Data(),i_sys,i_proj); // used for TF1s
		TString KEY_YieldKaon = Form("pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_phi_Psi_%d_%s_Kaon_SysError_%d_Proj_%d",i_pt,i_cent,i_charge,i_eta,i_phi,Order[mOrder].Data(),i_sys,i_proj);
		TString KEY_kaon = Form("f_1d_kaon_pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_phi_Psi_%d_%s_PiK_SysError_%d_Proj_%d",i_pt,i_cent,i_charge,i_eta,i_phi,Order[mOrder].Data(),i_sys,i_proj); // used for TF1s
		if(i_pt > 4)
		{
		  h_mPion[KEY_YieldPion]->Add(f_1d_kaon[KEY_kaon],-1.0);
		  Float_t Inte_pion_start = parfit[KEY_PiK_Tproj][1] - nSigma*parfit[KEY_PiK_Tproj][2];
		  Float_t Inte_pion_stop  = parfit[KEY_PiK_Tproj][1] + nSigma*parfit[KEY_PiK_Tproj][2];
		  Float_t bin_width_pion = h_mPion[KEY_YieldPion]->GetBinWidth(1);
		  TF1 *f_pion = new TF1(KEY_pion.Data(),student_t_1d_single,x_low[i_pt],x_up[i_pt],4);
		  for(Int_t i_par = 0; i_par < 4; ++i_par)
		  {
		    f_pion->ReleaseParameter(i_par);
		  }
		  f_pion->FixParameter(0,f_1d_pion[KEY_pion]->GetParameter(0));
		  f_pion->FixParameter(1,f_1d_pion[KEY_pion]->GetParameter(1));
		  f_pion->FixParameter(2,f_1d_pion[KEY_pion]->GetParameter(2));
		  f_pion->SetParameter(3,f_1d_pion[KEY_pion]->GetParameter(3));
		  h_mPion[KEY_YieldPion]->Fit(f_pion,"NR");
		  Float_t counts_pion_Inte = f_pion->Integral(Inte_pion_start,Inte_pion_stop)/bin_width_pion;
		  Float_t err_pion_Inte = f_pion->IntegralError(Inte_pion_start,Inte_pion_stop)/bin_width_pion;
		  h_mCounts_pion[KEY_Inte_Pion]->SetBinContent(i_phi+1,counts_pion_Inte);
		  h_mCounts_pion[KEY_Inte_Pion]->SetBinError(i_phi+1,err_pion_Inte);

		  Int_t Bin_pion_start = h_mPion[KEY_YieldPion]->FindBin(Inte_pion_start);
		  Int_t Bin_pion_stop = h_mPion[KEY_YieldPion]->FindBin(Inte_pion_stop);
		  Float_t counts_pion = 0.0;
		  Float_t err_pion = 0.0;
		  for(Int_t i_pi = Bin_pion_start; i_pi < Bin_pion_stop; ++i_pi)
		  {
		    counts_pion += h_mPion[KEY_YieldPion]->GetBinContent(i_pi);
		    err_pion += h_mPion[KEY_YieldPion]->GetBinError(i_pi)*h_mPion[KEY_YieldPion]->GetBinError(i_pi);
		  }
		  h_mCounts_pion[KEY_Counts_Pion]->SetBinContent(i_phi+1,counts_pion);
		  h_mCounts_pion[KEY_Counts_Pion]->SetBinError(i_phi+1,TMath::Sqrt(err_pion));

		  h_mKaon[KEY_YieldKaon]->Add(f_1d_pion[KEY_pion],-1.0);
		  Float_t Inte_kaon_start = parfit[KEY_PiK_Tproj][5] - nSigma*parfit[KEY_PiK_Tproj][6];
		  Float_t Inte_kaon_stop  = parfit[KEY_PiK_Tproj][5] + nSigma*parfit[KEY_PiK_Tproj][6];
		  Float_t bin_width_kaon = h_mKaon[KEY_YieldKaon]->GetBinWidth(1);
		  TF1 *f_kaon = new TF1(KEY_kaon.Data(),student_t_1d_single,x_low[i_pt],x_up[i_pt],4);
		  for(Int_t i_par = 0; i_par < 4; ++i_par)
		  {
		    f_kaon->ReleaseParameter(i_par);
		  }
		  f_kaon->FixParameter(0,f_1d_kaon[KEY_kaon]->GetParameter(0));
		  f_kaon->FixParameter(1,f_1d_kaon[KEY_kaon]->GetParameter(1));
		  f_kaon->FixParameter(2,f_1d_kaon[KEY_kaon]->GetParameter(2));
		  f_kaon->SetParameter(3,f_1d_kaon[KEY_kaon]->GetParameter(3));
		  h_mKaon[KEY_YieldKaon]->Fit(f_kaon,"NR");
		  Float_t counts_kaon_Inte = f_kaon->Integral(Inte_kaon_start,Inte_kaon_stop)/bin_width_kaon;
		  Float_t err_kaon_Inte = f_kaon->IntegralError(Inte_kaon_start,Inte_kaon_stop)/bin_width_kaon;
		  h_mCounts_kaon[KEY_Inte_Kaon]->SetBinContent(i_phi+1,counts_kaon_Inte);
		  h_mCounts_kaon[KEY_Inte_Kaon]->SetBinError(i_phi+1,err_kaon_Inte);

		  Int_t Bin_kaon_start = h_mKaon[KEY_YieldKaon]->FindBin(Inte_kaon_start);
		  Int_t Bin_kaon_stop = h_mKaon[KEY_YieldKaon]->FindBin(Inte_kaon_stop);
		  Float_t counts_kaon = 0.0;
		  Float_t err_kaon = 0.0;
		  for(Int_t i_k = Bin_kaon_start; i_k < Bin_kaon_stop; ++i_k)
		  {
		    counts_kaon += h_mKaon[KEY_YieldKaon]->GetBinContent(i_k);
		    err_kaon += h_mKaon[KEY_YieldKaon]->GetBinError(i_k)*h_mKaon[KEY_YieldKaon]->GetBinError(i_k);
		  }
		  h_mCounts_kaon[KEY_Counts_Kaon]->SetBinContent(i_phi+1,counts_kaon);
		  h_mCounts_kaon[KEY_Counts_Kaon]->SetBinError(i_phi+1,TMath::Sqrt(err_kaon));
		}
		else
		{
		  Float_t Inte_pion_start = -0.05;
		  Float_t Inte_pion_stop  = 0.05;
		  Int_t Bin_pion_start = h_mPion[KEY_YieldPion]->FindBin(Inte_pion_start);
		  Int_t Bin_pion_stop = h_mPion[KEY_YieldPion]->FindBin(Inte_pion_stop);
		  Float_t counts_pion = 0.0;
		  Float_t err_pion = 0.0;
		  for(Int_t i_pi = Bin_pion_start; i_pi < Bin_pion_stop; ++i_pi)
		  {
		    counts_pion += h_mPion[KEY_YieldPion]->GetBinContent(i_pi);
		    err_pion += h_mPion[KEY_YieldPion]->GetBinError(i_pi)*h_mPion[KEY_YieldPion]->GetBinError(i_pi);
		  }
		  h_mCounts_pion[KEY_Counts_Pion]->SetBinContent(i_phi+1,counts_pion);
		  h_mCounts_pion[KEY_Counts_Pion]->SetBinError(i_phi+1,TMath::Sqrt(err_pion));

		  Float_t Inte_kaon_start = 0.15;
		  Float_t Inte_kaon_stop = 0.30;
		  Int_t Bin_kaon_start = h_mKaon[KEY_YieldKaon]->FindBin(Inte_kaon_start);
		  Int_t Bin_kaon_stop = h_mKaon[KEY_YieldKaon]->FindBin(Inte_kaon_stop);
		  Float_t counts_kaon = 0.0;
		  Float_t err_kaon = 0.0;
		  for(Int_t i_k = Bin_kaon_start; i_k < Bin_kaon_stop; ++i_k)
		  {
		    counts_kaon += h_mKaon[KEY_YieldKaon]->GetBinContent(i_k);
		    err_kaon += h_mKaon[KEY_YieldKaon]->GetBinError(i_k)*h_mKaon[KEY_YieldKaon]->GetBinError(i_k);
		  }
		  h_mCounts_kaon[KEY_Counts_Kaon]->SetBinContent(i_phi+1,counts_kaon);
		  h_mCounts_kaon[KEY_Counts_Kaon]->SetBinError(i_phi+1,TMath::Sqrt(err_kaon));
		}
	      }
	    }
	  }
	}
      }
    }
  }

#if _PlotQA_
  TCanvas *c_Yield_Pion = new TCanvas("c_Yield_Pion","c_Yield_Pion",1200,1200);
  c_Yield_Pion->Divide(3,3);
  for(Int_t i_phi = 0; i_phi < 9; ++i_phi)
  {
    c_Yield_Pion->cd(i_phi+1);
    c_Yield_Pion->cd(i_phi+1)->SetLeftMargin(0.15);
    c_Yield_Pion->cd(i_phi+1)->SetBottomMargin(0.15);
    c_Yield_Pion->cd(i_phi+1)->SetTicks(1,1);
    c_Yield_Pion->cd(i_phi+1)->SetGrid(0,0);
    if(i_phi < phi_stop)
    {
      TString KEY_YieldPion_QA = Form("pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_phi_Psi_%d_%s_Pion_SysError_%d_Proj_%d",pt_QA,Cent_start,mCharge,Eta_start,i_phi,Order[mOrder].Data(),Sys_QA,Proj_QA);
      h_mPion[KEY_YieldPion_QA]->Draw("hE");

      TString KEY_pion_QA = Form("f_1d_pion_pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_phi_Psi_%d_%s_PiK_SysError_%d_Proj_%d",pt_QA,Cent_start,mCharge,Eta_start,i_phi,Order[mOrder].Data(),Sys_start,Proj_start);
      f_1d_pion[KEY_pion_QA]->Draw("l same");

      TString KEY_PiK_Tproj_QA = Form("pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_%s_PiK_SysError_%d_Proj_%d",pt_QA,Cent_start,mCharge,Eta_start,Order[mOrder].Data(),Sys_QA,Proj_QA);
      Float_t Inte_pion_start = parfit[KEY_PiK_Tproj_QA][1] - nSigma*parfit[KEY_PiK_Tproj_QA][2];
      Float_t Inte_pion_stop  = parfit[KEY_PiK_Tproj_QA][1] + nSigma*parfit[KEY_PiK_Tproj_QA][2];
      PlotLine(Inte_pion_start,Inte_pion_start,0,0.8*h_mPion[KEY_YieldPion_QA]->GetMaximum(),kAzure-2,2,2);
      PlotLine(Inte_pion_stop,Inte_pion_stop,0,0.8*h_mPion[KEY_YieldPion_QA]->GetMaximum(),kAzure-2,2,2);
    }
    if(i_phi == phi_stop)
    {
      TString KEY_Counts_Pion_QA = Form("Pion_Counts_pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_%s_SysError_%d_Proj_%d",pt_QA,Cent_start,mCharge,Eta_start,Order[mOrder].Data(),Sys_QA,Proj_QA); 
      h_mCounts_pion[KEY_Counts_Pion_QA]->SetMarkerStyle(24);
      h_mCounts_pion[KEY_Counts_Pion_QA]->SetMarkerColor(2);
      h_mCounts_pion[KEY_Counts_Pion_QA]->SetMarkerSize(1.4);
      h_mCounts_pion[KEY_Counts_Pion_QA]->Draw("pE");
      TString KEY_Inte_Pion_QA = Form("Pion_Inte_pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_%s_SysError_%d_Proj_%d",pt_QA,Cent_start,mCharge,Eta_start,Order[mOrder].Data(),Sys_QA,Proj_QA); 
      h_mCounts_pion[KEY_Inte_Pion_QA]->SetMarkerStyle(20);
      h_mCounts_pion[KEY_Inte_Pion_QA]->SetMarkerColor(kGray+2);
      h_mCounts_pion[KEY_Inte_Pion_QA]->SetMarkerSize(1.4);
      h_mCounts_pion[KEY_Inte_Pion_QA]->Draw("pE same");
    }
  }

  TCanvas *c_Yield_Kaon = new TCanvas("c_Yield_Kaon","c_Yield_Kaon",1200,1200);
  c_Yield_Kaon->Divide(3,3);
  for(Int_t i_phi = 0; i_phi < 9; ++i_phi)
  {
    c_Yield_Kaon->cd(i_phi+1);
    c_Yield_Kaon->cd(i_phi+1)->SetLeftMargin(0.15);
    c_Yield_Kaon->cd(i_phi+1)->SetBottomMargin(0.15);
    c_Yield_Kaon->cd(i_phi+1)->SetTicks(1,1);
    c_Yield_Kaon->cd(i_phi+1)->SetGrid(0,0);
    if(i_phi < phi_stop)
    {
      TString KEY_YieldKaon_QA = Form("pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_phi_Psi_%d_%s_Kaon_SysError_%d_Proj_%d",pt_QA,Cent_start,mCharge,Eta_start,i_phi,Order[mOrder].Data(),Sys_QA,Proj_QA);
      h_mKaon[KEY_YieldKaon_QA]->Draw("hE");

      TString KEY_kaon_QA = Form("f_1d_kaon_pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_phi_Psi_%d_%s_PiK_SysError_%d_Proj_%d",pt_QA,Cent_start,mCharge,Eta_start,i_phi,Order[mOrder].Data(),Sys_QA,Proj_QA);
      f_1d_kaon[KEY_kaon_QA]->Draw("l same");

      TString KEY_PiK_Tproj_QA = Form("pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_%s_PiK_SysError_%d_Proj_%d",pt_QA,Cent_start,mCharge,Eta_start,Order[mOrder].Data(),Sys_QA,Proj_QA);
      Float_t Inte_kaon_start = parfit[KEY_PiK_Tproj_QA][5] - nSigma*parfit[KEY_PiK_Tproj_QA][6];
      Float_t Inte_kaon_stop  = parfit[KEY_PiK_Tproj_QA][5] + nSigma*parfit[KEY_PiK_Tproj_QA][6];
      PlotLine(Inte_kaon_start,Inte_kaon_start,0,0.8*h_mKaon[KEY_YieldKaon_QA]->GetMaximum(),kGray+2,2,2);
      PlotLine(Inte_kaon_stop,Inte_kaon_stop,0,0.8*h_mKaon[KEY_YieldKaon_QA]->GetMaximum(),kGray+2,2,2);
    }
    if(i_phi == phi_stop)
    {
      TString KEY_Counts_Kaon_QA = Form("Kaon_Counts_pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_%s_SysError_%d_Proj_%d",pt_QA,Cent_start,mCharge,Eta_start,Order[mOrder].Data(),Sys_QA,Proj_QA); 
      h_mCounts_kaon[KEY_Counts_Kaon_QA]->SetMarkerStyle(24);
      h_mCounts_kaon[KEY_Counts_Kaon_QA]->SetMarkerColor(2);
      h_mCounts_kaon[KEY_Counts_Kaon_QA]->SetMarkerSize(1.4);
      h_mCounts_kaon[KEY_Counts_Kaon_QA]->Draw("pE");
      TString KEY_Inte_Kaon_QA = Form("Kaon_Inte_pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_%s_SysError_%d_Proj_%d",pt_QA,Cent_start,mCharge,Eta_start,Order[mOrder].Data(),Sys_QA,Proj_QA); 
      h_mCounts_kaon[KEY_Inte_Kaon_QA]->SetMarkerStyle(20);
      h_mCounts_kaon[KEY_Inte_Kaon_QA]->SetMarkerColor(kGray+2);
      h_mCounts_kaon[KEY_Inte_Kaon_QA]->SetMarkerSize(1.4);
      h_mCounts_kaon[KEY_Inte_Kaon_QA]->Draw("pE same");
    }
  }
#endif

  TString outputPion = Form("./OutPut/AuAu%s/Pion/Counts_%s_Charge_%d.root",Energy[mEnergy].Data(),Order[mOrder].Data(),mCharge);
  TFile *File_Pion = new TFile(outputPion.Data(),"RECREATE");
  File_Pion->cd();
  for(Int_t i_pt = pt_rebin_first; i_pt < pt_rebin_last; i_pt++) // pt bin
  {
    for(Int_t i_cent = Cent_start; i_cent < Cent_stop; i_cent++) // centrality bin
    {
      for(Int_t i_charge = mCharge; i_charge < mCharge+1; i_charge++) // charge bin
      {
	for(Int_t i_eta = Eta_start; i_eta < Eta_stop; i_eta++) // eta gap bin
	{
	  for(Int_t i_sys = Sys_start; i_sys < Sys_stop; i_sys++)
	  {
	    for(Int_t i_proj = Proj_start; i_proj < Proj_stop; ++i_proj)
	    { 
	      TString KEY_Counts_Pion = Form("Pion_Counts_pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_%s_SysError_%d_Proj_%d",i_pt,i_cent,i_charge,i_eta,Order[mOrder].Data(),i_sys,i_proj); 
	      h_mCounts_pion[KEY_Counts_Pion]->Write();
	      TString KEY_Inte_Pion = Form("Pion_Inte_pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_%s_SysError_%d_Proj_%d",i_pt,i_cent,i_charge,i_eta,Order[mOrder].Data(),i_sys,i_proj); 
	      h_mCounts_pion[KEY_Inte_Pion]->Write();
	    }
	  }
	}
      }
    }
  }
  File_Pion->Close();

  TString outputKaon = Form("./OutPut/AuAu%s/Kaon/Counts_%s_Charge_%d.root",Energy[mEnergy].Data(),Order[mOrder].Data(),mCharge);
  TFile *File_Kaon = new TFile(outputKaon.Data(),"RECREATE");
  File_Kaon->cd();
  for(Int_t i_pt = pt_rebin_first; i_pt < pt_rebin_last; i_pt++) // pt bin
  {
    for(Int_t i_cent = Cent_start; i_cent < Cent_stop; i_cent++) // centrality bin
    {
      for(Int_t i_charge = mCharge; i_charge < mCharge+1; i_charge++) // charge bin
      {
	for(Int_t i_eta = Eta_start; i_eta < Eta_stop; i_eta++) // eta gap bin
	{
	  for(Int_t i_sys = Sys_start; i_sys < Sys_stop; i_sys++)
	  {
	    for(Int_t i_proj = Proj_start; i_proj < Proj_stop; ++i_proj)
	    { 
	      TString KEY_Counts_Kaon = Form("Kaon_Counts_pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_%s_SysError_%d_Proj_%d",i_pt,i_cent,i_charge,i_eta,Order[mOrder].Data(),i_sys,i_proj); 
	      h_mCounts_kaon[KEY_Counts_Kaon]->Write();
	      TString KEY_Inte_Kaon = Form("Kaon_Inte_pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_%s_SysError_%d_Proj_%d",i_pt,i_cent,i_charge,i_eta,Order[mOrder].Data(),i_sys,i_proj); 
	      h_mCounts_kaon[KEY_Inte_Kaon]->Write();
	    }
	  }
	}
      }
    }
  }
  File_Kaon->Close();
}

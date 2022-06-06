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
static const Int_t pt_QA    = 7;

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

static const Int_t Sys_total = 6;
static const Int_t Sys_start = 0;
static const Int_t Sys_stop  = 6;

// Initial parameters 200 GeV
//Float_t order_pion[pt_total]   = {3.0,3.0,3.0,3.0,3.0,2.55,2.55,3.0,3.0,3.0,2.5,3.00,2.55,3.00,3.00,3.0}; // positive
//Float_t order_proton[pt_total] = {3.0,3.0,3.0,3.0,3.0,2.50,2.55,3.0,3.0,3.0,2.5,2.55,2.55,2.55,2.50,3.0};
Float_t order_pion[pt_total]   = {3.0,3.0,3.0,3.0,3.0,2.50,2.55,3.0,3.0,3.0,2.5,3.00,2.55,3.00,2.55,3.0}; //negative 
Float_t order_proton[pt_total] = {3.0,3.0,3.0,3.0,3.0,3.00,2.55,3.0,3.0,3.0,2.5,2.55,2.55,2.55,2.50,3.0};
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
//Float_t width_p_x[pt_total] = {0.03359,0.03034,0.03311,0.04130,0.05238,0.07000,0.08030,0.09100,0.11610,0.1384,0.158,0.1984,0.3567,0.3567,0.36,0.3567}; // positive
Float_t width_p_x[pt_total] = {0.03359,0.03034,0.03311,0.04130,0.05238,0.07200,0.08030,0.09100,0.11610,0.1384,0.158,0.1984,0.3567,0.3567,0.36,0.3567}; // negative
Float_t width_p_y[pt_total] = {0.01217,0.02191,0.03077,0.04085,0.05326,0.06540,0.07996,0.09776,0.11900,0.1419,0.169,0.2052,0.3111,0.3111,0.32,0.3111};
Float_t nu_p[pt_total] = {8.975,9.936,7.798,6.460,5.970,5.905,6.109,6.991,8.212,10.58,20.,47.92,49.,50.2,55,50.2};

typedef std::map<TString,TH1F*> TH1FMap;
typedef std::map<TString,TH2F*> TH2FMap;
typedef std::map<TString,TProfile*> TProMap;
typedef std::map<TString,TGraphAsymmErrors*> TGraMap;
typedef std::map<TString,std::vector<Float_t>> vecFMap;

// mMode:   0 for all pT bin, 1 for single pT bin
// mEnergy: 0 for 200GeV, 1 for 39GeV
// mCharge: 0 for positive, 1 for negative 
// mOrder:  0 for elliptic flow, 1 for triangular flow 
void CombPID_PiK(Int_t mMode = 0, Int_t mEnergy = 0, Int_t mCharge = 0, Int_t mOrder = 1)
{
  TGaxis::SetMaxDigits(4);

  TString InPutFile_Flow = Form("./Data/AuAu%s/nSigmaPion/Flow_%s.root",Energy[mEnergy].Data(),Energy[mEnergy].Data());
  TFile *File_Flow = TFile::Open(InPutFile_Flow.Data());

  TH2FMap h_mMass2_raw;
  for(Int_t i_pt = 0; i_pt < pt_total; i_pt++) // pt bin
  {
    for(Int_t i_cent = Cent_start; i_cent < Cent_stop; i_cent++) // centrality bin
    {
      for(Int_t i_charge = mCharge; i_charge < mCharge+1; i_charge++) // charge bin
      {
	for(Int_t i_eta = Eta_start; i_eta < Eta_stop; i_eta++) // eta gap bin
	{
	  for(Int_t i_phi = phi_start; i_phi < phi_stop; i_phi++) // phi-psi bin
	  {
	    for(Int_t i_sys = 0; i_sys < 6; i_sys++)
	    {
	      TString KEY_PiK_raw = Form("pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_phi_Psi_%d_%s_PiK_SysError_%d",i_pt,i_cent,i_charge,i_eta,i_phi,Order[mOrder].Data(),i_sys);
	      h_mMass2_raw[KEY_PiK_raw] = (TH2F*)File_Flow->Get(KEY_PiK_raw.Data());
	      /*
	      if(i_pt > 8)
	      {
		h_mMass2_raw[KEY_PiK_raw]->RebinX(2);
		h_mMass2_raw[KEY_PiK_raw]->RebinY(2);
		if(i_pt > 12)
		{
		  h_mMass2_raw[KEY_PiK_raw]->RebinX(2);
		  h_mMass2_raw[KEY_PiK_raw]->RebinY(2);
		}
	      }
	      */
	    }
	  }
	}
      }
    }
  }

  /*
  // QA plots Mass2 vs pT bins
  TCanvas *c_pT = new TCanvas("c_pT","c_pT",10,10,1600,1600);
  c_pT->Divide(4,4);
  for(Int_t i_pt = 0; i_pt < 16; i_pt++)
  {
    c_pT->cd(i_pt+1);
    c_pT->cd(i_pt+1)->SetLeftMargin(0.15);
    c_pT->cd(i_pt+1)->SetBottomMargin(0.15);
    c_pT->cd(i_pt+1)->SetTicks(1,1);
    c_pT->cd(i_pt+1)->SetGrid(0,0);
    c_pT->cd(i_pt+1)->SetLogz();
    TString KEY_PiK_raw_pT_QA = Form("pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_phi_Psi_%d_%s_PiK_SysError_%d",i_pt,Cent_start,mCharge,Eta_start,phi_start,Order[mOrder].Data(),Sys_start);
    h_mMass2_raw[KEY_PiK_raw_pT_QA]->SetStats(0);
    h_mMass2_raw[KEY_PiK_raw_pT_QA]->SetTitle("");
    h_mMass2_raw[KEY_PiK_raw_pT_QA]->GetXaxis()->SetTitle("x(n#sigma_{#pi},m^{2})");
    h_mMass2_raw[KEY_PiK_raw_pT_QA]->GetXaxis()->CenterTitle();
    h_mMass2_raw[KEY_PiK_raw_pT_QA]->GetXaxis()->SetTitleSize(0.06);
    h_mMass2_raw[KEY_PiK_raw_pT_QA]->GetYaxis()->SetTitle("y(n#sigma_{#pi},m^{2})");
    h_mMass2_raw[KEY_PiK_raw_pT_QA]->GetYaxis()->CenterTitle();
    h_mMass2_raw[KEY_PiK_raw_pT_QA]->GetYaxis()->SetTitleSize(0.06);
    h_mMass2_raw[KEY_PiK_raw_pT_QA]->Draw("colz");

    TString pT_range = Form("[%.2f,%.2f]",pt_low_raw[i_pt],pt_up_raw[i_pt]);
    plotTopLegend((char*)pT_range.Data(),0.65,0.7,0.08,1,0.0,42,1);
  }
  */

  // pT rebin
  TH2FMap h_mMass2; // rebinned m2 distribution
  for(Int_t i_cent = Cent_start; i_cent < Cent_stop; i_cent++) // Centrality loop
  {
    for(Int_t i_charge = mCharge; i_charge < mCharge+1; i_charge++) // charge bin
    {
      for(Int_t i_eta = Eta_start; i_eta < Eta_stop; i_eta++) // EtaGap loop
      {
	for(Int_t i_phi = phi_start; i_phi < phi_stop; i_phi++) // phi loop
	{
	  for(Int_t i_sys = Sys_start; i_sys < Sys_stop; i_sys++) // Systematic loop
	  {
	    for(Int_t pt_bin = pt_rebin_first; pt_bin < pt_rebin_last; pt_bin++) // pt loop
	    {
	      TString KEY_PiK = Form("pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_phi_Psi_%d_%s_PiK_SysError_%d",pt_bin,i_cent,i_charge,i_eta,i_phi,Order[mOrder].Data(),i_sys);
	      for(Int_t i_pt = pt_rebin_start[pt_bin]; i_pt <= pt_rebin_stop[pt_bin]; i_pt++)
	      {
		TString KEY_PiK_raw = Form("pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_phi_Psi_%d_%s_PiK_SysError_%d",i_pt,i_cent,i_charge,i_eta,i_phi,Order[mOrder].Data(),i_sys);
		if(i_pt == pt_rebin_start[pt_bin]) h_mMass2[KEY_PiK] = (TH2F*)h_mMass2_raw[KEY_PiK_raw]->Clone();
		else h_mMass2[KEY_PiK]->Add(h_mMass2_raw[KEY_PiK_raw],1.0);
	      }
	    }
	  }
	}
      }
    }
  }

  /*
  // QA plots Mass2 vs pT bins after pT rebin
  TCanvas *c_pT_rebin = new TCanvas("c_pT_rebin","c_pT_rebin",10,10,1600,1600);
  c_pT_rebin->Divide(4,4);
  for(Int_t i_pt = pt_rebin_first; i_pt < pt_rebin_last; i_pt++)
  {
    c_pT_rebin->cd(pt_rebin_start[i_pt]+1);
    c_pT_rebin->cd(pt_rebin_start[i_pt]+1)->SetLeftMargin(0.15);
    c_pT_rebin->cd(pt_rebin_start[i_pt]+1)->SetBottomMargin(0.15);
    c_pT_rebin->cd(pt_rebin_start[i_pt]+1)->SetTicks(1,1);
    c_pT_rebin->cd(pt_rebin_start[i_pt]+1)->SetGrid(0,0);
    c_pT_rebin->cd(pt_rebin_start[i_pt]+1)->SetLogz();
    TString KEY_PiK_raw_pT_QA = Form("pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_phi_Psi_%d_%s_PiK_SysError_%d",i_pt,Cent_start,mCharge,Eta_start,phi_start,Order[mOrder].Data(),Sys_start);
    h_mMass2[KEY_PiK_raw_pT_QA]->SetStats(0);
    h_mMass2[KEY_PiK_raw_pT_QA]->SetTitle("");
    h_mMass2[KEY_PiK_raw_pT_QA]->GetXaxis()->SetTitle("x(n#sigma_{#pi},m^{2})");
    h_mMass2[KEY_PiK_raw_pT_QA]->GetXaxis()->CenterTitle();
    h_mMass2[KEY_PiK_raw_pT_QA]->GetXaxis()->SetTitleSize(0.06);
    h_mMass2[KEY_PiK_raw_pT_QA]->GetYaxis()->SetTitle("y(n#sigma_{#pi},m^{2})");
    h_mMass2[KEY_PiK_raw_pT_QA]->GetYaxis()->CenterTitle();
    h_mMass2[KEY_PiK_raw_pT_QA]->GetYaxis()->SetTitleSize(0.06);
    h_mMass2[KEY_PiK_raw_pT_QA]->Draw("colz");

    TString pT_range = Form("[%.2f,%.2f]",pt_low[i_pt],pt_up[i_pt]);
    plotTopLegend((char*)pT_range.Data(),0.65,0.7,0.08,1,0.0,42,1);
  }
  */

  // integrate over phi to get first guess of initial parameters
  TH2FMap h_mMass2_total;
  vecFMap ParStudent_total_1st;
  vecFMap ParStudent_total_2nd; 
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
	    TString KEY_PiK_total = Form("pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_%s_PiK_SysError_%d",i_pt,i_cent,i_charge,i_eta,Order[mOrder].Data(),i_sys);
	    for(Int_t i_phi = phi_start; i_phi < phi_stop; i_phi++) // phi-psi bin
	    {
	      TString KEY_PiK = Form("pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_phi_Psi_%d_%s_PiK_SysError_%d",i_pt,i_cent,i_charge,i_eta,i_phi,Order[mOrder].Data(),i_sys);
	      if(i_phi == phi_start) h_mMass2_total[KEY_PiK_total] = (TH2F*)h_mMass2[KEY_PiK]->Clone();
	      else h_mMass2_total[KEY_PiK_total]->Add(h_mMass2[KEY_PiK],1.0);
	    }
	    TF2 *f_student_total_1st = new TF2("f_student_total_1st",student_t_2d_fit,x_low[i_pt],x_up[i_pt],y_low[i_pt],y_up[i_pt],19);
	    for(Int_t i_par = 0; i_par < 19; i_par++)
	    {
	      f_student_total_1st->ReleaseParameter(i_par);
	    }
	    f_student_total_1st->SetParLimits(0,0,nu_pi[i_pt]*10.0);
	    f_student_total_1st->SetParLimits(1,0,nu_pi[i_pt]*10.0);
	    f_student_total_1st->SetParLimits(2,x_pi[i_pt]-0.08,x_pi[i_pt]+0.08);
	    f_student_total_1st->SetParLimits(3,y_pi[i_pt]-0.08,y_pi[i_pt]+0.08);
	    f_student_total_1st->SetParLimits(4,width_pi[i_pt]*0.4,width_pi[i_pt]*1.6);
	    f_student_total_1st->SetParLimits(5,width_pi[i_pt]*0.4,width_pi[i_pt]*1.6);
	    f_student_total_1st->SetParLimits(7,x_k[i_pt]-0.18,x_k[i_pt]+0.18);
	    f_student_total_1st->SetParLimits(8,y_k[i_pt]-0.18,y_k[i_pt]+0.18);
	    f_student_total_1st->SetParLimits(9,width_k[i_pt]*0.4,width_k[i_pt]*1.6);
	    f_student_total_1st->SetParLimits(10,width_k[i_pt]*0.4,width_k[i_pt]*1.6);
	    f_student_total_1st->SetParLimits(12,0,nu_p[i_pt]*10.0);
//	    f_student_total_1st->SetParLimits(13,0,nu_p[i_pt]*10.0);
	    f_student_total_1st->SetParLimits(14,x_p[i_pt]-0.08,x_p[i_pt]+0.08);
	    f_student_total_1st->SetParLimits(15,y_p[i_pt]-0.08,y_p[i_pt]+0.08);
	    f_student_total_1st->SetParLimits(16,width_p_x[i_pt]*0.4,width_p_x[i_pt]*1.6);
	    f_student_total_1st->SetParLimits(17,width_p_y[i_pt]*0.4,width_p_y[i_pt]*1.6);

	    // parameter for pion start
	    f_student_total_1st->SetParameter(0,nu_pi[i_pt]);
	    f_student_total_1st->SetParameter(1,nu_pi[i_pt]);
	    f_student_total_1st->SetParameter(2,x_pi[i_pt]);
	    f_student_total_1st->SetParameter(3,y_pi[i_pt]);
	    f_student_total_1st->SetParameter(4,width_pi[i_pt]);
	    f_student_total_1st->SetParameter(5,width_pi[i_pt]);
	    f_student_total_1st->SetParameter(6,1);
	    // parameter for pion stop

	    // parameter for kaon start
	    f_student_total_1st->SetParameter(7,x_k[i_pt]);
	    f_student_total_1st->SetParameter(8,y_k[i_pt]);
	    f_student_total_1st->SetParameter(9,width_k[i_pt]);
	    f_student_total_1st->SetParameter(10,width_k[i_pt]);
	    f_student_total_1st->SetParameter(11,1);
	    // parameter for kaon stop

	    // parameter for proton start
	    f_student_total_1st->SetParameter(12,nu_p[i_pt]);
	    f_student_total_1st->SetParameter(13,nu_p[i_pt]);
	    f_student_total_1st->SetParameter(14,x_p[i_pt]);
	    f_student_total_1st->SetParameter(15,y_p[i_pt]);
	    f_student_total_1st->SetParameter(16,width_p_x[i_pt]);
	    f_student_total_1st->SetParameter(17,width_p_y[i_pt]);
	    f_student_total_1st->SetParameter(18,1);
	    // parameter for proton stop

	    Float_t pion_max = h_mMass2_total[KEY_PiK_total]->GetBinContent(h_mMass2_total[KEY_PiK_total]->GetXaxis()->FindBin(x_pi[i_pt]),h_mMass2_total[KEY_PiK_total]->GetYaxis()->FindBin(y_pi[i_pt]));
	    Float_t pion_fit = f_student_total_1st->Eval(x_pi[i_pt],y_pi[i_pt]);
	    Float_t pion_Nrom = pion_max/pion_fit;

	    Float_t kaon_max = h_mMass2_total[KEY_PiK_total]->GetBinContent(h_mMass2_total[KEY_PiK_total]->GetXaxis()->FindBin(x_k[i_pt]),h_mMass2_total[KEY_PiK_total]->GetYaxis()->FindBin(y_k[i_pt]));
	    Float_t kaon_fit = f_student_total_1st->Eval(x_k[i_pt],y_k[i_pt]);
	    Float_t kaon_Nrom = kaon_max/kaon_fit;

	    Float_t proton_max = h_mMass2_total[KEY_PiK_total]->GetBinContent(h_mMass2_total[KEY_PiK_total]->GetXaxis()->FindBin(x_p[i_pt]),h_mMass2_total[KEY_PiK_total]->GetYaxis()->FindBin(y_p[i_pt]));
	    Float_t proton_fit = f_student_total_1st->Eval(x_p[i_pt],y_p[i_pt]);
	    Float_t proton_Nrom = proton_max/proton_fit;

	    f_student_total_1st->SetParameter(6,pion_Nrom);
	    f_student_total_1st->SetParameter(11,kaon_Nrom);
	    f_student_total_1st->SetParameter(18,proton_Nrom);

	    f_student_total_1st->SetRange(x_pi[i_pt]-(order_pion[i_pt]+0.5)*width_pi[i_pt],y_p[i_pt],x_p[i_pt],y_pi[i_pt]+(order_pion[i_pt]+0.5)*width_pi[i_pt]);

	    f_student_total_1st->SetLineStyle(1);
	    f_student_total_1st->SetLineColor(2);

	    if((mMode == 0 && i_pt > 4) || (mMode == 1 && i_pt == pt_QA && i_pt > 4))
	    {
//	      cout << endl;
//	      cout<< "first fits:"  << " i_pt = " << i_pt << ", i_sys = " << i_sys << endl;
	      h_mMass2_total[KEY_PiK_total]->Fit(f_student_total_1st,"MQRN"); // first fits for total distribution
	      for(Int_t i_par = 0; i_par < 19; i_par++)
	      {
		ParStudent_total_1st[KEY_PiK_total].push_back(static_cast<Float_t>(f_student_total_1st->GetParameter(i_par)));
	      }

	      TF2 *f_student_total_2nd = new TF2("f_student_total_2nd",student_t_2d_fit,x_low[i_pt],x_up[i_pt],y_low[i_pt],y_up[i_pt],19);
	      for(Int_t i_par = 0; i_par < 19; i_par++)
	      {
		f_student_total_2nd->ReleaseParameter(i_par);
		f_student_total_2nd->SetParameter(i_par,ParStudent_total_1st[KEY_PiK_total][i_par]);
	      }
	      f_student_total_2nd->SetParLimits(0,0,342);
	      f_student_total_2nd->SetParLimits(1,0,342);
	      f_student_total_2nd->SetParLimits(2,x_pi[i_pt]-0.08,x_pi[i_pt]+0.08);
	      f_student_total_2nd->SetParLimits(3,y_pi[i_pt]-0.08,y_pi[i_pt]+0.08);
	      f_student_total_2nd->SetParLimits(7,x_k[i_pt]-0.18,x_k[i_pt]+0.18);
	      f_student_total_2nd->SetParLimits(8,y_k[i_pt]-0.18,y_k[i_pt]+0.18);
	      if(i_pt == 12) f_student_total_2nd->SetParLimits(8,y_k[i_pt]-0.36,y_k[i_pt]+0.36);
	      f_student_total_2nd->SetParLimits(12,0,342);
	      f_student_total_2nd->SetParLimits(13,0,342);
	      f_student_total_2nd->SetParLimits(14,x_p[i_pt]-0.08,x_p[i_pt]+0.08);
	      if(i_pt >= 14) f_student_total_2nd->SetParLimits(14,x_p[i_pt]-0.16,x_p[i_pt]+0.16);
	      f_student_total_2nd->SetParLimits(15,y_p[i_pt]-0.08,y_p[i_pt]+0.08);
	      if(i_pt >= 14) f_student_total_2nd->SetParLimits(15,y_p[i_pt]-0.16,y_p[i_pt]+0.16);
	      f_student_total_2nd->SetRange(x_pi[i_pt]-order_pion[i_pt]*ParStudent_total_1st[KEY_PiK_total][4],y_p[i_pt]-order_proton[i_pt]*ParStudent_total_1st[KEY_PiK_total][17],x_p[i_pt]+order_proton[i_pt]*ParStudent_total_1st[KEY_PiK_total][16],y_pi[i_pt]+order_pion[i_pt]*ParStudent_total_1st[KEY_PiK_total][5]);

	      f_student_total_2nd->SetLineStyle(3);
	      f_student_total_2nd->SetLineColor(4);

	      cout << endl;
	      cout<< "second fits:"  << " i_pt = " << i_pt << ", i_sys = " << i_sys << endl;
	      h_mMass2_total[KEY_PiK_total]->Fit(f_student_total_2nd,"MRN"); // second fits for total distribution
	      for(Int_t i_par = 0; i_par < 19; i_par++)
	      {
		ParStudent_total_2nd[KEY_PiK_total].push_back(static_cast<Float_t>(f_student_total_2nd->GetParameter(i_par)));
	      }
	    }
	  }
	}
      }
    }
  }

  /*
  // QA plots phi integrated distribution
  TCanvas *c_pT_rebin_total = new TCanvas("c_pT_rebin_total","c_pT_rebin_total",10,10,1600,1600);
  c_pT_rebin_total->Divide(4,4);
  for(Int_t i_pt = pt_rebin_first; i_pt < pt_rebin_last; i_pt++)
  {
    c_pT_rebin_total->cd(pt_rebin_start[i_pt]+1);
    c_pT_rebin_total->cd(pt_rebin_start[i_pt]+1)->SetLeftMargin(0.15);
    c_pT_rebin_total->cd(pt_rebin_start[i_pt]+1)->SetBottomMargin(0.15);
    c_pT_rebin_total->cd(pt_rebin_start[i_pt]+1)->SetTicks(1,1);
    c_pT_rebin_total->cd(pt_rebin_start[i_pt]+1)->SetGrid(0,0);
    c_pT_rebin_total->cd(pt_rebin_start[i_pt]+1)->SetLogz();
    TString KEY_PiK_total_pT_QA = Form("pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_%s_PiK_SysError_%d",i_pt,Cent_start,mCharge,Eta_start,Order[mOrder].Data(),Sys_start);
    h_mMass2_total[KEY_PiK_total_pT_QA]->SetStats(0);
    h_mMass2_total[KEY_PiK_total_pT_QA]->SetTitle("");
    h_mMass2_total[KEY_PiK_total_pT_QA]->GetXaxis()->SetTitle("x(n#sigma_{#pi},m^{2})");
    h_mMass2_total[KEY_PiK_total_pT_QA]->GetXaxis()->CenterTitle();
    h_mMass2_total[KEY_PiK_total_pT_QA]->GetXaxis()->SetTitleSize(0.06);
    h_mMass2_total[KEY_PiK_total_pT_QA]->GetYaxis()->SetTitle("y(n#sigma_{#pi},m^{2})");
    h_mMass2_total[KEY_PiK_total_pT_QA]->GetYaxis()->CenterTitle();
    h_mMass2_total[KEY_PiK_total_pT_QA]->GetYaxis()->SetTitleSize(0.06);
    h_mMass2_total[KEY_PiK_total_pT_QA]->SetMarkerStyle(24);
    h_mMass2_total[KEY_PiK_total_pT_QA]->SetMarkerColor(kGray+3);
    h_mMass2_total[KEY_PiK_total_pT_QA]->SetMarkerSize(0.4);
    h_mMass2_total[KEY_PiK_total_pT_QA]->Draw("colz");

    if((mMode == 0 && i_pt > 4) || (mMode == 1 && i_pt == pt_QA && i_pt > 4))
    {
      TF2 *f_student_total_1st_QA = new TF2("f_student_total_1st_QA",student_t_2d_fit,x_low[i_pt],x_up[i_pt],y_low[i_pt],y_up[i_pt],19);
      for(Int_t i_par = 0; i_par < 19; i_par++)
      {
	f_student_total_1st_QA->SetParameter(i_par,ParStudent_total_1st[KEY_PiK_total_pT_QA][i_par]);
      }
      f_student_total_1st_QA->SetRange(x_pi[i_pt]-(order_pion[i_pt]+0.5)*width_pi[i_pt],y_p[i_pt],x_p[i_pt],y_pi[i_pt]+(order_pion[i_pt]+0.5)*width_pi[i_pt]);
      f_student_total_1st_QA->SetLineStyle(1);
      f_student_total_1st_QA->SetLineColor(2);
      f_student_total_1st_QA->Draw("cont3 same");

      TF2 *f_student_total_2nd_QA = new TF2("f_student_total_2nd_QA",student_t_2d_fit,x_low[i_pt],x_up[i_pt],y_low[i_pt],y_up[i_pt],19);
      for(Int_t i_par = 0; i_par < 19; i_par++)
      {
	f_student_total_2nd_QA->SetParameter(i_par,ParStudent_total_2nd[KEY_PiK_total_pT_QA][i_par]);
      }
      f_student_total_2nd_QA->SetRange(x_pi[i_pt]-order_pion[i_pt]*ParStudent_total_1st[KEY_PiK_total_pT_QA][4],y_p[i_pt]-order_proton[i_pt]*ParStudent_total_1st[KEY_PiK_total_pT_QA][17],x_p[i_pt]+order_proton[i_pt]*ParStudent_total_1st[KEY_PiK_total_pT_QA][16],y_pi[i_pt]+order_pion[i_pt]*ParStudent_total_1st[KEY_PiK_total_pT_QA][5]);
      f_student_total_2nd_QA->SetLineStyle(3);
      f_student_total_2nd_QA->SetLineColor(4);
      f_student_total_2nd_QA->Draw("cont3 same");
    }

    TString pT_range = Form("[%.2f,%.2f]",pt_low[i_pt],pt_up[i_pt]);
    plotTopLegend((char*)pT_range.Data(),0.65,0.7,0.08,1,0.0,42,1);
  }
  */

  // fit for phi differential bin and subtract proton
  vecFMap ParStudent; 
  TH2FMap h_mMass2_sub;
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
	    TString KEY_PiK_total = Form("pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_%s_PiK_SysError_%d",i_pt,i_cent,i_charge,i_eta,Order[mOrder].Data(),i_sys);
	    for(Int_t i_phi = phi_start; i_phi < phi_stop; i_phi++) // phi-psi bin
	    {
	      TString KEY_PiK = Form("pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_phi_Psi_%d_%s_PiK_SysError_%d",i_pt,i_cent,i_charge,i_eta,i_phi,Order[mOrder].Data(),i_sys);
	      h_mMass2_sub[KEY_PiK] = (TH2F*)h_mMass2[KEY_PiK]->Clone();
	      if((mMode == 0 && i_pt > 4) || (mMode == 1 && i_pt == pt_QA && i_pt > 4))
	      {
		TF2 *f_student = new TF2("f_student",student_t_2d_fit,x_low[i_pt],x_up[i_pt],y_low[i_pt],y_up[i_pt],19);
		for(Int_t i_par = 0; i_par < 19; i_par++)
		{
		  f_student->ReleaseParameter(i_par);
		  f_student->FixParameter(i_par,ParStudent_total_2nd[KEY_PiK_total][i_par]);
		  if(i_par == 6 || i_par == 11 || i_par == 18)
		  {
		    f_student->ReleaseParameter(i_par);
		    f_student->SetParameter(i_par,1.0);
		    f_student->SetParError(i_par,0.0);
		  }
		}
		Float_t pion_max = h_mMass2[KEY_PiK]->GetBinContent(h_mMass2[KEY_PiK]->GetXaxis()->FindBin(ParStudent_total_2nd[KEY_PiK_total][2]),h_mMass2[KEY_PiK]->GetYaxis()->FindBin(ParStudent_total_2nd[KEY_PiK_total][3]));
		Float_t pion_fit = f_student->Eval(ParStudent_total_2nd[KEY_PiK_total][2],ParStudent_total_2nd[KEY_PiK_total][3]);
		Float_t pion_Nrom = pion_max/pion_fit;

		Float_t kaon_max = h_mMass2[KEY_PiK]->GetBinContent(h_mMass2[KEY_PiK]->GetXaxis()->FindBin(ParStudent_total_2nd[KEY_PiK_total][7]),h_mMass2[KEY_PiK]->GetYaxis()->FindBin(ParStudent_total_2nd[KEY_PiK_total][8]));
		Float_t kaon_fit = f_student->Eval(ParStudent_total_2nd[KEY_PiK_total][7],ParStudent_total_2nd[KEY_PiK_total][8]);
		Float_t kaon_Nrom = kaon_max/kaon_fit;

		Float_t proton_max = h_mMass2[KEY_PiK]->GetBinContent(h_mMass2[KEY_PiK]->GetXaxis()->FindBin(ParStudent_total_2nd[KEY_PiK_total][14]),h_mMass2[KEY_PiK]->GetYaxis()->FindBin(ParStudent_total_2nd[KEY_PiK_total][15]));
		Float_t proton_fit = f_student->Eval(ParStudent_total_2nd[KEY_PiK_total][14],ParStudent_total_2nd[KEY_PiK_total][15]);
		Float_t proton_Nrom = proton_max/proton_fit;

		f_student->SetParameter(6,pion_Nrom);
		f_student->SetParameter(11,kaon_Nrom);
		f_student->SetParameter(18,proton_Nrom);

		f_student->SetRange(x_pi[i_pt]-order_pion[i_pt]*ParStudent_total_2nd[KEY_PiK_total][4],y_p[i_pt]-order_proton[i_pt]*ParStudent_total_2nd[KEY_PiK_total][17],x_p[i_pt]+order_proton[i_pt]*ParStudent_total_2nd[KEY_PiK_total][16],y_pi[i_pt]+order_pion[i_pt]*ParStudent_total_2nd[KEY_PiK_total][5]);
		cout << "i_pt = " << i_pt << ", i_charge = " << i_charge << ", i_eta = " << i_eta << ", i_phi = " << i_phi << endl;
		h_mMass2[KEY_PiK]->Fit(f_student,"MRN");
		for(Int_t i_par = 0; i_par < 19; i_par++)
		{
		  ParStudent[KEY_PiK].push_back(static_cast<Float_t>(f_student->GetParameter(i_par)));
		}

		TF2 *f_proton = new TF2("f_proton",student_t_2d_single,x_low[i_pt],x_up[i_pt],y_low[i_pt],y_up[i_pt],7);
		for(Int_t i_par = 0; i_par < 7; i_par++)
		{
		  f_proton->ReleaseParameter(i_par);
		  f_proton->FixParameter(i_par,ParStudent[KEY_PiK][i_par+12]);
		}
		h_mMass2_sub[KEY_PiK]->Add(f_proton,-1.0);
	      }
	    }
	  }
	}
      }
    }
  }

  /*
  // QA plots phi differential distribution after proton subtraction
  TCanvas *c_pT_sub = new TCanvas("c_pT_sub","c_pT_sub",10,10,1600,1600);
  c_pT_sub->Divide(4,4);
  for(Int_t i_pt = pt_rebin_first; i_pt < pt_rebin_last; i_pt++)
  {
    c_pT_sub->cd(pt_rebin_start[i_pt]+1);
    c_pT_sub->cd(pt_rebin_start[i_pt]+1)->SetLeftMargin(0.15);
    c_pT_sub->cd(pt_rebin_start[i_pt]+1)->SetBottomMargin(0.15);
    c_pT_sub->cd(pt_rebin_start[i_pt]+1)->SetTicks(1,1);
    c_pT_sub->cd(pt_rebin_start[i_pt]+1)->SetGrid(0,0);
    c_pT_sub->cd(pt_rebin_start[i_pt]+1)->SetLogz();
    TString KEY_PiK_QA = Form("pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_phi_Psi_%d_%s_PiK_SysError_%d",i_pt,Cent_start,mCharge,Eta_start,phi_start,Order[mOrder].Data(),Sys_start);
    h_mMass2_sub[KEY_PiK_QA]->SetStats(0);
    h_mMass2_sub[KEY_PiK_QA]->SetTitle("");
    h_mMass2_sub[KEY_PiK_QA]->GetXaxis()->SetTitle("x(n#sigma_{#pi},m^{2})");
    h_mMass2_sub[KEY_PiK_QA]->GetXaxis()->CenterTitle();
    h_mMass2_sub[KEY_PiK_QA]->GetXaxis()->SetTitleSize(0.06);
    h_mMass2_sub[KEY_PiK_QA]->GetYaxis()->SetTitle("y(n#sigma_{#pi},m^{2})");
    h_mMass2_sub[KEY_PiK_QA]->GetYaxis()->CenterTitle();
    h_mMass2_sub[KEY_PiK_QA]->GetYaxis()->SetTitleSize(0.06);
    h_mMass2_sub[KEY_PiK_QA]->SetMarkerStyle(24);
    h_mMass2_sub[KEY_PiK_QA]->SetMarkerColor(kGray+3);
    h_mMass2_sub[KEY_PiK_QA]->SetMarkerSize(0.4);
    h_mMass2_sub[KEY_PiK_QA]->Draw("colz");

    TString pT_range = Form("[%.2f,%.2f]",pt_low[i_pt],pt_up[i_pt]);
    plotTopLegend((char*)pT_range.Data(),0.65,0.7,0.08,1,0.0,42,1);
  }
  */

  // Save fit parameters and proton subtracted distribution
  TString OutPutFile = Form("./OutPut/AuAu%s/nSigmaPion/nSigmaPion_Sub_Charge_%d.root",Energy[mEnergy].Data(),mCharge);
  TFile *File_OutPut = new TFile(OutPutFile.Data(),"RECREATE");
  File_OutPut->cd();
  TH1FMap h_ParStudnet_total_1st, h_ParStudnet_total_2nd;
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
	    if((mMode == 0 && i_pt > 4) || (mMode == 1 && i_pt == pt_QA && i_pt > 4))
	    {
	      TString KEY_PiK_total = Form("pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_%s_PiK_SysError_%d",i_pt,i_cent,i_charge,i_eta,Order[mOrder].Data(),i_sys);
	      h_ParStudnet_total_1st[KEY_PiK_total] = new TH1F((KEY_PiK_total+"_1st").Data(),(KEY_PiK_total+"_1st").Data(),19,-0.5,18.5);
	      h_ParStudnet_total_2nd[KEY_PiK_total] = new TH1F((KEY_PiK_total+"_2nd").Data(),(KEY_PiK_total+"_2nd").Data(),19,-0.5,18.5);
	      for(Int_t i_par = 0; i_par < 19; i_par++)
	      {
		h_ParStudnet_total_1st[KEY_PiK_total]->SetBinContent(i_par+1,ParStudent_total_1st[KEY_PiK_total][i_par]);
		h_ParStudnet_total_2nd[KEY_PiK_total]->SetBinContent(i_par+1,ParStudent_total_2nd[KEY_PiK_total][i_par]);
	      }
	      h_ParStudnet_total_1st[KEY_PiK_total]->Write();
	      h_ParStudnet_total_2nd[KEY_PiK_total]->Write();
	    }
	    for(Int_t i_phi = phi_start; i_phi < phi_stop; i_phi++) // phi-psi bin
	    {
	      TString KEY_PiK = Form("pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_phi_Psi_%d_%s_PiK_SysError_%d",i_pt,i_cent,i_charge,i_eta,i_phi,Order[mOrder].Data(),i_sys);
	      h_mMass2_sub[KEY_PiK]->Write();
	    }
	  }
	}
      }
    }
  }
  File_OutPut->Close();
}

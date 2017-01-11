#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "draw.h"
#include "TString.h"
#include "TMath.h"
#include "TGraphAsymmErrors.h"
#include "TGaxis.h"
#include "TLegend.h"
#include "TProfile.h"
#include <map>
#include <vector>

Double_t flow_3(Double_t *x_val, Double_t *par)
{
  Double_t x, y;
  Double_t Ampl, v3;
  x = x_val[0];
  Ampl = par[0];
  v3 = par[1];

  y = Ampl*(1.0 + 2.0*v3*TMath::Cos(3.0*x));

  return y;
}

Double_t ErrorAdd(Double_t x, Double_t y)
{
    return sqrt(x*x+y*y);
}

Double_t ErrTimes(Double_t x, Double_t y, Double_t dx, Double_t dy)
{
    return x*y*ErrorAdd(dx/x,dy/y);
}

Double_t ErrDiv(Double_t x, Double_t y, Double_t dx, Double_t dy)
{
    return x/y*ErrorAdd(dx/x,dy/y);
}

static const TString Energy[2] = {"200GeV","39GeV"};
static const TString Order[2] = {"2nd","3rd"};
static const Double_t PI_max[2] = {TMath::Pi()/2.0,TMath::Pi()/3.0};

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

static const Int_t Cent_total = 4; 
static const Int_t Cent_start = 0;
static const Int_t Cent_stop  = 1;
static const Int_t cent_low[4] = {0,7,4,0}; // 0 = 0-80%, 1 = 0-10%, 2 = 10-40%, 3 = 40-80%
static const Int_t cent_up[4]  = {8,8,6,3}; // 0 = 0-80%, 1 = 0-10%, 2 = 10-40%, 3 = 40-80%

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

static const TString PID[2] = {"Pion","Kaon"};
static const Float_t Inte_start[2] = {-0.05,0.15};
static const Float_t Inte_stop[2] = {0.05,0.30};

typedef std::map<TString,TH1F*> TH1FMap;
typedef std::map<TString,TProfile*> TProMap;
typedef std::map<TString,TGraphAsymmErrors*> TGraMap;
typedef std::map<TString,std::vector<Float_t>> vecFMap;

#ifndef _PlotQA_
#define _PlotQA_ 0
#endif

void PiKFlow(Int_t mEnergy = 0, Int_t mCharge = 0, Int_t mOrder = 1, Int_t mPID = 0)
{
  TGaxis::SetMaxDigits(6);

  TString InPutFile = Form("./OutPut/AuAu%s/%s/Counts_%s_Charge_%d.root",Energy[mEnergy].Data(),PID[mPID].Data(),Order[mOrder].Data(),mCharge);
  cout << "InPutFile set to: " << InPutFile.Data() << endl;
  TFile *File_InPut = TFile::Open(InPutFile.Data());

  TH1FMap h_mCounts;
  vecFMap ParFlow;
  TGraMap g_mRawFlow;
  for(Int_t i_cent = Cent_start; i_cent < Cent_stop; i_cent++) // centrality loop
  {
    for(Int_t i_charge = mCharge; i_charge < mCharge+1; i_charge++) // charge loop 
    {
      for(Int_t i_eta = Eta_start; i_eta < Eta_stop; i_eta++) // eta gap loop
      {
	for(Int_t i_sys = Sys_start; i_sys < Sys_stop; i_sys++) // sys cut loop
	{
	  for(Int_t i_proj = Proj_start; i_proj < Proj_stop; ++i_proj) // proj range loop
	  { 
	    TString KEY_Flow_Counts = Form("%s_Counts_Centrality_%d_Charge_%d_EtaGap_%d_%s_SysError_%d_Proj_%d",PID[mPID].Data(),i_cent,i_charge,i_eta,Order[mOrder].Data(),i_sys,i_proj); 
	    g_mRawFlow[KEY_Flow_Counts] = new TGraphAsymmErrors();
	    TString KEY_Flow_Inte = Form("%s_Inte_Centrality_%d_Charge_%d_EtaGap_%d_%s_SysError_%d_Proj_%d",PID[mPID].Data(),i_cent,i_charge,i_eta,Order[mOrder].Data(),i_sys,i_proj); 
	    g_mRawFlow[KEY_Flow_Inte] = new TGraphAsymmErrors();
	    for(Int_t i_pt = pt_rebin_first; i_pt < pt_rebin_last; i_pt++) // pt bin
	    {
	      Float_t pt_mean = (pt_low[i_pt]+pt_up[i_pt])/2.0;

	      TString KEY_Counts = Form("%s_Counts_pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_%s_SysError_%d_Proj_%d",PID[mPID].Data(),i_pt,i_cent,i_charge,i_eta,Order[mOrder].Data(),i_sys,i_proj); 
	      h_mCounts[KEY_Counts] = (TH1F*)File_InPut->Get(KEY_Counts.Data());
	      TF1 *f_flowCounts = new TF1("f_flowCounts",flow_3,0.0,PI_max[mOrder],2);
	      f_flowCounts->ReleaseParameter(0);
	      f_flowCounts->ReleaseParameter(1);
	      f_flowCounts->SetParameter(0,h_mCounts[KEY_Counts]->GetMaximum());
	      f_flowCounts->SetParameter(1,0.1);
	      h_mCounts[KEY_Counts]->Fit(f_flowCounts,"NQ");
	      ParFlow[KEY_Counts].clear();
	      ParFlow[KEY_Counts].push_back(static_cast<Float_t>(f_flowCounts->GetParameter(0)));
	      ParFlow[KEY_Counts].push_back(static_cast<Float_t>(f_flowCounts->GetParameter(1)));
	      g_mRawFlow[KEY_Flow_Counts]->SetPoint(i_pt,pt_mean,f_flowCounts->GetParameter(1));
	      g_mRawFlow[KEY_Flow_Counts]->SetPointError(i_pt,0.0,0.0,f_flowCounts->GetParError(1),f_flowCounts->GetParError(1));

	      TString KEY_Inte = Form("%s_Inte_pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_%s_SysError_%d_Proj_%d",PID[mPID].Data(),i_pt,i_cent,i_charge,i_eta,Order[mOrder].Data(),i_sys,i_proj); 
	      h_mCounts[KEY_Inte] = (TH1F*)File_InPut->Get(KEY_Inte.Data());
	      TF1 *f_flowInte = new TF1("f_flowInte",flow_3,0.0,PI_max[mOrder],2);
	      f_flowInte->ReleaseParameter(0);
	      f_flowInte->ReleaseParameter(1);
	      f_flowInte->SetParameter(0,h_mCounts[KEY_Inte]->GetMaximum());
	      f_flowInte->SetParameter(1,0.1);
	      if(i_pt > 4) h_mCounts[KEY_Inte]->Fit(f_flowInte,"NQ");
	      ParFlow[KEY_Inte].clear();
	      ParFlow[KEY_Inte].push_back(static_cast<Float_t>(f_flowInte->GetParameter(0)));
	      ParFlow[KEY_Inte].push_back(static_cast<Float_t>(f_flowInte->GetParameter(1)));
	      g_mRawFlow[KEY_Flow_Inte]->SetPoint(i_pt,pt_mean,f_flowInte->GetParameter(1));
	      g_mRawFlow[KEY_Flow_Inte]->SetPointError(i_pt,0.0,0.0,f_flowInte->GetParError(1),f_flowInte->GetParError(1));
	      if(i_pt <= 4)
	      {
		g_mRawFlow[KEY_Flow_Inte]->SetPoint(i_pt,pt_mean,-10);
		g_mRawFlow[KEY_Flow_Inte]->SetPointError(i_pt,0.0,0.0,0.0,0.0);
	      }
	    }
	  }
	}
      }
    }
  }

#if _PlotQA_
  TCanvas *c_Counts = new TCanvas("c_Counts","c_Counts",1200,1200);
  c_Counts->Divide(4,4);
  for(Int_t i_pt = pt_rebin_first; i_pt < pt_rebin_last; i_pt++) // pt bin
  {
    c_Counts->cd(i_pt+1);
    c_Counts->cd(i_pt+1)->SetLeftMargin(0.15);
    c_Counts->cd(i_pt+1)->SetBottomMargin(0.15);
    c_Counts->cd(i_pt+1)->SetTicks(1,1);
    c_Counts->cd(i_pt+1)->SetGrid(0,0);
    TString KEY_Counts_QA = Form("%s_Counts_pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_%s_SysError_%d_Proj_%d",PID[mPID].Data(),i_pt,Cent_start,mCharge,Eta_start,Order[mOrder].Data(),Sys_QA,Proj_QA); 
    h_mCounts[KEY_Counts_QA]->SetTitle("");
    h_mCounts[KEY_Counts_QA]->SetStats(0);

    h_mCounts[KEY_Counts_QA]->GetXaxis()->SetTitle("#phi-#Psi_{3}");
    h_mCounts[KEY_Counts_QA]->GetXaxis()->SetTitleSize(0.08);
    h_mCounts[KEY_Counts_QA]->GetXaxis()->CenterTitle();
    h_mCounts[KEY_Counts_QA]->GetXaxis()->SetTitleOffset(0.8);
    h_mCounts[KEY_Counts_QA]->SetNdivisions(505,"X");

    h_mCounts[KEY_Counts_QA]->GetYaxis()->SetTitle("counts/resolution");
    h_mCounts[KEY_Counts_QA]->GetYaxis()->SetTitleSize(0.08);
    h_mCounts[KEY_Counts_QA]->GetYaxis()->CenterTitle();
    h_mCounts[KEY_Counts_QA]->GetYaxis()->SetRangeUser(0.8*h_mCounts[KEY_Counts_QA]->GetMinimum(),1.2*h_mCounts[KEY_Counts_QA]->GetMaximum());
    h_mCounts[KEY_Counts_QA]->SetNdivisions(505,"Y");

    h_mCounts[KEY_Counts_QA]->SetMarkerStyle(24);
    h_mCounts[KEY_Counts_QA]->SetMarkerColor(kAzure-2);
    h_mCounts[KEY_Counts_QA]->SetMarkerSize(1.4);
    h_mCounts[KEY_Counts_QA]->Draw("pE");
    TF1 *f_flowCounts_QA = new TF1("f_flowCounts_QA",flow_3,0.0,PI_max[mOrder],2);
    f_flowCounts_QA->FixParameter(0,ParFlow[KEY_Counts_QA][0]);
    f_flowCounts_QA->FixParameter(1,ParFlow[KEY_Counts_QA][1]);
    f_flowCounts_QA->SetLineColor(kAzure-2);
    f_flowCounts_QA->SetLineStyle(2);
    f_flowCounts_QA->SetLineWidth(2);
    f_flowCounts_QA->Draw("l same");

    TString KEY_Inte_QA = Form("%s_Inte_pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_%s_SysError_%d_Proj_%d",PID[mPID].Data(),i_pt,Cent_start,mCharge,Eta_start,Order[mOrder].Data(),Sys_QA,Proj_QA); 
    h_mCounts[KEY_Inte_QA]->SetMarkerStyle(20);
    h_mCounts[KEY_Inte_QA]->SetMarkerColor(kGray+2);
    h_mCounts[KEY_Inte_QA]->SetMarkerSize(1.2);
    h_mCounts[KEY_Inte_QA]->Draw("pE same");
    TF1 *f_flowInte_QA = new TF1("f_flowInte_QA",flow_3,0.0,PI_max[mOrder],2);
    f_flowInte_QA->FixParameter(0,ParFlow[KEY_Inte_QA][0]);
    f_flowInte_QA->FixParameter(1,ParFlow[KEY_Inte_QA][1]);
    f_flowInte_QA->SetLineColor(kGray+2);
    f_flowInte_QA->SetLineStyle(1);
    f_flowInte_QA->SetLineWidth(2);
    f_flowInte_QA->Draw("l same");

    TLegend *leg_Counts_QA = new TLegend(0.3,0.6,0.8,0.8);
    leg_Counts_QA->SetBorderSize(0);
    leg_Counts_QA->SetFillColor(10);
    leg_Counts_QA->AddEntry(h_mCounts[KEY_Counts_QA],"bin counting","P");
    leg_Counts_QA->AddEntry(h_mCounts[KEY_Inte_QA],"student-t integrating","P");
    leg_Counts_QA->Draw("same");

    TString pt_range = Form("[%1.1f,%1.1f]",pt_low[i_pt],pt_up[i_pt]);
    plotTopLegend((char*)pt_range.Data(),0.6,0.3,0.08,1,0.0,42,1);
  }
#endif

#if _PlotQA_
  TCanvas *c_RawFlow = new TCanvas("c_RawFlow","c_RawFlow",10,10,800,800);
  c_RawFlow->cd();
  c_RawFlow->cd()->SetLeftMargin(0.15);
  c_RawFlow->cd()->SetBottomMargin(0.15);
  c_RawFlow->cd()->SetTicks(1,1);
  c_RawFlow->cd()->SetGrid(0,0);
  TH1F *h_frame_QA = new TH1F("h_frame_QA","h_frame_QA",100,-0.05,9.95);
  for(Int_t i_bin = 0; i_bin < 100; ++i_bin)
  {
    h_frame_QA->SetBinContent(i_bin+1,-10);
    h_frame_QA->SetBinError(i_bin+1,1.0);
  }
  h_frame_QA->SetTitle("");
  h_frame_QA->SetStats(0);

  h_frame_QA->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_frame_QA->GetXaxis()->SetTitleSize(0.06);
  h_frame_QA->GetXaxis()->CenterTitle();
  h_frame_QA->GetXaxis()->SetRangeUser(-0.05,5.05);
  h_frame_QA->SetNdivisions(505,"X");

  h_frame_QA->GetYaxis()->SetTitle("v_{3}^{raw}");
  h_frame_QA->GetYaxis()->SetTitleSize(0.06);
  h_frame_QA->GetYaxis()->CenterTitle();
  h_frame_QA->GetYaxis()->SetRangeUser(-0.005,0.05);
  h_frame_QA->SetNdivisions(505,"Y");
  h_frame_QA->Draw("pE");

  TString KEY_Flow_Counts_QA = Form("%s_Counts_Centrality_%d_Charge_%d_EtaGap_%d_%s_SysError_%d_Proj_%d",PID[mPID].Data(),Cent_start,mCharge,Eta_start,Order[mOrder].Data(),Sys_QA,Proj_QA); 
  Draw_TGAE_new_Symbol(g_mRawFlow[KEY_Flow_Counts_QA],24,kAzure-2,1.4);

  TString KEY_Flow_Inte_QA= Form("%s_Inte_Centrality_%d_Charge_%d_EtaGap_%d_%s_SysError_%d_Proj_%d",PID[mPID].Data(),Cent_start,mCharge,Eta_start,Order[mOrder].Data(),Sys_QA,Proj_QA); 
  Draw_TGAE_new_Symbol(g_mRawFlow[KEY_Flow_Inte_QA],24,kGray+2,1.2);
#endif

  // resolution correction
  TString InPutFile_Yield = Form("./Data/AuAu%s/nSigmaPion/Yield_%s.root",Energy[mEnergy].Data(),Energy[mEnergy].Data());
  TFile *File_Yield = TFile::Open(InPutFile_Yield.Data());

  TH1FMap h_mMass2_Yields;
  vecFMap Yields_Counts;
  for(Int_t i_cent = 0; i_cent < 9; i_cent++) // centrality bin
  {
    for(Int_t i_charge = mCharge; i_charge < mCharge+1; i_charge++) // charge bin
    {
      for(Int_t i_eta = Eta_start; i_eta < Eta_stop; i_eta++) // eta gap bin
      {
	for(Int_t i_sys = 0; i_sys < 6; i_sys++)
	{
	  TString KEY_PiK_Yield = Form("Centrality_%d_Charge_%d_EtaGap_%d_Yields_PiK_SysError_%d",i_cent,i_charge,i_eta,i_sys);
	  TString HistName = KEY_PiK_Yield+"_projX";
	  TH2F *h_Yields_temp = (TH2F*)File_Yield->Get(KEY_PiK_Yield.Data());
	  h_mMass2_Yields[KEY_PiK_Yield] = (TH1F*)h_Yields_temp->ProjectionX(HistName.Data());
	  Int_t bin_start = h_mMass2_Yields[KEY_PiK_Yield]->FindBin(Inte_start[mPID]);
	  Int_t bin_stop  = h_mMass2_Yields[KEY_PiK_Yield]->FindBin(Inte_stop[mPID]);
	  Float_t counts = 0.0;
	  Float_t errors = 0.0;
	  for(Int_t i_bin = bin_start; i_bin < bin_stop; i_bin++)
	  {
	    counts += h_mMass2_Yields[KEY_PiK_Yield]->GetBinContent(i_bin);
	    errors += h_mMass2_Yields[KEY_PiK_Yield]->GetBinError(i_bin)*h_mMass2_Yields[KEY_PiK_Yield]->GetBinError(i_bin);
	  }
	  Yields_Counts[KEY_PiK_Yield].clear();
	  Yields_Counts[KEY_PiK_Yield].push_back(static_cast<Float_t>(counts));
	  Yields_Counts[KEY_PiK_Yield].push_back(static_cast<Float_t>(TMath::Sqrt(errors)));
	}
      }
    }
  }

#if _PlotQA_
  TCanvas *c_Yields = new TCanvas("c_Yields","c_Yields",10,10,900,900);
  c_Yields->Divide(3,3);
  for(Int_t i_cent = 0; i_cent < 9; i_cent++)
  {
    c_Yields->cd(i_cent+1);
    c_Yields->cd(i_cent+1)->SetLeftMargin(0.15);
    c_Yields->cd(i_cent+1)->SetBottomMargin(0.15);
    c_Yields->cd(i_cent+1)->SetTicks(1,1);
    c_Yields->cd(i_cent+1)->SetGrid(0,0);
    TString KEY_PiK_Yield_QA = Form("Centrality_%d_Charge_%d_EtaGap_%d_Yields_PiK_SysError_%d",i_cent,mCharge,Eta_start,Sys_QA);
    h_mMass2_Yields[KEY_PiK_Yield_QA]->SetStats(0);
    h_mMass2_Yields[KEY_PiK_Yield_QA]->SetTitle("");
    h_mMass2_Yields[KEY_PiK_Yield_QA]->GetXaxis()->SetTitle("x(n#sigma#pion,m^{2})");
    h_mMass2_Yields[KEY_PiK_Yield_QA]->GetXaxis()->CenterTitle();
    h_mMass2_Yields[KEY_PiK_Yield_QA]->GetXaxis()->SetTitleSize(0.06);
    h_mMass2_Yields[KEY_PiK_Yield_QA]->GetYaxis()->SetTitle("counts");
    h_mMass2_Yields[KEY_PiK_Yield_QA]->GetYaxis()->CenterTitle();
    h_mMass2_Yields[KEY_PiK_Yield_QA]->GetYaxis()->SetTitleSize(0.06);
    h_mMass2_Yields[KEY_PiK_Yield_QA]->SetMarkerStyle(24);
    h_mMass2_Yields[KEY_PiK_Yield_QA]->SetMarkerSize(0.5);
    h_mMass2_Yields[KEY_PiK_Yield_QA]->SetMarkerColor(kGray+3);
    h_mMass2_Yields[KEY_PiK_Yield_QA]->SetLineColor(1);
    h_mMass2_Yields[KEY_PiK_Yield_QA]->Draw("pE");
    // cout << "i_cent = " << i_cent << ", Yields_Counts = " << Yields_Counts[KEY_PiK_Yield_QA][0] << endl;

    PlotLine(Inte_start[mPID],Inte_start[mPID],0.0,h_mMass2_Yields[KEY_PiK_Yield_QA]->GetMaximum()/2.0,4,2,2);
    PlotLine(Inte_stop[mPID],Inte_stop[mPID],0.0,h_mMass2_Yields[KEY_PiK_Yield_QA]->GetMaximum()/2.0,4,2,2);
  }
#endif

  // calculate final resolution correction factors and correct flow
  TString InPutFile_Res = Form("./Data/AuAu%s/file_%s_Resolution.root",Energy[mEnergy].Data(),Energy[mEnergy].Data());
  TFile *File_Res = TFile::Open(InPutFile_Res.Data());
  TString Res_Order[2] = {"Res2","Res3"};
  TProMap p_mRes;
  vecFMap ResValue;
  TGraMap g_mFlow;
  for(Int_t i_eta = Eta_start; i_eta < Eta_stop; i_eta++) // eta_gap loop
  {
    TString KEY_eta = Form("%s_EtaGap_%d_EP",Res_Order[mOrder].Data(),i_eta);
    p_mRes[KEY_eta] = (TProfile*)File_Res->Get(KEY_eta.Data()); // read in resolution file
    for(Int_t i_charge = mCharge; i_charge < mCharge+1; i_charge++) // charge bin
    {
      for(Int_t i_sys = Sys_start; i_sys < Sys_stop; i_sys++) // Systematic errors loop
      {
	for(Int_t i_cent = Cent_start; i_cent < Cent_stop; i_cent++) // centrality bin loop
	{
	  Float_t yields_total = 0.0;
	  for(Int_t cent = cent_low[i_cent]; cent <= cent_up[i_cent]; cent++) // calculate resolution and total yields in selected centrality bin
	  {
	    if(p_mRes[KEY_eta]->GetBinContent(cent+1) > 0) 
	    {
	      TString KEY_PiK_Yield = Form("Centrality_%d_Charge_%d_EtaGap_%d_Yields_PiK_SysError_%d",cent,i_charge,i_eta,i_sys);
	      ResValue[KEY_PiK_Yield].push_back(static_cast<Float_t>(TMath::Sqrt(p_mRes[KEY_eta]->GetBinContent(cent+1))));
	      yields_total += Yields_Counts[KEY_PiK_Yield][0];
	    }
	  }

	  Float_t mean_res = 0.0;
	  for(Int_t cent = cent_low[i_cent]; cent <= cent_up[i_cent]; cent++) // calculate final resolution correction factor <1/R(centrality)>
	  {
	    TString KEY_PiK_Yield = Form("Centrality_%d_Charge_%d_EtaGap_%d_Yields_PiK_SysError_%d",cent,i_charge,i_eta,i_sys);
	    mean_res += Yields_Counts[KEY_PiK_Yield][0]/(ResValue[KEY_PiK_Yield][0]*yields_total);
	  }
	  cout << "i_eta = " << i_eta << ", i_charge = " << i_charge << ", i_sys = " << i_sys << ", centrality_bin = " << i_cent << ", mean_res = " << mean_res << endl;

	  for(Int_t i_proj = Proj_start; i_proj < Proj_stop; ++i_proj) // proj range loop
	  { 
	    TString KEY_Flow_Counts = Form("%s_Counts_Centrality_%d_Charge_%d_EtaGap_%d_%s_SysError_%d_Proj_%d",PID[mPID].Data(),i_cent,i_charge,i_eta,Order[mOrder].Data(),i_sys,i_proj); 
	    g_mFlow[KEY_Flow_Counts] = new TGraphAsymmErrors();
	    for(Int_t i_point = 0; i_point < g_mRawFlow[KEY_Flow_Counts]->GetN(); ++i_point)
	    {
	      Double_t pt,v3Raw, err_v3Raw;
	      g_mRawFlow[KEY_Flow_Counts]->GetPoint(i_point,pt,v3Raw);
	      err_v3Raw = g_mRawFlow[KEY_Flow_Counts]->GetErrorYhigh(i_point);

	      Double_t v3 = v3Raw*mean_res;
	      Double_t err_v3 = ErrTimes(v3Raw,mean_res,err_v3Raw,0.0);
	      g_mFlow[KEY_Flow_Counts]->SetPoint(i_point,pt,v3);
	      g_mFlow[KEY_Flow_Counts]->SetPointError(i_point,0.0,0.0,err_v3,err_v3);
	    }

	    TString KEY_Flow_Inte = Form("%s_Inte_Centrality_%d_Charge_%d_EtaGap_%d_%s_SysError_%d_Proj_%d",PID[mPID].Data(),i_cent,i_charge,i_eta,Order[mOrder].Data(),i_sys,i_proj); 
	    g_mFlow[KEY_Flow_Inte] = new TGraphAsymmErrors();
	    for(Int_t i_point = 0; i_point < g_mRawFlow[KEY_Flow_Inte]->GetN(); ++i_point)
	    {
	      Double_t pt,v3Raw, err_v3Raw;
	      g_mRawFlow[KEY_Flow_Inte]->GetPoint(i_point,pt,v3Raw);
	      err_v3Raw = g_mRawFlow[KEY_Flow_Inte]->GetErrorYhigh(i_point);

	      Double_t v3 = v3Raw*mean_res;
	      Double_t err_v3 = ErrTimes(v3Raw,mean_res,err_v3Raw,0.0);
	      g_mFlow[KEY_Flow_Inte]->SetPoint(i_point,pt,v3);
	      g_mFlow[KEY_Flow_Inte]->SetPointError(i_point,0.0,0.0,err_v3,err_v3);
	    }
	  }
	}
      }
    }
  }

  TCanvas *c_v3 = new TCanvas("c_v3","c_v3",10,10,800,800);
  c_v3->cd();
  c_v3->cd()->SetLeftMargin(0.15);
  c_v3->cd()->SetBottomMargin(0.15);
  c_v3->cd()->SetTicks(1,1);
  c_v3->cd()->SetGrid(0,0);

  TH1F *h_frame = new TH1F("h_frame","h_frame",100,-0.05,9.95);
  for(Int_t i_bin = 0; i_bin < 100; ++i_bin)
  {
    h_frame->SetBinContent(i_bin+1,-10);
    h_frame->SetBinError(i_bin+1,1.0);
  }
  h_frame->SetTitle("");
  h_frame->SetStats(0);

  h_frame->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_frame->GetXaxis()->SetTitleSize(0.06);
  h_frame->GetXaxis()->CenterTitle();
  h_frame->GetXaxis()->SetRangeUser(-0.05,5.05);
  h_frame->SetNdivisions(505,"X");

  h_frame->GetYaxis()->SetTitle("v_{3}");
  h_frame->GetYaxis()->SetTitleSize(0.06);
  h_frame->GetYaxis()->CenterTitle();
  h_frame->GetYaxis()->SetRangeUser(-0.005,0.15);
  h_frame->SetNdivisions(505,"Y");
  h_frame->Draw("hE");

  for(Int_t i_cent = Cent_start; i_cent < Cent_stop; i_cent++) // centrality loop
  {
    for(Int_t i_charge = mCharge; i_charge < mCharge+1; i_charge++) // charge loop 
    {
      for(Int_t i_eta = Eta_start; i_eta < Eta_stop; i_eta++) // eta gap loop
      {
	for(Int_t i_sys = Sys_start; i_sys < Sys_stop; i_sys++) // sys cut loop
	{
	  for(Int_t i_proj = Proj_start; i_proj < Proj_stop; ++i_proj) // proj range loop
	  { 
	    TString KEY_Flow_Counts = Form("%s_Counts_Centrality_%d_Charge_%d_EtaGap_%d_%s_SysError_%d_Proj_%d",PID[mPID].Data(),i_cent,i_charge,i_eta,Order[mOrder].Data(),i_sys,i_proj); 
	    g_mFlow[KEY_Flow_Counts]->SetMarkerStyle(24);
	    g_mFlow[KEY_Flow_Counts]->SetMarkerColor(kAzure-2);
	    g_mFlow[KEY_Flow_Counts]->SetMarkerSize(1.4);
	    g_mFlow[KEY_Flow_Counts]->SetName(KEY_Flow_Counts.Data());
	    g_mFlow[KEY_Flow_Counts]->Draw("pE same");

	    TString KEY_Flow_Inte = Form("%s_Inte_Centrality_%d_Charge_%d_EtaGap_%d_%s_SysError_%d_Proj_%d",PID[mPID].Data(),i_cent,i_charge,i_eta,Order[mOrder].Data(),i_sys,i_proj); 
	    g_mFlow[KEY_Flow_Inte]->SetMarkerStyle(20);
	    g_mFlow[KEY_Flow_Inte]->SetMarkerColor(kGray+2);
	    g_mFlow[KEY_Flow_Inte]->SetMarkerSize(1.4);
	    g_mFlow[KEY_Flow_Inte]->SetName(KEY_Flow_Inte.Data());
	    g_mFlow[KEY_Flow_Inte]->Draw("pE same");
	  }
	}
      }
    }
  }

  TString OutPutFile = Form("./OutPut/AuAu%s/%s/flow_%s_Charge_%d.root",Energy[mEnergy].Data(),PID[mPID].Data(),Order[mOrder].Data(),mCharge);
  cout << "OutPutFile set to: " << OutPutFile.Data() << endl;
  TFile *File_OutPut = new TFile(OutPutFile.Data(),"RECREATE");
  File_OutPut->cd();
  for(Int_t i_cent = Cent_start; i_cent < Cent_stop; i_cent++) // centrality loop
  {
    for(Int_t i_charge = mCharge; i_charge < mCharge+1; i_charge++) // charge loop 
    {
      for(Int_t i_eta = Eta_start; i_eta < Eta_stop; i_eta++) // eta gap loop
      {
	for(Int_t i_sys = Sys_start; i_sys < Sys_stop; i_sys++) // sys cut loop
	{
	  for(Int_t i_proj = Proj_start; i_proj < Proj_stop; ++i_proj) // proj range loop
	  { 
	    TString KEY_Flow_Counts = Form("%s_Counts_Centrality_%d_Charge_%d_EtaGap_%d_%s_SysError_%d_Proj_%d",PID[mPID].Data(),i_cent,i_charge,i_eta,Order[mOrder].Data(),i_sys,i_proj); 
	    g_mFlow[KEY_Flow_Counts]->Write();

	    TString KEY_Flow_Inte = Form("%s_Inte_Centrality_%d_Charge_%d_EtaGap_%d_%s_SysError_%d_Proj_%d",PID[mPID].Data(),i_cent,i_charge,i_eta,Order[mOrder].Data(),i_sys,i_proj); 
	    g_mFlow[KEY_Flow_Inte]->Write();
	  }
	}
      }
    }
  }
  h_frame->Write();
  File_OutPut->Close();

}

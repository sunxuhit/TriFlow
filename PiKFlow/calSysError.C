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

static const TString PID[2] = {"Pion","Kaon"};
static const TString Energy[2] = {"200GeV","39GeV"};
static const TString Order[2] = {"2nd","3rd"};
static const TString Method[2] = {"Counts","Inte"};
static const TString ParType[2] = {"pion","kaon"};
static const TString Charge[2] = {"plus","minus"};

static const Int_t Cent_total = 4; 
static const Int_t Cent_start = 0;
static const Int_t Cent_stop  = 1;

static const Int_t Eta_total = 4;
static const Int_t Eta_start = 0;
static const Int_t Eta_stop = 1;

static const Int_t Sys_total = 6;
static const Int_t Sys_start = 0;
static const Int_t Sys_stop  = 6;

static const Int_t Proj_total = 3;
static const Int_t Proj_start = 0;
static const Int_t Proj_stop  = 3;

typedef std::map<TString,TH1F*> TH1FMap;
typedef std::map<TString,TProfile*> TProMap;
typedef std::map<TString,TGraphAsymmErrors*> TGraMap;
typedef std::map<TString,std::vector<Float_t>> vecFMap;

using namespace std;

#ifndef _PlotQA_
#define _PlotQA_ 1
#endif

void calSysError(Int_t mEnergy = 0, Int_t mCharge = 0, Int_t mOrder = 1, Int_t mPID = 0)
{
  TString InPutFile = Form("./OutPut/AuAu%s/%s/flow_%s_Charge_%d.root",Energy[mEnergy].Data(),PID[mPID].Data(),Order[mOrder].Data(),mCharge);
  cout << "InPutFile set to: " << InPutFile.Data() << endl;
  TFile *File_InPut = TFile::Open(InPutFile.Data());

  TGraMap g_mFlow;
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
	    for(Int_t i_method = 0; i_method < 2; ++i_method)
	    {
	      TString KEY_Flow = Form("%s_%s_Centrality_%d_Charge_%d_EtaGap_%d_%s_SysError_%d_Proj_%d",PID[mPID].Data(),Method[i_method].Data(),i_cent,i_charge,i_eta,Order[mOrder].Data(),i_sys,i_proj); 
	      g_mFlow[KEY_Flow] = (TGraphAsymmErrors*)File_InPut->Get(KEY_Flow.Data());
	    }
	  }
	}
      }
    }
  }
  TH1F *h_frame = (TH1F*)File_InPut->Get("h_frame");

#if _PlotQA_
  TCanvas *c_v3 = new TCanvas("c_v3","c_v3",10,10,800,800);
  c_v3->cd();
  c_v3->cd()->SetLeftMargin(0.15);
  c_v3->cd()->SetBottomMargin(0.15);
  c_v3->cd()->SetTicks(1,1);
  c_v3->cd()->SetGrid(0,0);
  Float_t pt_start = -0.05;
  Float_t pt_stop  = 3.8;
  h_frame->GetXaxis()->SetRangeUser(pt_start,pt_stop);
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
	    for(Int_t i_method = 0; i_method < 2; ++i_method)
	    {
	      TString KEY_Flow = Form("%s_%s_Centrality_%d_Charge_%d_EtaGap_%d_%s_SysError_%d_Proj_%d",PID[mPID].Data(),Method[i_method].Data(),i_cent,i_charge,i_eta,Order[mOrder].Data(),i_sys,i_proj); 
	      g_mFlow[KEY_Flow]->Draw("pE same");
	    }
	  }
	}
      }
    }
  }
  PlotLine(pt_start,pt_stop,0.0,0.0,1,2,2);
#endif

  TString KEY_Default = Form("%s_Counts_Centrality_0_Charge_%d_EtaGap_0_%s_SysError_0_Proj_0",PID[mPID].Data(),mCharge,Order[mOrder].Data()); 

  TGraphAsymmErrors *g_SysErrors = new TGraphAsymmErrors();
  for(Int_t i_point = 0; i_point < g_mFlow[KEY_Default]->GetN(); ++i_point)
  {
    Int_t counter = 0;
    TGraph *g_diff = new TGraph();
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
	      for(Int_t i_method = 0; i_method < 2; ++i_method)
	      {
		TString KEY_Flow = Form("%s_%s_Centrality_%d_Charge_%d_EtaGap_%d_%s_SysError_%d_Proj_%d",PID[mPID].Data(),Method[i_method].Data(),i_cent,i_charge,i_eta,Order[mOrder].Data(),i_sys,i_proj); 
		Double_t pt_val, v3_val;
		g_mFlow[KEY_Flow]->GetPoint(i_point,pt_val,v3_val);
		if(i_method == 0 || (i_method == 1 && i_point > 4))
		{
		  g_diff->SetPoint(counter,pt_val-0.1,v3_val);
		  counter++;
		}
	      }
	    }
	  }
	}
      }
    }
    Double_t v3_min = TMath::MinElement(g_diff->GetN(),g_diff->GetY());
    Double_t v3_max = TMath::MaxElement(g_diff->GetN(),g_diff->GetY());
    Double_t SysError_v3 = (v3_max-v3_min)/TMath::Sqrt(12.0);
    cout << "v3_min = " << v3_min << ", v3_max = " << v3_max << endl;
    Double_t pt, v3;
    g_mFlow[KEY_Default]->GetPoint(i_point,pt,v3);
    g_SysErrors->SetPoint(i_point,pt,v3);
    g_SysErrors->SetPointError(i_point,0.0,0.0,SysError_v3,SysError_v3);
    cout << "number of combinations is: " << counter << endl;
    delete g_diff;
  }

  TCanvas *c_v3_SysError = new TCanvas("c_v3_SysError","c_v3_SysError",600,10,800,800);
  c_v3_SysError->cd();
  c_v3_SysError->cd()->SetLeftMargin(0.15);
  c_v3_SysError->cd()->SetBottomMargin(0.15);
  c_v3_SysError->cd()->SetTicks(1,1);
  c_v3_SysError->cd()->SetGrid(0,0);
  h_frame->Draw("pE");
  g_mFlow[KEY_Default]->Draw("pE same");
  g_SysErrors->SetMarkerStyle(20);
  g_SysErrors->SetMarkerColor(kGray+2);
  g_SysErrors->Draw("pE3 same");

  TString OutPutFile = Form("./OutPut/AuAu%s/%s/TriFlow_%s_%s_SysError.root",Energy[mEnergy].Data(),PID[mPID].Data(),ParType[mPID].Data(),Charge[mCharge].Data());
  cout << "OutPutFile set to: " << OutPutFile.Data() << endl;
  TFile *File_OutPut = new TFile(OutPutFile.Data(),"RECREATE");
  File_OutPut->cd();
  h_frame->Write();
  TString StatErrorV3 = Form("g_TriFlow_%s_%s_%s_StatError",Energy[mEnergy].Data(),ParType[mPID].Data(),Charge[mCharge].Data());
  g_mFlow[KEY_Default]->SetName(StatErrorV3.Data());
  g_mFlow[KEY_Default]->Write();
  TString SysErrorV3 = Form("g_TriFlow_%s_%s_%s_SysError",Energy[mEnergy].Data(),ParType[mPID].Data(),Charge[mCharge].Data());
  g_SysErrors->SetName(SysErrorV3.Data());
  g_SysErrors->Write();
  File_OutPut->Close();
}

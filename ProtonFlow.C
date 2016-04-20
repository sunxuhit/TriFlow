#include "TString.h"
#include "TFile.h"
#include "TH1F.h"
#include <map>
#include "TCanvas.h"
#include <iostream>
#include "draw.h"
#include "TF1.h"
#include "functions.h"
#include <vector>
#include "TMath.h"
#include "TLegend.h"
#include "TGaxis.h"
#include "TProfile.h"
#include "TGraphAsymmErrors.h"
#include "./student_t_1d_single.h"

// mEnergy: 0, 200 GeV | 1, 39 GeV
static const TString Energy[2] = {"200GeV","39GeV"};
static const TString Order[2] = {"2nd","3rd"};
static const Double_t PI_max[2] = {TMath::Pi()/2.0,TMath::Pi()/3.0};
static const Float_t nSigProton = 3.0;
static const Float_t Flow_Order[2] = {2.0,3.0};

static const Int_t pt_total = 16; // pT loop
static const Int_t pt_start = 0;
static const Int_t pt_stop  = 16;

// pt bin
//                                      0 , 1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 ,10 ,11 ,12 ,13 ,14 ,15
static const Float_t pt_low_raw[16] = {0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4};
static const Float_t pt_up_raw[16]  = {0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,4.2};

// pt rebin
static const Int_t pt_rebin_total = 15;
static const Float_t pt_low[pt_rebin_total] = {0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.4};
static const Float_t pt_up[pt_rebin_total]  = {0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.4,4.2};
static const Int_t pt_rebin_start[pt_rebin_total] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,15};
static const Int_t pt_rebin_stop[pt_rebin_total]  = {0,1,2,3,4,5,6,7,8,9,10,11,12,14,15};
static const Int_t pt_rebin_first = 0;
static const Int_t pt_rebin_last  = 15;
static const Int_t pt_QA    = 4;
static const Float_t m2_cut_low[pt_rebin_total] = {0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.65,0.7,0.7,0.75,0.879,0.879,0.879,0.879};
static const Float_t m2_cut_up[pt_rebin_total]  = {1.2,1.2,1.2,1.2,1.2,1.2,1.2,1.20,1.2,1.2,1.20,1.400,1.400,1.400,1.400};
static const Float_t m2_yield_low  = 0.6;
static const Float_t m2_yield_high = 1.2;

static const Int_t Cent_total = 4; // Centrality loop
static const Int_t Cent_start = 0;
static const Int_t Cent_stop  = 1;
static const Int_t cent_low[4] = {0,7,4,0}; // 0 = 0-80%, 1 = 0-10%, 2 = 10-40%, 3 = 40-80%
static const Int_t cent_up[4]  = {8,8,6,3}; // 0 = 0-80%, 1 = 0-10%, 2 = 10-40%, 3 = 40-80%

static const Int_t Charge_total = 2;
static const Int_t Charge_start = 0;
static const Int_t Charge_stop  = 2;

static const Int_t Eta_total = 4; // Eta loop
static const Int_t Eta_start = 0;
static const Int_t Eta_stop  = 1;

static const Int_t phi_total = 7; // phi loop
static const Int_t phi_start = 0;
static const Int_t phi_stop  = 7;

static const Int_t Sys_total = 18; // Systematic loop
static const Int_t Sys_start = 0;
static const Int_t Sys_stop  = 18;

typedef std::map<TString,TH1F*> TH1FMap;
typedef std::map<TString,TProfile*> TProMap;
typedef std::map<TString,TGraphAsymmErrors*> TGraMap;
typedef std::map<TString,std::vector<Float_t>> vecFMap;

void ProtonFlow(Int_t mEnergy = 0, Int_t mOrder = 1)
{
  TGaxis::SetMaxDigits(4);

  TString InPutFile_Flow = Form("./Data/AuAu%s/Proton/Flow_%s.root",Energy[mEnergy].Data(),Energy[mEnergy].Data());
  TFile *File_Flow = TFile::Open(InPutFile_Flow.Data());

  // read in  histogram for flow calculation
  TH1FMap h_mMass2_raw;
  for(Int_t i_pt = 0; i_pt < pt_total; i_pt++) // pt bin
  {
    for(Int_t i_cent = Cent_start; i_cent < Cent_stop; i_cent++) // centrality bin
    {
      for(Int_t i_charge = Charge_start; i_charge < Charge_stop; i_charge++) // charge bin
      {
	for(Int_t i_eta = Eta_start; i_eta < Eta_stop; i_eta++) // eta gap bin
	{
	  for(Int_t i_phi = phi_start; i_phi < phi_stop; i_phi++) // phi-psi bin
	  {
	    for(Int_t i_sys = 0; i_sys < 18; i_sys++)
	    {
	      TString KEY_Proton_raw = Form("pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_phi_Psi_%d_%s_Proton_SysError_%d",i_pt,i_cent,i_charge,i_eta,i_phi,Order[mOrder].Data(),i_sys);
	      h_mMass2_raw[KEY_Proton_raw] = (TH1F*)File_Flow->Get(KEY_Proton_raw.Data());
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
    TString KEY_Proton_raw_pT_QA = Form("pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_phi_Psi_%d_%s_Proton_SysError_%d",i_pt,Cent_start,Charge_start,Eta_start,phi_start,Order[mOrder].Data(),Sys_start);
    h_mMass2_raw[KEY_Proton_raw_pT_QA]->SetStats(0);
    h_mMass2_raw[KEY_Proton_raw_pT_QA]->SetTitle("");
    h_mMass2_raw[KEY_Proton_raw_pT_QA]->GetXaxis()->SetTitle("m^{2} ((GeV/c^{2})^{2})");
    h_mMass2_raw[KEY_Proton_raw_pT_QA]->GetXaxis()->CenterTitle();
    h_mMass2_raw[KEY_Proton_raw_pT_QA]->GetXaxis()->SetTitleSize(0.06);
    h_mMass2_raw[KEY_Proton_raw_pT_QA]->GetYaxis()->SetTitle("Counts/Resolution");
    h_mMass2_raw[KEY_Proton_raw_pT_QA]->GetYaxis()->CenterTitle();
    h_mMass2_raw[KEY_Proton_raw_pT_QA]->GetYaxis()->SetTitleSize(0.06);
    h_mMass2_raw[KEY_Proton_raw_pT_QA]->Draw("pE");

    TString pT_range = Form("[%.2f,%.2f]",pt_low_raw[i_pt],pt_up_raw[i_pt]);
    plotTopLegend((char*)pT_range.Data(),0.65,0.7,0.08,1,0.0,42,1);
  }
  */

  // pT rebin
  TH1FMap h_mMass2; // rebinned m2 distribution
  for(Int_t i_cent = Cent_start; i_cent < Cent_stop; i_cent++) // Centrality loop
  {
    for(Int_t i_charge = Charge_start; i_charge < Charge_stop; i_charge++) // charge bin
    {
      for(Int_t i_eta = Eta_start; i_eta < Eta_stop; i_eta++) // EtaGap loop
      {
	for(Int_t i_phi = phi_start; i_phi < phi_stop; i_phi++) // phi loop
	{
	  for(Int_t i_sys = Sys_start; i_sys < Sys_stop; i_sys++) // Systematic loop
	  {
	    for(Int_t pt_bin = pt_rebin_first; pt_bin < pt_rebin_last; pt_bin++) // pt loop
	    {
	      TString KEY_Proton = Form("pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_phi_Psi_%d_%s_Proton_SysError_%d",pt_bin,i_cent,i_charge,i_eta,i_phi,Order[mOrder].Data(),i_sys);
	      for(Int_t i_pt = pt_rebin_start[pt_bin]; i_pt <= pt_rebin_stop[pt_bin]; i_pt++)
	      {
		TString KEY_Proton_raw = Form("pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_phi_Psi_%d_%s_Proton_SysError_%d",i_pt,i_cent,i_charge,i_eta,i_phi,Order[mOrder].Data(),i_sys);
		if(i_pt == pt_rebin_start[pt_bin]) h_mMass2[KEY_Proton] = (TH1F*)h_mMass2_raw[KEY_Proton_raw]->Clone();
		else h_mMass2[KEY_Proton]->Add(h_mMass2_raw[KEY_Proton_raw],1.0);
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
    TString KEY_Proton_raw_pT_QA = Form("pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_phi_Psi_%d_%s_Proton_SysError_%d",i_pt,Cent_start,Charge_start,Eta_start,phi_start,Order[mOrder].Data(),Sys_start);
    h_mMass2[KEY_Proton_raw_pT_QA]->SetStats(0);
    h_mMass2[KEY_Proton_raw_pT_QA]->SetTitle("");
    h_mMass2[KEY_Proton_raw_pT_QA]->GetXaxis()->SetTitle("m^{2} ((GeV/c^{2})^{2})");
    h_mMass2[KEY_Proton_raw_pT_QA]->GetXaxis()->CenterTitle();
    h_mMass2[KEY_Proton_raw_pT_QA]->GetXaxis()->SetTitleSize(0.06);
    h_mMass2[KEY_Proton_raw_pT_QA]->GetYaxis()->SetTitle("Counts/Resolution");
    h_mMass2[KEY_Proton_raw_pT_QA]->GetYaxis()->CenterTitle();
    h_mMass2[KEY_Proton_raw_pT_QA]->GetYaxis()->SetTitleSize(0.06);
    h_mMass2[KEY_Proton_raw_pT_QA]->Draw("pE");

    TString pT_range = Form("[%.2f,%.2f]",pt_low[i_pt],pt_up[i_pt]);
    plotTopLegend((char*)pT_range.Data(),0.65,0.7,0.08,1,0.0,42,1);
  }
  */

  /*
  //QA: h_mMass2_total vs rebinned pT and 1d student-t fits
  // integrate over phi and do the student-t fit to get fitting/counting range
  TH1FMap h_mMass2_total;
  vecFMap ParStudent;
  for(Int_t i_pt = pt_rebin_first; i_pt < pt_rebin_last; i_pt++) // pt bin
  {
    for(Int_t i_cent = Cent_start; i_cent < Cent_stop; i_cent++) // centrality bin
    {
      for(Int_t i_charge = Charge_start; i_charge < Charge_stop; i_charge++) // charge bin
      {
	for(Int_t i_eta = Eta_start; i_eta < Eta_stop; i_eta++) // eta gap bin
	{
	  for(Int_t i_sys = Sys_start; i_sys < Sys_stop; i_sys++)
	  {
	    TString KEY_Proton_total = Form("pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_%s_Proton_SysError_%d",i_pt,i_cent,i_charge,i_eta,Order[mOrder].Data(),i_sys);
	    for(Int_t i_phi = phi_start; i_phi < phi_stop; i_phi++) // phi-psi bin
	    {
	      TString KEY_Proton = Form("pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_phi_Psi_%d_%s_Proton_SysError_%d",i_pt,i_cent,i_charge,i_eta,i_phi,Order[mOrder].Data(),i_sys);
	      if(i_phi == phi_start) h_mMass2_total[KEY_Proton_total] = (TH1F*)h_mMass2[KEY_Proton]->Clone();
	      else h_mMass2_total[KEY_Proton_total]->Add(h_mMass2[KEY_Proton],1.0);
	    }

#if 0
	    // 1d student-t fits for yields
	    TF1 *f_student = new TF1("f_student",student_t_1d_single,0.5,1.5,4);
	    f_student->SetParameter(0,7.5);
	    f_student->SetParameter(1,0.938*0.938);
	    f_student->SetParameter(2,0.04);
	    f_student->SetParameter(3,100000);
	    f_student->SetRange(0.8,1.1);
	    h_mMass2_total[KEY_Proton_total]->Fit(f_student,"MQRN");

	    TF1 *f_student_2nd = new TF1("f_student_2nd",student_t_1d_single,0.5,1.5,4);
	    f_student_2nd->SetParameter(0,f_student->GetParameter(0));
	    f_student_2nd->SetParameter(1,f_student->GetParameter(1));
	    f_student_2nd->SetParameter(2,f_student->GetParameter(2));
	    f_student_2nd->SetParameter(3,f_student->GetParameter(3));
	    f_student_2nd->SetRange(f_student->GetParameter(1),f_student->GetParameter(1)+nSigProton*f_student->GetParameter(2));
	    cout << "i_pt = " << i_pt << ",i_charge = " << i_charge << ", i_sys = " << i_sys << endl;
	    h_mMass2_total[KEY_Proton_total]->Fit(f_student_2nd,"MRN");
	    for(Int_t i_par = 0; i_par < 4; i_par++) ParStudent[KEY_Proton_total].push_back(static_cast<Float_t>(f_student_2nd->GetParameter(i_par)));
#endif
	  }
	}
      }
    }
  }

  TCanvas *c_pT_rebin_total = new TCanvas("c_pT_rebin_total","c_pT_rebin_total",10,10,1600,1600);
  c_pT_rebin_total->Divide(4,4);
  for(Int_t i_pt = pt_rebin_first; i_pt < pt_rebin_last; i_pt++)
  {
    c_pT_rebin_total->cd(pt_rebin_start[i_pt]+1);
    c_pT_rebin_total->cd(pt_rebin_start[i_pt]+1)->SetLeftMargin(0.15);
    c_pT_rebin_total->cd(pt_rebin_start[i_pt]+1)->SetBottomMargin(0.15);
    c_pT_rebin_total->cd(pt_rebin_start[i_pt]+1)->SetTicks(1,1);
    c_pT_rebin_total->cd(pt_rebin_start[i_pt]+1)->SetGrid(0,0);
    TString KEY_Proton_total_pT_QA = Form("pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_%s_Proton_SysError_%d",i_pt,Cent_start,Charge_start,Eta_start,Order[mOrder].Data(),Sys_start);
    h_mMass2_total[KEY_Proton_total_pT_QA]->SetStats(0);
    h_mMass2_total[KEY_Proton_total_pT_QA]->SetTitle("");
    h_mMass2_total[KEY_Proton_total_pT_QA]->GetXaxis()->SetTitle("m^{2} ((GeV/c^{2})^{2})");
    h_mMass2_total[KEY_Proton_total_pT_QA]->GetXaxis()->CenterTitle();
    h_mMass2_total[KEY_Proton_total_pT_QA]->GetXaxis()->SetTitleSize(0.06);
    h_mMass2_total[KEY_Proton_total_pT_QA]->GetYaxis()->SetTitle("Counts/Resolution");
    h_mMass2_total[KEY_Proton_total_pT_QA]->GetYaxis()->CenterTitle();
    h_mMass2_total[KEY_Proton_total_pT_QA]->GetYaxis()->SetTitleSize(0.06);
    h_mMass2_total[KEY_Proton_total_pT_QA]->SetMarkerStyle(24);
    h_mMass2_total[KEY_Proton_total_pT_QA]->SetMarkerColor(kGray+3);
    h_mMass2_total[KEY_Proton_total_pT_QA]->SetMarkerSize(0.4);
    h_mMass2_total[KEY_Proton_total_pT_QA]->Draw("pE");

    PlotLine(m2_cut_low[i_pt],m2_cut_low[i_pt],0.0,h_mMass2_total[KEY_Proton_total_pT_QA]->GetMaximum()/2.0,4,2,2);
    PlotLine( m2_cut_up[i_pt], m2_cut_up[i_pt],0.0,h_mMass2_total[KEY_Proton_total_pT_QA]->GetMaximum()/2.0,4,2,2);

#if 0
    TF1 *f_student_QA = new TF1("f_student_QA",student_t_1d_single,0.5,1.5,4);
    f_student_QA->SetParameter(0,ParStudent[KEY_Proton_total_pT_QA][0]);
    f_student_QA->SetParameter(1,ParStudent[KEY_Proton_total_pT_QA][1]);
    f_student_QA->SetParameter(2,ParStudent[KEY_Proton_total_pT_QA][2]);
    f_student_QA->SetParameter(3,ParStudent[KEY_Proton_total_pT_QA][3]);
    f_student_QA->SetLineColor(2);
    f_student_QA->SetLineWidth(2);
    f_student_QA->SetRange(ParStudent[KEY_Proton_total_pT_QA][1],ParStudent[KEY_Proton_total_pT_QA][1]+nSigProton*ParStudent[KEY_Proton_total_pT_QA][2]);
    f_student_QA->Draw("l same");
#endif

    TString pT_range = Form("[%.2f,%.2f]",pt_low[i_pt],pt_up[i_pt]);
    plotTopLegend((char*)pT_range.Data(),0.65,0.7,0.08,1,0.0,42,1);
  }
  */

  // get counts for all phi bin
  TH1FMap h_mCounts, h_mRawFlow;
  vecFMap ParFlow;
  for(Int_t i_cent = Cent_start; i_cent < Cent_stop; i_cent++) // Centrality loop
  {
    for(Int_t i_charge = Charge_start; i_charge < Charge_stop; i_charge++) // charge bin
    {
      for(Int_t i_eta = Eta_start; i_eta < Eta_stop; i_eta++) // EtaGap loop
      {
	for(Int_t i_sys = Sys_start; i_sys < Sys_stop; i_sys++) // Systematic loop
	{
	  TString KEY_Proton_RawFlow = Form("RawFlow_Centrality_%d_Charge_%d_EtaGap_%d_%s_Proton_SysError_%d",i_cent,i_charge,i_eta,Order[mOrder].Data(),i_sys);
	  h_mRawFlow[KEY_Proton_RawFlow] = new TH1F(KEY_Proton_RawFlow.Data(),KEY_Proton_RawFlow.Data(),100,-0.05,9.95); // histogram for raw flow
	  for(Int_t i_pt = pt_rebin_first; i_pt < pt_rebin_last; i_pt++)
	  {
	    TString KEY_Proton_Counts = Form("pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_%s_Proton_SysError_%d",i_pt,i_cent,i_charge,i_eta,Order[mOrder].Data(),i_sys);
	    h_mCounts[KEY_Proton_Counts] = new TH1F(KEY_Proton_Counts.Data(),KEY_Proton_Counts.Data(),7,0.0,PI_max[mOrder]); // histogram for counts
	    for(Int_t i_phi = phi_start; i_phi < phi_stop; i_phi++) // phi loop
	    {
	      TString KEY_Proton = Form("pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_phi_Psi_%d_%s_Proton_SysError_%d",i_pt,i_cent,i_charge,i_eta,i_phi,Order[mOrder].Data(),i_sys);
	      Int_t proton_start = h_mMass2[KEY_Proton]->FindBin(m2_cut_low[i_pt]+10e-4);
	      Int_t proton_stop  = h_mMass2[KEY_Proton]->FindBin(m2_cut_up[i_pt]-10e-4);
	      Float_t counts_proton = 0.0;
	      Float_t errors_proton = 0.0;
	      for(Int_t i_bin = proton_start; i_bin < proton_stop; i_bin++)
	      {
		counts_proton += h_mMass2[KEY_Proton]->GetBinContent(i_bin);
		errors_proton += h_mMass2[KEY_Proton]->GetBinError(i_bin)*h_mMass2[KEY_Proton]->GetBinError(i_bin);
	      }
	      Float_t bin_center = PI_max[mOrder]/14.0+i_phi*PI_max[mOrder]/7.0;
	      Int_t bin_counts = h_mCounts[KEY_Proton_Counts]->FindBin(bin_center);
	      h_mCounts[KEY_Proton_Counts]->SetBinContent(bin_counts,counts_proton);
	      h_mCounts[KEY_Proton_Counts]->SetBinError(bin_counts,TMath::Sqrt(errors_proton));
	    }

	    TF1 *f_flow = new TF1("f_flow",flow,0.0,PI_max[mOrder],3); // flow fit
	    f_flow->SetParameter(0,1000);
	    f_flow->SetParameter(1,0.2);
	    f_flow->FixParameter(2,Flow_Order[mOrder]);
	    h_mCounts[KEY_Proton_Counts]->Fit(f_flow,"NQM");
	    ParFlow[KEY_Proton_Counts].clear();
	    ParFlow[KEY_Proton_Counts].push_back(static_cast<Float_t>(f_flow->GetParameter(0)));
	    ParFlow[KEY_Proton_Counts].push_back(static_cast<Float_t>(f_flow->GetParameter(1)));
	    ParFlow[KEY_Proton_Counts].push_back(static_cast<Float_t>(f_flow->GetParameter(2)));

	    Float_t pt_mean = (pt_low[i_pt]+pt_up[i_pt])/2.0;
	    h_mRawFlow[KEY_Proton_RawFlow]->SetBinContent(h_mRawFlow[KEY_Proton_RawFlow]->FindBin(pt_mean),f_flow->GetParameter(1));
	    h_mRawFlow[KEY_Proton_RawFlow]->SetBinError(h_mRawFlow[KEY_Proton_RawFlow]->FindBin(pt_mean),f_flow->GetParError(1));
	  }
	}
      }
    }
  }

  /*
  // QA plots: counts vs phi-Psi
  TCanvas *c_Counts = new TCanvas("c_Counts","c_Counts",10,10,900,900);
  c_Counts->Divide(3,3);
  for(Int_t i_phi = phi_start; i_phi < phi_stop; i_phi++)
  {
    c_Counts->cd(i_phi+1);
    c_Counts->cd(i_phi+1)->SetLeftMargin(0.15);
    c_Counts->cd(i_phi+1)->SetBottomMargin(0.15);
    c_Counts->cd(i_phi+1)->SetTicks(1,1);
    c_Counts->cd(i_phi+1)->SetGrid(0,0);
    TString KEY_Proton_QA = Form("pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_phi_Psi_%d_%s_Proton_SysError_%d",pt_QA,Cent_start,Charge_start,Eta_start,i_phi,Order[mOrder].Data(),Sys_start);
    h_mMass2[KEY_Proton_QA]->SetStats(0);
    h_mMass2[KEY_Proton_QA]->SetTitle("");
    h_mMass2[KEY_Proton_QA]->GetXaxis()->SetTitle("m^{2} ((GeV/c^{2})^{2})");
    h_mMass2[KEY_Proton_QA]->GetXaxis()->CenterTitle();
    h_mMass2[KEY_Proton_QA]->GetXaxis()->SetTitleSize(0.06);
    h_mMass2[KEY_Proton_QA]->GetYaxis()->SetTitle("Counts/Resolution");
    h_mMass2[KEY_Proton_QA]->GetYaxis()->CenterTitle();
    h_mMass2[KEY_Proton_QA]->GetYaxis()->SetTitleSize(0.06);
    h_mMass2[KEY_Proton_QA]->Draw("PE");
    PlotLine(m2_cut_low[pt_QA],m2_cut_low[pt_QA],0.0,h_mMass2[KEY_Proton_QA]->GetMaximum()/2.0,4,2,2);
    PlotLine( m2_cut_up[pt_QA], m2_cut_up[pt_QA],0.0,h_mMass2[KEY_Proton_QA]->GetMaximum()/2.0,4,2,2);
  }

  c_Counts->cd(8); // h_mCounts
  c_Counts->cd(8)->SetLeftMargin(0.15);
  c_Counts->cd(8)->SetBottomMargin(0.15);
  c_Counts->cd(8)->SetTicks(1,1);
  c_Counts->cd(8)->SetGrid(0,0);
  TString KEY_Proton_Counts_QA = Form("pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_%s_Proton_SysError_%d",pt_QA,Cent_start,Charge_start,Eta_start,Order[mOrder].Data(),Sys_start);
  h_mCounts[KEY_Proton_Counts_QA]->SetStats(0);
  h_mCounts[KEY_Proton_Counts_QA]->SetTitle("Counts/Resolution vs. #phi-#Psi_{3}");
  h_mCounts[KEY_Proton_Counts_QA]->GetXaxis()->SetTitle("#phi-#Psi_{3}");
  h_mCounts[KEY_Proton_Counts_QA]->GetXaxis()->CenterTitle();
  h_mCounts[KEY_Proton_Counts_QA]->GetXaxis()->SetTitleSize(0.06);
  h_mCounts[KEY_Proton_Counts_QA]->GetYaxis()->SetTitle("Counts/Resolution");
  h_mCounts[KEY_Proton_Counts_QA]->GetYaxis()->CenterTitle();
  h_mCounts[KEY_Proton_Counts_QA]->GetYaxis()->SetTitleSize(0.06);
  h_mCounts[KEY_Proton_Counts_QA]->SetMarkerStyle(24);
  h_mCounts[KEY_Proton_Counts_QA]->SetMarkerSize(0.5);
  h_mCounts[KEY_Proton_Counts_QA]->SetMarkerColor(kGray+3);
  h_mCounts[KEY_Proton_Counts_QA]->SetLineColor(1);
  h_mCounts[KEY_Proton_Counts_QA]->DrawCopy("PE");

  TF1 *f_flow_QA = new TF1("f_flow_QA",flow,0.0,PI_max[mOrder],3);
  f_flow_QA->FixParameter(0,ParFlow[KEY_Proton_Counts_QA][0]);
  f_flow_QA->FixParameter(1,ParFlow[KEY_Proton_Counts_QA][1]);
  f_flow_QA->FixParameter(2,ParFlow[KEY_Proton_Counts_QA][2]);
  f_flow_QA->SetLineColor(2);
  f_flow_QA->Draw("l same");

  TString pT_range_QA = Form("[%.2f,%.2f]",pt_low_raw[pt_QA],pt_up_raw[pt_QA]);
  plotTopLegend((char*)pT_range_QA.Data(),0.5,0.7,0.08,1,0.0,42,1);

  c_Counts->cd(9); // h_mCounts
  c_Counts->cd(9)->SetLeftMargin(0.15);
  c_Counts->cd(9)->SetBottomMargin(0.15);
  c_Counts->cd(9)->SetTicks(1,1);
  c_Counts->cd(9)->SetGrid(0,0);
  TString KEY_Proton_RawFlow_QA = Form("RawFlow_Centrality_%d_Charge_%d_EtaGap_%d_%s_Proton_SysError_%d",Cent_start,Charge_start,Eta_start,Order[mOrder].Data(),Sys_start);
  h_mRawFlow[KEY_Proton_RawFlow_QA]->SetStats(0);
  h_mRawFlow[KEY_Proton_RawFlow_QA]->SetTitle("raw v_{3} vs. p_{T}");
  h_mRawFlow[KEY_Proton_RawFlow_QA]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_mRawFlow[KEY_Proton_RawFlow_QA]->GetXaxis()->CenterTitle();
  h_mRawFlow[KEY_Proton_RawFlow_QA]->GetXaxis()->SetTitleSize(0.06);
  h_mRawFlow[KEY_Proton_RawFlow_QA]->GetYaxis()->SetTitle("raw v_{3}");
  h_mRawFlow[KEY_Proton_RawFlow_QA]->GetYaxis()->CenterTitle();
  h_mRawFlow[KEY_Proton_RawFlow_QA]->GetYaxis()->SetTitleSize(0.06);
  h_mRawFlow[KEY_Proton_RawFlow_QA]->SetMarkerStyle(24);
  h_mRawFlow[KEY_Proton_RawFlow_QA]->SetMarkerSize(0.5);
  h_mRawFlow[KEY_Proton_RawFlow_QA]->SetMarkerColor(kGray+3);
  h_mRawFlow[KEY_Proton_RawFlow_QA]->SetLineColor(1);
  h_mRawFlow[KEY_Proton_RawFlow_QA]->GetXaxis()->SetRangeUser(-0.05,3.4);
  h_mRawFlow[KEY_Proton_RawFlow_QA]->DrawCopy("PE");
  PlotLine(-0.05,3.4,0.0,0.0,1,2,2);
  */

  // resolution correction
  // read in yields of protons
  TString InPutFile_Yield = Form("./Data/AuAu%s/Proton/Yield_%s.root",Energy[mEnergy].Data(),Energy[mEnergy].Data());
  TFile *File_Yield = TFile::Open(InPutFile_Yield.Data());

  TH1FMap h_mMass2_Yields;
  vecFMap Yields_Counts;
  for(Int_t i_cent = 0; i_cent < 9; i_cent++) // centrality bin
  {
    for(Int_t i_charge = Charge_start; i_charge < Charge_stop; i_charge++) // charge bin
    {
      for(Int_t i_eta = Eta_start; i_eta < Eta_stop; i_eta++) // eta gap bin
      {
	for(Int_t i_sys = 0; i_sys < 18; i_sys++)
	{
	  TString KEY_Proton_Yield = Form("Centrality_%d_Charge_%d_EtaGap_%d_Yields_Proton_SysError_%d",i_cent,i_charge,i_eta,i_sys);
	  h_mMass2_Yields[KEY_Proton_Yield] = (TH1F*)File_Yield->Get(KEY_Proton_Yield.Data());
	  Int_t yields_start = h_mMass2_Yields[KEY_Proton_Yield]->FindBin(m2_yield_low);
	  Int_t yields_stop  = h_mMass2_Yields[KEY_Proton_Yield]->FindBin(m2_yield_high);
	  Float_t counts = 0.0;
	  Float_t errors = 0.0;
	  for(Int_t i_bin = yields_start; i_bin < yields_stop; i_bin++)
	  {
	    counts += h_mMass2_Yields[KEY_Proton_Yield]->GetBinContent(i_bin);
	    errors += h_mMass2_Yields[KEY_Proton_Yield]->GetBinError(i_bin)*h_mMass2_Yields[KEY_Proton_Yield]->GetBinError(i_bin);
	  }
	  Yields_Counts[KEY_Proton_Yield].clear();
	  Yields_Counts[KEY_Proton_Yield].push_back(static_cast<Float_t>(counts));
	  Yields_Counts[KEY_Proton_Yield].push_back(static_cast<Float_t>(TMath::Sqrt(errors)));
	}
      }
    }
  }

  /*
  // QA: h_mMass2_Yields vs Centrality
  TCanvas *c_Yields = new TCanvas("c_Yields","c_Yields",10,10,900,900);
  c_Yields->Divide(3,3);
  for(Int_t i_cent = 0; i_cent < 9; i_cent++)
  {
    c_Yields->cd(i_cent+1);
    c_Yields->cd(i_cent+1)->SetLeftMargin(0.15);
    c_Yields->cd(i_cent+1)->SetBottomMargin(0.15);
    c_Yields->cd(i_cent+1)->SetTicks(1,1);
    c_Yields->cd(i_cent+1)->SetGrid(0,0);
    TString KEY_Proton_Yield_QA = Form("Centrality_%d_Charge_%d_EtaGap_%d_Yields_Proton_SysError_%d",i_cent,Charge_start,Eta_start,Sys_start);
    h_mMass2_Yields[KEY_Proton_Yield_QA]->SetStats(0);
    h_mMass2_Yields[KEY_Proton_Yield_QA]->SetTitle("");
    h_mMass2_Yields[KEY_Proton_Yield_QA]->GetXaxis()->SetTitle("m^{2} ((GeV/c^{2})^{2})");
    h_mMass2_Yields[KEY_Proton_Yield_QA]->GetXaxis()->CenterTitle();
    h_mMass2_Yields[KEY_Proton_Yield_QA]->GetXaxis()->SetTitleSize(0.06);
    h_mMass2_Yields[KEY_Proton_Yield_QA]->GetYaxis()->SetTitle("Counts");
    h_mMass2_Yields[KEY_Proton_Yield_QA]->GetYaxis()->CenterTitle();
    h_mMass2_Yields[KEY_Proton_Yield_QA]->GetYaxis()->SetTitleSize(0.06);
    h_mMass2_Yields[KEY_Proton_Yield_QA]->SetMarkerStyle(24);
    h_mMass2_Yields[KEY_Proton_Yield_QA]->SetMarkerSize(0.5);
    h_mMass2_Yields[KEY_Proton_Yield_QA]->SetMarkerColor(kGray+3);
    h_mMass2_Yields[KEY_Proton_Yield_QA]->SetLineColor(1);
    h_mMass2_Yields[KEY_Proton_Yield_QA]->Draw("pE");

    KEY_Proton_Yield_QA = Form("Centrality_%d_Charge_%d_EtaGap_%d_Yields_Proton_SysError_%d",i_cent,Charge_start+1,Eta_start,Sys_start);
    h_mMass2_Yields[KEY_Proton_Yield_QA]->SetMarkerStyle(24);
    h_mMass2_Yields[KEY_Proton_Yield_QA]->SetMarkerSize(0.5);
    h_mMass2_Yields[KEY_Proton_Yield_QA]->SetMarkerColor(2);
    h_mMass2_Yields[KEY_Proton_Yield_QA]->SetLineColor(1);
    h_mMass2_Yields[KEY_Proton_Yield_QA]->Draw("pE same");

    PlotLine( m2_yield_low, m2_yield_low,0.0,h_mMass2_Yields[KEY_Proton_Yield_QA]->GetMaximum()/2.0,4,2,2);
    PlotLine(m2_yield_high,m2_yield_high,0.0,h_mMass2_Yields[KEY_Proton_Yield_QA]->GetMaximum()/2.0,4,2,2);
  }
  */

  // calculate final resolution correction factors and correct flow
  TString InPutFile_Res = Form("./Data/AuAu%s/file_%s_Resolution.root",Energy[mEnergy].Data(),Energy[mEnergy].Data());
  TFile *File_Res = TFile::Open(InPutFile_Res.Data());
  TString Res_Order[2] = {"Res2","Res3"};
  TProMap p_mRes;
  vecFMap ResValue;
  for(Int_t i_eta = Eta_start; i_eta < Eta_stop; i_eta++) // eta_gap loop
  {
    TString KEY_eta = Form("%s_EtaGap_%d_EP",Res_Order[mOrder].Data(),i_eta);
    p_mRes[KEY_eta] = (TProfile*)File_Res->Get(KEY_eta.Data()); // read in resolution file
    for(Int_t i_charge = Charge_start; i_charge < Charge_stop; i_charge++) // charge bin
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
	      TString KEY_Proton_Yield = Form("Centrality_%d_Charge_%d_EtaGap_%d_Yields_Proton_SysError_%d",cent,i_charge,i_eta,i_sys);
	      ResValue[KEY_Proton_Yield].push_back(static_cast<Float_t>(TMath::Sqrt(p_mRes[KEY_eta]->GetBinContent(cent+1))));
	      yields_total += Yields_Counts[KEY_Proton_Yield][0];
	    }
	  }

	  Float_t mean_res = 0.0;
	  for(Int_t cent = cent_low[i_cent]; cent <= cent_up[i_cent]; cent++) // calculate final resolution correction factor <1/R(centrality)>
	  {
	    TString KEY_Proton_Yield = Form("Centrality_%d_Charge_%d_EtaGap_%d_Yields_Proton_SysError_%d",cent,i_charge,i_eta,i_sys);
	    mean_res += Yields_Counts[KEY_Proton_Yield][0]/(ResValue[KEY_Proton_Yield][0]*yields_total);
	  }
//	  cout << "i_eta = " << i_eta << ", i_charge = " << i_charge << ", i_sys = " << i_sys << ", centrality_bin = " << i_cent << ", mean_res = " << mean_res << endl;

	  TString KEY_Proton_RawFlow = Form("RawFlow_Centrality_%d_Charge_%d_EtaGap_%d_%s_Proton_SysError_%d",i_cent,i_charge,i_eta,Order[mOrder].Data(),i_sys);
	  h_mRawFlow[KEY_Proton_RawFlow]->Scale(mean_res);
	}
      }
    }
  }

  // read in pT spectra
  TString InPutFile_Pt = Form("./OutPut/AuAu%s/Proton/h_pt_%s.root",Energy[mEnergy].Data(),Order[mOrder].Data());
  TFile *File_Spec = TFile::Open(InPutFile_Pt.Data());
  TH1FMap h_mPt;
  for(Int_t i_cent = Cent_start; i_cent < Cent_stop; i_cent++) // centrality bin
  {
    for(Int_t i_charge = Charge_start; i_charge < Charge_stop; i_charge++) // charge bin
    {
      for(Int_t i_eta = Eta_start; i_eta < Eta_stop; i_eta++) // eta gap bin
      {
	for(Int_t i_sys = 0; i_sys < 18; i_sys++)
	{
	  TString KEY_Proton_pT = Form("Spec_Centrality_%d_Charge_%d_EtaGap_%d_Proton_%s_SysErrors_%d",i_cent,i_charge,i_eta,Order[mOrder].Data(),i_sys);
	  h_mPt[KEY_Proton_pT] = (TH1F*)File_Spec->Get(KEY_Proton_pT.Data());
	}
      }
    }
  }

  // mean pT calculations
  Float_t pt_low_spec[32], pt_up_spec[32], pt_width[32], pt_center[32];
  for(Int_t i_pt = pt_start; i_pt < pt_stop; i_pt++)
  {
    pt_low_spec[2*i_pt] = pt_low_raw[i_pt];
    pt_up_spec[2*i_pt]  = 0.5*(pt_low_raw[i_pt]+pt_up_raw[i_pt]);
    pt_width[2*i_pt]    = pt_up_spec[2*i_pt]-pt_low_spec[2*i_pt];
    pt_center[2*i_pt]   = 0.5*(pt_up_spec[2*i_pt]+pt_low_spec[2*i_pt]);

    pt_low_spec[2*i_pt+1] = 0.5*(pt_low_raw[i_pt]+pt_up_raw[i_pt]);
    pt_up_spec[2*i_pt+1]  = pt_up_raw[i_pt];
    pt_width[2*i_pt+1]    = pt_up_spec[2*i_pt+1]-pt_low_spec[2*i_pt+1];
    pt_center[2*i_pt+1]   = 0.5*(pt_up_spec[2*i_pt+1]+pt_low_spec[2*i_pt+1]);
  }

  vecFMap mean_pt; // mean pt with systematic errors
  for(Int_t i_cent = Cent_start; i_cent < Cent_stop; i_cent++) // Centrality loop
  {
    for(Int_t i_charge = Charge_start; i_charge < Charge_stop; i_charge++) // charge bin
    {
      for(Int_t i_eta = Eta_start; i_eta < Eta_stop; i_eta++) // EtaGap loop
      {
	for(Int_t i_sys = Sys_start; i_sys < Sys_stop; i_sys++) // Systematic loop
	{
	  TString KEY_Proton_RawFlow = Form("RawFlow_Centrality_%d_Charge_%d_EtaGap_%d_%s_Proton_SysError_%d",i_cent,i_charge,i_eta,Order[mOrder].Data(),i_sys);
	  mean_pt[KEY_Proton_RawFlow].clear();
	  for(Int_t i_pt = pt_rebin_first; i_pt < pt_rebin_last; i_pt++) // loop over rebinned pT
	  {
//	    cout << "i_pt = " << i_pt << endl;
	    Float_t mean_pt_counts = 0.0;
	    Float_t spec_counts    = 0.0;
	    for(Int_t i_pt_raw = pt_rebin_start[i_pt]; i_pt_raw <= pt_rebin_stop[i_pt]; i_pt_raw++) // loop over raw pT bin
	    {
//	      cout << "i_pt_raw = " << i_pt_raw << endl;
	      for(Int_t i_pt_spec = 2*i_pt_raw; i_pt_spec <= 2*i_pt_raw+1; i_pt_spec++)
	      {
//		cout << "i_pt_spec = " << i_pt_spec << endl;
		TString KEY_Proton_pT = Form("Spec_Centrality_%d_Charge_%d_EtaGap_%d_Proton_%s_SysErrors_%d",i_cent,i_charge,i_eta,Order[mOrder].Data(),i_sys);
		mean_pt_counts += pt_center[i_pt_spec]*h_mPt[KEY_Proton_pT]->GetBinContent(i_pt_spec+1)*pt_width[i_pt_spec];
		spec_counts    += h_mPt[KEY_Proton_pT]->GetBinContent(i_pt_spec+1)*pt_width[i_pt_spec];

	      }
	    }
	    mean_pt[KEY_Proton_RawFlow].push_back(static_cast<Float_t>(mean_pt_counts/spec_counts));
	  }
	}
      }
    }
  }

  // set final pt and flow to one TGraphAsymmErrors
  TGraMap g_mFlow;
  for(Int_t i_cent = Cent_start; i_cent < Cent_stop; i_cent++) // Centrality loop
  {
    for(Int_t i_charge = Charge_start; i_charge < Charge_stop; i_charge++) // charge bin
    {
      for(Int_t i_eta = Eta_start; i_eta < Eta_stop; i_eta++) // EtaGap loop
      {
	for(Int_t i_sys = Sys_start; i_sys < Sys_stop; i_sys++) // Systematic loop
	{
	  TString KEY_Proton_RawFlow = Form("RawFlow_Centrality_%d_Charge_%d_EtaGap_%d_%s_Proton_SysError_%d",i_cent,i_charge,i_eta,Order[mOrder].Data(),i_sys);
	  g_mFlow[KEY_Proton_RawFlow] = new TGraphAsymmErrors();
	  for(Int_t i_pt = pt_rebin_first; i_pt < pt_rebin_last; i_pt++) // loop over rebinned pT
	  {
	    Float_t pt_mean = (pt_low[i_pt]+pt_up[i_pt])/2.0;

	    Float_t content = h_mRawFlow[KEY_Proton_RawFlow]->GetBinContent(h_mRawFlow[KEY_Proton_RawFlow]->FindBin(pt_mean)); // bin counting
	    Float_t error   = h_mRawFlow[KEY_Proton_RawFlow]->GetBinError(h_mRawFlow[KEY_Proton_RawFlow]->FindBin(pt_mean));
	    g_mFlow[KEY_Proton_RawFlow]->SetPoint(i_pt,mean_pt[KEY_Proton_RawFlow][i_pt],content);
	    g_mFlow[KEY_Proton_RawFlow]->SetPointError(i_pt,0.0,0.0,error,error);
	    TString Name = Form("Flow_%s_Proton_Centrality_%d_Charge_%d_EtaGap_%d_SysError_%d",Order[mOrder].Data(),i_cent,i_charge,i_eta,i_sys);
	    g_mFlow[KEY_Proton_RawFlow]->SetName(Name.Data());
	  }
	}
      }
    }
  }

  // QA flow vs. pt for gaussian fits and breit wigner fits with TGraphAsymmErrors
  TCanvas *c_flow = new TCanvas("c_flow","c_flow",10,10,800,800);
  c_flow->SetLeftMargin(0.15);
  c_flow->SetBottomMargin(0.15);
  c_flow->SetTicks(1,1);
  c_flow->SetGrid(0,0);
  TH1F *h_play = new TH1F("h_play","h_play",100,0.0,4.0);
  for(Int_t i_bin = 0; i_bin < 100; i_bin++)
  {
    h_play->SetBinContent(i_bin+1,-10.0);
    h_play->SetBinError(i_bin+1,1.0);
  }
  h_play->SetTitle("");
  h_play->SetStats(0);
  h_play->GetXaxis()->SetRangeUser(0.0,4.0);
  h_play->GetXaxis()->SetNdivisions(505,'N');
  h_play->GetXaxis()->SetLabelSize(0.03);
  h_play->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  h_play->GetXaxis()->SetTitleSize(0.05);
  h_play->GetXaxis()->SetTitleOffset(1.2);
  h_play->GetXaxis()->CenterTitle();

  h_play->GetYaxis()->SetRangeUser(-0.01,0.20);
  h_play->GetYaxis()->SetNdivisions(505,'N');
  h_play->GetYaxis()->SetTitle("v_{3}");
  h_play->GetYaxis()->SetTitleSize(0.05);
  h_play->GetYaxis()->SetLabelSize(0.03);
  h_play->GetYaxis()->CenterTitle();
  h_play->DrawCopy("pE");
  PlotLine(0.0,4.0,0.0,0.0,1,2,2);
  for(Int_t i_sys = Sys_start; i_sys < Sys_stop; i_sys++) // Systematic loop
  {
    TString KEY_Proton_RawFlow = Form("RawFlow_Centrality_%d_Charge_%d_EtaGap_%d_%s_Proton_SysError_%d",Cent_start,Charge_start,Eta_start,Order[mOrder].Data(),i_sys);
    g_mFlow[KEY_Proton_RawFlow]->SetMarkerStyle(24);
    g_mFlow[KEY_Proton_RawFlow]->SetMarkerColor(4);
    g_mFlow[KEY_Proton_RawFlow]->Draw("PE same");
  }

  TString OutPutFile = Form("./OutPut/AuAu%s/Proton/flow_%s.root",Energy[mEnergy].Data(),Order[mOrder].Data());
  TFile *File_OutPut = new TFile(OutPutFile.Data(),"RECREATE");
  File_OutPut->cd();
  for(Int_t i_cent = Cent_start; i_cent < Cent_stop; i_cent++) // Centrality loop
  {
    for(Int_t i_charge = Charge_start; i_charge < Charge_stop; i_charge++) // charge bin
    {
      for(Int_t i_eta = Eta_start; i_eta < Eta_stop; i_eta++) // EtaGap loop
      {
	for(Int_t i_sys = Sys_start; i_sys < Sys_stop; i_sys++) // Systematic loop
	{
	  TString KEY_Proton_RawFlow = Form("RawFlow_Centrality_%d_Charge_%d_EtaGap_%d_%s_Proton_SysError_%d",i_cent,i_charge,i_eta,Order[mOrder].Data(),i_sys);
	  g_mFlow[KEY_Proton_RawFlow]->SetMarkerStyle(24);
	  g_mFlow[KEY_Proton_RawFlow]->SetMarkerColor(4);
	  g_mFlow[KEY_Proton_RawFlow]->Write();
	}
      }
    }
  }
  h_play->Write();
  File_OutPut->Close();
}

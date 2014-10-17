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

#include "/u/xusun/STAR/Analysis_v3/StRoot/CombPiKP/student_t_1d_single.h"
#include "/u/xusun/STAR/Analysis_v3/StRoot/CombPiKP/student_t_1d_double.h"
#include "/u/xusun/STAR/Analysis_v3/StRoot/CombPiKP/draw.h"

// mEnergy: 0 for 200GeV, 1 for 39GeV
// mCharge: 0 for positive, 1 for negative
// mCentrality: 0 for 00-80, 1 for 00-10, 2 for 10-40, 3 for 40-80
void CombPID_1d(Int_t mEnergy = 0, Int_t mCharge = 0, Int_t mCentrality = 0)
{
  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(50000);

  const Int_t pt_total = 16;
  const Int_t cent_total = 4; // 4
  const Int_t cent_start[4] = {0,1,2,3};
  const Int_t cent_stop[4]  = {1,2,3,4};
  const Int_t charge_total = 2; // 2
  const Int_t eta_total = 4; // 4
  const Int_t eta_start = 0;
  const Int_t eta_stop  = 1;
  const Int_t phi_psi_total = 7; // 7
  const Int_t flow_total = 2; // 2
  TString Energy[2] = {"200GeV","39GeV"};
  TString Centrality1[4] = {"0-80%","0-10%","10-40%","40-80%"};
  TString Centrality2[4] = {"0-70%","0-10%","10-40%","40-70%"};
  TString Centrality[4] = {"0080","0010","1040","4080"};
  TString Charge[2] = {"pos","neg"};
  TString EtaGap[4] = {"#eta_{gap} = 0.05","#eta_{gap} = 0.10","#eta_{gap} = 0.20","#eta_{gap} = 0.50"};
  TString Order[2] = {"2nd","3rd"};
  Double_t parfit[pt_total][cent_total][charge_total][eta_total][flow_total][8];
  Double_t PI_max[2] = {TMath::Pi()/2.0,TMath::Pi()/3.0};

  /*
  // Initial parameters 200 GeV
  Float_t order_pion[pt_total]   = {3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,2.50,2.00,3.00,3.0};
  Float_t order_kaon[pt_total]   = {3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,2.55,2.50,3.00,3.0};
  Float_t order_proton[pt_total] = {3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,2.55,3.00,2.55,3.0};
  */

  /*
  // Initial parameters 39 GeV for 00-80
  Float_t order_pion[pt_total]   = {3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,2.50,3.0,3.00,2.5,3.00,3.0};
  Float_t order_kaon[pt_total]   = {3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,2.55,2.5,2.50,2.0,2.50,3.0};
  Float_t order_proton[pt_total] = {3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.00,3.0,2.50,2.5,2.55,3.0};
  */

  // Initial parameters 39 GeV for 10-40
  Float_t order_pion[pt_total]   = {3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,2.55,2.55,3.00,2.5,3.00,3.0};
  Float_t order_kaon[pt_total]   = {3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.00,3.00,3.00,2.0,2.00,3.0};
  Float_t order_proton[pt_total] = {3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,2.55,2.55,3.00,2.5,3.00,3.0};

  // pion counts
  Double_t sigma_pion[pt_total][cent_total][charge_total][eta_total][flow_total];
  Int_t bin_pion_start[pt_total][cent_total][charge_total][eta_total][phi_psi_total][flow_total];
  Int_t bin_pion_stop[pt_total][cent_total][charge_total][eta_total][phi_psi_total][flow_total];
  Double_t counts_pion[pt_total][cent_total][charge_total][eta_total][phi_psi_total][flow_total];
  Double_t errors_pion[pt_total][cent_total][charge_total][eta_total][phi_psi_total][flow_total];
  // kaon counts
  Double_t sigma_kaon[pt_total][cent_total][charge_total][eta_total][flow_total];
  Int_t bin_kaon_start[pt_total][cent_total][charge_total][eta_total][phi_psi_total][flow_total];
  Int_t bin_kaon_stop[pt_total][cent_total][charge_total][eta_total][phi_psi_total][flow_total];
  Double_t counts_kaon[pt_total][cent_total][charge_total][eta_total][phi_psi_total][flow_total];
  Double_t errors_kaon[pt_total][cent_total][charge_total][eta_total][phi_psi_total][flow_total];

  //-----------------------------------------------------------------------------------------
  // Initial parameters
  // pion
  Float_t y_pi[pt_total] = {0.0005,0.00125,0.00125,-0.00025,0.002,0.00825,0.01075,0.02175,0.02825,0.03137,0.02625,0.03375,0.0114,0.0236,0.0236,0.0825};

  // proton
  Float_t y_p[pt_total] = {0.0375,0.02075,-0.01125,-0.05975,-0.126,-0.2033,-0.2862,-0.3553,-0.4268,-0.4614,-0.5062,-0.5587,-0.5015,-0.5071,-0.4981,-0.5295};

  // pt bin
  Float_t pt_low[pt_total] = {0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4};
  Float_t pt_up[pt_total]  = {0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,4.2};
  TString Title[pt_total];
  for(Int_t i = 0; i < pt_total; i++)
  {
    Title[i] = Form("%1.1f-%1.1f GeV/c",pt_low[i],pt_up[i]);
  }

  // x and y range
  Float_t x_low[pt_total] = {-0.4, -0.6,-0.6, -0.6,-0.6,-0.6,-0.6,-0.8,-0.8,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0};
  Float_t x_up[pt_total]  = { 1.4,  1.4, 1.4,  1.4, 1.4, 1.4, 1.4, 1.6, 1.6, 2.0, 2.0, 2.0, 2.4, 2.4, 2.4, 2.4};

  Float_t y_low[pt_total] = {-0.15,-0.2,-0.3,-0.45,-0.6,-0.8,-1.0,-1.3,-1.5,-1.7,-1.8,-1.8,-2.0,-2.0,-2.0,-2.0};
  Float_t y_up[pt_total]  = { 0.15, 0.2, 0.2, 0.20, 0.3, 0.4, 0.6, 1.0, 1.0, 1.0, 1.2, 1.2, 1.4, 1.4, 1.4, 1.4};
  //-----------------------------------------------------------------------------------------
  TString inputname_hist = Form("/project/projectdirs/star/xusun/OutPut/AuAu%s/Mass2_nSigmaPion/merged_file/flow_pik/Mass2_nSigmaPion_cut_%s_charge_%d_etagap_%d.root",Energy[mEnergy].Data(),Centrality[mCentrality].Data(),mCharge,eta_start);
  TFile *File_hist = TFile::Open(inputname_hist.Data());
  TH2F *h_XY_cut[pt_total][cent_total][charge_total][eta_total][phi_psi_total][flow_total];
  for(Int_t pt_bin = 0; pt_bin < pt_total; pt_bin ++)
  {
    for(Int_t cent = cent_start[mCentrality]; cent < cent_stop[mCentrality]; cent++)
    {
      for(Int_t charge = mCharge; charge < mCharge+1; charge++)
      {
	for(Int_t eta_bin = eta_start; eta_bin < eta_stop; eta_bin++)
	{
	  for(Int_t phi_psi_bin = 0; phi_psi_bin < phi_psi_total; phi_psi_bin++)
	  {
	    for(Int_t flow_bin = 0; flow_bin < flow_total; flow_bin++)
	    {
	      TString HistName = Form("pt_%d_Centrality_%d_Charge_%d_etagap_%d_phi_psi_%d_%s_cut",pt_bin,cent,charge,eta_bin,phi_psi_bin,Order[flow_bin].Data());
	      h_XY_cut[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin] = (TH2F*)File_hist->Get(HistName.Data());
	    }
	  }
	}
      }
    }
  }
  TH2F *h_XY_total[pt_total][cent_total][charge_total][eta_total][flow_total];
  for(Int_t pt_bin = 0; pt_bin < pt_total; pt_bin ++)
  {
    for(Int_t cent = cent_start[mCentrality]; cent < cent_stop[mCentrality]; cent++)
    {
      for(Int_t charge = mCharge; charge < mCharge+1; charge++)
      {
	for(Int_t eta_bin = eta_start; eta_bin < eta_stop; eta_bin++)
	{
	  for(Int_t flow_bin = 0; flow_bin < flow_total; flow_bin++)
	  {
	    for(Int_t phi_psi_bin = 0; phi_psi_bin < phi_psi_total; phi_psi_bin++)
	    {
	      if(phi_psi_bin == 0)
	      {
		TString HistName = Form("pt_%d_Centrality_%d_Charge_%d_etagap_%d_phi_psi_%d_%s_total",pt_bin,cent,charge,eta_bin,phi_psi_bin,Order[flow_bin].Data());
		h_XY_total[pt_bin][cent][charge][eta_bin][flow_bin] = (TH2F*)h_XY_cut[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->Clone(HistName.Data());
	      }
	      if(phi_psi_bin > 0)
	      {
		h_XY_total[pt_bin][cent][charge][eta_bin][flow_bin]->Add(h_XY_cut[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin],1.0);
	      }
	    }
	  }
	}
      }
    }
  }

  TString inputname_par = Form("/project/projectdirs/star/xusun/OutPut/AuAu%s/Mass2_nSigmaPion/merged_file/flow_pik/parfit_%s_charge_%d_etagap_%d.root",Energy[mEnergy].Data(),Centrality[mCentrality].Data(),mCharge,eta_start);
  TFile *File_par = TFile::Open(inputname_par.Data());
  Double_t parfit_1st[pt_total][cent_total][charge_total][eta_total][flow_total][19];
  Double_t parfit_2nd[pt_total][cent_total][charge_total][eta_total][flow_total][19];
  TH1F *h_parfit_1st[pt_total][cent_total][charge_total][eta_total][flow_total];
  TH1F *h_parfit_2nd[pt_total][cent_total][charge_total][eta_total][flow_total];
  for(Int_t pt_bin = 0; pt_bin < pt_total; pt_bin ++)
  {
    for(Int_t cent = cent_start[mCentrality]; cent < cent_stop[mCentrality]; cent++)
    {
      for(Int_t charge = mCharge; charge < mCharge+1; charge++)
      {
	for(Int_t eta_bin = eta_start; eta_bin < eta_stop; eta_bin++)
	{
	  for(Int_t flow_bin = 0; flow_bin < flow_total; flow_bin++)
	  {
	    TString ParName_1st = Form("parfit_pt_%d_Centrality_%d_charge_%d_etagap_%d_%s_total_1st",pt_bin,cent,charge,eta_bin,Order[flow_bin].Data());
	    h_parfit_1st[pt_bin][cent][charge][eta_bin][flow_bin] = (TH1F*)File_par->Get(ParName_1st.Data());
	    TString ParName_2nd = Form("parfit_pt_%d_Centrality_%d_charge_%d_etagap_%d_%s_total_2nd",pt_bin,cent,charge,eta_bin,Order[flow_bin].Data());
	    h_parfit_2nd[pt_bin][cent][charge][eta_bin][flow_bin] = (TH1F*)File_par->Get(ParName_2nd.Data());
	    for(Int_t i = 0; i < 19; i++)
	    {
	      parfit_1st[pt_bin][cent][charge][eta_bin][flow_bin][i] = 0.0;
	      parfit_2nd[pt_bin][cent][charge][eta_bin][flow_bin][i] = 0.0;
	      parfit_1st[pt_bin][cent][charge][eta_bin][flow_bin][i] = h_parfit_1st[pt_bin][cent][charge][eta_bin][flow_bin]->GetBinContent(i+1);
	      parfit_2nd[pt_bin][cent][charge][eta_bin][flow_bin][i] = h_parfit_2nd[pt_bin][cent][charge][eta_bin][flow_bin]->GetBinContent(i+1);
	    }
	  }
	}
      }
    }
  }

  //----------------------------------------------------------------------------------------------------------
  // fit for integral phi-psi bin
  // ProjectionX
  TH1F *h_X_total[pt_total][cent_total][charge_total][eta_total][flow_total];
  for(Int_t pt_bin = 0; pt_bin < pt_total; pt_bin++)
  {
    for(Int_t cent = cent_start[mCentrality]; cent < cent_stop[mCentrality]; cent++)
    {
      for(Int_t charge = mCharge; charge < mCharge+1; charge++)
      {
	for(Int_t eta_bin = eta_start; eta_bin < eta_stop; eta_bin++)
	{
	  for(Int_t flow_bin = 0; flow_bin < flow_total; flow_bin++)
	  {
	    if(pt_bin < 5)
	    {
	      TString HistName = Form("pt_%d_Centrality_%d_charge_%d_etagap_%d_%s_X_total",pt_bin,cent,charge,eta_bin,Order[flow_bin].Data());

	      h_X_total[pt_bin][cent][charge][eta_bin][flow_bin] = (TH1F*)h_XY_total[pt_bin][cent][charge][eta_bin][flow_bin]->ProjectionX(HistName.Data(),h_XY_total[pt_bin][cent][charge][eta_bin][flow_bin]->GetYaxis()->FindBin(y_low[pt_bin]),h_XY_total[pt_bin][cent][charge][eta_bin][flow_bin]->GetYaxis()->FindBin(y_up[pt_bin]));
	    }
	    if(pt_bin >= 5)
	    {
	      TString HistName = Form("pt_%d_Centrality_%d_charge_%d_etagap_%d_%s_X_total",pt_bin,cent,charge,eta_bin,Order[flow_bin].Data());
	      h_X_total[pt_bin][cent][charge][eta_bin][flow_bin] = (TH1F*)h_XY_total[pt_bin][cent][charge][eta_bin][flow_bin]->ProjectionX(HistName.Data(),h_XY_total[pt_bin][cent][charge][eta_bin][flow_bin]->GetYaxis()->FindBin(y_p[pt_bin]-3.0*parfit_1st[pt_bin][cent][charge][eta_bin][flow_bin][17]),h_XY_total[pt_bin][cent][charge][eta_bin][flow_bin]->GetYaxis()->FindBin(y_pi[pt_bin]+3.0*parfit_1st[pt_bin][cent][charge][eta_bin][flow_bin][5]));
	    }
	    h_X_total[pt_bin][cent][charge][eta_bin][flow_bin]->SetLineStyle(2);
	    h_X_total[pt_bin][cent][charge][eta_bin][flow_bin]->SetLineColor(1);
	    h_X_total[pt_bin][cent][charge][eta_bin][flow_bin]->SetLineWidth(1);

	    h_X_total[pt_bin][cent][charge][eta_bin][flow_bin]->SetStats(0);
	    h_X_total[pt_bin][cent][charge][eta_bin][flow_bin]->SetTitle(Title[pt_bin].Data());

	    h_X_total[pt_bin][cent][charge][eta_bin][flow_bin]->GetXaxis()->SetNdivisions(505,'N');
	    h_X_total[pt_bin][cent][charge][eta_bin][flow_bin]->GetXaxis()->SetTitle("x (n#sigma_{#pi},m^{2})");
	    h_X_total[pt_bin][cent][charge][eta_bin][flow_bin]->GetXaxis()->SetTitleSize(0.05);
	    h_X_total[pt_bin][cent][charge][eta_bin][flow_bin]->GetXaxis()->CenterTitle();
	    h_X_total[pt_bin][cent][charge][eta_bin][flow_bin]->GetXaxis()->SetTitleOffset(0.9);

	    h_X_total[pt_bin][cent][charge][eta_bin][flow_bin]->GetYaxis()->SetNdivisions(505,'N');
	    h_X_total[pt_bin][cent][charge][eta_bin][flow_bin]->GetYaxis()->SetTitle("Counts");
	    h_X_total[pt_bin][cent][charge][eta_bin][flow_bin]->GetYaxis()->SetTitleSize(0.05);
	    h_X_total[pt_bin][cent][charge][eta_bin][flow_bin]->GetYaxis()->CenterTitle();
	    h_X_total[pt_bin][cent][charge][eta_bin][flow_bin]->GetYaxis()->SetTitleOffset(1.5);
	  }
	}
      }
    }
  }

  // fit integral distribution
  TF1 *f_1d_double_total[pt_total][cent_total][charge_total][eta_total][flow_total];
  TF1 *f_1d_pion_total[pt_total][cent_total][charge_total][eta_total][flow_total];
  TF1 *f_1d_kaon_total[pt_total][cent_total][charge_total][eta_total][flow_total];
  for(Int_t pt_bin = 5; pt_bin < pt_total; pt_bin ++)
  {
    for(Int_t cent = cent_start[mCentrality]; cent < cent_stop[mCentrality]; cent++)
    {
      for(Int_t charge = mCharge; charge < mCharge+1; charge++)
      {
	for(Int_t eta_bin = eta_start; eta_bin < eta_stop; eta_bin++)
	{
	  for(Int_t flow_bin = 0; flow_bin < flow_total; flow_bin++)
	  {
	    TString FuncName = Form("f_1d_double_total_pt_%d_cent_%d_charge_%d_etagap_%d_%s",pt_bin,cent,charge,eta_bin,Order[flow_bin].Data());
	    f_1d_double_total[pt_bin][cent][charge][eta_bin][flow_bin] = new TF1(FuncName.Data(),student_t_1d_double,x_low[pt_bin],x_up[pt_bin],8);
	    f_1d_double_total[pt_bin][cent][charge][eta_bin][flow_bin]->SetNpx(400);
	    /*
	    if(pt_bin > 8)
	    {
	      f_1d_double_total[pt_bin][cent][charge][eta_bin][flow_bin]->SetNpx(200);
	      if(pt_bin == 12)
	      {
		f_1d_double_total[pt_bin][cent][charge][eta_bin][flow_bin]->SetNpx(100);
	      }
	    }
	    */

	    // pion
	    f_1d_double_total[pt_bin][cent][charge][eta_bin][flow_bin]->SetParameter(0,parfit_2nd[pt_bin][cent][charge][eta_bin][flow_bin][0]);
	    f_1d_double_total[pt_bin][cent][charge][eta_bin][flow_bin]->SetParameter(1,parfit_2nd[pt_bin][cent][charge][eta_bin][flow_bin][2]);
	    f_1d_double_total[pt_bin][cent][charge][eta_bin][flow_bin]->SetParameter(2,parfit_2nd[pt_bin][cent][charge][eta_bin][flow_bin][4]);
//	    f_1d_double_total[pt_bin][cent][charge][eta_bin][flow_bin]->SetParameter(3,parfit_2nd[pt_bin][cent][charge][eta_bin][flow_bin][6]);
	    f_1d_double_total[pt_bin][cent][charge][eta_bin][flow_bin]->SetParameter(3,1.0);

	    // kaon
	    f_1d_double_total[pt_bin][cent][charge][eta_bin][flow_bin]->SetParameter(4,(2.0*parfit_2nd[pt_bin][cent][charge][eta_bin][flow_bin][0]+parfit_2nd[pt_bin][cent][charge][eta_bin][flow_bin][12])/3.0);
	    f_1d_double_total[pt_bin][cent][charge][eta_bin][flow_bin]->SetParameter(5,parfit_2nd[pt_bin][cent][charge][eta_bin][flow_bin][7]);
	    f_1d_double_total[pt_bin][cent][charge][eta_bin][flow_bin]->SetParameter(6,parfit_2nd[pt_bin][cent][charge][eta_bin][flow_bin][9]);
//	    f_1d_double_total[pt_bin][cent][charge][eta_bin][flow_bin]->SetParameter(7,parfit_2nd[pt_bin][cent][charge][eta_bin][flow_bin][11]);
	    f_1d_double_total[pt_bin][cent][charge][eta_bin][flow_bin]->SetParameter(7,1.0);

	    Float_t pi_max = h_X_total[pt_bin][cent][charge][eta_bin][flow_bin]->GetBinContent(h_X_total[pt_bin][cent][charge][eta_bin][flow_bin]->FindBin(parfit_2nd[pt_bin][cent][charge][eta_bin][flow_bin][2]));
	    Float_t pi_func = f_1d_double_total[pt_bin][cent][charge][eta_bin][flow_bin]->Eval(parfit_2nd[pt_bin][cent][charge][eta_bin][flow_bin][2]);
	    Float_t Norm_pi = pi_max/pi_func;

	    Float_t k_max = h_X_total[pt_bin][cent][charge][eta_bin][flow_bin]->GetBinContent(h_X_total[pt_bin][cent][charge][eta_bin][flow_bin]->FindBin(parfit_2nd[pt_bin][cent][charge][eta_bin][flow_bin][7]));
	    Float_t k_func = f_1d_double_total[pt_bin][cent][charge][eta_bin][flow_bin]->Eval(parfit_2nd[pt_bin][cent][charge][eta_bin][flow_bin][7]);
	    Float_t Norm_k = k_max/k_func;

	    f_1d_double_total[pt_bin][cent][charge][eta_bin][flow_bin]->SetParameter(3,Norm_pi);
	    f_1d_double_total[pt_bin][cent][charge][eta_bin][flow_bin]->SetParameter(7,Norm_k);

	    f_1d_double_total[pt_bin][cent][charge][eta_bin][flow_bin]->SetLineStyle(1);
	    f_1d_double_total[pt_bin][cent][charge][eta_bin][flow_bin]->SetLineColor(2);
	    f_1d_double_total[pt_bin][cent][charge][eta_bin][flow_bin]->SetFillStyle(0);

	    f_1d_double_total[pt_bin][cent][charge][eta_bin][flow_bin]->SetRange(parfit_2nd[pt_bin][cent][charge][eta_bin][flow_bin][2]-order_pion[pt_bin]*parfit_2nd[pt_bin][cent][charge][eta_bin][flow_bin][4],parfit_2nd[pt_bin][cent][charge][eta_bin][flow_bin][7]+order_kaon[pt_bin]*parfit_2nd[pt_bin][cent][charge][eta_bin][flow_bin][9]);

	    // Fit
	    cout << "pt_bin = " << pt_bin << ", charge = " << charge << ", eta_bin = " << eta_bin <<  ", " << Order[flow_bin].Data() << endl;
	    h_X_total[pt_bin][cent][charge][eta_bin][flow_bin]->Fit(f_1d_double_total[pt_bin][cent][charge][eta_bin][flow_bin],"MRN");

	    Float_t chi2 = f_1d_double_total[pt_bin][cent][charge][eta_bin][flow_bin]->GetChisquare();
	    Float_t Ndf = f_1d_double_total[pt_bin][cent][charge][eta_bin][flow_bin]->GetNDF();
	    cout << "chi2/Ndf = " << chi2 << "/" << Ndf << " = " << chi2/Ndf << endl;
	    for(Int_t i = 0; i < 8; i++)
	    {
	      parfit[pt_bin][cent][charge][eta_bin][flow_bin][i] = f_1d_double_total[pt_bin][cent][charge][eta_bin][flow_bin]->GetParameter(i);
	    }

	    f_1d_double_total[pt_bin][cent][charge][eta_bin][flow_bin]->SetRange(x_low[pt_bin],x_up[pt_bin]);

	    TString PionName = Form("f_1d_pion_total_pt_%d_cent_%d_charge_%d_etagap_%d__%s",pt_bin,cent,charge,eta_bin,Order[flow_bin].Data());
	    f_1d_pion_total[pt_bin][cent][charge][eta_bin][flow_bin] = new TF1(PionName.Data(),student_t_1d_single,x_low[pt_bin],x_up[pt_bin],4);
	    f_1d_pion_total[pt_bin][cent][charge][eta_bin][flow_bin]->SetNpx(400);
	    f_1d_pion_total[pt_bin][cent][charge][eta_bin][flow_bin]->FixParameter(0,f_1d_double_total[pt_bin][cent][charge][eta_bin][flow_bin]->GetParameter(0));
	    f_1d_pion_total[pt_bin][cent][charge][eta_bin][flow_bin]->FixParameter(1,f_1d_double_total[pt_bin][cent][charge][eta_bin][flow_bin]->GetParameter(1));
	    f_1d_pion_total[pt_bin][cent][charge][eta_bin][flow_bin]->FixParameter(2,f_1d_double_total[pt_bin][cent][charge][eta_bin][flow_bin]->GetParameter(2));
	    f_1d_pion_total[pt_bin][cent][charge][eta_bin][flow_bin]->FixParameter(3,f_1d_double_total[pt_bin][cent][charge][eta_bin][flow_bin]->GetParameter(3));
	    f_1d_pion_total[pt_bin][cent][charge][eta_bin][flow_bin]->SetFillColor(kAzure-2);
	    f_1d_pion_total[pt_bin][cent][charge][eta_bin][flow_bin]->SetLineColor(kAzure-2);
	    f_1d_pion_total[pt_bin][cent][charge][eta_bin][flow_bin]->SetLineStyle(1);
	    f_1d_pion_total[pt_bin][cent][charge][eta_bin][flow_bin]->SetLineWidth(2);
	    f_1d_pion_total[pt_bin][cent][charge][eta_bin][flow_bin]->SetFillStyle(3002);
	    Double_t max_pion = f_1d_pion_total[pt_bin][cent][charge][eta_bin][flow_bin]->GetMaximum(x_low[pt_bin],x_up[pt_bin]);
	    sigma_pion[pt_bin][cent][charge][eta_bin][flow_bin] = parfit[pt_bin][cent][charge][eta_bin][flow_bin][1] - (Double_t)f_1d_pion_total[pt_bin][cent][charge][eta_bin][flow_bin]->GetX(max_pion/2.0,-1.0,parfit[pt_bin][cent][charge][eta_bin][flow_bin][1]);

	    TString KaonName = Form("f_1d_kaon_total_pt_%d_cent_%d_charge_%d_etagap_%d_%s",pt_bin,cent,charge,eta_bin,Order[flow_bin].Data());
	    f_1d_kaon_total[pt_bin][cent][charge][eta_bin][flow_bin] = new TF1(KaonName.Data(),student_t_1d_single,x_low[pt_bin],x_up[pt_bin],4);
	    f_1d_kaon_total[pt_bin][cent][charge][eta_bin][flow_bin]->SetNpx(400);
	    f_1d_kaon_total[pt_bin][cent][charge][eta_bin][flow_bin]->FixParameter(0,f_1d_double_total[pt_bin][cent][charge][eta_bin][flow_bin]->GetParameter(4));
	    f_1d_kaon_total[pt_bin][cent][charge][eta_bin][flow_bin]->FixParameter(1,f_1d_double_total[pt_bin][cent][charge][eta_bin][flow_bin]->GetParameter(5));
	    f_1d_kaon_total[pt_bin][cent][charge][eta_bin][flow_bin]->FixParameter(2,f_1d_double_total[pt_bin][cent][charge][eta_bin][flow_bin]->GetParameter(6));
	    f_1d_kaon_total[pt_bin][cent][charge][eta_bin][flow_bin]->FixParameter(3,f_1d_double_total[pt_bin][cent][charge][eta_bin][flow_bin]->GetParameter(7));
	    f_1d_kaon_total[pt_bin][cent][charge][eta_bin][flow_bin]->SetFillColor(kGray+2);
	    f_1d_kaon_total[pt_bin][cent][charge][eta_bin][flow_bin]->SetLineColor(kGray+2);
	    f_1d_kaon_total[pt_bin][cent][charge][eta_bin][flow_bin]->SetLineStyle(1);
	    f_1d_kaon_total[pt_bin][cent][charge][eta_bin][flow_bin]->SetLineWidth(2);
	    f_1d_kaon_total[pt_bin][cent][charge][eta_bin][flow_bin]->SetFillStyle(3003);
	    Double_t max_kaon = f_1d_kaon_total[pt_bin][cent][charge][eta_bin][flow_bin]->GetMaximum(x_low[pt_bin],x_up[pt_bin]);
	    sigma_kaon[pt_bin][cent][charge][eta_bin][flow_bin] = parfit[pt_bin][cent][charge][eta_bin][flow_bin][5] - (Double_t)f_1d_kaon_total[pt_bin][cent][charge][eta_bin][flow_bin]->GetX(max_kaon/2.0,0.0,parfit[pt_bin][cent][charge][eta_bin][flow_bin][5]);
	  }
	}
      }
    }
  }

  // Canvas
  TCanvas *c_X_pik_total[cent_total][charge_total][eta_total][flow_total];
  for(Int_t cent = cent_start[mCentrality]; cent < cent_stop[mCentrality]; cent++)
  {
    for(Int_t charge = mCharge; charge < mCharge+1; charge++)
    {
      for(Int_t eta_bin = eta_start; eta_bin < eta_stop; eta_bin++)
      {
	for(Int_t flow_bin = 0; flow_bin < flow_total; flow_bin++)
	{
	  TString CanName = Form("c_X_cent_%d_charge_%d_etagap_%d_%s_pik",cent,charge,eta_bin,Order[flow_bin].Data());
	  c_X_pik_total[cent][charge][eta_bin][flow_bin] = new TCanvas(CanName.Data(),CanName.Data(),1400,10,1500,900);
	  c_X_pik_total[cent][charge][eta_bin][flow_bin]->Divide(5,3);
	  for(Int_t pt_bin = 0; pt_bin < pt_total; pt_bin++)
	  {
	    c_X_pik_total[cent][charge][eta_bin][flow_bin]->cd(pt_bin+1);
	    c_X_pik_total[cent][charge][eta_bin][flow_bin]->cd(pt_bin+1)->SetLeftMargin(0.15);
	    c_X_pik_total[cent][charge][eta_bin][flow_bin]->cd(pt_bin+1)->SetRightMargin(0.05);
	    c_X_pik_total[cent][charge][eta_bin][flow_bin]->cd(pt_bin+1)->SetTopMargin(0.1);
	    c_X_pik_total[cent][charge][eta_bin][flow_bin]->cd(pt_bin+1)->SetBottomMargin(0.15);
	    c_X_pik_total[cent][charge][eta_bin][flow_bin]->cd(pt_bin+1)->SetLogz();
	    c_X_pik_total[cent][charge][eta_bin][flow_bin]->cd(pt_bin+1)->SetTicks(1,1);

	    h_X_total[pt_bin][cent][charge][eta_bin][flow_bin]->DrawCopy("l");
	    if(pt_bin >= 5)
	    {
	      f_1d_pion_total[pt_bin][cent][charge][eta_bin][flow_bin]->Draw("l same");
	      f_1d_kaon_total[pt_bin][cent][charge][eta_bin][flow_bin]->Draw("l same");
	      f_1d_double_total[pt_bin][cent][charge][eta_bin][flow_bin]->Draw("l same");
	      h_X_total[pt_bin][cent][charge][eta_bin][flow_bin]->DrawCopy("l same");
	      PlotLine(parfit_2nd[pt_bin][cent][charge][eta_bin][flow_bin][2]-order_pion[pt_bin]*parfit_2nd[pt_bin][cent][charge][eta_bin][flow_bin][4],parfit_2nd[pt_bin][cent][charge][eta_bin][flow_bin][2]-order_pion[pt_bin]*parfit_2nd[pt_bin][cent][charge][eta_bin][flow_bin][4],0,h_X_total[pt_bin][cent][charge][eta_bin][flow_bin]->GetMaximum()/1.5,1,2,2);
	      PlotLine(parfit_2nd[pt_bin][cent][charge][eta_bin][flow_bin][7]+order_kaon[pt_bin]*parfit_2nd[pt_bin][cent][charge][eta_bin][flow_bin][9],parfit_2nd[pt_bin][cent][charge][eta_bin][flow_bin][7]+order_kaon[pt_bin]*parfit_2nd[pt_bin][cent][charge][eta_bin][flow_bin][9],0,h_X_total[pt_bin][cent][charge][eta_bin][flow_bin]->GetMaximum()/1.5,1,2,2);
	    }
	  }
	}
      }
    }
  }
  //----------------------------------------------------------------------------------------------------------

  //----------------------------------------------------------------------------------------------------------
  // fit for differential phi-psi bin
  // Projection X
  TH1F *h_X_cut[pt_total][cent_total][charge_total][eta_total][phi_psi_total][flow_total];
  for(Int_t pt_bin = 0; pt_bin < pt_total; pt_bin++)
  {
    for(Int_t cent = cent_start[mCentrality]; cent < cent_stop[mCentrality]; cent++)
    {
      for(Int_t charge = mCharge; charge < mCharge+1; charge++)
      {
	for(Int_t eta_bin = eta_start; eta_bin < eta_stop; eta_bin++)
	{
	  for(Int_t phi_psi_bin = 0; phi_psi_bin < phi_psi_total; phi_psi_bin++)
	  {
	    for(Int_t flow_bin = 0; flow_bin < flow_total; flow_bin++)
	    {
	      if(pt_bin < 5)
	      {
		TString HistName = Form("pt_%d_Centrality_%d_charge_%d_etagap_%d_phi_psi_%d_%s_X_cut",pt_bin,cent,charge,eta_bin,phi_psi_bin,Order[flow_bin].Data());

		h_X_cut[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin] = (TH1F*)h_XY_cut[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->ProjectionX(HistName.Data(),h_XY_cut[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->GetYaxis()->FindBin(y_low[pt_bin]),h_XY_cut[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->GetYaxis()->FindBin(y_up[pt_bin]));
	      }
	      if(pt_bin >= 5)
	      {
		TString HistName = Form("pt_%d_Centrality_%d_charge_%d_etagap_%d_phi_psi_%d_%s_X_cut",pt_bin,cent,charge,eta_bin,phi_psi_bin,Order[flow_bin].Data());
		h_X_cut[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin] = (TH1F*)h_XY_cut[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->ProjectionX(HistName.Data(),h_XY_cut[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->GetYaxis()->FindBin(y_p[pt_bin]-order_proton[pt_bin]*parfit_2nd[pt_bin][cent][charge][eta_bin][flow_bin][17]),h_XY_cut[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->GetYaxis()->FindBin(y_pi[pt_bin]+order_pion[pt_bin]*parfit_2nd[pt_bin][cent][charge][eta_bin][flow_bin][5]));
	      }
	      h_X_cut[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->SetLineStyle(2);
	      h_X_cut[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->SetLineColor(1);
	      h_X_cut[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->SetLineWidth(1);

	      h_X_cut[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->SetStats(0);
	      h_X_cut[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->SetTitle(Title[pt_bin].Data());

	      h_X_cut[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->GetXaxis()->SetNdivisions(505,'N');
	      h_X_cut[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->GetXaxis()->SetTitle("x (n#sigma_{#pi},m^{2})");
	      h_X_cut[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->GetXaxis()->SetTitleSize(0.05);
	      h_X_cut[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->GetXaxis()->CenterTitle();
	      h_X_cut[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->GetXaxis()->SetTitleOffset(0.9);

	      h_X_cut[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->GetYaxis()->SetNdivisions(505,'N');
	      h_X_cut[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->GetYaxis()->SetTitle("Counts");
	      h_X_cut[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->GetYaxis()->SetTitleSize(0.05);
	      h_X_cut[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->GetYaxis()->CenterTitle();
	      h_X_cut[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->GetYaxis()->SetTitleOffset(1.5);
	    }
	  }
	}
      }
    }
  }

  TF1 *f_1d_double[pt_total][cent_total][charge_total][eta_total][phi_psi_total][flow_total];
  TF1 *f_1d_pion[pt_total][cent_total][charge_total][eta_total][phi_psi_total][flow_total];
  TF1 *f_1d_kaon[pt_total][cent_total][charge_total][eta_total][phi_psi_total][flow_total];
  for(Int_t pt_bin = 5; pt_bin < pt_total; pt_bin ++)
  {
    for(Int_t cent = cent_start[mCentrality]; cent < cent_stop[mCentrality]; cent++)
    {
      for(Int_t charge = mCharge; charge < mCharge+1; charge++)
      {
	for(Int_t eta_bin = eta_start; eta_bin < eta_stop; eta_bin++)
	{
	  for(Int_t phi_psi_bin = 0; phi_psi_bin < phi_psi_total; phi_psi_bin++)
	  {
	    for(Int_t flow_bin = 0; flow_bin < flow_total; flow_bin++)
	    {
	      TString FuncName = Form("f_1d_double_pt_%d_cent_%d_charge_%d_etagap_%d_phi_psi_%d_%s",pt_bin,cent,charge,eta_bin,phi_psi_bin,Order[flow_bin].Data());
	      f_1d_double[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin] = new TF1(FuncName.Data(),student_t_1d_double,x_low[pt_bin],x_up[pt_bin],8);
	      f_1d_double[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->SetNpx(400);
	      /*
	      if(pt_bin > 8)
	      {
		f_1d_double[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->SetNpx(200);
		if(pt_bin == 12)
		{
		  f_1d_double[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->SetNpx(100);
		}
	      }
	      */

	      // pion
	      f_1d_double[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->FixParameter(0,parfit[pt_bin][cent][charge][eta_bin][flow_bin][0]);
	      f_1d_double[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->FixParameter(1,parfit[pt_bin][cent][charge][eta_bin][flow_bin][1]);
	      f_1d_double[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->FixParameter(2,parfit[pt_bin][cent][charge][eta_bin][flow_bin][2]);
	      f_1d_double[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->SetParameter(3,parfit[pt_bin][cent][charge][eta_bin][flow_bin][3]);
//	      f_1d_double[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->SetParameter(3,1.0);

	      // kaon
	      f_1d_double[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->FixParameter(4,parfit[pt_bin][cent][charge][eta_bin][flow_bin][4]);
	      f_1d_double[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->FixParameter(5,parfit[pt_bin][cent][charge][eta_bin][flow_bin][5]);
	      f_1d_double[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->FixParameter(6,parfit[pt_bin][cent][charge][eta_bin][flow_bin][6]);
	      f_1d_double[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->SetParameter(7,parfit[pt_bin][cent][charge][eta_bin][flow_bin][7]);
//	      f_1d_double[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->SetParameter(7,1.0);

	      Float_t pi_max = h_X_cut[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->GetBinContent(h_X_cut[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->FindBin(parfit[pt_bin][cent][charge][eta_bin][flow_bin][1]));
	      Float_t pi_func = f_1d_double[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->Eval(parfit[pt_bin][cent][charge][eta_bin][flow_bin][1]);
	      Float_t Norm_pi = pi_max/pi_func;

	      Float_t k_max = h_X_cut[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->GetBinContent(h_X_cut[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->FindBin(parfit[pt_bin][cent][charge][eta_bin][flow_bin][5]));
	      Float_t k_func = f_1d_double[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->Eval(parfit[pt_bin][cent][charge][eta_bin][flow_bin][5]);
	      Float_t Norm_k = k_max/k_func;

	      f_1d_double[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->SetParameter(3,Norm_pi);
	      f_1d_double[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->SetParameter(7,Norm_k);

	      f_1d_double[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->SetLineStyle(1);
	      f_1d_double[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->SetLineColor(2);
	      f_1d_double[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->SetFillStyle(0);

	      f_1d_double[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->SetRange(parfit_2nd[pt_bin][cent][charge][eta_bin][flow_bin][2]-order_pion[pt_bin]*parfit_2nd[pt_bin][cent][charge][eta_bin][flow_bin][4],parfit_2nd[pt_bin][cent][charge][eta_bin][flow_bin][7]+order_kaon[pt_bin]*parfit_2nd[pt_bin][cent][charge][eta_bin][flow_bin][9]);

	      // Fit
//	      cout << "pt_bin = " << pt_bin << ", charge = " << charge << ", eta_bin = " << eta_bin << ", phi_psi_bin = " << phi_psi_bin << ", " << Order[flow_bin].Data() << endl;
	      h_X_cut[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->Fit(f_1d_double[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin],"MQRN");

	      Float_t chi2 = f_1d_double[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->GetChisquare();
	      Float_t Ndf = f_1d_double[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->GetNDF();
//	      cout << "chi2/Ndf = " << chi2 << "/" << Ndf << " = " << chi2/Ndf << endl;

	      f_1d_double[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->SetRange(x_low[pt_bin],x_up[pt_bin]);

	      TString PionName = Form("f_1d_pion_pt_%d_cent_%d_charge_%d_etagap_%d_phi_psi_%d_%s",pt_bin,cent,charge,eta_bin,phi_psi_bin,Order[flow_bin].Data());
	      f_1d_pion[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin] = new TF1(PionName.Data(),student_t_1d_single,x_low[pt_bin],x_up[pt_bin],4);
	      f_1d_pion[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->SetNpx(400);
	      f_1d_pion[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->FixParameter(0,f_1d_double[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->GetParameter(0));
	      f_1d_pion[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->FixParameter(1,f_1d_double[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->GetParameter(1));
	      f_1d_pion[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->FixParameter(2,f_1d_double[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->GetParameter(2));
	      f_1d_pion[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->FixParameter(3,f_1d_double[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->GetParameter(3));
	      f_1d_pion[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->SetFillColor(kAzure-2);
	      f_1d_pion[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->SetLineColor(kAzure-2);
	      f_1d_pion[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->SetLineStyle(1);
	      f_1d_pion[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->SetLineWidth(2);
	      f_1d_pion[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->SetFillStyle(3002);

	      TString KaonName = Form("f_1d_kaon_pt_%d_cent_%d_charge_%d_etagap_%d_phi_psi_%d_%s",pt_bin,cent,charge,eta_bin,phi_psi_bin,Order[flow_bin].Data());
	      f_1d_kaon[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin] = new TF1(KaonName.Data(),student_t_1d_single,x_low[pt_bin],x_up[pt_bin],4);
	      f_1d_kaon[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->SetNpx(400);
	      f_1d_kaon[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->FixParameter(0,f_1d_double[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->GetParameter(4));
	      f_1d_kaon[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->FixParameter(1,f_1d_double[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->GetParameter(5));
	      f_1d_kaon[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->FixParameter(2,f_1d_double[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->GetParameter(6));
	      f_1d_kaon[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->FixParameter(3,f_1d_double[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->GetParameter(7));
	      f_1d_kaon[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->SetFillColor(kGray+2);
	      f_1d_kaon[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->SetLineColor(kGray+2);
	      f_1d_kaon[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->SetLineStyle(1);
	      f_1d_kaon[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->SetLineWidth(2);
	      f_1d_kaon[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->SetFillStyle(3003);
	    }
	  }
	}
      }
    }
  }

  // Canvas
  TCanvas *c_X_pik[cent_total][charge_total][eta_total][phi_psi_total][flow_total];
  for(Int_t cent = cent_start[mCentrality]; cent < cent_stop[mCentrality]; cent++)
  {
    for(Int_t charge = mCharge; charge < mCharge+1; charge++)
    {
      for(Int_t eta_bin = eta_start; eta_bin < eta_stop; eta_bin++)
      {
	for(Int_t phi_psi_bin = 0; phi_psi_bin < phi_psi_total; phi_psi_bin++)
	{
	  for(Int_t flow_bin = 0; flow_bin < flow_total; flow_bin++)
	  {
	    TString CanName = Form("c_X_cent_%d_charge_%d_etagap_%d_phi_psi_%d_%s_pik",cent,charge,eta_bin,phi_psi_bin,Order[flow_bin].Data());
	    c_X_pik[cent][charge][eta_bin][phi_psi_bin][flow_bin] = new TCanvas(CanName.Data(),CanName.Data(),1400,10,1500,900);
	    c_X_pik[cent][charge][eta_bin][phi_psi_bin][flow_bin]->Divide(5,3);
	    for(Int_t pt_bin = 0; pt_bin < pt_total; pt_bin++)
	    {
	      c_X_pik[cent][charge][eta_bin][phi_psi_bin][flow_bin]->cd(pt_bin+1);
	      c_X_pik[cent][charge][eta_bin][phi_psi_bin][flow_bin]->cd(pt_bin+1)->SetLeftMargin(0.15);
	      c_X_pik[cent][charge][eta_bin][phi_psi_bin][flow_bin]->cd(pt_bin+1)->SetRightMargin(0.05);
	      c_X_pik[cent][charge][eta_bin][phi_psi_bin][flow_bin]->cd(pt_bin+1)->SetTopMargin(0.1);
	      c_X_pik[cent][charge][eta_bin][phi_psi_bin][flow_bin]->cd(pt_bin+1)->SetBottomMargin(0.15);
	      c_X_pik[cent][charge][eta_bin][phi_psi_bin][flow_bin]->cd(pt_bin+1)->SetLogz();
	      c_X_pik[cent][charge][eta_bin][phi_psi_bin][flow_bin]->cd(pt_bin+1)->SetTicks(1,1);

	      h_X_cut[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->DrawCopy("l");
	      if(pt_bin >= 5)
	      {
		f_1d_pion[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->Draw("l same");
		f_1d_kaon[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->Draw("l same");
		f_1d_double[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->Draw("l same");
		h_X_cut[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->DrawCopy("l same");
		PlotLine(parfit[pt_bin][cent][charge][eta_bin][flow_bin][1]-order_pion[pt_bin]*parfit[pt_bin][cent][charge][eta_bin][flow_bin][2],parfit[pt_bin][cent][charge][eta_bin][flow_bin][1]-order_pion[pt_bin]*parfit[pt_bin][cent][charge][eta_bin][flow_bin][2],0,h_X_cut[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->GetMaximum()/1.5,1,2,2);
		PlotLine(parfit[pt_bin][cent][charge][eta_bin][flow_bin][5]+order_pion[pt_bin]*parfit[pt_bin][cent][charge][eta_bin][flow_bin][6],parfit[pt_bin][cent][charge][eta_bin][flow_bin][5]+order_pion[pt_bin]*parfit[pt_bin][cent][charge][eta_bin][flow_bin][6],0,h_X_cut[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->GetMaximum()/1.5,1,2,2);
	      }
	    }
	  }
	}
      }
    }
  }

  // extract the bins of pion and kaon
  TH1F *h_X_pion[pt_total][cent_total][charge_total][eta_total][phi_psi_total][flow_total]; // pion after subtract kaon
  TH1F *h_X_kaon[pt_total][cent_total][charge_total][eta_total][phi_psi_total][flow_total]; // kaon after subtract pion
  for(Int_t pt_bin = 0; pt_bin < pt_total; pt_bin++)
  {
    for(Int_t cent = cent_start[mCentrality]; cent < cent_stop[mCentrality]; cent++)
    {
      for(Int_t charge = mCharge; charge < mCharge+1; charge++)
      {
	for(Int_t eta_bin = eta_start; eta_bin < eta_stop; eta_bin++)
	{
	  for(Int_t phi_psi_bin = 0; phi_psi_bin < phi_psi_total; phi_psi_bin++)
	  {
	    for(Int_t flow_bin = 0; flow_bin < flow_total; flow_bin++)
	    {
	      if(pt_bin <= 4)
	      {
		bin_pion_start[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin] = h_X_cut[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->FindBin(-0.05);
		bin_pion_stop[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]  = h_X_cut[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->FindBin(0.05);
		bin_kaon_start[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin] = h_X_cut[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->FindBin(0.15);
		bin_kaon_stop[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]  = h_X_cut[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->FindBin(0.30);
	      }
	      if(pt_bin >= 5)
	      {
		TString HistName_pion = Form("pt_%d_Centrality_%d_charge_%d_etagap_%d_phi_psi_%d_%s_X_pion",pt_bin,cent,charge,eta_bin,phi_psi_bin,Order[flow_bin].Data());
		h_X_pion[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin] = (TH1F*)h_X_cut[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->Clone(HistName_pion.Data());
		TString HistName_kaon = Form("pt_%d_Centrality_%d_charge_%d_etagap_%d_phi_psi_%d_%s_X_kaon",pt_bin,cent,charge,eta_bin,phi_psi_bin,Order[flow_bin].Data());
		h_X_kaon[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin] = (TH1F*)h_X_cut[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->Clone(HistName_kaon.Data());


		h_X_pion[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->Add(f_1d_kaon[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin],-1.0); // pion after subtract kaon
		h_X_kaon[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->Add(f_1d_pion[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin],-1.0); // kaon after subtract pion

		bin_pion_start[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin] = h_X_pion[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->FindBin(parfit[pt_bin][cent][charge][eta_bin][flow_bin][1]-2.0*sigma_pion[pt_bin][cent][charge][eta_bin][flow_bin]);
		bin_pion_stop[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]  = h_X_pion[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->FindBin(parfit[pt_bin][cent][charge][eta_bin][flow_bin][1]+2.0*sigma_pion[pt_bin][cent][charge][eta_bin][flow_bin]);

		bin_kaon_start[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin] = h_X_kaon[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->FindBin(parfit[pt_bin][cent][charge][eta_bin][flow_bin][5]-2.0*sigma_kaon[pt_bin][cent][charge][eta_bin][flow_bin]);
		bin_kaon_stop[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]  = h_X_kaon[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->FindBin(parfit[pt_bin][cent][charge][eta_bin][flow_bin][5]+2.0*sigma_kaon[pt_bin][cent][charge][eta_bin][flow_bin]);
	      }
//	      cout << "pt_bin = " << pt_bin << endl;
//	      cout << "pion_start = " << bin_pion_start[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin] << ", pion_stop = " << bin_pion_stop[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin] << endl;
//	      cout << "kaon_start = " << bin_kaon_start[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin] << ", kaon_stop = " << bin_kaon_stop[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin] << endl;
	    }
	  }
	}
      }
    }
  }

  /*
  // Canvas
  TCanvas *c_X_pi[cent_total][charge_total][eta_total][phi_psi_total][flow_total];
  TCanvas *c_X_k[cent_total][charge_total][eta_total][phi_psi_total][flow_total];
  for(Int_t cent = cent_start[mCentrality]; cent < cent_stop[mCentrality]; cent++)
  {
    for(Int_t charge = mCharge; charge < mCharge+1; charge++)
    {
      for(Int_t eta_bin = eta_start; eta_bin < eta_stop; eta_bin++)
      {
	for(Int_t phi_psi_bin = 0; phi_psi_bin < phi_psi_total; phi_psi_bin++)
	{
	  for(Int_t flow_bin = 0; flow_bin < flow_total; flow_bin++)
	  {
	    TString CanName_pi = Form("c_X_cent_%d_charge_%d_etagap_%d_phi_psi_%d_%s_pi",cent,charge,eta_bin,phi_psi_bin,Order[flow_bin].Data());
	    c_X_pi[cent][charge][eta_bin][phi_psi_bin][flow_bin] = new TCanvas(CanName_pi.Data(),CanName_pi.Data(),1400,10,1500,900);
	    c_X_pi[cent][charge][eta_bin][phi_psi_bin][flow_bin]->Divide(5,3);
	    TString CanName_k = Form("c_X_cent_%d_charge_%d_etagap_%d_phi_psi_%d_%s_k",cent,charge,eta_bin,phi_psi_bin,Order[flow_bin].Data());
	    c_X_k[cent][charge][eta_bin][phi_psi_bin][flow_bin] = new TCanvas(CanName_k.Data(),CanName_k.Data(),1400,10,1500,900);
	    c_X_k[cent][charge][eta_bin][phi_psi_bin][flow_bin]->Divide(5,3);
	    for(Int_t pt_bin = 5; pt_bin < pt_total; pt_bin++)
	    {
	      c_X_pi[cent][charge][eta_bin][phi_psi_bin][flow_bin]->cd(pt_bin+1);
	      c_X_pi[cent][charge][eta_bin][phi_psi_bin][flow_bin]->cd(pt_bin+1)->SetLeftMargin(0.15);
	      c_X_pi[cent][charge][eta_bin][phi_psi_bin][flow_bin]->cd(pt_bin+1)->SetRightMargin(0.05);
	      c_X_pi[cent][charge][eta_bin][phi_psi_bin][flow_bin]->cd(pt_bin+1)->SetTopMargin(0.1);
	      c_X_pi[cent][charge][eta_bin][phi_psi_bin][flow_bin]->cd(pt_bin+1)->SetBottomMargin(0.15);
	      c_X_pi[cent][charge][eta_bin][phi_psi_bin][flow_bin]->cd(pt_bin+1)->SetLogz();
	      c_X_pi[cent][charge][eta_bin][phi_psi_bin][flow_bin]->cd(pt_bin+1)->SetTicks(1,1);

	      h_X_pion[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->DrawCopy("h");
	      PlotLine(parfit[pt_bin][cent][charge][eta_bin][flow_bin][1]-2.0*sigma_pion[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin],parfit[pt_bin][cent][charge][eta_bin][flow_bin][1]-2.0*sigma_pion[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin],0,h_X_pion[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->GetMaximum()/1,1,2,2);
	      PlotLine(parfit[pt_bin][cent][charge][eta_bin][flow_bin][1]+2.0*sigma_pion[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin],parfit[pt_bin][cent][charge][eta_bin][flow_bin][1]+2.0*sigma_pion[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin],0,h_X_pion[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->GetMaximum()/1,1,2,2);

	      c_X_k[cent][charge][eta_bin][phi_psi_bin][flow_bin]->cd(pt_bin+1);
	      c_X_k[cent][charge][eta_bin][phi_psi_bin][flow_bin]->cd(pt_bin+1)->SetLeftMargin(0.15);
	      c_X_k[cent][charge][eta_bin][phi_psi_bin][flow_bin]->cd(pt_bin+1)->SetRightMargin(0.05);
	      c_X_k[cent][charge][eta_bin][phi_psi_bin][flow_bin]->cd(pt_bin+1)->SetTopMargin(0.1);
	      c_X_k[cent][charge][eta_bin][phi_psi_bin][flow_bin]->cd(pt_bin+1)->SetBottomMargin(0.15);
	      c_X_k[cent][charge][eta_bin][phi_psi_bin][flow_bin]->cd(pt_bin+1)->SetLogz();
	      c_X_k[cent][charge][eta_bin][phi_psi_bin][flow_bin]->cd(pt_bin+1)->SetTicks(1,1);

	      h_X_kaon[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->DrawCopy("h");
	      PlotLine(parfit[pt_bin][cent][charge][eta_bin][flow_bin][5]-2.0*sigma_kaon[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin],parfit[pt_bin][cent][charge][eta_bin][flow_bin][5]-2.0*sigma_kaon[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin],0,h_X_kaon[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->GetMaximum()/1,1,2,2);
	      PlotLine(parfit[pt_bin][cent][charge][eta_bin][flow_bin][5]+2.0*sigma_kaon[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin],parfit[pt_bin][cent][charge][eta_bin][flow_bin][5]+2.0*sigma_kaon[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin],0,h_X_kaon[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->GetMaximum()/1,1,2,2);
	    }
	  }
	}
      }
    }
  }
  */

  // extract the yield and error of pion and kaon
  for(Int_t pt_bin = 0; pt_bin < pt_total; pt_bin++)
  {
    for(Int_t cent = cent_start[mCentrality]; cent < cent_stop[mCentrality]; cent++)
    {
      for(Int_t charge = mCharge; charge < mCharge+1; charge++)
      {
	for(Int_t eta_bin = eta_start; eta_bin < eta_stop; eta_bin++)
	{
	  for(Int_t phi_psi_bin = 0; phi_psi_bin < phi_psi_total; phi_psi_bin++)
	  {
	    for(Int_t flow_bin = 0; flow_bin < flow_total; flow_bin++)
	    {
	      // pion
	      Int_t pion_start = bin_pion_start[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin];
	      Int_t pion_stop  = bin_pion_stop[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin] + 1;
	      counts_pion[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin] = 0.0;
	      Float_t Err_pion  = 0.0;
	      if(pt_bin <= 4)
	      {
		for(Int_t bin_pion = pion_start; bin_pion < pion_stop; bin_pion++)
		{
		  counts_pion[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin] += h_X_cut[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->GetBinContent(bin_pion);
		  Err_pion += h_X_cut[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->GetBinError(bin_pion)*h_X_cut[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->GetBinError(bin_pion);
		}
	      }
	      if(pt_bin >= 5)
	      {
		for(Int_t bin_pion = pion_start; bin_pion < pion_stop; bin_pion++)
		{
		  counts_pion[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin] += h_X_pion[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->GetBinContent(bin_pion);
		  Err_pion += h_X_pion[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->GetBinError(bin_pion)*h_X_pion[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->GetBinError(bin_pion);
		}
	      }
	      errors_pion[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin] = TMath::Sqrt(Err_pion);

	      // kaon
	      Int_t kaon_start = bin_kaon_start[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin];
	      Int_t kaon_stop  = bin_kaon_stop[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin] + 1;
	      counts_kaon[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin] = 0.0;
	      Float_t Err_kaon  = 0.0;
	      if(pt_bin <= 4)
	      {
		for(Int_t bin_kaon = kaon_start; bin_kaon < kaon_stop; bin_kaon++)
		{
		  counts_kaon[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin] += h_X_cut[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->GetBinContent(bin_kaon);
		  Err_kaon += h_X_cut[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->GetBinError(bin_kaon)*h_X_cut[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->GetBinError(bin_kaon);
		}
	      }
	      if(pt_bin >= 5)
	      {
		for(Int_t bin_kaon = kaon_start; bin_kaon < kaon_stop; bin_kaon++)
		{
		  counts_kaon[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin] += h_X_kaon[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->GetBinContent(bin_kaon);
		  Err_kaon += h_X_kaon[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->GetBinError(bin_kaon)*h_X_kaon[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->GetBinError(bin_kaon);
		}
	      }
	      errors_kaon[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin] = TMath::Sqrt(Err_kaon);
	    }
	  }
	}
      }
    }
  }

  TH1F *h_counts_pion[pt_total][cent_total][charge_total][eta_total][flow_total]; // counts of pion
  TH1F *h_counts_kaon[pt_total][cent_total][charge_total][eta_total][flow_total]; // counts of kaon
  for(Int_t pt_bin = 0; pt_bin < pt_total; pt_bin++)
  {
    for(Int_t cent = cent_start[mCentrality]; cent < cent_stop[mCentrality]; cent++)
    {
      for(Int_t charge = mCharge; charge < mCharge+1; charge++)
      {
	for(Int_t eta_bin = eta_start; eta_bin < eta_stop; eta_bin++)
	{
	  for(Int_t flow_bin = 0; flow_bin < flow_total; flow_bin++)
	  {
	    TString HistName_pion = Form("pt_%d_Centrality_%d_charge_%d_etagap_%d_%s_counts_pion",pt_bin,cent,charge,eta_bin,Order[flow_bin].Data());
	    h_counts_pion[pt_bin][cent][charge][eta_bin][flow_bin] = new TH1F(HistName_pion.Data(),HistName_pion.Data(),7,0,PI_max[flow_bin]); // counts of pion
	    TString HistName_kaon = Form("pt_%d_Centrality_%d_charge_%d_etagap_%d_%s_counts_kaon",pt_bin,cent,charge,eta_bin,Order[flow_bin].Data());
	    h_counts_kaon[pt_bin][cent][charge][eta_bin][flow_bin] = new TH1F(HistName_kaon.Data(),HistName_kaon.Data(),7,0,PI_max[flow_bin]); // counts of kaon
	    for(Int_t phi_psi_bin = 0; phi_psi_bin < phi_psi_total; phi_psi_bin++)
	    {
	      h_counts_pion[pt_bin][cent][charge][eta_bin][flow_bin]->SetBinContent(phi_psi_bin+1,counts_pion[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]);
	      h_counts_pion[pt_bin][cent][charge][eta_bin][flow_bin]->SetBinError(phi_psi_bin+1,errors_pion[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]);
	      h_counts_kaon[pt_bin][cent][charge][eta_bin][flow_bin]->SetBinContent(phi_psi_bin+1,counts_kaon[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]);
	      h_counts_kaon[pt_bin][cent][charge][eta_bin][flow_bin]->SetBinError(phi_psi_bin+1,errors_kaon[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]);
	    }
	  }
	}
      }
    }
  }
  TString output = Form("/project/projectdirs/star/xusun/OutPut/AuAu%s/Mass2_nSigmaPion/merged_file/flow_pik/Counts_%s_charge_%d_etagap_%d.root",Energy[mEnergy].Data(),Centrality[mCentrality].Data(),mCharge,eta_start);
  TFile *File_counts = new TFile(output.Data(),"RECREATE");
  File_counts->cd();
  for(Int_t pt_bin = 0; pt_bin < pt_total; pt_bin++)
  {
    for(Int_t cent = cent_start[mCentrality]; cent < cent_stop[mCentrality]; cent++)
    {
      for(Int_t charge = mCharge; charge < mCharge+1; charge++)
      {
	for(Int_t eta_bin = eta_start; eta_bin < eta_stop; eta_bin++)
	{
	  for(Int_t flow_bin = 0; flow_bin < flow_total; flow_bin++)
	  {
	    h_counts_pion[pt_bin][cent][charge][eta_bin][flow_bin]->Write();
	    h_counts_kaon[pt_bin][cent][charge][eta_bin][flow_bin]->Write();
	  }
	}
      }
    }
  }
  File_counts->Close();
  File_hist->Close();
  File_par->Close();
}

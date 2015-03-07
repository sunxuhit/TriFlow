#include "TFile.h"
#include "TString.h"
#include "TProfile.h"
#include "TH1F.h"
#include "TF1.h"
#include "Math/MinimizerOptions.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"

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

Double_t PtFitFunc2_mod_x(Double_t* x_val, Double_t* par)
{
  Double_t x, y, m0, Temp, Ampl, shift;
  m0    = par[0];
  Temp  = par[1];
  Ampl  = par[2];
  shift = par[3];
  x = x_val[0];
  y = x*(Ampl*(x-shift)*sqrt((x-shift)*(x-shift)+m0*m0)*TMath::Exp(-(sqrt((x-shift)*(x-shift)+m0*m0)-m0)/Temp)) + par[4]*x + par[5]*x*x;
  return y;
}

//Blast Wave Fit for spectra
Double_t SpectraFunc(Double_t* x_val, Double_t* par)
{
  Double_t pt, alpha, beta, rho, rho_0, rho_a, phi,phi_p, T, s2, R, Inte, Norm;
  pt = x_val[0];
  Int_t PID = (Int_t)par[5];
  Float_t mass[10] = {0.1396,0.1396,0.4937,0.4937,0.9382,0.9382,1.116,1.116,0.4976,1.019};
  Double_t mt = TMath::Sqrt(pt*pt + mass[PID]*mass[PID]);
  Int_t nbins_phi = 10;
  Double_t phi_start = 0.0;
  Double_t phi_stop = 2.0*TMath::Pi();
  Double_t delta_phi = (phi_stop - phi_start)/((Double_t)nbins_phi);
  T = par[0];
  rho_0 = par[1];
  rho_a = par[2];
  rho_a = 0.0;
  s2 = par[3];
  s2 = 0.0;
  R = par[4];
  Norm = par[6];
  Double_t r;
  Int_t nbins_r = 10;
  Double_t delta_r = R/((Double_t)nbins_r);
  Double_t n_power = 2./3.;
  
  Inte = 0.0;

  for(Int_t i = 0; i < nbins_phi + 1; i++)
  {
    phi = phi_start + i*delta_phi;
    for(Int_t j = 0; j < nbins_r + 1 ; j++)
    {
      r = j*delta_r;
      rho = TMath::ATanH(TMath::TanH(rho_0)*TMath::Power((r/R),n_power)) + TMath::ATanH(TMath::TanH(rho_a)*TMath::Power((r/R),n_power))*TMath::Cos(2.0*phi);
      //    rho = TMath::ATanH(TMath::TanH(rho_0)*(r/R)) + TMath::ATanH(TMath::TanH(rho_a)*(r/R))*TMath::Cos(2.0*phi);
      alpha = (pt/T)*TMath::SinH(rho);
      beta = (mt/T)*TMath::CosH(rho);

      Inte += Norm*pt*r*mt*delta_phi*delta_r*TMath::BesselK1(beta)*(1.0 + 2.0*s2*TMath::Cos(2.0*phi))*TMath::BesselI0(alpha);
    }    
  }
  return Inte;
}

static TString Mode[2] = {"Default","StringMelting"};
static TString ScreenMass[3] = {"1mb","3mb","6mb"};
static TString Energy[7] = {"7GeV","11GeV","19GeV","27GeV","39GeV","62GeV","200GeV"};
static TString Order[2] = {"2nd","3rd"};
static TString Centrality[4] = {"0080","0010","1040","4080"};
static TString ParType[10] = {"pi_plus","pi_minus","K_plus","K_minus","p","pbar","Lambda","Lambdabar","K0s","phi"};
static Float_t ParMass[10] = {0.1396,0.1396,0.4937,0.4937,0.9382,0.9382,1.116,1.116,0.4976,1.019};
static Int_t ParOrder[5] = {1,3,5,7,9};
static Int_t ParStyle[5] = {24,26,32,30,28};
static Int_t ParColor[5] = {1,kRed,kAzure+4,kOrange+7,kGray+3};
static TString ParName[5] = {"#pi^{-}","K^{-}","#bar{p}","#bar{#Lambda}","#phi"};

static const Int_t mOrder_Total = 2;
static Int_t mOrder_Start = 0;
static Int_t mOrder_Stop  = 1;

static const Int_t mCentrality_Total = 4;
static Int_t mCentrality_Start = 0;
static Int_t mCentrality_Stop  = 1;

static const Int_t mParType_Total = 10;
static Int_t mParType_Start = 0;
static Int_t mParType_Stop  = 3;

// Calculate integrated v2 and v3
void InteFlow(Int_t mEnergy = 4, Int_t mMode = 0, Int_t mScreen = 0) // 0: 7.7 GeV, 1: 11.5 GeV, 2: 19.6 GeV, 3: 27 GeV, 4: 39 GeV, 5: 62.4 GeV, 6: 200 GeV | 0: Default, 1: String Melting | 0: 1mb, 1: 3mb, 2: 6mb
{
//  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(50000);
  //---------------------Constant------------------------
  const Float_t pt_start = 1.5; // start point for pt integration
  const Float_t pt_stop  = 3.0; // stop point for pt integration
  //---------------------Constant------------------------

  TString inputfile;
  if(mMode == 0)
  {
    inputfile = Form("/home/xusun/Data/AMPT_%s/Flow/%s_%s/Flow_%s.root",Mode[mMode].Data(),Energy[mEnergy].Data(),Mode[mMode].Data(),Energy[mEnergy].Data());
  }
  if(mMode == 1)
  {
    inputfile = Form("/home/xusun/Data/AMPT_%s/Flow/%s_%s/%s/Flow_%s.root",Mode[mMode].Data(),Energy[mEnergy].Data(),Mode[mMode].Data(),ScreenMass[mScreen].Data(),Energy[mEnergy].Data());
  }
  cout << "Input File: " << inputfile.Data() << endl;
  TFile *File_input = TFile::Open(inputfile.Data());


  // flow for pi, K, p, Lambda, K0s by using eta_sub event plane method
  // pi_plus,pi_minus,K_plus,K_minus,p,pbar,Lambda,Lambdabar,K0s,phi
  TProfile *p_mFlow[mParType_Total][mOrder_Total][mCentrality_Total]; // Particle types | 0 for 2nd, 1 for 3rd | 0 for 0-80%, 1 for 0-10%, 2 for 10-40%, 3 for 40-80%
  TH1F *h_mFlow[mParType_Total][mOrder_Total][mCentrality_Total]; // Particle types | 0 for 2nd, 1 for 3rd | 0 for 0-80%, 1 for 0-10%, 2 for 10-40%, 3 for 40-80%

  // pt spectra 
  TH1F *h_mPt[mParType_Total][mOrder_Total][mCentrality_Total]; // Particle types | 0 for 2nd, 1 for 3rd | 0 for 0-80%, 1 for 0-10%, 2 for 10-40%, 3 for 40-80%
//  TF1 *f_mPt[mParType_Total][mOrder_Total][mCentrality_Total];

  // Integrated v2 and v3
  TH1F *h_mInteFlow[mOrder_Total][mCentrality_Total]; // 0 for 2nd, 1 for 3rd | 0 for 0-80%, 1 for 0-10%, 2 for 10-40%, 3 for 40-80%


  // Read in flow TProfile => transfer to TH1F | Read in pt Spectra
  for(Int_t i_par = 0; i_par < 10; i_par++)
  {
    for(Int_t i_order = 0; i_order < 2; i_order++)
    {
      for(Int_t i_cent = 0; i_cent < 4; i_cent++)
      {
	TString ProName;
	TString HistName;

	// v2 and v3 relative to event plane
	//--------------------------------------------------------------------------------------------------
	ProName = Form("Flow_%s_%s_%s",ParType[i_par].Data(),Order[i_order].Data(),Centrality[i_cent].Data()); // pi_plus
	p_mFlow[i_par][i_order][i_cent] = (TProfile*)File_input->Get(ProName.Data());
	HistName = Form("h_Flow_%s_%s_%s",ParType[i_par].Data(),Order[i_order].Data(),Centrality[i_cent].Data());
	h_mFlow[i_par][i_order][i_cent] = new TH1F(HistName.Data(),HistName.Data(),25,0.1,5.1);
	for(Int_t i_bin = 1; i_bin < 26; i_bin++)
	{
	  h_mFlow[i_par][i_order][i_cent]->SetBinContent(i_bin,p_mFlow[i_par][i_order][i_cent]->GetBinContent(i_bin));
	  h_mFlow[i_par][i_order][i_cent]->SetBinError(i_bin,p_mFlow[i_par][i_order][i_cent]->GetBinError(i_bin));
	}

	HistName = Form("Pt_%s_%s_%s",ParType[i_par].Data(),Order[i_order].Data(),Centrality[i_cent].Data());
	h_mPt[i_par][i_order][i_cent] = (TH1F*)File_input->Get(HistName.Data());

	h_mFlow[i_par][i_order][i_cent]->Multiply(h_mPt[i_par][i_order][i_cent]);
	//--------------------------------------------------------------------------------------------------
      }
    }
  }

  // Calculate Integrated v2 and v3 => could be merged into loop above => using individual loop to be more clear
  for(Int_t i_order = 0; i_order < 2; i_order++)
  {
    for(Int_t i_cent = 0; i_cent < 4; i_cent++)
    {
      TString HistName = Form("h_mInteFlow_%s_%s",Order[i_order].Data(),Centrality[i_cent].Data());
      h_mInteFlow[i_order][i_cent] = new TH1F(HistName.Data(),HistName.Data(),10,-0.5,9.5);
      for(Int_t i_par = 0; i_par < 10; i_par++)
      {
	Double_t Inte_flow = 0.0;
	Double_t Err_flow  = 0.0;
	Double_t Inte_spec = 0.0;
	Double_t Err_spec  = 0.0;
	Inte_flow = h_mFlow[i_par][i_order][i_cent]->IntegralAndError(h_mFlow[i_par][i_order][i_cent]->FindBin(pt_start),h_mFlow[i_par][i_order][i_cent]->FindBin(pt_stop),Err_flow,"width");
	Inte_spec = h_mPt[i_par][i_order][i_cent]->IntegralAndError(h_mPt[i_par][i_order][i_cent]->FindBin(pt_start),h_mPt[i_par][i_order][i_cent]->FindBin(pt_stop),Err_spec,"width");
	if(Inte_spec != 0.0)
	{
	  Double_t Err_ratio = ErrDiv(Inte_flow,Inte_spec,Err_flow,Err_spec);
	  h_mInteFlow[i_order][i_cent]->SetBinContent(i_par+1,Inte_flow/Inte_spec);
	  h_mInteFlow[i_order][i_cent]->SetBinError(i_par+1,Err_ratio);
	}
      }
    }
  }

  // Calculate ratio of integrated v2 and v3
  TH1F *h_ratio[4];
  for(Int_t i_cent = 0; i_cent < 4; i_cent++)
  {
    TString HistName = Form("h_ratio_%s",Centrality[i_cent].Data());
    h_ratio[i_cent] = (TH1F*)h_mInteFlow[0][i_cent]->Clone(HistName.Data());
    h_ratio[i_cent]->Divide(h_mInteFlow[1][i_cent]);
  }

  // write Histogram into Output file
  TString outputfile;
  if(mMode == 0)
  {
    outputfile = Form("/home/xusun/Data/AMPT_%s/InteFlow/%s_%s/InteFlow.root",Mode[mMode].Data(),Energy[mEnergy].Data(),Mode[mMode].Data());
  }
  if(mMode == 1)
  {
    outputfile = Form("/home/xusun/Data/AMPT_%s/InteFlow/%s_%s/%s/InteFlow.root",Mode[mMode].Data(),Energy[mEnergy].Data(),Mode[mMode].Data(),ScreenMass[mScreen].Data());
  }
  cout << "OutPut File: " << outputfile.Data() << endl;
  TFile *File_output = new TFile(outputfile.Data(),"RECREATE");
  File_output->cd();
  for(Int_t i_cent = 0; i_cent < 4; i_cent++)
  {
    h_ratio[i_cent]->Write();
    for(Int_t i_order = 0; i_order < 2; i_order++)
    {
      h_mInteFlow[i_order][i_cent]->Write();
    }
  }

  File_output->Close();
}

#include "TString.h"
#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"

//----------------------------------------------------------------------------------------
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

// Boltzmann Distribution
Double_t Boltzmann(Double_t* x_val, Double_t* par)
{
  Double_t x, y, m0, Temp, Ampl, pol1, pol2;
  m0    = par[0];
  Temp  = par[1];
  Ampl  = par[2];
  pol1  = par[3];
  pol2  = par[4];

  x = x_val[0];
  Double_t alpha = TMath::Power((m0/(2*TMath::Pi()*Temp)),1.5);
  Double_t Boltz = x*Ampl*alpha*4*TMath::Pi()*x*x*TMath::Exp(-0.5*m0*x*x/Temp);
  Double_t Poly3 = pol1*x + pol2*x*x;

  y = Boltz + Poly3;

  return y;
}

Double_t Gaussian(Double_t *x_val, Double_t *par)
{
  Double_t mu    = par[0];
  Double_t sigma = par[1];
  Double_t Norm  = par[2];
  Double_t x     = x_val[0];
  
  Double_t y = Norm*TMath::Exp(-0.5*(x-mu)*(x-mu)/(sigma*sigma))/sigma;

  return y;
}

Double_t Poly(Double_t *x_val, Double_t *par)
{
  Double_t x = x_val[0];
  Double_t y = par[0] + par[1]*TMath::Power(x,1) + par[2]*TMath::Power(x,2) + par[3]*TMath::Power(x,3) + par[4]*TMath::Power(x,4) + par[5]*TMath::Power(x,5) + par[6]*TMath::Power(x,6) + par[7]*TMath::Power(x,7);

  return y;
}
//----------------------------------------------------------------------------------------

static TString Mode[2] = {"Default","StringMelting"};
static TString Energy[7] = {"7GeV","11GeV","19GeV","27GeV","39GeV","62GeV","200GeV"};
static TString Order[2] = {"2nd","3rd"};
static TString Centrality[4] = {"0080","0010","1040","4080"};
static TString ParType[10] = {"pi_plus","pi_minus","K_plus","K_minus","p","pbar","Lambda","Lambdabar","K0s","phi"};
static Float_t ParMass[10] = {0.1396,0.1396,0.4937,0.4937,0.9382,0.9382,1.116,1.116,0.4976,1.019};

static const Int_t mOrder_Total = 2;
static Int_t mOrder_Start = 0;
static Int_t mOrder_Stop  = 1;

static const Int_t mCentrality_Total = 4;
static Int_t mCentrality_Start = 2;
static Int_t mCentrality_Stop  = 3;

static const Int_t mParType_Total = 10;
static Int_t mParType_Start = 7;
static Int_t mParType_Stop  = 8;

void SpecFit(Int_t mEnergy = 4, Int_t mMode = 0) // 0: 7.7 GeV, 1: 11.5 GeV, 2: 19.6 GeV, 3: 27 GeV, 4: 39 GeV, 5: 62.4 GeV, 6: 200 GeV | 0: Default, 1: String Melting
{
  TString inputfile = Form("/home/xusun/Data/AMPT_%s/Flow/%s_%s/Flow_%s.root",Mode[mMode].Data(),Energy[mEnergy].Data(),Mode[mMode].Data(),Energy[mEnergy].Data());
  cout << "Input File: " << inputfile.Data() << endl;
  TFile *File_input = TFile::Open(inputfile.Data());

  TH1D *h_spec[mOrder_Total][mCentrality_Total][mParType_Total];
  TF1 *f_spec[mOrder_Total][mCentrality_Total][mParType_Total];
  for(Int_t i_order = mOrder_Start; i_order < mOrder_Stop; i_order++)
  {
    for(Int_t i_cent = mCentrality_Start; i_cent < mCentrality_Stop; i_cent++)
    {
      for(Int_t i_par = mParType_Start; i_par < mParType_Stop; i_par++)
      {
	TString HistName = Form("Pt_%s_%s_%s",ParType[i_par].Data(),Order[i_order].Data(),Centrality[i_cent].Data());
	h_spec[i_order][i_cent][i_par] = (TH1D*)File_input->Get(HistName.Data());
	TString FuncName = Form("f_Pt_%s_%s_%s",ParType[i_par].Data(),Order[i_order].Data(),Centrality[i_cent].Data());
	f_spec[i_order][i_cent][i_par] = new TF1(FuncName.Data(),PtFitFunc2_mod_x,0,5.0,4);
	f_spec[i_order][i_cent][i_par]->SetParameter(0,ParMass[i_par]);
	f_spec[i_order][i_cent][i_par]->SetParameter(1,0.2);
	f_spec[i_order][i_cent][i_par]->SetParameter(2,10000);
	f_spec[i_order][i_cent][i_par]->SetParameter(3,1.0);
	f_spec[i_order][i_cent][i_par]->SetRange(2.0,5.0);
	h_spec[i_order][i_cent][i_par]->Fit(f_spec[i_order][i_cent][i_par],"NRMI");
      }
    }
  }
  h_spec[0][2][7]->SetLineColor(4);
  h_spec[0][2][7]->Draw("pE");
  f_spec[0][2][7]->SetLineColor(2);
  f_spec[0][2][7]->Draw("l same");
}

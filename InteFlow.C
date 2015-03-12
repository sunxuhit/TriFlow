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

static TString Mode[2] = {"Default","StringMelting"};
static TString ScreenMass[3] = {"1mb","3mb","6mb"};
static TString Energy[7] = {"7GeV","11GeV","19GeV","27GeV","39GeV","62GeV","200GeV"};
static TString Order[2] = {"2nd","3rd"};
static TString Centrality[4] = {"0080","0010","1040","4080"};
static TString ParType[10] = {"pi_plus","pi_minus","K_plus","K_minus","p","pbar","Lambda","Lambdabar","K0s","phi"};
static Float_t ParMass[10] = {0.1396,0.1396,0.4937,0.4937,0.9382,0.9382,1.116,1.116,0.4976,1.019};

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
	h_mFlow[i_par][i_order][i_cent] = new TH1F(HistName.Data(),HistName.Data(),25,0.0,5.0);
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

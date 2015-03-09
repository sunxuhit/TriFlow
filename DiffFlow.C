#include "TFile.h"
#include "TString.h"
#include "TProfile.h"
#include "TH1F.h"

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

// Calculate integrated v2 and v3
void DiffFlow(Int_t mEnergy = 4, Int_t mMode = 0, Int_t mScreen = 0) // 0: 7.7 GeV, 1: 11.5 GeV, 2: 19.6 GeV, 3: 27 GeV, 4: 39 GeV, 5: 62.4 GeV, 6: 200 GeV | 0: Default, 1: String Melting | 0: 1mb, 1: 3mb, 2: 6mb
{
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
  TProfile *p_mFlow[10][2][4]; // Particle types | 0 for 2nd, 1 for 3rd | 0 for 0-80%, 1 for 0-10%, 2 for 10-40%, 3 for 40-80%
  TH1F *h_mFlow[10][2][4]; // Particle types | 0 for 2nd, 1 for 3rd | 0 for 0-80%, 1 for 0-10%, 2 for 10-40%, 3 for 40-80%

  // pt spectra 
  TH1F *h_mPt[10][2][4]; // Particle types | 0 for 2nd, 1 for 3rd | 0 for 0-80%, 1 for 0-10%, 2 for 10-40%, 3 for 40-80%

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
	for(Int_t i_bin = 1; i_bin < 101; i_bin++)
	{
	  h_mFlow[i_par][i_order][i_cent]->SetBinContent(i_bin,p_mFlow[i_par][i_order][i_cent]->GetBinContent(i_bin));
	  h_mFlow[i_par][i_order][i_cent]->SetBinError(i_bin,p_mFlow[i_par][i_order][i_cent]->GetBinError(i_bin));
	}
	//--------------------------------------------------------------------------------------------------
      }
    }
  }

  // Calculate differential ratio v2/v3 => could be merged into loop above => using individual loop to be more clear
  TH1F *h_ratio[10][4];
  for(Int_t i_par = 0; i_par < 10; i_par++)
  {
    for(Int_t i_cent = 0; i_cent < 4; i_cent++)
    {
      TString HistName = Form("h_ratio_%s_%s",ParType[i_par].Data(),Centrality[i_cent].Data());
      h_ratio[i_par][i_cent] = (TH1F*)h_mFlow[i_par][0][i_cent]->Clone(HistName.Data());
      h_ratio[i_par][i_cent]->Divide(h_mFlow[i_par][1][i_cent]);
    }
  }

  // write Histogram into Output file
  TString outputfile;
  if(mMode == 0)
  {
    outputfile = Form("/home/xusun/Data/AMPT_%s/InteFlow/%s_%s/DiffFlow.root",Mode[mMode].Data(),Energy[mEnergy].Data(),Mode[mMode].Data());
  }
  if(mMode == 1)
  {
    outputfile = Form("/home/xusun/Data/AMPT_%s/InteFlow/%s_%s/%s/DiffFlow.root",Mode[mMode].Data(),Energy[mEnergy].Data(),Mode[mMode].Data(),ScreenMass[mScreen].Data());
  }
  cout << "OutPut File: " << outputfile.Data() << endl;
  TFile *File_output = new TFile(outputfile.Data(),"RECREATE");
  File_output->cd();
  for(Int_t i_par = 0; i_par < 10; i_par++)
  {
    for(Int_t i_cent = 0; i_cent < 4; i_cent++)
    {
      h_ratio[i_par][i_cent]->Write();
    }
  }

  File_output->Close();
}

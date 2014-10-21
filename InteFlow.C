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
static TString Energy[7] = {"7GeV","11GeV","19GeV","27GeV","39GeV","62GeV","200GeV"};
static TString Order[2] = {"2nd","3rd"};
static TString Centrality[4] = {"0080","0010","1040","4080"};
static TString ParType[10] = {"pi_plus","pi_minus","K_plus","K_minus","p","pbar","Lambda","Lambdabar","K0s","phi"};

// Calculate integrated v2 and v3
void InteFlow(Int_t mEnergy = 4, Int_t mMode = 0) // 0: 7.7 GeV, 1: 11.5 GeV, 2: 19.6 GeV, 3: 27 GeV, 4: 39 GeV, 5: 62.4 GeV, 6: 200 GeV | 0: Default, 1: String Melting
{
  //---------------------Constant------------------------
  const Float_t pt_start = 0.2; // start point for pt integration
  const Float_t pt_stop  = 3.0; // stop point for pt integration
  //---------------------Constant------------------------

//  TString inputfile = Form("/home/xusun/Data/AMPT_%s/Flow/%s_%s/Flow_%s.root",Mode[mMode].Data(),Energy[mEnergy].Data(),Mode[mMode].Data(),Energy[mEnergy].Data());
  TString inputfile = Form("/home/xusun/Data/AMPT_%s/Flow/%s_%s/Flow_%s_1_100.root",Mode[mMode].Data(),Energy[mEnergy].Data(),Mode[mMode].Data(),Energy[mEnergy].Data()); // temperory file
  cout << "Input File: " << inputfile.Data() << endl;
  TFile *File_input = TFile::Open(inputfile.Data());


  // flow for pi, K, p, Lambda, K0s by using eta_sub event plane method
  // pi_plus,pi_minus,K_plus,K_minus,p,pbar,Lambda,Lambdabar,K0s,phi
  TProfile *p_mFlow[10][2][4]; // Particle types | 0 for 2nd, 1 for 3rd | 0 for 0-80%, 1 for 0-10%, 2 for 10-40%, 3 for 40-80%
  TH1F *h_mFlow[10][2][4]; // Particle types | 0 for 2nd, 1 for 3rd | 0 for 0-80%, 1 for 0-10%, 2 for 10-40%, 3 for 40-80%

  // pt spectra 
  TH1F *h_mPt[10][2][4]; // Particle types | 0 for 2nd, 1 for 3rd | 0 for 0-80%, 1 for 0-10%, 2 for 10-40%, 3 for 40-80%

  // Integrated v2 and v3
  TH1F *h_mInteFlow[2][4]; // 0 for 2nd, 1 for 3rd | 0 for 0-80%, 1 for 0-10%, 2 for 10-40%, 3 for 40-80%


  // Read in flow TProfile => transfer to TH1F | Read in pt Spectra
  for(Int_t i_par = 0; i_par < 10; i_par++)
  {
    for(Int_t i_order = 0; i_order < 2; i_order++)
    {
      for(Int_t i_cent = 0; i_cent < 4; i_cent++)
      {
	TString ProName;
	TString HistName;

	// v2 relative to event plane
	//--------------------------------------------------------------------------------------------------
	ProName = Form("Flow_%s_%s_%s",ParType[i_par].Data(),Order[i_order].Data(),Centrality[i_cent].Data()); // pi_plus
	p_mFlow[i_par][i_order][i_cent] = (TProfile*)File_input->Get(ProName.Data());
	HistName = Form("h_Flow_%s_%s_%s",ParType[i_par].Data(),Order[i_order].Data(),Centrality[i_cent].Data());
	h_mFlow[i_par][i_order][i_cent] = new TH1F(HistName.Data(),HistName.Data(),100,0.0,10.0);
	for(Int_t i_bin = 1; i_bin < 101; i_bin++)
	{
	  h_mFlow[i_par][i_order][i_cent]->SetBinContent(i_bin,p_mFlow[i_par][i_order][i_cent]->GetBinContent(i_bin));
	  h_mFlow[i_par][i_order][i_cent]->SetBinError(i_bin,p_mFlow[i_par][i_order][i_cent]->GetBinError(i_bin));
	}
	HistName = Form("Pt_%s_%s_%s",ParType[i_par].Data(),Order[i_order].Data(),Centrality[i_cent].Data());
	h_mPt[i_par][i_order][i_cent] = (TH1F*)File_input->Get(HistName.Data());
	h_mPt[i_par][i_order][i_cent]->Rebin(5);
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

  // write Histogram into Output file
  TString outputfile = Form("/home/xusun/Data/AMPT_%s/InteFlow/%s_%s/InteFlow.root",Mode[mMode].Data(),Energy[mEnergy].Data(),Mode[mMode].Data(),Energy[mEnergy].Data());
  TFile *File_output = new TFile(outputfile.Data(),"RECREATE");
  File_output->cd();
  for(Int_t i_order = 0; i_order < 2; i_order++)
  {
    for(Int_t i_cent = 0; i_cent < 4; i_cent++)
    {
      h_mInteFlow[i_order][i_cent]->Write();
    }
  }

  File_output->Close();
}

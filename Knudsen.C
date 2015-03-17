#include "TString.h"
#include "TFile.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TMath.h"
#include "TF1.h"

static TString Mode[2] = {"Default","StringMelting"};
static TString ScreenMass[3] = {"1mb","3mb","6mb"};
static TString Energy[7] = {"7GeV","11GeV","19GeV","27GeV","39GeV","62GeV","200GeV"};
static Double_t mu[3] = {3.2,3.2264,2.2814}; // 1.5mb, 3mb, 6mb

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


Double_t sigma(Double_t *x_val, Double_t *par)
{
  Double_t mu = x_val[0];
  Double_t s = par[0];
  Double_t alpha_s = par[1];

  Double_t z = mu*mu/s;
  Double_t sigma_0 = 9*TMath::Pi()*alpha_s*alpha_s/(2*mu*mu);
  Double_t sigma_t = sigma_0*4*z*(1+z)*((2*z+1)*TMath::Log(1+1/z)-2);
  Double_t sigma_Cs = 1.5*sigma_t/TMath::Sqrt(3);
  
  return sigma_Cs;
}

void Knudsen(Int_t mEnergy = 6, Int_t mMode = 1, Int_t mScreen = 0) // 0: 7.7 GeV, 1: 11.5 GeV, 2: 19.6 GeV, 3: 27 GeV, 4: 39 GeV, 5: 62.4 GeV, 6: 200 GeV | 0: Default, 1: String Melting | 0: 1mb, 1: 3mb, 2: 6mb
{
  TString inputfile;
  if(mMode == 0)
  {
    inputfile  = Form("/home/xusun/Data/AMPT_%s/Epsilon/%s_%s/Epsilon_%s.root",Mode[mMode].Data(),Energy[mEnergy].Data(),Mode[mMode].Data(),Energy[mEnergy].Data());
  }
  if(mMode == 1)
  {
    inputfile = Form("/home/xusun/Data/AMPT_%s/Epsilon/%s_%s/%s/Epsilon_%s.root",Mode[mMode].Data(),Energy[mEnergy].Data(),Mode[mMode].Data(),ScreenMass[mScreen].Data(),Energy[mEnergy].Data());
  }
  cout << "Input File: " << inputfile.Data() << endl;

  TFile *File_input = TFile::Open(inputfile.Data());

  // dNdy calculation
  TH1F *h_mEventCounter4 = (TH1F*)File_input->Get("h_mEventCounter4");
  TH1F *h_mRapWide[4];
  Double_t dNdy[4] = {0.0,0.0,0.0,0.0};
  Double_t err_dNdy[4] = {0.0,0.0,0.0,0.0};
  Double_t y_start = -0.5;
  Double_t y_stop  = 0.5;
  for(Int_t i_cent = 0; i_cent < 4; i_cent++)
  {
    Float_t Event_Counter = h_mEventCounter4->GetBinContent(i_cent+1);
    TString HistName = Form("h_mRapWide_%d",i_cent);
    h_mRapWide[i_cent] = (TH1F*)File_input->Get(HistName.Data());
    h_mRapWide[i_cent]->Scale(1/Event_Counter);
    Int_t bin_start = h_mRapWide[i_cent]->FindBin(y_start);
    Int_t bin_stop  = h_mRapWide[i_cent]->FindBin(y_stop);
    dNdy[i_cent] = h_mRapWide[i_cent]->IntegralAndError(bin_start,bin_stop,err_dNdy[i_cent],"width");
    cout << "dNdy = " << dNdy[i_cent] << " +/- " << err_dNdy[i_cent] << endl;
  }

  // AreaT calculation
  TProfile *p_mAreaT4;
  p_mAreaT4 = (TProfile*)File_input->Get("p_mAreaT4");
  Double_t AreaT[4] = {0.0,0.0,0.0,0.0};
  Double_t err_AreaT[4] = {0.0,0.0,0.0,0.0};
  for(Int_t i_cent = 0; i_cent < 4; i_cent++)
  {
    AreaT[i_cent] = (Double_t)p_mAreaT4->GetBinContent(i_cent+1);
    err_AreaT[i_cent] = (Double_t)p_mAreaT4->GetBinError(i_cent+1);
    cout << "AreaT = " << AreaT[i_cent] << " +/- " << err_AreaT[i_cent] << endl;
  }

  // sigma calculation
  TF1 *f_sigma = new TF1("f_sigma",sigma,0.0,5.0,2);
  f_sigma->FixParameter(0,0.35); // s
  f_sigma->FixParameter(1,0.4714); // alpha_s
  Double_t Sigma_Cs = f_sigma->Eval(mu[mScreen]);
  cout << "sigma*Cs = " << Sigma_Cs << endl;

  TH1F *h_Knudsen = new TH1F("h_Knudsen","h_Knudsen",4,-0.5,3.5);
  for(Int_t i_cent = 0; i_cent < 4; i_cent++)
  {
    Double_t Knudsen = Sigma_Cs*dNdy[i_cent]/AreaT[i_cent];
    Double_t err_Knudsen = ErrDiv(dNdy[i_cent],AreaT[i_cent],err_dNdy[i_cent],err_AreaT[i_cent]);
    h_Knudsen->SetBinContent(i_cent+1,Knudsen);
    h_Knudsen->SetBinError(i_cent+1,err_Knudsen);
  }
  h_Knudsen->Draw();

  /*
  // write Histogram into Output file
  TString outputfile;
  if(mMode == 0)
  {
    outputfile = Form("/home/xusun/Data/AMPT_%s/InteFlow/%s_%s/Knudsen.root",Mode[mMode].Data(),Energy[mEnergy].Data(),Mode[mMode].Data());
  }
  if(mMode == 1)
  {
    outputfile = Form("/home/xusun/Data/AMPT_%s/InteFlow/%s_%s/%s/Knudsen.root",Mode[mMode].Data(),Energy[mEnergy].Data(),Mode[mMode].Data(),ScreenMass[mScreen].Data());
  }
  cout << "OutPut File: " << outputfile.Data() << endl;
  TFile *File_output = new TFile(outputfile.Data(),"RECREATE");
  File_output->cd();
  File_output->Close();
  */
}

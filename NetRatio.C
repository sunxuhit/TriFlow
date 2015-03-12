#include "TString.h"
#include "TFile.h"
#include "TH1F.h"

static TString Mode[2] = {"Default","StringMelting"};
static TString ScreenMass[3] = {"1mb","3mb","6mb"};
static TString Energy[7] = {"7GeV","11GeV","19GeV","27GeV","39GeV","62GeV","200GeV"};
static TString Order[2] = {"2nd","3rd"};
static TString Centrality[4] = {"0080","0010","1040","4080"};
static TString ParType[10] = {"pi_plus","pi_minus","K_plus","K_minus","p","pbar","Lambda","Lambdabar","K0s","phi"};
static Float_t ParMass[10] = {0.1396,0.1396,0.4937,0.4937,0.9382,0.9382,1.116,1.116,0.4976,1.019};

void NetRatio(Int_t mEnergy = 6, Int_t mMode = 0, Int_t mScreen = 0) // 0: 7.7 GeV, 1: 11.5 GeV, 2: 19.6 GeV, 3: 27 GeV, 4: 39 GeV, 5: 62.4 GeV, 6: 200 GeV | 0: Default, 1: String Melting | 0: 1mb, 1: 3mb, 2: 6mb
{
  TString inputfile_flow, inputfile_eps;
  if(mMode == 0)
  {
    inputfile_flow = Form("/home/xusun/Data/AMPT_%s/InteFlow/%s_%s/InteFlow.root",Mode[mMode].Data(),Energy[mEnergy].Data(),Mode[mMode].Data());
    inputfile_eps  = Form("/home/xusun/Data/AMPT_%s/Epsilon/%s_%s/Epsilon_%s.root",Mode[mMode].Data(),Energy[mEnergy].Data(),Mode[mMode].Data(),Energy[mEnergy].Data());
  }
  if(mMode == 1)
  {
    inputfile_flow = Form("/home/xusun/Data/AMPT_%s/InteFlow/%s_%s/%s/InteFlow.root",Mode[mMode].Data(),Energy[mEnergy].Data(),Mode[mMode].Data(),ScreenMass[mScreen].Data());
    inputfile_eps  = Form("/home/xusun/Data/AMPT_%s/Epsilon/%s_%s/%s/Epsilon_%s.root",Mode[mMode].Data(),Energy[mEnergy].Data(),Mode[mMode].Data(),ScreenMass[mScreen].Data(),Energy[mEnergy].Data());
  }
  cout << "Input File flow: " << inputfile_flow.Data() << endl;
  cout << "Input File epsilon: " << inputfile_eps.Data() << endl;

  TFile *File_input_flow = TFile::Open(inputfile_flow.Data());
  TFile *File_input_eps  = TFile::Open(inputfile_eps.Data());

  TH1F *h_ratio[4]; // 0: 0-80%, 1: 0-10%, 2: 10-40%, 3: 40-80%
  TH1F *h_Inte_v2[4]; // 0: 0-80%, 1: 0-10%, 2: 10-40%, 3: 40-80%
  TH1F *h_Eps2;
  TH1F *h_Inte_v3[4]; // 0: 0-80%, 1: 0-10%, 2: 10-40%, 3: 40-80%
  TH1F *h_Eps3;

  for(Int_t i_cent = 0; i_cent < 4; i_cent++)
  {
    TString HistName_ratio = Form("h_ratio_%s",Centrality[i_cent].Data());
    h_ratio[i_cent] = (TH1F*)File_input_flow->Get(HistName_ratio.Data());

    TString HistName_v2 = Form("h_mInteFlow_2nd_%s",Centrality[i_cent].Data());
    h_Inte_v2[i_cent] = (TH1F*)File_input_flow->Get(HistName_v2.Data());

    TString HistName_v3 = Form("h_mInteFlow_3rd_%s",Centrality[i_cent].Data());
    h_Inte_v3[i_cent] = (TH1F*)File_input_flow->Get(HistName_v3.Data());
  }
  h_Eps2 = (TH1F*)File_input_eps->Get("p_mEpsilon4_2nd");
  h_Eps3 = (TH1F*)File_input_eps->Get("p_mEpsilon4_3rd");

  for(Int_t i_cent = 0; i_cent < 4; i_cent++) // calculate flow/epsilon
  {
    Float_t epsilon_2 = h_Eps2->GetBinContent(i_cent+1);
    h_Inte_v2[i_cent]->Scale(1/epsilon_2);
    TString HistName_v2 = Form("h_v2_e2_%s",Centrality[i_cent].Data());
    h_Inte_v2[i_cent]->SetName(HistName_v2.Data());

    Float_t epsilon_3 = h_Eps3->GetBinContent(i_cent+1);
    h_Inte_v3[i_cent]->Scale(1/epsilon_3);
    TString HistName_v3 = Form("h_v3_e3_%s",Centrality[i_cent].Data());
    h_Inte_v3[i_cent]->SetName(HistName_v3.Data());

    Float_t ratio = epsilon_2/epsilon_3;
    h_ratio[i_cent]->Scale(1/ratio);
    TString HistName_ratio = Form("h_mNetRatio_%s",Centrality[i_cent].Data());
    h_ratio[i_cent]->SetName(HistName_ratio.Data());
  }

  // write Histogram into Output file
  TString outputfile;
  if(mMode == 0)
  {
    outputfile = Form("/home/xusun/Data/AMPT_%s/InteFlow/%s_%s/NetRatio.root",Mode[mMode].Data(),Energy[mEnergy].Data(),Mode[mMode].Data());
  }
  if(mMode == 1)
  {
    outputfile = Form("/home/xusun/Data/AMPT_%s/InteFlow/%s_%s/%s/NetRatio.root",Mode[mMode].Data(),Energy[mEnergy].Data(),Mode[mMode].Data(),ScreenMass[mScreen].Data());
  }
  cout << "OutPut File: " << outputfile.Data() << endl;
  TFile *File_output = new TFile(outputfile.Data(),"RECREATE");
  File_output->cd();

  for(Int_t i_cent = 0; i_cent < 4; i_cent++)
  {
    h_ratio[i_cent]->Write();
    h_Inte_v2[i_cent]->Write();
    h_Inte_v3[i_cent]->Write();
  }

  File_output->Close();

}

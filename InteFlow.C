#include "TFile.h"
#include "TString.h"
#include "TProfile.h"
#include "TH1F.h"

static TString Mode[2] = {"Default","StringMelting"};
static TString Energy[7] = {"7GeV","11GeV","19GeV","27GeV","39GeV","62GeV","200GeV"};
static TString Order[2] = {"2nd","3rd"};
static TString Centrality[4] = {"0080","0010","1040","4080"};

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
  TProfile *p_mFlow_pi_plus[2][4]; // 0 for 2nd, 1 for 3rd | 0 for 0-80%, 1 for 0-10%, 2 for 10-40%, 3 for 40-80%
  TProfile *p_mFlow_pi_minus[2][4];
  TProfile *p_mFlow_K_plus[2][4];
  TProfile *p_mFlow_K_minus[2][4];
  TProfile *p_mFlow_p[2][4];
  TProfile *p_mFlow_pbar[2][4];
  TProfile *p_mFlow_Lambda[2][4];
  TProfile *p_mFlow_Lambdabar[2][4];
  TProfile *p_mFlow_K0s[2][4];
  TProfile *p_mFlow_phi[2][4];

  // pt spectra for pi, K, p, Lambda, K0s
  TH1F *h_mPt_pi_plus[2][4]; // 0 for 2nd, 1 for 3rd | 0 for 0-80%, 1 for 0-10%, 2 for 10-40%, 3 for 40-80%
  TH1F *h_mPt_pi_minus[2][4];
  TH1F *h_mPt_K_plus[2][4];
  TH1F *h_mPt_K_minus[2][4];
  TH1F *h_mPt_p[2][4];
  TH1F *h_mPt_pbar[2][4];
  TH1F *h_mPt_Lambda[2][4];
  TH1F *h_mPt_Lambdabar[2][4];
  TH1F *h_mPt_K0s[2][4];
  TH1F *h_mPt_phi[2][4];

  // Read flow TProfile
  for(Int_t i_order = 0; i_order < 2; i_order++)
  {
    for(Int_t i_cent = 0; i_cent < 4; i_cent++)
    {
      TString ProName;
      TString HistName;

      // v2 relative to event plane
      ProName = Form("Flow_pi_plus_%s_%s",Order[i_order].Data(),Centrality[i_cent].Data()); // pi_plus
      p_mFlow_pi_plus[i_order][i_cent] = (TProfile*)File_input->Get(ProName.Data());
      HistName = Form("Pt_pi_plus_%s_%s",Order[i_order].Data(),Centrality[i_cent].Data());
      h_mPt_pi_plus[i_order][i_cent] = (TH1F*)File_input->Get(HistName.Data());

      ProName = Form("Flow_pi_minus_%s_%s",Order[i_order].Data(),Centrality[i_cent].Data()); // pi_minus
      p_mFlow_pi_minus[i_order][i_cent] = (TProfile*)File_input->Get(ProName.Data());
      HistName = Form("Pt_pi_minus_%s_%s",Order[i_order].Data(),Centrality[i_cent].Data());
      h_mPt_pi_minus[i_order][i_cent] = (TH1F*)File_input->Get(HistName.Data());

      ProName = Form("Flow_K_plus_%s_%s",Order[i_order].Data(),Centrality[i_cent].Data()); // K_plus
      p_mFlow_K_plus[i_order][i_cent] = (TProfile*)File_input->Get(ProName.Data());
      HistName = Form("Pt_K_plus_%s_%s",Order[i_order].Data(),Centrality[i_cent].Data());
      h_mPt_K_plus[i_order][i_cent] = (TH1F*)File_input->Get(HistName.Data());

      ProName = Form("Flow_K_minus_%s_%s",Order[i_order].Data(),Centrality[i_cent].Data()); // K_minus
      p_mFlow_K_minus[i_order][i_cent] = (TProfile*)File_input->Get(ProName.Data());
      HistName = Form("Pt_K_minus_%s_%s",Order[i_order].Data(),Centrality[i_cent].Data());
      h_mPt_K_minus[i_order][i_cent] = (TH1F*)File_input->Get(HistName.Data());

      ProName = Form("Flow_p_%s_%s",Order[i_order].Data(),Centrality[i_cent].Data()); // p
      p_mFlow_p[i_order][i_cent] = (TProfile*)File_input->Get(ProName.Data());
      HistName = Form("Pt_p_%s_%s",Order[i_order].Data(),Centrality[i_cent].Data());
      h_mPt_p[i_order][i_cent] = (TH1F*)File_input->Get(HistName.Data());

      ProName = Form("Flow_pbar_%s_%s",Order[i_order].Data(),Centrality[i_cent].Data()); // pbar
      p_mFlow_pbar[i_order][i_cent] = (TProfile*)File_input->Get(ProName.Data());
      HistName = Form("Pt_pbar_%s_%s",Order[i_order].Data(),Centrality[i_cent].Data());
      h_mPt_pbar[i_order][i_cent] = (TH1F*)File_input->Get(HistName.Data());

      ProName = Form("Flow_Lambda_%s_%s",Order[i_order].Data(),Centrality[i_cent].Data()); // Lambda
      p_mFlow_Lambda[i_order][i_cent] = (TProfile*)File_input->Get(ProName.Data());
      HistName = Form("Pt_Lambda_%s_%s",Order[i_order].Data(),Centrality[i_cent].Data());
      h_mPt_Lambda[i_order][i_cent] = (TH1F*)File_input->Get(HistName.Data());

      ProName = Form("Flow_Lambdabar_%s_%s",Order[i_order].Data(),Centrality[i_cent].Data()); // Lambdabar
      p_mFlow_Lambdabar[i_order][i_cent] = (TProfile*)File_input->Get(ProName.Data());
      HistName = Form("Pt_Lambdabar_%s_%s",Order[i_order].Data(),Centrality[i_cent].Data());
      h_mPt_Lambdabar[i_order][i_cent] = (TH1F*)File_input->Get(HistName.Data());

      ProName = Form("Flow_K0s_%s_%s",Order[i_order].Data(),Centrality[i_cent].Data()); // K0s
      p_mFlow_K0s[i_order][i_cent] = (TProfile*)File_input->Get(ProName.Data());
      HistName = Form("Pt_K0s_%s_%s",Order[i_order].Data(),Centrality[i_cent].Data());
      h_mPt_K0s[i_order][i_cent] = (TH1F*)File_input->Get(HistName.Data());

      ProName = Form("Flow_phi_%s_%s",Order[i_order].Data(),Centrality[i_cent].Data()); // phi
      p_mFlow_phi[i_order][i_cent] = (TProfile*)File_input->Get(ProName.Data());
      HistName = Form("Pt_phi_%s_%s",Order[i_order].Data(),Centrality[i_cent].Data());
      h_mPt_phi[i_order][i_cent] = (TH1F*)File_input->Get(HistName.Data());
    }
  }

  TH1F *h_test = new TH1F("h_test","h_test",100,0.0,10.0);
  for(Int_t i_bin = 1; i_bin < 101; i_bin++)
  {
    h_test->SetBinContent(i_bin,p_mFlow_pi_plus[0][1]->GetBinContent(i_bin));
    h_test->SetBinError(i_bin,p_mFlow_pi_plus[0][1]->GetBinError(i_bin));
  }

  TH1F *h_mult = h_test->Clone("h_mult");;
  h_mPt_pi_plus[0][1]->Rebin(5);
  h_mult->Multiply(h_mPt_pi_plus[0][1]);
  h_mult->SetLineColor(2);
  h_mult->Draw("pE");
  h_mPt_pi_plus[0][1]->Draw("pE same");

  h_test->SetMarkerColor(2);
  h_test->SetMarkerStyle(24);
  h_test->SetMarkerSize(1.0);
  h_test->Draw("pE same");
}

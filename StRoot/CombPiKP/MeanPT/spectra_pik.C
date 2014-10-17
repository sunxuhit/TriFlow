#include "TFile.h"
#include "TH1F.h"
#include "TString.h"
#include "TF1.h"
#include "TMath.h"
#include "TCanvas.h"
#include "/u/xusun/STAR/Analysis_v3/StRoot/CombPiKP/draw.h"
#include "TF1.h"
#include "Math/MinimizerOptions.h"


//----------------------------------------------------------------------------------------
Double_t PtFitFunc2_mod_x(Double_t* x_val, Double_t* par)
{
  Double_t x, y, m0, Temp, Ampl, shift;
  m0    = par[0];
  Temp  = par[1];
  Ampl  = par[2];
  shift = par[3];
  Double_t pol0 = par[4];
  Double_t pol1 = par[5];
  Double_t pol2 = par[6];

  x = x_val[0];
  Double_t poly = pol0 + pol1*x + pol2*x*x;
  y = x*(Ampl*(x-shift)*sqrt((x-shift)*(x-shift)+m0*m0)*TMath::Exp(-(sqrt((x-shift)*(x-shift)+m0*m0)-m0)/Temp)) + poly;

  return y;
}

Double_t MeanPt(Double_t* x_val, Double_t* par)
{
  Double_t x, y, m0, Temp, Ampl, shift;
  m0    = par[0];
  Temp  = par[1];
  Ampl  = par[2];
  shift = par[3];
  Double_t pol0 = par[4];
  Double_t pol1 = par[5];
  Double_t pol2 = par[6];

  x = x_val[0];
  Double_t poly = pol0 + pol1*x + pol2*x*x;
  Double_t dNdpT = x*(Ampl*(x-shift)*sqrt((x-shift)*(x-shift)+m0*m0)*TMath::Exp(-(sqrt((x-shift)*(x-shift)+m0*m0)-m0)/Temp)) + poly;

  y = x*dNdpT;

  return y;
}

// mEnergy: 0 for 200GeV, 1 for 39GeV
// mCentrality: 0 for 00-80, 1 for 00-10, 2 for 10-40, 3 for 40-80
void spectra_pik(Int_t mEnergy = 0, Int_t mCentrality = 0)
{
  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(50000);

  const Int_t pt_total = 16;
  const Int_t cent_total = 4; // 4
  const Int_t cent_start = 0;
  const Int_t cent_stop  = 1;
  const Int_t charge_total = 2; // 2
  const Int_t charge_start = 0;
  const Int_t charge_stop  = 2;
  const Int_t eta_total = 4; // 4
  const Int_t eta_start = 0;
  const Int_t eta_stop = 1;
  const Int_t energy_total = 2; // 2

  TString Energy[2] = {"200GeV","39GeV"};
  TString Centrality[4] = {"0080","0010","1040","4080"};

  // pt bin
  //                            0 , 1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 ,10 ,11 ,12 ,13 ,14 ,15
  Float_t pt_low[pt_total] = {0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4};
  Float_t pt_up[pt_total]  = {0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,4.2};

  TFile *File_input[2];
  TH1F *h_pion[2];
  TH1F *h_kaon[2];

  for(Int_t i_char = 0; i_char < 2; i_char++)
  {
    TString inputfile = Form("/project/projectdirs/star/xusun/OutPut/AuAu%s/Mass2_nSigmaPion/merged_file/pt_spectra/Counts_%s_charge_%d_etagap_0.root",Energy[mEnergy].Data(),Centrality[mCentrality].Data(),i_char);
    File_input[i_char] = TFile::Open(inputfile.Data());
    TString HistName;
    HistName = Form("Spec_Centrality_0_charge_%d_etagap_0_counts_pion",i_char);
    h_pion[i_char] = (TH1F*)File_input[i_char]->Get(HistName.Data());
    HistName = Form("Spec_Centrality_0_charge_%d_etagap_0_counts_kaon",i_char);
    h_kaon[i_char] = (TH1F*)File_input[i_char]->Get(HistName.Data());
  }

  // Fit the distribution with Boltzmann Distritbution
  TF1 *f_Yield_pion[charge_total];
  TCanvas *c_Yields_pion[charge_total];
  for(Int_t charge = charge_start; charge < charge_stop; charge++) // charge bin
  {
    TString CanName = Form("c_Yileds_Charge_%d_pion",charge);
    c_Yields_pion[charge] = new TCanvas(CanName.Data(),CanName.Data(),1400,10,800,800);
    c_Yields_pion[charge]->cd();
    c_Yields_pion[charge]->cd()->SetLeftMargin(0.15);
    c_Yields_pion[charge]->cd()->SetBottomMargin(0.15);
    c_Yields_pion[charge]->cd()->SetGrid(0,0);
    c_Yields_pion[charge]->cd()->SetTicks(1,1);
    c_Yields_pion[charge]->cd()->SetLogy();
    h_pion[charge]->Draw("pE");
    TString FuncName = Form("f_Poly_Charge_%d_pion",charge);
    f_Yield_pion[charge] = new TF1(FuncName.Data(),PtFitFunc2_mod_x,0.2,3.4,7);
    f_Yield_pion[charge]->SetParameter(0,0.138);
    f_Yield_pion[charge]->SetParameter(1,0.160);
    f_Yield_pion[charge]->SetParLimits(1,0,10000);
    f_Yield_pion[charge]->SetParameter(2,2.87956e+13);
    f_Yield_pion[charge]->SetParameter(3,0.001);
    f_Yield_pion[charge]->SetParameter(4, 3.05968e+09);
    f_Yield_pion[charge]->SetParameter(5,-1.99156e+09);
    f_Yield_pion[charge]->SetParameter(6,0.0);
    f_Yield_pion[charge]->SetRange(0.2,3.4);
    h_pion[charge]->Fit(f_Yield_pion[charge],"NRI");
    f_Yield_pion[charge]->Draw("l same");
  }

  TF1 *f_MeanPT_pion[charge_total];
  for(Int_t charge = charge_start; charge < charge_stop; charge++) // charge bin
  {
    TString FuncName = Form("f_MeanPT_Charge_%d_pion",charge);
    f_MeanPT_pion[charge] = new TF1(FuncName.Data(),MeanPt,0.2,3.4,7);
    for(Int_t i_par = 0; i_par < 7; i_par++)
    {
      f_MeanPT_pion[charge]->SetParameter(i_par,f_Yield_pion[charge]->GetParameter(i_par));
    }
    f_MeanPT_pion[charge]->SetRange(0.2,3.4);
  }

  TF1 *f_Yield_kaon[charge_total];
  TCanvas *c_Yields_kaon[charge_total];
  for(Int_t charge = charge_start; charge < charge_stop; charge++) // charge bin
  {
    TString CanName = Form("c_Yileds_Charge_%d_kaon",charge);
    c_Yields_kaon[charge] = new TCanvas(CanName.Data(),CanName.Data(),1400,10,800,800);
    c_Yields_kaon[charge]->cd();
    c_Yields_kaon[charge]->cd()->SetLeftMargin(0.15);
    c_Yields_kaon[charge]->cd()->SetBottomMargin(0.15);
    c_Yields_kaon[charge]->cd()->SetGrid(0,0);
    c_Yields_kaon[charge]->cd()->SetTicks(1,1);
    c_Yields_kaon[charge]->cd()->SetLogy();
    h_kaon[charge]->Draw("pE");
    TString FuncName = Form("f_Poly_Charge_%d_kaon",charge);
    f_Yield_kaon[charge] = new TF1(FuncName.Data(),PtFitFunc2_mod_x,0.2,3.4,7);
    f_Yield_kaon[charge]->SetParameter(0,0.494);
    f_Yield_kaon[charge]->SetParameter(1,0.250);
    f_Yield_kaon[charge]->SetParLimits(1,0,10000);
    f_Yield_kaon[charge]->SetParameter(2,6.17689e+10);
    f_Yield_kaon[charge]->SetParameter(3,0.001);
    f_Yield_kaon[charge]->SetParameter(4,-7.65846e+08);
    f_Yield_kaon[charge]->SetParameter(5, 5.30772e+08);
    f_Yield_kaon[charge]->SetParameter(6,-5.28569e+07);
    f_Yield_kaon[charge]->SetRange(0.2,3.4);
    h_kaon[charge]->Fit(f_Yield_kaon[charge],"NRI");
    f_Yield_kaon[charge]->Draw("l same");
  }

  TF1 *f_MeanPT_kaon[charge_total];
  for(Int_t charge = charge_start; charge < charge_stop; charge++) // charge bin
  {
    TString FuncName = Form("f_MeanPT_Charge_%d_kaon",charge);
    f_MeanPT_kaon[charge] = new TF1(FuncName.Data(),MeanPt,0.2,3.4,7);
    for(Int_t i_par = 0; i_par < 7; i_par++)
    {
      f_MeanPT_kaon[charge]->SetParameter(i_par,f_Yield_kaon[charge]->GetParameter(i_par));
    }
    f_MeanPT_kaon[charge]->SetRange(0.2,3.4);
  }

  Float_t meanPT_pion[charge_total][pt_total-1];
  Float_t meanPT_kaon[charge_total][pt_total-1];
  for(Int_t charge = charge_start; charge < charge_stop; charge++) // charge bin
  {
    for(Int_t i_pt = 0; i_pt < pt_total-1; i_pt++)
    {
      meanPT_pion[charge][i_pt] = f_MeanPT_pion[charge]->Integral(pt_low[i_pt],pt_up[i_pt])/f_Yield_pion[charge]->Integral(pt_low[i_pt],pt_up[i_pt]);
      meanPT_kaon[charge][i_pt] = f_MeanPT_kaon[charge]->Integral(pt_low[i_pt],pt_up[i_pt])/f_Yield_kaon[charge]->Integral(pt_low[i_pt],pt_up[i_pt]);
//      if(charge == 0) cout << "pt_bin_low = " << pt_low[i_pt] << ", pt_bin_up = " << pt_up[i_pt] << ", Mean pT Kaon plus = "  << meanPT_kaon[charge][i_pt] << endl;
//      if(charge == 1) cout << "pt_bin_low = " << pt_low[i_pt] << ", pt_bin_up = " << pt_up[i_pt] << ", Mean pT Kaon minus = " << meanPT_kaon[charge][i_pt] << endl;
    }
  }

  // pion
  for(Int_t charge = charge_start; charge < charge_stop; charge++) // charge bin
  {
    for(Int_t i_pt = 0; i_pt < pt_total-1; i_pt++)
    {
      if(charge == 0 && i_pt == 0) cout << "MeanPT_pionPlus[14] = {";
      if(charge == 0 && i_pt < pt_total-3) cout << meanPT_pion[charge][i_pt] << ", ";
      if(charge == 0 && i_pt == pt_total-3) cout << meanPT_pion[charge][i_pt] << "} ";
      if(charge == 0 && i_pt == pt_total-2) cout << endl;

      if(charge == 1 && i_pt == 0) cout << "MeanPT_pionMinus[14] = {";
      if(charge == 1 && i_pt < pt_total-3) cout << meanPT_pion[charge][i_pt] << ", ";
      if(charge == 1 && i_pt == pt_total-3) cout << meanPT_pion[charge][i_pt] << "} ";
      if(charge == 1 && i_pt == pt_total-2) cout << endl;
    }
  }

  // kaon
  for(Int_t charge = charge_start; charge < charge_stop; charge++) // charge bin
  {
    for(Int_t i_pt = 0; i_pt < pt_total-1; i_pt++)
    {
      if(charge == 0 && i_pt == 0) cout << "MeanPT_KPlus[14] = {";
      if(charge == 0 && i_pt < pt_total-3) cout << meanPT_kaon[charge][i_pt] << ", ";
      if(charge == 0 && i_pt == pt_total-3) cout << meanPT_kaon[charge][i_pt] << "} ";
      if(charge == 0 && i_pt == pt_total-2) cout << endl;

      if(charge == 1 && i_pt == 0) cout << "MeanPT_KMinus[14] = {";
      if(charge == 1 && i_pt < pt_total-3) cout << meanPT_kaon[charge][i_pt] << ", ";
      if(charge == 1 && i_pt == pt_total-3) cout << meanPT_kaon[charge][i_pt] << "} ";
      if(charge == 1 && i_pt == pt_total-2) cout << endl;
    }
  }
}

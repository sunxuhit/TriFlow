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
Double_t Gaussian(Double_t *x_val, Double_t *par)
{
  Double_t x, y, mu, sigma, Norm;
  x = x_val[0];
  mu = par[0];
  sigma = par[1];
  Norm = par[2];

  y = Norm*TMath::Exp(-0.5*(x-mu)*(x-mu)/(sigma*sigma))/(sigma*TMath::Sqrt(2.0*TMath::Pi()));

  return y;
}

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
//----------------------------------------------------------------------------------------

// mEnergy: 0 for 200GeV, 1 for 39GeV
// mCentrality: 0 for 00-80, 1 for 00-10, 2 for 10-40, 3 for 40-80
void spectra_proton(Int_t mEnergy = 0, Int_t mCentrality = 0)
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
  TString Order[2] = {"2nd","3rd"};

  Double_t counts_proton[pt_total][cent_total][charge_total][eta_total];
  Double_t errors_proton[pt_total][cent_total][charge_total][eta_total];

  // pt bin
  Float_t pt_low[pt_total] = {0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4};
  Float_t pt_up[pt_total]  = {0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,4.2};


  Float_t bin_width[pt_total];
  for(Int_t i = 0; i < pt_total; i++)
  {
    bin_width[i] = pt_up[i] - pt_low[i];
  }

  Float_t pt_cut_low[pt_total] = {0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.65,0.7,0.7,0.75,0.879,0.879,0.879,0.879,0.879};
  Float_t pt_cut_up[pt_total]  = {1.2,1.2,1.2,1.2,1.2,1.2,1.2,1.20,1.2,1.2,1.20,1.400,1.400,1.400,1.400,1.400};

  TString Energy[2] = {"200GeV","39GeV"};
  TString Centrality[4] = {"0080","0010","1040","4080"};
  TString inputfile = Form("/project/projectdirs/star/xusun/OutPut/AuAu%s/Mass2_Proton/merged_file/merged_file_%s_M2_Proton_%s_etagap_00.root",Energy[mEnergy].Data(),Energy[mEnergy].Data(),Centrality[mCentrality].Data());
  cout << inputfile.Data() << endl;
  TFile *file_input = TFile::Open(inputfile.Data());

  // raw pt spectra
  TH1F *h_pt_spectra_proton[pt_total][cent_total][charge_total][eta_total];
  TF1  *f_pt_spectra_proton[pt_total][cent_total][charge_total][eta_total];
  for(Int_t pt_bin = 0; pt_bin < pt_total; pt_bin++) // pt bin
  {
    for(Int_t cent = cent_start; cent < cent_stop; cent++) // centrality bin
    {
      for(Int_t charge = charge_start; charge < charge_stop; charge++) // charge bin
      {
	for(Int_t eta_bin = eta_start; eta_bin < eta_stop; eta_bin++) // eta gap bin
	{
	  TString SpecName = Form("Spectra_pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_Proton_2nd",pt_bin,cent,charge,eta_bin);
	  h_pt_spectra_proton[pt_bin][cent][charge][eta_bin] = (TH1F*)file_input->Get(SpecName.Data());
	  TString FuncName = "f_" + SpecName;
	  f_pt_spectra_proton[pt_bin][cent][charge][eta_bin] = new TF1(FuncName.Data(),Gaussian,0.6,1.4,3);
	  f_pt_spectra_proton[pt_bin][cent][charge][eta_bin]->SetParameter(0,0.938*0.938);
	  f_pt_spectra_proton[pt_bin][cent][charge][eta_bin]->SetParameter(1,0.1);
	  f_pt_spectra_proton[pt_bin][cent][charge][eta_bin]->SetParameter(2,10000);
	}
      }
    }
  }

  // Draw the cut range
  TCanvas *c_cut[cent_total][charge_total][eta_total];
  for(Int_t cent = cent_start; cent < cent_stop; cent++) // centrality bin
  {
    for(Int_t charge = charge_start; charge < charge_stop; charge++) // charge bin
    {
      for(Int_t eta_bin = eta_start; eta_bin < eta_stop; eta_bin++) // eta gap bin
      {
	TString CanName = Form("c_cut_Centrality_%d_Charge_%d_EtaGap_%d",cent,charge,eta_bin);
	c_cut[cent][charge][eta_bin] = new TCanvas(CanName.Data(),CanName.Data(),1400,10,1200,1200);
	c_cut[cent][charge][eta_bin]->Divide(4,4);
	for(Int_t pt_bin = 0; pt_bin < pt_total; pt_bin++) // pt bin
	{
	  c_cut[cent][charge][eta_bin]->cd(pt_bin+1);
	  c_cut[cent][charge][eta_bin]->cd(pt_bin+1)->SetLeftMargin(0.1);
	  c_cut[cent][charge][eta_bin]->cd(pt_bin+1)->SetBottomMargin(0.1);
	  c_cut[cent][charge][eta_bin]->cd(pt_bin+1)->SetGrid(0,0);
	  c_cut[cent][charge][eta_bin]->cd(pt_bin+1)->SetTicks(1,1);
	  h_pt_spectra_proton[pt_bin][cent][charge][eta_bin]->DrawCopy("pE");
	  if(pt_bin < 8)
	  {
	    PlotLine(0.6,0.6,0.0,h_pt_spectra_proton[pt_bin][cent][charge][eta_bin]->GetMaximum()/2.0,1,2,2);
	    PlotLine(1.2,1.2,0.0,h_pt_spectra_proton[pt_bin][cent][charge][eta_bin]->GetMaximum()/2.0,1,2,2);
	  }
	  if(pt_bin >= 8)
	  {
	    f_pt_spectra_proton[pt_bin][cent][charge][eta_bin]->SetRange(0.75,1.00);
	    h_pt_spectra_proton[pt_bin][cent][charge][eta_bin]->Fit(f_pt_spectra_proton[pt_bin][cent][charge][eta_bin],"RQ");
	    Float_t mean = f_pt_spectra_proton[pt_bin][cent][charge][eta_bin]->GetParameter(0);
	    PlotLine(mean,mean,0.0,h_pt_spectra_proton[pt_bin][cent][charge][eta_bin]->GetMaximum()/2.0,1,2,2);
	  }
	}
      }
    }
  }

  // Get Yields in the cut range
  TH1F *h_Yield[cent_total][charge_total][eta_total];
  for(Int_t cent = cent_start; cent < cent_stop; cent++) // centrality bin
  {
    for(Int_t charge = charge_start; charge < charge_stop; charge++) // charge bin
    {
      for(Int_t eta_bin = eta_start; eta_bin < eta_stop; eta_bin++) // eta gap bin
      {
	TString HistName = Form("h_Yield_Centrality_%d_Charge_%d_EtaGap_%d",cent,charge,eta_bin);
	h_Yield[cent][charge][eta_bin] = new TH1F(HistName.Data(),HistName.Data(),pt_total-1,pt_low);
	for(Int_t pt_bin = 0; pt_bin < pt_total-1; pt_bin++)
	{
	  if(pt_bin < 8)
	  {
	    Int_t bin_start = h_pt_spectra_proton[pt_bin][cent][charge][eta_bin]->FindBin(0.6);
	    Int_t bin_stop  = h_pt_spectra_proton[pt_bin][cent][charge][eta_bin]->FindBin(1.2);
	    Float_t Counts = 0.0;
	    Float_t Errors = 0.0;
	    for(Int_t bin = bin_start; bin <= bin_stop; bin++)
	    {
	      Counts += h_pt_spectra_proton[pt_bin][cent][charge][eta_bin]->GetBinContent(bin)/bin_width[pt_bin];
	      Errors += h_pt_spectra_proton[pt_bin][cent][charge][eta_bin]->GetBinError(bin)*h_pt_spectra_proton[pt_bin][cent][charge][eta_bin]->GetBinError(bin)/bin_width[pt_bin];
	    }
	    Float_t pt_center = (pt_low[pt_bin]+pt_up[pt_bin])/2.0;
	    h_Yield[cent][charge][eta_bin]->SetBinContent(h_Yield[cent][charge][eta_bin]->FindBin(pt_center),Counts);
	    h_Yield[cent][charge][eta_bin]->SetBinError(h_Yield[cent][charge][eta_bin]->FindBin(pt_center),TMath::Sqrt(Errors));
	  }
	  if(pt_bin >= 8)
	  {
	    Float_t mean = f_pt_spectra_proton[pt_bin][cent][charge][eta_bin]->GetParameter(0);
	    Int_t bin_start = h_pt_spectra_proton[pt_bin][cent][charge][eta_bin]->FindBin(mean);
	    Int_t bin_stop  = h_pt_spectra_proton[pt_bin][cent][charge][eta_bin]->FindBin(1.6);
	    Float_t Counts = 0.0;
	    Float_t Errors = 0.0;
	    for(Int_t bin = bin_start; bin <= bin_stop; bin++)
	    {
	      Counts += 2.0*h_pt_spectra_proton[pt_bin][cent][charge][eta_bin]->GetBinContent(bin)/bin_width[pt_bin];
	      Errors += 2.0*h_pt_spectra_proton[pt_bin][cent][charge][eta_bin]->GetBinError(bin)*h_pt_spectra_proton[pt_bin][cent][charge][eta_bin]->GetBinError(bin)/bin_width[pt_bin];
	    }
	    Float_t pt_center = (pt_low[pt_bin]+pt_up[pt_bin])/2.0;
	    h_Yield[cent][charge][eta_bin]->SetBinContent(h_Yield[cent][charge][eta_bin]->FindBin(pt_center),Counts);
	    h_Yield[cent][charge][eta_bin]->SetBinError(h_Yield[cent][charge][eta_bin]->FindBin(pt_center),TMath::Sqrt(Errors));
	  }
	}
      }
    }
  }

  // Fit the distribution with Boltzmann Distritbution
  TF1 *f_Yield[cent_total][charge_total][eta_total];
  TCanvas *c_Yields[cent_total][charge_total][eta_total];
  for(Int_t cent = cent_start; cent < cent_stop; cent++) // centrality bin
  {
    for(Int_t charge = charge_start; charge < charge_stop; charge++) // charge bin
    {
      for(Int_t eta_bin = eta_start; eta_bin < eta_stop; eta_bin++) // eta gap bin
      {
	TString CanName = Form("c_Yileds_Centrality_%d_Charge_%d_EtaGap_%d",cent,charge,eta_bin);
	c_Yields[cent][charge][eta_bin] = new TCanvas(CanName.Data(),CanName.Data(),1400,10,800,800);
	c_Yields[cent][charge][eta_bin]->cd();
	c_Yields[cent][charge][eta_bin]->cd()->SetLeftMargin(0.15);
	c_Yields[cent][charge][eta_bin]->cd()->SetBottomMargin(0.15);
	c_Yields[cent][charge][eta_bin]->cd()->SetGrid(0,0);
	c_Yields[cent][charge][eta_bin]->cd()->SetTicks(1,1);
	h_Yield[cent][charge][eta_bin]->Draw("pE");
	TString FuncName = Form("f_Poly_Centrality_%d_Charge_%d_EtaGap_%d",cent,charge,eta_bin);
	f_Yield[cent][charge][eta_bin] = new TF1(FuncName.Data(),PtFitFunc2_mod_x,0.2,3.4,7);
	f_Yield[cent][charge][eta_bin]->SetParameter(0,0.938);
	f_Yield[cent][charge][eta_bin]->SetParameter(1,0.250);
	f_Yield[cent][charge][eta_bin]->SetParLimits(1,0,10000);
	f_Yield[cent][charge][eta_bin]->SetParameter(2,6.17689e+10);
	f_Yield[cent][charge][eta_bin]->SetParameter(3,0.001);
	f_Yield[cent][charge][eta_bin]->SetParameter(4,-7.65846e+08);
	f_Yield[cent][charge][eta_bin]->SetParameter(5, 5.30772e+08);
	f_Yield[cent][charge][eta_bin]->SetParameter(6,-5.28569e+07);
	f_Yield[cent][charge][eta_bin]->SetRange(0.2,3.4);
	h_Yield[cent][charge][eta_bin]->Fit(f_Yield[cent][charge][eta_bin],"NRI");
	f_Yield[cent][charge][eta_bin]->Draw("l same");
      }
    }
  }

  TF1 *f_MeanPT_proton[cent_total][charge_total][eta_total];
  for(Int_t cent = cent_start; cent < cent_stop; cent++) // centrality bin
  {
    for(Int_t charge = charge_start; charge < charge_stop; charge++) // charge bin
    {
      for(Int_t eta_bin = eta_start; eta_bin < eta_stop; eta_bin++) // eta gap bin
      {
	TString FuncName = Form("f_MeanPT_Centrality_%d_Charge_%d_EtaGap_%d_proton",cent,charge,eta_bin);
	f_MeanPT_proton[cent][charge][eta_bin] = new TF1(FuncName.Data(),MeanPt,0.2,3.4,7);
	for(Int_t i_par = 0; i_par < 7; i_par++)
	{
	  f_MeanPT_proton[cent][charge][eta_bin]->SetParameter(i_par,f_Yield[cent][charge][eta_bin]->GetParameter(i_par));
	}
	f_MeanPT_proton[cent][charge][eta_bin]->SetRange(0.2,3.4);
      }
    }
  }

  Float_t meanPT_proton[cent_total][charge_total][eta_total][pt_total-1];
  for(Int_t cent = cent_start; cent < cent_stop; cent++) // centrality bin
  {
    for(Int_t charge = charge_start; charge < charge_stop; charge++) // charge bin
    {
      for(Int_t eta_bin = eta_start; eta_bin < eta_stop; eta_bin++) // eta gap bin
      {
	for(Int_t i_pt = 0; i_pt < pt_total-1; i_pt++)
	{
	  meanPT_proton[cent][charge][eta_bin][i_pt] = f_MeanPT_proton[cent][charge][eta_bin]->Integral(pt_low[i_pt],pt_up[i_pt])/f_Yield[cent][charge][eta_bin]->Integral(pt_low[i_pt],pt_up[i_pt]);
//	  if(charge == 0) cout << "pt_bin_low = " << pt_low[i_pt] << ", pt_bin_up = " << pt_up[i_pt] << ", Mean pT proton plus = "  << meanPT_proton[cent][charge][eta_bin][i_pt] << endl;
//	  if(charge == 1) cout << "pt_bin_low = " << pt_low[i_pt] << ", pt_bin_up = " << pt_up[i_pt] << ", Mean pT proton minus = " << meanPT_proton[cent][charge][eta_bin][i_pt] << endl;
	}
      }
    }
  }

  // proton
  for(Int_t cent = cent_start; cent < cent_stop; cent++) // centrality bin
  {
    for(Int_t charge = charge_start; charge < charge_stop; charge++) // charge bin
    {
      for(Int_t eta_bin = eta_start; eta_bin < eta_stop; eta_bin++) // eta gap bin
      {
	for(Int_t i_pt = 1; i_pt < pt_total-1; i_pt++)
	{
	  if(charge == 0 && i_pt == 1) cout << "MeanPT_proton[14] = {";
	  if(charge == 0 && i_pt < pt_total-2) cout << meanPT_proton[cent][charge][eta_bin][i_pt] << ", ";
	  if(charge == 0 && i_pt == pt_total-2) cout << meanPT_proton[cent][charge][eta_bin][i_pt] << "} ";
	  if(charge == 0 && i_pt == pt_total-2) cout << endl;

	  if(charge == 1 && i_pt == 1) cout << "MeanPT_antiproton[14] = {";
	  if(charge == 1 && i_pt < pt_total-2) cout << meanPT_proton[cent][charge][eta_bin][i_pt] << ", ";
	  if(charge == 1 && i_pt == pt_total-2) cout << meanPT_proton[cent][charge][eta_bin][i_pt] << "} ";
	  if(charge == 1 && i_pt == pt_total-2) cout << endl;
	}
      }
    }
  }
}

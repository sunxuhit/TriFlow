#include "TFile.h"
#include "TH1F.h"
#include "TString.h"
#include "TF1.h"
#include "TMath.h"
#include "TCanvas.h"
#include "draw.h"
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

//Blast Wave Fit for spectra
Double_t SpectraFunc(Double_t* x_val, Double_t* par)
{
  Double_t pt, alpha, beta, beta_s, rho, m0, T, R, Inte, Norm;
  pt = x_val[0];
  T = par[0];
  beta_s = par[1];
  m0 = par[2];
  R = par[3];
  Norm = par[4];
  Double_t mt = TMath::Sqrt(pt*pt + m0*m0);
  Double_t r;
  Int_t nbins_r = 100;
  Double_t delta_r = R/((Double_t)nbins_r);
  
  Inte = 0.0;

  for(Int_t j = 0; j < nbins_r+1 ; j++)
  {
    r = j*delta_r;
    rho = TMath::ATanH(beta_s*r/R);
    alpha = (pt/T)*TMath::SinH(rho);
    beta = (mt/T)*TMath::CosH(rho);

    Inte += Norm*pt*r*mt*delta_r*TMath::BesselI0(alpha)*TMath::BesselK1(beta);
  }    

  return Inte;
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

// mEnergy: 0 for 200GeV, 1 for 39GeV
void spectra_proton(Int_t mEnergy = 0)
{
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
  TString inputfile = Form("/project/projectdirs/star/xusun/OutPut/AuAu%s/Mass2_Proton/merged_file/merged_file_%s_M2_Proton_0080_etagap_00.root",Energy[mEnergy].Data(),Energy[mEnergy].Data());
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
	c_Yields[cent][charge][eta_bin]->cd()->SetLogy();
	h_Yield[cent][charge][eta_bin]->Draw("pE");
	TString FuncName = Form("f_Poly_Centrality_%d_Charge_%d_EtaGap_%d",cent,charge,eta_bin);
	f_Yield[cent][charge][eta_bin] = new TF1(FuncName.Data(),PtFitFunc2_mod_x,0.2,3.4,6);
	f_Yield[cent][charge][eta_bin]->SetParameter(0,0.938);
	f_Yield[cent][charge][eta_bin]->SetParLimits(0,0,10000);
	f_Yield[cent][charge][eta_bin]->SetParameter(1,0.125);
	f_Yield[cent][charge][eta_bin]->SetParLimits(1,0,10000);
	f_Yield[cent][charge][eta_bin]->SetParameter(2,7.02565e+08);
	f_Yield[cent][charge][eta_bin]->FixParameter(3,0.0);
	f_Yield[cent][charge][eta_bin]->SetParameter(4,-3.25184e+07);
	f_Yield[cent][charge][eta_bin]->SetParameter(5, 3.01376e+06);
	/*
	f_Yield[cent][charge][eta_bin]->FixParameter(0,0.0);
	f_Yield[cent][charge][eta_bin]->SetParameter(1,3.53541e+08);
	f_Yield[cent][charge][eta_bin]->SetParameter(2,1.37825e+10);
	f_Yield[cent][charge][eta_bin]->SetParameter(3,-2.48249e+10);
	f_Yield[cent][charge][eta_bin]->SetParameter(4,1.82776e+10);
	f_Yield[cent][charge][eta_bin]->SetParameter(5,-6.83682e+09);
	f_Yield[cent][charge][eta_bin]->SetParameter(6,1.29065e+09);
	f_Yield[cent][charge][eta_bin]->SetParameter(7,-9.79036e+07);
	*/
	f_Yield[cent][charge][eta_bin]->SetRange(0.3,3.0);
	h_Yield[cent][charge][eta_bin]->Fit(f_Yield[cent][charge][eta_bin],"MNRI");
//	f_Yield[cent][charge][eta_bin]->SetRange(0.3,3.4);
	f_Yield[cent][charge][eta_bin]->Draw("l same");
      }
    }
  }
}

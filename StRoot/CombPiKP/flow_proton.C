#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "draw.h"
#include "TString.h"
#include "TMath.h"
#include "TProfile.h"
#include "TGraphAsymmErrors.h"

Double_t flow_2(Double_t *x_val, Double_t *par)
{
  Double_t x, y;
  Double_t Ampl, v2;
  x = x_val[0];
  Ampl = par[0];
  v2 = par[1];

  y = Ampl*(1.0 + 2.0*v2*TMath::Cos(2.0*x));

  return y;
}

Double_t flow_3(Double_t *x_val, Double_t *par)
{
  Double_t x, y;
  Double_t Ampl, v3;
  x = x_val[0];
  Ampl = par[0];
  v3 = par[1];

  y = Ampl*(1.0 + 2.0*v3*TMath::Cos(3.0*x));

  return y;
}

// Gauss Function
Double_t GaussPolyFitFunc(Double_t* x_val, Double_t* par)
{
  Double_t x, y, par0, par1, par2, pol0, pol1;
  par0  = par[0];
  par1  = par[1];
  par2  = par[2];
  pol0  = par[3];
  pol1  = par[4];
  x = x_val[0];
  y = par0*TMath::Gaus(x,par1,par2,0) + pol0 + pol1*x;
  return y;
}

// mEnergy: 0 for 200GeV, 1 for 39GeV, 2 for 27GeV
// mCentrality: 0 for 0-80, 1 for 0-10, 2 for 10-40, 3 for 40-80
void flow_proton(Int_t mEnergy = 0, Int_t mCentrality = 0)
{
  const Int_t pt_total = 16;
  const Int_t cent_total = 4; // 4
  const Int_t cent_start[4] = {0,1,2,3};
  const Int_t cent_stop[4]  = {1,2,3,4};
  const Int_t charge_total = 2; // 2
  const Int_t eta_total = 4; // 4
  const Int_t eta_start = 0;
  const Int_t eta_stop = 1;
  const Int_t phi_psi_total = 7; // 7
  const Int_t flow_total = 2; // 2
  const Int_t energy_total = 3; // 2
  TString Order[2] = {"2nd","3rd"};
  Double_t PI_max[2] = {TMath::Pi()/2.0,TMath::Pi()/3.0};
  TString Title[2] = {"#phi-#Psi_{2}","#phi-#Psi_{2}"};
  TString Title_Y[2] = {"v_{2}","v_{3}"};
  TString Centrality[4] = {"0080","0010","1040","4080"};

  Double_t counts_proton[pt_total][cent_total][charge_total][eta_total][phi_psi_total][flow_total];
  Double_t errors_proton[pt_total][cent_total][charge_total][eta_total][phi_psi_total][flow_total];

  Double_t flow_proton[pt_total][cent_total][charge_total][eta_total][flow_total];
  Double_t Errs_proton[pt_total][cent_total][charge_total][eta_total][flow_total];

  // pt bin
  Float_t pt_low[pt_total] = {0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4};
  Float_t pt_up[pt_total]  = {0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,4.2};

  Float_t pt_cut_low[pt_total] = {0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.65,0.7,0.7,0.75,0.879,0.879,0.879,0.879,0.879};
  Float_t pt_cut_up[pt_total]  = {1.2,1.2,1.2,1.2,1.2,1.2,1.2,1.20,1.2,1.2,1.20,1.400,1.400,1.400,1.400,1.400};

  TString Energy[3] = {"200GeV","39GeV","27GeV"};
  TString inputfile = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/Mass2_Proton/merged_file/merged_file_%s_M2_Proton_%s_etagap_00.root",Energy[mEnergy].Data(),Energy[mEnergy].Data(),Centrality[mCentrality].Data());
  cout << inputfile.Data() << endl;
  TFile *file_input = TFile::Open(inputfile.Data());

  // get the histogram
  TH1F *h_Mass2_Proton_EP[pt_total][cent_total][charge_total][eta_total][phi_psi_total][flow_total];
  for(Int_t pt_bin = 0; pt_bin < pt_total; pt_bin++) // pt bin
  {
    for(Int_t cent = cent_start[mCentrality]; cent < cent_stop[mCentrality]; cent++) // centrality bin
    {
      for(Int_t charge = 0; charge < charge_total; charge++) // charge bin
      {
	for(Int_t eta_bin = eta_start; eta_bin < eta_stop; eta_bin++) // eta gap bin
	{
	  for(Int_t phi_psi_bin = 0; phi_psi_bin < phi_psi_total; phi_psi_bin++) // phi-psi bin
	  {
	    for(Int_t flow_bin = 0; flow_bin < flow_total; flow_bin++)
	    {
	      TString HistName = Form("pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_phi_Psi_%d_%s_Proton_EP",pt_bin,cent,charge,eta_bin,phi_psi_bin,Order[flow_bin].Data());
	      h_Mass2_Proton_EP[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin] = (TH1F*)file_input->Get(HistName.Data());
	    }
	  }
	}
      }
    }
  }

  // merge the histogram
  TH1F *h_Mass2_Proton_total[pt_total][cent_total][charge_total][eta_total][flow_total];
  for(Int_t pt_bin = 0; pt_bin < pt_total; pt_bin++) // pt bin
  {
    for(Int_t cent = cent_start[mCentrality]; cent < cent_stop[mCentrality]; cent++) // centrality bin
    {
      for(Int_t charge = 0; charge < charge_total; charge++) // charge bin
      {
	for(Int_t eta_bin = eta_start; eta_bin < eta_stop; eta_bin++) // eta gap bin
	{
	  for(Int_t flow_bin = 0; flow_bin < flow_total; flow_bin++)
	  {
	    for(Int_t phi_psi_bin = 0; phi_psi_bin < phi_psi_total; phi_psi_bin++) // phi-psi bin
	    {
	      if(phi_psi_bin == 0)
	      {
		TString HistName = Form("pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_phi_Psi_%d_%s_Proton_total",pt_bin,cent,charge,eta_bin,phi_psi_bin,Order[flow_bin].Data());
		h_Mass2_Proton_total[pt_bin][cent][charge][eta_bin][flow_bin] = (TH1F*)h_Mass2_Proton_EP[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->Clone(HistName.Data());
	      }
	      if(phi_psi_bin > 0)
	      {
		h_Mass2_Proton_total[pt_bin][cent][charge][eta_bin][flow_bin]->Add(h_Mass2_Proton_EP[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin],1.0); 
	      }
	    }
	  }
	}
      }
    }
  }

  // get the cut range and the counts
  TCanvas *c_Mass2_Proton_total[cent_total][charge_total][eta_total][flow_total];
//  TF1 *GaussFit = new TF1("GaussFit",GaussPolyFitFunc,0.0,1.7,5);
  for(Int_t cent = cent_start[mCentrality]; cent < cent_stop[mCentrality]; cent++) // centrality bin
  {
    for(Int_t charge = 0; charge < charge_total; charge++) // charge bin
    {
      for(Int_t eta_bin = eta_start; eta_bin < eta_stop; eta_bin++) // eta gap bin
      {
	for(Int_t flow_bin = 0; flow_bin < flow_total; flow_bin++)
	{
	  TString CanName_total = Form("Centrality_%d_Charge_%d_EtaGap_%d_%s_total",cent,charge,eta_bin,Order[flow_bin].Data());
	  c_Mass2_Proton_total[cent][charge][eta_bin][flow_bin] = new TCanvas(CanName_total.Data(),CanName_total.Data(),1400,10,1200,1200);
	  c_Mass2_Proton_total[cent][charge][eta_bin][flow_bin]->Divide(4,4); 
	  for(Int_t pt_bin = 0; pt_bin < pt_total; pt_bin++) // pt bin
	  {
	    c_Mass2_Proton_total[cent][charge][eta_bin][flow_bin]->cd(pt_bin+1);
	    c_Mass2_Proton_total[cent][charge][eta_bin][flow_bin]->cd(pt_bin+1)->SetLeftMargin(0.15);
	    c_Mass2_Proton_total[cent][charge][eta_bin][flow_bin]->cd(pt_bin+1)->SetBottomMargin(0.15);
	    c_Mass2_Proton_total[cent][charge][eta_bin][flow_bin]->cd(pt_bin+1)->SetTicks(1,1);
	    c_Mass2_Proton_total[cent][charge][eta_bin][flow_bin]->cd(pt_bin+1)->SetGrid(0,0);
	    h_Mass2_Proton_total[pt_bin][cent][charge][eta_bin][flow_bin]->Draw();
	    /*
	    if(pt_bin >= 5)
	    {
	      // Gauss Fit
	      GaussFit->SetRange(0.5,1.25);
	      for(Int_t i = 0; i < 5; i++)
	      {
		GaussFit->ReleaseParameter(i);
		GaussFit->SetParameter(i,0.0);
		GaussFit->SetParError(i,0.0);
	      }
	      GaussFit->SetParameter(0,100);
	      GaussFit->SetParameter(1,0.85);
	      GaussFit->SetParameter(2,0.1);
	      GaussFit->SetParameter(3,1.0);
	      GaussFit->SetParameter(4,-1.0);
//	      h_Mass2_Proton_total[pt_bin][cent][charge][eta_bin][flow_bin]->Fit(GaussFit,"MR");
	    }
	    */
	    PlotLine(pt_cut_low[pt_bin],pt_cut_low[pt_bin],0.0,h_Mass2_Proton_total[pt_bin][cent][charge][eta_bin][flow_bin]->GetMaximum()/2.0,1,2,2);
	    PlotLine( pt_cut_up[pt_bin], pt_cut_up[pt_bin],0.0,h_Mass2_Proton_total[pt_bin][cent][charge][eta_bin][flow_bin]->GetMaximum()/2.0,1,2,2);
	    Int_t bin_proton_start = h_Mass2_Proton_total[pt_bin][cent][charge][eta_bin][flow_bin]->FindBin(pt_cut_low[pt_bin]+10e-4);
//	    if(pt_bin == 10 && cent == 0 && charge == 1 && eta_bin == 0 && flow_bin == 0)
//	    {
//	      cout << bin_proton_start << endl;
//	    }
	    Int_t bin_proton_stop  = h_Mass2_Proton_total[pt_bin][cent][charge][eta_bin][flow_bin]->FindBin(pt_cut_up[pt_bin]-10e-4);
	    for(Int_t phi_psi_bin = 0; phi_psi_bin < phi_psi_total; phi_psi_bin++)
	    {
	      counts_proton[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin] = 0.0;
	      errors_proton[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin] = 0.0;
	      Float_t Err_proton = 0.0;
	      for(Int_t binx = bin_proton_start; binx < bin_proton_stop+1; binx++)
	      {
		counts_proton[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin] += h_Mass2_Proton_EP[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->GetBinContent(binx);
		Err_proton += h_Mass2_Proton_EP[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->GetBinError(binx)*h_Mass2_Proton_EP[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin]->GetBinError(binx);
	      }
	      errors_proton[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin] = TMath::Sqrt(Err_proton);
	    }
	  }
	}
      }
    }
  }

  // fill TH1F
  TH1F *h_counts_proton[pt_total][cent_total][charge_total][eta_total][flow_total];
  for(Int_t pt_bin = 0; pt_bin < pt_total; pt_bin++)
  {
    for(Int_t cent = cent_start[mCentrality]; cent < cent_stop[mCentrality]; cent++)
    {
      for(Int_t charge = 0; charge < charge_total; charge++)
      {
	for(Int_t eta_bin = eta_start; eta_bin < eta_stop; eta_bin++)
	{
	  for(Int_t flow_bin = 0; flow_bin < flow_total; flow_bin++)
	  {
	    TString HistName = Form("pt_%d_Centrality_%d_Charge_%d_EtaGap_%d_Flow_%d",pt_bin,cent,charge,eta_bin,flow_bin);
	    h_counts_proton[pt_bin][cent][charge][eta_bin][flow_bin] = new TH1F(HistName.Data(),HistName.Data(),7,0,PI_max[flow_bin]);
	    for(Int_t phi_psi_bin = 0; phi_psi_bin < phi_psi_total; phi_psi_bin++)
	    {
	      Float_t bin_center = PI_max[flow_bin]/14.0+phi_psi_bin*PI_max[flow_bin]/7.0;
	      Float_t bin_content = counts_proton[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin];
	      Float_t bin_errors = errors_proton[pt_bin][cent][charge][eta_bin][phi_psi_bin][flow_bin];
//	      cout << "bin_center = " << bin_center << ", bin_content = " << bin_content << endl;
	      Int_t num_bin_center = h_counts_proton[pt_bin][cent][charge][eta_bin][flow_bin]->FindBin(bin_center);
	      h_counts_proton[pt_bin][cent][charge][eta_bin][flow_bin]->SetBinContent(num_bin_center,bin_content);
	      h_counts_proton[pt_bin][cent][charge][eta_bin][flow_bin]->SetBinError(num_bin_center,bin_errors);
	    }
	  }
	}
      }
    }
  }

  // set frame
  TH1F* pos_neg_dummy = new TH1F("pos_neg_dummy","pos_neg_dummy",500,0.,3.4);
  for(Int_t bin_x = 1; bin_x < pos_neg_dummy->GetNbinsX(); bin_x++)
  {
//     pos_neg_dummy->SetBinContent(bin_x,-100.0);
  }
  pos_neg_dummy->SetStats(0);
  pos_neg_dummy->SetTitle("");
  pos_neg_dummy->GetXaxis()->SetTitleOffset(1.0);
  pos_neg_dummy->GetYaxis()->SetTitleOffset(1.2);
  pos_neg_dummy->GetXaxis()->SetLabelSize(0.05);
  pos_neg_dummy->GetYaxis()->SetLabelSize(0.05);
  pos_neg_dummy->GetXaxis()->SetTitleSize(0.06);
  pos_neg_dummy->GetYaxis()->SetTitleSize(0.06);
  pos_neg_dummy->GetXaxis()->SetNdivisions(505,'N');
  pos_neg_dummy->GetYaxis()->SetNdivisions(505,'N');
  pos_neg_dummy->GetXaxis()->CenterTitle();
  pos_neg_dummy->GetYaxis()->CenterTitle();

  // declare TF1
  TF1 *f_proton[pt_total][cent_total][charge_total][eta_total][flow_total];
  for(Int_t pt_bin = 0; pt_bin < pt_total; pt_bin++)
  {
    for(Int_t cent = cent_start[mCentrality]; cent < cent_stop[mCentrality]; cent++)
    {
      for(Int_t charge = 0; charge < charge_total; charge++)
      {
	for(Int_t eta_bin = eta_start; eta_bin < eta_stop; eta_bin++)
	{
	  for(Int_t flow_bin = 0; flow_bin < flow_total; flow_bin++)
	  {
	    if(flow_bin == 0)
	    {
	      TString Flow_proton = Form("flow_pt_%d_Centrality_%d_charge_%d_etagap_%d_%s_proton",pt_bin,cent,charge,eta_bin,Order[flow_bin].Data());
	      f_proton[pt_bin][cent][charge][eta_bin][flow_bin] = new TF1(Flow_proton.Data(),flow_2,0.0,PI_max[flow_bin],2);
	      f_proton[pt_bin][cent][charge][eta_bin][flow_bin]->SetParameter(0,10000.0);
	      f_proton[pt_bin][cent][charge][eta_bin][flow_bin]->SetParameter(1,0.1);
	    }
	    if(flow_bin == 1)
	    {
	      TString Flow_proton = Form("flow_pt_%d_Centrality_%d_charge_%d_etagap_%d_%s_proton",pt_bin,cent,charge,eta_bin,Order[flow_bin].Data());
	      f_proton[pt_bin][cent][charge][eta_bin][flow_bin] = new TF1(Flow_proton.Data(),flow_3,0.0,PI_max[flow_bin],2);
	      f_proton[pt_bin][cent][charge][eta_bin][flow_bin]->SetParameter(0,10000.0);
	      f_proton[pt_bin][cent][charge][eta_bin][flow_bin]->SetParameter(1,0.2);
	    }
	  }
	}
      }
    }
  }

  // fit and draw counts
  TCanvas *c_Counts_Proton[cent_total][charge_total][eta_total][flow_total];
  for(Int_t cent = cent_start[mCentrality]; cent < cent_stop[mCentrality]; cent++) // centrality bin
  {
    for(Int_t charge = 0; charge < charge_total; charge++) // charge bin
    {
      for(Int_t eta_bin = eta_start; eta_bin < eta_stop; eta_bin++) // eta gap bin
      {
	for(Int_t flow_bin = 0; flow_bin < flow_total; flow_bin++)
	{
	  TString CanName = Form("Centrality_%d_Charge_%d_EtaGap_%d_%s",cent,charge,eta_bin,Order[flow_bin].Data());
	  c_Counts_Proton[cent][charge][eta_bin][flow_bin] = new TCanvas(CanName.Data(),CanName.Data(),1400,10,1200,1200);
	  c_Counts_Proton[cent][charge][eta_bin][flow_bin]->Divide(4,4); 
	  for(Int_t pt_bin = 0; pt_bin < pt_total; pt_bin++) // pt bin
	  {
	    c_Counts_Proton[cent][charge][eta_bin][flow_bin]->cd(pt_bin+1);
	    c_Counts_Proton[cent][charge][eta_bin][flow_bin]->cd(pt_bin+1)->SetLeftMargin(0.15);
	    c_Counts_Proton[cent][charge][eta_bin][flow_bin]->cd(pt_bin+1)->SetBottomMargin(0.15);
	    c_Counts_Proton[cent][charge][eta_bin][flow_bin]->cd(pt_bin+1)->SetTicks(1,1);
	    c_Counts_Proton[cent][charge][eta_bin][flow_bin]->cd(pt_bin+1)->SetGrid(0,0);
	    pos_neg_dummy->GetXaxis()->SetTitle(Title[flow_bin].Data());
	    pos_neg_dummy->GetXaxis()->SetRangeUser(0.0,PI_max[flow_bin]);
	    pos_neg_dummy->GetYaxis()->SetTitle("Counts");
	    pos_neg_dummy->GetYaxis()->SetRangeUser(counts_proton[pt_bin][cent][charge][eta_bin][6][flow_bin]*0.99,counts_proton[pt_bin][cent][charge][eta_bin][0][flow_bin]*1.01);
	    pos_neg_dummy->DrawCopy("pE");

	    flow_proton[pt_bin][cent][charge][eta_bin][flow_bin] = 0.0;
	    Errs_proton[pt_bin][cent][charge][eta_bin][flow_bin] = 0.0;

	    h_counts_proton[pt_bin][cent][charge][eta_bin][flow_bin]->Fit(f_proton[pt_bin][cent][charge][eta_bin][flow_bin],"QI");
	    flow_proton[pt_bin][cent][charge][eta_bin][flow_bin] = f_proton[pt_bin][cent][charge][eta_bin][flow_bin]->GetParameter(1);
	    Errs_proton[pt_bin][cent][charge][eta_bin][flow_bin] = f_proton[pt_bin][cent][charge][eta_bin][flow_bin]->GetParError(1);
	    Float_t chi2 = f_proton[pt_bin][cent][charge][eta_bin][flow_bin]->GetChisquare();
	    Float_t ndf  = f_proton[pt_bin][cent][charge][eta_bin][flow_bin]->GetNDF();
	    h_counts_proton[pt_bin][cent][charge][eta_bin][flow_bin]->Draw("PE");
	  }
	}
      }
    }
  }

  // get flow TGraphAsymmErrors

  Float_t mean_res_proton[energy_total][charge_total][cent_total][flow_total];

  /*
  // 200 GeV
  mean_res_proton[0][0][0] = 1.91524;
  mean_res_proton[0][0][1] = 4.13591;
                    
  mean_res_proton[0][1][0] = 1.91616;
  mean_res_proton[0][1][1] = 4.155;
  */

  // 200 GeV
  mean_res_proton[0][0][0][0] = 1.83556; // 200 GeV, pos, 00-80, 2nd
  mean_res_proton[0][0][0][1] = 4.00321; // 200 GeV, pos, 00-80, 3rd
                    
  mean_res_proton[0][1][0][0] = 1.83723; // 200 GeV, neg, 00-80, 2nd
  mean_res_proton[0][1][0][1] = 4.02360; // 200 GeV, neg, 00-80, 3rd

  mean_res_proton[0][0][2][0] = 1.53782; // 200 GeV, pos, 10-40, 2nd
  mean_res_proton[0][0][2][1] = 3.49166; // 200 GeV, pos, 10-40, 3rd
                    
  mean_res_proton[0][1][2][0] = 1.53783; // 200 GeV, neg, 10-40, 2nd
  mean_res_proton[0][1][2][1] = 3.49404; // 200 GeV, neg, 10-40, 3rd

  // 39 GeV
  mean_res_proton[1][0][0][0] = 2.50517; // 39  GeV, pos, 00-80, 2nd
  mean_res_proton[1][0][0][1] = 6.84600; // 39  GeV, pos, 00-80, 3rd
                                                            
  mean_res_proton[1][1][0][0] = 2.52602; // 39  GeV, neg, 00-80, 2nd
  mean_res_proton[1][1][0][1] = 7.20994; // 39  GeV, neg, 00-80, 3rd

  mean_res_proton[1][0][2][0] = 2.02102; // 39  GeV, pos, 10-40, 2nd
  mean_res_proton[1][0][2][1] = 5.76580; // 39  GeV, pos, 10-40, 3rd

  mean_res_proton[1][1][2][0] = 2.02141; // 39  GeV, neg, 10-40, 2nd
  mean_res_proton[1][1][2][1] = 5.79576; // 39  GeV, neg, 10-40, 3rd

  // 27 GeV
  mean_res_proton[2][0][0][0] = 2.73086; // 39  GeV, pos, 00-80, 2nd
  mean_res_proton[2][0][0][1] = 7.36466; // 39  GeV, pos, 00-80, 3rd
                                                            
  mean_res_proton[2][1][0][0] = 2.76464; // 39  GeV, neg, 00-80, 2nd
  mean_res_proton[2][1][0][1] = 7.73661; // 39  GeV, neg, 00-80, 3rd

  /*
  mean_res_proton[1][0][2][0] = 2.02102; // 39  GeV, pos, 10-40, 2nd
  mean_res_proton[1][0][2][1] = 5.76580; // 39  GeV, pos, 10-40, 3rd

  mean_res_proton[1][1][2][0] = 2.02141; // 39  GeV, neg, 10-40, 2nd
  mean_res_proton[1][1][2][1] = 5.79576; // 39  GeV, neg, 10-40, 3rd
  */

  TGraphAsymmErrors *g_flow_proton[cent_total][charge_total][eta_total][flow_total];
  TCanvas *c_Flow_Proton[cent_total][charge_total][eta_total][flow_total];
  for(Int_t cent = cent_start[mCentrality]; cent < cent_stop[mCentrality]; cent++) // centrality bin
  {
    for(Int_t charge = 0; charge < charge_total; charge++) // charge bin
    {
      for(Int_t eta_bin = eta_start; eta_bin < eta_stop; eta_bin++) // eta gap bin
      {
	for(Int_t flow_bin = 0; flow_bin < flow_total; flow_bin++)
	{
	  TString g_Name_proton = Form("g_Centrality_%d_charge_%d_etagap_%d_%s_flow_proton",cent,charge,eta_bin,Order[flow_bin].Data());
	  g_flow_proton[cent][charge][eta_bin][flow_bin] = new TGraphAsymmErrors();
	  g_flow_proton[cent][charge][eta_bin][flow_bin]->SetName(g_Name_proton.Data());
	  for(Int_t pt_bin = 0; pt_bin < pt_total; pt_bin++) // pt bin
	  {
	    g_flow_proton[cent][charge][eta_bin][flow_bin]->SetPoint(pt_bin,(pt_low[pt_bin]+pt_up[pt_bin])/2.0,flow_proton[pt_bin][cent][charge][eta_bin][flow_bin]*mean_res_proton[mEnergy][charge][mCentrality][flow_bin]);
	    g_flow_proton[cent][charge][eta_bin][flow_bin]->SetPointError(pt_bin,0.0,0.0,Errs_proton[pt_bin][cent][charge][eta_bin][flow_bin]*mean_res_proton[mEnergy][charge][mCentrality][flow_bin],Errs_proton[pt_bin][cent][charge][eta_bin][flow_bin]*mean_res_proton[mEnergy][charge][mCentrality][flow_bin]);
	  }
	  TString CanName = Form("Centrality_%d_Charge_%d_EtaGap_%d_%s_flow",cent,charge,eta_bin,Order[flow_bin].Data());
	  c_Flow_Proton[cent][charge][eta_bin][flow_bin] = new TCanvas(CanName.Data(),CanName.Data(),1400,10,800,800);
	  c_Flow_Proton[cent][charge][eta_bin][flow_bin]->cd();
	  c_Flow_Proton[cent][charge][eta_bin][flow_bin]->cd()->SetLeftMargin(0.15);
	  c_Flow_Proton[cent][charge][eta_bin][flow_bin]->cd()->SetBottomMargin(0.15);
	  c_Flow_Proton[cent][charge][eta_bin][flow_bin]->cd()->SetTicks(1,1);
	  c_Flow_Proton[cent][charge][eta_bin][flow_bin]->cd()->SetGrid(0,0);
	  pos_neg_dummy->GetXaxis()->SetTitle("p_{T} (GeV/c)");
	  pos_neg_dummy->GetXaxis()->SetRangeUser(0.0,3.4);
	  pos_neg_dummy->GetYaxis()->SetTitle(Title_Y[flow_bin].Data());
	  pos_neg_dummy->GetYaxis()->SetRangeUser(-0.05,0.25);
	  pos_neg_dummy->DrawCopy("pE");
	  Draw_TGAE_new_Symbol((TGraphAsymmErrors*)g_flow_proton[cent][charge][eta_bin][flow_bin],20,2,0.8);
	  PlotLine(0.0,3.4,0.0,0.0,1,2,2);
	}
      }
    }
  }
  // save flow
  TString output_proton = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/Mass2_Proton/merged_file/flow_proton/Flow_proton_%s_etagap_%d.root",Energy[mEnergy].Data(),Centrality[mCentrality].Data(),eta_start);
  TFile *file_output = new TFile(output_proton.Data(),"RECREATE");
  file_output->cd();
  for(Int_t cent = cent_start[mCentrality]; cent < cent_stop[mCentrality]; cent++)
  {
    for(Int_t charge = 0; charge < charge_total; charge++)
    {
      for(Int_t eta_bin = eta_start; eta_bin < eta_stop; eta_bin++)
      {
	for(Int_t flow_bin = 0; flow_bin < flow_total; flow_bin++)
	{
	  g_flow_proton[cent][charge][eta_bin][flow_bin]->Write();
	}
      }
    }
  }
  file_output->Close();
}

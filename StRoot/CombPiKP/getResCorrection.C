#include "TFile.h"
#include "TString.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TProfile.h"
#include <iostream>
#include "draw.h"

// mEnergy: 0 for 200GeV, 1 for 39GeV, 2 for 27GeV
// mCentrality: 0 for 0-80, 1 for 0-10, 2 for 10-40, 3 for 40-80

void getResCorrection(Int_t mEnergy = 0, Int_t mCentrality = 0)
{
  const Int_t charge_total = 2; // 2
  const Int_t eta_total = 4; // 4
  const Int_t eta_start = 0;
  const Int_t eta_stop = 1;
  const Int_t flow_total = 2; // 2
  const Float_t pt_pion_low = -0.1;
  const Float_t pt_pion_up = 0.1;
  const Float_t pt_kaon_low = 0.15;
  const Float_t pt_kaon_up = 0.3;
  const Float_t pt_proton_low = 0.6;
  const Float_t pt_proton_up = 1.2;

  Float_t Counts_Pion[9][charge_total][eta_total];
  Float_t Counts_Kaon[9][charge_total][eta_total];
  Float_t Counts_Proton[9][charge_total][eta_total];

  TString Energy[3] = {"200GeV","39GeV","27GeV"};
  TString Centrality[4] = {"0080","0010","1040","4080"};
  TString inputfile = Form("/project/projectdirs/star/xusun/OutPut/AuAu%s/Yields/merged_file/merged_file_%s_Yields_0070_etagap_00.root",Energy[mEnergy].Data(),Energy[mEnergy].Data());
  TFile *file_input = TFile::Open(inputfile.Data());


  Int_t cent_low[4] = {0,7,4,0}; // 0 = 0-80%, 1 = 0-10%, 2 = 10-40%, 3 = 40-80%
  Int_t cent_up[4]  = {8,8,6,3}; // 0 = 0-80%, 1 = 0-10%, 2 = 10-40%, 3 = 40-80%
  // particle yields
  // 0 = centrality: 0 = 70-80%, 1 = 60-70%, 2 = 50-60%, 3 = 40-50%, 4 = 30-40%, 5 = 20-30%, 6 = 10-20%, 7 = 5-10%, 8 = 0-5%
  // 1 = charge: 0 = pos, 1 = neg
  TH2F *h_mMass2_nSigmaPion_Yields_PiK_EP[9][2][4];
  TH1F *h_mMass2_Yields_Proton_EP[9][2][4];
  for(Int_t cent = cent_low[mCentrality]; cent <= cent_up[mCentrality]; cent++) // centrality bin
  {
    for(Int_t charge = 0; charge < charge_total; charge++) // charge bin
    {
      for(Int_t eta_bin = eta_start; eta_bin < eta_stop; eta_bin++) // eta gap bin
      {
	TString HistName;
	HistName = Form("Centrality_%d_Charge_%d_EtaGap_%d_Yields_PiK_EP",cent,charge,eta_bin);
	h_mMass2_nSigmaPion_Yields_PiK_EP[cent][charge][eta_bin] = (TH2F*)file_input->Get(HistName.Data());
	Counts_Pion[cent][charge][eta_bin] = 0.0;
	Counts_Kaon[cent][charge][eta_bin] = 0.0;
	HistName = Form("Centrality_%d_Charge_%d_EtaGap_%d_Yields_Proton_EP",cent,charge,eta_bin);
	h_mMass2_Yields_Proton_EP[cent][charge][eta_bin] = (TH1F*)file_input->Get(HistName.Data());
	Counts_Proton[cent][charge][eta_bin] = 0.0;
      }
    }
  }

  // get range
  TCanvas *c_Yields_PiK[charge_total][eta_total];
  for(Int_t charge = 0; charge < charge_total; charge++) // charge bin
  {
    for(Int_t eta_bin = eta_start; eta_bin < eta_stop; eta_bin++) // eta gap bin
    {
      TString CanName = Form("charge_%d_etagap_%d_pik",charge,eta_bin);
      c_Yields_PiK[charge][eta_bin] = new TCanvas(CanName.Data(),CanName.Data(),1400,10,900,900);
      c_Yields_PiK[charge][eta_bin]->Divide(3,3);
//      cout << CanName.Data() << endl;
      for(Int_t cent = cent_low[mCentrality]; cent <= cent_up[mCentrality]; cent++)
      {
	c_Yields_PiK[charge][eta_bin]->cd(cent+1);
	c_Yields_PiK[charge][eta_bin]->cd(cent+1)->SetLeftMargin(0.15);
	c_Yields_PiK[charge][eta_bin]->cd(cent+1)->SetBottomMargin(0.15);
	c_Yields_PiK[charge][eta_bin]->cd(cent+1)->SetTicks(1,1);
	c_Yields_PiK[charge][eta_bin]->cd(cent+1)->SetGrid(0,0);
	h_mMass2_nSigmaPion_Yields_PiK_EP[cent][charge][eta_bin]->GetXaxis()->SetTitle("x (m^{2},n#sigma_{#pi})"); 
	h_mMass2_nSigmaPion_Yields_PiK_EP[cent][charge][eta_bin]->GetXaxis()->SetRangeUser(-0.3,0.6);
	h_mMass2_nSigmaPion_Yields_PiK_EP[cent][charge][eta_bin]->GetYaxis()->SetTitle("Counts"); 
	h_mMass2_nSigmaPion_Yields_PiK_EP[cent][charge][eta_bin]->ProjectionX()->Draw("l");
	PlotLine(pt_pion_low,pt_pion_low,0.0,h_mMass2_nSigmaPion_Yields_PiK_EP[cent][charge][eta_bin]->ProjectionX()->GetMaximum()/2.0,1,2,2);
	PlotLine(pt_pion_up,pt_pion_up,0.0,h_mMass2_nSigmaPion_Yields_PiK_EP[cent][charge][eta_bin]->ProjectionX()->GetMaximum()/2.0,1,2,2);
	PlotLine(pt_kaon_low,pt_kaon_low,0.0,h_mMass2_nSigmaPion_Yields_PiK_EP[cent][charge][eta_bin]->ProjectionX()->GetMaximum()/4.0,4,2,2);
	PlotLine(pt_kaon_up,pt_kaon_up,0.0,h_mMass2_nSigmaPion_Yields_PiK_EP[cent][charge][eta_bin]->ProjectionX()->GetMaximum()/4.0,4,2,2);

	Int_t pion_start = h_mMass2_nSigmaPion_Yields_PiK_EP[cent][charge][eta_bin]->ProjectionX()->FindBin(pt_pion_low);
	Int_t pion_stop = h_mMass2_nSigmaPion_Yields_PiK_EP[cent][charge][eta_bin]->ProjectionX()->FindBin(pt_pion_up);
	for(Int_t binx = pion_start; binx < pion_stop+1; binx++)
	{
	  Counts_Pion[cent][charge][eta_bin] += h_mMass2_nSigmaPion_Yields_PiK_EP[cent][charge][eta_bin]->ProjectionX()->GetBinContent(binx);
	}
//	cout << "pion = " << Counts_Pion[cent][charge][eta_bin] << endl; 

	Int_t kaon_start = h_mMass2_nSigmaPion_Yields_PiK_EP[cent][charge][eta_bin]->ProjectionX()->FindBin(pt_kaon_low);
	Int_t kaon_stop = h_mMass2_nSigmaPion_Yields_PiK_EP[cent][charge][eta_bin]->ProjectionX()->FindBin(pt_kaon_up);
	for(Int_t binx = kaon_start; binx < kaon_stop+1; binx++)
	{
	  Counts_Kaon[cent][charge][eta_bin] += h_mMass2_nSigmaPion_Yields_PiK_EP[cent][charge][eta_bin]->ProjectionX()->GetBinContent(binx);
	}
//	cout << "kaon = " << Counts_Kaon[cent][charge][eta_bin] << endl; 
      }
    }
  }

  TCanvas *c_Yields_Proton[charge_total][eta_total];
  for(Int_t charge = 0; charge < charge_total; charge++) // charge bin
  {
    for(Int_t eta_bin = eta_start; eta_bin < eta_stop; eta_bin++) // eta gap bin
    {
      TString CanName = Form("charge_%d_etagap_%d_proton",charge,eta_bin);
      c_Yields_Proton[charge][eta_bin] = new TCanvas(CanName.Data(),CanName.Data(),1400,10,900,900);
      c_Yields_Proton[charge][eta_bin]->Divide(3,3);
//      cout << CanName.Data() << endl;
      for(Int_t cent = cent_low[mCentrality]; cent <= cent_up[mCentrality]; cent++)
      {
	c_Yields_Proton[charge][eta_bin]->cd(cent+1);
	c_Yields_Proton[charge][eta_bin]->cd(cent+1)->SetLeftMargin(0.15);
	c_Yields_Proton[charge][eta_bin]->cd(cent+1)->SetBottomMargin(0.15);
	c_Yields_Proton[charge][eta_bin]->cd(cent+1)->SetTicks(1,1);
	c_Yields_Proton[charge][eta_bin]->cd(cent+1)->SetGrid(0,0);
	h_mMass2_Yields_Proton_EP[cent][charge][eta_bin]->GetXaxis()->SetTitle("mass^{2}"); 
	h_mMass2_Yields_Proton_EP[cent][charge][eta_bin]->GetXaxis()->SetRangeUser(-0.3,1.7);
	h_mMass2_Yields_Proton_EP[cent][charge][eta_bin]->GetYaxis()->SetTitle("Counts"); 
	h_mMass2_Yields_Proton_EP[cent][charge][eta_bin]->Draw("l");
	PlotLine(pt_proton_low,pt_proton_low,0.0,h_mMass2_Yields_Proton_EP[cent][charge][eta_bin]->GetMaximum()/2.0,1,2,2);
	PlotLine(pt_proton_up,pt_proton_up,0.0,h_mMass2_Yields_Proton_EP[cent][charge][eta_bin]->GetMaximum()/2.0,1,2,2);

	Int_t proton_start = h_mMass2_Yields_Proton_EP[cent][charge][eta_bin]->FindBin(pt_proton_low);
	Int_t proton_stop = h_mMass2_Yields_Proton_EP[cent][charge][eta_bin]->FindBin(pt_proton_up);
	for(Int_t binx = proton_start; binx < proton_stop+1; binx++)
	{
	  Counts_Proton[cent][charge][eta_bin] += h_mMass2_Yields_Proton_EP[cent][charge][eta_bin]->GetBinContent(binx);
	}
//	cout << "charge = " << charge << ", eta_bin = " << eta_bin << ", proton = " << Counts_Proton[cent][charge][eta_bin] << endl; 
      }
    }
  }

  // Resolution Correction
  TString input_res = Form("/project/projectdirs/star/xusun/OutPut/AuAu%s/Resolution/file_%s_Resolution.root",Energy[mEnergy].Data(),Energy[mEnergy].Data());
  TFile *input = TFile::Open(input_res.Data());
  TProfile *p_res2[eta_total];
  TProfile *p_res3[eta_total];
  Double_t mean_res_2_pion[charge_total][eta_total];
  Double_t mean_res_2_kaon[charge_total][eta_total];
  Double_t mean_res_2_proton[charge_total][eta_total];
  Double_t mean_res_3_pion[charge_total][eta_total];
  Double_t mean_res_3_kaon[charge_total][eta_total];
  Double_t mean_res_3_proton[charge_total][eta_total];
  for(Int_t eta_bin = eta_start; eta_bin < eta_stop; eta_bin++)
  {
    TString ProName;
    ProName = Form("Res2_EtaGap_%d_EP",eta_bin);
    p_res2[eta_bin] = (TProfile*)input->Get(ProName.Data());
    ProName = Form("Res3_EtaGap_%d_EP",eta_bin);
    p_res3[eta_bin] = (TProfile*)input->Get(ProName.Data());
  }
  for(Int_t charge = 0; charge < charge_total; charge++) // charge bin
  {
    for(Int_t eta_bin = eta_start; eta_bin < eta_stop; eta_bin++) // eta gap bin
    {
      mean_res_2_pion[charge][eta_bin] = 0.0;
      mean_res_2_kaon[charge][eta_bin] = 0.0;
      mean_res_2_proton[charge][eta_bin] = 0.0;
      mean_res_3_pion[charge][eta_bin] = 0.0;
      mean_res_3_kaon[charge][eta_bin] = 0.0;
      mean_res_3_proton[charge][eta_bin] = 0.0;

      Double_t res_2[9];
      Float_t pion_total = 0.0;
      Float_t kaon_total = 0.0;
      Float_t proton_total = 0.0;
      for(Int_t cent = cent_low[mCentrality]; cent <= cent_up[mCentrality]; cent++)
      {
	res_2[cent] = TMath::Sqrt(p_res2[eta_bin]->GetBinContent(cent+1));
	pion_total += Counts_Pion[cent][charge][eta_bin];
	kaon_total += Counts_Kaon[cent][charge][eta_bin];
	proton_total += Counts_Proton[cent][charge][eta_bin];
      }
      for(Int_t cent = cent_low[mCentrality]; cent <= cent_up[mCentrality]; cent++)
      {
	mean_res_2_pion[charge][eta_bin] += Counts_Pion[cent][charge][eta_bin]/(res_2[cent]*pion_total);
	mean_res_2_kaon[charge][eta_bin] += Counts_Kaon[cent][charge][eta_bin]/(res_2[cent]*kaon_total);
	mean_res_2_proton[charge][eta_bin] += Counts_Proton[cent][charge][eta_bin]/(res_2[cent]*proton_total);
      }
      cout << "charge = " << charge << ", eta_bin = " << eta_bin << ", mean_res_2_pion = " << mean_res_2_pion[charge][eta_bin] << endl;
      cout << "charge = " << charge << ", eta_bin = " << eta_bin << ", mean_res_2_kaon = " << mean_res_2_kaon[charge][eta_bin] << endl;
      cout << "charge = " << charge << ", eta_bin = " << eta_bin << ", mean_res_2_proton = " << mean_res_2_proton[charge][eta_bin] << endl;

      Double_t res_3[9];
      for(Int_t i = 0; i < 9; i++)
      {
	res_3[i] = -999.9;
      }

      for(Int_t cent = cent_low[mCentrality]; cent <= cent_up[mCentrality]; cent++)
      {
	if(p_res3[eta_bin]->GetBinContent(cent+1) > 0)
	{
	  res_3[cent] = TMath::Sqrt(p_res3[eta_bin]->GetBinContent(cent+1));
	}
      }
      for(Int_t cent = cent_low[mCentrality]; cent <= cent_up[mCentrality]; cent++)
      {
	mean_res_3_pion[charge][eta_bin] += Counts_Pion[cent][charge][eta_bin]/(res_3[cent]*pion_total);
	mean_res_3_kaon[charge][eta_bin] += Counts_Kaon[cent][charge][eta_bin]/(res_3[cent]*kaon_total);
	mean_res_3_proton[charge][eta_bin] += Counts_Proton[cent][charge][eta_bin]/(res_3[cent]*proton_total);
      }
      cout << "charge = " << charge << ", eta_bin = " << eta_bin << ", mean_res_3_pion = " << mean_res_3_pion[charge][eta_bin] << endl;
      cout << "charge = " << charge << ", eta_bin = " << eta_bin << ", mean_res_3_kaon = " << mean_res_3_kaon[charge][eta_bin] << endl;
      cout << "charge = " << charge << ", eta_bin = " << eta_bin << ", mean_res_3_proton = " << mean_res_3_proton[charge][eta_bin] << endl;
    }
  }
}

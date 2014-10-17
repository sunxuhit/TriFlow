#ifndef StStrangenessHistoManger_h
#define StStrangenessHistoManger_h

#include "StMessMgr.h"

class TH1F;

class StStrangenessHistoManger
{
  public:
    StStrangenessHistoManger();
    ~StStrangenessHistoManger();

    void Init(Int_t X_flag, Int_t mode);
    void Fill(Float_t pt, Int_t Cent9, Int_t eta_gap, Float_t phi_psi2, Float_t Res2, Float_t phi_psi3, Float_t Res3, Float_t Mass2, Double_t reweight);
    void Fill_sub(Float_t pt, Int_t Cent9, Int_t eta_gap, Float_t phi_psi2, Float_t Res2, Float_t phi_psi3, Float_t Res3, Float_t Mass2, Double_t reweight);
    void Write();

  private:
    // flow analysis
    // 0 = pt bin
    // 1 = centrality: 0 = 0-80%(0-70%), 1 = 0-10%, 2 = 10-40%, 3 = 40-80%(40-70%)
    // 2 = eta_gap
    // 3 = phi - Psi
    TH1F *h_mMass2_EP[23][4][4][7];
    TH1F *h_mMass3_EP[23][4][4][7];
    // subtract k0s
    TH1F *h_mMass2_EP_sub[23][4][4][7];
    TH1F *h_mMass3_EP_sub[23][4][4][7];

    // raw pt spectra
    // 0 = pt bin
    // 1 = centrality: 0 = 0-80%(0-70%), 1 = 0-10%, 2 = 10-40%, 3 = 40-80%(40-70%)
    // 2 = eta_gap
    TH1F *h_mMass_Spec[23][4][4];
    // subtract k0s
    TH1F *h_mMass_Spec_sub[23][4][4];

    // event plane resolution correction
    // 0 = centrality
    // 1 = eta_gap
    TH1F *h_mMass_Yields[9][4];
    // subtract k0s
    TH1F *h_mMass_Yields_sub[9][4];

  ClassDef(StStrangenessHistoManger,1)
};
#endif

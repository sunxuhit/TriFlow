#ifndef StTriFlowHistoManger_h
#define StTriFlowHistoManger_h

#include "StMessMgr.h"

class TH2F;
class TH1F;
class StPicoTrack;
class StTriFlowCut;

class StTriFlowHistoManger
{
  public:
    StTriFlowHistoManger(Int_t energy);
    ~StTriFlowHistoManger();

    void InitHist();
    void InitProton();
    void InitYields_nSigPion();
    void InitYields_Proton();
    void InitQA_Detector();

    void FillHist(Float_t pt, Int_t Cent9, Int_t charge_bin, Int_t eta_gap, Float_t phi_psi2, Float_t Res2, Float_t phi_psi3, Float_t Res3, Float_t New_X, Float_t New_Y, Double_t reweight);
    void FillProton(Float_t pt, Int_t Cent9, Int_t charge_bin, Int_t eta_gap, Float_t phi_psi2, Float_t Res2, Float_t phi_psi3, Float_t Res3, Float_t Mass2, Double_t reweight);
    void FillYields_PiK(Int_t Cent9, Int_t charge_bin, Int_t eta_gap, Float_t New_X, Float_t New_Y, Double_t reweight);
    void FillYields_Proton(Int_t Cent9, Int_t charge_bin, Int_t eta_gap, Float_t Mass2, Double_t reweight);
    void FillQA_before(Int_t eta_gap, Float_t Mass2, Float_t dEdx, Float_t pq);
    void FillQA_after(Int_t eta_gap, Float_t Mass2, Float_t dEdx, Float_t pq);
    void FillToFLocal(StPicoTrack*);
    void FillQA_Event(Int_t RefMult, Float_t Vz);
    void FillQA_Detector(Float_t dEdx, Float_t Mass2, Float_t p);

    void WriteHist();
    void WriteProton();
    void WriteYileds_nSigPion();
    void WriteYileds_Proton();
    void WriteQA();
    void WriteQA_Detector();
    
  private:
    // 0 = pt bin
    // 1 = centrality: 0 = 0-80%(0-70%), 1 = 0-10%, 2 = 10-40%, 3 = 40-80%(40-70%)
    // 2 = charge: 0 = pos, 1 = neg
    // 3 = eta_gap
    // 4 = phi - Psi
    TH2F *h_mMass2_nSigmaPion2_EP[16][4][2][4][7];
    TH2F *h_mMass2_nSigmaPion3_EP[16][4][2][4][7];
//    TH2F *h_mMass2_nSigmaPion2_SP[14][4][2][4][7];
//    TH2F *h_mMass2_nSigmaPion3_SP[14][4][2][4][7];
//    TH2F *h_mMass2_nSigmaPion_Test[13][4][2];
    TH2F *h_mToFYLocal_Mass2;
    TH2F *h_mToFZLocal_Mass2;
    TH1F *h_mMass2_Proton2_EP[16][4][2][4][7];
    TH1F *h_mMass2_Proton3_EP[16][4][2][4][7];

    // raw pt spectra
    // 0 = pt bin
    // 1 = centrality: 0 = 0-80%(0-70%), 1 = 0-10%, 2 = 10-40%, 3 = 40-80%(40-70%)
    // 2 = charge: 0 = pos, 1 = neg
    // 3 = eta_gap
    TH2F *h_pt_spectra2_pik[16][4][2][4];
    TH2F *h_pt_spectra3_pik[16][4][2][4];
    TH1F *h_pt_spectra2_proton[16][4][2][4];
    TH1F *h_pt_spectra3_proton[16][4][2][4];

    // particle yields
    // 0 = centrality: 0 = 70-80%, 1 = 60-70%, 2 = 50-60%, 3 = 40-50%, 4 = 30-40%, 5 = 20-30%, 6 = 10-20%, 7 = 5-10%, 8 = 0-5%
    // 1 = charge: 0 = pos, 1 = neg
    // 2 = eta_gap
    TH2F *h_mMass2_nSigmaPion_Yields_PiK_EP[9][2][4];
    TH1F *h_mMass2_Yields_Proton_EP[9][2][4];

    TH1F *h_phi_psi2;
    TH1F *h_phi_psi3;
    TH1F *h_yield_phi_psi2;
    TH1F *h_yield_phi_psi3;

    // 0 = eta_gap
    TH2F *h_mDEdx_pq_before[4];
    TH2F *h_mDEdx_pq_after[4];
    TH2F *h_mMass2_pq_before[4];
    TH2F *h_mMass2_pq_after[4];

    // QA plots
    TH2F *h_mDEdx;
    TH2F *h_mMass2;

    TH2F *h_mMass2_pt;
    TH1F *h_mRefMult;
    TH1F *h_mVz;

    Int_t mEnergy;

    StTriFlowCut *mTriFlowCut;

  ClassDef(StTriFlowHistoManger,1)
};
#endif

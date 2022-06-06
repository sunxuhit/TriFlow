#ifndef StTriFlowHistManger_h
#define StTriFlowHistManger_h

#include "StMessMgr.h"

class TH1F;

class StTriFlowHistManger
{
  public:
    StTriFlowHistManger();
    ~StTriFlowHistManger();

    void InitEventPlane();

    void FillEventPlane_East_EP(Float_t, Float_t, Float_t, Float_t, Int_t, Int_t);
    void FillEventPlane_West_EP(Float_t, Float_t, Float_t, Float_t, Int_t, Int_t);
    void FillEventPlane_A_EP(Float_t, Float_t, Float_t, Float_t, Int_t);
    void FillEventPlane_B_EP(Float_t, Float_t, Float_t, Float_t, Int_t);
    void FillEventPlane_Full_EP(Float_t, Float_t, Float_t, Float_t, Int_t);

    void FillEventPlane_East_SP(Float_t, Float_t, Float_t, Float_t, Int_t, Int_t);
    void FillEventPlane_West_SP(Float_t, Float_t, Float_t, Float_t, Int_t, Int_t);
    void FillEventPlane_A_SP(Float_t, Float_t, Float_t, Float_t, Int_t);
    void FillEventPlane_B_SP(Float_t, Float_t, Float_t, Float_t, Int_t);
    void FillEventPlane_Full_SP(Float_t, Float_t, Float_t, Float_t, Int_t);

    void WriteEventPlane();

  private:
    // Event Plane method
    // ReCenter Correction
    TH1F *h_mPsi2_East_ReCenter_EP[10][4]; // 0: 0-8 = Centrality, 9 = minBias | 1: eta_gap
    TH1F *h_mPsi2_West_ReCenter_EP[10][4]; // 0: 0-8 = Centrality, 9 = minBias | 1: eta_gap
    TH1F *h_mPsi2_A_ReCenter_EP[10];       // 0: 0-8 = Centrality, 9 = minBias
    TH1F *h_mPsi2_B_ReCenter_EP[10];       // 0: 0-8 = Centrality, 9 = minBias
    TH1F *h_mPsi2_Full_ReCenter_EP[10];    // 0: 0-8 = Centrality, 9 = minBias

    TH1F *h_mPsi3_East_ReCenter_EP[10][4]; // 0: 0-8 = Centrality, 9 = minBias | 1: eta_gap
    TH1F *h_mPsi3_West_ReCenter_EP[10][4]; // 0: 0-8 = Centrality, 9 = minBias | 1: eta_gap
    TH1F *h_mPsi3_A_ReCenter_EP[10];       // 0: 0-8 = Centrality, 9 = minBias
    TH1F *h_mPsi3_B_ReCenter_EP[10];       // 0: 0-8 = Centrality, 9 = minBias
    TH1F *h_mPsi3_Full_ReCenter_EP[10];    // 0: 0-8 = Centrality, 9 = minBias
    // Shift Correction
    TH1F *h_mPsi2_East_Shift_EP[10][4]; // 0: 0-8 = Centrality, 9 = minBias | 1: eta_gap
    TH1F *h_mPsi2_West_Shift_EP[10][4]; // 0: 0-8 = Centrality, 9 = minBias | 1: eta_gap
    TH1F *h_mPsi2_A_Shift_EP[10];       // 0: 0-8 = Centrality, 9 = minBias
    TH1F *h_mPsi2_B_Shift_EP[10];       // 0: 0-8 = Centrality, 9 = minBias
    TH1F *h_mPsi2_Full_Shift_EP[10];    // 0: 0-8 = Centrality, 9 = minBias

    TH1F *h_mPsi3_East_Shift_EP[10][4]; // 0: 0-8 = Centrality, 9 = minBias | 1: eta_gap
    TH1F *h_mPsi3_West_Shift_EP[10][4]; // 0: 0-8 = Centrality, 9 = minBias | 1: eta_gap
    TH1F *h_mPsi3_A_Shift_EP[10];       // 0: 0-8 = Centrality, 9 = minBias
    TH1F *h_mPsi3_B_Shift_EP[10];       // 0: 0-8 = Centrality, 9 = minBias
    TH1F *h_mPsi3_Full_Shift_EP[10];    // 0: 0-8 = Centrality, 9 = minBias

    // Scalor Product method
    // ReCenter Correction
    TH1F *h_mPsi2_East_ReCenter_SP[10][4]; // 0: 0-8 = Centrality, 9 = minBias | 1: eta_gap
    TH1F *h_mPsi2_West_ReCenter_SP[10][4]; // 0: 0-8 = Centrality, 9 = minBias | 1: eta_gap
    TH1F *h_mPsi2_A_ReCenter_SP[10];       // 0: 0-8 = Centrality, 9 = minBias
    TH1F *h_mPsi2_B_ReCenter_SP[10];       // 0: 0-8 = Centrality, 9 = minBias
    TH1F *h_mPsi2_Full_ReCenter_SP[10];    // 0: 0-8 = Centrality, 9 = minBias

    TH1F *h_mPsi3_East_ReCenter_SP[10][4]; // 0: 0-8 = Centrality, 9 = minBias | 1: eta_gap
    TH1F *h_mPsi3_West_ReCenter_SP[10][4]; // 0: 0-8 = Centrality, 9 = minBias | 1: eta_gap
    TH1F *h_mPsi3_A_ReCenter_SP[10];       // 0: 0-8 = Centrality, 9 = minBias
    TH1F *h_mPsi3_B_ReCenter_SP[10];       // 0: 0-8 = Centrality, 9 = minBias
    TH1F *h_mPsi3_Full_ReCenter_SP[10];    // 0: 0-8 = Centrality, 9 = minBias
    // Shift Correction
    TH1F *h_mPsi2_East_Shift_SP[10][4]; // 0: 0-8 = Centrality, 9 = minBias | 1: eta_gap
    TH1F *h_mPsi2_West_Shift_SP[10][4]; // 0: 0-8 = Centrality, 9 = minBias | 1: eta_gap
    TH1F *h_mPsi2_A_Shift_SP[10];       // 0: 0-8 = Centrality, 9 = minBias
    TH1F *h_mPsi2_B_Shift_SP[10];       // 0: 0-8 = Centrality, 9 = minBias
    TH1F *h_mPsi2_Full_Shift_SP[10];    // 0: 0-8 = Centrality, 9 = minBias

    TH1F *h_mPsi3_East_Shift_SP[10][4]; // 0: 0-8 = Centrality, 9 = minBias | 1: eta_gap
    TH1F *h_mPsi3_West_Shift_SP[10][4]; // 0: 0-8 = Centrality, 9 = minBias | 1: eta_gap
    TH1F *h_mPsi3_A_Shift_SP[10];       // 0: 0-8 = Centrality, 9 = minBias
    TH1F *h_mPsi3_B_Shift_SP[10];       // 0: 0-8 = Centrality, 9 = minBias
    TH1F *h_mPsi3_Full_Shift_SP[10];    // 0: 0-8 = Centrality, 9 = minBias

  ClassDef(StTriFlowHistManger,1)
};

#endif

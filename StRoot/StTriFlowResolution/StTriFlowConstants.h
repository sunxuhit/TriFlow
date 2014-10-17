#ifndef StTriFolwConstants_h
#define StTriFolwConstants_h
#include "TString.h"
#include "Rtypes.h"
#include <map>

class TriFlow
{
  public:

  // event cut
  static std::map<float,float> mVzMaxMap;
  static Float_t mVrMax; // 2.0
  static Int_t mMatchedToFMin; // 2

  // track cut
  static Float_t mDcaMax; // 1.0
  static Int_t mHitsDedxMin; // 15
  static Int_t mHitsFitTPCMin; // 15
  static Int_t mHitsMaxTPCMin; // 0
  static Float_t mHitsRatioTPCMin; // 0.52
  static Float_t mEtaMax; // 1.0
  static Float_t mPrimPtMin; // 0.2
  static Float_t mPrimPtMax; // 2.0
  static Float_t mPrimPtWeight; // 2.0
  static Float_t mPrimMomMax; // 10.0
  static Float_t mMass2Min; // -10.0
  static Double_t MAGFIELDFACTOR; // kilogauss  
  static Int_t mTrackMin; // 2
  static Int_t mTrackMin_Full; // 4
  static Float_t mToFYLocalMax; // 1.7
  static Float_t mToFZLocalMax; // 1.7
  static Float_t mNSigmaProtonMax; // 2.5
  static std::map<float,float> mSigScaleMap;

  // used constant
  static Float_t mShiftOrder2[5];
  static Float_t mShiftOrder3[5];

  static Float_t mEta_Gap[4]; // 0 = 0.05, 1 = 0.1, 2 = 0.2, 3 = 0.5

  static Float_t pt_low[16];
  static Float_t pt_up[16];

  static Float_t x_low[16];
  static Float_t x_up[16];
  static Float_t y_low[16];
  static Float_t y_up[16];

  static Int_t cent_low[4];
  static Int_t cent_up[4];
  static TString Centrality_01[4];
  static TString Centrality_23[4];

  static Float_t phi_Psi2_low[7];
  static Float_t phi_Psi2_up[7];
  static Float_t phi_Psi3_low[7];
  static Float_t phi_Psi3_up[7];
  static Float_t Psi2_low[3];
  static Float_t Psi2_up[3];
  static Float_t Psi3_low[5];
  static Float_t Psi3_up[5];

  static Int_t pt_total;

  static Int_t Centrality_total;
  static Int_t Centrality_start;
  static Int_t Centrality_stop;

  static Int_t EtaGap_total;
  static Int_t EtaGap_start;
  static Int_t EtaGap_stop;

  static Int_t Charge_total;
  static Int_t Charge_start;
  static Int_t Charge_stop;

  static TString Energy[3];


  private:

  static std::map<float,float> createVzMaxMap()
  {
    std::map<float,float> m;

    m[7.7] = 70.0;
    m[11.5] = 50.0;
    m[19.6] = 70.0;
    m[27.0] = 70.0;
    m[39.0] = 40.0;
    m[62.4] = 40.0;
    m[200.0] = 30.0;

    return m;
  }

  static std::map<float,float> createSigScaleMap()
  {
    std::map<float,float> m;

    m[7.7] = 1.0;
    m[11.5] = 1.0;
    m[19.6] = 1.0;
    m[27.0] = 1.9;
    m[39.0] = 1.0;
    m[62.4] = 1.0;
    m[200.0] = 1.0;

    return m;
  }

  ClassDef(TriFlow, 1)
};

#endif

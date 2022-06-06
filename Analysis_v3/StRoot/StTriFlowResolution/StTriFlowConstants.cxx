#include "StTriFlowConstants.h"
#include "TMath.h"
#include "StarClassLibrary/SystemOfUnits.h"

ClassImp(TriFlow)

// event cut
std::map<float,float> TriFlow::mVzMaxMap = TriFlow::createVzMaxMap();
Float_t TriFlow::mVrMax = 2.0;
Int_t TriFlow::mMatchedToFMin = 2;

// track cut
std::map<float,float> TriFlow::mSigScaleMap = TriFlow::createSigScaleMap();
Float_t TriFlow::mDcaMax = 1.0;
Int_t TriFlow::mHitsDedxMin = 5;
Int_t TriFlow::mHitsFitTPCMin = 15;
Int_t TriFlow::mHitsMaxTPCMin = 0;
Float_t TriFlow::mHitsRatioTPCMin = 0.52;
Float_t TriFlow::mEtaMax = 1.0;
Float_t TriFlow::mPrimPtMin = 0.2;
Float_t TriFlow::mPrimPtMax = 2.0;
Float_t TriFlow::mPrimPtWeight = 2.0;
Float_t TriFlow::mPrimMomMax = 10.0;
Float_t TriFlow::mMass2Min = -10.0;
Double_t TriFlow::MAGFIELDFACTOR = kilogauss;
Int_t TriFlow::mTrackMin = 2;
Int_t TriFlow::mTrackMin_Full = 4;
Float_t TriFlow::mToFYLocalMax = 1.7;
Float_t TriFlow::mToFZLocalMax = 1.7;
Float_t TriFlow::mNSigmaProtonMax = 2.5;
Float_t TriFlow::mEta_Gap[4] = {0.05,0.10,0.20,0.50};

// used constant
Float_t TriFlow::mShiftOrder2[5] = {2.0, 4.0, 6.0, 8.0, 10.0};
Float_t TriFlow::mShiftOrder3[5] = {3.0, 6.0, 9.0, 12.0, 15.0};

// pt bin
//                              0 , 1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 ,10 ,11 ,12 ,13 ,14 ,15
Float_t TriFlow::pt_low[16] = {0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.2,3.6};
Float_t TriFlow::pt_up[16]  = {0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.2,3.6,4.2};

// x and y range
Float_t TriFlow::x_low[16] = {-0.4,-0.6,-0.6,-0.6,-0.6,-0.6,-0.6,-0.8,-0.8,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0};
Float_t TriFlow::x_up[16]  = { 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.6, 1.6, 2.0, 2.0, 2.0, 2.4, 2.4, 2.4, 2.4};

Float_t TriFlow::y_low[16] = {-0.15,-0.2,-0.3,-0.45,-0.6,-0.8,-1.0,-1.3,-1.5,-1.7,-1.8,-1.8,-2.0,-2.0,-2.0,-2.0};
Float_t TriFlow::y_up[16]  = { 0.15, 0.2, 0.2, 0.20, 0.3, 0.4, 0.6, 1.0, 1.0, 1.0, 1.2, 1.2, 1.4, 1.4, 1.4, 1.4};

// Centrality bin
Int_t TriFlow::cent_low[4] = {0,7,4,0}; // 0 = 0-80%, 1 = 0-10%, 2 = 10-40%, 3 = 40-80%
Int_t TriFlow::cent_up[4]  = {8,8,6,3}; // 0 = 0-80%, 1 = 0-10%, 2 = 10-40%, 3 = 40-80%
TString TriFlow::Centrality_01[4] = {"0080","0010","1040","4080"};
TString TriFlow::Centrality_23[4] = {"0070","0010","1040","4070"};

// phi-Psi bin
Double_t TriFlow::phi_Psi2_low[7] = {0.0,TMath::Pi()/14.0,2.0*TMath::Pi()/14.0,3.0*TMath::Pi()/14.0,4.0*TMath::Pi()/14.0,5.0*TMath::Pi()/14.0,6.0*TMath::Pi()/14.0};
Double_t TriFlow::phi_Psi2_up[7]  = {TMath::Pi()/14.0,2.0*TMath::Pi()/14.0,3.0*TMath::Pi()/14.0,4.0*TMath::Pi()/14.0,5.0*TMath::Pi()/14.0,6.0*TMath::Pi()/14.0,7.0*TMath::Pi()/14.0};
Double_t TriFlow::phi_Psi3_low[7] = {0.0,TMath::Pi()/21.0,2.0*TMath::Pi()/21.0,3.0*TMath::Pi()/21.0,4.0*TMath::Pi()/21.0,5.0*TMath::Pi()/21.0,6.0*TMath::Pi()/21.0};
Double_t TriFlow::phi_Psi3_up[7]  = {TMath::Pi()/21.0,2.0*TMath::Pi()/21.0,3.0*TMath::Pi()/21.0,4.0*TMath::Pi()/21.0,5.0*TMath::Pi()/21.0,6.0*TMath::Pi()/21.0,7.0*TMath::Pi()/21.0};

Double_t TriFlow::Psi2_low[3] = {-3.0*TMath::Pi()/2.0,-1.0*TMath::Pi()/2.0,1.0*TMath::Pi()/2.0};
Double_t TriFlow::Psi2_up[3]  = {-1.0*TMath::Pi()/2.0, 1.0*TMath::Pi()/2.0,3.0*TMath::Pi()/2.0};
Double_t TriFlow::Psi3_low[5] = {-4.0*TMath::Pi()/3.0,-3.0*TMath::Pi()/3.0,-1.0*TMath::Pi()/3.0,1.0*TMath::Pi()/3.0,3.0*TMath::Pi()/3.0};
Double_t TriFlow::Psi3_up[5]  = {-3.0*TMath::Pi()/3.0,-1.0*TMath::Pi()/3.0, 1.0*TMath::Pi()/3.0,3.0*TMath::Pi()/3.0,4.0*TMath::Pi()/3.0};

Int_t TriFlow::pt_total = 16;

Int_t TriFlow::Centrality_total = 4;
Int_t TriFlow::Centrality_start = 0;
Int_t TriFlow::Centrality_stop  = 1;

Int_t TriFlow::EtaGap_total = 4;
Int_t TriFlow::EtaGap_start = 0;
Int_t TriFlow::EtaGap_stop  = 2;

Int_t TriFlow::Charge_total = 2;
Int_t TriFlow::Charge_start = 0;
Int_t TriFlow::Charge_stop  = 2;

TString TriFlow::Energy[3] = {"200GeV","39GeV","27GeV"};

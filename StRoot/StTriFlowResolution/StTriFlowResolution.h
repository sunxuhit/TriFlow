#ifndef StTriFlowResolution_h
#define StTriFlowResolution_h

#include "StMessMgr.h"
#include "TString.h"
#include "TVector2.h"

class TNtuple;
class TFile;
class TProfile;
class StTriFlowHistManger;

class StTriFlowResolution
{
  public:
    StTriFlowResolution(const Int_t jobCounter, const Int_t energy);
    ~StTriFlowResolution();

    Int_t Init();

    void InitReCenterNutple(TString);
    void InitShiftCorrection(Int_t);

    // Event Plane method
    TVector2 calShiftAngle2East_EP(Int_t vz_sign, Int_t eta_gap); // x is Psi_ReCenter, y is delta_Psi
    TVector2 calShiftAngle2West_EP(Int_t vz_sign, Int_t eta_gap); // x is Psi_ReCenter, y is delta_Psi
    TVector2 calShiftAngle2A_EP(Int_t vz_sign); // x is Psi_ReCenter, y is delta_Psi
    TVector2 calShiftAngle2B_EP(Int_t vz_sign); // x is Psi_ReCenter, y is delta_Psi
    TVector2 calShiftAngle2Full_EP(Int_t vz_sign); // x is Psi_ReCenter, y is delta_Psi

    TVector2 calShiftAngle3East_EP(Int_t vz_sign, Int_t eta_gap); // x is Psi_ReCenter, y is delta_Psi
    TVector2 calShiftAngle3West_EP(Int_t vz_sign, Int_t eta_gap); // x is Psi_ReCenter, y is delta_Psi
    TVector2 calShiftAngle3A_EP(Int_t vz_sign); // x is Psi_ReCenter, y is delta_Psi
    TVector2 calShiftAngle3B_EP(Int_t vz_sign); // x is Psi_ReCenter, y is delta_Psi
    TVector2 calShiftAngle3Full_EP(Int_t vz_sign); // x is Psi_ReCenter, y is delta_Psi

    // Scalor Product method
    TVector2 calShiftAngle2East_SP(Int_t vz_sign, Int_t eta_gap); // x is Psi_ReCenter, y is delta_Psi
    TVector2 calShiftAngle2West_SP(Int_t vz_sign, Int_t eta_gap); // x is Psi_ReCenter, y is delta_Psi
    TVector2 calShiftAngle2A_SP(Int_t vz_sign); // x is Psi_ReCenter, y is delta_Psi
    TVector2 calShiftAngle2B_SP(Int_t vz_sign); // x is Psi_ReCenter, y is delta_Psi
    TVector2 calShiftAngle2Full_SP(Int_t vz_sign); // x is Psi_ReCenter, y is delta_Psi

    TVector2 calShiftAngle3East_SP(Int_t vz_sign, Int_t eta_gap); // x is Psi_ReCenter, y is delta_Psi
    TVector2 calShiftAngle3West_SP(Int_t vz_sign, Int_t eta_gap); // x is Psi_ReCenter, y is delta_Psi
    TVector2 calShiftAngle3A_SP(Int_t vz_sign); // x is Psi_ReCenter, y is delta_Psi
    TVector2 calShiftAngle3B_SP(Int_t vz_sign); // x is Psi_ReCenter, y is delta_Psi
    TVector2 calShiftAngle3Full_SP(Int_t vz_sign); // x is Psi_ReCenter, y is delta_Psi
    
    void InitResolution();
    Float_t Psi_Shift(Float_t Psi, Float_t order);
    void getResolution();
    void WriteResolution();

    Int_t Finish();

  private:
    TNtuple *mNtuple;
    Float_t mRunId, mEventId, mRefMult, mZDCx, mBBCx, mVzVpd, mCentrality9;
    Float_t mVx, mVy, mVz, mNToFMatched, mRunIndex;

    // Event Plane method
    Float_t mQ2X_East_EP[4], mQ2Y_East_EP[4],mQ2X_West_EP[4], mQ2Y_West_EP[4];
    Float_t mQ3X_East_EP[4], mQ3Y_East_EP[4],mQ3X_West_EP[4], mQ3Y_West_EP[4];
    Float_t mQ2X_Full_EP, mQ2Y_Full_EP, mQ3X_Full_EP, mQ3Y_Full_EP;
    Float_t mQ2X_A_EP, mQ2Y_A_EP, mQ2X_B_EP, mQ2Y_B_EP;
    Float_t mQ3X_A_EP, mQ3Y_A_EP, mQ3X_B_EP, mQ3Y_B_EP;

    // Scalor Product method
    Float_t mQ2X_East_SP[4], mQ2Y_East_SP[4],mQ2X_West_SP[4], mQ2Y_West_SP[4];
    Float_t mQ3X_East_SP[4], mQ3Y_East_SP[4],mQ3X_West_SP[4], mQ3Y_West_SP[4];
    Float_t mQ2X_Full_SP, mQ2Y_Full_SP, mQ3X_Full_SP, mQ3Y_Full_SP;
    Float_t mQ2X_A_SP, mQ2Y_A_SP, mQ2X_B_SP, mQ2Y_B_SP;
    Float_t mQ3X_A_SP, mQ3Y_A_SP, mQ3X_B_SP, mQ3Y_B_SP;

    // Counter
    Float_t mQCounter_East[4], mQCounter_West[4], mQCounter_Full, mQCounter_Full_East, mQCounter_Full_West;
    Float_t mQCounter_A, mQCounter_B;

    Int_t mN_Entires;

    TString mInPut_Corr_ReCenter;
    TString mOutPut_Resolution;
    TString mOutPut_EventPlane;

    TFile *mInPutFile_Shift;
    TFile *mFile_Resolution;
    TFile *mFile_EventPlane;

    static TString mVStr[2];

    TProfile *p_mRes2_EP[4]; // eta_gap
    TProfile *p_mRes3_EP[4]; // eta_gap
    TProfile *p_mRes2_SP[4]; // eta_gap
    TProfile *p_mRes3_SP[4]; // eta_gap

    TProfile *p_mRes2_ran_EP; // full event plane
    TProfile *p_mRes3_ran_EP; // full event plane
    TProfile *p_mRes2_ran_SP; // full event plane
    TProfile *p_mRes3_ran_SP; // full event plane

    StTriFlowHistManger *mTriFlowHistManger;

    Int_t mEnergy;


  ClassDef(StTriFlowResolution,1)
};
#endif

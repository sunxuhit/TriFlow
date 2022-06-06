#ifndef StTriFlowCorrection_h
#define StTriFlowCorrection_h

#include "TObject.h"
#include "TVector2.h"
#include "TString.h"

class StPicoDst;
class StPicoTrack;
class TProfile2D;
class TFile;
class TNtuple;

class StTriFlowCorrection : public TObject
{
  public:
    StTriFlowCorrection(Int_t energy);
    virtual ~StTriFlowCorrection();


    // ReCenter Correction
    bool passTrackEtaEast(StPicoTrack*,Int_t,Int_t);
    bool passTrackEtaWest(StPicoTrack*,Int_t,Int_t);
    bool passTrackFull(StPicoTrack*);

    TVector2 calq2Vector(StPicoTrack*);
    TVector2 calq3Vector(StPicoTrack*);
    Float_t getWeight(StPicoTrack*);

    void InitReCenterCorrection(Int_t);
    void addTrack_East(StPicoTrack* track, Int_t Cent9, Int_t RunIndex, Int_t i, Int_t j); // i = vz_sign, j = eta_gap
    void addTrack_West(StPicoTrack* track, Int_t Cent9, Int_t RunIndex, Int_t i, Int_t j); // i = vz_sign, j = eta_gap
    void addTrack_Full(StPicoTrack* track, Int_t Cent9, Int_t RunIndex, Int_t i); // i = vz_sign, j = eta_gap
    void addTrack_A(StPicoTrack* track, Int_t Cent9, Int_t RunIndex, Int_t i); // i = vz_sign, j = eta_gap || random sub A
    void addTrack_B(StPicoTrack* track, Int_t Cent9, Int_t RunIndex, Int_t i); // i = vz_sign, j = eta_gap || random sub B
    void Randomization();

    TVector2 getReCenterPar_East(Int_t order, Int_t Cent9, Int_t RunIndex, Int_t vz_sign, Int_t eta_gap, Int_t method); // order: 0 = 2nd, 1 = 3rd || method: 0 = EP, 1 = SP
    TVector2 getReCenterPar_West(Int_t order, Int_t Cent9, Int_t RunIndex, Int_t vz_sign, Int_t eta_gap, Int_t method); // order: 0 = 2nd, 1 = 3rd || method: 0 = EP, 1 = SP
    TVector2 getReCenterPar_Full(Int_t order, Int_t Cent9, Int_t RunIndex, Int_t vz_sign, Int_t method); // order: 0 = 2nd, 1 = 3rd || method: 0 = EP, 1 = SP

    void print(TVector2);
    void clear();

    void InitNtuple();
    Int_t fillNtuple(StPicoDst*, Int_t, Int_t, Int_t);
    void writeNtuple();

    // Shift Correction
    bool passTrackEtaNumCut(Int_t); // eta_gap
    bool passTrackFullNumCut();

    // Event Plane method
    TVector2 calPsi2_East_EP(Int_t, Int_t); // 0 = eta_gap, 1 = ShiftOrder: 2, 4, 6, 8, 10
    TVector2 calPsi2_West_EP(Int_t, Int_t); // 0 = eta_gap, 1 = ShiftOrder: 2, 4, 6, 8, 10
    TVector2 calPsi2_Full_EP(Int_t); // 1 = ShiftOrder: 2, 4, 6, 8, 10

    TVector2 calPsi3_East_EP(Int_t, Int_t); // 0 = eta_gap, 1 = ShiftOrder: 3, 6, 9, 12, 15
    TVector2 calPsi3_West_EP(Int_t, Int_t); // 0 = eta_gap, 1 = ShiftOrder: 3, 6, 9, 12, 15
    TVector2 calPsi3_Full_EP(Int_t); // 1 = ShiftOrder: 3, 6, 9, 12, 15

    // Scalor Product method
    TVector2 calPsi2_East_SP(Int_t, Int_t); // 0 = eta_gap, 1 = ShiftOrder: 2, 4, 6, 8, 10
    TVector2 calPsi2_West_SP(Int_t, Int_t); // 0 = eta_gap, 1 = ShiftOrder: 2, 4, 6, 8, 10
    TVector2 calPsi2_Full_SP(Int_t); // 1 = ShiftOrder: 2, 4, 6, 8, 10

    TVector2 calPsi3_East_SP(Int_t, Int_t); // 0 = eta_gap, 1 = ShiftOrder: 3, 6, 9, 12, 15
    TVector2 calPsi3_West_SP(Int_t, Int_t); // 0 = eta_gap, 1 = ShiftOrder: 3, 6, 9, 12, 15
    TVector2 calPsi3_Full_SP(Int_t); // 1 = ShiftOrder: 3, 6, 9, 12, 15

    void InitShiftCorrection(Int_t);
    Float_t AngleShift(Float_t Psi_raw, Float_t order);

    // Event Plane method
    Float_t calShiftAngle2East_EP(Int_t runIndex, Int_t Cent9, Int_t vz_sign, Int_t eta_gap);
    Float_t calShiftAngle2West_EP(Int_t runIndex, Int_t Cent9, Int_t vz_sign, Int_t eta_gap);
    Float_t calShiftAngle2Full_EP(Int_t runIndex, Int_t Cent9, Int_t vz_sign, StPicoTrack *track);

    Float_t calShiftAngle3East_EP(Int_t runIndex, Int_t Cent9, Int_t vz_sign, Int_t eta_gap);
    Float_t calShiftAngle3West_EP(Int_t runIndex, Int_t Cent9, Int_t vz_sign, Int_t eta_gap);
    Float_t calShiftAngle3Full_EP(Int_t runIndex, Int_t Cent9, Int_t vz_sign, StPicoTrack *track);

    // Scalor Product method
    TVector2 calQVector2East_SP(Int_t runIndex, Int_t Cent9, Int_t vz_sign, Int_t eta_gap);
    TVector2 calQVector2West_SP(Int_t runIndex, Int_t Cent9, Int_t vz_sign, Int_t eta_gap);
    TVector2 calQVector2Full_SP(Int_t runIndex, Int_t Cent9, Int_t vz_sign, StPicoTrack *track);

    TVector2 calQVector3East_SP(Int_t runIndex, Int_t Cent9, Int_t vz_sign, Int_t eta_gap);
    TVector2 calQVector3West_SP(Int_t runIndex, Int_t Cent9, Int_t vz_sign, Int_t eta_gap);
    TVector2 calQVector3Full_SP(Int_t runIndex, Int_t Cent9, Int_t vz_sign, StPicoTrack *track);

    Float_t getResolution2_EP(Int_t Cent9, Int_t eta_gap);
    Float_t getResolution3_EP(Int_t Cent9, Int_t eta_gap);
    Float_t getResolution2_SP(Int_t Cent9, Int_t eta_gap);
    Float_t getResolution3_SP(Int_t Cent9, Int_t eta_gap);

    Float_t getResolution2_Full_EP(Int_t Cent9);
    Float_t getResolution3_Full_EP(Int_t Cent9);
    Float_t getResolution2_Full_SP(Int_t Cent9);
    Float_t getResolution3_Full_SP(Int_t Cent9);

    TVector2 getQVector(Int_t j, Int_t k, Int_t l); // 0 = eta_gap, 1 = flow type, 2 = east/west

    Int_t getNumTrack(Int_t, Int_t); // 0 = eta_gap, 1 = east/west


  private:
    //Event Plane method
    TVector2 mQ2Vector_East_EP[4], mQ2Vector_West_EP[4], mQ2Vector_Full_EP, mQ2Vector_A_EP, mQ2Vector_B_EP;
    TVector2 mQ3Vector_East_EP[4], mQ3Vector_West_EP[4], mQ3Vector_Full_EP, mQ3Vector_A_EP, mQ3Vector_B_EP;

    // Scalor Product method
    TVector2 mQ2Vector_East_SP[4], mQ2Vector_West_SP[4], mQ2Vector_Full_SP, mQ2Vector_A_SP, mQ2Vector_B_SP;
    TVector2 mQ3Vector_East_SP[4], mQ3Vector_West_SP[4], mQ3Vector_Full_SP, mQ3Vector_A_SP, mQ3Vector_B_SP;

    Int_t    mQCounter_East[4], mQCounter_West[4], mQCounter_Full;
    Int_t    mQCounter_Full_East, mQCounter_Full_West;
    Int_t    mQCounter_A, mQCounter_B;
    Int_t    mEnergy;

    TFile *mInPutFile;
    TNtuple *mNtuple;
    Float_t mFillNtuple[113];

    TFile *mInPutFile_Shift;
    TFile *mInPutFile_Res;

    static TString mVStr[2];
    static TString mOrder[2];
    static TString mMethod[2];

  ClassDef(StTriFlowCorrection,1)
};

#endif

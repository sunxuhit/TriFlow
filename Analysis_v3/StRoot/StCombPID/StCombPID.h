//----------------------------------------------------------------------------------------------------
//  authors: Alexander Schmah, Hiroshi Masui
//----------------------------------------------------------------------------------------------------

#ifndef __StCombPID_h__
#define __StCombPID_h__

#include <vector>
#include "Rtypes.h"

//____________________________________________________________________________________________________
// Class to correct z-vertex dependence of refmult
class StCombPID {
  public:
    StCombPID();
    virtual ~StCombPID(); /// Default destructor

    /// Get rotated and scaled x and y values
    void     setInitValues(const Double_t P_t, const Double_t NSigma_pi, const Double_t M2, const Double_t in_ext_scale_factor); // in_ext_scale_factor = 1.0 for BES, in_ext_scale_factor = 0.8273
    Double_t getNewX() const;
    Double_t getNewY() const;

private:
    // Functions
    //Bool_t isIndexOk() const ; /// 0 <= mParameterIndex < maxArraySize

    Double_t getMeanNSigma(const Int_t PID) const;
    Double_t getWidthM2(const Int_t PID) const;
    Double_t getScaleFactor() const;
    Double_t getRotAngle(Double_t ScaleFactor) const;
    void     calcNewXY();

    Double_t mP_t;
    Double_t mNSigma_pi;
    Double_t mM2;
    Double_t mRotAngle;
    Double_t mScaleFactor;
    Double_t mNewX;
    Double_t mNewY;

    ClassDef(StCombPID, 0)
};
#endif


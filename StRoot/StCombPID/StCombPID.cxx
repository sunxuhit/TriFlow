//----------------------------------------------------------------------------------------------------
//
// Revision 1.0 aschmah
// First version of StCombPID class
//
//----------------------------------------------------------------------------------------------------

#include <fstream>
#include <iostream>
#include <string>
#include "StCombPID.h"
#include "TError.h"
#include "TRandom.h"
#include "TString.h"
#include "TMath.h"

ClassImp(StCombPID)

using std::cout ;
using std::endl ;
using std::ifstream ;
using std::string ;
using std::vector ;

static const Double_t mMeanM2_Pi      = 0.019;
static const Double_t mMeanM2_K       = 0.241;
static const Double_t mMeanM2_P       = 0.895;
static const Double_t mMeanNSigma_Pi  = 0.0;
static const Double_t mWidthNSigma_Pi = 1.0;

//____________________________________________________________________________________________________
// Default constructor
StCombPID::StCombPID()
{
    mP_t = 0.0;
    mNSigma_pi = 0.0;
    mM2 = 0.0;
    mRotAngle = 0.0;
    mScaleFactor = 0.0;
    mNewX = 0.0;
    mNewY = 0.0;
}

//____________________________________________________________________________________________________
// Default destructor
StCombPID::~StCombPID()
{
}

//____________________________________________________________________________________________________
Double_t StCombPID::getMeanNSigma(const Int_t PID) const
{
    // Returns the mean nSigma_Pi value for pions, kaons or protons for the set P_t value
    // PID: 0 = pions, 1 = kaons, 2 = protons
    Double_t x, meanNSigma, par0, par1, par2, par3;
    meanNSigma = 0.0;

    if(PID == 0) return 0.0; // pions
    if(PID == 1) // kaons
    {
        par0  = 5.84263e+00;
        par1  = -6.65912e+00;
        par2  = 8.32204e-01;
        par3  = 8.28440e-01;
    }
    if(PID == 2) // protons
    {
        par0  = 1.09685e+01;
        par1  = -8.61499e+00;
        par2  = 8.92938e-01;
        par3  = 8.08397e-01;
    }

    x = mP_t;
    if(x != 0.0) meanNSigma = par0*TMath::Power(1.0/x,par2) + par1 + par3*x;
    return meanNSigma;
}

//____________________________________________________________________________________________________
Double_t StCombPID::getWidthM2(const Int_t PID) const
{
    // Returns the width in mass2 for pions, kaons or protons for the set P_t value
    // PID: 0 = pions, 1 = kaons, 2 = protons

    Double_t x, widthM2, pol0, pol1, pol2;
    widthM2 = 0.0;
    if(PID == 0 || PID == 1 || PID == 2) // only implemented for pions so far
    {
        pol0  = 8.74975e-04;
        pol1  = -1.62659e-03;
        pol2  = 2.89828e-02;
    }

    x = mP_t;
    widthM2 = pol0 + pol1*x + pol2*x*x;
    return widthM2;
}




//____________________________________________________________________________________________________
Double_t StCombPID::getScaleFactor() const
{
    Double_t ScaleFactor = 1.0;

    Double_t WidthM2_Pi   = getWidthM2(0);

    if(WidthM2_Pi != 0) ScaleFactor = mWidthNSigma_Pi/WidthM2_Pi;

    return ScaleFactor;
}

//____________________________________________________________________________________________________
Double_t StCombPID::getRotAngle(Double_t ScaleFactor) const
{
    Double_t RotAngle = 0.0;

    Double_t MeanNSigma_K = getMeanNSigma(1);

    Double_t New_x_scale     = (MeanNSigma_K-mMeanNSigma_Pi)/(ScaleFactor); // nSigma: (kaon_mean - pion_mean) / scale_factor
    Double_t New_y_scale     = mMeanM2_K-mMeanM2_Pi;                  // m2: (kaon_mean - pion_mean)
    RotAngle = -TMath::ATan2(New_y_scale,New_x_scale); // angle between new kaon center and x-axis

    return RotAngle;
}

//____________________________________________________________________________________________________
void StCombPID::calcNewXY()
{
    Double_t NewX = 0.0;
    Double_t NewY = 0.0;

    Double_t New_x_scale     = (mNSigma_pi-mMeanNSigma_Pi)/(mScaleFactor); // nSigma: (kaon_mean - pion_mean) / scale_factor
    Double_t New_y_scale     = mM2-mMeanM2_Pi;                  // m2: (kaon_mean - pion_mean)

    NewX = TMath::Cos(mRotAngle)*New_x_scale - TMath::Sin(mRotAngle)*New_y_scale;  // rotate the scaled kaon center
    NewY = TMath::Sin(mRotAngle)*New_x_scale + TMath::Cos(mRotAngle)*New_y_scale;

    mNewX = NewX;
    mNewY = NewY;
}

//____________________________________________________________________________________________________
void StCombPID::setInitValues(const Double_t P_t, const Double_t NSigma_pi, const Double_t M2, const Double_t in_ext_scale_factor)
{
    mP_t         = P_t;
    mNSigma_pi   = NSigma_pi;
    mM2          = M2;
    mScaleFactor = getScaleFactor()/in_ext_scale_factor;
    mRotAngle    = getRotAngle(mScaleFactor);
    calcNewXY();

    //cout << "mP_t = " << mP_t << ", mNSigma_pi = " << mNSigma_pi << ", mM2 = " << mM2 << ", mScaleFactor = " << mScaleFactor << ", mRotAngle = " << mRotAngle << endl;
}

//____________________________________________________________________________________________________
Double_t StCombPID::getNewX() const
{
    return mNewX;
}

//____________________________________________________________________________________________________
Double_t StCombPID::getNewY() const
{
    return mNewY;
}


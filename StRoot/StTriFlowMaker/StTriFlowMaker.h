#ifndef StTriFlowMaker_h
#define StTriFlowMaker_h

#include "StMaker.h"
#include "TString.h"

class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;
class StRefMultCorr;
class StRunIdEventsDb;
class StCombPID;
class StTriFlowCut;
class StTriFlowProManger;
class StTriFlowCorrection;
class StTriFlowHistoManger;
class StTriFlowV0;

class StTriFlowMaker : public StMaker {
  public:
    StTriFlowMaker(const char *name, StPicoDstMaker *picoMaker, const Int_t jobCounter, const Int_t Mode, const Int_t Energy, const Int_t Flag_ME);
    virtual ~StTriFlowMaker();
    
    virtual Int_t Init();
    virtual Int_t Make();
    virtual void  Clear(Option_t *opt="");
    virtual Int_t Finish();
    
  private:
    StPicoDstMaker *mPicoDstMaker;
    StPicoDst      *mPicoDst;
    StPicoEvent *mPicoEvent;
    static StRefMultCorr *mRefMultCorr;
    StRunIdEventsDb *mRunIdEventsDb;
    StCombPID *mCombPID;
    StTriFlowCut *mTriFlowCut;
    StTriFlowProManger *mTriFlowProManger;
    StTriFlowCorrection *mTriFlowCorrection;
    StTriFlowHistoManger *mTriFlowHistoManger;
    StTriFlowV0 *mTriFlowV0;
    
    Int_t mMode;
    Int_t mEnergy;
    Int_t mFlag_ME;

    TString mInPut_Corr_ReCenter;

    TString mOutPut_ReCenterPar;
    TString mOutPut_Corr_ReCenter;
    TString mOutPut_Corr_Shift;
    TString mOutPut_ChargedFLow;
    TString mOutPut_M2_nSigPion;
    TString mOutPut_M2_Proton;
    TString mOutPut_Yields;
    TString mOutPut_Phi;
    TString mOutPut_Lambda;
    TString mOutPut_AntiLambda;

    TFile *mFile_ReCenterPar;
    TFile *mFile_Corr_ReCenter;
    TFile *mFile_Corr_Shift;
    TFile *mFile_ChargedFlow;
    TFile *mFile_M2_nSigPion;
    TFile *mFile_M2_Proton;
    TFile *mFile_Yields;
    TFile *mFile_Phi;
    TFile *mFile_Lambda;
    TFile *mFile_AntiLambda;

    Int_t mUsedTrackCounter;

    ClassDef(StTriFlowMaker, 1)
};

#endif

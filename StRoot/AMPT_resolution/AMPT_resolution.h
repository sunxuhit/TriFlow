#ifndef AMPT_resolution_h
#define AMPT_resolution_h
#include "StMessMgr.h"
#include "TString.h"

class TTree;
class TChain;
class TBranch;
class TProfile;
class TFile;
class TH1F;
class TH2F;

class AMPT_resolution // TODO: change the Tree structure
{
  public:
    AMPT_resolution(Int_t Energy, Int_t Mode, Int_t List, Long64_t StartEvent, Long64_t StopEvent); // read in energy, AMPT mode, data list, StartEvent and StopEvent
    ~AMPT_resolution();

    void SetInPutList(const TString inputlist);
    void SetOutPutFile(const TString outputfile);
    void SetStartEvent(Long64_t StartEvent);
    void SetStopEvent(Long64_t StopEvent);

    void Init(); // initialize tree structure
    void Make(); // main loop: calculate resolution with eta_sub method
    void Finish(); // save resolution 

  private:
    Int_t   mEnergy;
    Int_t mMode; // 0 for default, 1 for string melting
    TString mInPutList;
    TString mOutPutFile;
    Long64_t mStartEvent;
    Long64_t mStopEvent;
    static Int_t mInput_flag;
    static TString mBeamEnergy[7];
    static TString mMode_AMPT[2];
    static Int_t mCentrality[2][7][10]; // centrality definition
    static Int_t mList_start[10];
    static Int_t mList_stop[10];
    TProfile *p_mRes2;
    TProfile *p_mRes3;
    TFile *mFile_OutPut;

    // QA Plot
    TH1F *h_mPart;
    TH1F *h_mMult;
    TH1F *h_mRefMult;
    TH1F *h_mEta;
    TH1F *h_mPsi2_East;
    TH1F *h_mPsi2_West;
    TH1F *h_mPsi3_East;
    TH1F *h_mPsi3_West;
    TH2F *h_mPsi2;
    TH2F *h_mPsi3;

    //---------------------------------------------------------------
    TChain         *mChain_Input;
    Int_t           fCurrent; //!current Tree number in a TChain

    // Declaration of leaf types
    Int_t           Event;
    Int_t           Mult;
    Int_t           Npartp;
    Int_t           Npartt;
    Int_t           Nesp;
    Int_t           Ninesp;
    Int_t           Nest;
    Int_t           Ninest;
    Float_t         Imp;
    Int_t           Na;
    Int_t           Nb;
    Int_t           Nab;
    Float_t         Psi;
    Float_t         Nx[394];   //[Nab]
    Float_t         Ny[394];   //[Nab]
    Float_t         Nz[394];   //[Nab]
    Int_t           Stat[394];   //[Nab]
    Int_t           PID[38570];   //[Mult]
    Float_t         Px[38570];   //[Mult]
    Float_t         Py[38570];   //[Mult]
    Float_t         Pz[38570];   //[Mult]
    Float_t         Mass[38570];   //[Mult]
    Float_t         XX[38570];   //[Mult]
    Float_t         YY[38570];   //[Mult]
    Float_t         ZZ[38570];   //[Mult]
    Float_t         TT[38570];   //[Mult]

   // List of branches
    TBranch        *b_Event;   //!
    TBranch        *b_Mult;   //!
    TBranch        *b_Npartp;   //!
    TBranch        *b_Npartt;   //!
    TBranch        *b_Nesp;   //!
    TBranch        *b_Ninesp;   //!
    TBranch        *b_Nest;   //!
    TBranch        *b_Ninest;   //!
    TBranch        *b_Imp;   //!
    TBranch        *b_Na;   //!
    TBranch        *b_Nb;   //!
    TBranch        *b_Nab;   //!
    TBranch        *b_Psi;   //!
    TBranch        *b_Nx;   //!
    TBranch        *b_Ny;   //!
    TBranch        *b_Nz;   //!
    TBranch        *b_Stat;   //!
    TBranch        *b_PID;   //!
    TBranch        *b_Px;   //!
    TBranch        *b_Py;   //!
    TBranch        *b_Pz;   //!
    TBranch        *b_Mass;   //!
    TBranch        *b_XX;   //!
    TBranch        *b_YY;   //!
    TBranch        *b_ZZ;   //!
    TBranch        *b_TT;   //!

  ClassDef(AMPT_resolution,1)
};
#endif

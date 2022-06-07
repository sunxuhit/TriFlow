#ifndef AMPT_epsilon_h
#define AMPT_epsilon_h
#include "StMessMgr.h"
#include "TString.h"

class TTree;
class TChain;
class TBranch;
class TProfile;
class TFile;
class TH1D;

class AMPT_epsilon // TODO: change the Tree structure
{
  public:
    AMPT_epsilon(Int_t Energy, Int_t Mode, Int_t Screen, Int_t List, Long64_t StartEvent, Long64_t StopEvent); // read in energy, AMPT mode, data list, StartEvent and StopEvent
    ~AMPT_epsilon();

    void SetInPutList(const TString inputlist);
    void SetOutPutFile(const TString outputfile);
    void SetInPutRes(const TString inputres);
    void SetStartEvent(Long64_t StartEvent);
    void SetStopEvent(Long64_t StopEvent);

    void Init(); // initialize tree structure
    void Make(); // main loop: calculate resolution with eta_sub method
    void Finish(); // save resolution 

    Float_t getResolution(Int_t order, Int_t i_cent); // get Resolution 
    Int_t getCentrality(Int_t refMult); // get Centrality

  private:
    Int_t   mEnergy;
    Int_t mMode; // 0 for default, 1 for string melting
    Int_t mScreen; // 0 for 3mb, 1 for 6mb
    TString mInPutList;
    TString mOutPutFile;
    TString mInPutRes;
    Long64_t mStartEvent;
    Long64_t mStopEvent;
    static Int_t mInput_flag;
    static TString mBeamEnergy[7];
    static TString mMode_AMPT[2];
    static TString mScreenMass_AMPT[3];
    static Int_t mRefMult[2][7][10]; // centrality definition
    static Int_t mList_start[25];
    static Int_t mList_stop[25];
    static Int_t cent_low[4];
    static Int_t cent_up[4];
    static Int_t Centrality_start;
    static Int_t Centrality_stop;
    TFile *mFile_OutPut;
    TFile *mFile_Res;

    // resolution
    TProfile *p_mRes[2];

    // dN/dy|[-0.5,0.5]
    TH1D *h_mRapNarrow[9]; // 9 narrow centrality bins
    TH1D *h_mRapWide[4]; // 4 wide centrality bins

    // Transverse overlap area
    TProfile *p_mAreaT9;
    TProfile *p_mAreaT4;

    // Epsilon
    TProfile *p_mEpsilon9[2]; // 0 for 2nd, 1 for 3rd
    TProfile *p_mEpsilon4[2]; // 0 for 2nd, 1 for 3rd

    // Event Counter
    TH1D *h_mEventCounter9;
    TH1D *h_mEventCounter4;

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

  ClassDef(AMPT_epsilon,1)
};
#endif

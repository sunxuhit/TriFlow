#ifndef AMPT_phi_h
#define AMPT_phi_h
#include "StMessMgr.h"
#include "TString.h"
#include <vector>
#include "TLorentzVector.h"
#include "TVector2.h"

class TTree;
class TChain;
class TBranch;
class TProfile;
class TFile;
class TH1F;

class AMPT_phi
{
  public:
    AMPT_phi(Int_t Energy, Int_t Mode, Int_t List, Long64_t StartEvent, Long64_t StopEvent, Int_t Flag_ME); // read in energy, AMPT mode, data list, StartEvent and StopEvent
    ~AMPT_phi();

    void SetInPutList(const TString inputlist);
    void SetOutPutFile(const TString outputfile);
    void SetInPutRes(const TString inputres);
    void SetStartEvent(Long64_t StartEvent);
    void SetStopEvent(Long64_t StopEvent);

    void Init(); // initialize tree structure
    void Make(); // main loop: calculate resolution with eta_sub method
    void FillHist2nd(Float_t pt, Int_t cent9, Float_t phi, Float_t res, Float_t InvMass); // fill Histogram for flow calculation
    void FillHist3rd(Float_t pt, Int_t cent9, Float_t phi, Float_t res, Float_t InvMass); // fill Histogram for flow calculation
    void Finish(); // save resolution 
    void clear_phi(Int_t cent9); // clear everything used for phi reconstruction
    void doPhi(Int_t cent9); // reconstruct phi meson

    Float_t getResolution(Int_t order, Int_t i_cent); // get Resolution 
    Int_t getCentrality(Int_t refMult); // get Centrality

  private:
    Int_t mEnergy;
    Int_t mMode; // 0 for default, 1 for string melting
    Int_t mFlag_ME;
    TString mInPutList;
    TString mOutPutFile;
    TString mInPutRes;
    Long64_t mStartEvent;
    Long64_t mStopEvent;
    static Int_t mInput_flag;
    static TString mBeamEnergy[7];
    static TString mMode_AMPT[2];
    static Int_t mRefMult[2][7][10]; // centrality definition
    static Int_t mList_start[20];
    static Int_t mList_stop[20];
    static Int_t cent_low[4];
    static Int_t cent_up[4];
    static Int_t Centrality_start;
    static Int_t Centrality_stop;
    static Float_t mMassKaon;
    TFile *mFile_OutPut;
    TFile *mFile_Res;


    // pt bin
    static Float_t pt_low_phi[23];
    static Float_t pt_up_phi[23];

    // phi-Psi bin
    static Float_t phi_Psi2_low[7];
    static Float_t phi_Psi2_up[7];
    static Float_t phi_Psi3_low[7];
    static Float_t phi_Psi3_up[7];

    static Float_t Psi2_low[3];
    static Float_t Psi2_up[3];
    static Float_t Psi3_low[5];
    static Float_t Psi3_up[5];

    static Int_t pt_total_phi;
    static Int_t Phi_Psi_total;

    // mixed event
    static Int_t Buffer_depth;
    Int_t mEventCounter[9]; // centrality bin
    std::vector<TVector2> mQ2East[9]; // 0 = centrality bin | push_back->event
    std::vector<TVector2> mQ2West[9];
    std::vector<TVector2> mQ3East[9];
    std::vector<TVector2> mQ3West[9];
//    std::vector<Int_t> mRefMult[9];
    std::vector<Int_t> mCentrality[9];
    // store daughter particles of phi
    std::vector<TLorentzVector> mKplus[9][5]; // 0 = centrality bin, 1 = event bin | push_back->track(4-Vector)
    std::vector<TLorentzVector> mKminus[9][5];
    std::vector<Int_t> mFlag_Kplus[9][5]; // 0 = centrality bin, 1 = event bin | push_back->track(mEventCounter)
    std::vector<Int_t> mFlag_Kminus[9][5];

    // QA Plot
    TH1F *h_mPart;
    TH1F *h_mMult;
    TH1F *h_mRefMult;
    TH1F *h_mEta;
    TH1F *h_mPsi2_East;
    TH1F *h_mPsi2_West;
    TH1F *h_mPsi3_East;
    TH1F *h_mPsi3_West;
    TH1F *h_mCentrality;
    // invariant mass distribution for resolution correction
    TH1F *h_mPhi[9]; 

    // resolution
    TProfile *p_mRes[2];

    // flow for phi by using eta_sub event plane method
    TH1F *h_mFlow_phi[2][4][23][7]; // 0 for 2nd, 1 for 3rd | 0 for 0-80%, 1 for 0-10%, 2 for 10-40%, 3 for 40-80% | pt bin | phi-Psi bin

    // pt spectra for phi
    TH1F *h_mPt_phi[2][4][23]; // 0 for 2nd, 1 for 3rd | 0 for 0-80%, 1 for 0-10%, 2 for 10-40%, 3 for 40-80% | pt bin

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

  ClassDef(AMPT_phi,1)
};
#endif

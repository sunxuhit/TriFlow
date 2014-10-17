#ifndef __STPICOALEXEVENT_H__
#define __STPICOALEXEVENT_H__

#include "TMath.h"
#include "TObject.h"
#include "TClonesArray.h"
#include "StThreeVectorF.hh"
#include "TVector2.h"

// A. Schmah 12.10.2011

class StPicoAlexTrack : public TObject
{
private:
    // Track properties
    Int_t    mId;               // track Id
    Float_t  mChi2;             // chi2*1000
    Float_t  mChi2Prob;         // chi2prob*1000
    StThreeVectorF mGMomentum;  // Global momentum
    StThreeVectorF mPMomentum;  // primary momentum, (0.,0.,0.) if none
    Int_t    mFlowFlag;         // 1 - tpc EP, 2 - ftpc EP, 0 - none
    Float_t  mQXi;              //
    Float_t  mQYi;              // Q-vector for this track
    Float_t  mOriginX;          // global helix origin X * 100
    Float_t  mOriginY;          // global helix origin Y * 100
    Float_t  mOriginZ;          // global helix origin Z * 100
    Float_t  mGDca;             // global dca*1000
    Float_t  mDedx;             // dEdx*1000
    Int_t    mNHitsFit;         // q*nHitsFit
    Int_t    mNHitsMax;         // nHitsMax
    Int_t    mNHitsDedx;        // nHitsDedx
    Float_t  mNSigmaPion;       // nsigmaPi * 100
    Float_t  mNSigmaKaon;       // nsigmaK * 100
    Float_t  mNSigmaProton;     // nsigmaP * 100
    Float_t  mNSigmaElectron;   // nsigmaE * 100

    // BTOF variables
    Int_t    mBTofCellId;       // (tray-1)*192+(module-1)*6+(cell-1): -1 - no match
    Int_t    mBTofMatchFlag;    // 0 - no match, 1 - one-to-one, 2 - one-to-multiple
    Float_t  mBTof;             // time-Of-Flight * 1000 in ns
    Float_t  mBTofBeta;         // beta * 20000
    Float_t  mBTofYLocal;       // ylocal * 1000
    Float_t  mBTofZLocal;       // zlocal * 1000
    Float_t  mBTofHitPosX;      // projected hit position X * 100
    Float_t  mBTofHitPosY;      // projected hit position Y * 100
    Float_t  mBTofHitPosZ;      // projected hit position Z * 100

    // these variables are extracted from the standard BEMC cluster algorithm
    Int_t    mBEMCId;           // index in bemcPoint array
    Int_t    mBTOWADC0;         // adc0 higest adc in the cluster
    Float_t  mBTOWE0;           // E0*1000 highest tower in the cluster
    Float_t  mBTOWE;            // EMC point E*1000
    Float_t  mBEMCDistZ;        // z*100
    Float_t  mBEMCDistPhi;      // phi*10000
    Int_t    mBSMDNEta;         // # of hits in eta
    Int_t    mBSMDNPhi;         // # of hits in phi

    // these variables are purely from single tower or nearby towers
    Int_t    mBTOWId;           // projected tower Id 1-4800
    Int_t    mBTOWId23;         // emc 2nd and 3rd closest tower local id  ( 2nd X 10 + 3rd), each id 0-8
    Float_t  mBTOWE1;           // E1*1000 matched (closest) tower E
    Float_t  mBTOWE2;           // E2*1000 2nd closest tower E
    Float_t  mBTOWE3;           // E3*1000 3rd closest tower E
    Float_t  mBTOWDistEta;      // eta*10000 distance between track and matched tower center
    Float_t  mBTOWDistPhi;      // phi*10000 distance between track and matched tower center

public:
    StPicoAlexTrack() :

        mId(0),mChi2(0),mChi2Prob(0),mGMomentum(0,0,0),mPMomentum(0,0,0),mFlowFlag(0),mQXi(0),mQYi(0),mOriginX(0),mOriginY(0),mOriginZ(0),mGDca(0),
        mDedx(0),mNHitsFit(0),mNHitsMax(0),mNHitsDedx(0),mNSigmaPion(0),mNSigmaKaon(0),mNSigmaProton(0),mNSigmaElectron(0),mBTofCellId(0),
        mBTofMatchFlag(0),mBTof(0),mBTofBeta(0),mBTofYLocal(0),mBTofZLocal(0),mBTofHitPosX(0),mBTofHitPosY(0),mBTofHitPosZ(0),mBEMCId(0),
        mBTOWADC0(0),mBTOWE0(0),mBTOWE(0),mBEMCDistZ(0),mBEMCDistPhi(0),mBSMDNEta(0),mBSMDNPhi(0),mBTOWId(0),mBTOWId23(0),mBTOWE1(0),mBTOWE2(0),
        mBTOWE3(0),mBTOWDistEta(0),mBTOWDistPhi(0)
    {
    }
        ~StPicoAlexTrack() {}

        // getters
        Int_t          id() const             { return (Int_t)mId; }
        Float_t        chi2() const           { return (Float_t)mChi2; }
        Float_t        chi2prob() const       { return (Float_t)mChi2Prob;}
        StThreeVectorF gMom() const           { return mGMomentum; }
        StThreeVectorF pMom() const           { return mPMomentum; }
        StThreeVectorF origin() const         { return StThreeVectorF(mOriginX,mOriginY,mOriginZ); }
        Int_t          flowFlag() const       { return (Int_t)mFlowFlag; }
        TVector2       Qi() const             { return TVector2(mQXi, mQYi); }
        Float_t        dca() const            { return (Float_t)mGDca; }
        Short_t        charge() const         { return (mNHitsFit>0) ? +1 : -1; }
        Int_t          nHitsFit() const       { return (mNHitsFit>0) ? (Int_t)mNHitsFit : (Int_t)(-1*mNHitsFit); }
        Int_t          nHitsMax() const       { return (Int_t)mNHitsMax; }
        Int_t          nHitsDedx() const      { return (Int_t)mNHitsDedx; }
        Float_t        dEdx() const           { return (Float_t)mDedx; }
        Float_t        nSigmaPion() const     { return (Float_t)mNSigmaPion; }
        Float_t        nSigmaKaon() const     { return (Float_t)mNSigmaKaon; }
        Float_t        nSigmaProton() const   { return (Float_t)mNSigmaProton; }
        Float_t        nSigmaElectron() const { return (Float_t)mNSigmaElectron; }
        Int_t          btofCellId() const     { return (Int_t)mBTofCellId; }
        Int_t          btofMatchFlag() const  { return (Int_t)mBTofMatchFlag; }
        Float_t        btof() const           { return (Float_t)mBTof; }
        Float_t        btofBeta() const       { return (Float_t)mBTofBeta; }
        Float_t        btofYLocal() const     { return (Float_t)mBTofYLocal; }
        Float_t        btofZLocal() const     { return (Float_t)mBTofZLocal; }
        StThreeVectorF btofHisPos() const     { return StThreeVectorF(mBTofHitPosX, mBTofHitPosY, mBTofHitPosZ); }
        Int_t          bemcId() const         { return (Int_t)mBEMCId; }
        Int_t          adc0() const           { return (Int_t)mBTOWADC0; }
        Float_t        e0() const             { return (Float_t)mBTOWE0; }
        Float_t        e() const              { return (Float_t)mBTOWE; }
        Float_t        zDist() const          { return (Float_t)mBEMCDistZ; }
        Float_t        phiDist() const        { return (Float_t)mBEMCDistPhi; }
        Int_t          nEta() const           { return (Int_t)mBSMDNEta; }
        Int_t          nPhi() const           { return (Int_t)mBSMDNPhi; }
        Int_t          btowId() const         { return (Int_t)mBTOWId; }
        Int_t          btowId2() const        { return (Int_t)mBTOWId23; }
        Int_t          btowId3() const        { return (Int_t)mBTOWId23; }
        Float_t        e1() const             { return (Float_t)mBTOWE1; }
        Float_t        e2() const             { return (Float_t)mBTOWE2; }
        Float_t        e3() const             { return (Float_t)mBTOWE3; }
        Float_t        etaTowDist() const     { return (Float_t)mBTOWDistEta; }
        Float_t        phiTowDist() const     { return (Float_t)mBTOWDistPhi; }


        // setters
        void setid(Int_t i)                    { mId = i; }
        void setchi2(Float_t f)                { mChi2 = f; }
        void setchi2prob(Float_t f)            { mChi2Prob = f;}
        void setgMom(StThreeVectorF thv)       { mGMomentum = thv; }
        void setpMom(StThreeVectorF thv)       { mPMomentum = thv; }
        void setorigin(Float_t fx, Float_t fy, Float_t fz) {mOriginX = fx; mOriginY = fy; mOriginZ = fz; }
        void setflowFlag(Int_t i)            { mFlowFlag = i; }
        void setQi( Float_t fx, Float_t fy)    { mQXi = fx; mQYi = fy; }
        void setdca(Float_t f)                 { mGDca = f; }
        void setnHitsFit(Int_t i)            { mNHitsFit = i; }
        void setnHitsMax(Int_t i)            { mNHitsMax = i; }
        void setnHitsDedx(Int_t i)           { mNHitsDedx = i; }
        void setdEdx(Float_t f)                { mDedx = f; }
        void setnSigmaPion(Float_t f)          { mNSigmaPion = f; }
        void setnSigmaKaon(Float_t f)          { mNSigmaKaon = f; }
        void setnSigmaProton(Float_t f)        { mNSigmaProton = f; }
        void setnSigmaElectron(Float_t f)      { mNSigmaElectron = f; }
        void setbtofCellId(Int_t i)          { mBTofCellId = i; }
        void setbtofMatchFlag(Int_t i)       { mBTofMatchFlag = i; }
        void setbtof(Float_t f)                { mBTof = f; }
        void setbtofBeta(Float_t f)            { mBTofBeta = f; }
        void setbtofYLocal(Float_t f)          { mBTofYLocal = f; }
        void setbtofZLocal(Float_t f)          { mBTofZLocal = f; }
        void setbtofHisPos(Float_t fx, Float_t fy, Float_t fz) { mBTofHitPosX = fx; mBTofHitPosY = fy; mBTofHitPosZ = fz; }
        void setbemcId(Int_t i)              { mBEMCId = i; }
        void setadc0(Int_t i)                { mBTOWADC0 = i; }
        void sete0(Float_t f)                  { mBTOWE0 = f; }
        void sete(Float_t f)                   { mBTOWE = f; }
        void setzDist(Float_t f)               { mBEMCDistZ = f; }
        void setphiDist(Float_t f)             { mBEMCDistPhi = f; }
        void setnEta(Int_t i)                { mBSMDNEta = i; }
        void setnPhi(Int_t i)                { mBSMDNPhi = i; }
        void setbtowId(Int_t i)              { mBTOWId = i; }
        void setbtowId2(Int_t i)             { mBTOWId23 = i; }
        void setbtowId3(Int_t i)             { mBTOWId23 = i; }
        void sete1(Float_t f)                  { mBTOWE1 = f; }
        void sete2(Float_t f)                  { mBTOWE2 = f; }
        void sete3(Float_t f)                  { mBTOWE3 = f; }
        void setetaTowDist(Float_t f)          { mBTOWDistEta = f; }
        void setphiTowDist(Float_t f)          { mBTOWDistPhi = f; }

        ClassDef(StPicoAlexTrack,1)  // A simple track of a particle
};

class StPicoAlexEvent : public TObject
{
private:
    Int_t          mRunId;           // run number
    Int_t          mEventId;         // event number
    Int_t          mFillId;          // fill number
    Float_t        mBField;          // B field in kilogauss
    StThreeVectorF mPrimaryVertex;   // primary Vertex (1st)
    StThreeVectorF mSecondVertex;    // second Vertex position (for study)
    Int_t          mTriggerWord;     // self-defined trigger word - see code for details
    Int_t          mRefMultFtpcEast; // FTPC refMult east
    Int_t          mRefMultFtpcWest; // FTPC refMult west
    Int_t          mRefMultNeg;      // TPC refMult neg
    Int_t          mRefMultPos;      // TPC refMult pos

    Int_t          mNVpdHitsEast;    // Vpd Hits east;
    Int_t          mNVpdHitsWest;    // vpd hits west;
    Int_t          mNT0;             // number of T0 particles in BTOF self calibration
    Float_t        mVzVpd;           // VzVpd*100.

    Float_t        mZDCx;           // zdcX
    Float_t        mBBCx;
    Float_t mBackgroundRate;
    Float_t mBbcBlueBackgroundRate;
    Float_t mBbcYellowBackgroundRate;
    Float_t mBbcEastRate;
    Float_t mBbcWestRate;
    Float_t mZdcEastRate;
    Float_t mZdcWestRate;
    //Nov.10, 2008, Na
    Float_t mVpd[64];
    Float_t mZdcSumAdcEast;
    Float_t mZdcSumAdcWest;
    Float_t mZdcSmdEastHorizontal[8];
    Float_t mZdcSmdEastVertical[8];
    Float_t mZdcSmdWestHorizontal[8];
    Float_t mZdcSmdWestVertical[8];
    Float_t mSpaceCharge;

    UShort_t mbTofTrayMultiplicity; // BTOF tray multiplicity
    UShort_t mNumberOfGlobalTracks; // # of global tracks

    // From StMuPrimaryVertex
    Float_t mRanking;
    UShort_t mNBEMCMatch;

    // BBC ADC for q-vectors (Hiroshi)
    UShort_t mBbcAdcEast[24]; /// BBC East ADC: 0-23
    UShort_t mBbcAdcWest[24]; /// BBC West ADC: 24-47

    Int_t mYear;
    Int_t mDay;
    Float_t mEnergy;
    Bool_t mIsMinBias;
    Bool_t mIsMBSlow;
    Bool_t mIsCentral;
    Bool_t mIsHT;
    Bool_t mIsHT11;
    Bool_t mIsHT15;

    UShort_t      fNumTracks;
    TClonesArray* fTracks;      //->

public:
    StPicoAlexEvent() :
        mRunId(0),mEventId(0),mFillId(0),mBField(0),mPrimaryVertex(0),mSecondVertex(0),mTriggerWord(0),mRefMultFtpcEast(0),mRefMultFtpcWest(0),
        mRefMultNeg(0),mRefMultPos(0),mNVpdHitsEast(0),mNVpdHitsWest(0),mNT0(0),mVzVpd(0),mZDCx(0),mBBCx(0),mBackgroundRate(0),
        mBbcBlueBackgroundRate(0),mBbcYellowBackgroundRate(0),mBbcEastRate(0),mBbcWestRate(0),mZdcEastRate(0),mZdcWestRate(0),
        mVpd(),
        mZdcSumAdcEast(0),mZdcSumAdcWest(0),mZdcSmdEastHorizontal(),mZdcSmdEastVertical(),mZdcSmdWestHorizontal(),mZdcSmdWestVertical(),
        mSpaceCharge(0),mbTofTrayMultiplicity(0),mNumberOfGlobalTracks(0),mRanking(0),mNBEMCMatch(0),mBbcAdcEast(),mBbcAdcWest(),
        mYear(0),mDay(0),mEnergy(0),mIsMinBias(0),mIsMBSlow(0),mIsCentral(0),mIsHT(0),mIsHT11(0),mIsHT15(0)
    {
        fTracks      = new TClonesArray( "StPicoAlexTrack", 10 );
    }
        ~StPicoAlexEvent()
        {
            delete fTracks;
            fTracks = NULL;
        }


        // getters
        Int_t          runId() const                       { return mRunId; }
        Int_t          eventId() const                     { return mEventId; }
        Int_t          fillId() const                      { return (Int_t)mFillId; }
        Float_t        bField() const                      { return mBField; }
        StThreeVectorF primaryVertex() const               { return mPrimaryVertex; }
        Int_t          triggerWord() const                 { return mTriggerWord; }
        Int_t          refMultPos() const                  { return (Int_t)mRefMultPos; }
        Int_t          refMultNeg() const                  { return (Int_t)mRefMultNeg; }
        Int_t          refMultFtpcEast() const             { return (Int_t)mRefMultFtpcEast; }
        Int_t          refMultFtpcWest() const             { return (Int_t)mRefMultFtpcWest; }
        Int_t          refMult() const                     { return (Int_t)(mRefMultPos+mRefMultNeg); }
        Int_t          refMultFtpc() const                 { return (Int_t)(mRefMultFtpcEast+mRefMultFtpcWest); }
        Int_t          nVpdHitsEast() const                { return (Int_t)mNVpdHitsEast; }
        Int_t          nVpdHitsWest() const                { return (Int_t)mNVpdHitsWest; }
        Int_t          nT0() const                         { return (Int_t)mNT0; }
        Float_t        vzVpd() const                       { return (Float_t)mVzVpd; }
        Float_t        ZDCx() const                        { return mZDCx; }
        Float_t        BBCx() const                        { return mBBCx; }
        Float_t        Vpd(Int_t i) const                  { return (Float_t)mVpd[i]; }
        Float_t        ZdcSumAdcEast() const               { return (Float_t)mZdcSumAdcEast; }
        Float_t        ZdcSumAdcWest() const               { return (Float_t)mZdcSumAdcWest; }
        Float_t        ZdcSmdEastHorizontal(Int_t i) const { return (Float_t)mZdcSmdEastHorizontal[i]; }
        Float_t        ZdcSmdEastVertical(Int_t i) const   { return (Float_t)mZdcSmdEastVertical[i]; }
        Float_t        ZdcSmdWestHorizontal(Int_t i) const { return (Float_t)mZdcSmdWestHorizontal[i]; }
        Float_t        ZdcSmdWestVertical(Int_t i) const   { return (Float_t)mZdcSmdWestVertical[i]; }
        Float_t        backgroundRate() const              { return mBackgroundRate; }
        Float_t        bbcBlueBackgroundRate() const       { return mBbcBlueBackgroundRate; }
        Float_t        bbcYellowBackgroundRate() const     { return mBbcYellowBackgroundRate; }
        Float_t        bbcEastRate() const                 { return mBbcEastRate; }
        Float_t        bbcWestRate() const                 { return mBbcWestRate; }
        Float_t        zdcEastRate() const                 { return mZdcEastRate; }
        Float_t        zdcWestRate() const                 { return mZdcWestRate; }
        Float_t        spaceCharge() const                 { return mSpaceCharge; }
        UShort_t       btofTrayMultiplicity() const        { return mbTofTrayMultiplicity; }
        UShort_t       numberOfGlobalTracks() const        { return mNumberOfGlobalTracks; }
        Float_t        ranking() const                     { return mRanking; }
        UShort_t       nBEMCMatch() const                  { return mNBEMCMatch ; }
        UShort_t       bbcAdcEast(const Int_t i)           { return mBbcAdcEast[i]; }
        UShort_t       bbcAdcWest(const Int_t i)           { return mBbcAdcWest[i]; }

        // other user's functions
        Int_t          year() const                        { return mYear; }
        Int_t          day() const                         { return mDay; }
        Float_t        energy() const                      { return mEnergy; }
        Bool_t         isMinBias() const                   { return mIsMinBias; }
        Bool_t         isMBSlow() const                    { return mIsMBSlow; }
        Bool_t         isCentral() const                   { return mIsCentral; }
        Bool_t         isHT() const                        { return mIsHT; }
        Bool_t         isHT11() const                      { return mIsHT11; }
        Bool_t         isHT15() const                      { return mIsHT15; }



        // setters
        void setrunId(Int_t i)                           { mRunId = i; }
        void seteventId(Int_t i)                         { mEventId = i; }
        void setfillId(Int_t i)                          { mFillId = i; }
        void setbField(Float_t f)                        { mBField = f; }
        void setprimaryVertex(StThreeVectorF thv)        { mPrimaryVertex = thv; }
        void settriggerWord(Int_t i)                     { mTriggerWord = i; }
        void setrefMultPos(Int_t i)                      { mRefMultPos = i; }
        void setrefMultNeg(Int_t i)                      { mRefMultNeg = i; }
        void setrefMultFtpcEast(Int_t i)                 { mRefMultFtpcEast = i; }
        void setrefMultFtpcWest(Int_t i)                 { mRefMultFtpcWest = i; }
        void setnVpdHitsEast(Int_t i)                    { mNVpdHitsEast = i; }
        void setnVpdHitsWest(Int_t i)                    { mNVpdHitsWest = i; }
        void setnT0(Int_t i)                             { mNT0 = i; }
        void setvzVpd(Float_t f)                         { mVzVpd = f; }
        void setZDCx(Float_t f)                          { mZDCx = f; }
        void setBBCx(Float_t f)                          { mBBCx = f; }
        void setVpd(Int_t i, Float_t f)                  { mVpd[i] = f; }
        void setZdcSumAdcEast(Float_t f)                 { mZdcSumAdcEast = f; }
        void setZdcSumAdcWest(Float_t f)                 { mZdcSumAdcWest = f; }
        void setZdcSmdEastHorizontal(Int_t i, Float_t f) { mZdcSmdEastHorizontal[i] = f; }
        void setZdcSmdEastVertical(Int_t i, Float_t f)   { mZdcSmdEastVertical[i] = f; }
        void setZdcSmdWestHorizontal(Int_t i, Float_t f) { mZdcSmdWestHorizontal[i] = f; }
        void setZdcSmdWestVertical(Int_t i, Float_t f)   { mZdcSmdWestVertical[i] = f; }
        void setbackgroundRate(Float_t f)                { mBackgroundRate = f; }
        void setbbcBlueBackgroundRate(Float_t f)         { mBbcBlueBackgroundRate = f; }
        void setbbcYellowBackgroundRate(Float_t f)       { mBbcYellowBackgroundRate = f; }
        void setbbcEastRate(Float_t f)                   { mBbcEastRate = f; }
        void setbbcWestRate(Float_t f)                   { mBbcWestRate = f; }
        void setzdcEastRate(Float_t f)                   { mZdcEastRate = f; }
        void setzdcWestRate(Float_t f)                   { mZdcWestRate = f; }
        void setspaceCharge(Float_t f)                   { mSpaceCharge = f; }
        void setbtofTrayMultiplicity(UShort_t us)        { mbTofTrayMultiplicity = us; }
        void setnumberOfGlobalTracks(UShort_t us)        { mNumberOfGlobalTracks = us; }
        void setranking(Float_t f)                       { mRanking = f; }
        void setnBEMCMatch(UShort_t us)                  { mNBEMCMatch  = us; }
        void setbbcAdcEast(Int_t i, UShort_t us)         { mBbcAdcEast[i] = us; }
        void setbbcAdcWest(Int_t i, UShort_t us)         { mBbcAdcWest[i] = us; }

        // other user's functions
        void setyear(Int_t i)                            { mYear = i; }
        void setday(Int_t i)                             { mDay = i; }
        void setenergy(Float_t f)                        { mEnergy = f; }
        void setisMinBias(Bool_t b)                      { mIsMinBias = b; }
        void setisMBSlow(Bool_t b)                       { mIsMBSlow = b; }
        void setisCentral(Bool_t b)                      { mIsCentral = b; }
        void setisHT(Bool_t b)                           { mIsHT = b; }
        void setisHT11(Bool_t b)                         { mIsHT11 = b; }
        void setisHT15(Bool_t b)                         { mIsHT15 = b; }


        StPicoAlexTrack* createTrack()
        {
            if (fNumTracks == fTracks->GetSize())
                fTracks->Expand( fNumTracks + 10 );
            if (fNumTracks >= 10000)
            {
                Fatal( "StPicoAlexEvent::createTrack()", "ERROR: Too many tracks (>10000)!" );
                exit( 2 );
            }

            new((*fTracks)[fNumTracks++]) StPicoAlexTrack;
            return (StPicoAlexTrack*)((*fTracks)[fNumTracks - 1]);
        }
        void clearTrackList()
        {
            fNumTracks   = 0;
            fTracks      ->Clear();
        }
        UShort_t numberOfTracks() const
        {
            return fNumTracks;
        }
        StPicoAlexTrack* track(UShort_t i) const
        {
            return i < fNumTracks ? (StPicoAlexTrack*)((*fTracks)[i]) : NULL;
        }

        ClassDef(StPicoAlexEvent,1)  // A simple event compiled of tracks
};


#endif // __STPICOALEXEVENT_H__

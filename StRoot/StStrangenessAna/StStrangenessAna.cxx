#include "StStrangenessAna.h"
#include "StStrangenessCons.h"
#include "StStrangenessCut.h"
#include "StStrangenessCorr.h"
#include "StStrangenessHistoManger.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"
#include "StRoot/StAlexPhiMesonEvent/StAlexPhiMesonEvent.h"
#include "StRoot/StV0Event/StV0Event.h"
#include "StRoot/StRunIdEventsDb/StRunIdEventsDb.h"
#include "StThreeVectorF.hh"
#include "StMessMgr.h"
#include "TFile.h"
#include "TChain.h"
#include "TVector2.h"
#include "TLorentzVector.h"
#include <fstream>

ClassImp(StStrangenessAna)

StRefMultCorr* StStrangenessAna::mRefMultCorr = NULL;
Int_t StStrangenessAna::mSE_input_flag = 1;
Int_t StStrangenessAna::mME_input_flag = 1;
char* StStrangenessAna::XUV0_EVENT_TREE = NULL;
char* StStrangenessAna::XUV0_EVENT_BRANCH = NULL;

// Systematic Errors estimation
// phi meson
Float_t StStrangenessAna::pt_add[20]        = {0.0,0.3,-0.2,0.2,0.2,0.1,-0.2,0.2,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-0.1,0.1,-0.05,0.05};
Float_t StStrangenessAna::p_add[20]         = {0.0,0.15,-0.2,0.2,-0.2,0.2,-0.2,0.2,0.0,0.0,0.4,-0.1,0.2,-0.2,0.0,0.0,-0.1,0.1,-0.05,0.05};
Float_t StStrangenessAna::nsLow_add[20]     = {0.0,0.15,-0.2,0.2,0.1,-0.1,0.0,0.0,-0.1,-0.2,0.0,0.1,0.0,0.0,-0.4,0.4,0.0,0.0,0.0,0.0};
Float_t StStrangenessAna::nsHigh_add[20]    = {0.0,-0.15,-0.2,0.2,0.0,0.0,0.0,0.0,0.1,0.2,0.0,0.1,0.0,0.0,-0.4,0.4,0.0,0.0,0.0,0.0};

Int_t StStrangenessAna::n_cuts = 14; // temporary for Lambda, anti-Lambda and K0S => Code developing
//----------------------------------------------------
StStrangenessAna::StStrangenessAna(Int_t energy, Int_t X_flag, Int_t List, Long64_t start_event, Long64_t stop_event, Int_t mode)
{
  mEnergy = energy;
  mX_flag = X_flag;
  mList = List;
  mStart_Event = start_event;
  mStop_Event = stop_event;
  mMode = mode;
  if(!mRefMultCorr)
  {
    mRefMultCorr = CentralityMaker::instance()->getRefMultCorr();
  }
  mStrangenessCorr = new StStrangenessCorr();
  mStrangenessCut = new StStrangenessCut(mEnergy);
  mStrangenessHistoManger = new StStrangenessHistoManger();
}

StStrangenessAna::~StStrangenessAna()
{
}
//----------------------------------------------------
// set Input/Output
void StStrangenessAna::setInputDir(const TString inputdir)
{
    mInputdir = inputdir.Copy();
    cout << "Input directory was set to: " << mInputdir.Data() << endl;
}
void StStrangenessAna::setOutputfile(const TString outputfile)
{
    mOutputfile = outputfile.Copy();
    cout << "Output file was set to: " << mOutputfile.Data() << endl;
}
void StStrangenessAna::setSEList(const TString iSEList)
{
    mSEList = iSEList.Copy();
    cout << "Same event list was set to: " << mSEList.Data() << endl;
}
void StStrangenessAna::setMEList(const TString iMEList)
{
    mMEList = iMEList.Copy();
    cout << "Mixed event list was set to: " << mMEList.Data() << endl;
}
void StStrangenessAna::setStopEvent_SE(const Long64_t StopEvent_SE)
{
    mStopEvent_SE = StopEvent_SE;
    cout << "nStopEvent_SE = " << mStopEvent_SE << endl;
}
void StStrangenessAna::setStartEvent_SE(const Long64_t StartEvent_SE)
{
    mStartEvent_SE = StartEvent_SE;
    cout << "nStartEvent_SE = " << mStartEvent_SE << endl;
}
void StStrangenessAna::setStopEvent_ME(const Long64_t StopEvent_ME)
{
    mStopEvent_ME = StopEvent_ME;
    cout << "nStopEvent_ME = " << mStopEvent_ME << endl;
}
void StStrangenessAna::setStartEvent_ME(const Long64_t StartEvent_ME)
{
    mStartEvent_ME = StartEvent_ME;
    cout << "nStartEvent_ME = " << mStartEvent_ME << endl;
}
//----------------------------------------------------
// initial functions
void StStrangenessAna::Init()
{
  mStrangenessCorr->InitReCenterCorrection(mEnergy);
  mStrangenessCorr->InitShiftCorrection(mEnergy);
  mStrangenessHistoManger->Init(mX_flag,mMode);

  TString inputdir = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/%s/",Strangeness::Energy[mEnergy].Data(),Strangeness::Partype[mMode].Data());
  setInputDir(inputdir);

  const Int_t list_start = Strangeness::mList_Delta*mList + 1; // start list
  const Int_t list_stop  = Strangeness::mList_Delta*(mList+1); // stop list

  if(mX_flag == 0)
  {
    TString SEList = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/List/%s_list/Split_SE_%s_%d_%d.list",Strangeness::Energy[mEnergy].Data(),Strangeness::Partype[mMode].Data(),Strangeness::Energy[mEnergy].Data(),list_start,list_stop);
    setSEList(SEList);
    cout << SEList << endl;

    TString outputfile = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/%s/flow_%s/Yields_SE_%s_%d_%d.root",Strangeness::Energy[mEnergy].Data(),Strangeness::Partype[mMode].Data(),Strangeness::Partype[mMode].Data(),Strangeness::Energy[mEnergy].Data(),list_start,list_stop);
    setOutputfile(outputfile);

    setStartEvent_SE(Long64_t(mStart_Event));
    setStopEvent_SE(Long64_t(mStop_Event));

    InitSE();
  }

  if(mX_flag == 1)
  {
    TString MEList = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/List/%s_list/Split_ME_%s_%d_%d.list",Strangeness::Energy[mEnergy].Data(),Strangeness::Partype[mMode].Data(),Strangeness::Energy[mEnergy].Data(),list_start,list_stop);
    setMEList(MEList);

    TString outputfile = Form("/global/project/projectdirs/starprod/rnc/xusun/OutPut/AuAu%s/%s/flow_%s/Yields_ME_%s_%d_%d.root",Strangeness::Energy[mEnergy].Data(),Strangeness::Partype[mMode].Data(),Strangeness::Partype[mMode].Data(),Strangeness::Energy[mEnergy].Data(),list_start,list_stop);
    setOutputfile(outputfile);

    setStartEvent_ME(Long64_t(mStart_Event));
    setStopEvent_ME(Long64_t(mStop_Event));

    InitME();
  }
}

// Initialize Same Event
void StStrangenessAna::InitSE()
{
  TString Notification = Form("Initializing parameters and input/output for %s Same Event",Strangeness::Partype[mMode].Data());
  cout << Notification.Data() << endl;
  mFile_OutPut = new TFile(mOutputfile.Data(),"RECREATE");

  XUV0_EVENT_TREE       = (char*)Strangeness::v0_tree[mMode].Data();
  XUV0_EVENT_BRANCH     = (char*)Strangeness::v0_branch[mMode].Data();

  //----------------------------------------------------------------------------------------------------
  // Same event input
  if (!mSEList.IsNull())   // if input file is ok
  {
    cout << "Open same event file list " << mSEList << endl;
    ifstream in(mSEList);  // input stream
    if(in)
    {
      cout << "file list is ok" << endl;
      mInPut_SE  = new TChain( XUV0_EVENT_TREE, XUV0_EVENT_TREE );
      char str[255];       // char array for each file name
      Long64_t entries_save = 0;
      while(in)
      {
	in.getline(str,255);  // take the lines of the file list
	if(str[0] != 0)
	{
	  TString addfile;
	  addfile = str;
	  addfile = mInputdir+addfile;
	  mInPut_SE->AddFile(addfile.Data(),-1, XUV0_EVENT_TREE );
	  Long64_t file_entries = mInPut_SE->GetEntries();
	  cout << "File added to data chain: " << addfile.Data() << " with " << (file_entries-entries_save) << " entries" << endl;
	  entries_save = file_entries;
	}
      }
    }
    else
    {
      cout << "WARNING: SE file input is problemtic" << endl;
      mSE_input_flag = 0;
    }
  }

  // Set the input tree
  if (mSE_input_flag == 1 && !mInPut_SE->GetBranch( XUV0_EVENT_BRANCH ))
  {
    cerr << "ERROR: Could not find branch '"
      << XUV0_EVENT_BRANCH << "'in tree!" << endl;
  }

  if(mMode == 0) mXuPhiMeson_event = new StAlexPhiMesonEvent();
  if(mMode == 1) mLambda_event = new StV0Event();
  if(mMode == 2) mLambda_event = new StV0Event();
  if(mMode == 3) mK0S_event = new StV0Event();

  if(mSE_input_flag == 1)
  {
    if(mMode == 0) mInPut_SE->SetBranchAddress( XUV0_EVENT_BRANCH, &mXuPhiMeson_event );
    if(mMode == 1) mInPut_SE->SetBranchAddress( XUV0_EVENT_BRANCH, &mLambda_event);
    if(mMode == 2) mInPut_SE->SetBranchAddress( XUV0_EVENT_BRANCH, &mLambda_event);
    if(mMode == 3) mInPut_SE->SetBranchAddress( XUV0_EVENT_BRANCH, &mK0S_event);

    Int_t num_events_SE = mInPut_SE->GetEntriesFast();
    cout << "Number of events in file(s) = " << num_events_SE << endl;
    if(mStartEvent_SE > num_events_SE) mStartEvent_SE = num_events_SE;
    if(mStopEvent_SE  > num_events_SE) mStopEvent_SE  = num_events_SE;

    cout << "New nStartEvent_SE = " << mStartEvent_SE << ", new nStopEvent_SE = " << mStopEvent_SE << endl;
  }
}

// Initialize Mixed Event
void StStrangenessAna::InitME()
{
  TString Notification = Form("Initializing parameters and input/output for %s Mixed Event",Strangeness::Partype[mMode].Data());
  cout << Notification.Data() << endl;

  mFile_OutPut = new TFile(mOutputfile.Data(),"RECREATE");

  XUV0_EVENT_TREE       = (char*)Strangeness::v0_tree[mMode].Data();
  XUV0_EVENT_BRANCH     = (char*)Strangeness::v0_branch[mMode].Data();

  //----------------------------------------------------------------------------------------------------
  // Same event input
  if (!mMEList.IsNull())   // if input file is ok
  {
    cout << "Open mixed event file list " << mMEList << endl;
    ifstream in(mMEList);  // input stream
    if(in)
    {
      cout << "file list is ok" << endl;
      mInPut_ME  = new TChain( XUV0_EVENT_TREE, XUV0_EVENT_TREE );
      char str[255];       // char array for each file name
      Long64_t entries_save = 0;
      while(in)
      {
	in.getline(str,255);  // take the lines of the file list
	if(str[0] != 0)
	{
	  TString addfile;
	  addfile = str;
	  addfile = mInputdir+addfile;
	  mInPut_ME->AddFile(addfile.Data(),-1, XUV0_EVENT_TREE );
	  Long64_t file_entries = mInPut_ME->GetEntries();
	  cout << "File added to data chain: " << addfile.Data() << " with " << (file_entries-entries_save) << " entries" << endl;
	  entries_save = file_entries;
	}
      }
    }
    else
    {
      cout << "WARNING: ME file input is problemtic" << endl;
      mME_input_flag = 0;
    }
  }

  // Set the input tree
  if (mME_input_flag == 1 && !mInPut_ME->GetBranch( XUV0_EVENT_BRANCH ))
  {
    cerr << "ERROR: Could not find branch '"
      << XUV0_EVENT_BRANCH << "'in tree!" << endl;
  }

  if(mMode == 0) mXuPhiMeson_event = new StAlexPhiMesonEvent();
  if(mMode == 1) mLambda_event = new StV0Event();
  if(mMode == 2) mLambda_event = new StV0Event();
  if(mMode == 3) mK0S_event = new StV0Event();

  if(mME_input_flag == 1)
  {
    if(mMode == 0) mInPut_ME->SetBranchAddress( XUV0_EVENT_BRANCH, &mXuPhiMeson_event );
    if(mMode == 1) mInPut_ME->SetBranchAddress( XUV0_EVENT_BRANCH, &mLambda_event);
    if(mMode == 2) mInPut_ME->SetBranchAddress( XUV0_EVENT_BRANCH, &mLambda_event);
    if(mMode == 3) mInPut_ME->SetBranchAddress( XUV0_EVENT_BRANCH, &mK0S_event);

    Int_t num_events_ME = mInPut_ME->GetEntriesFast();
    cout << "Number of events in file(s) = " << num_events_ME << endl;
    if(mStartEvent_ME > num_events_ME) mStartEvent_ME = num_events_ME;
    if(mStopEvent_ME  > num_events_ME) mStopEvent_ME  = num_events_ME;

    cout << "New nStartEvent_ME = " << mStartEvent_ME << ", new nStopEvent_ME = " << mStopEvent_ME << endl;
  }
}
//----------------------------------------------------
void StStrangenessAna::Make()
{
  if(mX_flag == 0)
  {
    if(mMode == 0) MakePhiSE();
    if(mMode == 1) MakeLambdaSE();
    if(mMode == 2) MakeLambdaSE();
    if(mMode == 3) MakeK0S(mX_flag);
  }

  if(mX_flag == 1)
  {
    if(mMode == 0) MakePhiME();
    if(mMode == 1) MakeLambdaME();
    if(mMode == 2) MakeLambdaME();
    if(mMode == 3) MakeK0S(mX_flag);
  }
}


// loop phi meson Same Event
void StStrangenessAna::MakePhiSE()
{
  Long64_t start_event_use = mStartEvent_SE;
  Long64_t stop_event_use  = mStopEvent_SE;

  mInPut_SE->SetBranchAddress( XUV0_EVENT_BRANCH, &mXuPhiMeson_event );
  mInPut_SE->GetEntry(0); // For unknown reasons root doesn't like it if someone starts to read a file not from the 0 entry

  // Initialise Event Head
  StThreeVectorF PrimaryVertex(0.0,0.0,0.0);
  Int_t          RunId = 0;
  Int_t          RefMult = 0;
  Int_t          Centrality = 0;
  Int_t          N_prim = 0;
  Int_t          N_non_prim = 0;
  Int_t          N_Tof_match = 0;
  Float_t        ZDCx = 0.0; 
  Float_t        BBCx = 0.0; 
  Float_t        VzVpd = 0.0;
  Int_t          NumTrackUsed = 0;
  // ---------------------------------------QVector---------------------------------------------
  TVector2       Q2East[4];
  TVector2       Q2West[4];
  TVector2       Q3East[4];
  TVector2       Q3West[4];
  // -----------------------------------Number of Tracks----------------------------------------
  Int_t          NumTrackEast[4];
  Int_t          NumTrackWest[4];
  for(Int_t j = 0; j < 4; j++)
  {
    Q2East[j].Set(0.0,0.0);
    Q2West[j].Set(0.0,0.0);
    Q3East[j].Set(0.0,0.0);
    Q3West[j].Set(0.0,0.0);
    NumTrackEast[j] = 0;
    NumTrackWest[j] = 0;
  }

  for(Long64_t counter = start_event_use; counter < stop_event_use; counter++)
  {
    if (!mInPut_SE->GetEntry( counter )) // take the event -> information is stored in event
      break;  // end of data chunk

    // get Event Header
    PrimaryVertex     = mXuPhiMeson_event->getPrimaryVertex();
    RunId             = mXuPhiMeson_event->getRunId();
    RefMult           = mXuPhiMeson_event->getRefMult();
    Centrality        = mXuPhiMeson_event->getCentrality();
    N_prim            = mXuPhiMeson_event->getN_prim();
    N_non_prim        = mXuPhiMeson_event->getN_non_prim();
    N_Tof_match       = mXuPhiMeson_event->getN_Tof_match();
    ZDCx              = mXuPhiMeson_event->getZDCx(); 
    BBCx              = mXuPhiMeson_event->getBBCx(); 
    VzVpd             = mXuPhiMeson_event->getVzVpd();
    NumTrackUsed      = mXuPhiMeson_event->getNumTracks();

    for(Int_t j = 0; j < 4; j++)
    {
      Q2East[j]       = mXuPhiMeson_event->getQ2East(j);
      Q2West[j]       = mXuPhiMeson_event->getQ2West(j);
      Q3East[j]       = mXuPhiMeson_event->getQ3East(j);
      Q3West[j]       = mXuPhiMeson_event->getQ3West(j);
      NumTrackEast[j] = mXuPhiMeson_event->getNumTrackEast(j);
      NumTrackWest[j] = mXuPhiMeson_event->getNumTrackWest(j);
    }
 
    // Initialise Track 
    Float_t m2A = -999.9;
    Float_t m2B = -999.9;
    Float_t nsA = -999.9;
    Float_t nsB = -999.9;
    Float_t dcaA = -999.9;
    Float_t dcaB = -999.9;
    TLorentzVector lTrackA(0.0,0.0,0.0,0.0);
    TLorentzVector lTrackB(0.0,0.0,0.0,0.0);
    Int_t flagA = -1;
    Int_t flagB = -1;

    // vz sign
    Int_t vz_sign;
    if(PrimaryVertex.z() > 0.0)
    {
      vz_sign = 0;
    }
    else
    {
      vz_sign = 1;
    }

    // Centrality
    mRefMultCorr->init(RunId);
    if(mEnergy == 0) mRefMultCorr->initEvent(RefMult,PrimaryVertex.z(),ZDCx);
    if(mEnergy != 0) mRefMultCorr->initEvent(RefMult,PrimaryVertex.z());
    const Int_t cent9 = Centrality;
    const Double_t reweight = mRefMultCorr->getWeight();

    // runIndex
    mRunIdEventsDb = StRunIdEventsDb::Instance(Strangeness::mBeamEnergy[mEnergy],Strangeness::mBeamYear[mEnergy]);
    const Int_t runIndex = mRunIdEventsDb->getRunIdIndex(RunId); // expensive
    // cout << runIndex << endl;

    if (counter != 0  &&  counter % 1000 == 0)
      cout << "." << flush;
    if (counter != 0  &&  counter % 10000 == 0)
    {
      if((stop_event_use-start_event_use) > 0)
      {
	Double_t event_percent = 100.0*((Double_t)(counter-start_event_use))/((Double_t)(stop_event_use-start_event_use));
	cout << " " << counter-start_event_use << " (" << event_percent << "%) " << "\n" << "==> Processing data (strangeness_flow) " << flush;
      }
    }

    // get Track Information
    for(Int_t j = 0; j < Strangeness::mEtaGap_total; j++)
    {
      if(mStrangenessCorr->passTrackNumCut(NumTrackEast[j],NumTrackWest[j]))
      {
	for(UShort_t nTracks = 0; nTracks < NumTrackUsed; nTracks++) // loop over all tracks of the actual event
	{
	  mXuPhiMeson_track = mXuPhiMeson_event->getTrack(nTracks);
	  m2A = mXuPhiMeson_track->getMass2A();
	  m2B = mXuPhiMeson_track->getMass2B();
	  nsA = mXuPhiMeson_track->getNSigKaonA();
	  nsB = mXuPhiMeson_track->getNSigKaonB();
	  dcaA = mXuPhiMeson_track->getDcaA();
	  dcaB = mXuPhiMeson_track->getDcaB();
	  lTrackA = mXuPhiMeson_track->getTrackA();
	  lTrackB = mXuPhiMeson_track->getTrackB();
	  flagA = mXuPhiMeson_track->getFlagA();
	  flagB = mXuPhiMeson_track->getFlagB();

	  Float_t pA = lTrackA.P();
	  Float_t pB = lTrackB.P();
	  TLorentzVector lTrack = lTrackA + lTrackB;
	  Float_t pt_lTrack = lTrack.Perp();

	  for(Int_t i_cut = 0; i_cut < Strangeness::N_cuts; i_cut++)
	  {
	    // apply additional PID cut to increase significance
	    if(
		((fabs(pA) <= 0.65+p_add[i_cut] && m2A < -10) || (m2A > 0 && ((fabs(pA) < 1.5 && m2A > 0.16 && m2A < 0.36) || (fabs(pA) >= 1.5 && m2A > 0.125 && m2A < 0.36)) )) &&
		((fabs(pB) <= 0.65+p_add[i_cut] && m2B < -10) || (m2B > 0 && ((fabs(pB) < 1.5 && m2B > 0.16 && m2B < 0.36) || (fabs(pB) >= 1.5 && m2B > 0.125 && m2B < 0.36)) )) &&
		((pt_lTrack+pt_add[i_cut]) < 0.8 || ((pt_lTrack+pt_add[i_cut]) >= 0.8 && ( (m2A > 0.16 && m2A < 0.36) || (m2B > 0.16 && m2B < 0.36)))) &&
		(
		 ((m2A < -10 && nsA < 2.5+nsHigh_add[i_cut] && nsA > -(1.5+nsLow_add[i_cut])) || (m2A > 0.16 && m2A < 0.36)) &&
		 ((m2B < -10 && nsB < 2.5+nsHigh_add[i_cut] && nsB > -(1.5+nsLow_add[i_cut])) || (m2B > 0.16 && m2B < 0.36))
		)
	      )
	    {
	      Float_t phi_lTrack = lTrack.Phi();
	      Float_t InvMass_lTrack = lTrack.M();

	      if(mStrangenessCut->passPhiEtaEast(lTrack)) // neg eta(east)
	      { // Below is West Only
		TVector2 Q2Vector = Q2West[j];
		TVector2 Q3Vector = Q3West[j];
		// subtract auto-correlation from pos eta(west) event plane
		if(flagA == 0 && mStrangenessCut->passTrackEP(lTrackA,dcaA) && mStrangenessCut->passTrackEtaWest(lTrackA,j)) // trackA
		{
		  Float_t  w = mStrangenessCorr->getWeight(lTrackA);

		  TVector2 q2VectorA = mStrangenessCorr->calq2Vector(lTrackA);
		  TVector2 q2CorrA   = mStrangenessCorr->getReCenterPar_West(0,cent9,runIndex,vz_sign,j,0); // 2nd
		  Q2Vector = Q2Vector - w*(q2VectorA-q2CorrA);

		  TVector2 q3VectorA = mStrangenessCorr->calq3Vector(lTrackA);
		  TVector2 q3CorrA   = mStrangenessCorr->getReCenterPar_West(1,cent9,runIndex,vz_sign,j,0); // 3rd
		  Q3Vector = Q3Vector - w*(q3VectorA-q3CorrA);
		}
		if(flagB == 0 && mStrangenessCut->passTrackEP(lTrackB,dcaB) && mStrangenessCut->passTrackEtaWest(lTrackB,j)) // trackB
		{
		  Float_t  w = mStrangenessCorr->getWeight(lTrackB);

		  TVector2 q2VectorB = mStrangenessCorr->calq2Vector(lTrackB);
		  TVector2 q2CorrB   = mStrangenessCorr->getReCenterPar_West(0,cent9,runIndex,vz_sign,j,0); // 2nd
		  Q2Vector = Q2Vector - w*(q2VectorB-q2CorrB);

		  TVector2 q3VectorB = mStrangenessCorr->calq3Vector(lTrackB);
		  TVector2 q3CorrB   = mStrangenessCorr->getReCenterPar_West(1,cent9,runIndex,vz_sign,j,0); // 3rd
		  Q3Vector = Q3Vector - w*(q3VectorB-q3CorrB);
		}
		Float_t Res2 = mStrangenessCorr->getResolution2_EP(cent9,j);
		Float_t Psi2_west = mStrangenessCorr->calShiftAngle2West_EP(Q2Vector,runIndex,cent9,vz_sign,j);
		Float_t Res3 = mStrangenessCorr->getResolution3_EP(cent9,j);
		Float_t Psi3_west = mStrangenessCorr->calShiftAngle3West_EP(Q3Vector,runIndex,cent9,vz_sign,j);
		Float_t phi_Psi2 = phi_lTrack - Psi2_west;
		Float_t phi_Psi3 = phi_lTrack - Psi3_west;

		mStrangenessHistoManger->Fill(pt_lTrack,cent9,j,phi_Psi2,Res2,phi_Psi3,Res3,InvMass_lTrack,reweight,i_cut);
	      }

	      if(mStrangenessCut->passPhiEtaWest(lTrack)) // pos eta
	      { // Below is East Only
		TVector2 Q2Vector = Q2East[j];
		TVector2 Q3Vector = Q3East[j];
		// subtract auto-correlation from pos eta(west) event plane
		if(flagA == 0 && mStrangenessCut->passTrackEP(lTrackA,dcaA) && mStrangenessCut->passTrackEtaEast(lTrackA,j)) // trackA
		{
		  Float_t  w = mStrangenessCorr->getWeight(lTrackA);

		  TVector2 q2VectorA = mStrangenessCorr->calq2Vector(lTrackA);
		  TVector2 q2CorrA   = mStrangenessCorr->getReCenterPar_East(0,cent9,runIndex,vz_sign,j,0); // 2nd
		  Q2Vector = Q2Vector - w*(q2VectorA-q2CorrA);

		  TVector2 q3VectorA = mStrangenessCorr->calq3Vector(lTrackA);
		  TVector2 q3CorrA   = mStrangenessCorr->getReCenterPar_East(1,cent9,runIndex,vz_sign,j,0); // 3rd
		  Q3Vector = Q3Vector - w*(q3VectorA-q3CorrA);
		}
		if(flagB == 0 && mStrangenessCut->passTrackEP(lTrackB,dcaB) && mStrangenessCut->passTrackEtaEast(lTrackB,j)) // trackB
		{
		  Float_t  w = mStrangenessCorr->getWeight(lTrackB);

		  TVector2 q2VectorB = mStrangenessCorr->calq2Vector(lTrackB);
		  TVector2 q2CorrB   = mStrangenessCorr->getReCenterPar_East(0,cent9,runIndex,vz_sign,j,0); // 2nd
		  Q2Vector = Q2Vector - w*(q2VectorB-q2CorrB);

		  TVector2 q3VectorB = mStrangenessCorr->calq3Vector(lTrackB);
		  TVector2 q3CorrB   = mStrangenessCorr->getReCenterPar_East(1,cent9,runIndex,vz_sign,j,0); // 3rd
		  Q3Vector = Q3Vector - w*(q3VectorB-q3CorrB);
		}
		Float_t Res2 = mStrangenessCorr->getResolution2_EP(cent9,j);
		Float_t Psi2_east = mStrangenessCorr->calShiftAngle2East_EP(Q2Vector,runIndex,cent9,vz_sign,j);
		Float_t Res3 = mStrangenessCorr->getResolution3_EP(cent9,j);
		Float_t Psi3_east = mStrangenessCorr->calShiftAngle3East_EP(Q3Vector,runIndex,cent9,vz_sign,j);
		Float_t phi_Psi2 = phi_lTrack - Psi2_east;
		Float_t phi_Psi3 = phi_lTrack - Psi3_east;

		mStrangenessHistoManger->Fill(pt_lTrack,cent9,j,phi_Psi2,Res2,phi_Psi3,Res3,InvMass_lTrack,reweight,i_cut);
	      }
	    }
	  }
	}
      }
    }
  }

  cout << "." << flush;
  cout << " " << stop_event_use-start_event_use << "(" << 100 << "%)";
  cout << endl;
}

// loop phi meson Mixed Event
void StStrangenessAna::MakePhiME()
{
  Long64_t start_event_use = mStartEvent_ME;
  Long64_t stop_event_use  = mStopEvent_ME;

  mInPut_ME->SetBranchAddress( XUV0_EVENT_BRANCH, &mXuPhiMeson_event );
  mInPut_ME->GetEntry(0); // For unknown reasons root doesn't like it if someone starts to read a file not from the 0 entry

  // Initialise Event Head
  StThreeVectorF PrimaryVertex(0.0,0.0,0.0);
  Int_t          RunId = 0;
  Int_t          RefMult = 0;
  Int_t          Centrality = 0;
  Int_t          N_prim = 0;
  Int_t          N_non_prim = 0;
  Int_t          N_Tof_match = 0;
  Float_t        ZDCx = 0.0; 
  Float_t        BBCx = 0.0; 
  Float_t        VzVpd = 0.0;
  Int_t          NumTrackUsed = 0;
  // ---------------------------------------QVector---------------------------------------------
  TVector2       Q2East[4];
  TVector2       Q2West[4];
  TVector2       Q3East[4];
  TVector2       Q3West[4];
  // -----------------------------------Number of Tracks----------------------------------------
  Int_t          NumTrackEast[4];
  Int_t          NumTrackWest[4];
  for(Int_t j = 0; j < 4; j++)
  {
    Q2East[j].Set(0.0,0.0);
    Q2West[j].Set(0.0,0.0);
    Q3East[j].Set(0.0,0.0);
    Q3West[j].Set(0.0,0.0);
    NumTrackEast[j] = 0;
    NumTrackWest[j] = 0;
  }

  for(Long64_t counter = start_event_use; counter < stop_event_use; counter++)
  {
    if (!mInPut_ME->GetEntry( counter )) // take the event -> information is stored in event
      break;  // end of data chunk

    // get Event Header
    PrimaryVertex     = mXuPhiMeson_event->getPrimaryVertex();
    RunId             = mXuPhiMeson_event->getRunId();
    RefMult           = mXuPhiMeson_event->getRefMult();
    Centrality        = mXuPhiMeson_event->getCentrality();
    N_prim            = mXuPhiMeson_event->getN_prim();
    N_non_prim        = mXuPhiMeson_event->getN_non_prim();
    N_Tof_match       = mXuPhiMeson_event->getN_Tof_match();
    ZDCx              = mXuPhiMeson_event->getZDCx(); 
    BBCx              = mXuPhiMeson_event->getBBCx(); 
    VzVpd             = mXuPhiMeson_event->getVzVpd();
    NumTrackUsed      = mXuPhiMeson_event->getNumTracks();

    for(Int_t j = 0; j < 4; j++)
    {
      Q2East[j]       = mXuPhiMeson_event->getQ2East(j);
      Q2West[j]       = mXuPhiMeson_event->getQ2West(j);
      Q3East[j]       = mXuPhiMeson_event->getQ3East(j);
      Q3West[j]       = mXuPhiMeson_event->getQ3West(j);
      NumTrackEast[j] = mXuPhiMeson_event->getNumTrackEast(j);
      NumTrackWest[j] = mXuPhiMeson_event->getNumTrackWest(j);
    }
 
    // Initialise Track 
    Float_t m2A = -999.9;
    Float_t m2B = -999.9;
    Float_t nsA = -999.9;
    Float_t nsB = -999.9;
    Float_t dcaA = -999.9;
    Float_t dcaB = -999.9;
    TLorentzVector lTrackA(0.0,0.0,0.0,0.0);
    TLorentzVector lTrackB(0.0,0.0,0.0,0.0);
    Int_t flagA = -1;
    Int_t flagB = -1;

    // vz sign
    Int_t vz_sign;
    if(PrimaryVertex.z() > 0.0)
    {
      vz_sign = 0;
    }
    else
    {
      vz_sign = 1;
    }

    // Centrality
    mRefMultCorr->init(RunId);
    if(mEnergy == 0) mRefMultCorr->initEvent(RefMult,PrimaryVertex.z(),ZDCx);
    if(mEnergy != 0) mRefMultCorr->initEvent(RefMult,PrimaryVertex.z());
    const Int_t cent9 = Centrality;
    const Double_t reweight = mRefMultCorr->getWeight();

    // runIndex
    mRunIdEventsDb = StRunIdEventsDb::Instance(Strangeness::mBeamEnergy[mEnergy],Strangeness::mBeamYear[mEnergy]);
    const Int_t runIndex = mRunIdEventsDb->getRunIdIndex(RunId); // expensive
    // cout << runIndex << endl;

    if (counter != 0  &&  counter % 1000 == 0)
      cout << "." << flush;
    if (counter != 0  &&  counter % 10000 == 0)
    {
      if((stop_event_use-start_event_use) > 0)
      {
	Double_t event_percent = 100.0*((Double_t)(counter-start_event_use))/((Double_t)(stop_event_use-start_event_use));
	cout << " " << counter-start_event_use << " (" << event_percent << "%) " << "\n" << "==> Processing data (strangeness_flow) " << flush;
      }
    }

    // get Track Information
    for(Int_t j = 0; j < Strangeness::mEtaGap_total; j++)
    {
      if(mStrangenessCorr->passTrackNumCut(NumTrackEast[j],NumTrackWest[j]))
      {
	for(UShort_t nTracks = 0; nTracks < NumTrackUsed; nTracks++) // loop over all tracks of the actual event
	{
	  mXuPhiMeson_track = mXuPhiMeson_event->getTrack(nTracks);
	  m2A = mXuPhiMeson_track->getMass2A();
	  m2B = mXuPhiMeson_track->getMass2B();
	  nsA = mXuPhiMeson_track->getNSigKaonA();
	  nsB = mXuPhiMeson_track->getNSigKaonB();
	  dcaA = mXuPhiMeson_track->getDcaA();
	  dcaB = mXuPhiMeson_track->getDcaB();
	  lTrackA = mXuPhiMeson_track->getTrackA();
	  lTrackB = mXuPhiMeson_track->getTrackB();
	  Float_t pA = lTrackA.P();
	  Float_t pB = lTrackB.P();
	  flagA = mXuPhiMeson_track->getFlagA();
	  flagB = mXuPhiMeson_track->getFlagB();
	  TLorentzVector lTrack = lTrackA + lTrackB;
	  Float_t pt_lTrack = lTrack.Perp();

	  for(Int_t i_cut = 0; i_cut < Strangeness::N_cuts; i_cut++)
	  {
	    // apply additional PID cut to increase significance
	    if(
		((fabs(pA) <= 0.65+p_add[i_cut] && m2A < -10) || (m2A > 0 && ((fabs(pA) < 1.5 && m2A > 0.16 && m2A < 0.36) || (fabs(pA) >= 1.5 && m2A > 0.125 && m2A < 0.36)) )) &&
		((fabs(pB) <= 0.65+p_add[i_cut] && m2B < -10) || (m2B > 0 && ((fabs(pB) < 1.5 && m2B > 0.16 && m2B < 0.36) || (fabs(pB) >= 1.5 && m2B > 0.125 && m2B < 0.36)) )) &&
		((pt_lTrack+pt_add[i_cut]) < 0.8 || ((pt_lTrack+pt_add[i_cut]) >= 0.8 && ( (m2A > 0.16 && m2A < 0.36) || (m2B > 0.16 && m2B < 0.36)))) &&
		(
		 ((m2A < -10 && nsA < 2.5+nsHigh_add[i_cut] && nsA > -(1.5+nsLow_add[i_cut])) || (m2A > 0.16 && m2A < 0.36)) &&
		 ((m2B < -10 && nsB < 2.5+nsHigh_add[i_cut] && nsB > -(1.5+nsLow_add[i_cut])) || (m2B > 0.16 && m2B < 0.36))
		)
	      )
	    {
	      Float_t phi_lTrack = lTrack.Phi();
	      Float_t InvMass_lTrack = lTrack.M();

	      if(mStrangenessCut->passPhiEtaEast(lTrack)) // neg eta(east)
	      { // Below is West Only
		TVector2 Q2Vector = Q2West[j];
		TVector2 Q3Vector = Q3West[j];
		// subtract auto-correlation from pos eta(west) event plane
		if(flagA == 0 && mStrangenessCut->passTrackEP(lTrackA,dcaA) && mStrangenessCut->passTrackEtaWest(lTrackA,j)) // trackA
		{
		  Float_t  w = mStrangenessCorr->getWeight(lTrackA);

		  TVector2 q2VectorA = mStrangenessCorr->calq2Vector(lTrackA);
		  TVector2 q2CorrA   = mStrangenessCorr->getReCenterPar_West(0,cent9,runIndex,vz_sign,j,0); // 2nd
		  Q2Vector = Q2Vector - w*(q2VectorA-q2CorrA);

		  TVector2 q3VectorA = mStrangenessCorr->calq3Vector(lTrackA);
		  TVector2 q3CorrA   = mStrangenessCorr->getReCenterPar_West(1,cent9,runIndex,vz_sign,j,0); // 3rd
		  Q3Vector = Q3Vector - w*(q3VectorA-q3CorrA);
		}
		if(flagB == 0 && mStrangenessCut->passTrackEP(lTrackB,dcaB) && mStrangenessCut->passTrackEtaWest(lTrackB,j)) // trackB
		{
		  Float_t  w = mStrangenessCorr->getWeight(lTrackB);

		  TVector2 q2VectorB = mStrangenessCorr->calq2Vector(lTrackB);
		  TVector2 q2CorrB   = mStrangenessCorr->getReCenterPar_West(0,cent9,runIndex,vz_sign,j,0); // 2nd
		  Q2Vector = Q2Vector - w*(q2VectorB-q2CorrB);

		  TVector2 q3VectorB = mStrangenessCorr->calq3Vector(lTrackB);
		  TVector2 q3CorrB   = mStrangenessCorr->getReCenterPar_West(1,cent9,runIndex,vz_sign,j,0); // 3rd
		  Q3Vector = Q3Vector - w*(q3VectorB-q3CorrB);
		}
		Float_t Res2 = mStrangenessCorr->getResolution2_EP(cent9,j);
		Float_t Psi2_west = mStrangenessCorr->calShiftAngle2West_EP(Q2Vector,runIndex,cent9,vz_sign,j);
		Float_t Res3 = mStrangenessCorr->getResolution3_EP(cent9,j);
		Float_t Psi3_west = mStrangenessCorr->calShiftAngle3West_EP(Q3Vector,runIndex,cent9,vz_sign,j);
		Float_t phi_Psi2 = phi_lTrack - Psi2_west;
		Float_t phi_Psi3 = phi_lTrack - Psi3_west;

		mStrangenessHistoManger->Fill(pt_lTrack,cent9,j,phi_Psi2,Res2,phi_Psi3,Res3,InvMass_lTrack,reweight,i_cut);
	      }

	      if(mStrangenessCut->passPhiEtaWest(lTrack)) // pos eta
	      { // Below is East Only
		TVector2 Q2Vector = Q2East[j];
		TVector2 Q3Vector = Q3East[j];
		// subtract auto-correlation from neg eta(east) event plane
		if(flagA == 0 && mStrangenessCut->passTrackEP(lTrackA,dcaA) && mStrangenessCut->passTrackEtaEast(lTrackA,j)) // trackA
		{
		  Float_t  w = mStrangenessCorr->getWeight(lTrackA);

		  TVector2 q2VectorA = mStrangenessCorr->calq2Vector(lTrackA);
		  TVector2 q2CorrA   = mStrangenessCorr->getReCenterPar_East(0,cent9,runIndex,vz_sign,j,0); // 2nd
		  Q2Vector = Q2Vector - w*(q2VectorA-q2CorrA);

		  TVector2 q3VectorA = mStrangenessCorr->calq3Vector(lTrackA);
		  TVector2 q3CorrA   = mStrangenessCorr->getReCenterPar_East(1,cent9,runIndex,vz_sign,j,0); // 3rd
		  Q3Vector = Q3Vector - w*(q3VectorA-q3CorrA);
		}
		if(flagB == 0 && mStrangenessCut->passTrackEP(lTrackB,dcaB) && mStrangenessCut->passTrackEtaEast(lTrackB,j)) // trackB
		{
		  Float_t  w = mStrangenessCorr->getWeight(lTrackB);

		  TVector2 q2VectorB = mStrangenessCorr->calq2Vector(lTrackB);
		  TVector2 q2CorrB   = mStrangenessCorr->getReCenterPar_East(0,cent9,runIndex,vz_sign,j,0); // 2nd
		  Q2Vector = Q2Vector - w*(q2VectorB-q2CorrB);

		  TVector2 q3VectorB = mStrangenessCorr->calq3Vector(lTrackB);
		  TVector2 q3CorrB   = mStrangenessCorr->getReCenterPar_East(1,cent9,runIndex,vz_sign,j,0); // 3rd
		  Q3Vector = Q3Vector - w*(q3VectorB-q3CorrB);
		}
		Float_t Res2 = mStrangenessCorr->getResolution2_EP(cent9,j);
		Float_t Psi2_east = mStrangenessCorr->calShiftAngle2East_EP(Q2Vector,runIndex,cent9,vz_sign,j);
		Float_t Res3 = mStrangenessCorr->getResolution3_EP(cent9,j);
		Float_t Psi3_east = mStrangenessCorr->calShiftAngle3East_EP(Q3Vector,runIndex,cent9,vz_sign,j);
		Float_t phi_Psi2 = phi_lTrack - Psi2_east;
		Float_t phi_Psi3 = phi_lTrack - Psi3_east;

		mStrangenessHistoManger->Fill(pt_lTrack,cent9,j,phi_Psi2,Res2,phi_Psi3,Res3,InvMass_lTrack,reweight,i_cut);
	      }
	    }
	  }
	}
      }
    }
  }

  cout << "." << flush;
  cout << " " << stop_event_use-start_event_use << "(" << 100 << "%)";
  cout << endl;
}

// loop Lambda Same Event
void StStrangenessAna::MakeLambdaSE()
{
  Long64_t start_event_use = mStartEvent_SE;
  Long64_t stop_event_use  = mStopEvent_SE;

  mInPut_SE->SetBranchAddress( XUV0_EVENT_BRANCH, &mLambda_event);
  mInPut_SE->GetEntry(0); // For unknown reasons root doesn't like it if someone starts to read a file not from the 0 entry

  // Initialise Event Head
  StThreeVectorF PrimaryVertex(0.0,0.0,0.0);
  Int_t          RunId = 0;
  Int_t          RefMult = 0;
  Int_t          Centrality = 0;
  Int_t          N_prim = 0;
  Int_t          N_non_prim = 0;
  Int_t          N_Tof_match = 0;
  Float_t        ZDCx = 0.0; 
  Float_t        BBCx = 0.0; 
  Float_t        VzVpd = 0.0;
  Int_t          NumTrackUsed = 0;
  // ---------------------------------------QVector---------------------------------------------
  TVector2       Q2East[4];
  TVector2       Q2West[4];
  TVector2       Q3East[4];
  TVector2       Q3West[4];
  // -----------------------------------Number of Tracks----------------------------------------
  Int_t          NumTrackEast[4];
  Int_t          NumTrackWest[4];
  for(Int_t j = 0; j < 4; j++)
  {
    Q2East[j].Set(0.0,0.0);
    Q2West[j].Set(0.0,0.0);
    Q3East[j].Set(0.0,0.0);
    Q3West[j].Set(0.0,0.0);
    NumTrackEast[j] = 0;
    NumTrackWest[j] = 0;
  }

  for(Long64_t counter = start_event_use; counter < stop_event_use; counter++)
  {
    if (!mInPut_SE->GetEntry( counter )) // take the event -> information is stored in event
      break;  // end of data chunk

    // get Event Header
    PrimaryVertex     = mLambda_event->getPrimaryVertex();
    RunId             = mLambda_event->getRunId();
    RefMult           = mLambda_event->getRefMult();
    Centrality        = mLambda_event->getCentrality();
    N_prim            = mLambda_event->getN_prim();
    N_non_prim        = mLambda_event->getN_non_prim();
    N_Tof_match       = mLambda_event->getN_Tof_match();
    ZDCx              = mLambda_event->getZDCx(); 
    BBCx              = mLambda_event->getBBCx(); 
    VzVpd             = mLambda_event->getVzVpd();
    NumTrackUsed      = mLambda_event->getNumTracks();

    for(Int_t j = 0; j < 4; j++)
    {
      Q2East[j]       = mLambda_event->getQ2East(j);
      Q2West[j]       = mLambda_event->getQ2West(j);
      Q3East[j]       = mLambda_event->getQ3East(j);
      Q3West[j]       = mLambda_event->getQ3West(j);
      NumTrackEast[j] = mLambda_event->getNumTrackEast(j);
      NumTrackWest[j] = mLambda_event->getNumTrackWest(j);
    }
 
    // Initialise Track 
    Float_t m2A = -999.9;
    Float_t m2B = -999.9;
    Float_t nsA = -999.9;
    Float_t nsB = -999.9;
    Float_t dcaA = -999.9;
    Float_t dcaB = -999.9;
    TLorentzVector lGTrackA(0.0,0.0,0.0,0.0);
    TLorentzVector lGTrackB(0.0,0.0,0.0,0.0);
    TLorentzVector lPTrackA(0.0,0.0,0.0,0.0);
    TLorentzVector lPTrackB(0.0,0.0,0.0,0.0);
    Int_t flagA = -1;
    Int_t flagB = -1;
    Float_t dcaAB = -999.9;
    Float_t decaylength = -999.9; 
    Float_t dcaV0 = -999.9;

    // vz sign
    Int_t vz_sign;
    if(PrimaryVertex.z() > 0.0)
    {
      vz_sign = 0;
    }
    else
    {
      vz_sign = 1;
    }

    // Centrality
    mRefMultCorr->init(RunId);
    if(mEnergy == 0) mRefMultCorr->initEvent(RefMult,PrimaryVertex.z(),ZDCx);
    if(mEnergy != 0) mRefMultCorr->initEvent(RefMult,PrimaryVertex.z());
    const Int_t cent9 = Centrality;
    const Double_t reweight = mRefMultCorr->getWeight();

    // runIndex
    mRunIdEventsDb = StRunIdEventsDb::Instance(Strangeness::mBeamEnergy[mEnergy],Strangeness::mBeamYear[mEnergy]);
    const Int_t runIndex = mRunIdEventsDb->getRunIdIndex(RunId); // expensive
    //cout << runIndex << endl;

    if (counter != 0  &&  counter % 1000 == 0)
      cout << "." << flush;
    if (counter != 0  &&  counter % 10000 == 0)
    {
      if((stop_event_use-start_event_use) > 0)
      {
	Double_t event_percent = 100.0*((Double_t)(counter-start_event_use))/((Double_t)(stop_event_use-start_event_use));
	cout << " " << counter-start_event_use << " (" << event_percent << "%) " << "\n" << "==> Processing data (strangeness_flow) " << flush;
      }
    }

    // get Track Information
    for(Int_t j = 0; j < Strangeness::mEtaGap_total; j++)
    {
      if(mStrangenessCorr->passTrackNumCut(NumTrackEast[j],NumTrackWest[j]))
      {
	for(UShort_t nTracks = 0; nTracks < NumTrackUsed; nTracks++) // loop over all tracks of the actual event
	{
	  mLambda_track = mLambda_event->getTrack(nTracks);
	  m2A = mLambda_track->getMass2A();
	  m2B = mLambda_track->getMass2B();
	  nsA = mLambda_track->getNSigA();
	  nsB = mLambda_track->getNSigB();
	  dcaA = mLambda_track->getDcaA();
	  dcaB = mLambda_track->getDcaB();
	  lGTrackA = mLambda_track->getGTrackA();
	  lGTrackB = mLambda_track->getGTrackB();
	  lPTrackA = mLambda_track->getPTrackA();
	  lPTrackB = mLambda_track->getPTrackB();
	  flagA = mLambda_track->getFlagA();
	  flagB = mLambda_track->getFlagB();
	  dcaAB = mLambda_track->getDcaAB();
	  decaylength = mLambda_track->getDecayLength(); 
	  dcaV0 = mLambda_track->getDcaV0();

	  Float_t GpA = lGTrackA.P();
	  Float_t GpB = lGTrackB.P();
	  TLorentzVector lGTrack = lGTrackA + lGTrackB; // Lorentz vector for mother particle
	  Float_t pt_lGTrack = lGTrack.Perp();

	  // apply additional PID cut to increase significance
	  // final mass2 cut
	  Float_t Mass2_low_proton = 0.7;
	  Float_t Mass2_up_proton  = 1.1;
	  Float_t Mass2_low_pi_minus;
	  Float_t Mass2_up_pi_minus;
	  if(GpB < 0.75)
	  {
	    Mass2_low_pi_minus = -0.013;
	    Mass2_up_pi_minus  =  0.047;
	  }
	  if(GpB >= 0.75)
	  {
	    Mass2_low_pi_minus =  0.053 - 0.088*GpB;
	    Mass2_up_pi_minus  = -0.004 + 0.068*GpB;
	  }

	  if(
	      ( // no particle has a mass
		(m2A < -10.0) && (m2B < -10.0)
		&& fabs(dcaA)    > 0.6
		&& fabs(dcaB)    > 1.7
		&& dcaAB         < 1.0
		&& decaylength   > 4.0
		&& dcaV0         < 0.75
	      )
	      ||
	      ( // only pion has a mass
		(m2A < -10.0) && (m2B > Mass2_low_pi_minus && m2B < Mass2_up_pi_minus)
		&& fabs(dcaA)    > 0.5
		&& fabs(dcaB)    > 1.5
		&& dcaAB         < 1.0
		&& decaylength   > 3.5
		&& dcaV0         < 0.75
	      )
	      ||
	      ( // only proton has a mass
		(m2A > Mass2_low_proton && m2A < Mass2_up_proton) && (m2B < -10.0)
		&& fabs(dcaA)    > 0.15
		&& fabs(dcaB)    > 0.8
		&& dcaAB         < 1.0
		&& decaylength   > 2.5
		&& dcaV0         < 1.2
	      )
	      ||
	      ( // both particles have a mass
		(m2A > Mass2_low_proton && m2A < Mass2_up_proton) && (m2B > Mass2_low_pi_minus && m2B < Mass2_up_pi_minus)
		&& fabs(dcaA)    > 0.1
		&& fabs(dcaB)    > 0.7
		&& dcaAB         < 1.0
		&& decaylength   > 2.0
		&& dcaV0         < 1.3
	      )
	    )
	  {
	    Float_t phi_lGTrack = lGTrack.Phi();
	    Float_t InvMass_lGTrack = lGTrack.M();

	    if(mStrangenessCut->passPhiEtaEast(lGTrack)) // neg eta(east)
	    { // Below is West Only
	      TVector2 Q2Vector = Q2West[j];
	      TVector2 Q3Vector = Q3West[j];
	      // subtract auto-correlation from pos eta(west) event plane
	      if(flagA == 0 && mStrangenessCut->passTrackEP(lPTrackA,dcaA) && mStrangenessCut->passTrackEtaWest(lPTrackA,j)) // trackA
	      {
		Float_t  w = mStrangenessCorr->getWeight(lPTrackA);

		TVector2 q2VectorA = mStrangenessCorr->calq2Vector(lPTrackA);
		TVector2 q2CorrA   = mStrangenessCorr->getReCenterPar_West(0,cent9,runIndex,vz_sign,j,0); // 2nd
		Q2Vector = Q2Vector - w*(q2VectorA-q2CorrA);

		TVector2 q3VectorA = mStrangenessCorr->calq3Vector(lPTrackA);
		TVector2 q3CorrA   = mStrangenessCorr->getReCenterPar_West(1,cent9,runIndex,vz_sign,j,0); // 3rd
		Q3Vector = Q3Vector - w*(q3VectorA-q3CorrA);
	      }
	      if(flagB == 0 && mStrangenessCut->passTrackEP(lPTrackB,dcaB) && mStrangenessCut->passTrackEtaWest(lPTrackB,j)) // trackB
	      {
		Float_t  w = mStrangenessCorr->getWeight(lPTrackB);

		TVector2 q2VectorB = mStrangenessCorr->calq2Vector(lPTrackB);
		TVector2 q2CorrB   = mStrangenessCorr->getReCenterPar_West(0,cent9,runIndex,vz_sign,j,0); // 2nd
		Q2Vector = Q2Vector - w*(q2VectorB-q2CorrB);

		TVector2 q3VectorB = mStrangenessCorr->calq3Vector(lPTrackB);
		TVector2 q3CorrB   = mStrangenessCorr->getReCenterPar_West(1,cent9,runIndex,vz_sign,j,0); // 3rd
		Q3Vector = Q3Vector - w*(q3VectorB-q3CorrB);
	      }
	      Float_t Res2 = mStrangenessCorr->getResolution2_EP(cent9,j);
	      Float_t Psi2_west = mStrangenessCorr->calShiftAngle2West_EP(Q2Vector,runIndex,cent9,vz_sign,j);
	      Float_t Res3 = mStrangenessCorr->getResolution3_EP(cent9,j);
	      Float_t Psi3_west = mStrangenessCorr->calShiftAngle3West_EP(Q3Vector,runIndex,cent9,vz_sign,j);
	      Float_t phi_Psi2 = phi_lGTrack - Psi2_west;
	      Float_t phi_Psi3 = phi_lGTrack - Psi3_west;

	      mStrangenessHistoManger->Fill(pt_lGTrack,cent9,j,phi_Psi2,Res2,phi_Psi3,Res3,InvMass_lGTrack,reweight,n_cuts);

	      TLorentzVector lGTrackC;
	      lGTrackC.SetXYZM(lGTrackA.X(),lGTrackA.Y(),lGTrackA.Z(),0.49368);
	      TLorentzVector lGTrackCB = lGTrackC + lGTrackB; // misidentify K* -> Lambda
	      Double_t InvMassCB = lGTrackCB.M(); // Inv Mass of K*
	      if(!(InvMassCB > 0.89594-3*0.0487 && InvMassCB < 0.89594+3*0.0487))
	      {
		mStrangenessHistoManger->Fill_sub(pt_lGTrack,cent9,j,phi_Psi2,Res2,phi_Psi3,Res3,InvMass_lGTrack,reweight,n_cuts);
	      }
	    }

	    if(mStrangenessCut->passPhiEtaWest(lGTrack)) // pos eta
	    { // Below is East Only
	      TVector2 Q2Vector = Q2East[j];
	      TVector2 Q3Vector = Q3East[j];
	      // subtract auto-correlation from neg eta(east) event plane
	      if(flagA == 0 && mStrangenessCut->passTrackEP(lPTrackA,dcaA) && mStrangenessCut->passTrackEtaEast(lPTrackA,j)) // trackA
	      {
		Float_t  w = mStrangenessCorr->getWeight(lPTrackA);

		TVector2 q2VectorA = mStrangenessCorr->calq2Vector(lPTrackA);
		TVector2 q2CorrA   = mStrangenessCorr->getReCenterPar_East(0,cent9,runIndex,vz_sign,j,0); // 2nd
		Q2Vector = Q2Vector - w*(q2VectorA-q2CorrA);

		TVector2 q3VectorA = mStrangenessCorr->calq3Vector(lPTrackA);
		TVector2 q3CorrA   = mStrangenessCorr->getReCenterPar_East(1,cent9,runIndex,vz_sign,j,0); // 3rd
		Q3Vector = Q3Vector - w*(q3VectorA-q3CorrA);
	      }
	      if(flagB == 0 && mStrangenessCut->passTrackEP(lPTrackB,dcaB) && mStrangenessCut->passTrackEtaEast(lPTrackB,j)) // trackB
	      {
		Float_t  w = mStrangenessCorr->getWeight(lPTrackB);

		TVector2 q2VectorB = mStrangenessCorr->calq2Vector(lPTrackB);
		TVector2 q2CorrB   = mStrangenessCorr->getReCenterPar_East(0,cent9,runIndex,vz_sign,j,0); // 2nd
		Q2Vector = Q2Vector - w*(q2VectorB-q2CorrB);

		TVector2 q3VectorB = mStrangenessCorr->calq3Vector(lPTrackB);
		TVector2 q3CorrB   = mStrangenessCorr->getReCenterPar_East(1,cent9,runIndex,vz_sign,j,0); // 3rd
		Q3Vector = Q3Vector - w*(q3VectorB-q3CorrB);
	      }
	      Float_t Res2 = mStrangenessCorr->getResolution2_EP(cent9,j);
	      Float_t Psi2_east = mStrangenessCorr->calShiftAngle2East_EP(Q2Vector,runIndex,cent9,vz_sign,j);
	      Float_t Res3 = mStrangenessCorr->getResolution3_EP(cent9,j);
	      Float_t Psi3_east = mStrangenessCorr->calShiftAngle3East_EP(Q3Vector,runIndex,cent9,vz_sign,j);
	      Float_t phi_Psi2 = phi_lGTrack - Psi2_east;
	      Float_t phi_Psi3 = phi_lGTrack - Psi3_east;

	      mStrangenessHistoManger->Fill(pt_lGTrack,cent9,j,phi_Psi2,Res2,phi_Psi3,Res3,InvMass_lGTrack,reweight,n_cuts);

	      TLorentzVector lGTrackC;
	      lGTrackC.SetXYZM(lGTrackA.X(),lGTrackA.Y(),lGTrackA.Z(),0.49368);
	      TLorentzVector lGTrackCB = lGTrackC + lGTrackB; // misidentify K* -> Lambda
	      Double_t InvMassCB = lGTrackCB.M(); // Inv Mass of K*
	      if(!(InvMassCB > 0.89594-3*0.0487 && InvMassCB < 0.89594+3*0.0487))
	      {
		mStrangenessHistoManger->Fill_sub(pt_lGTrack,cent9,j,phi_Psi2,Res2,phi_Psi3,Res3,InvMass_lGTrack,reweight,n_cuts);
	      }
	    }
	  }
	}
      }
    }
  }

  cout << "." << flush;
  cout << " " << stop_event_use-start_event_use << "(" << 100 << "%)";
  cout << endl;
}

// loop Lambda Mixed Event
void StStrangenessAna::MakeLambdaME()
{
  Long64_t start_event_use = mStartEvent_ME;
  Long64_t stop_event_use  = mStopEvent_ME;

  mInPut_ME->SetBranchAddress( XUV0_EVENT_BRANCH, &mLambda_event);
  mInPut_ME->GetEntry(0); // For unknown reasons root doesn't like it if someone starts to read a file not from the 0 entry

  // Initialise Event Head
  StThreeVectorF PrimaryVertex(0.0,0.0,0.0);
  Int_t          RunId = 0;
  Int_t          RefMult = 0;
  Int_t          Centrality = 0;
  Int_t          N_prim = 0;
  Int_t          N_non_prim = 0;
  Int_t          N_Tof_match = 0;
  Float_t        ZDCx = 0.0; 
  Float_t        BBCx = 0.0; 
  Float_t        VzVpd = 0.0;
  Int_t          NumTrackUsed = 0;
  // ---------------------------------------QVector---------------------------------------------
  TVector2       Q2East[4];
  TVector2       Q2West[4];
  TVector2       Q3East[4];
  TVector2       Q3West[4];
  // -----------------------------------Number of Tracks----------------------------------------
  Int_t          NumTrackEast[4];
  Int_t          NumTrackWest[4];
  for(Int_t j = 0; j < 4; j++)
  {
    Q2East[j].Set(0.0,0.0);
    Q2West[j].Set(0.0,0.0);
    Q3East[j].Set(0.0,0.0);
    Q3West[j].Set(0.0,0.0);
    NumTrackEast[j] = 0;
    NumTrackWest[j] = 0;
  }

  for(Long64_t counter = start_event_use; counter < stop_event_use; counter++)
  {
    if (!mInPut_ME->GetEntry( counter )) // take the event -> information is stored in event
      break;  // end of data chunk

    // get Event Header
    PrimaryVertex     = mLambda_event->getPrimaryVertex();
    RunId             = mLambda_event->getRunId();
    RefMult           = mLambda_event->getRefMult();
    Centrality        = mLambda_event->getCentrality();
    N_prim            = mLambda_event->getN_prim();
    N_non_prim        = mLambda_event->getN_non_prim();
    N_Tof_match       = mLambda_event->getN_Tof_match();
    ZDCx              = mLambda_event->getZDCx(); 
    BBCx              = mLambda_event->getBBCx(); 
    VzVpd             = mLambda_event->getVzVpd();
    NumTrackUsed      = mLambda_event->getNumTracks();

    for(Int_t j = 0; j < 4; j++)
    {
      Q2East[j]       = mLambda_event->getQ2East(j);
      Q2West[j]       = mLambda_event->getQ2West(j);
      Q3East[j]       = mLambda_event->getQ3East(j);
      Q3West[j]       = mLambda_event->getQ3West(j);
      NumTrackEast[j] = mLambda_event->getNumTrackEast(j);
      NumTrackWest[j] = mLambda_event->getNumTrackWest(j);
    }
 
    // Initialise Track 
    Float_t m2A = -999.9;
    Float_t m2B = -999.9;
    Float_t nsA = -999.9;
    Float_t nsB = -999.9;
    Float_t dcaA = -999.9;
    Float_t dcaB = -999.9;
    TLorentzVector lGTrackA(0.0,0.0,0.0,0.0);
    TLorentzVector lGTrackB(0.0,0.0,0.0,0.0);
    TLorentzVector lPTrackA(0.0,0.0,0.0,0.0);
    TLorentzVector lPTrackB(0.0,0.0,0.0,0.0);
    Int_t flagA = -1;
    Int_t flagB = -1;
    Float_t dcaAB = -999.9;
    Float_t decaylength = -999.9; 
    Float_t dcaV0 = -999.9;

    // vz sign
    Int_t vz_sign;
    if(PrimaryVertex.z() > 0.0)
    {
      vz_sign = 0;
    }
    else
    {
      vz_sign = 1;
    }

    // Centrality
    mRefMultCorr->init(RunId);
    if(mEnergy == 0) mRefMultCorr->initEvent(RefMult,PrimaryVertex.z(),ZDCx);
    if(mEnergy != 0) mRefMultCorr->initEvent(RefMult,PrimaryVertex.z());
    const Int_t cent9 = Centrality;
    const Double_t reweight = mRefMultCorr->getWeight();

    // runIndex
    mRunIdEventsDb = StRunIdEventsDb::Instance(Strangeness::mBeamEnergy[mEnergy],Strangeness::mBeamYear[mEnergy]);
    const Int_t runIndex = mRunIdEventsDb->getRunIdIndex(RunId); // expensive
    // cout << runIndex << endl;

    if (counter != 0  &&  counter % 1000 == 0)
      cout << "." << flush;
    if (counter != 0  &&  counter % 10000 == 0)
    {
      if((stop_event_use-start_event_use) > 0)
      {
	Double_t event_percent = 100.0*((Double_t)(counter-start_event_use))/((Double_t)(stop_event_use-start_event_use));
	cout << " " << counter-start_event_use << " (" << event_percent << "%) " << "\n" << "==> Processing data (strangeness_flow) " << flush;
      }
    }

    // get Track Information
    for(Int_t j = 0; j < Strangeness::mEtaGap_total; j++)
    {
      if(mStrangenessCorr->passTrackNumCut(NumTrackEast[j],NumTrackWest[j]))
      {
	for(UShort_t nTracks = 0; nTracks < NumTrackUsed; nTracks++) // loop over all tracks of the actual event
	{
	  mLambda_track = mLambda_event->getTrack(nTracks);
	  m2A = mLambda_track->getMass2A();
	  m2B = mLambda_track->getMass2B();
	  nsA = mLambda_track->getNSigA();
	  nsB = mLambda_track->getNSigB();
	  dcaA = mLambda_track->getDcaA();
	  dcaB = mLambda_track->getDcaB();
	  lGTrackA = mLambda_track->getGTrackA();
	  lGTrackB = mLambda_track->getGTrackB();
	  lPTrackA = mLambda_track->getPTrackA();
	  lPTrackB = mLambda_track->getPTrackB();
	  flagA = mLambda_track->getFlagA();
	  flagB = mLambda_track->getFlagB();
	  dcaAB = mLambda_track->getDcaAB();
	  decaylength = mLambda_track->getDecayLength(); 
	  dcaV0 = mLambda_track->getDcaV0();

	  Float_t GpA = lGTrackA.P();
	  Float_t GpB = lGTrackB.P();
	  TLorentzVector lGTrack = lGTrackA + lGTrackB;
	  Float_t pt_lGTrack = lGTrack.Perp();

	  // apply additional PID cut to increase significance
	  // final mass2 cut
	  Float_t Mass2_low_proton = 0.7;
	  Float_t Mass2_up_proton  = 1.1;
	  Float_t Mass2_low_pi_minus;
	  Float_t Mass2_up_pi_minus;
	  if(GpB < 0.75)
	  {
	    Mass2_low_pi_minus = -0.013;
	    Mass2_up_pi_minus  =  0.047;
	  }
	  if(GpB >= 0.75)
	  {
	    Mass2_low_pi_minus =  0.053 - 0.088*GpB;
	    Mass2_up_pi_minus  = -0.004 + 0.068*GpB;
	  }

	  if(
	      ( // no particle has a mass
		(m2A < -10.0) && (m2B < -10.0)
		&& fabs(dcaA)    > 0.6
		&& fabs(dcaB)    > 1.7
		&& dcaAB         < 1.0
		&& decaylength   > 4.0
		&& dcaV0         < 0.75
	      )
	      ||
	      ( // only pion has a mass
		(m2A < -10) && (m2B > Mass2_low_pi_minus && m2B < Mass2_up_pi_minus)
		&& fabs(dcaA)    > 0.5
		&& fabs(dcaB)    > 1.5
		&& dcaAB         < 1.0
		&& decaylength   > 3.5
		&& dcaV0         < 0.75
	      )
	      ||
	      ( // only proton has a mass
		(m2B < -10) && (m2A > Mass2_low_proton && m2A < Mass2_up_proton)
		&& fabs(dcaA)    > 0.15
		&& fabs(dcaB)    > 0.8
		&& dcaAB         < 1.0
		&& decaylength   > 2.5
		&& dcaV0         < 1.2
	      )
	      ||
	      ( // both particles have a mass
		(m2A > Mass2_low_proton && m2A < Mass2_up_proton) && (m2B > Mass2_low_pi_minus && m2B < Mass2_up_pi_minus)
		&& fabs(dcaA)    > 0.1
		&& fabs(dcaB)    > 0.7
		&& dcaAB         < 1.0
		&& decaylength   > 2.0
		&& dcaV0         < 1.3
	      )
	    )
	  {
	    Float_t phi_lGTrack = lGTrack.Phi();
	    Float_t InvMass_lGTrack = lGTrack.M();

	    if(mStrangenessCut->passPhiEtaEast(lGTrack)) // neg eta(east)
	    { // Below is West Only
	      TVector2 Q2Vector = Q2West[j];
	      TVector2 Q3Vector = Q3West[j];
	      // subtract auto-correlation from pos eta(west) event plane
	      if(flagA == 0 && mStrangenessCut->passTrackEP(lPTrackA,dcaA) && mStrangenessCut->passTrackEtaWest(lPTrackA,j)) // trackA
	      {
		Float_t  w = mStrangenessCorr->getWeight(lPTrackA);

		TVector2 q2VectorA = mStrangenessCorr->calq2Vector(lPTrackA);
		TVector2 q2CorrA   = mStrangenessCorr->getReCenterPar_West(0,cent9,runIndex,vz_sign,j,0); // 2nd
		Q2Vector = Q2Vector - w*(q2VectorA-q2CorrA);

		TVector2 q3VectorA = mStrangenessCorr->calq3Vector(lPTrackA);
		TVector2 q3CorrA   = mStrangenessCorr->getReCenterPar_West(1,cent9,runIndex,vz_sign,j,0); // 3rd
		Q3Vector = Q3Vector - w*(q3VectorA-q3CorrA);
	      }
	      if(flagB == 0 && mStrangenessCut->passTrackEP(lPTrackB,dcaB) && mStrangenessCut->passTrackEtaWest(lPTrackB,j)) // trackB
	      {
		Float_t  w = mStrangenessCorr->getWeight(lPTrackB);

		TVector2 q2VectorB = mStrangenessCorr->calq2Vector(lPTrackB);
		TVector2 q2CorrB   = mStrangenessCorr->getReCenterPar_West(0,cent9,runIndex,vz_sign,j,0); // 2nd
		Q2Vector = Q2Vector - w*(q2VectorB-q2CorrB);

		TVector2 q3VectorB = mStrangenessCorr->calq3Vector(lPTrackB);
		TVector2 q3CorrB   = mStrangenessCorr->getReCenterPar_West(1,cent9,runIndex,vz_sign,j,0); // 3rd
		Q3Vector = Q3Vector - w*(q3VectorB-q3CorrB);
	      }
	      Float_t Res2 = mStrangenessCorr->getResolution2_EP(cent9,j);
	      Float_t Psi2_west = mStrangenessCorr->calShiftAngle2West_EP(Q2Vector,runIndex,cent9,vz_sign,j);
	      Float_t Res3 = mStrangenessCorr->getResolution3_EP(cent9,j);
	      Float_t Psi3_west = mStrangenessCorr->calShiftAngle3West_EP(Q3Vector,runIndex,cent9,vz_sign,j);
	      Float_t phi_Psi2 = phi_lGTrack - Psi2_west;
	      Float_t phi_Psi3 = phi_lGTrack - Psi3_west;

	      mStrangenessHistoManger->Fill(pt_lGTrack,cent9,j,phi_Psi2,Res2,phi_Psi3,Res3,InvMass_lGTrack,reweight,n_cuts);

	      TLorentzVector lGTrackC;
	      lGTrackC.SetXYZM(lGTrackA.X(),lGTrackA.Y(),lGTrackA.Z(),0.49368);
	      TLorentzVector lGTrackCB = lGTrackC + lGTrackB; // misidentify K* -> Lambda
	      Double_t InvMassCB = lGTrackCB.M(); // Inv Mass of K*
	      if(!(InvMassCB > 0.89594-3*0.0487 && InvMassCB < 0.89594+3*0.0487))
	      {
		mStrangenessHistoManger->Fill_sub(pt_lGTrack,cent9,j,phi_Psi2,Res2,phi_Psi3,Res3,InvMass_lGTrack,reweight,n_cuts);
	      }
	    }

	    if(mStrangenessCut->passPhiEtaWest(lGTrack)) // pos eta
	    { // Below is East Only
	      TVector2 Q2Vector = Q2East[j];
	      TVector2 Q3Vector = Q3East[j];
	      // subtract auto-correlation from neg eta(east) event plane
	      if(flagA == 0 && mStrangenessCut->passTrackEP(lPTrackA,dcaA) && mStrangenessCut->passTrackEtaEast(lPTrackA,j)) // trackA
	      {
		Float_t  w = mStrangenessCorr->getWeight(lPTrackA);

		TVector2 q2VectorA = mStrangenessCorr->calq2Vector(lPTrackA);
		TVector2 q2CorrA   = mStrangenessCorr->getReCenterPar_East(0,cent9,runIndex,vz_sign,j,0); // 2nd
		Q2Vector = Q2Vector - w*(q2VectorA-q2CorrA);

		TVector2 q3VectorA = mStrangenessCorr->calq3Vector(lPTrackA);
		TVector2 q3CorrA   = mStrangenessCorr->getReCenterPar_East(1,cent9,runIndex,vz_sign,j,0); // 3rd
		Q3Vector = Q3Vector - w*(q3VectorA-q3CorrA);
	      }
	      if(flagB == 0 && mStrangenessCut->passTrackEP(lPTrackB,dcaB) && mStrangenessCut->passTrackEtaEast(lPTrackB,j)) // trackB
	      {
		Float_t  w = mStrangenessCorr->getWeight(lPTrackB);

		TVector2 q2VectorB = mStrangenessCorr->calq2Vector(lPTrackB);
		TVector2 q2CorrB   = mStrangenessCorr->getReCenterPar_East(0,cent9,runIndex,vz_sign,j,0); // 2nd
		Q2Vector = Q2Vector - w*(q2VectorB-q2CorrB);

		TVector2 q3VectorB = mStrangenessCorr->calq3Vector(lPTrackB);
		TVector2 q3CorrB   = mStrangenessCorr->getReCenterPar_East(1,cent9,runIndex,vz_sign,j,0); // 3rd
		Q3Vector = Q3Vector - w*(q3VectorB-q3CorrB);
	      }
	      Float_t Res2 = mStrangenessCorr->getResolution2_EP(cent9,j);
	      Float_t Psi2_east = mStrangenessCorr->calShiftAngle2East_EP(Q2Vector,runIndex,cent9,vz_sign,j);
	      Float_t Res3 = mStrangenessCorr->getResolution3_EP(cent9,j);
	      Float_t Psi3_east = mStrangenessCorr->calShiftAngle3East_EP(Q3Vector,runIndex,cent9,vz_sign,j);
	      Float_t phi_Psi2 = phi_lGTrack - Psi2_east;
	      Float_t phi_Psi3 = phi_lGTrack - Psi3_east;

	      mStrangenessHistoManger->Fill(pt_lGTrack,cent9,j,phi_Psi2,Res2,phi_Psi3,Res3,InvMass_lGTrack,reweight,n_cuts);

	      TLorentzVector lGTrackC;
	      lGTrackC.SetXYZM(lGTrackA.X(),lGTrackA.Y(),lGTrackA.Z(),0.49368);
	      TLorentzVector lGTrackCB = lGTrackC + lGTrackB; // misidentify K* -> Lambda
	      Double_t InvMassCB = lGTrackCB.M(); // Inv Mass of K*
	      if(!(InvMassCB > 0.89594-3*0.0487 && InvMassCB < 0.89594+3*0.0487))
	      {
		mStrangenessHistoManger->Fill_sub(pt_lGTrack,cent9,j,phi_Psi2,Res2,phi_Psi3,Res3,InvMass_lGTrack,reweight,n_cuts);
	      }
	    }
	  }
	}
      }
    }
  }

  cout << "." << flush;
  cout << " " << stop_event_use-start_event_use << "(" << 100 << "%)";
  cout << endl;
}

// K0S loop for Same Event and Mixed Event
void StStrangenessAna::MakeK0S(Int_t X_flag)
{
  Long64_t start_event_use;
  Long64_t stop_event_use;
  if(X_flag == 0)
  {
    start_event_use = mStartEvent_SE;
    stop_event_use  = mStopEvent_SE;
    mInPut_SE->SetBranchAddress( XUV0_EVENT_BRANCH, &mK0S_event);
    mInPut_SE->GetEntry(0); // For unknown reasons root doesn't like it if someone starts to read a file not from the 0 entry
  }
  if(X_flag == 1)
  {
    start_event_use = mStartEvent_ME;
    stop_event_use  = mStopEvent_ME;
    mInPut_ME->SetBranchAddress( XUV0_EVENT_BRANCH, &mK0S_event);
    mInPut_ME->GetEntry(0); // For unknown reasons root doesn't like it if someone starts to read a file not from the 0 entry
  }

  // Initialise Event Head
  StThreeVectorF PrimaryVertex(0.0,0.0,0.0);
  Int_t          RunId = 0;
  Int_t          RefMult = 0;
  Int_t          Centrality = 0;
  Int_t          N_prim = 0;
  Int_t          N_non_prim = 0;
  Int_t          N_Tof_match = 0;
  Float_t        ZDCx = 0.0; 
  Float_t        BBCx = 0.0; 
  Float_t        VzVpd = 0.0;
  Int_t          NumTrackUsed = 0;
  // ---------------------------------------QVector---------------------------------------------
  TVector2       Q2East[4];
  TVector2       Q2West[4];
  TVector2       Q3East[4];
  TVector2       Q3West[4];
  // -----------------------------------Number of Tracks----------------------------------------
  Int_t          NumTrackEast[4];
  Int_t          NumTrackWest[4];
  for(Int_t j = 0; j < 4; j++)
  {
    Q2East[j].Set(0.0,0.0);
    Q2West[j].Set(0.0,0.0);
    Q3East[j].Set(0.0,0.0);
    Q3West[j].Set(0.0,0.0);
    NumTrackEast[j] = 0;
    NumTrackWest[j] = 0;
  }

  for(Long64_t counter = start_event_use; counter < stop_event_use; counter++)
  {
    if (X_flag == 0 && !mInPut_SE->GetEntry( counter )) // take the event -> information is stored in event
      break;  // end of data chunk

    if (X_flag == 1 && !mInPut_ME->GetEntry( counter )) // take the event -> information is stored in event
      break;  // end of data chunk

    // get Event Header
    PrimaryVertex     = mK0S_event->getPrimaryVertex();
    RunId             = mK0S_event->getRunId();
    RefMult           = mK0S_event->getRefMult();
    Centrality        = mK0S_event->getCentrality();
    N_prim            = mK0S_event->getN_prim();
    N_non_prim        = mK0S_event->getN_non_prim();
    N_Tof_match       = mK0S_event->getN_Tof_match();
    ZDCx              = mK0S_event->getZDCx(); 
    BBCx              = mK0S_event->getBBCx(); 
    VzVpd             = mK0S_event->getVzVpd();
    NumTrackUsed      = mK0S_event->getNumTracks();

    for(Int_t j = 0; j < 4; j++)
    {
      Q2East[j]       = mK0S_event->getQ2East(j);
      Q2West[j]       = mK0S_event->getQ2West(j);
      Q3East[j]       = mK0S_event->getQ3East(j);
      Q3West[j]       = mK0S_event->getQ3West(j);
      NumTrackEast[j] = mK0S_event->getNumTrackEast(j);
      NumTrackWest[j] = mK0S_event->getNumTrackWest(j);
    }
 
    // Initialise Track 
    Float_t m2A = -999.9;
    Float_t m2B = -999.9;
    Float_t nsA = -999.9;
    Float_t nsB = -999.9;
    Float_t dcaA = -999.9;
    Float_t dcaB = -999.9;
    TLorentzVector lGTrackA(0.0,0.0,0.0,0.0);
    TLorentzVector lGTrackB(0.0,0.0,0.0,0.0);
    TLorentzVector lPTrackA(0.0,0.0,0.0,0.0);
    TLorentzVector lPTrackB(0.0,0.0,0.0,0.0);
    Int_t flagA = -1;
    Int_t flagB = -1;
    Float_t dcaAB = -999.9;
    Float_t decaylength = -999.9; 
    Float_t dcaV0 = -999.9;

    // vz sign
    Int_t vz_sign;
    if(PrimaryVertex.z() > 0.0)
    {
      vz_sign = 0;
    }
    else
    {
      vz_sign = 1;
    }

    // Centrality
    mRefMultCorr->init(RunId);
    if(mEnergy == 0) mRefMultCorr->initEvent(RefMult,PrimaryVertex.z(),ZDCx);
    if(mEnergy != 0) mRefMultCorr->initEvent(RefMult,PrimaryVertex.z());
    const Int_t cent9 = Centrality;
    const Double_t reweight = mRefMultCorr->getWeight();

    // runIndex
    mRunIdEventsDb = StRunIdEventsDb::Instance(Strangeness::mBeamEnergy[mEnergy],Strangeness::mBeamYear[mEnergy]);
    const Int_t runIndex = mRunIdEventsDb->getRunIdIndex(RunId); // expensive
    //cout << runIndex << endl;

    if (counter != 0  &&  counter % 1000 == 0)
      cout << "." << flush;
    if (counter != 0  &&  counter % 10000 == 0)
    {
      if((stop_event_use-start_event_use) > 0)
      {
	Double_t event_percent = 100.0*((Double_t)(counter-start_event_use))/((Double_t)(stop_event_use-start_event_use));
	cout << " " << counter-start_event_use << " (" << event_percent << "%) " << "\n" << "==> Processing data (strangeness_flow) " << flush;
      }
    }

    // get Track Information
    for(Int_t j = 0; j < Strangeness::mEtaGap_total; j++)
    {
      if(mStrangenessCorr->passTrackNumCut(NumTrackEast[j],NumTrackWest[j]))
      {
	for(UShort_t nTracks = 0; nTracks < NumTrackUsed; nTracks++) // loop over all tracks of the actual event
	{
	  mK0S_track = mK0S_event->getTrack(nTracks);
	  m2A = mK0S_track->getMass2A();
	  m2B = mK0S_track->getMass2B();
	  nsA = mK0S_track->getNSigA();
	  nsB = mK0S_track->getNSigB();
	  dcaA = mK0S_track->getDcaA();
	  dcaB = mK0S_track->getDcaB();
	  lGTrackA = mK0S_track->getGTrackA();
	  lGTrackB = mK0S_track->getGTrackB();
	  lPTrackA = mK0S_track->getPTrackA();
	  lPTrackB = mK0S_track->getPTrackB();
	  flagA = mK0S_track->getFlagA();
	  flagB = mK0S_track->getFlagB();
	  dcaAB = mK0S_track->getDcaAB();
	  decaylength = mK0S_track->getDecayLength(); 
	  dcaV0 = mK0S_track->getDcaV0();

	  Float_t GpA = lGTrackA.P();
	  Float_t GpB = lGTrackB.P();
	  TLorentzVector lGTrack = lGTrackA + lGTrackB; // Lorentz vector for mother particle
	  Float_t pt_lGTrack = lGTrack.Perp();

	  // apply additional PID cut to increase significance
	  // final mass2 cut
	  Float_t Mass2_low_pi_plus;
	  Float_t Mass2_up_pi_plus;
	  if(GpA < 0.75)
	  {
	    Mass2_low_pi_plus = -0.013;
	    Mass2_up_pi_plus =  0.047;
	  }
	  if(GpA >= 0.75)
	  {
	    Mass2_low_pi_plus =  0.053 - 0.088*GpB;
	    Mass2_up_pi_plus = -0.004 + 0.068*GpB;
	  }

	  Float_t Mass2_low_pi_minus;
	  Float_t Mass2_up_pi_minus;
	  if(GpB < 0.75)
	  {
	    Mass2_low_pi_minus = -0.013;
	    Mass2_up_pi_minus  =  0.047;
	  }
	  if(GpB >= 0.75)
	  {
	    Mass2_low_pi_minus =  0.053 - 0.088*GpB;
	    Mass2_up_pi_minus  = -0.004 + 0.068*GpB;
	  }

	  if(
	       fabs(dcaA)    > 0.7
	    && fabs(dcaB)    > 0.7
	    && dcaAB         < 1.0
	    && decaylength   > 3.0
	    && dcaV0         < 0.8
	    && ((m2A > Mass2_low_pi_plus && m2A < Mass2_up_pi_plus) || (m2A < -10))
	    && ((m2B > Mass2_low_pi_minus && m2B < Mass2_up_pi_minus) || (m2B < -10))
	    )
	  {
	    Float_t phi_lGTrack = lGTrack.Phi();
	    Float_t InvMass_lGTrack = lGTrack.M();

	    if(mStrangenessCut->passPhiEtaEast(lGTrack)) // neg eta(east)
	    { // Below is West Only
	      TVector2 Q2Vector = Q2West[j];
	      TVector2 Q3Vector = Q3West[j];
	      // subtract auto-correlation from pos eta(west) event plane
	      if(flagA == 0 && mStrangenessCut->passTrackEP(lPTrackA,dcaA) && mStrangenessCut->passTrackEtaWest(lPTrackA,j)) // trackA
	      {
		Float_t  w = mStrangenessCorr->getWeight(lPTrackA);

		TVector2 q2VectorA = mStrangenessCorr->calq2Vector(lPTrackA);
		TVector2 q2CorrA   = mStrangenessCorr->getReCenterPar_West(0,cent9,runIndex,vz_sign,j,0); // 2nd
		Q2Vector = Q2Vector - w*(q2VectorA-q2CorrA);

		TVector2 q3VectorA = mStrangenessCorr->calq3Vector(lPTrackA);
		TVector2 q3CorrA   = mStrangenessCorr->getReCenterPar_West(1,cent9,runIndex,vz_sign,j,0); // 3rd
		Q3Vector = Q3Vector - w*(q3VectorA-q3CorrA);
	      }
	      if(flagB == 0 && mStrangenessCut->passTrackEP(lPTrackB,dcaB) && mStrangenessCut->passTrackEtaWest(lPTrackB,j)) // trackB
	      {
		Float_t  w = mStrangenessCorr->getWeight(lPTrackB);

		TVector2 q2VectorB = mStrangenessCorr->calq2Vector(lPTrackB);
		TVector2 q2CorrB   = mStrangenessCorr->getReCenterPar_West(0,cent9,runIndex,vz_sign,j,0); // 2nd
		Q2Vector = Q2Vector - w*(q2VectorB-q2CorrB);

		TVector2 q3VectorB = mStrangenessCorr->calq3Vector(lPTrackB);
		TVector2 q3CorrB   = mStrangenessCorr->getReCenterPar_West(1,cent9,runIndex,vz_sign,j,0); // 3rd
		Q3Vector = Q3Vector - w*(q3VectorB-q3CorrB);
	      }
	      Float_t Res2 = mStrangenessCorr->getResolution2_EP(cent9,j);
	      Float_t Psi2_west = mStrangenessCorr->calShiftAngle2West_EP(Q2Vector,runIndex,cent9,vz_sign,j);
	      Float_t Res3 = mStrangenessCorr->getResolution3_EP(cent9,j);
	      Float_t Psi3_west = mStrangenessCorr->calShiftAngle3West_EP(Q3Vector,runIndex,cent9,vz_sign,j);
	      Float_t phi_Psi2 = phi_lGTrack - Psi2_west;
	      Float_t phi_Psi3 = phi_lGTrack - Psi3_west;

	      mStrangenessHistoManger->Fill(pt_lGTrack,cent9,j,phi_Psi2,Res2,phi_Psi3,Res3,InvMass_lGTrack,reweight,n_cuts);

	      //-----------------------------------------------------------------------------
	      // get Lambda by misidentification
	      TLorentzVector  lGTrackP, lGTrackPbar;
	      lGTrackP.SetXYZM(lGTrackA.X(),lGTrackA.Y(),lGTrackA.Z(),0.93827); // set Lorentz vector of pi+ to p
	      lGTrackPbar.SetXYZM(lGTrackB.X(),lGTrackB.Y(),lGTrackB.Z(),0.93827); // set Lorentz vector of pi- to pbar

	      // Invariant mass calculations K0s
	      TLorentzVector trackLambda      = lGTrackP+lGTrackB; // p + pi-
	      TLorentzVector trackAntiLambda      = lGTrackPbar+lGTrackA; // pbar + pi+
	      Double_t InvMassLambda = trackLambda.M(); // invariant mass of Lambda
	      Double_t InvMassAntiLambda = trackAntiLambda.M(); // invariant mass of Lambda
	      //-----------------------------------------------------------------------------

	      if(!((InvMassLambda > 1.1157-0.006*3 && InvMassLambda < 1.1157+0.006*3) || (InvMassAntiLambda > 1.1157-0.006*3 && InvMassAntiLambda < 1.1157+0.006*3)))
	      {
		mStrangenessHistoManger->Fill_sub(pt_lGTrack,cent9,j,phi_Psi2,Res2,phi_Psi3,Res3,InvMass_lGTrack,reweight,n_cuts);
	      }
	    }

	    if(mStrangenessCut->passPhiEtaWest(lGTrack)) // pos eta
	    { // Below is East Only
	      TVector2 Q2Vector = Q2East[j];
	      TVector2 Q3Vector = Q3East[j];
	      // subtract auto-correlation from neg eta(east) event plane
	      if(flagA == 0 && mStrangenessCut->passTrackEP(lPTrackA,dcaA) && mStrangenessCut->passTrackEtaEast(lPTrackA,j)) // trackA
	      {
		Float_t  w = mStrangenessCorr->getWeight(lPTrackA);

		TVector2 q2VectorA = mStrangenessCorr->calq2Vector(lPTrackA);
		TVector2 q2CorrA   = mStrangenessCorr->getReCenterPar_East(0,cent9,runIndex,vz_sign,j,0); // 2nd
		Q2Vector = Q2Vector - w*(q2VectorA-q2CorrA);

		TVector2 q3VectorA = mStrangenessCorr->calq3Vector(lPTrackA);
		TVector2 q3CorrA   = mStrangenessCorr->getReCenterPar_East(1,cent9,runIndex,vz_sign,j,0); // 3rd
		Q3Vector = Q3Vector - w*(q3VectorA-q3CorrA);
	      }
	      if(flagB == 0 && mStrangenessCut->passTrackEP(lPTrackB,dcaB) && mStrangenessCut->passTrackEtaEast(lPTrackB,j)) // trackB
	      {
		Float_t  w = mStrangenessCorr->getWeight(lPTrackB);

		TVector2 q2VectorB = mStrangenessCorr->calq2Vector(lPTrackB);
		TVector2 q2CorrB   = mStrangenessCorr->getReCenterPar_East(0,cent9,runIndex,vz_sign,j,0); // 2nd
		Q2Vector = Q2Vector - w*(q2VectorB-q2CorrB);

		TVector2 q3VectorB = mStrangenessCorr->calq3Vector(lPTrackB);
		TVector2 q3CorrB   = mStrangenessCorr->getReCenterPar_East(1,cent9,runIndex,vz_sign,j,0); // 3rd
		Q3Vector = Q3Vector - w*(q3VectorB-q3CorrB);
	      }
	      Float_t Res2 = mStrangenessCorr->getResolution2_EP(cent9,j);
	      Float_t Psi2_east = mStrangenessCorr->calShiftAngle2East_EP(Q2Vector,runIndex,cent9,vz_sign,j);
	      Float_t Res3 = mStrangenessCorr->getResolution3_EP(cent9,j);
	      Float_t Psi3_east = mStrangenessCorr->calShiftAngle3East_EP(Q3Vector,runIndex,cent9,vz_sign,j);
	      Float_t phi_Psi2 = phi_lGTrack - Psi2_east;
	      Float_t phi_Psi3 = phi_lGTrack - Psi3_east;

	      mStrangenessHistoManger->Fill(pt_lGTrack,cent9,j,phi_Psi2,Res2,phi_Psi3,Res3,InvMass_lGTrack,reweight,n_cuts);

	      //-----------------------------------------------------------------------------
	      // get Lambda by misidentification
	      TLorentzVector  lGTrackP, lGTrackPbar;
	      lGTrackP.SetXYZM(lGTrackA.X(),lGTrackA.Y(),lGTrackA.Z(),0.93827); // set Lorentz vector of pi+ to p
	      lGTrackPbar.SetXYZM(lGTrackB.X(),lGTrackB.Y(),lGTrackB.Z(),0.93827); // set Lorentz vector of pi- to pbar

	      // Invariant mass calculations K0s
	      TLorentzVector trackLambda      = lGTrackP+lGTrackB; // p + pi-
	      TLorentzVector trackAntiLambda      = lGTrackPbar+lGTrackA; // pbar + pi+
	      Double_t InvMassLambda = trackLambda.M(); // invariant mass of Lambda
	      Double_t InvMassAntiLambda = trackAntiLambda.M(); // invariant mass of Lambda
	      //-----------------------------------------------------------------------------

	      // subtract Lambda and anti-Lambda contaminations
	      if(!((InvMassLambda > 1.1157-0.006*3 && InvMassLambda < 1.1157+0.006*3) || (InvMassAntiLambda > 1.1157-0.006*3 && InvMassAntiLambda < 1.1157+0.006*3)))
	      {
		mStrangenessHistoManger->Fill_sub(pt_lGTrack,cent9,j,phi_Psi2,Res2,phi_Psi3,Res3,InvMass_lGTrack,reweight,n_cuts);
	      }
	    }
	  }
	}
      }
    }
  }

  cout << "." << flush;
  cout << " " << stop_event_use-start_event_use << "(" << 100 << "%)";
  cout << endl;
}
//-------------------------------------------------------------------
void StStrangenessAna::Finish()
{
  mFile_OutPut->cd();
  mStrangenessHistoManger->Write();
  mFile_OutPut->Close();
}

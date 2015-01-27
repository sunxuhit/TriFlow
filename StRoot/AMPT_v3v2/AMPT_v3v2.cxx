#include "AMPT_v3v2.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"
#include "TProfile.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TMath.h"

ClassImp(AMPT_v3v2)

Int_t AMPT_v3v2::mInput_flag = 1;
TString AMPT_v3v2::mBeamEnergy[7] = {"7GeV","11GeV","19GeV","27GeV","39GeV","62GeV","200GeV"};
TString AMPT_v3v2::mMode_AMPT[2] = {"Default","StringMelting"};
TString AMPT_v3v2::mScreenMass_AMPT[2] = {"3mb","6mb"};
Int_t AMPT_v3v2::mRefMult[2][7][10] = {
					   { // Default
					     {9,16,27,42,63,92,130,182,219,329}, //  7GeV
					     {10,18,31,49,74,106,148,207,246,352}, // 11GeV
					     {12,21,35,56,83,119,166,229,270,385}, // 19GeV
					     {13,23,39,60,90,128,178,246,289,422}, // 27GeV
					     {14,26,43,67,99,141,196,270,318,478}, // 39GeV | with phi decay
					     {16,28,48,75,112,161,225,312,368,521}, // 62GeV
					     {23,42,72,115,174,254,361,507,603,900}  // 200GeV
					   },
					   { // String Melting
					     {9,16,27,41,61,86,119,163,193,291}, //  7GeV
					     {11,19,32,51,75,106,147,201,237,345}, // 11GeV
					     {12,22,38,59,88,126,175,240,281,414}, // 19GeV
					     {13,24,40,64,95,137,190,261,307,443}, // 27GeV
//					     {15,26,44,70,105,150,209,288,338,510}, // 39GeV | with phi decay | 3 mb
					     {15,26,43,68,101,144,201,277,325,492}, // 39GeV | with phi decay | 6 mb
					     {16,29,49,77,117,169,238,331,390,628},  // 62GeV
//					     {22,41,72,118,182,269,387,550,657,1043}  // 200GeV | with phi decay | 3 mb
					     {21,42,72,117,181,269,388,551,659,978}  // 200GeV | with phi decay | 3 mb v1.21/v2.21
//					     {22,41,71,114,176,261,376,536,641,1088}  // 200GeV | with phi decay | 6 mb
					   }
					 }; // 80%,70%,60%,50%,40%,30%,20%,10%,5%,0%
// Centrality bin
Int_t AMPT_v3v2::cent_low[4] = {0,7,4,0}; // 0 = 0-80%, 1 = 0-10%, 2 = 10-40%, 3 = 40-80%
Int_t AMPT_v3v2::cent_up[4]  = {8,8,6,3}; // 0 = 0-80%, 1 = 0-10%, 2 = 10-40%, 3 = 40-80%
Int_t AMPT_v3v2::mList_start[25] = {  1,101,201,301,401,501,601,701,801, 901,1001,1101,1201,1301,1401,1501,1601,1701,1801,1901,2001,2101,2201,2301,2401};
Int_t AMPT_v3v2::mList_stop[25]  = {100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2100,2200,2300,2400,2500};
Int_t AMPT_v3v2::Centrality_start = 0;
Int_t AMPT_v3v2::Centrality_stop  = 4;
//------------------------------------------------------------
AMPT_v3v2::AMPT_v3v2(Int_t Energy, Int_t Mode, Int_t Screen, Int_t List, Long64_t StartEvent, Long64_t StopEvent)
{
  mEnergy = Energy;
  mMode = Mode;
  mScreen = Screen;

  TString InPutList;
  if(mMode == 0)
  {
    InPutList = Form("/project/projectdirs/star/xusun/OutPut/AMPT_%s/List/%s_List/Split_%s_%d_%d.list",mMode_AMPT[Mode].Data(),mBeamEnergy[mEnergy].Data(),mBeamEnergy[mEnergy].Data(),mList_start[List],mList_stop[List]);
  }
  if(mMode == 1)
  {
    InPutList = Form("/project/projectdirs/star/xusun/OutPut/AMPT_%s/List/%s_List/%s/Split_%s_%d_%d.list",mMode_AMPT[Mode].Data(),mBeamEnergy[mEnergy].Data(),mScreenMass_AMPT[mScreen].Data(),mBeamEnergy[mEnergy].Data(),mList_start[List],mList_stop[List]);
  }
  SetInPutList(InPutList); // set input list

  SetStartEvent(StartEvent); // set start event
  SetStopEvent(StopEvent); // set stop event

  TString InPutRes;
  if(mMode == 0)
  {
    InPutRes = Form("/project/projectdirs/star/xusun/OutPut/AMPT_%s/Resolution/%s_Resolution/Resolution_%s.root",mMode_AMPT[Mode].Data(),mBeamEnergy[Energy].Data(),mBeamEnergy[Energy].Data());
  }
  if(mMode == 1)
  {
    InPutRes = Form("/project/projectdirs/star/xusun/OutPut/AMPT_%s/Resolution/%s_Resolution/%s/Resolution_%s.root",mMode_AMPT[Mode].Data(),mBeamEnergy[Energy].Data(),mScreenMass_AMPT[mScreen].Data(),mBeamEnergy[Energy].Data());
  }
  SetInPutRes(InPutRes); // set input resolution

  TString OutPutFile;
  if(mMode == 0)
  {
    OutPutFile = Form("/project/projectdirs/star/xusun/OutPut/AMPT_%s/Flow/%s_%s/Flow_%s_%d_%d.root",mMode_AMPT[Mode].Data(),mBeamEnergy[mEnergy].Data(),mMode_AMPT[Mode].Data(),mBeamEnergy[mEnergy].Data(),mList_start[List],mList_stop[List]);
  }
  if(mMode == 1)
  {
    OutPutFile = Form("/project/projectdirs/star/xusun/OutPut/AMPT_%s/Flow/%s_%s/%s/Flow_%s_%d_%d.root",mMode_AMPT[Mode].Data(),mBeamEnergy[mEnergy].Data(),mMode_AMPT[Mode].Data(),mScreenMass_AMPT[mScreen].Data(),mBeamEnergy[mEnergy].Data(),mList_start[List],mList_stop[List]);
  }
  SetOutPutFile(OutPutFile); // set output file
}

AMPT_v3v2::~AMPT_v3v2()
{
}
//------------------------------------------------------------
void AMPT_v3v2::SetInPutList(const TString inputlist)
{
  mInPutList= inputlist.Copy();
  cout << "Input list was set to: " << mInPutList.Data() << endl;
}

void AMPT_v3v2::SetOutPutFile(const TString outputfile)
{
  mOutPutFile = outputfile.Copy();
  cout << "Output file was set to: " << mOutPutFile.Data() << endl;
}

void AMPT_v3v2::SetInPutRes(const TString inputres)
{
  mInPutRes = inputres.Copy();
  cout << "Input resolution file was set to: " << mInPutRes.Data() << endl;
}

void AMPT_v3v2::SetStartEvent(const Long64_t StartEvent)
{
  mStartEvent = StartEvent;
  cout << "nStartEvent = " << mStartEvent << endl;
}

void AMPT_v3v2::SetStopEvent(const Long64_t StopEvent)
{
  mStopEvent = StopEvent;
  cout << "nStopEvent = " << mStopEvent << endl;
}
//------------------------------------------------------------
Float_t AMPT_v3v2::getResolution(Int_t order, Int_t i_cent) // 0 for 2nd, 1 for 3rd
{
  Float_t res = -999.9;
  if(p_mRes[order]->GetBinContent(p_mRes[order]->FindBin(i_cent)) > 0.0)
  {
    res = TMath::Sqrt(p_mRes[order]->GetBinContent(p_mRes[order]->FindBin(i_cent)));
  }
  return res;
}

Int_t AMPT_v3v2::getCentrality(Int_t refMult)
{
  Int_t cent9 = -1;

  // Centrality defination
  for(Int_t i_cent = 0; i_cent < 9; i_cent++)
  {
    if(refMult >= mRefMult[mMode][mEnergy][i_cent] && refMult < mRefMult[mMode][mEnergy][i_cent+1])
    {
      cent9 = i_cent;
    }
  }

  return cent9;
}
//------------------------------------------------------------
void AMPT_v3v2::Init()
{
  // initialize the resolution
  mFile_Res = TFile::Open(mInPutRes.Data());
  p_mRes[0] = (TProfile*)mFile_Res->Get("p_mRes2"); // 2nd event plane resolution
  p_mRes[1] = (TProfile*)mFile_Res->Get("p_mRes3"); // 3rd event plane resolution

  mFile_OutPut = new TFile(mOutPutFile.Data(),"RECREATE");

  // flow TProfile
  for(Int_t i_order = 0; i_order < 2; i_order++)
  {
    for(Int_t i_cent = 0; i_cent < 4; i_cent++)
    {
      TString Order[2] = {"2nd","3rd"};
      TString Centrality[4] = {"0080","0010","1040","4080"};
      TString ProName;
      TString HistName;

      // flow relative to event plane
      ProName = Form("Flow_pi_plus_%s_%s",Order[i_order].Data(),Centrality[i_cent].Data()); // pi_plus
      p_mFlow_pi_plus[i_order][i_cent] = new TProfile(ProName.Data(),ProName.Data(),25,0.1,5.1);
      HistName = Form("Pt_pi_plus_%s_%s",Order[i_order].Data(),Centrality[i_cent].Data());
      h_mPt_pi_plus[i_order][i_cent] = new TH1F(HistName.Data(),HistName.Data(),25,0.1,5.1);

      ProName = Form("Flow_pi_minus_%s_%s",Order[i_order].Data(),Centrality[i_cent].Data()); // pi_minus
      p_mFlow_pi_minus[i_order][i_cent] = new TProfile(ProName.Data(),ProName.Data(),25,0.1,5.1);
      HistName = Form("Pt_pi_minus_%s_%s",Order[i_order].Data(),Centrality[i_cent].Data());
      h_mPt_pi_minus[i_order][i_cent] = new TH1F(HistName.Data(),HistName.Data(),25,0.1,5.1);

      ProName = Form("Flow_K_plus_%s_%s",Order[i_order].Data(),Centrality[i_cent].Data()); // K_plus
      p_mFlow_K_plus[i_order][i_cent] = new TProfile(ProName.Data(),ProName.Data(),25,0.1,5.1);
      HistName = Form("Pt_K_plus_%s_%s",Order[i_order].Data(),Centrality[i_cent].Data());
      h_mPt_K_plus[i_order][i_cent] = new TH1F(HistName.Data(),HistName.Data(),25,0.1,5.1);

      ProName = Form("Flow_K_minus_%s_%s",Order[i_order].Data(),Centrality[i_cent].Data()); // K_minus
      p_mFlow_K_minus[i_order][i_cent] = new TProfile(ProName.Data(),ProName.Data(),25,0.1,5.1);
      HistName = Form("Pt_K_minus_%s_%s",Order[i_order].Data(),Centrality[i_cent].Data());
      h_mPt_K_minus[i_order][i_cent] = new TH1F(HistName.Data(),HistName.Data(),25,0.1,5.1);

      ProName = Form("Flow_p_%s_%s",Order[i_order].Data(),Centrality[i_cent].Data()); // p
      p_mFlow_p[i_order][i_cent] = new TProfile(ProName.Data(),ProName.Data(),25,0.1,5.1);
      HistName = Form("Pt_p_%s_%s",Order[i_order].Data(),Centrality[i_cent].Data());
      h_mPt_p[i_order][i_cent] = new TH1F(HistName.Data(),HistName.Data(),25,0.1,5.1);

      ProName = Form("Flow_pbar_%s_%s",Order[i_order].Data(),Centrality[i_cent].Data()); // pbar
      p_mFlow_pbar[i_order][i_cent] = new TProfile(ProName.Data(),ProName.Data(),25,0.1,5.1);
      HistName = Form("Pt_pbar_%s_%s",Order[i_order].Data(),Centrality[i_cent].Data());
      h_mPt_pbar[i_order][i_cent] = new TH1F(HistName.Data(),HistName.Data(),25,0.1,5.1);

      ProName = Form("Flow_Lambda_%s_%s",Order[i_order].Data(),Centrality[i_cent].Data()); // Lambda
      p_mFlow_Lambda[i_order][i_cent] = new TProfile(ProName.Data(),ProName.Data(),25,0.1,5.1);
      HistName = Form("Pt_Lambda_%s_%s",Order[i_order].Data(),Centrality[i_cent].Data());
      h_mPt_Lambda[i_order][i_cent] = new TH1F(HistName.Data(),HistName.Data(),25,0.1,5.1);

      ProName = Form("Flow_Lambdabar_%s_%s",Order[i_order].Data(),Centrality[i_cent].Data()); // Lambdabar
      p_mFlow_Lambdabar[i_order][i_cent] = new TProfile(ProName.Data(),ProName.Data(),25,0.1,5.1);
      HistName = Form("Pt_Lambdabar_%s_%s",Order[i_order].Data(),Centrality[i_cent].Data());
      h_mPt_Lambdabar[i_order][i_cent] = new TH1F(HistName.Data(),HistName.Data(),25,0.1,5.1);

      ProName = Form("Flow_K0s_%s_%s",Order[i_order].Data(),Centrality[i_cent].Data()); // K0s
      p_mFlow_K0s[i_order][i_cent] = new TProfile(ProName.Data(),ProName.Data(),25,0.1,5.1);
      HistName = Form("Pt_K0s_%s_%s",Order[i_order].Data(),Centrality[i_cent].Data());
      h_mPt_K0s[i_order][i_cent] = new TH1F(HistName.Data(),HistName.Data(),25,0.1,5.1);

      ProName = Form("Flow_phi_%s_%s",Order[i_order].Data(),Centrality[i_cent].Data()); // phi
      p_mFlow_phi[i_order][i_cent] = new TProfile(ProName.Data(),ProName.Data(),25,0.1,5.1);
      HistName = Form("Pt_phi_%s_%s",Order[i_order].Data(),Centrality[i_cent].Data());
      h_mPt_phi[i_order][i_cent] = new TH1F(HistName.Data(),HistName.Data(),25,0.1,5.1);

      // flow relative to Reaction Plane
      ProName = Form("Flow_pi_plus_%s_%s_RP",Order[i_order].Data(),Centrality[i_cent].Data()); // pi_plus
      p_mFlow_pi_plus_RP[i_order][i_cent] = new TProfile(ProName.Data(),ProName.Data(),25,0.1,5.1);

      ProName = Form("Flow_pi_minus_%s_%s_RP",Order[i_order].Data(),Centrality[i_cent].Data()); // pi_minus
      p_mFlow_pi_minus_RP[i_order][i_cent] = new TProfile(ProName.Data(),ProName.Data(),25,0.1,5.1);

      ProName = Form("Flow_K_plus_%s_%s_RP",Order[i_order].Data(),Centrality[i_cent].Data()); // K_plus
      p_mFlow_K_plus_RP[i_order][i_cent] = new TProfile(ProName.Data(),ProName.Data(),25,0.1,5.1);

      ProName = Form("Flow_K_minus_%s_%s_RP",Order[i_order].Data(),Centrality[i_cent].Data()); // K_minus
      p_mFlow_K_minus_RP[i_order][i_cent] = new TProfile(ProName.Data(),ProName.Data(),25,0.1,5.1);

      ProName = Form("Flow_p_%s_%s_RP",Order[i_order].Data(),Centrality[i_cent].Data()); // p
      p_mFlow_p_RP[i_order][i_cent] = new TProfile(ProName.Data(),ProName.Data(),25,0.1,5.1);

      ProName = Form("Flow_pbar_%s_%s_RP",Order[i_order].Data(),Centrality[i_cent].Data()); // pbar
      p_mFlow_pbar_RP[i_order][i_cent] = new TProfile(ProName.Data(),ProName.Data(),25,0.1,5.1);

      ProName = Form("Flow_Lambda_%s_%s_RP",Order[i_order].Data(),Centrality[i_cent].Data()); // Lambda
      p_mFlow_Lambda_RP[i_order][i_cent] = new TProfile(ProName.Data(),ProName.Data(),25,0.1,5.1);

      ProName = Form("Flow_Lambdabar_%s_%s_RP",Order[i_order].Data(),Centrality[i_cent].Data()); // Lambdabar
      p_mFlow_Lambdabar_RP[i_order][i_cent] = new TProfile(ProName.Data(),ProName.Data(),25,0.1,5.1);

      ProName = Form("Flow_K0s_%s_%s_RP",Order[i_order].Data(),Centrality[i_cent].Data()); // K0s
      p_mFlow_K0s_RP[i_order][i_cent] = new TProfile(ProName.Data(),ProName.Data(),25,0.1,5.1);

      ProName = Form("Flow_phi_%s_%s_RP",Order[i_order].Data(),Centrality[i_cent].Data()); // phi 
      p_mFlow_phi_RP[i_order][i_cent] = new TProfile(ProName.Data(),ProName.Data(),25,0.1,5.1);
    }
  }

  // QA Plot
  h_mPart = new TH1F("h_mPart","h_mPart",2000,0,2000.0);
  h_mMult = new TH1F("h_mMult","h_mMult",10000,0,10000.0);
  h_mRefMult = new TH1F("h_mRefMult","h_mRefMult",10000,0,10000.0);
  h_mEta = new TH1F("h_mEta","h_mEta",1001,-10.01,10.01);
  for(Int_t i_cent = 0; i_cent < 9; i_cent++)
  {
    TString HistName;
    HistName = Form("h_mPsi2_East_%d",i_cent);
    h_mPsi2_East[i_cent] = new TH1F(HistName.Data(),HistName.Data(),100,-TMath::Pi(),TMath::Pi());
    HistName = Form("h_mPsi2_West_%d",i_cent);
    h_mPsi2_West[i_cent] = new TH1F(HistName.Data(),HistName.Data(),100,-TMath::Pi(),TMath::Pi());
    HistName = Form("h_mPsi3_East_%d",i_cent);
    h_mPsi3_East[i_cent] = new TH1F(HistName.Data(),HistName.Data(),100,-TMath::Pi(),TMath::Pi());
    HistName = Form("h_mPsi3_West_%d",i_cent);
    h_mPsi3_West[i_cent] = new TH1F(HistName.Data(),HistName.Data(),100,-TMath::Pi(),TMath::Pi());
  }
  h_mCentrality = new TH1F("h_mCentrality","h_mCentrality",9,-0.5,8.5);

  // initialize the TChain
  if (!mInPutList.IsNull())   // if input file is ok
  {
    TString COUT = Form("Open AMPT %s file list ",mMode_AMPT[mMode].Data());
    cout << COUT.Data() << mInPutList.Data() << endl;
    ifstream in(mInPutList);  // input stream
    if(in)
    {
      cout << "file list is ok" << endl;
      mChain_Input = new TChain("tr");
      char str[255];       // char array for each file name
      Long64_t entries_save = 0;
      while(in)
      {
	in.getline(str,255);  // take the lines of the file list
	if(str[0] != 0)
	{
	  TString addfile;
	  addfile = str;
	  mChain_Input->AddFile(addfile.Data(),-1,"tr");
	  Long64_t file_entries = mChain_Input->GetEntries();
	  cout << "File added to data chain: " << addfile.Data() << " with " << (file_entries-entries_save) << " entries" << endl;
	  entries_save = file_entries;
	}
      }
    }
    else
    {
      cout << "WARNING: SE file input is problemtic" << endl;
      mInput_flag = 0;
    }
  }

  // set input Tree
  if(mInput_flag == 1 && !mChain_Input->GetBranch("Event"))
  {
    cerr << "ERROR: Could not find branch 'tr' in tree!" << endl;
  }

  if(mInput_flag == 1)
  {
    mChain_Input->SetBranchAddress("Event", &Event, &b_Event);
    mChain_Input->SetBranchAddress("Mult", &Mult, &b_Mult);
    mChain_Input->SetBranchAddress("Npartp", &Npartp, &b_Npartp);
    mChain_Input->SetBranchAddress("Npartt", &Npartt, &b_Npartt);
    mChain_Input->SetBranchAddress("Nesp", &Nesp, &b_Nesp);
    mChain_Input->SetBranchAddress("Ninesp", &Ninesp, &b_Ninesp);
    mChain_Input->SetBranchAddress("Nest", &Nest, &b_Nest);
    mChain_Input->SetBranchAddress("Ninest", &Ninest, &b_Ninest);
    mChain_Input->SetBranchAddress("Imp", &Imp, &b_Imp);
    mChain_Input->SetBranchAddress("Na", &Na, &b_Na);
    mChain_Input->SetBranchAddress("Nb", &Nb, &b_Nb);
    mChain_Input->SetBranchAddress("Nab", &Nab, &b_Nab);
    mChain_Input->SetBranchAddress("Psi", &Psi, &b_Psi);
//    mChain_Input->SetBranchAddress("Nx", Nx, &b_Nx);
//    mChain_Input->SetBranchAddress("Ny", Ny, &b_Ny);
//    mChain_Input->SetBranchAddress("Nz", Nz, &b_Nz);
//    mChain_Input->SetBranchAddress("Stat", Stat, &b_Stat);
    mChain_Input->SetBranchAddress("PID", PID, &b_PID);
    mChain_Input->SetBranchAddress("Px", Px, &b_Px);
    mChain_Input->SetBranchAddress("Py", Py, &b_Py);
    mChain_Input->SetBranchAddress("Pz", Pz, &b_Pz);
    mChain_Input->SetBranchAddress("Mass", Mass, &b_Mass);
    mChain_Input->SetBranchAddress("XX", XX, &b_XX);
    mChain_Input->SetBranchAddress("YY", YY, &b_YY);
    mChain_Input->SetBranchAddress("ZZ", ZZ, &b_ZZ);
    mChain_Input->SetBranchAddress("TT", TT, &b_TT);

    Int_t num_events = mChain_Input->GetEntriesFast();
    cout << "Number of events in file(s) = " << num_events << endl;

    if(mStartEvent > num_events) mStartEvent = num_events;
    if(mStopEvent  > num_events) mStopEvent  = num_events;

    cout << "New nStartEvent = " << mStartEvent << ", new nStopEvent = " << mStopEvent << endl;
  }
}
//------------------------------------------------------------
void AMPT_v3v2::Make()
{
  Long64_t start_event_use = mStartEvent;
  Long64_t stop_event_use = mStopEvent;
  mChain_Input->GetEntry(0);

  for(Long64_t i_event = start_event_use; i_event < stop_event_use; i_event++) // event loop
  {
    if (i_event != 0  &&  i_event % 1000 == 0)
      cout << "." << flush;
    if (i_event != 0  &&  i_event % 10000 == 0)
    {
      if((stop_event_use-start_event_use) > 0)
      {
	Double_t event_percent = 100.0*((Double_t)(i_event-start_event_use))/((Double_t)(stop_event_use-start_event_use));
	cout << " " << i_event-start_event_use << " (" << event_percent << "%) " << "\n" << "==> Processing data (AMPT_v3v2) " << flush;
      }
    }

    mChain_Input->GetEntry(i_event);

    TVector3 track;
    Float_t Q2x_east = 0.0, Q2y_east = 0.0, Q2x_west = 0.0, Q2y_west = 0.0;
    Float_t Q3x_east = 0.0, Q3y_east = 0.0, Q3x_west = 0.0, Q3y_west = 0.0;
    Int_t refMult = 0;

    h_mPart->Fill(Npartp+Npartt);
    h_mMult->Fill(Mult);

    for(Int_t i_track = 0; i_track < Mult; i_track++) // 1st track loop for event plane reconstruction
    {
      if(Px[i_track] == 0. && Py[i_track] == 0.) continue;
      track.SetXYZ(Px[i_track],Py[i_track],Pz[i_track]);
      Float_t p_track = track.Mag();
      Float_t pt_track = track.Perp();
      Float_t phi_track = track.Phi(); // -pi to pi
      Float_t eta_track = track.Eta();
      h_mEta->Fill(eta_track);
      // track selection
      if(TMath::Abs(PID[i_track]) == 211 || TMath::Abs(PID[i_track]) == 321 || TMath::Abs(PID[i_track]) == 2212) // pi^{+/-}, K^{+/-}, p and pbar
      {
	if(TMath::Abs(eta_track) < 0.5) refMult++; // refMult calculation
	if(pt_track > 0.2 && pt_track < 2.0 && p_track < 10.0) // pt and p cut
	{
	  if(TMath::Abs(eta_track) <= 1.0) // eta cut
	  {
	    if(eta_track < -0.05) // east
	    {
	      Q2x_east += pt_track*TMath::Cos(2.0*phi_track);
	      Q2y_east += pt_track*TMath::Sin(2.0*phi_track);
	      Q3x_east += pt_track*TMath::Cos(3.0*phi_track);
	      Q3y_east += pt_track*TMath::Sin(3.0*phi_track);
	    }
	    if(eta_track > 0.05) // west
	    {
	      Q2x_west += pt_track*TMath::Cos(2.0*phi_track);
	      Q2y_west += pt_track*TMath::Sin(2.0*phi_track);
	      Q3x_west += pt_track*TMath::Cos(3.0*phi_track);
	      Q3y_west += pt_track*TMath::Sin(3.0*phi_track);
	    }
	  }
	}
      }
    }
    h_mRefMult->Fill(refMult); // fill refMult distribution

    // calculate event plane angle 
    Float_t Psi2_East = TMath::ATan2(Q2y_east,Q2x_east)/2.0; // 2nd: -pi/2 to pi/2
    Float_t Psi2_West = TMath::ATan2(Q2y_west,Q2x_west)/2.0;
    Float_t Psi3_East = TMath::ATan2(Q3y_east,Q3x_east)/3.0; // 3rd: -pi/3 to pi/3
    Float_t Psi3_West = TMath::ATan2(Q3y_west,Q3x_west)/3.0;

    Int_t cent9 = getCentrality(refMult);
    Float_t res2 = getResolution(0,cent9);
    Float_t res3 = getResolution(1,cent9);

//    cout << "refMult = " << refMult << ", Centrality = " << cent9 << endl;

    if(cent9 > -1.0)
    {
      if( // 2nd track loop for v2 calculation
	  !(Q2x_east == 0.0 && Q2y_east == 0.0) 
       && !(Q2x_west == 0.0 && Q2y_west == 0.0)
       && res2 > 0.0
	)
      {
	track.SetXYZ(-999.9,-999.9,-999.9);
	for(Int_t i_track = 0; i_track < Mult; i_track++) // 2nd track loop for flow calculation
	{
	  if(Px[i_track] == 0. && Py[i_track] == 0.) continue;
	  track.SetXYZ(Px[i_track],Py[i_track],Pz[i_track]);
	  // Float_t p_track = track.Mag();
	  Float_t pt_track = track.Perp();
	  Float_t phi_track = track.Phi(); // -pi to pi
	  Float_t eta_track = track.Eta();

	  // track selection
	  if(TMath::Abs(eta_track) < 1.0) // eta cut
	  {
	    if(pt_track < 0.1) continue;
	    Float_t v2;
	    Float_t v2_RP = TMath::Cos(2.0*phi_track);
	    if(eta_track <= 0.0) // east track => west event plane
	    {
//	      v2 = TMath::Cos(2.0*(phi_track-Psi2_West))/res2;
	      v2 = TMath::Cos(2.0*(phi_track-Psi2_West));
	    }
	    if(eta_track > 0.0) // west track => east event plane
	    {
//	      v2 = TMath::Cos(2.0*(phi_track-Psi2_East))/res2;
	      v2 = TMath::Cos(2.0*(phi_track-Psi2_East));
	    }
	    // Centrality bin selection
	    for(Int_t i_cent = AMPT_v3v2::Centrality_start; i_cent < AMPT_v3v2::Centrality_stop; i_cent++)
	    {
	      if(cent9 >= AMPT_v3v2::cent_low[i_cent] && cent9 <= AMPT_v3v2::cent_up[i_cent])
	      {
		// PID
		if(PID[i_track] == 211) // pi_plus
		{
		  p_mFlow_pi_plus[0][i_cent]->Fill(pt_track,v2);
		  p_mFlow_pi_plus_RP[0][i_cent]->Fill(pt_track,v2_RP);
		  h_mPt_pi_plus[0][i_cent]->Fill(pt_track);
		}
		if(PID[i_track] == -211) // pi_minus
		{
		  p_mFlow_pi_minus[0][i_cent]->Fill(pt_track,v2);
		  p_mFlow_pi_minus_RP[0][i_cent]->Fill(pt_track,v2_RP);
		  h_mPt_pi_minus[0][i_cent]->Fill(pt_track);
		}
		if(PID[i_track] == 321) // K_plus
		{
		  p_mFlow_K_plus[0][i_cent]->Fill(pt_track,v2);
		  p_mFlow_K_plus_RP[0][i_cent]->Fill(pt_track,v2_RP);
		  h_mPt_K_plus[0][i_cent]->Fill(pt_track);
		}
		if(PID[i_track] == -321) // K_minus
		{
		  p_mFlow_K_minus[0][i_cent]->Fill(pt_track,v2);
		  p_mFlow_K_minus_RP[0][i_cent]->Fill(pt_track,v2_RP);
		  h_mPt_K_minus[0][i_cent]->Fill(pt_track);
		}
		if(PID[i_track] == 2212) // p
		{
		  p_mFlow_p[0][i_cent]->Fill(pt_track,v2);
		  p_mFlow_p_RP[0][i_cent]->Fill(pt_track,v2_RP);
		  h_mPt_p[0][i_cent]->Fill(pt_track);
		}
		if(PID[i_track] == -2212) // pbar
		{
		  p_mFlow_pbar[0][i_cent]->Fill(pt_track,v2);
		  p_mFlow_pbar_RP[0][i_cent]->Fill(pt_track,v2_RP);
		  h_mPt_pbar[0][i_cent]->Fill(pt_track);
		}
		if(PID[i_track] == 3122) // Lambda
		{
		  p_mFlow_Lambda[0][i_cent]->Fill(pt_track,v2);
		  p_mFlow_Lambda_RP[0][i_cent]->Fill(pt_track,v2_RP);
		  h_mPt_Lambda[0][i_cent]->Fill(pt_track);
		}
		if(PID[i_track] == -3122) // Lambdabar
		{
		  p_mFlow_Lambdabar[0][i_cent]->Fill(pt_track,v2);
		  p_mFlow_Lambdabar_RP[0][i_cent]->Fill(pt_track,v2_RP);
		  h_mPt_Lambdabar[0][i_cent]->Fill(pt_track);
		}
		if(PID[i_track] == 310) // K0s
		{
		  p_mFlow_K0s[0][i_cent]->Fill(pt_track,v2);
		  p_mFlow_K0s_RP[0][i_cent]->Fill(pt_track,v2_RP);
		  h_mPt_K0s[0][i_cent]->Fill(pt_track);
		}
		if(PID[i_track] == 333) // phi
		{
		  p_mFlow_phi[0][i_cent]->Fill(pt_track,v2);
		  p_mFlow_phi_RP[0][i_cent]->Fill(pt_track,v2_RP);
		  h_mPt_phi[0][i_cent]->Fill(pt_track);
		}
	      }
	    }
	  }
	}

	// QA: 2nd event plane and centrality
	h_mPsi2_East[cent9]->Fill(Psi2_East);
	h_mPsi2_West[cent9]->Fill(Psi2_West);
	h_mCentrality->Fill(cent9);
      }
      if( // 3rd track loop for v3 calculation
	  !(Q3x_east == 0.0 && Q3y_east == 0.0) 
       && !(Q3x_west == 0.0 && Q3y_west == 0.0)
       && res3 > 0.0
	)
      {
	track.SetXYZ(-999.9,-999.9,-999.9);
	for(Int_t i_track = 0; i_track < Mult; i_track++) // 2nd track loop for flow calculation
	{
	  if(Px[i_track] == 0. && Py[i_track] == 0.) continue;
	  track.SetXYZ(Px[i_track],Py[i_track],Pz[i_track]);
//	  Float_t p_track = track.Mag();
	  Float_t pt_track = track.Perp();
	  Float_t phi_track = track.Phi(); // -pi to pi
	  Float_t eta_track = track.Eta();

	  // track selection
	  if(TMath::Abs(eta_track) < 1.0) // eta cut
	  {
	    if(pt_track < 0.1) continue;
	    Float_t v3;
	    Float_t v3_RP = TMath::Cos(3.0*phi_track);
	    if(eta_track <= 0.0) // east track => west event plane
	    {
//	      v3 = TMath::Cos(3.0*(phi_track-Psi3_West))/res3;
	      v3 = TMath::Cos(3.0*(phi_track-Psi3_West));
	    }
	    if(eta_track > 0.0) // west track => east event plane
	    {
//	      v3 = TMath::Cos(3.0*(phi_track-Psi3_East))/res3;
	      v3 = TMath::Cos(3.0*(phi_track-Psi3_East));
	    }
	    // Centrality bin selection
	    for(Int_t i_cent = AMPT_v3v2::Centrality_start; i_cent < AMPT_v3v2::Centrality_stop; i_cent++)
	    {
	      if(cent9 >= AMPT_v3v2::cent_low[i_cent] && cent9 <= AMPT_v3v2::cent_up[i_cent])
	      {
		// PID
		if(PID[i_track] == 211) // pi_plus
		{
		  p_mFlow_pi_plus[1][i_cent]->Fill(pt_track,v3);
		  p_mFlow_pi_plus_RP[1][i_cent]->Fill(pt_track,v3_RP);
		  h_mPt_pi_plus[1][i_cent]->Fill(pt_track);
		}
		if(PID[i_track] == -211) // pi_minus
		{
		  p_mFlow_pi_minus[1][i_cent]->Fill(pt_track,v3);
		  p_mFlow_pi_minus_RP[1][i_cent]->Fill(pt_track,v3_RP);
		  h_mPt_pi_minus[1][i_cent]->Fill(pt_track);
		}
		if(PID[i_track] == 321) // K_plus
		{
		  p_mFlow_K_plus[1][i_cent]->Fill(pt_track,v3);
		  p_mFlow_K_plus_RP[1][i_cent]->Fill(pt_track,v3_RP);
		  h_mPt_K_plus[1][i_cent]->Fill(pt_track);
		}
		if(PID[i_track] == -321) // K_minus
		{
		  p_mFlow_K_minus[1][i_cent]->Fill(pt_track,v3);
		  p_mFlow_K_minus_RP[1][i_cent]->Fill(pt_track,v3_RP);
		  h_mPt_K_minus[1][i_cent]->Fill(pt_track);
		}
		if(PID[i_track] == 2212) // p
		{
		  p_mFlow_p[1][i_cent]->Fill(pt_track,v3);
		  p_mFlow_p_RP[1][i_cent]->Fill(pt_track,v3_RP);
		  h_mPt_p[1][i_cent]->Fill(pt_track);
		}
		if(PID[i_track] == -2212) // pbar
		{
		  p_mFlow_pbar[1][i_cent]->Fill(pt_track,v3);
		  p_mFlow_pbar_RP[1][i_cent]->Fill(pt_track,v3_RP);
		  h_mPt_pbar[1][i_cent]->Fill(pt_track);
		}
		if(PID[i_track] == 3122) // Lambda
		{
		  p_mFlow_Lambda[1][i_cent]->Fill(pt_track,v3);
		  p_mFlow_Lambda_RP[1][i_cent]->Fill(pt_track,v3_RP);
		  h_mPt_Lambda[1][i_cent]->Fill(pt_track);
		}
		if(PID[i_track] == -3122) // Lambdabar
		{
		  p_mFlow_Lambdabar[1][i_cent]->Fill(pt_track,v3);
		  p_mFlow_Lambdabar_RP[1][i_cent]->Fill(pt_track,v3_RP);
		  h_mPt_Lambdabar[1][i_cent]->Fill(pt_track);
		}
		if(PID[i_track] == 310) // K0s
		{
		  p_mFlow_K0s[1][i_cent]->Fill(pt_track,v3);
		  p_mFlow_K0s_RP[1][i_cent]->Fill(pt_track,v3_RP);
		  h_mPt_K0s[1][i_cent]->Fill(pt_track);
		}
		if(PID[i_track] == 333) // phi
		{
		  p_mFlow_phi[1][i_cent]->Fill(pt_track,v3);
		  p_mFlow_phi_RP[1][i_cent]->Fill(pt_track,v3_RP);
		  h_mPt_phi[1][i_cent]->Fill(pt_track);
		}
	      }
	    }
	  }
	}

	// QA: 3rd event plane
	h_mPsi3_East[cent9]->Fill(Psi3_East);
	h_mPsi3_West[cent9]->Fill(Psi3_West);
      }
    }
  }

  cout << "." << flush;
  cout << " " << stop_event_use-start_event_use << "(" << 100 << "%)";
  cout << endl;
}

void AMPT_v3v2::Finish()
{
  mFile_Res->Close();
  mFile_OutPut->cd();
  h_mPart->Write();
  h_mMult->Write();
  h_mEta->Write();
  h_mRefMult->Write();
  for(Int_t i_cent = 0; i_cent < 9; i_cent++)
  {
    h_mPsi2_East[i_cent]->Write();
    h_mPsi2_West[i_cent]->Write();
    h_mPsi3_East[i_cent]->Write();
    h_mPsi3_West[i_cent]->Write();
  }
  h_mCentrality->Write();
  for(Int_t i_order = 0; i_order < 2; i_order++)
  {
    for(Int_t i_cent = 0; i_cent < 4; i_cent++)
    {
      // v2 relative to Event Plane
      p_mFlow_pi_plus[i_order][i_cent]->Write();
      p_mFlow_pi_minus[i_order][i_cent]->Write();
      p_mFlow_K_plus[i_order][i_cent]->Write();
      p_mFlow_K_minus[i_order][i_cent]->Write();
      p_mFlow_p[i_order][i_cent]->Write();
      p_mFlow_pbar[i_order][i_cent]->Write();
      p_mFlow_Lambda[i_order][i_cent]->Write();
      p_mFlow_Lambdabar[i_order][i_cent]->Write();
      p_mFlow_K0s[i_order][i_cent]->Write();
      p_mFlow_phi[i_order][i_cent]->Write();

      // v2 relative to Reaction Plane
      p_mFlow_pi_plus_RP[i_order][i_cent]->Write();
      p_mFlow_pi_minus_RP[i_order][i_cent]->Write();
      p_mFlow_K_plus_RP[i_order][i_cent]->Write();
      p_mFlow_K_minus_RP[i_order][i_cent]->Write();
      p_mFlow_p_RP[i_order][i_cent]->Write();
      p_mFlow_pbar_RP[i_order][i_cent]->Write();
      p_mFlow_Lambda_RP[i_order][i_cent]->Write();
      p_mFlow_Lambdabar_RP[i_order][i_cent]->Write();
      p_mFlow_K0s_RP[i_order][i_cent]->Write();
      p_mFlow_phi_RP[i_order][i_cent]->Write();

      h_mPt_pi_plus[i_order][i_cent]->Write();
      h_mPt_pi_minus[i_order][i_cent]->Write();
      h_mPt_K_plus[i_order][i_cent]->Write();
      h_mPt_K_minus[i_order][i_cent]->Write();
      h_mPt_p[i_order][i_cent]->Write();
      h_mPt_pbar[i_order][i_cent]->Write();
      h_mPt_Lambda[i_order][i_cent]->Write();
      h_mPt_Lambdabar[i_order][i_cent]->Write();
      h_mPt_K0s[i_order][i_cent]->Write();
      h_mPt_phi[i_order][i_cent]->Write();
    }
  }
  mFile_OutPut->Close();
}

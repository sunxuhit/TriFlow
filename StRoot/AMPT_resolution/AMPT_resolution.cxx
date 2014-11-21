#include "AMPT_resolution.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"
#include "TProfile.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"

ClassImp(AMPT_resolution)

Int_t AMPT_resolution::mInput_flag = 1;
TString AMPT_resolution::mBeamEnergy[7] = {"7GeV","11GeV","19GeV","27GeV","39GeV","62GeV","200GeV"};
TString AMPT_resolution::mMode_AMPT[2] = {"Default","StringMelting"};
Int_t AMPT_resolution::mCentrality[2][7][10] = {
						 { // Default
						   {9,16,27,42,63,92,130,182,219,329}, //  7GeV
						   {10,18,31,49,74,106,148,207,246,352}, // 11GeV
						   {12,21,35,56,83,119,166,229,270,385}, // 19GeV
						   {13,23,39,60,90,128,178,246,289,422}, // 27GeV
						   {15,26,43,68,100,143,198,273,321,478}, // 39GeV | with phi decay
//						   {14,26,43,67,99,142,197,271,318,462}, // 39GeV | with phi partly decay
						   {16,28,48,75,112,161,225,312,368,521}, // 62GeV
						   {23,42,72,116,176,257,364,513,610,890}  // 200GeV
						 },
						 { // String Melting
						   {9,16,27,41,61,86,119,163,193,291}, //  7GeV
						   {11,19,32,51,75,106,147,201,237,345}, // 11GeV
						   {12,22,38,59,88,126,175,240,281,414}, // 19GeV
						   {13,24,40,64,95,137,190,261,307,443}, // 27GeV
						   {15,26,43,68,100,143,198,273,321,491}, // 39GeV | with phi decay
//						   {15,26,45,70,105,151,210,289,339,503}, // 39GeV | with phi partly decay
						   {16,29,49,77,117,169,238,331,390,628},  // 62GeV
						   {23,42,73,119,184,273,392,557,666,1020}  // 200GeV
						 }
					       }; // 80%,70%,60%,50%,40%,30%,20%,10%,5%,0%
Int_t AMPT_resolution::mList_start[20] = {  1,101,201,301,401,501,601,701,801, 901,1001,1101,1201,1301,1401,1501,1601,1701,1801,1901};
Int_t AMPT_resolution::mList_stop[20]  = {100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000};
//------------------------------------------------------------
AMPT_resolution::AMPT_resolution(Int_t Energy, Int_t Mode, Int_t List, Long64_t StartEvent, Long64_t StopEvent)
{
  mEnergy = Energy;
  mMode = Mode;
  TString InPutList = Form("/project/projectdirs/star/xusun/OutPut/AMPT_%s/List/%s_List/Split_%s_%d_%d.list",mMode_AMPT[Mode].Data(),mBeamEnergy[mEnergy].Data(),mBeamEnergy[mEnergy].Data(),mList_start[List],mList_stop[List]);
  SetInPutList(InPutList);
  SetStartEvent(StartEvent);
  SetStopEvent(StopEvent);
  TString OutPutFile = Form("/project/projectdirs/star/xusun/OutPut/AMPT_%s/Resolution/%s_Resolution/Resolution_%s_%d_%d.root",mMode_AMPT[Mode].Data(),mBeamEnergy[mEnergy].Data(),mBeamEnergy[mEnergy].Data(),mList_start[List],mList_stop[List]);
  SetOutPutFile(OutPutFile);
}

AMPT_resolution::~AMPT_resolution()
{
}
//------------------------------------------------------------
void AMPT_resolution::SetInPutList(const TString inputlist)
{
  mInPutList= inputlist.Copy();
  cout << "Input list was set to: " << mInPutList.Data() << endl;
}

void AMPT_resolution::SetOutPutFile(const TString outputfile)
{
  mOutPutFile = outputfile.Copy();
  cout << "Output file was set to: " << mOutPutFile.Data() << endl;
}

void AMPT_resolution::SetStartEvent(const Long64_t StartEvent)
{
  mStartEvent = StartEvent;
  cout << "nStartEvent = " << mStartEvent << endl;
}

void AMPT_resolution::SetStopEvent(const Long64_t StopEvent)
{
  mStopEvent = StopEvent;
  cout << "nStopEvent = " << mStopEvent << endl;
}
//------------------------------------------------------------
void AMPT_resolution::Init()
{
  mFile_OutPut = new TFile(mOutPutFile.Data(),"RECREATE");

  p_mRes2 = new TProfile("p_mRes2","p_mRes2",9,-0.5,8.5);
  p_mRes3 = new TProfile("p_mRes3","p_mRes3",9,-0.5,8.5);

  // QA Plot
  h_mPart = new TH1F("h_mPart","h_mPart",2000,0,2000.0);
  h_mMult = new TH1F("h_mMult","h_mMult",10000,0,10000.0);
  h_mRefMult = new TH1F("h_mRefMult","h_mRefMult",10000,0,10000.0);
  h_mEta = new TH1F("h_mEta","h_mEta",800,-10.0,10.0);
  h_mPsi2_East = new TH1F("h_mPsi2_East","h_mPsi2_East",100,-TMath::Pi()/2.0,TMath::Pi()/2.0);
  h_mPsi2_West = new TH1F("h_mPsi2_West","h_mPsi2_West",100,-TMath::Pi()/2.0,TMath::Pi()/2.0);
  h_mPsi3_East = new TH1F("h_mPsi3_East","h_mPsi3_East",100,-TMath::Pi()/3.0,TMath::Pi()/3.0);
  h_mPsi3_West = new TH1F("h_mPsi3_West","h_mPsi3_West",100,-TMath::Pi()/3.0,TMath::Pi()/3.0);
  h_mPsi2 = new TH2F("h_mPsi2","h_mPsi2",100,-TMath::Pi()/2.0,TMath::Pi()/2.0,100,-TMath::Pi()/2.0,TMath::Pi()/2.0);
  h_mPsi3 = new TH2F("h_mPsi3","h_mPsi3",100,-TMath::Pi()/3.0,TMath::Pi()/3.0,100,-TMath::Pi()/3.0,TMath::Pi()/3.0);

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
    mChain_Input->SetBranchAddress("Nx", Nx, &b_Nx);
    mChain_Input->SetBranchAddress("Ny", Ny, &b_Ny);
    mChain_Input->SetBranchAddress("Nz", Nz, &b_Nz);
    mChain_Input->SetBranchAddress("Stat", Stat, &b_Stat);
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
void AMPT_resolution::Make()
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
	cout << " " << i_event-start_event_use << " (" << event_percent << "%) " << "\n" << "==> Processing data (AMPT_resolution) " << flush;
      }
    }

    // start to calculate Resolution
    mChain_Input->GetEntry(i_event);

    TVector3 track;
    Float_t Q2x_full = 0.0, Q2y_full = 0.0, Q2x_east = 0.0, Q2y_east = 0.0, Q2x_west = 0.0, Q2y_west = 0.0;
    Float_t Q3x_full = 0.0, Q3y_full = 0.0, Q3x_east = 0.0, Q3y_east = 0.0, Q3x_west = 0.0, Q3y_west = 0.0;
    Int_t Num_esat = 0, Num_west = 0;
    Int_t refMult = 0;

    h_mPart->Fill(Npartp+Npartt);
    h_mMult->Fill(Mult);

    for(Int_t i_track = 0; i_track < Mult; i_track++) // track loop
    {
      if(Px[i_track] == 0. && Py[i_track] == 0.) continue;
      track.SetXYZ(Px[i_track],Py[i_track],Pz[i_track]);
      Float_t p_track = track.Mag();
      Float_t pt_track = track.Perp();
      Float_t phi_track = track.Phi(); // -pi to pi
      Float_t eta_track = track.Eta();
      h_mEta->Fill(eta_track);
//      cout << "px = " << Px[i_track] << ", py = " << Py[i_track] << ", pz = " << Pz[i_track] << endl;
//      cout << "p_track = " << p_track << ", pt_track = " << pt_track << endl;
      // track selection
      if(TMath::Abs(PID[i_track]) == 211 || TMath::Abs(PID[i_track]) == 321 || TMath::Abs(PID[i_track]) == 2212) // pi^{+/-}, K^{+/-}, p and pbar
      {
	if(TMath::Abs(eta_track) < 0.5) refMult++; // refMult calculation
	if(pt_track > 0.2 && pt_track < 2.0 && p_track < 10.0) // pt and p cut
//	if(pt_track >= 0.2 && pt_track < 2.0) // pt and p cut
	{
	  if(TMath::Abs(eta_track) <= 1.0) // eta cut
	  {
	    // full event
	    Q2x_full += TMath::Cos(2.0*phi_track);
	    Q2y_full += TMath::Sin(2.0*phi_track);
	    Q3x_full += TMath::Cos(3.0*phi_track);
	    Q3y_full += TMath::Sin(3.0*phi_track);

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
    // Centrality defination
    Int_t cent9 = -1;;
    for(Int_t i_cent = 0; i_cent < 10; i_cent++)
    {
      if(refMult >= mCentrality[mMode][mEnergy][i_cent] && refMult < mCentrality[mMode][mEnergy][i_cent+1])
      {
	cent9 = i_cent;
      }
    }

    if(cent9 > -1.0)
    {
      // eta_sub method event plane resolution calculation
      if( 
	  !(Q2x_east == 0.0 && Q2y_east == 0.0) 
       && !(Q2x_west == 0.0 && Q2y_west == 0.0)
	)
      {
	Float_t Psi2_east = TMath::ATan2(Q2y_east,Q2x_east)/2.0; // -pi/2 to pi/2
	Float_t Psi2_west = TMath::ATan2(Q2y_west,Q2x_west)/2.0;
	p_mRes2->Fill(cent9,TMath::Cos(2.0*(Psi2_east-Psi2_west)));
	h_mPsi2_East->Fill(Psi2_east);
	h_mPsi2_West->Fill(Psi2_west);
	h_mPsi2->Fill(Psi2_east,Psi2_west);
      }
      if( 
	  !(Q3x_east == 0.0 && Q3y_east == 0.0) 
       && !(Q3x_west == 0.0 && Q3y_west == 0.0)
	)
      {
	Float_t Psi3_east = TMath::ATan2(Q3y_east,Q3x_east)/3.0; // -pi/3 to pi/3
	Float_t Psi3_west = TMath::ATan2(Q3y_west,Q3x_west)/3.0;
	p_mRes3->Fill(cent9,TMath::Cos(3.0*(Psi3_east-Psi3_west)));
	h_mPsi3_East->Fill(Psi3_east);
	h_mPsi3_West->Fill(Psi3_west);
	h_mPsi3->Fill(Psi3_east,Psi3_west);
      }
    }
  }

  cout << "." << flush;
  cout << " " << stop_event_use-start_event_use << "(" << 100 << "%)";
  cout << endl;
}

void AMPT_resolution::Finish()
{
  mFile_OutPut->cd();
  p_mRes2->Write();
  p_mRes3->Write();
  h_mPart->Write();
  h_mMult->Write();
  h_mEta->Write();
  h_mPsi2_East->Write();
  h_mPsi2_West->Write();
  h_mPsi3_East->Write();
  h_mPsi3_West->Write();
  h_mPsi2->Write();
  h_mPsi3->Write();
  h_mRefMult->Write();
  mFile_OutPut->Close();
}

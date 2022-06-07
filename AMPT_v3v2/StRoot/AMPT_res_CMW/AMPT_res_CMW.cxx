#include "AMPT_res_CMW.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"
#include "TProfile.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"

ClassImp(AMPT_res_CMW)

Int_t AMPT_res_CMW::mInput_flag = 1;
TString AMPT_res_CMW::mBeamEnergy[7] = {"7GeV","11GeV","19GeV","27GeV","39GeV","62GeV","200GeV"};
TString AMPT_res_CMW::mMode_AMPT[2] = {"Default","StringMelting"};
TString AMPT_res_CMW::mScreenMass_AMPT[3] = {"1mb","3mb","6mb"};
Int_t AMPT_res_CMW::mCentrality[2][7][10] = {
						 { // Default
						   {9,16,27,42,63,92,130,182,219,329}, //  7GeV
						   {9,17,30,48,72,104,147,205,243,378}, // 11GeV
						   {12,21,35,56,83,119,166,229,270,385}, // 19GeV
						   {12,22,37,59,88,126,176,242,286,432}, // 27GeV
//						   {10,19,32,52,77,111,155,214,252,377}, // 39GeV 1.5mb
						   {13,25,42,66,98,140,195,269,316,474}, // 39GeV | with phi decay
						   {15,27,46,73,110,158,222,307,362,521}, // 62GeV
//						   {15,29,51,84,129,190,273,387,462,693}  // 200GeV new parameters
						   {22,41,71,114,173,253,360,506,602,899}  // 200GeV
						 },
						 { // String Melting
						   {9,16,27,41,61,86,119,163,193,291}, //  7GeV
//						   {8,15,27,43,66,95,134,186,220,338}, // 11GeV | 1.5mb
//						   {10,18,31,49,73,105,145,200,235,365}, // 11GeV | 3mb
						   {10,18,31,48,71,100,139,190,223,340}, // 11GeV | 6mb
						   {12,22,38,59,88,126,175,240,281,414}, // 19GeV
						   {9,18,31,50,76,111,156,217,256,385}, // 27GeV | 1.5 mb
//						   {12,23,39,62,94,134,187,258,303,442}, // 27GeV | 3 mb
//						   {11,20,34,55,84,122,172,240,283,433}, // 39GeV | with phi decay | 1.5 mb
//						   {14,25,43,69,104,149,208,287,337,504}, // 39GeV | with phi decay | 3 mb
						   {14,25,42,67,100,143,200,276,324,491}, // 39GeV | with phi decay | 6 mb
						   {11,21,36,59,92,135,193,271,321,477},  // 62GeV | 1.5 mb
//						   {15,27,47,75,114,166,234,326,384,559},  // 62GeV | 3 mb
//						   {15,29,51,86,137,206,300,431,517,802}  // 200GeV | with phi decay | 1.5 mb
//						   {21,40,71,117,181,268,386,549,656,1042}  // 200GeV | with phi decay | 3 mb
//						   {21,42,72,117,181,269,388,551,659,978}  // 200GeV | with phi decay | 3 mb v1.21/v2.21
						   {21,40,70,113,175,260,375,535,640,1087}  // 200GeV | with phi decay | 6 mb
						 }
					       }; // 80%,70%,60%,50%,40%,30%,20%,10%,5%,0%
Int_t AMPT_res_CMW::mList_start[25] = {  1,101,201,301,401,501,601,701,801, 901,1001,1101,1201,1301,1401,1501,1601,1701,1801,1901,2001,2101,2201,2301,2401};
Int_t AMPT_res_CMW::mList_stop[25]  = {100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2100,2200,2300,2400,2500};
//------------------------------------------------------------
AMPT_res_CMW::AMPT_res_CMW(Int_t Energy, Int_t Mode, Int_t Screen, Int_t List, Long64_t StartEvent, Long64_t StopEvent)
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
  SetInPutList(InPutList);
  SetStartEvent(StartEvent);
  SetStopEvent(StopEvent);

  TString OutPutFile;
  if(mMode == 0)
  {
    OutPutFile = Form("/project/projectdirs/star/xusun/OutPut/AMPT_%s/Res_CMW/%s_Resolution/Resolution_%s_%d_%d.root",mMode_AMPT[Mode].Data(),mBeamEnergy[mEnergy].Data(),mBeamEnergy[mEnergy].Data(),mList_start[List],mList_stop[List]);
  }
  if(mMode == 1)
  {
    OutPutFile = Form("/project/projectdirs/star/xusun/OutPut/AMPT_%s/Res_CMW/%s_Resolution/%s/Resolution_%s_%d_%d.root",mMode_AMPT[Mode].Data(),mBeamEnergy[mEnergy].Data(),mScreenMass_AMPT[mScreen].Data(),mBeamEnergy[mEnergy].Data(),mList_start[List],mList_stop[List]);
  }
  SetOutPutFile(OutPutFile);
}

AMPT_res_CMW::~AMPT_res_CMW()
{
}
//------------------------------------------------------------
void AMPT_res_CMW::SetInPutList(const TString inputlist)
{
  mInPutList= inputlist.Copy();
  cout << "Input list was set to: " << mInPutList.Data() << endl;
}

void AMPT_res_CMW::SetOutPutFile(const TString outputfile)
{
  mOutPutFile = outputfile.Copy();
  cout << "Output file was set to: " << mOutPutFile.Data() << endl;
}

void AMPT_res_CMW::SetStartEvent(const Long64_t StartEvent)
{
  mStartEvent = StartEvent;
  cout << "nStartEvent = " << mStartEvent << endl;
}

void AMPT_res_CMW::SetStopEvent(const Long64_t StopEvent)
{
  mStopEvent = StopEvent;
  cout << "nStopEvent = " << mStopEvent << endl;
}
//------------------------------------------------------------
void AMPT_res_CMW::Init()
{
  mFile_OutPut = new TFile(mOutPutFile.Data(),"RECREATE");

  p_mRes2 = new TProfile("p_mRes2","p_mRes2",9,-0.5,8.5);
  p_mRes3 = new TProfile("p_mRes3","p_mRes3",9,-0.5,8.5);

  // QA Plot
  h_mPart = new TH1D("h_mPart","h_mPart",2000,0,2000.0);
  h_mMult = new TH1D("h_mMult","h_mMult",10000,0,10000.0);
  h_mRefMult = new TH1D("h_mRefMult","h_mRefMult",10000,-0.5,9999.5);
  h_mEta = new TH1D("h_mEta","h_mEta",1001,-10.01,10.01);
  for(Int_t i_cent = 0; i_cent < 9; i_cent++)
  {
    TString HistName;
    HistName = Form("h_mPsi2_East_%d",i_cent);
    h_mPsi2_East[i_cent] = new TH1D(HistName.Data(),HistName.Data(),100,-TMath::Pi(),TMath::Pi());
    HistName = Form("h_mPsi2_West_%d",i_cent);
    h_mPsi2_West[i_cent] = new TH1D(HistName.Data(),HistName.Data(),100,-TMath::Pi(),TMath::Pi());
    HistName = Form("h_mPsi3_East_%d",i_cent);
    h_mPsi3_East[i_cent] = new TH1D(HistName.Data(),HistName.Data(),100,-TMath::Pi(),TMath::Pi());
    HistName = Form("h_mPsi3_West_%d",i_cent);
    h_mPsi3_West[i_cent] = new TH1D(HistName.Data(),HistName.Data(),100,-TMath::Pi(),TMath::Pi());
    HistName = Form("h_mPsi2_%d",i_cent);
    h_mPsi2[i_cent] = new TH2D(HistName.Data(),HistName.Data(),100,-TMath::Pi(),TMath::Pi(),100,-TMath::Pi(),TMath::Pi());
    HistName = Form("h_mPsi3_%d",i_cent);
    h_mPsi3[i_cent] = new TH2D(HistName.Data(),HistName.Data(),100,-TMath::Pi(),TMath::Pi(),100,-TMath::Pi(),TMath::Pi());
  }

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
void AMPT_res_CMW::Make()
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
	cout << " " << i_event-start_event_use << " (" << event_percent << "%) " << "\n" << "==> Processing data (AMPT_res_CMW) " << flush;
      }
    }

    // start to calculate Resolution
    mChain_Input->GetEntry(i_event);

    TVector3 track;
    Float_t Q2x_east = 0.0, Q2y_east = 0.0, Q2x_west = 0.0, Q2y_west = 0.0;
    Float_t Q3x_east = 0.0, Q3y_east = 0.0, Q3x_west = 0.0, Q3y_west = 0.0;
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
    for(Int_t i_cent = 0; i_cent < 9; i_cent++)
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
	h_mPsi2_East[cent9]->Fill(Psi2_east);
	h_mPsi2_West[cent9]->Fill(Psi2_west);
	h_mPsi2[cent9]->Fill(Psi2_east,Psi2_west);
      }
      if( 
	  !(Q3x_east == 0.0 && Q3y_east == 0.0) 
       && !(Q3x_west == 0.0 && Q3y_west == 0.0)
	)
      {
	Float_t Psi3_east = TMath::ATan2(Q3y_east,Q3x_east)/3.0; // -pi/3 to pi/3
	Float_t Psi3_west = TMath::ATan2(Q3y_west,Q3x_west)/3.0;
	p_mRes3->Fill(cent9,TMath::Cos(3.0*(Psi3_east-Psi3_west)));
	h_mPsi3_East[cent9]->Fill(Psi3_east);
	h_mPsi3_West[cent9]->Fill(Psi3_west);
	h_mPsi3[cent9]->Fill(Psi3_east,Psi3_west);
      }
    }
  }

  cout << "." << flush;
  cout << " " << stop_event_use-start_event_use << "(" << 100 << "%)";
  cout << endl;
}

void AMPT_res_CMW::Finish()
{
  mFile_OutPut->cd();
  p_mRes2->Write();
  p_mRes3->Write();
  h_mPart->Write();
  h_mMult->Write();
  h_mEta->Write();
  for(Int_t i_cent = 0; i_cent < 9; i_cent++)
  {
    h_mPsi2_East[i_cent]->Write();
    h_mPsi2_West[i_cent]->Write();
    h_mPsi3_East[i_cent]->Write();
    h_mPsi3_West[i_cent]->Write();
    h_mPsi2[i_cent]->Write();
    h_mPsi3[i_cent]->Write();
  }
  h_mRefMult->Write();
  mFile_OutPut->Close();
}

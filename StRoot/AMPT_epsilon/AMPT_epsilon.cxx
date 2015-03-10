#include "AMPT_epsilon.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"
#include "TProfile.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TMath.h"

ClassImp(AMPT_epsilon)

Int_t AMPT_epsilon::mInput_flag = 1;
TString AMPT_epsilon::mBeamEnergy[7] = {"7GeV","11GeV","19GeV","27GeV","39GeV","62GeV","200GeV"};
TString AMPT_epsilon::mMode_AMPT[2] = {"Default","StringMelting"};
TString AMPT_epsilon::mScreenMass_AMPT[3] = {"1mb","3mb","6mb"};
Int_t AMPT_epsilon::mRefMult[2][7][10] = {
	    				   { // Default
					     {9,16,27,42,63,92,130,182,219,329}, //  7GeV
					     {9,17,30,48,72,104,147,205,243,378}, // 11GeV
					     {12,21,35,56,83,119,166,229,270,385}, // 19GeV
					     {13,23,39,60,90,128,178,246,289,422}, // 27GeV
//					     {10,19,32,52,77,111,155,214,252,377}, // 39GeV 1.5mb
					     {13,25,42,66,98,140,195,269,316,474}, // 39GeV | with phi decay
					     {16,28,48,75,112,161,225,312,368,521}, // 62GeV
//					     {15,29,51,84,129,190,273,387,462,693}  // 200GeV new parameters
					     {22,41,71,114,173,253,360,506,602,899}  // 200GeV
					   },
					   { // String Melting
					     {9,16,27,41,61,86,119,163,193,291}, //  7GeV
//					     {8,15,27,43,66,95,134,186,220,338}, // 11GeV | 1.5mb
//					     {10,18,31,49,73,105,145,200,235,365}, // 11GeV | 3mb
					     {10,18,31,48,71,100,139,190,223,340}, // 11GeV | 6mb
					     {12,22,38,59,88,126,175,240,281,414}, // 19GeV
					     {13,24,40,64,95,137,190,261,307,443}, // 27GeV
//					     {11,20,34,55,84,122,172,240,283,433}, // 39GeV | with phi decay | 1.5 mb
//					     {14,25,43,69,104,149,208,287,337,504}, // 39GeV | with phi decay | 3 mb
					     {14,25,42,67,100,143,200,276,324,491}, // 39GeV | with phi decay | 6 mb
					     {16,29,49,77,117,169,238,331,390,628},  // 62GeV
//					     {15,29,51,86,137,206,300,431,517,802}  // 200GeV | with phi decay | 1.5 mb
//					     {21,40,71,117,181,268,386,549,656,1042}  // 200GeV | with phi decay | 3 mb
//					     {21,42,72,117,181,269,388,551,659,978}  // 200GeV | with phi decay | 3 mb v1.21/v2.21
					     {21,40,70,113,175,260,375,535,640,1087}  // 200GeV | with phi decay | 6 mb
					   }
				         }; // 80%,70%,60%,50%,40%,30%,20%,10%,5%,0%
// Centrality bin
Int_t AMPT_epsilon::cent_low[4] = {0,7,4,0}; // 0 = 0-80%, 1 = 0-10%, 2 = 10-40%, 3 = 40-80%
Int_t AMPT_epsilon::cent_up[4]  = {8,8,6,3}; // 0 = 0-80%, 1 = 0-10%, 2 = 10-40%, 3 = 40-80%
Int_t AMPT_epsilon::mList_start[25] = {  1,101,201,301,401,501,601,701,801, 901,1001,1101,1201,1301,1401,1501,1601,1701,1801,1901,2001,2101,2201,2301,2401};
Int_t AMPT_epsilon::mList_stop[25]  = {100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2100,2200,2300,2400,2500};
Int_t AMPT_epsilon::Centrality_start = 0;
Int_t AMPT_epsilon::Centrality_stop  = 4;
//------------------------------------------------------------
AMPT_epsilon::AMPT_epsilon(Int_t Energy, Int_t Mode, Int_t Screen, Int_t List, Long64_t StartEvent, Long64_t StopEvent)
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
    OutPutFile = Form("/project/projectdirs/star/xusun/OutPut/AMPT_%s/Epsilon/%s_%s/Epsilon_%s_%d_%d.root",mMode_AMPT[Mode].Data(),mBeamEnergy[mEnergy].Data(),mMode_AMPT[Mode].Data(),mBeamEnergy[mEnergy].Data(),mList_start[List],mList_stop[List]);
  }
  if(mMode == 1)
  {
    OutPutFile = Form("/project/projectdirs/star/xusun/OutPut/AMPT_%s/Epsilon/%s_%s/%s/Epsilon_%s_%d_%d.root",mMode_AMPT[Mode].Data(),mBeamEnergy[mEnergy].Data(),mMode_AMPT[Mode].Data(),mScreenMass_AMPT[mScreen].Data(),mBeamEnergy[mEnergy].Data(),mList_start[List],mList_stop[List]);
  }
  SetOutPutFile(OutPutFile); // set output file
}

AMPT_epsilon::~AMPT_epsilon()
{
}
//------------------------------------------------------------
void AMPT_epsilon::SetInPutList(const TString inputlist)
{
  mInPutList= inputlist.Copy();
  cout << "Input list was set to: " << mInPutList.Data() << endl;
}

void AMPT_epsilon::SetOutPutFile(const TString outputfile)
{
  mOutPutFile = outputfile.Copy();
  cout << "Output file was set to: " << mOutPutFile.Data() << endl;
}

void AMPT_epsilon::SetInPutRes(const TString inputres)
{
  mInPutRes = inputres.Copy();
  cout << "Input resolution file was set to: " << mInPutRes.Data() << endl;
}

void AMPT_epsilon::SetStartEvent(const Long64_t StartEvent)
{
  mStartEvent = StartEvent;
  cout << "nStartEvent = " << mStartEvent << endl;
}

void AMPT_epsilon::SetStopEvent(const Long64_t StopEvent)
{
  mStopEvent = StopEvent;
  cout << "nStopEvent = " << mStopEvent << endl;
}
//------------------------------------------------------------
Float_t AMPT_epsilon::getResolution(Int_t order, Int_t i_cent) // 0 for 2nd, 1 for 3rd
{
  Float_t res = -999.9;
  if(p_mRes[order]->GetBinContent(p_mRes[order]->FindBin(i_cent)) > 0.0)
  {
    res = TMath::Sqrt(p_mRes[order]->GetBinContent(p_mRes[order]->FindBin(i_cent)));
  }
  return res;
}

Int_t AMPT_epsilon::getCentrality(Int_t refMult)
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
void AMPT_epsilon::Init()
{
  // initialize the resolution
  mFile_Res = TFile::Open(mInPutRes.Data());
  p_mRes[0] = (TProfile*)mFile_Res->Get("p_mRes2"); // 2nd event plane resolution
  p_mRes[1] = (TProfile*)mFile_Res->Get("p_mRes3"); // 3rd event plane resolution

  mFile_OutPut = new TFile(mOutPutFile.Data(),"RECREATE");

  // initialize Epsilon9 and Epsilon4
  p_mEpsilon9[0] = new TProfile("p_mEpsilon9_2nd","p_mEpsilon9_2nd",9,-0.5,8.5);
  p_mEpsilon9[1] = new TProfile("p_mEpsilon9_3rd","p_mEpsilon9_3rd",9,-0.5,8.5);

  p_mEpsilon4[0] = new TProfile("p_mEpsilon4_2nd","p_mEpsilon4_2nd",4,-0.5,3.5);
  p_mEpsilon4[1] = new TProfile("p_mEpsilon4_3rd","p_mEpsilon4_3rd",4,-0.5,3.5);

  // QA Plot

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
void AMPT_epsilon::Make()
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
	cout << " " << i_event-start_event_use << " (" << event_percent << "%) " << "\n" << "==> Processing data (AMPT_epsilon) " << flush;
      }
    }

    mChain_Input->GetEntry(i_event);

    Int_t refMult = 0;
    TVector3 track;
    for(Int_t i_track = 0; i_track < Mult; i_track++) // refMult calculation
    {
      if(Px[i_track] == 0. && Py[i_track] == 0.) continue;
      track.SetXYZ(Px[i_track],Py[i_track],Pz[i_track]);
      Float_t eta_track = track.Eta();
      // track selection
      if(TMath::Abs(PID[i_track]) == 211 || TMath::Abs(PID[i_track]) == 321 || TMath::Abs(PID[i_track]) == 2212) // pi^{+/-}, K^{+/-}, p and pbar
      {
	if(TMath::Abs(eta_track) < 0.5) refMult++; // refMult calculation
      }
    }

    Int_t cent9 = getCentrality(refMult);
    Float_t res[2]; // 0 for 2nd, 1 for 3rd
    res[0] = getResolution(0,cent9);
    res[1] = getResolution(1,cent9);

    TVector3 particle;
    Float_t epsilon[2]; // 0: 2nd epsilon, 1: 3rd epsilon
    Float_t Order[2] = {2.0,3.0};
    Float_t mean_R2X[2] = {0.0,0.0};
    Float_t mean_R2Y[2] = {0.0,0.0};
    Float_t mean_R2 = 0.0;
    Int_t nPart = Npartp+Npartt;

    for(Int_t i_par = 0; i_par < Nab; i_par++)
    {
      if(Stat[i_par] > 0) // participant
      {
	particle.SetXYZ(Nx[i_par],Ny[i_par],Nz[i_par]);
	Float_t R = particle.Perp(); // return the transverse component (R in cylindrical coordinate system) 
	Float_t phi = particle.Phi(); // return the  azimuth angle. returns phi from -pi to pi

	mean_R2 += R*R/nPart;
	for(Int_t i_order = 0; i_order < 2; i_order++)
	{
	  mean_R2X[i_order] += R*R*TMath::Cos(Order[i_order]*phi)/nPart;
	  mean_R2Y[i_order] += R*R*TMath::Sin(Order[i_order]*phi)/nPart;
	}
      }
    }
    for(Int_t i_order = 0; i_order < 2; i_order++)
    {
      epsilon[i_order] = TMath::Sqrt(mean_R2X[i_order]*mean_R2X[i_order]+mean_R2Y[i_order]*mean_R2Y[i_order])/mean_R2;
    }
    if(cent9 > -1)
    {
      for(Int_t i_order = 0; i_order < 2; i_order++)
      {
	if(res[i_order] > 0.0) 
	{
	  p_mEpsilon9[i_order]->Fill(cent9,epsilon[i_order],refMult); // calculate epsilon for narrow centrality bin

	  for(Int_t i_cent = AMPT_epsilon::Centrality_start; i_cent < AMPT_epsilon::Centrality_stop; i_cent++) // calculate epsilon for wide centrality bin
	  {
	    if(cent9 >= AMPT_epsilon::cent_low[i_cent] && cent9 <= AMPT_epsilon::cent_up[i_cent])
	    {
	      p_mEpsilon4[i_order]->Fill(i_cent,epsilon[i_order],refMult); // calculate epsilon for narrow centrality bin
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

void AMPT_epsilon::Finish()
{
  mFile_Res->Close();

  mFile_OutPut->cd();
  for(Int_t i_order = 0; i_order < 2; i_order++)
  {
    p_mEpsilon9[i_order]->Write();
    p_mEpsilon4[i_order]->Write();
  }
  mFile_OutPut->Close();
}

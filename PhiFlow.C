#include <TSystem>

void PhiFlow(const Int_t energy = 1, const Int_t X_flag = 0, const Int_t List = 0, const Long64_t start_event = 1000, const Long64_t stop_event = 100024, const Int_t mode = 3)
{
  // energy: 0 for 200GeV, 1 for 39GeV
  // X_flag: 0 for Same Event, 1 for Mixed Event
  // List: different number for different TTree list
  // mode: 0 for phi meson, 1 for Lambda, 2 for anti-Lambda, 3 for K0s
//  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
//  loadSharedLibraries();
  // root4star -b -q PhiFlow.C\(1,0,0,1000,1024,1\)

  gSystem->Load("StRefMultCorr");
  gSystem->Load("StAlexPhiMesonEvent");
  gSystem->Load("StV0Event");
  gSystem->Load("StStrangenessAna");
  gSystem->Load("StRunIdEventsDb");

  cout << "All libraries are loaded!!!!" << endl;

  //**************************** Set graphic style ***************************************
  gStyle->SetPalette(1);
  gStyle->SetCanvasColor(10);
  gStyle->SetFrameFillColor(10);
  gStyle->SetGridColor(0);
  gStyle->SetGridStyle(0);
  gStyle->SetGridWidth(0);
  //gStyle->SetFillColor(4);
  TGaxis::SetMaxDigits(4);
  gStyle->SetPadTopMargin(0.15);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadRightMargin(0.14);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetLabelSize(0.07,"X");
  gStyle->SetLabelSize(0.07,"Y");
  gStyle->SetTitleSize(0.07,"X");
  gStyle->SetTitleSize(0.07,"Y");

  const Int_t NRGBs = 5;
  const Int_t NCont = 255;
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t reds[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t greens[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blues[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  Int_t  FI = TColor::CreateGradientColorTable(NRGBs, stops, reds,
      greens, blues, NCont);
  gStyle->SetNumberContours(NCont);
  //**************************************************************************************

  cout << "Start to Read Trees!" << endl;

  StStrangenessAna *mStrangenessAna = new StStrangenessAna(energy,X_flag,List,start_event,stop_event,mode);
  mStrangenessAna->Init();
  mStrangenessAna->Make();
  mStrangenessAna->Finish();

//  cout << "End of the Calculation!!" << endl;
}

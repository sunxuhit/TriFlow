#include <TSystem>

void Resolution(const Int_t job_start = 123, const Int_t job_stop = 124, const Int_t energy = 2)
{
//  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
//  loadSharedLibraries();

  gSystem->Load("StTriFlowResolution");

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

  cout << "Start to Calculate Resolution!!" << endl;

  for(Int_t jobCounter = job_start; jobCounter < job_stop; jobCounter++)
  {
    cout << "job No. " << jobCounter << " is running!!" << endl;
    StTriFlowResolution *mTriFlowResolution = new StTriFlowResolution(jobCounter,energy);
    mTriFlowResolution->Init();
    mTriFlowResolution->getResolution();
    mTriFlowResolution->Finish();
    delete mTriFlowResolution;
  }

  cout << "End of the Calculation!!" << endl;
}

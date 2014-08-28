#include <TSystem>

void Resolution(const Int_t Energy = 4, const Int_t List = 4, const Long64_t StartEvent = 1000, const Long64_t StopEvent = 10024) // Energy = 0: 7GeV, 1: 11GeV, 2: 19GeV, 3: 27GeV, 4: 39GeV, 5: 62GeV, 6: 200 GeV
{
  gSystem->Load("AMPT_resolution");

  cout << "All libraries are loaded!!!!" << endl;
  cout << "Start to Calculate Resolution!!" << endl;

  AMPT_resolution *mRes = new AMPT_resolution(Energy,List,StartEvent,StopEvent);
  mRes->Init();
  mRes->Make();
  mRes->Finish();
  delete mRes;

  cout << "End of the Calculation!!" << endl;
}

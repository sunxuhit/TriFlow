#include <TSystem>

// root4star -b -q Flow.C\(4,0,0,1,1024\)
void Flow(const Int_t Energy = 4, const Int_t Mode = 0, const Int_t Screen = 0, const Int_t List = 0, const Long64_t StartEvent = 1000, const Long64_t StopEvent = 10024) // Energy = 0: 7GeV, 1: 11GeV, 2: 19GeV, 3: 27GeV, 4: 39GeV, 5: 62GeV | Mode = 0: Default, 1: String Melting | Screen = 0: 3mb, 1: 6mb
{
  gSystem->Load("AMPT_v3v2");

  cout << "All libraries are loaded!!!!" << endl;
  cout << "Start to Calculate Flow!!" << endl;

  AMPT_v3v2 *mFlow = new AMPT_v3v2(Energy,Mode,List,StartEvent,StopEvent);
  mFlow->Init();
  mFlow->Make();
  mFlow->Finish();
  delete mFlow;

  cout << "End of the Calculation!!" << endl;
}

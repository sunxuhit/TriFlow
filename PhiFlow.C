#include <TSystem>

// root4star -b -q PhiFlow.C\(4,0,0,1,1024\)
void PhiFlow(const Int_t Energy = 4, const Int_t Mode = 0, const Int_t List = 0, const Long64_t StartEvent = 1000, const Long64_t StopEvent = 10024, const Int_t Flag_ME = 0) // Energy = 0: 7GeV, 1: 11GeV, 2: 19GeV, 3: 27GeV, 4: 39GeV, 5: 62GeV | Mode = 0: Default, 1: String Melting | Flag_ME = 0: Same Event, 1: Mixed Event
{
  gSystem->Load("AMPT_phi");

  cout << "All libraries are loaded!!!!" << endl;
  cout << "Start to Calculate Flow!!" << endl;

  AMPT_phi *mPhiFlow = new AMPT_phi(Energy,Mode,List,StartEvent,StopEvent,Flag_ME);
  mPhiFlow->Init();
  mPhiFlow->Make();
  mPhiFlow->Finish();
  delete mPhiFlow;

  cout << "End of the Calculation!!" << endl;
}

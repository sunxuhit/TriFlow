
#include <TSystem>

class StMaker;
class StChain;
class StPicoDstMaker;


StChain *chain;
//void TriFlow(const Char_t *inputFile="/project/projectdirs/star/xusun/OutPut/AuAu200GeV/List/run_list/200GeV_511.list", const Int_t jobCounter = 511, const Int_t Mode = 7, const Int_t energy = 0)
void TriFlow(const Char_t *inputFile="/project/projectdirs/star/xusun/OutPut/AuAu39GeV/List/run_list/39GeV_511.list", const Int_t jobCounter = 511, const Int_t Mode = 7, const Int_t energy = 1)
//void TriFlow(const Char_t *inputFile="./List/27GeV/run_list/27GeV_134.list", const Int_t jobCounter = 134, const Int_t Mode = 3, const Int_t energy = 2)
{
//        Int_t nEvents = 10000000;
	Int_t nEvents = 50000;
	
        gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
	loadSharedLibraries();

	gSystem->Load("StRefMultCorr");
	gSystem->Load("StPicoDstMaker");
	gSystem->Load("StAlexPhiMesonEvent");
	gSystem->Load("StV0Event");
	gSystem->Load("StV0TofCorrection");
        gSystem->Load("StTriFlowMaker");
	gSystem->Load("StRunIdEventsDb");
	gSystem->Load("StCombPID");

	chain = new StChain();

	StPicoDstMaker *picoMaker = new StPicoDstMaker(0,inputFile,"picoDst");

        StTriFlowMaker *TriFlowMaker = new StTriFlowMaker("TriFlow",picoMaker,jobCounter,Mode,energy);

	chain->Init();
	cout<<"chain->Init();"<<endl;
	int total = picoMaker->chain()->GetEntries();
        cout << " Total entries = " << total << endl;
        if(nEvents>total) nEvents = total;
	for (Int_t i=0; i<nEvents; i++)
	{
	  if(i != 0 && i%50 == 0)
	  {
	    cout << "." << flush;
	  }
	  if(i != 0 && i%500 == 0)
	  {
	    Float_t event_percent = 100.0*i/nEvents;
	    cout << " " << i << "(" << event_percent << "%)" << "\n" << " ==> Processing data " << flush;
	  }
		
	  chain->Clear();
	  int iret = chain->Make(i);
		
	  if (iret)
	  { 
	    cout << "Bad return code!" << iret << endl;
	    break;
	  }

	  total++;
		
	}
	
	cout << "." << flush;
	cout << " " << nEvents << "(" << 100 << "%)";
	cout << endl;
	cout << "****************************************** " << endl;
	cout << "Work done... now its time to close up shop!"<< endl;
	cout << "****************************************** " << endl;
	chain->Finish();
	cout << "****************************************** " << endl;
	cout << "total number of events  " << nEvents << endl;
	cout << "****************************************** " << endl;
	
	delete chain;
	
	
}

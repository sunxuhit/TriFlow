#include "TF1.h"
#include "TCanvas.h"
#include "./student_t_combPID.h"

void calKaonPosition()
{
  Double_t x1, x2, y1, y2;
  setInitValues(2.6,-3.0,0.494*0.494);
  x1 = getNewX();
  y1 = getNewY();
  setInitValues(2.6,3.0,0.494*0.494);
  x2 = getNewX();
  y2 = getNewY();
  TF1 *cutline = new TF1("cutline",getLine,-3.0,3.0,4);
  cutline->FixParameter(0,x1);
  cutline->FixParameter(1,y1);
  cutline->FixParameter(2,x2);
  cutline->FixParameter(3,y2);

  Float_t x_val = cutline->GetX(0.3375);
  cout << x_val << endl;
  TCanvas *c1 = new TCanvas("c1","c1",10,10,800,800);
  cutline->Draw("L");
}

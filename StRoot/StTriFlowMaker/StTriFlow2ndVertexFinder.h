#ifndef StTriFlow2ndVertexFinder_h
#define StTriFlow2ndVertexFinder_h

#include "StPhysicalHelixD.hh"
#include "StThreeVectorF.hh"
#include "TLorentzVector.h"

class StTriFlow2ndVertexFinder
{
  public:
    StTriFlow2ndVertexFinder();
    ~StTriFlow2ndVertexFinder();

    void Find2ndVertex(StPhysicalHelixD,StPhysicalHelixD,StThreeVectorF,Float_t,Float_t,TLorentzVector &,TLorentzVector &,Float_t &,Float_t &, StThreeVectorF &, Float_t &,Int_t);

  private:
    // Alex's Secondary Vertex Finder
    Int_t fDCA_Helix_Estimate(StPhysicalHelixD helixA, StPhysicalHelixD helixB, Float_t &pathA, Float_t &pathB, Float_t &dcaAB);
    Int_t fCross_points_Circles(Double_t x1, Double_t y1, Double_t r1, Double_t x2, Double_t y2, Double_t r2, Double_t &x1_c, Double_t &y1_c, Double_t &x2_c, Double_t &y2_c);
    StThreeVectorF calcVertexAnalytical(StThreeVectorF &base1, StThreeVectorF &dir1, StThreeVectorF &base2, StThreeVectorF &dir2);
    StThreeVectorF calculatePointOfClosestApproach(StThreeVectorF &base1, StThreeVectorF &dir1, StThreeVectorF &base2, StThreeVectorF &dir2);
    StThreeVectorF calculateCrossPoint(StThreeVectorF &base1, StThreeVectorF &dir1, StThreeVectorF &base2, StThreeVectorF &dir2);
    Double_t calcDeterminant(StThreeVectorF& v1,StThreeVectorF& v2,StThreeVectorF& v3);
    Double_t calculateMinimumDistance(StThreeVectorF &base1, StThreeVectorF &dir1, StThreeVectorF &base2, StThreeVectorF &dir2);
    Double_t calculateMinimumDistanceStraightToPoint(StThreeVectorF &base, StThreeVectorF &dir, StThreeVectorF &point);
    void fHelixABdca_start_params(StPhysicalHelixD helixA, StPhysicalHelixD helixB, Float_t &pathA, Float_t &pathB, Float_t &dcaAB, Float_t path_in_A, Float_t path_in_B);
    void fHelixABdca(StPhysicalHelixD helixA, StPhysicalHelixD helixB, Float_t &pathA, Float_t &pathB, Float_t &dcaAB);
    void fHelixAtoPointdca_start_params(StThreeVectorF space_vec, StPhysicalHelixD helixA, Float_t &pathA, Float_t &dcaAB, Float_t path_in_A);
    void fHelixAtoPointdca(StThreeVectorF space_vec, StPhysicalHelixD helixA, Float_t &pathA, Float_t &dcaAB);

  ClassDef(StTriFlow2ndVertexFinder,1)
};
#endif

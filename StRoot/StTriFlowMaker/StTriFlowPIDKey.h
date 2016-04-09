#ifndef StTriFlowPIDKey_h
#define StTriFlowPIDKey_h
#include "Rtypes.h"
#include <map>

class TH2F;

class PIDKey
{
  public:
    Int_t Centrality;
    Int_t charge;
    Int_t eta_gap;
    Int_t eta_pos_neg;
    Int_t pt;

    PIDKey() {}
    PIDKey(Int_t i, Int_t j, Int_t k, Int_t l, Int_t m) : Centrality(i), charge(j), eta_gap(k), eta_pos_neg(l), pt(m) {}

    bool operator<(const PIDKey &r) const
    { // "r == right"
      if(Centrality == r.Centrality)
      {
        if(charge == r.charge)
	{
	  if(eta_gap == r.eta_gap)
	  {
	    if(eta_pos_neg == r.eta_pos_neg)
	    {
	      return pt < r.pt;
	    }
	    else
	    {
	      return eta_pos_neg < r.eta_pos_neg;
	    }
	  }
	  else
	  {
	    return eta_gap < r.eta_gap;
	  }
	}
	else
	{
	  return charge < r.charge;
	}
      }
      else
      {
        return Centrality < r.Centrality;
      }
    }

    void print()
    {
      cout << "Centrality = " << PIDKey::Centrality << endl;
      cout << "charge = " << PIDKey::charge << endl;
      cout << "eta_gap = " << PIDKey::eta_gap << endl;
      cout << "eta_pos_neg = " << PIDKey::eta_pos_neg << endl;
      cout << "pt = " << PIDKey::pt << endl;
    }

    ClassDef(PIDKey, 1)
};

// typedef std::map<PIDKey, TH2F*> TH2FMap;

#endif

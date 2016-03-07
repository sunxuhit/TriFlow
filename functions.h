#include "TMath.h"

Double_t PolyBreitWigner(Double_t *x_val, Double_t *par) 
{
  Double_t x = x_val[0];
  Double_t m0 = par[0];
  Double_t Gamma = par[1];
  Double_t Norm = par[2];

  Double_t denom = 2.0*TMath::Pi()*((x-m0)*(x-m0)+Gamma*Gamma/4.0);
  Double_t BW = Norm*Gamma/denom;

  Double_t Poly = par[3] + par[4]*x;

  Double_t y = BW + Poly;

  return y;
}

Double_t Poly(Double_t *x_val, Double_t *par) 
{
  Double_t x = x_val[0];
  Double_t y = par[0] + par[1]*x;

  return y;
}

Double_t BreitWigner(Double_t *x_val, Double_t *par) 
{
  Double_t x = x_val[0];
  Double_t m0 = par[0];
  Double_t Gamma = par[1];
  Double_t Norm = par[2];

  Double_t denom = 2.0*TMath::Pi()*((x-m0)*(x-m0)+Gamma*Gamma/4.0);
  Double_t BW = Norm*Gamma/denom;

  return BW;
}

Double_t flow(Double_t *x_val, Double_t *par)
{
  Double_t x, y;
  Double_t Ampl, v3, order;
  x = x_val[0];
  Ampl = par[0];
  v3 = par[1];
  order = par[2];

  y = Ampl*(1.0 + 2.0*v3*TMath::Cos(order*x));

  return y;
}

Double_t Gaussion(Double_t *x_val, Double_t *par)
{
  Double_t x, y;
  x = x_val[0];

  Double_t mu = par[0];
  Double_t sigma = par[1];
  Double_t norm = par[2];

  y = norm*TMath::Exp(-1.0*(x-mu)*(x-mu)/(2.0*sigma*sigma))/(sigma*TMath::Sqrt(2*TMath::Pi()));

  return y;
}

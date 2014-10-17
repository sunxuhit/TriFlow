Double_t student_t_1d_single(Double_t *x_val, Double_t*par)
{
  Double_t x, z;
  Double_t z_x;

  Double_t    nu_x;
  Double_t  mean_x;
  Double_t sigma_x;
  Double_t    norm;
  Double_t      pi;

  pi = TMath::Pi();

  x = x_val[0];

  nu_x = par[0];
  mean_x = par[1];
  sigma_x = par[2];
  norm = par[3];

//------------------------------------------------------------------------------------------------
// student_t function of pion start

// x component start
  z_x = (TMath::Gamma((nu_x+1.0)/2.0)/(TMath::Gamma(nu_x/2.0)*TMath::Sqrt(pi*nu_x)*sigma_x))*TMath::Power((1.0+TMath::Power((x-mean_x)/sigma_x,2.0)/nu_x),(-(nu_x+1.0)/2.0));
// x component stop

// student_t function of pion

  z= norm*z_x;
  
// student_t function of pion stop
//------------------------------------------------------------------------------------------------

  return z;
}

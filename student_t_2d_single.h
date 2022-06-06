
Double_t student_t_2d_single(Double_t *x_val, Double_t*par)
{
  Double_t x, y, z;
  Double_t z_x;
  Double_t z_y;

  Double_t    nu_x;
  Double_t    nu_y;
  Double_t  mean_x;
  Double_t  mean_y;
  Double_t sigma_x;
  Double_t sigma_y;
  Double_t    norm;
  Double_t      pi;

  pi = TMath::Pi();

  x = x_val[0];
  y = x_val[1];

  nu_x= par[0];
  nu_y= par[1];
  mean_x= par[2];
  mean_y= par[3];
  sigma_x= par[4];
  sigma_y= par[5];
  norm= par[6];

//------------------------------------------------------------------------------------------------
// student_t function start

// x component start
  z_x = (TMath::Gamma((nu_x+1.0)/2.0)/(TMath::Gamma(nu_x/2.0)*TMath::Sqrt(pi*nu_x)*sigma_x))*TMath::Power((1.0+TMath::Power((x-mean_x)/sigma_x,2.0)/nu_x),(-(nu_x+1.0)/2.0));
// x component stop

// y component start
  z_y= (TMath::Gamma((nu_y+1.0)/2.0)/(TMath::Gamma(nu_y/2.0)*TMath::Sqrt(pi*nu_y)*sigma_y))*TMath::Power((1.0+TMath::Power((y-mean_y)/sigma_y,2.0)/nu_y),(-(nu_y+1.0)/2.0));
// y component stop

// student_t function of pion

  z= norm*z_x*z_y;
  
// student_t function of pion stop
//------------------------------------------------------------------------------------------------

  return z;
}

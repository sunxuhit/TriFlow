Double_t student_t_1d_double(Double_t *x_val, Double_t*par)
{
  Double_t x, y, z, z_pion, z_kaon;
  Double_t z_pion_x, z_kaon_x;

  Double_t    nu_x_pion,    nu_x_kaon;
  Double_t  mean_x_pion,  mean_x_kaon;
  Double_t sigma_x_pion, sigma_x_kaon;
  Double_t    norm_pion,    norm_kaon;
  Double_t pi;

  pi = TMath::Pi();

  x = x_val[0];

// pion start
  nu_x_pion = par[0];
  mean_x_pion = par[1];
  sigma_x_pion = par[2];
  norm_pion = par[3];
// pion stop

// kaon start
  nu_x_kaon = par[4];
  mean_x_kaon = par[5];
  sigma_x_kaon = par[6];
  norm_kaon = par[7];
// kaon stop

//------------------------------------------------------------------------------------------------
// student_t function of pion start

// x component start
  z_pion_x = (TMath::Gamma((nu_x_pion+1.0)/2.0)/(TMath::Gamma(nu_x_pion/2.0)*TMath::Sqrt(pi*nu_x_pion)*sigma_x_pion))*TMath::Power((1.0+TMath::Power((x-mean_x_pion)/sigma_x_pion,2.0)/nu_x_pion),(-(nu_x_pion+1.0)/2.0));
// x component stop

// student_t function of pion

  z_pion = norm_pion*z_pion_x;
  
// student_t function of pion stop
//------------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------------
// student_t function of kaon start

// x component start
  z_kaon_x = (TMath::Gamma((nu_x_kaon+1.0)/2.0)/(TMath::Gamma(nu_x_kaon/2.0)*TMath::Sqrt(pi*nu_x_kaon)*sigma_x_kaon))*TMath::Power((1.0+TMath::Power((x-mean_x_kaon)/sigma_x_kaon,2.0)/nu_x_kaon),(-(nu_x_kaon+1.0)/2.0));
// x component stop
  
// student_t function of kaon 

  z_kaon = norm_kaon*z_kaon_x;
  
// student_t function of kaon stop
//------------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------------
// student_t function start
  z = z_pion + z_kaon;
// student_t function stop
//------------------------------------------------------------------------------------------------

  return z;
}

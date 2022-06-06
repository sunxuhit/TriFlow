Double_t student_t_1d_fit(Double_t *x_val, Double_t*par)
{
  Double_t x, z, z_pion, z_kaon, z_proton;
  Double_t z_pion_x, z_kaon_x, z_proton_x;

  Double_t    nu_x_pion,    nu_x_kaon,    nu_x_proton;
  Double_t  mean_x_pion,  mean_x_kaon,  mean_x_proton;
  Double_t sigma_x_pion, sigma_x_kaon, sigma_x_proton;
  Double_t    norm_pion,    norm_kaon,    norm_proton;
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

// proton start
  nu_x_proton = par[8];
  mean_x_proton = par[9];
  sigma_x_proton = par[10];
  norm_proton = par[11];
// proton stop

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
// student_t function of proton start

// x component start
  z_proton_x = (TMath::Gamma((nu_x_proton+1.0)/2.0)/(TMath::Gamma(nu_x_proton/2.0)*TMath::Sqrt(pi*nu_x_proton)*sigma_x_proton))*TMath::Power((1.0+TMath::Power((x-mean_x_proton)/sigma_x_proton,2.0)/nu_x_proton),(-(nu_x_proton+1.0)/2.0));
// x component stop
  
// student_t function of proton 

  z_proton = norm_proton*z_proton_x;
  
// student_t function of proton stop
//------------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------------
// student_t function start
  z = z_pion + z_kaon + z_proton;
//  z = z_proton;
// student_t function stop
//------------------------------------------------------------------------------------------------

  return z;
}

Double_t student_t_2d_fit(Double_t *x_val, Double_t*par)
{
  Double_t x, y, z, z_pion, z_kaon, z_proton;
  Double_t z_pion_x, z_kaon_x, z_proton_x;
  Double_t z_pion_y, z_kaon_y, z_proton_y;

  Double_t    nu_x_pion,    nu_x_kaon,    nu_x_proton;
  Double_t    nu_y_pion,    nu_y_kaon,    nu_y_proton;
  Double_t  mean_x_pion,  mean_x_kaon,  mean_x_proton;
  Double_t  mean_y_pion,  mean_y_kaon,  mean_y_proton;
  Double_t sigma_x_pion, sigma_x_kaon, sigma_x_proton;
  Double_t sigma_y_pion, sigma_y_kaon, sigma_y_proton;
  Double_t    norm_pion,    norm_kaon,    norm_proton;
  Double_t pi;

  pi = TMath::Pi();

  x = x_val[0];
  y = x_val[1];

// pion start
  nu_x_pion = par[0];
  nu_y_pion = par[1];
  mean_x_pion = par[2];
  mean_y_pion = par[3];
  sigma_x_pion = par[4];
  sigma_y_pion = par[5];
  norm_pion = par[6];
// pion stop

// kaon start
  mean_x_kaon = par[7];
  mean_y_kaon = par[8];
  sigma_x_kaon = par[9];
  sigma_y_kaon = par[10];
  norm_kaon = par[11];
// kaon stop

// proton start
  nu_x_proton = par[12];
  nu_y_proton = par[13];
  mean_x_proton = par[14];
  mean_y_proton = par[15];
  sigma_x_proton = par[16];
  sigma_y_proton = par[17];
  norm_proton = par[18];
// proton stop

  // kaon 
  nu_x_kaon = (2.0*par[0]+par[12])/3.0;
  nu_y_kaon = (2.0*par[1]+par[13])/3.0;

//------------------------------------------------------------------------------------------------
// student_t function of pion start

// x component start
  z_pion_x = (TMath::Gamma((nu_x_pion+1.0)/2.0)/(TMath::Gamma(nu_x_pion/2.0)*TMath::Sqrt(pi*nu_x_pion)*sigma_x_pion))*TMath::Power((1.0+TMath::Power((x-mean_x_pion)/sigma_x_pion,2.0)/nu_x_pion),(-(nu_x_pion+1.0)/2.0));
// x component stop

// y component start
  z_pion_y = (TMath::Gamma((nu_y_pion+1.0)/2.0)/(TMath::Gamma(nu_y_pion/2.0)*TMath::Sqrt(pi*nu_y_pion)*sigma_y_pion))*TMath::Power((1.0+TMath::Power((y-mean_y_pion)/sigma_y_pion,2.0)/nu_y_pion),(-(nu_y_pion+1.0)/2.0));
// y component stop

// student_t function of pion

  z_pion = norm_pion*z_pion_x*z_pion_y;
  
// student_t function of pion stop
//------------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------------
// student_t function of kaon start

// x component start
  z_kaon_x = (TMath::Gamma((nu_x_kaon+1.0)/2.0)/(TMath::Gamma(nu_x_kaon/2.0)*TMath::Sqrt(pi*nu_x_kaon)*sigma_x_kaon))*TMath::Power((1.0+TMath::Power((x-mean_x_kaon)/sigma_x_kaon,2.0)/nu_x_kaon),(-(nu_x_kaon+1.0)/2.0));
// x component stop

// y component start
  z_kaon_y = (TMath::Gamma((nu_y_kaon+1.0)/2.0)/(TMath::Gamma(nu_y_kaon/2.0)*TMath::Sqrt(pi*nu_y_kaon)*sigma_y_kaon))*TMath::Power((1.0+TMath::Power((y-mean_y_kaon)/sigma_y_kaon,2.0)/nu_y_kaon),(-(nu_y_kaon+1.0)/2.0));
// y component stop
  
// student_t function of kaon 

  z_kaon = norm_kaon*z_kaon_x*z_kaon_y;
  
// student_t function of kaon stop
//------------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------------
// student_t function of proton start

// x component start
  z_proton_x = (TMath::Gamma((nu_x_proton+1.0)/2.0)/(TMath::Gamma(nu_x_proton/2.0)*TMath::Sqrt(pi*nu_x_proton)*sigma_x_proton))*TMath::Power((1.0+TMath::Power((x-mean_x_proton)/sigma_x_proton,2.0)/nu_x_proton),(-(nu_x_proton+1.0)/2.0));
// x component stop

// y component start
  z_proton_y = (TMath::Gamma((nu_y_proton+1.0)/2.0)/(TMath::Gamma(nu_y_proton/2.0)*TMath::Sqrt(pi*nu_y_proton)*sigma_y_proton))*TMath::Power((1.0+TMath::Power((y-mean_y_proton)/sigma_y_proton,2.0)/nu_y_proton),(-(nu_y_proton+1.0)/2.0));
// y component stop
  
// student_t function of proton 

  z_proton = norm_proton*z_proton_x*z_proton_y;
  
// student_t function of proton stop
//------------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------------
// student_t function start
  z = z_pion + z_kaon + z_proton;
//  z = z_pion;
// student_t function stop
//------------------------------------------------------------------------------------------------

  return z;
}

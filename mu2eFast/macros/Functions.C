Double_t splitgaus(Double_t *x, Double_t *par) {
  Double_t retval;
  Double_t core;
  Double_t tail;
  Float_t xval = x[0];
  if(xval > par[1]) {
    core = exp(-0.5*pow((xval-par[1])/par[2],2))/par[2];
    tail = par[4]*exp(-0.5*pow((xval-par[1])/par[5],2))/par[5];
  } else {
    core = exp(-0.5*pow((xval-par[1])/par[3],2))/par[3];
    tail = (1/par[2]-1/par[3]+par[4]/par[5])*exp(-0.5*pow((xval-par[1])/par[6],2));
  }
  retval = par[0]*0.398942*(core+tail);
// add a tail Gaussian
  return retval;
}

Double_t doublegaus(Double_t *x, Double_t *par) {
  Double_t retval;
  Float_t xval = x[0];
  retval = par[0]*0.398942*( par[4]*exp(-0.5*pow((xval-par[1])/par[2],2))/par[2] + 
    (1.0-par[4])*exp(-0.5*pow((xval-par[1])/par[3],2))/par[3]);
  return retval;
}

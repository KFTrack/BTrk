#include <Rtypes.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TTree.h>
#include <TCut.h>
#include <TLegend.h>
#include <TString.h>
#include <TStyle.h>

Double_t doublegaus(Double_t *x, Double_t *par) {
  Float_t xval = x[0];
  double gaus1 = exp(-0.5*pow(xval/par[1],2))/par[1];
  double gaus2 = exp(-0.5*pow(xval/par[2],2))/par[2];
  Double_t retval = par[0]*0.39894*( (1.0-par[3])*gaus1 + par[3]*gaus2);
  return retval;
}

Double_t tripplegaus(Double_t *x, Double_t *par) {
  Float_t xval = x[0];
  double gaus1 = exp(-0.5*pow(xval/par[1],2))/par[1];
  double gaus2 = exp(-0.5*pow(xval/par[2],2))/par[2];
  double gaus3 = exp(-0.5*pow(xval/par[3],2))/par[3];
  Double_t retval = par[0]*0.39894*( (1.0-par[4]-par[5])*gaus1 + par[4]*gaus2 + par[5]*gaus3 );
  return retval;
}

Double_t Rudi_scatter(Double_t* x, Double_t *par) {
// formula from R. Fruhwirth
// 0=overall normalization
// 1=Gaussian sigma
// 2,3 = tail min,max
// 4=tail fraction
  Double_t retval;
  Float_t xval = x[0];
  Double_t gaus = 0.398942*exp(-0.5*pow(xval/par[1],2))/par[1];
  double asq = par[2]*par[2];
  double bsq = par[3]*par[3];
  double xsq = xval*xval;
  Double_t tail = 0.318310*asq*(sqrt(bsq-xsq)/((asq+bsq)*(xsq+asq)) + acos(sqrt((xsq+asq)/(bsq+asq)))/pow(xsq+asq,1.5) );
  retval = par[0]*( (1.0-par[4])*gaus + par[4]*tail);
  return retval;
}

void mu2e_EGS(TCanvas* can, TTree* tree, const char* cpage="eloss" ) {
  TString page(cpage);
  TF1* dgau = new TF1("dgau",doublegaus,-1.,1.,4);
  TF1* tgau = new TF1("tgau",tripplegaus,-1.,1.,6);
  TF1* scat = new TF1("scat",Rudi_scatter,-1.,1.,5);
  
  TCut one_e("np==1");
  TCut noloss("E1-E0>-0.02");
  if( page == "eloss"){
    TH1F* eloss = new TH1F("eloss","energy loss",200,-0.1, 0.0);
    TH1F* eloss2 = new TH1F("eloss2","energy loss",200,-1, 0.0);

    tree->Project("eloss","E1-E0",one_e);
    tree->Project("eloss2","E1-E0",one_e);
    can->Clear();
    can->Divide(1,2);
    can->cd(1);
    gPad->SetLogy();
    eloss->Draw();
    can->cd(2);
    gPad->SetLogy();
    eloss2->Draw();
  } else if (page == "scatter"){
    
    gStyle->SetOptFit(1111);
    TH1F* du = new TH1F("du","U angle change",200,-0.02,0.02);
    tree->Project("du","Ux",one_e+noloss);
    can->Clear();
    can->Divide(1,3);
    can->cd(1);
    gPad->SetLogy();
    double integral = du->GetEntries()*du->GetBinWidth(1);
    dgau->SetParameters(integral,du->GetRMS()/1.5,du->GetRMS()*2.0,0.05);
    du->Fit("dgau");
    can->cd(2);
    TH1F* du2 = new TH1F(*du);
    gPad->SetLogy();
    scat->SetParameters(integral,du2->GetRMS()/1.5,0.005,0.02,0.05);
    du2->Fit("scat");
    can->cd(3);
    TH1F* du3 = new TH1F(*du);
    gPad->SetLogy();
    tgau->SetParameters(integral,du3->GetRMS()/1.5,du3->GetRMS()*2,du3->GetRMS()*5,0.1,0.01);
    du3->Fit("tgau");
    
  }
}
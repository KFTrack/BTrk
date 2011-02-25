#include <Rtypes.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TTree.h>
#include <TCut.h>
#include <TLegend.h>
#include <TString.h>
#include <TStyle.h>
#include <TRandom3.h>
#include <TMath.h>
#include "mu2eFast/include/VavRan.h"

Double_t VavilovFun(Double_t *x, Double_t *par ) {
  return par[0]*TMath::Vavilov(x[0],par[1],par[2]);
}

void testVavilov(double rkappa=2.0, double beta2=0.5,unsigned ntrials=100000) {
  unsigned nbins(200);
  double range[2] = {-5,5};
  TH1F* vavh = new TH1F("vav","VavRan test",nbins,range[0],range[1]);
  TF1* vfun = new TF1("Vavilov",VavilovFun,range[0],range[1],3);
  vfun->SetParameter(0,ntrials*(range[1]-range[0])/nbins);
  vfun->SetParameter(1,rkappa);
  vfun->SetParameter(2,beta2);
  TRandom3* tr = new TRandom3();
  for(unsigned itry=0;itry<ntrials;itry++){
    double vran = VavRan::Instance()->gen(rkappa,beta2,tr->Rndm());
    vavh->Fill(vran);
  }  
  TCanvas* can = new TCanvas("can");
  can->Clear();
//  can->Divide(1,2);
//  can->cd(1);
//  vavh->Draw();
//  vfun->Draw("same");
//  can->cd(2);
  vavh->Fit(vfun);
  can->SaveAs("testVavilov.png");
}
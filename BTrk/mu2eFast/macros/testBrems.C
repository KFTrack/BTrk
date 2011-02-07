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

static const double ft(4.0/3.0);
static const double ot(1.0/3.0);

// spectrum of photons (dN/dy) in terms of the energy fraction of the photon
Double_t dNdy(Double_t* x, Double_t* par){
  return par[0]*( ft*(1.0/x[0] - 1.0) + x[0]);
}

// total photon energy fraction below a cut on the fractional energy.  This is the
// integral of the above times the photon energy
Double_t ELow(Double_t* x, Double_t* par) {
  return par[0]*( ft*(x[0]-0.5*x[0]*x[0]) + ot*x[0]*x[0]*x[0]);
}

// # of photons emitted with fractional energy above a cut.  This is the integral
// of the above
Double_t NGam(Double_t* x, Double_t* par) {
  return par[0]*( ft*(x[0]-log(x[0])-1.0) + 0.5*(1-x[0]*x[0]) );
}

void testBrems(double rfrac,unsigned ntrials=100000) {
  unsigned nbins(500);
  TH1F* elowh = new TH1F("elowh","Low-E Brems energy fraction vs cut",nbins,0,1);
  TF1* elow = new TF1("Elow",ELow,1e-3,1,1);
  elow->SetParameter(0,rfrac);
  elowh->SetMaximum(rfrac);

  TH1F* dndyh = new TH1F("dndyh","Brems photon energy spectrum",nbins,0,1);
  TF1* dndy = new TF1("dNdy",dNdy,1e-3,1,1);
  dndy->SetParameter(0,rfrac);
  dndyh->SetMaximum(10);

  TH1F* ngamh = new TH1F("ngamh","Average # high-E Brems photons vs cut",nbins,0,1);
  TF1* ngam = new TF1("NGam",NGam,1e-3,1,1);
  ngam->SetParameter(0,rfrac);
  ngamh->SetMaximum(10*rfrac);
  
  TCanvas* can1 = new TCanvas("can1");
  can1->Divide(2,2);
  can1->cd(1);
  elowh->Draw();
  elow->Draw("same");

  can1->cd(2);
  gPad->SetLogy();
  dndyh->Draw();
  dndy->Draw("same");

  can1->cd(3);
  gPad->SetLogy();
  ngamh->Draw();
  ngam->Draw("same");
  
  
// simulate a bunch of brems
  TH1F* efrac = new TH1F("efrac","Fractional Brems energy loss",nbins,0.0,1.0);
  TH1F* efracc = new TH1F("efracc","Fractional Brems energy loss, with cut",nbins,0.0,1.0);
  TH1F* efracc2 = new TH1F("efracc2","Fractional Brems energy loss, sampled, with cut",nbins,0.0,1.0);
  TProfile* efracp = new TProfile("efracp","N Brems",nbins,0.0,1.0,0,1e5);
  TH1F* ngh = new TH1F("ngh","N gammas above cut",10,-0.5,9.5);
  TRandom3* tr = new TRandom3();
  double ybin= 1.0/nbins;
  double ycut=0.05;
// set integration range
  dndy->SetRange(ycut,1.0);
  for(unsigned itrial=0;itrial<ntrials;itrial++){
//  First simulate each bin of the photon energy spectrm
    double eloss(0.0);
    double elossc = ELow(&ycut,&rfrac);
    double elossc2 = elossc;
    double ng = NGam(&ycut,&rfrac);
    unsigned ngn = tr->Poisson(ng);
    ngh->Fill((float)ngn);
    for(unsigned igh=0;igh<ngn;igh++){
      double y = dndy->GetRandom();
      elossc2 += y;
    }
    efracc2->Fill(elossc2);
// bin-by-bin test
    for(unsigned ibin=0;ibin<nbins;ibin++){
      double y = (ibin+0.5)*ybin;
      double dngam = ybin*dNdy(&y,&rfrac);
      efracp->Fill(y,dngam/ybin);
      unsigned nga = tr->Poisson(dngam);
      eloss += nga*y;
      if(y>ycut)elossc += nga*y;
//      if(ng > 0) std::cout << "ng =" << ng << " y=" << y << " eloss = " << eloss << std::endl;
    }
    if(eloss>1.0)eloss=1.0;
    efrac->Fill(eloss);
    efracc->Fill(elossc);
  }
  
  can1->cd(4);
  gPad->SetLogy();
  efracp->Draw();
  
  
  TCanvas* can2 = new TCanvas("can2");
  can2->Divide(2,2);
  can2->cd(1);
  gPad->SetLogy();
  efrac->Draw();
  
  can2->cd(2);
  gPad->SetLogy();
  efracc->Draw();

  can2->cd(3);
  gPad->SetLogy();
  ngh->Draw();
  
  can2->cd(4);
  gPad->SetLogy();
  efracc2->Draw();
  
}
  

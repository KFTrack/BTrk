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
#include <TGraph.h>
#include <TStyle.h>
#include <TSpline.h>
#include <TRandom3.h>
#include <TMath.h>
#include "mu2eFast/include/VavRan.h"

static const double ft(4.0/3.0);
static const double ot(1.0/3.0);

// spectrum of photons (dN/dy) in terms of the energy fraction of the photon
Double_t dNdy(Double_t* x, Double_t* par){
  return par[0]*( ft*(1.0/x[0] - 1.0) + x[0]);
}
// # of photons emitted with fractional energy above a given fraction.  This is the integral
// of the above from the faction to 1.0.
Double_t NGam(Double_t* x, Double_t* par) {
  return par[0]*( ft*(x[0]-log(x[0])-1.0) + 0.5*(1-x[0]*x[0]) );
}


// total photon energy fraction below a cut on the fractional energy.  This is the
// integral of the above times the photon energy
Double_t ELow(Double_t* x, Double_t* par) {
  return par[0]*( ft*(x[0]-0.5*x[0]*x[0]) + ot*x[0]*x[0]*x[0]);
}


void testBrems(double rfrac,double ycut=0.05,unsigned ntrials=100000) {
  unsigned nbins(200);
  TH1F* elowh = new TH1F("elowh","Low-E Brems energy fraction vs cut",nbins,0,1);
  TF1* elow = new TF1("Elow",ELow,1e-3,1,1);
  elow->SetParameter(0,rfrac);
  elowh->SetMaximum(rfrac);

  TH1F* dndyh = new TH1F("dndyh","Brems photon number spectrum",nbins,0,1);
  TF1* dndy = new TF1("dNdy",dNdy,1e-3,1,1);
  dndy->SetParameter(0,rfrac);
  dndyh->SetMaximum(10);

  TH1F* ngamh = new TH1F("ngamh","Average # high-E Brems photons vs cut",nbins,0,1);
  TF1* ngam = new TF1("NGam",NGam,1e-3,1,1);
  ngam->SetParameter(0,rfrac);
  ngamh->SetMaximum(10*rfrac);
  
  // sample the normalized integrated N gamma distribution in order to invert it
  double ybin= 1.0/nbins;
  std::vector<Double_t> xv;
  std::vector<Double_t> yv;
  xv.reserve(nbins);
  yv.reserve(nbins);
  double norm = NGam(&ycut,&rfrac);
  double yval(0.0);
  for(unsigned ibin=0;ibin<nbins;ibin++){
    if(yval>=ycut){
      xv.push_back(yval);
      yv.push_back((norm-NGam(&yval,&rfrac))/norm);
    }
    yval += ybin;
  }
//  TSpline3* sp3 = new TSpline3("spline",&xv.front(),&yv.front(),xv.size(),"",dNdy(&ycut,&rfrac)/norm,rfrac/norm);
//  TSpline3* sp3_inv = new TSpline3("invspline",&yv.front(),&xv.front(),xv.size(),"",norm/dNdy(&ycut,&rfrac),norm/rfrac);

  TGraph* g = new TGraph(xv.size(),&xv.front(),&yv.front());
  TGraph* ginv = new TGraph(xv.size(),&yv.front(),&xv.front());
  TSpline3* sp3_inv = new TSpline3("invspline",ginv);
  TSpline3* sp3 = new TSpline3("spline",g);
  g->SetLineColor(kRed);
  ginv->SetLineColor(kRed);
//  sp3_inv->SetXmin(ycut);
//  sp3_inv->SetXmax(yv[0]);
  TH1F* hsp = new TH1F("hsp","spline fit to NGam",nbins,0,1);
  TH1F* hspinv = new TH1F("hspinv","spline fit to Inverse NGam",nbins,0,1);
  hsp->SetMaximum(1.0);
  hsp->SetMinimum(0.0);
  hspinv->SetMaximum(1.0);
  hspinv->SetMinimum(0.0);
  
  // convert the inverse spline into an evenly-binned spline, which is much
  // faster to compute
  std::vector<Double_t> ixv;
  std::vector<Double_t> iyv;
  ixv.reserve(nbins);
  iyv.reserve(nbins);
  double dx = 1.0/(nbins);
  double xval = (0.0);
  for(unsigned ibin=0;ibin<nbins;ibin++){
    ixv.push_back(xval);
    iyv.push_back(sp3_inv->Eval(xval));
    xval += dx;
  }
  TGraph* ginv2 = new TGraph(ixv.size(),&ixv.front(),&iyv.front());
  TSpline3* sp3_inv2 = new TSpline3("invspline2",ginv2);

  
// generate a random distribution
  TRandom3* tr = new TRandom3();
  TH1F* hspect = new TH1F("spect","N spectrum",nbins,0,1);
  for(unsigned itrial=0;itrial<ntrials;itrial++){
    double ran = tr->Rndm();
    double sx = sp3_inv2->Eval(ran);
    hspect->Fill(sx);
  }
  
  
  TCanvas* can1 = new TCanvas("can1");
  can1->Divide(2,3);
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
  
  can1->cd(4);
  hsp->Draw();
  sp3->Draw("same");
  g->Draw("same");
  
  can1->cd(5);
  hspinv->Draw();
  sp3_inv2->Draw("same");
  ginv2->Draw("same");

  can1->cd(6);
  gPad->SetLogy();
  TF1* newfun = new TF1(*dndy);
  newfun->SetName("newdNdy");
  cout << newfun->GetName() << endl;
  hspect->Fit(newfun);
  
  sp3_inv2->SaveAs("invngam.cc");
  
  return;
  
// simulate a bunch of brems
  TH1F* efrac = new TH1F("efrac","Fractional Brems energy loss",nbins,0.0,1.0);
  TH1F* efracc = new TH1F("efracc","Fractional Brems energy loss, with cut",nbins,0.0,1.0);
  TH1F* efracc2 = new TH1F("efracc2","Fractional Brems energy loss, sampled, with cut",nbins,0.0,1.0);
  TProfile* efracp = new TProfile("efracp","N Brems",nbins,0.0,1.0,0,1e5);
  TH1F* ngh = new TH1F("ngh","N gammas above cut",10,-0.5,9.5);
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
  
  
  TCanvas* can2 = new TCanvas("can2");
  can2->Divide(2,2);
  can2->cd(1);
  gPad->SetLogy();
  efrac->Draw();
  
  can2->cd(2);
  gPad->SetLogy();
//  efracc->Draw();
  efracp->Draw();
  
  can2->cd(3);
  gPad->SetLogy();
  ngh->Draw();
  
  can2->cd(4);
  gPad->SetLogy();
  efracc2->Draw();
  
}
  

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



void testAmbiguity(TCanvas* can,double tres) {
  TF1* fwrong = new TF1("wrong","1.0-TMath::Freq(x/[0])",0.0,0.25);
  char title[100];
  snprintf(title,100,"Ambiguity error, interpolation resolution=%f cm",tres);
  TH1F* hwrong = new TH1F("hwrong",title,200,0,0.25);
  hwrong->GetXaxis()->SetTitle("distance to wire (cm)");
  hwrong->GetYaxis()->SetTitle("Wrong Ambiguity Fraction");
  hwrong->SetMaximum(0.5);
  hwrong->SetStats(0);
  
  TF1* frms = new TF1("frms","2*x*sqrt((1.0-TMath::Freq(x/[0]))*TMath::Freq(x/[0]))",0,0.25);
  TH1F* hrms = new TH1F("hrms","RMS due to Ambiguity error",100,0,0.25);
  hrms->GetXaxis()->SetTitle("distance to wire (cm)");
  hrms->GetYaxis()->SetTitle("resolution RMS");
  hrms->SetMaximum(0.02);
  hrms->SetStats(0);
  
  can->Clear();
  can->Divide(1,2);
  can->cd(1);
  hwrong->Draw();
  fwrong->SetParameter(0,tres);
  fwrong->Draw("same");
  
  can->cd(2);
  hrms->Draw();
  frms->SetParameter(0,tres);
  frms->Draw("same");
  
}

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


void testDIO(TCanvas* can,const char* file="~/Downloads/dio_spectrum.csv") {
  TGraph* spectrum = new TGraph(file, "%lg,%lg");
  TSpline3* spspect = new TSpline3("spspect",spectrum);
  Int_t npt = 1000;
  double range=105;
  double dx = range/npt;
  std::vector<Double_t> x;
  std::vector<Double_t> y;
  double xval = dx/2.0;
  double yval(0.0);  
  for(Int_t ipt=0;ipt<npt;ipt++){
    x.push_back(xval);
    y.push_back(yval);
    xval += dx;
    yval += spspect->Eval(xval);
  }
// normalize
  for(Int_t ipt=0;ipt<npt;ipt++){
    y[ipt] /= yval;
  }
// integral spectrum
  TGraph* intspect = new TGraph(x.size(),&x.front(),&y.front());
// inverse spectrum
  TGraph* invintspect = new TGraph(x.size(),&y.front(),&x.front());  
  
  TH1D* dio = new TH1D("dio","dio spectrum",200,0,106);
  dio->SetMaximum(4e-2);
  dio->SetMinimum(1e-6);
  dio->GetXaxis()->SetTitle("MeV");
  
  TH1D* intdio = new TH1D("intdio","integrated dio spectrum",200,0,106);
  intdio->SetMaximum(1.01);
  intdio->GetXaxis()->SetTitle("MeV");
  
  TH1D* invdio = new TH1D("invdio","inverse integral dio spectrum",201,0,1.01);
  invdio->SetMaximum(105);
  invdio->GetYaxis()->SetTitle("MeV");
  
  TH1D* sdio = new TH1D("dio","sampled dio spectrum",200,0,106);
  sdio->GetXaxis()->SetTitle("MeV");
  sdio->SetMinimum(1);
// sample the spectrum
  int ntrials=1000000;
  TRandom3* tr = new TRandom3();
  
  for(int itry=0;itry<ntrials;itry++){
    double mom = invintspect->Eval(tr->Rndm());
    sdio->Fill(mom);
  }
  
  
  can->Clear();
  can->Divide(2,2);
  can->cd(1);
  dio->Draw();
//  spectrum->Draw("same");
  spspect->Draw("same");
  can->cd(2);
  intdio->Draw();
  
  intspect->Draw("same");
  cout << "Integral 0-40 = " << intspect->Eval(40.0) << endl;
  
  can->cd(3);
  invdio->Draw();
  invintspect->Draw("same");
  
  can->cd(4);
  sdio->Draw();
}

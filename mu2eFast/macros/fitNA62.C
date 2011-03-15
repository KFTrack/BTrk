#include <Rtypes.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TStyle.h>

// fit the hit resolution function from NA62

void FitNA62(TCanvas* can,const char* resfile) {
  TGraph* star = new TGraph(resfile);
  TH1F* resfun= new TH1F("resfun","NA62 straw hit resolution, NIMA 604(2009) 307-309",100,0,2);
  resfun->SetMaximum(0.025);
  star->SetMarkerStyle(3);
  resfun->GetXaxis()->SetTitle("radius/2.5mm");
  resfun->GetYaxis()->SetTitle("average resolution (cm)");
  resfun->SetStats(0);
  can->Clear();
  resfun->Draw();
  star->Fit("pol3");
  star->Draw("same,*");
}
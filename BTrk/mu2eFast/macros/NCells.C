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

void NCells(double maxpath=5.0,double threshold=0.2,unsigned ntrials=100000) {
  TProfile* mync = new TProfile("mync","<N cells> vs relative pathlength, my formula",200,0,maxpath,0,100);
  mync->SetLineColor(kBlue);
  TProfile* gfnc = new TProfile("gfnc","<N cells> vs relative pathlength, Gianfranco",200,0,maxpath,0,100);
  gfnc->SetLineColor(kRed);
  TRandom3* tr = new TRandom3();
  for(unsigned itry=0;itry<ntrials;itry++){
    double cpath = tr->Rndm()*maxpath;
    unsigned myncells = max((unsigned)1,(unsigned)ceil(cpath-threshold));
    double mycfrac = cpath-myncells+1;
    if(mycfrac - threshold > tr->Rndm())myncells++;
    mync->Fill(cpath,myncells);
    
    
    unsigned gfncells = (unsigned)floor(cpath);
    double gfcfrac = cpath-gfncells;
    gfncells++;
    if(gfcfrac - threshold > tr->Rndm())gfncells++;
    gfnc->Fill(cpath,gfncells);
  }  
  TCanvas* can = new TCanvas("can");
  can->Clear();
  can->Divide(1,2);
  can->cd(1);
  mync->Draw();
  can->cd(2);
  gfnc->Draw();
  can->SaveAs("ncells.png");
}
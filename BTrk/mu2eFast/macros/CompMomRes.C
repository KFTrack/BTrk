#include "../include/PacSimHitInfo.hh"
#include "../include/PacSimTrkSummary.hh"
#ifdef __MAKECINT__
#pragma link C++ class PacSimHitInfo;
#pragma link C++ class PacSimTrkSummary;
#pragma link C++ class vector<PacSimHitInfo>;
#endif


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
#include <TProfile.h>
#include <TObjArray.h>
#include <TFile.h>
#include <stdio.h>
#include "Riostream.h"

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

void CompMomRes(TCanvas* can, std::vector<std::string>& files) {
  int colors[6] = {kRed,kBlue,kGreen,kMagenta,kCyan,kBlack};
  std::vector<TFile*> tfiles;
  std::vector<TTree*> trees;
  std::vector<TH1F*> fitprob;
  std::vector<TH1F*> momres;
  TCut rec("rec_ndof>15");
  TCut goodradius("abs(rec_d0)<10.0 && abs(2.0/rec_omega - rec_d0)<68.0");
  TCut gooddip("rec_tandip>0.5773&&rec_tandip<1.0");
  TCut goodfit("rec_fitprob>0.05&&rec_mom_err<0.0005");
  TCut goodrec = rec+goodradius+gooddip+goodfit;
  
  for(unsigned ifile=0;ifile<files.size();ifile++){
    TFile* tfile = new TFile(files[ifile].c_str());
    if(tfile != 0){
      tfiles.push_back(tfile);
      TTree* tree = (TTree*)tfile->Get("tracks");
      if(tree != 0){
        trees.push_back(tree);
        char name[20];

        snprintf(name,20,"fitprob%i",ifile);
        TH1F* fprob = new TH1F(name,"Kalman fit consistency",200,0.0,1.0);
        fprob->SetLineColor(colors[ifile]);
        fprob->SetStats(0);
        fprob->SetMinimum(1);
        fprob->SetMaximum(tree->GetEntriesFast());
        
        fitprob.push_back(fprob);
        tree->Project(name,"rec_fitprob",rec+gooddip+goodradius);
        
        snprintf(name,20,"momres%i+",ifile);
        TH1F* momr = new TH1F(name,"momentum resolution",200,-2,2);
        momr->SetLineColor(colors[ifile]);
        momr->SetStats(0);
        momr->GetXaxis()->SetTitle("MeV");
        
        momres.push_back(momr);
        tree->Project(name,"1000*(rec_mom_mag-sim_mom_mag)",goodrec);
        
      } else {
        std::cerr << "Can find tree 'tracks' in file " << files[ifile] << std::endl;
        return;
      }
    } else {
      std::cerr << "Can't find file " << files[ifile] << std::endl;
      return;
    }
  }
  can->Divide(1,2);
  can->cd(1);
  gPad->SetLogy(1);
  TLegend* leg = new TLegend(0.6,0.6,0.9,0.9);
  for(unsigned ifile=0;ifile<files.size();ifile++){
    fitprob[ifile]->SetLineColor(colors[ifile]);
    if(ifile > 0)
      fitprob[ifile]->Draw("same");
    else
      fitprob[ifile]->Draw();

    TString fullname(files[ifile].c_str());
    TObjArray* sstrings = fullname.Tokenize(TString("/."));
    
    TString fname = ((TObjString*)(*sstrings)[sstrings->GetEntries()-2])->String();
    double eff = momres[ifile]->GetEntries()/trees[ifile]->GetEntriesFast();
    cout << "Selection efficiency for file " << files[ifile] << " = " << eff << endl;
    char effs[20];
    snprintf(effs,20," efficiency = %g",eff);
    fname += TString(effs);
    leg->AddEntry(fitprob[ifile],fname,"L");
  }
  leg->Draw();
  
  can->cd(2);
  gPad->SetLogy(1);
  for(unsigned ifile=0;ifile<files.size();ifile++){
    if(ifile > 0)
      momres[ifile]->Draw("same");
    else
      momres[ifile]->Draw();
  }
// compute efficiency of final cuts
}

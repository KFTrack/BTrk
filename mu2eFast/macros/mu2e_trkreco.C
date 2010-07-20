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


void mu2e_trkreco(TCanvas* can, TTree* tree, const char* cpage="rec" ) {
  TString page(cpage);
  TCut rec("rec_ndof>0");
  TCut goodfit("rec_fitprob>0.05&&rec_ndof>=10&&rec_mom_err<0.0005&&abs(rec_d0)<10.0");
  
  TF1* sgau = new TF1("sgau",splitgaus,-1.,1.,7);
  sgau->SetParName(0,"Norm");
  sgau->SetParName(1,"Mean");
  sgau->SetParName(2,"SigH");
  sgau->SetParName(3,"SigL");
  sgau->SetParName(4,"TFH");
  sgau->SetParName(5,"TSigH");
  sgau->SetParName(6,"TSigL");
  TF1* dgau = new TF1("dgau",doublegaus,-1.,1.,5);
  if( page == "sim"){
    TH1F* mom = new TH1F("mom","momentum",100,0.09,0.11);
    TH1F* td = new TH1F("td","tandip",100,0,1.4);
    TH1F* z0 = new TH1F("z0","z position",100,-330,-240);
    TH1F* d0 = new TH1F("d0","transverse position",100,-25,25);
    TH2F* nvtd = new TH2F("nvtd","N planes vs tandip",51,0.5,1.1,51,-0.5,50.5);
    TH2F* nvd0 = new TH2F("nvd0","N planes vs transverse position",51,-25.5,25.5,51,-0.5,50.5);
    

    tree->Project("mom","sim_mom_mag");
    tree->Project("td","sim_tandip");
    tree->Project("z0","sim_z0");
    tree->Project("d0","sim_d0");
    tree->Project("nvtd","simtrk.ngas:sim_tandip");
    tree->Project("nvd0","simtrk.ngas:sim_d0");
    can->Clear();
    can->Divide(3,2);
    can->cd(1);
    mom->Draw();
    can->cd(2);
    td->Draw();
    can->cd(3);
    z0->Draw();
    can->cd(4);
    d0->Draw();
    can->cd(5);
    nvtd->Draw("box");
    can->cd(6);
    nvd0->Draw("box");
    
  } else if(page == "rec"){

    TH1F* ndof = new TH1F("ndof","N DOF",40,-0.5,39.5);
    tree->Project("ndof","rec_ndof");
    
    TH1F* nhit = new TH1F("nhit","N hits",50,-0.5,49.5);
    tree->Project("nhit","rec_nhit");

    TH1F* chindof = new TH1F("chindof","Chisquare/NDOF", 100,0.0, 10.0);
    tree->Project("chindof","rec_chisqr/rec_ndof",rec);
    
    TH1F* fitprob = new TH1F("fitprob","fit consistency",200,0.0,1.0);
    tree->Project("fitprob","rec_fitprob",rec);
    
    can->Clear();
    can->Divide(2,2);
    can->cd(1);
    ndof->Draw();
    can->cd(2);
    nhit->Draw();
    can->cd(3);
    chindof->Draw();
    can->cd(4);
    gPad->SetLogy();
    fitprob->Fit("pol1","","",0.05,1);
    
  } else if(page == "eff"){

    TH1F* td_s = new TH1F("td_s","TanDip",100,0,1.4);
    TH1F* td_r = new TH1F("td_r","TanDip",100,0,1.4);
    TH1F* td_g = new TH1F("td_g","TanDip",100,0,1.4);
    tree->Project("td_s","sim_tandip");
    tree->Project("td_r","sim_tandip",rec);
    tree->Project("td_g","sim_tandip",rec+goodfit);
    td_r->Divide(td_s);
    td_g->Divide(td_s);
    td_r->SetLineColor(kRed);
    td_g->SetLineColor(kBlue);
    
    TH1F* z0_s = new TH1F("z0_s","z0",100,-350,-230);
    TH1F* z0_r = new TH1F("z0_r","z0",100,-350,-230);
    TH1F* z0_g = new TH1F("z0_g","z0",100,-350,-230);
    tree->Project("z0_s","sim_z0");
    tree->Project("z0_r","sim_z0",rec);
    tree->Project("z0_g","sim_z0",rec+goodfit);
    z0_r->Divide(z0_s);
    z0_g->Divide(z0_s);
    z0_r->SetLineColor(kRed);
    z0_g->SetLineColor(kBlue);
    
    
    TH2F* pvtd_s = new TH2F("pvtd_s","Phi vs TanDip",50,0,1.4,50,-3.15,3.15);
    TH2F* pvtd_g = new TH2F("pvtd_g","Phi vs TanDip",50,0,1.4,50,-3.15,3.15);
    tree->Project("pvtd_s","sim_phi0:sim_tandip");
    tree->Project("pvtd_g","sim_phi0:sim_tandip",rec+goodfit);
    pvtd_g->Divide(pvtd_s);
    pvtd_g->SetLineColor(kBlue);
    
    
    
    can->Clear();
    can->Divide(2,2);

    can->cd(1);
    td_r->Draw();
    td_g->Draw("same");
    TLegend* leg = new TLegend(0.1,0.7,0.4,0.9);
    leg->AddEntry(td_r,"All Fits","L");
    leg->AddEntry(td_g,"Fit con>0.01","L");
    leg->Draw();

    can->cd(2);
    z0_r->Draw();
    z0_g->Draw("same");

    can->cd(4);
    pvtd_g->Draw("box");
    
  } else if (page == "pull"){
    gStyle->SetOptFit(1111);
    TH1F* d0p = new TH1F("d0p","d0 pull",100,-10,10);
    TH1F* p0p = new TH1F("p0p","phi0 pull",100,-10,10);
    TH1F* omp = new TH1F("omp","omega pull",100,-10,10);
    TH1F* z0p = new TH1F("z0p","z0 pull",100,-10,10);
    TH1F* tdp = new TH1F("tdp","tandip pull",100,-10,10);
    TH1F* momp = new TH1F("momp","momentum pull",100,-10,10);
    d0p->SetLineColor(kRed);
    p0p->SetLineColor(kRed);
    omp->SetLineColor(kRed);
    z0p->SetLineColor(kRed);
    tdp->SetLineColor(kRed);
    momp->SetLineColor(kRed);
    
    TH1F* d0pg = new TH1F("d0pg","d0 pull",100,-10,10);
    TH1F* p0pg = new TH1F("p0pg","phi0 pull",100,-10,10);
    TH1F* ompg = new TH1F("ompg","omega pull",100,-10,10);
    TH1F* z0pg = new TH1F("z0pg","z0 pull",100,-10,10);
    TH1F* tdpg = new TH1F("tdpg","tandip pull",100,-10,10);
    TH1F* mompg = new TH1F("mompg","momentum pull",100,-10,10);
    d0pg->SetLineColor(kBlue);
    p0pg->SetLineColor(kBlue);
    ompg->SetLineColor(kBlue);
    z0pg->SetLineColor(kBlue);
    tdpg->SetLineColor(kBlue);
    mompg->SetLineColor(kBlue);
    
    tree->Project("d0p","pull_d0",rec);
    tree->Project("p0p","pull_phi0",rec);
    tree->Project("omp","pull_omega",rec);
    tree->Project("z0p","pull_z0",rec);
    tree->Project("tdp","pull_tandip",rec);
    tree->Project("momp","(rec_mom_mag-sim_mom_mag)/rec_mom_err",rec);

    
    tree->Project("d0pg","pull_d0",rec+goodfit);
    tree->Project("p0pg","pull_phi0",rec+goodfit);
    tree->Project("ompg","pull_omega",rec+goodfit);
    tree->Project("z0pg","pull_z0",rec+goodfit);
    tree->Project("tdpg","pull_tandip",rec+goodfit);
    tree->Project("mompg","(rec_mom_mag-sim_mom_mag)/rec_mom_err",rec+goodfit);
    
    can->Clear();
    can->Divide(3,2);
    can->cd(1);
    gPad->SetLogy();
    d0p->Draw();
    d0pg->Fit("gaus","","sames");
    TLegend* leg = new TLegend(0.1,0.7,0.4,0.9);
    leg->AddEntry(d0p,"All Fits","L");
    leg->AddEntry(d0pg,"Fit con>0.01","L");
    leg->Draw();
    
    can->cd(2);
    gPad->SetLogy();
    p0p->Draw();
    p0pg->Fit("gaus","","sames");
    can->cd(3);
    gPad->SetLogy();
    omp->Draw();
    ompg->Fit("gaus","","sames");
    can->cd(4);
    gPad->SetLogy();
    z0p->Draw();
    z0pg->Fit("gaus","","sames");
    can->cd(5);
    gPad->SetLogy();
    tdp->Draw();
    tdpg->Fit("gaus","","sames");
    can->cd(6);
    gPad->SetLogy();
    momp->Draw();
    mompg->Fit("gaus","","sames");
        
  } else if (page=="res"){
    gStyle->SetOptFit(1111);
    TH1F* d0r = new TH1F("d0p","d0 resolution",100,-5,5);
    TH1F* p0r = new TH1F("p0p","phi0 resolution",100,-0.1,0.1);
    TH1F* omr = new TH1F("omp","omega resolution",100,-0.002,0.002);
    TH1F* z0r = new TH1F("z0p","z0 resolution",100,-10,10);
    TH1F* tdr = new TH1F("tdp","tandip resolution",100,-0.1,0.1);
    TH1F* momr = new TH1F("momp","momentum resolution",200,-0.0025,0.0025);
    
    tree->Project("d0p","rec_d0-sim_d0",rec+goodfit);
    tree->Project("p0p","rec_phi0-sim_phi0",rec+goodfit);
    tree->Project("omp","rec_omega-sim_omega",rec+goodfit);
    tree->Project("z0p","rec_z0-sim_z0",rec+goodfit);
    tree->Project("tdp","rec_tandip-sim_tandip",rec+goodfit);
    tree->Project("momp","rec_mom_mag-sim_mom_mag",rec+goodfit);

    
    can->Clear();
    can->Divide(3,2);
    can->cd(1);
    gPad->SetLogy();
    d0r->Fit("gaus");
    can->cd(2);
    gPad->SetLogy();
    p0r->Fit("gaus");
    can->cd(3);
    gPad->SetLogy();
    omr->Fit("gaus");
    can->cd(4);
    gPad->SetLogy();
    z0r->Fit("gaus");
    can->cd(5);
    gPad->SetLogy();
    tdr->Fit("gaus");
    can->cd(6);
    gPad->SetLogy();
    double integral = momr->GetEntries()*momr->GetBinWidth(1);
    sgau->SetParameters(integral,0.0,momr->GetRMS(),momr->GetRMS(),0.01,2*momr->GetRMS(),2*momr->GetRMS());
    sgau->SetParLimits(5,1.0*momr->GetRMS(),1.0);
    sgau->SetParLimits(6,1.0*momr->GetRMS(),1.0);
    sgau->SetParLimits(4,0.0,0.1);
    momr->Fit("sgau","L");
//    momr->Fit("sgau","M");

  } else if (page == "mom"){
    gStyle->SetOptFit(1111);
    TH1F* nhit = new TH1F("nhit","N hits",50,-0.5,49.5);
    tree->Project("nhit","rec_nhit",rec+goodfit);
    
    TH1F* mome = new TH1F("mome","estimated fit mom error",100,0.00001,0.0006);
    tree->Project("mome","rec_mom_err",rec+goodfit);
    
    TH1F* mompg = new TH1F("mompg","momentum pull",100,-10,10);
    tree->Project("mompg","(rec_mom_mag-sim_mom_mag)/rec_mom_err",rec+goodfit);
    
    TH1F* momr = new TH1F("momp","momentum resolution",200,-0.002,0.002);
    tree->Project("momp","rec_mom_mag-sim_mom_mag",rec+goodfit);
    
    can->Clear();
    can->Divide(2,2);
    can->cd(1);
    nhit->Draw();
    can->cd(2);
    mome->Draw();
    can->cd(3);
    gPad->SetLogy();
    mompg->Fit("gaus");
    can->cd(4);
    gPad->SetLogy();
    double integral = momr->GetEntries()*momr->GetBinWidth(1);
    sgau->SetParameters(integral,0.0,momr->GetRMS(),momr->GetRMS(),0.01,2*momr->GetRMS(),2*momr->GetRMS());
    sgau->SetParLimits(5,1.0*momr->GetRMS(),1.0);
    sgau->SetParLimits(6,1.0*momr->GetRMS(),1.0);
    sgau->SetParLimits(4,0.0,0.1);
    momr->Fit("sgau","L");
//    momr->Fit("sgau","M");
    
  } else if(page == "mat") {
    gStyle->SetOptFit(1111);
    
    TH1F* dmom = new TH1F("dmom","momentum mag change",200,-0.00015,0.00005);
    TH1F* dang = new TH1F("dang","momentum dir change",200,-0.05,0.05);
    tree->Project("dmom","simtrk.dmom");
    tree->Project("dang","simtrk.ddir*cos(simtrk.ddirphi)");
    tree->Project("dang","simtrk.ddir*sin(simtrk.ddirphi)");
    can->Clear();
    can->Divide(2,1);
    can->cd(1);
    gPad->SetLogy();
    dmom->Draw();
    can->cd(2);
    gPad->SetLogy();
    double integral = dang->GetEntries()*dang->GetBinWidth(1);
    dgau->SetParameters(integral,0.0,dang->GetRMS()*5.0,dang->GetRMS()/2.0,0.05);
    dang->Fit("dgau");
  } else if(page == "trajdiff") {
    TCut adjacent = ("trajdiff.endglen-trajdiff.startglen>20 && trajdiff.endglen-trajdiff.startglen<40");
    TCut middle = ("trajdiff.endglen>700 && trajdiff.endglen < 800 && trajdiff.startglen> 550 && trajdiff.startglen< 650");
    TCut early = ("trajdiff.startglen<600");
    
    TH2F* glen = new TH2F("glen","end vs start flightlen",100,400,1300,100,400,1300);
    TH1F* dg = new TH1F("dg","global length difference",100,0,100);
    TH1F* adiff = new TH1F("adiff","average traj diff, adjacent stations",100,0.0,0.5);
    TH1F* mdiff = new TH1F("mdiff","average traj diff, middle of tracker",100,0.0,0.5);
    tree->Project("glen","trajdiff.endglen:trajdiff.startglen");
    tree->Project("dg","trajdiff.endglen-trajdiff.startglen");
    tree->Project("adiff","trajdiff.ddiff",adjacent+rec+goodfit+early);
    tree->Project("mdiff","trajdiff.ddiff",middle+rec+goodfit+early);
    can->Clear();
    can->Divide(2,2);
    can->cd(1);
    glen->Draw();
//    tree->Draw("trajdiff.endglen:trajdiff.startglen>>glen");
    can->cd(2);
    dg->Draw();
    can->cd(3);
    gPad->SetLogy();
    adiff->Draw();
    can->cd(4);
    gPad->SetLogy();
    mdiff->Draw();
  } else if(page == "bintdiff") {
    TCut adjacent = ("trajdiff.endglen-trajdiff.startglen>20 && trajdiff.endglen-trajdiff.startglen<40");
    TCut middle = ("trajdiff.endglen>700 && trajdiff.endglen < 800 && trajdiff.startglen> 550 && trajdiff.startglen< 650");
    TCut early = ("trajdiff.startglen<600");

    TH1F* tadiff = new TH1F("tadiff","True Delta P, adjacent stations",100,0,2e-3);
    TH1F* tmdiff = new TH1F("tmdiff","True Delta P, middle of tracker",100,0,2e-3);
    TH1F* radiff = new TH1F("radiff","Reco Delta P, adjacent stations",100,0,2e-3);
    TH1F* rmdiff = new TH1F("rmdiff","Reco Delta P, middle of tracker",100,0,2e-3);
    TH1F* dadiff = new TH1F("dadiff","Delta Delta P, adjacent stations",100,0,1e-4);
    TH1F* dmdiff = new TH1F("dmdiff","Delta Delta P, middle of tracker",100,0,1e-4);
    
    tree->Project("tadiff","trajdiff.truedp",adjacent+rec+goodfit+early);
    tree->Project("tmdiff","trajdiff.truedp",middle+rec+goodfit+early);
    tree->Project("radiff","trajdiff.recodp",adjacent+rec+goodfit+early);
    tree->Project("rmdiff","trajdiff.recodp",middle+rec+goodfit+early);
    tree->Project("dadiff","trajdiff.deltadp",adjacent+rec+goodfit+early);
    tree->Project("dmdiff","trajdiff.deltadp",middle+rec+goodfit+early);
    can->Clear();
    can->Divide(2,3);
    can->cd(1);
    tadiff->Draw();
    can->cd(2);
    tmdiff->Draw();
    can->cd(3);
    radiff->Draw();
    can->cd(4);
    rmdiff->Draw();    
    can->cd(5);
    gPad->SetLogy();
    dadiff->Draw();
    can->cd(6);
    gPad->SetLogy();
    dmdiff->Draw();    
  } else if(page == "caloresid") {
    gStyle->SetOptFit(1111);
    TCut calor("simhit.shelemnum<100&&simhit.shtypenum>600&&simhit.shnhot>0");
    TCut trans("abs(simhit.shx)<1.0");
    TH1F* zpos = new TH1F("zpos","Z position",100,310,480);
    TH1F* rpos = new TH1F("rpos","R position",100,35,80);
    TH1F* zres = new TH1F("zres","Z calo resid",100,-3,3);
    TH1F* rres = new TH1F("rres","R calo resid",100,-1,1);
    tree->Project("zpos","simhit.shz",calor);
    tree->Project("rpos","sqrt(simhit.shx^2+simhit.shy^2)",calor+trans);
    tree->Project("zres","simhit.shdz",calor);
    tree->Project("rres","simhit.shdy",calor+trans);
    can->Clear();
    can->Divide(2,2);
    can->cd(1);
    zpos->Draw();
    can->cd(2);
    rpos->Draw();
    can->cd(3);
    zres->Fit("gaus");
    can->cd(4);
    rres->Fit("gaus");
  }
}

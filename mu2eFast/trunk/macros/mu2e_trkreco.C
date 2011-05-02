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
#include <TMath.h>
#include <Math.h>
#include <TH2F.h>
#include <TTree.h>
#include <TCut.h>
#include <TLegend.h>
#include <TString.h>
#include <TStyle.h>
#include <TProfile.h>
#include <TLine.h>
#include <TRandom3.h>
#include <stdio.h>
#include "Riostream.h"
#include "Math/QuantFuncMathCore.h"

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

void mu2e_trkreco(TCanvas* can,TTree* tree, const char* cpage="rec" ) {
  TString page(cpage);
  TCut rec("rec_ndof>0");
  TCut goodradius("abs(rec_d0)<10.0 && abs(2.0/rec_omega + rec_d0)<65.0");
  TCut gooddip("rec_tandip>0.5773&&rec_tandip<1.0");
  TCut goodfit("rec_fitprob>0.01&&rec_ndof>=20");
  TCut goodres("rec_mom_err<0.001");
  TCut goodhits("rec_nhit-rec_nactive<10");
  TCut gen("sim_inipos_z<-300 && sim_pdgid==11");
  TCut goodrec = goodhits+goodradius+gooddip+goodfit+goodres;
  TCut goodfitp("rec_fitprob>0.01");
  TCut goodndof("rec_ndof>=20");
  TCut goodmerr("rec_mom_err<0.001");
  TCut goodd0("abs(rec_d0)<10.0");
  TCut goodrmax("abs(2.0/rec_omega + rec_d0)<68.0");
  TCut signalbox("rec_mom_mag>0.104");
  
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
    TH1F* mom = new TH1F("mom","momentum",100,0.0,0.11);
    mom->GetXaxis()->SetTitle("GeV");
    TH1F* cost = new TH1F("cost","Cos(#theta)",100,-1,1);
    cost->SetMinimum(0);
    TH1F* z0 = new TH1F("z0","z position",91,-475.5,-384.5);
    z0->GetXaxis()->SetTitle("cm from tracker center");
    z0->SetStats(0);
    TH1F* rprof = new TH1F("rprof","radial profile",100,-11,11);
    rprof->GetXaxis()->SetTitle("cm");
//    rprof->SetStats(0);
    TH2F* xyprof = new TH2F("xyprof","xy profile",50,-11,11,50,-11,11);
    xyprof->GetXaxis()->SetTitle("cm");
    xyprof->GetYaxis()->SetTitle("cm");
    xyprof->SetStats(0);
//    TH1F* d0 = new TH1F("d0","transverse position",100,-25,25);
//    TH2F* nvtd = new TH2F("nvtd","N planes vs tandip",51,0.5,1.1,51,-0.5,50.5);
//    TH2F* nvd0 = new TH2F("nvd0","N planes vs transverse position",51,-25.5,25.5,51,-0.5,50.5);
    

    tree->Project("mom","sim_mom_mag",gen);
    tree->Project("cost","sim_mom_cost",gen);
    tree->Project("z0","sim_inipos_z",gen);
    tree->Project("xyprof","sim_inipos_y:sim_inipos_x",gen);
    tree->Project("rprof","sim_inipos_x",gen);
    tree->Project("rprof+","sim_inipos_y",gen);
    
//    tree->Project("d0","sim_d0");
//    tree->Project("nvtd","simtrk.ngas:sim_tandip");
//    tree->Project("nvd0","simtrk.ngas:sim_d0");
    can->Clear();
    can->Divide(2,2);
    can->cd(1);
    mom->Draw();
    can->cd(2);
    cost->Draw();
    can->cd(3);
    z0->Draw();
    can->cd(4);
//    xyprof->Draw("box");
    rprof->Draw();
//    can->cd(5);
//    nvtd->Draw("box");
//    can->cd(6);
//    nvd0->Draw("box");
    
  } else if(page == "rec"){

    TH1F* ndof = new TH1F("ndof","N DOF",80,-0.5,79.5);
    tree->Project("ndof","rec_ndof");
    
    TH1F* nhit = new TH1F("nhit","N hits",80,-0.5,79.5);
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
  } else if(page == "selection1")  {
    TH1F* d0 = new TH1F("d0","DOCA to Z axis",100,0,20);
    TH1F* sd0 = new TH1F("sd0","DOCA to Z axis",100,0,20);
    TH1F* gd0 = new TH1F("gd0","DOCA to Z axis",100,0,20);
    d0->SetLineColor(kBlue);
    sd0->SetLineColor(kGreen);
    gd0->SetLineColor(kRed);
    d0->GetXaxis()->SetTitle("cm");
    
    TH1F* rmax = new TH1F("rmax","Maximum radius",100,40,70);
    TH1F* srmax = new TH1F("srmax","Maximum radius",100,40,70);
    TH1F* grmax = new TH1F("grmax","Maximum radius",100,40,70);
    rmax->SetLineColor(kBlue);
    srmax->SetLineColor(kGreen);
    grmax->SetLineColor(kRed);
    rmax->GetXaxis()->SetTitle("cm");

    TH1F* td = new TH1F("td","tanDip",100,0.25,1.5);
    TH1F* std = new TH1F("std","tanDip",100,0.25,1.5);
    TH1F* gtd = new TH1F("gtd","tanDip",100,0.25,1.5);
    td->SetLineColor(kBlue);
    std->SetLineColor(kGreen);
    gtd->SetLineColor(kRed);

    TH1F* ndof = new TH1F("ndof","Fit N DOF",100,-0.5,99.5);
    TH1F* sndof = new TH1F("sndof","Fit N DOF",100,-0.5,99.5);
    TH1F* gndof = new TH1F("gndof","Fit N DOF",100,-0.5,99.5);
    ndof->SetLineColor(kBlue);
    sndof->SetLineColor(kGreen);
    gndof->SetLineColor(kRed);
    
    
    tree->Project("d0","abs(rec_d0)",rec);
    tree->Project("sd0","abs(rec_d0)",rec+gen);
    tree->Project("gd0","abs(rec_d0)",rec+goodrmax+gooddip+goodfit+goodres+goodhits);
    
    tree->Project("rmax","abs(2.0/rec_omega - rec_d0)",rec);
    tree->Project("srmax","abs(2.0/rec_omega - rec_d0)",rec+gen);
    tree->Project("grmax","abs(2.0/rec_omega - rec_d0)",rec+goodd0+gooddip+goodfit+goodres+goodhits);
    
    tree->Project("td","rec_tandip",rec);
    tree->Project("std","rec_tandip",rec+gen);
    tree->Project("gtd","rec_tandip",rec+goodradius+goodfit+goodres+goodhits);

    tree->Project("ndof","rec_ndof",rec);
    tree->Project("sndof","rec_ndof",rec+gen);
    tree->Project("gndof","rec_ndof",rec+goodradius+gooddip+goodfitp+goodres+goodmerr+goodhits);
    
    can->Clear();
    can->Divide(2,2);
    can->cd(1);
//    gPad->SetLogy();
    d0->SetMinimum(1);
    d0->Draw();
    sd0->Draw("same");
    gd0->Draw("same");
    TLine* d0cut = new TLine(10.0,0.0,10.0,0.5*d0->GetMaximum());
    d0cut->Draw("same");
    can->cd(2);
    rmax->Draw();
    srmax->Draw("same");
    grmax->Draw("same");
    TLine* rmaxcut = new TLine(65.0,0.0,65.0,0.5*rmax->GetMaximum());
    rmaxcut->Draw("same");
    can->cd(3);
    td->Draw();
    std->Draw("same");
    gtd->Draw("same");
    TLine* tdcut1 = new TLine(0.5774,0.0,0.5774,0.5*td->GetMaximum());
    TLine* tdcut2 = new TLine(1.0,0.0,1.0,0.5*td->GetMaximum());
    tdcut1->Draw("same");
    tdcut2->Draw("same");
    
    can->cd(4);
//    gPad->SetLogy();
    ndof->SetMinimum(1);
    ndof->Draw();
    sndof->Draw("same");
    gndof->Draw("same");
    TLine* ndofcut = new TLine(20,0.0,20.0,0.5*ndof->GetMaximum());
    ndofcut->Draw("same");
    TLegend* leg = new TLegend(0.5,0.7,0.9,0.9);
//    leg->AddEntry(fitp,"Conversion + DIO","L");
    leg->AddEntry(sndof,"Conversion","L");
    leg->AddEntry(gndof,"All other cuts applied","L");
    leg->Draw();
    
    
  } else if(page == "selection2")  {
    
    
    TH1F* nmiss = new TH1F("nmiss","N missing hits",51,-0.5,50.5);
    TH1F* snmiss = new TH1F("snmiss","N missing hits",51,-0.5,50.5);
    TH1F* gnmiss = new TH1F("gnmiss","N missing hits",51,-0.5,50.5);
    nmiss->SetLineColor(kBlue);
    snmiss->SetLineColor(kGreen);
    gnmiss->SetLineColor(kRed);
    
    TH1F* fitp = new TH1F("fitp","Fit consistency",501,-0.001,1.001);
    TH1F* sfitp = new TH1F("sfitp","Fit consistency",501,-0.001,1.001);
    TH1F* gfitp = new TH1F("gfitp","Fit consistency",501,-0.001,1.001);
    fitp->SetLineColor(kBlue);
    sfitp->SetLineColor(kGreen);
    gfitp->SetLineColor(kRed);

    TH1F* merr = new TH1F("merr","Estimated mom error",100,0,2.5);
    TH1F* smerr = new TH1F("smerr","Estimated mom error",100,0,2.5);
    TH1F* gmerr = new TH1F("gmerr","Estimated mom error",100,0,2.5);
    merr->SetLineColor(kBlue);
    smerr->SetLineColor(kGreen);
    gmerr->SetLineColor(kRed);
    merr->GetXaxis()->SetTitle("MeV");
    
    TH1F* mom = new TH1F("mom","momentum",100,103,107);
    TH1F* smom = new TH1F("smom","momentum",100,103,107);
    TH1F* gmom = new TH1F("gmom","momentum",100,103,107);
    mom->SetLineColor(kBlue);
    smom->SetLineColor(kGreen);
    gmom->SetLineColor(kRed);
    mom->GetXaxis()->SetTitle("MeV");
    

    tree->Project("fitp","rec_fitprob",rec);
    tree->Project("sfitp","rec_fitprob",rec+gen);
    tree->Project("gfitp","rec_fitprob",rec+goodrmax+gooddip+goodndof+goodmerr+goodhits);

    tree->Project("merr","1000*rec_mom_err",rec);
    tree->Project("smerr","1000*rec_mom_err",rec+gen);
    tree->Project("gmerr","1000*rec_mom_err",rec+goodrmax+gooddip+goodndof+goodfitp+goodhits);

    tree->Project("nmiss","rec_nhit-rec_nactive",rec);
    tree->Project("snmiss","rec_nhit-rec_nactive",rec+gen);
    tree->Project("gnmiss","rec_nhit-rec_nactive",rec+goodrmax+gooddip+goodfit+goodres);

    tree->Project("mom","1000*rec_mom_mag",rec);
    tree->Project("smom","1000*rec_mom_mag",rec+gen);
    tree->Project("gmom","1000*rec_mom_mag",rec+goodrec);
    
    can->Clear();
    can->Divide(2,2);
    
    can->cd(1);
    fitp->Draw();
    sfitp->Draw("same");
    gfitp->Draw("same");
    TLine* fitpcut = new TLine(0.01,0.0,0.01,0.5*fitp->GetMaximum());
    fitpcut->Draw("same");
    
    
    TLegend* leg = new TLegend(0.5,0.7,0.9,0.9);
//    leg->AddEntry(fitp,"Conversion + DIO","L");
    leg->AddEntry(sfitp,"Conversion","L");
    leg->AddEntry(gfitp,"All other cuts applied","L");
    leg->Draw();
    
    
    can->cd(2);
    gPad->SetLogy();
    merr->SetMinimum(1);
    merr->Draw();
    smerr->Draw("same");
    gmerr->Draw("same");
    TLine* merrcut = new TLine(1,0.0,1,0.5*merr->GetMaximum());
    merrcut->Draw("same");
    
    can->cd(3);
    gPad->SetLogy();
    nmiss->SetMinimum(1);
    nmiss->Draw();
    snmiss->Draw("same");
    gnmiss->Draw("same");
    TLine* nmisscut = new TLine(10.0,0.0,10.0,0.5*nmiss->GetMaximum());
    nmisscut->Draw("same");
    
    can->cd(4);
    gPad->SetLogy();
    mom->SetMinimum(1);
    mom->Draw();
    smom->Draw("same");
    gmom->Draw("same");

    
  } else if(page == "eff"){
    TCut goodrec_noa = goodhits+goodradius+goodfit+goodres;

    TH1F* td_s = new TH1F("td_s","TanDip",100,0.0,2.0);
    TH1F* td_r = new TH1F("td_r","TanDip",100,0.0,2.0);
    TH1F* td_g = new TH1F("td_g","TanDip",100,0.0,2.0);
    tree->Project("td_s","sim_tandip");
    tree->Project("td_r","sim_tandip",rec);
    tree->Project("td_g","sim_tandip",goodrec_noa);
    td_r->Divide(td_s);
    td_g->Divide(td_s);
    td_r->SetLineColor(kRed);
    td_g->SetLineColor(kBlue);
    td_r->SetStats(0);
    
    TH1F* ct_s = new TH1F("ct_s","Cos(#theta)",100,0.0,1.0);
    TH1F* ct_r = new TH1F("ct_r","Cos(#theta)",100,0.0,1.0);
    TH1F* ct_g = new TH1F("ct_g","Cos(#theta)",100,0.0,1.0);
    tree->Project("ct_s","sim_mom_cost");
    tree->Project("ct_r","sim_mom_cost",rec);
    tree->Project("ct_g","sim_mom_cost",goodrec_noa);
    ct_r->Divide(ct_s);
    ct_g->Divide(ct_s);
    ct_r->SetLineColor(kRed);
    ct_g->SetLineColor(kBlue);
    ct_r->SetStats(0);
    
    
    TH1F* z0_s = new TH1F("z0_s","Production z",100,-480,-380);
    TH1F* z0_r = new TH1F("z0_r","Production z",100,-480,-380);
    TH1F* z0_g = new TH1F("z0_g","Production z",100,-480,-380);
    tree->Project("z0_s","sim_inipos_z");
    tree->Project("z0_r","sim_inipos_z",rec);
    tree->Project("z0_g","sim_inipos_z",goodrec);
    z0_r->Divide(z0_s);
    z0_g->Divide(z0_s);
    z0_r->SetLineColor(kRed);
    z0_g->SetLineColor(kBlue);
    z0_r->SetStats(0);
    TH1F* mom_s = new TH1F("mom_s","momentum",100,50,160);
    TH1F* mom_r = new TH1F("mom_r","momentum",100,50,160);
    TH1F* mom_g = new TH1F("mom_g","momentum",100,50,160);
    tree->Project("mom_s","1000*sim_mom_mag");
    tree->Project("mom_r","1000*sim_mom_mag",rec);
    tree->Project("mom_g","1000*sim_mom_mag",goodrec);
    mom_r->Divide(mom_s);
    mom_g->Divide(mom_s);
    mom_r->SetLineColor(kRed);
    mom_g->SetLineColor(kBlue);
    mom_r->GetXaxis()->SetTitle("MeV");
    mom_r->GetYaxis()->SetTitle("efficiency");
    mom_r->SetStats(0);
    
    
    TH1F* pt_s = new TH1F("pt_s","transverse momentum",100,40,150);
    TH1F* pt_r = new TH1F("pt_r","transverse momentum",100,40,150);
    TH1F* pt_g = new TH1F("pt_g","transverse momentum",100,40,150);
    tree->Project("pt_s","1000*sim_mom_pt");
    tree->Project("pt_r","1000*sim_mom_pt",rec);
    tree->Project("pt_g","1000*sim_mom_pt",goodrec);
    pt_r->Divide(pt_s);
    pt_g->Divide(pt_s);
    pt_r->SetLineColor(kRed);
    pt_g->SetLineColor(kBlue);
    pt_r->GetXaxis()->SetTitle("MeV");
    pt_r->GetYaxis()->SetTitle("efficiency");
    pt_r->SetStats(0);
    
    
    can->Clear();
    can->Divide(2,2);

    can->cd(1);
    td_r->Draw();
    td_g->Draw("same");
    TLine* tdcut1 = new TLine(0.5774,0.0,0.5774,0.9*td_r->GetMaximum());
    TLine* tdcut2 = new TLine(1.0,0.0,1.0,0.9*td_r->GetMaximum());
    tdcut1->Draw("same");
    tdcut2->Draw("same");
    

    can->cd(2);
//    z0_r->Draw();
//    z0_g->Draw("same");
    ct_r->Draw();
    ct_g->Draw("same");
    TLine* ctcut1 = new TLine(0.5,0.0,0.5,0.9*ct_r->GetMaximum());
    TLine* ctcut2 = new TLine(0.7071,0.0,0.7071,0.9*ct_r->GetMaximum());
    ctcut1->Draw("same");
    ctcut2->Draw("same");

    can->cd(3);
    mom_r->Draw();
    mom_g->Draw("same");
    TLegend* leg = new TLegend(0.1,0.7,0.4,0.9);
    leg->AddEntry(td_r,"All Fits","L");
    leg->AddEntry(td_g,"Good Fits","L");
    leg->Draw();
    
    can->cd(4);
    pt_r->Draw();
    pt_g->Draw("same");
    
  } else if(page =="count"){
    
    TH1F* nhit = new TH1F("nhit","N layer hits",81,-0.5,80.5);
    TH1F* nlay = new TH1F("nlay","N panel hits",45,-0.5,44.5);
    TH1F* nslay = new TH1F("nslay","N panel hits",45,-0.5,44.5);
    TH1F* ndlay = new TH1F("ndlay","N panel hits",45,-0.5,44.5);
    TH1F* nstation = new TH1F("nstation","N station hits",18,-0.5,17.5);
    TH1F* nsingle = new TH1F("nsingle","N station hits",18,-0.5,17.5);
    TH1F* ndouble = new TH1F("ndouble","N station hits",18,-0.5,17.5);
    TH1F* ntriple = new TH1F("ntriple","N plane hits",18,-0.5,17.5);
    TH2F* n3v2 = new TH2F("n3v2","N triple vs double",15,-0.5,14.5,15,-0.5,14.5);
    
    nhit->SetLineColor(kBlack);
    nlay->SetLineColor(kBlack);
    nslay->SetLineColor(kBlue);
    ndlay->SetLineColor(kRed);
    nstation->SetLineColor(kBlack);
    nsingle->SetLineColor(kGreen);
    ndouble->SetLineColor(kBlue);
    ntriple->SetLineColor(kRed);
    
    nlay->SetStats(0);
    nslay->SetStats(0);
    ndlay->SetStats(0);
    
    nstation->SetStats(0);
    nsingle->SetStats(0);
    ndouble->SetStats(0);
    ntriple->SetStats(0);
  
    nhit->SetLineWidth(2);
    nlay->SetLineWidth(2);
    nslay->SetLineWidth(2);
    ndlay->SetLineWidth(2);
    nstation->SetLineWidth(2);
    nsingle->SetLineWidth(2);
    ndouble->SetLineWidth(2);
    ntriple->SetLineWidth(2);
    n3v2->SetLineWidth(2);

    tree->Project("nhit","simtrk.nwiremeas",goodrec);
    tree->Project("nlay","simtrk.nwiremeas-sim_ndlayer",goodrec);
    tree->Project("nslay","simtrk.nwiremeas-2*sim_ndlayer",goodrec);
    tree->Project("ndlay","sim_ndlayer",goodrec);

    tree->Project("nstation","sim_nstation",goodrec);
    tree->Project("nsingle","sim_nsingle",goodrec);
    tree->Project("ndouble","sim_ndouble",goodrec);
    tree->Project("ntriple","sim_ntriple_ge",goodrec);

    tree->Project("n3v2","sim_ntriple_ge:sim_ndouble",goodrec);
    
    TLegend* legl = new TLegend(0.6,0.7,1.0,0.9);
    legl->AddEntry(nlay,"All","L");
    legl->AddEntry(nslay,"Single Layer","L");
    legl->AddEntry(ndlay,"Double Layer","L");
    
    TLegend* legs = new TLegend(0.6,0.6,1.0,0.9);
    legs->AddEntry(nstation,"All","L");
    legs->AddEntry(nsingle,"Single Plane","L");
    legs->AddEntry(ndouble,"Double Plane","L");
    legs->AddEntry(ntriple,">=Triple Plane","L");
    
    can->Clear();
    can->Divide(2,2);
    can->cd(1);
    nhit->Draw();
    can->cd(2);
    nslay->Draw();
    nlay->Draw("same");
    ndlay->Draw("same");
    legl->Draw();
    can->cd(3);
    nsingle->Draw();
    nstation->Draw("same");
    ndouble->Draw("same");
    ntriple->Draw("same");
    legs->Draw();
    can->cd(4);
    n3v2->Draw("box");
    
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

    
    tree->Project("d0pg","pull_d0",goodrec);
    tree->Project("p0pg","pull_phi0",goodrec);
    tree->Project("ompg","pull_omega",goodrec);
    tree->Project("z0pg","pull_z0",goodrec);
    tree->Project("tdpg","pull_tandip",goodrec);
    tree->Project("mompg","(rec_mom_mag-sim_mom_mag)/rec_mom_err",goodrec);
    
    can->Clear();
    can->Divide(3,2);
    can->cd(1);
    gPad->SetLogy();
    d0p->Draw();
    d0pg->Fit("gaus","","sames");
    TLegend* leg = new TLegend(0.1,0.7,0.4,0.9);
    leg->AddEntry(d0p,"All Fits","L");
    leg->AddEntry(d0pg,"Fit con>0.05","L");
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
    
    tree->Project("d0p","rec_d0-sim_d0",goodrec);
    tree->Project("p0p","rec_phi0-sim_phi0",goodrec);
    tree->Project("omp","rec_omega-sim_omega",goodrec);
    tree->Project("z0p","rec_z0-sim_z0",goodrec);
    tree->Project("tdp","rec_tandip-sim_tandip",goodrec);
    tree->Project("momp","rec_mom_mag-sim_mom_mag",goodrec);

    
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
    sgau->SetParLimits(4,0.0,0.8);
    momr->Fit("sgau","L");
//    momr->Fit("sgau","M");

  } else if (page == "mom"){
    gStyle->SetOptFit(1111);
//    TH1F* nhit = new TH1F("nhit","N hits",200,-0.5,199.5);
//    tree->Project("nhit","rec_nhit",goodrec);
    
    TH1F* mome = new TH1F("mome","estimated fit mom error",100,0.01,1.0);
    tree->Project("mome","1000*rec_mom_err",goodrec);
    mome->GetXaxis()->SetTitle("MeV");
    
    TH1F* mompg = new TH1F("mompg","momentum pull",100,-10,10);
    tree->Project("mompg","(rec_mom_mag-sim_mom_mag)/rec_mom_err",goodrec);
    
    TH1F* momr = new TH1F("momr","momentum resolution",200,-2,2);
    tree->Project("momr","1000*(rec_mom_mag-sim_mom_mag)",goodrec);
    momr->GetXaxis()->SetTitle("MeV");
    
    TH1F* mom = new TH1F("mom","Reconstructed Momentum Magnitude",200,100,110);
    tree->Project("mom","1000*(rec_mom_mag)",goodrec);
    mom->GetXaxis()->SetTitle("MeV");
    
    
    can->Clear();
    can->Divide(2,2);
    can->cd(1);
//    nhit->Draw();
//    can->cd(2);
    mome->Draw();
    can->cd(2);
    gPad->SetLogy();
    mompg->Fit("gaus","","",-2,10);
    can->cd(3);
    gPad->SetLogy();
    double integral = momr->GetEntries()*momr->GetBinWidth(1);
    sgau->SetParameters(integral,0.0,0.8*momr->GetRMS(),0.8*momr->GetRMS(),0.01,1.5*momr->GetRMS(),1.5*momr->GetRMS());
    sgau->SetParLimits(5,1.0*momr->GetRMS(),1.0);
    sgau->SetParLimits(6,1.0*momr->GetRMS(),1.0);
    sgau->SetParLimits(4,0.0,0.49);
    momr->Fit("sgau","L","",-0.5,2.0);
    can->cd(4);
    mom->Draw();
//    momr->Fit("sgau","M");

  } else if (page == "momeff"){
    const unsigned npcut(5);
    double pcutval[npcut] = {0.0,1e-3,1e-2,2e-2,5e-2};
    const unsigned necut(5);
    double ecutval[necut] = {1.0,0.25,0.2,0.18,0.16};
    TH1F* momresp[npcut];
    TH1F* momrese[necut];
    TLegend* pleg = new TLegend(0.5,0.7,0.9,0.9);
    TLegend* eleg = new TLegend(0.5,0.7,0.9,0.9);
    int colors[6] = {kRed,kBlue,kGreen,kMagenta,kCyan,kOrange};
    for(unsigned imom=0;imom<npcut;imom++){
      char name[100];
      snprintf(name,100,"momp_%d",imom);
      char pcut[100];
      snprintf(pcut,100,"rec_fitprob>%f",pcutval[imom]);
      TCut probcut(pcut);
      momresp[imom] = new TH1F(name,"momentum resolution",200,-2,3);
      tree->Project(name,"1000*(rec_mom_mag-sim_mom_mag)",rec+goodrmax+gooddip+goodndof+goodhits+probcut);
      momresp[imom]->SetMinimum(0.1);
      momresp[imom]->SetStats(0);
      momresp[imom]->GetXaxis()->SetTitle("MeV");
      momresp[imom]->SetLineColor(colors[imom]);
      double eff = momresp[imom]->GetEntries()/momresp[0]->GetEntries();
      char label[100];
      snprintf(label,100,"fit con>%5.3f, eff = %3.2f",pcutval[imom],eff);
      pleg->AddEntry(momresp[imom],label,"L");
    }
    
    for(unsigned imom=0;imom<necut;imom++){
      char name[100];
      snprintf(name,100,"mome_%d",imom);
      char ecut[100];
      snprintf(ecut,100,"1000*rec_mom_err<%f",ecutval[imom]);
      TCut errcut(ecut);
      momrese[imom] = new TH1F(name,"momentum resolution",200,-2,3);
      tree->Project(name,"1000*(rec_mom_mag-sim_mom_mag)",rec+goodrmax+gooddip+goodndof+goodhits+errcut);
      momrese[imom]->SetMinimum(0.1);
      momrese[imom]->SetStats(0);
      momrese[imom]->GetXaxis()->SetTitle("MeV");
      momrese[imom]->SetLineColor(colors[imom]);
      double eff = momrese[imom]->GetEntries()/momrese[0]->GetEntries();
      char label[100];
      snprintf(label,100,"mom err<%5.3f, eff = %3.2f",ecutval[imom],eff);
      eleg->AddEntry(momrese[imom],label,"L");
    }
    
    
    can->Clear();
    can->Divide(2,1);
    can->cd(1);
    gPad->SetLogy();
    for(unsigned imom=0;imom<npcut;imom++){
      if(imom==0)
        momresp[imom]->Draw();
      else
        momresp[imom]->Draw("same");
    }
    pleg->Draw();
    
    can->cd(2);
    gPad->SetLogy();
    for(unsigned imom=0;imom<necut;imom++){
      if(imom==0)
        momrese[imom]->Draw();
      else
        momrese[imom]->Draw("same");
    }
    eleg->Draw();
    
    
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
    TProfile* dprof = new TProfile("dprof","average reco, true transverse separation vs Z (WRT tracker center)",40,-170,170,0,0.3);
    TH2F* d2d = new TH2F("d2d","reco, true transverse separation vs Z (WRT tracker center)",50,-170,170,50,0,1.0);
    tree->Project("glen","trajdiff.endglen:trajdiff.startglen",goodrec);
    tree->Project("dg","trajdiff.endglen-trajdiff.startglen",goodrec);
    tree->Project("adiff","trajdiff.ddiff",adjacent+goodrec);
    tree->Project("mdiff","trajdiff.ddiff",middle+goodrec);
    tree->Project("dprof","trajdiff.ddiff:0.5*(trajdiff.startz+trajdiff.endz)",goodrec);
    tree->Project("d2d","trajdiff.ddiff:0.5*(trajdiff.startz+trajdiff.endz)",goodrec);


    dprof->GetXaxis()->SetTitle("cm");
    dprof->GetYaxis()->SetTitle("cm");
    dprof->SetStats(false);
    
    can->Clear();
    can->Divide(2,3);
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
//    can->Divide(1,3);
    can->cd(5);
    d2d->Draw();
    can->cd(6);
    dprof->Draw();
  } else if(page == "bintdiff") {
    TCut adjacent = ("trajdiff.endglen-trajdiff.startglen>20 && trajdiff.endglen-trajdiff.startglen<40");
    TCut middle = ("trajdiff.endglen>700 && trajdiff.endglen < 800 && trajdiff.startglen> 550 && trajdiff.startglen< 650");
    TCut early = ("trajdiff.startglen<1200");

    TH1F* tadiff = new TH1F("tadiff","True Delta P perp, adjacent stations",100,0,2e-3);
    TH1F* tmdiff = new TH1F("tmdiff","True Delta P perp, middle of tracker",100,0,2e-3);
    TH1F* radiff = new TH1F("radiff","Reco Delta P perp, adjacent stations",100,0,2e-3);
    TH1F* rmdiff = new TH1F("rmdiff","Reco Delta P perp, middle of tracker",100,0,2e-3);
    TH1F* dadiff = new TH1F("dadiff","Reco-True Delta P perp, adjacent stations",200,-2e-4,2e-4);
    TH1F* dmdiff = new TH1F("dmdiff","Reco-True Delta P perp, middle of tracker",200,-1e-4,1e-4);
    
    tadiff->GetXaxis()->SetTitle("GeV");
    tmdiff->GetXaxis()->SetTitle("GeV");
    radiff->GetXaxis()->SetTitle("GeV");
    rmdiff->GetXaxis()->SetTitle("GeV");
    dadiff->GetXaxis()->SetTitle("GeV");
    dmdiff->GetXaxis()->SetTitle("GeV");
    
    tree->Project("tadiff","trajdiff.truedpp",adjacent+goodrec+early);
    tree->Project("tmdiff","trajdiff.truedpp",middle+goodrec+early);
    tree->Project("radiff","trajdiff.recodpp",adjacent+goodrec+early);
    tree->Project("rmdiff","trajdiff.recodpp",middle+goodrec+early);
    tree->Project("dadiff","trajdiff.recodpp-trajdiff.truedpp",adjacent+goodrec+early);
    tree->Project("dmdiff","trajdiff.recodpp-trajdiff.truedpp",middle+goodrec+early);
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
  } else if(page == "binttrk") {
    TH1F* tbint = new TH1F("tbint","Field correction integral, entire track, true trajectory",100,0,0.005);
    TH1F* rbint = new TH1F("rbint","Field correction integral, entire track, reco trajectory",100,0,0.005);
    TH2F* cint = new TH2F("cint","Field correction integeral, entire track, reco vs true",50,0,0.005,50,0,0.005);
    TH1F* dint = new TH1F("dint","Field correction integral, entire track, reco - true trajectory",200,-1e-4,1e-4);
    tree->Project("tbint","simtrk.binttru",goodrec);
    tree->Project("rbint","simtrk.bintrec",goodrec);
    tree->Project("cint","simtrk.bintrec:simtrk.binttru",goodrec);
    tree->Project("dint","simtrk.bintrec-simtrk.binttru",goodrec);
    can->Clear();
    can->Divide(2,2);
    can->cd(1);
    cint->Draw();
    can->cd(2);
    rbint->Draw();
    can->cd(3);
    tbint->Draw();
    can->cd(4);
    dint->Draw();
    
  
  } else if(page == "caloresid") {
    gStyle->SetOptFit(1111);
    TCut calor("shelemnum>=10010&&shnhot>0&&shmeastype==1");
    TCut nopreshower("shmomin>0.10");
    TCut xyhit("hview==0");
    TCut zhit("hview==1");
    TH1F* zpos = new TH1F("zpos","Z track position at calo",100,160,340);
    TH1F* rpos = new TH1F("rpos","R track position at calo",100,30,75);
    TH1F* zres = new TH1F("zres","Z track residual at calo",100,-3,3);
    TH1F* rres = new TH1F("rres","R track residual at calo",100,-1,1);
    
    zpos->GetXaxis()->SetTitle("cm");
    rpos->GetXaxis()->SetTitle("cm");
    
    zres->GetXaxis()->SetTitle("cm");
    rres->GetXaxis()->SetTitle("cm");
    
    tree->Project("zpos","simhit.shz",goodrec+calor+nopreshower+xyhit);
    tree->Project("rpos","sqrt(simhit.shx^2+simhit.shy^2)",goodrec+calor+nopreshower+xyhit);
    // correct for sign convention (normal points in +phi)
//    tree->Project("zdir","simhit.mdot*simhit.sdot/sqrt(1.0-simhit.mdot^2)",goodrec+calor+nopreshower+zhit);
//    tree->Project("rdir","simhit.mdot*simhit.sdot/sqrt(1.0-simhit.mdot^2)",goodrec+calor+nopreshower+xyhit);
    // correct for projection effect of residual (which is distance in space)
    tree->Project("zres","simhit.tresid",goodrec+calor+nopreshower+zhit);
    tree->Project("rres","simhit.tresid",goodrec+calor+nopreshower+xyhit);
    
    can->Clear();
    can->Divide(2,2);
    can->cd(1);
    zpos->Draw();
    can->cd(2);
    rpos->Draw();
    can->cd(3);
//    zdir->Draw();
//    can->cd(4);
//    rdir->Draw();
//    can->cd(5);
    zres->Fit("gaus");
    can->cd(4);
    rres->Fit("gaus");
  } else if(page == "calores") {
    gStyle->SetOptFit(1111);
    TCut calor("simhit.shelemnum>10008&&simhit.shnhot>0");
    TCut nopreshower("simhit.shmomin>0.103");
    TCut xyhit("simhit.hview==0");
    TCut zhit("simhit.hview==1");

    TH2F* zrvz = new TH2F("zrvz","Z track calo res. vs Z",50,160,340,50,-3,3);
    TH2F* rrvz = new TH2F("rrvz","R track calo res. vs Z",50,160,340,50,-1,1);
    TH2F* zrvr = new TH2F("zrvr","Z track calo res. vs R",50,30,75,50,-3,3);
    TH2F* rrvr = new TH2F("rrvr","R track calo res. vs R",50,30,75,50,-1,1);

    zrvz->GetXaxis()->SetTitle("cm");
    rrvz->GetXaxis()->SetTitle("cm");

    tree->Project("zrvz","tresid:shz",goodrec+calor+nopreshower+zhit);
    tree->Project("rrvz","tresid:shz",goodrec+calor+nopreshower+xyhit);
    tree->Project("zrvr","tresid:sqrt(shx^2+shy^2)",goodrec+calor+nopreshower+zhit);
    tree->Project("rrvr","tresid:sqrt(shx^2+shy^2)",goodrec+calor+nopreshower+xyhit);

    can->Clear();
    can->Divide(2,2);
    can->cd(1);
    zrvz->Draw("box");
    can->cd(2);
    rrvz->Draw("box");
    can->cd(3);
    zrvr->Draw("box");
    can->cd(4);
    rrvr->Draw("box");
  } else if(page == "merged") {
    TH1F* nbkgevt = new TH1F("nbkgevt","# in-time DIO tracks/signal track",100,-0.5,99.5);
    TH1F* nbkghit = new TH1F("nbkghit","# in-time DIO hits/signal track",100,0,1000);
    TH1F* nmerged = new TH1F("nmerged","# replaced hits/signal track",10,-0.5,9.5);
    TH1F* nmergeda = new TH1F("nmergeda","# replaced hits/signal track",10,-0.5,9.5);
    nmerged->SetLineColor(kRed);
    nmergeda->SetLineColor(kBlue);
    
    TH1F* nshadowed = new TH1F("nshadowed","# shadowed hits/signal track",10,-0.5,9.5);
//    TH1F* mfprob = new TH1F("mfprob","Fit consistency",200,0.0,1.0);
//    TH1F* nmfprob = new TH1F("nmfprob","Fit consistency",200,0.0,1.0);
    
//    nmfprob->SetLineColor(kBlue);
//    mfprob->SetLineColor(kRed);
    tree->Project("nbkgevt","bkg_ntrks",gen);
    tree->Project("nbkghit","bkg_nhits",gen);
    tree->Project("nmerged","rec_nmerged",gen);
    tree->Project("nmergeda","rec_nmergeda",gen);
    tree->Project("nshadowed","rec_nshadowed",gen);
    tree->Project("mfprob","rec_fitprob",gen+"rec_nmergeda>0");
    tree->Project("nmfprob","rec_fitprob",gen+"rec_nmergeda==0");
    can->Clear();
    can->Divide(2,2);
    can->cd(1);
//    gPad->SetLogy();
    nbkgevt->Draw();
    can->cd(2);
//    gPad->SetLogy();
    nbkghit->Draw();
    can->cd(3);
    gPad->SetLogy();
//    nshadowed->SetLineColor(kGreen);
    nmerged->Draw();
    nmergeda->Draw("same");
    TLegend* leg = new TLegend(0.3,0.5,0.7,0.7);
    leg->AddEntry(nmerged,"replaced","L");
    leg->AddEntry(nmergeda,"replaced(active)","L");
//    leg->AddEntry(nshadowed,"shadowed","L");
    leg->Draw();
    can->cd(4);
    gPad->SetLogy();
    nshadowed->Draw();
//    can->cd(4);
//    gPad->SetLogy();
//    mfprob->Scale(10.0);
//    mfprob->Draw();
//    nmfprob->Draw("same");
   
//    TLegend* leg2 = new TLegend(0.2,0.5,0.7,0.8);
//    leg2->AddEntry(nmfprob,"#replaced==0","L");
//    leg2->AddEntry(mfprob,"#replaced>0 (X10)","L");
//    leg2->Draw();
    
  } else if(page == "bkg") {
    TH1F* nhit = new TH1F("nhit","# track hits",60,-0.5,59.5);
    TH1F* rad = new TH1F("rad","Origin transverse radius",200,0,11);
    TH1F* radh = new TH1F("radh","Origin transverse radius",200,0,11);
    TH1F* mom = new TH1F("mom","Momentum",200,0,100.);
    TH1F* momh = new TH1F("momh","Momentum",200,0,100.);
//    TH1F* pt = new TH1F("pt","Transverse momentum",200,0,100.);
//    TH1F* pth = new TH1F("pth","Transverse momentum",200,0,100.);
    TH1F* cost = new TH1F("cost","cos(#theta)",100,-1.0,1.0);
    TH1F* costh = new TH1F("costh","cos(#theta)",100,-1.0,1.0);

    rad->SetLineColor(kBlue);
    rad->GetXaxis()->SetTitle("cm");
    rad->SetMinimum(0.5);
    rad->SetStats(0);
    radh->SetLineColor(kRed);
    radh->GetXaxis()->SetTitle("cm");
    mom->GetXaxis()->SetTitle("MeV");
    mom->SetMinimum(0.5);
    mom->SetLineColor(kBlue);
    mom->SetStats(0);
    momh->GetXaxis()->SetTitle("MeV");
    momh->SetLineColor(kRed);
//    pt->SetLineColor(kBlue);
//    pt->SetMinimum(0.5);
//    pth->SetLineColor(kRed);
//    pt->GetXaxis()->SetTitle("MeV");
//    pth->GetXaxis()->SetTitle("MeV");
    cost->SetLineColor(kBlue);
    cost->SetMinimum(0.5);
    cost->SetStats(0);
    costh->SetLineColor(kRed);

    int max=100000000;
    tree->Project("nhit","nwiremeas",gen,"",max);
    tree->Project("rad","sqrt(sim_inipos_x^2+sim_inipos_y^2)",gen,"",max);
    tree->Project("radh","sqrt(sim_inipos_x^2+sim_inipos_y^2)",gen+"nwiremeas>0","",max);
    tree->Project("mom","1000*sim_mom_mag",gen,"",max);
    tree->Project("momh","1000*sim_mom_mag",gen+"nwiremeas>0","",max);
//    tree->Project("pt","1000*sim_mom_pt","","",max);
//    tree->Project("pth","1000*sim_mom_pt","nwiremeas>0","",max);
    tree->Project("cost","sim_mom_cost",gen,"",max);
    tree->Project("costh","sim_mom_cost",gen+"nwiremeas>0","",max);

    can->Clear();
    can->Divide(2,2);
    can->cd(1);
    gPad->SetLogy();
    nhit->Draw();
    can->cd(2);
    gPad->SetLogy();
    rad->Draw();
    radh->Draw("sames");
    can->cd(3);
    gPad->SetLogy();
    mom->Draw();
    momh->Draw("sames");
    can->cd(4);
    gPad->SetLogy();
    cost->Draw();
    costh->Draw("sames");
//    pt->Draw();
//    pth->Draw("same");
    TLegend* leg = new TLegend(0.3,0.5,0.7,0.7);
    leg->AddEntry(rad,"All","L");
    leg->AddEntry(radh,"# Hits>0","L");
    leg->Draw();
  } else if(page=="eloss"){
    TCut target("shelemnum==1");
    TCut absorber("shelemnum==0");    
    
//    TH1F* recmom = new TH1F("recmom","e- momentum, target+absorber",201,80,108);
//    TH1F* trumom = new TH1F("trumom","e- momentum, target+absorber",201,80,108);

    TH1F* grecmom = new TH1F("grecmom","e- momentum, target+absorber, good tracks",201,100,111);
    TH1F* gtrumom = new TH1F("gtrumom","e- momentum, target+absorber, good tracks",201,100,111);
    TH1F* dmom = new TH1F("dmom","e- energy loss to first hit",201,-0.01,20);
    
    TH1F* tdmom = new TH1F("tdmom","Energy loss/Intersection",201,0,50);
    TH1F* admom = new TH1F("admom","Energy loss/intersection",201,0,50);
    
    TH1F* ntar = new TH1F("ntar","# intersections/track",10,-0.5,9.5);
    TH1F* nabs = new TH1F("nabs","# intersections/track",10,-0.5,9.5);
    
//    recmom->SetStats(0);
//    trumom->SetStats(0);
    grecmom->SetStats(0);
    gtrumom->SetStats(0);
    tdmom->SetStats(0);
    admom->SetStats(0);
    ntar->SetStats(0);
    nabs->SetStats(0);
    
    tree->Project("tdmom","1000*(shmomin-shmomout)",target+gen);
    tree->Project("admom","1000*(shmomin-shmomout)",absorber+gen);
    
    tree->Project("grecmom","1000*rec_mom_mag",goodrec);
    tree->Project("gtrumom","1000*shmomin[ifirsthit]",goodrec);
    tree->Project("dmom","1000*(shmomin[0]-shmomin[ifirsthit])",goodrec);
    
    tree->Project("ntar","sim_ntarget",goodrec);
    tree->Project("nabs","sim_nabsorber",goodrec);

//    tree->Project("recmom","1000*rec_mom_mag",rec);
///    tree->Project("trumom","1000*shmomin[ifirsthit]",rec);

//    recmom->SetLineColor(kRed);
//    trumom->SetLineColor(kBlue);
//    recmom->GetXaxis()->SetTitle("MeV");
//    trumom->GetXaxis()->SetTitle("MeV");

    grecmom->SetLineColor(kRed);
    gtrumom->SetLineColor(kBlue);
    grecmom->GetXaxis()->SetTitle("MeV");
    gtrumom->GetXaxis()->SetTitle("MeV");
    tdmom->GetXaxis()->SetTitle("MeV");
    admom->GetXaxis()->SetTitle("MeV");
    dmom->GetXaxis()->SetTitle("MeV");
    
    tdmom->SetLineColor(kRed);
    admom->SetLineColor(kBlue);

    ntar->SetLineColor(kRed);
    nabs->SetLineColor(kBlue);
    

//    TLegend* leg = new TLegend(0.2,0.65,0.85,0.8);
    char title[100];
//    sprintf(title,"True mom at 1st track hit, mean=%4.1f MeV",trumom->GetMean());    
//    leg->AddEntry(trumom,title,"L");
//    sprintf(title,"Reconstructed momentum at origin, mean=%4.1f MeV",recmom->GetMean());    
//    leg->AddEntry(recmom,title,"L");
    
    TLegend* gleg = new TLegend(0.1,0.65,0.75,0.8);
    sprintf(title,"True mom at 1st hit, mean=%4.1f MeV",gtrumom->GetMean());    
    gleg->AddEntry(gtrumom,title,"L");
    sprintf(title,"Reco mom at origin, mean=%4.1f MeV",grecmom->GetMean());    
    gleg->AddEntry(grecmom,title,"L");
    
    TLegend* nleg = new TLegend(0.5,0.65,0.95,0.8);
    sprintf(title,"target, mean=%4.1f",ntar->GetMean());    
    nleg->AddEntry(ntar,title,"L");
    sprintf(title,"absorber, mean=%4.1f",nabs->GetMean());    
    nleg->AddEntry(nabs,title,"L");

    TLegend* deleg = new TLegend(0.2,0.65,0.85,0.8);
    sprintf(title,"target, mean=%4.2f MeV",tdmom->GetMean());    
    deleg->AddEntry(tdmom,title,"L");
    sprintf(title,"absorber, mean=%4.2f MeV",admom->GetMean());    
    deleg->AddEntry(admom,title,"L");
    
    can->Clear();
    can->Divide(2,2);
    can->cd(1);
    gPad->SetLogy();
    grecmom->Draw();
    gtrumom->Draw("same");
    gleg->Draw();

    can->cd(2);
    gPad->SetLogy();
    dmom->Draw();
    
    can->cd(3);
    ntar->Draw();
    nabs->Draw("same");
    nleg->Draw();
    
    can->cd(4);
    gPad->SetLogy();
    admom->Draw();
    tdmom->Draw("same");
    deleg->Draw();
    
  } else if(page == "hits") {
    gStyle->SetOptStat(0);
    TCut bkg("sim_mom_mag<0.1");
    TCut sig("sim_mom_mag>0.1");
//    TCut intime("shtime[ifirsthit]<0.05e-6");
    TCut trkhit("simhit.shmeastype>0&&simhit.shz<170");

    TH1F* bmom = new TH1F("bmom","track momentum",100,30,110);
    bmom->GetXaxis()->SetTitle("MeV");
    bmom->SetLineColor(kRed);
    TH1F* smom = new TH1F("smom","track momentum",100,30,110);
    smom->GetXaxis()->SetTitle("MeV");
    smom->SetLineColor(kBlue);
    smom->SetStats(0);

    TH1F* brad = new TH1F("brad","Tracker hit transverse radius",100,37,69);
    brad->GetXaxis()->SetTitle("cm");
    brad->SetLineColor(kRed);
    TH1F* srad = new TH1F("srad","Tracker hit transverse radius",100,37,69);
    srad->GetXaxis()->SetTitle("cm");
    srad->SetLineColor(kBlue);
    srad->SetStats(0);
    
    TH1F* btime = new TH1F("btime","Tracker hit time",100,-0.01,0.1);
    btime->GetXaxis()->SetTitle("#mu seconds");
    btime->SetLineColor(kRed);
    TH1F* stime = new TH1F("stime","Tracker hit time",100,-0.01,0.1);
    stime->GetXaxis()->SetTitle("#mu seconds");
    stime->SetLineColor(kBlue);
    stime->SetStats(0);
    
    TH1F* bplen = new TH1F("bplen","particle pathlength in cell",150,0,10);
    bplen->GetXaxis()->SetTitle("cm");
    bplen->SetLineColor(kRed);
    TH1F* splen = new TH1F("splen","particle pathlength in cell",150,0,10);
    splen->GetXaxis()->SetTitle("cm");
    splen->SetLineColor(kBlue);
    splen->SetStats(0);
    
    tree->Project("bmom","1000*sim_mom_mag",bkg+gen);
    tree->Project("smom","1000*sim_mom_mag",sig+gen);
    tree->Project("brad","sqrt(simhit.shy^2+simhit.shx^2)",bkg+trkhit+gen);
    tree->Project("srad","sqrt(simhit.shy^2+simhit.shx^2)",sig+trkhit+gen);
    tree->Project("btime","1e6*(shtime-shtime[ifirsthit])",bkg+trkhit+gen);
    tree->Project("stime","1e6*(shtime-shtime[ifirsthit])",sig+trkhit+gen);
    tree->Project("bplen","simhit.shpathlen",bkg+trkhit+gen);
    tree->Project("splen","simhit.shpathlen",sig+trkhit+gen);
    can->Clear();
    can->Divide(2,2);
    can->cd(1);
    brad->Draw();
    srad->Draw("same");
    
    TLegend* leg = new TLegend(0.5,0.65,0.9,0.9);
    leg->AddEntry(bmom,"DIO hits","L");
    leg->AddEntry(smom,"Conversion hits","L");
    leg->Draw();
    
    can->cd(2);
    btime->Draw();
    stime->Draw("same");
    can->cd(3);
    bmom->Draw();
    smom->Draw("same");
    can->cd(4);
    bplen->Draw();
    splen->Draw("same");
  } else if (page == "cost") {
    
    TH2F* cttvct = new TH2F("cttvct",";cos(#theta prod);cos(#theta rec)",400,-0.7,0.85,400,0.45,0.95);
    cttvct->SetStats(0);
    TH1F* cost = new TH1F("cost","#epsilon vs cos(#theta);prod cos(#theta);#epsilon",100,-1,1);
    TH1F* cose = new TH1F("cose","#epsilon vs cos(#theta);prod cos(#theta);#epsilon",100,-1,1);
    tree->Project("cost","sim_mom_cost");
    tree->Project("cose","sim_mom_cost","rec_nhit>19");
    cose->Divide(cost);
    tree->Project("cttvct","sqrt(rec_mom_mag^2-rec_mom_pt^2)/rec_mom_mag:sim_mom_cost","rec_ndof>15");
    can->Clear();
    can->Divide(1,2);
    can->cd(1);
    cttvct->Draw();
    can->cd(2);
    cose->Draw();
  } else if(page == "grad") {
    TProfile* angle = new TProfile("angle","angle ratio vs Z_{prod};Z_{prod};sin(#theta)_{tracker}/sin(#theta)_{prod}",100,-490,-380,0.0,5);
    angle->SetMinimum(0.75);
    angle->SetMaximum(0.9);
    angle->SetStats(0);
    angle->SetMarkerStyle(30);
    angle->SetMarkerColor(kRed);
    tree->Project("angle","sqrt(1.0-simt_mom_cost^2)/sqrt(1.0-sim_mom_cost^2):sim_inipos_z","simt_mom_mag>0");
    can->Clear();
    angle->Draw();
    TF1* bgrad = new TF1("bgrad","sqrt(1.0/(1.0-0.003*(x-[0])))",-500,-300);
    bgrad->SetParameter(0,-270.0);
    bgrad->SetLineColor(kBlue);
    bgrad->Draw("same");
    TLegend* leg = new TLegend(0.1,0.65,0.5,0.9);
    leg->AddEntry(angle,"conversion electrons","P");
    leg->AddEntry(bgrad,"#sqrt{B_{tracker}/B_{prod}}","L");
    leg->Draw();
  } else if(page == "dioeff"){
    TH1F* trmax = new TH1F("trmax","Maximum radius at tracker",100,0,60);
    TH1F* trmaxw = new TH1F("trmaxw","Maximum radius at tracker",100,0,60);
    TH1F* prmax = new TH1F("prmax","Maximum radius at production",100,0,60);
    TH1F* prmaxw = new TH1F("prmaxw","Maximum radius at production",100,0,60);
    
    tree->Project("trmax","2.0/simt_omega+simt_d0","simt_mom_mag>0");
    tree->Project("trmaxw","2.0/simt_omega+simt_d0","simt_mom_mag>0 && nwiremeas>0");
    
    tree->Project("prmax","2.0/sim_omega+sim_d0","simt_mom_mag>0");
    tree->Project("prmaxw","2.0/sim_omega+sim_d0","simt_mom_mag>0 && nwiremeas>0");
    
    trmax->SetLineColor(kBlue);
    trmaxw->SetLineColor(kRed);
    trmaxw->SetMinimum(1);
    trmax->SetStats(0);
    trmaxw->SetStats(0);
    trmax->GetXaxis()->SetTitle("cm");
    
    prmax->SetLineColor(kBlue);
    prmaxw->SetLineColor(kRed);
    prmaxw->SetMinimum(1);
    prmax->SetStats(0);
    prmaxw->SetStats(0);
    prmax->GetXaxis()->SetTitle("cm");
    
    can->Clear();
    can->Divide(1,2);
    can->cd(1);
    gPad->SetLogy();
    trmax->Draw();
    trmaxw->Draw("same");
    TLegend* leg = new TLegend(0.6,0.6,0.9,0.9);
    leg->AddEntry(trmax,"All DIO","L");
    leg->AddEntry(trmaxw,"DIO Nhits>0","L");
    leg->Draw();
    can->cd(2);
    gPad->SetLogy();
    prmax->Draw();
    prmaxw->Draw("same");
        
  } else if (page=="origin2d"){
    TH2F* orig = new TH2F("orig","Dio origin;Z (cm);Radius (cm)",17,-472.5,-387.5,50,0,10.5);
    TH2F* origh = new TH2F("origh","Dio origin;Z (cm);Radius (cm)",17,-472.5,-387.5,50,0,10.5);
    tree->Project("orig","sqrt(sim_inipos_x^2+sim_inipos_y^2):sim_inipos_z");
    tree->Project("origh","sqrt(sim_inipos_x^2+sim_inipos_y^2):sim_inipos_z","nwiremeas>0");
    origh->SetLineColor(kRed);
    orig->SetStats(0);
    origh->SetStats(0);
    can->Clear();
    can->Divide(1,2);
    can->cd(1);
    orig->Draw("box");
    can->cd(2);
    origh->Draw("box");
    TLegend* leg = new TLegend(.5,.7,.9,.9);
    leg->AddEntry(orig,"All DIO","L");
    leg->AddEntry(origh,"DIO, # hits>0","L");
    leg->Draw();
  } else if(page == "acc"){
    TCut tracker("nwiremeas>0");
    TH1F* cutf = new TH1F("cutf","Acceptance",8,-0.5,7.5);
    TH1F* cutfs = new TH1F("cutfs","Acceptance",8,-0.5,7.5);
    TH1F* ncutf = new TH1F("ncutf","Acceptance",8,-0.5,7.5);
    cutf->GetXaxis()->SetBinLabel(1,"All Conversions");
    cutf->GetXaxis()->SetBinLabel(2,"Reaches Tracker");
    cutf->GetXaxis()->SetBinLabel(3,"Track Reconstructed");
    cutf->GetXaxis()->SetBinLabel(4,"Pitch Angle");
    cutf->GetXaxis()->SetBinLabel(5,"Fiducial Geometry");
    cutf->GetXaxis()->SetBinLabel(6,"Fit Quality");
    cutf->GetXaxis()->SetBinLabel(7,"Mom. Resolution");
    cutf->GetXaxis()->SetBinLabel(8,"Signal Box");
    
    tree->Project("cutf","0",gen);
    tree->Project("+cutf","1",gen+tracker);
    tree->Project("+cutf","2",gen+tracker+rec);
    tree->Project("+cutf","3",gen+tracker+rec+gooddip);
    tree->Project("+cutf","4",gen+tracker+rec+gooddip+goodradius);
    tree->Project("+cutf","5",gen+tracker+rec+goodradius+gooddip+goodfit);
    tree->Project("+cutf","6",gen+tracker+rec+goodradius+gooddip+goodfit+goodres);
    tree->Project("+cutf","7",gen+tracker+rec+goodradius+gooddip+goodfit+goodres+signalbox);
    
    tree->Project("cutfs","0",gen);
    tree->Project("+cutfs","1",gen+tracker);
    tree->Project("+cutfs","2",gen+rec);
    tree->Project("+cutfs","3",gen+rec+gooddip);
    tree->Project("+cutfs","4",gen+rec+goodradius);
    tree->Project("+cutfs","5",gen+rec+goodfit);
    tree->Project("+cutfs","6",gen+rec+goodres);
    tree->Project("+cutfs","7",gen+rec+signalbox);
    
    tree->Project("ncutf","0",gen);
    tree->Project("+ncutf","1",gen);
    tree->Project("+ncutf","2",gen);
    tree->Project("+ncutf","3",gen);
    tree->Project("+ncutf","4",gen);
    tree->Project("+ncutf","5",gen);
    tree->Project("+ncutf","6",gen);
    tree->Project("+ncutf","7",gen);

    cutf->Divide(ncutf);
    cutf->SetStats(0);
    cutf->SetLineColor(kRed);
    cutf->SetMinimum(0.0);
    cutf->SetMaximum(1.01);
    cutf->SetLineWidth(2);

    cutfs->Divide(ncutf);
    cutfs->SetStats(0);
    cutfs->SetLineColor(kBlue);
    cutfs->SetMinimum(0.0);
    cutfs->SetMaximum(1.01);
    cutfs->SetLineWidth(2);

    can->Clear();
    cutf->Draw();
    cutfs->Draw("same");
    
    TLegend* leg = new TLegend(.5,.7,.9,.9);
    leg->AddEntry(cutfs,"Individual","L");
    leg->AddEntry(cutf,"Cumulative","L");
    leg->Draw();
    
    
  } else if(page == "accres"){
    TH1F* mom = new TH1F("mom","Reconstructed Momentum Magnitude",200,100,110);
    TH1F* momp = new TH1F("momp","Reconstructed Momentum Magnitude",200,100,110);
    TH1F* momq = new TH1F("momq","Reconstructed Momentum Magnitude",200,100,110);
    TH1F* momr = new TH1F("momr","Reconstructed Momentum Magnitude",200,100,110);
    
    tree->Project("mom","1000*(rec_mom_mag)",rec);
    tree->Project("momp","1000*(rec_mom_mag)",rec+gooddip+goodradius);
    tree->Project("momq","1000*(rec_mom_mag)",rec+gooddip+goodradius+goodfit+goodhits);
    tree->Project("momr","1000*(rec_mom_mag)",rec+gooddip+goodradius+goodfit+goodhits+goodres);
    
    mom->GetXaxis()->SetTitle("MeV");
    mom->SetStats(0);
    mom->SetLineColor(kBlack);
    momp->SetLineColor(kBlue);
    momq->SetLineColor(kGreen);
    momr->SetLineColor(kRed);
    can->Clear();
    mom->Draw();
    momp->Draw("same");
    momq->Draw("same");
    momr->Draw("same");
    TLegend* leg = new TLegend(.65,.7,.9,.9);
    leg->AddEntry(mom,"All Reco","L");
    leg->AddEntry(momp,"Geometry","L");
    leg->AddEntry(momq,"+ Fit Quality","L");
    leg->AddEntry(momr,"+ Momentum Error","L");
    leg->Draw();

  } else if(page=="mustop") {
    TCut muon("sim_pdgid==13");
    TH1F* stime = new TH1F("stime","Muon stop time",100,0,1e-6);
    TH1F* stimes = new TH1F("stimes","Muon stop time",100,0,1e-6);
    tree->Project("stime","shtime[nsimhit-1]",muon+"sheffect[nsimhit-1]==7");
    TF1* lognorm = new TF1("lognorm","[0]*TMath::LogNormal(x,[1],[2],[3])");
    lognorm->SetParName(0,"norm");
    lognorm->SetParName(1,"sigma");
    lognorm->SetParName(2,"theta");
    lognorm->SetParName(3,"scale");
    lognorm->SetParameters(1.0e-8*stime->GetEntries(),0.4,3.6e-8,2.2e-7);
    
    can->Clear();
    can->Divide(2);
    can->cd(1);
    gStyle->SetOptFit(111111);
    stime->Fit("lognorm");
    can->cd(2);
    
    double sigma = lognorm->GetParameter(1);
    double theta = lognorm->GetParameter(2);
    double scale = lognorm->GetParameter(3);
    TRandom3 myrand;
    for(unsigned iln=0;iln<stime->GetEntries();iln++){
      double y = myrand.Uniform();
      double x = ROOT::Math::lognormal_quantile(y,theta,sigma)*scale + theta;
//      double ynew = lognorm->Eval(x)/stime->GetEntries();
//      cout << "y = " << y << " x = " << x << endl;
      stimes->Fill(x);
    }  
    stimes->Draw();
  }
}

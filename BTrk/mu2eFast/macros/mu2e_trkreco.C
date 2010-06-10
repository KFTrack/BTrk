#include <rtypes.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TTree.h>
#include <TCut.h>
#include <TLegend.h>
#include <TString.h>
#include <TStyle.h>


Double_t splitgaus(Double_t *x, Double_t *par) {
  Double_t retval;
  Float_t xval = x[0];
  if(xval > par[1]) {
    retval = exp(-0.5*pow((xval-par[1])/par[2],2));
  } else {
    retval = exp(-0.5*pow((xval-par[1])/par[3],2));  
  }
  retval *= par[0]*0.398942/par[2];
  return retval;
}

void mu2e_trkreco(TCanvas* can, TTree* tree, const char* cpage="rec" ) {
  TString page(cpage);
  TCut rec("rec_ndof>0");
  TCut goodfit("rec_fitprob>0.05&&rec_ndof>=15&&rec_mom_err<0.0005");
  
  TF1* sgau = new TF1("sgau",splitgaus,-100.,100.,4);
  if( page == "sim"){
    TH1F* mom = new TH1F("mom","momentum",100,0.09,0.11);
    TH1F* td = new TH1F("td","tandip",100,0,1.4);
    TH1F* z0 = new TH1F("z0","z position",100,-330,-240);
    TH1F* d0 = new TH1F("d0","transverse position",100,-25,25);

    tree->Project("mom","sim_mom_mag");
    tree->Project("td","sim_tandip");
    tree->Project("z0","sim_z0");
    tree->Project("d0","sim_d0");
    can->Clear();
    can->Divide(2,2);
    can->cd(1);
    mom->Draw();
    can->cd(2);
    td->Draw();
    can->cd(3);
    z0->Draw();
    can->cd(4);
    d0->Draw();
  } else if(page == "rec"){

    TH1F* ndof = new TH1F("ndof","N DOF",40,-0.5,39.5);
    tree->Project("ndof","rec_ndof");
    
    TH1F* nhit = new TH1F("nhit","N hits",50,-0.5,49.5);
    tree->Project("nhit","rec_nhit");

    TH1F* nactive = new TH1F("nactive","N active hits",50,-0.5,49.5);
    tree->Project("nactive","rec_nactive");    

    TH1F* chindof = new TH1F("chindof","Chisquare/NDOF", 100,0.0, 10.0);
    tree->Project("chindof","rec_chisqr/rec_ndof",rec);
    
    TH1F* fitprob = new TH1F("fitprob","fit consistency",200,0.0,1.0);
    tree->Project("fitprob","rec_fitprob",rec);
    
    can->Clear();
    can->Divide(3,2);
    can->cd(1);
    ndof->Draw();
    can->cd(2);
    nhit->Draw();
    can->cd(3);
    nactive->Draw();

    can->cd(4);
    chindof->Draw();
    can->cd(5);
    fitprob->Fit("pol1","","",0.05,1);
    can->cd(6);
    gPad->SetLogy();
    fitprob->Draw();
    
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
    TH1F* momr = new TH1F("momp","momentum resolution",100,-0.0025,0.0025);
    
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
    sgau->SetParameters(integral,-0.0004,momr->GetRMS(),momr->GetRMS()/2.0);
    
    momr->Fit("sgau");

  }
}
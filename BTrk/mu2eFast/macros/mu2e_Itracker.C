{
  gROOT->LoadMacro("../mu2eFast/macros/functions.C+");
  
  TCut ITgood("rec_nactive>40");
  TCut goodfit("rec_fitprob>0.05&&rec_ndof>=10&&rec_mom_err<0.0005&&abs(rec_d0)<10.0");
  TCut endplate("shelemnum[1]==232");
  TCut innercyl("shelemnum[1]==230");
  TCut twowalls("nstraw==2");
  TCut fourwalls("nstraw>=4");
  
  TF1* sgau = new TF1("sgau",splitgaus,-1.,1.,7);
  sgau->SetParName(0,"Norm");
  sgau->SetParName(1,"Mean");
  sgau->SetParName(2,"SigH");
  sgau->SetParName(3,"SigL");
  sgau->SetParName(4,"TFH");
  sgau->SetParName(5,"TSigH");
  sgau->SetParName(6,"TSigL");
  
  TFile I("/Volumes/brownd/SuperB/FastSim/V0.2.6/workdir/Itracker_test.root");
  TFile IC_CFwalls("/Volumes/brownd/SuperB/FastSim/V0.2.6/workdir/IC_CFwalls.root");
  TFile IC_nogas("/Volumes/brownd/SuperB/FastSim/V0.2.6/workdir/IC_nogas.root");
  TFile IC_nowalls("/Volumes/brownd/SuperB/FastSim/V0.2.6/workdir/IC_nowalls.root");
  TFile IC("/Volumes/brownd/SuperB/FastSim/V0.2.6/workdir/IC_test.root");
  TFile IC_vacuum("/Volumes/brownd/SuperB/FastSim/V0.2.6/workdir/IC_vacuum.root");
  
  TTree* I_T = (TTree*)I.Get("tracks");
  TTree* IC_T = (TTree*)IC.Get("tracks");
  TTree* IC_CFwalls_T = (TTree*)IC_CFwalls.Get("tracks");
  TTree* IC_nowalls_T  = (TTree*)IC_nowalls.Get("tracks");
  TTree* IC_nogas_T = (TTree*)IC_nogas.Get("tracks");
  TTree* IC_vacuum_T = (TTree*)IC_vacuum.Get("tracks");
  
  if (false){
  TH1F* I_radlen = new TH1F("I_radlen","Int. Rad. Len. from start to last hit",100,0,0.02);
  TH1F* IC_radlen = new TH1F("IC_radlen","Int. Rad. Len. from start to last hit",100,0,0.02);
  TH1F* IC_CFwalls_radlen = new TH1F("IC_CFwalls_radlen","Int. Rad. Len. from start to last hit",100,0,0.02);

  I_radlen->SetLineColor(kBlack);
  IC_radlen->SetLineColor(kBlue);
  IC_CFwalls_radlen->SetLineColor(kRed);

  I_T->Project("I_radlen","shradlenint[ilasthit]",goodfit+ITgood);
  IC_T->Project("IC_radlen","shradlenint[ilasthit]",goodfit+ITgood);
  IC_CFwalls_T->Project("IC_CFwalls_radlen","shradlenint[ilasthit]",goodfit+ITgood);


  TH1F* I_nhits = new TH1F("I_nhits","N hits",201,-0.5,200.5);
  TH1F* IC_nhits = new TH1F("IC_nhits","N hits",201,-0.5,200.5);
  TH1F* IC_CFwalls_nhits = new TH1F("IC_CFwalls_nhits","N hits",201,-0.5,200.5);

  I_nhits->SetLineColor(kBlack);
  IC_nhits->SetLineColor(kBlue);
  IC_CFwalls_nhits->SetLineColor(kRed);

  I_T->Project("I_nhits","rec_nactive",goodfit);
  IC_T->Project("IC_nhits","rec_nactive",goodfit);
  IC_CFwalls_T->Project("IC_CFwalls_nhits","rec_nactive",goodfit);

  TH1F* I_nwalls = new TH1F("I_nwalls","N Walls hit",9,-0.5,8.5);
  TH1F* IC_nwalls = new TH1F("IC_nwalls","N Walls hit",9,-0.5,8.5);
  TH1F* IC_CFwalls_nwalls = new TH1F("IC_CFwalls_nwalls","N Walls hit",9,-0.5,8.5);

  I_nwalls->SetLineColor(kBlack);
  IC_nwalls->SetLineColor(kBlue);
  IC_CFwalls_nwalls->SetLineColor(kRed);

  I_T->Project("I_nwalls","nstraw",ITgood);
  IC_T->Project("IC_nwalls","nstraw",ITgood);
  IC_CFwalls_T->Project("IC_CFwalls_nwalls","nstraw",ITgood);
  
  TLegend* compleg = new TLegend(0.6,0.6,0.9,0.9);
  compleg->AddEntry(I_radlen,"Detailed Itracker","L");
  compleg->AddEntry(IC_radlen,"Crude Itracker, Al walls","L");
  compleg->AddEntry(IC_CFwalls_radlen,"Crude Itracker, CF walls","L");
  
  TH1F* end_2_radlen = new TH1F("end_2_radlen","Int. Rad. Len. from start to last hit, detailed Itracker",100,0,0.02);
  TH1F* end_4_radlen = new TH1F("end_4_radlen","Int. Rad. Len. from start to last hit, detailed Itracker",100,0,0.02);
  TH1F* cyl_2_radlen = new TH1F("cyl_2_radlen","Int. Rad. Len. from start to last hit, detailed Itracker",100,0,0.02);
  TH1F* cyl_4_radlen = new TH1F("cyl_4_radlen","Int. Rad. Len. from start to last hit, detailed Itracker",100,0,0.02);

  end_2_radlen->SetLineColor(kBlack);
  end_4_radlen->SetLineColor(kBlue);
  cyl_2_radlen->SetLineColor(kRed);
  cyl_4_radlen->SetLineColor(kGreen);

  I_T->Project("end_2_radlen","shradlenint[ilasthit]",goodfit+ITgood+endplate+twowalls);
  I_T->Project("end_4_radlen","shradlenint[ilasthit]",goodfit+ITgood+endplate+fourwalls);
  I_T->Project("cyl_2_radlen","shradlenint[ilasthit]",goodfit+ITgood+innercyl+twowalls);
  I_T->Project("cyl_4_radlen","shradlenint[ilasthit]",goodfit+ITgood+innercyl+fourwalls);

  TLegend* radleg = new TLegend(0.6,0.6,0.9,0.9);
  radleg->AddEntry(end_2_radlen,"starts endplate, hits 2 walls","L");
  radleg->AddEntry(end_4_radlen,"starts endplate, hits 4 walls","L");
  radleg->AddEntry(cyl_2_radlen,"starts inner cyl., hits 2 walls","L");
  radleg->AddEntry(cyl_4_radlen,"starts inner cyl., hits 4 walls","L");

  TCanvas* basiccomp = new TCanvas("basiccomp","Basic comparison",1200,800);
  basiccomp->Divide(2,2);
  basiccomp->cd(1);
  IC_CFwalls_radlen->Draw();
  I_radlen->Draw("same");
  IC_radlen->Draw("same");
  compleg->Draw();
  basiccomp->cd(2);
  IC_CFwalls_nhits->Draw();
  I_nhits->Draw("same");
  IC_nhits->Draw("same");
  basiccomp->cd(3);
  cyl_2_radlen->Draw();
  end_2_radlen->Draw("same");
  end_4_radlen->Draw("same");
  cyl_4_radlen->Draw("same");
  radleg->Draw();
  basiccomp->cd(4);
  IC_nwalls->Draw();
  IC_CFwalls_nwalls->Draw("same");
  I_nwalls->Draw("same");
  IC_nwalls->Draw("same");
  
  }
  gStyle->SetOptStat(0);
  
  TH1F* I_mome = new TH1F("I_mome","est. mom error, Detailed Itracker",100,0.01,0.6);
  I_T->Project("I_mome","1000*rec_mom_err",ITgood+goodfit);
  I_mome->GetXaxis()->SetTitle("MeV");
  I_mome->SetLineColor(kBlue);
  
  TH1F* I_mompg = new TH1F("I_mompg","mom pull, Detailed Itracker",100,-10,10);
  I_T->Project("I_mompg","(rec_mom_mag-sim_mom_mag)/rec_mom_err",ITgood+goodfit);
  I_mompg->SetLineColor(kBlue);
  
  TH1F* I_momr = new TH1F("I_momp","mom res, Detailed Itracker",400,-2,2);
  I_T->Project("I_momp","1000*(rec_mom_mag-sim_mom_mag)",ITgood+goodfit);
  I_momr->GetXaxis()->SetTitle("MeV");
  I_momr->SetLineColor(kBlue);
  
  
  TH1F* IC_mome = new TH1F("IC_mome","est. mom error, Crude Itracker",100,0.01,0.6);
  IC_T->Project("IC_mome","1000*rec_mom_err",ITgood+goodfit);
  IC_mome->GetXaxis()->SetTitle("MeV");
  IC_mome->SetLineColor(kBlue);
  
  TH1F* IC_mompg = new TH1F("IC_mompg","mom pull, Crude Itracker",100,-10,10);
  IC_T->Project("IC_mompg","(rec_mom_mag-sim_mom_mag)/rec_mom_err",ITgood+goodfit);
  IC_mompg->SetLineColor(kBlue);
  
  TH1F* IC_momr = new TH1F("IC_momp","mom res, Crude Itracker",400,-2,2);
  IC_T->Project("IC_momp","1000*(rec_mom_mag-sim_mom_mag)",ITgood+goodfit);
  IC_momr->GetXaxis()->SetTitle("MeV");
  IC_momr->SetLineColor(kBlue);
  
  
  TH1F* IC_CFwalls_mome = new TH1F("IC_CFwalls_mome","est. mom error, CF walls",100,0.01,0.6);
  IC_CFwalls_T->Project("IC_CFwalls_mome","1000*rec_mom_err",ITgood+goodfit);
  IC_CFwalls_mome->GetXaxis()->SetTitle("MeV");
  IC_CFwalls_mome->SetLineColor(kBlue);
  
  TH1F* IC_CFwalls_mompg = new TH1F("IC_CFwalls_mompg","mom pull, CF walls",100,-10,10);
  IC_CFwalls_T->Project("IC_CFwalls_mompg","(rec_mom_mag-sim_mom_mag)/rec_mom_err",ITgood+goodfit);
  IC_CFwalls_mompg->SetLineColor(kBlue);
  
  TH1F* IC_CFwalls_momr = new TH1F("IC_CFwalls_momp","mom res, CF walls",400,-2,2);
  IC_CFwalls_T->Project("IC_CFwalls_momp","1000*(rec_mom_mag-sim_mom_mag)",ITgood+goodfit);
  IC_CFwalls_momr->GetXaxis()->SetTitle("MeV");
  IC_CFwalls_momr->SetLineColor(kBlue);
  
  
  TH1F* IC_nowalls_mome = new TH1F("IC_nowalls_mome","est. mom error, no walls",100,0.01,0.6);
  IC_nowalls_T->Project("IC_nowalls_mome","1000*rec_mom_err",ITgood+goodfit);
  IC_nowalls_mome->GetXaxis()->SetTitle("MeV");
  IC_nowalls_mome->SetLineColor(kBlue);
  
  TH1F* IC_nowalls_mompg = new TH1F("IC_nowalls_mompg","mom pull, no walls",100,-10,10);
  IC_nowalls_T->Project("IC_nowalls_mompg","(rec_mom_mag-sim_mom_mag)/rec_mom_err",ITgood+goodfit);
  IC_nowalls_mompg->SetLineColor(kBlue);
  
  TH1F* IC_nowalls_momr = new TH1F("IC_nowalls_momp","mom res, no walls",400,-2,2);
  IC_nowalls_T->Project("IC_nowalls_momp","1000*(rec_mom_mag-sim_mom_mag)",ITgood+goodfit);
  IC_nowalls_momr->GetXaxis()->SetTitle("MeV");
  IC_nowalls_momr->SetLineColor(kBlue);


  TH1F* IC_nogas_mome = new TH1F("IC_nogas_mome","est. mom error, no walls or gas",100,0.01,0.6);
  IC_nogas_T->Project("IC_nogas_mome","1000*rec_mom_err",ITgood+goodfit);
  IC_nogas_mome->GetXaxis()->SetTitle("MeV");
  IC_nogas_mome->SetLineColor(kBlue);
  
  TH1F* IC_nogas_mompg = new TH1F("IC_nogas_mompg","mom pull, no walls or gas",100,-10,10);
  IC_nogas_T->Project("IC_nogas_mompg","(rec_mom_mag-sim_mom_mag)/rec_mom_err",ITgood+goodfit);
  IC_nogas_mompg->SetLineColor(kBlue);
  
  TH1F* IC_nogas_momr = new TH1F("IC_nogas_momp","mom res, no walls or gas",400,-2,2);
  IC_nogas_T->Project("IC_nogas_momp","1000*(rec_mom_mag-sim_mom_mag)",ITgood+goodfit);
  IC_nogas_momr->GetXaxis()->SetTitle("MeV");
  IC_nogas_momr->SetLineColor(kBlue);

  
  TH1F* IC_vacuum_mome = new TH1F("IC_vacuum_mome","est. mom error, no material",100,0.01,0.6);
  IC_vacuum_T->Project("IC_vacuum_mome","1000*rec_mom_err",ITgood+goodfit);
  IC_vacuum_mome->GetXaxis()->SetTitle("MeV");
  IC_vacuum_mome->SetLineColor(kBlue);
  
  TH1F* IC_vacuum_mompg = new TH1F("IC_vacuum_mompg","mom pull, no material",100,-10,10);
  IC_vacuum_T->Project("IC_vacuum_mompg","(rec_mom_mag-sim_mom_mag)/rec_mom_err",ITgood+goodfit);
  IC_vacuum_mompg->SetLineColor(kBlue);
  
  TH1F* IC_vacuum_momr = new TH1F("IC_vacuum_momp","mom res, no material",400,-2,2);
  IC_vacuum_T->Project("IC_vacuum_momp","1000*(rec_mom_mag-sim_mom_mag)",ITgood+goodfit);
  IC_vacuum_momr->GetXaxis()->SetTitle("MeV");
  IC_vacuum_momr->SetLineColor(kBlue);

  TCanvas* momcomp = new TCanvas("momcomp","Momentum Comparison",1500,1000);
  momcomp->Divide(3,5);
  gStyle->SetOptFit(1);
  
  momcomp->cd(1);
  I_mome->Draw();
  momcomp->cd(2);
  I_mompg->Fit("gaus");
  momcomp->cd(3);
  double integral = I_momr->GetEntries()*I_momr->GetBinWidth(1);
  sgau->SetParameters(integral,0.0,0.8*I_momr->GetRMS(),0.8*I_momr->GetRMS(),0.01,1.5*I_momr->GetRMS(),1.5*I_momr->GetRMS());
  sgau->SetParLimits(5,1.0*I_momr->GetRMS(),1.0);
  sgau->SetParLimits(6,1.0*I_momr->GetRMS(),1.0);
  sgau->SetParLimits(4,0.0,0.3);
  I_momr->Fit("sgau","L");
  
  momcomp->cd(4);
  IC_mome->Draw();
  momcomp->cd(5);
  IC_mompg->Fit("gaus");
  momcomp->cd(6);
  double integral = IC_momr->GetEntries()*IC_momr->GetBinWidth(1);
  sgau->SetParameters(integral,0.0,0.8*IC_momr->GetRMS(),0.8*IC_momr->GetRMS(),0.01,1.5*IC_momr->GetRMS(),1.5*IC_momr->GetRMS());
  sgau->SetParLimits(5,1.0*IC_momr->GetRMS(),1.0);
  sgau->SetParLimits(6,1.0*IC_momr->GetRMS(),1.0);
  sgau->SetParLimits(4,0.0,0.3);
  IC_momr->Fit("sgau","L");
  
//  momcomp->cd(4);
//  IC_CFwalls_mome->Draw();
//  momcomp->cd(5);
//  IC_CFwalls_mompg->Fit("gaus");
//  momcomp->cd(6);
//  double integral = IC_CFwalls_momr->GetEntries()*IC_CFwalls_momr->GetBinWidth(1);
//  sgau->SetParameters(integral,0.0,0.8*IC_CFwalls_momr->GetRMS(),0.8*IC_CFwalls_momr->GetRMS(),0.01,1.5*IC_CFwalls_momr->GetRMS(),1.5*IC_CFwalls_momr->GetRMS());
//  sgau->SetParLimits(5,1.0*IC_CFwalls_momr->GetRMS(),1.0);
//  sgau->SetParLimits(6,1.0*IC_CFwalls_momr->GetRMS(),1.0);
//  sgau->SetParLimits(4,0.0,0.3);
//  IC_CFwalls_momr->Fit("sgau","L");
  
  momcomp->cd(7);
  IC_nowalls_mome->Draw();
  momcomp->cd(8);
  IC_nowalls_mompg->Fit("gaus");
  momcomp->cd(9);
  double integral = IC_nowalls_momr->GetEntries()*IC_nowalls_momr->GetBinWidth(1);
  sgau->SetParameters(integral,0.0,0.8*IC_nowalls_momr->GetRMS(),0.8*IC_nowalls_momr->GetRMS(),0.01,1.5*IC_nowalls_momr->GetRMS(),1.5*IC_nowalls_momr->GetRMS());
  sgau->SetParLimits(5,1.0*IC_nowalls_momr->GetRMS(),1.0);
  sgau->SetParLimits(6,1.0*IC_nowalls_momr->GetRMS(),1.0);
  sgau->SetParLimits(4,0.0,0.3);
  IC_nowalls_momr->Fit("sgau","L");
  
  momcomp->cd(10);
  IC_nogas_mome->Draw();
  momcomp->cd(11);
  IC_nogas_mompg->Fit("gaus");
  momcomp->cd(12);
  double integral = IC_nogas_momr->GetEntries()*IC_nogas_momr->GetBinWidth(1);
  sgau->SetParameters(integral,0.0,0.8*IC_nogas_momr->GetRMS(),0.8*IC_nogas_momr->GetRMS(),0.01,1.5*IC_nogas_momr->GetRMS(),1.5*IC_nogas_momr->GetRMS());
  sgau->SetParLimits(5,1.0*IC_nogas_momr->GetRMS(),1.0);
  sgau->SetParLimits(6,1.0*IC_nogas_momr->GetRMS(),1.0);
  sgau->SetParLimits(4,0.0,0.3);
  IC_nogas_momr->Fit("sgau","L");


  momcomp->cd(13);
  IC_vacuum_mome->Draw();
  momcomp->cd(14);
  IC_vacuum_mompg->Fit("gaus");
  momcomp->cd(15);
  double integral = IC_vacuum_momr->GetEntries()*IC_vacuum_momr->GetBinWidth(1);
  sgau->SetParameters(integral,0.0,0.8*IC_vacuum_momr->GetRMS(),0.8*IC_vacuum_momr->GetRMS(),0.01,1.5*IC_vacuum_momr->GetRMS(),1.5*IC_vacuum_momr->GetRMS());
  sgau->SetParLimits(5,1.0*IC_vacuum_momr->GetRMS(),1.0);
  sgau->SetParLimits(6,1.0*IC_vacuum_momr->GetRMS(),1.0);
  sgau->SetParLimits(4,0.0,0.3);
  IC_vacuum_momr->Fit("sgau","L");
}


{
  TCut calo("simhit.shelemnum>=9990&&simhit.shnhot>0&&simhit.hview==1");
  TCut calor("simhit.shelemnum>=10010&&simhit.shnhot>0&&simhit.hview==1");
  TCut calorf("simhit.shelemnum>=10000&&simhit.shelemnum<10010&&simhit.shnhot>0&&simhit.hview==1");
  TCut calorr("simhit.shelemnum>=9990&&simhit.shelemnum<10000&&simhit.shnhot>0&&simhit.hview==1");
  TCut nopreshower("simhit.shmomin>0.1");
  
  TCut goodradius("abs(rec_d0)<10.0 && abs(2.0/rec_omega - rec_d0)<68.0");
  TCut gooddip("rec_tandip>0.5773&&rec_tandip<1.0");
  TCut goodfit("rec_fitprob>0.01&&rec_ndof>=20&&rec_mom_err<0.0005");
  TCut goodhits("rec_nhit-rec_nactive<10");
  TCut gen("sim_inipos_z<-300");
  TCut goodrec = goodhits+goodradius+gooddip+goodfit;
  
  TFile* calo3 = new TFile("/Volumes/brownd/SuperB/FastSim/V0.2.6_new/workdir/test_calo3.root");
  TFile* calo4 = new TFile("/Volumes/brownd/SuperB/FastSim/V0.2.6_new/workdir/test_calo4.root");
  TFile* calo6 = new TFile("/Volumes/brownd/SuperB/FastSim/V0.2.6_new/workdir/test_calo6.root");
  TFile* calo8 = new TFile("/Volumes/brownd/SuperB/FastSim/V0.2.6_new/workdir/test_calo8.root");
  TTree* c3 = (TTree*)calo3->Get("tracks");
  TTree* c4 = (TTree*)calo4->Get("tracks");
  TTree* c6 = (TTree*)calo6->Get("tracks");
  TTree* c8 = (TTree*)calo8->Get("tracks");
  if(c3 == 0) cout << "Cant find c3"<< endl;
  if(c4 == 0) cout << "Cant find c4"<< endl;
  if(c6 == 0) cout << "Cant find c6"<< endl;
  if(c8 == 0) cout << "Cant find c8"<< endl;

  TH1F* calz3 = new TH1F("calz3","Calorimeter hit Z, 3 vanes",200,170,340);
  TH1F* calz4 = new TH1F("calz4","Calorimeter hit Z, 4 vanes",200,170,340);
  TH1F* calz6 = new TH1F("calz6","Calorimeter hit Z, 6 vanes",200,170,340);
  TH1F* calz8 = new TH1F("calz8","Calorimeter hit Z, 8 vanes",200,170,340);
  
  TH1F* calfz3 = new TH1F("calfz3","Calorimeter hit Z, 3 vanes",200,170,340);
  TH1F* calfz4 = new TH1F("calfz4","Calorimeter hit Z, 4 vanes",200,170,340);
  TH1F* calfz6 = new TH1F("calfz6","Calorimeter hit Z, 6 vanes",200,170,340);
  TH1F* calfz8 = new TH1F("calfz8","Calorimeter hit Z, 8 vanes",200,170,340);

  TH1F* calrz3 = new TH1F("calrz3","Calorimeter hit Z, 3 vanes",200,170,340);
  TH1F* calrz4 = new TH1F("calrz4","Calorimeter hit Z, 4 vanes",200,170,340);
  TH1F* calrz6 = new TH1F("calrz6","Calorimeter hit Z, 6 vanes",200,170,340);
  TH1F* calrz8 = new TH1F("calrz8","Calorimeter hit Z, 8 vanes",200,170,340);
  
  calz3->SetLineColor(kBlack);
  calz4->SetLineColor(kBlack);
  calz6->SetLineColor(kBlack);
  calz8->SetLineColor(kBlack);
  
  calfz3->SetLineColor(kBlue);
  calfz4->SetLineColor(kBlue);
  calfz6->SetLineColor(kBlue);
  calfz8->SetLineColor(kBlue);
  
  calrz3->SetLineColor(kRed);
  calrz4->SetLineColor(kRed);
  calrz6->SetLineColor(kRed);
  calrz8->SetLineColor(kRed);
  
  calz3->GetXaxis()->SetTitle("cm WRT tracker center");
  calz4->GetXaxis()->SetTitle("cm WRT tracker center");
  calz6->GetXaxis()->SetTitle("cm WRT tracker center");
  calz8->GetXaxis()->SetTitle("cm WRT tracker center");
  
  c3->Project("calz3","shz",calor+nopreshower+goodrec);
  c4->Project("calz4","shz",calor+nopreshower+goodrec);
  c6->Project("calz6","shz",calor+nopreshower+goodrec);
  c8->Project("calz8","shz",calor+nopreshower+goodrec);

  c3->Project("calfz3","shz",calorf+nopreshower+goodrec);
  c4->Project("calfz4","shz",calorf+nopreshower+goodrec);
  c6->Project("calfz6","shz",calorf+nopreshower+goodrec);
  c8->Project("calfz8","shz",calorf+nopreshower+goodrec);

  c3->Project("calrz3","shz",calorr+nopreshower+goodrec);
  c4->Project("calrz4","shz",calorr+nopreshower+goodrec);
  c6->Project("calrz6","shz",calorr+nopreshower+goodrec);
  c8->Project("calrz8","shz",calorr+nopreshower+goodrec);
  
  TLegend* leg = new TLegend(0.5,0.4,0.8,0.7);
  leg->AddEntry(calz3,"Face","L");
  leg->AddEntry(calfz3,"Front edge","L");
  leg->AddEntry(calrz3,"Radial edge","L");

  TCanvas* can1 = new TCanvas("can1");
  can1->Divide(2,2);
  can1->cd(1);
  calz3->Draw();
  calfz3->Draw("same");
  calrz3->Draw("same");
  can1->cd(2);
  calz4->Draw();
  calfz4->Draw("same");
  calrz4->Draw("same");
  can1->cd(3);
  calz6->Draw();
  calfz6->Draw("same");
  calrz6->Draw("same");
  can1->cd(4);
  calz8->Draw();
  calfz8->Draw("same");
  calrz8->Draw("same");
  leg->Draw();
  
  TH2F* cal3 = new TH2F("cal3","Calorimeter hit R vs Z, 3 vanes",100,165,340,100,35,75);
  TH2F* cal4 = new TH2F("cal4","Calorimeter hit R vs Z, 4 vanes",100,165,340,100,35,75);
  TH2F* cal6 = new TH2F("cal6","Calorimeter hit R vs Z, 6 vanes",100,165,340,100,35,75);
  TH2F* cal8 = new TH2F("cal8","Calorimeter hit R vs Z, 8 vanes",100,165,340,100,35,75);
  
  cal3->GetXaxis()->SetTitle("cm WRT tracker center");
  cal4->GetXaxis()->SetTitle("cm WRT tracker center");
  cal6->GetXaxis()->SetTitle("cm WRT tracker center");
  cal8->GetXaxis()->SetTitle("cm WRT tracker center");
  
  cal3->GetYaxis()->SetTitle("cm WRT beamline");
  cal4->GetYaxis()->SetTitle("cm WRT beamline");
  cal6->GetYaxis()->SetTitle("cm WRT beamline");
  cal8->GetYaxis()->SetTitle("cm WRT beamline");

  TCanvas* can2 = new TCanvas("can2");
  can2->Divide(2,2);
  can2->cd(1);
  c3->Draw("sqrt(shx^2+shy^2):shz>>cal3",calo+nopreshower+goodrec);
  cout << "3 vanes, front edge : inner radius : face " << calfz3->GetEntries() << " : " 
    << calrz3->GetEntries() << " : " << calz3->GetEntries() << endl;
  can2->cd(2);
  c4->Draw("sqrt(shx^2+shy^2):shz>>cal4",calo+nopreshower+goodrec);
  cout << "4 vanes, front edge : inner radius : face " << calfz4->GetEntries() << " : " 
    << calrz4->GetEntries() << " : " << calz4->GetEntries() << endl;
  can2->cd(3);
  c6->Draw("sqrt(shx^2+shy^2):shz>>cal6",calo+nopreshower+goodrec);
  cout << "6 vanes, front edge : inner radius : face " << calfz6->GetEntries() << " : " 
    << calrz6->GetEntries() << " : " << calz6->GetEntries() << endl;
  can2->cd(4);
  c8->Draw("sqrt(shx^2+shy^2):shz>>cal8",calo+nopreshower+goodrec);  
  cout << "8 vanes, front edge : inner radius : face " << calfz8->GetEntries() << " : " 
    << calrz8->GetEntries() << " : " << calz8->GetEntries() << endl;
  
}

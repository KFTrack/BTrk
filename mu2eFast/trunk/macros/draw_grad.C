{
  gStyle->SetOptStat(0);
  TFile g1("/Volumes/HD2/groups/mu2e/data/FastSim/test_grad_210.root");
  TFile g2("/Volumes/HD2/groups/mu2e/data/FastSim/test_grad_240.root");
  TFile g3("/Volumes/HD2/groups/mu2e/data/FastSim/test_grad_270.root");
  TFile g4("/Volumes/HD2/groups/mu2e/data/FastSim/test_grad_300.root");
  TFile g5("/Volumes/HD2/groups/mu2e/data/FastSim/test_grad_330.root");
  TFile g6("/Volumes/HD2/groups/mu2e/data/FastSim/test_grad_360.root");
  TFile g7("/Volumes/HD2/groups/mu2e/data/FastSim/test_grad_390.root");
  TTree* t1 = (TTree*)g1.Get("tracks");
  TTree* t2 = (TTree*)g2.Get("tracks");
  TTree* t3 = (TTree*)g3.Get("tracks");
  TTree* t4 = (TTree*)g4.Get("tracks");
  TTree* t5 = (TTree*)g5.Get("tracks");
  TTree* t6 = (TTree*)g6.Get("tracks");
  TTree* t7 = (TTree*)g7.Get("tracks");

  TH1F* sint1 = new TH1F("sint1","reconstructed sin(#theta), Nhit>=20",100,0.5,1.01);
  TH1F* sint2 = new TH1F("sint2","reconstructed sin(#theta), Nhit>=20",100,0.5,1.01);
  TH1F* sint3 = new TH1F("sint3","reconstructed sin(#theta), Nhit>=20",100,0.5,1.01);
  TH1F* sint4 = new TH1F("sint4","reconstructed sin(#theta), Nhit>=20",100,0.5,1.01);
  TH1F* sint5 = new TH1F("sint5","reconstructed sin(#theta), Nhit>=20",100,0.5,1.01);
  TH1F* sint6 = new TH1F("sint6","reconstructed sin(#theta), Nhit>=20",100,0.5,1.01);
  TH1F* sint7 = new TH1F("sint7","reconstructed sin(#theta), Nhit>=20",100,0.5,1.01);

  TH1F* sint01 = new TH1F("sint01","produced sin(#theta), Nhit>=20",100,0.5,1.02);
  TH1F* sint02 = new TH1F("sint02","produced sin(#theta), Nhit>=20",100,0.5,1.02);
  TH1F* sint03 = new TH1F("sint03","produced sin(#theta), Nhit>=20",100,0.5,1.02);
  TH1F* sint04 = new TH1F("sint04","produced sin(#theta), Nhit>=20",100,0.5,1.02);
  TH1F* sint05 = new TH1F("sint05","produced sin(#theta), Nhit>=20",100,0.5,1.02);
  TH1F* sint06 = new TH1F("sint06","produced sin(#theta), Nhit>=20",100,0.5,1.02);
  TH1F* sint07 = new TH1F("sint07","produced sin(#theta), Nhit>=20",100,0.5,1.02);  
  
  t1->Project("sint1","rec_mom_pt/rec_mom_mag","rec_ndof>14");
  t2->Project("sint2","rec_mom_pt/rec_mom_mag","rec_ndof>14");
  t3->Project("sint3","rec_mom_pt/rec_mom_mag","rec_ndof>14");
  t4->Project("sint4","rec_mom_pt/rec_mom_mag","rec_ndof>14");
  t5->Project("sint5","rec_mom_pt/rec_mom_mag","rec_ndof>14");
  t6->Project("sint6","rec_mom_pt/rec_mom_mag","rec_ndof>14");
  t7->Project("sint7","rec_mom_pt/rec_mom_mag","rec_ndof>14");
  
  sint1->SetLineColor(kRed);
  sint2->SetLineColor(kBlue);
  sint3->SetLineColor(kGreen);
  sint4->SetLineColor(kMagenta);
  sint5->SetLineColor(kCyan);
  sint6->SetLineColor(kOrange);
  sint7->SetLineColor(kBlack);
  
  t1->Project("sint01","sqrt(1.0-sim_mom_cost^2)","rec_ndof>14");
  t2->Project("sint02","sqrt(1.0-sim_mom_cost^2)","rec_ndof>14");
  t3->Project("sint03","sqrt(1.0-sim_mom_cost^2)","rec_ndof>14");
  t4->Project("sint04","sqrt(1.0-sim_mom_cost^2)","rec_ndof>14");
  t5->Project("sint05","sqrt(1.0-sim_mom_cost^2)","rec_ndof>14");
  t6->Project("sint06","sqrt(1.0-sim_mom_cost^2)","rec_ndof>14");
  t7->Project("sint07","sqrt(1.0-sim_mom_cost^2)","rec_ndof>14");
    
  sint01->SetLineColor(kRed);
  sint02->SetLineColor(kBlue);
  sint03->SetLineColor(kGreen);
  sint04->SetLineColor(kMagenta);
  sint05->SetLineColor(kCyan);
  sint06->SetLineColor(kOrange);
  sint07->SetLineColor(kBlack);
  TLegend* leg = new TLegend(0.1,0.4,0.4,0.9);
  
  double n1 = sint1->Integral(sint1->FindBin(0.70711),sint1->FindBin(0.86603));
  double n2 = sint2->Integral(sint2->FindBin(0.70711),sint2->FindBin(0.86603));
  double n3 = sint3->Integral(sint3->FindBin(0.70711),sint3->FindBin(0.86603));
  double n4 = sint4->Integral(sint4->FindBin(0.70711),sint4->FindBin(0.86603));
  double n5 = sint5->Integral(sint5->FindBin(0.70711),sint5->FindBin(0.86603));
  double n6 = sint6->Integral(sint6->FindBin(0.70711),sint6->FindBin(0.86603));
  double n7 = sint7->Integral(sint7->FindBin(0.70711),sint7->FindBin(0.86603));
  
  
  double r[7];
  double p[7];
  double rerr[7];
  double perr[7];
  r[0] = sint1->GetMean();
  r[1] = sint2->GetMean();
  r[2] = sint3->GetMean();
  r[3] = sint4->GetMean();
  r[4] = sint5->GetMean();
  r[5] = sint6->GetMean();
  r[6] = sint7->GetMean();

  rerr[0] = sint1->GetMeanError();
  rerr[1] = sint2->GetMeanError();
  rerr[2] = sint3->GetMeanError();
  rerr[3] = sint4->GetMeanError();
  rerr[4] = sint5->GetMeanError();
  rerr[5] = sint6->GetMeanError();
  rerr[6] = sint7->GetMeanError();

  p[0] = sint01->GetMean();
  p[1] = sint02->GetMean();
  p[2] = sint03->GetMean();
  p[3] = sint04->GetMean();
  p[4] = sint05->GetMean();
  p[5] = sint06->GetMean();
  p[6] = sint07->GetMean();

  perr[0] = sint01->GetMeanError();
  perr[1] = sint02->GetMeanError();
  perr[2] = sint03->GetMeanError();
  perr[3] = sint04->GetMeanError();
  perr[4] = sint05->GetMeanError();
  perr[5] = sint06->GetMeanError();
  perr[6] = sint07->GetMeanError();



  cout << "Z = -210 cm rec <sin(t)> = " <<  sint1->GetMean() << " +- " << sint1->GetMeanError()
  << " prod <sin(t)> = " <<  sint01->GetMean() << " +- " << sint01->GetMeanError() << endl;
  cout << "Z = -240 cm rec <sin(t)> = " <<  sint2->GetMean() << " +- " << sint2->GetMeanError()
  << " prod <sin(t)> = " <<  sint02->GetMean() << " +- " << sint02->GetMeanError() << endl;
  cout << "Z = -270 cm rec <sin(t)> = " <<  sint3->GetMean() << " +- " << sint3->GetMeanError()
  << " prod <sin(t)> = " <<  sint03->GetMean() << " +- " << sint03->GetMeanError() << endl;
  cout << "Z = -300 cm rec <sin(t)> = " <<  sint4->GetMean() << " +- " << sint4->GetMeanError()
  << " prod <sin(t)> = " <<  sint04->GetMean() << " +- " << sint04->GetMeanError() << endl;
  cout << "Z = -330 cm rec <sin(t)> = " <<  sint5->GetMean() << " +- " << sint5->GetMeanError()
  << " prod <sin(t)> = " <<  sint05->GetMean() << " +- " << sint05->GetMeanError() << endl;
  cout << "Z = -360 cm rec <sin(t)> = " <<  sint6->GetMean() << " +- " << sint6->GetMeanError()
  << " prod <sin(t)> = " <<  sint06->GetMean() << " +- " << sint06->GetMeanError() << endl;
  cout << "Z = -390 cm rec <sin(t)> = " <<  sint7->GetMean() << " +- " << sint7->GetMeanError()
  << " prod <sin(t)> = " <<  sint07->GetMean() << " +- " << sint07->GetMeanError() << endl;
  
  char label[100];
  snprintf(label,100,"Zgrad= -210cm N=%4.0f",n1);  
  leg->AddEntry(sint1,label,"L");
  snprintf(label,100,"Zgrad= -240cm N=%4.0f",n2);  
  leg->AddEntry(sint2,label,"L");
  snprintf(label,100,"Zgrad= -270cm N=%4.0f",n3);  
  leg->AddEntry(sint3,label,"L");
  snprintf(label,100,"Zgrad= -300cm N=%4.0f",n4);  
  leg->AddEntry(sint4,label,"L");
  snprintf(label,100,"Zgrad= -330cm N=%4.0f",n5);  
  leg->AddEntry(sint5,label,"L");
  snprintf(label,100,"Zgrad= -360cm N=%4.0f",n6);  
  leg->AddEntry(sint6,label,"L");
  snprintf(label,100,"Zgrad= -390cm N=%4.0f",n7);  
  leg->AddEntry(sint7,label,"L");
  
  TCanvas* can = new TCanvas("grad");
  can->Clear();
  can->Divide(1,2);
  can->cd(1);
  
  sint03->Draw();
  sint02->Draw("same");
  sint05->Draw("same");
  sint04->Draw("same");
  sint01->Draw("same");
  sint06->Draw("same");
  sint07->Draw("same");

  can->cd(2);

  sint3->Draw();
  sint2->Draw("same");
  sint5->Draw("same");
  sint4->Draw("same");
  sint1->Draw("same");
  sint6->Draw("same");
  sint7->Draw("same");
  leg->Draw();
  
  TLine* tdcut1 = new TLine(0.86603,0.0,0.86603,1.0*sint5->GetMaximum());
  TLine* tdcut2 = new TLine(0.70711,0.0,0.70711,1.0*sint5->GetMaximum());
  tdcut1->Draw("same");
  tdcut2->Draw("same");
  
  double zstop = -440;
  double bgrad = -0.003;
  double zgrad[7] = {-210,-240,-270,-300,-330,-360,-390};
  double zerr[7] = {10,10,10,10,10,10,10};
  double ratio[7];
  double ratioerr[7];
  for(int iz=0;iz<7;iz++){
    ratio[iz] = r[iz]/(p[iz]);
    ratioerr[iz] = rerr[iz]/p[iz];
  }
  
  
  TCanvas* can2 = new TCanvas("comp");
  TGraphErrors* tg = new TGraphErrors(7,zgrad,ratio,zerr,ratioerr);
  tg->SetTitle("<reco sin(#theta)>/<produced sin(#theta)> vs Zgrad");
  tg->SetMarkerStyle(21);
  tg->SetMarkerColor(kRed);
  tg->Draw("AP");
  
  TF1* flip = new TF1("flip","sqrt([0]/([0]+([1]-x)*[2]))",-450,-200);
  flip->SetParameters(1.0,zstop,bgrad);
  flip->Draw("same");
  
  TLegend* leg2 = new TLegend(0.5,0.5,0.9,0.9);
  leg2->AddEntry(tg,"FastSim scan","P");
  leg2->AddEntry(flip,"Simple Model","L");
  leg2->Draw();
  
  
  
}
#include "TFile.h"
#include "TTree.h"
#include "TParticle.h"
#include "TClonesArray.h"
#include "TParticlePDG.h"
#include "TTree.h"
#include "Riostream.h"
#include <iostream>

// convert Rick's G4beamline ascii file dump into particles
void convertAsciiBeam(const char* asciifile, const char* rootfile) {
// configuration parameters: these should be read as input
  double xoffset = -7808;
  double zoffset = 7929 + 10200;

  ifstream in;
  in.open(asciifile);
   // configure the output tree
  TFile* file = new TFile(rootfile,"RECREATE");
  const char treename[] = {"Mu2eParticles"};
  TTree* _tree = new TTree(treename,treename);
  TClonesArray _particles("TParticle",10000);
  Int_t _evtnum;
  Float_t _evtwt;
  UInt_t _nsum;
  UInt_t _npar;
  _tree->Branch("EvtNr",&_evtnum);
  _tree->Branch("EvtW",&_evtwt);
  _tree->Branch("NEvt",&_nsum);
  _tree->Branch("NPar",&_npar);
  _tree->Branch("Particles",&_particles);
// read the input file: skip any comments   
// format is: 
// x y z Px Py Pz t PDGid EventID TrackID ParentID Weight
// mm mm mm MeV/c MeV/c MeV/c ns - - - - -
  char commentline[256];  
  unsigned npar(0);
  double xpos,ypos,zpos;
  double xmom,ymom,zmom,ptime;
  int pdgid, eventid, trackid, parentid;
  double weight;
  const char pound[] = {"#"};
  while (true) {
    char first =in.peek();
//    cout << "peeking found " << first << endl;
    if(strncmp(&first,pound,1)==0) {
      in.getline(commentline,256);
//      cout << " found comment " << commentline << endl;
    } else {
      in >> xpos >> ypos >> zpos >> xmom >> ymom >> zmom >> ptime >> pdgid >> eventid >> trackid >> parentid >> weight;

      
// shift to tracker coordinates and convert to cm GeV
      xpos = 0.1*(xpos-xoffset);
      ypos = 0.1*ypos;
      zpos = 0.1*(zpos-zoffset);
      
      xmom *= 1e-3;
      ymom *= 1e-3;
      zmom *= 1e-3;

//      cout << "found particle vertex = " << xpos <<","<< ypos<<","<< zpos 
//      << " momentum =  " << xmom <<","<< ymom <<","<< zmom
//      << " time = " << ptime << " pdgid = " <<  pdgid  << endl;
      
      
// only 1 particle/tree
      TParticle* part = new(_particles[0]) TParticle();
      part->SetPdgCode(pdgid);
      if(part->GetPDG() != 0){
        double mass = part->GetPDG()->Mass();
//      cout << " mass = " << mass << endl;
        double energy = sqrt(mass*mass + xmom*xmom + ymom*ymom + zmom*zmom);
        part->SetMomentum(xmom,ymom,zmom,energy);
        part->SetProductionVertex(xpos,ypos,zpos,ptime);
        part->SetWeight(weight); // all particles have same weight
        part->SetStatusCode(eventid);
        _tree->Fill();
        _particles.Clear();        
        npar++;
      }
    }

    if (!in.good()) break;
  }
  cout << " found " << npar << " particles" << endl;
// cleanup
  in.close();
  file->Write();
  file->Close();
}

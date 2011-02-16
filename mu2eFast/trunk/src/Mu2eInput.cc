//
//  Define the input file structure for reading ROOT data particles, and setup
// a structure to read/rewind etc.  David Brown, LBNL 7 Jan. 2011
//
#include "mu2eFast/Mu2eInput.hh"

Mu2eInput::Mu2eInput(){
}

Mu2eInput::~Mu2eInput(){
}

void
Mu2eInput::clear(Mu2eEvent& event,bool del) {
  if(del){
    for(std::vector<TParticle*>::iterator ipart = event._particles.begin();
    ipart != event._particles.end(); ipart++){
      delete *ipart;
    }
  }
  event._particles.clear();
}


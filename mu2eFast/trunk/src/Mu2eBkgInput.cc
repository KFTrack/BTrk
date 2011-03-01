//
//  Define the input file structure for reading ROOT data particles, and setup
// a structure to read/rewind etc.  David Brown, LBNL 7 Jan. 2011
//
#include "mu2eFast/Mu2eBkgInput.hh"
#include "PacEnv/PacConfig.hh"
#include <iostream>
using std::cerr;
using std::endl;

Mu2eBkgInput::Mu2eBkgInput(PacConfig& config) : Mu2eRootInput(config){
// normalize over the bunch time window
  _bunchtime = config.getfloat("bunchtime",1.0e-6);
  _lambda = config.getfloat("signaldecay",0.8e-6);
  _halfwindow = config.getfloat("sensitivehalfwindow",0.1e-6);
  _nbkg = config.getfloat("nbunchbkg",4e4);
  _bkgeff = config.getfloat("bkgeff",3.2e-3);
  _norm = 1.0-exp(-_bunchtime/_lambda);
  _bkgtime = config.getbool("bkgtime",true);
  unsigned rndseed = config.getint("rndseed", 123872342);
  _rng.SetSeed(rndseed);
}

Mu2eBkgInput::~Mu2eBkgInput(){
}

bool
Mu2eBkgInput::nextEvent(Mu2eEvent& event) {
  bool retval(true);
  clear(event,false);
// sample time of signal event within the bunch
  double stime = -_lambda*log(1.0 - _norm*_rng.Uniform());
// sample the # of background events in the sensitive time window around this,
// and read that many events from the bkg 
  double nsen =(_bkgeff*_nbkg/_norm)*(exp((-stime+_halfwindow)/_lambda) - exp((-stime-_halfwindow)/_lambda));
  unsigned nbkg = _rng.Poisson(nsen);
  Mu2eEvent temp;
  for(unsigned ibkg=0;ibkg<nbkg;ibkg++){
    if(Mu2eRootInput::nextEvent(temp) && temp._particles.size() > 0){
// compute a random time within the window, and move the particles to that.  By Dfn. t=0 is the signal production time
      Double_t btime = _rng.Uniform(-_halfwindow,_halfwindow);
      for(std::vector<TParticle*>::iterator ipart=temp._particles.begin();ipart!=temp._particles.end();ipart++){
        TParticle* part = *ipart;
        if(_bkgtime)btime += part->T();
        part->SetProductionVertex(part->Vx(),part->Vy(),part->Vz(),btime);
      }
      event.append(temp);
    }
// check if the file needs rewinding
    if(_nread >= _nevents)rewind();
  }
  return retval;
}


//
//  Define the input file structure for reading ROOT data particles, and setup
// a structure to read/rewind etc.  David Brown, LBNL 7 Jan. 2011
//
#include "mu2eFast/Mu2eBkgInput.hh"
#include "PacEnv/PacConfig.hh"
#include <iostream>
#include <assert.h>
using std::cerr;
using std::endl;

Mu2eBkgInput::Mu2eBkgInput(PacConfig& config) : Mu2eRootInput(config){
// normalize over the bunch time window
  _bunchtime = config.getfloat("bunchtime",1.7e-6);
// average # of background events is the # stopped * BF * efficiency
  _nbkg = config.getfloat("nstopped",0) * config.getfloat("bkgBF",0) * config.getfloat("bkgeff",0); 
  if(_nbkg == 0){
    cerr << "# stopped muons, background BF or background efficiency not specified: aborting" << endl;
    assert(0);
  }
// parameters for generating bunch fluctuations
  _nspread = config.getfloat("nbunchspread",0);
  _ymin = 0.25*pow(1-_nspread,2)/_nspread;
  _ymax = 0.25*pow(1+_nspread,2)/_nspread;
  assert(_ymax > _ymin);
  unsigned rndseed = config.getint("rndseed", 123872342);
  _rng.SetSeed(rndseed);
}

Mu2eBkgInput::~Mu2eBkgInput(){
}

bool
Mu2eBkgInput::nextEvent(Mu2eEvent& event) {
  bool retval(true);
  clear(event,true);
// sample the # of events in this bunch assuming the production of signal and bkg are proportional
  double nbunch = 2.0*_nbkg*sqrt(_nspread*_rng.Uniform(_ymin,_ymax));
// Poisson fluctuation
  unsigned nbkg = _rng.Poisson(nbunch);
  Mu2eEvent temp;
  for(unsigned ibkg=0;ibkg<nbkg;ibkg++){
    if(Mu2eRootInput::nextEvent(temp) && temp._particles.size() > 0){
      for(std::vector<TParticle*>::iterator ipart=temp._particles.begin();ipart!=temp._particles.end();ipart++){
// synchronize to the current bunch crossing.  This isn't quite right, as the # of particles/bunch are
// uncorrelated, but...
        double btime = fmod((*ipart)->T(),_bunchtime);
        (*ipart)->SetProductionVertex((*ipart)->Vx(),(*ipart)->Vy(),(*ipart)->Vz(),btime);
// must create new particles as root placement new overwrites old events
        TParticle* part = new TParticle(**ipart);
        event._particles.push_back(part);
      }
    }
// check if the file needs rewinding
    if(_nread >= _nevents)rewind();
  }
// update general event information
  event._evtnum = _evtnum;
  event._evtwt = _evtwt;
  event._nevt = _nevt;
  event._npar = _npar;
  return retval;
}


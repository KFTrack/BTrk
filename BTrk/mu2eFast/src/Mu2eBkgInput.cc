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
  _bunchtime = config.getfloat("bunchtime",0.975e-6);
  _lambda = config.getfloat("signaldecay",0.86e-6);
  _halfwindow = config.getfloat("bkghalfwindow",0.1e-6);
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
// normalization factor for the lifetime sampling
  _norm = 1.0-exp(-_bunchtime/_lambda);
// if true, add the input particle time to the background frame time.  This normally accounts for transit
// from the production to tyhe tracker
  _bkgtime = config.getbool("bkgtime",true);
  unsigned rndseed = config.getint("rndseed", 123872342);
  _rng.SetSeed(rndseed);
}

Mu2eBkgInput::~Mu2eBkgInput(){
}

bool
Mu2eBkgInput::nextEvent(Mu2eEvent& event) {
  bool retval(true);
  clear(event,true);
// sample time of signal event within the bunch
  double stime = -_lambda*log(1.0 - _norm*_rng.Uniform());
// sample the # of events in this bunch assuming the production of signal and bkg
// both are proportional to this
  double nbunch = 2.0*_nbkg*sqrt(_nspread*_rng.Uniform(_ymin,_ymax));
// sample the # of background events in the sensitive time window around this,
// weighted by the muon decay probability, with Poisson fluctuations
// and read that many events from the bkg 
  double nsen =(nbunch/_norm)*(exp((-stime+_halfwindow)/_lambda) - exp((-stime-_halfwindow)/_lambda));
  unsigned nbkg = _rng.Poisson(nsen);
  Mu2eEvent temp;
  for(unsigned ibkg=0;ibkg<nbkg;ibkg++){
    if(Mu2eRootInput::nextEvent(temp) && temp._particles.size() > 0){
// compute a random time within the window, and move the particles to that.  By Dfn. t=0 is the signal production time
      Double_t btime = _rng.Uniform(-_halfwindow,_halfwindow);
      for(std::vector<TParticle*>::iterator ipart=temp._particles.begin();ipart!=temp._particles.end();ipart++){
// must create new particles as root placement new overwrites old events
        TParticle* part = new TParticle(**ipart);
        if(_bkgtime)btime += part->T();
        part->SetProductionVertex(part->Vx(),part->Vy(),part->Vz(),btime);
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


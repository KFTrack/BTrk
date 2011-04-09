#include "mu2eFast/Mu2eSimpleInput.hh"
#include "PDT/PdtPdg.hh"
#include "PDT/Pdt.hh"
#include "PacEnv/PacConfig.hh"
#include <TRandom.h>
#include <TParticle.h>
#include "G3Data/GVertex.hh"
#include "G3Data/GTrack.hh"

Mu2eSimpleInput::Mu2eSimpleInput(PacConfig& config) : _ievt(0) {
  PdtPdg::PdgType pdgid = (PdtPdg::PdgType)config.getint("PdtPdg",11);
  _pdt = Pdt::lookup(pdgid);
  // Read helix generation parameters 
  _p_min = config.getdouble("p_min");
  _p_max = config.getdouble("p_max");
  _cost_min = config.getdouble("cost_min");
  _cost_max = config.getdouble("cost_max");
  _r0_mean = config.getdouble("r0_mean");
  _r0_sigma = config.getdouble("r0_sigma");
  _z0_mean = config.getdouble("z0_mean");
  _z0_sigma = config.getdouble("z0_sigma");
  _nevents = config.getint("nevents");
  // initialize random number
  unsigned rndseed = config.getint("rndseed", 1238783);
  _rng.SetSeed(rndseed);
}

bool
Mu2eSimpleInput::nextEvent(Mu2eEvent& event) {
// clear existing event
  clear(event,true);
  if(_ievt < _nevents){
// create a particle
    TParticle* part = create();
    if(part != 0)
      event._particles.push_back(part);
    event._evtnum = _ievt;
    event._evtwt = 1.0;
    event._nevt = 1;
    event._npar = 1;
    _ievt++;
    return true;
  } else {
    return false;
  }
}

void
Mu2eSimpleInput::rewind(){
// reset counter
  _ievt = 0;
}

TParticle*
Mu2eSimpleInput::create() {
// Generate track parameters; first, origin vertex
  double posphi =  _rng.Uniform(0, 2*M_PI);
  double dx = _rng.Gaus(0, _r0_sigma);
  double dy = _rng.Gaus(0, _r0_sigma);
  double x = _r0_mean*cos(posphi)+dx;
  double y = _r0_mean*sin(posphi)+dy;
  double z = _rng.Gaus(_z0_mean, _z0_sigma);
  // now momentum
  double mom	= fabs(_rng.Uniform(_p_min, _p_max));                // transverse momentum
  double phi	= _rng.Uniform(0, 2*M_PI);
  double cost = _rng.Uniform(_cost_min,_cost_max);
  double pz	= mom*cost;                // longitudinal momentum
  double pt = mom*sqrt(1.0-cost*cost);
  Hep3Vector momvec(pt*cos(phi), pt*sin(phi), pz);
//Create particle
  TParticle* part = new TParticle();
  part->SetPdgCode((Int_t)_pdt->pdgId());
  double mass = _pdt->mass();
  double energy = sqrt(mass*mass+mom*mom);
  part->SetMomentum(pt*cos(phi),pt*sin(phi),pz,energy);
  part->SetProductionVertex(x,y,z,0.0);
  part->SetWeight(1.0); // all particles have same weight
  part->SetStatusCode(1);
  return part;
}


#include "mu2eFast/Mu2eSimpleInput.hh"
#include "PDT/PdtPdg.hh"
#include "PDT/Pdt.hh"
#include "PacEnv/PacConfig.hh"
#include <TRandom.h>
#include <TParticle.h>
#include "G3Data/GVertex.hh"
#include "G3Data/GTrack.hh"

Mu2eSimpleInput::Mu2eSimpleInput(PacConfig& config) : _ievt(0) {
  PdtPdg::PdgType pdgid = (PdtPdg::PdgType)gconfig.getint("SG.PdtPdg",11);
  _pdt = Pdt::lookup(pdgid);
  // Read helix generation parameters 
  _p_min = gconfig.getdouble("SG.p_min");
  _p_max = gconfig.getdouble("SG.p_max");
  _cost_min = gconfig.getdouble("SG.cost_min");
  _cost_max = gconfig.getdouble("SG.cost_max");
  _r0_mean = gconfig.getdouble("SG.r0_mean");
  _r0_sigma = gconfig.getdouble("SG.r0_sigma");
  _z0_mean = gconfig.getdouble("SG.z0_mean");
  _z0_sigma = gconfig.getdouble("SG.z0_sigma");
  _nevents = gconfig.getint("SG.nevents");
  // initialize random number
  unsigned rndseed = gconfig.getint("SG.rndseed", 1238783);
  _rng.SetSeed(rndseed);
}

bool
Mu2eSimpleInput::nextEvent(Mu2eEvent& event) {
// clear existing event
  clear(event);
  if(_ievt < _nevents){
// create a particle
    TParticle* part = create();
    if(part != 0)
      event._particles.push_back(part);
// set event
    Int_t _evtnum;
    Float_t _evtwt;
    UInt_t _nevt, _npar;

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

void
Mu2eSimpleInput::clear(Mu2eEvent& event) {
  for(std::vector<TParticle*>::iterator ipart = event._particles.begin();
  ipart != event._particles.end(); ipart++){
    delete *ipart;
  }
  event._particles.clear();
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

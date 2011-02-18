#include "mu2eFast/Mu2eTargetInput.hh"
#include "PDT/PdtPdg.hh"
#include "PDT/Pdt.hh"
#include "PacEnv/PacConfig.hh"
#include <TRandom.h>
#include <TParticle.h>
#include "G3Data/GVertex.hh"
#include "G3Data/GTrack.hh"

Mu2eTargetInput::Mu2eTargetInput(PacConfig& config) : _ievt(0) {
// conversion to electron
  _pdt = Pdt::lookup((PdtPdg::PdgType)11);
// Read beam configuration 
  _beamxsig = gconfig.getdouble("TargetInput.beamxsigma");
  _beamtsig = gconfig.getdouble("TargetInput.beamthetasigma");
  _beamzlambda = gconfig.getdouble("TargetInput.beamzlambda");
  _p_min = gconfig.getdouble("TargetInput.p_min");
  _p_max = gconfig.getdouble("TargetInput.p_max");
  _cost_min = gconfig.getdouble("TargetInput.cost_min");
  _cost_max = gconfig.getdouble("TargetInput.cost_max");
  _nevents = gconfig.getint("TargetInput.nevents");
  // initialize random number
  unsigned rndseed = gconfig.getint("TargetInput.rndseed", 1238783);
  _rng.SetSeed(rndseed);
// find target geometry
  _diskradii = gconfig.getvector("TargetInput.diskradii");
  _diskz = gconfig.getvector("TargetInput.diskz");
  
}

bool
Mu2eTargetInput::nextEvent(Mu2eEvent& event) {
// clear existing event
  clear(event,true);
  if(_ievt < _nevents){
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
Mu2eTargetInput::rewind(){
// reset counter
  _ievt = 0;
}

TParticle*
Mu2eTargetInput::create() {
// compute z position using truncated exponential
  double deltaz = (_diskz.back() - _diskz.front());
  double zspace = deltaz/(_diskz.size()-1);
// extend the space to account for first and last foils
  deltaz += zspace;
  double z0 = _diskz[0] - zspace/2.0;
  static double norm = 1.0-exp(-deltaz/_beamzlambda);
  double zpos = z0 -_beamzlambda*log(1.0 - norm*_rng.Uniform());
// find nearest foil
  double dist = 2*deltaz;
  unsigned ifoil(0);
  for(unsigned jfoil=0;jfoil<_diskz.size();jfoil++){
    double dz = fabs(zpos - _diskz[jfoil]);
    if(dz < dist){
      dist = dz;
      ifoil=jfoil;
    }
  }
  double z = _diskz[ifoil];
// generate position.  It must be inside the target radius
  double x,y;
  double radius(100.0);
  while(radius > _diskradii[ifoil]){
    x = _rng.Gaus(0, _beamxsig);
    y = _rng.Gaus(0, _beamxsig);
    radius = sqrt(x*x + y*y);
  }
//  nom momentum
  double mom	= fabs(_rng.Uniform(_p_min, _p_max));
  double phi	= _rng.Uniform(0, 2*M_PI);
  double cost = _rng.Uniform(_cost_min,_cost_max);
  double pz	= mom*cost;                // longitudinal momentum
  double pt = mom*sqrt(1.0-cost*cost);
// create the particle
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

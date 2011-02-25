#include "mu2eFast/Mu2eTargetInput.hh"
#include "PDT/PdtPdg.hh"
#include "PDT/Pdt.hh"
#include "PacEnv/PacConfig.hh"
#include <TRandom.h>
#include <TParticle.h>
#include <TGraph.h>
#include <TSpline.h>

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
  int ispect = gconfig.getint("TargetInput.spectrumtype",0);
  switch(ispect) {
    case flat: default:
    _stype = flat;
    break;
    case file:
    _stype = file;
  }
  _cost_min = gconfig.getdouble("TargetInput.cost_min");
  _cost_max = gconfig.getdouble("TargetInput.cost_max");
  _nevents = gconfig.getint("TargetInput.nevents");
  // initialize random number
  unsigned rndseed = gconfig.getint("TargetInput.rndseed", 1238783);
  _rng.SetSeed(rndseed);
// find target geometry
  _diskradii = gconfig.getvector("TargetInput.diskradii");
  _diskz = gconfig.getvector("TargetInput.diskz");
  _halfthickness = 0.5*gconfig.getdouble("TargetInput.thickness",0.01);
  
// prepare the dio spectrum generator if necessary
  if(_stype == file){
    const char* sfile = gconfig.getcstr("TargetInput.spectrumfile");
    TGraph spectrum(sfile, "%lg,%lg");
    if(spectrum.GetN()>0){
      double lowedge = spectrum.GetX()[0];
      double hiedge = spectrum.GetX()[spectrum.GetN()-1];
      double range = hiedge-lowedge;
      TSpline3 spspect("spspect",&spectrum);
// integrate
      Int_t npt = 1000;
      double dx = range/npt;
      std::vector<Double_t> x;
      std::vector<Double_t> y;
      double xval = lowedge+dx/2.0;
      double yval(0.0);
      for(Int_t ipt=0;ipt<npt;ipt++){
        x.push_back(xval);
        y.push_back(yval);
        xval += dx;
        yval += spspect.Eval(xval);
      }
    // normalize
      for(Int_t ipt=0;ipt<npt;ipt++){
        y[ipt] /= yval;
      }
    // inverse integral spectrum
      _invintspect = new TGraph(x.size(),&y.front(),&x.front());
    }
  }
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
// randomize the position within the thickness
  double z = _diskz[ifoil] + _rng.Uniform(-_halfthickness,_halfthickness);
// generate position.  It must be inside the target radius
  double x,y;
  double radius(100.0);
  while(radius > _diskradii[ifoil]){
    x = _rng.Gaus(0, _beamxsig);
    y = _rng.Gaus(0, _beamxsig);
    radius = sqrt(x*x + y*y);
  }
//  nom momentum
  double mom;
  if(_stype == flat)
    mom	= fabs(_rng.Uniform(_p_min, _p_max));
  else
    mom	= fabs(_invintspect->Eval(_rng.Uniform()));
  
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

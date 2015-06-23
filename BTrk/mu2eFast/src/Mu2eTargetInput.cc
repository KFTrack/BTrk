#include "BaBar/BaBar.hh"
#include "mu2eFast/Mu2eTargetInput.hh"
#include "PacGeom/PacRingDetType.hh"
#include "PacGeom/PacDetector.hh"
#include "DetectorModel/DetSet.hh"
#include "PDT/PdtPdg.hh"
#include "PDT/Pdt.hh"
#include "PacEnv/PacConfig.hh"
#include <TRandom.h>
#include <TParticle.h>
#include <TGraph.h>
#include <TSpline.h>
#include <Math/QuantFuncMathCore.h>

#include <assert.h>
#include "PacGeom/PacPlaneDetElem.hh"
#include <iostream>
using namespace std;

Mu2eTargetInput::Mu2eTargetInput(PacConfig& config,const PacDetector* detector) : _ievt(0) {
  prepareBeam(config);
// find the target in the detector
  const std::vector<DetSet*>& sets = detector->detectorModel()->setList();
  for(std::vector<DetSet*>::const_iterator iset=sets.begin();iset!=sets.end();iset++){
    if((*iset)->setName()=="Target"){
      const std::vector<DetElem*>& elems = (*iset)->elementList();
      for(std::vector<DetElem*>::const_iterator ielem=elems.begin(); ielem!=elems.end(); ielem++){
        const DetElem* elem = *ielem;
        const PacPlaneDetElem* pelem = dynamic_cast<const PacPlaneDetElem*>(elem);
        if(pelem != 0){
          const PacRingDetType* rtype = dynamic_cast<const PacRingDetType*>(pelem->planeType());
          if(rtype !=0){
            cout << "found Target disk element " << pelem->elementName() 
            << " z = " << pelem->midpoint().z() << " radius = " << rtype->highrad() 
            << " thickness = " << rtype->thick() << endl;
            _halfthickness.push_back(rtype->thick()/2.0);
            _diskz.push_back(pelem->midpoint().z());
            _diskradii.push_back(rtype->highrad());
          }
        }
      }
    }
  }
}

void
Mu2eTargetInput::prepareBeam(PacConfig& config){
// particle type; by default, negative electrons
  _pdt = Pdt::lookup((PdtPdg::PdgType)config.getint("PdtPdg",11));
  _mass2 = _pdt->mass()*_pdt->mass();
// Read beam configuration 
  _beamxsig = config.getdouble("beamxsigma");
  _beamtsig = config.getdouble("beamthetasigma");
  _beamzlambda = config.getdouble("beamzlambda");
  _p_min = config.getdouble("p_min");
  _p_max = config.getdouble("p_max");
  int ispect = config.getint("spectrumtype",0);
  switch(ispect) {
    case flat: default:
    _stype = flat;
    break;
    case file:
    _stype = file;
  }
  _cost_min = config.getdouble("cost_min");
  _cost_max = config.getdouble("cost_max");
  _nevents = config.getint("nevents");
  // initialize random number
  unsigned rndseed = config.getint("rndseed", 1238783);
  _rng.SetSeed(rndseed);
  // prepare the spectrum generator if necessary
  if(_stype == file){
    const char* sfile = config.getcstr("spectrumfile");
    double scale = config.getdouble("spectrumscale",1.0);
    TGraph spectrum(sfile, "%lg,%lg");
    if(spectrum.GetN()>0){
      double lowedge = scale*spectrum.GetX()[0];
      double hiedge = scale*spectrum.GetX()[spectrum.GetN()-1];
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
        yval += spspect.Eval(xval/scale);
      }
    // normalize
      for(Int_t ipt=0;ipt<npt;ipt++){
        y[ipt] /= yval;
      }
    // inverse integral spectrum
      _invintspect = new TGraph(x.size(),&y.front(),&x.front());
    // integral spectrum
      TGraph intspect(x.size(),&x.front(),&y.front());
  // compute integral from _pmin to _pmax
      double spectint = intspect.Eval(_p_max) - intspect.Eval(_p_min);
      double eff = (_cost_max-_cost_min)*spectint/2.0;
      std::cout << "efficiency between " << _cost_min << " < cos(theta) < " <<_cost_max 
        << " and " << _p_min << " < P < " << _p_max << " = " << eff << std::endl;
    }
  }
  // timing parameters
  _bunchtime = config.getfloat("bunchtime",1.7e-6);
  _lifetime = config.getfloat("lifetime",0.86e-6);
  // lognormal distribution parameters (for muon stops)
  _lnsigma = config.getfloat("lnsigma",3.791);
  _lntheta = config.getfloat("lntheta",3.581e-8);
  _lnscale = config.getfloat("lnscale",2.191e-7);  
}

bool
Mu2eTargetInput::nextEvent(Mu2eEvent& event) {
// clear existing event
  clear(event,true);
  if(_ievt < _nevents){
    static TLorentzVector pos;
    static TLorentzVector mom;
    createPosition(pos);
    createMomentum(mom);
    TParticle* part = create(pos,mom);
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

Mu2eTargetInput::~Mu2eTargetInput(){
  std::cout << "Mu2eTargetInput processed " << _ievt << " events." << std::endl;
}


void
Mu2eTargetInput::rewind(){
// reset counter
  _ievt = 0;
}

TParticle*
Mu2eTargetInput::create(const TLorentzVector& pos, const TLorentzVector& mom) const {
  // create the particle
  TParticle* part = new TParticle();
  part->SetPdgCode((Int_t)_pdt->pdgId());
  part->SetMomentum(mom);
  part->SetProductionVertex(pos);
  part->SetWeight(1.0); // all particles have same weight
  part->SetStatusCode(1);
  return part;
}

void
Mu2eTargetInput::createPosition(TLorentzVector& pos) const {
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
  double z = _diskz[ifoil] + _rng.Uniform(-_halfthickness[ifoil],_halfthickness[ifoil]);
// generate position.  It must be inside the target radius
  double x,y;
  double radius(100.0);
  while(radius > _diskradii[ifoil]){
    x = _rng.Gaus(0, _beamxsig);
    y = _rng.Gaus(0, _beamxsig);
    radius = sqrt(x*x + y*y);
  }
// set the time according to the log-normal distribution for stopping muons, plus the muon decay
  double stime = ROOT::Math::lognormal_quantile(_rng.Uniform(),_lntheta,_lnsigma)*_lnscale + _lntheta + _rng.Exp(_lifetime);
// Synchronize this to the nearest bunch: this assumes infinite bunch trains
  double btime = fmod(stime,_bunchtime);
  pos = TLorentzVector(x,y,z,btime);
}
 
void
Mu2eTargetInput::createMomentum(TLorentzVector& momvec) const {
  double mom(0.0);
  if(_stype == flat)
    mom	= fabs(_rng.Uniform(_p_min, _p_max));
  else {
// loop until we are within range
    while(mom < _p_min || mom > _p_max){
      mom	= fabs(_invintspect->Eval(_rng.Uniform()));
    }
  }
  double phi	= _rng.Uniform(0, 2*M_PI);
  double cost = _rng.Uniform(_cost_min,_cost_max);
  double pz	= mom*cost;                // longitudinal momentum
  double pt = mom*sqrt(1.0-cost*cost);
  double energy = sqrt(_mass2+mom*mom);
  momvec = TLorentzVector(pt*cos(phi),pt*sin(phi),pz,energy);
}

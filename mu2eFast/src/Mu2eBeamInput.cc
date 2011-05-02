#include "BaBar/BaBar.hh"
#include "mu2eFast/Mu2eBeamInput.hh"
#include "PacEnv/PacConfig.hh"
#include <TParticle.h>
#include "PacSim/PacSimulate.hh"
#include "PacSim/PacSimTrack.hh"
using namespace std;


Mu2eBeamInput::Mu2eBeamInput(PacConfig& config,const PacSimulate* sim) :
  Mu2eTargetInput(config,sim->getDetector()),
  _sim(sim)
{
  PacConfig bconfig = config.getconfig("Beam.");
  _rinput = new Mu2eRootInput(bconfig);
  _lifetime = config.getfloat("lifetime",0.86e-6);
  _bunchtime = config.getfloat("bunchtime",1.7e-6);
}

Mu2eBeamInput::~Mu2eBeamInput(){
  std::cout << "Mu2eBeamInput read " << _rinput->nread() << " beam events." << std::endl;
}


bool
Mu2eBeamInput::nextEvent(Mu2eEvent& event) {
// clear existing event
  clear(event,true);  
  if(_ievt < _nevents){
// simulate the beam particle stopping
    bool stops(false);
    Mu2eEvent beamevt;
    std::vector<PacSimTrack*> strks;
    while(!stops){
// cleanup
      for(std::vector<PacSimTrack*>::iterator istrk=strks.begin();istrk!=strks.end();istrk++)
        delete *istrk;
      strks.clear();
// get a beam muon
      bool goodevent = _rinput->nextEvent(beamevt);
      if(goodevent){
// simulate the beam particles, requiring at least 1 to stop in the target
        for(std::vector<TParticle*>::const_iterator ipar = beamevt._particles.begin();ipar != beamevt._particles.end();ipar++){
          TParticle* part = *ipar;
          PacSimTrack* simtrk = _sim->simulateParticle(part);
          stops |= stopsInTarget(simtrk);
          strks.push_back(simtrk);
        }
      } else {
        std::cout << "end of beam input file " << std::endl;
        break;
      }
    }
// use the stopping point(s) to generate daughters
    for(std::vector<PacSimTrack*>::iterator istrk=strks.begin();istrk!=strks.end();istrk++){
      if(stopsInTarget(*istrk)){
        // randomize the position along the trajector within this element
        const PacSimHit* sthit = (*istrk)->lastHit();
        Hep3Vector dir = sthit->momentumIn().unit();
        double flen = _rng.Uniform(sthit->detIntersection().pathrange[0]-sthit->detIntersection().pathlen,
          sthit->detIntersection().pathrange[1]-sthit->detIntersection().pathlen);
        HepPoint spos = sthit->position() + dir*flen;
        // add a random decay time to the stopping time
        double stime = (*istrk)->lastHit()->time() + _rng.Exp(_lifetime);
        // Synchronize this to the nearest bunch: this assumes infinite bunch trains
        double btime = fmod(stime,_bunchtime);
        TLorentzVector pos(spos.x(),spos.y(),spos.z(),btime);
        static TLorentzVector mom;
        createMomentum(mom);
        TParticle* part = create(pos,mom);
        if(part != 0){
          event._particles.push_back(part);
          event._strks.push_back(*istrk);
        } else {
          delete *istrk;          
        }
      } else {
        delete *istrk;
      }
    }
    event._evtnum = _ievt;
    event._evtwt = 1.0;
    event._nevt = 1;
    event._npar = 1;
    _ievt++;
    return stops;
  } else {
    return false;
  }
}

void 
Mu2eBeamInput::rewind() {
  _rinput->rewind();
}


bool
Mu2eBeamInput::stopsInTarget(const PacSimTrack* strk) const {
  return strk->lastHit()->detEffect() == PacSimHit::stop &&
    strk->lastHit()->position().z() >= _diskz.front() &&
    strk->lastHit()->position().z() <= _diskz.back();
}
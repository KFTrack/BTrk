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
        HepPoint spos = (*istrk)->lastHit()->position();
        double stime = (*istrk)->lastHit()->time();
        TLorentzVector pos(spos.x(),spos.y(),spos.z(),stime);
        static TLorentzVector mom;
        createMomentum(mom);
        TParticle* part = create(pos,mom);
        if(part != 0)
          event._particles.push_back(part); 
      }
// cleanup; I would like to keep these simtracks!
      delete *istrk;
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
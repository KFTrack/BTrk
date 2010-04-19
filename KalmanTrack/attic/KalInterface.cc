//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: KalInterface.cc,v 1.11 2002/04/24 00:41:17 brownd Exp $
//
// Description:
//     
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author(s): Steve Schaffner
//
//------------------------------------------------------------------------

#include "BaBar/BaBar.hh"
#include "KalmanTrack/KalInterface.hh"
#include "KalmanTrack/KalRep.hh"
#include "KalmanTrack/KalMiniRep.hh"
#include "KalmanTrack/KalPairRep.hh"
#include "ProxyDict/IfdIntKey.hh"

KalInterface::KalInterface() : _mini(false) {}

//------------------------------------------------------------------------
KalInterface::~KalInterface() {
//------------------------------------------------------------------------
}

//------------------------------------------------------------------------
KalRep* 
KalInterface::kalRep() {
//------------------------------------------------------------------------
  return (KalRep*) myRep();
}

const KalRep* 
KalInterface::kalmanRep() const {
//------------------------------------------------------------------------
  return (const KalRep*) myConstRep();
}

//------------------------------------------------------------------------
TrkErrCode 
KalInterface::extendThrough(const TrkVolume& trkVol, trkDirection trkDir) {
//------------------------------------------------------------------------
  return kalRep()->extendThrough(trkVol, trkDir);
}

//------------------------------------------------------------------------
const IfdKey& 
KalInterface::myKey() const {
//------------------------------------------------------------------------
  static IfdIntKey _theKey(7261);
  return _theKey;
}

//------------------------------------------------------------------------
double
KalInterface::refMomentum()const {
//------------------------------------------------------------------------
  return kalmanRep()->refMomentum();
}

//------------------------------------------------------------------------
double
KalInterface::refMomFltLen()const {
//------------------------------------------------------------------------
  return kalmanRep()->refMomFltLen();
}

//------------------------------------------------------------------------
const TrkSimpTraj* 
KalInterface::seed() const {
//------------------------------------------------------------------------
  return kalmanRep()->seed();
}

bool
KalInterface::attach(TrkRep* rep) {
  bool retval = TrkExtInterface::attach(rep);
  if(retval)
    _mini = false;
  else if(_miniiface.attach(rep) &&
	  _miniiface.kalMiniRep()->active()){
    retval = true;
    _mini = true;
    setRep(_miniiface.kalMiniRep()->kalRep());
  } else if(_pairiface.attach(rep)){
    _mini = false;
    retval = true;
    setRep(_pairiface.kalPairRep());
  }
  return retval;
}

bool
KalInterface::attach(const TrkRep* rep) {
  bool retval = TrkExtInterface::attach(rep);
  if(retval)
    _mini = false;
  else if(_miniiface.attach(rep) ){
// be explicit with constness here to help the compiler figure it out
    const KalMiniRep* minirep = _miniiface.kalmanMiniRep();
    if(minirep != 0 && minirep->active()){
      retval = true;
      _mini = true;
      setRep(minirep->kalRep());
    }
  } else if(_pairiface.attach(rep)){
    _mini = false;
    retval = true;
    setRep(_pairiface.kalmanPairRep());
  }
  return retval;
}











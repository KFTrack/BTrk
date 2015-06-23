//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: KalPairInterface.cc,v 1.3 2000/12/05 19:56:22 brownd Exp $
//
// Description:
//     
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author(s): Doug Roberts
//
//------------------------------------------------------------------------

#include "BaBar/BaBar.hh"
#include "KalmanTrack/KalPairInterface.hh"
#include "KalmanTrack/KalPairRep.hh"
#include "ProxyDict/IfdIntKey.hh"

//------------------------------------------------------------------------
KalPairInterface::~KalPairInterface() {
//------------------------------------------------------------------------
}

//------------------------------------------------------------------------
KalPairRep* 
KalPairInterface::kalPairRep() {
//------------------------------------------------------------------------
  return (KalPairRep*) myRep();
}

const KalPairRep* 
KalPairInterface::kalmanPairRep() const {
//------------------------------------------------------------------------
  return (const KalPairRep*) myConstRep();
}

//------------------------------------------------------------------------
//TrkErrCode 
//KalInterface::extendThrough(const TrkVolume& trkVol, trkDirection trkDir) {
//------------------------------------------------------------------------
//  return kalRep()->extendThrough(trkVol, trkDir);
//}

//------------------------------------------------------------------------
const IfdKey& 
KalPairInterface::myKey() const {
//------------------------------------------------------------------------
  static IfdIntKey _theKey(4952957);
  return _theKey;
}

//------------------------------------------------------------------------
double
KalPairInterface::refMomentum()const {
//------------------------------------------------------------------------
  return kalmanPairRep()->refMomentum();
}

//------------------------------------------------------------------------
double
KalPairInterface::refMomFltLen()const {
//------------------------------------------------------------------------
  return kalmanPairRep()->refMomFltLen();
}

//------------------------------------------------------------------------
const TrkSimpTraj* 
KalPairInterface::seed() const {
//------------------------------------------------------------------------
  return kalmanPairRep()->seed();
}

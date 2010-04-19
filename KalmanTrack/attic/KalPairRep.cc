//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: KalPairRep.cc,v 1.24 2003/01/24 06:15:21 brownd Exp $
//
// Description:
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Copyright Information:
//	Copyright (C) 1998	Lawrence Berkeley Laboratory
//
// Author List:
//      Dave Brown 3/21/98
//      Doug Roberts
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"
#include "BbrGeom/BbrPointErr.hh"
#include "KalmanTrack/KalPairRep.hh"
#include "KalmanTrack/KalPairConduit.hh"
#include "KalmanTrack/KalPairSite.hh"
#include "ProxyDict/IfdIntKey.hh"
#include "TrkBase/TrkParams.hh"

KalPairRep::KalPairRep(TrkRecoTrk* newtrk,const KalRep* seedRep,
		       KalPairConduit* conduit) :
  KalRep(*seedRep,newtrk), // start by copying the plus branch
  _conduit(conduit)
{
  // Tell the conduit that I exist
  assert(0 != _conduit);
  _conduit->addRep(this);

  // Now add a PairSite, giving the constraints calculated by the conduit
  TrkParams params = _conduit->getConstraint(this).trackParameters();
  double fltlen = _conduit->getFltLen(this);
  bool constrainAll[5] = {true, true, true, true, true};
  KalPairSite* psite = new KalPairSite(referenceTraj(), params, constrainAll,
				       fltlen, _conduit);

  // Tell the conduit about this site
  _conduit->addSite(psite, this);

  // Insert PairSite in the rep's site list
  kalSites().push_back(psite);

  // Do we have to invalidate sites, too?
  // Make fit uncurrent
  setCurrent(false);
  resetFit();

  // update  _hitsites
  findHitSites();
}

TrkErrCode
KalPairRep::fit()
{
  assert( 0 != _conduit);  
  return _conduit->coordinateFit();
}

// int
// KalPairRep::iterations() const
// {
//   assert (0 != _conduit);
//   return _conduit->iterations();
// }

double
KalPairRep::pairFltLen()
{
  assert (0 != _conduit);
  return _conduit->getFltLen(this);
}

bool
KalPairRep::converged() const
{
  assert (0 != _conduit);
  return _conduit->converged();
}

BbrPointErr
KalPairRep::prodPoint() const
{
  assert (0 != _conduit);
  return _conduit->prodPoint();
}

double
KalPairRep::prodPointChi2() const
{
  assert (0 != _conduit);
  return _conduit->prodPointChi2();
}

KalPairRep*
KalPairRep::otherPairRep()
{
  return _conduit->otherPairRep(this);
}

KalPairRep::~KalPairRep()
{
  // Responsible for deleting the conduit
  if (0 != _conduit){
    // Tell the conduit that it is about to be deleted
    _conduit->killedBy(this);
    delete _conduit;
  }
}


//Returns arbitrary key, distinguishing KalPairRep from all other reps
const IfdKey&
KalPairRep::myKey() const {
  static IfdIntKey _theKey(4952957);
  return _theKey;
}

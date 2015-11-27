#ifndef TRKHITONTRKUPDATER_HH
#define TRKHITONTRKUPDATER_HH

//--------------------------------------------------------------------------
//
// Environment:
//      This software was developed for the BaBar collaboration.  If you
//      use all or part of it, please give an appropriate acknowledgement.
//
// Copyright Information:
//      Copyright (C) 2000      University of California, San Diego
//
//------------------------------------------------------------------------

#include "BTrk/TrkBase/TrkHit.hh"
#include "BTrk/TrkBase/TrkErrCode.hh"

/**
 *  TrkHitUpdater. 
 *     this class regulates access to some protected functions
 *     in TrkHit: only classes inheriting from this one can 
 *     modify the status of a TrkHit
 *
 *
 *  This software was developed for the BaBar collaboration.  If you
 *  use all or part of it, please give an appropriate acknowledgement.
 *
 *  Copyright (C) 2000 University of California, San Diego
 *
 *  @version $Id: TrkHitUpdater.hh,v 1.5 2002/04/22 00:44:23 raven Exp $
 *
 *  @author (Gerhard Raven)           (based on an idea of Steve Schaffner)
 */
#include "BTrk/TrkBase/TrkFunctors.hh"

class TrkHitUpdater
{
public:
  virtual ~TrkHitUpdater() = 0;
protected:
  TrkErrCode updateMeasurement(TrkHit &thit, const TrkDifTraj* traj=0) const
    { return thit.updateMeasurement(traj);}
// allow subclasses (essentially TrkReps) to set thit activity directly
  void setActivity(TrkHit& thit,bool active) const {
    thit.setActive(active); }
// allow changing the parent to which a thit is assigned
  void setParent(TrkHit& thit,TrkRep* parent) const {
    thit._parentRep = parent;
  }

  TrkBase::Functors::updateMeasurement updateMeasurement( const TrkDifTraj* traj=0) const 
  { return TrkBase::Functors::updateMeasurement(traj); }
  TrkBase::Functors::setParent setParent(TrkRep* parent) const 
  { return TrkBase::Functors::setParent(parent); }
  TrkBase::Functors::setActive setActive(bool active) const 
  { return TrkBase::Functors::setActive(active); }


};

#endif

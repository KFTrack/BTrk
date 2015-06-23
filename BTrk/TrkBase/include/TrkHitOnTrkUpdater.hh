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

#include "TrkBase/TrkHitOnTrk.hh"
#include "TrkBase/TrkErrCode.hh"

/**
 *  TrkHitOnTrkUpdater. 
 *     this class regulates access to some protected functions
 *     in TrkHitOnTrk: only classes inheriting from this one can 
 *     modify the status of a TrkHitOnTrk
 *
 *
 *  This software was developed for the BaBar collaboration.  If you
 *  use all or part of it, please give an appropriate acknowledgement.
 *
 *  Copyright (C) 2000 University of California, San Diego
 *
 *  @version $Id: TrkHitOnTrkUpdater.hh,v 1.5 2002/04/22 00:44:23 raven Exp $
 *
 *  @author (Gerhard Raven)           (based on an idea of Steve Schaffner)
 */
#include "TrkBase/TrkFunctors.hh"

class TrkHitOnTrkUpdater
{
public:
  virtual ~TrkHitOnTrkUpdater() = 0;
protected:
  TrkErrCode updateMeasurement(TrkHitOnTrk &hot, const TrkDifTraj* traj=0) const
    { return hot.updateMeasurement(traj);}
// allow subclasses (essentially TrkReps) to set hot activity directly
  void setActivity(TrkHitOnTrk& hot,bool active) const {
    hot.setActive(active); }
// allow changing the parent to which a hot is assigned
  void setParent(TrkHitOnTrk& hot,TrkRep* parent) const {
    hot._parentRep = parent;
  }

  TrkBase::Functors::updateMeasurement updateMeasurement( const TrkDifTraj* traj=0) const 
  { return TrkBase::Functors::updateMeasurement(traj); }
  TrkBase::Functors::setParent setParent(TrkRep* parent) const 
  { return TrkBase::Functors::setParent(parent); }
  TrkBase::Functors::setActive setActive(bool active) const 
  { return TrkBase::Functors::setActive(active); }


};

#endif

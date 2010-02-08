//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkTimeCalculator.hh,v 1.2 2001/07/23 22:25:25 raven Exp $
//
// Description:
//   class TrkTimeCalculator.  An abstract base class for computing a new
//   track time based on the hots in the track.
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Copyright Information:
//	Copyright (C) 2001	Lawrence Berkeley Laboratory
//
// Author List:
//      Dave Brown 7/17/01
//------------------------------------------------------------------------

#ifndef TRKTIMECALCULATOR_HH
#define TRKTIMECALCULATOR_HH

class TrkRecoTrk;
#include "TrkBase/TrkHotSelector.hh"

class TrkTimeCalculator {
public:
// only one constructor
  TrkTimeCalculator(const TrkHotSelector& selector) :
    _selector(selector){};
  virtual ~TrkTimeCalculator(){};
// The Function
  virtual bool trackTime(const TrkRecoTrk& trk,
                         double& time, double& timeerr,
                         int& nHotsUsed) const = 0;
  bool trackTime(const TrkRecoTrk& trk,
                 double& time, double& timeerr) const
  { int dummy; return trackTime(trk,time,timeerr,dummy);}
protected:
  bool useHot(const TrkHitOnTrk& x) const { return _selector.useHot(x);}
private:
  const TrkHotSelector& _selector;
// disallow
  TrkTimeCalculator(const TrkTimeCalculator&);
  TrkTimeCalculator& operator = (const TrkTimeCalculator&);
};
#endif

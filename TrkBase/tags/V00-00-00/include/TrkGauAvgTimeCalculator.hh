//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkGauAvgTimeCalculator.hh,v 1.2 2001/07/23 22:25:25 raven Exp $
//
// Description:
//   class TrkGauAvgTimeCalculator.  A simple implementation of TrkTimeCalculator
//   that takes the weighted average of all selected hot's times assuming a
//   Gaussian error.
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

#ifndef TRKGAUAVGTIMECALCULATOR_HH
#define TRKGAUAVGTIMECALCULATOR_HH

#include "TrkBase/TrkTimeCalculator.hh"

class TrkGauAvgTimeCalculator : public TrkTimeCalculator {
public:
// only a default constructor
  TrkGauAvgTimeCalculator(const TrkHotSelector& selector);
  virtual ~TrkGauAvgTimeCalculator();
// The Function
  virtual bool trackTime(const TrkRecoTrk& trk,
                         double& time, double& timeerr,
                         int &nHotsUsed) const;
private:
// disallow
  TrkGauAvgTimeCalculator(const TrkGauAvgTimeCalculator&);
  TrkGauAvgTimeCalculator& operator = (const TrkGauAvgTimeCalculator&);
};
#endif

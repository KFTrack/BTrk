//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkKalTrk.hh,v 1.5 2004/08/06 06:31:42 bartoldu Exp $
//
// Description:
//      class TrkKalTrk.  Transient interface base class for persistent
//      form of Kalman track (mini-track).
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Copyright Information:
//	Copyright (C) 2000	Lawrence Berkeley Laboratory
//
// Author List:
//      Dave Brown 10/31/00
//------------------------------------------------------------------------
#ifndef TRKKALTRK_HH
#define TRKKALTRK_HH

class TrkHitOnTrk;
class TrkHotList;
class TrkSimpTraj;
class TrkErrCode;
#include <iosfwd>

#include "PDT/PdtPid.hh"
#include "TrkBase/TrkEnums.hh"

class TrkKalTrk {
public:
// only default constructor, as there can be NO DATA MEMBERS in this class
  TrkKalTrk();
// pure virtual destructor, as this is a pure interface class
  virtual ~TrkKalTrk() = 0;
// a TrkKalTrk knows what its default PID Is
  virtual PdtPid::PidType defaultHypo() const = 0;
// it can tell you which fit a given hypo points to
  virtual PdtPid::PidType fitHypo(PdtPid::PidType hypo) const = 0;
// it knows what its (persistent) trk id is
  virtual unsigned long trackId() const = 0;
// it can return a seed trajectory.  This function _RETURNS OWNERSHIP_
  virtual TrkSimpTraj* seedTrajectory() const = 0;
// it can return the T0 used when the track was fit
  virtual double usedT0() const = 0;
// it can report the fit status of a paritcular hypo
  virtual bool isValid(PdtPid::PidType hypo) const = 0;
  virtual bool isCurrent(PdtPid::PidType hypo) const = 0;
  virtual TrkErrCode fitStatus(PdtPid::PidType hypo) const = 0;
// A KalTrk should know how many svt, dch hots and fits it has
  virtual unsigned nSvt() const = 0;
  virtual unsigned nDch() const = 0;
  virtual unsigned nFit() const = 0;
// one implemented function; printout.
  void print(std::ostream& os) const;
private:
};

std::ostream& operator<<(std::ostream& os, const TrkKalTrk&);

#endif

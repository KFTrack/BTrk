//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkFitSummary.hh,v 1.2 2003/03/27 19:58:21 brownd Exp $
//
//  Description:
//  Class TrkFitSummary: summarize the results of a track fit
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Copyright Infomation;
//	Copyright (C) 2002	Lawrence Berkeley Laboratory
//
// Author(s): Dave Brown 10/202/02
//
//------------------------------------------------------------------------

#ifndef TRKFITSUMMARY_HH
#define TRKFITSUMMARY_HH

#include "BaBar/PdtPid.hh"
#include "TrkBase/TrkEnums.hh"

class TrkSimpTraj;

class TrkFitSummary {
public:
// construct from input information
  TrkFitSummary();
  TrkFitSummary(PdtPid::PidType type,
		TrkEnums::PackFlag flag,
		double chisqprob,
		TrkSimpTraj* traj);
// copy and equivalence are OK
  TrkFitSummary(const TrkFitSummary&);
  TrkFitSummary& operator =(const TrkFitSummary&);
// needed for rw
  bool operator == (const TrkFitSummary& other) const;
  virtual ~TrkFitSummary();
// obvious accessors
  PdtPid::PidType pidType() const { return _type; }
  TrkEnums::PackFlag fitFlag() const { return _flag; }
  double chisqProb() const { return _chisqprob; }
  const TrkSimpTraj* traj() const { return _traj; }
  TrkSimpTraj* traj() { return _traj; }
private:
  PdtPid::PidType _type;
  TrkEnums::PackFlag _flag;
  double _chisqprob;
  TrkSimpTraj* _traj;
};

#endif

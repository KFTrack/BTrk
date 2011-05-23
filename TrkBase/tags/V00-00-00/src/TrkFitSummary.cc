//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkFitSummary.cc,v 1.1 2002/10/08 18:42:33 brownd Exp $
//
//  Description:
//  Class TrkFitSummary
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
#include "BaBar/BaBar.hh"
#include "TrkBase/TrkFitSummary.hh"

TrkFitSummary::TrkFitSummary() :
  _type(PdtPid::null),_flag(TrkEnums::KalFit),_chisqprob(-1.0),_traj(0)
{}

TrkFitSummary::TrkFitSummary(PdtPid::PidType type,
		TrkEnums::PackFlag flag,
		double chisqprob,
			       TrkSimpTraj* traj) :
  _type(type),_flag(flag),_chisqprob(chisqprob),_traj(traj)
{}

TrkFitSummary::TrkFitSummary(const TrkFitSummary& other) :
  _type(other._type),
  _flag(other._flag),
  _chisqprob(other._chisqprob),
  _traj(other._traj)
{}

TrkFitSummary&
TrkFitSummary::operator =(const TrkFitSummary& other) {
  if(this != &other){
    _type = other._type;
    _flag = other._flag;
    _chisqprob = other._chisqprob;
    _traj = other._traj;
  }
  return *this;
}

TrkFitSummary::~TrkFitSummary()
{}


bool
TrkFitSummary::operator == (const TrkFitSummary& other) const {
  return _type == other._type &&
    _flag == other._flag &&
    _traj == other._traj;
}


// ------------------------------------------------------------------------------
// File and Version Information:
//      $Id: DetExternalSet.cc,v 1.7 2002/12/17 04:03:57 dbrown Exp $
//
// Trivial implementations of intersection methods for a set of elements
// outside the tracking volume
//
//
// Copyright Information:
//	Copyright (C) 1997	Lawrence Berkeley Laboratory
//
// Author : Gautier Hamel de Monchenault
//
//------------------------------------------------------------------------------
#include "BaBar/BaBar.hh"
#include "ErrLogger/ErrLog.hh"
#include "DetectorModel/DetIntersection.hh"
#include "DetectorModel/DetExternalSet.hh"
#include <vector>

DetExternalSet::DetExternalSet(const char* name,int IDNum) :
  DetSet(name,IDNum)
{
}

void
DetExternalSet::intersection(std::vector<DetIntersection>& divec,
			     const Trajectory* traj,
			     double* myrange, bool clear) const 
{
  if(clear)divec.clear();
}

bool
DetExternalSet::nextIntersection(const Trajectory* traj,
			      DetIntersection& next,
			      DetIntersection* lastinter) const
{
  return false;
}


bool
DetExternalSet::firstIntersection(const Trajectory* traj,
				  DetIntersection& next,
				  double* myrange) const 
{
  return false;
}



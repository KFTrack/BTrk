// ------------------------------------------------------------------------------
// File and Version Information:
//      $Id: DetIntersection.cc,v 1.4 2000/05/19 23:56:58 brownd Exp $
//
//  Description:
//  Structure for intersection return values.  This has to obey all the rules
//  for ValSortedVector
//
// Copyright Information:
//	Copyright (C) 1996	Lawrence Berkeley Laboratory
//
//  Authors: Dave Brown, 1/8/97
//------------------------------------------------------------------------------
#include "BaBar/BaBar.hh"
#include "ErrLogger/ErrLog.hh"
#include "DetectorModel/DetIntersection.hh"
//
//  Constructors
//
DetIntersection::DetIntersection() :
  delem(0),pathlen(0.0),trajet(0)
{
  pathrange[0] = pathrange[1] = 0.0; 
  flag[0] = flag[1] = 0;
}
DetIntersection::DetIntersection(const DetIntersection& other) :
  delem(other.delem),pathlen(other.pathlen),trajet(other.trajet)
{
  pathrange[0] = other.pathrange[0];
  pathrange[1] = other.pathrange[1];
  flag[0]      = other.flag[0];
  flag[1]      = other.flag[1];
}

DetIntersection::DetIntersection(const DetElem* elem,const Trajectory* traj,double len):
  delem(elem),pathlen(len),trajet(traj)
{
  pathrange[0] = pathrange[1] = len;
  flag[0] = flag[1] = 0;
}

DetIntersection::DetIntersection(const DetElem* elem,const Trajectory* traj,double range[2]):
  delem(elem),trajet(traj)
{
  pathlen = pathrange[0] = range[0];
  pathrange[1] = range[1];
  flag[0] = flag[1] = 0;
}

DetIntersection::DetIntersection(const DetElem* elem,const Trajectory* traj,
					   double lowrange,double highrange):
  delem(elem),pathlen(lowrange),trajet(traj)
{
  pathrange[0] = lowrange;
  pathrange[1] = highrange;
  flag[0] = flag[1] = 0;
}

DetIntersection::DetIntersection(const DetElem* elem,const Trajectory* traj,double path,
					   double lowrange,double highrange):
  delem(elem),pathlen(path),trajet(traj)
{
  pathrange[0] = lowrange;
  pathrange[1] = highrange;
  flag[0] = flag[1] = 0;
}

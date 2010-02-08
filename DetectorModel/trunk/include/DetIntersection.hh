// ------------------------------------------------------------------------------
// File and Version Information:
//      $Id: DetIntersection.hh,v 1.7 2001/09/06 17:20:21 steinke Exp $
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
#ifndef DETECTORINTERSECTION_HH
#define DETECTORINTERSECTION_HH
#define EPSILON 1.0e-3 // test value for intersections
#include "DetectorModel/DetElem.hh"
#include "BbrGeom/Trajectory.hh"
#include "CLHEP/config/CLHEP.h"
//
struct DetIntersection{
  const DetElem* delem; // the concerned detector element
  const Trajectory* trajet; // the concerned trajectory
  double pathlen; // path length at intersection (typically the average of the range)
  double pathrange[2]; // path length at entrance and exit (can be equal, = pathlen)
  int flag[2]; // flag describing entrance and exit.
  DetIntersection();
  DetIntersection(const DetIntersection&);
  DetIntersection(const DetElem*,const Trajectory*,double);
  DetIntersection(const DetElem*,const Trajectory*,double[2]);
  DetIntersection(const DetElem*,const Trajectory*,double,double);
  DetIntersection(const DetElem*,const Trajectory*,double,double,double);
  ~DetIntersection(){;}
  bool operator == (const DetIntersection& other) const {
    return (delem == other.delem) && (trajet == other.trajet) &&
      fabs(pathlen - other.pathlen) < EPSILON; }
//  Comparison use only path length, to enable sorts of heterogeneous lists
  int operator < (const DetIntersection& other) const {
    return pathlen < other.pathlen; }
  int operator > (const DetIntersection& other) const {
    return pathlen > other.pathlen; }
  DetIntersection& operator = (const DetIntersection& other){
    delem = other.delem;
    trajet = other.trajet;
    pathlen = other.pathlen;
    pathrange[0] = other.pathrange[0];
    pathrange[1] = other.pathrange[1];
    flag[0] = other.flag[0];
    flag[1] = other.flag[1];
    return *this;
  }
  double pathLength() const { return pathrange[1]-pathrange[0]; }
};
#endif

// ------------------------------------------------------------------------------
// File and Version Information:
// $Id: Intersection.hh,v 1.8 2004/10/11 19:08:05 brownd Exp $
//
//
// Intersection is a special class that handles the association of a Ray and a
// Surface object to find the intersection if any.
//
// Both Trajectory and Surface are polymorphic
//
// Coded stolen from gismo 2 
//    and modified for BaBar by Gautier Hamel de Monchenault
//
// Copyright Information:
//	Copyright (C) 1997		Lawrence Berkeley Laboratory
//	Copyright (C) 1997	       CEA - Centre d'Etude de Saclay
//
//------------------------------------------------------------------------
#ifndef INTERSECTION_HH
#define INTERSECTION_HH

#include "BTrk/TrkBase/TrkErrCode.hh"
#include "BTrk/BbrGeom/Trajectory.hh"
#include "BTrk/TrkBase/TrkDirection.hh"
#include "BTrk/DetectorModel/DetSurface.hh"

class SurfacePoint;

class Intersection
{
 public:

  Intersection(const Trajectory& t, const DetSurface & s )
  : traj(t), surf(s)
  {}
  
  TrkErrCode intersect( double& f, SurfacePoint& uv,
			trkDirection tdir = trkOut,double* frange = 0 ) const ;
// same, without surfave point
  TrkErrCode intersect( double& f,
			trkDirection tdir = trkOut,double* frange = 0 ) const ;
// old code, for verification
  TrkErrCode oldintersect( double& f, SurfacePoint& uv,
			   trkDirection tdir = trkOut,double* frange = 0 ) const ;
  
  // return the actual intersection f, with a condition flag

  void funcd(double d, double& f, double& df) const;
  // for distance along the trajectory d, sets f to distance to the surface,
  // and df to its derivate.
  TrkErrCode root(const double& x1, const double& x2, double& f) const;
  // return special Newton-Raphson solution designed for this case
 private:
// find approximate intersection using nonlinear techniques
  DetSurface::intertype findProximity(double&,double,double,trkDirection) const;
  const Trajectory& traj;
  const DetSurface& surf;
};
#endif

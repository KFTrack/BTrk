//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: TrkCentralVolume.hh 105 2010-01-15 12:14:11Z stroili $
//
// Description:
//	TrkCentralVolume Class - Implements a TrkVolume as a single cylinder
//      centered on the Z axis.  The inner flight distance is determined as
//      the distance of closest approach to a straight line trajectory on
//      the Z axis
//
// Author List:
//      David Brown, LBL 5/13/97
//
//      
//
// Copyright Information:
//	Copyright (C) 1997		Lawrence Berkeley Laboratory
//------------------------------------------------------------------------

#ifndef TRKCENTRALVOLUME_HH
#define TRKCENTRALVOLUME_HH

#include "TrkBase/TrkVolume.hh"
#include "DetectorModel/DetCylinder.hh"
#include "TrajGeom/TrkLineTraj.hh"


class TrkCentralVolume : public TrkVolume {

public:
// Constructors
  TrkCentralVolume(const char* name,double rmax,double zmin, double zmax);
// Destructor
  virtual ~TrkCentralVolume();
// Extend the flight range of the trajectory through the volume
  virtual bool extendThrough( const Trajectory* theTraj, 
			      double& theFlightDist,
			      trkDirection theDirection=trkOut,
			      double* theStartingFlightDist=0) const;
// returns false if the point is outside the volume 
  virtual bool isInside( const HepPoint& ) const;

private:
  DetCylinder _outercyl; // outer cylinder defining the volume
  TrkLineTraj _zaxis; // line traj on Z axis
  double _zmin;
  double _zmax; // Z limits
};


#endif


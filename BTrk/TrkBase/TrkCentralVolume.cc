//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: TrkCentralVolume.cc 105 2010-01-15 12:14:11Z stroili $
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

#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/BaBar/ErrLog.hh"
#include "BTrk/TrkBase/TrkCentralVolume.hh"
#include <assert.h>
#include "BTrk/BaBar/Constants.hh"
#include "CLHEP/Vector/ThreeVector.h"
#include "BTrk/BbrGeom/Transformation.h"
#include "BTrk/BbrGeom/HepPoint.h"
#include "BTrk/BbrGeom/Trajectory.hh"
#include "BTrk/TrkBase/TrkErrCode.hh"
#include "BTrk/TrkBase/TrkPoca.hh"
#include "BTrk/DetectorModel/Intersection.hh"

TrkCentralVolume::TrkCentralVolume(const char* name,double rmax,double zmin,double zmax) :
  TrkVolume(name),
  _outercyl(HepTransformation(),rmax),_zaxis(HepPoint(0,0,zmin),HepPoint(0,0,zmax)),
  _zmin(zmin),_zmax(zmax)
{;}
//
//  Destructor
//
TrkCentralVolume::~TrkCentralVolume()
{;}

//
// returns true if the onput point lies inside the tracking volume
//
bool
TrkCentralVolume::isInside( const HepPoint& point )  const 
{
  const double r ( sqrt( sqr(point.x()) + sqr(point.y()) ) );
  const double z ( point.z() );
  return ( r <= _outercyl.radius()  && z >= _zmin && z <= _zmax );
}
//
//  Extension function
//
bool 
TrkCentralVolume::extendThrough( const Trajectory* traj, 
				 double& flight, 
				 trkDirection trkdir,
				 double* startflight) const 
{
//
//  Switch on direction
//
  TrkErrCode iflag;
  SurfacePoint uv;
  double f;
//
//  Set the starting flight length according to the specified initial value
//
  if(startflight != 0)
    f = *startflight;
  else
    f = traj->lowRange();
  bool hitit = false;
//
//  find POCA to the line traj if we're going inward.  This will be the
//  starting point for an intersection test agains the cylinder.
  if(trkdir == trkIn){
    double zlen = -_zmin;
    TrkPoca zpoca(*traj,f,_zaxis,zlen,1.0e-3);
    iflag = zpoca.status();
    f = zpoca.flt1();
    zlen = zpoca.flt2();
    if(iflag.success()){
      flight = f;
      hitit = true;
    }
  }
// determine a reasonable range to extend the flight distance :
// The range increment must be smaller than the length of a full loop,
// assuming that the trajectory is described locally by an helix
  HepPoint point( traj->position( f ) );
  Hep3Vector vect( traj->direction( f ) );
  double cosDip( vect.perp()/vect.mag() );
  Hep3Vector delDir( traj->delDirect( f ) );
  double curv( delDir.mag() );
  double maxDeltaF = curv > 0 ? 2*Constants::pi*cosDip/curv: 500;
//
// set the range for intersections
//
  double frange[2];
  Intersection inter( *traj,  _outercyl);
  if(trkdir == trkOut) {
    frange[0] = f;
    frange[1] = f+maxDeltaF;
  } else {
    frange[0] = f-maxDeltaF;
    frange[1] = f;
  }
//
// try intersecting with the outer cylinder
//    
  iflag = inter.intersect(f,uv,trkdir,frange);
  if( iflag.success() && uv[1] >= _zmin &&  uv[1] <= _zmax ){
    flight = f;
    hitit = true;
  }
  return hitit;
}

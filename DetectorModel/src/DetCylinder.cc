//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DetCylinder.cc,v 1.20 2004/08/06 05:58:29 bartoldu Exp $
//
// Description:
//	DetCylinder Class 
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      David N Brown                - Lawrence Berkeley Lab
//	Gautier Hamel de Monchenault - CEN Saclay & Lawrence Berkeley Lab
//
// Copyright Information:
//	Copyright (C) 1996	Lawrence Berkeley Laboratory
//	Copyright (C) 1996	CEA - Centre d'Etude de Saclay
//
//------------------------------------------------------------------------

//----------------
// BaBar Header --
//----------------
#include "BaBar/BaBar.hh"

//-----------------------
// This Class's Header --
//-----------------------
#include "DetectorModel/DetCylinder.hh"

//---------------
// C++ Headers --
//---------------
#include <assert.h>
#include <iostream>

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "CLHEP/Geometry/HepPoint.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Geometry/Transformation.h"
#include "CLHEP/Utilities/CLHEP.h"
using std::endl;
using std::ostream;
//-----------------------------------------------------------------------
// Local Macros, Typedefs, Structures, Unions and Forward Declarations --
//-----------------------------------------------------------------------
//----------------
// Constructors --
//----------------
//
DetCylinder::DetCylinder(const HepTransformation& transf, double radius)
  :   DetSurface(transf),_r(radius)
{;}
//
// copy constructor
DetCylinder::DetCylinder( const DetCylinder& c ) 
  :  DetSurface( c ), _r( c._r )
{;}
//--------------
// Destructor --
//--------------
DetCylinder::~DetCylinder( ) {
}
//
//  copyOf function
//
DetSurface*
DetCylinder::copyOf(double fparam) const{
  if(fparam == 0.0)
    return (DetSurface*) new DetCylinder(*this);
  else {
//
//  Family is defined by changing the radius
//
    assert(_r+fparam > 0.0); // only positive radii
    return (DetSurface*) new DetCylinder(*transform(),_r+fparam);
  }
}
// 
// Member functions
//
HepPoint
DetCylinder::spacePoint(const SurfacePoint& uv) const {
  assert ( _r!=0. );
  double x = _r*cos( uv[0] );
  double y = _r*sin( uv[0] );
  double z = uv[1];
  return gotoGlobal( HepPoint( x, y, z ) );
}

Hep3Vector
DetCylinder::normal(const SurfacePoint& uv) const {
  double x = cos( uv[0] );
  double y = sin( uv[0] );
  return gotoGlobal( Hep3Vector( x, y, 0.) );
}

Hep3Vector
DetCylinder::surfaceDirection(const SurfacePoint& uv,int idir) const {
  switch(idir){
  case 1:
    return axis();
  case 0: default:
    double x = -sin( uv[0] );
    double y = cos( uv[0] );
    return gotoGlobal( Hep3Vector( x, y, 0.) );
  }    
}

int 
DetCylinder::surfacePoint( const HepPoint& x, SurfacePoint& uv, double tol )const  {
  HepPoint p( gotoLocal(x) );
  uv[0] = p.phi();
  uv[1] = p.z();
  double r = sqrt( p.x()*p.x() + p.y()*p.y() );
  if( fabs( r - _r ) > tol ) return 1;
  return 0;
}

double
DetCylinder::normTo( const HepPoint& x, Hep3Vector& n ) const {
  Hep3Vector r = x - centerPoint();
  double    rl = r*axis();  
  Hep3Vector radial = ( r - rl*axis() );
  double rt = radial.mag();
//  Protection against points on the axis
  n = rt > 0 ? radial*(1.0/rt) : basis(0);
  return _r-rt;    // signed quantity !
}

double
DetCylinder::firstDeriv( const HepPoint& x, const Hep3Vector& v ) const {
  Hep3Vector r = x - centerPoint();
  double    rl = r*axis();  
  Hep3Vector n = ( r - rl*axis() ).unit();
  return v.unit()*n;
}

DetSurface::intertype
DetCylinder::distTo( const HepPoint& p, const Hep3Vector& d, 
			 double& dist,intermode mode) const {
// normalize the direction
  Hep3Vector dir(d.unit());
// simple vectors
  Hep3Vector axe(axis());
  Hep3Vector diff(centerPoint()-p);
// compute cross-product of directions; this is radial to POCA
  Hep3Vector radial(dir.cross(axe));
// check that the directions aren't parallel, in which case there's no
// crossing
  double cmag = radial.mag();
  DetSurface::intertype retval = DetSurface::nointersect;
  if(cmag > 0.0){
// compute transverse distance
    radial *= (1.0/cmag);
// compute the other unit vector, perp to POCA
    Hep3Vector raddir(axe.cross(radial));
// projections
    double dtrans = fabs(diff * radial);
    double dlong = diff * raddir;
// projection factor to get distance along the original direction
    double cost = axe*dir;
    double proj = 1.0/sqrt(1.0-cost*cost);
// check to see if there's a real intersection
    if(dtrans <= _r){
// 2 solutions; choose the one specified by the mode
      double delta = sqrt(_r*_r - dtrans*dtrans);
      double long1=dlong-delta;
      double long2=dlong+delta;
      retval = DetSurface::intersect;
      switch (mode) {
      case closest:
	dlong = fabs(long1)<fabs(long2) ? long1 : long2;
	break;
      case forward:
	dlong = long1>0 ? long1 : long2;
	if(dlong<0.0) retval = DetSurface::outofrange;
	break;
      case backward:
	dlong  = long2<0 ? long2 : long1;
	if(dlong>0.0) retval = DetSurface::outofrange;
	break;
      }
// correct for the projection
      dist =dlong*proj;
      return retval;
    } else {
// no intersection; return POCA
      retval = DetSurface::localmin;
      dist =dlong*proj;
      if( (mode==forward && dist<0.0) ||
	  (mode==backward && dist>0.0) )
	retval = DetSurface::outofrange;
    }      
  }
  return retval;
}
//
double
DetCylinder::curvature(const SurfacePoint& uv) const {
  return 1/_r;
}
//
bool
DetCylinder::wrappedCoordinate(int icoord) const {
  switch(icoord){
  case 0:
    return true;
  case 1: default:
    return false;
  }
}
//
void
DetCylinder::segmentMinMax(const HepPoint& p1,const HepPoint& p2,
			   double& mindist,double& maxdist) const {
// find the radii from the cylinder axis to the points
  Hep3Vector radial1 = axis().cross(p1 - centerPoint());
  Hep3Vector radial2 = axis().cross(p2 - centerPoint());
  double mag1 = radial1.mag();
  double mag2 = radial2.mag();
// maximum (signed!) distance is given by the maximum radius
  maxdist = max(mag1,mag2) - radius();
// minimum distance can be at a point or at the segment tangent
  if(mag1>0.0 && mag2>0.0){
    double theta = radial1.angle(radial2);
    double cost = cos(theta);
    if(mag1/mag2 <= cost || mag2/mag1 <= cost)
      mindist = min(mag1,mag2) - radius();
    else {
      double dist = (radial1-radial2).mag();
      if(dist>0.0)
	mindist = mag1*mag2*sin(theta)/dist - radius();
      else
	mindist = mag1-radius();
    }
  } else
    mindist = -radius();
}

void
DetCylinder::print(ostream& os) const {
  os << " Cylinder Surface with radius " << _r << endl;
}

void
DetCylinder::printAll(ostream& os) const {
  HepPoint origin = transform()->origin();
  Hep3Vector axis = transform()->unit();
  Hep3Vector xaxis = transform()->unit(0);
  os << " Cylinder Surface with radius " << _r << endl;
  os << "Origin = " << origin.x() << " " << origin.y() << " " << origin.z() << endl
     << "Axis direction = " 
     << axis.x() << " " << axis.y() << " " << axis.z() << endl
     << "Local X direction = " 
     << xaxis.x() << " " << xaxis.y() << " " << xaxis.z() << endl;
}

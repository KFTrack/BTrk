//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DetCone.cc,v 1.17 2004/08/06 05:58:29 bartoldu Exp $
//
// Description:
//	DetCone Class 
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      David N Brown                - Lawrence Berkeley Lab
//	Gautier Hamel de Monchenault - CEN Saclay & Lawrence Berkeley Lab
//      Phil Strother                - Imperial College
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
#include "DetectorModel/DetCone.hh"

//---------------
// C++ Headers --
//---------------
#include <assert.h>
#include <iostream>

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "BaBar/Constants.hh"
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
DetCone::DetCone(const HepTransformation& transf, double tanTheta, double vertex)
  :   DetSurface(transf), _tanTheta(tanTheta),_z0(vertex)
{
//  tanTheta = fabs(tanTheta);
  // Constrain theta to lie between 0 and pi/2 - if theta>pi/2 then theta==>pi-theta

  assert(tanTheta!=0); 
  // Make sure we haven't got a plane
  
  // assign the quadrant by hand; we wan the sine to be always positive!
  _cosTheta = copysign(1.0/sqrt(1.0 + tanTheta*tanTheta),tanTheta);
  _sinTheta = tanTheta*_cosTheta;
}

DetCone::DetCone(const HepTransformation& transf, double z1, double z2, double r1, double r2)
  :   DetSurface(transf)
{
  assert(z1!=z2);
  assert(r1!=r2);
  // Make sure we haven't got a plane or a cylinder

  r1=fabs(r1);
  r2=fabs(r2);
  
  _tanTheta = (r2-r1)/(z2-z1);

  _z0 = z1 - r1/_tanTheta;
//  _tanTheta = fabs(_tanTheta);
  double theta = atan2(r2-r1,z2-z1);
  _cosTheta = cos(theta);
  _sinTheta = sin(theta);

 
}
  
//
// copy constructor
DetCone::DetCone( const DetCone& c ) 
  : DetSurface( c ),
    _tanTheta( c._tanTheta ),
    _z0(c._z0),
    _cosTheta( c._cosTheta ),
    _sinTheta( c._sinTheta )
{
}

//--------------
// Destructor --
//--------------

DetCone::~DetCone( ) {
}

//
//  copyOf function
//
DetSurface* DetCone::copyOf(double fparam) const{
  if (fparam==0.0){
    return (DetSurface*) new DetCone(*this);
  } else {
    // Assume that a family of cones means a bunch of cones with the same angle but
    // different vertices.
    return (DetSurface*) new DetCone(*transform(),_tanTheta,
				     _z0-fparam/_sinTheta);
  }
}
// 
// Member functions
//

HepPoint DetCone::spacePoint(const SurfacePoint& uv) const {

  double x = uv[1]*_sinTheta*cos( uv[0] );
  double y = uv[1]*_sinTheta*sin( uv[0] );
  double z = uv[1]*_cosTheta + _z0;
//  int signOfR = (int)copysign( 1, uv[1] );
  int signOfR = (uv[1] < 0.0) ? -1 : 1;
  return gotoGlobal( HepPoint( signOfR*x, signOfR*y, z ) );
}

Hep3Vector DetCone::normal(const SurfacePoint& uv) const {
  // Defined as the outward normal to the surface 
  double x = _cosTheta * cos( uv[0] );
  double y = _cosTheta *sin( uv[0] );
  double z = -_sinTheta;
//  int signOfNormal = (int)copysign( 1, uv[1] );
  int signOfNormal = (uv[1] < 0.0) ? -1 : 1;
  return gotoGlobal(   Hep3Vector(x,y,signOfNormal*z));
}

Hep3Vector DetCone::surfaceDirection(const SurfacePoint& uv,int idir) const {
  double x,y,z;
  int signOfR;
  switch(idir){
  case 1:
    x = uv[1]*_sinTheta*cos(uv[0]);
    y = uv[1]*_sinTheta*sin(uv[0]);
    z = uv[1]*_cosTheta;
//    signOfR = (int)copysign( 1, uv[1] );
    signOfR = (uv[1] < 0.0) ? -1 : 1;
    return gotoGlobal(Hep3Vector(signOfR*x,signOfR*y,z));
  case 0: default:
    x = -uv[1]*_sinTheta*sin( uv[0] );
    y = uv[1]*_sinTheta*cos( uv[0] );
//    signOfR = (int)copysign( 1, uv[1] );
    signOfR = (uv[1] < 0.0) ? -1 : 1;
    return gotoGlobal( Hep3Vector( signOfR*x, signOfR*y, 0.) );
    // Unit vector in the direction of +ve phi
  }    
}

int DetCone::surfacePoint( const HepPoint& thePoint, SurfacePoint& uv, double tolerance )const  {
  HepPoint theLocalPoint( gotoLocal(thePoint) );
  uv[0] = theLocalPoint.phi();
  if (uv[0] < 0) uv[0]+=Constants::twoPi;
  // Ensure phi lies between 0 and 2 pi

  uv[1] = (theLocalPoint.z() - _z0)/_cosTheta;
  double radius = sqrt( theLocalPoint.x()*theLocalPoint.x() + theLocalPoint.y()*theLocalPoint.y() );
  double actualRadius = fabs(uv[1]) *_sinTheta;
  if( fabs( radius - actualRadius ) > tolerance ) return 1;
  return 0;
}

double DetCone::normTo( const HepPoint& thePoint, Hep3Vector& theNormal) const {
  Hep3Vector localPositionVector = thePoint - centerPoint();

  // The plane with normal along the surface at the phi given by the requested point
  // contains the point in question and its nearest point on the cone
//  int isign = (int)copysign( 1, localPositionVector.z()-_z0);
  int isign = (localPositionVector.z()-_z0 < 0.0) ? -1 : 1;
  // This sign needed since theta lies between 0 and pi/2
  
  Hep3Vector normalToPlane(_sinTheta * cos(localPositionVector.phi()), 
			   _sinTheta * sin(localPositionVector.phi()), 
			   isign*_cosTheta);


  double thePlanesDefiningDistance(localPositionVector*normalToPlane);
  // No apostrophes in C++......
  double theClosestPointsZCoord(((thePlanesDefiningDistance - _z0*isign*_cosTheta)*normalToPlane).z()+_z0);  
  double radiusAtClosestZ(fabs((theClosestPointsZCoord - _z0) * _tanTheta));
  double radiusOfGivenPoint(sqrt(localPositionVector.x()*localPositionVector.x() + 
				 localPositionVector.y()*localPositionVector.y()));
  double difference(radiusAtClosestZ/_cosTheta - radiusOfGivenPoint/_cosTheta  );
  // Not sure about the sign here.  I think all it is is a relative one i.e. that the
  // sign changes inside to outside.  Have to check aginst Intersection::intersect();
  double x(((thePlanesDefiningDistance - _z0*isign*_cosTheta)*normalToPlane).x());
  double y(((thePlanesDefiningDistance - _z0*isign*_cosTheta)*normalToPlane).y());
  double z(((thePlanesDefiningDistance - _z0*isign*_cosTheta)*normalToPlane).z()+_z0);
  HepPoint localpoint(x,y,z);
  SurfacePoint uv;
  int ierr = surfacePoint( gotoGlobal(localpoint), uv );
  assert( ierr == 0);
  theNormal = normal(uv);
  return difference;
}

double DetCone::firstDeriv( const HepPoint& thePoint, const Hep3Vector& theDirection ) const {
  Hep3Vector theNormal;
  normTo(thePoint, theNormal);
  return theDirection.unit()*theNormal;
}

DetSurface::intertype
DetCone::distTo( const HepPoint& thePoint, const Hep3Vector& theDirection,
		 double& dist,intermode mode) const {
  HepPoint   localPoint = gotoLocal(thePoint);
  Hep3Vector localDirection = (gotoLocal(theDirection)).unit();
// Store some useful quantities that get used a lot. 
  // Hep3Vector inlines thes functions, so we don't gain much speed but it makes the code 
  // look slightly less of a complete mess than it would otherwise
  double lpx(localPoint.x()),lpy(localPoint.y()),lpz(localPoint.z());
  double ldx(localDirection.x()),ldy(localDirection.y()),ldz(localDirection.z());
  
  // Now set up and solve our quadratic

  double a = ldx*ldx + ldy*ldy - ldz*ldz*_tanTheta*_tanTheta;

  double b = lpx*ldx + lpy*ldy - ldz*(lpz-_z0)*_tanTheta*_tanTheta;

  double c = lpx*lpx + lpy*lpy - (lpz - _z0)*(lpz - _z0)*(_tanTheta*_tanTheta);

  double disc = b*b - a*c;

  if( disc < 0. ) return DetSurface::nointersect;
  DetSurface::intertype retval = DetSurface::intersect;
  if( a==0. ) { 
    dist = -0.5*c/b;
    if( (mode==forward && dist<0.0) ||
	(mode==backward && dist>0.0) )retval = DetSurface::outofrange;
  } else {
// multiple roots; take the one specified by the mode
    double delta=sqrt(disc);
    double r1 =  a>0.0 ? (-b+delta)/a : (-b-delta)/a;
    double r2 =  a>0.0 ? (-b-delta)/a : (-b+delta)/a;;
    switch (mode) {
    case closest:
      dist = ( fabs(r1) < fabs(r2) ) ? r1:r2;
      break;
    case forward:
      dist = r2 > 0 ? r2 : r1;
      if(dist<0.0) retval = DetSurface::outofrange;
      break;
    case backward:
      dist = r1 < 0 ? r1 : r2;
      if(dist>0.0) retval = DetSurface::outofrange;
      break;
    }
  }
  return retval;
}


double DetCone::curvature(const SurfacePoint& uv) const {
  double radius = uv[1]*_sinTheta;
  return 1.0/radius;
}

bool DetCone::wrappedCoordinate(int theCoordNumber) const{
  switch(theCoordNumber){
  case 1: 
    return false;
  case 0: default:
    return true;
  }
}

void DetCone::segmentMinMax(const HepPoint &pointOne, const HepPoint &pointTwo, 
			    double &mindist, double &maxdist) const{

  Hep3Vector firstVect, secondVect;

  double dist1 = normTo(pointOne, firstVect);
  double dist2 = normTo(pointTwo, secondVect);

  mindist = min(dist1, dist2);
  maxdist = max(dist1,dist2);

}

void
DetCone::print(ostream& os) const {
  os << " Cone Surface with opening angle = " << atan(_tanTheta) << endl;
}

void
DetCone::printAll(ostream& os) const {
  HepPoint origin = transform()->origin();
  Hep3Vector axis = transform()->unit();
  Hep3Vector xaxis = transform()->unit(0);
  os << " Cone Surface with opening angle = " << atan2(_sinTheta,_cosTheta) << endl;
  os << "Origin = " << origin.x() << " " << origin.y() << " " << origin.z() << endl
     << "Axis direction = " 
     << axis.x() << " " << axis.y() << " " << axis.z() << endl
     << "Local X direction = " 
     << xaxis.x() << " " << xaxis.y() << " " << xaxis.z() << endl;
}

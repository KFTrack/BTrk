//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DetPlane.cc,v 1.14 2004/08/06 05:58:31 bartoldu Exp $
//
// Description:
//	DetPlane Class 
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
#include "DetectorModel/DetPlane.hh"

//---------------
// C++ Headers --
//---------------
#include <assert.h>
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
DetPlane::DetPlane(const HepTransformation& transf)
  : DetSurface(transf)
{;}
//
// copy constructor
DetPlane::DetPlane( const DetPlane& p ) 
  :  DetSurface( (DetSurface&)p )
{;}
//--------------
// Destructor --
//--------------
DetPlane::~DetPlane( ) {
}
//
//  copyOf function
//
DetSurface*
DetPlane::copyOf(double fparam) const {
//
//  Make sure the simple copy is identical
//
  if(fparam == 0.0)
    return new DetPlane(*this);
  else {
//
//  Build a new transform based on extending the translation
//
    Hep3Vector transvec = transform()->tra_vec().trans_vec();
    HepTranslation newtranslate(transvec+fparam*axis());
    HepTransformation newtransform(transform()->rot_mat(),newtranslate);
    DetPlane* newDetPlane = new DetPlane(newtransform);
    return (DetSurface*) newDetPlane;
  }
}
// 
// Member functions
//
HepPoint
DetPlane::spacePoint(const SurfacePoint& uv) const {
  return gotoGlobal( HepPoint( uv[0], uv[1], 0. ) );
}

Hep3Vector
DetPlane::normal(const SurfacePoint& ) const {
  return axis();
}

Hep3Vector
DetPlane::surfaceDirection(const SurfacePoint& ,int idir) const {
  switch(idir){
  case 1:
    return gotoGlobal( Hep3Vector( 0.0, 1.0, 0. ) );
  case 0: default:
    return gotoGlobal( Hep3Vector( 1.0, 0.0, 0. ) );
  }    
}

int 
DetPlane::surfacePoint( const HepPoint& x, SurfacePoint& uv, double tol ) const  {
  HepPoint p( gotoLocal(x) );
  uv[0] = p.x();
  uv[1] = p.y();
  if( fabs( p.z() ) > tol ) return 1;
  return 0;
}

double
DetPlane::normTo( const HepPoint& x, Hep3Vector& n ) const {
  n = axis();  
  double d = (centerPoint()-x)*n;    // signed quantity !
  return d;
}

double
DetPlane::firstDeriv( const HepPoint& x, const Hep3Vector& v ) const {
  return v.unit()*axis();
}

DetSurface::intertype
DetPlane::distTo( const HepPoint& x, const Hep3Vector& v,
		  double& dist,intermode mode) const {
  double d = (x-centerPoint())*axis();
  double s = -(v.unit())*axis();
  if( s != 0.0){
    dist = d/s;
    if( mode == closest ||
	(mode==forward && dist>0.0) ||
	(mode==backward && dist<0.0) )
      return DetSurface::intersect;
    else
      return DetSurface::outofrange;
  } else
    return DetSurface::nointersect;  // direction parallel to the DetPlane !
}

void
DetPlane::segmentMinMax(const HepPoint& p1,const HepPoint& p2,
			double& mindist,double& maxdist) const{

  Hep3Vector point1 = p1 - centerPoint();
  Hep3Vector point2 = p2 - centerPoint();
  double dist1 = point1*axis();
  double dist2 = point2*axis();
  mindist = min(dist1,dist2);
  maxdist = max(dist1,dist2);
}

void
DetPlane::print(ostream& os) const {
  os << " Planar Surface " << endl;
}

void
DetPlane::printAll(ostream& os) const {
  HepPoint origin = transform()->origin();
  Hep3Vector axis = transform()->unit();
  Hep3Vector xaxis = transform()->unit(0);
  os << " Planar Surface with origin = " 
     << origin.x() << " " << origin.y() << " " << origin.z() << endl
     << "Normal direction = " 
     << axis.x() << " " << axis.y() << " " << axis.z() << endl
     << "Local X direction = " 
     << xaxis.x() << " " << xaxis.y() << " " << xaxis.z() << endl;
}

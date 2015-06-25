//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DetSurface.cc,v 1.10 2004/08/06 05:58:32 bartoldu Exp $
//
// Description:
//	DetSurface virtual base Class 
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
#include "BTrk/BaBar/BaBar.hh"

//-----------------------
// This Class's Header --
//-----------------------
#include "BTrk/DetectorModel/DetSurface.hh"

//---------------
// C++ Headers --
//---------------
#include <assert.h>


#include "BTrk/BbrGeom/HepPoint.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "BTrk/BbrGeom/Transformation.h"
using std::endl;
using std::ostream;

//-----------------------------------------------------------------------
// Local Macros, Typedefs, Structures, Unions and Forward Declarations --
//-----------------------------------------------------------------------

//----------------
// Constructors --
//----------------
//
DetSurface::DetSurface(const HepTransformation& transf) :
  _transf( new HepTransformation( transf ) )
{ 
  assert ( _transf != 0 );
}
//
// copy constructor
DetSurface::DetSurface( const DetSurface& p ) 
  :  _transf( new HepTransformation( *p._transf ) )
{ 
  assert ( _transf != 0 );
}

//--------------
// Destructor --
//--------------
DetSurface::~DetSurface( ) {
  if( _transf ) delete _transf;
}

HepPoint   
DetSurface::centerPoint() const {
  return transform()->origin();
}

Hep3Vector 
DetSurface::basis(int ax) const {
  return transform()->unit(ax);
}

Hep3Vector 
DetSurface::axis() const {
  return basis(); 
}

HepPoint   
DetSurface::gotoLocal ( const HepPoint& x )   const {
  return transform()->transTo(x) ; 
}

HepPoint   
DetSurface::gotoGlobal( const HepPoint& x )   const {
  return transform()->transFrom(x) ; 
}

Hep3Vector 
DetSurface::gotoLocal ( const Hep3Vector& p ) const {
  return transform()->transTo(p) ; 
}

Hep3Vector 
DetSurface::gotoGlobal( const Hep3Vector& p ) const {
  return transform()->transFrom(p) ; 
}

// default implementation of print
void
DetSurface::print(ostream& os) const {
  os << "Detector Surface" << endl; }
//
void
DetSurface::printAll(ostream& os) const {
  HepPoint origin = transform()->origin();
  Hep3Vector axis = transform()->unit();
  Hep3Vector xaxis = transform()->unit(0);
  os << "Detector Surface with transform origin = "
     << origin.x() << " " << origin.y() << " " << origin.z() << endl
     << "localZ direction = " 
     << axis.x() << " " << axis.y() << " " << axis.z() << endl
     << "Local X direction = " 
     << xaxis.x() << " " << xaxis.y() << " " << xaxis.z() << endl;
}

double
DetSurface::normalTo( const HepPoint& point,Hep3Vector& norm,
		      SurfacePoint& uv) const {
  double dist = normTo(point,norm);
  HepPoint onsurf = point+dist*norm;
  int isonsurf = surfacePoint(onsurf,uv);
  assert(isonsurf == 0);
  return dist;
}

DetSurface::intertype
DetSurface::distanceTo( const HepPoint& point,const Hep3Vector& dir,
			double& dist, SurfacePoint& uv,
			intermode mode) const {
  intertype retval = distTo(point,dir,dist,mode);
  if(retval == intersect){
    Hep3Vector direction = dir.unit();
    HepPoint onsurf = point+dist*direction;
    int isonsurf = surfacePoint(onsurf,uv);
    assert(isonsurf == 0);
  }
  return retval;
}

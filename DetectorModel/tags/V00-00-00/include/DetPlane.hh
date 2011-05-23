//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DetPlane.hh,v 1.10 2004/08/06 05:58:31 bartoldu Exp $
//
// Description:
//	Description of a DetPlane for BaBar reconstruction
//      A DetPlane is a Surface defined by equation z=0 in its 
//      local coordinate system.
//      The 2D system of local coordinates is u=x and v=y. 
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      David N Brown                - Lawrence Berkeley Lab
//	Gautier Hamel de Monchenault - CEN Saclay & Lawrence Berkeley Lab
//
// History (add to end):
//
// Copyright Information:
//	Copyright (C) 1996	Lawrence Berkeley Laboratory
//	Copyright (C) 1996	CEA - Centre d'Etude de Saclay
//
//------------------------------------------------------------------------

#ifndef DETPLANE_HH
#define DETPLANE_HH
//
//  other includes
//
#include "DetectorModel/DetSurface.hh"

class DetPlane : public DetSurface {

public:
//
//  Transform constructor;
//
  DetPlane(const HepTransformation&);
//
// Copy Constructor
//
  DetPlane( const DetPlane& );
//
// Destructor
//
  virtual ~DetPlane( );
//
//  copyOf function
//
  DetSurface* copyOf(double fparam=0.0) const;
//  
// methods which define the geometry of this surface
//
//  Closest distance to the DetPlane from the given point
//  Also output is the normal to the DetPlane  and the
//  2D-vector in 'local coordinates' at the point
//  on the surface closest to the given point
  virtual double normTo( const HepPoint&, // input point
			 Hep3Vector&) const;     // ouput normal to the surface
//
// first derivative of the distance to the surface from the given point
// as a function of position along the given direction
  virtual double firstDeriv( const HepPoint&,           // input point
			     const Hep3Vector& ) const; // input direction
//
// distance to the plane from a givenm point along a given direction
  virtual
  DetSurface::intertype distTo(const HepPoint&,   // input point
			       const Hep3Vector&, // input direction
			       double&,
			       intermode mode=closest)  const;
//
  bool wrappedCoordinate(int icoord) const { return false; }
//  min/max distance to a line segment defined by 2 points
  void segmentMinMax(const HepPoint&,const HepPoint&,double&,double&) const;
//
// Access functions
//
  HepPoint spacePoint( const SurfacePoint& ) const ;
  Hep3Vector normal( const SurfacePoint& ) const;
  Hep3Vector surfaceDirection( const SurfacePoint&, int ) const;
  int surfacePoint( const HepPoint&, SurfacePoint&, double tol=0.01 ) const;
  double curvature(const SurfacePoint& uv) const {return 0;}
// printout
  void print(std::ostream& os) const;
  void printAll(std::ostream& os) const;
};
#endif


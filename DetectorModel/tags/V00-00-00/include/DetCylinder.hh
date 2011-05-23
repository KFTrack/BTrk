//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DetCylinder.hh,v 1.11 2004/08/06 05:58:29 bartoldu Exp $
//
// Description:
//	Description of a DetCylinder for BaBar reconstruction
//      A DetCylinder is a Surface defined by equation x^2+y^2=R in its 
//      local coordinate system.
//      The 2D system of local coordinates is u=phi and v=z. 
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      David N Brown                - Lawrence Berkeley Lab
//	Gautier Hamel de Monchenault - CEN Saclay & Lawrence Berkeley Lab
//
// History (add to end):
//      Gautier Sep 26, 1996 - creation
//      BobJ    Feb 4, 1997 - change from Cylinder to DetCylinder
//
// Copyright Information:
//	Copyright (C) 1996	Lawrence Berkeley Laboratory
//	Copyright (C) 1996	CEA - Centre d'Etude de Saclay
//
//------------------------------------------------------------------------

#ifndef DETCYLINDER_HH
#define DETCYLINDER_HH
#include "DetectorModel/DetSurface.hh"

class DetCylinder : public DetSurface {

public:
//
//  Construct from a transform;
//
  DetCylinder(const HepTransformation&,double radius);
//
// Copy Constructor
//
  DetCylinder( const DetCylinder& );
//
// Destructor
//
  virtual ~DetCylinder( );
//
//  copyOf
//
  DetSurface* copyOf(double fparam = 0.0) const;
//  
// methods which define the geometry of this surface
//
//  Closest distance to the DetCylinder from the given point
//  Also output is the normal to the DetCylinder  and the
//  2D-vector in 'local coordinates' at the point
//  on the surface closest to the given point
  virtual double normTo( const HepPoint&, // input point
			 Hep3Vector&) const;
//
// first derivative of the distance to the surface from the given point
// as a function of position along the given direction
  virtual double firstDeriv( const HepPoint&,           // input point
			      const Hep3Vector& ) const; // input direction
//
// distance to the DetCylinder from a givenm point along a given direction
  virtual
  DetSurface::intertype distTo(const HepPoint&,   // input point
			       const Hep3Vector&, // input direction
			       double&,
			       intermode mode=closest)  const;
//  Curvature at a surface point
//
  double curvature(const SurfacePoint&) const ;
//  Coordinate wrap: phi (0th coordinate) wraps
  bool wrappedCoordinate(int) const;
//  min/max distance to a line segment defined by 2 points
  void segmentMinMax(const HepPoint&,const HepPoint&,double&,double&) const;
//
// Access functions
//
  HepPoint spacePoint( const SurfacePoint& ) const ;
  Hep3Vector normal( const SurfacePoint& ) const;
  Hep3Vector surfaceDirection( const SurfacePoint&, int ) const;
  int surfacePoint( const HepPoint&, SurfacePoint&, double tol=0.01 ) const;
  double radius() const { return _r; } 
// printout
  void print(std::ostream& os) const;
  void printAll(std::ostream& os ) const;
private:
  double _r;
};
#endif


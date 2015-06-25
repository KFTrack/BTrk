//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DetCone.hh,v 1.7 2004/08/06 05:58:29 bartoldu Exp $
//
// Description:
//	Description of a DetCone for BaBar reconstruction
//      A DetCone is an infinite Surface defined by equation
//              R = |z-z0| * tan(theta) 
//      in its local coordinate system where theta is its half angle
//      and z is its vertex. The 2D system of local coordinates is u=R*phi 
//      and v=z. 
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      David N Brown                - Lawrence Berkeley Lab
//	Gautier Hamel de Monchenault - CEN Saclay & Lawrence Berkeley Lab
//      Phil Strother                - Imperial College 
//
// History (add to end):
//      Gautier Sep 26, 1996 - creation
//      Phil    March 1997   - adapt for cone
//
// Copyright Information:
//	Copyright (C) 1996	Lawrence Berkeley Laboratory
//	Copyright (C) 1996	CEA - Centre d'Etude de Saclay
//
//------------------------------------------------------------------------

#ifndef DETCONE_HH
#define DETCONE_HH
#include "BTrk/DetectorModel/DetSurface.hh"

class DetCone : public DetSurface {

public:
//
//  Construct from a transform;
//
  DetCone(const HepTransformation&,double tanTheta, double vertex);
  // tanTheta is the slope of the cone w.r.t the positive z axis and is constrained to 
  // lie in the range 0 < theta < pi/2 

  DetCone(const HepTransformation&,double z1, double z2, double r1, double r2=0.0);
  // The absolute values of r1 and r2 are taken (i.e. the z axis, in the local
  // coordinate system,  is assumed to be the axis of symmetry). 


  //
  // Copy Constructor
  //
  DetCone( const DetCone& );

  //
  // Destructor
  //
  virtual ~DetCone( );

  //
  //  copyOf
  //
  virtual DetSurface* copyOf(double fparam=0.0) const;

  //  
  // methods which define the geometry of this surface
  //
  //  Closest distance to the DetCone from the given point
  //  Also output is the normal to the DetCone  and the
  //  2D-vector in 'local coordinates' at the point
  //  on the surface closest to the given point

  virtual double normTo( const HepPoint&, // input point
			 CLHEP::Hep3Vector&) const;     // ouput normal to the surface

  //
  // first derivative of the distance to the surface from the given point
  // as a function of position along the given direction
  virtual double firstDeriv( const HepPoint&,           // input point
			     const CLHEP::Hep3Vector& ) const; // input direction
  
  //
  // distance to the DetCone from a givenm point along a given direction
  virtual
  DetSurface::intertype distTo(const HepPoint&,   // input point
			       const CLHEP::Hep3Vector&, // input direction
			       double&,
			       intermode=closest) const;           // output distance

  //
  //  Curvature at a surface point
  //
  virtual double curvature(const SurfacePoint&) const ;

  //
  //  Special function to find the min/max perpendicular distance from the surface
  //  to a line segment defined by 2 points.  Used in DetSurfaceSet. The returned 
  //  arguments are min, max in that order.
  //
  
  virtual void segmentMinMax(const HepPoint &, const HepPoint &, double &, double&) const;
  
  // Z shift
  void applyZShift(double zshift) { _z0 += zshift; }

  //
  // Access functions
  //
  virtual HepPoint spacePoint( const SurfacePoint& ) const ;
  virtual CLHEP::Hep3Vector normal( const SurfacePoint& ) const;
  virtual CLHEP::Hep3Vector surfaceDirection( const SurfacePoint&, int ) const;
  virtual int surfacePoint( const HepPoint&, SurfacePoint&, double tol=0.01 ) const;
  virtual bool wrappedCoordinate(int) const;

  //
  // data members
  //

  inline double tanTheta() const { return _tanTheta; }
  inline double cosTheta() const { return _cosTheta;}
  inline double sinTheta() const { return _sinTheta;}
  inline double vertex() const { return _z0;}
// printout
  void print(std::ostream& os) const;
  void printAll(std::ostream& os) const;
private:
  double _tanTheta;
  double _z0;
  
  // Useful qauntities not necessary for the representation
  double _cosTheta;
  double _sinTheta;
};
#endif


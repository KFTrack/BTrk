//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DetSCone.hh,v 1.4 2004/08/06 05:58:31 bartoldu Exp $
//
// Description:
//	Description of a DetSCone for BaBar reconstruction
//
//      A DetSCone is an semi-infinite DetSurface defined by equation
//
//        (z-z0)^2-tan(alpha)^2*rho^2=0 ; z0!=0, tan(alpha)!=0 
//
//      in its local coordinate system where alpha is its half angle
//      and z0 is its vertex.
//
//      The 2D system of local coordinates are u=phi and v=|r-vertex|/scale.
//      For "scones" created with copyOf, scale proportional to |z0| so 
//      coordinates of the central projection are the same for all clones.
//
//      Cordinate v positive. The surface divide space on inside and
//      outside parts. Inside part is sharp angle (convex) one. Normal vector
//      is directed from inside to outside.
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      David N Brown                - Lawrence Berkeley Lab
//	Gautier Hamel de Monchenault - CEN Saclay & Lawrence Berkeley Lab
//      Phil Strother                - Imperial College 
//      Alexander Korol              - BINP, Novosibirsk
//
// History (add to end):
//      Gautier Sep 26, 1996 - creation
//      Phil    March 1997   - adapt for cone
//      Alex    May 1999     - make "scone" surface
//
// Copyright Information:
//	Copyright (C) 1996	Lawrence Berkeley Laboratory
//	Copyright (C) 1996	CEA - Centre d'Etude de Saclay
//      Copyright (C) 1999      BINP, Novosibirsk
//
//------------------------------------------------------------------------

#ifndef DETSCONE_HH
#define DETSCONE_HH

#include "BTrk/DetectorModel/DetSurface.hh"

class DetSCone : public DetSurface {

  friend class DetSConeEquation; // helper class described in implementation

public:
  // returns the region of cone and closest (signed) distance to it
  enum coneRegion { nearVertex, onVertex, onAxis, insideCone, outsideCone };

//
//  Construct from a transform;
//
  DetSCone(const HepTransformation&, // (DetSurface) transformation
	   double,                   // tan(theta), signed to positive z
	   double,                   // z0, coordinate of vertex
	   double=1.0);              // scale, used in "copyOf" method

  DetSCone(const HepTransformation&, // (DetSurface) transformation
	   double,                   // z (local) of the first point
	   double,                   // z (local) of the second point
	   double,                   // distance to axis for the first point
	   double,                   // distance to axis for the second point
	   double=1.0);              // scale, used in "copyOf" method

  //
  // Copy Constructor
  //
  DetSCone( const DetSCone& );

  //
  // Destructor
  //
  virtual ~DetSCone( );

  //
  //  copyOf clones cone so that distance between new and old surface (signed)
  //  is its argument. Must be consistent with normTo and segmentMaxMin.
  //
  virtual DetSurface* copyOf(double=0.0) const;

  //
  // methods which define the geometry of this surface
  //

  //  Returns closest distance (signed!) to the DetSCone from the given point
  //  Also output is the normal to the DetSCone. Must be consistent
  //  with segmentMinMax and with copyOf (reverse sign).
  virtual double normTo( const HepPoint&,        // input point
			 CLHEP::Hep3Vector&) const;     // ouput normal to the surface

  // first derivative of the distance to the surface from the given point
  // as a function of position along the given direction
  virtual double firstDeriv( const HepPoint&,           // input point
			     const CLHEP::Hep3Vector&) const;  // input direction

  // distance to the DetSCone from a givenm point along a given direction
  virtual
  DetSurface::intertype distTo(const HepPoint&,   // input point
			       const CLHEP::Hep3Vector&, // input direction
			       double&,           // output distance (signed)
			       intermode=closest) const;
  //
  //  Special function to find the min/max perpendicular distance from the 
  //  surface to a line segment defined by 2 points.  Used in DetSurfaceSet.
  //  The returned arguments are min, max in that order. Must be consistent
  //  with copyOf and normTo.
  virtual void
  segmentMinMax(const HepPoint&,  // input, first point
		const HepPoint&,  // input, second point
		double&,          // output, minimal (signed) distance
		double&) const;   // output, maximal (signed) distance

  // returns true if the point is close to surface and its "surface" coords.
  virtual int surfacePoint( const HepPoint&,     // input
			    SurfacePoint&,
			    double=0.01 ) const;

  //
  // geometry available from the surface point view
  //
  virtual double curvature(const SurfacePoint&) const ;
  virtual HepPoint spacePoint( const SurfacePoint&) const ;
  virtual CLHEP::Hep3Vector normal( const SurfacePoint&) const;
  virtual CLHEP::Hep3Vector surfaceDirection( const SurfacePoint& sp, int indx) const;

  // is the specified coordinate cyclic ?
  virtual bool wrappedCoordinate(int) const;

  // printout
  virtual void print(std::ostream& os) const;
  virtual void printAll(std::ostream& os) const;

  //
  // data accessors
  //
  inline double tanAlpha() const { return _tanAlpha; }
  inline double z0() const { return _z0; }
  inline double scale() const { return _scale; }

private:
  double _tanAlpha;
  double _z0;
  double _scale; // to scale r-coordinate for "cloned" cones

  // not definitive but frequently used values cached
  void updateCache(); // use only in constructors
  double _cosAlpha;
  double _sinAlpha;

  // helper functions

  // some alternative angle accessors (partially cached)
  // make them private to have  "inlined" in implementation
  double cosA() const;  // positive or negative (angle in [0,pi])
  double sinA() const;  // positive             (angle in [0,pi])
  double sSinA() const; // positive or negative (angle in [-pi/2,pi/2])
  double uCosA() const; // positive             (angle in [-pi/2,pi/2])

  inline coneRegion theRegion(const CLHEP::Hep3Vector&, // input point
			  double&) const;    // output, distance to it

};
#endif

//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DetSurface.hh,v 1.15 2007/09/20 23:01:42 gapon Exp $
//
// Description:
//	DetSurface virtual base Class 
//      Base class for describing a surface in 3-space.  A local
//      2-d coordinate system is assumed for the subclasses, which
//      can be of any type (cartesian, polar,  etc.).
//      the so-called 'surface coordinates'. 
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
//  Dave Brown, 9/19/96
//  modified by Gautier, 9/23/96
// ------------------------------------------------------------------------------
#ifndef DETSURFACE_HH
#define DETSURFACE_HH

//----------------
// BaBar header --
//----------------
#include <iostream>
class HepPoint;
#include "CLHEP/Vector/ThreeVector.h"
class HepTransformation;
enum Axis  { I_x, I_y, I_z };
//
//  Structure to define a point on the surface (IE surface coordinates).  This
//  is interpreted locally according to the surface geometry; planes use
//  cartesian coordinates, cylinders use cylindrical.
//
struct SurfacePoint{
  double _surfpoint[2];
  SurfacePoint(){;}
  SurfacePoint(const double& a,const double& b){
    _surfpoint[0] = a;
    _surfpoint[1] = b; }
  SurfacePoint(const double* ab){
    _surfpoint[0] = ab[0];
    _surfpoint[1] = ab[1]; }
  double& operator [](int icoord) { return _surfpoint[icoord]; }
  double operator [](int icoord) const { return _surfpoint[icoord]; }
  SurfacePoint& operator = (const SurfacePoint& other) {
    _surfpoint[0] = other._surfpoint[0];
    _surfpoint[1] = other._surfpoint[1];
    return *this;
  }
  double* array(){ return _surfpoint;}
//
//  Functions for specific interpretations
//
  double& phi() { return _surfpoint[0]; }
  double& zp() {return _surfpoint[1]; }
  double& xp() {return _surfpoint[0]; }
  double& yp() {return _surfpoint[1]; }
};
//
//  Define the main class
//
class DetSurface{
public:
// distanceTo return values
  enum intertype {intersect=0,outofrange,localmin,nointersect };
// distT0 modes
  enum intermode {forward=0,backward=1,closest=2 };
//
//  Construct from a transform.
//
  DetSurface(const HepTransformation&);
//  copy constructor
  DetSurface(const DetSurface& other);
//
//  Destructor
//
  virtual ~DetSurface();
//
//  Copy function.  This includes the possibilty of generating
//  a family of surfaces from a given one, according to a parameter.
//  The parameter represents the perpendicular 'size' of the new
//  surface relative to the old.
//
  virtual DetSurface* copyOf(double fparam=0.0) const = 0;
//
//  Return the direction normal to the surface closest to the
//  input point (returned value is the closest distance to the
//  surface).  The sign convention on the distance is that the HepPoint on the
//  surface can be calculated as surfpoint = testpoint + dist*norm.
  virtual 
  double normTo( const HepPoint& testpoint ,CLHEP::Hep3Vector& norm) const = 0;
//
// same, including computation of the associated surface point
  double normalTo( const HepPoint&,CLHEP::Hep3Vector&,SurfacePoint&) const;
//
// first derivative of the distance to the surface from the given point
// as a function of position along the given direction
//
  virtual 
  double firstDeriv( const HepPoint&,
		     const CLHEP::Hep3Vector& ) const = 0; 
//
//  Return the distance to the surface from a given point along
//  the given direction.  The behavior is now explicit regarding whether
//  the intersection forward, backward, or closest to the surface is desired
//  (the old convention of forward is preserved as default).  The return
//  value is now an enum, with success the same value as before, but with
//  failure now specifying a variety of conditions.
//  The sign convention of the returned distance is such that the HepPoint
//  on the surface can be computed as
//  surfpoint = testpoint + distance*dir.norm()
//
  virtual 
  intertype distTo(const HepPoint& testpoint ,const CLHEP::Hep3Vector& dir,
		   double& dist,
		   intermode mode=closest) const = 0;
// same, returning surface point.  This onl
  intertype distanceTo( const HepPoint&,const CLHEP::Hep3Vector&, double&,
			SurfacePoint&,intermode mode=closest) const;
//
// scalar curvature at a surface point.  This is the maximum of the 2-d curvature
// for all points
//
  virtual double curvature( const SurfacePoint& ) const = 0;
//
//  Special function to find the min/max perpendicular distance from the surface
//  to a line segment defined by 2 points.  Used in DetSurfaceSet
//
  virtual void segmentMinMax(const HepPoint&,const HepPoint&,double&,double&) const = 0;
//
//  Access functions
//
  HepTransformation* transform(){ return _transf;}
  const HepTransformation* transform() const { return _transf;}
  HepPoint   centerPoint() const;
  CLHEP::Hep3Vector basis(int ax=I_z) const;
  CLHEP::Hep3Vector axis() const;
  HepPoint   gotoLocal ( const HepPoint& )   const;
  HepPoint   gotoGlobal( const HepPoint& )   const;
  CLHEP::Hep3Vector gotoLocal ( const CLHEP::Hep3Vector& ) const;
  CLHEP::Hep3Vector gotoGlobal( const CLHEP::Hep3Vector& ) const;
//
// Functions that define the 'surface system of coordinates'
//  space point for a point on the surface for the given surface coordinates 
  virtual HepPoint   spacePoint(const SurfacePoint&) const = 0;
//  normal to the surface at a surface point for the given surface coordinates
  virtual CLHEP::Hep3Vector normal(const SurfacePoint&) const = 0;
//  Surface coordinate directions at a surface point
  virtual CLHEP::Hep3Vector surfaceDirection(const SurfacePoint&,int) const = 0;
//  returns the surface coordinates for a given space point, with a condition
//  flag if the point is not on the surface within tolerances.
  virtual int surfacePoint(const HepPoint&, SurfacePoint&,double tol=0.01 ) const = 0;
//  Defines the coordinates as infinite or wrapped.  Wrapped coordinates must
//  go from 0 to 2pi.
  virtual bool wrappedCoordinate(int) const = 0;
// superseed operator ==
  bool operator==( const DetSurface& otherSurface ) const { return false; }
// printout
  virtual void print(std::ostream& os) const;
  virtual void printAll(std::ostream& os) const;
private:
  HepTransformation* _transf;  // transformation describing the 
                              // local system of coordinates 
                              // inside the global system
// prohibit
  DetSurface& operator = (const DetSurface&);
};

#endif

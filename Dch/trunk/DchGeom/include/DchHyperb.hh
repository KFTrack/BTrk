//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DchHyperb.hh 91 2010-01-14 12:37:23Z stroili $
//
// Description:
//	Class DchHyperb: "simple" with circular cross section hyperboloid 
//      surface of equation
//
//          x^2   y^2   z^2
//          --- + --- - --- = 1
//          a^2   a^2   c^2
//
//          a^2*x^2 + a^2*y^2 - c^2*z^2 = 1
//
//      in the local reference frame
//      The surface coordinates are given by [phi,z]  
//          
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	R. Stroili		originator
//	
//
// Copyright Information:
//	Copyright (C) 1997	INFN-Pd
//
//------------------------------------------------------------------------

#ifndef DCHHYPERB_HH
#define DCHHYPERB_HH

//-------------
// C Headers --
//-------------
#include <math.h>

//---------------
// C++ Headers --
//---------------

//----------------------
// Base Class Headers --
//----------------------
#include "DetectorModel/DetSurface.hh"

//-------------------------------
// Collaborating Class Headers --
//-------------------------------

//------------------------------------
// Collaborating Class Declarations --
//------------------------------------

//		---------------------
// 		-- Class Interface --
//		---------------------

//DchHyperb<class X> 
class DchHyperb : public DetSurface {

//--------------------
// Declarations     --
//--------------------

  // Typedefs, consts, and enums

//--------------------
// Instance Members --
//--------------------

public:

  // Constructors
//
//  Construct from a transform;
//
  DchHyperb(const HepTransformation&,double a, double c, double z,
	    double twist);
  DchHyperb(const HepTransformation&,double rend, double zlength, 
	    double twist);
//
// Copy Constructor
//
//  DchHyperb( const DchHyperb& );
//
// Destructor
//
  virtual ~DchHyperb( );
//
//  Copy function.  This includes the possibilty of generating
//  a family of surfaces from a given one, according to a parameter.
//  The parameter represents the perpendicular 'size' of the new
//  surface relative to the old.
//
  DetSurface* copyOf(double fparam=0.0) const;
  
//  
// methods which define the geometry of this surface
//
//  Closest distance to the DchHyperboloid from the given point
//  Also output is the normal to the DchHyperb  and the
//  2D-vector in 'local coordinates' at the point
//  on the surface closest to the given point
  double normTo( const HepPoint&,         // input point
		 Hep3Vector& ) const;     // ouput normal to the surface

//   double normalTo( const HepPoint&, // input point
// 		   Hep3Vector&,     // ouput normal to the surface
// 		   SurfacePoint&) const;  // local coordinates at closest app
//
// first derivative of the distance to the surface from the given point
// as a function of position along the given direction
// Similar to Surface:dHowNear in gismo 2.
//
  double firstDeriv( const HepPoint&,          
		     const Hep3Vector& ) const; 
//
//  Return the distance to the surface from a given point along
//  the given direction (this distance might be negative if the point
//  is 'in front' of the surface).  A 2-vector in 'surface coordinates'
//  at the intersection point is also returned.
//  Similar to Surface::distance in gismo 2.
//
  DetSurface::intertype distTo( const HepPoint&, const Hep3Vector&,
				double&, intermode mode=forward ) const;
//   int distanceTo( const HepPoint&,
// 		  const Hep3Vector&, double&, 
// 		  SurfacePoint& ) const;
//
//
// scalar curvature at a surface point.  This is the maximum of the 2-d 
// curvature for all points
//
  double curvature( const SurfacePoint& uv ) const;
//
//  Special function to find the min/max perpendicular distance from the surface
//  to a line segment defined by 2 points.  Used in DetSurfaceSet
//
  void segmentMinMax( const HepPoint&, const HepPoint&, double&,
		      double&) const;
//
//  Accessor functions
// normal to surface at given point in space on the surface
  Hep3Vector normal( const SurfacePoint& ) const;
  int surfacePoint( const HepPoint&, SurfacePoint&, 
		    double tol=0.01 ) const;
  HepPoint spacePoint( const SurfacePoint& ) const;
  Hep3Vector surfaceDirection(const SurfacePoint&,int) const;
  bool wrappedCoordinate(int) const;

  // Operators
  //            --->   all these work in the local reference frame
  double radFromZ( const double& z ) const { 
    return sqrt( (1.+z*z*(_parC*_parC)) )/_parA; }
//   double zFromRad( const double& ) const;
  double radFromSurf( const SurfacePoint& uv ) const { 
    return radFromZ(uv[1]); }
  double rMax(void) const { return sqrt(1+_zmax*_zmax*(_parC*_parC))/_parA; }
  double rMin(void) const { return 1./_parA; }

  // it gets a point in the global reference point and returns an angle in
  // the global reference frame
  double projectPhiToRear( const HepPoint& point ) const;
  double twist( void ) const { return _twist; }

//    DchHyperb&       operator= ( const DchHyperb& );
//    virtual Boolean operator==( const DchHyperb& ) const;
//            Boolean operator!=( const DchHyperb& ) const;

  // Selectors (const)

private:

  // Data members
  double _parA;
  double _parC;
  double _zmax;  // max z
  double _twist;

};


#endif

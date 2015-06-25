//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DetSCone.cc,v 1.11 2004/08/06 05:58:31 bartoldu Exp $
//
// Description:
//	DetSCone Class.
//      Details in DetSCone.hh
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      David N Brown                - Lawrence Berkeley Lab
//	Gautier Hamel de Monchenault - CEN Saclay & Lawrence Berkeley Lab
//      Alexander Korol              - BINP, Novosibirsk
//
// Copyright Information:
//	Copyright (C) 1996	Lawrence Berkeley Laboratory
//	Copyright (C) 1996	CEA - Centre d'Etude de Saclay
//      Copyright (C) 1999      BINP, Novosibirsk
//
//------------------------------------------------------------------------

//----------------
// BaBar Header --
//----------------
#include "BTrk/BaBar/BaBar.hh"

//-----------------------
// This Class's Header --
//-----------------------
#include "BTrk/DetectorModel/DetSCone.hh"

//---------------
// C++ Headers --
//---------------
#include <assert.h>
#include <cfloat>
#include <math.h>
#include <iostream>

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "BTrk/BaBar/Constants.hh"
#include "BTrk/BbrGeom/HepPoint.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "BTrk/BbrGeom/Transformation.h"
using std::endl;
using std::ostream;
//-----------------------------------------------------------------------
// Local Macros, Typedefs, Structures, Unions and Forward Declarations --
//-----------------------------------------------------------------------

static inline double 
sconeMax(double x, double y) { return x>y?x:y; } 

static inline double 
sconeMin(double x, double y) { return x<y?x:y; } 

//
// Helper class which describes and solve square equation of the cone surface
//
class DetSConeEquation {

public:
  DetSConeEquation( const DetSCone* cone, const Hep3Vector& r, const Hep3Vector& vu);
  int findRoots(double x[2]) const;
  int findExtrema(double x[2], double dist[2]) const;
  double a() const { return _a; }
  double b() const { return _b; }
  double c() const { return _c; }

private:
  const double _thelimit;
  const DetSCone *const _cone;
  const Hep3Vector& _r;
  const Hep3Vector& _vu;
  double _r_vu_prod; // scalar production r*vu
  double _a, _b, _c; // square root coeffs

};

//
// build equation coefficients and store cone/vector information
//
inline
DetSConeEquation::DetSConeEquation(const DetSCone* cone, const Hep3Vector& r, const Hep3Vector& vu)
  : _thelimit(1e7),_cone(cone), _r(r), _vu(vu)
{
  double vuZ = _vu.z()/_cone->cosA();
  double rZ = _r.z()/_cone->cosA();
  
  _r_vu_prod = _r*_vu;
  _a = 1.0-vuZ*vuZ;       // 
  _b = _r_vu_prod-vuZ*rZ; // 
  _c = _r.mag2()-rZ*rZ;   // ==0 means we are on the surface, >0 outside, <0 inside
}


//
// find roots of the square equation
//
inline int
DetSConeEquation::findRoots(double x[2]) const
{
  // check is it really square?
  if( fabs(_a)<=fabs(_b)/_thelimit ) {
    // no, it is (almost) linear
    if( 2.0*fabs(_b)<=fabs(_c)/_thelimit ) {
      return 0; // no solutions in reasonable region
    } else {
      // linear and solvable
      x[0] = -0.5*_c/_b;
      return 1;
    }
  } else {
    // classic square equation
    double D = _b*_b-_a*_c;
    if( D<0.0 ) {
      return 0; // no solutions
    } else if( D==0.0 ) {
      // one solution
      x[0] = -_b/_a;
      return 1;
    } else {
      // two solutions, order NOT guaranteed
      x[0] = (-_b-sqrt(D))/_a;
      x[1] = (-_b+sqrt(D))/_a;
      return 2;
    }
  }
}

//
// find points with extremal distance to the surface and vertex
//
inline int
DetSConeEquation::findExtrema( double x[2], double dist[2]) const
{
  int roots = 0;

  // find the local minimum for distance to vertex
  x[0] = -_r_vu_prod;
  const Hep3Vector rvertex = _r+x[0]*_vu;
  if( _cone->theRegion( rvertex, dist[0] )==DetSCone::nearVertex ) {
    // and in the "vertex" region
    roots++;
  }

  // find the local extremum for distance to cone surface
  if( fabs(_a)>fabs(_b)/_thelimit ) {
    // it can exist inside the reasonable range
    x[roots] = -_b/_a;

    const Hep3Vector rcone = _r+x[roots]*_vu;
    if( _cone->theRegion(rcone, dist[roots])!=DetSCone::nearVertex ) {
      // and "cone surface" region
      roots++;
    }
  }

  return roots;
}

//----------------
// Constructors --
//----------------
//

void
DetSCone::updateCache()
{
  // alpha in (0,pi)
//  _cosAlpha = copysign(1/sqrt(1+_tanAlpha*_tanAlpha),_tanAlpha); // + or -
  _cosAlpha = (_tanAlpha < 0.0) ? -1/sqrt(1+_tanAlpha*_tanAlpha) : 1/sqrt(1+_tanAlpha*_tanAlpha); 
  _sinAlpha = _tanAlpha*_cosAlpha;                               // positive
}

DetSCone::DetSCone(const HepTransformation& t, 
		   double tanAlpha, double vertex, double scale )
  : DetSurface(t), _tanAlpha(tanAlpha), _z0(vertex), _scale(scale)
{
  // _tanAlpha is the slope of the cone w.r.t the positive z axis 
  // _tanAlpha can positive or negative.
  assert(_tanAlpha!=0.0);
  assert(_z0!=0);
  assert(_scale>0);

  updateCache();
}

DetSCone::DetSCone(const HepTransformation& t,
		   double z1, double z2, double r1, double r2, double scale)
  : DetSurface(t), _scale(scale)
{
  // The absolute values of r1 and r2 are taken (i.e. the z axis, in the local
  // coordinate system,  is assumed to be the axis of symmetry). 
  assert(r1!=r2);
  assert(z1!=z2);
  assert(scale>0);

  _tanAlpha = (r2-r1)/(z2-z1);
  double r, z;
  if( r2>r1 ) {
    z = z2;
    r = r2;
  } else {
    z = z1;
    r = r1;
  }
  _z0 = z-r/_tanAlpha;
  assert( _z0!=0.0 );
  updateCache();
}

//
// Copy Constructor
//

DetSCone::DetSCone( const DetSCone& c )
  : DetSurface(c), _tanAlpha(c._tanAlpha), _z0(c._z0), _scale(c._scale)
{
  updateCache();
}

//--------------
// Destructor --
//--------------
DetSCone::~DetSCone( )
{
}

inline double
DetSCone::cosA() const
{
  return _cosAlpha;
}

inline double
DetSCone::sinA() const
{
  return _sinAlpha;
}

inline double 
DetSCone::sSinA() const
{
//  return copysign(_sinAlpha,_tanAlpha);
  return (_tanAlpha < 0.0) ? -_sinAlpha : _sinAlpha;
}

inline double 
DetSCone::uCosA() const
{
  return fabs(_cosAlpha);
}

DetSCone::coneRegion
DetSCone::theRegion(const Hep3Vector& point, double& dist) const
{
  const double z   = point.z();
  const double rho = point.perp();
  assert(!(rho<0));
  if( rho>0 ) {
    const double along = z*cosA()+rho*sinA();
    if( along<0 ) {
      dist = sqrt(rho*rho+z*z);
      return nearVertex;
    } else {
      dist = -(z*sSinA())+rho*uCosA();
      return (dist>0)?outsideCone:insideCone;
    }
  } else { // rho==0, axis
    if( z==0 ) { // vertex
      dist = 0;
      return onVertex;
    } else if( (z*sSinA())>0 ) {
      dist = -(z*sSinA()); // inside cone, so dist(to surface)<0
      return onAxis;
    } else {
      dist = fabs(z);     // outside cone, dist(to vertex)>0
      return nearVertex;
    }
  }
}

//
//  copyOf function
//
DetSurface*
DetSCone::copyOf(double fparam) const
{
  if(fparam == 0.0) {
    return (DetSurface*) new DetSCone(*this);
  } else {
    //
    //  Moving the vertex for fabs(fparam) to be distance
    //  between surfaces and out of sharp angle direction
    //
    double newz = _z0-fparam/sSinA();
    double newscale = _scale*newz/_z0;
    assert(newscale>0); // make sure that we have not crossed the centre
    return (DetSurface*) new DetSCone(*transform(), _tanAlpha, newz, newscale );
  }
}

// 
// Member functions
//
HepPoint
DetSCone::spacePoint(const SurfacePoint& uv) const
{
  assert(uv[1]>=0); // r>=0

  if( uv[1]>0 ) {
    double rho = _scale*uv[1]*sinA();
    double x   = rho*cos(uv[0]);
    double y   = rho*sin(uv[0]);
    double z   = _z0+_scale*uv[1]*cosA();
    return gotoGlobal( HepPoint( x, y, z ) );
  } else { // vertex
    return gotoGlobal( HepPoint( 0.0, 0.0, _z0 ) );
  }
}

Hep3Vector
DetSCone::normal(const SurfacePoint& uv) const
{
  assert(uv[1]>0); // r>0, singularity when r==0

  double x =  uCosA()*cos( uv[0] );
  double y =  uCosA()*sin( uv[0] );
  double z = -sSinA();
  return gotoGlobal( Hep3Vector( x, y, z ) );
}

Hep3Vector
DetSCone::surfaceDirection(const SurfacePoint& uv,int idir) const
{
  assert(uv[1]>0); // r>0, singularity when r==0

  double x, y, z;
  switch(idir){
  case 1:
    x = sinA()*cos( uv[0] );
    y = sinA()*sin( uv[0] );
    z = cosA();
    break;
  case 0: default:
    x = -sin( uv[0] );
    y =  cos( uv[0] );
    z =  0.0;
    break;
  }    
  return gotoGlobal( Hep3Vector( x, y, z) );
}

int 
DetSCone::surfacePoint( const HepPoint& x, SurfacePoint& uv, double tol ) const
{
  assert(_scale>0.0);

  Hep3Vector r( gotoLocal(x)-HepPoint(0.0,0.0,_z0) );

  uv[0] = r.phi();
  if( uv[0]<0 ) uv[0] += Constants::twoPi;

  uv[1] = r.mag()/_scale; // 1 coord is scaled distance from vertex along the surface

  double delta = (r.z()-_scale*uv[1]*cosA())/sinA(); // delta=>distance when distance=>0
  if( fabs( delta ) > tol ) return 1;

  return 0;
}

// have to check on math correctness though
double
DetSCone::curvature(const SurfacePoint& uv) const
{
  assert(uv[1]>0); // infinity when r==0, undefined when r<0
  return 1/(_scale*uv[1]*uCosA());
}

// 0: phi, 1: scaled r
bool
DetSCone::wrappedCoordinate(int icoord) const
{
  switch(icoord){
  case 0:
    return true;
  case 1: default:
    return false;
  }
}

//
// Here we use extended definition of normal as a direction to/from nearest point on the surface
// The direction is from inside to outside cone, or from vertex to point, what is closer.
// The sign is positive when inside cone and negative when outside (which opposite
// to segmentMinMax and copyOf usage)
//
double
DetSCone::normTo( const HepPoint& x, Hep3Vector& n ) const
{
  const Hep3Vector r = gotoLocal(x) - HepPoint(0.0,0.0,_z0);
  double dist;
  double phi;
  double temp;

  switch( theRegion(r, dist) ) {
  case insideCone:
  case outsideCone:
    // continuous surface place
    phi = r.phi();
    n = gotoGlobal(Hep3Vector(uCosA()*cos(phi), uCosA()*sin(phi), -sSinA()));
    break;
  case onAxis:
    // any phi, my choise =0.0
    n = gotoGlobal(Hep3Vector(uCosA(), 0.0, -sSinA()));
    break;
  case nearVertex:
    // normal is a direction from vertex to point
    n = gotoGlobal( r.unit() );
    break;
  case onVertex:
    // any possible outside-cone direction, why not along z
//    n = gotoGlobal(Hep3Vector(0.0, 0.0, -copysign(1.0,_tanAlpha)));
    temp = (_tanAlpha < 0.0) ? 1.0 : -1.0;
    n = gotoGlobal(Hep3Vector(0.0, 0.0, temp));
    break;
  default:
    ::abort(); // impossible value
  }
  return -dist; // negative sign, so (point+dist*n==surface_point)
}

double
DetSCone::firstDeriv( const HepPoint& x, const Hep3Vector& v ) const
{
  const Hep3Vector r = gotoLocal(x) - HepPoint(0.0,0.0,_z0);
  const Hep3Vector vu = gotoLocal(v).unit();
  double dist;

  switch( theRegion(r, dist) ) {
  case insideCone:
  case outsideCone:
    // continuous surface place: (-z*sin(const)+rho*cos(const))'=
    return -vu.z()*sSinA()+(vu.x()*r.x()+vu.y()*r.y())*uCosA()/r.perp();
  case onAxis:
    // on axis, ru.perp()==0, same, but rho'=vperp on the axis
    return -vu.z()*sSinA()+vu.perp()*uCosA();
  case onVertex:
    // on vertex, depends on where are we moving
    if( vu.x()*cosA()+vu.y()*sinA()<=0 ) {
      // move toward "near vertex" region, dr/dlen=1
      return 1.0;
    } else {
      // move toward "contigous" region, and still on the axis, like "onAxis"
      return -vu.z()*sSinA()+vu.perp()*uCosA();
    }
  case nearVertex:
    // derivative of distance to vertex
    return vu*r.unit();
  default:
    ::abort(); // impossible value
    return FLT_MAX; // for compilers
  }
}

// the same choose procedure for intersect and localmin
static inline DetSurface::intertype
chooseRoot(const int roots,
	   double x[2],
	   const DetSurface::intermode mode,
	   const DetSurface::intertype goodOne,
	   double& dist)
{
  switch( roots ) {
  case 2: // both solutions
    // ensure x[0]<x[1]
    if( x[0]>x[1] ) {
      double tmp = x[0];
      x[0] = x[1];
      x[1] = tmp;
    }
    switch( mode ) {
    case DetSurface::forward:
      if( x[0]>=0 ) {
	dist = x[0];
	return goodOne;
      } else {
	dist = x[1];
	return (dist>=0)? goodOne : DetSurface::outofrange;
      }
    case DetSurface::backward:
      if( x[1]<=0 ) {
	dist = x[1];
	return goodOne;
      } else {
	dist = x[0];
	return (dist<=0)? goodOne : DetSurface::outofrange;
      }
    case DetSurface::closest:
    default:
      dist = (fabs(x[0])<fabs(x[1]))?x[0]:x[1];
      return goodOne; 
    }
  case 1: // only one solution
      dist = x[0];
      switch( mode ) {
      case DetSurface::forward:
	return (dist>=0)? goodOne : DetSurface::outofrange;
      case DetSurface::backward:
	return (dist<=0)? goodOne : DetSurface::outofrange;
      case DetSurface::closest:
      default:
	return goodOne;    
      }
  case 0:
    return DetSurface::nointersect;
  default:
    ::abort(); // impossible
    return DetSurface::nointersect; // for compiler
  }
}

//
// well, it looks complicated, but every step is very simple actually
//
DetSurface::intertype
DetSCone::distTo( const HepPoint& p, const Hep3Vector& v, 
		  double& dist, intermode mode) const
{
  const Hep3Vector r = gotoLocal(p) - HepPoint(0.0,0.0,_z0);
  const Hep3Vector vu = gotoLocal(v).unit();

  DetSConeEquation eq(this, r, vu );

  // check for "any"/"no" solutions
  if( eq.b()==0.0 && eq.a()==0.0 ) {
    dist = 0.0;
    if( eq.c()==0.0 ) { // any solution
      // here (and everywhere) intersection
      return intersect;
    } else { // no solutions
      return nointersect;
    }
  }

  // solve the square equation
  double x[2];
  int roots = eq.findRoots(x);

  // check for fictive "cone surface".
  switch( roots ) {
  case 2:
    if( (r.z()+vu.z()*x[1])*_tanAlpha < 0 ) {
      roots--; // second solution is in fictive region
    }
    // no break !!
  case 1:
    if( (r.z()+vu.z()*x[0])*_tanAlpha < 0 ) {
      if( roots>1 ) x[0] = x[1]; // shift
      roots--;                   // first solution is in fictive region
    }
  case 0:
    break;
  default:
    ::abort(); // impossible value
  }

  // choose root and test with direction for intersections
  DetSurface::intertype inter = chooseRoot( roots, x, mode, intersect, dist );
  if( inter!=nointersect ) return inter;

  // no intersections, try to find local minima
  double fromsurf[2];
  roots = eq.findExtrema( x, fromsurf );
  switch( roots ) {
  case  2: // both "vertex" and "surface" extrema present (vertex first)
    if( fromsurf[1]<0 ) {
      // inside, can be maximum only
      roots--;
    }
    break;
  case  1: // only "surface" or vertex extremum presents
    if( fromsurf[0]<0 ) {
      // inside, can be maximum only
      roots--;
    }
    break;
  case  0: // no extrema
    // roots = 0;
    return nointersect;
  default:
    ::abort(); // impossible value
  }

  // choose root and test with direction for local minima
  return chooseRoot( roots, x, mode, localmin, dist );
}

// well it can be really difficult
void
DetSCone::segmentMinMax(const HepPoint& p1,const HepPoint& p2,
			double& mindist,double& maxdist) const
{
  double tmax = (p1-p2).mag();
  if( tmax==0.0 ) {
    // the same point, have to check, because we use unit on p1-p2
    const Hep3Vector r = gotoLocal(p1) - HepPoint(0.0,0.0,_z0); 

    DetSCone::coneRegion where = theRegion( r, mindist );
    maxdist = mindist;
  } else {
    const Hep3Vector r1 = gotoLocal(p1) - HepPoint(0.0,0.0,_z0); 
    const Hep3Vector r2 = gotoLocal(p2) - HepPoint(0.0,0.0,_z0);
    const Hep3Vector vu = (r2-r1).unit();

    DetSCone::coneRegion where1 = theRegion( r1, mindist );
    maxdist = mindist;

    double thedist;
    DetSCone::coneRegion where2 = theRegion( r2, thedist );
    maxdist = sconeMax(maxdist, thedist);
    mindist = sconeMin(mindist, thedist);

    // find the local extremum for distances to vertex and surface 
    DetSConeEquation eq(this, r1, vu);
    double x[2];
    double exdist[2];
    switch( eq.findExtrema( x, exdist ) ) {
    case 2: // both vertex and surface
      if( 0.0<x[1] && x[1]<tmax ) {
	// inside the segment
	maxdist = sconeMax(maxdist, exdist[1]);
	mindist = sconeMin(mindist, exdist[1]);
      }
      // no break!!
    case 1: // either vertex or surface
      if( 0.0<x[0] && x[0]<tmax ) {
	// inside the segment
	maxdist = sconeMax(maxdist, exdist[0]);
	mindist = sconeMin(mindist, exdist[0]);
      }
    case 0: // no extrema found, do nothing
      return;
    default:
      ::abort(); // imposiible value
    }
  }
} // end segmentMinMax

void
DetSCone::print(ostream& os) const {
  os << " S(ingle)Cone Surface with" 
     << " tan="   << _tanAlpha
     << " z0="    << _z0
     << " scale=" << _scale
     << endl;
} // end print

void
DetSCone::printAll(ostream& os) const {
  HepPoint   origin = transform()->origin();
  Hep3Vector axis = transform()->unit();
  Hep3Vector xaxis = transform()->unit(0);

  print(os);
  os << "Origin = " 
     << origin.x() << " " << origin.y() << " " << origin.z()
     << endl
     << "Axis direction = " 
     << axis.x() << " " << axis.y() << " " << axis.z() 
     << endl
     << "Local X direction = "
     << xaxis.x() << " " << xaxis.y() << " " << xaxis.z() 
     << endl;
} // end printAll

//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DchCell.cc 91 2010-01-14 12:37:23Z stroili $
//
// Description:
//	Class DchCell
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
//	NOTE:  Compile with "-DDCHGEOM_DEBUG" to get extensive debugging
//	       messages.
//
// Author List:
//	R. Stroili		originator
//
// Copyright Information:
//	Copyright (C) 1997	INFN-Pd
//
// Revision History:
//	20020409  M. Kelsey -- Collapse duplicate helixPath() functions
//		  into one.  Replace HelixTraj& argument with HelixParams
//		  and HepPoint.
//	20020411  Add DchSWire data member.
//		  bug fix -- operator<< ought to take const DchCell&
//	20020419  BUG FIX -- Undo part of changes made 20020409 above:
//		  helixPath MUST take HelixTraj&, and not HelixParams.
//	20020425  Protect debugging messages and calculations with #ifdef
//	20020428  Modify helixPath() to compute arc lengths from hit position
//		  directly, rather than from helix base, protected by #ifdef.
//	20020506  Remove #ifdef protection -- use hit-based length always
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"

//-----------------------
// This Class's Header --
//-----------------------
#include "DchGeom/DchCell.hh"

//-------------
// C Headers --
//-------------

//---------------
// C++ Headers --
//---------------
#include <algorithm>
#include <functional>
#include <functional>
#include <iterator>
#include <iostream>

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "BaBar/Constants.hh"
#include "BbrGeom/BbrAngle.hh"
#include "BbrGeom/Trajectory.hh"
#include "BaBar/BbrCollectionUtils.hh"
#include "CLHEP/Geometry/HepPoint.h"
#include "CLHEP/Geometry/Transformation.h"
#include "CLHEP/config/CLHEP.h"	// *** To get ATAN2
#include "CLHEP/Vector/ThreeVector.h"
#include "DchGeom/DchCellPlaneType.hh"
#include "DchGeom/DchFWire.hh"
#include "DchGeom/DchSWire.hh"
#include "DetectorModel/DetIntersection.hh"
#include "DetectorModel/DetPlane.hh"
#include "DetectorModel/DetSurfaceElem.hh"
#include "ErrLogger/ErrLog.hh"
#include "TrajGeom/TrkLineTraj.hh"
#include "TrkBase/HelixTraj.hh"
#include "TrkBase/TrkDifTraj.hh"
#include "TrkBase/TrkErrCode.hh"
#include "TrkBase/HelixParams.hh"
#include "TrkBase/TrkFit.hh"
#include "TrkBase/TrkParams.hh"
#include "TrkBase/TrkPoca.hh"
#include "TrkBase/TrkSimpTraj.hh"
using std::endl;
using std::ostream;

static const double INTERS_TOL = 0.1;
static const double Z_INTERVAL = 10.;
static const double POCA_PREC = 1.e-2;


//----------------
// Constructors --
//----------------
DchCell::DchCell(std::vector<DchFWire*>& fieldWires)
  : _cellWires(fieldWires),_sense(0), _cellNum(0)
{
}

DchCell::DchCell( const DchCell& copy )
  :_sense(copy._sense), _cellNum(copy._cellNum)
{
  for (int index=0; index<copy._cellWires.size(); index++) {
    DchFWire* wire = copy._cellWires[index];
    HepPoint newR = *wire->getRearPoint();
    HepPoint newF = *wire->getForwPoint();
    DchFWire* fWire = new DchFWire(newR,newF);
    _cellWires.push_back(fWire);
  }
}

DchCell::DchCell(const DchCell& copy, double phi, int cellnum,
		 const DchSWire* sense)
  : _sense(sense), _cellNum(cellnum)
{
  for (int index=0; index<copy._cellWires.size(); index++) {
    DchFWire* wire = copy._cellWires[index];
    HepPoint newR = *wire->getRearPoint();
    newR.rotateZ(phi);
    HepPoint newF = *wire->getForwPoint();
    newF.rotateZ(phi);
    DchFWire* fWire = new DchFWire(newR,newF);
    _cellWires.push_back(fWire);
  }
}    

DchCell::~DchCell()
{
  std::for_each(_cellWires.begin(),_cellWires.end(),babar::Collection::DeleteObject());
  _cellWires.clear();
}

int
DchCell::operator == ( const DchCell& otherCell ) const
{
  return _cellNum == otherCell._cellNum;
}

float
DchCell::intersect( const Trajectory* traj, double range[2] ) const 
{
  static std::vector<HepPoint> vec;
  vec.clear();
  return intersect(traj, range, vec);
}

float
DchCell::intersect( const Trajectory* traj, double range[2],
                    std::vector<HepPoint>& intersections ) const
{
  intersections.clear();  // clear intersections vector
  HepPoint start = traj->position((range[0]+range[1])*0.5);
  double z0r = start.z();
  // loop over couples of wires
  int index2=0;
  //  double zInterv=100.;

  int intersectionCounter = 0;
  std::vector< double > intersFLT;
  std::vector< TwoDCoord > theCorners;
  double firstZ = z0r-Z_INTERVAL;
  double lastZ = z0r+Z_INTERVAL;
  DetIntersection inters;

  for (int index=0; index<_cellWires.size(); index++) {
    index2 = (index+1)%_cellWires.size();
    // build the transform necessary to build the plane
    DchFWire* wire1 = _cellWires[index];
    DchFWire* wire2 = _cellWires[index2];

    // space points on wires at trajectory z
    HepPoint p1 = wire1->getPoint(z0r);
    HepPoint p2 = wire2->getPoint(z0r);
    // vectors to define the plane
    Hep3Vector v1 = p2 - p1;
    Hep3Vector v2 = wire1->getTraj()->direction(z0r);
    Hep3Vector norm = v1.cross(wire1->getTraj()->direction(z0r));
    //    Hep3Vector norm = wire1->getTraj()->direction(z0r).cross(v1);

    Hep3Vector vec(p1.x(), p1.y(), p1.z());

    theCorners.clear();
    // build the plane
//     double firstZ = traj->position(range[0]).z();
//     double lastZ = traj->position(range[1]).z();
    // first point: first wire lower z
    HepPoint point = wire1->getPoint(firstZ);
    TwoDCoord thePoint(0.,point.z()-z0r);
    theCorners.push_back(thePoint);
    // second point: first wire higher z
    point = wire1->getPoint(lastZ);
    thePoint = TwoDCoord(0.,point.z()-z0r);
    theCorners.push_back(thePoint);
    // third point: second wire higher z
    point = wire2->getPoint(lastZ);
    thePoint = TwoDCoord(v1.mag(),point.z()-z0r);
    theCorners.push_back(thePoint);
    // fourth point: second wire lower z
    point = wire2->getPoint(firstZ);
    thePoint = TwoDCoord(v1.mag(),point.z()-z0r);
    theCorners.push_back(thePoint);

    DchCellPlaneType plType("DchCellPlane",100000+index,theCorners);

    HepTransformation trans(vec,norm,v1);
    DetPlane plane(trans);
    DetSurfaceElem planeEl((DetSurfaceType*)&plType,"Cell Plane",200000+index,
                           plane);

    // intersect it

    inters.pathrange[0] = range[0];
    inters.pathrange[1] = range[1];
    inters.pathlen = range[0];
    int good = planeEl.intersect(traj,inters);
    if ( good ) {
//        cout<<"\t---\t"<<index<<"\t-\t"<<index2<<endl;
      HepPoint pcross = traj->position(inters.pathrange[0]);
      intersections.push_back(pcross);
      intersectionCounter++;
      intersFLT.push_back(inters.pathrange[0]);
    }
  } //end loop

  double path = 0.;
  if ( intersFLT.size() == 2 ) {
    range[0] = intersFLT[0]+0.5*(intersFLT[1]-intersFLT[0]);
    range[1] = range[0];
  } else if ( intersFLT.size() > 2 ) {
    double dist = 100.;
    bool found = true;
    while ( intersFLT.size() > 2 && found ) {
      found = false;
      std::vector<double>::iterator couple0 = intersFLT.end();
      std::vector<double>::iterator couple1 = intersFLT.end();
      std::vector<double>::iterator j;
      for (j=intersFLT.begin();j!=intersFLT.end();++j) {
        std::vector<double>::iterator jj;
        for ( jj=j+1; jj!=intersFLT.end(); ++jj ) {
          dist = fabs(*jj-*j);
          if ( dist < INTERS_TOL ) {
            couple0 = j;
            couple1 = jj;
            found = true;
          }
        }
      }
      if ( found ) {
        std::vector<HepPoint>::iterator c1=intersections.begin();
        c1 += couple1-intersFLT.begin();
        std::vector<HepPoint>::iterator c0=intersections.begin();
        c0 += couple0-intersFLT.begin();
        HepPoint p2 = *c1 ; intersections.erase(c1);
        HepPoint p1 = *c0 ; intersections.erase(c0);
        double int1 = *couple1 ; intersFLT.erase(couple1);
        double int2 = *couple0 ; intersFLT.erase(couple0);
        HepPoint midpoint = traj->position( (int1+int2)*0.5 );
        intersections.push_back(midpoint);
        intersFLT.push_back( (int1+int2)*0.5 );
      }
    }
  } else if ( intersections.size() != 2 &&
              intersections.size() != 0 ) {
    if ( ErrLogging(debugging) ) {
      ErrMsg(debugging) << "DchCell::intersect: # of intersections: "
			<< intersections.size() << endmsg;
    }
  }
  if ( intersFLT.size() == 2 ) {
    path = fabs(intersFLT[1]-intersFLT[0]);
  }
  if ( path == 0. ) {
    if ( ErrLogging(debugging) ) {
      ErrMsg(debugging) << "DchCell::intersect: null path in cell #"
			<< cellNum() << endmsg;
    }
  }
  return path;
}

float
DchCell::path( const Trajectory* traj, double trkRange[2], 
               double wireFlt ) const 
{
  static std::vector<HepPoint> vec;
  vec.clear();
  return path(traj, trkRange, wireFlt, vec);
}

float
DchCell::path( const Trajectory* traj, double trkRange[2], 
	       double wireFlt, 
               std::vector<HepPoint>& intersections ) const 
{
  // NOTE: ::intersect() does exactly the same thing, without any POCAs
  return intersect(traj, trkRange, intersections);
}

double
DchCell::helixPath(const TrkFit* traj, double fltlen) const
{
  if ( traj == 0 ) {
    ErrMsg(error) << "a non valid trajectory was given" << endmsg;
    return 0.;
  }

  // Return path length from local helix referenced to origin
  return helixPath(traj->helix(fltlen), HepPoint(0.,0.,0.),
		   traj->position(fltlen), fltlen);
}

double
DchCell::helixPath(const HelixTraj& theHelix, double fltlen) const
{
  // NOTE:  Must convert TrkParams to HelixParams by hand
  HelixParams pars(theHelix.parameters()->parameter(),
		      theHelix.parameters()->covariance());

  // Return path length from helix at built-in reference point
  return helixPath(pars, theHelix.referencePoint(),
		   theHelix.position(fltlen), fltlen);
}


// calculate the intersection of a helix with the cell walls and
// return the flight length distance between the entry and exit point
// FFW 25-APR-99

// This version of helixPath is INTERNAL to DchCell, and is called from
// both of the public interfaces above.

double
DchCell::helixPath(const HelixParams& helix, const HepPoint& ref,
		   const HepPoint& hit, double fltlen) const
{
#ifdef DCHGEOM_DEBUG
  ErrMsg(trace) << "DchCell::helixPath(): " << endl
		<< " " << helix << " referenced to " << ref << endl
		<< " through " << hit << endmsg;
#endif

  // z coordinate at hit position
  double z = hit.z();
  double lowZ  = z - Z_INTERVAL;
  double highZ = z + Z_INTERVAL; 

#ifdef DCHGEOM_DEBUG
  ErrMsg(trace) << " Z range is " << lowZ << " to " << highZ << endmsg;
#endif

  // Reference point of helix, used for intersections
  double rad = 1./helix.omega();
  double r2 = rad*rad;
  double phi0 = helix.phi0();
  double x = ref.x() - (rad + helix.d0()) * sin(phi0);
  double y = ref.y() + (rad + helix.d0()) * cos(phi0);
  double raddip = rad*sqrt(1.+helix.tanDip()*helix.tanDip());

  // Hit direction is perpendicular to pos. vector, parallel to phi0
  double xh = hit.x()-x, yh = hit.y()-y;
  double phih = atan2(xh/rad,-yh/rad);

#ifdef DCHGEOM_DEBUG
  // Detailed computations to verify that hit point lies on helix
  ErrMsg(trace) << " center=(" << x << "," << y << ")"
		<< ", radius=" << rad << ", raddip=" << raddip << endmsg;

  double r2hit = xh*xh + yh*yh;
  
  ErrMsg(trace) << " helix r^2 = " << r2 << ", hit r^2 = " << r2hit
		<< ", diff = " << r2hit-r2 << endmsg;
  
  // NOTE:  Implicit assumption that BaBar tracks only loop once!
  double delphi = phih - phi0;
  if (delphi*rad<0.) delphi += Constants::twoPi*rad/fabs(rad);

  double alen = delphi*rad;		// Arc length around circle
  
  ErrMsg(trace) << " phi dir at hit = " << phih << ", delta = " << delphi 
		<< ", arc length = " << alen << endmsg;
  
  double zlen = alen * helix.tanDip();	// Flight distance in z
  double zpos = ref.z() + helix.z0() + zlen;
  
  ErrMsg(trace) << " helix z = " << zpos << ", hit z = " << z
		<< ", diff = " << z-zpos << endmsg;
#endif	// DCHGEOM_DEBUG

  int index2(0);
  double intersects[40];
  int ninters(0);

  static HepPoint p1;
  static Hep3Vector v12, v13, vplane;

  // loop over field wires of this cell -- relies on particular ordering
  unsigned nWires = _cellWires.size();
  for (int index=0; index<nWires; index++) {
    index2 = (index+1)%nWires;

    // wires
    DchFWire* wire1 = _cellWires[index];  assert(0 != wire1);
    DchFWire* wire2 = _cellWires[index2]; assert(0 != wire2);

#ifdef DCHGEOM_DEBUG
    ErrMsg(trace) << " Looking for intersection in plane of field wires "
		  << index << " and " << index2 << ":" << endl
		  << " " << wire1->getPoint(lowZ) << " to"
		  << " " << wire2->getPoint(lowZ) << " to"
		  << " " << wire2->getPoint(highZ) << endmsg;
#endif

    p1 = wire1->getPoint(lowZ);    	  // space points on wire
    Hep3Vector v1(p1.x(),p1.y(),p1.z());  // vector to first wire

    v12 = wire2->getPoint(lowZ) - p1;    // wire 1 to wire 2 (plane's height)
    v13 = wire2->getPoint(highZ) - p1;   // wire 1 to wire 2 (plane's length)
    vplane = v12.cross(v13);		 // vector normal to plane

#ifdef DCHGEOM_DEBUG
    ErrMsg(trace) << " v12 = " << v12 << endl << " v13 = " << v13 << endl
		  << " vplane = " << vplane << endmsg;
#endif

    double a = vplane.x();		// components of the normal
    double b = vplane.y();
    double c = vplane.z();

#ifdef DCHGEOM_DEBUG
    ErrMsg(trace) << " Cross check: a,b,c = " << a << "," << b << "," << c
		  << endmsg;
#endif

    // now calculate the discriminant part
    double bigC  = a*x + b*y + c*z  - vplane.dot(v1);
    double bigD2 = r2*(a*a + b*b) - bigC*bigC;

#ifdef DCHGEOM_DEBUG
    ErrMsg(trace) << " bigC = " << bigC << ", bigD2 = " << bigD2 << endmsg;
#endif

    // calculate other terms
    if (bigD2<0) continue;    		// real roots only
    double bigD  = sqrt(bigD2);
    double denom = bigC + b*rad;

    double bigZ1 =  (bigD - a*rad) / denom;
    double bigZ2 = -(bigD + a*rad) / denom;

#ifdef DCHGEOM_DEBUG
    ErrMsg(trace) << " bigD = " << bigD << ", denom = " << denom
		  << ", bigZ1 = " << bigZ1 << ", bigZ2 = " << bigZ2 << endmsg;
#endif

    // Phi directions at intersection points
    double dphiZ1 = 2.*atan(bigZ1) - phih;
    if (dphiZ1 < -Constants::pi) dphiZ1 += Constants::twoPi;
    if (dphiZ1 >  Constants::pi) dphiZ1 -= Constants::twoPi;

    double dphiZ2 = 2.*atan(bigZ2) - phih;
    if (dphiZ2 < -Constants::pi) dphiZ2 += Constants::twoPi;
    if (dphiZ2 >  Constants::pi) dphiZ2 -= Constants::twoPi;

#ifdef DCHGEOM_DEBUG
    ErrMsg(trace)
      << " phiZ1-phih = " << dphiZ1 << ", phiZ2-phih = " << dphiZ2 << endmsg;
#endif

    // Intersections computed relative to hit location
    intersects[ninters++] = dphiZ1 * raddip;
    intersects[ninters++] = dphiZ2 * raddip;

#ifdef DCHGEOM_DEBUG
    ErrMsg(trace)
      << " int[" << ninters-2 << "]: " << intersects[ninters-2]
      << " int[" << ninters-1 << "]: " << intersects[ninters-1]
      << endmsg;
#endif
  }  // for (int index=0

  // loop over found intersects and choose one with smallest path length
  // (closest intersections: one before input hit point, one after)

  double path = 0.;
  if (ninters > 1) {
    int mint = -1, pint = -1;
    double mflt = -50000.0, pflt = 50000.0;
    for (int i=0; i<ninters; i++) {
      if (intersects[i] > 0.) {
	if (intersects[i] < pflt) {
	  pflt = intersects[i];
	  pint = i;
	}
      } else {
	if (intersects[i] > mflt) {
	  mflt = intersects[i];
	  mint = i;
	}
      }
    } // for (int i...

#ifdef DCHGEOM_DEBUG
    ostream& debs = ErrMsg(trace);
    debs << "Closest intersections are " << mint;
    if (mint>=0) debs << " (" << intersects[mint] << ")";
    debs << " to " << pint;
    if (pint>=0) debs << " (" << intersects[pint] << ")";
    debs << endmsg;
#endif

    if (pint>=0 && mint>=0) path = intersects[pint] - intersects[mint];
  } // if (ninters > 1)

  return path;
}


// This computation walks counterclockwise along the edges of the cell, as
// defined by the lines connecting the field wires; It then checks that
// the point is on the left of each of the edges. As soon as it is on the
// right for one of the edges, it knows that the point is _not_ within
// the cell. 
//
// WARNING: This routine relies on the ordering of the field wires in
//          _cellWires (as it is this ordering that defines the edges
//          and their orientation!)

bool
DchCell::inCell(const HepPoint& p) const
{
  // ugly statics, but they _do_ speed up things...
  static HepPoint current,next;
  static Hep3Vector v1,v2;

  unsigned n=_cellWires.size();
  double z = p.z();
  current = _cellWires[0]->getPoint(z);
  for (unsigned i=0; i<n; ++i) {
    next = _cellWires[(i+1)%n]->getPoint(z);
    v1 = next - current;
    v2 = p - current;
    if ( v1.x()*v2.y() < v1.y()*v2.x() ) return false;
    current = next;
  }
  return true;
}

void
DchCell::print(ostream& o) const 
{
  o << " DchCell #: " << _cellNum << " layer # "<< layerNum() 
    << " cell # " << wireNum() << endl;
}

void
DchCell::printAll(ostream& o) const 
{
  o << " DchCell #: " << _cellNum << " layer # "<< layerNum() 
    << " cell # " << wireNum() << endl;
  for (int index=0; index<_cellWires.size(); index++) {
    o << index;
    _cellWires[index]->printAll(o);
  }
  if (0 != _sense) o << *_sense;
}

ostream&  
operator<<(ostream& o, const DchCell& cell) 
{
  cell.print(o);
  return o;
}

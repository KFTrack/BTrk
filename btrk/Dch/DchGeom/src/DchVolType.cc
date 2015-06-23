//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DchVolType.cc 123 2010-04-29 14:41:45Z stroili $
//
// Description:
//	Class DchVolType
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	R. Stroili		originator
//	
//
// Copyright Information:
//	Copyright (C) 1997	INFN - Pd
//
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"

//-----------------------
// This Class's Header --
//-----------------------
#include "DchGeom/DchVolType.hh"

//-------------
// C Headers --
//-------------

//---------------
// C++ Headers --
//---------------

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "ErrLogger/ErrLog.hh"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Geometry/Translation.h"
#include "CLHEP/Geometry/Transformation.h"
#include "DetectorModel/DetCylinder.hh"
#include "DetectorModel/DetPlane.hh"
#include "DetectorModel/DetSurface.hh"

//-----------------------------------------------------------------------
// Local Macros, Typedefs, Structures, Unions and Forward Declarations --
//-----------------------------------------------------------------------

static const char rscid[] =
    "$Id: DchVolType.cc 123 2010-04-29 14:41:45Z stroili $";
const double _tolerance = 1.0e-10;

//		----------------------------------------
// 		-- Public Function Member Definitions --
//		----------------------------------------

//----------------
// Constructors --
//----------------
DchVolType::DchVolType(const char* name, double rmin, double rmax, double zmin,
    double zmax, const DetMaterial* theMat, int idnum) :
  DetVolumeType(name, idnum), _rmin(rmin), _rmax(rmax), _zmin(zmin),
      _zmax(zmax), _theMaterial(theMat), _debug(false), _isPhiSegmented(false)
{
  HepTransformation nullTransf;
  DetCylinder* inCyl = 0;
  inCyl = new DetCylinder(nullTransf, _rmin);
  mySides()->push_back(inCyl);
  DetCylinder* outCyl = 0;
  outCyl = new DetCylinder(nullTransf, _rmax);
  mySides()->push_back(outCyl);
  HepTransformation tform(Hep3Vector(0., 0., _zmin), Hep3Vector(0, 0, 1.));
  DetPlane* rEP = 0;
  rEP = new DetPlane(tform);
  mySides()->push_back(rEP);
  tform = HepTransformation(Hep3Vector(0., 0., _zmax), Hep3Vector(0, 0, 1.));
  DetPlane* fEP = 0;
  fEP = new DetPlane(tform);
  mySides()->push_back(fEP);
}

//--------------
// Destructor --
//--------------
DchVolType::~DchVolType()
{
}

//-------------
// Selectors --
//-------------
bool
DchVolType::insideLimitsOf(int side, const SurfacePoint& thisPoint) const
{
  bool answer = true;
  size_t i = side;
  double z;

  if (side >= sides()->size()) {
    ErrMsg(error) << " DchVolType: wrong surface " << side << endmsg;
    return false;
  }
  HepPoint point = ((*sides())[i])->spacePoint(thisPoint);
  double radius = sqrt(point.x() * point.x() + point.y() * point.y());

  if (ErrLogging(debugging)) {
    ErrMsg(debugging) << "DchVolType: " << point << "\n" << " radius: "
        << radius << " r diff: " << radius - _rmin << " " << radius - _rmax
        << " z: " << point.z() << " z diff: " << point.z() - _zmin << " "
        << point.z() - _zmax << endmsg;
  }

  z = point.z();
  if (((fabs(z - _zmin) < _tolerance || fabs(z - _zmax) < _tolerance)
      && (radius > _rmin && radius < _rmax)) || ((fabs(radius - _rmin)
      < _tolerance || fabs(radius - _rmax) < _tolerance) && (z > _zmin && z
      < _zmax))) {
//if ( (z - _zmin) > -_tolerance && (z - _zmax) < _tolerance
//    && (radius - _rmin) > -_tolerance && (radius - _rmax) < _tolerance) {
    answer = true;
  } else {
    if (ErrLogging(debugging)) {
      ErrMsg(debugging) << "DchVolType: " << (fabs(z - _zmin) < _tolerance)
          << "\t" << (fabs(z - _zmax) < _tolerance) << "\t" << (radius > _rmin)
          << "\t" << (radius < _rmax) << "\t" << (fabs(radius - _rmin)
          < _tolerance) << "\t" << (fabs(radius - _rmax) < _tolerance) << "\t"
          << (z > _zmin) << "\t" << (z < _zmax) << endmsg;
    }

    answer = false;
  }
  if (ErrLogging(debugging)) {
    ErrMsg(debugging) << " DchVolType ANSWER: " << (int) answer << endmsg;
  }
  return answer;
}

const DetMaterial&
DchVolType::material(const TypeCoord* here) const
{
  return *_theMaterial;
}

bool
DchVolType::insideLine(const SurfacePoint& toTest, const SurfacePoint& p1,
    const SurfacePoint& p2) const
{
  return true;
}


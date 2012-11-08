//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id:  $
//
// Description:
//	Class DchPhiSegmVolType
//
// Environment:
//
//
// Author List:
//			originator
//	
//
// Copyright Information:
//
//
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"

//-----------------------
// This Class's Header --
//-----------------------
#include "DchGeom/DchPhiSegmVolType.hh"

//-------------
// C Headers --
//-------------

//---------------
// C++ Headers --
//---------------
#include <cmath>

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
    "$Id: DchPhiSegmVolType.cc 123 2010-04-29 14:41:45Z stroili $";
const double _tolerance = 1.0e-10;

//		----------------------------------------
// 		-- Public Function Member Definitions --
//		----------------------------------------

//----------------
// Constructors --
//----------------
DchPhiSegmVolType::DchPhiSegmVolType(const char* name, double rmin, double rmax, double zmin,
    double zmax, double phi0, double solidDeltaPhi, double hollowDeltaPhi,
    const DetMaterial* theMat, int idnum) :
    DchVolType(name, rmin, rmax, zmin, zmax, theMat, idnum),
      _phi0(phi0), _solidDeltaPhi(solidDeltaPhi), _hollowDeltaPhi(hollowDeltaPhi)
{
  _deltaPhi = _solidDeltaPhi + _hollowDeltaPhi;
  double nPhi = CLHEP::twopi/_deltaPhi;
  _nRot = ceil(nPhi);
  _isPhiSegmented=true;
}

//--------------
// Destructor --
//--------------
DchPhiSegmVolType::~DchPhiSegmVolType()
{
}

////-------------
//// Selectors --
////-------------
//bool
//DchPhiSegmVolType::insideLimitsOf(int side, const SurfacePoint& thisPoint) const
//{
//  bool answer = true;
//  size_t i = side;
//  double z;
//
//  if (side >= sides()->size()) {
//    ErrMsg(error) << " DchPhiSegmVolType: wrong surface " << side << endmsg;
//    return false;
//  }
//  HepPoint point = ((*sides())[i])->spacePoint(thisPoint);
//  double radius = sqrt(point.x() * point.x() + point.y() * point.y());
//
//  if (ErrLogging(debugging)) {
//    ErrMsg(debugging) << "DchPhiSegmVolType: " << point << "\n" << " radius: "
//        << radius << " r diff: " << radius - _rmin << " " << radius - _rmax
//        << " z: " << point.z() << " z diff: " << point.z() - _zmin << " "
//        << point.z() - _zmax << endmsg;
//  }
//
//  z = point.z();
//  double phi=atan2(point.y(),point.x());
//  phi+=(phi<0.0) ? CLHEP::twopi : 0.0;
//  //phi-=_phi0
//  double relativePhi = fmod(phi,_deltaPhi);
//  //std::cout<<"point "<<point<<std::endl;
//  //std::cout<<"point phi "<<phi<<" relativePhi "<<relativePhi<<" solidDeltaPhi "<<_solidDeltaPhi<<std::endl;
//  if (((fabs(z - _zmin) < _tolerance || fabs(z - _zmax) < _tolerance)
//      && (radius > _rmin && radius < _rmax)) || ((fabs(radius - _rmin)
//      < _tolerance || fabs(radius - _rmax) < _tolerance) && (z > _zmin && z
//      < _zmax)) && ((relativePhi-_solidDeltaPhi)<_tolerance) ) {
//          //std::cout<<"Discrete material hit "<<std::endl;
//    answer = true;
//  } else {
//    if (ErrLogging(debugging)) {
//      ErrMsg(debugging) << "DchPhiSegmVolType: " << (fabs(z - _zmin) < _tolerance)
//          << "\t" << (fabs(z - _zmax) < _tolerance) << "\t" << (radius > _rmin)
//          << "\t" << (radius < _rmax) << "\t" << (fabs(radius - _rmin)
//          < _tolerance) << "\t" << (fabs(radius - _rmax) < _tolerance) << "\t"
//          << (z > _zmin) << "\t" << (z < _zmax) << "\t"<<relativePhi-_solidDeltaPhi<< endmsg;
//    }
//
//    answer = false;
//  }
//  if (ErrLogging(debugging)) {
//    ErrMsg(debugging) << " DchPhiSegmVolType ANSWER: " << (int) answer << endmsg;
//  }
//  return answer;
//}

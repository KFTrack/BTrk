//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id:  $
//
// Description:
//	Class DchPhiSegmVolElem
//
// Environment:
//
//
// Author List:
//	
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
#include "DchGeom/DchPhiSegmVolElem.hh"

//-------------
// C Headers --
//-------------

//---------------
// C++ Headers --
//---------------
#include <assert.h>
#include <algorithm>

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "CLHEP/Geometry/HepPoint.h"
#include "CLHEP/Geometry/Transformation.h"
#include "CLHEP/Geometry/Translation.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"
#include "DetectorModel/DetCylinder.hh"
#include "DetectorModel/DetPlane.hh"
#include "DetectorModel/DetIntersection.hh"
#include "DetectorModel/DetSurface.hh"
#include "DetectorModel/DetVolSideIntersection.hh"
#include "DetectorModel/DetVolumeType.hh"
#include "DetectorModel/Intersection.hh"
//#include "DchGeom/DchVolType.hh"
#include "DchGeom/DchPhiSegmVolType.hh"
#include "ErrLogger/ErrLog.hh"
using std::endl;
using std::ends;
using std::ostream;

//-----------------------------------------------------------------------
// Local Macros, Typedefs, Structures, Unions and Forward Declarations --
//-----------------------------------------------------------------------

//		----------------------------------------
// 		-- Public Function Member Definitions --
//		----------------------------------------

//----------------
// Constructors --
//----------------
DchPhiSegmVolElem::DchPhiSegmVolElem(DchVolType* itsType, const char* name, int id,
    const HepTransformation& theAlignment) :
  DetVolumeElem(itsType, name, id, theAlignment), /*_extSettedStrtStp(false),*/ _debug(false)
{
  ;
}

//--------------
// Destructor --
//--------------
DchPhiSegmVolElem::~DchPhiSegmVolElem()
{
}

int
DchPhiSegmVolElem::intersect(const Trajectory* traj, DetIntersection& dinter) const
{
        if (ErrLogging(debugging)) {
                ErrMsg(debugging) << elementName()<<" DchPhiSegmVolElem: intersect \t" << traj->lowRange() << "\t"
                                << traj->hiRange() << "\t" << traj->range() << "\t" << dinter.pathlen
                                << "\t" << dinter.pathrange[0] << "\t" << dinter.pathrange[1] << endl;
                traj->printAll(std::cout);
                ErrMsg(debugging) << endmsg;
        }
        bool success = false;
        double sEntrance, sExit, sPath;
        sEntrance = sExit = sPath = 0.;
        double path = 0.;

        if ( DetVolumeElem::intersect(traj,dinter) ) {
                if ( insideVolume(traj->position(dinter.pathlen)) ) {

                        if (ErrLogging(debugging)) {
                                ErrMsg(debugging) << "intersection found :" << dinter.pathrange[0] <<" = "<<traj->position(dinter.pathrange[0])<< "\t"
                                                << dinter.pathrange[1] <<" = "<<traj->position(dinter.pathrange[1]) << "\t" << dinter.pathlen <<"\t"<< dinter.pathLength() << endmsg;
                        }
                        success = true;
                }
        }

        return success;
}

bool
DchPhiSegmVolElem::insideVolume(const HepPoint& point) const
{
  DchPhiSegmVolType* volType = (DchPhiSegmVolType*) _dtype;
  //  it does not compile on SUN ....
  //   const DchVolType* volType = static_cast<const DchVolType*> (_dtype);

  //   transform to local coordinate system
  HepPoint locPoint = transform().transTo(point);

  double pointRad = sqrt(locPoint.x() * locPoint.x() + locPoint.y()
      * locPoint.y());

  if (ErrLogging(debugging)) {
          ErrMsg(debugging) <<"insideVolume for vol "<<volType->typeName()<<" is ssegm "<<volType->isPhiSegmented()<<std::endl
                          <<"point "<<point<<" local point "<<locPoint<<std::endl
                          <<"vol rmin "<<volType->rmin()<<" vol rmax "<<volType->rmax()<<" vol zmin "<<volType->zmin()
                          <<" vol zmax "<<volType->zmax()<<endmsg;
  }
  if (pointRad > volType->rmin() && pointRad < volType->rmax() && locPoint.z()
      > volType->zmin() && locPoint.z() < volType->zmax()) {
          if (volType->isPhiSegmented()) {
                    double phi=atan2(locPoint.y(),locPoint.x());
                    phi+=(phi<0.0) ? CLHEP::twopi : 0.0;
                    phi-=volType->phi0();
                    double relativePhi = fmod(phi,volType->totDeltaPhi());
                    if (ErrLogging(debugging)) {
                            ErrMsg(debugging) <<"point phi "<<phi<<" relativePhi "<<relativePhi<<" solidDeltaPhi "<<volType->solidDeltaPhi()<<" Pphi0 "<<volType->phi0()<<endmsg;
                    }
                    if ( fabs(relativePhi)>volType->solidDeltaPhi() ) {
                            return false;
                    }
          }
          return true;
  }

  return false;
}

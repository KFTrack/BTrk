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
//#include "TrkBase/HelixTraj.hh"
#include "ErrLogger/ErrLog.hh"
using std::endl;
using std::ends;
using std::ostream;

const double tolerance = 1.0e-6;

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
        //double sEntrance, sExit, sPath;
        //sEntrance = sExit = sPath = 0.;
        //double path = 0.;

        if ( DetVolumeElem::intersect(traj,dinter) ) {
                //double extraAngle(0.0);
                //if ( insideVolume(traj->position(dinter.pathlen),extraAngle) ) {
                if ( chckIntrsctInAglLmts(traj, dinter) ) {

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
DchPhiSegmVolElem::insideVolume(const HepPoint& point, double &extraDelta, int nRot) const
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
//  std::cout <<"insideVolume for vol "<<volType->typeName()<<" is ssegm "<<volType->isPhiSegmented()<<std::endl
//                  <<"point "<<point<<" local point "<<locPoint<<std::endl
//                  <<"vol rmin "<<volType->rmin()<<" vol rmax "<<volType->rmax()<<" vol zmin "<<volType->zmin()
//                  <<" vol zmax "<<volType->zmax()<<std::endl;

  if (pointRad > (volType->rmin()-tolerance) && pointRad < (volType->rmax()+tolerance)
      && locPoint.z() > (volType->zmin()-tolerance) && locPoint.z() < (volType->zmax()+tolerance) ) {
          if (volType->isPhiSegmented()) {
                    double phi=locPoint.phi();
                    phi+=(phi<0.0) ? CLHEP::twopi : 0.0;
                    phi-=volType->phi0();
                    if (nRot>=volType->maxNRotation()) { nRot=volType->maxNRotation()-1; }
                    else if (nRot<0) { nRot = std::floor( phi/volType->totDeltaPhi() + 0.500001); }
                    double relativePhi = phi-((double)nRot)*volType->totDeltaPhi();
                    if (ErrLogging(debugging)) {
                            ErrMsg(debugging) <<"point phi "<<phi<<" relativePhi "<<relativePhi<<" totDeltaPhi "<<volType->totDeltaPhi()
                                              <<" solidDeltaPhi "<<volType->solidDeltaPhi()<<" Pphi0 "<<volType->phi0()
                                              <<" maxNRotation "<<volType->maxNRotation()<<endmsg;
                    }
//                    std::cout <<"point phi "<<phi<<" relativePhi "<<relativePhi<<" totDeltaPhi "<<volType->totDeltaPhi()
//                              <<" solidDeltaPhi "<<volType->solidDeltaPhi()<<" Pphi0 "<<volType->phi0()
//                              <<" maxNRotation "<<volType->maxNRotation()<<std::endl;
                    extraDelta=0.0;
                    if ( relativePhi>volType->solidDeltaPhi() ) {
                            extraDelta = relativePhi-volType->solidDeltaPhi();
                            return false;
                    }
                    if ( relativePhi<0 ) {
                            extraDelta = relativePhi;
                            return false;
                    }
          }
          return true;
  }
  extraDelta=-999;
  return false;
}

bool
DchPhiSegmVolElem::chckIntrsctInAglLmts(const Trajectory* traj, DetIntersection& dinter) const {
        double extraAngle0(0.0), extraAngle1(0.0);
        DchPhiSegmVolType* volType = (DchPhiSegmVolType*) _dtype;
        if (volType->isPhiSegmented()) {
                double minpath(0.0), maxpath(0.0);
                //double k = fabs(traj->curvature(dinter.pathlen));
                //double vrho = (traj->direction(dinter.pathlen)).perp();
                //std::cout<<"trak dir perp "<<vrho<<" cosDip "<<((HelixTraj*)traj)->cosDip()<<std::endl;
                //double angleTpoPath = 1.0/(vrho*k);
                //double angleTpoPath = vrho/k;
                double phi0=traj->position(dinter.pathrange[0]).phi();
                phi0+=(phi0<0.0) ? CLHEP::twopi : 0.0;
                double phi1=traj->position(dinter.pathrange[1]).phi();
                phi1+=(phi1<0.0) ? CLHEP::twopi : 0.0;
                double phi=phi0-volType->phi0();
                int nRot = phi/volType->totDeltaPhi();
                //int nRot = std::floor( phi/volType->totDeltaPhi() + 0.500001);
                double angleTpoPath = (phi1-phi0);
                if (angleTpoPath<0) { angleTpoPath += CLHEP::twopi; }
                angleTpoPath = (dinter.pathrange[1]-dinter.pathrange[0])/angleTpoPath;
                if (ErrLogging(debugging)) {
                        ErrMsg(debugging) << "intersection before check :" << dinter.pathrange[0] <<" = "<<traj->position(dinter.pathrange[0])<< "\t"
                                        << dinter.pathrange[1] <<" = "<<traj->position(dinter.pathrange[1]) << "\t" << dinter.pathlen <<"\t"<< dinter.pathLength() <<endmsg;
                        ErrMsg(debugging) /*<<"track rad "<<1.0/k*/<<" theta "<<traj->direction(dinter.pathlen).theta()<<" dir mag "<<traj->direction(dinter.pathlen).mag()
                                        <<" angleTpoPath "<<angleTpoPath<<" nRot "<<nRot<<endmsg;
                }
//                std::cout << "intersection before check :" << dinter.pathrange[0] <<" = "<<traj->position(dinter.pathrange[0])<< "\t"
//                                << dinter.pathrange[1] <<" = "<<traj->position(dinter.pathrange[1]) << "\t" << dinter.pathlen <<"\t"<< dinter.pathLength() <<std::endl;
//                std::cout /*<<"track rad "<<1.0/k*/<<" theta "<<traj->direction(dinter.pathlen).theta()<<" dir mag "<<traj->direction(dinter.pathlen).mag()
//                                <<" angleTpoPath "<<angleTpoPath<<" nRot "<<nRot<<std::endl;
                bool inside0(false), inside1(false);
                inside0 = insideVolume(traj->position(dinter.pathrange[0]),extraAngle0,nRot);
                if (!inside0 && extraAngle0>0) {
                        ++nRot;
                        inside0 = insideVolume(traj->position(dinter.pathrange[0]),extraAngle0,nRot);
                }
                inside1 = insideVolume(traj->position(dinter.pathrange[1]),extraAngle1,nRot);

                if (inside0 && inside1) { return true; }
                if (extraAngle0<-998 || extraAngle1<-998 ||
                    fabs(extraAngle0)>volType->solidDeltaPhi() || fabs(extraAngle1)>volType->solidDeltaPhi()) { return false; }
                if (extraAngle0*extraAngle1>0.0) { return false; }
                minpath = dinter.pathrange[0]+fabs(extraAngle0)*angleTpoPath;
                maxpath = dinter.pathrange[1]-fabs(extraAngle1)*angleTpoPath;

                dinter.pathrange[0]=minpath;
                dinter.pathrange[1]=maxpath;
                dinter.pathlen=(minpath+maxpath)*0.5;

                if (ErrLogging(debugging)) {
                        ErrMsg(debugging)<< "intersection after check :" << dinter.pathrange[0] <<" = "<<traj->position(dinter.pathrange[0])<< "\t"
                                        << dinter.pathrange[1] <<" = "<<traj->position(dinter.pathrange[1]) << "\t" << dinter.pathlen <<"\t"<< dinter.pathLength() << endmsg;
                }
//                std::cout<< "intersection after check :" << dinter.pathrange[0] <<" = "<<traj->position(dinter.pathrange[0])<< "\t"
//                                << dinter.pathrange[1] <<" = "<<traj->position(dinter.pathrange[1]) << "\t" << dinter.pathlen <<"\t"<< dinter.pathLength() << std::endl;
                return true;
        } else {
                return insideVolume(traj->position(dinter.pathlen),extraAngle0);
        }
}

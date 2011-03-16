//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DchVolElem.cc 123 2010-04-29 14:41:45Z stroili $
//
// Description:
//	Class DchVolElem
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
#include "BaBar/BaBar.hh"

//-----------------------
// This Class's Header --
//-----------------------
#include "DchGeom/DchVolElem.hh"

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
#include "DchGeom/DchVolType.hh"
#include "ErrLogger/ErrLog.hh"
using std::endl;
using std::ends;
using std::ostream;

//-----------------------------------------------------------------------
// Local Macros, Typedefs, Structures, Unions and Forward Declarations --
//-----------------------------------------------------------------------

#define STEP 4.0  // define gross step size
static const double _eps = 1.0e-5;
const double epsilon = 1.0e-5;
const Trajectory* _oldTraj = 0;
int _counter = 0;
int _samples = 0;
double _entrance = 0.;
double _exitp = 0.;
double _stepSize = 0.;
//		----------------------------------------
// 		-- Public Function Member Definitions --
//		----------------------------------------

//----------------
// Constructors --
//----------------
DchVolElem::DchVolElem(DchVolType* itsType, const char* name, int id,
    const HepTransformation& theAlignment) :
  DetVolumeElem(itsType, name, id, theAlignment), _debug(false)
{
  ;
}

//--------------
// Destructor --
//--------------
DchVolElem::~DchVolElem()
{
}

int
DchVolElem::intersect(const Trajectory* traj, DetIntersection& dinter) const
{
  if (ErrLogging(debugging)) {
    ostream& debs = ErrMsg(debugging);
    debs << "DchVolElem: intersect \t" << traj->lowRange() << "\t"
        << traj->hiRange() << "\t" << traj->range() << "\t" << dinter.pathlen
        << "\t" << dinter.pathrange[0] << "\t" << dinter.pathrange[0] << endl;
    traj->printAll(debs);
    debs << endmsg;
  }
  bool success = false;
  double sEntrance, sExit, sPath;
  sEntrance = sExit = sPath = 0.;
  double path = 0.;

  if (ErrLogging(debugging)) {
    ErrMsg(debugging) << "DchVolElem: trajecories: " << _oldTraj << "\t"
        << traj << endmsg;
  }
  if (_oldTraj != traj) { // new track, calculate new path in detector
    // first some initializations
    _oldTraj = traj;
    //     setValidPath(0.);
    _counter = 0;

    std::vector<DetVolSideIntersection> ilist;

    sideIntersect(traj, ilist, dinter.pathrange[0] + epsilon, dinter.pathrange);

    if (ErrLogging(debugging)) {
      ErrMsg(debugging) << "DchVolElem: # of intersections: " << ilist.size()
          << endmsg;
    }

    //  found some intersection
    if (ilist.size() != 0) {

      //  maybe track is generated or dies inside Dch gas volume 
      //  
      if (ilist.size() < 2) {

        if (ErrLogging(debugging)) {
          ostream& outst = ErrMsg(debugging);
          outst << " DchVolElem: NUMBER OF INTERSECTIONS != 2 " << ilist.size()
              << endmsg;
          outst << " intersections:" << endl;
          for (int j = 0; j < ilist.size(); j++) {
            outst << j << ":  entry: " << ilist[j].pathLength() << endl;
          }
          outst << endmsg;
        }

        //  check if it's generated inside the Dch gas volume
        if (insideVolume(traj->position(dinter.pathrange[0]))) {
          _entrance = dinter.pathrange[0];
          _exitp = ilist[0].pathLength();
          if (ErrLogging(debugging)) {
            ErrMsg(debugging) << "track starts in the chamber\t" << _entrance
                << "\t" << _exitp << endl << "___\t" << traj->lowRange()
                << "\t" << traj->hiRange() << "\t" << ilist[0].pathLength()
                << "\t" << dinter.pathrange[1] << endmsg;
          }
          setStep();
          path = 0.5 * (_entrance + _exitp);
          dinter = DetIntersection((DetElem*) this, traj, path, _entrance,
              _exitp);
          dinter.flag[0] = -1;
          dinter.flag[1] = ilist[0].side();
          success = true;
          if (ErrLogging(debugging)) {
            ErrMsg(debugging) << "found track origin inside Dch gas volume"
                << endmsg;
          }
        } else if (insideVolume(traj->position(dinter.pathrange[1]))) {
          _entrance = ilist[0].pathLength();
          _exitp = dinter.pathrange[1];
          if (ErrLogging(debugging)) {
            ErrMsg(debugging) << "track dies in the chamber\t" << _entrance
                << "\t" << _exitp << endl << "...\t" << traj->lowRange()
                << "\t" << traj->hiRange() << "\t" << ilist[0].pathLength()
                << "\t" << dinter.pathrange[0] << endmsg;
          }
          setStep();
          path = 0.5 * (_entrance + _exitp);
          dinter = DetIntersection((DetElem*) this, traj, path, _entrance,
              _exitp);
          dinter.flag[0] = ilist[0].side();
          dinter.flag[1] = -1;
          success = true;
          ErrMsg(debugging) << "found track end inside Dch gas volume"
              << endmsg;
        } else {
          _oldTraj = 0;
          if (ErrLogging(debugging)) {
            ErrMsg(debugging)
                << "the track doesn't intersect the Dch gas volume" << endl
                << "\t--->\t" << dinter.pathrange[0] << "\t"
                << dinter.pathrange[1] << "\t" << dinter.pathlen << "\tfalse"
                << endmsg;
          }
          return false;
        }
      } else if (ilist.size() > 2) {
        ErrMsg(warning) << "more than 2 intersections " << ilist.size()
            << endmsg;
        return false;
      } else {
        if (ErrLogging(debugging)) {
          ErrMsg(debugging) << "standard track" << endmsg;
        }
        unsigned ienter = (ilist[0].pathLength() < ilist[1].pathLength()) ? 0
            : 1;
        unsigned iexit = ienter == 0 ? 1 : 0;
        _entrance = ilist[ienter].pathLength();
        _exitp = ilist[iexit].pathLength();
        setStep();
        path = 0.5 * (_entrance + _exitp);
        dinter
            = DetIntersection((DetElem*) this, traj, path, _entrance, _exitp);
        dinter.flag[0] = ilist[ienter].side();
        dinter.flag[1] = ilist[iexit].side();
        success = true;
      }

      if (ErrLogging(debugging)) {
        ErrMsg(debugging) << " DchVolElem: PATH " << path << " ENTRANCE: "
            << _entrance << " EXIT: " << _exitp << " step: " << _stepSize
            << endmsg;
      }
    } else {

      if (ErrLogging(debugging)) {
        ErrMsg(debugging) << "DchVolElem: no intersections found!!!" << endl
            << "\ttrack starts at position " << traj->position(
            dinter.pathrange[0]) << endl << "\tand ends at position "
            << traj->position(dinter.pathrange[1]) << endmsg;
      }
      if (insideVolume(traj->position(dinter.pathrange[0])) && insideVolume(
          traj->position(dinter.pathrange[1]))) {
        if (ErrLogging(debugging)) {
          ErrMsg(debugging) << " Track originates and dies inside the chamber"
              << endmsg;
        }
        _entrance = dinter.pathrange[0];
        _exitp = dinter.pathrange[1];
        unsigned ienter = (dinter.pathrange[0] < dinter.pathrange[1]) ? 0 : 1;
        unsigned iexit = ienter == 0 ? 1 : 0;
        _entrance = dinter.pathrange[ienter];
        _exitp = dinter.pathrange[iexit];
        setStep();
        path = 0.5 * (_entrance + _exitp);
        dinter
            = DetIntersection((DetElem*) this, traj, path, _entrance, _exitp);
        dinter.flag[0] = -1;
        dinter.flag[1] = -1;
        if (ErrLogging(debugging)) {
          ErrMsg(debugging) << "... track starts and ends in Dch gas volume"
              << endmsg;
        }
        success = true;
      } else {
        if (ErrLogging(debugging)) {
          ErrMsg(debugging) << "DchVolElem: REALLY NO INTERSECTIONS FOUND!!!"
              << endmsg;
        }
        _oldTraj = 0;
        if (ErrLogging(debugging)) {
          ErrMsg(debugging) << "\t--->\t" << dinter.pathrange[0] << "\t"
              << dinter.pathrange[1] << "\t" << dinter.pathlen << "\tfalse"
              << endmsg;
        }
        return false;
      }
    }
  } else {
    if (ErrLogging(debugging)) {
      ErrMsg(debugging) << "DchVolElem: OLD TRACK!" << endmsg;
    }
  }

  if (_counter < _samples) {
    sEntrance = _entrance + _counter * _stepSize;
    sExit = sEntrance + _stepSize;
    sPath = 0.5 * (sEntrance + sExit);

    dinter = DetIntersection((DetElem*) this, traj, sPath, sEntrance, sExit);
    _counter++;
    success = true;
    if (ErrLogging(debugging)) {
      ErrMsg(debugging) << " DchVolElem: STEP # " << _counter << " path: "
          << sPath << " entrance: " << sEntrance << " exit: " << sExit
          << endmsg;
    }
  } else {
    if (ErrLogging(debugging)) {
      ErrMsg(debugging) << " DchVolElem: start " << _counter << "\t"
          << _samples << endmsg;
    }
    _oldTraj = 0;
    dinter = DetIntersection((DetElem*) this, traj, sPath, sEntrance, sExit
        + _eps);

    success = false;
  }
  if (ErrLogging(debugging)) {
    ErrMsg(debugging) << "\tend\t->\t" << dinter.pathrange[0] << "\t"
        << dinter.pathrange[1] << "\t" << dinter.pathlen << "\ttrue\t"
        << success << "\t" << _counter << "\t" << _samples << endmsg;

  }

  if (_counter == _samples) {
    _oldTraj = 0;
  }

  return success;
}

void
DchVolElem::sideIntersect(const Trajectory* traj, std::vector<
    DetVolSideIntersection>& ilist, double distance, double* range) const
{
  if (ErrLogging(debugging)) {
    ErrMsg(debugging) << " sideIntersect:\t" << distance << "\t" << range[0]
        << "\t" << range[1] << endmsg;
  }
  DchVolType* volType = (DchVolType*) _dtype;
  SurfacePoint surfcoord;
  size_t side = 0, nSides = sides()->size();
  for (side = 0; side < nSides; side++) {
    if (ErrLogging(debugging)) {
      ErrMsg(debugging) << "\t___SIDE:\t" << side << endmsg;
    }
    //
    //  Use the initial values to define the starting point and the 
    //  search range
    //
    double flightdist = distance;
    double flightrange[2];
    flightrange[0] = distance;
    flightrange[1] = range[1];
    Intersection intersection(*traj, *(*sides())[side]);
    TrkErrCode iflag = intersection.intersect(flightdist, surfcoord, trkOut,
        flightrange);
    if (ErrLogging(debugging)) {
      ErrMsg(debugging) << "\t...success:\t" << iflag << "\t" << flightrange[0]
          << "\t" << flightrange[1] << "\t" << flightdist << endmsg;
    }
    while (iflag.success()) {
      if (ErrLogging(debugging)) {
        ErrMsg(debugging) << "\tintersection found" << endmsg;
      }
      TwoDCoord detcoord(surfcoord.array());
      if (!(volType->insideLimitsOf(side, surfcoord))) {
        if (ErrLogging(debugging)) {
          ErrMsg(debugging) << "\tNOT inside limits\t->\t" << flightrange[0]
              << "\t" << flightrange[1] << "\t" << flightdist << "\t" << iflag
              << endmsg;
        }
        // not inside the physical boundaries for this plane,
        // go to the next intersection
        flightrange[0] = flightdist + epsilon;
        iflag = intersection.intersect(flightdist, surfcoord, trkOut,
            flightrange);
        if (ErrLogging(debugging)) {
          ErrMsg(debugging) << "\t....success:\t" << iflag << "\t"
              << flightrange[0] << "\t" << flightrange[1] << "\t" << flightdist
              << endmsg;
        }
      } else {
        ilist.push_back(DetVolSideIntersection(flightdist, side, detcoord));
        std::sort(ilist.begin(), ilist.end());
        if (ErrLogging(debugging)) {
          ErrMsg(debugging) << "\tfound intersection\t->\t" << flightdist
              << "\t" << side << endmsg;
        }
        break;
      }
    }
  }

  if (ilist.size() == 1) { // looper, try to find next intersection
    if (ErrLogging(debugging)) {
      ErrMsg(debugging) << "\tfound ONE intersection" << endmsg;
    }
    for (side = 0; side < nSides; side++) {
      if (ErrLogging(debugging)) {
        ErrMsg(debugging) << "\t\tSIDE:\t" << side << endmsg;
      }
      //
      //  Use the initial values to define the starting point and the 
      //  search range
      //
      double flightdist = ilist[0].pathLength() + epsilon;
      double flightrange[2];
      flightrange[0] = flightdist;
      flightrange[1] = range[1];
      Intersection intersection(*traj, *(*sides())[side]);
      TrkErrCode iflag = intersection.intersect(flightdist, surfcoord, trkOut,
          flightrange);
      if (ErrLogging(debugging)) {
        ErrMsg(debugging) << "\t+++success:\t" << iflag << "\t"
            << flightrange[0] << "\t" << flightrange[1] << "\t" << flightdist
            << endmsg;
      }
      while (iflag.success()) {
        TwoDCoord detcoord(surfcoord.array());
        if (!(volType->insideLimitsOf(side, surfcoord))) {
          if (ErrLogging(debugging)) {
            ErrMsg(debugging) << "\tNOT inside limits\t->\t" << flightrange[0]
                << "\t" << flightrange[1] << "\t" << flightdist << "\t"
                << iflag << endmsg;
          }
          // not inside the physical boundaries for this plane,
          // go to the next intersection
          flightrange[0] = flightdist + epsilon;
          iflag = intersection.intersect(flightdist, surfcoord, trkOut,
              flightrange);
          if (ErrLogging(debugging)) {
            ErrMsg(debugging) << "\t+++success:\t" << iflag << "\t"
                << flightrange[0] << "\t" << flightrange[1] << "\t"
                << flightdist << endmsg;
          }
        } else {
          ilist.push_back(DetVolSideIntersection(flightdist, side, detcoord));
          std::sort(ilist.begin(), ilist.end());
          if (ErrLogging(debugging)) {
            ErrMsg(debugging) << "\tfound intersection\t->\t" << flightdist
                << "\t" << side << endmsg;
          }
          break;
        }
      }
    }
  }

  while (ilist.size() > 2) { // too many intersections ... remove some
    if (ErrLogging(debugging)) {
      ErrMsg(debugging) << "\tfound too many intersection\t->\t"
          << ilist.size() << "\t" << endmsg;
    }
    // loop over intersections end remove the one with the bigger pathlength
    ilist.erase(std::max_element(ilist.begin(), ilist.end()));
  }
  std::sort(ilist.begin(), ilist.end());
}

bool
DchVolElem::insideVolume(const HepPoint& point) const
{
  DchVolType* volType = (DchVolType*) _dtype;
  //  it does not compile on SUN ....
  //   const DchVolType* volType = static_cast<const DchVolType*> (_dtype);

  //   transform to local coordinate system
  HepPoint locPoint = transform().transTo(point);

  double pointRad = sqrt(locPoint.x() * locPoint.x() + locPoint.y()
      * locPoint.y());

  if (pointRad > volType->rmin() && pointRad < volType->rmax() && locPoint.z()
      > volType->zmin() && locPoint.z() < volType->zmax()) return true;

  return false;
}

void
DchVolElem::setStep() const
{
  _samples = (int) ((_exitp - _entrance) / STEP);
  if (_samples < 0) {
    ErrMsg(warning) << "negative number of samples... something must be wrong"
        << endmsg;
  }
  if (_samples == 0) _samples = 1;
  _stepSize = (_exitp - _entrance) / _samples;

  if (ErrLogging(debugging)) {
    ErrMsg(debugging) << "DchVolElem: stepsize " << _stepSize << "   samples "
        << _samples << endmsg;
  }

}


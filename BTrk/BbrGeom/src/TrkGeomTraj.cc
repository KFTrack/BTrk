//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkGeomTraj.cc 491 2010-01-13 16:59:16Z stroili $
//
// Description:
//     
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author(s): Steve Schaffner
//
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"
#include "BbrGeom/TrkGeomTraj.hh"
#include "ErrLogger/ErrLog.hh"

TrkGeomTraj::TrkGeomTraj(double lowlim, double hilim) : 
  Trajectory(lowlim, hilim)
{
}

TrkGeomTraj::~TrkGeomTraj()
{
}

void 
TrkGeomTraj::accept(TrkGeomTrajVisitor&) const
{
  ErrMsg(warning) << 
    "TrkGeomTraj: accept() invoked for derived class that has\n"
		  << "not overridden it.  No action taken." << endmsg;
}

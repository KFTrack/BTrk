//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkGeomTraj.hh 491 2010-01-13 16:59:16Z stroili $
//
// Description:
//     Base class for all simple, geometric trajectories (i.e. trajs that 
//     don't describe tracks).  Inherits from Trajectory.  Supports 
//     (via accept() function) Visitor pattern for adding functionality.
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author(s): Steve Schaffner
//
//------------------------------------------------------------------------
#ifndef TRKGEOMTRAJ_HH
#define TRKGEOMTRAJ_HH
#include "BTrk/BbrGeom/Trajectory.hh"

class TrkGeomTrajVisitor;

// Class interface //
class TrkGeomTraj : public Trajectory {

public:
  TrkGeomTraj(double lowlim, double hilim);
  virtual ~TrkGeomTraj();

  virtual void accept(TrkGeomTrajVisitor&) const;
  virtual TrkGeomTraj* clone() const = 0;

private:
  // Preempt
  TrkGeomTraj&   operator= (const TrkGeomTraj&);
  TrkGeomTraj(const TrkGeomTraj &);
};

#endif

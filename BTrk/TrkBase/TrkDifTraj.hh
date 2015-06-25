//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkDifTraj.hh,v 1.18 2003/01/21 12:55:08 raven Exp $
//
// Description:  Abstract base class for Trajectory objects designed to 
//     live inside TrkReps as descriptions of tracks.  Adds to the base 
//     class functions that provide derivatives related to the underlying 
//     independent parameters.  The interface for getting the parameters 
//     is deferred to more derived classes.
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author(s): Steve Schaffner, Dave Brown
//
//------------------------------------------------------------------------

#ifndef TRKDIFTRAJ_HH
#define TRKDIFTRAJ_HH
#include "BTrk/BbrGeom/Trajectory.hh"

class DifPoint;
class DifVector;
class TrkSimpTraj;

// Class interface //
class TrkDifTraj : public Trajectory {
public:
  //**************
  // Constructors, etc.
  //**************
  //       By default, the valid flight distance range is really big
  TrkDifTraj(const double lowlim = -99999.,const double hilim = 99999.);
  virtual ~TrkDifTraj();

  //**************
  // Access 
  //**************
  // DifNumber version of position, direction, 2nd deriv information
  virtual void  getDFInfo(double fltLen, DifPoint& pos, DifVector& direction,
                          DifVector& delDirect) const = 0;
  virtual void  getDFInfo2(double fltLen, DifPoint& pos, DifVector& direction) const;
  // Return locally-valid simple trajectory (complete with parameters) and the
  //    equivalent flight length along it.  Trivial except for piece-wise trajs
  virtual const TrkSimpTraj* localTrajectory(double fltLen, double& localFlt)
    const = 0;

private:
// Preempt
  TrkDifTraj& operator= (const TrkDifTraj&);
  TrkDifTraj(const TrkDifTraj &);
};

#endif

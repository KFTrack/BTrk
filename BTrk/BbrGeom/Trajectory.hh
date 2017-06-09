//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: Trajectory.hh 491 2010-01-13 16:59:16Z stroili $
//
// Description:
//  Abstract base class to describe 3-dimensional trajectories in space.
//  Defines an interface, and provides one data member -- the pathlength
//  range for which a trejectory object is valid (infinite by default).
//  In all cases, "flightlength" = 3-d pathlength along traj.
//
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author(s): Dave Brown, Steve Schaffner
//
//------------------------------------------------------------------------
#ifndef TRAJECTORY_HH
#define TRAJECTORY_HH

//----------------
// BaBar header --
//----------------

#include "BTrk/BbrGeom/HepPoint.h"
// class HepPoint;
#include <iosfwd>
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Vector/ThreeVector.h"

class Trajectory {
   public:
    //*******************
    // Constructors, etc.
    //*******************
    //       By default, the valid flight distance range is really big
    Trajectory(double lowlim, double hilim);
    virtual ~Trajectory();

    //**********
    //  Access
    //**********
    // As a function of the flight distance
    virtual HepPoint position(double) const = 0;
    virtual CLHEP::Hep3Vector direction(double) const = 0;
    virtual CLHEP::Hep3Vector delDirect(double) const = 0;  // 2nd deriv. WRT pathlen
    virtual double curvature(double) const = 0;             // = |delDirect|
                                                            // For more efficient calling:
    virtual void getInfo(double fltLen, HepPoint& pos, CLHEP::Hep3Vector& direction) const = 0;
    virtual void getInfo(double fltLen,
                         HepPoint& pos,
                         CLHEP::Hep3Vector& direction,
                         CLHEP::Hep3Vector& delDirect) const = 0;

    // How far can you go using given approximation in direction (+/- 1) pathDir
    //      before error > tolerance?  Only the sign of pathDir matters; magnitude
    //      must be 1.  Returned distance is >= 0.
    double distTo0thError(double s, double tol, int pathDir) const;
    virtual double distTo1stError(double s, double tol, int pathDir) const = 0;
    virtual double distTo2ndError(double s, double tol, int pathDir) const = 0;

    //  CopyOf function
    virtual Trajectory* clone() const = 0;

    // Range of valid flight distances:
    bool validFlightDistance(double f, double tolerance = 0.0) const;
    virtual void setFlightRange(double newrange[2]);
    double lowRange() const;
    double hiRange() const;
    double range() const;
    //  Print functions

    virtual void print(std::ostream& os) const;
    virtual void printAll(std::ostream& os) const;
    //**************
    // End interface
    //**************
   protected:
    Trajectory& operator=(const Trajectory&);
    double flightrange[2];  // validity range for the flight distance parameter
};

// inline functions

inline bool Trajectory::validFlightDistance(double f, double tol) const {
    return f >= flightrange[0] - tol && f <= flightrange[1] + tol;
}
inline double Trajectory::lowRange() const {
    return flightrange[0];
}
inline double Trajectory::hiRange() const {
    return flightrange[1];
}
inline double Trajectory::range() const {
    return hiRange() - lowRange();
}
#endif

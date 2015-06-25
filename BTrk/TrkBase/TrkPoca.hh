//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkPoca.hh,v 1.15 2006/03/25 15:15:55 brownd Exp $
//
// Description:
//     Calculate the point of closest approach between two trajectories, 
//     or between a trajectory and a point.    
//     Calculates (in ctor) the distance and the flight lengths alongs 
//     the trajectory or trajectories; calculated values are obtained 
//     through accessors.  "Precision" is maximum allowed error on distance 
//     (in cm).  The input flightlengths are used as a starting point; the 
//     code will find the point-of-closest-approach that is closest to that
//     point.  (A good starting point also reduces CPU time.)
//     Note that distance is a signed quantity for two trajectories.
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author(s): Steve Schaffner; core algorithm stolen from Art Snyder.
//
//------------------------------------------------------------------------
#ifndef TRKPOCA_HH
#define TRKPOCA_HH
class TrkDifTraj;
#include "BTrk/TrkBase/TrkPocaBase.hh"

// Class interface //
class TrkPoca : public TrkPocaBase {
public:
  TrkPoca(const Trajectory& traj1, double flt1,
          const Trajectory& traj2, double flt2,
          double precision=1.e-5);
// same, allowing optionally to restrict POCA to the allowed flight
// range of either or both input trajectories.
  TrkPoca(const Trajectory& traj1, double flt1, bool restrictfltlen1,
          const Trajectory& traj2, double flt2, bool restrictfltlen2,
          double precision=1.e-5);

  TrkPoca(const Trajectory& traj,  double flt,
          const HepPoint& pt, double precision=1.e-5);

  TrkPoca(const Trajectory& traj,  double flt, bool restrictfltlen,
          const HepPoint& pt, double precision=1.e-5);

  TrkPoca();

  TrkPoca(const TrkPoca& other);

  TrkPoca& operator = (const TrkPoca& other);

  virtual ~TrkPoca();

  virtual TrkPoca* clone() const;
  
  virtual double doca() const;                   // distance of closest approach
  /*
  // The following inherited functions are also available:
    const TrkErrCode& status() const;      // did the calculation succeed?
    double flt1() const;                   // path length on traj 1 @ poca
    double flt2() const;
    double precision();                    // In case anyone wants to know:
  */

private:
  double _doca;

  // private functions
  void calcDist(const Trajectory& traj1, const Trajectory& traj2);
};


#endif

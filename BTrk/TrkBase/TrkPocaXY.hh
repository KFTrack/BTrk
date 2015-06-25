//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkPocaXY.hh,v 1.8 2006/03/25 15:15:56 brownd Exp $
//
// Description:
//     Calculate the point of closest approach between two trajectories 
//     (TrkDifTraj only) or between a trajectory and a point in the XY plane. 
//     Calculates (in ctor) the distance and the flight lengths along
//     the trajectory or trajectories;
//     The input flightlengths are used as a starting point; 
//     (A good starting point also reduces CPU time.)
//     Note that distance is a signed quantity for two trajectories.
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author(s): Selenia Dittongo                    Univ. Ferrara
//
//     Dave Brown 3/31/03 
//     Re-implemented to use the TrkPoca interface and implementation
//------------------------------------------------------------------------
#ifndef TRKPOCAXY_HH
#define TRKPOCAXY_HH

#include "BTrk/TrkBase/TrkPocaBase.hh"
#include "BTrk/TrkBase/TrkDifTraj.hh"

class TrkPocaXY : public TrkPocaBase {
public:
  
//  Replace the original constructors with the equivalent, standard TrkPoca interface
//
//  TrkPocaXY(const TrkDifTraj& traj1, const double& fltl1, 
//	    const TrkDifTraj& traj2, const double& fltl2);
//
//  TrkPocaXY(const TrkDifTraj& traj, const double& fltl, const HepPoint& pt);
  TrkPocaXY(const Trajectory& traj,  double flt,
	    const HepPoint& pt, double precision=1.0e-4);

  TrkPocaXY(const Trajectory& traj1,  double flt1,
	    const Trajectory& traj2,  double flt2,
	    double precision=1.0e-4);

  TrkPocaXY(const TrkPocaXY& other);


  virtual ~TrkPocaXY();

  virtual TrkPocaXY* clone() const; // returns ownership

  virtual double doca() const;
 
  inline double docaXY() const;             // distance of closest approach in XY plane
// provide the following for backwards compatibility.  Code copying is really an ugly thing.

  double fltl1() const{ return flt1(); }
  double fltl2() const{ return flt2(); }
                  
private:

  double _docaxy;

  void interLineCircle
  (const double& m, const double& q,
   const double& xc, const double& yc, const double& radius,
   double& xint1,  double& yint1, double& xint2,  double& yint2);

  void interTwoLines
  (const double& m1, const double& q1, const double& m2, const double& q2, 
   double& xint,  double& yint);

  void interTwoCircles
  (const double& xc1, const double& yc1, const double& r1,
   const double& xc2, const double& yc2, const double& r2,
   double& xint1,  double& yint1, double& xint2,  double& yint2);

};

// Inlined functions
double TrkPocaXY::docaXY() const {return _docaxy;}

#endif


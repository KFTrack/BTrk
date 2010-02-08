//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkDifPoca.hh,v 1.9 2006/03/25 15:15:55 brownd Exp $
//
// Description:
//     Calculate the point of closest approach between two trajectories.  
//     Calculates (in ctor) the distance and the flight lengths alongs 
//     the two trajectories; calculated values are obtained through accessors.
//     "Precision" is maximum allowed error on distance (in cm).
//     The input flightlengths are used as a starting point; the 
//     code will find the point-of-closest-approach that is closest to that
//     point.  (A good starting point also reduces CPU time.)
//     Note that distance is a signed quantity for two trajectories.
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author(s): Steve Schaffner; core algorithm stolen from Art Snyder
//
//------------------------------------------------------------------------
#ifndef TRKDIFPOCA_HH
#define TRKDIFPOCA_HH
#include "TrkBase/TrkPocaBase.hh"
#include "difAlgebra/DifNumber.hh"

class TrkDifTraj;

// Class interface //
class TrkDifPoca : public TrkPocaBase {

public:
  TrkDifPoca(const TrkDifTraj& traj1, double flt1, 
	     const Trajectory& traj2, double flt2, 
	     double precision=1.e-5);
  TrkDifPoca(const TrkDifTraj& traj,  double flt, 
             const HepPoint& pt, double precision=1.e-5);

  TrkDifPoca(const TrkDifPoca& other);
  virtual ~TrkDifPoca();

  virtual TrkDifPoca* clone() const; // returns ownership

  inline const DifNumber& difDoca() const;       // distance of closest approach
  virtual double doca() const;
  inline const HepVector derivs() const;     //deriv of doca w/r/t track params
  inline void fetchDerivs(HepVector&) const;   // faster version

  // The following inherited functions are also available:
  //const TrkErrCode& status() const;      // did the calculation succeed?
  //double flt1() const;                   // path length on traj 1 @ poca
  //double flt2() const;
  //double precision();                    // In case anyone wants to know:
  
private:	
  DifNumber _doca;

  void calcDist(const TrkDifTraj& traj1, const Trajectory& traj2);
  void calcDist(const TrkDifTraj& , const HepPoint& );

  // Preempt 
  // TrkDifPoca&   operator= (const TrkDifPoca&);
  //  TrkDifPoca(const TrkDifPoca &);
};

// Inlined functions
const DifNumber& TrkDifPoca::difDoca() const                     {return _doca;}
const HepVector TrkDifPoca::derivs() const         {return _doca.derivatives();}
void TrkDifPoca::fetchDerivs(HepVector& dv) const {_doca.fetchDerivatives(dv);}
#endif

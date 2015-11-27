//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkHit.hh,v 1.50 2007/09/24 21:56:27 gapon Exp $
//
// Description:
//  Abstract base class for reconstruction.  Provides common interface 
//   for any tracking hit for fitters, etc.   
//   Stable and semi-stable quantities are cached: weight, flight and 
//   length.  Residual, deltaChi (normalized resid), and derivs of delChi 
//   w/r/t track params.  The flight and hit lengths
//   are ONLY updated whenever updateMeasurement is called
//
//   Also contains flags that tell whether is (or can be) active (i.e. 
//   actually used in fit).
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author(s): Steve Schaffner
//
// Modified 3-july-97 by Leon Rochester to add '<' & '==' operators for 
// Roguewave Sorted Vector
//------------------------------------------------------------------------

#ifndef TRKHITONTRK_HH
#define TRKHITONTRK_HH
#include "BTrk/TrkBase/TrkParticle.hh"
#include "BTrk/TrkBase/TrkPoca.hh"
#include <iostream>
#include <functional>

#include "CLHEP/Matrix/Vector.h"
class TrkRep;
class TrkDifTraj;
class Trajectory;
class TrkErrCode;
class TrkDifPoca;
class TrkPoca;

class TrkHitUpdater;
namespace TrkBase { namespace Functors {
   class updateMeasurement;
   class setActive;
   class setParent;
} };

class TrkHit {
// allow TrkRep to set activity
  friend class TrkRep;
public:
  typedef std::unary_function<TrkHit,bool> predicate_type;
  //****************
  // Constructors and such
  //****************
  TrkHit(); 
  virtual ~TrkHit();
public:

  //****************
  // Accessors -- current state of hit
  //****************
  const TrkRep* getParentRep() const                       {return _parentRep;}
  TrkParticle const& particleType() const;
  const TrkDifTraj* trkTraj() const { return _trkTraj;}

  inline bool isActive() const;    // false => leave out of current fit calc
  double hitRms() const                                    {return _hitRms;}
  double weight() const;  
  double fltLen() const                                     {return _trkLen;}
  double hitLen() const                                     {return _hitLen;}
// ambiguity functions.  These are implemented here as no-ops, and
  virtual int ambig() const;
  virtual void setAmbig(int newambig);

  bool operator==(const TrkHit&) const;
  bool operator< (const TrkHit& rhs) const { return fltLen()<rhs.fltLen();}

  virtual const Trajectory* hitTraj() const = 0;

// test whether residual information is present
  bool hasResidual() const { return _poca.status().success(); }
  TrkPoca const& poca() const { return _poca; }

  //getFitStuff: returns derivs and deltaChi (based on current state of track)
  //updateMeasurement: update internal representation, weight/sigma

  TrkErrCode getFitStuff(CLHEP::HepVector& derivs, double& deltaChi) const;
  TrkErrCode getFitStuff(double& deltaChi) const;

// allow updating POCA, as this can be used in other contexts (ie material model)

  TrkErrCode updatePoca(const TrkDifTraj *trkTraj);

  // return the *external* residual (this calls the Rep, which may call
  // down to the internal residual -- the Rep is responsible for computing
  // this quantity. In the case of a KalRep, this could be the unbiased 
  // residual (i.e. the one wrt to the track 'without' this hit included in
  // the track). This is implemented fully in this baseclass
  //   NOTE: this form of 'resid' is here for backwards compatibilty 
  //         please use the 'bool resid(double&,double&,bool) const'
  //         version for new code...
  double resid(bool exclude=false) const;
  // This version of 'resid' will also return the 'full' error on the
  // residual. 'full' implies that the error due to the uncertainty in the
  // track parameters is included. In addition, it is capable of returning
  // 'false' if the computation failed (basically because a TrkPoca failed).
  // Only trust the answer if 'true' is returned.
  bool   resid(double &resid, double &residErr, bool exclude=false) const;

  // return the *internal* residual (used to satisfy getFitStuff)
  double residual() const;

  //****************
  // Set values
  //****************
  // 
  void setActivity(bool turnOn);   // this is the other function that directly calls
                                   // non-const members of TrkRep, and as such the
                                   // reason we need a non-const TrkRep *
  void setFltLen(double f)                                      {_trkLen = f;}
  void setHitLen(double h)                               {_hitLen = h;}

  //****************
  // Set values that shouldn't normally be set
  //****************
  // Use *only* if you want to circumvent the standard 
  //   calculation of the quantity in question and set it by hand
  void setHitRms(double newRms)                            {_hitRms = newRms;}

  //****************
  // Printing
  //****************
  virtual void print(std::ostream& ) const;
  virtual void printAll(std::ostream& ) const;

protected:
  TrkRep* _parentRep;
  bool _isActive;
  double _hitRms;
  double _trkLen;
  double _hitLen;
  double _resid;
  const TrkDifTraj *_trkTraj;
  TrkPoca _poca;
  // define tolerance for POCA
  static double _tolerance;
  void setTolerance(double newtol);

protected:
  void setHitResid(double newResid)                        {_resid = newResid;}
  TrkRep* parentRep() const { return _parentRep;}
  virtual TrkErrCode updateMeasurement(const TrkDifTraj* traj) = 0;
private:
  TrkHit&   operator= (const TrkHit&);    // Preempt 
  TrkHit(const TrkHit& hit);  // preempt; use 1st protected ctor
  // FIXME: have special 'friend' functors for each operation
  //        that requires friendship which are friends, and then
  //        arrange it such that only the "allowed" classes can
  //        create one of those functors. 
  friend class TrkHitUpdater;
  friend class TrkBase::Functors::updateMeasurement;
  friend class TrkBase::Functors::setActive;
  friend class TrkBase::Functors::setParent;
// allow friends (essentially reps) to change the activity directly
  TrkHit *setActive(bool active) { _isActive = active; return this; }
  TrkHit *setParent(TrkRep* rep) { _parentRep = rep; return this; }
};

// Inline functions
inline bool TrkHit::isActive() const {return _isActive;}

std::ostream& operator<<(std::ostream& o, const TrkHit& x) ;
#endif

//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkHitOnTrk.hh,v 1.50 2007/09/24 21:56:27 gapon Exp $
//
// Description:
//  Abstract base class for reconstruction.  Provides common interface 
//   for any tracking hit (e.g. SVT or DCH) for fitters, etc.   
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
#include "TrkBase/TrkParticle.hh"
#include "TrkBase/TrkEnums.hh"
#include <iostream>
#include <functional>

#include "CLHEP/Matrix/Vector.h"
class TrkRep;
class TrkDifTraj;
class Trajectory;
class TrkErrCode;
class TrkDifPoca;
class TrkPoca;

class TrkHitOnTrkUpdater;
namespace TrkBase { namespace Functors {
   class updateMeasurement;
   class setActive;
   class setParent;
} };

class TrkHitOnTrk {
  friend class TrkHotList;
// allow TrkRep to set activity
  friend class TrkRep;
public:
  typedef std::unary_function<TrkHitOnTrk,bool> predicate_type;
  //****************
  // Constructors and such
  //****************
  TrkHitOnTrk(double tolerance); 
  virtual ~TrkHitOnTrk();
  virtual TrkHitOnTrk* clone(TrkRep* parentRep, const TrkDifTraj* trkTraj=0) const = 0;
protected:
  TrkHitOnTrk(const TrkHitOnTrk& hitToBeCopied, TrkRep* newRep, const TrkDifTraj* trkTraj=0 );
public:

  //****************
  // Accessors -- current state of hit
  //****************
  const TrkRep* getParentRep() const                       {return _parentRep;}
  TrkParticle const& particleType() const;
  const TrkDifTraj* trkTraj() const { return _trkTraj;}



  inline bool isActive() const;    // false => leave out of current fit calc
  inline bool isUsable() const;    // false => cannot be made active
  inline int usability() const { return _isUsable;}
  inline bool mustUse() const;     // true => cannot be made inactive
  virtual TrkEnums::TrkViewInfo whatView() const = 0;
  virtual unsigned layerNumber() const = 0;
  double hitRms() const                                    {return _hitRms;}
  double weight() const;  
  double fltLen() const                                     {return _trkLen;}
  double hitLen() const                                     {return _hitLen;}
// ambiguity functions.  These are implemented here as no-ops, and
// are overridden where necessary in DchHot
  virtual int ambig() const;
  virtual void setAmbig(int newambig);

  bool operator==(const TrkHitOnTrk&) const;
  bool operator< (const TrkHitOnTrk& rhs) const { return fltLen()<rhs.fltLen();}

  virtual const Trajectory* hitTraj() const = 0;

// test whether residual information is present
  bool hasResidual() const { return _poca != 0; }
// poca STATUS
  TrkErrCode pocaStatus() const;
  const TrkPoca* poca() const { return _poca; }

  //getFitStuff: returns derivs and deltaChi (based on current state of track)
  //updateMeasurement: update internal representation, weight/sigma

  TrkErrCode getFitStuff(HepVector& derivs, double& deltaChi) const;
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

  // timing information; note that this returns in units of SECONDS!!!
  // First, relative to the track time
  virtual bool timeResid(double& resid, double& error) const = 0; 
  // then, in 'absolute' units (relative to the trigger time)
  virtual bool timeAbsolute(double& time,double& error) const = 0;

  //****************
  // Set values
  //****************
  // 
  void setActivity(bool turnOn);   // this is the other function that directly calls
                                   // non-const members of TrkRep, and as such the
                                   // reason we need a non-const TrkRep *
  void setUsability(int usability);         // 0=unusable; 1=usable; 2=must use
                                            // setUsability will call setActivity
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
  int  _isUsable;
  double _hitRms;
  double _trkLen;
  double _hitLen;
  double _resid;
  const TrkDifTraj *_trkTraj;
  TrkPoca *_poca;
  double _tolerance;

protected:
  void setHitResid(double newResid)                        {_resid = newResid;}
  TrkRep* parentRep() const { return _parentRep;}
  void setUsedHit();               // tell underlying hit 
  void setUnusedHit(); 
  virtual TrkErrCode updateMeasurement(const TrkDifTraj* traj) = 0;
private:
  TrkHitOnTrk&   operator= (const TrkHitOnTrk&);    // Preempt 
  TrkHitOnTrk(const TrkHitOnTrk& hit);  // preempt; use 1st protected ctor
  // FIXME: have special 'friend' functors for each operation
  //        that requires friendship which are friends, and then
  //        arrange it such that only the "allowed" classes can
  //        create one of those functors. 
  friend class TrkHitOnTrkUpdater;
  friend class TrkBase::Functors::updateMeasurement;
  friend class TrkBase::Functors::setActive;
  friend class TrkBase::Functors::setParent;
// allow friends (essentially reps) to change the activity directly
  TrkHitOnTrk *setActive(bool active) { _isActive = active; return this; }
  TrkHitOnTrk *setParent(TrkRep* rep) { _parentRep = rep; return this; }
};

// Inline functions
inline bool TrkHitOnTrk::isActive() const {return _isActive;}
inline bool TrkHitOnTrk::isUsable() const {return (_isUsable > 0);}
inline bool TrkHitOnTrk::mustUse()  const {return (_isUsable > 1);}

std::ostream& operator<<(std::ostream& o, const TrkHitOnTrk& x) ;
#endif

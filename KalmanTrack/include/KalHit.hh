// ----------------------------------------------------------------------------
// File and Version Information:
//      $Id: KalHit.hh,v 1.51 2005/10/04 19:34:10 brownd Exp $
//
//  Description:
//  Class to describe a kalman filter hit site.  This puts measurement information
//  into the fit
//
// Copyright Information:
//	Copyright (C) 1996	Lawrence Berkeley Laboratory
//
//  Authors: Dave Brown, 3/7/97
//-----------------------------------------------------------------------------
#ifndef KALHIT_HH
#define KALHIT_HH
#include <assert.h>
#include "KalmanTrack/KalSite.hh"
#include "TrkBase/TrkHitOnTrk.hh"
#include "TrkBase/TrkHitOnTrkUpdater.hh"
#include "TrkBase/TrkEnums.hh"
#include "KalmanTrack/KalWeight.hh"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Matrix/Matrix.h"
//
class TrkDifPieceTraj;
class TrkRecoTrk;
class TrkParams;
class TrkRep;

//
//  Define the class
//
class KalHit : public KalSite, TrkHitOnTrkUpdater {
public:
//  Constructors
  KalHit(const TrkDifPieceTraj*,TrkHitOnTrk*);
// copy constructor
  KalHit(const KalHit&);
// clone function
  KalSite* clone(const KalRep*) const;
//
  virtual ~KalHit();
//
//  Fit functions
//
  bool process(const KalSite*,trkDirection idir);
  bool update(const TrkDifPieceTraj*,double);
// override needsFit, to be able to follow HOT state changes
  virtual bool needsFit(trkDirection idir) const;
//
//  Access
//
  void printAll(std::ostream& os = std::cout) const;
  const TrkHitOnTrk* hitOnTrack() const { return _hot; }
  TrkHitOnTrk* hitOnTrack() { return _hot; }
// the following are now deprecated
  double chisquared(const KalSite*, trkDirection) const;
  double chisquared(const KalParams& params) const;
// chisq including boolean value
  bool chisquared(double& chisq,const KalSite*, trkDirection) const;
// can also compute the chi (signed) of the hit, and the total error (squared!!!)
// on this residual (including the parameter errors). 
// return value is true for OK, false if the calculation fails.
//  The activity state of the hot can optionally be ignored.
  bool chi(const KalParams& params,double& chival,double& chierr2,
	   bool ignoreactive=true) const;
// same as above, except just computing chisquared.
  bool chisquared(double& chisq,const KalParams& params,
		  bool ignoreactive=true) const;
// compute the uncorrected residual (chi) relative to the reference trajectory
  double referenceChi() const { return _residual; }
  virtual unsigned nDof(TrkEnums::TrkViewInfo view=TrkEnums::bothView) const;
  virtual bool isActive() const { return _hot!=0 ? _hot->isActive() : false; }
// delete the hot; normally this is owned outside
  void deleteHOT() { delete _hot; }
  void setActivity(bool turnOn);
// return the flight length change since last update.  This can be useful
// for identifying POCA failures
  double flightLengthChange() const { return _fltdif; }
// override invert
  void invert();
private:
  TrkHitOnTrk* _hot; // hot for thqis site
  KalWeight _hitweight; // hit parameters/weight
  double _residual; // residual vector
  HepVector _linrel; // hit derrivatives
  bool _hotstate; // record the state of the HOT
  double _fltdif; // flightlength difference since last iteration
// private functions
  bool updateCache(const TrkDifPieceTraj*); // update the above elements when necessary
  const HepVector& refvec() const { return
				      localTrajectory()->parameters()->parameter(); }
// clone the hot
  void cloneHot(const KalRep* newrep);
};
#endif

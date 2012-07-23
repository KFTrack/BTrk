// ----------------------------------------------------------------------------
// File and Version Information:
//      $Id: KalConstraint.hh,v 1.14 2005/10/04 19:34:09 brownd Exp $
//
//  Description:
//  Class to describe a direct constraint on the parameters in a Kalman
//  fit.  This site is identical to a hit site, except that the 'residual'
//  is measured directly WRT the track parameters.
//
// Copyright Information:
//	Copyright (C) 1999	Lawrence Berkeley Laboratory
//
//  Authors: Dave Brown, 7/8/99
//-----------------------------------------------------------------------------
#ifndef KALCONSTRAINT_HH
#define KALCONSTRAINT_HH
#include <assert.h>
#include "KalmanTrack/KalSite.hh"
#include "KalmanTrack/KalWeight.hh"
#include "CLHEP/Matrix/SymMatrix.h"

class TrkDifPieceTraj;
class TrkParams;
class TrkRep;

//
//  Define the class
//
class KalConstraint : public KalSite {
public:
//
//  Constructors
//
  KalConstraint(const TrkDifPieceTraj*,const TrkParams&,
		bool* constrainParams, double fltlen);
// copy constructor
  KalConstraint(const KalConstraint&);
// clone function
  KalConstraint* clone(const KalRep*) const;
//
  virtual ~KalConstraint();
//
//  Fit functions
//
  bool process(const KalSite*,trkDirection idir);
// update does nothing for this site
  bool update(const TrkDifPieceTraj*,double);
  bool chisquared(double& chisq,const KalSite*,trkDirection) const;
// direct chisquared calculation
  bool chisquared(double& chisq,const KalParams& params) const;
  unsigned nDof(TrkEnums::TrkViewInfo view=TrkEnums::bothView) const;
//  Access
  void printAll(std::ostream& os = std::cout) const;
// KalConstraint specific
  const KalWeight& constraintWeight() const { return _cweight; }
  const KalParams& constraintParams() const { return _cparams; }
  bool isConstrained(unsigned iparam) const { return _constrain[iparam]; }
// allow the constraint to be manipulated.  Return value tells if anything changed
  bool setConstraint(bool* constrain);
// override isactive
  virtual bool isActive() const;
// can provide the effective trajectory this constraint represents
  const TrkSimpTraj* constraintTrajectory() const;
// override invert
  virtual void invert();
protected:
// subclass constructor
  KalConstraint(const TrkDifPieceTraj*,const TrkParams&,
		bool* constrainParams, double fltlen,
		KalSite::siteType type);
  KalWeight _cweight; // constraint weight (masked)
  KalParams _cparams; // constraint parameters (unmasked)
private:
// both weight and params are needed, since masking makes them no longer
// invertible
  bool* _constrain; // which parameters to constrain
// utility function to compute the weight
  void maskWeight();
// cache of simptraj representing this constraint
  mutable TrkSimpTraj* _traj;
};
#endif


// ------------------------------------------------------------------------------
// File and Version Information:
//      $Id: KalSite.hh,v 1.62 2006/05/16 18:18:31 brownd Exp $
//
//  Description:
//  Class to describe a generic Kalman filter 'site'.  A site is
//  defined as somewhere along the trajectory where something happens.
//  Specific daughter classes inherit from this to describe either
//  material interactions or measurements.  This abstract base class serves
//  as a placeholder in a list of sites used to compute the Kalman
//  filter track fit.  It also holds information used in the non-gaussian
//  kalman extension.
//
// Copyright Information:
//	Copyright (C) 1996	Lawrence Berkeley Laboratory
//
//  Authors: Dave Brown, 12/18/96
//------------------------------------------------------------------------------
//
#ifndef KALSITE_HH
#define KALSITE_HH
#include <iostream>
#include "BaBar/PdtPid.hh"
#include "TrkBase/TrkDirection.hh"
#include "TrkBase/TrkEnums.hh"
#include "KalmanTrack/KalParams.hh"
#include "KalmanTrack/KalWeight.hh"
//
class TrkSimpTraj;
class TrkDifPieceTraj;
class HepPoint;
class KalRep;
class KalHit;
class KalBend;
class KalMaterial;
class KalPairSite;
class KalConstraint;
class KalSmear;
class KalScatter;
class KalBrems;
class KalBetaCons;
//
class KalSite{
//  friend class KalPairRep;

public:
//  Enums
  enum siteType{matSite=0,hitSite,bendSite,endSite,pairSite,constraintSite,smearSite,scatterSite,bremsSite,betaConsSite,unknown};
//  Constructors; 'default'
  KalSite(const siteType);
// copy constructor. 
  KalSite(const KalSite&);
// clone onto a new rep
  virtual KalSite* clone(const KalRep*) const = 0;
// destructor
  virtual ~KalSite();
//  Fit functions
  virtual bool process(const KalSite*,trkDirection) = 0;
// updating
  virtual bool update(const TrkDifPieceTraj*,double momentum) = 0;
//  Basic information
  virtual void printAll(std::ostream& os=std::cout) const;
  void print(std::ostream& os=std::cout) const;
  bool hasFit(trkDirection idir) const {return !needsFit(idir); }
  virtual bool needsFit(trkDirection idir) const {return !_siteflag[idir];}
  const siteType& type() const { return _stype; }
// trajectory information
  double globalLength() const { return _globlen; }
  double localLength() const { return _loclen; }
  const TrkSimpTraj* localTrajectory() const { return _loctraj;}
  const HepVector& localParameters() const { return _lparams; }
// does the site carry information?  If so, how much.  Boolean indicates
// whether the site contributed at all to chisq in the fit
  virtual bool chisquared(double& chsiq,const KalSite*,trkDirection) const;
  virtual unsigned nDof(TrkEnums::TrkViewInfo view=TrkEnums::bothView) const;
// activity of the site (ie, does it have any impact on the fit)
  virtual bool isActive() const;
//  Filtering results: these just return the data member if it's been built,
//  otherwise it tries to build them (a kind of lazy evaluation).
  const KalParams& filterParameters(trkDirection idir) const;
  const KalWeight& filterWeights(trkDirection idir) const;
// Force a site to be invalid in a particular direction
  void invalidateSite(trkDirection idir) {  _siteflag[idir] = false; }
// set the state (parameters and vectors) of TrkSimpTraj
// according to the state of this site (in the specified direction)
  bool setTrajState(trkDirection tdir,TrkSimpTraj* straj) const;
//  Functions for rw vectors
  bool operator == (const KalSite& other) const {
    return _globlen == other._globlen; }
  bool operator < (const KalSite& other) const {
    return _globlen < other._globlen; }
// momentum change on traversing this site in the specified direction
  virtual double momentumChange(trkDirection idir) const { return 0.0; }
// gap incured by this site
  double gap() const { return _gap; }
// merge parameters with another sites: this takes care of (potential) 
// reference point differences.  The matrixOK flag of the output KalParams
// object describes the status.
  void mergeParams(const KalSite* other,KalParams& addparams) const;
// add and subtract parameters, accounting for reference point differences
  void addParams(const KalSite* other,KalParams& addparams) const;
// chisq wrt parameters of another site, optionally only for a subset of parameters
  double parameterChisq(const KalSite* other,bool* tparams=0) const;
// add weights (weight difference has no meaning)
  void addWeights(const KalSite* other,KalWeight& addweights) const;
// safe downcasting
  KalMaterial* kalMaterial() {
    return (_stype == matSite) ? (KalMaterial*) this : 0;
  }
  KalHit* kalHit() {
    return (_stype == hitSite) ? (KalHit*) this : 0;
  }
  KalBend* kalBend() {
    return (_stype == bendSite) ? (KalBend*) this : 0;
  }
  KalPairSite* kalPairSite() {
    return (_stype == pairSite) ? (KalPairSite*) this : 0;
  }
  KalConstraint* kalConstraint() {
    return (_stype == constraintSite) ? (KalConstraint*) this : 0;
  }
//
  const KalMaterial* kalMaterial() const {
    return (_stype == matSite) ? (const KalMaterial*) this : 0;
  }
  const KalHit* kalHit() const {
    return (_stype == hitSite) ? (const KalHit*) this : 0;
  }
  const KalBend* kalBend() const {
    return (_stype == bendSite) ? (const KalBend*) this : 0;
  }
  const KalPairSite* kalPairSite() const {
    return (_stype == pairSite) ? (const KalPairSite*) this : 0;
  }
  const KalConstraint* kalConstraint() const {
    return (_stype == constraintSite) ? (const KalConstraint*) this : 0;
  }
  const KalSmear* kalSmear() const {
    return (_stype == smearSite) ? (const KalSmear*) this : 0;
  }
  const KalScatter* kalScatter() const {
    return (_stype == scatterSite) ? (const KalScatter*) this : 0;
  }
  const KalBrems* kalBrems() const {
    return (_stype == bremsSite) ? (const KalBrems*) this : 0;
  }
  const KalBetaCons* kalBetaCons() const {
    return (_stype == betaConsSite) ? (const KalBetaCons*) this : 0;
  }
// simple utility to decide if a site is outwards or inwards of myself
  trkDirection relativeDirection(const KalSite* other) const {
    return (other->globalLength() > globalLength()) ? trkOut : trkIn; }
// reset the fit flags
  void reset();
  void reset(trkDirection);
  virtual void invert(); // invert the direction (sign of flightlength) for this site
// set the gap
  void setGap(double gap) { _gap = gap; }
private:
//  Filtering results; these are lazy-evaluated caches
  KalParams _params[2];
  KalWeight _weight[2];
  bool _siteflag[2]; // flag to keep track of the condition of the site
  siteType _stype;
  double _globlen; // global trajectory flight length of this site
  double _loclen; // local (simple traj) flight length of this site
  const TrkSimpTraj* _loctraj; // local trajectory
  HepVector _lparams; // local trajectory parameters
  double _gap; // gap incured in trajectory by this site
protected:
// equivalence
  KalSite& operator = (const KalSite& other);
// default constructor
  KalSite();
//  Access the filtering results by reference
  KalParams& params(trkDirection tdir) const {
    return (KalParams&)_params[tdir]; }
  KalWeight& weight(trkDirection tdir) const {
    return (KalWeight&)_weight[tdir]; }
//  Helpful functions for internal manipulations
  void setFit(trkDirection idir,bool sval = true) {_siteflag[idir] = sval;}
// set the trajectory parameters
  virtual bool setTraj(const TrkDifPieceTraj*,double globlen);
// update the parameters: this can handle the case that the reference points
// aren't the same
  void copyParams(const KalSite* other,trkDirection idir);
// same for weights
  void copyWeights(const KalSite* other,trkDirection idir);
// copy whatever is active in the site
  void copySite(const KalSite* other,trkDirection idir);
// test if the site has valid weight or parameters
  bool validSite(trkDirection idir) const;
// generic parameter-space processing.
  bool processParams(const KalSite* prevsite,
		     trkDirection tdir,
		     const KalParams& transport);
// generic weight-space processing
  bool processWeight(const KalSite* prevsite,
		     trkDirection tdir,
		     const KalWeight& hweight);
// update parameters for a given (large) momentum fraction change
  void processDeltaP(KalParams& params,double dpfract) const;
};
#endif


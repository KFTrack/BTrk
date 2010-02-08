// ------------------------------------------------------------------------------
// File and Version Information:
//      $Id: KalMaterial.hh,v 1.50 2006/04/24 18:53:06 brownd Exp $
//
//  Description:
//  Class to describe a kalman filter material interaction site.
//  This models all material effects which do not change the particle
//  count (essentially multiple scattering and energy loss).
//
// Copyright Information:
//	Copyright (C) 1996	Lawrence Berkeley Laboratory
//
//  Authors: Dave Brown, 12/18/96
//------------------------------------------------------------------------------
#ifndef KALMATERIAL_HH
#define KALMATERIAL_HH
#include <assert.h>
#include "KalmanTrack/KalSite.hh"
#include "DetectorModel/DetIntersection.hh"
#include "CLHEP/Matrix/Matrix.h"
class DetElem;
class TrkDifPieceTraj;
class HepPoint;
//
//  Define the class
//
class KalMaterial : public KalSite {
public:
//
//  Constructors; note that the intersection trajectory MUST be a TrkDifPieceTraj
//
  KalMaterial(const DetIntersection&,const TrkDifPieceTraj*,
	      double momentum,PdtPid::PidType pid);
// copy constructor
  KalMaterial(const KalMaterial&);
// clone operator
  KalMaterial* clone(const KalRep*) const;
//  destructor
  virtual ~KalMaterial();
//  Fit function
  bool process(const KalSite*,trkDirection);
// update
  virtual bool update(const TrkDifPieceTraj*,double momentum);
// momentum change
  virtual double momentumChange(trkDirection idir) const;
  double energyChange(trkDirection idir) const;
// materials can be inactive
  virtual bool isActive() const { return _active; }
  void setActivity(bool active);
//
//  Access
//
  void printAll(std::ostream& os = std::cout) const;
  const DetElem* detElem() const { return _dinter.delem; }
  const KalParams& transport() const { return _transport; }
  const DetIntersection& detIntersection() const { return _dinter; }
  double momFraction() const { return _pfract; }
  double momFractionRMS() const { return _pfractrms; }
  double deflectRMS() const { return _deflectrms; }
// momentum on the specificed _side_ of the material
  double momentum(trkDirection idir) const {
    return idir == trkOut ? _momentum : _momentum + momentumChange(trkIn); }
  PdtPid::PidType particle() const { return _partid; }
// number of radation lengths through this site
  double radiationFraction() const;
// material doesn't contribute to chisquared
// override invert
  virtual void invert();
private:
  PdtPid::PidType _partid; // define the particle type (and hence mass)
  DetIntersection _dinter; // intersection with this element
  KalParams _transport; // parameter transport (column matrix, P' = P + DP) for small dE/dx
//  HepMatrix _covrot; // covariance similarity transform for energy loss
  void updateCache(const TrkDifPieceTraj*); // update the above elements when necessary
  double _momentum; // approximate momentum of track in the middle of the site
  bool _active; // is the site active or not
  HepSymMatrix _scatter; // scattering effects, normalized to unit angle
  HepSymMatrix _eloss; // energy loss effects, normalized to unit dedx
  HepMatrix _pderiv; //derivative of momentum (energy) loss
  double _pfract;// cache the fractional momentum change
  double _pfractrms; // same for rms
  double _deflectrms; // cache the scattering sigma
// override setting the trajectory to deal with the detectorintersection
  virtual bool setTraj(const TrkDifPieceTraj*,double globlen);
// reset the PID
  void setPID(PdtPid::PidType newtype) { _partid = newtype; }
};
#endif

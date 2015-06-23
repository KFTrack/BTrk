// ------------------------------------------------------------------------------
// File and Version Information:
//      $Id: KalBrems.hh,v 1.3 2006/04/27 17:59:15 brownd Exp $
//
//  Description:
//  KalSite subclass representing a bremsstrahlung, where the photon is explicitly
//  reconstructed as a neutral cluster.
//
// Copyright Information:
//	Copyright (C) 2006	Lawrence Berkeley Laboratory
//
//  Authors: Dave Brown, 3/7/2006
//------------------------------------------------------------------------------
#ifndef KALBREMS_HH
#define KALBREMS_HH
#include <assert.h>
#include "KalmanTrack/KalSite.hh"
#include "DetectorModel/DetIntersection.hh"
#include "CLHEP/Matrix/Matrix.h"
class DetElem;
class TrkDifPieceTraj;
class HepPoint;
class AbsRecoCalo;
//
//  Define the class
//
class KalBrems : public KalSite {
public:

// construct from cluster
  KalBrems(const TrkDifPieceTraj*,
	   const AbsRecoCalo* cluster,
	   double momentum,
	   PdtPid::PidType pid);
// copy constructor
  KalBrems(const KalBrems&);
// clone operator
  KalBrems* clone(const KalRep*) const;
//  destructor
  virtual ~KalBrems();
//  Fit function
  bool process(const KalSite*,trkDirection);
// update
  virtual bool update(const TrkDifPieceTraj*,double momentum);
// cluster 
  double clusterEnergy() const;
  double clusterEnergyErr() const;
  virtual double momentumChange(trkDirection idir) const;
// bremss will by default be inactive except on electron sites, but allow overriding this
  virtual bool isActive() const { return _active; }
  void setActivity(bool active);
//
//  Access
//
  void printAll(std::ostream& os = std::cout) const;
  const AbsRecoCalo* recoCalo() const { return _cluster; }
// momentum _after_ the brems loss
  double momentum() const { return _momentum; }
  PdtPid::PidType particle() const { return _partid; }
// brems doesn't contribute to chisquared
private:
  const AbsRecoCalo* _cluster; // reference my cluster
  PdtPid::PidType _partid; // define the particle type (and hence mass)
  double _momentum; // approximate momentum of track in the middle of the site
  bool _active; // is the site active or not
  HepMatrix _pderiv; // derivative of momentum change effect on parameters
  double _gamphi; // phi direction of the gamma
  double _gamE; // cluster energy
  double _gamEerr; // cluster energy error

  void updateCache(const TrkDifPieceTraj*); // update the above elements when necessary
// reset the PID
  void setPID(PdtPid::PidType newtype);
// locate tangent point
  double tangent(const TrkDifPieceTraj*) const;
};
#endif

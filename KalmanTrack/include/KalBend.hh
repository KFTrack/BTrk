// ------------------------------------------------------------------------------
// File and Version Information:
//      $Id: KalBend.hh,v 1.23 2006/04/24 18:53:06 brownd Exp $
//
//  Description: KalBend
//  KalBend is a KalSite subclass for describing the effect of
//  non-uniform magnetic field upon a track. The hard work (the
//  path integral of the non-solenoidal part of the field) must
//  be done outside this class.
//
// Copyright Information:
//	Copyright (C) 1998	Lawrence Berkeley Laboratory
//
//  Authors: Dave Brown, 3/30/98
//------------------------------------------------------------------------------
#ifndef KALBEND_HH
#define KALBEND_HH

#include "BField/BFieldIntegrator.hh"
#include "KalmanTrack/KalSite.hh"
#include "CLHEP/Vector/ThreeVector.h"

class KalBend : public KalSite {
public:
// construct from the field integrator, the trajectory, range,
// and momentum (direction from the middle of the range)
  KalBend(const BFieldIntegrator&,
	  const TrkDifPieceTraj*,BFieldIntRange const& range,
	  double momentum, int charge);
// copy constructor
  KalBend(const KalBend& other);
// clone operator
  KalBend* clone(const KalRep*) const;
// destructor
  virtual ~KalBend();
// fit processing
  bool process(const KalSite*,trkDirection);
// updating
  bool update(const TrkDifPieceTraj*,double momentum);
// bend adds nothing to chisquared
// KalBend specific; fractional components of the momentum change
  double deltaTheta() const {
    return _charge*_delmom.dot(_thetahat)/_momentum; }
  double deltaPhi() const {
    return _charge*_delmom.dot(_phihat)/_momentum; }
  double deltaMomentum() const { return _delmom.mag(); }
  double momentum() const { return _momentum; }
  double lowRange() const { return _range._slo; }
  double hiRange() const { return _range._shi; }
  double range() const { return _range.range(); }
// informative print
  virtual void printAll(std::ostream& os=std::cout) const;
// override invert
  virtual void invert();
private:
  const BFieldIntegrator& _integrator; // integrator
  BFieldIntRange _range; // integral range
  double _momentum; // momentum magnitude (unsigned)
  int _charge; // particle charge
  Hep3Vector _delmom; // momentum change
  KalParams _transport; // transport vector (no process noise)
// update the cache
  void updateCache(const TrkDifPieceTraj*);
// unit vectors, derived from the momentum direction
  Hep3Vector _momhat;
  Hep3Vector _thetahat;
  Hep3Vector _phihat;
  double midpoint() const { return _range._smid; }
};
#endif

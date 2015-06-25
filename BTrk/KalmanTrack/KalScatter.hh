// ------------------------------------------------------------------------------
// File and Version Information:
//      $Id: KalScatter.hh,v 1.3 2006/04/24 18:53:07 brownd Exp $
//
//  Description:
//  Class to describe a hard scatter at a material interaction site.
//
// Copyright Information:
//	Copyright (C) 2005     Lawrence Berkeley Laboratory
//
//  Authors: Dave Brown, 2/14/05
//------------------------------------------------------------------------------
#ifndef KALSCATTER_HH
#define KALSCATTER_HH
#include <assert.h>
#include "BTrk/KalmanTrack/KalSite.hh"
#include "BTrk/DetectorModel/DetIntersection.hh"
#include "CLHEP/Matrix/Matrix.h"
class DetElem;
class KalMaterial;
//
//  Define the class
//
class KalScatter : public KalSite {
public:
//  Construct from KalMaterial, with scattering factor
  KalScatter(const KalMaterial& km,double factor);
// copy constructor
  KalScatter(const KalScatter&);
// clone operator
  KalScatter* clone(const KalRep*) const;
//  destructor
  virtual ~KalScatter();
//  Fit function
  virtual bool process(const KalSite*,trkDirection);
  virtual bool update(const TrkDifPieceTraj* newtraj,double newmom);
//
//  Access
//
  void printAll(std::ostream& os = std::cout) const;
  const KalParams& transport() const { return _transport; }
  double deflectRMS() const { return _deflectrms; }
// scatter doesn't contribute to chisquared
private:
  KalParams _transport; // parameter transport 
  double _deflectrms; // scattering sigma
  double _pfractrms; // same for energy loss
  CLHEP::HepSymMatrix _scatter; // scattering effects, normalized to unit angle
//
  void updateCache(const TrkDifPieceTraj*,double fltlen); // update the above elements when necessary
};
#endif


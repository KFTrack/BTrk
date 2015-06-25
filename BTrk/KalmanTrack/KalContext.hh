//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: KalContext.hh 103 2010-01-15 12:12:27Z stroili $
//
// Description:
//      class KalContext
//      This class collects the configuration and conditions data
//	needed to perform the Kalman fit
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Copyright Information:
//	Copyright (C) 1997	Lawrence Berkeley Laboratory
//
// Author List:
//      Dave Brown 3/15/97
//------------------------------------------------------------------------
#ifndef KALCONTEXT_HH
#define KALCONTEXT_HH

#include <math.h>
#include "BTrk/TrkBase/TrkEnums.hh"
#include "BTrk/TrkBase/TrkParticle.hh"
#include "BTrk/TrkBase/TrkDirection.hh"
#include "BTrk/BField/BFieldIntegrator.hh"
#include "BTrk/DetectorModel/DetSet.hh"

#include <math.h>
class TrkVolume;

class KalContext {
public:
  KalContext(); // only a default constructor, subclasses should provide configuration method
  virtual ~KalContext();
// Accessors
  double minGap() const { return _mingap; }
  double maxParamDiff(trkDirection trkdir) const {
    return _maxpardif[trkdir]; }
  double distanceTolerance() const {return _disttol; }
  unsigned maxIterations() const { return _maxiter; }
  unsigned maxIntersections() const { return _maxinter; }
  double intersectionTolerance() const {return _intertol; }
  bool materialSites() const {return _matcorr; }
  bool bendSites() const { return _fieldcorr; }
  double smearFactor() const { return _smearfactor; }
  double maxSiteDMom() const { return _sitethresh; }
  double maxDMom() const { return _momthresh; }
// minimum DOFs required for each view.  bothview means 'overall DOFs'
  unsigned minDOF(TrkEnums::TrkViewInfo view=TrkEnums::bothView) const { 
    return _mindof[view]; }
  double maxMomDiff() const { return _maxmomdiff; }
  double momUpdateFactor() const { return _momfac; }
  const BFieldIntConfig& bFieldIntConfig() const { return _bintconfig; }
// conditions; these must be supplied by subclassing
  virtual const TrkVolume* trkVolume(trkDirection trkdir) const = 0;
  virtual BField const& bField() const = 0;
  BFieldIntegrator const& bFieldIntegrator() const;
  double minFltLen() const { return _minfltlen; }
  double minMom() const { return _minmom; }
  double fltEpsilon() const { return _fltepsilon; }
  double divergeFlt() const { return _divergeflt; }
  double minDot() const { return _mindot; }
  const DetSet *getDetModel () const { return _trkmodel; }
protected:
// configuration data
  double _disttol; // tolerance on the maximum distance for iteration
  double _intertol; // tolerance on the maximum distance for re-intersection
  double _maxpardif[2]; // tolerance on parameter difference (each end)
  ulong _maxiter; // maximum number of iterations allowed
  ulong _maxinter; // maximum number of intersections allowed
  bool _matcorr; // use material sites
  bool _fieldcorr; // use bend sites
  double _smearfactor; // initial covariance smearing factor
  double _sitethresh; // single site maximum momentum fraction before stopping a track
  double _momthresh; // minimum momentum fraction before stopping a track
  double _mingap; // minimum gap between adjacent sites to build a new traj piece
  double _minfltlen; // minimum flight length allowed for a KalRep
  double _minmom; // minimum momentum
  double _fltepsilon; // small flight length buffer
  double _divergeflt; // flight length change to signify a diverging fit
  double _mindot;  // minimum direction dot product change for a traj to be 'reasonable'
  ulong _mindof[3]; // minimum number of DOFs to allow fit to succeed (can be 0)
  double _maxmomdiff; // maximum momentum difference before forcing iteration
  double _momfac; // factor for updating momentum on iteration; 0=full update,
// infinity = don't update at all.  Scale is set by track momentum
  BFieldIntConfig _bintconfig;
  mutable BFieldIntegrator* _bint; // field integrator, created on first access

  const DetSet* _trkmodel;//detector material description

//disallow
  KalContext(const KalContext& other);
  KalContext& operator = (const KalContext& other);
};

#endif

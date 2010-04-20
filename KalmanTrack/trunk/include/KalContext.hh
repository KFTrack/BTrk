//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: KalContext.hh 103 2010-01-15 12:12:27Z stroili $
//
// Description:
//      class KalContext
//      This class describes the basic parameters used in the Kalman fit.
//      It is implemented using odmg types so that it can be stored in the
//      config database (eventually)
//
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
#include "BaBar/BaBarODMGTypes.h"
#include "TrkEnv/TrkVolumeHandle.hh"
#include "TrkBase/TrkDirection.hh"
#include "TrkBase/TrkEnums.hh"
#include "BaBar/PdtPid.hh"
#include <math.h>
class TrkVolume;

class KalContext {
public:
  KalContext(); // only a default constructor, sets default values
  ~KalContext();
// Accessors
  double minGap() const { return _mingap; }
  double tBuffer() const { return _trajbuff; }
  double maxParamDiff(trkDirection trkdir) const {
    return _maxpardif[trkdir]; }
  double distanceTolerance() const {return _disttol; }
  unsigned maxIterations() const { return _maxiter; }
  unsigned maxIntersections() const { return _maxinter; }
  double intersectionTolerance() const {return _intertol; }
  d_Boolean materialSites() const {return _matsites; }
  d_Boolean bendSites() const { return _bends; }
  double smearFactor() const { return _smearfactor; }
  double maxSiteDMom() const { return _sitethresh; }
  double maxDMom() const { return _momthresh; }
  double localSiteDMom() const { return _sitepfrac; }
  double localSiteDeflect() const { return _sitedflct; }
// This function uses the _environment_ to provide the requested
// Tracking Volume
  const TrkVolume* trkVolume(trkDirection trkdir) const;
  double bFieldIntMinStep() const { return _bintminstep; }
  double bFieldIntMaxStep() const { return _bintmaxstep; }
  double bFieldIntMaxFraction() const { return _bintmaxfrac; }
  double bFieldIntTolerance() const { return _binttolerance; }
  double bFieldDivMinStep() const { return _bdivminstep; }
  double bFieldDivMaxStep() const { return _bdivmaxstep; }
  double bFieldDivMaxFraction() const { return _bdivmaxfrac; }
  double bFieldDivTolerance() const { return _bdivtolerance; }
  PdtPid::PidType defaultType() const { return (PdtPid::PidType)_defpid; }
// minimum DOFs required for each view.  bothview means 'overall DOFs'
  unsigned minDOF(TrkEnums::TrkViewInfo view=TrkEnums::bothView) const { 
    return _mindof[view]; }
  double maxMomDiff() const { return _maxmomdiff; }
  d_Boolean stopHots() const { return _stophots; }
  d_Boolean forbidAmbigFlips() const { return _ambigflip; }
  double momUpdateFactor() const { return _momfac; }
// Set functions
  void setMinGap(double mingap) { _mingap = mingap; }
  void setTrajBuffer(double tbuf) { _trajbuff = tbuf; }
  void setMaxParamDiff(trkDirection trkdir,double mpdif) {
    _maxpardif[trkdir] = mpdif; }
  void setDistanceTolerance(double disttol)  {_disttol = disttol; }
  void setIntersectionTolerance(double intertol)  {_intertol = intertol; }
  void setMaxIterations(unsigned maxiter)  { _maxiter = maxiter; }
  void setMaxIntersections(unsigned maxinter)  { _maxinter = maxinter; }
  void setMaterialSites(d_Boolean matsites)  {_matsites = matsites; }
  void setBendSites(d_Boolean bends)  { _bends = bends; }
  void setStopHots(d_Boolean stop)  { _stophots = stop; }
  void setForbidAmbigFlips(d_Boolean flip)  { _ambigflip = flip; }
  void setMomUpdateFactor(double factor) { _momfac = fabs(factor); }
  void setSmearFactor(double smearfactor)  { _smearfactor = smearfactor; }
  void setMaxSiteDMom(double sitethresh)  { _sitethresh = sitethresh; }
  void setMaxDMom(double momthresh)  { _momthresh = momthresh; }
  void setLocalSiteDMom(double sitepfrac)  { _sitepfrac = sitepfrac; }
  void setLocalSiteDeflect(double sitedflct)  { _sitedflct = sitedflct; }
  void setVolume(TrkVolumeHandle::trkvolumes tvol,trkDirection trkdir){
    _volumes[trkdir] = tvol; }
  void setBFieldIntegration(double minstep,double maxstep,
			    double maxfrac,double tolerance,
			    double divminstep,double divmaxstep,
			    double divmaxfrac,double divtolerance) {
    _bintminstep = minstep; // limits for the BField integrator
    _bintmaxstep = maxstep;
    _bintmaxfrac = maxfrac;
    _binttolerance = tolerance;
    _bdivminstep = divminstep; // limits for the track divider
    _bdivmaxstep = divmaxstep;
    _bdivmaxfrac = divmaxfrac;
    _bdivtolerance = divtolerance;
  }
  void setDefaultType(PdtPid::PidType newdefaulttype) {
    _defpid = newdefaulttype; }
  void setMinDOF(unsigned mindof,TrkEnums::TrkViewInfo view) {
    _mindof[view] = mindof; }
  void setMaxMomDiff(double momdiff) {
    _maxmomdiff = momdiff; }
private:
  d_Double _disttol; // tolerance on the maximum distance for iteration
  d_Double _intertol; // tolerance on the maximum distance for re-intersection
  d_Double _maxpardif[2]; // tolerance on parameter difference (each end)
  d_ULong _maxiter; // maximum number of iterations allowed
  d_ULong _maxinter; // maximum number of intersections allowed
  d_Boolean _matsites; // use material sites
  d_Boolean _bends; // use bend sites
  d_Double _smearfactor; // initial covariance smearing factor
  d_Double _sitethresh; // single site maximum momentum fraction before stopping a track
  d_Double _momthresh; // minimum momentum fraction before stopping a track
  d_Double _sitepfrac; // threshold on p-frac change to use local reference
  d_Double _sitedflct; // threshold on deflection to use local reference
  d_Double _mingap; // minimum gap between adjacent sites to build a new traj piece
  d_Double _trajbuff; // trajectory piece merging buffer size
  d_Long _volumes[2];  // Inner and out tracking volumes
  d_Double _bintminstep; // BField integration parameters
  d_Double _bintmaxstep;
  d_Double _bintmaxfrac;
  d_Double _binttolerance;
  d_Double _bdivminstep; // BField track divider parameters
  d_Double _bdivmaxstep;
  d_Double _bdivmaxfrac;
  d_Double _bdivtolerance;
  d_ULong _defpid; // default PID to use in Kalman fit
  d_ULong _mindof[3]; // minimum number of DOFs to allow fit to succeed (can be 0)
  d_Double _maxmomdiff; // maximum momentum difference before forcing iteration
  d_Boolean _stophots; // deactivate hots beyond the dE/dx stopping point
  d_Boolean _ambigflip; // allow ambiguity flips when updating HOTs
  d_Double _momfac; // factor for updating momentum on iteration; 0=full update,
// infinity = don't update at all.  Scale is set by track momentum
//disallow
  KalContext(const KalContext& other);
  KalContext& operator = (const KalContext& other);

};

#endif

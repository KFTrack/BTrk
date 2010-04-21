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
  bool materialSites() const {return _matsites; }
  bool bendSites() const { return _bends; }
  double smearFactor() const { return _smearfactor; }
  double maxSiteDMom() const { return _sitethresh; }
  double maxDMom() const { return _momthresh; }
  double localSiteDMom() const { return _sitepfrac; }
  double localSiteDeflect() const { return _sitedflct; }
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
  bool stopHots() const { return _stophots; }
  bool forbidAmbigFlips() const { return _ambigflip; }
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
  void setMaterialSites(bool matsites)  {_matsites = matsites; }
  void setBendSites(bool bends)  { _bends = bends; }
  void setStopHots(bool stop)  { _stophots = stop; }
  void setForbidAmbigFlips(bool flip)  { _ambigflip = flip; }
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
  double _disttol; // tolerance on the maximum distance for iteration
  double _intertol; // tolerance on the maximum distance for re-intersection
  double _maxpardif[2]; // tolerance on parameter difference (each end)
  ulong _maxiter; // maximum number of iterations allowed
  ulong _maxinter; // maximum number of intersections allowed
  bool _matsites; // use material sites
  bool _bends; // use bend sites
  double _smearfactor; // initial covariance smearing factor
  double _sitethresh; // single site maximum momentum fraction before stopping a track
  double _momthresh; // minimum momentum fraction before stopping a track
  double _sitepfrac; // threshold on p-frac change to use local reference
  double _sitedflct; // threshold on deflection to use local reference
  double _mingap; // minimum gap between adjacent sites to build a new traj piece
  double _trajbuff; // trajectory piece merging buffer size
  d_Long _volumes[2];  // Inner and out tracking volumes
  double _bintminstep; // BField integration parameters
  double _bintmaxstep;
  double _bintmaxfrac;
  double _binttolerance;
  double _bdivminstep; // BField track divider parameters
  double _bdivmaxstep;
  double _bdivmaxfrac;
  double _bdivtolerance;
  ulong _defpid; // default PID to use in Kalman fit
  ulong _mindof[3]; // minimum number of DOFs to allow fit to succeed (can be 0)
  double _maxmomdiff; // maximum momentum difference before forcing iteration
  bool _stophots; // deactivate hots beyond the dE/dx stopping point
  bool _ambigflip; // allow ambiguity flips when updating HOTs
  double _momfac; // factor for updating momentum on iteration; 0=full update,
// infinity = don't update at all.  Scale is set by track momentum
//disallow
  KalContext(const KalContext& other);
  KalContext& operator = (const KalContext& other);

};

#endif

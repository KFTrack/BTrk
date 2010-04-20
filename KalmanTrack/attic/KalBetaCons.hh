// ----------------------------------------------------------------------------
// File and Version Information:
//      $Id: KalBetaCons.hh,v 1.9 2007/09/06 23:12:38 brownd Exp $
//
//  Description:
//  Class to describe a constraint on the particle beta, coming from PID
//  into the fit
//
// Copyright Information:
//	Copyright (C) 2005	Lawrence Berkeley Laboratory
//
//  Authors: Marco Battaglia, David Brown
//-----------------------------------------------------------------------------
#ifndef KALBETACONS_HH
#define KALBETACONS_HH
#include "KalmanTrack/KalSite.hh"
#include "TrkBase/TrkEnums.hh"
#include "KalmanTrack/KalWeight.hh"
#include "BbrGeom/BbrDoubleErr.hh"
#include "SvtPid/SvtPidCalib.hh"
//
class TrkDifPieceTraj;
class TrkRecoTrk;
class TrkParams;
class TrkRep;
class BField;
class SvtPidInfo;
class SvtHitOnTrack;
class PdtEntry;
//
//  Define the class
//
class KalBetaCons : public KalSite {
public:
//  Constructors
  KalBetaCons(const TrkDifPieceTraj* reftraj, const BField*, double fltlen,
	      const BbrDoubleErr& momcons,PdtPid::PidType pid);
// build from SvtPidInfo
  KalBetaCons(const TrkDifPieceTraj* reftraj, const BField*,
	      const SvtPidInfo* svtpid,PdtPid::PidType pid);
// build from SvtHOTs
  KalBetaCons(const TrkDifPieceTraj* reftraj, const BField*,
              const std::vector<const SvtHitOnTrack*>& svthots,PdtPid::PidType pid);
// copy constructor
  KalBetaCons(const KalBetaCons&);
// clone function
  KalSite* clone(const KalRep*) const;
//
  virtual ~KalBetaCons();
//
//  Fit functions
//
  bool process(const KalSite*,trkDirection idir);
  bool update(const TrkDifPieceTraj*,double);
//
//  Access
//
  void printAll(std::ostream& os = std::cout) const;
// chisq including boolean value
  bool chisquared(double& chisq,const KalSite*, trkDirection) const;
// same as above, except just computing chisquared.
  bool chisquared(double& chisq, bool ignoreactive=true) const;
  virtual unsigned nDof(TrkEnums::TrkViewInfo view=TrkEnums::bothView) const;
  virtual bool isActive() const { return _active; }
  bool setActivity(bool active);
  PdtPid::PidType particleType() const { return _partid; }
  const BbrDoubleErr dedxMom() const { return _dedxmom; }
  const SvtPidInfo* svtPid() const { return _svtpid; }
  int usedHits() const { return _nused; }
  SvtPidCalib::errcode errorCode() const { return _errcode;}
private:
  PdtPid::PidType _partid; // define the particle type (and hence mass)
  const PdtEntry* _pdt;
  BbrDoubleErr _dedxmom;  // dE/dx momentum and covariance
  double _refmom;  // reference momentum, from fit
  KalWeight _weight; // hit parameters/weight
  bool _active; // activity of this site
  const BField* _bfield; // cache of bfield, needed to compute momentum
  const SvtPidInfo* _svtpid; // reference to SvtPidInfo
  int _nused; // # of svt hits
  SvtPidCalib::errcode _errcode; // status of dEdx momentum function
// private functions
  bool updateCache(const TrkDifPieceTraj*);
};
#endif


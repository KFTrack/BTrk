//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkHitOnTrk.cc,v 1.60 2006/03/04 19:52:08 brownd Exp $
//
// Description:
//
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author(s): Steve Schaffner
//
//
// Modified 3-july-97 by Leon Rochester to add '<' & '==' operators for
// Roguewave Sorted Vector
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"
#include "TrkBase/TrkHitOnTrk.hh"
#include <assert.h>
#include "TrkBase/TrkRecoTrk.hh"
#include "TrkBase/TrkFundHit.hh"
#include "TrkBase/TrkRep.hh"
#include "TrkBase/TrkDifTraj.hh"
#include "TrkBase/TrkPoca.hh"
#include "TrkBase/TrkDifPoca.hh"
#include "TrkBase/TrkSimpTraj.hh"
#include "ErrLogger/ErrLog.hh"
using std::endl;
using std::ostream;

const DchHitOnTrack* TrkHitOnTrk::dchHitOnTrack() const {return 0;}
const SvtHitOnTrack* TrkHitOnTrk::svtHitOnTrack() const {return 0;}

TrkHitOnTrk::TrkHitOnTrk(const TrkFundHit* hit,double tolerance) :
  _parentRep(0),
  _theHit(const_cast<TrkFundHit*>(hit)),
  _isActive(true),
  _isUsable(true),
  //make caches invalid
  _hitRms(-9.e50),
  _trkLen(0.0),
  _hitLen(0.0),
  _trkTraj(0),
  _poca(0),
  _tolerance(tolerance)
{  }

// This is effectively a copy ctor:
TrkHitOnTrk::TrkHitOnTrk(const TrkHitOnTrk &oldHit, TrkRep *newRep,
                         const TrkDifTraj *trkTraj)
  : _parentRep(newRep),
    _theHit(oldHit._theHit),
    _isActive(oldHit._isActive),
    _isUsable(oldHit._isUsable),
    _hitRms(oldHit._hitRms),
    _trkLen(oldHit._trkLen),
    _hitLen(oldHit._hitLen),
    _resid(9999.9),
    _trkTraj(0),
    _poca(0),
    _tolerance(oldHit._tolerance)
{
  assert (0 != newRep);
  if( oldHit._trkTraj!=0 && trkTraj!=0 && oldHit._trkTraj == trkTraj ) {
                              // re-use cache as traj are the same
        _resid=oldHit._resid;
        _trkTraj=trkTraj;
        _poca = (oldHit._poca==0 ? 0 :new TrkPoca(*oldHit._poca));
  } else {
       double fl = oldHit.fltLen();
       double dum;
       const TrkSimpTraj *t1= (oldHit._trkTraj==0?0:oldHit._trkTraj->localTrajectory(fl,dum));
       const TrkSimpTraj *t2= (trkTraj==0?0:trkTraj->localTrajectory(fl,dum));
       if( t1 != 0 && t2 != 0 && t1->parameters()->parameter() == t2->parameters()->parameter() ) {
                              // re-use cache as traj are sufficiently equiv
           _resid=oldHit._resid;
           _trkTraj=trkTraj;
           _poca = (oldHit._poca==0 ? 0: new TrkPoca(*oldHit._poca));

       }
  }
  // Only record hots if on default TrkRep
  const TrkRep* rep = getParentRep();
  if (rep->particleType() == rep->parentTrack()->defaultType()) setUsedHit();
}

TrkHitOnTrk::~TrkHitOnTrk()
{
  delete _poca;

  const TrkRecoTrk *parentTrack = _parentRep==0?0: _parentRep->parentTrack();

  if (_parentRep!=0 && parentTrack!=0 && _parentRep->particleType() == parentTrack->defaultType() ) {
    setUnusedHit();
  }
}

void
TrkHitOnTrk::setActivity(bool turnOn)
{
  if (!isUsable() || isActive()==turnOn) return;
  if (getParentRep() != 0) {    // needed until Rep-less HoTs go away
    turnOn ? parentRep()->activateHot(this)
           : parentRep()->deactivateHot(this);
  } else {
    _isActive = turnOn;
  }
}

void
TrkHitOnTrk::setUsability(int usability)
{
  _isUsable = usability;
  if (isActive() && !isUsable()) {
    _isActive = false;
    if(getParentRep() != 0)parentRep()->deactivateHot(this);
  }
  if (!isActive() && mustUse()) {
    _isActive = true;
    if(getParentRep() != 0)parentRep()->activateHot(this);
  }
}

double
TrkHitOnTrk::weight() const
{
  // could be cached
  double rms=hitRms();
  assert(rms > 0);
  return double(1) / ( rms * rms );
}

void
TrkHitOnTrk::print(ostream& o) const
{
  o << (dchHitOnTrack()?"Dch":svtHitOnTrack()?"Svt":"Unknown")
      << " HOT: hitlength: " << hitLen() << "  flightlength: "
      << fltLen() << "  active: " << (isActive() != 0) << endl;
}

void
TrkHitOnTrk::printAll(ostream& o) const
{
  print(o);
}

TrkRecoTrk*
TrkHitOnTrk::parentTrack() const
{
  return parentRep()->parentTrack();
}

const TrkRecoTrk*
TrkHitOnTrk::getParentTrack() const
{
  return getParentRep()->parentTrack();
}

PdtPid::PidType
TrkHitOnTrk::particleType() const
{
  return getParentRep()->particleType();
}

void
TrkHitOnTrk::setUsedHit()
{
  if (hit() != 0) hit()->setUsedHit(this);
}

void
TrkHitOnTrk::setUnusedHit()
{
  if (hit() != 0) hit()->setUnusedHit(this);
}

bool TrkHitOnTrk::operator==(const TrkHitOnTrk &rhs) const
{
  return this == &rhs;
}

int
TrkHitOnTrk::ambig() const
{
  return 0;// by default no ambiguity
}

void
TrkHitOnTrk::setAmbig(int newambig)
{} // by default nothing to set

double
TrkHitOnTrk::resid(bool exclude) const
{
    double r(-99999.9),re(-9999.9);
    bool s=getParentRep()->resid(this,r,re,exclude);
    if (!s && r<-99999.8) {
      ErrMsg(routine) << "error calling parentRep()->residual()" << endmsg;
    }
    return r;
}

bool
TrkHitOnTrk::resid(double &resid, double &residErr, bool exclude) const
{
    assert(getParentRep()!=0);
    return getParentRep()->resid(this,resid,residErr,exclude);
}

double
TrkHitOnTrk::residual() const
{
    assert (_trkTraj == &(getParentRep()->traj()));
    return _resid;
}

TrkErrCode
TrkHitOnTrk::updatePoca(const TrkDifTraj* trkTraj, bool maintainAmb)
{
    _trkTraj = (trkTraj!=0?trkTraj:&(getParentRep()->traj()));
    if (_poca==0) {
        _poca = new TrkPoca(*_trkTraj,fltLen(),
                            *hitTraj(), hitLen(),_tolerance);
    } else {
        *_poca = TrkPoca(*_trkTraj,fltLen(),
                         *hitTraj(), hitLen(),_tolerance);
    }
    if(_poca->status().failure()) {
        if(isActive())ErrMsg(routine) << " TrkPoca failed in TrkHitOnTrk::updatePoca"
                                      << endmsg;
        delete _poca; _poca=0;
        return TrkErrCode(TrkErrCode::fail,4);
    }
    _trkLen = _poca->flt1();
    _hitLen = _poca->flt2();
    double dca=_poca->doca();
    if (!maintainAmb) setAmbig(dca>0?+1:-1);
    return TrkErrCode(TrkErrCode::succeed);
}

TrkErrCode
TrkHitOnTrk::getFitStuff(HepVector &derivs, double &deltaChi ) const
{
    if (_poca==0 || _poca->status().failure()) {
        return TrkErrCode(TrkErrCode::fail);
    }
    // FIXME: I wish I could tell poca to NOT iterate
    //        and ONLY compute the distance & derivatives...
    TrkDifPoca poca(*_trkTraj,fltLen(),*hitTraj(), hitLen(),_tolerance);
    if (poca.status().failure()) {
        return TrkErrCode(TrkErrCode::fail);
    }
    if (derivs.num_row() != 0) {
        poca.fetchDerivs(derivs);
    } else {
        derivs = poca.derivs();
    }
    double sigInv = 1. / hitRms();
    deltaChi = _resid * sigInv; // NOTE: use _INTERNAL_ residual
    derivs *= sigInv;
    return TrkErrCode(TrkErrCode::succeed);
}

TrkErrCode
TrkHitOnTrk::getFitStuff(double &deltaChi)  const
{
    assert (_trkTraj == &(getParentRep()->traj()));
    deltaChi=_resid/hitRms(); // NOTE: use _INTERNAL_ residual
    return TrkErrCode(TrkErrCode::succeed);
}


ostream&
operator<<(ostream& o, const TrkHitOnTrk& x)
{
        x.print(o); return o;
}

TrkErrCode
TrkHitOnTrk::pocaStatus() const {
  if(_poca != 0)
    return _poca->status();
  else
    return TrkErrCode(TrkErrCode::fail);
}

//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkHit.cc,v 1.60 2006/03/04 19:52:08 brownd Exp $
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
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/TrkBase/TrkHit.hh"
#include <assert.h>
#include "BTrk/TrkBase/TrkRep.hh"
#include "BTrk/TrkBase/TrkDifTraj.hh"
#include "BTrk/TrkBase/TrkPoca.hh"
#include "BTrk/TrkBase/TrkDifPoca.hh"
#include "BTrk/TrkBase/TrkSimpTraj.hh"
#include "BTrk/BaBar/ErrLog.hh"
using std::endl;
using std::ostream;
using namespace CLHEP;

double TrkHit::_tolerance(1e-5);
void TrkHit::setTolerance(double tol) { _tolerance = tol; }

TrkHit::TrkHit() :
  _parentRep(0),
  _trkTraj(0),
  _isActive(false),
  //make caches invalid
  _hitRms(-1.0),
  _trkLen(0.0),
  _hitLen(0.0),
  _resid(-1.0)
{  }

TrkHit::~TrkHit()
{
}

void
TrkHit::setActivity(bool turnOn)
{
  if (isActive()==turnOn ) return;
  if (getParentRep() != 0) {    // needed until Rep-less THits go away
    turnOn ? parentRep()->activateHit(this)
           : parentRep()->deactivateHit(this);
  } else {
    _isActive = turnOn;
  }
}

double
TrkHit::weight() const
{
  // could be cached
  double rms=hitRms();
  assert(rms > 0);
  return double(1) / ( rms * rms );
}

void
TrkHit::print(ostream& o) const
{
}

void
TrkHit::printAll(ostream& o) const
{
  print(o);
}

TrkParticle const&
TrkHit::particleType() const
{
  return getParentRep()->particleType();
}

bool TrkHit::operator==(const TrkHit &rhs) const
{
  return this == &rhs;
}

int
TrkHit::ambig() const
{
  return 0;// by default no ambiguity
}

void
TrkHit::setAmbig(int newambig)
{} // by default nothing to set

double
TrkHit::resid(bool exclude) const
{
    double r(-99999.9),re(-9999.9);
    bool s=getParentRep()->resid(this,r,re,exclude);
    if (!s && r<-99999.8) {
      ErrMsg(routine) << "error calling parentRep()->residual()" << endmsg;
    }
    return r;
}

bool
TrkHit::resid(double &resid, double &residErr, bool exclude) const
{
    assert(getParentRep()!=0);
    return getParentRep()->resid(this,resid,residErr,exclude);
}

double
TrkHit::residual() const
{
    // cachching resid separate from poca risks inconsistency, FIXME!!!
    return _resid;
}

TrkErrCode
TrkHit::updatePoca(const TrkDifTraj* trkTraj)
{
    if (trkTraj==0)trkTraj = &getParentRep()->traj();
    _poca = TrkPoca(*trkTraj,fltLen(),*hitTraj(), hitLen(),_tolerance);
    if(_poca.status().success()){
    // this copying of flightlens between poca and TrkHit is error-prone
    //  and risks cache inconsistency, FIXME!!!
      _trkLen = _poca.flt1();
      _hitLen = _poca.flt2();
    }
    return _poca.status();
}

TrkErrCode
TrkHit::getFitStuff(HepVector &derivs, double &deltaChi ) const
{
    if (_poca.status().failure()) {
        return TrkErrCode(TrkErrCode::fail);
    }
    // This copy risks having inconsistent values in poca, difpoca and TrkHit
    // It is also inefficient.
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
TrkHit::getFitStuff(double &deltaChi)  const
{
    assert (_trkTraj == &(getParentRep()->traj()));
    deltaChi=_resid/hitRms(); // NOTE: use _INTERNAL_ residual
    return TrkErrCode(TrkErrCode::succeed);
}


ostream&
operator<<(ostream& o, const TrkHit& x)
{
        x.print(o); return o;
}



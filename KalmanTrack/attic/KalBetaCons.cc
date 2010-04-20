// ----------------------------------------------------------------------------
// File and Version Information:
//      $Id: KalBetaCons.cc,v 1.16 2007/09/06 23:12:30 brownd Exp $
//
//  Description:
//  Class to describe a constraint on the momentum
//  into the fit
//
// Copyright Information:
//	Copyright (C) 2005	Lawrence Berkeley Laboratory
//
//  Authors: David Brown
//-----------------------------------------------------------------------------
#include "BaBar/BaBar.hh"
#include "KalmanTrack/KalBetaCons.hh"
#include "KalmanTrack/KalRep.hh"
#include "TrkBase/TrkDifPieceTraj.hh"
#include "TrkBase/TrkParams.hh"
#include "TrkBase/TrkErrCode.hh"
#include "TrkBase/TrkMomCalculator.hh"
#include "TrkBase/HelixTraj.hh"
#include "CLHEP/Vector/ThreeVector.h"
#include "difAlgebra/DifNumber.hh"
#include "PDT/Pdt.hh"
#include "SvtPid/SvtPidInfo.hh"
#include "BbrGeom/BbrDoubleErr.hh"
#include "PDT/Pdt.hh"
#include <math.h>
#include <assert.h>
#include <algorithm>
using std::endl;
using std::ostream;

KalBetaCons::KalBetaCons(const TrkDifPieceTraj* reftraj,const BField* bfield,double fltlen,const BbrDoubleErr& momcons,PdtPid::PidType pid) : 
  KalSite(betaConsSite),_partid(pid),_pdt(0),_dedxmom(momcons),
  _refmom(0.0),_active(true),_bfield(bfield),_svtpid(0),_nused(5),_errcode(SvtPidCalib::Fail)
{
  _pdt = Pdt::lookup(_partid,-1);
  setTraj(reftraj,fltlen);
  updateCache(reftraj);
}

KalBetaCons::KalBetaCons(const TrkDifPieceTraj* reftraj,const BField* bfield,
                         const SvtPidInfo* svtpid,PdtPid::PidType pid) : 
  KalSite(betaConsSite),_partid(pid),_pdt(0),
  _dedxmom(0.0),_refmom(0.0),_active(false),_bfield(bfield),_svtpid(svtpid),_nused(0),
  _errcode(SvtPidCalib::Fail)
{
  _pdt = Pdt::lookup(_partid,-1);
  if(_svtpid != 0 && _pdt != 0){
    _dedxmom = _svtpid->getMomentum(_errcode,_pdt);
// temporary (till we call truncation algorithm explicitly)
    _nused = _svtpid->getNPoints();
    if(_errcode <= SvtPidCalib::OK) {
      if(_dedxmom.value() > 0.0) {
        _active = true;
	setTraj(reftraj,_svtpid->getCalibrationfltlen());
        updateCache(reftraj);
      } else
        ErrMsg(error) << "Negative momentum from SvtPidInfo" << endmsg;
    }
  } else
    ErrMsg(error) <<"Undefined constraint!" << endmsg;
}

KalBetaCons::KalBetaCons(const TrkDifPieceTraj* reftraj, const BField* bfield,
                         const std::vector<const SvtHitOnTrack*>& svthots,PdtPid::PidType pid) :
  KalSite(betaConsSite),_partid(pid),_pdt(0),
  _dedxmom(0.0),_refmom(0.0),_active(false),_bfield(bfield),_svtpid(0),_nused(0),_errcode(SvtPidCalib::Fail)
{
  _pdt = Pdt::lookup(_partid,-1);
  if(_pdt != 0){
// need to call getMomentum on each hot and truncate: FIXME!!!!!
    _nused = svthots.size();
  } else
    ErrMsg(error) <<"Undefined PdtPid!" << endmsg;
}


KalBetaCons::KalBetaCons(const KalBetaCons& other) :
  KalSite(other),
  _partid(other._partid),_pdt(other._pdt),_dedxmom(other._dedxmom),
  _refmom(other._refmom),
  _weight(other._weight),
  _active(other._active),
  _bfield(other._bfield),
  _svtpid(other._svtpid),
  _errcode(other._errcode)
{
}

KalSite*
KalBetaCons::clone(const KalRep* rep) const {
  if(_svtpid != 0)
    return new KalBetaCons(rep->referenceTraj(),_bfield,
			   _svtpid, rep->particleType());
  else
    return new KalBetaCons(rep->referenceTraj(),_bfield,
			   globalLength(),_dedxmom,
			   rep->particleType());
}

KalBetaCons::~KalBetaCons() {}

bool
KalBetaCons::process(const KalSite* prevsite,trkDirection idir) {
  bool status(false);
  if(isActive())
    status = processWeight(prevsite,idir,_weight);
  else {
// copy the previous site
    copySite(prevsite,idir);
    if(validSite(idir)){
// set the site processed in this direction
      setFit(idir);
      status = true;
    }
  }
  return status;
}


bool
KalBetaCons::update(const TrkDifPieceTraj* reftraj,double){
  return updateCache(reftraj);
}

void
KalBetaCons::printAll(std::ostream& os) const {
  os << "KalBetaCons ";
  KalSite::printAll(os);
  os << "momentum constraint value = " << _dedxmom <<  endl;
}

bool
KalBetaCons::chisquared(double& chisq,const KalSite* prevsite,trkDirection tdir) const {
  return chisquared(chisq,false);
}

bool
KalBetaCons::chisquared(double& chisq, bool ignoreactive) const {
  bool retval(false);
  if(_active || ignoreactive) {
    retval = true;
    chisq = sqr(_dedxmom.value() - _refmom)/_dedxmom.covariance();
  }
  return retval;
}

unsigned
KalBetaCons::nDof(TrkEnums::TrkViewInfo) const {
  return _active ? 1 : 0;
}

bool
KalBetaCons::updateCache(const TrkDifPieceTraj* reftraj) {
  bool retval(true);
// get the local trajectory
  double loclen;
  const TrkSimpTraj* loctraj = reftraj->localTrajectory(globalLength(),loclen);
// make sure this is a helix traj
  const HelixTraj* htraj = static_cast<const HelixTraj*>(loctraj);
  if(htraj != 0){
//  update reference momentum
    DifNumber mommag = TrkMomCalculator::momMag(*htraj,*_bfield);
// take the new momentum value from this
    _refmom = mommag.number();
// also the derivatives
    HepVector linrel = mommag.derivatives();
// expand the measurement error into parameter space, invert to get weight
    HepSymMatrix hwt = vT_times_v(linrel)/_dedxmom.covariance();
// compute the residual: always measurement - prediction;
    double resid = _dedxmom.value() - _refmom;
// get the reference parameters
    const HepVector& refvec = loctraj->parameters()->parameter();
// now compute the information vector
    HepVector deltapar = linrel*resid;
    _weight = KalWeight(hwt*(refvec + deltapar),hwt);
  } else 
    retval = false;
  return retval;
}

bool
KalBetaCons::setActivity(bool active) {
// we must reset the site and the fit to be not current when the state changes
  if(active != _active){
    setFit(trkIn, false);
    setFit(trkOut, false);
    _active = active;
    return true;
  } else
    return false;
}

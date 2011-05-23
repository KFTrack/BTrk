// ------------------------------------------------------------------------------
// File and Version Information:
//      $Id: KalHit.cc,v 1.68 2007/10/17 21:41:11 brownd Exp $
//
//  Description:
//  Class to describe a kalman filter hit site.  This puts measurement information
//  into the fit
//
// Copyright Information:
//	Copyright (C) 1996	Lawrence Berkeley Laboratory
//
//  Authors: Dave Brown, 3/7/97
//------------------------------------------------------------------------------
#include "BaBar/BaBar.hh"
#include "KalmanTrack/KalHit.hh"
#include "KalmanTrack/KalRep.hh"
#include "TrkBase/TrkDifPieceTraj.hh"
#include "TrkBase/TrkParams.hh"
#include "TrkBase/TrkErrCode.hh"
#include "CLHEP/Vector/ThreeVector.h"
#include <math.h>
#include <assert.h>
using std::endl;
using std::ostream;

//
//  Constructor
//
KalHit::KalHit(const TrkDifPieceTraj* straj,TrkHitOnTrk* hot,
	       bool ambigflip) :
  KalSite(hitSite),_hot(hot),_residual(0.0),_hotstate(_hot->isActive()),
  _ambigflip(ambigflip),_fltdif(0.0)
{ updateCache(straj); }
//
//  copy constructor
//
KalHit::KalHit(const KalHit& other) :
  KalSite(other),_hot(other._hot),
  _hitweight(other._hitweight),
  _residual(other._residual),
  _linrel(other._linrel),
  _hotstate(other._hotstate),
  _ambigflip(other._ambigflip),
  _fltdif(other._fltdif)
{;}
// clone function
KalSite*
KalHit::clone(const KalRep* krep) const {
  KalHit* newhit = new KalHit(*this);
  assert(newhit != 0);
// clone a new HOT to point to the new rep
  newhit->cloneHot(krep);
// update the hit site
  if(newhit->hitOnTrack()->isActive())
    newhit->updateCache(krep->referenceTraj());
  return newhit;
}
//
KalHit::~KalHit(){;}
//
//  Update for a new reference trajectory
//
bool
KalHit::update(const TrkDifPieceTraj* newtraj,double mom) {
  return updateCache(newtraj);
}
//
//  Process the hit
//
bool
KalHit::process(const KalSite* prevsite,trkDirection idir){
// record the state of the hot when processed
  _hotstate = _hot->isActive();
  return processWeight(prevsite,idir,_hitweight);
}
//
//  access
//
void
KalHit::printAll(ostream& os) const {
  os << "Hit ";
  KalSite::printAll(os);
  os << "Associated HOT = ";
  _hot->print(os);
//  os << "Weight parameters and info matrix " << endl;
//  _hitweight.printAll(os);
  os << endl << "Residual = " << _residual << endl;
  os << "Derivative vector = " << _linrel << endl;
}
//
bool
KalHit::updateCache(const TrkDifPieceTraj* reftraj){
// test against the old flightlength, and only update
// if the parameters have changed
//
//  DNB12/06/05.  We Must call updateMeasurement whenever reftraj changes, otherwise
// the hot can become self-inconsistent
//
  double oldlen=_hot->fltLen();
  bool retval(true);
  if(setTraj(reftraj,oldlen) || reftraj != _hot->trkTraj()){
    TrkErrCode residerr = updateMeasurement(*_hot,reftraj,_ambigflip);
    if(residerr.success())
      residerr = _hot->getFitStuff(_linrel,_residual);
    retval = residerr.success();
    if (retval) {
      setTraj(reftraj,_hot->fltLen());
//  compute information quantities
// the sign here defines the convention residual = measurement - expectation
//, together with the geometric definition of residual = hit.cross(traj)
      HepSymMatrix herr = vT_times_v(_linrel);
      HepVector rvec = herr*refvec();
      HepVector dvec = _linrel*_residual;
      _hitweight = KalWeight(rvec-dvec ,herr);
// record the flightlength change
      _fltdif = fabs(_hot->fltLen() - oldlen);
    } else { // turn off the hot
      _hot->setUsability(false);
    }
  }
// reset the site
  reset();
  return retval;
}
//  
double
KalHit::chisquared(const KalSite* site, trkDirection tdir) const {
  double chisq(-1.0);
  chisquared(chisq,site,tdir);
  return chisq;
}
double
KalHit::chisquared(const KalParams& params) const {
  double chisq(-1.0);
  chisquared(chisq,params);
  return chisq;
}
//
//  Chisquared wrt a previous site
//
bool
KalHit::chisquared(double& chisq,const KalSite* prevsite,trkDirection tdir) const {
  return chisquared(chisq,prevsite->filterParameters(tdir),false);
}
// real chisquared calculation, wrt some parameters
bool
KalHit::chisquared(double& chisq,const KalParams& params,bool ignoreactive) const {
  double chival,chierr2;
  if(chi(params,chival,chierr2,ignoreactive)){
    chisq = chival*chival/chierr2;
    return true;
  } else
    return false;
}

bool
KalHit::chi(const KalParams& params,
	    double& chival,
	    double& chierr2,
	    bool ignoreactive) const {
  bool retval(false);
// check basic things
  if(_hot->isUsable() && params.matrixOK()){
// only active hots contribute to chisquared
    if(_hot->isActive() || ignoreactive){
// compute the track parameters effective residual 'error'.  This is
// actually normalized by the hit error, as we work in chi space instead
// of residual space
      chierr2 = params.covarianceMatrix().similarity(_linrel);
//  add the unit matrix to account for the hit error
      chierr2 += 1.0;
// compute the (linearly) corrected chi
      chival = _residual + dot(_linrel,(params.parameterVector()-refvec()));
      retval = true;
    } else {
// no info here
      chival = 0.0;
      chierr2 = 1.0;
    }
  }
  return retval;
}
//
//  Activate/deactivate the HOT
//
void
KalHit::setActivity(bool newstate){
//
//  Only set the 'needs fit' flag if the state changed
//
  if(_hotstate != _hot->isActive()){
    setFit(trkIn, false);
    setFit(trkOut, false);
  }
}

void
KalHit::cloneHot(const KalRep* newrep) {
// cast off const
  _hot = _hot->clone(const_cast<KalRep*>(newrep),newrep->referenceTraj());
  assert(_hot != 0);
}

bool
KalHit::needsFit(trkDirection idir) const {
// hit needs fit if it hasn't been processed, or
// the HOT state has changed since it was processed
  return KalSite::needsFit(idir) ||
    _hot->isActive() != _hotstate;
}

unsigned
KalHit::nDof( TrkEnums::TrkViewInfo view) const {
  if(view == TrkEnums::bothView ||
     _hot->whatView() == TrkEnums::bothView ||
     _hot->whatView() == view)
    return _hot->isActive() ? 1 : 0;
  else
    return 0;
}

void
KalHit::invert() {
// Need to invert lengths for Hot, too.
  _hot->setFltLen(-1*_hot->fltLen());
// we also need to flip the ambiguity state.  This only affects
// DCH hots
  _hot->setAmbig( -1*_hot->ambig() );
 KalSite::invert();
}


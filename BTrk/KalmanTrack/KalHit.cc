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
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/KalmanTrack/KalHit.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrk/TrkBase/TrkDifPieceTraj.hh"
#include "BTrk/TrkBase/TrkParams.hh"
#include "BTrk/TrkBase/TrkErrCode.hh"
#include "CLHEP/Vector/ThreeVector.h"
#include <math.h>
#include <assert.h>
using std::endl;
using std::ostream;
using namespace CLHEP;

//
//  Constructor
//
KalHit::KalHit(const TrkDifPieceTraj* straj,TrkHit* hit) :
  KalSite(hitSite),_hit(hit),_residual(0.0),_hitstate(_hit->isActive()),
  _fltdif(0.0)
{ updateCache(straj); }

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
// record the state of the hit when processed
  _hitstate = _hit->isActive();
  return processWeight(prevsite,idir,_hitweight);
}
//
//  access
//
void
KalHit::printAll(ostream& os) const {
  os << "Hit ";
  KalSite::printAll(os);
  os << "Associated TrkHit = ";
  _hit->print(os);
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
// the hit can become self-inconsistent
//
  double oldlen=_hit->fltLen();
  bool retval(true);
  if(setTraj(reftraj,oldlen) || reftraj != _hit->trkTraj()){
    TrkErrCode residerr = updateMeasurement(*_hit,reftraj);
    if(residerr.success())
      residerr = _hit->getFitStuff(_linrel,_residual);
    retval = residerr.success();
    if (retval) {
      setTraj(reftraj,_hit->fltLen());
//  compute information quantities
// the sign here defines the convention residual = measurement - expectation
//, together with the geometric definition of residual = hit.cross(traj)
      HepSymMatrix herr = vT_times_v(_linrel);
      HepVector rvec = herr*refvec();
      HepVector dvec = _linrel*_residual;
      _hitweight = KalWeight(rvec-dvec ,herr);
// record the flightlength change
      _fltdif = fabs(_hit->fltLen() - oldlen);
    } else { // turn off the hit
      _hit->setActivity(false);
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
  if((_hit->isActive() || ignoreactive) && params.matrixOK()){
    // only active hits contribute to chisquared
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
  return retval;
}
//
//  Activate/deactivate the TrkHit
//
void
KalHit::setActivity(bool newstate){
//
//  Only set the 'needs fit' flag if the state changed
//
  if(_hitstate != _hit->isActive()){
    setFit(trkIn, false);
    setFit(trkOut, false);
  }
}

bool
KalHit::needsFit(trkDirection idir) const {
// hit needs fit if it hasn't been processed, or
// the TrkHit state has changed since it was processed
  return KalSite::needsFit(idir) ||
    _hit->isActive() != _hitstate;
}

unsigned
KalHit::nDof( ) const {
  return _hit->isActive() ? 1 : 0;
}

void
KalHit::invert() {
// Need to invert lengths for Hit, too.
  _hit->setFltLen(-1*_hit->fltLen());
// we also need to flip the ambiguity state
  _hit->setAmbig( -1*_hit->ambig() );
 KalSite::invert();
}


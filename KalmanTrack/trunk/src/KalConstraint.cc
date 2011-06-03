// ------------------------------------------------------------------------------
// File and Version Information:
//      $Id: KalConstraint.cc,v 1.26 2005/10/04 19:34:09 brownd Exp $
//
//  Description:
//  Class to describe a kalman filter constraint site.  This puts measurement information
//  into the fit
//
// Copyright Information:
//	Copyright (C) 1996	Lawrence Berkeley Laboratory
//
//  Authors: Dave Brown, 3/7/97
//------------------------------------------------------------------------------
#include "BaBar/BaBar.hh"
#include "KalmanTrack/KalConstraint.hh"
#include "KalmanTrack/KalRep.hh"
#include "TrkBase/TrkDifPieceTraj.hh"
#include "TrkBase/TrkParams.hh"
#include "TrkBase/TrkErrCode.hh"
#include "TrkBase/TrkSimpTraj.hh"
#include "CLHEP/Vector/ThreeVector.h"
#include <math.h>
#include <assert.h>
using std::endl;
using std::ostream;
//
//  Constructor
//
KalConstraint::KalConstraint(const TrkDifPieceTraj* straj,
			     const TrkParams& params,
			     bool* constrain,
			     double fltlen):
  KalSite(constraintSite),
  _cweight(params.parameter().num_row()),
  _cparams(params),
  _constrain(new bool[params.parameter().num_row()]),
  _traj(0)
{
// Initialize the array: this simply avoids uninitialized memory read
  for(unsigned ibool=0;ibool<params.parameter().num_row();ibool++)
    _constrain[ibool]=false;
// apply the masking.  This creates the weight.
  setConstraint(constrain);
// set the trajectory
  setTraj(straj,fltlen);
}
// subclass constructor
KalConstraint::KalConstraint(const TrkDifPieceTraj* straj,
			     const TrkParams& params,
			     bool* constrain,
			     double fltlen,
			     KalSite::siteType type):
  KalSite(type),
  _cweight(params.parameter().num_row()),
  _cparams(params),
  _constrain(new bool[params.parameter().num_row()]),
  _traj(0)
{
// initialize the array
  for(unsigned ibool=0;ibool<params.parameter().num_row();ibool++)
    _constrain[ibool]=false;
// apply the masking.  This creates the weight.
  setConstraint(constrain);
// set the trajectory
  setTraj(straj,fltlen);
}
//
//  copy constructor
//
KalConstraint::KalConstraint(const KalConstraint& other) :
  KalSite(other),_cweight(other._cweight),_cparams(other._cparams),
  _constrain(new bool[other._cparams.parameterVector().num_row()]),
  _traj(other._traj == 0 ? 0 : other._traj->clone())
{
// copy the constraint vector
  for(unsigned iparam=0;iparam<_cparams.parameterVector().num_row();iparam++)
    _constrain[iparam] = other._constrain[iparam];
}
// clone function
KalConstraint*
KalConstraint::clone(const KalRep* krep) const {
  KalConstraint* newconstraint = new KalConstraint(*this);
// set the trajectory to the new rep's reference trajectory
// For now, assume the length doesn't change
  newconstraint->setTraj(krep->referenceTraj(),globalLength());
  return newconstraint;
}
//
KalConstraint::~KalConstraint(){
  delete[] _constrain;
  delete _traj;
}
//
//  Update for a new reference trajectory.  The constraint doesn't
//  change, as it's in 'absolute' parameter space
//
bool
KalConstraint::update(const TrkDifPieceTraj* reftraj,double) {
  setTraj(reftraj,globalLength());
  reset();
  return true;
}
//
//  Process the constraint; this is identical to KalHit
//
bool
KalConstraint::process(const KalSite* prevsite,trkDirection idir){
  bool status(false);
  if(nDof()>0)
    status = processWeight(prevsite,idir,_cweight);
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
//
//  access
//
void
KalConstraint::printAll(ostream& os) const {
  os << "KalConstraint ";
  KalSite::printAll(os);
  os << "constraint parameters: " << endl;
  _cparams.printAll(os);
}
// chisquared wrt the previous site (calls down to parameter-based chisquared)
bool
KalConstraint::chisquared(double& chisq,const KalSite* prevsite,trkDirection tdir) const {
  if(prevsite != 0)
    return chisquared(chisq,prevsite->filterParameters(tdir));
  else
    return false;
}
// chisquared WRT parameters
bool
KalConstraint::chisquared(double& chisq,const KalParams& params) const {
  bool retval(false);
  if(params.matrixOK()) {
    chisq = _cparams.chisq(params,_constrain);
    retval = true;
  }
  return retval;
}

unsigned
KalConstraint::nDof(TrkEnums::TrkViewInfo) const {
//
// view info can't be used since TrksimpTraj doesn't provide it
  unsigned ndof(0);
  unsigned npar = _cparams.parameterVector().num_row();
  for(unsigned ipar=0;ipar<npar;++ipar)
    if(_constrain[ipar]) ndof++;
  return ndof;
}

bool
KalConstraint::setConstraint(bool* constrain) {
  bool changed(false);
  if(constrain != 0){
    unsigned npar = _cparams.parameterVector().num_row();
    for(unsigned ipar=0;ipar<npar;ipar++){
      changed |= _constrain[ipar] != constrain[ipar];
      _constrain[ipar] = constrain[ipar];
    }
  } else {
// no input; turn on all constraints
    unsigned npar = _cparams.parameterVector().num_row();
    for(unsigned ipar=0;ipar<npar;ipar++){
      changed |= _constrain[ipar] != true;
      _constrain[ipar] = true;
    }
  }
// re-computed the weight and reset the status if changed
  if(changed){
    maskWeight();
    setFit(trkIn, false);
    setFit(trkOut, false);
  }
  return changed;
}
  
void
KalConstraint::maskWeight() {
// by default, set the weight invalid
  _cweight = KalWeight();
// we must mask off the unconstrained DOF _BEFORE INVERTING_ !!!!
// mask off parameters by reducing the dimensionality
// First, compute a reduction matrix
  unsigned ndof = nDof();
  if(ndof > 0){
    unsigned npar = _cparams.nPar();
    HepMatrix reduce(ndof,npar,0);
    unsigned idof(0);
    for(unsigned ipar=0;ipar<npar;ipar++){
      if(_constrain[ipar]){
        reduce[idof][ipar] = 1.0;
        idof++;
      }
    }
// project out the unconstrained parameters
    HepSymMatrix rmat = _cparams.covarianceMatrix().similarity(reduce);
// invert the reduced covariance to get the reduced weight
    int ierr;
    rmat.invert(ierr);
    if(ierr == 0){
// expand back out the masked parameters with 0s to recover the dimensionality
      HepMatrix expand = reduce.T();
      HepSymMatrix wmat = rmat.similarity(expand);
// construct the weight.  Note that this automatically
// zeros any masked components.
      HepVector wvec(wmat*_cparams.parameterVector());
      _cweight = KalWeight(wvec,wmat);
    }
  }
}

bool
KalConstraint::isActive() const {
  return nDof() > 0;
}

const TrkSimpTraj*
KalConstraint::constraintTrajectory() const {
  if(_traj == 0){
    const KalParams& cparams = constraintParams();
    _traj = localTrajectory()->clone();
    assert(_traj != 0);
    _traj->parameters()->parameter() = cparams.parameterVector();
    _traj->parameters()->covariance() = cparams.covarianceMatrix();
  }
  return _traj;
}
      
void
KalConstraint::invert() {
// invert the original constraint
  _cparams.invert(localTrajectory());
// re-compute the weight
  maskWeight();
// call-down to base
  KalSite::invert();
}

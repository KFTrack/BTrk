// ------------------------------------------------------------------------------
// File and Version Information:
//      $Id: KalSite.cc,v 1.50 2006/04/25 04:29:53 brownd Exp $
//
//  Description:
//  Class to describe a generic Kalman filter 'site'.
//
// Copyright Information:
//	Copyright (C) 1996	Lawrence Berkeley Laboratory
//
//  Authors: Dave Brown, 12/18/96
//------------------------------------------------------------------------------
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/BaBar/Constants.hh"
#include "BTrk/KalmanTrack/KalSite.hh"
#include "BTrk/TrkBase/TrkDifPieceTraj.hh"
#include "BTrk/TrkBase/TrkSimpTraj.hh"
#include "BTrk/TrkBase/TrkParams.hh"
#include "BTrk/TrkBase/HelixTraj.hh"
#include "CLHEP/Matrix/Vector.h"
using std::endl;
using std::ostream;
using namespace CLHEP;
//
KalSite::KalSite(const siteType stype):
  _stype(stype),_globlen(0.0),_loclen(0.0),_loctraj(0),_gap(0.0)
{
  _siteflag[trkOut] = _siteflag[trkIn] = false;
}
//
//  Copy constructor
//
KalSite::KalSite(const KalSite& other) :
  _stype(other._stype),_globlen(other._globlen),_loclen(other._loclen),
  _loctraj(other._loctraj),_gap(other._gap)
{
  for(int idir=0;idir<2;idir++){
    _siteflag[idir] = other._siteflag[idir];
    trkDirection tdir = idir == trkIn ? trkIn : trkOut;
    params(tdir) = other.params(tdir);
    weight(tdir) = other.weight(tdir);
  }
}
//
KalSite::~KalSite(){;}

// reset the site (delete the cache, reset flags)
void
KalSite::reset(){
  for(int idir=0;idir<2;idir++){
    static KalParams emptyPar;
    static KalWeight emptyWeight;
    _params[idir] = emptyPar;
    _weight[idir] = emptyWeight;
    _siteflag[idir] = false;
  }
}

void
KalSite::reset(trkDirection idir){
  static KalParams emptyPar;
  static KalWeight emptyWeight;
  _params[idir] = emptyPar;
  _weight[idir] = emptyWeight;
  _siteflag[idir] = false;
}

//
//  simple printout
//
void
KalSite::printAll(ostream& os) const {
  os << "Kalman Site at flightlength " << _globlen << endl;
//  double olen = localLength() - globalLength();
//  _loctraj->printAll(os);
  if(needsFit(trkOut)) os << "Outer result needs fit " << endl;
  else {
    os << "Outer result" << endl;
// report results at the origin
    filterParameters(trkOut).printAll(os);
  }
  if(needsFit(trkIn)) os << "Inner result needs fit " << endl;
  else {
    os << "Inner result" << endl;
    filterParameters(trkIn).printAll(os);
  }
}
void
KalSite::print(ostream& os) const {
  if(type() == hitSite)
    os << "Kalman Hit Site ";
  else if (type() == matSite)
    os << "Kalman Mat Site ";
  else if (type() == bendSite)
    os << "Kalman Bnd Site ";
  os <<" at flightlength " << _globlen;
  if(needsFit(trkOut)) os << " Outer needs fit ";
  else os << " Outer  has  fit ";
  if(needsFit(trkIn)) os << " Inner needs fit" << endl;
  else os << " Inner  has  fit" << endl;
}
// set the trajectory parameters given the piece traj and length
bool
KalSite::setTraj(const TrkDifPieceTraj* reftraj,double globlen) {
  _globlen = globlen;
  _loctraj = reftraj->localTrajectory(_globlen,_loclen);
  assert(_loctraj != 0);
// set return value before overwriting parameters
  bool retval = ! (_lparams == _loctraj->parameters()->parameter());
  _lparams = _loctraj->parameters()->parameter();
  return retval;
}


void
KalSite::copyParams(const KalSite* other,trkDirection idir) {
// copy the other sites parameters
  params(idir) = other->filterParameters(idir);
}

void
KalSite::copyWeights(const KalSite* other,trkDirection idir) {
  weight(idir) = other->filterWeights(idir);
}

void
KalSite::copySite(const KalSite* other,trkDirection idir) {
  if(other->params(idir).matrixOK())
    params(idir) = other->params(idir);
  if(other->weight(idir).matrixOK())
    weight(idir) = other->weight(idir);
}

bool
KalSite::validSite(trkDirection idir) const {
  return params(idir).matrixOK() || weight(idir).matrixOK();
}

const KalWeight& 
KalSite::filterWeights(trkDirection idir) const {
  if( (! weight(idir).matrixOK()) && params(idir).matrixOK())
    weight(idir) = params(idir);
  return _weight[idir];
}

const KalParams&
KalSite::filterParameters(trkDirection idir) const {
  if( (! params(idir).matrixOK()) && weight(idir).matrixOK())
    params(idir) = weight(idir);
  return _params[idir];
}


bool
KalSite::processParams(const KalSite* prevsite,
		       trkDirection tdir,
		       const KalParams& transport) {
  bool status(false);
  invalidateSite(tdir);
// previous site must have been processed
  if(prevsite != 0 && prevsite->hasFit(tdir)){
    if(isActive()){
// copy over the previous site's parameters 
      copyParams(prevsite,tdir);
      if(params(tdir).matrixOK() && transport.matrixOK()){
// sign of transport is given by direction; by convention, inwards is additive, outwards subtractive (like momentum)
        if(tdir == trkIn)
          params(tdir).addEffect(transport);
        else
          params(tdir).subEffect(transport);
// invalidate weights to force re-computation (if necessary)
        weight(tdir).invalidate();
// set fit flag
        setFit(tdir);
        status = true;
      }
    } else {
// copy the previous site
      copySite(prevsite,tdir);
      if(validSite(tdir)){
// set the site processed in this direction
	setFit(tdir);
	status = true;
      }
    }
  }
  return status;
}

bool
KalSite::processWeight(const KalSite* prevsite,
		       trkDirection tdir,
		       const KalWeight& hweight) {
// Assume failure
  bool status(false);
  invalidateSite(tdir);
// previous site must have been processed
  if(prevsite != 0 && prevsite->hasFit(tdir)){
    if(isActive()){
// copy the weights
      copyWeights(prevsite,tdir);
//  Update the information parameters and matrix
      if(weight(tdir).matrixOK() && hweight.matrixOK()){
	weight(tdir) +=  hweight;
// invalidate parameters to force recomputation
	params(tdir).invalidate();
// set the site processed in this direction
	setFit(tdir);
// set the fit status
	status = true;
      }
    } else {
// copy the previous site
      copySite(prevsite,tdir);
      if(validSite(tdir)){
// set the site processed in this direction
	setFit(tdir);
	status = true;
      }
    }
  }
  return status;
}


void
KalSite::mergeParams(const KalSite* other,KalParams& mergeparams) const {
  KalWeight weight;
  addWeights(other,weight);
  mergeparams = KalParams(weight);
}

void
KalSite::addParams(const KalSite* other,KalParams& addparams) const {
  trkDirection tdir = relativeDirection(other);
  trkDirection otherdir = tdir==trkIn ? trkOut : trkIn;
  addparams = filterParameters(tdir);
  addparams += other->filterParameters(otherdir);
}

double
KalSite::parameterChisq(const KalSite* other,bool* tparams) const {
  trkDirection tdir = relativeDirection(other);
  trkDirection otherdir = tdir==trkIn ? trkOut : trkIn;
  const KalParams& params = filterParameters(tdir);
// have to translate the parameters in order to add weights
  const KalParams& nparams(other->filterParameters(otherdir));
  return params.chisq(nparams,tparams);
}

void
KalSite::addWeights(const KalSite* other,KalWeight& addweights) const {
  trkDirection tdir = relativeDirection(other);
  trkDirection otherdir = tdir==trkIn ? trkOut : trkIn;
  addweights = filterWeights(tdir);
  addweights += other->filterWeights(otherdir);
}

bool
KalSite::setTrajState(trkDirection tdir, TrkSimpTraj* straj) const {
// find the parameters for the specified direction
  const KalParams& params = filterParameters(tdir);
  bool retval(params.matrixOK());
  if(retval){
// this is really ugly: this should be done through TrkSimpTraj
// operators
    *(straj->parameters()) = params.trackParameters();
  }
  return retval;
}

unsigned
KalSite::nDof( TrkEnums::TrkViewInfo view) const {
  return 0;
}

bool
KalSite::chisquared(double& ,const KalSite*, trkDirection) const {
  return false;
}

bool
KalSite::isActive() const {
  return true;
}

KalSite&
KalSite::operator = (const KalSite& other ) {
  if(&other != this){
    _stype = other._stype;
    _globlen = other._globlen;
    _loclen = other._loclen;
    _loctraj = other._loctraj;
    _gap = other._gap;
    for(int idir=0;idir<2;idir++){
      _siteflag[idir] = other._siteflag[idir];
      trkDirection tdir = idir == trkIn ? trkIn : trkOut;
      params(tdir) = other.params(tdir);
      weight(tdir) = other.weight(tdir);
    }
  }
  return *this;
}

KalSite::KalSite() :
  _stype(unknown),_globlen(0.0),_loclen(0.0),
  _loctraj(0) {
   for(int idir=0;idir<2;idir++){
     _siteflag[idir] = false;
   }
}

void
KalSite::invert() {
  _globlen = -_globlen;
  _loclen = -_loclen;
// none of the cache is valid: reset it
  reset(trkIn);
  reset(trkOut);
}

void
KalSite::processDeltaP(KalParams& params,double dpfract) const {
 // here we have to process the parameter change analytically as the brem photons can be so large
// that linear approximation fails.
// first, bring in the previous site's parameters
  HepVector& pvec = params.parameterVector();
// remember the old paramters for future calculations
  HepVector oldpvec(pvec);
// I Must interpret the parameters physically.  This needs to be separated out into
// a dedicated service in future.
// first, remember a few of the input parameters
// tandip is unaffected
// omega is scaled by the energy.  Apply the convention that the reference momentum is after the brem
  double omfactor(1.0);
  omfactor = 1.0/(1.0 + dpfract);
  double oldomega = oldpvec[HelixTraj::omegaIndex];
  double omega = oldomega*omfactor;
  double oldd0 = oldpvec[HelixTraj::d0Index];
  pvec[HelixTraj::omegaIndex] = omega;
// compute new phi0
  double dominv = 1.0/oldomega - 1.0/omega;
  double offset = 1.0/oldomega + oldd0;
  double oldphi0 = oldpvec[HelixTraj::phi0Index];
  double tandip =  oldpvec[HelixTraj::tanDipIndex];
  double cosdip =  1./sqrt(1.+sqr(tandip));
  double gamphi = oldphi0 + localLength()*oldomega*cosdip;
  double yprime = sin(gamphi)*dominv - sin(oldphi0)*offset;
  double xprime = cos(gamphi)*dominv - cos(oldphi0)*offset;
  double phi0 = atan2(yprime,xprime);
// first, check for wrapping of phi
  if(fabs(phi0-oldphi0)>Constants::pi)
    phi0 += phi0 < oldphi0 ? Constants::twoPi : -Constants::twoPi;
// choose the new phi0 angle phase to be close to the old.  This avoids
// problems due to the pi-periodicity of the tan function
  if(fabs(phi0-oldphi0)>Constants::pi/2.0)
    phi0 += phi0 < oldphi0 ? Constants::pi : -Constants::pi;
  pvec[HelixTraj::phi0Index] = phi0;
// compute d0; be careful to choose an appropriate equation to avoid numerical problems
  double dphi = gamphi-phi0;
  double olddphi = gamphi-oldphi0;
  double d0;
  if(fabs(dphi)>0.01 && fabs(olddphi)>0.01){
    d0 = (1.0/oldomega + oldd0)*sin(olddphi)/sin(dphi) - 1.0/omega;
  } else {
    double cdphi = cos(dphi);
    double oldcdphi = cos(olddphi);
    d0 = ( (1.0-cdphi)/omega - (1.0-oldcdphi)/oldomega + oldd0*oldcdphi)/cdphi;
  }
  pvec[HelixTraj::d0Index] = d0;
// z0 is derived from the flight lengths of the local trajectory
  double newflt = (gamphi - phi0)/omega;
  double oldflt = (gamphi - oldphi0)/oldomega;
  pvec[HelixTraj::z0Index] = oldpvec[HelixTraj::z0Index] + tandip*(oldflt-newflt);
}

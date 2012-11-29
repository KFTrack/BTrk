// File and Version Information:
//      $Id: KalRep.cc,v 1.303 2008/09/08 22:05:48 brownd Exp $
//// Description:
//      Implementation class for TrkRep for a Kalman fit
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

#include "BaBar/BaBar.hh"
#include "BaBar/BbrCollectionUtils.hh"
#include <math.h>
#include <algorithm>
#include "KalmanTrack/KalRep.hh"
#include "KalmanTrack/KalHit.hh"
#include "KalmanTrack/KalMaterial.hh"
#include "KalmanTrack/KalEndSite.hh"
#include "KalmanTrack/KalBend.hh"
#include "KalmanTrack/KalSmear.hh"
#include "KalmanTrack/KalConstraint.hh"
#include "KalmanTrack/KalStub.hh"
#include "KalmanTrack/KalCodes.hh"
#include "DetectorModel/DetSet.hh"
#include "DetectorModel/DetMaterial.hh"
#include "TrkBase/TrkHitOnTrk.hh"
#include "TrkBase/TrkHotList.hh"
#include "TrkBase/TrkDifTraj.hh"
#include "TrkBase/HelixParams.hh"
#include "TrkBase/TrkVolume.hh"
#include "TrkBase/TrkSimpTraj.hh"
#include "TrkBase/HelixTraj.hh"
#include "CLHEP/Vector/ThreeVector.h"
#include "TrkBase/TrkMomCalculator.hh"
#include "TrkBase/TrkSimpTraj.hh"
#include "TrkBase/TrkHotListFull.hh"
#include "TrkBase/TrkHotList.hh"
#include "TrkBase/TrkHitUse.hh"
#include "ErrLogger/ErrLog.hh"
//#include "DchGeom/DchDetector.hh"

#include <algorithm>
#include <vector>
#include <deque>
using std::endl;
using std::ostream;
using std::min;
using std::max;

// predicates for site updating

struct FindMatBendSites {
  bool operator () ( KalSite* site ) const {
    return site->type() != KalSite::bendSite &&
      site->type() != KalSite::matSite; }
};


void
KalRep::init(const HelixParams& inPar)
{
//  Construct the seed traj as a helix from the input parameters, using the origin
  _seedtraj = new HelixTraj(inPar);
  initFromSeed();
}

// initialize from the seed traj;
void
KalRep::initFromSeed() {
	initRep();
	initSites();
}

void
KalRep::initRep() {
// sort the hots
  hotList()->sort();
// reserve a nominal size for sites
  _sites.reserve(128);
  setMultScat(true);   // in TrkFitStatus
  _siteflag[trkIn] = _siteflag[trkOut] = false;
  _chisq[trkIn] = _chisq[trkOut] = -1.;
// set the initial fitrange according to the hots
  if(hotList()->nActive() > 0){
    _fitrange[0] = hotList()->startFoundRange()-_kalcon.fltEpsilon();
    _fitrange[1] = hotList()->endFoundRange()+_kalcon.fltEpsilon();
  } else {
// not hots! set a small but non-null range
    _fitrange[0] = 0.0;
    _fitrange[1] = 1.0;
  }
  _extendable[trkIn] = _extendable[trkOut] = true; // by default, can be extended
//  set the seed flight range
  _seedtraj->setFlightRange(_fitrange);
// set the seed mom to the begining of the hit range
  _refmomfltlen = _fitrange[0];
// compute the seed momentum from this (+ the bfield)
  Hep3Vector momvec = TrkMomCalculator::vecMom(*_seedtraj,_kalcon.bField(),_refmomfltlen);
  _refmom = momvec.mag();
  _charge = TrkMomCalculator::charge(*_seedtraj,_kalcon.bField(),_refmomfltlen);
// build the initial reference trajectory
  if(_reftraj == 0)buildRefTraj();
// Initialize _ptraj to _reftraj. This is needed
// by UpdateFitStuff, which calls down to momentum(), which needs the traj
  _ptraj = _reftraj;
}
void
KalRep::initSites() {
// build the hit sites
  buildHitSites();
// build the KalMaterial sites
  if(_kalcon.materialSites()){
		std::vector<DetIntersection> tlist; tlist.reserve(64);
		buildMaterialSites(_fitrange,tlist);
  }
// build the bend sites
  if(_kalcon.bendSites()){
    buildBendSites(_fitrange);
  }
// re-find hit sites
  findHitSites();
// create the end sites from the helix.  increase the smearing
  updateEndSites(_kalcon.smearFactor(),true);
//  On construction, the track is neither current or valid
}

// construct from hots and intersections
KalRep::KalRep(const TrkSimpTraj& seed, TrkHotList* hotl,
	const std::vector<DetIntersection>& dinters,
	const KalContext& context,
	TrkParticle const& tpart) : TrkRep(hotl,tpart,true), _maxdist(0),_maxfltdif(0),
	_niter(0),_ninter(0),_ptraj(0),_reftraj(0),
	_kalcon(context),
	_seedtraj((TrkSimpTraj*)(seed.clone())),
	_stopsite(0)
{
// basic initialization
  initRep();
// build the hit sites
  buildHitSites();
// build the KalMaterial sites
  if(_kalcon.materialSites()){
// convert the intersections to use the new reftraj;
		std::vector<DetIntersection> tlist(dinters);
		for(std::vector<DetIntersection>::iterator iinter = tlist.begin();iinter!= tlist.end();iinter++)
			iinter->trajet = _reftraj;
		buildMaterialSites(_fitrange,tlist);
  }
// build the bend sites
  if(_kalcon.bendSites()){
    buildBendSites(_fitrange);
  }
// re-find hit sites
  findHitSites();
// create the end sites from the helix.  increase the smearing
  updateEndSites(_kalcon.smearFactor(),true);
//  On construction, the track is neither current or valid
}

// construct from hots and intersections and reference trajectory
KalRep::KalRep(const TrkDifPieceTraj* rtraj, TrkHotList* hotl,
	const std::vector<DetIntersection>& dinters,
	const KalContext& context,
	TrkParticle const& tpart) : TrkRep(hotl,tpart,true), _maxdist(0),_maxfltdif(0),
	_niter(0),_ninter(0),_ptraj(0),_reftraj(0),
	_kalcon(context),
	_seedtraj(0),
	_stopsite(0)
{
// clone the ref traj
  _reftraj = rtraj->clone();
  assert(_reftraj != 0);
// create the seed from a piece of the piecetraj
  double locflt;
  double midflt = 0.5*(_reftraj->lowRange() + _reftraj->hiRange());
  const TrkSimpTraj* loctraj = _reftraj->localTrajectory(midflt, locflt);
  if(loctraj != 0)_seedtraj = loctraj->clone();
  assert(_seedtraj != 0);
// basic initialization
  initRep();
// build the hit sites
  buildHitSites();
// build the KalMaterial sites
  if(_kalcon.materialSites()){
// convert the intersections to use the new reftraj;
    std::vector<DetIntersection> tlist(dinters);
    for(std::vector<DetIntersection>::iterator iinter = tlist.begin();iinter!= tlist.end();iinter++)
      iinter->trajet = _reftraj;
    buildMaterialSites(_fitrange,tlist);
  }
  // build the bend sites
  if(_kalcon.bendSites()){
    buildBendSites(_fitrange);
  }
// re-find hit sites
  findHitSites();
// create the end sites from the helix.  increase the smearing
  updateEndSites(_kalcon.smearFactor(),true);
//  On construction, the track is neither current or valid
}



//
//  Copy constructor to change mass tpart. Note that the copy is ALWAYS invalid and
//  must be fit.  No protection is provided against copying to the same mass tpart,
//  this is provided by the 'cloneNewHypo' function.
//
KalRep::KalRep(const KalRep& other,TrkParticle const& tpart) :
  TrkRep(other,tpart),
  _maxdist(0),_maxfltdif(0),  _niter(0), _ninter(0),
  _ptraj(0),
  _reftraj(other._reftraj->clone()),
  _kalcon(other._kalcon),
  _seedtraj(other._seedtraj->clone()),
  _stopsite(0),
  _refmom(1.0),
  _refmomfltlen(other._refmomfltlen),_charge(other._charge)
{
  _sites.reserve(128);
// create a new hotlist
  _hotList.reset( new TrkHotListFull );
  setValid(false);
  setMultScat(true);   // in TrkFitStatus
//  set direction-dependent information 
  for(int idir=0;idir<2;idir++){
    _siteflag[idir] = false;
    _extendable[idir] = other._extendable[idir];
    _chisq[idir] = -1.0;
    _fitrange[idir] = other._fitrange[idir];
  }
// remove any effect of truncation due to energy loss in the previous fit
  if(other._stopsite != 0)
    _extendable[trkOut] = true;
// Initialize _ptraj to _reftraj
  _ptraj = _reftraj;
// reset the initial momentum to the seed value: this makes sure lower-mass states aren't
// found to 'stop'
  Hep3Vector momvec = TrkMomCalculator::vecMom(*_seedtraj,_kalcon.bField(),_refmomfltlen);
  _refmom = momvec.mag();
// clone-Copy the sites
  for(unsigned isite=0;isite<other._sites.size();isite++){
    KalSite* newsite = other._sites[isite]->clone(this);
    assert(newsite != 0);
    _sites.push_back(newsite);
// append the HOTS to the hotlist
    if(newsite->kalHit() != 0)
      hotList()->append(newsite->kalHit()->hitOnTrack());
  }
// Locate the first and last hit site
  findHitSites();
  updateEndSites(_kalcon.smearFactor(),true);
}
//
//  Copy constructor to move rep to a new track.
//
KalRep::KalRep(const KalRep& other ) : 
  TrkRep(other,other.particleType()),
  _maxdist(other._maxdist),
  _maxfltdif(other._maxfltdif),
  _niter(other._niter),
  _ninter(other._ninter),
  _ptraj(0),
  _reftraj(other._reftraj->clone()),
  _kalcon(other._kalcon),
  _seedtraj( other._seedtraj->clone()),
  _stopsite(0),
  _refmom(other._refmom),
  _refmomfltlen(other._refmomfltlen),_charge(other._charge)
{
  setMultScat(true);   // in TrkFitStatus
// if a piecetraj exists, clone it
  if(other._ptraj == other._reftraj)
    _ptraj = _reftraj;
  else if(other._ptraj != 0)
    _ptraj = other._ptraj->clone();
// create a new hotlist
  _hotList.reset( new TrkHotListFull );
//  Copy direction-dependent information 
  for(int idir=0;idir<2;idir++){
    _siteflag[idir] = other._siteflag[idir];
    _extendable[idir] = other._extendable[idir];
    _chisq[idir] = other._chisq[idir];
    _hitrange[idir] = other._hitrange[idir];
    _fitrange[idir] = other._fitrange[idir];
  }
// clone-Copy the sites
  for(unsigned isite=0;isite<other._sites.size();isite++){
    KalSite* oldsite = other._sites[isite];
    KalSite* newsite = oldsite->clone(this);
    assert(newsite != 0);
    _sites.push_back(newsite);
// append the HOTS to the hotlist
    if(newsite->kalHit() != 0)
      hotList()->append(newsite->kalHit()->hitOnTrack());
  }
// Locate the first and last hit site
  findHitSites();
// create end sites
  updateEndSites(_kalcon.smearFactor());
}
//  Destructor
KalRep::~KalRep(){
//  Delete all the KalSites
  std::for_each(_sites.begin(),_sites.end(),babar::Collection::DeleteObject());
// Delete the trajectories
  delete _seedtraj;
  delete _reftraj;
  if(_reftraj != _ptraj)
    delete _ptraj;
// delete integrator
}
//
//  Clone operator
//
KalRep*
KalRep::clone() const {
  return new KalRep(*this);
}

TrkRep*
KalRep::cloneNewHypo(TrkParticle const& tpart) {
  if(tpart != particleType()) {
    KalRep* newrep = new KalRep(*this, tpart);
    newrep->setValid(false);
    return newrep;
  }
  else
    return this;
}

//
//  Return exchange par of rep.
//
HelixParams
KalRep::helix(double fltlen) const {
  double locflight;
  const TrkSimpTraj* ltraj = localTrajectory(fltlen,locflight);
// make sure the reference point is at the origin
// set origin to a nominal value DNB_RKK
  static const HepPoint origin(0.0,0.0,0.0);
  if(ltraj->referencePoint() == origin){
    return HelixParams(ltraj->parameters()->parameter(), 
			  ltraj->parameters()->covariance());
  } else {
// clone the trajectory, and change its reference point to the origin
    TrkSimpTraj* otraj = ltraj->clone();
    double olen = locflight - fltlen;
    otraj->changePoint(origin,olen);
    HelixParams outpar(otraj->parameters()->parameter(), 
			  otraj->parameters()->covariance());
    delete otraj;
    return outpar;
  }
}
//
//  Chisquared function; allow this to work even with 1-directional fits
double
KalRep::chisq() const {
  if(hasFit(trkIn))
    return chisquared(trkIn);
  else
    return chisquared(trkOut);
}
//
double
KalRep::chisquared(trkDirection tdir) const {
//  See if the cache is valid (it's self-flagging)
  if(_chisq[tdir] < 0.0) {
    if(hasFit(tdir)){
      int nsites = _sites.size();
      int startsite = tdir==trkOut ? 0 : nsites-1;
//    must cast-off const to use the cache
      double& chisum = (double&)_chisq[tdir];
      chisum = chisquared(startsite,nsites,tdir);
    }
  }
  return _chisq[tdir];
}

double
KalRep::chisquared(int startsite,int nsites,trkDirection tdir) const {
  int sitestep = tdir==trkOut ? 1 : -1;
  int isite(startsite);
  const KalSite* prevsite = _sites[isite];
  nsites--; // 1st site doesn't contribute to chisquared
  double chisum(0.0); // zero chisquared
  double chisq(0.0);
  while(nsites>0){
    isite += sitestep;
    KalSite* thesite = _sites[isite];
    if(thesite->chisquared(chisq,prevsite,tdir))
      chisum += chisq;
    prevsite = thesite;
    nsites--;
  }
  return chisum;
}

//
double
KalRep::chisquared(double fltlen,trkDirection tdir) const {
// first, find the sites; I really want 1 past the one found by this function
  int nearest = findNearestSite(fltlen) + 1 ;
//  Determine the start for this direction
  int startsite = tdir==trkOut ? 0 : _sites.size()-1;
  int nsites = tdir==trkOut ? nearest : _sites.size() - nearest;
  return chisquared(startsite,nsites,tdir);
}

int
KalRep::nDof() const{
  return nDof(TrkEnums::bothView);
}

int
KalRep::nDof(double fltlen,trkDirection tdir) const {
// first, find the sites; I really want 1 past the one found by this function
  int nearest = findNearestSite(fltlen) + 1;
//  Determine the start for this direction
  int isite = tdir==trkOut ? 0 : _sites.size()-1;
  int nsites = tdir==trkOut ? nearest : _sites.size() - nearest;
  int sitestep = tdir==trkOut ? 1 : -1;
  int dof(0);
  while(nsites>0){
    dof += _sites[isite]->nDof();
    isite += sitestep;
    nsites--;
  }
// correct for the parameter DOFs
  dof -= _seedtraj->parameters()->parameter().num_row();
  return dof;
}

int 
KalRep::nDof(TrkEnums::TrkViewInfo view) const {
  int dof(0);
  for(unsigned isite= 0;isite<_sites.size();isite++)
    dof += _sites[isite]->nDof(view);
// subtract the parameters (for overall DOFs)
  if(view == TrkEnums::bothView)
    dof -= _seedtraj->parameters()->parameter().num_row();
  return dof;
}

bool
KalRep::enoughDofs() const {
// initialize separate counting in each view
  int dof(-_seedtraj->parameters()->parameter().num_row());
// TrkSimpTraj doesn't know how to provide this for different
// views; too bad
  int xydof(0);
  int zdof(0);
  unsigned isite(0);
  unsigned nsites(_sites.size());
  bool retval(false);
  while(isite<nsites && !retval) {
    dof += _sites[isite]->nDof(TrkEnums::bothView);
    xydof += _sites[isite]->nDof(TrkEnums::xyView);
    zdof += _sites[isite]->nDof(TrkEnums::zView);
    isite++;
    retval = 
      dof >= (int)_kalcon.minDOF(TrkEnums::bothView) &&
      xydof >= (int)_kalcon.minDOF(TrkEnums::xyView) &&
      zdof >= (int)_kalcon.minDOF(TrkEnums::zView);
  }
  return retval;
}

double
KalRep::hitChisq(const TrkHitUse& hituse) const {
  double chisq = -1.0;
  if(fitValid()){
// create an 'end site' from the trajectory
    KalEndSite trajsite(_ptraj,hituse.fltLen(),trkIn,1.0,false);
// create a HOT for this hit usage
    TrkHitOnTrk* hot = hituse.createHitOnTrk(*this);
    if(hot != 0) {
      KalHit kalhit(_ptraj,hot);
// compute the chisquared using this hit
      kalhit.chisquared(chisq,&trajsite,trkIn);
      kalhit.deleteHOT();
    }
  }
  return chisq;
}

double
KalRep::hitChisq(const TrkHitOnTrk* hot,bool exclude) const {
  double chisq = -1.0;
  double chival,chierr2;
  if(hitChi(hot,chival,chierr2,exclude))
    chisq = chival*chival/chierr2;
  return chisq;
}

bool
KalRep::hitChi(const TrkHitOnTrk* hot,
	       double& chival,
	       double& chierr2,
	       bool exclude) const {
  bool retval(false);
  if(fitValid()){
// find the hot site
    unsigned index;
    const KalHit* hitsite = findHotSite(hot,index);
    if(hitsite != 0){
      const KalSite* refsite(hitsite);
      KalParams smoothed;
// merge the sites on either side of the hot (if necessary)
      if(index > _hitrange[0] && index < _hitrange[1] ){
	const KalSite* prevsite(0);
	const KalSite* nextsite(0);
	if(exclude){
	  prevsite = nextActive(index,trkIn);
	  nextsite = nextActive(index,trkOut);
	  refsite = prevsite;
	} else {
	  prevsite = hitsite;
	  nextsite = nextActive(index,trkOut);
	}
// merge the inner and outer parameters
	if(prevsite != 0 && nextsite != 0)
	  prevsite->mergeParams(nextsite,smoothed);
      } else if(index <= _hitrange[0] ) {
	if(exclude)
	  refsite = nextActive(index,trkOut);
	if(refsite != 0)
	  smoothed = refsite->filterParameters(trkIn);
      } else if(index >= _hitrange[1] ) {
	if(exclude)
	  refsite = nextActive(index,trkIn);
	if(refsite != 0)
	  smoothed = refsite->filterParameters(trkOut);
      }
// index outside of hitrange is implicitly a (soft) error
      if(smoothed.matrixOK()){
	retval = hitsite->chi(smoothed,chival,chierr2);
      }
    }
  }
  return retval;
}

double
KalRep::hitChisq(const TrkHitOnTrk* hot,trkDirection tdir) const {
  double chisq = -1.0;
  if(!needsFit(tdir)){
// find the hot site
    unsigned index;
    const KalHit* hitsite = findHotSite(hot,index);
    if(hitsite != 0 && ( (tdir == trkIn && index < _sites.size()-1) ||
			 (tdir == trkOut && index > 0 ))){
      trkDirection odir = tdir == trkIn ? trkOut : trkIn;
      const KalSite* prevsite = nextActive(index,odir);
// compute the chisquared
      if(prevsite != 0)
	hitsite->chisquared(chisq,prevsite,tdir);
    }
  }
  return chisq;
}
double
KalRep::matchChisq(double fltlen,bool* tparams) const {
  double chisqFB = -1.0;
//  Find the sites which bound this flightlength
  const KalSite* insite;
  const KalSite* outsite;
  findBoundingSites(fltlen,insite,outsite);
  if(fitValid() && insite!=0 && outsite!= 0)
    chisqFB = insite->parameterChisq(outsite,tparams);
  return chisqFB;
}

bool
KalRep::resid(const TrkHitOnTrk *hot, 
              double& resval, double& reserr,
              bool exclude) const
{
   double chival,chierr2;
   bool b=hitChi(hot,chival,chierr2,exclude); 
   if (b) {
     resval=chival*hot->hitRms();
     reserr=sqrt(chierr2)*hot->hitRms();
   }
   return b;
}

//
//  Print functions;
//
void
KalRep::printAll(ostream& ostr) const {
  ostr << "Kalman track Rep, seeded from ";
  _seedtraj->printAll(ostr);
  ostr << "Reference trajectory = " << endl;
  _reftraj->printAll(ostr);
  if(needsFit(trkOut))
    ostr << " needs outwards fit ";
  else
    ostr << " has   outwards fit";
  if(needsFit(trkIn))
    ostr << " and needs inwards fit ";
  else
    ostr << " and has   inwards fit";
  ostr << endl;
  ostr << "There are " << _sites.size()
    << " KalSites on this rep, as follows: " << endl;
  for(unsigned isite=0;isite<_sites.size();isite++){
    ostr << "Site " << isite << " : ";
    _sites[isite]->printAll(ostr);
  }
  if(_ptraj != 0) {
    ostr << "The piecetraj for this rep is as follows:" << endl;
    _ptraj->printAll(ostr);
  } else
    ostr << "No piecetraj for this rep" << endl;
}
//
void
KalRep::print(ostream& ostr) const {
  ostr << "Kalman track Rep has " << _sites.size() << " KalSites and ";
  if(needsFit(trkOut))
    ostr << "needs outwards fit ";
  else
    ostr << "has   outwards fit";
  if(needsFit(trkIn))
    ostr << " and needs inwards fit ";
  else
    ostr << " and has   inwards fit";
}
//
//  Full fit function.  This will fit the rep from scratch, clearing out any
//  previous results.  It uses the existing trajectory.
//
TrkErrCode
KalRep::fit(){
// if the fit is already current, no need to do anything
  if(fitCurrent())
    return TrkErrCode(TrkErrCode::succeed,KalCodes::current,
		      "Fit is already current");
//
  TrkErrCode fiterr = iterateFit();
  if(fiterr.success()){
// Extend the fit before and after the last hit
    fiterr = extendSites(_hitrange[0],trkIn);
    if(fiterr.success())
      fiterr = extendSites(_hitrange[1],trkOut);
//  set fit as valid and current
    if(fiterr.success()){
      setValid(true);
      setCurrent(true);
// The fit iterator leaves _niter at one more than maxIterations 
// when there is no convergence
      if(_niter == _kalcon.maxIterations() + 1)
	fiterr = TrkErrCode(TrkErrCode::succeed,2);
    }
// re-compute the charge using the fit result.  Rarely this will differ from the seed charge
    double loclen(0.0);
    const TrkSimpTraj* loctraj = localTrajectory(_refmomfltlen,loclen);
    if(loctraj == 0){
      ErrMsg(error) << "Can't find local trajectory for charge measurement!" << endmsg;
      loctraj = _seedtraj;
    }
    _charge = TrkMomCalculator::charge(*loctraj,_kalcon.bField(),loclen);
  }  else {
   // make sure a failed fit is neither valid nor current
       setValid(false);
       setCurrent(false);
  }
  // update the track itself for the new status: this is only relevant in mini-refit jobs
  setValid(fitValid());
  setCurrent(fitCurrent());
  return fiterr;
}
// iteration, used to converge in the fit
TrkErrCode
KalRep::iterateFit(){
  TrkErrCode fiterr(TrkErrCode::fail);
  while(_niter <= _kalcon.maxIterations() && !converged()){
// iteration sequence: update, fit, build the trajectory
// Check status after every step
    if (ErrLogging(debugging))
      std::cout<<"fit iteration "<<_niter<<std::endl;
    _niter++;
    if(_niter>1){
      fiterr = updateSites();
      if (ErrLogging(debugging))
	std::cout<<"update Sites err="<<fiterr<<std::endl;
      if(fiterr.failure())return fiterr;
    }
    fiterr = fitSites();
    if (ErrLogging(debugging)){
      std::cout<<"fit Sites err="<<fiterr<<std::endl;
      printAll(std::cout);
    }
    if(fiterr.failure())return fiterr;
    fiterr = buildTraj();
    if (ErrLogging(debugging))
      std::cout<<"build Traj err="<<fiterr<<std::endl;
    if(fiterr.failure())return fiterr;
  }
  return fiterr;
}
// 'fit' a single direction.  No trajectory is produced.  This function is only useful
// for creating KalStubs
TrkErrCode
KalRep::fit(trkDirection tdir){
  TrkErrCode retval;
// if the fit is already current, no need to do anything
  if(needsFit(tdir)){
// check that the track is fitable
    if(isFitable(retval)){
// make sure sites are OK
      fixupSites();
// process _all_ the sites
      bool status = tdir == trkOut? process(&_endsites[trkIn],0,_sites.size()-1,tdir)
    	: process(&_endsites[trkOut],_sites.size()-1,0,tdir);
      if(status) {
	setFit(tdir);
// build the 1-d traj
	retval = buildTraj(tdir);
	if(retval.success())
	  setValid(true); // the track now has a 'valid' fit
      } else
	retval = TrkErrCode(TrkErrCode::fail);
    }
  } else
    retval = TrkErrCode(TrkErrCode::succeed,KalCodes::valid,
			"Fit in specified direction is already valid");
  return retval;
}
//
//  allow the fit range to be changed
//
void
KalRep::setFitRange(double newrange[2]){
  _fitrange[0] = newrange[0];
  _fitrange[1] = newrange[1];
}
void
KalRep::setFitRange(double lowrange, double hirange){
  if( lowrange < hirange ) {
    _fitrange[0] = lowrange;
    _fitrange[1] = hirange;
  } else {
    _fitrange[0] = hirange;
    _fitrange[1] = lowrange;
  }
}
//
//  build the piecewise trajectory from the processed sites.  This constructs
//  the smoothed parameters from between the sites.  The sites must have been
//  already fit in both directions.
//
TrkErrCode
KalRep::buildTraj(){
  TrkErrCode retval; // preset to success
// set its parameters according the inner site parameters
  KalSite* firstsite = _sites[_hitrange[0]];
// clone the seed traj: cast it back to a simptraj
  TrkSimpTraj* straj = _seedtraj->clone();
  if(!firstsite->setTrajState(trkIn,straj)){
    delete straj;
    return TrkErrCode(TrkErrCode::fail,KalCodes::matrix,"Matrix inversion error");
  }
// correct fitrange for site's local range
  double range[2];
  double sitedist = firstsite->localLength() - firstsite->globalLength();
  range[0] = _fitrange[0] + sitedist;
  range[1] = _fitrange[1] + sitedist;
  straj->setFlightRange(range);
// build the new piece traj starting with this piece
  TrkDifPieceTraj* newtraj = 
    new TrkDifPieceTraj(straj,_fitrange[0],_fitrange[1]);
  _maxdist = 0.0;
  _maxfltdif = 0.0;
  double oldlen = firstsite->globalLength();
  unsigned isite = _hitrange[0];
  while(true){
// find the first active non-hit site that's far enough away from the last one
    do
      isite++;
    while(isite<_hitrange[1] && ( _sites[isite]->nDof()> 0 || 
				  (!_sites[isite]->isActive()) ||
				  (_sites[isite]->globalLength()-oldlen<_kalcon.minGap())));
// test for loop termination
    if(isite>=_hitrange[1])break;
    KalSite* thesite = _sites[isite];
//  Compute the (transverse) distance between new and old trajectories at 
//  this site for testing convergence
    double sitelen = thesite->globalLength() ;
    //  offset to avoid a mismatch at trajectory discontinuities
    double testlen = sitelen - 0.1;
    HepPoint oldpoint = _reftraj->position(testlen);
    HepPoint newpoint = newtraj->position(testlen);
    Hep3Vector dire = newtraj->direction(testlen);
    Hep3Vector diff = oldpoint - newpoint;
    double dist = diff.perp(dire);
    if(dist > _maxdist)_maxdist = dist;    
//  Smoothed result is the statistical combination of the inner site's outward
//  weight with the outer site's inward weight
    KalParams smoothed;
    thesite->mergeParams(_sites[isite+1],smoothed);
    if(smoothed.matrixOK()){
// build a piece out of the smoothed result.
      TrkSimpTraj* straj = _seedtraj->clone();
      *(straj->parameters()) = smoothed.trackParameters();
// set the range on this trajectory
      double fltrng[2];
      fltrng[0] = thesite->localLength();
      fltrng[1] = fltrng[0] + std::max(_fitrange[1] - sitelen,_kalcon.minFltLen());
      straj->setFlightRange(fltrng);
//  append this to the piecetraj
      double gap(0);
      const TrkErrCode& appenderror = newtraj->append(sitelen,straj,gap);
      if(appenderror.success()){
        thesite->setGap(gap);
      } else {
// cleanup and exit
        delete straj;
        delete newtraj;
        return appenderror;
      }
    } else {
      ErrMsg(error) << "Failed to construct merged traj segment" << endmsg;
    }
// mark the site's position as the starting point for the next segment
    oldlen = thesite->globalLength();
  }
// test the new trajectory before making it official
  double hflt = _sites[_hitrange[0]]->globalLength();
  double dot = newtraj->direction(hflt).dot(_ptraj->direction(hflt));
  if( dot > _kalcon.minDot() || _niter <= 1){
//  Delete the old trajectory, and set the new one
    if(_ptraj != _reftraj)
      delete _ptraj;
    _ptraj = newtraj;
  } else {
    delete newtraj;
    retval = TrkErrCode(TrkErrCode::fail,KalCodes::badtraj,"Unreasonable trajectory");
  }
  return retval;
}
//  Process the sites in a given direction between the specified limits
bool
KalRep::process(KalSite* psite,int startindex,int endindex,trkDirection fdir) {
//
  int step = fdir==trkOut ? 1 : -1;
  int nstep = fdir==trkOut ? endindex-startindex+1 : startindex-endindex+1;
  const KalSite* prevsite = psite;
//  Filter through the sites and process them for this direction. The starting
//  parameters are the inflated-error intial parameters.  Only process sites
//  which need it.  Note that the processing latches to force processing any
//  sites appearing after a site which needed processing.
  bool status(true);
  bool needsfit(false);
  _chisq[fdir] = -1.0; // previous chisq no longer valid
  int iindex(startindex);
  while( nstep>0) {
    KalSite* thesite = _sites[iindex];
    needsfit = thesite->needsFit(fdir) || needsfit;
    if(needsfit){
      if(!thesite->process(prevsite,fdir)){
	status = false;
	break;
      }
    }
    prevsite = thesite;
    iindex += step;
    nstep--;
  }
  if(status)setFit(fdir);
  return status;
}
//
//  Hot functions
//
KalHit*
KalRep::findHotSite(const TrkHitOnTrk* thehot,unsigned& siteindex) const {
//  Find the hot; exhaustive search for now
  KalHit* returnhot(0);
  for(unsigned isite=0;isite<_sites.size();isite++)
    if(_sites[isite]->kalHit() != 0) {
      if(_sites[isite]->kalHit()->hitOnTrack() == thehot){
	returnhot = _sites[isite]->kalHit();
	siteindex = isite;
	break;
      }
    }
  return returnhot;
}

const KalHit*
KalRep::findHotSite(const TrkHitOnTrk* thehot) const {
  KalHit* returnhot(0);
  for(unsigned isite=0;isite<_sites.size();isite++)
    if(_sites[isite]->kalHit() != 0) {
      if(_sites[isite]->kalHit()->hitOnTrack() == thehot){
	returnhot = _sites[isite]->kalHit();
	break;
      }
    }
  return returnhot;
}

//
void
KalRep::activateHot(TrkHitOnTrk* thehot) {
// action only necessary if the hot isn't already active
  if( (! thehot->isActive()) && thehot->isUsable()){
    TrkRep::activateHot(thehot);
// Require the HOT be useable and associated with a KalHit site
    unsigned index;
    KalHit* thesite = findHotSite(thehot,index);
    if(thesite != 0 ){
      resetFit(); // force refitting
      thesite->setActivity(true);
      findHitSites();
    }
  }
}
//
void
KalRep::deactivateHot(TrkHitOnTrk* thehot) {
  if(thehot->isActive()){
    TrkRep::deactivateHot(thehot);
    unsigned index;
    KalHit* thesite = findHotSite(thehot,index);
    if(thesite != 0){
      resetFit(); // force refitting
      thesite->setActivity(false);
      findHitSites();
    }
  }
}
//
void
KalRep::updateHot(TrkHitOnTrk* thehot) {
// Require the HOT be useable and associated with a KalHit site
  unsigned index;
  KalHit* thesite = findHotSite(thehot,index);
  if(thesite != 0 && thehot->isUsable()){
    resetFit(); // force refitting
    thesite->update(_reftraj,_refmom); // note, hit sites don't care about momentum,
  }
}
//
void
KalRep::addHot(TrkHitOnTrk* thehot) {
// try to find the HOT in the existing track
  unsigned index;
  KalHit* thesite = findHotSite(thehot,index);
  if(thesite == 0){
    TrkRep::addHot(thehot);
// create a new KalHit from this hot
    thesite = new KalHit(_reftraj,thehot);
    _sites.push_back(thesite);
// refind the hit sites
    findHitSites();
// update hit site limits
    if(thehot->isActive()){
      resetFit(); // force refitting
    } else {
// find the hot we just inserted, so we have its index
      KalHit* thesite = findHotSite(thehot,index);
      assert(thesite != 0);
// if there's a valid outward fit, process the inactive hot.  This doesn't change
// the fit, but it makes sure the sitelist remains internally self-consistent (no gaps)
      if(hasFit(trkOut)){
	if(index > 0)
	  thesite->process(_sites[index-1],trkOut);
	else {
// have to use end site to handle an inactive hot on the end of a rep
	  thesite->process(&_endsites[trkIn],trkOut);
	}
      }
// same for inwards fit
      if(hasFit(trkIn)){
	if(index< _sites.size()-1)
	  thesite->process(_sites[index+1],trkIn);
	else {
	  thesite->process(&_endsites[trkOut],trkIn);
	}
      }
    }
// add the hot to the rep.  This will reset the current flag as needed
//    TrkRep::addHot(thehot);
  } else
    ErrMsg(error) << "cannot add HOT already on this rep" << endmsg;
}
//
void
KalRep::removeHot(TrkHitOnTrk* thehot) {
// try to find the HOT in the existing track
  unsigned index;
  KalHit* thesite = findHotSite(thehot,index);
  if(thesite != 0){
// Go to the sites on either side of the one to be removed and invalidate them
    if(index > 0)_sites[index-1]->invalidateSite(trkOut);
    if(index < _sites.size()-2)_sites[index+1]->invalidateSite(trkIn);
// remove and delete the site; note that the HOT must be deleted elsewhere!
    KalSite* rmsite = _sites[index];
    _sites.erase(_sites.begin()+index);
    resetFit(); // force refitting
    if ( rmsite != 0 ){
      delete thesite;
// update hit sites
      findHitSites();
// if the hot wasn't active, it doesn't change the fit
      if(thehot->isActive())
	resetFit();
    } else
      ErrMsg(error) << "Error removing HOT" << endmsg;
  } else
    ErrMsg(error) << "Error, requested HOT not found" << endmsg;
// remove the hot from the rep hotlist, no matter what
  TrkRep::removeHot(thehot);
}
// add intersection to the track
void
KalRep::addInter(DetIntersection const& detinter){
// use the reference  momentum; should be site-specific, FIXME!!!!
  double sitemom = _refmom;
//  Build the KalMaterial site on this new piece
  KalMaterial* newmat = new KalMaterial(detinter,_reftraj,sitemom,particleType());
  assert(newmat != 0);
// insert this into the sites and update
  _sites.push_back(newmat);
// refind the hit sites; this also reorders material sites
  findHitSites();
  resetFit(); // force refitting
}
//
//  Extend the track
//
TrkErrCode
KalRep::extendThrough(const TrkVolume& tvol, trkDirection trkdir){
  TrkErrCode retval(TrkErrCode::succeed,KalCodes::alreadyextended,"KalRep is already extended");
  if(fitValid() && _ptraj != 0){
//  Extend the trajectory through the volume
    double newf(0.0); // initial value here is unimportant, extension uses the ptraj range
    bool isextended = tvol.extendThrough( _ptraj, newf, trkdir );
    if(isextended)
      retval = extendThrough(newf);
  } else
    retval = TrkErrCode(TrkErrCode::fail,KalCodes::notready,
			"KalRep not ready for extending");
  return retval;
}

TrkErrCode
KalRep::extendThrough(double newf) {
  TrkErrCode retval(TrkErrCode::succeed,KalCodes::alreadyextended,"KalRep is already extended");
// test the obvious
  if(newf < startValidRange() || newf > endValidRange()){
// see if we're extending inwards or outwards
    trkDirection trkdir = newf<startValidRange()? trkIn : trkOut;
// basic tests
    if( _extendable[trkdir] && fitValid() && _ptraj != 0 ){
      double piecerange[2];
      KalSite* endsite;
      switch(trkdir){
      case trkOut: default:
	piecerange[0] = _ptraj->hiRange();
	piecerange[1] = newf;
	endsite = _sites.back();
	break;
      case trkIn:
	piecerange[0] = newf;
	piecerange[1] = _ptraj->lowRange();
	endsite = _sites.front();
	break;
      }
//  Don't bother doing anything if the traj wasn't really extended
      if(piecerange[1] > piecerange[0]){
	if(_kalcon.materialSites()){
//  Get the top-level detector set
// For now, use an empty detector set.  In future, this should come from the geometry service DNB_RKK
//    static const DetSet* trkmodel = new DetSet("dummy",1);
//	  const DetSet* trkmodel= &(DchDetector::GetInstance()->dchSet());
//  Find the material intersections, and create KalMaterial sites for them
	  std::vector<DetIntersection> tlist; tlist.reserve(16);
//	  trkmodel->intersection(tlist,_reftraj,piecerange);
	  if (_kalcon.getDetModel()!=0) {
	    _kalcon.getDetModel()->intersection(tlist,_reftraj,piecerange);
	  }

// get the momentum at the current end of the track
	  double loclen;
	  double extendmom;
	  const TrkSimpTraj* reftraj = 
	    trkdir == trkOut ? localTrajectory(_sites.back()->globalLength(),loclen) :
	    localTrajectory(_sites.front()->globalLength(),loclen);
	  if(reftraj != 0){
// get the momentum from this using the mom calculator
	    Hep3Vector momvec = TrkMomCalculator::vecMom(*reftraj,_kalcon.bField(),loclen);
	    extendmom = momvec.mag();
	  } else
	    return TrkErrCode(TrkErrCode::fail,KalCodes::momentum,
			      "Track momentum cannot be calculated");
// Add the material sites
	  KalMaterial* stopsite = trkdir == trkOut ?
	    createMaterialSites( tlist,_sites,0,tlist.size()-1,extendmom,trkdir) :
	    createMaterialSites( tlist,_sites,tlist.size()-1,0,extendmom,trkdir);
	  if(stopsite != 0 && trkdir == trkOut){
// limit the track in this direction
	    _stopsite = stopsite;
	    _fitrange[1] = std::min(_fitrange[1],stopsite->globalLength());
	    _extendable[trkOut] = false;
	    retval = TrkErrCode(TrkErrCode::fail,KalCodes::cannotfullyextend,"cannot fully extend");
	  }
	}
// Add the bend sites
	if(_kalcon.bendSites()){
	  createBendSites(piecerange,_sites);
	}
// refind the hit sites (in case they've moved)
	findHitSites();
//  process the new sites, and extend the trajectory for them (if the fit is current)
	if(fitCurrent()){
	  unsigned isite=0;
	  for(;isite<_sites.size();isite++){
	    if(endsite == _sites[isite]) break;
	  }
	  if (isite < _sites.size() ){
	    retval = extendSites(isite,trkdir);  // endsite still exists
	  } else
	    retval = TrkErrCode(TrkErrCode::fail,KalCodes::endextend,
				"Can't find endsite");
	}
// if the extension was successfull, update the fit range
// update fitrange
	if(retval.success()){
	  _fitrange[0] = std::min(_fitrange[0],newf-_kalcon.fltEpsilon());
	  _fitrange[1] = std::max(_fitrange[1],newf+_kalcon.fltEpsilon());
// update the range of the trajectory
	  _ptraj->setFlightRange(_fitrange);
	  _reftraj->setFlightRange(_fitrange);
	}
      }
    } else {
      if(!_extendable[trkdir] )
	retval =  trkdir==trkOut ?
	  TrkErrCode(TrkErrCode::fail,KalCodes::extendable,
		     "Track is not extendable in direction trkOut") : 
	  TrkErrCode(TrkErrCode::fail,KalCodes::extendable,
		     "Track is not extendable in direction trkIn");
      else if( (!fitValid()) ||  _ptraj == 0)
	retval = TrkErrCode(TrkErrCode::fail,KalCodes::notready,
			    "KalRep not ready for extending");
    }
  }
  return retval;
}



//
//  Update all the sites after building the trajectory.  This is used to iterate
//  the fit.
TrkErrCode
KalRep::updateSites() {
//  The output trajectory has to have built
  if(_ptraj != 0 && _ptraj != _reftraj ){
//  Switch around the trajs
    delete _reftraj;
    _reftraj = _ptraj;
//  Set flags: the fit may still be valid, but is not current
    setCurrent(false);
// update the refence momentum
    updateRefMom();
// update the sites
// split the sites at the reference momentum flight length
    int losite = findNearestSite(_refmomfltlen);
// reset flightlength change
    _maxfltdif = 0.0;
    updateSites(losite,0,_refmom,trkIn);
    updateSites(losite+1,_sites.size()-1,_refmom,trkOut);
// update the end sites
    updateEndSites(_kalcon.smearFactor());
// if the flightlength of the hit sites has changed a lot, this means the
// intersection and field integrals may no longer be correct, so rebuild the
// material and bend sites
    if( _kalcon.materialSites() && _maxfltdif > _kalcon.intersectionTolerance() && 
	_ninter < _kalcon.maxIntersections())
      reIntersect();
    return TrkErrCode(TrkErrCode::succeed);
  } else{
    if(_ptraj == 0)
      return TrkErrCode(TrkErrCode::fail,KalCodes::notready,
		      "KalRep not ready for updating");
		else
      return TrkErrCode(TrkErrCode::succeed,KalCodes::valid,
		      "KalRep already updated");
  }
}
//
//  fit existing sites in both directions
//
TrkErrCode
KalRep::fitSites() {
// perform any fixups necessary for the hots moving around
  fixupSites();
// check that the track is fitable
  TrkErrCode retval;
  if(isFitable(retval)){
// process all the sites between the first and last hit
    bool status = process(&_endsites[trkIn],_hitrange[0],_hitrange[1],trkOut);
    if(status)
      status = process(&_endsites[trkOut],_hitrange[1],_hitrange[0],trkIn);
    if(!status)
      retval =  TrkErrCode(TrkErrCode::fail,KalCodes::processing,
			   "Processing Failed");
  }
  return retval;
}

void
KalRep::setFitRange(const Trajectory* traj) {
//  If initial volumes have been specified, use them to extend these limits
  double range;
  const TrkVolume* invol = _kalcon.trkVolume(trkIn);
  const TrkVolume* outvol = _kalcon.trkVolume(trkOut);
  if( invol != 0 && invol->extendThrough( traj,range,trkIn))
    _fitrange[0] = std::min(_fitrange[0],range);
  if(outvol != 0 && outvol->extendThrough( traj,range,trkOut))
    _fitrange[1] = std::max(_fitrange[1],range);
// use the seed to start the reference piece traj (initialy 1 piece)
// require a non-zero flight range
  if( (_fitrange[1] - _fitrange[0]) < _kalcon.minFltLen() ){
// this is really bad: the track has <=1 hit and _no_ extension volumes.  This
// track will be unfittable no matter what, but I don't want it crash.  Thus
// we create a phony range.
    double midrange = (_fitrange[0] + _fitrange[1] )/2.0;
    _fitrange[0] = midrange - _kalcon.minFltLen()/2.0;
    _fitrange[1] = midrange + _kalcon.minFltLen()/2.0;
  }
}

//
// build the reference trajectory from the seed
//
void
KalRep::buildRefTraj() {
  setFitRange(_seedtraj);
// reset the seed trajectory range
  _seedtraj->setFlightRange(_fitrange);
  _reftraj = new TrkDifPieceTraj(*_seedtraj,_fitrange[0],_fitrange[1]);
  assert(_reftraj != 0);
}
// find first and last 'hit' sites.  We really want all the sites
// which add information to the fit.
void
KalRep::findHitSites() {
// sort the vector
  std::sort(_sites.begin(),_sites.end(),babar::Collection::PtrLess());
// find the first and last hot sites
  _hitrange[0] = _sites.size()+1;
  _hitrange[1] = -1;
  for(unsigned isite=0;isite<_sites.size();isite++)
    if(_sites[isite]->nDof() > 0) {
      _hitrange[0] = isite;
      break;
    }
  for(int isite=_sites.size()-1;isite>=0;isite--)
    if(_sites[isite]->nDof() > 0) {
      _hitrange[1] = isite;
      break;
    }
}
//
//  Process the sites after/before the hits
//
TrkErrCode
KalRep::extendSites(int startsite,trkDirection tdir) {
  TrkErrCode retval;
// compute some index limits
  int endsite,nsites,sitestep;
  switch(tdir){
  case trkIn:
    endsite = 0;
    sitestep = -1;
    nsites = startsite-endsite;
    break;
  case trkOut:
    endsite = _sites.size()-1;
    sitestep = 1;
    nsites = endsite-startsite;
    break;
  }
// process the sites
  if(nsites>0){
    bool status(false);
    KalSite* psite(_sites[startsite]);
// make sure to start processing the site _after_ the first one
    status = process(psite,startsite+sitestep,endsite,tdir);
// extend the trajectory
    if(status)
      retval = extendTraj(startsite+sitestep,tdir);
    else
      retval =  TrkErrCode(TrkErrCode::fail,KalCodes::processing,
			   "Processing Failed");
  }
  return retval;
}
// create KalMaterial sites.
void
KalRep::buildMaterialSites(double range[2],std::vector<DetIntersection>& tlist) {
// avoid a detmodel error
  if(range[1]>range[0]){
//  Get the Tracking DetectorModel tree from TrkEnv
// For now, use an empty detector set.  In future, this should come from the geometry service DNB_RKK
//    static const DetSet* trkmodel = new DetSet("dummy",1);
//    const DetSet* trkmodel= &(DchDetector::GetInstance()->dchSet());
// find the intersections with the reference trajectory, if no intersections provided
    if (ErrLogging(debugging)){
      if(tlist.size() == 0)
        std::cout << "No detector " << std::endl;
    }
//      trkmodel->intersection(tlist,_reftraj,range);
//    std::cout<<"Intersections n="<<tlist.size()<<std::endl;
    if (_kalcon.getDetModel()!=0) {
      _kalcon.getDetModel()->intersection(tlist,_reftraj,range);
      if (ErrLogging(debugging)){ std::cout<<"Intersections n="<<tlist.size()<<std::endl; }
    }
// split the intersection list into those before and after the reference momentum
    unsigned below(0);
    while(below<tlist.size() && tlist[below].pathlen<_refmomfltlen)
      below++;
// build material sites on either side of the first inwards
    if(below>0)
      createMaterialSites(tlist,_sites,below-1,0,_refmom,trkIn);
    _stopsite = 
      createMaterialSites(tlist,_sites,below,tlist.size()-1,_refmom,trkOut);
// if the track stops before the end of this range (going outwards), limit the
// rep and cleanup
    if(_stopsite != 0){
// limit the track in this direction
      _fitrange[1] = std::min(_fitrange[1],_stopsite->globalLength());
      _extendable[trkOut] = false;
    }
  }
}
// same thing for hits
void
KalRep::buildHitSites() {
// intialize some HOT stuff
//  Loop over the hots, and create KalHit sites for them, and add them to the site list
  TrkHotList* hots = hotList();
  TrkHotList::nc_hot_iterator end = hots->end();
  for(TrkHotList::nc_hot_iterator ihit=hots->begin();ihit!=end;ihit++){
//  Only 'useable' hots should be made into sites
    if(ihit->isUsable()){
      KalHit* newhit = new KalHit(_reftraj,ihit.get());
      assert(newhit != 0);
      _sites.push_back(newhit);
    }
  }
// Locate the first and last hit site
  findHitSites();
// see if any hots are active
// sometimes the last hit will move outside the fit range.  This is not really an
// error, so just adjust the fit range appropriately
  if( _hitrange[0]<_sites.size() && _hitrange[1]>=0){
    if(_extendable[trkIn] && _fitrange[0] > _sites[_hitrange[0]]->globalLength())
      _fitrange[0] = _sites[_hitrange[0]]->globalLength() - _kalcon.fltEpsilon();
    if(_extendable[trkOut] && _fitrange[1]< _sites[_hitrange[1]]->globalLength())
      _fitrange[1] = _sites[_hitrange[1]]->globalLength() + _kalcon.fltEpsilon();
  }
}
//
// same thing for bend sites
void
KalRep::buildBendSites(double frange[2]) {
// create the bend sites all at once
  createBendSites(frange,_sites);
  // refind the hit sites
  findHitSites();
}

void
KalRep::createBendSites(double range[2], std::vector<KalSite*>& sites) const {
// divide the trajectory range; split it at 0.0 if possible to minimize end
// effects at the origin
  BFieldIntRange brange(range[0],range[1]);
  if(range[0]>range[1]){
    ErrMsg(warning) << "Reversed integration range, inverting" <<  endmsg;
    brange.invert();
  }
  std::vector<BFieldIntRange> ranges;
  _kalcon.bFieldIntegrator().divideRange(_reftraj,brange,ranges);
// loop over the divisions, creating the bend sites.
  for(size_t irange=0;irange<ranges.size();++irange){
// reference momentum is good enough for bend corrections
    KalBend* bendsite = new KalBend(_kalcon.bFieldIntegrator(),_reftraj,ranges[irange],_refmom,_charge);
    assert(bendsite != 0);
    sites.push_back(bendsite);
  }
}

//----------------------------------------------------------------------
Hep3Vector 
KalRep::momentum(double fltL) const {
//----------------------------------------------------------------------
//  const BField& theField = _kalcon.bField();
  // kludge DNB_RKK
//  static const BField* theField = new BFieldFixed(0.0,0.0,1.0);
  double localFlt = 0.;
  const TrkSimpTraj* locTraj = localTrajectory(fltL,localFlt);
  return TrkMomCalculator::vecMom(*locTraj, _kalcon.bField(), localFlt);
}

//----------------------------------------------------------------------
double 
KalRep::pt(double fltL) const {
//----------------------------------------------------------------------
  Hep3Vector mom = momentum(fltL);
  return mom.perp();
}

//----------------------------------------------------------------------
BbrVectorErr 
KalRep::momentumErr(double fltL) const {
//----------------------------------------------------------------------
//  const BField& theField = _kalcon.bField();
  // kludge DNB_RKK
  
  double localFlt = 0.;
  const TrkSimpTraj* locTraj = localTrajectory(fltL,localFlt);
  return TrkMomCalculator::errMom(*locTraj, _kalcon.bField(), localFlt);
}


void
KalRep::updateHots() {
//  update the KalHit sites
  for(unsigned isite=0;isite<_sites.size();isite++){
    KalSite* thesite = _sites[isite];
    if(thesite->kalHit() != 0)
      thesite->update(_reftraj,_refmom); // note, hit sites don't care about momentum,
// so the seed value is OK here
  }
// perform any fixups necessary for the hots moving around
  fixupSites();
// reset the fit and make sure everyone knows about the change
  resetFit();
}

double
KalRep::parameterDifference(trkDirection tdir) const {
// compare the parameters at the specified end of the track,
// computing a 'chisquared' as if the current fit results and
// the reference were independent (which they aren't). This
// is used for convergence testing.  Take the first/last hit site as reference

// protect against invalid fits
  if(_hitrange[0]>_hitrange[1])
    return 1.0e8;
  double fltlen(0.0);
  switch(tdir){
  case trkIn:
    fltlen = _sites[_hitrange[0]]->globalLength();
    break;
  case trkOut:
    fltlen = _sites[_hitrange[1]]->globalLength();
  }
  double loclen;
// find the parameters at the specified flight length.
  KalParams params(*_ptraj->localTrajectory(fltlen,loclen)->parameters());
// subtract the previous iteration's parameters (= reference trajectory)
  KalParams ref_params(*_reftraj->localTrajectory(fltlen,loclen)->parameters());
  return params.chisq(ref_params);
}

//------------------------------------------------------------------------
HepMatrix 
KalRep::posmomCov(double fltL) const {
//------------------------------------------------------------------------
  const BField& theField = _kalcon.bField();
  double localFlt = 0.;
  const TrkSimpTraj* locTraj = localTrajectory(fltL,localFlt);
  return TrkMomCalculator::posmomCov(*locTraj, theField, localFlt);
}

//------------------------------------------------------------------------
void 
KalRep::getAllCovs(double fltL,
			HepSymMatrix& xxCov,
			HepSymMatrix& ppCov,
			HepMatrix&    xpCov)      const {
//------------------------------------------------------------------------
  const BField& theField = _kalcon.bField();
  double localFlt = 0.;
  const TrkSimpTraj* locTraj = localTrajectory(fltL,localFlt);
  TrkMomCalculator::getAllCovs(*locTraj, theField, localFlt,
			      xxCov,ppCov,xpCov); 
}

//------------------------------------------------------------------------
void 
KalRep::getAllWeights(double fltL,
				 HepVector& pos,
				 HepVector& mom,
				 HepSymMatrix& xxWeight,
				 HepSymMatrix& ppWeight,
				 HepMatrix&    xpWeight) const {
//------------------------------------------------------------------------
  const BField& theField = _kalcon.bField();
  double localFlt = 0.;
  const TrkSimpTraj* locTraj = localTrajectory(fltL,localFlt);
  TrkMomCalculator::getAllWeights(*locTraj, theField, localFlt,
			      pos,mom,xxWeight,ppWeight,xpWeight); 

}

TrkErrCode
KalRep::append(KalStub& stub) {
// make sure this stub belongs to this rep!!!
  if(stub.sameRep(*this) ){
// set fit not current
    setCurrent(false);
// reset iterations (if necessary)
    _niter = 0;
// Append or prepend the stub sites, as appropriate
    std::deque<KalSite*>& stubsites = stub.sites();
    for(std::deque<KalSite*>::iterator siter = stubsites.begin();
	siter!= stubsites.end();siter++){
// NB: we should be using the 'insert' command here to put all the stub sites
// in at the same time, however the Sun compiler doesn't currently support this,
// so the sites are added one at a time (out of order).
      _sites.push_back(*siter);
// we also have to append the HOTs to the hot list
      if((*siter)->kalHit() != 0)
	TrkRep::addHot((*siter)->kalHit()->hitOnTrack());
    }
// sort the hotlist
    hotList()->sort();
// now remove the sites from the stub (so they won't be deleted when it is)
    stubsites.clear();
// extend the track range
    if(stub.direction() == trkOut)
      _fitrange[1] = std::max(_fitrange[1],stub.hitRange(trkOut)+_kalcon.fltEpsilon());
    else
      _fitrange[0] = std::min(_fitrange[0],stub.hitRange(trkIn)-_kalcon.fltEpsilon());
// find the new hits
    findHitSites();
// User is responsable for the refit
    return TrkErrCode();
  } else
    return TrkErrCode(TrkErrCode::fail,KalCodes::stubmatch,
		      "KalRep error: KalStub doesn't match KalRep or is not a primary stub!");
}

double
KalRep::radiationFraction(double* range) const {
  double radfrac(0.0);
  unsigned nsites = _sites.size();
  KalSite* site;
  double minlen(_fitrange[0]);
  double maxlen(_fitrange[1]);
  if(range != 0){
    minlen = range[0];
    maxlen = range[1];
  }
  for(unsigned isite= 0;isite<nsites;isite++){
    site = _sites[isite];
    if(site->kalMaterial() != 0 && 
       site->globalLength() >= minlen && 
       site->globalLength() <= maxlen)
      radfrac += site->kalMaterial()->radiationFraction();
  }
  return radfrac;
}
//
double
KalRep::energyLoss(double* range) const {
  double eloss(0.0);
  unsigned nsites = _sites.size();
  const KalSite* site;
  double minlen(_fitrange[0]);
  double maxlen(_fitrange[1]);
  if(range != 0){
    minlen = range[0];
    maxlen = range[1];
  }
  for(unsigned isite= 0;isite<nsites;isite++){
    site = _sites[isite];
    if(site->kalMaterial() != 0 && 
       site->globalLength() >= minlen && 
       site->globalLength() <= maxlen)
      eloss += site->kalMaterial()->energyChange(trkOut);
  }
  return eloss;
}

double
KalRep::thetaBend(double* range) const {
  double bend(0.0);
  unsigned nsites = _sites.size();
  KalSite* site;
  double minlen(_fitrange[0]);
  double maxlen(_fitrange[1]);
  if(range != 0){
    minlen = range[0];
    maxlen = range[1];
  }
  for(unsigned isite= 0;isite<nsites;isite++){
    site = _sites[isite];
    if(site->kalBend() != 0 &&
       site->globalLength() >= minlen && 
       site->globalLength() <= maxlen)
      bend += site->kalBend()->deltaTheta();
  }
  return bend;
}

double
KalRep::phiBend(double* range) const {
  double bend(0.0);
  unsigned nsites = _sites.size();
  KalSite* site;
  double minlen(_fitrange[0]);
  double maxlen(_fitrange[1]);
  if(range != 0){
    minlen = range[0];
    maxlen = range[1];
  }
  for(unsigned isite= 0;isite<nsites;isite++){
    site = _sites[isite];
    if(site->kalBend() != 0 &&
       site->globalLength() >= minlen && 
       site->globalLength() <= maxlen)
      bend += site->kalBend()->deltaPhi();
  }
  return bend;
}

const TrkSimpTraj*
KalRep::localTrajectory(double fltlen,double& loclen) const {
  return pieceTraj().localTrajectory(fltlen, loclen);
}


const TrkDifPieceTraj&
KalRep::pieceTraj() const
{ return *_ptraj; }

TrkDifPieceTraj& 
KalRep::pieceTraj() 
{ return *_ptraj; }

bool
KalRep::converged() const {
// convergence requres both directions to have been processed
  bool converged = hasFit(trkIn) && hasFit(trkOut);
// only reps with hits benifit from iteration
  if(hotList()->hitCapable()){
// several convergence tests, do them roughly in order of importance.
// trajectory convergence  
    if(converged) converged = _maxdist < _kalcon.distanceTolerance();
// momentum convergence
    if(converged) converged = fabs(estimatedMomDiff()) < _kalcon.maxMomDiff();
//  Parameter convergence; if the parameter differences tolerance is set < 0.0  skip the test
    if(converged) converged = _kalcon.maxParamDiff(trkOut) < 0.0 ||
		    parameterDifference(trkOut) < _kalcon.maxParamDiff(trkOut);
    if(converged) converged = _kalcon.maxParamDiff(trkIn) < 0.0 ||
		    parameterDifference(trkIn) < _kalcon.maxParamDiff(trkIn);
  }
  if(_niter<=1) converged=false;
  //std::cout<<"converged "<<converged<<" "<<_niter<<std::endl;

  return converged;
}

  
int
KalRep::findNearestSite(double fltlen) const {
// binary search in index space, testing against flight length
  int retsite(-1);
  int maxsite = _sites.size()-1;
  if(maxsite > 0 ){
    if(fltlen >= _sites[0]->globalLength() &&
       fltlen < _sites[maxsite]->globalLength() ) {
      unsigned step = std::max(1,(int)maxsite/4);
      int losite = (maxsite+1)/2;
      while(true){
	if(fltlen < _sites[losite]->globalLength())
	  losite -= step;
	else if(losite<maxsite && 
		fltlen> _sites[losite+1]->globalLength())
	  losite += step;
	else
	  break;
	step = std::max(1,(int)std::min(step/2,step-1));
      }
      retsite = losite;
    } else if( fltlen >= _sites[maxsite]->globalLength() ) {
      retsite = maxsite;
    }
  }
  return retsite;
}

void
KalRep::findBoundingSites(double fltlen,
			  const KalSite*& insite,
			  const KalSite*& outsite) const {
// preset to failure
  insite = outsite = 0;
  int losite = findNearestSite(fltlen);
  if(losite >= 0 && losite < _sites.size()-1){
    insite = _sites[losite];
    outsite = _sites[losite+1];
  }
}

bool
KalRep::filterTrajs(double fltlen,
		    TrkSimpTraj* intraj,
		    TrkSimpTraj* outtraj) const {
  bool retval(false);
// find the bounding sites
  const KalSite* insite(0);
  const KalSite* outsite(0);
  findBoundingSites(fltlen,insite,outsite);
  if(insite != 0 && outsite != 0){
    retval = insite->setTrajState(trkOut,intraj) &&
      outsite->setTrajState(trkIn,outtraj);
// must set ranges too; damn this stupid double[2] interface
    double newrange[2];
    newrange[0] = _fitrange[0];
    newrange[1] = fltlen;
    intraj->setFlightRange(newrange);
    newrange[1] = _fitrange[1];
    newrange[0] = fltlen;
    outtraj->setFlightRange(newrange);
  }
  return retval;
}

bool
KalRep::filterTraj(double fltlen,
		   trkDirection tdir,
		   TrkSimpTraj* traj) const {
  bool retval(false);
// find the bounding sites
  const KalSite* insite(0);
  const KalSite* outsite(0);
  findBoundingSites(fltlen,insite,outsite);
  double newrange[2];
  if(tdir == trkIn){
// going inwards we want the inward fit of the _outer_ site
    retval = outsite != 0 &&
      outsite->setTrajState(trkIn,traj);
    if(retval){
      newrange[1] = _fitrange[1];
      newrange[0] = fltlen;
      traj->setFlightRange(newrange);
    }
  } else {
    retval = insite != 0 &&
      insite->setTrajState(trkOut,traj);
    if(retval){
      newrange[0] = _fitrange[0];
      newrange[1] = fltlen;
      traj->setFlightRange(newrange);
    }
  }
  return retval;
}

bool
KalRep::smoothedTraj(const KalHit* hitsite,TrkSimpTraj* traj) const {
  // find this site
  std::vector<KalSite*>::const_iterator ifnd = std::find(_sites.begin(),_sites.end(),hitsite);
  if(ifnd != _sites.end()) {
   return smoothedTraj(ifnd,ifnd,traj);
  } else
    return false;
}

bool
KalRep::smoothedTraj(std::vector<KalSite*>::const_iterator ifirst,
    std::vector<KalSite*>::const_iterator isecond, TrkSimpTraj* traj) const {
  bool retval(false);
    // find the sites on either side: if none, mark them as the end of the container
    std::vector<KalSite*>::const_iterator iprev = ifirst;
    if(iprev != _sites.begin())
      --iprev;
    else
      iprev = _sites.end();
    std::vector<KalSite*>::const_iterator inext = isecond;
    ++inext;
   if(iprev != _sites.end() && inext != _sites.end()){
    // both bounding sites exist; merge their parameters
    KalParams smoothed;
    (*iprev)->mergeParams(*inext,smoothed);
    if(smoothed.matrixOK()){
      retval = true;      
      *(traj->parameters()) = smoothed.trackParameters();
    }
  } else if(iprev != _sites.end()){
    retval = (*iprev)->setTrajState(trkOut,traj);      
  } else if(inext != _sites.end()) {
    retval = (*inext)->setTrajState(trkIn,traj);
  }
  return retval;
}

void
KalRep::updateRefMom() {
// get the local trajectory at the reference flight length
  double loclen;
  const TrkSimpTraj* reftraj = localTrajectory(_refmomfltlen,loclen);
  if(reftraj != 0){
// get the momentum from this using the mom calculator
    Hep3Vector momvec = TrkMomCalculator::vecMom(*reftraj,_kalcon.bField(),loclen);
    double delmom = momvec.mag()-_refmom;
    double momfac = 1.0;
    double floor = 0.5 ;
// damp momentum update after 3rd iteration
      if(_niter>3) {
         double arguement = (_niter-3) * _kalcon.momUpdateFactor();
         momfac = floor + (1-floor)*exp(-arguement);
      }
   _refmom += momfac*delmom;
  }
}

KalMaterial *
KalRep::createMaterialSites( std::vector<DetIntersection> tlist,
			     std::vector<KalSite*>& sites,
			     int startindex,int endindex,
			     double initialmom,trkDirection dedxdir) const {
//  Loop over the intersections creating the KalMaterial sites
  double sitemom(initialmom);
  int step = dedxdir==trkOut ? 1 : -1;
  int nstep = dedxdir==trkOut ? endindex-startindex+1 : startindex-endindex+1;
  int iindex(startindex);
  KalMaterial* retval(0);
  while(nstep>0) {
// ignore intersections without a material specified: these are placeholders
// here we have to deal with a DetectorModel design flaw, that the material is returned by reference.
		if( &(tlist[iindex].delem->material(tlist[iindex])) != 0){
//  Build the KalMaterial site on this new piece
			KalMaterial* newmat = new KalMaterial(tlist[iindex],_reftraj,
				sitemom,particleType());
			assert(newmat != 0);
// if going outward, test momentum loss
// check if too much momentum has been lost: if so, stop the sitelist
// and prevent the track from being extended
			if(retval==0 && dedxdir == trkOut && stopsIn(newmat))
				retval = newmat;
// deactivate sites after stopping
			if(retval != 0)newmat->setActivity(false);
			sites.push_back(newmat);
// increment
			if(newmat->isActive())sitemom += newmat->momentumChange(dedxdir); // modify momentum for the next site
		}
    iindex += step;
    nstep--;
  }
  return retval;
}

void
KalRep::updateSites( int startindex,int endindex,
		     double initialmom,trkDirection dedxdir) {
// require a minimum momentum
  double sitemom = std::max(initialmom,_kalcon.minMom());
  int step = dedxdir==trkOut ? 1 : -1;
  int nstep = dedxdir==trkOut ? endindex-startindex+1 : startindex-endindex+1;
  int iindex(startindex);
// reset extendability (outwards)
  if(dedxdir == trkOut && !_extendable[dedxdir]){
    _fitrange[1] = std::max(_fitrange[1],_sites.back()->globalLength()+_kalcon.fltEpsilon());
    _extendable[dedxdir] = true;
  }
  while(nstep>0) {
    KalSite* thesite = _sites[iindex];
// if we're going outwards, update the site momentum _before_ the actual site.  This
// ensures consistency in how momentum is interpreted (sitemom is always outwards
// of the effect of this site)
    if(dedxdir == trkOut && thesite->isActive()){
      sitemom += thesite->momentumChange(dedxdir);
      sitemom = std::max(sitemom,_kalcon.minMom());
    }
    if( thesite->update(_reftraj,sitemom)){
// note update resets the activity flag of the material
// reset material sites
      if(thesite->kalMaterial() != 0 && thesite->isActive()){
        KalMaterial* mat = thesite->kalMaterial();
// if we're going out (losing energy), check for the track stopping
        if(dedxdir==trkOut && stopsIn(mat)){
          _stopsite = mat;
          _extendable[dedxdir] = false;
          _fitrange[1] = std::max(thesite->globalLength()+_kalcon.fltEpsilon(),_fitrange[0]+ _kalcon.minFltLen());
        }
// deactivate the materials past the stopping site
        if(_stopsite != 0 && _stopsite->globalLength() <= mat->globalLength())
          mat->setActivity(false);
      }
// update momentum for inwards processing
      if(dedxdir == trkIn && thesite->isActive()){
        sitemom += thesite->momentumChange(dedxdir);
        sitemom = std::max(sitemom,_kalcon.minMom());
    }
// check hot sites for flightlenght changes
      if(thesite->kalHit() != 0 && thesite->isActive())
	_maxfltdif = std::max(_maxfltdif,thesite->kalHit()->flightLengthChange());
    }
// move to the next site
    iindex +=step;
    nstep--;
  }
}

void
KalRep::fixupSites() {
// sites can move during update: re-order the ones that have
  std::sort(_sites.begin(),
	    _sites.end(),
	    babar::Collection::PtrLess());
// sort the hotlist too
  hotList()->sort();
// find the extreme hit sites again, in case they've moved
  findHitSites();
// sometimes the last hit will move outside the fit range.  This is not really an
// error, so just adjust the fit range appropriately
  if( _hitrange[0]<_sites.size() && _hitrange[1]>0){
    if(_extendable[trkIn] && _fitrange[0] > _sites[_hitrange[0]]->globalLength())
      _fitrange[0] = _sites[_hitrange[0]]->globalLength()+_kalcon.fltEpsilon();
    if(_extendable[trkOut] && _fitrange[1]< _sites[_hitrange[1]]->globalLength())
      _fitrange[1] = _sites[_hitrange[1]]->globalLength()+_kalcon.fltEpsilon();
  }
}


KalStub*
KalRep::createStub(const TrkVolume& extendvolume,
		   trkDirection extenddir,double tolerance,
		   const KalContext* kalcon) {
// preset to failure
  KalStub* retval(0);
// the rep must have been fit in the direction specified and be extendable in that direction
  if( hasFit(extenddir)){
// sites in both directions  need re-processing
    trkDirection odir = extenddir==trkOut ? trkIn : trkOut;
    _siteflag[odir] = false;
// invalidate all the sites in the other direction
    unsigned nsites = _sites.size();
    for(unsigned isite=0;isite<nsites;isite++)
      _sites[isite]->invalidateSite(odir);
// create a list for holding the sites
    std::deque<KalSite*> sites;
// Extend the trajectory through volume
    extendThrough(extendvolume,extenddir);
// set the range
    double xrange[2];
    xrange[0] = xrange[1] = extenddir == trkIn ?
      _ptraj->lowRange() : _ptraj->hiRange();
// give away the non-hit sites before the 'first' hit to the stub
// be careful of inactive hots and constraints!
    KalSite* firsthit = extenddir==trkIn ? _sites[_hitrange[0]] :
      _sites[_hitrange[1]];
    KalSite* willremove = extenddir==trkIn ? _sites.front() : _sites.back();
    while(willremove != firsthit && willremove->kalHit()==0 
	  && willremove->nDof() == 0){
// update the range
      xrange[0] = std::min(xrange[0],willremove->globalLength());
      xrange[1] = std::max(xrange[1],willremove->globalLength());
// check if we've given away the stopping site, and if so, unset it
      if(willremove == _stopsite)
	_stopsite = 0;

      if(extenddir==trkIn){
	_sites.erase(_sites.begin());
	sites.push_front(willremove);
	willremove = _sites.front();
      } else {
	_sites.pop_back();
	sites.push_back(willremove);
	willremove = _sites.back();
      }
    }
// explicitly sort the sites
    std::sort(sites.begin(),sites.end(),babar::Collection::PtrLess());
// create the KalStub from these
    retval = new KalStub(*this,extenddir,tolerance,xrange,sites,_kalcon);
// cleanup
    findHitSites();
  }
  return retval;
}

double
KalRep::estimatedMomDiff() const {
// kludge value; this estimates the total amount of material seen on average by
// tracks in terms of dE/dx.  Particle mass dependence (beta dependence) is already
// taken out in this formulation.  Estimated from data.
  static double dedxfactor(0.006); // GeV, assumes dE/dx ~ 1/beta^2 (dP/dx ~ 1/beta^3)
// compare the current momentum with the reference momentum at the reference flight
// length; this sets the scale for how much the true momentum might change next iteration
  double curmom = momentum(refMomFltLen()).mag();
  double dmom = curmom - refMomentum();
  double avgmom = (curmom + refMomentum())/2.0;
// estimated change in momentum scales as mass^2/momentum^3
  double energy = DetMaterial::particleEnergy(refMomentum(),particleType());
  double ddmom = dedxfactor*dmom*energy*sqr(particleType().mass())/
    sqr(sqr(avgmom));
  return ddmom;
}


bool 
KalRep::isFitable(TrkErrCode& err) const {
  bool fitable(true);
//  count the degrees of freedom, and abort if this is < the minimum
  if(enoughDofs()){
// Check for stopping in material before the last hit; use a tolerance, since the intersections aren't perfect
    if(_stopsite == 0 ||
       (_sites[_hitrange[1]]->globalLength() <= _stopsite->globalLength() + _kalcon.intersectionTolerance() )){
// check for parameters diverging
      if(_maxfltdif > _kalcon.divergeFlt()){
	if (ErrLogging(debugging)) {
	  std::cout<<"max fligth difference "<<_maxfltdif<<std::endl;
	}
	fitable = false;
	err = TrkErrCode(TrkErrCode::fail,KalCodes::diverge,"Iterations diverge");
      }
    } else {
      if (ErrLogging(debugging)) {
	if(_stopsite!=0)
	  std::cout<<"Will stop must be on "<<_sites[_hitrange[1]]->globalLength()<<" "
		   <<" but stop on length "<<_stopsite->globalLength()
		   <<" eps "<<_kalcon.intersectionTolerance()<<std::endl;
      }
      fitable = false;
      err = TrkErrCode(TrkErrCode::fail,KalCodes::stops,
		       "Track stops due to energy loss before last hit");
    }
  } else {
    fitable = false;
    err = TrkErrCode(TrkErrCode::fail,KalCodes::dof,"Insufficient degress of freedom");
  }
  return fitable;
}
// reset the state to allow 'fit' to work properly when the
// hots have changed
void
KalRep::resetFit() {
  _niter = 0; // allow iterations to start again
// sites in both directions  need re-processing
  _siteflag[trkIn] = false;
  _siteflag[trkOut] = false;
  setCurrent(false);
// invalidate all the sites
  unsigned nsites = _sites.size();
  for(unsigned isite=0;isite<nsites;isite++){
    _sites[isite]->invalidateSite(trkIn);
    _sites[isite]->invalidateSite(trkOut);
  }
// update the end sites
  updateEndSites(_kalcon.smearFactor());
}

// reset everything about the fit to construction level quantities
TrkErrCode
KalRep::resetAll(bool invert) {
  if(invert){
// loop over all sites and invert them.  This also inverts the hots
    unsigned nsites = _sites.size();
    for(unsigned isite=0;isite<nsites;isite++)
      _sites[isite]->invert();
// invert the seed
    _seedtraj->invert();
// re-order the sites
    std::sort(_sites.begin(),_sites.end(),babar::Collection::PtrLess());
// reverse the fit range
    double temp = _fitrange[0];
    _fitrange[0] = - _fitrange[1];
    _fitrange[1] = -temp;
  }
// reset the fit
  resetFit();
// delete existing trajectory
  if(_reftraj != _ptraj)
    delete _ptraj;
// recalculate initial momentum and charge from the seed.
// remake the trajectory by hand
  _ptraj = new TrkDifPieceTraj(*_seedtraj,_fitrange[0],_fitrange[1]);
  assert(_ptraj != 0);
// update the sites
  return updateSites();
}

void
KalRep::addConstraint(const TrkSimpTraj* traj,double cfltlen,bool* cparams) {
// construct a constraint site from this info
  KalConstraint* csite = new KalConstraint(_reftraj,
					   *traj->parameters(),
					   cparams,
					   cfltlen);
  if(csite != 0){
    _sites.push_back(csite);
    findHitSites(); // find the new range (constraint is counted as a HOT)
    resetFit(); // force refitting
  } else
    ErrMsg(error) << "cannot create constraint site" << endmsg;
}

void
KalRep::addSite(KalSite* site) {
  if(site != 0){
    _sites.push_back(site);
    findHitSites(); // find the new range (constraint is counted as a HOT)
    resetFit(); // force refitting
  } else
    ErrMsg(error) << "cannot add null site" << endmsg;
}

void
KalRep::setConstraints(bool* cparams) {
// loop over all sites and find the constraints
  bool changed(false);
  for(unsigned isite=0;isite<_sites.size();isite++){
    if(_sites[isite]->kalConstraint() != 0)
      changed |= _sites[isite]->kalConstraint()->setConstraint(cparams);
  }
  if(changed){
    setCurrent(false);
    findHitSites();
    resetFit(); // force refitting
  }
}

TrkErrCode
KalRep::buildTraj(trkDirection tdir) {
// clone the seed traj: cast it back to a simptraj
  TrkSimpTraj* straj = (TrkSimpTraj*)(_seedtraj->clone());
// set its parameters according the 'last' site
  int ilastsite = tdir == trkOut ? _hitrange[1] : _hitrange[0];
  KalSite* lastsite(_sites[ilastsite]);
  if(!lastsite->setTrajState(tdir,straj))
    return TrkErrCode(TrkErrCode::fail,KalCodes::matrix,
		      "Matrix inversion error");
// the traj will only be valid from the 'last' site out to the 'end'
// buffer the range slightly to avoid problems with degenerate pathlength sites 
  double range[2];
  if(tdir == trkOut){
    range[0] = lastsite->globalLength()- _kalcon.minFltLen();
    range[1] = std::max(_fitrange[1],range[0]+ _kalcon.minFltLen());
  } else {
    range[1] = lastsite->globalLength()+ _kalcon.minFltLen();
    range[0] = std::min(_fitrange[0],range[1]- _kalcon.minFltLen());
  }
  straj->setFlightRange(range);
// build the new piece traj starting with this piece
//  Delete the old trajectory, and set the new one
  if(_ptraj != _reftraj)
    delete _ptraj;
  _ptraj = new TrkDifPieceTraj(straj,range[0],range[1]);
  _maxdist = 0.0;
// extend the traj
  return extendTraj(ilastsite,tdir);
}

TrkErrCode
KalRep::extendTraj(int startsite,trkDirection tdir) {
  TrkErrCode retval;
// compute some index limits
  int endsite(0),nsites(0),sitestep(0);
  switch(tdir){
  case trkIn:
    endsite = 0;
    sitestep = -1;
    nsites = startsite-endsite+1;
    break;
  case trkOut:
    endsite = _sites.size()-1;
    sitestep = 1;
    nsites = endsite-startsite+1;
    break;
  }
  if(nsites > 0){
// kludge DNB_RKK
    static const HepPoint origin(0.0,0.0,0.0);
// loop over the new sites and build trajectory pieces from them
    int isite=startsite;
    while(nsites>0){
      KalSite* thesite = _sites[isite];
// exclude hit sites.  Yes, these can be here, if they are _inactive_
      if(thesite->nDof() == 0 && thesite->isActive()){
        KalSite* nextsite = nsites>1 ? _sites[isite+sitestep] : 0;
        if( nextsite == 0 ||
            fabs(thesite->globalLength() -nextsite->globalLength())> _kalcon.minGap()){
          double sitelen = thesite->globalLength();
// build the extra trajectory pieces for these sites
// clone the seed traj: cast it back to a simptraj
          TrkSimpTraj* straj = (TrkSimpTraj*)(_seedtraj->clone());
          if(thesite->setTrajState(tdir,straj)){
// mover inner parameters to reference the origin
            if(tdir == trkIn){
              double olen = thesite->localLength() - sitelen;
              straj->changePoint(origin,olen);
            }
// set the range on this trajectory; this helps the append function
            double fltrng[2];
            switch(tdir){
            case trkIn:
              fltrng[0] = _fitrange[0];
              fltrng[1] = std::max(sitelen,_fitrange[0]);
              break;
              case trkOut:
                fltrng[0] = thesite->localLength();
                fltrng[1] = std::max(fltrng[0] + _fitrange[1] - sitelen,fltrng[0]);
                break;
            }
            straj->setFlightRange(fltrng);
            double gap(0.0);
            const TrkErrCode& adderror = 
              (tdir==trkIn  ? _ptraj->prepend(sitelen,straj,gap): _ptraj->append(sitelen,straj,gap));
            if(adderror.success()){
              thesite->setGap(gap);
            } else {
// cleanup and abort
              retval = adderror;
              delete straj;
              break;;
            }
          } else {
            retval = TrkErrCode(TrkErrCode::fail,KalCodes::matrix,"Matrix inversion error");
            delete straj;
            break;
          }
        }
      }
      nsites--;
      isite += sitestep;
    }
  }
  return retval;
}


void
KalRep::setExtendable(trkDirection idir) {
  _extendable[idir] = true;
// make sure the fit range covers the sites
  switch(idir) {
  case trkIn:
    _fitrange[0] = std::min(_fitrange[0],_sites.front()->globalLength()-_kalcon.fltEpsilon());
    break;
  case trkOut:
    _fitrange[1] = std::max(_fitrange[1],_sites.back()->globalLength()+_kalcon.fltEpsilon());
    break;
  }
}

bool
KalRep::stopsIn(const KalMaterial* mat) const {
  return fabs(mat->momFraction()) > _kalcon.maxSiteDMom() ||
    fabs(mat->momentumChange(trkIn)) > _kalcon.maxDMom();
}

TrkErrCode
KalRep::removeStopSite() {
// only invoke if necessary
  if(_stopsite == 0 )
    return TrkErrCode(TrkErrCode::fail,KalCodes::nostop,"Track doesn't stop");
// find the stop site in the site list
  size_t stopper;
  for(stopper =0;stopper<_sites.size();stopper++){
    if(_sites[stopper] == _stopsite)
      break;
  }
// reset the sites on either side, to force re-processing
  if(stopper < _sites.size()){
    unsigned isite;
    unsigned nsites = _sites.size();
    for(isite=stopper+1;isite<nsites;isite++)
      _sites[isite]->reset(trkOut);
    for(isite=0;isite<stopper;isite++)
      _sites[isite]->reset(trkIn);
// remove and delete the stopping site
    if(_sites[stopper]  == _stopsite){
      _sites.erase(_sites.begin()+stopper);
      delete _stopsite;
      _stopsite = 0;
    } else {
      return TrkErrCode(TrkErrCode::fail,KalCodes::inconsistent,"KalRep state is inconsistent");
    }
// refind the hits
    findHitSites();
// reset the fit
    resetFit();
  } else
    return TrkErrCode(TrkErrCode::fail,KalCodes::inconsistent,"KalRep state is inconsistent");
  return TrkErrCode(TrkErrCode::succeed);
}

bool
KalRep::validFlightLength(double fltL,double tolerance) const {
// start with basic information
  bool retval = _ptraj != 0 && 
    fltL >= _ptraj->lowRange()-tolerance &&
    fltL <= _ptraj->hiRange()+tolerance;
// if there are no hots, include test on found range
  return retval;
}

void
KalRep::constraintTrajectories(std::vector<TrkSimpTraj*>& cons) const {
  unsigned nsites = _sites.size();
  for(unsigned isite=0;isite<nsites;isite++){
    if(_sites[isite]->kalConstraint() != 0){
// create a TrkSimpTraj from the site.  give it a very small range
      const KalConstraint* csite = _sites[isite]->kalConstraint();
      cons.push_back(const_cast<TrkSimpTraj*>(csite->constraintTrajectory()));
    }
  }
}

void
KalRep::reIntersect() {
// increment # of intersections
  ++_ninter;
// find all material and bend sites
  FindMatBendSites matbend;
  std::vector<KalSite*>::iterator mid =
    std::partition(_sites.begin(),_sites.end(),matbend);
// delete and erase removed sites.
  std::for_each(mid,_sites.end(),babar::Collection::DeleteObject());
  _sites.erase(mid,_sites.end());
// re-calculate the fitrange
  _fitrange[0] = hotList()->startFoundRange()-_kalcon.fltEpsilon();
  _fitrange[1] = hotList()->endFoundRange()+_kalcon.fltEpsilon();
  setFitRange(_reftraj);
// re-make materials and bends using current ref traj
	std::vector<DetIntersection> tlist; tlist.reserve(64);
  buildMaterialSites(_fitrange,tlist);
  if(_kalcon.bendSites() )
    buildBendSites(_fitrange);
}


void
KalRep::updateEndSites(double smear,bool diagonly) {
// diagonalization doesn't seem to work
  diagonly = false;
  double firstlen = _fitrange[0];
  double lastlen = _fitrange[1];
// use the hits if possible
  if(_hitrange[0] < _sites.size() && _hitrange[1]> 0){
    firstlen = _sites[_hitrange[0]]->globalLength();
    lastlen = _sites[_hitrange[1]]->globalLength();
  }
// only update on the first N iterations
  if(_niter <= 100){
    _endsites[trkIn] = KalEndSite(_reftraj,firstlen,trkOut,smear,diagonly);
    _endsites[trkOut] = KalEndSite(_reftraj,lastlen,trkIn,smear,diagonly);
  } else {
    _endsites[trkIn].update(_reftraj,firstlen);
    _endsites[trkOut].update(_reftraj,lastlen);
  }
}

KalSite* 
KalRep::nextActive(unsigned index,trkDirection tdir) const {
  KalSite* retval(0);
  int step = tdir == trkOut ? 1 : -1;
  int jndex = index+step;
  while(jndex >= 0 && jndex < _sites.size()){
    if( _sites[jndex]->isActive()){
      retval = _sites[jndex];
      break;
    }
    jndex += step;
  }
  return retval;
}

const KalMaterial*
KalRep::findMaterial(double fltlen,double minrad) const {
  unsigned nsites = _sites.size();
  double dist(1.0e12);
  const KalMaterial* retval(0);
  for(unsigned isite=0;isite<nsites;isite++){
    const KalSite* site = _sites[isite];
    if(site->kalMaterial() != 0 &&
       site->kalMaterial()->radiationFraction() >= minrad &&
       fabs(site->globalLength()-fltlen) < dist){
      dist = fabs(site->globalLength()-fltlen);
      retval = site->kalMaterial();
    }
  }
  return retval;
}

const KalSite*
KalRep::findSite(KalSite::siteType stype) const {
  const KalSite* retval(0);
  for(unsigned isite=0; isite < _sites.size(); isite++){
    if(_sites[isite]->type() == stype){
      retval = _sites[isite];
      break;
    }
  }
  return retval;
}

//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: KalMiniRep.cc,v 1.74 2007/09/06 05:21:20 brownd Exp $
//
//   Description: class KalMiniRep.
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Copyright Information:
//	Copyright (C) 2000	Lawrence Berkeley Laboratory
//
// Author List:
//      Dave Brown 11/6/00
//------------------------------------------------------------------------

#include "BaBar/BaBar.hh"
#include "KalmanTrack/KalMiniRep.hh"
#include "KalmanTrack/KalRep.hh"
#include "KalmanTrack/KalCodes.hh"
#include "KalmanTrack/KalConstraint.hh"
#include "KalmanTrack/KalMiniInterface.hh"
#include "KalmanTrack/KalSite.hh"
#include "KalmanTrack/KalBetaCons.hh"
#include "TrkBase/TrkSimpTraj.hh"
#include "TrkBase/TrkDifPieceTraj.hh"
#include "TrkBase/TrkHotListFull.hh"
#include "TrkBase/TrkHotListUnowned.hh"
#include "TrkBase/TrkHotListEmpty.hh"
#include "TrkBase/TrkExchangePar.hh"
#include "TrkBase/TrkRecoTrk.hh"
#include "TrkBase/TrkDifTraj.hh"
#include "TrkBase/TrkMomCalculator.hh"
#include "TrkBase/TrkVolume.hh"
#include "ProxyDict/IfdIntKey.hh"
#include "ErrLogger/ErrLog.hh"
#include "BbrStdUtils/BbrCollectionUtils.hh"
#include <iostream>
using std::endl;
using std::ostream;

bool KalMiniRep::_useKalRep(true);
bool KalMiniRep::_autoExtend(false);

void
KalMiniRep::setAutoExtend(bool autoex) {
  _autoExtend = autoex;
}

KalMiniRep::KalMiniRep(TrkRecoTrk* trk,PdtPid::PidType hypo,
                       TrkSimpTraj* seed,
                       const KalContext& context,
                       DnaPtr<TrkHotList> fullhotlist,
                       TrkHotList* emptyhotlist,
                       const std::vector<TrkSimpTraj*>& constraints,
                       const std::vector<TrkSimpTraj*>& fitlist,
                       double chiprob) :
  TrkRep(trk,hypo),
  _kalcon(context),
  _hotrep(0),
  _xrep(0),
  _fullhotlist(fullhotlist),
  _emptyhotlist(emptyhotlist),
  _seed(seed),
  _fittrajs(fitlist),
  _contrajs(constraints),
  _state(cache),
  _ptraj(0),
  _consistency(1.0,1),
  _usesvtdedx(false)
{
// sanity checks
   bool fullhotlist_transok = fullhotlist.valid() &&
       (fullhotlist.rawPtr() != 0);
    if(fullhotlist_transok) assert(fullhotlist->hitCapable());
  assert(fullhotlist_transok || emptyhotlist != 0);
  setMultScat(true);   // in TrkFitStatus
// setup trajectory
  if(!_fittrajs.empty())
// if there are fits, turn them into a trajectory, otherwise use the seed
    _ptraj.reset( new TrkDifPieceTraj(_fittrajs) );
  else 
    _ptraj.reset( new TrkDifPieceTraj(*_seed,_seed->lowRange(),_seed->hiRange()) );
// make the hots on the full list  point to this rep
  if(fullhotlist_transok)
    std::for_each(_fullhotlist->begin(),_fullhotlist->end(),setParent(this));
// set the consistency
  int ndof = nDof();
  if(ndof<0) ndof=0;
  _consistency = ChisqConsistency((unsigned)ndof,chiprob);
// default setting of status, validity, and history are done in the calling code
}

// constructor for non-default reps on a track, including constraints, fits, and chisq cache.

KalMiniRep::KalMiniRep(KalMiniRep* defrep, PdtPid::PidType hypo,
		       const std::vector<TrkSimpTraj*>& constraints,
		       const std::vector<TrkSimpTraj*>& fitlist,
		       double chiprob) :
  TrkRep(defrep->parentTrack(),hypo),
  _kalcon(defrep->kalContext()),
  _hotrep(0),
  _xrep(0),
  _fullhotlist(0),
  _emptyhotlist(0),
  _seed(defrep->seedTrajectory()->clone()),
  _fittrajs(fitlist),
  _contrajs(constraints),
  _state(cache),
  _ptraj(0),
  _consistency(1.0,1),
  _usesvtdedx(false)
{
  setMultScat(true);   // in TrkFitStatus
// setup trajectory
  if(!_fittrajs.empty())
// if there are fits, turn them into a trajectory, and consider the fit current
    _ptraj.reset( new TrkDifPieceTraj(_fittrajs) );
  else
    _ptraj.reset( new TrkDifPieceTraj(*_seed,_seed->lowRange(),_seed->hiRange()) );
// if the default rep has an empty hotlist, clone it.  Don't clone the full hot list
  if(defrep->emptyHotList() != 0){
    TrkBase::Functors::cloneHot cloner(this);
    _emptyhotlist.reset(defrep->emptyHotList()->clone(cloner));
  } else
// create an empty hotlist from the def rep's full hot list
    _emptyhotlist.reset(new TrkHotListEmpty(*defrep->hotList()));
// set the consistency: must use the default rep to get ndof as this rep does not
// necessarily have a hotlist
  int ndof = defrep->nDof();
  if(ndof<0) ndof=0;
  _consistency = ChisqConsistency((unsigned)ndof,chiprob);
// default setting of status, validity, and history are done in the calling code
}


// copy constructor for cloning to a new track.  In this case, copy the full state
KalMiniRep::KalMiniRep(const KalMiniRep& rep,TrkRecoTrk* trk) :
  TrkRep(rep,trk,rep.particleType()),
  _kalcon(rep._kalcon),
  _hotrep(rep._hotrep.get() != 0 ? rep._hotrep->clone(trk) : 0),
  _xrep(rep._xrep.get() != 0 ? rep._xrep->clone(trk) : 0),
  _fullhotlist(0),
  _emptyhotlist(rep._emptyhotlist.get() != 0 ? rep._emptyhotlist->clone(this) : 0),
  _seed(((TrkSimpTraj*)rep._seed->clone())),
  _state(rep._state),
  _ptraj(0),
  _consistency(rep._consistency),
  _usesvtdedx(false)
{
  _contrajs.reserve(rep._contrajs.size());
  _fittrajs.reserve(rep._fittrajs.size());
// clone the mini-reps own trajectory
  _ptraj.reset( (TrkDifPieceTraj*)rep._ptraj->clone() );
// some care is needed cloning the full hotlists, as the hots may be unowned by
// this rep and the pointers need to be reset
  if(_hotrep.get() != 0)
// simply copy a reference to the real list owned by the hotrep
    _fullhotlist.assign( new TrkHotListUnowned((_hotrep->hotList())));
  else if (rep._fullhotlist.rawPtr() != 0)
    _fullhotlist.assign(rep._fullhotlist->clone(this));
// clone the constraint and fit trajectories
  for(TrkSimpTrajs::const_iterator icon=rep._contrajs.begin();
      icon!=rep._contrajs.end();++icon)
    _contrajs.push_back((TrkSimpTraj*)((*icon)->clone()));
  for(TrkSimpTrajs::const_iterator ifit=rep._fittrajs.begin();
      ifit!=rep._fittrajs.end();++ifit)
    _fittrajs.push_back((TrkSimpTraj*)((*ifit)->clone()));
  setValid(rep.fitValid());
  setCurrent(rep.fitCurrent());
}

// copy constructor for cloning new hypos.  This only works for a hot-base
// mini-rep
KalMiniRep::KalMiniRep(const KalMiniRep& rep,PdtPid::PidType hypo) :
  TrkRep((TrkRecoTrk*)rep.parentTrack(),hypo),
  _kalcon(rep._kalcon),
  _hotrep(0),
  _xrep(0),
  _fullhotlist(0),
  _emptyhotlist(0),
  _seed(((TrkSimpTraj*)rep._seed->clone())),
  _state(KalMiniRep::hots),
  _ptraj(0),
  _consistency(rep._consistency),
  _usesvtdedx(false)
{
//  check that there's really hots here
  if (rep.fullHotList() ==0){
    ErrMsg(error) << "Copying KalMiniRep without hots! state willl be set to seed" << endmsg;
    _state = KalMiniRep::seed;
  } else
// create a copy of the hotlist
    _fullhotlist.assign( new TrkHotListFull(*rep.fullHotList(),this) );
// create a new piecetraj
  _ptraj.reset( new TrkDifPieceTraj(*_seed,_seed->lowRange(),_seed->hiRange()) );
// create the hots rep
  TrkErrCode err = createHotKalRep();
  if(err.success() && _hotrep.get() != 0){
    setCurrent(_hotrep->fitCurrent());
    setValid(_hotrep->fitValid());
  } else {
// complete failure
    setValid(true);
    setCurrent(false);
    _state = KalMiniRep::seed;
    ErrMsg(error) << "Failed to clone hot-based mini rep" << endmsg;
  }
}

KalMiniRep::~KalMiniRep()
{
//
// delete my hots first, to avoid problems with the Hot destructor
 _fullhotlist.destroy();
 // clear out the fit and constraint trajectories
  std::for_each(_fittrajs.begin(),_fittrajs.end(),babar::Collection::DeleteObject());
  std::for_each(_contrajs.begin(),_contrajs.end(),babar::Collection::DeleteObject());
}



KalMiniRep*
KalMiniRep::clone(TrkRecoTrk* newTrack) const {
  return new KalMiniRep(*this,newTrack);
}


TrkRep*
KalMiniRep::cloneNewHypo(PdtPid::PidType hypo) {
  return new KalMiniRep(*this,hypo);
}

// all the tracking functions are implemented as call-down to the
// appropriate KalRep.  If the state is cache, I return null or
// something fake

TrkDifTraj& 
KalMiniRep::traj() {
  if(active())
    return kalRep()->traj();
  else
    return *_ptraj;
}

const TrkDifTraj& 
KalMiniRep::traj() const {
  if(active())
    return kalRep()->traj();
  else
    return *_ptraj;
}

int
KalMiniRep::nDof() const {
// extended (constrained) KalReps don't have chisq or NDOF values; use the
// MiniRep for that
  if(_state == hots)
    return kalRep()->nDof();
  else
    return hotList()->nActive()-_seed->parameters()->parameter().num_row();
}

double
KalMiniRep::chisq() const {
  if(_state == hots)
    return kalRep()->chisq();
  else
    return _consistency.chisqValue();
}

int               
KalMiniRep::charge() const {
  if(active())
    return kalRep()->charge();
  else {
    double midpoint = (_seed->lowRange() + _seed->hiRange())/2.0;
// if fit trajectories exist, use that, otherwise the seed
    const TrkSimpTraj* traj = fitTraj(0.0);
    if(traj == 0) traj = _seed.get();
    return TrkMomCalculator::charge(*traj,parentTrack()->bField(),midpoint);
  }
}

HepPoint
KalMiniRep::position(double fltL)         const {
  extendIfNeeded(fltL);
  return TrkRep::position(fltL);
}

Hep3Vector
KalMiniRep::direction(double fltL)        const {
  extendIfNeeded(fltL);
  return TrkRep::direction(fltL);
}

BbrPointErr
KalMiniRep::positionErr(double fltL)      const {
  extendIfNeeded(fltL);
  return TrkRep::positionErr(fltL);
}

BbrVectorErr
KalMiniRep::directionErr(double fltL)     const {
  extendIfNeeded(fltL);
  return TrkRep::directionErr(fltL);
}

Hep3Vector
KalMiniRep::momentum(double fltL) const {
  extendIfNeeded(fltL);
  if(active())
    return kalRep()->momentum(fltL);
  else {
    double localflt;
    const TrkSimpTraj* locTraj = _ptraj->localTrajectory(fltL,localflt);
    if(locTraj != 0){
      return TrkMomCalculator::vecMom(*locTraj, parentTrack()->bField(), localflt);
    } else{
      ErrMsg(error) << "Can't find local trajectory: using seed" << endmsg;
      return TrkMomCalculator::vecMom(*_seed, parentTrack()->bField(), 0.0);
    }
  }
}

double
KalMiniRep::pt(double fltL) const {
  extendIfNeeded(fltL);
  if(active())
    return kalRep()->pt(fltL);
  else
    return momentum(fltL).perp();
}

BbrVectorErr
KalMiniRep::momentumErr(double fltL) const {
  extendIfNeeded(fltL);
  if(active())
    return kalRep()->momentumErr(fltL);
  else {
    double localflt;
    const TrkSimpTraj* locTraj = _ptraj->localTrajectory(fltL,localflt);
    if(locTraj != 0){
      return TrkMomCalculator::errMom(*locTraj, parentTrack()->bField(), localflt);
    } else {
      ErrMsg(error) << "Can't find local trajectory: using seed" << endmsg;
      return TrkMomCalculator::errMom(*_seed, parentTrack()->bField(), 0.0);
    }
  }
}

TrkExchangePar
KalMiniRep::helix(double fltL) const {
  extendIfNeeded(fltL);
  if(active())
    return kalRep()->helix(fltL);
  else {
    double localflt;
    const TrkSimpTraj* ltraj = _ptraj->localTrajectory(fltL,localflt);
    return TrkExchangePar(ltraj->parameters()->parameter(), 
			  ltraj->parameters()->covariance());
  }
}

HepMatrix
KalMiniRep::posmomCov(double fltL) const {
  extendIfNeeded(fltL);
  if(active())
    return kalRep()->posmomCov(fltL);
  else {
    double localflt;
    const TrkSimpTraj* locTraj = _ptraj->localTrajectory(fltL,localflt);
    if(locTraj != 0){
      return TrkMomCalculator::posmomCov(*locTraj, parentTrack()->bField(), localflt);
    } else {
      ErrMsg(error) << "Can't find local trajectory: using seed" << endmsg;
      return TrkMomCalculator::posmomCov(*_seed, parentTrack()->bField(), 0.0);
    }
  }
}

void
KalMiniRep::getAllCovs(double fltL,
		       HepSymMatrix& xxCov,
		       HepSymMatrix& ppCov,
		       HepMatrix& xpCov) const {
  extendIfNeeded(fltL);
  if(active())
    kalRep()->getAllCovs(fltL,xxCov,ppCov,xpCov);
  else {
    double localflt;
    const TrkSimpTraj* locTraj = _ptraj->localTrajectory(fltL,localflt);
// check if this local flightlength is really on the local trajectory
    if(locTraj != 0)
      TrkMomCalculator::getAllCovs(*locTraj, parentTrack()->bField(), localflt,
				   xxCov,ppCov,xpCov);
    else {
      ErrMsg(error) << "Can't find local trajectory: using seed" << endmsg;
      TrkMomCalculator::getAllCovs(*_seed, parentTrack()->bField(), 0.0,
				   xxCov,ppCov,xpCov);
    }
  }
}

void
KalMiniRep::getAllWeights(double fltL,
			  HepVector& pos,
			  HepVector& mom,
			  HepSymMatrix& xxWeight,
			  HepSymMatrix& ppWeight,
			  HepMatrix&    xpWeight) const {
  extendIfNeeded(fltL);
  if(active())
    kalRep()->getAllWeights(fltL,pos,mom,xxWeight,ppWeight,xpWeight);
  else {
    double localflt;
    const TrkSimpTraj* locTraj = _ptraj->localTrajectory(fltL,localflt);
    if(locTraj != 0)
      TrkMomCalculator::getAllWeights(*locTraj, parentTrack()->bField(), localflt,
				      pos,mom,xxWeight,ppWeight,xpWeight);
    else {
      ErrMsg(error) << "Can't find local trajectory: using seed" << endmsg;
      TrkMomCalculator::getAllWeights(*_seed, parentTrack()->bField(), 0.0,
				      pos,mom,xxWeight,ppWeight,xpWeight);
    }
  }
}

const IfdKey&     
KalMiniRep::myKey() const {
  static IfdIntKey _theKey(2154);
  return _theKey;
}

// hot-based operations.  Defer to the hot rep if that exists, otherwise perform the operation
// (it will have no effect on the cached result, but it will propogate to the hot rep if that
// is subsequently created.
void
KalMiniRep::addHot(TrkHitOnTrk *theHot) {
  if(_state == KalMiniRep::hots){
// repoint the hot to point to the KalRep underneath before adding
    setParent(*theHot,_hotrep.get());
    _hotrep->addHot(theHot);
  } else
    TrkRep::addHot(theHot);
}

void
KalMiniRep::removeHot(TrkHitOnTrk *theHot) {
  if(_state == KalMiniRep::hots){
    _hotrep->removeHot(theHot);
  } else
    TrkRep::removeHot(theHot);
}

void
KalMiniRep::activateHot(TrkHitOnTrk *theHot) {
  if(_state == KalMiniRep::hots){
    _hotrep->activateHot(theHot);
  } else
    TrkRep::activateHot(theHot);
}

void
KalMiniRep::deactivateHot(TrkHitOnTrk *theHot) {
  if(_state == KalMiniRep::hots){
    _hotrep->deactivateHot(theHot);
  } else
    TrkRep::deactivateHot(theHot);
}

bool
KalMiniRep::resid(const TrkHitOnTrk *theHot, 
		  double &residual, double &residErr,
		  bool exclude) const {
  if(_state == KalMiniRep::hots)
    return _hotrep->resid(theHot,residual,residErr,exclude);
  else
    return false;
}

TrkErrCode
KalMiniRep::fit() {
  TrkErrCode retval(TrkErrCode::succeed);
// always defer to the sub-rep if that exists
  if(active()){
    retval = kalRep()->fit();
// propogate the state
    setValid(kalRep()->fitValid());
    setCurrent(kalRep()->fitCurrent());
  } else {
// nothing to do if the fit is already current
    if(fitCurrent())
      retval =  TrkErrCode(TrkErrCode::succeed,KalCodes::current,
			 "Fit is already current");
    else if (fitValid()){
      setCurrent(true);
    }
  }
  return retval;
}

void
KalMiniRep::printAll(ostream& ostr) const {
// print what's here
  static char* names[3] = {"Cache","Hot Rep Active","Fit Rep Active"};
  ostr << "KalMiniRep for ";
  TrkRep::printType(ostr);
  ostr << " in state " << names[_state] <<" with fit probability " 
       << _consistency.significanceLevel() << endl;
  ostr << "Seed trajectory = ";_seed->printAll(ostr);
  ostr << _fittrajs.size() << " fit results, as follows:" << endl;
  for(unsigned ifit=0;ifit<_fittrajs.size();ifit++){
    ostr << "Fit " <<ifit << " : ";
    _fittrajs[ifit]->printAll(ostr);
  }
  ostr << _contrajs.size() << " constraints, as follows:" << endl;
  for(unsigned ihot=0;ihot<_contrajs.size();ihot++){
    ostr << "Constraint " <<ihot << " : ";
    _contrajs[ihot]->printAll(ostr);
  }
  if(_hotrep.get() != 0){
    ostr << "Hot Rep ";
    _hotrep->printAll(ostr);
  }
  if(_xrep.get() != 0){
    ostr << "Fit Rep ";
    _xrep->printAll(ostr);
  }
}

void
KalMiniRep::print(ostream& ostr) const {
// print what's here
  static char* names[3] = {"Cache","Hot Rep Active","Fit Rep Active"};
  ostr << "KalMiniRep for ";
  TrkRep::printType(ostr);
  ostr << " in state " << names[_state] <<" with fit probability " 
       << _consistency.significanceLevel() << endl;
  ostr << "Seed trajectory = ";_seed->print(ostr);
  ostr << _fittrajs.size() << " fit results, as follows:" << endl;
  for(unsigned ifit=0;ifit<_fittrajs.size();ifit++){
    ostr << "Fit " <<ifit << " : ";
    _fittrajs[ifit]->print(ostr);
  }
  ostr << _contrajs.size() << " constraints, as follows:" << endl;
  for(unsigned icon=0;icon<_contrajs.size();icon++){
    ostr << "Constraint " <<icon << " : ";
    _contrajs[icon]->print(ostr);
  }
  if(_hotrep.get() != 0){
    ostr << "Hot Rep ";
    _hotrep->print(ostr);
  }
  if(_xrep.get() != 0){
    ostr << "Fit Rep ";
    _xrep->print(ostr);
  }
}
// mini-rep specific functions


const KalRep*
KalMiniRep::kalRep(KalMiniRep::miniState state) const {
  switch (state ) {
  case KalMiniRep::hots:
    return _hotrep.get();
  case KalMiniRep::extendedcache:
    return _xrep.get();
  case KalMiniRep::cache: default:
    return 0;
  }
}

KalRep*
KalMiniRep::kalRep(KalMiniRep::miniState state) {
  switch (state ) {
  case KalMiniRep::hots:
    return _hotrep.get();
  case KalMiniRep::extendedcache:
    return _xrep.get();
  case KalMiniRep::cache: default:
    return 0;
  }
}

TrkErrCode
KalMiniRep::changeState(KalMiniRep::miniState newstate) {
  TrkErrCode retval(TrkErrCode::succeed,KalCodes::unchanged,"KalMiniRep state unchanged");
  if(_state != newstate){
    switch(newstate){
    case (KalMiniRep::hots):
      if( _hotrep.get() == 0)
	retval = createHotKalRep();
      else
	retval = TrkErrCode(TrkErrCode::succeed);
// update status of this rep to reflect the hots rep
      if(retval.success() && _hotrep.get() != 0){
	setCurrent(_hotrep->fitCurrent());
	setValid(_hotrep->fitValid());
	_state = newstate;
      }
      break;
    case (KalMiniRep::extendedcache):
      if(_xrep.get() == 0)
	retval = createXKalRep();
      else
	retval = TrkErrCode(TrkErrCode::succeed);
// update status of this rep to reflect the hots rep
      if(retval.success() && _xrep.get() != 0){
	setCurrent(_xrep->fitCurrent());
	setValid(_xrep->fitValid());
	_state = newstate;
      }
      break;
    case (KalMiniRep::cache):
// restore the constructor settings.  The hots will be valid if
// either I've updated them, or the hotrep exists.
      if(_fittrajs.size() > 0){
	setFoundRange();
	setCurrent(true);
	setValid(true);
	retval = TrkErrCode(TrkErrCode::succeed);
	_state = newstate;
      } else {
	setCurrent(false);
	setValid(false);
	retval = TrkErrCode(TrkErrCode::fail,KalCodes::nocache,
			    "No cache exists for this PID");
      }
      break;
    case (KalMiniRep::seed):
// use the seed for everything
      _ptraj.reset( new TrkDifPieceTraj(*_seed,_seed->lowRange(),_seed->hiRange()) );
      setCurrent(true);
      setValid(true);      
      _state = newstate;
      break;
    }
  }
  return retval;
}

TrkErrCode
KalMiniRep::createHotKalRep() {
  TrkErrCode retval(TrkErrCode::succeed);
// if we're a non-default rep, we want to clone the KalRep from the default rep (if possible)
  static KalMiniInterface kiface;
  KalMiniRep* defrep(0);
  TrkRecoTrk* parent = parentTrack();
  if(parent->defaultType() == particleType())
    defrep = this;
  else if(parent->attach(kiface,parent->defaultType()))
    defrep = kiface.kalMiniRep();
  if(defrep != 0 && defrep != this) {
// first, see if the default rep has a kalrep.  If not, try to force it to make one
    const KalRep* krep = defrep->kalRep(hots);
    if(krep == 0){
      retval = defrep->createHotKalRep();
      if(!retval.success())return retval;
      krep = defrep->kalRep(hots);
    }
// clone the kalrep to my own hotrep
    _hotrep.reset(new KalRep(*krep,particleType()));
    assert(_hotrep.get() != 0);
// create an unowned hotlist as my new full hotlist
    _fullhotlist.destroy();
    _fullhotlist.assign(new TrkHotListUnowned(_hotrep->hotList()));
// otherwise, create the KalRep from the hot list.  The hotrep will take
// ownership of the hots from the mini (if it has them), otherwise
// it clones new hots itself
  } else if(_fullhotlist.rawPtr() != 0){
//  Use the cache'd origin trajectory if it exists as seed, otherwise use the
// original seed (Backed out 9/3/02 by DNB: this causes inconsistencies when
// re-persisting mini tracks)
    const TrkSimpTraj* seedtraj(_seed.get());
// create an unowned hot list to give the kalrep
    DnaPtr<TrkHotList> hotl(0); 
    hotl.assign( new TrkHotListUnowned(_fullhotlist.rawPtr()));
// create the new rep
    _hotrep.reset( new KalRep(*seedtraj,
                              hotl.rawPtr(), // DnaPtr is unowned
                              parentTrack(),
                              _kalcon,particleType(),0.0,0,true) );
    assert(_hotrep.get() != 0);
// add the constraints.  By default, all parameters will be constrained,
// though the set of parameters to be constrained can be modified before calling 'fit'
    unsigned ncons = _contrajs.size();
    for(unsigned icon=0;icon<ncons;icon++){
      double cfltlen = (_contrajs[icon]->lowRange()+_contrajs[icon]->hiRange())/2.0;
      _hotrep->addConstraint(_contrajs[icon],cfltlen,0);
    }
// create history
    _hotrep->addHistory(retval,"Created by KalMiniRep");
  } else 
    retval = TrkErrCode(TrkErrCode::fail,KalCodes::nohots);
// create Svt dedx constraint if necessary
  setSvtdEdx(_usesvtdedx);
  return retval;
}

TrkErrCode
KalMiniRep::createXKalRep() {
  TrkErrCode retval(TrkErrCode::succeed);
// create the KalRep using the cached fit results
  if(!_fittrajs.empty()){
// create an empty hot list
    TrkHotListEmpty emptyhots(*hotList());
// create the rep
    _xrep.reset( new KalRep(parentTrack(),particleType(),
                            *_fittrajs[0],_kalcon,emptyhots,
                            _fittrajs) );
    assert(_xrep.get() != 0);
// create history
    _xrep->addHistory(retval,"Created by KalMiniRep");
// fitting is left for the user
  } else
     retval = TrkErrCode(TrkErrCode::fail,KalCodes::nocache,
                            "No cache exists for this PID");
  return retval;
}

ChisqConsistency
KalMiniRep::chisqConsistency() const {
  if(_state == hots)
    return kalRep()->chisqConsistency();
  else
    return _consistency;
}

bool
KalMiniRep::validFlightLength(double fltL,double tolerance) const {
  bool retval(false);
  if(active())
    retval =  kalRep()->validFlightLength(fltL,tolerance);
  else
    retval = 
      fltL >= _ptraj->lowRange()-tolerance &&
      fltL <= _ptraj->hiRange()+tolerance;
  return retval;
}

TrkErrCode
KalMiniRep::extendThrough(const TrkVolume& vol,trkDirection trkDir) const {
  if(active()){
    KalRep* krep = const_cast<KalRep*>(kalRep());
    return krep->extendThrough(vol,trkDir);
  } else
    return TrkErrCode(TrkErrCode::fail,KalCodes::mustbeactive,
		      "KalMiniRep in cache state cannot perform this operation");
}

TrkErrCode
KalMiniRep::extendThrough(double fltlen) const {
// only extend if necessary
  if(active()){
    KalRep* krep = const_cast<KalRep*>(kalRep());
    return krep->extendThrough(fltlen);
  } else
    return TrkErrCode(TrkErrCode::fail,KalCodes::mustbeactive,
		      "KalMiniRep in cache state cannot perform this operation");
}

void
KalMiniRep::updateHots() {
  if(_hotrep.get() != 0)
    _hotrep->updateHots();
  else if(_fullhotlist.rawPtr() != 0){
    TrkHotList::nc_hot_iterator end = _fullhotlist->end();
    for(TrkHotList::nc_hot_iterator ihot=_fullhotlist->begin();ihot!=end;++ihot) {
      if(ihot->hasResidual()) updateMeasurement(*ihot,_ptraj.get(),true);
    }
  }
}

void
KalMiniRep::setFoundRange() {
// the first and last hots must be up-to-date to establish the start and end found range
  if(_fullhotlist.rawPtr() != 0){
    TrkHotList::nc_hot_iterator fhot = _fullhotlist->begin();
    if(fhot.get() != 0 && (!fhot->hasResidual()))
      updateMeasurement(*fhot,_ptraj.get(),true);
    TrkHotList::nc_hot_iterator lhot = _fullhotlist->end();
    --lhot;
    if(lhot.get() != 0 && (!lhot->hasResidual()))
      updateMeasurement(*lhot,_ptraj.get(),true);
  }
}

const TrkSimpTraj*
KalMiniRep::localTrajectory(double fltlen,double& loclen) const {
  extendIfNeeded(fltlen);
  if(active())
    return kalRep()->localTrajectory(fltlen,loclen);
  else
    return pieceTraj().localTrajectory(fltlen, loclen);
}

const char*
KalMiniRep::stateName(KalMiniRep::miniState state) {
  switch(state) {
  case KalMiniRep::seed:
    return "Seed";
    break;
  case KalMiniRep::cache:
    return "Cache";
    break;
  case KalMiniRep::extendedcache:
    return "ExtendedCache";
    break;
  case KalMiniRep::hots:
    return "Hots";
    break;
  default:
    return "Unknown";
    break;
  }
}


const TrkSimpTraj* 
KalMiniRep::fitTraj(double fltlen) const {
  static double tol(2.0); // tolerance for finding trajectories
  const TrkSimpTraj* retval(0);
  if(_hotrep.get() != 0 && _hotrep->fitCurrent()){
// use the fit result
    double localflt;
    retval = _hotrep->localTrajectory(fltlen,localflt);
  } else {
// loop over the fits to find one that matches.  Allow some tolerance
    for(TrkSimpTrajs::const_iterator ifit=_fittrajs.begin();
                               ifit!=_fittrajs.end();++ifit){
      if((*ifit)->validFlightDistance(fltlen,tol)){
	retval = *ifit;
	break;
      }
    }
// failsafe
    if(retval == 0){
      ErrMsg(debugging) << "Requested fit not cached, using local trajectory" << endmsg;
      double localflt(0.0);
      retval = localTrajectory(fltlen,localflt);
    }
  }
  return retval;
}

void
KalMiniRep::fitTrajectories(std::vector<TrkSimpTraj*>& fits,const char* stream) const {
// If we're using the hotrep and we are configured to use it, take the fits from the
// KalRep
  if(_state == KalMiniRep::hots && _hotrep.get() != 0 && useKalRep()){
// make sure the hot rep is -CURRENT- before doing this
    if(_hotrep->fitCurrent()){
      _hotrep->fitTrajectories(fits,stream);
      return;
    } else 
      ErrMsg(error) << "Refit KalRep is not current!  Using MiniRep cached fits " << endmsg;
  }
// otherwise use the cached fits.  Someday we should filter these to make sure they
// agree with the requests.  FIXME!!!!!
  std::vector<TrkSimpTraj*>::const_iterator ifit = _fittrajs.begin();
  while(ifit != _fittrajs.end()){
    fits.push_back(*ifit);
    ifit++;
  }
}

void
KalMiniRep::constraintTrajectories(std::vector<TrkSimpTraj*>& cons) const {
// If we're using the hotrep and we are configured to use it, take the constraints
// KalRep
  if(_state == KalMiniRep::hots && _hotrep.get() != 0 && useKalRep()){
// make sure the hot rep is -CURRENT- before doing this
    if(_hotrep->fitCurrent()){
      _hotrep->constraintTrajectories(cons);
      return;
    } else 
      ErrMsg(error) << "Refit KalRep is not current!  Using MiniRep cached fits " << endmsg;
  }
// otherwise use the cached constraints.
  std::vector<TrkSimpTraj*>::const_iterator icon = _contrajs.begin();
  while(icon != _contrajs.end()){
    cons.push_back(*icon);
    icon++;
  }
}

TrkHotList*
KalMiniRep::hotList() {
// if we've fit the kalrep, use the full hot list
  if(_hotrep.get() != 0)
    return _fullhotlist.rawPtr();
  else {
// otherwise, return the empty first (if it exists)
    if(_emptyhotlist.get() != 0)
      return _emptyhotlist.get();
    else if(_fullhotlist.rawPtr() != 0)
      return _fullhotlist.rawPtr();
    else
      return 0;
  }
}

const TrkHotList*
KalMiniRep::hotList() const {
  return const_cast<KalMiniRep*>(this)->hotList();
}

void
KalMiniRep::extendIfNeeded(double fltlen) const {
  static const double tol(0.5); // 1/2 cm tolerance
  if(_autoExtend ) {
// see if we're already in range
    if (!validFlightLength(fltlen,tol)) {
// if we're in cache, try to upgrade to extendedcache
      if(currentState() == KalMiniRep::cache) {
// we're outside the range.  promote this fit to extend mode
// cast-off const
	KalMiniRep* ncthis = const_cast<KalMiniRep*>(this);
	TrkErrCode exerr = ncthis->createXKalRep();
	if(exerr.success() ) {
	  exerr = _xrep->fit();
	  if(exerr.success()){
// update the state
	    ncthis->setValid(_xrep->fitValid());
	    ncthis->setCurrent(_xrep->fitCurrent());
	    ncthis->_state = KalMiniRep::extendedcache;
	  }
	}
      }
// extend tracks in refit and extend mode
      TrkErrCode exerr = extendThrough(fltlen);
    }
  }
}


void
KalMiniRep::addHistory(const TrkErrCode& status,const char* modulename) {
  if(active())
    kalRep()->addHistory(status,modulename);
// call the base class too
  TrkFitStatus::addHistory(status,modulename);
}

bool
KalMiniRep::svtdEdx() const {
  if(active()){
    const KalSite* bcsite = kalRep()->findSite(KalSite::betaConsSite);
    return (bcsite != 0 && bcsite->isActive());
  } else
    return _usesvtdedx;
}

void 
KalMiniRep::setSvtdEdx(bool usesvtdedx) {
// if we have a full Kalman fit, we have to add the KalBetaCons site
  KalRep* krep = kalRep(hots);
  if(krep != 0){
    const KalSite* bcsite = krep->findSite(KalSite::betaConsSite);
    if(bcsite != 0){
// activate it if necessary
      if(bcsite->kalBetaCons()->isActive() != usesvtdedx){
        const_cast<KalBetaCons*>(bcsite->kalBetaCons())->setActivity(usesvtdedx);
        krep->setCurrent(false);
      }
    } else if(usesvtdedx){
// no existing KalBetaCons site: create one and add it to the fit
      std::vector<const SvtHitOnTrack*> svthots;
      fillHots(svthots,krep->hotList());
      KalBetaCons* kbc = new KalBetaCons(&krep->pieceTraj(),&krep->parentTrack()->bField(),svthots,
                            krep->particleType());
      krep->addSite(kbc);
    }
  }
// update cached flag
  _usesvtdedx = usesvtdedx;
}

void
KalMiniRep::fillHots(std::vector<const SvtHitOnTrack*>& svthots,const TrkHotList* hotlist) const {
  svthots.clear();
  TrkHotList::hot_iterator ihot = hotlist->begin();
  while(ihot != hotlist->end()){
    const TrkHitOnTrk* hot = ihot.get();
    if(hot->svtHitOnTrack() != 0)
      svthots.push_back(hot->svtHitOnTrack());
    ihot++;
  }
}

 


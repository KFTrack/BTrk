//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: KalStub.cc,v 1.35 2004/08/06 06:12:49 bartoldu Exp $
//
// Description:
//      class KalStub.
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Copyright Information:
//	Copyright (C) 1998	Lawrence Berkeley Laboratory
//
// Author List:
//      Dave Brown 9/6/98
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"

#include "KalmanTrack/KalStub.hh"
#include "KalmanTrack/KalRep.hh"
#include "KalmanTrack/KalMaterial.hh"
#include "KalmanTrack/KalBend.hh"
#include "TrkBase/TrkHitUse.hh"
#include "TrkBase/TrkVolume.hh"
#include "DetectorModel/DetSet.hh"
#include "ErrLogger/ErrLog.hh"
#include "BaBar/BbrCollectionUtils.hh"


KalStub::KalStub(const KalRep& seedrep,trkDirection extenddir,
		 double tolerance,double* range,
		 std::deque<KalSite*> sites,
		 const KalContext& kalcon) :
  _kalcon(kalcon),
  _kalrep(seedrep),
  _ptraj(_kalrep.referenceTraj()),
  _tdir(extenddir),
  _deltachi(0.0),
  _tolerance(tolerance),
  _sites(sites.begin(),sites.end()),
  _primary(true)
{
  _hitcache.reserve(16);
// set the range
  _xrange[0] = range[0];
  _xrange[1] = range[1];
// find the 'last' site of the rep
  _lastsite = lastSite(_kalrep);
}

KalStub::KalStub(const KalStub& other) :
  _kalcon(other._kalcon),
  _kalrep(other._kalrep),
  _lastsite(other._lastsite),
  _ptraj(other._ptraj),
  _tdir(other._tdir),
  _status(other._status),
  _deltachi(other._deltachi),
  _tolerance(other._tolerance),
  _primary(other._primary)
{
  _hitcache.reserve(other._hitcache.size());
// deep-copy the sites
  unsigned nsites = other._sites.size();
  for(unsigned isite=0;isite<nsites;isite++){
    KalSite* newsite = other._sites[isite]->clone(&_kalrep);
    assert(newsite != 0);
    _sites.push_back(newsite);
  }
// copy the range
  _xrange[0]= other._xrange[0];
  _xrange[1]= other._xrange[1];
}

KalStub::KalStub(KalStub& other,
		 const TrkHitUse* hitToAdd) :
  _kalcon(other._kalcon),
  _kalrep(other._kalrep),
  _ptraj(other._ptraj),
  _tdir(other._tdir),
  _status(other._status),
  _deltachi(other._deltachi),
  _tolerance(other._tolerance),
  _primary(other._primary)
{
  _hitcache.reserve(other._hitcache.size());
  // find the site corresponding to the Hot
  int startsite,endsite,sitestep;
  double epsilon(0.);
  switch(_tdir){
  case trkIn:
    startsite = 0;
    endsite = other._sites.size();
    sitestep = 1;
    epsilon = -0.0001;
    break;
  case trkOut:
    startsite = other._sites.size()-1;
    endsite = -1;
    sitestep = -1;
    epsilon = 0.0001;
    break;
  }

  int hitsite = -1;
  unsigned isite(startsite);
  for(isite=startsite;isite!=endsite;isite+=sitestep){
    if( other._sites[isite]->kalHit() == 0 ) continue;
    _lastsite = other._sites[isite]->kalHit();
    hitsite = isite;
    _primary = false;
    break;
  }
// fall-back to a regular KalStub if we couldn't find the hot
  if(hitsite<0){
    ErrMsg(debugging) << "Could not find requested hot in stub " << endmsg;
    _lastsite = other._lastsite;
    hitsite = 0;
  }

// make sure that the old KalStub is valid up to _lastsite
  other.process(_lastsite->globalLength() + epsilon);

// deep-copy the remaining non-hit sites
  for(isite=0;isite<other._sites.size();isite++){
    if(other._sites[isite]->kalHit() == 0){
      KalSite* newsite = other._sites[isite]->clone(&_kalrep);
      assert(newsite != 0);
      _sites.push_back(newsite);
    }
  }

// copy the range
  _xrange[0]= other._xrange[0];
  _xrange[1]= other._xrange[1];
  
  if ( hitToAdd ) {
    int addHitIndex = other.findHitInCache(*hitToAdd);
    if ( addHitIndex < 0 ) {
      ErrMsg(warning) << "Couldn't find hit to add to stub in cache" << endmsg;
    } else {
      _hitcache.push_back(KalHitCache(other._hitcache[addHitIndex],&_kalrep));
      TrkErrCode errC = addHit(*hitToAdd);
      if ( errC.failure() ) {
	ErrMsg(warning) << "Failed to add requested hit to sub-stub " << errC << endmsg;
      }
    }
  }
}

KalStub*
KalStub::clone(bool flushCache) const {
  KalStub* newstub = new KalStub(*this);
  //cout << "KalStub (copy) " << newstub << endl;
// if necessary, clone the cache
  if(!flushCache){
    unsigned ncache = _hitcache.size();
    for(unsigned icache=0;icache<ncache;icache++)
      newstub->_hitcache.push_back(KalHitCache(_hitcache[icache],&_kalrep));
  }
  return newstub;
}

KalStub::~KalStub() {
  // debug printout
  //  cout << "kill KalStub " << this << endl;
// clean out the cache
  unsigned ncache = _hitcache.size();
  for(unsigned icache=0;icache<ncache;icache++)
    _hitcache[icache].deleteAll();
// delete the HOTs
  unsigned nsites = _sites.size();
  for(unsigned isite=0;isite<nsites;isite++)
    if(_sites[isite]->kalHit() != 0)
      delete _sites[isite]->kalHit()->hitOnTrack();
// clean out the sites
  std::for_each(_sites.begin(),_sites.end(),babar::Collection::DeleteObject());
}

TrkErrCode
KalStub::chisquared(const TrkHitUse& hituse,double& chisq,const TrkHitOnTrk*& rethot) {
  rethot = 0;
  chisq = -1.0;
  if(!_status.success())return _status;
// process up to the hituse flightlength
  process(hituse.fltLen());
// see if the hituse exists in cache
  int index = findHitInCache(hituse);

  TrkHitOnTrk* hot(0);
  KalHit* kalhit(0);
  KalHitCache* cachedHit(0);
  if(index < 0){
// create a HOT for this hit usage
    hot = hituse.createHitOnTrk(_kalrep);
    if(hot == 0) {
      return TrkErrCode(TrkErrCode::fail,63,"Couldn't create HOT");      
    }    
    kalhit = new KalHit(_ptraj,hot);
    if (kalhit == 0 ) {
      delete hot;
      return TrkErrCode(TrkErrCode::fail,62,"Couldn't create KalHit");      
    }
  } else {
// hit is already in cache
// otherwise get the kal hit for the chi2 calculation and access to the hot
    cachedHit = &(_hitcache[index]);
    assert(cachedHit != 0);
    kalhit = cachedHit->kalHit();
    assert(kalhit != 0);
  }   
  if ( cachedHit && cachedHit->valid() ) {
// cache exists, and is valid, return it
    chisq = cachedHit->chisquared();      
  } else {
// cache doesn't exist, or isn't valid, so we (re)-compute the chisq
// make sure the hit is in range
    if(! inRange(kalhit->globalLength()) ){
// delete everything!
      kalhit->deleteHOT();	  
      delete kalhit;
      if ( cachedHit ) {
	_hitcache.erase(_hitcache.begin()+index);
      }
      return TrkErrCode(TrkErrCode::fail,64,"HitUse not in extension range");
    }
// find the site closest to this flight length
    const KalSite* prevsite = findSite(kalhit->globalLength());
// compute chisquared relative to this site
    kalhit->chisquared(chisq,prevsite,_tdir);
// if needed, build a cache for this item, and insert it in the cache list
    if ( cachedHit == 0 ) {
      _hitcache.push_back(KalHitCache(hituse,kalhit,chisq));
    } else {
// update the cache (which also sets it valid)
      cachedHit->setChisq(chisq);
    }
  }
// pass back the relevent hot and return
  rethot = kalhit->hitOnTrack();
  return _status;
}

int
KalStub::findHitInCache(const TrkHitUse& hit) const {
// see if this hit is already part of the cache
  unsigned ncache = _hitcache.size();
  for(unsigned icache=0;icache<ncache;icache++)
    if(_hitcache[icache] == hit)
      return icache;
  return -1;
}

const KalSite*
KalStub::findSite(double glen) const {
  int nsites = _sites.size();
  int jsite=-1;
  int isite;
  if(_tdir == trkOut){
    for(isite=0;isite<nsites;isite++)
      if(glen < _sites[isite]->globalLength())
	break;
    jsite = isite-1;
  } else {
    for(isite=nsites-1;isite>=0;isite--)
      if(glen > _sites[isite]->globalLength())
	break;
    jsite = isite+1;
  }
// return the appropriate site in the stub list
  if(jsite>=0 && jsite <= nsites-1)
    return _sites[jsite];
  else
// return the last/first site on the rep
    return _lastsite;
}

bool
KalStub::inRange(double fltlen) const {

  bool val(false);
  switch ( _tdir ) {
  case trkOut:
    val = fltlen-_xrange[0] > -_tolerance;
    break;
  case trkIn:
    val = fltlen-_xrange[1] < _tolerance;
    break;
  }

  if ( !val ) {
    ErrMsg(debugging) << "Candidate hit not in Range " << fltlen << ' ' << _xrange[0] << " - " << _xrange[1] << endmsg;
  }
  return val;
}

TrkErrCode
KalStub::addHit(const TrkHitUse& hituse,bool flushCache) {
  if(!_status.success())return _status;
  TrkErrCode status;
// see if the hituse exists in cache
  int index = findHitInCache(hituse);
  if(index >= 0){
// update the chisquared count
    _deltachi += _hitcache[index].chisquared();
// insert the hit (remove it from the cache first)
    KalHitCache hitc = _hitcache[index];
    _hitcache.erase(_hitcache.begin()+index);
    KalHit* khit = hitc.kalHit();
// if we're adding outpwards, push back, otherwards push forwards
    bool needssort(false);
    if(_tdir == trkOut){
// check explicitly for sorting
      needssort = _sites.size() != 0 &&
	*khit < *_sites.back();
      _sites.push_back(khit);
    } else {
      needssort = _sites.size() == 0 &&
	*_sites.front() < *khit;
      _sites.push_front(khit);
    }
// for now, always sort.  This class will need to be re-designed to 
// deal with stl efficiently.
    std::sort(_sites.begin(),_sites.end(),babar::Collection::PtrLess());
// delete or clear the cache
    unsigned ncache = _hitcache.size();
    for(unsigned icache=0;icache<ncache;icache++) {
      if ( flushCache ) { _hitcache[icache].deleteAll(); }
      else { _hitcache[icache].setValid(false); }
    }
    if ( flushCache ) _hitcache.clear();
// update the range
    unsigned idir = _tdir == trkIn ? 1 : 0;
    _xrange[idir] = khit->globalLength();
// flush the cache for any material/bend sites beyond the hit
// which may have been processed in evaluating chisquared for
// other hits
    double glen = khit->globalLength();
    int isite;
    int nsites=_sites.size();
    switch (_tdir) {
    case trkIn:
      for(isite=0;isite<nsites;isite++){
	if(_sites[isite]->globalLength()<glen && _sites[isite]->hasFit(_tdir))
	  _sites[isite]->invalidateSite(_tdir);
      }
      break;
    case trkOut:
      for(isite=nsites-1;isite>=0;isite--){
	if(_sites[isite]->globalLength()>glen && _sites[isite]->hasFit(_tdir))
	  _sites[isite]->invalidateSite(_tdir);
      }
      break;
    }
  } else
    status = TrkErrCode(TrkErrCode::fail,65,"Cannot add hit not already in cache");
  return status;
}


//  Process the sites
bool
KalStub::process(double fltlen) {
  int startsite,endsite,sitestep;
  switch(_tdir){
  case trkOut:
    startsite = 0;
    endsite = _sites.size()-1;
    sitestep = 1;
    break;
  case trkIn:
    startsite = _sites.size()-1;
    endsite = 0;
    sitestep = -1;
    break;
  }
//  start processing from the last site of the rep
  const KalSite* prevsite = _lastsite;
//  Filter through the sites and process them for this direction. The starting
//  parameters are the inflated-error intial parameters.  Only process sites
//  which need it.  Note that the processing latches to force processing any
//  sites appearing after a site which needed processing.
  bool status(true);
  bool needsfit(false);
  for(int isite=startsite;isite!=endsite+sitestep;isite+= sitestep){
    KalSite* thesite = _sites[isite];
    bool beyond = _tdir==trkOut ? thesite->globalLength() > fltlen :
      thesite->globalLength() < fltlen;
    if(beyond)break;
    needsfit = thesite->needsFit(_tdir) || needsfit;
    if(needsfit){
      if(!thesite->process(prevsite,_tdir)){
	status = false;
	break;
      }
    }
    prevsite = thesite;
  }
  return status;
}

const TrkHitOnTrk*
KalStub::findHOT(const TrkHitUse& hituse) const {
// find the hit in cache
  const TrkHitOnTrk* retval(0);
  int index = findHitInCache(hituse);
  if(index >= 0){
    retval = _hitcache[index].kalHit()->hitOnTrack();
  }
  return retval;
}

bool
KalStub::sameRep(const KalRep& other) const {
// require the KalRep to be the same _and_ not to
// have been altered since the KalStub was created
  const KalSite* otherlast = lastSite(other);
  return &_kalrep == &other &&
    _ptraj == other.referenceTraj() &&
    _lastsite == otherlast;
}

unsigned
KalStub::nDOF() const {
  unsigned nhits(0);
  unsigned nsites = _sites.size();
  for(unsigned isite=0;isite<nsites;isite++)
    if(_sites[isite]->kalHit() != 0)nhits++;
  return nhits;
}

const KalSite* 
KalStub::lastSite(const KalRep& rep) const {
  const KalSite* last = _tdir == trkOut ? rep.siteList().back() :
    rep.siteList().front();
  return last;
}

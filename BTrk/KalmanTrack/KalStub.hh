//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: KalStub.hh,v 1.23 2004/08/06 06:12:49 bartoldu Exp $
//
// Description:
//      class KalStub. This class is a pattern-recognition tool which allows
//      extending a (prefit) KalRep in a _single_ direction.  It allows testing
//      the effect of adding new hits to the rep in an efficient manner,
//      automatically adding material and bend sites as needed.  Once
//      the PatRec phase is over, this class can also be used to update the
//      original KalRep, or it can be passed to a KalMaker to construct a new track.
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Copyright Information:
//	Copyright (C) 1998	Lawrence Berkeley Laboratory
//
// Author List:
//      Dave Brown 9/6/98
//------------------------------------------------------------------------
#ifndef KALSTUB_HH
#define KALSTUB_HH

#include "BTrk/KalmanTrack/KalHit.hh"
#include "BTrk/TrkBase/TrkErrCode.hh"
#include "BTrk/TrkBase/TrkDirection.hh"
#include "BTrk/DetectorModel/DetIntersection.hh"
#include "BTrk/KalmanTrack/KalHitCache.hh"
#include <iostream>
#include <vector>
#include <deque>

class KalRep;
class KalContext;
class TrkHitUse;
class TrkVolume;

class KalStub {
public:
// copy constructor; does NOT copy cache
  KalStub(const KalStub& other);
// sub-stub constructor: may provide a single hit use to copy from cache
  KalStub(KalStub& other,
	  const TrkHitUse* hitToAdd);
// clone; optional cache copy
  virtual KalStub* clone(bool flushCache=true) const;
  virtual ~KalStub();
// The internal state of this object may not allow all possible operations: they
// user can determine which are possible by checking the status
  const TrkErrCode& status() const { return _status; }
// Flightlength range for which hits can be added to this stub
  double hitRange(trkDirection tdir) const { 
    return tdir == trkIn ? _xrange[0] : _xrange[1]; }
// access the rep
  const KalRep& kalRep() const { return _kalrep; }
// direction
  trkDirection direction() const { return _tdir; }
// compute the chisquared penalty for adding a hit onto the end.  This
// adds the hit to the consideration cache.
  TrkErrCode chisquared(const TrkHitUse& hituse,double& chisq,const TrkHitOnTrk*& hot);
//  Add the hit to the stub.  This uses a cache if the provided hit was
//  previously considered.
  TrkErrCode addHit(const TrkHitUse& hituse,bool flushCache=true);
// access the HOT behind a cache object
  const TrkHitOnTrk* findHOT(const TrkHitUse& hituse) const;
// stubs current delta-chisquared and degrees-of-freedom
  double deltaChisquared() const { return _deltachi; }
  unsigned nDOF() const;
  unsigned nSites() const { return _sites.size(); }
// access to the sites
  const std::deque<KalSite*>& sites() const { return _sites;}
// verify this stub against a KalRep
  bool sameRep(const KalRep&) const;
// inquire as to hit adding status
  unsigned nHitsAdded() const { return nDOF(); }
// is this a primary KalStub?
  bool isPrimary() const { return _primary; }
// config information
  const KalContext& kalContext() const { return _kalcon; }
private:
// only KalRep can call; this provides the material and bend sites
// construct from a KalRep.  Optionally override the KalContext
  KalStub(const KalRep& rep,trkDirection extenddir,
	  double tolerance,double* range,std::deque<KalSite*> sites,
	  const KalContext& kalcon);
// non-const access to the sites, again for KalRep
  std::deque<KalSite*>& sites() { return _sites;}
//
  const KalContext& _kalcon; // context for this rep
  const KalRep& _kalrep; // what we started with
  const KalSite* _lastsite; // last site of the rep
  const TrkDifPieceTraj* _ptraj; // trajectory of the rep
  trkDirection _tdir; // what direction we're going in
  TrkErrCode _status; // current status, defines what operations are legal
  double _xrange[2]; // range in which hits can be added
  double _deltachi; // current chisquared of stub
  double _tolerance; // pathlength tolerance for adding new hits on the 'end'
  std::deque<KalSite*> _sites; // sorted list of sites added to the stub
  std::vector<KalHitCache> _hitcache; // cache of hits under consideration
  bool _primary; // whether or not this stub was made directly from a rep
// utility functions
  int findHitInCache(const TrkHitUse& hituse) const;
  const KalSite* findSite(double len) const;
  const KalSite* lastSite(const KalRep& rep) const;
  bool inRange(double) const;
  bool process(double); // process the sites
// preempt
  KalStub& operator = (const KalStub& other);
// allow KalRep to get to private parts
  friend class KalRep;
};

#endif

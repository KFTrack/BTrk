//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: KalHitCache.hh,v 1.6 2002/07/15 20:31:49 echarles Exp $
//
// Description:
//      class KalHitCache.  Wrapper class to hold together pat. rec. hypos.
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Copyright Information:
//	Copyright (C) 1998	Lawrence Berkeley Laboratory
//
// Author List:
//      Dave Brown10/20/98
//------------------------------------------------------------------------
#ifndef KALHITCACHE_HH
#define KALHITCACHE_HH

class TrkHitUse;
class KalHit;
class KalRep;

class KalHitCache {
public:
  KalHitCache();
  KalHitCache(const TrkHitUse& hituse,KalHit*,double chisq);
// the next constructor _share hits_
  KalHitCache(const KalHitCache& other);
// create with a new KalHit
  KalHitCache(const KalHitCache& other,const KalRep* );
// this next function
  KalHitCache& operator = (const KalHitCache& other);
  virtual ~KalHitCache();
  KalHit* kalHit() { return _kalhit; }
  const KalHit* kalHit() const { return _kalhit; }
  double chisquared() const { return _chisquared; }
  bool valid() const { return _chi2Valid; }
  bool operator == (const KalHitCache& other) const;
  bool operator == (const TrkHitUse& hituse) const;

// operations that extend to the KalHit
  void deleteAll();
private:

  void setValid(bool isValid = false) { _chi2Valid = isValid; }
  void setChisq(double chisq) { _chisquared = chisq; _chi2Valid = true; }

  friend class KalStub;  

  KalHit* _kalhit;
  const TrkHitUse* _hituse; // has to be a pointer, to support nul constructor
  double _chisquared;
  bool _chi2Valid;
};

#endif

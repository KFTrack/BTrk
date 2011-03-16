//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: DchHitUse.hh 89 2010-01-14 12:34:14Z stroili $
//
// Description:
//     Minimal description of how a hit is to be used on a track; used for 
//     adding hits to tracks.
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author(s): Steve Schaffner
//
//------------------------------------------------------------------------
#ifndef DCHHITUSE_HH
#define DCHHITUSE_HH

#include "TrkBase/TrkHitUse.hh"
#include "TrkBase/TrkHitOnTrkUpdater.hh"

class DchHit;

// Class interface //
class DchHitUse : public TrkHitUse, TrkHitOnTrkUpdater {

public:
  DchHitUse(const DchHit&, double fltLen, int ambig,
            bool active=true, int usable=1);
  virtual ~DchHitUse();
  virtual bool operator==(const TrkHitUse&) const;

  int ambig() const                                          {return _ambig;}
  virtual TrkHitOnTrk* createHitOnTrk(const TrkRep&) const;

  // Cast (safe) of base hit
  const DchHit* dchHit() const;

private:
  int _ambig;
  double _t0;

  // Preempt 
  DchHitUse& operator= (const DchHitUse&);
  DchHitUse(const DchHitUse &);
};

#endif

//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: DchRecoHitOnTrack.hh 89 2010-01-14 12:34:14Z stroili $
//
// Description:
//   Contains drift chamber hit info, as hit is used on a particular track
//   Inherits from TrkHitOnTrk.  The drift distance is stored as an 
//   absolute value, but returned as |drift|*ambig.  Ambiguity stored as +/- 1; 
//   ambiguity = 0 => pick better value @ first call to updateMeasurement.
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Authors: Steve Schaffner
//------------------------------------------------------------------------

#ifndef DCHRECOHITONTRACK_HH
#define DCHRECOHITONTRACK_HH

#include "DchData/DchHitOnTrack.hh"

class DchLayer;
#include "CLHEP/Matrix/Vector.h"
class TrkHitOnTrk;
class TrkRep;
class TrkDifTraj;
class Trajectory;
class TrkRecoTrk;

// Class definition//
class DchRecoHitOnTrack : public DchHitOnTrack {
public:
  DchRecoHitOnTrack(const DchHit& baseHit, int ambig=0, double bunchTime=0,bool isProtoII=false);
  virtual ~DchRecoHitOnTrack();
  virtual TrkHitOnTrk* clone(TrkRep *, const TrkDifTraj *trkTraj=0) const;

protected:
  DchRecoHitOnTrack(const DchRecoHitOnTrack &hitToBeCopied, TrkRep *newRep,
                const TrkDifTraj* trkTraj=0 );

public:
  unsigned                      status()          const;
  const DchHit*                 dchHit()          const;
  unsigned                      tdcIndex()        const;

  int                           whichView()       const;

private:

  DchRecoHitOnTrack(const DchRecoHitOnTrack&);
  DchRecoHitOnTrack& operator=(const DchRecoHitOnTrack&);
};

#endif

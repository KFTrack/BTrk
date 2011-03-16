//  DchRecoHitOnTrack.cc
//
#include "BaBar/BaBar.hh"

#include "DchData/DchRecoHitOnTrack.hh"
#include "DchData/DchHit.hh"

// DchGeom needed to verify if hit is inside of chamber...
// and to find the trajectory describing the hit, i.e. wire
#include "DchData/DchHitOnTrack.hh"
#include "TrkBase/TrkEnums.hh"
#include "TrkBase/TrkSimpTraj.hh"
#include "ErrLogger/ErrLog.hh"


//-------------
// Constructors
//-------------

DchRecoHitOnTrack::DchRecoHitOnTrack(const DchHit &baseHit, int ambig, double bunchTime,bool isProtoII)
  : DchHitOnTrack(baseHit,baseHit, ambig,  bunchTime, isProtoII)
{ }


DchRecoHitOnTrack::DchRecoHitOnTrack(const DchRecoHitOnTrack &hot,
                                     TrkRep *newRep,
                                     const TrkDifTraj* trkTraj)
  : DchHitOnTrack(hot,newRep,trkTraj)
{ }

DchRecoHitOnTrack::~DchRecoHitOnTrack()
{ ; }

TrkHitOnTrk*
DchRecoHitOnTrack::clone(TrkRep *rep, const TrkDifTraj *trkTraj) const
{
  return new DchRecoHitOnTrack(*this,rep,trkTraj);
}

const DchHit*
DchRecoHitOnTrack::dchHit() const 
{
  return static_cast<const DchHit*>(hit());
}

unsigned
DchRecoHitOnTrack::status() const
{
  return dchHit()->status();
}

unsigned
DchRecoHitOnTrack::tdcIndex() const
{
  return dchHit()->tdcIndex();
}

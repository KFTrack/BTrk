//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: DchHitUse.cc 89 2010-01-14 12:34:14Z stroili $
//
// Description:
//     
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author(s): Steve Schaffner
//
//------------------------------------------------------------------------

#include "BaBar/BaBar.hh"
#include "DchData/DchHitUse.hh"
#include "DchData/DchRecoHitOnTrack.hh"
#include "DchData/DchHit.hh"


DchHitUse::DchHitUse(const DchHit& thehit, double flt, int amb,
                     bool active, int usable) :
  TrkHitUse(thehit, flt, active, usable) 
{
    _ambig = amb;
}

DchHitUse::~DchHitUse() 
{ }

TrkHitOnTrk* 
DchHitUse::createHitOnTrk(const TrkRep& rep) const 
{
  // This is a kludge -- fix once DchHitOnTrack gets more cleaned up
  //   (time should not be used for anything (I hope)
  // Also note that I temporarily have to cast off const from TrkRep
  const DchHit *h=dchHit(); assert(h!=0);
  DchRecoHitOnTrack tempHot(*h, ambig(), 0.);
  TrkHitOnTrk* newHot = tempHot.clone(&const_cast<TrkRep&>(rep));
  newHot->setFltLen( fltLen() );
  updateMeasurement(*newHot);
  return newHot;
}

bool 
DchHitUse::operator==(const TrkHitUse& rhs) const 
{
  // This is not going to win any design prizes:
  if (dchHit() == 0 || rhs.dchHit() ==0) return false;
  const DchHitUse& x = static_cast<const DchHitUse&>(rhs);
  return ( ambig() == x.ambig() && TrkHitUse::operator==(x) );
}

const DchHit* 
DchHitUse::dchHit() const 
{
  return static_cast<const DchHit*>(&(hit()));
}

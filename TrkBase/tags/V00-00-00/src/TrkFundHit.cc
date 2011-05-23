//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkFundHit.cc,v 1.19 2004/08/06 06:31:41 bartoldu Exp $
//
// Description:
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author(s): Steve Schaffner
//
//  We may want to inline some of these functions some day.
//
// Revision History:
//  20000523  M. Kelsey -- Add concrete printAll() implementation which
//		calls through to subclass print(), then dumps HOT list.
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"
#include "TrkBase/TrkFundHit.hh"
#include "TrkBase/TrkHitOnTrk.hh"
#include "TrkBase/TrkHitOnTrkIter.hh"
#include "TrkBase/TrkPredicates.hh"
#include <assert.h>
#include <algorithm>
#include <iostream>
// Remove!:
// MHK Restore, for use in printAll()
#include "TrkBase/TrkRecoTrk.hh"
using std::endl;
using std::ostream;


TrkFundHit::~TrkFundHit()
{
  // hitlist cleanup done in derived classes
}

TrkFundHit::TrkFundHit(const TrkFundHit& )
{
  assert(_hitList.empty());
}

TrkFundHit&
TrkFundHit::operator= (const TrkFundHit& x)
{
  assert(_hitList.empty()); 
  return *this;
}

const TrkHitOnTrk*
TrkFundHit::setUsedHit(const TrkHitOnTrk *hit)
{
  //  if (hitList->contains(hit)) {
  //    return;
  //  }
  //FIXME: check hot corresponds to this hit??
  _hitList.push_back(hit);
  return hit;
}

// return hit if on list, return 0 if not on list...
const TrkHitOnTrk *
TrkFundHit::setUnusedHit(const TrkHitOnTrk *hit)
{
  if (_hitList.empty()) return 0;
  std::vector<const TrkHitOnTrk*>::iterator i=std::find(_hitList.begin(),_hitList.end(),hit);
  if (i==_hitList.end()) return 0;
  assert(*i==hit);
  _hitList.erase(i);
  return hit;
}

int
TrkFundHit::nUsedHits() const
{
  return _hitList.size();
}

const TrkHitOnTrk*
TrkFundHit::getHitOnTrack(const TrkRecoTrk *trk) const
{
  hot_iterator i = std::find_if(begin(), end(),
                                std::bind2nd(TrkBase::Predicates::isHotOnTrack(),trk));
  return (i==end()?0:i.get());
}

void
TrkFundHit::printAll(ostream& os) const
{
  print(os);                            // Call through to get subclass info
  os << " used by " << nUsedHits() << " HOTs" << endl;
  if (usedHit()) {
    for (hot_iterator i=begin(); i != end(); ++i) {
      os << "\ttrack " << i->parentTrack()->id() << ": ";
      i->print(os);                // NOTE:  includes endl!
    }
    os << endl;
  }
}

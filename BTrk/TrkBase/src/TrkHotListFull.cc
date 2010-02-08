//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkHotListFull.cc,v 1.32 2005/02/11 16:27:11 raven Exp $
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
#include "TrkBase/TrkHotListFull.hh"
#include "TrkBase/TrkPredicates.hh"
#include "TrkBase/TrkHitOnTrk.hh"
#include "TrkBase/TrkView.hh"
#include "ErrLogger/ErrLog.hh"
#include "BbrStdUtils/BbrCollectionUtils.hh"
#include <assert.h>
#include <algorithm>
#include <functional>
#include <iterator>

// Default ctor
// Special case of ctor below:

TrkHotListFull::TrkHotListFull()
{
}


TrkHotListFull::TrkHotListFull(const TrkHotList& inHots,TrkBase::Functors::cloneHot f)
{
    _hotlist.reserve(dfltCapac());
//Clones Hots, and makes each point at the new track.
    std::transform(inHots.begin(),inHots.end(),std::back_inserter(_hotlist),f);
}


TrkHotListFull::TrkHotListFull(TrkHotList& inHots, TrkBase::Functors::setParent f)
{
    _hotlist.reserve(dfltCapac());
// shallow copy the hots and remove from input hotList
    std::transform(inHots.begin(),inHots.end(),std::back_inserter(_hotlist),f);
    inHots.hotlist().clear();
}


TrkHotList*
TrkHotListFull::clone(TrkBase::Functors::cloneHot f) const
{
  return new TrkHotListFull(*this, f);
}

TrkHotListFull::~TrkHotListFull()
{
// turn the parents off before deleting hots.  This avoids a cyclic delete error
// when deleting a track
//   std::for_each(begin(),end(),setParent(0));
  std::for_each(_hotlist.begin(),_hotlist.end(),babar::Collection::DeleteObject());
}

size_t
TrkHotListFull::dfltCapac()
{
  static size_t _dfltCapac = 75;
  return _dfltCapac;
}

void
TrkHotListFull::updateHots()
{
  std::for_each(begin(),end(),updateMeasurement(0,false));
  sort();
}

void
TrkHotListFull::append(TrkHitOnTrk* newHot)
{
  _hotlist.push_back(newHot);
}

void
TrkHotListFull::remove(TrkHitOnTrk* deadHot)
{
  typedef std::vector<TrkHitOnTrk*>::iterator iter_t;
  iter_t i = std::find(_hotlist.begin(),_hotlist.end(),deadHot);
  if (i!=_hotlist.end()) {
    delete *i;
    _hotlist.erase(i);
  } else
    ErrMsg(error) << " you asked to remove a hit which I don't have! " << endmsg;
}

TrkHitOnTrk*
TrkHotListFull::findHot(const TrkFundHit* theHit) const
{
  TrkBase::Predicates::hotMatchesFundHit match(theHit);
  TrkHotList::hot_iterator i = std::find_if(begin(),end(),match);
  return i==end()?0:const_cast<TrkHitOnTrk*>( i.get() );  // FIXME: return (non)const TrkHitOnTrk from (non)const
}

int
TrkHotListFull::nActive(TrkEnums::TrkViewInfo view) const
{
  int nAct = 0;
  for (TrkHotList::hot_iterator i = begin();i!=end();++i) {
    if (i->isActive())
      if(view == TrkEnums::bothView || i->whatView() == view)++nAct;
  }
  return nAct;
}

int
TrkHotListFull::nHit(TrkEnums::TrkViewInfo view) const
{
  if(view == TrkEnums::bothView)
    return end()-begin();
  else {
    int nAct = 0;
    for (TrkHotList::hot_iterator i = begin();i!=end();++i) {
      if(i->whatView() == view)++nAct;
    }
    return nAct;
  }
}

bool
TrkHotListFull::hitCapable() const
{
  return true;
}

int
TrkHotListFull::nDch(TrkEnums::TrkViewInfo view) const
{
  bool activeOnly(true);
  TrkBase::Predicates::isDchHitOnTrack dch(activeOnly);
  TrkBase::Predicates::hasView v(view);
// FIXME: W6U1 doesn't have std::count (at least, not in the config used by BaBar!)
// FIXME: std::compose2 is an SGI extension...
//  return std::count(begin(),end(),std::compose2(std::logical_and<bool>(),dch,v);
  int n = 0;
  for (TrkHotList::hot_iterator i = begin();i!=end();++i)
      if(dch(*i)&&v(*i)) ++n;
  return n;
}

int
TrkHotListFull::nSvt(TrkEnums::TrkViewInfo view) const
{
  bool activeOnly(true);
  TrkBase::Predicates::isSvtHitOnTrack svt(activeOnly);
  TrkBase::Predicates::hasView v(view);
// FIXME: W6U1 doesn't have std::count (at least, not in the config used by BaBar!)
// FIXME: std::compose2 is an SGI extension...
//  return std::count(begin(),end(),std::compose2(std::logical_and<bool>(),svt,v);
  int n = 0;
  for (TrkHotList::hot_iterator i = begin();i!=end();++i)
          if (svt(*i)&&v(*i)) ++n;
  return n;
}

double
TrkHotListFull::startFoundRange() const
{
  TrkBase::Predicates::isHotActive active;
  TrkHotList::hot_iterator i = std::find_if(begin(),end(),active);
  return i == end() ? 9999 : i->fltLen();
}

double
TrkHotListFull::endFoundRange() const
{
  TrkBase::Predicates::isHotActive predicate;
#if 1
  double maxFlt = -9999;
  TrkHotList::hot_iterator i = end();
  TrkHotList::hot_iterator b = begin();
  while (i-- != b) {
    if (predicate(*i)) {
      maxFlt = i->fltLen();
      break;
    }
  }
  return maxFlt;
#else
  typedef std::reverse_iterator<TrkHotList::hot_iterator> riter_t;
  riter_t rbegin = std::reverse_iterator(end());
  riter_t rend = std::reverse_iterator(begin());
  riter_t i = std::find_if(rbegin,rend,predicate);
  return i == rend ? -9999 : i->fltLen();
#endif

}


TrkView
TrkHotListFull::svtView(int layer) const
{
  TrkView retval;
  bool activeOnly(true);
  TrkBase::Predicates::isSvtHitOnTrack svt(activeOnly);
  TrkBase::Predicates::isLayer l(layer);
  // FIXME: std::compose2 is an SGI extension...
  for(TrkHotList::hot_iterator i=begin();i!=end();++i)
        if (svt(*i)&&l(*i)) retval.addView(i->whatView());
  return retval;
}

unsigned
TrkHotListFull::firstDchLayer() const
{
  unsigned firstlay(50);
  bool activeOnly(true);
  TrkBase::Predicates::isDchHitOnTrack dch(activeOnly);
  for (TrkHotList::hot_iterator i = begin();i!=end();++i) {
    if (dch(*i)) firstlay = std::min(firstlay,i->layerNumber());
  }
  return firstlay<50?firstlay:0;
}

unsigned
TrkHotListFull::lastDchLayer() const
{
  unsigned lastlay(0);
  bool activeOnly(true);
  TrkBase::Predicates::isDchHitOnTrack dch(activeOnly);
  for (TrkHotList::hot_iterator i = begin();i!=end();++i) {
    if (dch(*i)) lastlay = std::max(lastlay,i->layerNumber());
  }
  return lastlay;
}

const std::vector<TrkHitOnTrk*>&
TrkHotListFull::hotlist() const
{
  return _hotlist;
}

std::vector<TrkHitOnTrk*>&
TrkHotListFull::hotlist()
{
  return _hotlist;
}

bool
TrkHotListFull::isActive(unsigned ihot) const
{
  return ihot<_hotlist.size() && _hotlist[ihot]->isActive();
}

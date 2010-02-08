//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkHitList.cc,v 1.27 2003/02/19 21:42:21 raven Exp $
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
#include "TrkBase/TrkHitList.hh"
#include "TrkBase/TrkRecoTrk.hh"
#include "TrkBase/TrkHotList.hh"
#include "TrkBase/TrkHitOnTrk.hh"
#include "TrkBase/TrkHitOnTrkIter.hh"
#include "TrkBase/TrkErrCode.hh"
#include "TrkBase/TrkFundHit.hh"
#include "TrkBase/TrkRep.hh"
#include "TrkBase/TrkHitUse.hh"
#include "TrkBase/TrkRepIter.hh"
#include "ErrLogger/ErrLog.hh"


TrkHitList::TrkHitList(TrkRecoTrk* trk, PdtPid::PidType hypo)
        : _theTrack(trk)
        , _myHypo(hypo)
{
}

TrkHitList::~TrkHitList()
{
}

TrkRep*
TrkHitList::theRep()
{
  return _theTrack->getRep(_myHypo);
}

const TrkRep*
TrkHitList::theRep() const
{
  return _theTrack->getRep(_myHypo);
}

const TrkHotList&
TrkHitList::hotList() const
{
  const TrkHotList* x = theRep()->hotList();
  assert(x!=0);
  return *x;
}

TrkErrCode
TrkHitList::fit()
{
  _theTrack->_fitNumber[_myHypo]++;
  TrkErrCode err = theRep()->fit();
  theRep()->hotList()->sort();
  return err;
}

bool
TrkHitList::removeHit(const TrkFundHit *theHit)
{
  // This would be more efficient if the Rep did the finding of the Hot
  //   (save one search through each hotlist).
  if (theHit == 0) return false;

  if (!theHit->usedOnTrack(_theTrack)) {
    ErrMsg(warning) <<
      "TrkHitList: you just deleted a hit that was not on the track." << endmsg;
    return false;
  }
  std::pair<TrkRepIter,TrkRepIter> reps = _theTrack->uniqueReps();
  for (TrkRepIter i= reps.first; i != reps.second; ++i) {
    // Find the Hot and ask the Rep to remove it
    TrkHitOnTrk* h = i->hotList()->findHot(theHit);
    if (h != 0) i->removeHot(h);
  }
  return true;
}

TrkHitOnTrk*
TrkHitList::appendHot(const TrkHitOnTrk *theHot)
{
  if (theHot == 0) return 0;
  // Test whether hit already on this track
  if (theHot->hit() != 0 && theHot->hit()->usedOnTrack(_theTrack)) {
      ErrMsg(warning) << "You just tried to add a hit to a track twice.  Don't do that. "
                      << endmsg;
      return 0;
  }
  TrkHitOnTrk* defaultHot = 0;
  std::pair<TrkRepIter,TrkRepIter> reps = _theTrack->uniqueReps();
  for (TrkRepIter i=reps.first; i != reps.second; ++i) {
    TrkHitOnTrk* h = theHot->clone(i.get());
    i->addHot(h);
    if (i->particleType() == _theTrack->defaultType()) defaultHot = h;
  }
  return defaultHot;
}

TrkHitOnTrk*
TrkHitList::appendHit(const TrkHitUse& theHit)
{
  // Test whether hit already on this track
  if (theHit.hit().usedOnTrack(_theTrack)) {
    ErrMsg(warning) <<
      "You just tried to add a hit to a track twice.  Don't do that. " << endmsg;
    return 0;
  }
  TrkHitOnTrk* defaultHot = 0;
  std::pair<TrkRepIter,TrkRepIter> reps = _theTrack->uniqueReps();
  for (TrkRepIter i= reps.first; i != reps.second; ++i) {
    TrkHitOnTrk* h = theHit.createHitOnTrk(*i);
    i->addHot(h);
    if (i->particleType() == _theTrack->defaultType()) defaultHot = h;
  }
  return defaultHot;
}

bool
TrkHitList::append(const TrkHitList& list)
{
   bool x(true);
   for(TrkHitList::hot_iterator i = list.begin(); i!=list.end();++i) {
        TrkHitOnTrk* h = appendHot(i.get());
        x = ( x && h!=0);
   }
   return x;
}

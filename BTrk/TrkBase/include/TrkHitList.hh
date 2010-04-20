#ifndef TRKHITLIST_HH
#define TRKHITLIST_HH
//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkHitList.hh,v 1.16 2003/02/19 21:42:21 raven Exp $
//
// Description:
//     Provides track users with an interface to the track's list of hits.
//     It's really just a split-off piece of the TrkRecoTrk interface, and 
//     remains intimately connected to TrkRecoTrk.
//     Comments: 
//      Adding or removing a hit adds it to all hypotheses
//      Removing a hit destroys any Hots associated with it
//      Appending Hots clones them; the track does *not* take 
//             ownership of the hots passed in.
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author(s): Steve Schaffner
//------------------------------------------------------------------------

#include "BaBar/PdtPid.hh"
#include "TrkBase/TrkHitOnTrkIter.hh"
#include "TrkBase/TrkHotList.hh"

class TrkRecoTrk;
class TrkHitOnTrk;
class TrkFundHit;
class TrkRep;
class TrkErrCode;
class TrkHitUse;

// Class interface //
class TrkHitList {

public:
  typedef TrkHotList::hot_iterator hot_iterator;
  virtual ~TrkHitList();

  //**********************************
  // Access to the hits
  //**********************************
  unsigned      nHit()  const { return end()-begin();};
  hot_iterator  begin() const { return hotList().begin(); }
  hot_iterator  end()   const { return hotList().end(); }

  const TrkHotList& hotList() const;

  //**********************************
  // Modify
  //**********************************
  bool             removeHit(const TrkFundHit *theHit);  //ret. false if not found
  TrkHitOnTrk*     appendHot(const TrkHitOnTrk *theHot);
  TrkHitOnTrk*     appendHit(const TrkHitUse& theHit);
  bool             append(const TrkHitList& list);
  void             setActivity(const TrkHitOnTrk&);   

  //**********************************
  // Fit the track
  //**********************************
  // Perform the fit; should have no effect on an up-to-date (fitCurrent=true)
  //    fit 
  TrkErrCode          fit();

  bool operator==(const TrkHitList& right) const {return this == &right;}
private:
  // Data
  TrkRecoTrk* _theTrack;
  PdtPid::PidType _myHypo;

  // ctor
  TrkHitList(TrkRecoTrk*, PdtPid::PidType hypo);

  // Functions
  const TrkRep* theRep() const;
  TrkRep* theRep();

  // Preempt
  TrkHitList&   operator= (const TrkHitList&);
  TrkHitList(const TrkHitList &);

  friend class TrkRecoTrk;

};

#endif

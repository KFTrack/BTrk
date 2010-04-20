//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkFundHit.hh,v 1.23 2004/08/06 06:31:41 bartoldu Exp $
//
// Description:  Abstract base class.  Derived classes describe detector 
//  (Svt & Dch) hits, used as input to tracking routines.  All 
//  the base class does is provide list of pointers to HitOnTrack objects 
//  that are currently pointing at the underlying hit.  
//  Note: if you copy the contents 
//  of a FundHit into another (with either the copy ctor or operator=), 
//  the list of HOT pointers does not get copied (since the HOT points 
//  back to the original, not the copy).  
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Authors: Steve Schaffner
//
// Revision History:
//  20000523  M. Kelsey -- Add concrete printAll() implementation which
//		calls through to subclass print(), then dumps HOT list.
//------------------------------------------------------------------------

#ifndef TRKFUNDHIT_HH
#define TRKFUNDHIT_HH
 #include "TrkBase/TrkEnums.hh"
#include "TrkBase/TrkHitOnTrkIter.hh"

#include <vector>
#include <utility>
class TrkRecoTrk;
class TrkHitOnTrk;
class TrkDetElemId;
class GTrack;
#include <iosfwd>

// Class interface //
class TrkFundHit  {

public:
  typedef TrkHitOnTrkIter<TrkFundHit> hot_iterator;

  //**********************
  // Constructors and such
  //**********************
  TrkFundHit() {;}
  virtual ~TrkFundHit();

  //**********************
  // Access list of HitOnTrack objects this hit is associated with
  //**********************
  int nUsedHits() const;
  bool usedHit(void) const                       {return !_hitList.empty();}

  std::pair<TrkFundHit::hot_iterator,
            TrkFundHit::hot_iterator > getUsedHits() const {
        return std::pair<TrkFundHit::hot_iterator,TrkFundHit::hot_iterator >(begin(),end());
  }
  inline TrkFundHit::hot_iterator begin() const;
  inline TrkFundHit::hot_iterator end() const;

  // Is this hit used on track trk?
  bool usedOnTrack(const TrkRecoTrk *t) const {return getHitOnTrack(t) != 0;}
  // return HOT connecting this hit to track trk (return 0 if none)
  //  TrkHitOnTrk* getHitOnTrack(TrkRecoTrk *trk, PdtPid::PidType) const;
  const TrkHitOnTrk* getHitOnTrack(const TrkRecoTrk *trk) const;

  //**********************
  // Modify list of HitOnTrack objects
  //**********************
  const TrkHitOnTrk* setUsedHit(const TrkHitOnTrk *hit); // return hit if OK, 0 if not
  const TrkHitOnTrk* setUnusedHit(const TrkHitOnTrk *hit); // return hit if OK, 0 if not

  //**********************
  // Pattern-recognition functions
  //**********************
  virtual TrkEnums::TrkViewInfo whatView() const = 0;
  virtual TrkDetElemId elemId() const = 0;

  // MC truth (this may not survive until data-taking)
  virtual const GTrack* getGTrack() const = 0;

  //**********************
  // Dump list of HOTs (for debugging)
  //**********************
  virtual void printAll(std::ostream& os) const;

protected:
  friend class TrkHitOnTrkIter<TrkFundHit>;

  typedef std::vector<const TrkHitOnTrk*>::iterator  iterator_implementation;
  typedef const TrkHitOnTrk                    iterator_value_type;

  std::vector<const TrkHitOnTrk*> _hitList;
  // Operators
  TrkFundHit&   operator= (const TrkFundHit&);

private:
  // Copy ctor
  TrkFundHit(const TrkFundHit &);

};

// Might need this again someday:
  // Create a HitOnTrk object of the correct type
  //  virtual TrkHitOnTrk* makeHot(TrkRecoTrk *) = 0;


TrkFundHit::hot_iterator
TrkFundHit::begin() const
{
        return TrkFundHit::hot_iterator(const_cast<std::vector<const TrkHitOnTrk*>&>(_hitList).begin());
}

TrkFundHit::hot_iterator
TrkFundHit::end() const
{
        return TrkFundHit::hot_iterator(const_cast<std::vector<const TrkHitOnTrk*>&>(_hitList).end());
}

#endif

//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkHotList.hh,v 1.39 2004/08/06 06:31:42 bartoldu Exp $
//
// Description: List of hits (as HitOnTrk objects) associated with a 
//  reconstructed track.  Abstract base class.
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author(s): Steve Schaffner
//
//------------------------------------------------------------------------
#ifndef TRKHOTLIST_HH
#define TRKHOTLIST_HH

class TrkHitOnTrk;
class TrkRep;
class TrkView;
#include <iosfwd>
#include <vector>
#include "BTrk/TrkBase/TrkHitOnTrkIter.hh"
#include "BTrk/TrkBase/TrkEnums.hh"
#include "BTrk/TrkBase/TrkFunctors.hh"

// Class interface //
class TrkHotList  {
private:
  struct iterator_traits;
  struct const_iterator_traits;
public:

  // constructors and such
  TrkHotList();
  virtual TrkHotList* clone(TrkBase::Functors::cloneHot) const = 0;
  // this Hotlist is about to be usurped by a new TrkRep...
  virtual TrkHotList* resetParent(TrkBase::Functors::setParent);
  virtual ~TrkHotList();

  typedef TrkHitOnTrkIter<TrkHotList::const_iterator_traits> hot_iterator;
  hot_iterator  begin() const { return hot_iterator(hotlist().begin()); }
  hot_iterator  end() const   { return hot_iterator(hotlist().end()); }

  typedef TrkHitOnTrkIter<TrkHotList::iterator_traits> nc_hot_iterator;
  nc_hot_iterator  begin()    { return nc_hot_iterator(hotlist().begin()); }
  nc_hot_iterator  end()      { return nc_hot_iterator(hotlist().end()); }

  virtual bool         hitCapable()       const = 0;
  virtual int          nActive(TrkEnums::TrkViewInfo view=TrkEnums::bothView) const = 0;
  virtual int          nHit(TrkEnums::TrkViewInfo view=TrkEnums::bothView) const = 0;
  virtual double       startFoundRange()  const = 0;
  virtual double       endFoundRange()    const = 0;
  virtual bool         isActive(unsigned ihot) const =0;

  virtual void         append(TrkHitOnTrk* ) = 0;
  virtual void         remove(TrkHitOnTrk* ) = 0;
  virtual void         updateHots() = 0;
  virtual void         sort();

  void                 print(std::ostream &o) const;
  void                 printAll(std::ostream &o) const;

protected:
  friend struct iterator_traits;  // these two want access to hotlist_t...
  friend struct const_iterator_traits;
  typedef std::vector<TrkHitOnTrk*> hotlist_t;

private:
  // Preempt: Hots have to be cloned, not naively copied
  TrkHotList(const TrkHotList& rhs);
  TrkHotList&   operator= (const TrkHotList&);

  struct iterator_traits {
    typedef TrkHitOnTrk         iterator_value_type;
    typedef hotlist_t::iterator iterator_implementation;
  };
  struct const_iterator_traits {
    typedef const TrkHitOnTrk         iterator_value_type;
    typedef hotlist_t::const_iterator iterator_implementation;
  };

  friend class TrkHotListUnowned;  // these two want access to hotlist()
  friend class TrkHotListFull;
  virtual const hotlist_t&  hotlist() const = 0;
  virtual       hotlist_t&  hotlist() = 0;
};

#endif

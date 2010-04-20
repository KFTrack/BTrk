//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkHotListUnowned.hh,v 1.10 2003/10/25 21:58:49 brownd Exp $
//
// Description: class TrkHotListUnowned.  A copy of TrkHotListFull where the
// hots are not owned by this class.  This is useful in nested reps like the
//  mini-rep, where only 1 can truely own the hots but both may need access to
// them.
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
//// Copyright Information:
//	Copyright (C) 2000	Lawrence Berkeley Laboratory
//
// Author List:
//      Dave Brown 11/6/00
//------------------------------------------------------------------------

#ifndef TRKHOTLISTUNOWNED_HH
#define TRKHOTLISTUNOWNED_HH
#include "TrkBase/TrkHotList.hh"
#include "TrkBase/TrkHitOnTrkUpdater.hh"

class TrkFundHit;
class TrkHitOnTrk;
class TrkRep;

// Class interface //
class TrkHotListUnowned : public TrkHotList, TrkHitOnTrkUpdater {
public:
// create new HotList from another
  TrkHotListUnowned(TrkHotList* other,bool takehotlist=false);
// copy and assignment are OK for this class
  TrkHotListUnowned(const TrkHotListUnowned& rhs);
  TrkHotListUnowned&   operator= (const TrkHotListUnowned&);
// cloning
  TrkHotList* clone(TrkBase::Functors::cloneHot) const;
  virtual ~TrkHotListUnowned();

  virtual bool      hitCapable()      const;
  virtual int        nActive(TrkEnums::TrkViewInfo view=TrkEnums::bothView) const;
  virtual int        nDch(TrkEnums::TrkViewInfo view=TrkEnums::bothView) const;
  virtual int        nSvt(TrkEnums::TrkViewInfo view=TrkEnums::bothView) const;
  virtual int        nHit(TrkEnums::TrkViewInfo view=TrkEnums::bothView) const;
  virtual TrkView svtView(int layer) const;
  virtual unsigned  firstDchLayer() const;
  virtual unsigned  lastDchLayer() const;
  virtual double    startFoundRange() const;
  virtual double    endFoundRange()   const;
  virtual bool         isActive(unsigned ihot) const;

  void              append(TrkHitOnTrk* );
  void              remove(TrkHitOnTrk* );
  TrkHitOnTrk*      findHot(const TrkFundHit*) const;
  virtual void      updateHots();
// special function for this class
  const TrkHotList* myHotList() const { return _hotl; }
// take ownership
  TrkHotList* takeHotList();
protected:
  virtual const std::vector<TrkHitOnTrk*>&   hotlist()   const;
  virtual       std::vector<TrkHitOnTrk*>&   hotlist();
private:
// reference to the real HotList
  TrkHotList* _hotl;
  bool _ownhots; // do I own these or not?
};

#endif

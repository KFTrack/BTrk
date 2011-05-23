//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: DchDummyHit.hh 89 2010-01-14 12:34:14Z stroili $
//
// Description:
//      class DchDummyHit.  Dummy implementation of TrkFundHit, used with
//      DchMiniHitOnTrack.
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Authors:
//      16/10/00 R. Stroili
//
//      Copyright (C) INFN & Padova University
//------------------------------------------------------------------------

#ifndef DCHDUMMYHIT_HH
#define DCHDUMMYHIT_HH

#include "TrkBase/TrkFundHit.hh"
#include "TrkBase/TrkDetElemId.hh"

class DchDummyHit : public TrkFundHit {
public:
// pass in the view; this is the only function which this class can
// really implement of TrkFundHit (and it's completely irrelevant, however...)
  DchDummyHit(TrkEnums::TrkViewInfo view)
          : _view(view),_eid(0,TrkDetElemId::dch),_refcount(0)
  { addReference(); }
  virtual ~DchDummyHit(){}
// TrkFundHit required interface
  virtual TrkEnums::TrkViewInfo whatView() const { return _view; }
  virtual TrkDetElemId elemId() const { return _eid; }
  virtual const GTrack* getGTrack() const { return 0; } // deprecated function
  unsigned referenceCount() const { return _refcount; }
// functions to manipulate reference count
  void addReference() { ++_refcount; }
  void delReference() { --_refcount; }
private:
  TrkEnums::TrkViewInfo _view;
  TrkDetElemId _eid; // completely bogus element id.
  unsigned _refcount;

};

#endif

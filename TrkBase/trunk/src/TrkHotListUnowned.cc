//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkHotListUnowned.cc,v 1.10 2003/10/25 21:58:49 brownd Exp $
//
// Description: class TrkHotListUnowned.  A copy of TrkHotListUnowned where the
// hots are not owned by this class.  This is useful in nested reps like the
//  mini-rep.
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

#include "BaBar/BaBar.hh"
#include "TrkBase/TrkHotListUnowned.hh"
#include "TrkBase/TrkView.hh"
#include <assert.h>
#include "ErrLogger/ErrLog.hh"


TrkHotListUnowned::TrkHotListUnowned(TrkHotList* other,bool take) :
  _hotl(other),_ownhots(take)
{
}

TrkHotList*
TrkHotListUnowned::clone(TrkBase::Functors::cloneHot) const
{
  return new TrkHotListUnowned(_hotl);
}

TrkHotListUnowned::~TrkHotListUnowned()
{
  if(_ownhots)delete _hotl;
}

void
TrkHotListUnowned::append(TrkHitOnTrk* newHot)
{
  _hotl->append(newHot);
}

void
TrkHotListUnowned::remove(TrkHitOnTrk* deadHot)
{
  _hotl->remove(deadHot);
}

TrkHitOnTrk*
TrkHotListUnowned::findHot(const TrkFundHit* theHit) const
{
  return _hotl->findHot(theHit);
}

int
TrkHotListUnowned::nActive(TrkEnums::TrkViewInfo view) const
{
  return _hotl->nActive(view);
}

int
TrkHotListUnowned::nHit(TrkEnums::TrkViewInfo view) const
{
  return _hotl->nHit(view);
}


bool
TrkHotListUnowned::hitCapable() const
{
  return _hotl->hitCapable();
}

int
TrkHotListUnowned::nDch(TrkEnums::TrkViewInfo view) const
{
  return _hotl->nDch(view);
}

int
TrkHotListUnowned::nSvt(TrkEnums::TrkViewInfo view) const
{
  return _hotl->nSvt(view);
}

double
TrkHotListUnowned::startFoundRange() const
{
  return _hotl->startFoundRange();
}

double
TrkHotListUnowned::endFoundRange() const
{
  return _hotl->endFoundRange();
}

void
TrkHotListUnowned::updateHots()
{
  _hotl->updateHots();
}

TrkView
TrkHotListUnowned::svtView(int layer) const {
  return _hotl->svtView(layer);
}

unsigned
TrkHotListUnowned::firstDchLayer() const {
  return _hotl->firstDchLayer();
}

unsigned
TrkHotListUnowned::lastDchLayer() const
{
  return _hotl->lastDchLayer();
}

const std::vector<TrkHitOnTrk*>&
TrkHotListUnowned::hotlist() const
{
  return _hotl->hotlist();
}

std::vector<TrkHitOnTrk*>&
TrkHotListUnowned::hotlist()
{
  return _hotl->hotlist();
}

TrkHotList*
TrkHotListUnowned::takeHotList() {
  if(_ownhots){
    _ownhots = false;
    return _hotl;
  } else
    return 0;
}

bool
TrkHotListUnowned::isActive(unsigned ihot) const {
  return _hotl->isActive(ihot);
}


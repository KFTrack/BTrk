//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkView.cc,v 1.1 2002/10/26 00:17:47 brownd Exp $
//
//  Description:
//  Class TrkView.
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Copyright Infomation;
//	Copyright (C) 2002	Lawrence Berkeley Laboratory
//
// Author(s): Dave Brown 10/25/02
//
//------------------------------------------------------------------------

#include "BaBar/BaBar.hh"
#include "TrkBase/TrkView.hh"

const unsigned short
TrkView::_xyview(0x1);
const unsigned short
TrkView::_zview(0x2);
const unsigned short
TrkView::_bothview(0x3);

TrkView::TrkView(TrkEnums::TrkViewInfo view) :
  _view(0)
{
  addView(view);
}

TrkView::TrkView(unsigned short& pat) :
  _view( (pat & _bothview))
{
}

TrkView::TrkView(const TrkView& other) : _view(other._view)
{}

TrkView& 
TrkView::operator =(const TrkView& other) {
  if(this != &other){
    _view = other._view;
  }
  return *this;
}

TrkView::~TrkView()
{}

void 
TrkView::addView(TrkEnums::TrkViewInfo view) {
  switch (view) {
  case TrkEnums::noView: default:
    break;
  case TrkEnums::xyView:
    _view |= _xyview;
    break;
  case TrkEnums::zView:
    _view |= _zview;
    break;
  case TrkEnums::bothView:
    _view |= _bothview;
    break;
  }
}

TrkEnums::TrkViewInfo 
TrkView::view() const {
  switch (_view) {
  case 0: default:
    return TrkEnums::noView;
  case _xyview:
    return TrkEnums::xyView;
  case _zview:
    return TrkEnums::zView;
  case _bothview:
    return TrkEnums::bothView;
  } 
}

bool
TrkView::contains(const TrkView& other) const {
  return (other.viewData() & _view) == other.viewData();
}

bool
TrkView::operator ==(const TrkView& other) const {
  return other.viewData() == _view;
}

bool
TrkView::operator != (const TrkView& other) const {
  return other.viewData() !=  _view;
}

bool
TrkView::contains(TrkEnums::TrkViewInfo view ) const {
  TrkView other(view);
  return contains(other);
}


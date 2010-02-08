// ------------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkStoreHypo.cc,v 1.2 2002/10/23 21:12:16 brownd Exp $
//
//  Description: TrkStoreHypo
//
// Copyright Information:
//	Copyright (C) 2000	Lawrence Berkeley Laboratory
//
//  Authors: Dave Brown, 11/20/00
//------------------------------------------------------------------------------

#include "BaBar/BaBar.hh"
#include "TrkBase/TrkStoreHypo.hh"
#include <math.h>

double
TrkStoreHypo::_flttol(1.0); // call requests within this equivalent

bool 
TrkStoreHypo::operator == (const TrkStoreHypo& other) const {
  return _hypo == other._hypo && fabs(_fltlen-other._fltlen) < _flttol;
}
 

bool
TrkStoreHypo::operator < (const TrkStoreHypo& other) const {
  if(!operator ==(other)){
    if(_hypo != other._hypo)
      return _hypo < other._hypo;
    else
      return _fltlen < other._fltlen;
  } else
    return false;
}

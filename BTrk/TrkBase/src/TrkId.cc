//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkId.cc,v 1.9 2002/03/27 22:02:57 steinke Exp $
//
// Description:
//     Implementation of TrkId class.  Not a whole lot here.
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Authors:
//
// Copyright (C)  1996  The Board of Trustees of  
// The Leland Stanford Junior University.  All Rights Reserved.
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"
#include "TrkBase/TrkId.hh"
#include <assert.h>

unsigned TrkId::_maxval(0);

// Ctors
//------------------------------------------------------------------------
TrkId::TrkId(long myval) : _value(myval) {
//------------------------------------------------------------------------
}

//------------------------------------------------------------------------
TrkId::TrkId() : _value(++_maxval)  {
//------------------------------------------------------------------------
}

//------------------------------------------------------------------------
TrkId::~TrkId() { 
//------------------------------------------------------------------------
}

// Copy ctor
//------------------------------------------------------------------------
TrkId::TrkId(const TrkId &rhs) {
//------------------------------------------------------------------------
  _value = rhs._value;
}

//------------------------------------------------------------------------
TrkId& 
TrkId::operator= (const TrkId& rhs) {
//------------------------------------------------------------------------
  _value = rhs._value;
  return *this;
}

//------------------------------------------------------------------------
bool
TrkId::operator<(const TrkId& other) const {
  return _value < other._value;
}
//------------------------------------------------------------------------


//------------------------------------------------------------------------
void 
TrkId::setNewValue(const TrkId& source) {
//------------------------------------------------------------------------
  _value = source._value;
}

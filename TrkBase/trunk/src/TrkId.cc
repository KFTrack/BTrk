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
#include "TrkBase/TrkIdManager.hh"


// Ctors
//------------------------------------------------------------------------
TrkId::TrkId(long myval, TrkIdManager* man) : _value(myval), _idman(man) {
//------------------------------------------------------------------------
}

//------------------------------------------------------------------------
TrkId::TrkId(TrkIdManager* man) : _value(man->nextId()), _idman(man) {
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
  _idman = rhs.idManager();
  _value = rhs._value;
}

//------------------------------------------------------------------------
TrkId& 
TrkId::operator= (const TrkId& rhs) {
//------------------------------------------------------------------------
  _idman = rhs.idManager();
  _value = rhs._value;
  return *this;
}

//------------------------------------------------------------------------
bool
TrkId::operator<(const TrkId& other) const {
  if (*idManager() == *(other.idManager()) && _value < other._value) {
    return true;
  }
  else {
    if (idManager() < other.idManager()) {
      return true;
    }
    else
      return false;
  }
}
//------------------------------------------------------------------------


//------------------------------------------------------------------------
void 
TrkId::setNewValue(const TrkId& source) {
//------------------------------------------------------------------------
  _idman = source.idManager();
  assert (_idman != 0);
  _value = idManager()->nextId();
}

//------------------------------------------------------------------------
TrkIdManager* 
TrkId::idManager() const {
//------------------------------------------------------------------------
  return _idman;
}

//------------------------------------------------------------------------
void 
TrkId::setIdManager(TrkIdManager* idMan){
//------------------------------------------------------------------------
  _idman = idMan;
}

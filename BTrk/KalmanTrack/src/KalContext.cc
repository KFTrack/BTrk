//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: KalContext.cc 103 2010-01-15 12:12:27Z stroili $
//
// Description:
//      class KalContext
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Copyright Information:
//	Copyright (C) 1997	Lawrence Berkeley Laboratory
//
// Author List:
//      Dave Brown 3/15/97
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"
#include "KalmanTrack/KalContext.hh"
#include "BField/BFieldIntegrator.hh"
#include <assert.h>

KalContext::KalContext() : _bint(0), _trkmodel(0)
{}

KalContext::~KalContext(){
  delete _bint;
}

BFieldIntegrator const&
KalContext::bFieldIntegrator() const {
  if(_bint == 0){
    _bint = new BFieldIntegrator(bField(),_bintconfig);
    assert(_bint != 0);
  }
  return *_bint;
}

//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkExchangeData.cc,v 1.2 2004/07/05 05:43:26 hulsberg Exp $
//
//  Description:
//  Class TrkExchangeData;
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Copyright Infomation;
//	Copyright (C) 2000	Lawrence Berkeley Laboratory
//
// Author(s): Dave Brown
//
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"
#include "TrkBase/TrkExchangeData.hh"
#include "CommonUtils/ComPackFlatFloat.hh"
#include "CommonUtils/ComPackSignedExpFloat.hh"
#include "CommonUtils/ComPackBase.hh"
#include "BaBar/Constants.hh"
#include <assert.h>
#include <math.h>
#define BITSPERBYTE 8

// packing statics.  Note that these define the meaning of the data
// and should NEVER BE CHANGED.  Or, if they must be changed, this
// class will need to be subclassed and any persistent data based on it
// have its schema 'migrated'.

const ComPackSignedExpFloat TrkExchangeData::_packd0(16,6,128.0);
const ComPackFlatFloat TrkExchangeData::_packphi0(-Constants::pi,Constants::pi,16);
const ComPackSignedExpFloat TrkExchangeData::_packomega(16,6,0.5);
const ComPackSignedExpFloat TrkExchangeData::_packz0(16,6,256.0);
const ComPackFlatFloat TrkExchangeData::_packlambda(-Constants::pi/2.0,Constants::pi/2.0,16);

const ComPackBase<double>&
TrkExchangeData::paramPacker(int index) {
  switch(index) {
  case TrkExchangePar::ex_d0:
    return _packd0;
  case TrkExchangePar::ex_phi0:
    return _packphi0;
  case TrkExchangePar::ex_omega:
    return _packomega;
  case TrkExchangePar::ex_z0:
    return _packz0;
  case TrkExchangePar::ex_tanDip:
    return _packlambda;
  default:
    assert(0);
    return _packd0;
  }
}

TrkExchangeData::TrkExchangeData() {
  for(unsigned ipar=0;ipar<TrkExchangePar::nParam;ipar++)
    _params[ipar] = 0;
}

TrkExchangeData::TrkExchangeData(const TrkExchangePar* traj) {
// pack the parameters and the diagonal errors
  assert(traj != 0);
  for(int ipar=TrkExchangePar::ex_d0;ipar<TrkExchangePar::nParam;ipar++){
    unsigned parpack(0);
    const ComPackBase<double>& packpar = paramPacker(ipar);
// convert from tandip to dip
    double pval(0) ;
    switch (ipar) {
    case TrkExchangePar::ex_d0:
      pval = traj->d0();
      break;
    case TrkExchangePar::ex_phi0:
      pval = traj->phi0();
      break;
    case TrkExchangePar::ex_omega:
      pval = traj->omega();
      break;
    case TrkExchangePar::ex_z0:
      pval = traj->z0();
      break;
    case TrkExchangePar::ex_tanDip:
      pval = atan(traj->tanDip());
      break;
    }
    ComPackBaseBase::StatusCode pcode = packpar.pack(pval,parpack);
    assert(pcode !=ComPackBaseBase::TAG_BAD);
    _params[ipar] = parpack;
  }
}

TrkExchangeData::TrkExchangeData(const TrkExchangeData& other){
  *this = other;
}

TrkExchangeData& 
TrkExchangeData::operator =(const TrkExchangeData& other) {
  if(this != &other){
    for (int ipar=TrkExchangePar::ex_d0;ipar<TrkExchangePar::nParam;ipar++)
      _params[ipar] = other._params[ipar];
  }
  return *this;
}

TrkExchangeData::~TrkExchangeData(){}

TrkExchangePar*
TrkExchangeData::exchange() const {
// unpack the parameters and the errors
  HepVector pvec(5,0);
  for(int ipar=TrkExchangePar::ex_d0;ipar<TrkExchangePar::nParam;ipar++){
    const ComPackBase<double>& packpar = paramPacker(ipar);
    ComPackBaseBase::StatusCode pcode = packpar.unpack(_params[ipar],pvec[ipar]);
    assert(pcode != ComPackBaseBase::TAG_BAD );
  }
// convert back from dip to tandip
  pvec[TrkExchangePar::ex_tanDip] = tan(pvec[TrkExchangePar::ex_tanDip]);
  return new TrkExchangePar(pvec);
}

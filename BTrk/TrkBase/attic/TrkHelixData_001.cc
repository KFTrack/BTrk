//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkHelixData_001.cc,v 1.5 2004/07/05 05:43:27 hulsberg Exp $
//
//  Description:
//  Class TrkHelixData_001;
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
#include "TrkBase/TrkHelixData_001.hh"
#include "TrkBase/TrkMomCalculator.hh"
#include "CommonUtils/ComPackFlatFloat.hh"
#include "CommonUtils/ComPackInt.hh"
#include "CommonUtils/ComPackSignedExpFloat.hh"
#include "CommonUtils/ComPackBase.hh"
#include "BaBar/Constants.hh"
#include "ErrLogger/ErrLog.hh"
#include "BField/BField.hh"
#include "CLHEP/Vector/ThreeVector.h"
#include <assert.h>
#include <math.h>
#define BITSPERBYTE 8

// packing statics.  Note that these define the meaning of the data
// and should NEVER BE CHANGED.  Or, if they must be changed, this
// class will need to be subclassed and any persistent data based on it
// have its schema 'migrated'.

const ComPackSignedExpFloat TrkHelixData_001::_packd0(18,6,128.0);
const ComPackFlatFloat TrkHelixData_001::_packphi0(-Constants::pi,Constants::pi,18);
const ComPackSignedExpFloat TrkHelixData_001::_packomega(18,6,0.5);
const ComPackSignedExpFloat TrkHelixData_001::_packz0(18,6,256.0);
const ComPackFlatFloat TrkHelixData_001::_packlambda(-Constants::pi/2.0,Constants::pi/2.0,18);

const ComPackExpFloat TrkHelixData_001::_packd0err(14,3,2.0e-4,4.0);
const ComPackExpFloat TrkHelixData_001::_packphi0err(14,2,2.0e-5,Constants::pi/8.0);
const ComPackExpFloat TrkHelixData_001::_packomegaerr(14,6,2.0e-6,0.0625);
const ComPackExpFloat TrkHelixData_001::_packz0err(14,5,2.0e-4,16.0);
const ComPackExpFloat TrkHelixData_001::_packlambdaerr(14,2,2.0e-5,Constants::pi/8.0);

const ComPackExpFloat TrkHelixData_001::_packcorrfine(9,3,1.0,0.0);
const ComPackExpFloat TrkHelixData_001::_packcorrcoarse(8,3,0.0,1.0);
const ComPackFlatFloat TrkHelixData_001::_packcorrflat(0.0,1.0,8);
const ComPackInt TrkHelixData_001::_packsigns(0,(1<<10)-1);

const ComPackSignedExpFloat TrkHelixData_001::_packfltlen(16,3,256.0);
const ComPackFlatFloat TrkHelixData_001::_packfltrange(0.0,256.0,16);

const double TrkHelixData_001::_rescale(0.001); // rescale off-diags
const double TrkHelixData_001::_mindet(1.0e-8);// minimum unpacked determinant value

const ComPackBase<double>&
TrkHelixData_001::errorPacker(HelixTraj::ParIndex index) {
  switch(index) {
  case HelixTraj::d0Index:
    return _packd0err;
  case HelixTraj::phi0Index:
    return _packphi0err;
  case HelixTraj::omegaIndex:
    return _packomegaerr;
  case HelixTraj::z0Index:
    return _packz0err;
  case HelixTraj::tanDipIndex:
    return _packlambdaerr;
  default:
    assert(0);
    return _packd0err;
  }
}

const ComPackBase<double>&
TrkHelixData_001::paramPacker(HelixTraj::ParIndex index) {
  switch(index) {
  case HelixTraj::d0Index:
    return _packd0;
  case HelixTraj::phi0Index:
    return _packphi0;
  case HelixTraj::omegaIndex:
    return _packomega;
  case HelixTraj::z0Index:
    return _packz0;
  case HelixTraj::tanDipIndex:
    return _packlambda;
  default:
    assert(0);
    return _packd0;
  }
}

void
TrkHelixData_001::corrTerm(unsigned iterm,
		       HelixTraj::ParIndex& ipar,
		       HelixTraj::ParIndex& jpar) {
// this function defines the order in which correlation
// matrix terms are packed int the 96 bits.  First we pack
// the coarse terms, then the fine.
  switch (iterm) {
  case 0:
    jpar = HelixTraj::d0Index;
    ipar = HelixTraj::phi0Index;
    break;
  case 1:
    jpar = HelixTraj::z0Index;
    ipar = HelixTraj::tanDipIndex;
    break;
  case 2:
    jpar = HelixTraj::d0Index;
    ipar = HelixTraj::z0Index;
    break;
  case 3:
    jpar = HelixTraj::d0Index;
    ipar = HelixTraj::tanDipIndex;
    break;
  case 4:
    jpar = HelixTraj::phi0Index;
    ipar = HelixTraj::z0Index;
    break;
  case 5:
    jpar = HelixTraj::phi0Index;
    ipar = HelixTraj::tanDipIndex;
    break;
  case 6:
    jpar = HelixTraj::omegaIndex;
    ipar = HelixTraj::z0Index;
    break;
  case 7:
    jpar = HelixTraj::omegaIndex;
    ipar = HelixTraj::tanDipIndex;
    break;
  case 8:
    jpar = HelixTraj::d0Index;
    ipar = HelixTraj::omegaIndex;
    break;
  case 9:
    jpar = HelixTraj::phi0Index;
    ipar = HelixTraj::omegaIndex;
    break;
  }
}

const ComPackBase<double>&
TrkHelixData_001::corrPacker(HelixTraj::ParIndex ipar,
			 HelixTraj::ParIndex jpar) {
// i assume jpar<ipar
  assert(jpar<ipar);
  switch(jpar) {
  case HelixTraj::d0Index:
    switch(ipar) {
    case HelixTraj::phi0Index:
      return _packcorrfine;
    case HelixTraj::z0Index:  
    case HelixTraj::tanDipIndex:
      return _packcorrcoarse;
    case HelixTraj::omegaIndex:
      return _packcorrflat; 
    case HelixTraj::d0Index:
      assert(0);
    }
  case HelixTraj::phi0Index:
    switch(ipar) {
    case HelixTraj::omegaIndex:
      return _packcorrflat;
    case HelixTraj::z0Index:  
    case HelixTraj::tanDipIndex:
      return _packcorrcoarse; 
    case HelixTraj::phi0Index:
    case HelixTraj::d0Index:
      assert(0);
    }     
  case HelixTraj::omegaIndex:
    return _packcorrflat;
  case HelixTraj::z0Index:
    return _packcorrfine;
  default:
    assert(0);
    return _packd0;
  }
}

TrkHelixData_001::TrkHelixData_001() : _fltlen(0) {
  for(unsigned ipar=0;ipar<HelixTraj::NHLXPRM;ipar++)
    _params[ipar] = 0;
  for(unsigned iword=0;iword<3;iword++)
    _corr[iword] = 0;
}

TrkHelixData_001::TrkHelixData_001(const d_ULong params[5], const d_ULong corr[3],const d_ULong& fltlen)
{
  for (unsigned i=0; i<5; i++) _params[i]= params[i];
  for (unsigned i=0; i<3; i++) _corr[i]= corr[i];
  _fltlen= fltlen;
}

TrkHelixData_001::TrkHelixData_001(const HelixTraj* traj) {
// pack the parameters and the diagonal errors
  assert(traj != 0);
  const TrkParams* params= traj->parameters();
  assert(params != 0);
  const HepVector& pvec = params->parameter();
  HepVector perr(HelixTraj::NHLXPRM);
  const HepSymMatrix& pcov = params->covariance();
  unsigned ipar,jpar;
  for(ipar=HelixTraj::d0Index;ipar<HelixTraj::NHLXPRM;ipar++){
    unsigned parpack(0);
    unsigned errpack(0);
    const ComPackBase<double>& packpar = paramPacker((HelixTraj::ParIndex)ipar);
    const ComPackBase<double>& packerr = errorPacker((HelixTraj::ParIndex)ipar);
// convert from tandip to dip
    double pval =  (ipar != HelixTraj::tanDipIndex) ? pvec[ipar] : atan(pvec[ipar]);
    perr[ipar] = sqrt(max(pcov.fast(ipar+1,ipar+1),1.0e-30));
    ComPackBaseBase::StatusCode pcode = packpar.pack(pval,parpack);
    ComPackBaseBase::StatusCode perrcode =packerr.pack(perr[ipar],errpack);
    assert(pcode !=ComPackBaseBase::TAG_BAD &&
	   perrcode!=ComPackBaseBase::TAG_BAD);
    _params[ipar] = parpack | (errpack<<packpar.bitRange());
  }
// create the correlation matrix
  HepSymMatrix corr(5,1);
  for(ipar=HelixTraj::d0Index;ipar<HelixTraj::NHLXPRM;ipar++){
    for(jpar=HelixTraj::d0Index;jpar<ipar;jpar++){
      corr.fast(ipar+1,jpar+1) = pcov.fast(ipar+1,jpar+1)/(perr[ipar]*perr[jpar]);
    }
  }
  unsigned signs;
  double newdet;
  double rescale(_rescale);
  unsigned iused[3];
// pack the correlation matrix.  Loop till the unpacked determinant is OK
// note that the signs of the correlation terms are packed separately from their
// values, to avoid splitting data across longword boundaries
  do {
    unsigned iword(0);
// leave space for the signs in word 0.
    iused[0] =  _packsigns.bitRange();
    iused[1] = iused[2] = 0;
    _corr[0] = _corr[1] = _corr[2] = 0;
    signs = 0;
    for(unsigned icor=0;icor<10;icor++){
      HelixTraj::ParIndex ipar,jpar;
      corrTerm(icor,ipar,jpar);
// record the signs
      signs |= ((corr.fast(ipar+1,jpar+1)<0)<<icor);
      const ComPackBase<double>& packcorr = corrPacker(ipar,jpar);
      unsigned corrpack(0);
      ComPackBaseBase::StatusCode ccode = packcorr.pack(fabs(corr.fast(ipar+1,jpar+1)),corrpack);
      assert(ccode != ComPackBaseBase::TAG_BAD);
// advance to the next word if there's not enough space left in this one
      if(iused[iword] + packcorr.bitRange()>sizeof(d_ULong)*BITSPERBYTE){
	iword++;
	assert(iword<3);
      }
      _corr[iword] |= (corrpack<<iused[iword]);
// increment the used bit count
      iused[iword] += packcorr.bitRange();
    }
// pack the signs
    unsigned packsigns;
    ComPackBaseBase::StatusCode scode = _packsigns.pack(signs,packsigns);
    assert(scode != ComPackBase<int>::TAG_BAD);
    _corr[0] |= packsigns;
// unpack correlation matrix to test roundoff
    HepSymMatrix newcor;
    unsigned jused[3];
    correlation(newcor,jused);
    newdet = newcor.determinant();
// if the unpacked determininant is too small,
//  Reduce the off-diagonal terms to improve its mathematical behavior
    if(newdet<_mindet){
      double scale = max(1.0-rescale,0.0);
      for(unsigned ipar=HelixTraj::d0Index+1;ipar<=HelixTraj::NHLXPRM;ipar++)
	for(unsigned jpar=HelixTraj::d0Index+1;jpar<ipar;jpar++)
	  corr.fast(ipar,jpar) *= scale;
// next iteration, increase rescaling.  This speeds up convergence
      rescale *= 2.0;
    }
  } while (newdet<_mindet && rescale<2.0);
// pack the flightlength/range and flag
  double fltlen,fltrange;
  fltlen = traj->lowRange();
  fltrange = traj->hiRange()-fltlen;
  d_ULong fltlenpack,fltrangepack;
  ComPackBaseBase::StatusCode fltcode = _packfltlen.pack(fltlen,fltlenpack);
  ComPackBaseBase::StatusCode rngcode = _packfltrange.pack(fltrange,fltrangepack);

  assert(fltcode != ComPackBaseBase::TAG_BAD &&
	 rngcode != ComPackBaseBase::TAG_BAD);
  _fltlen =  fltlenpack | ( fltrangepack<<(_packfltlen.bitRange()));
}

TrkHelixData_001::TrkHelixData_001(const TrkHelixData_001& other){
  *this = other;
}

TrkHelixData_001& 
TrkHelixData_001::operator =(const TrkHelixData_001& other) {
  if(this != &other){
    for (unsigned ipar=0;ipar<HelixTraj::NHLXPRM;ipar++)
      _params[ipar] = other._params[ipar];
    for (unsigned icor=0;icor<3;icor++)
      _corr[icor] = other._corr[icor];
    _fltlen = other._fltlen;
  }
  return *this;
}

TrkHelixData_001::~TrkHelixData_001(){}

HelixTraj*
TrkHelixData_001::helix() const {
// unpack flightrange
  double range[2];
  flightRange(range);
// unpack the correlation matrix
  HepSymMatrix corr;
  unsigned iused[3];
  correlation(corr,iused);
// unpack the parameters and the errors
  HepVector pvec(5,0);
  HepVector perr(5,0);
  unsigned ipar,jpar;
  for(ipar=HelixTraj::d0Index;ipar<HelixTraj::NHLXPRM;ipar++){
    const ComPackBase<double>& packpar = paramPacker((HelixTraj::ParIndex)ipar);
    const ComPackBase<double>& packerr = errorPacker((HelixTraj::ParIndex)ipar);
    unsigned iperr = _params[ipar]>>packpar.bitRange();
    ComPackBaseBase::StatusCode pcode = packpar.unpack(_params[ipar],pvec[ipar]);
    ComPackBaseBase::StatusCode perrcode =  packerr.unpack(iperr,perr[ipar]);
    assert(pcode != ComPackBaseBase::TAG_BAD &&
	  perrcode != ComPackBaseBase::TAG_BAD);
  }
// convert back from dip to tandip
  pvec[HelixTraj::tanDipIndex] = tan(pvec[HelixTraj::tanDipIndex]);
// build back the covariance matrix; this also establishes the diagonal terms
  for(ipar=HelixTraj::d0Index;ipar<HelixTraj::NHLXPRM;ipar++){
    for(jpar=HelixTraj::d0Index;jpar<=ipar;jpar++)
      corr.fast(ipar+1,jpar+1) *= perr[ipar]*perr[jpar];
  }
// Finally, we're ready to construct the helix traj
  return new HelixTraj(pvec,corr,range[0],range[1]);
}

void
TrkHelixData_001::flightRange(double range[2]) const {
  double fltlen,fltrange;
  unsigned ifltlen = _fltlen;
  unsigned ifltrange = _fltlen>>(_packfltlen.bitRange());
  ComPackBaseBase::StatusCode fltcode = _packfltlen.unpack(ifltlen,fltlen);
  ComPackBaseBase::StatusCode rngcode = _packfltrange.unpack(ifltrange,fltrange);
  assert(fltcode != ComPackBaseBase::TAG_BAD &&
	 rngcode != ComPackBaseBase::TAG_BAD);
// convert the flight range information back to standard format
  range[0]=fltlen;
  range[1]=fltlen+fltrange;
}

void
TrkHelixData_001::correlation(HepSymMatrix& corr,unsigned* iused) const {
// dimension the array and initialize to unit
  corr = HepSymMatrix(5,1);
  iused[0] = iused[1] = iused[2] = 0;
// unpack the signs
  unsigned iword(0);
  assert(iused[iword]+_packsigns.bitRange() <= sizeof(d_ULong)*BITSPERBYTE);
  int signs(0);
  ComPackBaseBase::StatusCode scode = _packsigns.unpack(_corr[iword]>>iused[iword],signs);
  assert(scode != ComPackBase<int>::TAG_BAD);
  iused[iword] += _packsigns.bitRange();
// unpack the correlation magnitudes
  HelixTraj::ParIndex ipar,jpar;
  unsigned icor;
  for(icor=0;icor<10;icor++){
    corrTerm(icor,ipar,jpar);
    const ComPackBase<double>& packcorr = corrPacker(ipar,jpar);
// advance to the next word if there's not enough space left in this one
    if(iused[iword] + packcorr.bitRange()>sizeof(d_ULong)*BITSPERBYTE){
      iword++;
      assert(iword<3);
    }
    ComPackBaseBase::StatusCode ccode = packcorr.unpack(_corr[iword]>>iused[iword],corr.fast(ipar+1,jpar+1));
    assert(ccode != ComPackBaseBase::TAG_BAD);
// sign them
    if(signs & 1<<icor)corr.fast(ipar+1,jpar+1)*= -1.0;
    iused[iword] += packcorr.bitRange();
  }
}

int
TrkHelixData_001::charge(const BField* field ) const {
// get omega packer
  const ComPackBase<double>& omegapacker = paramPacker(HelixTraj::omegaIndex);
// unpack omega
  double omega;
  ComPackBaseBase::StatusCode pcode = omegapacker.unpack(_params[HelixTraj::omegaIndex],omega);
  if(pcode == ComPackBaseBase::TAG_BAD || field == 0 || omega == 0.0){
    ErrMsg(error) << "Can't determine charge" << endmsg;
    return 0;
  }
  double product = field->bFieldNominal() * omega;
  return product < 0.0 ? 1 : -1;
}

Hep3Vector
TrkHelixData_001::momentum(const BField* field,double fltlen) const {
  static HelixTraj traj(TrkParams(5));
  unpackHelix(traj);
  return TrkMomCalculator::vecMom(traj,*field,fltlen);
}

void
TrkHelixData_001::setData(d_ULong params[5],d_ULong corr[3],d_ULong& fltlen) {
  for(unsigned ipar=0;ipar<5;ipar++)
    params[ipar] = _params[ipar];
  for(unsigned icorr=0;icorr<3;icorr++)
    corr[icorr] = _corr[icorr];
  fltlen = _fltlen;
}

void
TrkHelixData_001::unpackHelix(HelixTraj& traj) const {
//  only need parameters, do don't unpack anything else
  static HepSymMatrix corr(5,1);
  static double range[2] = {0.0,1.0};
  HepVector pvec(5,0);
  for(unsigned ipar=HelixTraj::d0Index;ipar<HelixTraj::NHLXPRM;ipar++){
    const ComPackBase<double>& packpar = paramPacker((HelixTraj::ParIndex)ipar);
    ComPackBaseBase::StatusCode pcode = packpar.unpack(_params[ipar],pvec[ipar]);
    assert(pcode != ComPackBaseBase::TAG_BAD);
  }
// convert back from dip to tandip
  pvec[HelixTraj::tanDipIndex] = tan(pvec[HelixTraj::tanDipIndex]);
// create a helixtraj from this
  traj = HelixTraj(pvec,corr,range[0],range[1]);
}


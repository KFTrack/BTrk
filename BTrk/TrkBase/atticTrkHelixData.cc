//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkHelixData.cc,v 1.19 2004/07/05 05:43:27 hulsberg Exp $
//
//  Description:
//  Class TrkHelixData;
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
#include "TrkBase/TrkHelixData.hh"
#include "TrkBase/TrkFitSummary.hh"
#include "CommonUtils/ComPackFlatFloat.hh"
#include "CommonUtils/ComPackInt.hh"
#include "CommonUtils/ComPackSignedExpFloat.hh"
#include "CommonUtils/ComPackBase.hh"
#include "BaBar/Constants.hh"
#include "ErrLogger/ErrLog.hh"
#include "BField/BField.hh"
#include <assert.h>
#include <math.h>
#define BITSPERBYTE 8

// packing statics.  Note that these define the meaning of the data
// and should NEVER BE CHANGED.  Or, if they must be changed, this
// class will need to be subclassed and any persistent data based on it
// have its schema 'migrated'.

const ComPackSignedExpFloat TrkHelixData::_packd0(18,6,128.0);
const ComPackFlatFloat TrkHelixData::_packphi0(-Constants::pi,Constants::pi,18);
const ComPackSignedExpFloat TrkHelixData::_packomega(18,6,0.5);
const ComPackSignedExpFloat TrkHelixData::_packz0(18,6,256.0);
const ComPackFlatFloat TrkHelixData::_packlambda(-Constants::pi/2.0,Constants::pi/2.0,18);

const ComPackExpFloat TrkHelixData::_packd0err(14,3,2.0e-4,4.0);
const ComPackExpFloat TrkHelixData::_packphi0err(14,2,2.0e-5,Constants::pi/8.0);
const ComPackExpFloat TrkHelixData::_packomegaerr(14,6,2.0e-6,0.0625);
const ComPackExpFloat TrkHelixData::_packz0err(14,5,2.0e-4,16.0);
const ComPackExpFloat TrkHelixData::_packlambdaerr(14,2,2.0e-5,Constants::pi/8.0);

const ComPackExpFloat TrkHelixData::_packcorrfine(9,3,1.0,0.0);
const ComPackExpFloat TrkHelixData::_packcorrcoarse(8,3,0.0,1.0);
const ComPackFlatFloat TrkHelixData::_packcorrflat(0.0,1.0,8);
const ComPackInt TrkHelixData::_packpid(0,PdtPid::nPidType);
const ComPackInt TrkHelixData::_packsigns(0,(1<<10)-1);

const ComPackExpFloat TrkHelixData::_packchiprob(10,9,0.0,1.0);
const ComPackSignedExpFloat TrkHelixData::_packfltlen(11,3,256.0);
const ComPackFlatFloat TrkHelixData::_packfltrange(0.0,51.2,9);
const ComPackInt TrkHelixData::_packflag(0,3);

const double TrkHelixData::_rescale(0.001); // rescale off-diags
const double TrkHelixData::_mindet(1.0e-8);// minimum unpacked determinant value

const ComPackBase<double>&
TrkHelixData::errorPacker(HelixTraj::ParIndex index) {
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
TrkHelixData::paramPacker(HelixTraj::ParIndex index) {
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
TrkHelixData::corrTerm(unsigned iterm,
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
TrkHelixData::corrPacker(HelixTraj::ParIndex ipar,
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

TrkHelixData::TrkHelixData() : _fltlen(0) {
  for(unsigned ipar=0;ipar<HelixTraj::NHLXPRM;ipar++)
    _params[ipar] = 0;
  for(unsigned iword=0;iword<3;iword++)
    _corr[iword] = 0;
}

TrkHelixData::TrkHelixData(d_ULong params[5], d_ULong corr[3], d_ULong fltlen)
{
  for (unsigned i=0; i<5; i++) _params[i]= params[i];
  for (unsigned i=0; i<3; i++) _corr[i]= corr[i];
  _fltlen= fltlen;
}

TrkHelixData::TrkHelixData(const HelixTraj* traj,
			   double chiprob,PdtPid::PidType pid,
			   TrkEnums::PackFlag flag) {
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
// leave space for the pid and signs in word 0.
    iused[0] = _packpid.bitRange() + _packsigns.bitRange();
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
    _corr[0] |= ((packsigns)<<_packpid.bitRange());
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
// pack the pid 
  unsigned packpid(0);
  ComPackBaseBase::StatusCode pidcode =_packpid.pack(pid,packpid);
  assert(pidcode != ComPackBase<int>::TAG_BAD);
  _corr[0] |= packpid;
// pack the chisq prob and flightlength/range and flag
  double fltlen,fltrange;
  fltlen = traj->lowRange();
  fltrange = traj->hiRange()-fltlen;
  d_ULong chiprobpack,fltlenpack,fltrangepack,flagpack;
  ComPackBaseBase::StatusCode chicode = _packchiprob.pack(chiprob,chiprobpack);
  ComPackBaseBase::StatusCode fltcode = _packfltlen.pack(fltlen,fltlenpack);
  ComPackBaseBase::StatusCode rngcode = _packfltrange.pack(fltrange,fltrangepack);
  ComPackBaseBase::StatusCode flgcode = _packflag.pack(flag,flagpack);

  assert(chicode != ComPackBaseBase::TAG_BAD &&
	 fltcode != ComPackBaseBase::TAG_BAD &&
	 rngcode != ComPackBaseBase::TAG_BAD &&
	 flgcode != ComPackBaseBase::TAG_BAD);
  _fltlen = chiprobpack | (fltlenpack<<_packchiprob.bitRange()) |
    ( fltrangepack<<(_packchiprob.bitRange()+_packfltlen.bitRange())) |
    (flagpack<<(_packchiprob.bitRange()+_packfltlen.bitRange()+_packfltrange.bitRange()));
}

TrkHelixData::TrkHelixData(const TrkHelixData& other){
  *this = other;
}

TrkHelixData& 
TrkHelixData::operator =(const TrkHelixData& other) {
  if(this != &other){
    for (unsigned ipar=0;ipar<HelixTraj::NHLXPRM;ipar++)
      _params[ipar] = other._params[ipar];
    for (unsigned icor=0;icor<3;icor++)
      _corr[icor] = other._corr[icor];
    _fltlen = other._fltlen;
  }
  return *this;
}

TrkHelixData::~TrkHelixData(){}

HelixTraj*
TrkHelixData::helix() const {
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
		       
double
TrkHelixData::chisqProb() const{
  double chiprob;
  ComPackBaseBase::StatusCode chicode = _packchiprob.unpack(_fltlen,chiprob);
  assert(chicode != ComPackBaseBase::TAG_BAD);
  return chiprob;
}

void
TrkHelixData::flightRange(double range[2]) const {
  double fltlen,fltrange;
  unsigned ifltlen = _fltlen>>_packchiprob.bitRange();
  unsigned ifltrange = _fltlen>>(_packchiprob.bitRange() + _packfltlen.bitRange());
  ComPackBaseBase::StatusCode fltcode = _packfltlen.unpack(ifltlen,fltlen);
  ComPackBaseBase::StatusCode rngcode = _packfltrange.unpack(ifltrange,fltrange);
  assert(fltcode != ComPackBaseBase::TAG_BAD &&
	 rngcode != ComPackBaseBase::TAG_BAD);
// convert the flight range information back to standard format
  range[0]=fltlen;
  range[1]=fltlen+fltrange;
}

void
TrkHelixData::correlation(HepSymMatrix& corr,unsigned* iused) const {
// dimension the array and initialize to unit
  corr = HepSymMatrix(5,1);
  iused[0] = _packpid.bitRange(); // skip the pid work, which is stored first
  iused[1] = iused[2] = 0;
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

PdtPid::PidType
TrkHelixData::pidType() const {
  int ipid(0);
  ComPackBaseBase::StatusCode pidcode = _packpid.unpack(_corr[0],ipid);
  assert(pidcode != ComPackBaseBase::TAG_BAD);
  return (PdtPid::PidType)ipid;
}

TrkEnums::PackFlag
TrkHelixData::fitFlag() const {
  unsigned iflag = _fltlen>>(_packchiprob.bitRange() + _packfltlen.bitRange() +
			     _packfltrange.bitRange());
  int jflag;
  ComPackBaseBase::StatusCode fcode = _packflag.unpack(iflag,jflag);
  assert(fcode != ComPackBaseBase::TAG_BAD);
  return (TrkEnums::PackFlag)jflag;
}

TrkFitSummary
TrkHelixData::fitSummary() const {
  return TrkFitSummary(pidType(),
		       fitFlag(),
		       chisqProb(),
		       helix());
}

int
TrkHelixData::charge(const BField* field ) const {
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

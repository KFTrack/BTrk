//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: KalMiniTrk_001.cc,v 1.3 2003/10/02 23:23:02 brownd Exp $
//
// Description:
//      class KalMiniTrk_001.
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Copyright Information:
//	Copyright (C) 2000	Lawrence Berkeley Laboratory
//
// Author List:
//      Dave Brown 11/14/00
//------------------------------------------------------------------------

#include "BaBar/BaBar.hh"
#include "KalmanTrack/KalMiniTrk_001.hh"
#include "TrkBase/TrkRecoTrk.hh"
#include "ErrLogger/ErrLog.hh"
#include "PDT/PdtPid.hh"
#include "CommonUtils/ComPackFlatFloat.hh"

//#define TESTHYPO

const ComPackFlatFloat KalMiniTrk_001::_packt0(0.0,2.0e-6,16); // pack into 16 bits
const unsigned KalMiniTrk_001::_hypomsk(0x7); // mask for 1 hypo
const unsigned KalMiniTrk_001::_default(0x7); // default hypo flag
const unsigned KalMiniTrk_001::_failed(0x6); // failed fit flag
const unsigned KalMiniTrk_001::_hyposft(3); // shift (number of bits) for 1 hypo

unsigned
KalMiniTrk_001::hypoShift(unsigned ihypo) {
  return _hyposft*ihypo;
}

KalMiniTrk_001::KalMiniTrk_001() : _hypomap(0),_t0(0),_id(0)
{}

KalMiniTrk_001::KalMiniTrk_001(d_UShort hypomap,d_UShort t0,d_UShort id): 
  _hypomap(hypomap),_t0(t0),_id(id)
{}

KalMiniTrk_001&
KalMiniTrk_001::operator = (const KalMiniTrk_001& other ) {
  if(this != &other){
    _hypomap = other._hypomap;
    _t0 = other._t0;
    _id = other._id;
  }
  return *this;
}

KalMiniTrk_001::KalMiniTrk_001(const TrkRecoTrk* trk) : _hypomap(0),_t0(0),
						_id(trk->id())
{
// pack t0
  unsigned packt0;
  assert(_packt0.pack(trk->trackT0(),packt0) != ComPackBase<double>::TAG_BAD);
  _t0 = packt0;
// pack the hypo map.  Default is given a special value
  PdtPid::PidType defhypo = trk->defaultType(); 
// map the other hypos
  for(unsigned ihypo=0;ihypo<PdtPid::nPidType;ihypo++){
    PdtPid::PidType hypo = (PdtPid::PidType)ihypo;
    if(hypo == defhypo){
// the default fit MUST be valid
      assert(trk->status()->fitValid());
      _hypomap |= (_default<<hypoShift(hypo));
    } else {
      if(trk->status(hypo)->fitValid())
	_hypomap |= ((trk->whichFit(hypo)&_hypomsk)<<hypoShift(hypo));
      else
	_hypomap |= (_failed<<hypoShift(hypo));
    }
  }
}

KalMiniTrk_001::~KalMiniTrk_001()
{}

PdtPid::PidType
KalMiniTrk_001::defaultHypo() const {
  for (unsigned ihypo=0; ihypo<PdtPid::nPidType;ihypo++){
    if(((_hypomap>>hypoShift(ihypo))&_hypomsk) == _default)
      return (PdtPid::PidType)ihypo;
  }
  return PdtPid::null;
}

PdtPid::PidType
KalMiniTrk_001::fitHypo(PdtPid::PidType hypo) const {
  unsigned ihypo = (_hypomap>>hypoShift(hypo))&_hypomsk;
  if(ihypo == _default)
    return hypo;
  else if(ihypo == _failed)
    return PdtPid::null;
  else
    return (PdtPid::PidType)ihypo;
}

unsigned long
KalMiniTrk_001::trackId() const {
  return _id;
}

double
KalMiniTrk_001::usedT0() const {
  double t0val;
  assert(_packt0.unpack(_t0,t0val) != ComPackBase<double>::TAG_BAD);
  return t0val;
}

unsigned
KalMiniTrk_001::nFit() const {
  unsigned nfit(0);
  for (unsigned ihypo=0; ihypo<PdtPid::nPidType;ihypo++){
    PdtPid::PidType hypo = (PdtPid::PidType)ihypo;
    if(fitHypo(hypo) == hypo)nfit++;
  }
  return nfit;
}

unsigned
KalMiniTrk_001::nFailed() const {
  unsigned nfail(0);
  for (unsigned ihypo=0; ihypo<PdtPid::nPidType;ihypo++){
    PdtPid::PidType hypo = (PdtPid::PidType)ihypo;
    if(fitHypo(hypo) == PdtPid::null)nfail++;
  }
  return nfail;
}

unsigned
KalMiniTrk_001::nMapped() const {
  unsigned nmap(0);
  for (unsigned ihypo=0; ihypo<PdtPid::nPidType;ihypo++){
    PdtPid::PidType hypo = (PdtPid::PidType)ihypo;
    PdtPid::PidType fithypo = fitHypo(hypo);
    if(fithypo != PdtPid::null && fithypo != hypo)nmap++;
  }
  return nmap;
}

void
KalMiniTrk_001::fits(std::vector<PdtPid::PidType>& hypos) const {
  hypos.clear();
  hypos.reserve(PdtPid::nPidType);
  for (unsigned ihypo=0; ihypo<PdtPid::nPidType;ihypo++){
    PdtPid::PidType hypo = (PdtPid::PidType)ihypo;
    if(fitHypo(hypo) == hypo)
      hypos.push_back(hypo);
  }
}

int
KalMiniTrk_001::index(PdtPid::PidType myhypo) const {
  int retval(-1);
  unsigned icount(0);
  for (unsigned ihypo=0; ihypo<PdtPid::nPidType;ihypo++){
    PdtPid::PidType hypo = (PdtPid::PidType)ihypo;
    if(fitHypo(hypo) == hypo){
      if(hypo == myhypo){
	retval = icount;
	break;
      }
      icount++;
    }
  }
  return retval;
}


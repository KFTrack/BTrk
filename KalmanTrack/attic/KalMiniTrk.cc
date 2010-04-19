//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: KalMiniTrk.cc,v 1.21 2004/08/06 06:12:48 bartoldu Exp $
//
// Description:
//      class KalMiniTrk.
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
#include "KalmanTrack/KalMiniTrk.hh"
#include "KalmanTrack/KalMiniRep.hh"
#include "KalmanTrack/KalInterface.hh"
#include "KalmanTrack/KalMiniInterface.hh"
#include "KalmanTrack/KalRep.hh"
#include "TrkBase/TrkRecoTrk.hh"
#include "TrkBase/TrkFitSummary.hh"
#include "TrkBase/TrkErrCode.hh"
#include "TrkBase/TrkKalComposite.hh"
#include "TrkBase/TrkFitStatus.hh"
#include "ErrLogger/ErrLog.hh"
#include "PDT/PdtPid.hh"
#include "CommonUtils/ComPackFlatFloat.hh"
using std::cout;
using std::endl;

//#define TESTHYPO

const ComPackFlatFloat KalMiniTrk::_packt0(0.0,2.0e-6,16); // pack into 16 bits
const unsigned KalMiniTrk::_hypomsk(0x7); // mask for 1 hypo
const unsigned KalMiniTrk::_hyposft(3); // shift (number of bits) for 1 hypo
const unsigned KalMiniTrk::_validmsk(1); // mask for validity flag
const unsigned KalMiniTrk::_currentmsk(2); // mask for current flag
const unsigned KalMiniTrk::_fitstatmsk(4); // mask for status flag
const unsigned KalMiniTrk::_nsvtmsk(0x7f);
const unsigned KalMiniTrk::_nsvtsft(0);
const unsigned KalMiniTrk::_ndchmsk(0x3ff);
const unsigned KalMiniTrk::_ndchsft(_nsvtsft+7);
const unsigned KalMiniTrk::_nfitmsk(0x1f);
const unsigned KalMiniTrk::_nfitsft(_ndchsft+10);

unsigned
KalMiniTrk::hypoShift(unsigned ihypo) {
  return _hyposft*ihypo;
}

KalMiniTrk::KalMiniTrk() : _hypomap(0),_fitstat(0),_t0(0),
			   _id(0),_counters(0)
{}

KalMiniTrk&
KalMiniTrk::operator = (const KalMiniTrk& other ) {
  if(this != &other){
    _hypomap = other._hypomap;
    _fitstat = other._fitstat;
    _t0 = other._t0;
    _id = other._id;
    _counters = other._counters;
    _seed = other._seed;
  }
  return *this;
}

KalMiniTrk::KalMiniTrk(const TrkRecoTrk* trk) : _hypomap(0),_fitstat(0),_t0(0),
						_id(trk->id()),_counters(0)
{
// counters get packed later
// pack t0
  unsigned packt0;
  assert(_packt0.pack(trk->trackT0(),packt0) != ComPackBase<double>::TAG_BAD);
  _t0 = packt0;
// pack the hypo map.  Default is given a special value
  PdtPid::PidType defhypo = trk->defaultType(); 
// map the other hypos
  for(unsigned ihypo=0;ihypo<PdtPid::nPidType;ihypo++){
    PdtPid::PidType hypo = (PdtPid::PidType)ihypo;
// store the fit status for all reps
    const TrkFitStatus* fstat = trk->status(hypo);
    unsigned stat(0);
    if(fstat->fitCurrent()) stat |= _currentmsk;
    if(fstat->fitValid()) stat |= _validmsk;
    if(fstat->fitStatus().success()) stat |= _fitstatmsk;
    _fitstat |= ((stat&_hypomsk)<<hypoShift(hypo));
    if(ihypo == defhypo)
      _hypomap |= (_hypomsk<<hypoShift(hypo));
    else
      _hypomap |= ((trk->whichFit(hypo)&_hypomsk)<<hypoShift(hypo));
  }
// get the hotlist for this, and create the hot info from it
  KalInterface kface;
  KalMiniInterface kmface;
  if(trk->attach(kface,defhypo) && kface.kalmanRep() != 0){
// set the seed
    _seed = TrkHelixData(static_cast<const HelixTraj*>(kface.kalmanRep()->seed()),0.0,defhypo);
// try mini interface if that failed
  } else if(trk->attach(kmface,defhypo)&& kmface.kalmanMiniRep() != 0){
    _seed = TrkHelixData(static_cast<const HelixTraj*>(kmface.kalmanMiniRep()->seedTrajectory()),0.0,defhypo);
  } else
    ErrMsg(error) << " cannot persist non-Kalman fit" << endmsg;
// testing
#ifdef TESTHYPO
  for(unsigned ihypo=0;ihypo<PdtPid::nPidType;ihypo++){
    PdtPid::PidType hypo = (PdtPid::PidType)ihypo;
    if(trk->whichFit(hypo) != fitHypo(hypo))
      cout << "Hypos for " << hypo << " don't match, trk = "
	   << trk->whichFit(hypo) << " , KalMiniTrk = " << fitHypo(hypo) << endl;
  }
  if(trk->defaultType() != defaultHypo())
    cout << "Default hypos don't match, trk = " << trk->defaultType()
	 << " , KalMiniTrk = " << defaultHypo() << endl;
#endif
}

KalMiniTrk::~KalMiniTrk()
{}

TrkRecoTrk*
KalMiniTrk::createMiniTrack(const KalContext& kalcon,
                            const TrkSimpTraj& seed,
                            TrkHotList* hotlist,
                            std::vector<TrkFitSummary>& fits) const {
// create the TrkRecoTrk
  PdtPid::PidType defhypo = defaultHypo();
  assert(defhypo != PdtPid::null);
  if( !isCurrent(defhypo)){
    ErrMsg(error) << "Default hypo has invalid  fit" << endmsg;
    return 0;
  }
  if( !fitStatus(defhypo).success() ){
    ErrMsg(error) << "Default hypo has failed fit" << endmsg;
    //We got here twice during Kan conversion. Commenting out the following
    //line gets rid of the crash without any side effect. AMokhtarani 2/4/04 
    //   return 0;
  }
  unsigned long trkid = trackId();
  unsigned nfit = nFit();
  TrkRecoTrk* thetrk = createTrack(defhypo,trkid,usedT0());
  if(thetrk != 0){
// find the fits
    std::vector<TrkSimpTraj*> conlist; conlist.reserve(nfit);
    std::vector<TrkSimpTraj*> fitlist; fitlist.reserve(nfit);
    double fitprob(0.0);
    double dummy(0.0);
    fitResult(fits,defhypo,TrkEnums::KalConstraint,conlist,dummy);
    fitResult(fits,defhypo,TrkEnums::KalFit,fitlist,fitprob);
// this code is only invoked by the Objy mini, which has only the full hot list
// create the default miniRep, and add it to the track.
    KalMiniRep* defrep = new KalMiniRep(thetrk,defhypo,
                                        seed.clone(),kalcon,
					hotlist,0,
                                        conlist,
                                        fitlist,
                                        fitprob);
    if(defrep != 0){
      setRep(*thetrk,defrep);
// copy over initial status information
      setStatus(defrep);
// loop over the other reps
      for (unsigned ihypo=0; ihypo<PdtPid::nPidType;ihypo++){
	PdtPid::PidType newhypo = (PdtPid::PidType)ihypo;
// create reps for non-default hypos
	if(newhypo != defhypo){
// if this hypo has a unique fit, set it up
	  if(newhypo == fitHypo(newhypo)){
// get the fits for this hypo
	    fitResult(fits,newhypo,TrkEnums::KalConstraint,conlist,dummy);
	    fitResult(fits,newhypo,TrkEnums::KalFit,fitlist,fitprob);
// create the KalMiniRep
	    KalMiniRep* newrep = new KalMiniRep(defrep,newhypo,
						conlist,
						fitlist,
						fitprob);
	    if(newrep != 0){
	      addHypoTo(*thetrk,newrep,newhypo);
	      setStatus(newrep);
	    }
	  } else {
// no unique fit for this hypo; set the track to point to the correct real hypo
	    repointHypo(*thetrk,newhypo,fitHypo(newhypo));
	  }
	}
      }
    }
  }
#ifdef TESTHYPO
  for(unsigned ihypo=0;ihypo<PdtPid::nPidType;ihypo++){
    PdtPid::PidType hypo = (PdtPid::PidType)ihypo;
    if(thetrk->whichFit(hypo) != fitHypo(hypo))
      cout << "Hypos for " << hypo << " don't match, trk = "
           << thetrk->whichFit(hypo) << " , KalMiniTrk = " << fitHypo(hypo) << endl;
  }
  if(thetrk->defaultType() != defaultHypo())
    cout << "Default hypos don't match, trk = " << thetrk->defaultType()
         << " , KalMiniTrk = " << defaultHypo() << endl;
#endif

  return thetrk;
}

PdtPid::PidType
KalMiniTrk::defaultHypo() const {
  for (unsigned ihypo=0; ihypo<PdtPid::nPidType;ihypo++){
    if(((_hypomap>>hypoShift(ihypo))&_hypomsk) == _hypomsk)
      return (PdtPid::PidType)ihypo;
  }
  return PdtPid::null;
}

PdtPid::PidType
KalMiniTrk::fitHypo(PdtPid::PidType hypo) const {
  unsigned ihypo = (_hypomap>>hypoShift(hypo))&_hypomsk;
  if(ihypo != _hypomsk)
    return (PdtPid::PidType)ihypo;
  else
    return hypo;
}

unsigned long
KalMiniTrk::trackId() const {
  return _id;
}

double
KalMiniTrk::usedT0() const {
  double t0val;
  assert(_packt0.unpack(_t0,t0val) != ComPackBase<double>::TAG_BAD);
  return t0val;
}

bool
KalMiniTrk::isValid(PdtPid::PidType hypo) const {
  unsigned istat = (_fitstat>>hypoShift(hypo))&_hypomsk;
  return (istat & _validmsk) != 0;
}

bool
KalMiniTrk::isCurrent(PdtPid::PidType hypo) const {
  return ((_fitstat>>hypoShift(hypo)) & _currentmsk) != 0;
}


TrkErrCode
KalMiniTrk::fitStatus(PdtPid::PidType hypo) const {
  if(((_fitstat>>hypoShift(hypo)) & _fitstatmsk) != 0)
    return TrkErrCode(TrkErrCode::succeed);
  else
    return TrkErrCode(TrkErrCode::fail);
}

unsigned
KalMiniTrk::nSvt() const {
  return (_counters>>_nsvtsft)&_nsvtmsk;
}

unsigned
KalMiniTrk::nDch() const {
  return (_counters>>_ndchsft)&_ndchmsk;
}

unsigned
KalMiniTrk::nFit() const {
  return (_counters>>_nfitsft)&_nfitmsk;
}

void
KalMiniTrk::setCounters(unsigned nsvt,unsigned ndch,unsigned nfit) {
  _counters = 0;
  _counters |= ((nsvt&_nsvtmsk)<<_nsvtsft);
  _counters |= ((ndch&_ndchmsk)<<_ndchsft);
  _counters |= ((nfit&_nfitmsk)<<_nfitsft);
}

TrkSimpTraj*
KalMiniTrk::seedTrajectory() const {
  return _seed.helix();
}

void
KalMiniTrk::setStatus(KalMiniRep* kmrep) const {
  PdtPid::PidType hypo = kmrep->particleType();
  kmrep->setValid(isValid(hypo));
  kmrep->setCurrent(isCurrent(hypo));
  kmrep->addHistory(fitStatus(hypo),"KalMiniTrk");
}


void
KalMiniTrk::fitResult(std::vector<TrkFitSummary>& fits,
                      PdtPid::PidType hypo,TrkEnums::PackFlag flag,
                      std::vector<TrkSimpTraj*>& trajs,double& fitprob) const {
  if(!trajs.empty()) trajs.clear();
  fitprob = 0.0;
// loop over all the fit results, and unpack the matches into the traj vector
  typedef std::vector<TrkFitSummary>::iterator Sums_iter;
  Sums_iter send = fits.end();
  for(Sums_iter ifit=fits.begin();ifit!=send;++ifit){
    if(ifit->pidType() == hypo &&
       ifit->fitFlag() == flag){
// check for duplicates
      bool found(false);
      typedef std::vector<TrkSimpTraj*>::iterator Trajs_iter;
      Trajs_iter end = trajs.end();
      for(Trajs_iter itraj=trajs.begin();itraj!=end&&!found;++itraj)
        found = (**itraj == *(ifit->traj()));
      if(!found){
        trajs.push_back(ifit->traj());
        fitprob = ifit->chisqProb();
      }
    }
  }
}

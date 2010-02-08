//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: KalMiniListPromoter.cc,v 1.1 2005/04/04 19:54:54 brownd Exp $
//
// Description:
//	Module KalMiniListPromoter.
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//     David Brown  April 1, 2005
//
// Copyright Information:
//	Copyright (C) 2005		Lawrence Berkeley Laboratory
//
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"
#include "KalmanTrack/KalMiniListPromoter.hh"
#include "KalmanTrack/KalMiniPromoter.hh"
#include "TrkBase/TrkRecoTrk.hh"
#include "ErrLogger/ErrLog.hh"
#include <assert.h>

KalMiniListPromoter::KalMiniListPromoter( const char* const theName, const char* const theDescription ) :
  AppModule(theName,theDescription),
  _trklist("TrackList",this,"Default"),
  _hypos("HyposToPromote",this)
{
  commands( )->append( &_trklist );
  commands( )->append( &_hypos );
  _hypos.addItem("Default",defhypo);
  _hypos.addItem("AllExisting",allexisting);
}


KalMiniListPromoter::KalMiniListPromoter( const KalMiniListPromoter& other, const char* const theName ) :
  AppModule(theName,"KalMiniListPromoter clone"),
  _trklist(other._trklist,this),
  _hypos(other._hypos,this)
{
  commands( )->append( &_trklist );
  commands( )->append( &_hypos );
// copy basic state
  _verbose.set(other._verbose.value());
  _production.set(other._production.value());
  _enableFrames.set(other._enableFrames.value());
}

KalMiniListPromoter* 
KalMiniListPromoter::clone(const char* cloneName) {
  return new KalMiniListPromoter(*this,cloneName);
}

KalMiniListPromoter::~KalMiniListPromoter( )
{}

AppResult 
KalMiniListPromoter::beginJob( AbsEvent* anEvent ) {

  return AppResult::OK;
}

AppResult 
KalMiniListPromoter::event( AbsEvent* anEvent ) {
  HepAList<TrkRecoTrk>* trklist = Ifd<HepAList<TrkRecoTrk> >::get(anEvent,_trklist.value());
  if(trklist != 0){
    unsigned ntrk = trklist->length();
    for(unsigned itrk=0;itrk<ntrk;++itrk){
      TrkRecoTrk* trk = (*trklist)[itrk];
// always promote the default hypo
      TrkErrCode def =
	KalMiniPromoter::promoteToHots(trk,trk->defaultType());
// then, try the others if requested
      if(def.success() && _hypos.value()==allexisting){
	for(unsigned ihypo=0;ihypo<PdtPid::nPidType;++ihypo){
	  PdtPid::PidType hypo = (PdtPid::PidType)ihypo;
// only promote if this hypo isn't mapped to another
	  if(trk->whichFit(hypo)==hypo)
	    TrkErrCode hypoerr =
	      KalMiniPromoter::promoteToHots(trk,hypo);
	}
      }
    }
  } else
    ErrMsg(error) <<"No track list found to promote" << _trklist.value() << endmsg;
  return AppResult::OK;
}

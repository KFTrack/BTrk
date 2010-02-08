//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: KalMiniRX.cc,v 1.9 2007/11/16 23:54:03 brownd Exp $
//
// Description:
//	Class KalMiniRX
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      David Brown 11/7/01
//
// Copyright Information:
//	Copyright (C) 2001		Lawrence Berkeley Laboratory
//
//------------------------------------------------------------------------

#include "BaBar/BaBar.hh"
//-----------------------
// This Class's Header --
//-----------------------
#include "KalmanTrack/KalMiniRX.hh"
//-------------
// C Headers --
//-------------
#include <assert.h>
//---------------
// C++ Headers --
//---------------
#include <iostream>

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "CLHEP/Alist/AList.h"
#include "AbsEvent/AbsEvent.hh"
#include "TrkBase/TrkRecoTrk.hh"
#include "KalmanTrack/KalMiniRep.hh"
#include "KalmanTrack/KalRep.hh"
#include "KalmanTrack/KalCodes.hh"
#include "KalmanTrack/KalMaker.hh"
#include "ErrLogger/ErrLog.hh"

//
//----------------
// Constructors --
//----------------

KalMiniRX::KalMiniRX( const char* const theName, 
	      const char* const theDescription )
  : KalRX( theName, theDescription ),
    _downgrade("Downgrade",this,true)
{
  commands( )->append( &_downgrade );
}
//--------------
// Destructor --
//--------------

KalMiniRX::~KalMiniRX( )
{}
//--------------
// Operations --
//--------------
//
//

void
KalMiniRX::repairTrack(TrkRecoTrk* trk) const {
// invoke the base class version
  KalRX::repairTrack(trk);
// if requested, look for failed fits and downgrade them.
// loop over hypos and downgrade the ones which failed to fit
  for(int ihypo=0;ihypo<PdtPid::nPidType;ihypo++){
    PdtPid::PidType hypo = PdtPid::PidType(ihypo);
    if(trk->whichFit(hypo) == hypo &&
       trk->fitResult(hypo) == 0) {
// we _must_ have a valid default fit for all tracks in the mini
      if(_downgrade.value() || hypo == trk->defaultType()){
	if(trk->attach(_kmiface,hypo)){
	  KalMiniRep* krep = _kmiface.kalMiniRep();
	  if(krep != 0){
	    TrkErrCode repaired = downgrade(krep);
// append this to the history
	    krep->addHistory(repaired,name());
	  } else
	    ErrMsg(error) << "Can't find KalMiniRep " << endmsg;
	}
      }
    }
  }
  if(trk->fitResult(trk->defaultType()) == 0)
    ErrMsg(error) << "Track has no default rep" << endmsg;
}

TrkErrCode
KalMiniRX::downgrade(KalMiniRep* krep) const {
  TrkErrCode retval(TrkErrCode::fail,91,"Failed to repair KalMiniRep");
// get the current cache state
  KalMiniRep::miniState state = krep->currentState();
  KalMiniRep::miniState newstate = (KalMiniRep::miniState)(state-1);
  ErrMsg(routine) << "Found failed Fit for hypo" << krep->particleType()
		  << " in state "
		  << KalMiniRep::stateName(state)
		  << " , downgrading to state "
		  << KalMiniRep::stateName(newstate) << endmsg;
//  lower the state by 1
  if(state > KalMiniRep::seed){
    retval = krep->changeState(newstate);
// try refitting
    if(retval.success()){
      retval = krep->fit();
      if(!retval.success()) // try lowering the state again
	return downgrade(krep);
    } else {
      ErrMsg(error) << "Error setting KalMiniRep state to state" 
		    << KalMiniRep::stateName(newstate)
		    << " , setting state to seed" << endmsg;
// fall back all the way: this MUST work
      retval = krep->changeState(KalMiniRep::seed);
    }
  } else
    ErrMsg(error) << "Cannot downgrade from seed state" << endmsg;
  return retval;
}


void
KalMiniRX::changeToElectron(TrkRecoTrk* trk) const {
// get default KalMiniRep state
  static KalMiniInterface defkmiface;
  KalMiniRep::miniState defstate(KalMiniRep::seed);
  if(trk->attach(defkmiface,trk->defaultType()) &&
     defkmiface.kalmanMiniRep() != 0)
    defstate = defkmiface.kalmanMiniRep()->currentState();
  TrkErrCode repaired;
  KalMiniRep* elerep(0);
// find electron rep electron
  if(defstate != KalMiniRep::seed){
// see if the electron Kalman fit already exists, if not add one
    if(trk->attach(_kmiface,PdtPid::electron) &&
       _kmiface.kalMiniRep() != 0){
      elerep = _kmiface.kalMiniRep();
    } else {
// create a new KalMiniRep for electron hypo.  By dfn this has no
// fits, constraints, or chisq prob.
      elerep = new KalMiniRep(*defkmiface.kalmanMiniRep(),
			      PdtPid::electron);
      assert(elerep != 0);
// add it to the track
      _maker->addHypoTo(*trk,elerep,PdtPid::electron);
    }
// set the state of the electron fit to be the same as the default
    repaired = elerep->changeState(defstate);
    if(repaired.success())
      repaired = elerep->fit();
// set the history for this rep regardless
    elerep->addHistory(repaired,name());
    if(repaired.success())
// successful fit: change track default
      _maker->changeDefault(*trk,PdtPid::electron);
    else
// downgrade this fit
      downgrade(elerep);
  } else
    ErrMsg(error) << "No default KalMiniRep! cannot change to electron fit" << endmsg;
}

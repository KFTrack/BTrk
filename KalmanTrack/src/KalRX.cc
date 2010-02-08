//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: KalRX.cc,v 1.27 2007/11/16 23:54:03 brownd Exp $
//
// Description:
//	Class KalRX
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
// Revision History:
//	20030620  M. Kelsey -- Add TCL to remove or delete irreparable tracks
//------------------------------------------------------------------------

#include "BaBar/BaBar.hh"
//-----------------------
// This Class's Header --
//-----------------------
#include "KalmanTrack/KalRX.hh"
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
#include "AbsEvent/getTmpAList.hh"
#include "AbsEnv/AbsEnv.hh"
#include "TrkEnv/TrkEnv.hh"
#include "TrkBase/TrkRecoTrk.hh"
#include "TrkBase/TrkExchangePar.hh"
#include "KalmanTrack/KalCodes.hh"
#include "KalmanTrack/KalRep.hh"
#include "KalmanTrack/KalMaker.hh"
#include "ErrLogger/ErrLog.hh"
#include "PDT/PdtPid.hh"
using std::cout;
using std::endl;

// statics

TrkErrCode
KalRX::_repairstopping(TrkErrCode::succeed,50,"Repair stopping track");
TrkErrCode
KalRX::_refitting(TrkErrCode::succeed,51,"Refit not-current track");
TrkErrCode
KalRX::_electron(TrkErrCode::succeed,52,"Change default PID to electron");
TrkErrCode
KalRX::_grazer(TrkErrCode::succeed,53,"Remove grazing material");
TrkErrCode
KalRX::_deactivate(TrkErrCode::succeed,54,"Deactivate hots beyond stopping point");
TrkErrCode
KalRX::_unphysical(TrkErrCode::fail,55,"Fit invalidated because parameters unphysical");
TrkErrCode
KalRX::_diverged(TrkErrCode::fail,56,"Track repair diverging ");
TrkErrCode
KalRX::_unknown(TrkErrCode::fail,57,"Unknown failure, cannot repair ");

//
//----------------
// Constructors --
//----------------

KalRX::KalRX( const char* const theName, 
	      const char* const theDescription )
  : AppModule( theName, theDescription ),
    _inputkey("InputList",this,"Default"),
    _junkkey("JunkList",this,"RXJunk"),
    _makecurrent("MakeCurrent",this,true),
    _ele("MakeElectron",this,true),
    _inputaction("inputListAction",this),
    _repairhypos("HyposToRepair",this),
    _grazernormcut("GrazerCosIncident",this,0.1),
    _matlencut("GrazerMaterialLength",this,0.3),
    _omegacut("UnphysicalOmegaCut",this,100.0),
    _npastcut("NHitsPastStop",this,2),
    _maker(0)
{
  commands( )->append( &_inputkey );
  commands( )->append( &_junkkey );
  commands( )->append( &_makecurrent );
  commands( )->append( &_ele );
  commands( )->append( &_inputaction );
  commands( )->append( &_repairhypos );
  commands( )->append( &_grazernormcut );
  commands( )->append( &_matlencut );
  commands( )->append( &_omegacut ) ;
  commands( )->append( &_npastcut );

  // Define name-value pairs, (first is default) for failed-track bookeeping
  _inputaction.addItem("junktrk",KalRX::junktrk);
  _inputaction.addItem("noaction",KalRX::noaction);
  _inputaction.addItem("deletetrk",KalRX::deletetrk);
  _repairhypos.addItem("default",KalRX::defhypo);
  _repairhypos.addItem("all",KalRX::all);

}
//--------------
// Destructor --
//--------------

KalRX::~KalRX( )
{
  delete _maker;
}
//--------------
// Operations --
//--------------
AppResult
KalRX::beginJob( AbsEvent* anEvent )
{
  ErrMsg(trace)  << " begin Job" << endmsg;

  const KalContext* kalcon = gblEnv->getTrk()->kalContext();
//  Construct a maker object around this context
  if(kalcon != 0)_maker = new KalMaker(*kalcon);
  assert(_maker != 0);
  return AppResult::OK;
}
//
//
AppResult
KalRX::event( AbsEvent* anEvent )
{
//  Get the correct track list from the event
  HepAList<TrkRecoTrk>* trklist = Ifd< HepAList<TrkRecoTrk> >::get(anEvent,_inputkey.value());
  if(trklist == 0){
    ErrMsg(error) << ": no list found in event, aborting." << endmsg;
    return AppResult::OK;
  }
  HepAList<TrkRecoTrk> junk; // list of junk tracks
//  Loop over the tracks
  unsigned ntrks = trklist->length();
  for(unsigned itrk=0;itrk<ntrks;itrk++){
    TrkRecoTrk* trk = (*trklist)[itrk];
// repair this track if necessary
    repairTrack(trk);
// if the default fit is still invalid, junk the track
    if (trk->fitResult() == 0)
      junk.append(trk);
  }
// if the junk list has entries, deal with them
  if(junk.length() != 0 && _inputaction.value() != noaction){
    unsigned njunk = junk.length();
    for(unsigned ijunk=0;ijunk<njunk;++ijunk){
// remove the track from the input
      trklist->remove(junk[ijunk]);
      if(_inputaction.value()==deletetrk)
	delete junk[ijunk];
    }
// put the junk list into the event if requested
    if(_inputaction.value() == junktrk){
      HepAList<TrkRecoTrk>* junklist(0);
      getTmpAList(anEvent,junklist,_junkkey.value());
      *junklist = junk;
    }	  
  }
  return AppResult::OK;
}

void
KalRX::repairTrack(TrkRecoTrk* trk) const {
// loop over the reps
  for(int ihypo=0;ihypo<PdtPid::nPidType;ihypo++){
    PdtPid::PidType hypo = PdtPid::PidType(ihypo);
// only look at requested fits
    if(trk->whichFit(hypo) == hypo &&
       (_repairhypos.value() == all ||
	(_repairhypos.value() == defhypo && trk->defaultType() == hypo))) {
// some fixes are applied even to valid fits
      if(trk->fitResult(hypo)){
	if(_makecurrent.value() && (!trk->status(hypo)->fitCurrent())){
// look for tracks which aren't current
	  if(verbose()){
	    cout << "Found not-current track, history:" << endl;
	    trk->status()->printHistory(cout);
	    cout << endl;
	  }
// try to refit the track
	  trk->status()->addHistory(_refitting,name());
	  TrkErrCode refit = trk->hits(hypo)->fit();
	  trk->status()->addHistory(refit,name());
// if the fit failed, try to repair it
	  if(refit.failure())
	    repairFit(trk,hypo);
	}
      } else
	repairFit(trk,hypo);
// failed fit: try to repair it
// make sure the fit is 'physical', and if not invalidate it
      invalidateUnphysical(trk,hypo);
    }
  }
// if the default fit is invalid, try to change the default to the electron fit
  if(trk->fitResult() == 0  &&  trk->defaultType() != PdtPid::electron && _ele.value()) {
    changeToElectron(trk);
// make sure the electron hypo is physical
    invalidateUnphysical(trk,PdtPid::electron);
  }
}

void
KalRX::changeToElectron(TrkRecoTrk* trk) const {
// try refitting it as an electron
  if(trk->defaultType() != PdtPid::electron ) {
    TrkErrCode repaired;
// see if the electron Kalman fit already exists, if not add one
    if(trk->attach(_kiface,PdtPid::electron) &&
       _kiface.kalRep() != 0 && 
       _kiface.kalRep()->particleType()==PdtPid::electron){
// change track default
      _maker->changeDefault(*trk,PdtPid::electron);
// make sure it's current
      repaired = _kiface.kalRep()->fit();
    } else {
      repaired = _maker->addHypo(*trk,PdtPid::electron,true);
      if(repaired.success()) {
// try fitting the new hypo
	repaired = _maker->fitHypo(*trk,PdtPid::electron);
      }
    }
// if this fails, try to repair it
    if(repaired.failure())
      repairFit(trk,PdtPid::electron);
  }
}

void
KalRX::repairFit(TrkRecoTrk* trk,PdtPid::PidType hypo) const {
  TrkErrCode repaired = _unknown;
  KalRep* krep(0);
  if(trk->attach(_kiface,hypo))
    krep = _kiface.kalRep();
  if(krep != 0){
// see if the failure code of the last fit is one we recognize
    switch(krep->fitStatus().failure()){
    case KalCodes::stops:
      repaired = repairStopping(krep);
      break;
    case KalCodes::diverge:
      break;
    case KalCodes::dof:
      break;
    case KalCodes::matrix:
      break;
    case KalCodes::processing:
      break;
    default:
      break;
    }
// add the status to the history
    if(repaired.failure() != _unknown.failure())
      krep->addHistory(repaired,name());
    if(verbose()){
      cout << "Found failed track, history:" << endl;
      krep->printHistory(cout);
      cout << endl;
    }
  } else
    ErrMsg(error) << "Can't repair non-Kalman fit" << endmsg;
}

TrkErrCode
KalRX::repairStopping(KalRep* krep) const {
  TrkErrCode repaired( TrkErrCode::fail,KalCodes::stops,
		       "Track stops due to energy loss before last hit");
// first, decide if this track is a 'grazer'.  If so , there are 2 possible fixes:
// If the grazing intersection is spurious, try removing it.  Otherwise if we have just
// picked up a few hits after a real intersection, we can recover the fit by dropping those
// hits.  If the track isn't a grazer, there's nothing that can be done.
//
  return repaired;
}

void
KalRX::invalidateUnphysical(TrkRecoTrk* trk, PdtPid::PidType hypo) const {
  TrkFitStatus* status = trk->status(hypo);
  const TrkFit* fit = trk->fitResult(hypo);
// make a cut on omega to remove very 'unphysical' track fits
  if( fit != 0 && !(fit->helix(0).omega()<_omegacut.value()) ) {
    status->setValid(false) ;
    status->addHistory(_unphysical,name()) ;
    if(verbose()){
      cout << "Removing unphysical fit. Track history: " << endl ;
      trk->status()->printHistory(cout);
      cout << endl;
    }
  }
}

//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: KalFit.cc,v 1.95 2006/11/20 18:31:34 desilva Exp $
//
// Description:
//	Class KalFit. This module tests the kalman track fit
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      David Brown 4/17/97
//
// Copyright Information:
//	Copyright (C) 1997		Lawrence Berkeley Laboratory
//
//------------------------------------------------------------------------

#include "BaBar/BaBar.hh"
//-----------------------
// This Class's Header --
//-----------------------
#include "KalmanTrack/KalFit.hh"
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
#include "CLHEP/Alist/AIterator.h"
#include "AbsEvent/AbsEvent.hh"
#include "AbsEvent/getTmpAList.hh"
#include "AbsEnv/AbsEnv.hh"
#include "TrkEnv/TrkEnv.hh"
#include "TrkBase/TrkRecoTrk.hh"
#include "TrkBase/TrkExchangePar.hh"
#include "TrkBase/TrkErrCode.hh"
#include "TrkBase/TrkVolume.hh"
#include "TrkBase/TrkDirection.hh"
#include "RecoUtils/RecAListCopyParms.hh"
#include "BField/BFieldIntegrator.hh"
#include "KalmanTrack/KalMaker.hh"
#include "ErrLogger/ErrLog.hh"
#include "PDT/PdtPid.hh"
#include "CLHEP/String/Strings.h"
#include "TrkEnv/KalContext.hh"
#include "TrkFitter/TrkMergeMap.hh"
#include "AbsParm/AbsParmIfdStrKey.hh"
class TrkRecoTrk;
class KalMaker;
#define NCONS 10
//
//----------------
// Constructors --
//----------------

KalFit::KalFit( const char* const theName, 
		const char* const theDescription )
  : AppModule( theName, theDescription ),
    _copyParms(0),
    _fitdir("FitDirection",this),
    _inputtype("InputType",this),
    _ministate("MiniState", this),
    _fithypos("HyposToFit",this),
    _constraints("Constraints",this),
    _maptype("MapTypes",this,true),
    _maker(0), _cons(new bool[NCONS])
{
  _copyParms = new RecAListCopyParms<TrkRecoTrk>(this, "Default", "Default");

  commands( )->append( &_inputtype );
  commands( )->append( &_ministate );
  commands( )->append( &_fithypos );
  commands( )->append( &_constraints );
  commands( )->append( &_fitdir );
  commands( )->append( &_maptype );

// first item added becomes default.
  _ministate.addItem("Hots",KalMiniRep::hots);
  _ministate.addItem("Cache",KalMiniRep::cache);
  _ministate.addItem("ExtendedCache",KalMiniRep::extendedcache);
  _inputtype.addItem("Helices",KalFit::helix);
  _inputtype.addItem("Kalman",KalFit::kalman);
  _inputtype.addItem("KalMini",KalFit::kalmini);
  _fithypos.addItem("Default",KalFit::defhypo);
  _fithypos.addItem("None",KalFit::none);
  _fithypos.addItem("All",KalFit::all);
  _fithypos.addItem("AllExisting",KalFit::allexisting);
  _fithypos.addItem("electron",KalFit::electron);
  _fithypos.addItem("muon",KalFit::muon);
  _fithypos.addItem("pion",KalFit::pion);
  _fithypos.addItem("kaon",KalFit::kaon);
  _fithypos.addItem("proton",KalFit::proton);
  _fithypos.addItem("Electron",KalFit::electron);
  _fithypos.addItem("Muon",KalFit::muon);
  _fithypos.addItem("Pion",KalFit::pion);
  _fithypos.addItem("Kaon",KalFit::kaon);
  _fithypos.addItem("Proton",KalFit::proton);
  _constraints.addItem("D0",TrkExchangePar::ex_d0);
  _constraints.addItem("Phi0",TrkExchangePar::ex_phi0);
  _constraints.addItem("Omega",TrkExchangePar::ex_omega);
  _constraints.addItem("Z0",TrkExchangePar::ex_z0);
  _constraints.addItem("TanDip",TrkExchangePar::ex_tanDip);
  _fitdir.addItem("Both",KalFit::fitBoth);
  _fitdir.addItem("In",KalFit::fitIn);
  _fitdir.addItem("Out",KalFit::fitOut);
  for(unsigned icons=0;icons<NCONS;icons++)
    _cons[icons]=false;
}

// copy constructor, used in cloning
KalFit::KalFit(const KalFit& other,const char* name) :
  AppModule( name, other.description() ),
  _copyParms(0),
  _fitdir(other._fitdir,this),
  _inputtype(other._inputtype,this),
  _ministate(other._ministate, this),
  _fithypos(other._fithypos,this),
  _constraints(other._constraints,this),
  _maptype(other._maptype,this),
  _maker(0), _cons(new bool[NCONS])
{
// copy the copyParms object.  This is convoluted
  _copyParms = new RecAListCopyParms<TrkRecoTrk>(this,
						 other._copyParms->inputList()->valueString().c_str(),
						 other._copyParms->outputList()->valueString().c_str(),
						 other._copyParms->action());
  commands( )->append( &_inputtype );
  commands( )->append( &_ministate );
  commands( )->append( &_fithypos );
  commands( )->append( &_constraints );
  commands( )->append( &_fitdir );
  commands( )->append( &_maptype );
// set the state of base-class parameters
  _verbose.set(other._verbose.value());
  _production.set(other._production.value());
  _enableFrames.set(other._enableFrames.value());
  for(unsigned icons=0;icons<NCONS;icons++)
    _cons[icons]=other._cons[icons];
}

KalFit*
KalFit::clone(const char* clonename) {
  return new KalFit(*this,clonename);
}

//--------------
// Destructor --
//--------------

KalFit::~KalFit( )
{
  delete _maker;
  delete[] _cons;
  delete _copyParms;
}
//--------------
// Operations --
//--------------
AppResult
KalFit::beginJob( AbsEvent* anEvent )
{
  ErrMsg(trace)  << " begin Job" << endmsg;
// get KalContext from TrkEnv
  const KalContext* kalcon = gblEnv->getTrk()->kalContext();
//  Construct a maker object around this context
  if(kalcon != 0)_maker = new KalMaker(*kalcon);
  assert(_maker != 0);
// setup the constraints as desired
  unsigned ncons = _constraints.value().size();
  for(unsigned icons=0;icons<ncons;icons++){
    if(_constraints.value()[icons]>=0 &&
       _constraints.value()[icons]<NCONS)
      _cons[_constraints.value()[icons]] = true;
  }
  return AppResult::OK;
}
//
//
AppResult
KalFit::event( AbsEvent* anEvent )
{
// setup the lists
  _copyParms->prepare(anEvent);
//  Get the correct track list from the event
  HepAList<TrkRecoTrk>* trklist =  _copyParms->findInputList();
  if(trklist == 0){
    ErrMsg(routine) << ": no list found in event, aborting." << endmsg;
    return AppResult::OK;
  }
  TrkMergeMap* mMap(0);
// if we're copying, setup the merge map
  if(_copyParms->copyFromInput()){
    static const IfdStrKey defaultKey("Default");
    mMap = Ifd< TrkMergeMap >::get( anEvent, defaultKey );
    if (mMap == 0) { 
      mMap = new TrkMergeMap;
      IfdDataProxy<TrkMergeMap> *theProxy = new IfdDataProxy<TrkMergeMap>(mMap);
      
      if(!Ifd<TrkMergeMap>::put(anEvent, theProxy, defaultKey))
	ErrMsg(error) << ": Couldn't add merge map to event." << endmsg;
    }
  }
//  Loop over the tracks, and fit them using the Kalman fit
  TrkRecoTrk *thetrk(0);
  HepAListIterator<TrkRecoTrk> iterTrk(*trklist);
  while ( (thetrk = iterTrk())) {
// get the output
    TrkRecoTrk* fittrk = _copyParms->addOutput(*thetrk);
    if(fittrk != thetrk){
// add an entry in the merge map (lie about the 2nd parent)
// follow the hidden convention track1=svtonly,track2=dchonly
      if(thetrk->fitResult()->nDch()==0)
	mMap->insert(fittrk,thetrk,0);
      else
	mMap->insert(fittrk,0,thetrk);
    }
    unsigned ihypo;
    PdtPid::PidType hypo;
    TrkErrCode fitresult;
    KalMiniRep::miniState ministate = (KalMiniRep::miniState)_ministate.value();
// Setup the track for fitting
    switch ( _inputtype.value()) {
    case KalFit::helix: default:
// make sure the default hypo isn't already a Kalman track
      if(! fittrk->attach(_kiface,fittrk->defaultType())){
// setup the Kalman rep, for the specified particle type.
	fitresult = _maker->changeFit(*fittrk);
// add fits as required
	for (ihypo = 0; ihypo < PdtPid::nPidType; ihypo++) {
	  fitresult = TrkErrCode();
	  hypo = (PdtPid::PidType) ihypo;
	  if(fitThisHypo(fittrk,hypo)){
	    if(hypo != fittrk->defaultType())
	      fitresult = fittrk->addFit(hypo,false);
	    if(fitresult.success() && fittrk->attach(_kiface,hypo)){
	      fitresult = fitKalRep(_kiface.kalRep());
	      if(verbose() && fitresult.failure())
		ErrMsg(warning) << ": error fitting hypo " 
				<< hypo << " " << fitresult << endmsg;
	    } else
	      ErrMsg(error) << "Cannot attach kalman interface to hypo" << hypo << endmsg;
	  }
	}
      } else
	ErrMsg(error) << "Input track is already a kalman track" << endmsg;
      break;
    case KalFit::kalman:
// make sure the default hypo is really a Kalman track
      if(fittrk->attach(_kiface,fittrk->defaultType())){
// Fit all the mass hypos which exist explicitly
	for(ihypo=0;ihypo<PdtPid::nPidType;ihypo++){
	  fitresult = TrkErrCode();
	  hypo = PdtPid::PidType(ihypo);
	  if(fitThisHypo(fittrk,hypo)){
// add the fit if required
	    if(fittrk->whichFit(hypo) != hypo)
	      fitresult = fittrk->addFit(hypo,false);
	    if(fitresult.success() && fittrk->attach(_kiface,hypo)){
// set constraints
	      _kiface.kalRep()->setConstraints(_cons);
// fit it
	      fitresult = fitKalRep(_kiface.kalRep());
	      if(verbose() && fitresult.failure())
		ErrMsg(warning) << ": error fitting hypo " 
				<< hypo << " " << fitresult << endmsg;
	    } else
 	      ErrMsg(error) << name() << "Cannot attach kalman interface to hypo" << hypo << endmsg;
	  }
	}
      } else
	ErrMsg(error) << name() << "Input track is not a kalman track" << endmsg;
      break;
    case KalFit::kalmini:
// make sure the default hypo is a KalMiniRep
      if(fittrk->attach(_kminiiface,fittrk->defaultType())){
// fit the default hypo first.  This can improve the efficiency in refit mode
	std::vector<PdtPid::PidType> types;
	types.push_back(fittrk->defaultType());
	for(ihypo=0;ihypo<PdtPid::nPidType;ihypo++){
	  hypo = PdtPid::PidType(ihypo);
	  if(hypo != fittrk->defaultType())
	    types.push_back(hypo);
	}
	std::vector<PdtPid::PidType>::iterator ihypo = types.begin();
	while(ihypo != types.end()){
	  hypo = *ihypo;
	  ihypo++;
	  fitresult = TrkErrCode();
// change the mini-state regardless of the requested fit hypos
	  if(fitThisHypo(fittrk,hypo)){
// if we're in hots mode, see if we need to add hypos
	    if(ministate == KalMiniRep::hots){
	      if(fittrk->whichFit(hypo) != hypo)
		fitresult = fittrk->addFit(hypo,false);
	      if(! fittrk->attach(_kminiiface,hypo)){
// this is a non-Kalminrep rep (probably a dead rep).  replace it
		if(fittrk->attach(_kminiiface,fittrk->defaultType())){
		  static std::vector<TrkSimpTraj*> cache;
		  KalMiniRep* newrep = new KalMiniRep(_kminiiface.kalMiniRep(),
						      hypo,cache,cache,-1.0);
		  assert(newrep != 0);
		  _maker->addHypoTo(*fittrk,newrep,hypo);
		} else
		  ErrMsg(error) << "Default rep is not KalMiniRep" << endmsg;
	      }
	    }
	    if(fitresult.success() && fittrk->attach(_kminiiface,hypo)){
// set it to the desired state
	      KalMiniRep* minirep = _kminiiface.kalMiniRep();
	      fitresult = minirep->changeState(ministate);
	      if(fitresult.success()){
// If the state allows setting constraints, do so
		if(ministate == KalMiniRep::hots)
		  minirep->kalRep()->setConstraints(_cons);
// fit it
		fitresult = minirep->fit();;
		minirep->addHistory(fitresult,name());
		if(verbose() && fitresult.failure())
		  ErrMsg(warning) << ": error fitting hypo " 
				  << hypo << " " << fitresult << endmsg;
	      } else {
		ErrMsg(error) << "Error changing KalMini state" << endmsg;
		minirep->addHistory(fitresult,name());
	      }
	    } else
 	      ErrMsg(error) << "Cannot attach kalman interface to hypo" << hypo << endmsg;
	  }
	}
      } else
	ErrMsg(error) << "Input track is not a kalman mini track" << endmsg;
    }
  }
// cleanup
  _copyParms->finalize();
  return AppResult::OK;
}

AppResult
KalFit::endJob( AbsEvent* anEvent )
{
  ErrMsg(trace) << " end Job" << endmsg;
  delete _maker;
  _maker = 0;
  return AppResult::OK;
}

bool
KalFit::fitThisHypo(TrkRecoTrk* fittrk,PdtPid::PidType hypo) const {
  const std::vector<int>& hypos = _fithypos.value();
  bool retval = false;
  std::vector<int>::const_iterator ihypo = hypos.begin();
  while(!retval && ihypo != hypos.end()){
    fitHypos fhypo = (fitHypos)*ihypo;
    if((int)fhypo <= (int)PdtPid::proton) {
// explicit hypo requested. check if this hypo matches
      if(_maptype.value())
	retval |= fittrk->whichFit((PdtPid::PidType)fhypo) == hypo;
      else
	retval |= (PdtPid::PidType)fhypo == hypo;
    } else
// test logical values.  Must exclude dead reps from matching 'allexisting'
      retval |= fhypo == KalFit::all || 
	(fhypo == KalFit::allexisting && fittrk->whichFit(hypo) == hypo && 
	 (!fittrk->attach(_deadiface,hypo))) ||
	(fhypo == KalFit::defhypo && fittrk->defaultType() == hypo);
    ihypo++;
  }
  return retval;
}

TrkErrCode
KalFit::fitKalRep(KalRep* rep) const {
  TrkErrCode result;
  switch(_fitdir.value()){
  case fitBoth: default:
    result = rep->fit();
    break;
  case fitIn:
    result =  rep->fit(trkIn);
    break;
  case fitOut:
    result = rep->fit(trkOut);
    break;
  }
  rep->addHistory(result,name());
  return result;
}

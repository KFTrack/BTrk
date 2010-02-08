//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: KalPairFit.cc,v 1.23 2005/04/29 18:04:05 desilva Exp $
//
// Description:
//	Class KalPairFit
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      David Brown 4/1/98
//
// Copyright Information:
//	Copyright (C) 1998		Lawrence Berkeley Laboratory
//
//------------------------------------------------------------------------

#include "BaBar/BaBar.hh"
#include "KalmanTrack/KalPairFit.hh"
#include <assert.h>
#include <iostream>

#include "CLHEP/Alist/AList.h"
#include "CLHEP/Alist/AIterator.h"
#include "AbsEvent/AbsEvent.hh"
#include "AbsEvent/getTmpAList.hh"
#include "AbsEnv/AbsEnv.hh"
#include "AbsEvent/getTmpAList.hh"
#include "ErrLogger/ErrLog.hh"
#include "KalmanTrack/KalInterface.hh"
#include "KalmanTrack/KalPairConduit.hh"
#include "KalmanTrack/KalPairInterface.hh"
#include "KalmanTrack/KalPairMaker.hh"
#include "KalmanTrack/KalRep.hh"
#include "PepEnv/PepEnv.hh"
#include "PepEnvData/PepBeams.hh"
#include "TrkBase/TrkErrCode.hh"
#include "TrkBase/TrkRecoTrk.hh"
#include "TrkFitter/TrkMergeMap.hh"

//
//----------------
// Constructors --
//----------------

KalPairFit::KalPairFit( const char* const theName, 
		  const char* const theDescription )
  : AppModule( theName, theDescription ),
    _trklistkey("FitList",this,"Default"),
    _pairlistkey("PairList",this,"KalPair"),
    _mergeMapKey("MergeMapKey", this, "KalPair"),
    _makeMergeMap("MakeMergeMap", this, true),
    _deltaChi2Cut("DeltaChi2Cut", this, 100.),
    _pxbias("PxBias", this, 0.),
    _pybias("PyBias", this, 0.),
    _pzbias("PzBias", this, 0.),
    _pxerrscale("PxErrScale", this, 1.), 
    _pyerrscale("PyErrScale", this, 1.), 
    _pzerrscale("PzErrScale", this, 1.), 
    _maker(0)
{
  commands( )->append( &_trklistkey );
  commands( )->append( &_pairlistkey );
  commands( )->append( &_mergeMapKey );
  commands( )->append( &_makeMergeMap );
  commands( )->append( &_deltaChi2Cut );
  commands( )->append( &_pxbias );
  commands( )->append( &_pybias );
  commands( )->append( &_pzbias );
  commands( )->append( &_pxerrscale );
  commands( )->append( &_pyerrscale );
  commands( )->append( &_pzerrscale );
}

KalPairFit::~KalPairFit(){
}

AppModule*
KalPairFit::clone(const char* cloneName)
{
  return new KalPairFit(cloneName, this->description());
}

AppResult
KalPairFit::beginJob( AbsEvent* anEvent )
{
  ErrMsg(routine) << " begin Job" << endmsg;
  //
  //  Construct a maker object 
  //
  _maker = new KalPairMaker();
  
  return AppResult::OK;
}

//
AppResult
KalPairFit::event( AbsEvent* anEvent )
{
  // Get the list where the KalPair tracks will be stored.  This list should
  // already be in the event, but with no entries.
  HepAList<TrkRecoTrk> *kalpairlist;
  getTmpAList( anEvent, kalpairlist, _pairlistkey.value());
  if( kalpairlist == 0) {
    ErrMsg(error) << " Failed to find or put HepAList " << _pairlistkey.value()
		  << " into the event" << endmsg;
	return AppResult::OK;
  }

  // Get the input track list from the event
  HepAList<TrkRecoTrk>* trklist =
    Ifd<HepAList<TrkRecoTrk> >::get(anEvent, _trklistkey.value());

//   HepAList<TrkRecoTrk>* trklist =
//     Ifd<HepAList<TrkRecoTrk> >::
//     get(anEvent, _trklistkey.value());

  if( trklist == 0 ){
    ErrMsg(error) << " List " << _trklistkey.value() 
 		  << " not found in event, aborting." << endmsg;
    return AppResult::OK;
  }

  // Require there be exactly 2 tracks of opposite sign
  if ( trklist->length() != 2) {
    if( _verbose.value() )
      ErrMsg(routine) << " Tracklist " << _trklistkey.value() 
		      << " doesn't have 2 tracks, do nothing" << endmsg;
    return AppResult::OK;
  }

  TrkRecoTrk* trk0 = (*trklist)[0];
  TrkRecoTrk* trk1 = (*trklist)[1];
  TrkRecoTrk* pairtrack0(0);
  TrkRecoTrk* pairtrack1(0);
  if( 0 == trk0->fitResult() || 0 == trk1->fitResult() ) {
    ErrMsg(error) << "KalPairFit: track with no fitResult(), aborting."
		  << endmsg;
    return AppResult::OK;
  }

  // Make sure the tracks don't have the same charge
  if (trk0->fitResult()->charge() == trk1->fitResult()->charge()) {
    ErrMsg(error) << "KalPairFit: tracks with same charge, aborting."
		  << endmsg;
    return AppResult::OK;
  }

  // Make trk0 the positive track and trk1 the negative track
  if( trk0->fitResult()->charge() < trk1->fitResult()->charge() ) {
    TrkRecoTrk* temptrk = trk0;
    trk0 = trk1;
    trk1 = temptrk;
  }

  // Get the helix parameter at the origin of the two input tracks
  // to seed the constraints in the conduit
  KalInterface kalInter0;
  KalInterface kalInter1;
  if (!trk0->attach(kalInter0, trk0->defaultType()) ||
      !trk1->attach(kalInter1, trk1->defaultType())) {
    ErrMsg(error) << "Couldn't attach KalRep" << endmsg;
    return AppResult::OK;
  }

  // Get momentum constraint from PepEnv
  BbrVectorErr beamMomentum = 
    gblEnv->getPep()->pepBeams()->totalMomentumWErr();

  // Add option to tweak the beam momentum.  Mostly for doing studies.
  // First, add a possible bias.  Defaults are all zero.
  Hep3Vector bias(_pxbias.value(), _pybias.value(), _pzbias.value());
  BbrVectorErr biasBbr(bias);
  beamMomentum += biasBbr;

  // Next, scale covariance matrix. Default scales are 1.
  BbrError scaledError(beamMomentum.covMatrix());
  scaledError.fast(1,1) *= _pxerrscale.value()*_pxerrscale.value();
  scaledError.fast(2,1) *= _pyerrscale.value()*_pxerrscale.value();
  scaledError.fast(3,1) *= _pzerrscale.value()*_pxerrscale.value();
  scaledError.fast(2,2) *= _pyerrscale.value()*_pyerrscale.value();
  scaledError.fast(3,2) *= _pzerrscale.value()*_pyerrscale.value();
  scaledError.fast(3,3) *= _pzerrscale.value()*_pzerrscale.value();
  beamMomentum.setCovMatrix(scaledError);

  // Tell the conduit to update the constraints that it will pass on to the
  // KalPairSites.  These need to be loaded before the tracks are made, because
  // the KalPairRep ctor adds a PairSite, and it gets it's info from the
  // conduit
  const KalRep* kalRep0 = kalInter0.kalmanRep();
  const KalRep* kalRep1 = kalInter1.kalmanRep();

  // Make sure the original fits converged
  if (!kalRep0->converged() || !kalRep1->converged()) {
    ErrMsg(error) << "KalRep fits didn't converge" << endmsg;
    return AppResult::OK;
  }
  // Make a Conduit to control communication between the tracks
  // Deletion of the conduit is the responsibility of the PairReps
  KalPairConduit* conduit = new KalPairConduit(kalRep0->pieceTraj(),
					       kalRep1->pieceTraj(),
					       beamMomentum);

  // make new tracks
  pairtrack0 = _maker->makePairTrack(*trk0, conduit);
  pairtrack1 = _maker->makePairTrack(*trk1, conduit);
  if( pairtrack0 == 0 || pairtrack1 == 0) {
    ErrMsg(error) << "KalPairFit: failed to make tracks!" << endmsg;
    delete pairtrack0;
    delete pairtrack1;
    return AppResult::OK;
  }

  // make the fits (Should do this through the rep instead?
  KalPairInterface pairInter0;
  KalPairInterface pairInter1;
  if (!pairtrack0->attach(pairInter0, pairtrack0->defaultType()) ||
      !pairtrack1->attach(pairInter1, pairtrack1->defaultType())) {
    ErrMsg(error) << "Failed to attach KalPairInterface!" << endmsg;
    delete pairtrack0;
    delete pairtrack1;
    return AppResult::OK;
  }
  KalPairRep* pairRep0 = const_cast<KalPairRep*>(pairInter0.kalmanPairRep());
  KalPairRep* pairRep1 = const_cast<KalPairRep*>(pairInter1.kalmanPairRep());

  // Only need to explicitly fit one of these reps.  The other will
  // automatically get fit.
  TrkErrCode goodfit = pairRep0->fit();
  if( goodfit.failure() ) {
    ErrMsg(error) << "KalPairFit: error fitting, "
		  << goodfit << endmsg;
    delete pairtrack0;
    delete pairtrack1;
    return AppResult::OK;
  }

  // Compare chi2 of the pair fit with the single track fit
  if ( (pairRep0->chisq() - kalRep0->chisq()) > _deltaChi2Cut.value() ||
       (pairRep1->chisq() - kalRep1->chisq()) > _deltaChi2Cut.value()) {
    delete pairtrack0;
    delete pairtrack1;
    return AppResult::OK;
  }

  // Write out both HepAList and RW list
  kalpairlist->append( pairtrack0 );
  kalpairlist->append( pairtrack1 );

  // Make MergeMap if requested.  If not already there, make it.
  if (_makeMergeMap.value()) {
    TrkMergeMap* mergeMap = 0;
    mergeMap = Ifd<TrkMergeMap>::get(anEvent, _mergeMapKey.value());
    if ( mergeMap == 0 ) {
      mergeMap = new TrkMergeMap();
      if (!Ifd<TrkMergeMap>::put(anEvent, mergeMap, _mergeMapKey.value())) {
	ErrMsg(error) << " Failed to put MergeMap " << _mergeMapKey.value()
		      << " into the event" << endmsg;
	delete mergeMap;
      }
    }
    // Convention is that track 1 is the primary Rep the PairRep was made from
    if (0 != mergeMap) {
      mergeMap->insert(pairtrack0, trk0, trk1);
      mergeMap->insert(pairtrack1, trk1, trk0);
    }
  }

  return AppResult::OK;
}  
//
AppResult
KalPairFit::endJob( AbsEvent* anEvent )
{
  ErrMsg(routine) << name( ) << " end Job" << endmsg;
  if(_maker != 0) delete _maker;
  
  return AppResult::OK;
}

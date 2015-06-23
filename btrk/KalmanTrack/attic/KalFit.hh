//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: KalFit.hh,v 1.46 2005/04/29 18:04:05 desilva Exp $
//
// Description:
//	Class KalFit. This module takes an existing track list and
//      refits the tracks using a kalman filter fit.  The output can be
//      chosen either to overwrite the fit (track rep) of the input
//      list, or to create a (local) copy of the track list.
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//     David Brown  4/17/97
//
// Copyright Information:
//	Copyright (C) 1997		Lawrence Berkeley Laboratory
//
//------------------------------------------------------------------------

#ifndef KALFIT_HH
#define KALFIT_HH

//----------------------
// Base Class Headers --
//----------------------
#include "Framework/AppModule.hh"
#include "AbsParm/AbsParmNamedValue.hh"
#include "AbsParm/AbsParmVector.hh"
#include "Framework/AbsParmBool.hh"
#include "KalmanTrack/KalMiniRep.hh"
#include "KalmanTrack/KalMiniInterface.hh"
#include "KalmanTrack/KalInterface.hh"
#include "TrkBase/TrkDeadInterface.hh"

class KalMaker;
class TrkErrCode;
template <class A> class RecAListCopyParms;

//		---------------------
// 		-- Class Interface --
//		---------------------
class KalFit : public AppModule {

//--------------------
// Instance Members --
//--------------------

public:
// describe what we're fitting
  enum inputType { helix=0,kalman,kalmini};
  enum fitHypos { electron=PdtPid::electron,
		  muon=PdtPid::muon,
		  pion=PdtPid::pion,
		  kaon=PdtPid::kaon,
		  proton=PdtPid::proton,
		  none=100,defhypo,allexisting,all};
  enum fitDirection {fitIn =trkIn,fitOut=trkOut,fitBoth=3};
  
// Constructors
  KalFit( const char* const theName, const char* const theDescription );
  KalFit( const KalFit& other, const char* const theName );
// clone function (with covariant return)
  virtual KalFit* clone(const char* cloneName);
// Destructor
  virtual ~KalFit( );
// Operations
  virtual AppResult           beginJob( AbsEvent* anEvent );
  virtual AppResult      event( AbsEvent* anEvent );
  virtual AppResult           endJob  ( AbsEvent* anEvent );
  
private:
  RecAListCopyParms<TrkRecoTrk>* _copyParms;	// Lists and input/output track conversions
  AbsParmNamedValue<KalFit::fitDirection> _fitdir; // what direction to fit the tracks in
  AbsParmNamedValue<KalFit::inputType> _inputtype; // what's the input
  AbsParmNamedValue<KalMiniRep::miniState> _ministate; // state to build mini in
  AbsParmVector<int> _fithypos; // hypos to build/fit
  AbsParmVector<int> _constraints; // constraints on parameters
  AbsParmBool _maptype; // map requested fit types or not 
  KalMaker* _maker; // fit maker
  mutable KalInterface _kiface; // interfaces
  mutable KalMiniInterface _kminiiface;
  mutable TrkDeadInterface _deadiface;

  bool* _cons; // constraint vector
// utility functions
  bool fitThisHypo(TrkRecoTrk* trk,PdtPid::PidType hypo) const;
  TrkErrCode fitKalRep(KalRep* rep) const;
};
#endif

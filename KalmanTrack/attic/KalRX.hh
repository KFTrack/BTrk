//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: KalRX.hh,v 1.12 2004/08/16 00:58:56 bartoldu Exp $
//
// Description:
//	Class KalRX.  This module tries to 'repair' the tracks that it
//      finds for certain known problems.
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//     David Brown  11/7/01
//
// Copyright Information:
//	Copyright (C) 2001		Lawrence Berkeley Laboratory
//
// Revision History:
//	20030620  M. Kelsey -- Add TCL to remove or delete irreparable tracks
//------------------------------------------------------------------------

#ifndef KALRX_HH
#define KALRX_HH

//----------------------
// Base Class Headers --
//----------------------
#include "Framework/AppModule.hh"
#include "KalmanTrack/KalInterface.hh"
#include "Framework/AbsParmBool.hh"
#include "Framework/AbsParmDouble.hh"
#include "Framework/AbsParmGeneral.hh"
#include "AbsParm/AbsParmIfdStrKey.hh"
#include "AbsParm/AbsParmNamedValue.hh"

#include "TrkBase/TrkHitOnTrkUpdater.hh"
class TrkRecoTrk;
class KalMaker;
class KalRep;
//		---------------------
// 		-- Class Interface --
//		---------------------
class KalRX : public AppModule, TrkHitOnTrkUpdater {

//--------------------
// Instance Members --
//--------------------

public:
// Constructors
  KalRX( const char* const theName, const char* const theDescription );
// Destructor
  virtual ~KalRX( );
// Operations
  virtual AppResult           beginJob( AbsEvent* anEvent );
  virtual AppResult           event( AbsEvent* anEvent );
protected:
  // Value flags for list management
  // "delete" means actually delete the track pointer
  // "remove" means take pointer off list, but don't delete
  enum listAction { noaction=0, deletetrk=1, junktrk=2 };
  enum repairHypos { defhypo=0, all=1};
// repair a fit.  Return code is the fit status (history is added too)
  virtual void repairTrack(TrkRecoTrk* trk) const;
  virtual void changeToElectron(TrkRecoTrk* trk) const;
// utility functions

  void repairFit(TrkRecoTrk* trk,PdtPid::PidType hypo) const;
  TrkErrCode repairStopping(KalRep* rep) const;
  void invalidateUnphysical(TrkRecoTrk* trk, PdtPid::PidType hypo) const;

// statics
  static TrkErrCode _repairstopping;
  static TrkErrCode _refitting;
  static TrkErrCode _electron;
  static TrkErrCode _grazer;
  static TrkErrCode _deactivate;
  static TrkErrCode _unphysical;
  static TrkErrCode _diverged;
  static TrkErrCode _unknown;

private:
  mutable KalInterface _kiface; // interfaces

  AbsParmIfdStrKey _inputkey;		// list key for input tracks
  AbsParmIfdStrKey _junkkey;		// list key for junk tracks
  AbsParmBool _makecurrent; // make non-current tracks current?
  AbsParmBool _ele; // convert to electron types if failed?
  AbsParmNamedValue<listAction> _inputaction; // what to do with bad tracks
  AbsParmNamedValue<repairHypos> _repairhypos; // which hypos to repair
  AbsParmDouble _grazernormcut; // incident angle cut to define a grazer
  AbsParmDouble _matlencut; // length of material intersection to consider 'unphysical'
  AbsParmDouble _omegacut;  // value of omega to consider trackfits 'unphysical'
  AbsParmGeneral<int> _npastcut; // # of hits past a graze to consider it 'false'
protected:
  KalMaker* _maker; // fit maker

};
#endif

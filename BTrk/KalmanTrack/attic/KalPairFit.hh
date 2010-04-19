//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: KalPairFit.hh,v 1.13 2005/04/29 18:04:06 desilva Exp $
//
// Description:
//	Class KalPairFit. This module refits an event with exactly 2 tracks
//      in it to a single track fit, and adds that track to the event track
//      list.
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//     David Brown  4/1/98
//
// Copyright Information:
//	Copyright (C) 1998		Lawrence Berkeley Laboratory
//
//------------------------------------------------------------------------

#ifndef KALPAIRFIT_HH
#define KALPAIRFIT_HH

//----------------------
// Base Class Headers --
//----------------------
#include "Framework/AppModule.hh"
#include "Framework/AbsParmDouble.hh"
#include "Framework/AbsParmBool.hh"
#include "AbsParm/AbsParmIfdStrKey.hh"

class KalPairMaker;
//		---------------------
// 		-- Class Interface --
//		---------------------
class KalPairFit : public AppModule {

//--------------------
// Instance Members --
//--------------------

public:
  // Constructors
  KalPairFit( const char* const theName, const char* const theDescription );
  
  // Destructor
  virtual ~KalPairFit( );
  
  // Clone function
  virtual AppModule* clone(const char* cloneName);
  
  // Operations
  
  virtual AppResult           beginJob( AbsEvent* anEvent );
  virtual AppResult           event   ( AbsEvent* anEvent );
  virtual AppResult           endJob  ( AbsEvent* anEvent );
  //
private:
  AbsParmIfdStrKey _trklistkey; // list key
  AbsParmIfdStrKey _pairlistkey; // pair list key
  AbsParmIfdStrKey _mergeMapKey;
  AbsParmBool _makeMergeMap;
  AbsParmDouble _deltaChi2Cut;
  // Things to allow user input for beam momentum.  The user can specify a 
  // bias in px, py or pz in GeV.  These will be added to the values read
  // from the database.  Useful for studying systematic effects of mismeasured
  // beam parameters on the pair fit.  The user can also specify scale factors
  // for the errors used on the beam momenta.  Can either deweight (scale >1.)
  // or over-weight (scale <1.) the beam momenta information.
  AbsParmDouble _pxbias;  // in GeV, default is 0.
  AbsParmDouble _pybias;  // in GeV, default is 0.
  AbsParmDouble _pzbias;  // in GeV, default is 0.
  AbsParmDouble _pxerrscale;  // Default is 1.
  AbsParmDouble _pyerrscale;  // Default is 1.
  AbsParmDouble _pzerrscale;  // Default is 1.

  KalPairMaker* _maker; // fit maker
};
#endif

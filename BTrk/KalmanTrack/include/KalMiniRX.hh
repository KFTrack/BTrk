//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: KalMiniRX.hh,v 1.5 2004/08/16 00:58:56 bartoldu Exp $
//
// Description:
//	Class KalMiniRX.  This module tries to 'repair' mini tracks
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//     David Brown  4/23/02
//
// Copyright Information:
//	Copyright (C) 2002		Lawrence Berkeley Laboratory
//
//------------------------------------------------------------------------

#ifndef KALMINIRX_HH
#define KALMINIRX_HH

//----------------------
// Base Class Headers --
//----------------------
#include "Framework/AppModule.hh"
#include "KalmanTrack/KalMiniInterface.hh"
#include "KalmanTrack/KalRX.hh"

class KalMiniRep;
class TrkErrCode;
//		---------------------
// 		-- Class Interface --
//		---------------------
class KalMiniRX : public KalRX {

//--------------------
// Instance Members --
//--------------------

public:
// Constructors
  KalMiniRX( const char* const theName, const char* const theDescription );
// Destructor
  virtual ~KalMiniRX( );
// override the repairTrack functions
  virtual void repairTrack(TrkRecoTrk* trk) const;
  virtual void changeToElectron(TrkRecoTrk* trk) const;
private:
  mutable KalMiniInterface _kmiface; // interfaces
  AbsParmBool _downgrade; // downgrade failed fits
  TrkErrCode downgrade(KalMiniRep* krep) const; // re-entrant repair
};
#endif

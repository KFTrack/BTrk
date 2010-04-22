//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: TrkTypeUtil.hh,v 1.1 2000/01/31 20:12:58 echarles Exp $
//
// Description:
//	Class TrkTypeUtil
//
// Environment:
//	Software developed for BaBar expirment @ SLAC B-Factory
//
// Author List:
//      Eric A Charles
//
// Copyright Information:
//	Copyright (C) 1998	Univ. Wisconsin-Madsion
//
//------------------------------------------------------------------------

#ifndef TRKTYPEUTIL_HH
#define TRKTYPEUTIL_HH

//-------------
// C Headers --
//-------------
extern "C" {
}

//---------------
// C++ Headers --
//---------------

#include "BaBar/PdtPid.hh"
class TrkRecoTrk;
class TrkRep;
class TrkFitTypeKey;
class TrkHitList;

//		---------------------
// 		-- Class Interface --
//		---------------------

class TrkTypeUtil {

//------------------
// Static Members --
//------------------

public:

  static const TrkRep* getRep( const TrkRecoTrk& aTrk, 
			       const PdtPid::PidType& pid);

  static const TrkHitList* getHits( const TrkRecoTrk& aTrk, 
				    const PdtPid::PidType& pid);

  static const TrkRep* getRep( const TrkRecoTrk& aTrk,
			       const TrkFitTypeKey& key );

  static const TrkHitList* getHits( const TrkRecoTrk& aTrk, 
				    const TrkFitTypeKey& key);
  
  static PdtPid::PidType pidType( const TrkRecoTrk& aTrk, 
				  const TrkFitTypeKey& key );

//--------------------
// Instance Members --
//--------------------

private:

  // preempt everything
  ~TrkTypeUtil( );

};

// Inline implementations
//#include "TrkTypeUtil.icc"

#endif

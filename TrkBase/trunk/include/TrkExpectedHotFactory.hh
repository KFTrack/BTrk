//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: TrkExpectedHotFactory.hh,v 1.4 2003/01/21 12:55:08 raven Exp $
//
// Description:
//	Class Template
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

#ifndef TRKEXPECTEDHOTFACTORY_HH
#define TRKEXPECTEDHOTFACTORY_HH

//-------------
// C Headers --
//-------------
extern "C" {
}

//---------------
// C++ Headers --
//---------------
#include <memory>

//-------------------------------
// Collaborating Class Headers --
//-------------------------------

#include "TrkBase/TrkDetElemId.hh"

//------------------------------------
// Collaborating Class Declarations --
//------------------------------------

class AbsEvent;
class TrkExpectedHot;
class TrkExpectedHotSet;
class TrkFitTypeKey;
class TrkExpectedMap;
class TrkExpectedTrk;
class TrkRecoTrk;
class TrkHitOnTrk;
class TrkDetElemId;
class TrkFundHit;
class IfdKey;
class GTrack;

template <class K, class V> class AstSTLMap;

//		---------------------
// 		-- Class Interface --
//		---------------------

class TrkExpectedHotFactory {

protected:

  // helper function
  static const TrkExpectedHotSet* parseIntoExTrk(TrkExpectedTrk* exTrk,
						 const TrkHitOnTrk* hot,
						 bool fillGaps = true );

  static const TrkExpectedHotSet* parseSetIntoExTrk(TrkExpectedTrk*,
						    const TrkDetElemId&,
						    bool fillGaps=true);

  // Data members
  
  static const AstSTLMap<const TrkRecoTrk, const TrkFundHit>* _patRecMap;

//--------------------
// Instance Members --
//--------------------

public:

  // Constructors
  TrkExpectedHotFactory( const IfdKey& patRecKey );
  
  // Destructor
  virtual ~TrkExpectedHotFactory( );
  
  // Selectors (const)
  virtual bool makeExHotsFromTraj( TrkExpectedTrk* exTrk,
				   const TrkFitTypeKey& key,
				   const double* theRange) const = 0;
  virtual bool makeExHotsFromHots( TrkExpectedTrk* exTrk,
				   const TrkFitTypeKey& key,
				   const double* theRange,
				   bool redoInt=false,
				   bool fillGaps=true) const = 0;
  virtual bool makeExHotsFromGTrk( TrkExpectedTrk* exTrk,
				   const TrkFitTypeKey& key,
				   const GTrack* gTrk) const = 0;

  virtual bool fillExHotsAfterGTrk( TrkExpectedTrk* exTrk,
				    const TrkFitTypeKey& key ) const;

  virtual bool fillExHot( const TrkHitOnTrk* hot, 
			  TrkExpectedHot* exHot ) const;

  virtual const TrkExpectedHot* parseExHotIntoSet(const TrkHitOnTrk* hot,
						  const TrkFitTypeKey& key,
						  TrkExpectedHotSet& hotSet,
						  bool reCalcInt=false) const = 0;
  bool parseHotUsingMap( TrkExpectedTrk& exTrk,
			 const TrkExpectedMap* exMap,
			 const TrkHitOnTrk* hot, 
			 const bool fillGaps,
			 const bool redoInt) const;
  
  bool fillHotSetUsingMap( const TrkExpectedMap* exMap,
                           TrkExpectedTrk* exTrk,
                           TrkExpectedHotSet* hotSet ) const;

  bool parseHotsFromMap( const TrkExpectedMap* ) const;
  bool parseHotsFromMap( const TrkExpectedMap*, 
                         const TrkFitTypeKey& key ) const;

  virtual bool nextEvent( AbsEvent* anEvent );

  virtual bool usedHitInPatRec( const TrkRecoTrk*, const TrkFundHit* ) const;
  virtual bool usedInPatRec( const TrkRecoTrk*, const TrkExpectedHot* ) const = 0;

  virtual int nUsedInPatRec( const TrkRecoTrk*) const;

private:
  std::auto_ptr<IfdKey> _patRecKey;

};

#endif

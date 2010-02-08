//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: TrkExpectedMap.hh,v 1.4 2004/08/06 06:31:41 bartoldu Exp $
//
// Description:
//	Class TrkExpectedMap.  Class for mapping sets of 
//   TrkRecoTrk's into TrkExpectedTrks.  Also handles access to 
//   TrkExpectedTrks 
//
// Environment:
//	Software developed for BaBar expirment @ SLAC B-Factory
//
// Author List:
//      Eric A Charles
//
// Copyright Information:
//	Copyright (C) 1999	Univ. Wisconsin-Madsion
//
//------------------------------------------------------------------------

#ifndef TRKEXPECTEDMAP_HH
#define TRKEXPECTEDMAP_HH

//-------------
// C Headers --
//-------------
extern "C" {
}

//---------------
// C++ Headers --
//---------------

#include <iostream>

//-----------------
// BaBar Headers --
//-----------------

//----------------------
// Base Class Headers --
//----------------------

#include "AbsEvent/AbsEvtObj.hh"

//-------------------------------
// Collaborating Class Headers --
//-------------------------------

#include <map>
#include <vector>
#include "PDT/PdtPid.hh"

//------------------------------------
// Collaborating Class Declarations --
//------------------------------------

class AbsEvent;
class TrkRep;
class TrkExpectedTrk;
class TrkDetElemId;
class TrkHitOnTrk;
class TrkExpectedHotSet;
class TrkExpectedHot;
class TrkRecoTrk;
class TrkFitTypeKey;
class TrkExpectedHotFactory;
class TrkRep;
class GTrack;

//		---------------------
// 		-- Class Interface --
//		---------------------

class TrkExpectedMap : public AbsEvtObj {

  //--------------------
  // Static Members   --
  //--------------------

public:  

  // make the expected tracks
  static TrkExpectedTrk* makeExpectedTrk( const TrkRep* aRep,
                                          TrkExpectedMap* map,
                                          const PdtPid::PidType pid=PdtPid::null );
public:

  TrkExpectedMap( TrkExpectedHotFactory* factory,
                  AbsEvent* anEvent );

  virtual ~TrkExpectedMap( );

  bool parseInFitType( TrkExpectedTrk* exTrk,
                       const TrkFitTypeKey& key,
                       const double* theRange,
                       bool fromTraj=false,
                       bool redoInt=false);

  bool parseFitTypeFromHots( TrkExpectedTrk* exTrk,
                             const TrkFitTypeKey& key,
                             const double* theRange,
                             bool fillGaps=true,
                             bool redoInt=false);

  bool parseExHotsFromGTrk( TrkExpectedTrk* exTrk, const TrkFitTypeKey& key );

  // Selectors (const) 
  // return the exHot factory
  const TrkExpectedHotFactory* exHotFactory() const { return _hotFactory; }

  // return the expected track for a particular object,
  // all of these use std::map.
  const TrkExpectedTrk* getExpectedTrk( const TrkRecoTrk* aTrk ) const;
  const TrkExpectedTrk* getExpectedTrk( const TrkHitOnTrk* aHot ) const;
  const TrkExpectedTrk* getExpectedTrk( const TrkRep* aRep ) const;

  // get the associated gTrack
  const GTrack* getGTrack( const TrkRep* aRep ) const;

  // get the fit type for a particular object 
  TrkFitTypeKey getFitType( const TrkHitOnTrk* aHot ) const;
  TrkFitTypeKey getFitType( const TrkRecoTrk* aTrk,
                            const PdtPid::PidType pidType=PdtPid::null) const;
  TrkFitTypeKey getFitType( const TrkRep* aRep ) const;

  // return the list of expected tracks
  const std::vector<TrkExpectedTrk*>& exTracks() const { return _exTrks; };

  // fill a list of all the expected hot or hot sets
  void fillHotSetList( std::vector<TrkExpectedHotSet*>& hotSets ) const;
  void fillHotsList( const TrkFitTypeKey& key, 
                     std::vector<TrkExpectedHot*>& hotSets ) const;

  // Modifiers
  // add a rep to the map, all reps should be added before 
  // makeReferenceTracks and crossLink are called
  bool addRep( const TrkRecoTrk* aTrk, const TrkFitTypeKey& key, const GTrack* gTrk = 0);
  bool addRep( const TrkRep* aRep, const TrkFitTypeKey& key, const GTrack* gTrk = 0);

  bool addGTrkLink( const GTrack* gTk, const TrkExpectedTrk* exTrk );


  // build the original TrkExpectedHotSets based on all the
  // reps that are of type refKey.  
  // if doEff is true also make TrkExpectedHotSets for elements
  // that were hit, but don't have associated TrkHitOnTrks
  bool makeReferenceTracks( const TrkFitTypeKey& refKey,
                            const bool doIntersect = false,
                            const bool doEff = false );

   // build the original TrkExpectedHotSets based on all the
  // reps that are of type refKey.  
  // if doEff is true also make TrkExpectedHotSets for elements
  // that were hit, but don't have associated TrkHitOnTrks
  bool makeReferenceTracksFromMC( const TrkFitTypeKey& refKey );

  // parse in all the other (non-reference) TrkReps in the map
  bool crossLink(const bool fillGaps = false,
                 const bool doIntersect = false);

  bool crossLinkFromMC();

  // these two functions make expected tracks and set the reference type
  TrkExpectedTrk* makeExpectedTrk( const TrkRecoTrk* aTrk,
                                   const PdtPid::PidType pidType=PdtPid::null);

  TrkExpectedTrk* makeExpectedTrk( const TrkRep* aRep );

  // print operators
  void printAll( std::ostream& os = std::cout ) const;
  void printIdMap( std::ostream& os = std::cout ) const;
  void printRepMap( std::ostream& os = std::cout ) const;

private:

  // unowned pointer
  const TrkExpectedHotFactory* _hotFactory;

  // Data members    
  // List of the expected tracks
  std::vector<TrkExpectedTrk*> _exTrks;
  // Map from reps to fit types
  std::map< const TrkRep*,TrkFitTypeKey> _fitTypeMap;
  // Map from reps to GTracks
  std::map< const TrkRep*, const GTrack*> _gTkMap;
  // Map from GTracks to expected tracks
  std::map< const GTrack*, const TrkExpectedTrk*> _gTkExMap;
  // Map from track id's to expected tracks
  std::map< int, const TrkExpectedTrk*> _exTrkMap;
};

#endif

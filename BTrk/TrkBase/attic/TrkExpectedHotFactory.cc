//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: TrkExpectedHotFactory.cc,v 1.8 2004/09/10 18:00:17 bartoldu Exp $
//
// Description:
//	Class TrkExpectedHotFactory
//
// Environment:
//	Software developed for BaBar expirment @ SLAC B-Factory
//
// Author List:
//	Eric A Charles
//
// Copyright Information:
//	Copyright (C) 1998	Univ. Wisconsin-Madison
//
//------------------------------------------------------------------------

//----------------
// BaBar header
//----------------
#include "BaBar/BaBar.hh"

//-----------------------
// This Class's Header --
//-----------------------
#include "TrkBase/TrkExpectedHotFactory.hh"

//-------------
// C Headers --
//-------------
extern "C" {
}

//---------------
// C++ Headers --
//---------------

#include <map>
#include <vector>
#include <algorithm>

//-------------------------------
// Collaborating Class Headers --
//-------------------------------

#include "TrkBase/TrkExpectedTrk.hh"
#include "TrkBase/TrkExpectedHotSet.hh"
#include "TrkBase/TrkExpectedHot.hh"
#include "TrkBase/TrkExpectedMap.hh"
#include "TrkBase/TrkFitTypeKey.hh"
#include "TrkBase/TrkDetElemId.hh"
#include "TrkBase/TrkHitOnTrk.hh"
#include "TrkBase/TrkFundHit.hh"
#include "TrkBase/TrkRecoTrk.hh"

#include "ProxyDict/IfdKey.hh"
#include "ProxyDict/Ifd.hh"
#include "AbsEvent/AbsEvent.hh"
#include "AssocTools/AstSTLMap.hh"


//-----------------------------------------------------------------------
// Local Macros, Typedefs, Structures, Unions and Forward Declarations --
//-----------------------------------------------------------------------



//		-----------------------------------------------
// 		-- Static Data & Function Member Definitions --
//		-----------------------------------------------

const AstSTLMap<const TrkRecoTrk, const TrkFundHit>*
TrkExpectedHotFactory::_patRecMap(0);


TrkExpectedHotFactory::TrkExpectedHotFactory(const IfdKey& patRecKey)
  :_patRecKey(patRecKey.clone())
{
}

//--------------
// Destructor --
//--------------

TrkExpectedHotFactory::~TrkExpectedHotFactory( )
{
}

bool 
TrkExpectedHotFactory::parseHotUsingMap( TrkExpectedTrk& exTrk,
                                         const TrkExpectedMap* exMap,
                                         const TrkHitOnTrk* hot, 
                                         const bool fillGaps,
                                         const bool redoInt ) const
{
  // get the fit type key
  if ( exMap == 0 || hot == 0 ) return false;
  TrkFitTypeKey fitTypeKey = exMap->getFitType( hot );
  if ( fitTypeKey.value() < 0 ) return false;
  
  // get the hot set or make a new one
  const TrkExpectedHotSet* hotSet = exTrk.exHotSet(hot);
  TrkExpectedHotSet* theSet(0);
  if ( hotSet == 0 ) {
    if ( !fillGaps ) return false;
    const TrkDetElemId id = hot->hit()->elemId();
    TrkExpectedHotSet* theSet = 
      new TrkExpectedHotSet( id, TrkFitTypeKey::currentKey() ); 
    exTrk.addExHotSet( id, theSet );
  }
  else { 
    theSet = (TrkExpectedHotSet*)hotSet;
  }
  
  // make a new ex hot and put it in the set
  //const TrkExpectedHot* exHot = 
  parseExHotIntoSet( hot,fitTypeKey,*theSet,redoInt );

  exTrk._reps.push_back( hot->getParentRep() );
  return true;
}

bool 
TrkExpectedHotFactory::fillHotSetUsingMap( const TrkExpectedMap* exMap,
                                           TrkExpectedTrk* exTrk,
                                           TrkExpectedHotSet* hotSet ) const
{
  if ( hotSet == 0 || exMap == 0 || exTrk == 0 ) return false;
  // get the other hots associated with same fund hit
  std::vector<TrkHitOnTrk*> hots; 
  hotSet->getHots(hots);
  // loop on the hots  
  for ( int i(0); i < hots.size(); ++i ) {
    const TrkHitOnTrk* hot = hots[i];
    // get the fit type and see if it is in the map
    TrkFitTypeKey fitTypeKey =  exMap->getFitType( hot );
    if ( fitTypeKey.value() < 0 ) continue;
    // see if we already have exHot of this type
    if ( hotSet->hasFitType( fitTypeKey ) ) continue;
    // make the ex hot and put it in the set
    //const TrkExpectedHot* exHot = 
    parseExHotIntoSet( hot, fitTypeKey, *hotSet );    
    // put in the the rep list, even if allready in, 
    // because that is how we test the number of overlaps
    exTrk->_reps.push_back( hot->getParentRep() );
  }  
  return true;
}

bool 
TrkExpectedHotFactory::fillExHot( const TrkHitOnTrk* hot, 
                                  TrkExpectedHot* exHot ) const
{
  if ( hot == 0 || exHot == 0 ) return false;
  exHot->setHot(hot);
  return true;
}

bool
TrkExpectedHotFactory::fillExHotsAfterGTrk( TrkExpectedTrk* exTrk,
                                            const TrkFitTypeKey& key ) const
{
  if ( exTrk == 0 ) return false;
  const TrkRecoTrk* theTrk = exTrk->track();
  TrkHitList::hot_iterator end = theTrk->hits()->end();
  for (TrkHitList::hot_iterator i = theTrk->hits()->begin();
       i != end ;
       ++i ) {
    TrkExpectedHot* exHot = const_cast<TrkExpectedHot*>(exTrk->exHot(i.get()));
    if ( exHot == 0 ) {
      const TrkDetElemId id = i->hit()->elemId();
      TrkExpectedHotSet* newHotSet =
        const_cast<TrkExpectedHotSet*>(parseSetIntoExTrk(exTrk,id));
      //const TrkExpectedHot* exHot = 
      parseExHotIntoSet(i.get(),key,*newHotSet);
    } else {
      fillExHot(i.get(),exHot);
    }
  }
  return true;
}


const TrkExpectedHotSet*
TrkExpectedHotFactory::parseIntoExTrk(TrkExpectedTrk* exTrk,
                                      const TrkHitOnTrk* hot,
                                      bool fillGaps )
{
  if ( exTrk == 0 || hot == 0 ) return 0;
  TrkDetElemId id = hot->hit()->elemId();
  return parseSetIntoExTrk(exTrk,id,fillGaps);
}

const TrkExpectedHotSet* 
TrkExpectedHotFactory::parseSetIntoExTrk(TrkExpectedTrk* exTrk,
                                         const TrkDetElemId& id,
                                         bool fillGaps)
{
  // is there an exHot set for this element yet?
  // if so, get it, if not make it 
  TrkExpectedHotSet* anExHotSet = (TrkExpectedHotSet*)(exTrk->exHotSet(id));
  if ( anExHotSet == 0 && fillGaps ) {
    anExHotSet = new TrkExpectedHotSet( id, TrkFitTypeKey::currentKey() );
    exTrk->addExHotSet(id,anExHotSet);
  }
  return anExHotSet;
}

bool 
TrkExpectedHotFactory::parseHotsFromMap( const TrkExpectedMap* map) const
{
  if ( map == 0 ) return false;
  typedef std::vector<TrkExpectedTrk*> vec_t;
  const vec_t& exTrks = map->exTracks();
  for ( vec_t::const_iterator i=exTrks.begin(); i!=exTrks.end(); ++i ) {
    if ( ! (*i)->parseHotsFromMap() ) continue;
  }
  return true;
}

bool 
TrkExpectedHotFactory::parseHotsFromMap( const TrkExpectedMap* map, 
                                         const TrkFitTypeKey& key ) const
{
  if ( map == 0 ) return false;
  typedef std::vector<TrkExpectedTrk*> vec_t;
  const vec_t& exTrks = map->exTracks();
  for ( vec_t::const_iterator i=exTrks.begin(); i!=exTrks.end(); ++i ) {
    if ( ! (*i)->parseHotsFromMap(key) ) continue;
  }
  return true;
}

bool
TrkExpectedHotFactory::nextEvent( AbsEvent* anEvent )
{
  if ( _patRecKey.get() == 0 ) return false;
  _patRecMap = Ifd<AstSTLMap<const TrkRecoTrk, const TrkFundHit> >::get(anEvent,*_patRecKey);
  return ( _patRecMap != 0 );
}


bool 
TrkExpectedHotFactory::usedHitInPatRec( const TrkRecoTrk* trk, 
                                        const TrkFundHit* hit) const
{
  bool found = false;
  std::map<const TrkRecoTrk*, std::vector<const TrkFundHit*> >::const_iterator
    iter = _patRecMap->find(trk);
  if (iter != _patRecMap->end()) {
    const std::vector<const TrkFundHit*>* hits = & iter->second;
    std::vector<const TrkFundHit*>::const_iterator hitIter = hits->begin();
    while (! found && hitIter != hits->end() ) {
      if ( **hitIter == *hit ) found = true;
      hitIter++;
    }
  }

  return found;
}
 
int 
TrkExpectedHotFactory::nUsedInPatRec( const TrkRecoTrk* aTrk) const{
  if ( aTrk == 0 ) return 0;
  std::map<const TrkRecoTrk*, std::vector<const TrkFundHit*> >::const_iterator
    iter = _patRecMap->find(aTrk);
  if (iter != _patRecMap->end()) {
    const std::vector<const TrkFundHit*>* hits = & iter->second;
    return hits->size();
  }
  return 0;
}

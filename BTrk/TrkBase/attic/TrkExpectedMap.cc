//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: TrkExpectedMap.cc,v 1.6 2004/09/10 18:00:17 bartoldu Exp $
//
// Description:
//	Class TrkExpectedMap
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
#include "TrkBase/TrkExpectedMap.hh"

//-------------
// C Headers --
//-------------
extern "C" {
}

//---------------
// C++ Headers --
//---------------

//-------------------------------
// Collaborating Class Headers --
//-------------------------------

#include "BbrStdUtils/BbrCollectionUtils.hh"
#include "TrkBase/TrkExpectedTrk.hh"
#include "TrkBase/TrkHitOnTrk.hh"
#include "TrkBase/TrkRecoTrk.hh"
#include "TrkBase/TrkFitTypeKey.hh"
#include "TrkBase/TrkHitList.hh"
#include "TrkBase/TrkRep.hh"
#include "TrkBase/TrkTypeUtil.hh"
#include "TrkBase/TrkExpectedHotFactory.hh"
using std::endl;
using std::ostream;


//-----------------------------------------------------------------------
// Local Macros, Typedefs, Structures, Unions and Forward Declarations --
//-----------------------------------------------------------------------



//		-----------------------------------------------
// 		-- Static Data & Function Member Definitions --
//		-----------------------------------------------

TrkExpectedTrk*
TrkExpectedMap::makeExpectedTrk( const TrkRep* aRep,
				 TrkExpectedMap* map,
				 const PdtPid::PidType pid ){
  const GTrack* gTk(0);
  if ( map ) {
    gTk = map->getGTrack(aRep);
  }
  TrkExpectedTrk* exTrk = new TrkExpectedTrk(aRep->parentTrack(),pid,gTk,(TrkExpectedMap*)map);
  if ( gTk != 0 && map != 0 ) {
    map->addGTrkLink(gTk,exTrk);
  }
  return exTrk;
}

bool
TrkExpectedMap::parseInFitType( TrkExpectedTrk* exTrk,
				const TrkFitTypeKey& key,
				const double* theRange,
				bool fromTraj,
				bool redoInt) {
  if ( exHotFactory() == 0 ) return false;
  return fromTraj ?
    exHotFactory()->makeExHotsFromTraj(exTrk,key,theRange) :
    exHotFactory()->makeExHotsFromHots(exTrk,key,theRange,redoInt);
}



bool
TrkExpectedMap::parseFitTypeFromHots( TrkExpectedTrk* exTrk,
				      const TrkFitTypeKey& key,
				      const double* theRange,
				      bool fillGaps,
				      bool redoInt){
  if ( exHotFactory() == 0 ) return false;
  return  exHotFactory()->makeExHotsFromHots(exTrk,key,theRange,
					     fillGaps,redoInt);
}

bool
TrkExpectedMap::parseExHotsFromGTrk( TrkExpectedTrk* exTrk, const TrkFitTypeKey& key ){
  if ( exHotFactory() == 0 ) return false;
  if ( exTrk == 0 ) return false;
  const GTrack* gTrk = exTrk->gTrack();
  if ( gTrk == 0 ) return false;
  return exHotFactory()->makeExHotsFromGTrk(exTrk, key, gTrk);
}

//		----------------------------------------
// 		-- Public Function Member Definitions --
//		----------------------------------------

//----------------
// Constructors --
//----------------

TrkExpectedMap::TrkExpectedMap( TrkExpectedHotFactory* factory,
                                AbsEvent* anEvent )
  : _hotFactory(factory)
{
     _exTrks.reserve(10);
     assert(_hotFactory!=0);
     factory->nextEvent(anEvent);
}

//--------------
// Destructor --
//--------------

TrkExpectedMap::~TrkExpectedMap()
{
  // _exTrks are owned, so delete
  std::for_each(_exTrks.begin(),_exTrks.end(),babar::Collection::DeleteObject());
}

//-------------
// Operators --
//-------------

//-------------
// Selectors --
//-------------

const TrkExpectedTrk*
TrkExpectedMap::getExpectedTrk( const TrkRecoTrk* aTrk ) const
{
  typedef std::map< int, const TrkExpectedTrk*> map_t;
  const int id = aTrk->id();
  map_t::const_iterator i = _exTrkMap.find(id);
  return (i==_exTrkMap.end())?0:i->second;
}

const TrkExpectedTrk*
TrkExpectedMap::getExpectedTrk( const TrkRep* aRep ) const
{
  if ( aRep == 0 ) return 0;
  const TrkRecoTrk* aTrk = aRep->parentTrack();
  return getExpectedTrk( aTrk );
}

const TrkExpectedTrk*
TrkExpectedMap::getExpectedTrk( const TrkHitOnTrk* aHot ) const
{
  if ( aHot == 0 ) return 0;
  const TrkRecoTrk* aTrk = aHot->parentTrack();
  return getExpectedTrk( aTrk );
}

const GTrack*
TrkExpectedMap::getGTrack( const TrkRep* aRep ) const
{
  typedef std::map< const TrkRep*, const GTrack*> map_t;
  map_t::const_iterator i = _gTkMap.find(aRep);
  return i == _gTkMap.end() ? 0 : i->second;
}

TrkFitTypeKey
TrkExpectedMap::getFitType( const TrkHitOnTrk* aHot ) const
{
  const TrkRep* aRep = aHot->getParentRep();
  return aRep == 0 ? TrkFitTypeKey((char*)0) : getFitType(aRep);
}

TrkFitTypeKey
TrkExpectedMap::getFitType( const TrkRecoTrk* aTrk,
                            const PdtPid::PidType pidType ) const
{
  const TrkRep* theRep = TrkTypeUtil::getRep(*aTrk,pidType);
  return getFitType( theRep );
}

TrkFitTypeKey
TrkExpectedMap::getFitType( const TrkRep* aRep ) const
{
  typedef std::map< const TrkRep*,TrkFitTypeKey> map_t;
  map_t::const_iterator i = _fitTypeMap.find(aRep);
  return i == _fitTypeMap.end() ? 0 : i->second;
}

void
TrkExpectedMap::fillHotSetList( std::vector<TrkExpectedHotSet*>& hotSets ) const
{
  typedef std::vector<TrkExpectedTrk*>::const_iterator iter;
  for ( iter i = _exTrks.begin(); i < _exTrks.end(); ++i ) {
    if ( *i == 0 ) continue;
    (*i)->fillExHotSets(hotSets);
  }
}

void
TrkExpectedMap::fillHotsList( const TrkFitTypeKey& key,
                              std::vector<TrkExpectedHot*>& hots ) const
{
  typedef std::vector<TrkExpectedTrk*>::const_iterator iter;
  for ( iter i = _exTrks.begin(); i < _exTrks.end(); ++i ) {
    if ( *i == 0 ) continue;
    (*i)->fillExHots(key,hots);
  }
}

bool
TrkExpectedMap::addRep( const TrkRecoTrk* aTrk,
                        const TrkFitTypeKey& key,
                        const GTrack* gTrk )
{
  const TrkRep* aRep = TrkTypeUtil::getRep(*aTrk,key);
  return addRep(aRep,key,gTrk);
}

bool
TrkExpectedMap::addRep( const TrkRep* aRep,
                        const TrkFitTypeKey& key,
                        const GTrack* gTrk )
{
  if ( aRep == 0 ) return false;
  typedef std::map< const TrkRep*,TrkFitTypeKey> map_t;
  map_t::iterator i = _fitTypeMap.find(aRep);
  if ( i != _fitTypeMap.end()) return false;
  _fitTypeMap.insert(std::make_pair(aRep,key));
  _gTkMap.insert(std::make_pair(aRep,gTrk));
  return true;
}

bool
TrkExpectedMap::addGTrkLink( const GTrack* gTk, const TrkExpectedTrk* exTrk )
{
  _gTkExMap.insert(std::make_pair(gTk,exTrk));
  return true;
}


bool
TrkExpectedMap::crossLink(const bool fillGaps,
                          const bool doIntersect)
{
  // get all the extra hots
  exHotFactory()->parseHotsFromMap(this);

  // now make the cross listings
  typedef std::map< const TrkRep*,TrkFitTypeKey> map_t;
  for ( map_t::iterator itr = _fitTypeMap.begin(); itr!=_fitTypeMap.end();++itr ) {
    const TrkRecoTrk* aTrk = itr->first->parentTrack();
    if ( aTrk == 0 ) continue;
    int id = aTrk->id();
    if (_exTrkMap.find(id)!=_exTrkMap.end()) continue;
    int bestOverlap(-1);
    TrkExpectedTrk* exTrk(0);
    typedef std::vector<TrkExpectedTrk*>::iterator iter;
    for ( iter iEx = _exTrks.begin(); iEx < _exTrks.end(); ++iEx ) {
      const int test = (*iEx)->overLap(itr->first);
      if ( test <= bestOverlap ) continue;
      bestOverlap = test;
      exTrk = *iEx;
    }
    if ( exTrk == 0 ) continue;
    const TrkExpectedTrk *cexTrk = exTrk;
    _exTrkMap.insert(std::make_pair(id,cexTrk));
    const GTrack* gTk = exTrk->gTrack();
    if ( gTk != 0 ) {
      _gTkExMap.insert(std::make_pair(gTk,cexTrk));
    }

    // now link in all the hits from this rep
    const TrkHitList* hitList = TrkTypeUtil::getHits( *aTrk, itr->second );
    if ( hitList == 0 ) continue;
    if ( ! fillGaps ) continue;
    TrkHitList::hot_iterator end = hitList->end();
    for ( TrkHitList::hot_iterator iHit = hitList->begin(); iHit != end; ++iHit ) {
      if ( ! exHotFactory()->parseHotUsingMap(*exTrk,this,iHit.get(),
                                              fillGaps,doIntersect) ) continue;
    }
  }
  return true;
}

bool
TrkExpectedMap::crossLinkFromMC()
{
  // get all the extra hots
  exHotFactory()->parseHotsFromMap(this);
  typedef std::map< const TrkRep*,TrkFitTypeKey> map_t;
  for (map_t::iterator itr = _fitTypeMap.begin(); itr!=_fitTypeMap.end();++itr ) {
    const TrkRecoTrk* aTrk = itr->first->parentTrack();
    if ( aTrk == 0 ) continue;
    int id = aTrk->id();
    typedef std::map< int, const TrkExpectedTrk*> map2_t;
    map2_t::iterator i = _exTrkMap.find(id);
    if (i!=_exTrkMap.end()) continue;
    typedef std::map< const TrkRep*, const GTrack*> map3_t;
    map3_t::iterator j = _gTkMap.find(itr->first);
    const GTrack* gTk =  ( j==_gTkMap.end()?0:j->second);
    if ( gTk == 0 ) continue;
    typedef std::map< const GTrack*, const TrkExpectedTrk*> map4_t;
    map4_t::iterator k = _gTkExMap.find(gTk);
    const TrkExpectedTrk* exTrk = (k==_gTkExMap.end()?0:k->second);
    _exTrkMap.insert(std::make_pair(id,exTrk));
  }
  return true;
}

bool
TrkExpectedMap::makeReferenceTracks( const TrkFitTypeKey& refKey,
                                     const bool doIntersect,
                                     const bool doEff )
{
  typedef std::map< const TrkRep*,TrkFitTypeKey> map_t;
  for(map_t::iterator itr=_fitTypeMap.begin();itr!=_fitTypeMap.end();++itr){
    if ( itr->second.value() == refKey.value() ) {
      TrkExpectedTrk* newExTrk = makeExpectedTrk( itr->first );
      bool madeHots = parseInFitType(newExTrk,refKey,0,doEff,doIntersect);
      if ( !madeHots ) return false;
    }
  }
  return true;
}

bool
TrkExpectedMap::makeReferenceTracksFromMC( const TrkFitTypeKey& refKey )
{
  const TrkExpectedHotFactory* exHotFact = exHotFactory();
  typedef std::map< const TrkRep*,TrkFitTypeKey> map_t;
  for (map_t::iterator itr=_fitTypeMap.begin();itr!=_fitTypeMap.end();++itr){
    if ( itr->second.value() == refKey.value() ) {
      TrkExpectedTrk* newExTrk = makeExpectedTrk( itr->first );
      const GTrack* gTk = newExTrk->gTrack();
      if ( exHotFact && gTk ) {
        exHotFact->makeExHotsFromGTrk(newExTrk,refKey,gTk);
        exHotFact->fillExHotsAfterGTrk(newExTrk,refKey);
      }
    }
  }
  return true;
}

TrkExpectedTrk*
TrkExpectedMap::makeExpectedTrk( const TrkRecoTrk* aTrk,
                                 const PdtPid::PidType pidType)
{
  return makeExpectedTrk(TrkTypeUtil::getRep(*aTrk,pidType));
}

TrkExpectedTrk*
TrkExpectedMap::makeExpectedTrk( const TrkRep* aRep)
{
  if ( aRep == 0 ) return 0;
  const TrkFitTypeKey key = getFitType(aRep);
  if ( key.value() < 0 ) return 0;
  const int tkId = aRep->parentTrack()->id();
  typedef std::map< int, const TrkExpectedTrk*> map_t;
  map_t::iterator i = _exTrkMap.find(tkId);
  if ( i!=_exTrkMap.end() ) return 0;
  // pass this along to static member which calls TrkExpectedTrk c'tor
  TrkExpectedTrk* exTrk = makeExpectedTrk( aRep,
                                           this,
                                           key.pidType() );
  assert(exTrk!=0);

  // book-keeping
  _exTrkMap.insert(std::make_pair(tkId,(const TrkExpectedTrk*)exTrk));
  _exTrks.push_back(exTrk);
  return exTrk;

}

void
TrkExpectedMap::printAll( ostream& os ) const
{
  printIdMap(os);
  printRepMap(os);
  // loop on fit types
  for ( int i(0); i < TrkFitTypeKey::currentKey(); ++i ) {
    TrkFitTypeKey key(i);
    os << "TrkExpectedMap info for fit type: " << key << endl;
    // loop on tracks
    typedef std::vector<TrkExpectedTrk*>::const_iterator iter;
    for ( iter iTk = _exTrks.begin(); iTk < _exTrks.end(); ++iTk ) {
      os << "Starting track number " << (iTk-_exTrks.begin()) << " ptr = "
         << *iTk << endl;
      (*iTk)->printFit(key,os);
    }
    os << endl;
  }
}

void
TrkExpectedMap::printIdMap( ostream& os ) const
{
  typedef std::map<int, const TrkExpectedTrk*> map_t;
  for (map_t::const_iterator itr = _exTrkMap.begin();itr!=_exTrkMap.end();++itr) {
    int tkId = itr->first;
    const TrkExpectedTrk* exTrk = itr->second;
    os << "Track " << tkId << " goes to expected track " << exTrk << endl;
  }
}

void
TrkExpectedMap::printRepMap( ostream& os ) const
{
  typedef std::map<const TrkRep*,TrkFitTypeKey>  map_t;
  for(map_t::const_iterator itr = _fitTypeMap.begin(); itr!=_fitTypeMap.end();++itr) {
    const TrkRep* rep = itr->first;
    const TrkRecoTrk* trk = rep->parentTrack();
    const int& tkId = trk->id();
    TrkFitTypeKey key = itr->second;
    os << "Rep " << rep << " from track " << tkId
       << " is of type " << key << endl;
  }
}

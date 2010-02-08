//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: TrkExpectedTrk.cc,v 1.8 2004/09/10 18:00:17 bartoldu Exp $
//
// Description:
//	Class TrkExpectedTrack.  Stores lists of expected hits for a given
//    track in the Trk.
//
// Environment:
//	Software developed for Babar experiment @ Slac B-factory.
//
// Author List:
//	Eric A Charles
//
//------------------------------------------------------------------------

//----------------
// BaBar header
//----------------
#include "BaBar/BaBar.hh"

//-----------------------
// This Class's Header --
//-----------------------
#include "TrkBase/TrkExpectedTrk.hh"

//-------------------------------
// Collaborating Class Headers --
//-------------------------------

#include "BbrStdUtils/BbrPairUtils.hh"
#include "TrkBase/TrkFitTypeKey.hh"
#include "TrkBase/TrkExpectedHotSet.hh"
#include "TrkBase/TrkExpectedHot.hh"
#include "TrkBase/TrkExpectedMap.hh"
#include "TrkBase/TrkHitOnTrk.hh"
#include "TrkBase/TrkRecoTrk.hh"
#include "TrkBase/TrkRep.hh"
#include "TrkBase/TrkTypeUtil.hh"
#include "TrkBase/TrkFundHit.hh"
#include "TrkBase/TrkExpectedHotFactory.hh"
using std::ostream;

//-----------------------------------------------------------------------
// Local Macros, Typedefs, Structures, Unions and Forward Declarations --
//-----------------------------------------------------------------------



//		-----------------------------------------------
// 		-- Static Data & Function Member Definitions --
//		-----------------------------------------------


//		----------------------------------------
// 		-- Public Function Member Definitions --
//		----------------------------------------

//----------------
// Constructors --
//----------------


TrkExpectedTrk::TrkExpectedTrk ( const TrkRecoTrk* trk,
                                 const PdtPid::PidType pidType,
                                 const GTrack* gTrk,
                                 TrkExpectedMap* exMap)
  :_parentMap(exMap),
   _refRep(TrkTypeUtil::getRep(*trk,pidType)),
   _gTrk(gTrk)
{
}

TrkExpectedTrk::~TrkExpectedTrk( )
{
  // This class owns expected Hots.
  typedef std::map<TrkDetElemId, TrkExpectedHotSet*> map_t;
  for (map_t::const_iterator itr=_exHotTable.begin();itr!=_exHotTable.end();++itr) {
    delete itr->second;
  }
}

const TrkRecoTrk*
TrkExpectedTrk::track() const
{
  return _refRep == 0 ? 0 : _refRep->parentTrack();
}

const TrkExpectedHotSet*
TrkExpectedTrk::exHotSet( const TrkDetElemId& elemId) const
{
  typedef std::map<TrkDetElemId,TrkExpectedHotSet*> map_t;
  map_t::const_iterator i = _exHotTable.find(elemId);
  return (i==_exHotTable.end()?0:i->second);
}

const TrkExpectedHotSet*
TrkExpectedTrk::exHotSet( const TrkHitOnTrk* hot ) const
{
  if ( hot == 0 ) return 0;
  const TrkDetElemId id = hot->hit()->elemId();
  return exHotSet( id );
}

const TrkExpectedHotSet*
TrkExpectedTrk::exHotSet( const TrkFundHit* hit ) const
{
  if ( hit == 0 ) return 0;
  const TrkDetElemId id = hit->elemId();
  return exHotSet( id );
}


const TrkExpectedHot*
TrkExpectedTrk::exHot( const TrkDetElemId& elemId,
                       const TrkFitTypeKey key ) const
{
  const TrkExpectedHotSet* hotSet = exHotSet(elemId);
  return hotSet != 0 ? hotSet->exHot( key ): 0;
}

const TrkExpectedHot*
TrkExpectedTrk::exHot( const TrkHitOnTrk* hot ) const
{
  const TrkExpectedHotSet* hotSet = exHotSet(hot);
  const TrkFitTypeKey key =
    _parentMap == 0 ? TrkFitTypeKey() : _parentMap->getFitType( hot );
  return hotSet != 0 ? hotSet->exHot( key ): 0;
}

const TrkExpectedHot*
TrkExpectedTrk::exHot( const TrkFundHit* hit,
                       const TrkFitTypeKey key ) const
{
  const TrkExpectedHotSet* hotSet = exHotSet(hit);
  return hotSet != 0 ? hotSet->exHot( key ) : 0;
}

void
TrkExpectedTrk::fillExHotSets( std::vector<TrkExpectedHotSet*>& hotSets ) const
{
  typedef std::map<TrkDetElemId, TrkExpectedHotSet*> map_t;
  typedef std::back_insert_iterator<std::vector<TrkExpectedHotSet*> > insert_t;
  std::transform(_exHotTable.begin(),_exHotTable.end(), insert_t(hotSets), babar::Pair::select2nd<map_t::value_type>());
}

void
TrkExpectedTrk::fillExHots( const TrkFitTypeKey& key,
                            std::vector<TrkExpectedHot*>& hots ) const
{
  typedef std::map<TrkDetElemId, TrkExpectedHotSet*> map_t;
  for (map_t::const_iterator itr=_exHotTable.begin();itr!=_exHotTable.end();++itr) {
    TrkExpectedHotSet* set = itr->second;
    TrkExpectedHot* aHot = const_cast<TrkExpectedHot*>(set->exHot(key));
    if ( aHot == 0 ) continue;
    hots.push_back(aHot);
  }
}


void
TrkExpectedTrk::printFit( const TrkFitTypeKey& key,
                          ostream& os ) const
{
  const TrkRecoTrk* refTrack = track();
  if ( refTrack ) refTrack->print(os);
  typedef std::map<TrkDetElemId, TrkExpectedHotSet*> map_t;
  for (map_t::const_iterator itr=_exHotTable.begin();itr!=_exHotTable.end();++itr) {
    TrkExpectedHotSet* exHots = itr->second;
    exHots->printFit(key,os);
  }
}

bool
TrkExpectedTrk::addExHotSet( const TrkDetElemId& id,
                             TrkExpectedHotSet* val )
{
  _exHotTable.insert(std::make_pair(id,val));
  val->setExTrk(this);
  return true;
}

//		-----------------------------------------
// 		-- Private Function Member Definitions --
//		-----------------------------------------

bool
TrkExpectedTrk::parseHotSetFromMap( TrkExpectedHotSet* hotSet,
				    const TrkFitTypeKey& key ){

  if ( hotSet == 0 || _parentMap == 0 ) return false;
  // get the other hots associated with same fund hit
  std::vector<TrkHitOnTrk*> hots;
  hotSet->getHots(hots);
  // loop on the hots
  for ( int i(0); i < hots.size(); ++i ) {
    const TrkHitOnTrk* hot = hots[i];
    // get the fit type and see if it is in the map
    TrkFitTypeKey fitTypeKey =  _parentMap->getFitType( hot );
    if ( ! ( fitTypeKey == key ) ) continue;
    // see if we already have exHot of this type
    //const TrkExpectedHot* exHot =
    _parentMap->exHotFactory()->parseExHotIntoSet( hot, fitTypeKey,*hotSet, false );
    // put in the the rep list, even if allready in,
    // because that is how me test the number of overlaps
    _reps.push_back( hot->getParentRep() );
  }
  return true;
}

bool
TrkExpectedTrk::parseHotSetFromMap( TrkExpectedHotSet* hotSet){

  if ( hotSet == 0 || _parentMap == 0 ) return false;
  // get the other hots associated with same fund hit
  std::vector<TrkHitOnTrk*> hots;
  hotSet->getHots(hots);
  // loop on the hots
  for ( int i(0); i < hots.size(); ++i ) {
    const TrkHitOnTrk* hot = hots[i];
    // get the fit type and see if it is in the map
    TrkFitTypeKey fitTypeKey =  _parentMap->getFitType( hot );
    if ( fitTypeKey.value() < 0 ) continue;
    // see if we already have exHot of this type
    //const TrkExpectedHot* exHot =
    _parentMap->exHotFactory()->parseExHotIntoSet( hot, fitTypeKey,*hotSet, false );
    // put in the the rep list, even if allready in,
    // because that is how me test the number of overlaps
    _reps.push_back( hot->getParentRep() );
  }
  return true;
}

bool TrkExpectedTrk::parseHotFromMap( const TrkHitOnTrk* hot,
                                      const bool fillGaps )
{
  // check parent map and fit type
  if ( _parentMap == 0 || hot == 0 ) return false;
  TrkFitTypeKey fitTypeKey = _parentMap->getFitType( hot );

  if ( fitTypeKey.value() < 0 ) return false;

  const TrkExpectedHotSet* hotSet = exHotSet(hot);
  TrkExpectedHotSet* theSet(0);
  if ( hotSet == 0 ) {
    if ( !fillGaps ) return false;
    const TrkDetElemId id = hot->hit()->elemId();
    TrkExpectedHotSet* theSet =
      new TrkExpectedHotSet( id, TrkFitTypeKey::currentKey() );
    addExHotSet( id, theSet );
  } else {
    theSet = (TrkExpectedHotSet*)hotSet;
  }

  //const TrkExpectedHot* exHot =
  _parentMap->exHotFactory()->parseExHotIntoSet( hot, fitTypeKey,*theSet, false);

  _reps.push_back( hot->getParentRep() );

  return true;

}

bool
TrkExpectedTrk::parseHotsFromMap( const TrkFitTypeKey& key )
{
  typedef std::map<TrkDetElemId, TrkExpectedHotSet*> map_t;
  for (map_t::const_iterator itr=_exHotTable.begin();itr!=_exHotTable.end();++itr) {
    if ( !parseHotSetFromMap( itr->second, key ) ) continue;
  }
  return true;
}

bool
TrkExpectedTrk::parseHotsFromMap( )
{
  typedef std::map<TrkDetElemId, TrkExpectedHotSet*> map_t;
  for (map_t::const_iterator itr=_exHotTable.begin();itr!=_exHotTable.end();++itr) {
    if ( !parseHotSetFromMap( itr->second ) ) continue;
  }
  return true;
}

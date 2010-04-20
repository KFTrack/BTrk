//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: TrkExpectedHotSet.cc,v 1.7 2004/09/10 18:00:17 bartoldu Exp $
//
// Description:
//	Class TrkExpectedHotSet, stores information about projected location 
//      of Tracks
//
// Environment:
//	Software developed for B Factory.
//
// Author List:
//	Eric A Charles        UW-Madison   
//
//------------------------------------------------------------------------

//----------------
// BaBar header
//----------------
#include "BaBar/BaBar.hh"

//-----------------------
// This Class's Header --
//-----------------------
#include "TrkBase/TrkExpectedHotSet.hh"

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

#include "TrkBase/TrkFitTypeKey.hh"
#include "TrkBase/TrkExpectedHot.hh"
using std::endl;
using std::ostream;

//----------------------------------------------------------------------
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

TrkExpectedHotSet::TrkExpectedHotSet( const TrkDetElemId& elemId,
				      const int& maxFits )
  :_elemId(elemId),
   _nFit(-1),
   _fitTypeMap(maxFits,-1),
   _exTrk(0)
{
}

//--------------
// Destructor --
//--------------

TrkExpectedHotSet::~TrkExpectedHotSet(){
  // owns the exHost
  const int iEx = _theHots.size();
  for ( int iKill(0); iKill < iEx; ++iKill ) {
    TrkExpectedHot* exHot = _theHots[iKill];
    delete exHot;
  }  
}

//-------------
// Operators --
//-------------
    
//-------------
// Selectors --
//-------------

const TrkExpectedHot* 
TrkExpectedHotSet::exHot( const TrkFitTypeKey& key ) const{
  if ( key.value() < 0 || key.value() >= _fitTypeMap.size() ) return 0;
  const int index = _fitTypeMap[key.value()];
  return index < 0 ? 0 : _theHots[index];
}

const TrkExpectedHot* 
TrkExpectedHotSet::exHot( ) const{
  return _nFit < 0 ? 0 : _theHots[0];
}

void 
TrkExpectedHotSet::getHots( std::vector<TrkHitOnTrk*>& hots ) const {
  const TrkExpectedHot* anExHot(0);
  for ( int i(0); i < _theHots.size(); ++i ) {
    anExHot = _theHots[i];
    if ( anExHot->getHots(hots) ) break;
  }
}

bool 
TrkExpectedHotSet::hasFitType( const TrkFitTypeKey& key ) const{
  return ( exHot(key) != 0 );  
}


void
TrkExpectedHotSet::printFit( const TrkFitTypeKey& key,
			     ostream& os ) const{
  os << " Intersects element "; 
  _elemId.printAll(os);
  const TrkExpectedHot* theExHot = exHot(key);
  if ( theExHot == 0 ) {
    os << " _" << endl;
    return;
  }
  theExHot->printAll(os);
  os << endl;
}

//-------------
// Modifiers --
//-------------

bool
TrkExpectedHotSet::addInfo ( TrkExpectedHot* aHot,
			     const TrkFitTypeKey& key ){
  if ( key.value() < 0 ) return false;
  if ( key.value() >= _fitTypeMap.size() ) return false;
  if ( aHot == 0 ) return false;
  const int test = _fitTypeMap[key.value()];
  if ( test != -1 ) {
    TrkExpectedHot* myHot = _theHots[test];
    if ( aHot->hasHot() ) { myHot->setHots(aHot); }
    delete aHot;
    return true;
  }
  ++_nFit;
  _fitTypeMap[key.value()] = _nFit;
  aHot->setExTrk(_exTrk);
  _theHots.push_back(aHot);
  return true;
}

void
TrkExpectedHotSet::setExTrk( const TrkExpectedTrk* anExTrk ) {
  _exTrk = anExTrk;
  const int nHot = _theHots.size();
  for ( int iHot(0); iHot < nHot; ++iHot ) {
    TrkExpectedHot* aHot = _theHots[iHot];
    aHot->setExTrk(_exTrk);
  }
}

bool 
TrkExpectedHotSet::isSvt() const {
  return _elemId.sysInd() == TrkDetElemId::svt;
}

bool 
TrkExpectedHotSet::isDch() const {
  return _elemId.sysInd() == TrkDetElemId::dch;
}

//		-------------------------------------------
// 		-- Protected Function Member Definitions --
//		-------------------------------------------

//		-----------------------------------------
// 		-- Private Function Member Definitions --
//		-----------------------------------------

//		-----------------------------------
// 		-- Internal Function Definitions --
//		-----------------------------------













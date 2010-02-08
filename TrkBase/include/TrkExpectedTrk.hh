//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: TrkExpectedTrk.hh,v 1.5 2004/08/06 06:31:41 bartoldu Exp $
//
// Description:
//	Class TrkExpectedTrack Stores lists of expected hits for a given 
//    set of tracks. 
//
// Environment:
//	Software developed for Babar experiment @ Slac B-factory.
//
// Author List:
//	Eric A Charles         UW-Madison 
//
//------------------------------------------------------------------------

#ifndef TRKEXPECTEDTRACK_HH
#define TRKEXPECTEDTRACK_HH

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
#include <iostream>

//----------------------
// Base Class Headers --
//----------------------

//-------------------------------
// Collaborating Class Headers --
//-------------------------------

#include "PDT/PdtPid.hh"
#include "TrkBase/TrkFitTypeKey.hh"

//------------------------------------
// Collaborating Class Declarations --
//------------------------------------

class TrkRecoTrk;
class TrkRep;
class TrkFit;
class TrkExpectedHotSet;
class TrkExpectedHot;
class TrkHitOnTrk;
class TrkFundHit;
class TrkExpectedMap;
class TrkDetElemId;
class GTrack;

//		---------------------
// 		-- Class Interface --
//		---------------------

class TrkExpectedTrk {

//--------------------
// Instance Members --
//--------------------

public:

  // Constructor
  
  // c'tor from a track, uses the default fit unless overridden,
  // this also registers the expected track in the map,
  // finally, it sets the tracks fit type as the default type,
  // it does not, however, actully make the TrkExpectedHotSets
  TrkExpectedTrk( const TrkRecoTrk* aTrack,
		  const PdtPid::PidType = PdtPid::null,
		  const GTrack* gTrk = 0,
		  TrkExpectedMap* exMap = 0 );

  // Destructor
  virtual ~TrkExpectedTrk( );
  
  // Operators
  bool operator==( const TrkExpectedTrk& rhs ) const {
    return this == &rhs; 
  }
   
  // Selectors (const)
  const TrkHitOnTrk* getHot( const TrkFundHit* hit,
			     const TrkFitTypeKey& key ) const;

  const TrkRep* refRep() const{
    return _refRep;
  }

  const GTrack* gTrack() const {
    return _gTrk;
  }

  const TrkRecoTrk* track() const;

  // these are ways of getting the hots sets
  // by element id 
  const TrkExpectedHotSet* exHotSet( const TrkDetElemId& elemId ) const;
  // by hot and fund hit
  const TrkExpectedHotSet* exHotSet( const TrkHitOnTrk* hot ) const;
  const TrkExpectedHotSet* exHotSet( const TrkFundHit* hit ) const;

  // these are ways of getting the ex hots them-selves
  // by element id
  const TrkExpectedHot* exHot( const TrkDetElemId& elemId,
			       const TrkFitTypeKey key = TrkFitTypeKey(0) ) const;
  
  // by hot and fund hit
  const TrkExpectedHot* exHot( const TrkHitOnTrk* hot ) const;
  const TrkExpectedHot* exHot( const TrkFundHit* hit,
			       const TrkFitTypeKey key = TrkFitTypeKey(0)) const;

  void fillExHotSets( std::vector<TrkExpectedHotSet*>& hotSets ) const;
  void fillExHots( const TrkFitTypeKey& key, 
                  std::vector<TrkExpectedHot*>& hots ) const;

  // get the number of overlaps with a rep
  int overLap( const TrkRep* rep ) const
  {
          // FIXME: std::count doesn't work with the BaBar config of Sun WS6U1..
    // return std::count(_reps.begin(),_reps.end(),rep);
    typedef std::vector<const TrkRep*>::const_iterator i_t;
    int j=0;
    for (i_t i=_reps.begin(); i!=_reps.end();++i) {
            if (*i == rep) ++j;
    }
    return j;
  }

  //Accessing methods for the class private members:
  const std::map<TrkDetElemId, TrkExpectedHotSet*>& expectedHots() const 
  {
    return _exHotTable;
  }

  // modifiers
  void printFit( const TrkFitTypeKey& key, 
                 std::ostream& os = std::cout ) const;

protected:

  // Helper functions
  // add in the un-parsed hots from a rep
  bool parseHotFromMap( const TrkHitOnTrk* hot, const bool fillGaps );
  bool parseHotSetFromMap( TrkExpectedHotSet* hotSet,
                           const TrkFitTypeKey& key );
  bool parseHotSetFromMap( TrkExpectedHotSet* hotSet );
  bool parseHotsFromMap( const TrkFitTypeKey& key );
  bool parseHotsFromMap( );

  const TrkRep* getRep( const TrkRecoTrk* aTrk,
                        const PdtPid::PidType type ) const;

  bool addExHotSet( const TrkDetElemId& elemId,
                    TrkExpectedHotSet* val );

private:
  // friends
  friend class TrkExpectedHotFactory;

  // Data members
  std::vector<const TrkRep*> _reps;
  const TrkExpectedMap* _parentMap;
  std::map<TrkDetElemId,TrkExpectedHotSet*> _exHotTable;
  const TrkRep* _refRep;
  const GTrack* _gTrk;
};

#endif

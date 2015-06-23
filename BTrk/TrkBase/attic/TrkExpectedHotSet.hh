//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkExpectedHotSet.hh,v 1.4 2004/08/06 06:31:40 bartoldu Exp $
//
// Description:
//      Defines a set of Expected hit on a track.
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Authors: Eric Charles            UW Madison
//           
//------------------------------------------------------------------------

#ifndef TRKEXPECTEDHOTSET_HH
#define TRKEXPECTEDHOTSET_HH

//---------------
// C++ Headers --
//---------------

#include <iostream>
#include "TrkBase/TrkDetElemId.hh"

//----------------
//  BaBar Headers
//----------------

#include <vector>

class TrkExpectedTrk;
class TrkExpectedHot;
class TrkFitTypeKey;
class TrkFundHit;
class TrkHitOnTrk;
class TrkHitList;

// Class interface //
class TrkExpectedHotSet {

public:

  // Construct from the cell and the hit
  TrkExpectedHotSet( const TrkDetElemId& elemId,
		     const int& maxFits );

  // Destructor
  virtual ~TrkExpectedHotSet();

  // operator == for vectors, check identity
  bool operator== ( const TrkExpectedHotSet& rhs ) {
    return this == &rhs;
  }

  // Access to the cell id (offset by subsystem)
  const TrkDetElemId& elemId() const{ 
    return _elemId;
  };

  // access to the hot ( from the key )
  const TrkExpectedHot* exHot( const TrkFitTypeKey& key ) const;
  // the "reference" exHot
  const TrkExpectedHot* exHot( ) const;  

  // access to the other hots from fund hits
  void getHots( std::vector<TrkHitOnTrk*>& hots ) const;

  // is this fit type already in?
  bool hasFitType( const TrkFitTypeKey& key ) const;

  // other access
  const std::vector<int>& fitTypeMap() const { return _fitTypeMap; }  
  std::vector<TrkExpectedHot*> theHots() const { return _theHots; }

  // add info from a track
  bool addInfo( TrkExpectedHot* aHot,
	        const TrkFitTypeKey& key );
  
  void printFit( const TrkFitTypeKey& key,
		 std::ostream& os = std::cout ) const;
  
  bool isSvt() const;
  bool isDch() const;

  const TrkExpectedTrk* exTrk() const { return _exTrk; }

protected:

  void setExTrk( const TrkExpectedTrk* anExTrk );

private:

  const TrkDetElemId _elemId;
  int _nFit;
  
  std::vector<int> _fitTypeMap;
  std::vector<TrkExpectedHot*> _theHots;

  friend class TrkExpectedTrk;
  friend class TrkExpectedHotFactory;

  const TrkExpectedTrk* _exTrk;

};

#endif






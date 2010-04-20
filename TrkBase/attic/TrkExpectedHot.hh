//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkExpectedHot.hh,v 1.4 2004/08/06 06:31:40 bartoldu Exp $
//
// Description:
//      Defines a Expected hit on a track, and stores information about it.
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Authors: Eric Charles            UW Madison
//           
//
//------------------------------------------------------------------------

#ifndef TRKEXPECTEDHOT_HH
#define TRKEXPECTEDHOT_HH

//---------------
// C++ Headers --
//---------------

#include <iostream>
#include "TrkBase/TrkDetElemId.hh"

//----------------
//  BaBar Headers
//----------------

class TrkExpectedTrk;
class TrkHitOnTrk;
class TrkDifPoca;

#include <vector>

// Class interface //
class TrkExpectedHot {

public:
		 
  // copy c'tor
  TrkExpectedHot( const TrkExpectedHot& );

  // Destructor
  virtual ~TrkExpectedHot();

  // access to the internal id of the associated element 
  // ( offset by subsystem )
  const TrkDetElemId& id() const{ return _id; }

  // sorting
  // sorting uses the internal id
  bool operator == ( const TrkExpectedHot& rhs ) const;
  bool operator <  ( const TrkExpectedHot& rhs ) const;

  virtual bool setHot( const TrkHitOnTrk* aHot ) = 0;
  virtual bool setHots( const TrkExpectedHot* aHot ) = 0;

  virtual bool hasHot( ) const = 0;
  virtual bool hasGHit( ) const = 0;

  bool isSvt() const;
  bool isDch() const;

  virtual bool getHots( std::vector<TrkHitOnTrk*>& hots ) const=0;

  virtual void printAll( std::ostream& os ) const = 0;

  const TrkExpectedTrk* exTrk() const { return _exTrk; };
  
protected:

  // c'tors are protected, only called by concrete classes
  TrkExpectedHot( const TrkDetElemId& id );
  TrkExpectedHot( const int& id, TrkDetElemId::systemIndex sysInd );

  void setExTrk( const TrkExpectedTrk* anExTrk ) { _exTrk = anExTrk; }

private:

  // elem id
  const TrkDetElemId _id;
  const TrkExpectedTrk* _exTrk;

  friend class TrkExpectedHotSet;

};

#endif


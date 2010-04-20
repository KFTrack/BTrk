//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: TrkFitTypeKey.hh,v 1.3 2004/08/06 06:31:41 bartoldu Exp $
//
// Description:
//	Class TrkFitTypeKey
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

#ifndef TRKFITTYPEKEY_HH
#define TRKFITTYPEKEY_HH

//-------------
// C Headers --
//-------------
#include <vector>
#include <iostream>
#include "PDT/PdtPid.hh"

class IfdKey;

class TrkFitTypeKey {

  //--------------------
  // static Members --
  //--------------------

public:
  // methods
  static int currentKey() { return _keys.size(); }

  static void killFitType(const TrkFitTypeKey key);

protected:
  static int newFitType(const char* name=0);

private:

  // Data members
  static std::vector<IfdKey*> _keys;
  static std::vector<PdtPid::PidType> _pidTypes;
  static const int _defaultValue;

  //--------------------
  // Instance Members --
  //--------------------

public:

  // standard c'tor makes a new IfdKey from a char* and associates a PidType
  TrkFitTypeKey( const char* name = 0,
                 const PdtPid::PidType pid = PdtPid::null );

  // copy c'tor
  TrkFitTypeKey( const TrkFitTypeKey& );
  TrkFitTypeKey( const int& );

  // Destructor
  virtual ~TrkFitTypeKey( );

  // Operators
  bool operator==( const TrkFitTypeKey& rhs ) const{
    return _value == rhs.value();
  }
  bool operator< ( const TrkFitTypeKey& rhs ) const{
    return _value < rhs.value();
  }
  TrkFitTypeKey& operator= ( const TrkFitTypeKey& rhs ) {
    _value = rhs.value();
    return *this;
  }

  // Selectors (const)
  const int& value() const { return _value; }

  const IfdKey* ifdKey() const {
    if ( _value < 0 || _value >= currentKey() ) return 0;
    return _keys[_value];
  }

  PdtPid::PidType pidType() const {
    if ( _value < 0 || _value >= currentKey() ) return PdtPid::null;
    return _pidTypes[_value];
  }

  void printAll( std::ostream& os = std::cout ) const;

private:

  // Data members
  int _value;

};

std::ostream& operator<<(std::ostream& os, const TrkFitTypeKey& key);

#endif

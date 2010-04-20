//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: TrkDetElemId.hh,v 1.2 2004/08/06 06:31:39 bartoldu Exp $
//
// Description:
//	Class TrkDetElemId
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

#ifndef TRKDETELEMID_HH
#define TRKDETELEMID_HH

//-------------
// C Headers --
//-------------
extern "C" {
}

//---------------
// C++ Headers --
//---------------

#include <iostream>

//		---------------------
// 		-- Class Interface --
//		---------------------

class TrkDetElemId {

//------------------
// Static Members --
//------------------

public:

  // Typedefs, consts, and enums
  enum systemIndex { null=0};

  // functions
  static int calcValue( const int& id, 
			TrkDetElemId::systemIndex sysInd );

private:

  static const int nullElemID;

//--------------------
// Instance Members --
//--------------------

public:

  // Constructors
  TrkDetElemId( const int& id, TrkDetElemId::systemIndex sysInd );

  // Copy Constructor
  TrkDetElemId( const TrkDetElemId& );
  
  // Destructor
  virtual ~TrkDetElemId( );

  // Operators
    
  TrkDetElemId& operator= ( const TrkDetElemId& );

  bool operator==( const TrkDetElemId& rhs ) const{
    return elemId() == rhs.elemId(); 
  }
  bool operator<( const TrkDetElemId& rhs) const {
    return elemId() < rhs.elemId(); 
  }

  // Selectors (const)
  int elemId() const{ 
    return calcValue(_id,_sysInd);
  };

  const int& systemElemId() const{
    return _id;
  }

  const TrkDetElemId::systemIndex& sysInd() const {
    return _sysInd;
  }

  void printAll( std::ostream& os=std::cout ) const;

private:
  
  // Data members
  int _id;
  TrkDetElemId::systemIndex _sysInd;

};
 
std::ostream& operator<<(std::ostream& os, const TrkDetElemId& id);

// Inline implementations
//#include "TrkDetElemId.icc"

#endif

//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: TrkExpectedHot.cc,v 1.3 2004/09/10 18:00:17 bartoldu Exp $
//
// Description:
//	Class TrkExpectedHot
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
#include "TrkBase/TrkExpectedHot.hh"

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

TrkExpectedHot::TrkExpectedHot( const TrkDetElemId& id )
  :_id(id),
   _exTrk(0){
}

TrkExpectedHot::TrkExpectedHot(  const int& id, 
				 TrkDetElemId::systemIndex sysInd )
  :_id(id,sysInd),
   _exTrk(0){
}


// copy c'tor
TrkExpectedHot::TrkExpectedHot( const TrkExpectedHot& rhs )
  :_id(rhs.id()),
   _exTrk(0){
}

//--------------
// Destructor --
//--------------

TrkExpectedHot::~TrkExpectedHot(){}

//-------------
// Operators --
//-------------
    
bool
TrkExpectedHot::operator==(const TrkExpectedHot& rhs) const{
  return _id == rhs.id();
}

bool
TrkExpectedHot::operator< (const TrkExpectedHot& rhs) const{
  return _id < rhs.id();
}

bool 
TrkExpectedHot::isSvt() const {
  return _id.sysInd() == TrkDetElemId::svt;
}

bool 
TrkExpectedHot::isDch() const {
  return _id.sysInd() == TrkDetElemId::dch;
}



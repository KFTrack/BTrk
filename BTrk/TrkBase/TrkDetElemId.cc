//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: TrkDetElemId.cc,v 1.3 2004/09/10 18:00:17 bartoldu Exp $
//
// Description:
//	Class TrkDetElemId
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
#include "BTrk/BaBar/BaBar.hh"

//-----------------------
// This Class's Header --
//-----------------------
#include "BTrk/TrkBase/TrkDetElemId.hh"
using std::endl;
using std::ostream;

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

const int TrkDetElemId::nullElemID(-1);

int TrkDetElemId::calcValue( const int& id, 
			     TrkDetElemId::systemIndex sysInd ){
  int retVal(id);
  switch (sysInd) {
  default:
    break;
  }
  return TrkDetElemId::nullElemID;
}

//		----------------------------------------
// 		-- Public Function Member Definitions --
//		----------------------------------------

//----------------
// Constructors --
//----------------

TrkDetElemId::TrkDetElemId( const int& id, TrkDetElemId::systemIndex sysInd )
  :_id(id),_sysInd(sysInd){
}


TrkDetElemId::TrkDetElemId( const TrkDetElemId& rhs )
  :_id(rhs.systemElemId()),_sysInd(rhs.sysInd()){
}

//--------------
// Destructor --
//--------------

TrkDetElemId::~TrkDetElemId(){
}

//-------------
// Operators --
//-------------

TrkDetElemId& 
TrkDetElemId::operator= ( const TrkDetElemId& rhs ){
  _id = rhs.systemElemId();
  _sysInd = rhs.sysInd();
  return *this;
}
   
void 
TrkDetElemId::printAll( ostream& os ) const{
  os << _id << ' ' << _sysInd << ' ' << elemId() << endl;  
}

ostream& operator<<(ostream& os, const TrkDetElemId& id){
  id.printAll(os);
  return os;
}


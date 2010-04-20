//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: TrkFitTypeKey.cc,v 1.5 2004/09/10 18:00:18 bartoldu Exp $
//
// Description:
//	Class TrkFitTypeKey
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
#include "TrkBase/TrkFitTypeKey.hh"

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

#include "ProxyDict/IfdStrKey.hh"
using std::endl;
using std::ostream;

//-----------------------------------------------------------------------
// Local Macros, Typedefs, Structures, Unions and Forward Declarations --
//-----------------------------------------------------------------------





//		-----------------------------------------------
// 		-- Static Data & Function Member Definitions --
//		-----------------------------------------------

const int TrkFitTypeKey::_defaultValue(-1);
std::vector<IfdKey*> TrkFitTypeKey::_keys ;
std::vector<PdtPid::PidType> TrkFitTypeKey::_pidTypes ;

int
TrkFitTypeKey::newFitType( const char* name )
{
  if ( name == 0 ) return _defaultValue;
  _keys.push_back( new IfdStrKey(name) );
  return (currentKey()-1);
}

void
TrkFitTypeKey::killFitType( const TrkFitTypeKey key )
{
  const int i = key.value();
  if ( i < 0 || i >= currentKey() ) return;
  IfdKey* theKey = _keys[i];
  assert(theKey != 0 );
  delete theKey;
  theKey = 0;
}

//----------------
// Constructors --
//----------------

TrkFitTypeKey::TrkFitTypeKey( const char* name,
                              const PdtPid::PidType pid )
  :_value(newFitType(name))
{
    if ( _value >= 0 ) _pidTypes.push_back(pid);
}

TrkFitTypeKey::TrkFitTypeKey( const TrkFitTypeKey& rhs )
  :_value(rhs.value())
{
}

TrkFitTypeKey::TrkFitTypeKey( const int& val )
  :_value(val)
{
}

//--------------
// Destructor --
//--------------

TrkFitTypeKey::~TrkFitTypeKey( )
{
}


void
TrkFitTypeKey::printAll( ostream& os ) const
{
  char* pidName(0);
  switch ( pidType() ) {
  case PdtPid::electron:
    pidName = "Electron";
    break;
  case PdtPid::muon:
    pidName = "Muon";
    break;
  case PdtPid::pion:
    pidName = "Pion";
    break;
  case PdtPid::kaon:
    pidName = "Kaon";
    break;
  case PdtPid::proton:
    pidName = "Proton";
    break;
  case PdtPid::null:
  default:
    pidName = "Unknown";
    break;
  }
  os << ifdKey() << ' ' << pidName << " key: " << _value << endl;
}

ostream&
operator<<(ostream& os, const TrkFitTypeKey& key)
{
  key.printAll(os);
  return os;
}

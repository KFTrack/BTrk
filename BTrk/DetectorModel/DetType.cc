//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DetType.cc,v 1.16 2004/08/06 05:58:33 bartoldu Exp $
//
// Description:
//	Class DetType
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//
// Copyright Information:
//	Copyright (C) 1997	
//
//------------------------------------------------------------------------
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/BaBar/ErrLog.hh"

//-----------------------
// This Class's Header --
//-----------------------
#include "BTrk/DetectorModel/DetType.hh"

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
#include <vector>
using std::ostream;

//-----------------------------------------------------------------------
// Local Macros, Typedefs, Structures, Unions and Forward Declarations --
//-----------------------------------------------------------------------

//		----------------------------------------
// 		-- Public Function Member Definitions --
//		----------------------------------------

//----------------
// Constructors --
//----------------
DetType::DetType()
  : _tname("Unknown"),
    _itype(-1),
    _outline( new std::vector< TypeCoord* > ), 
    _pointStore( new std::vector< TypeCoord* > )
{
}

DetType::DetType(const char* name,int ityp)
  : _tname(name), 
    _itype(ityp), 
    _outline( new std::vector< TypeCoord* > ), 
    _pointStore( new std::vector< TypeCoord* > )
{
}

//--------------
// Destructor --
//--------------
DetType::~DetType()
{
  delete _outline;
  std::vector<TypeCoord*>::iterator piter = _pointStore->begin();
  TypeCoord* aTC = 0;
  while( piter != _pointStore->end() ) {
    aTC = *piter++;
    if( aTC != 0 ) {
      delete aTC;
      aTC = 0;
    }
  }
  _pointStore->clear();
  delete _pointStore;
}

//-------------
// Methods   --
//-------------
    
//-------------
// Operators --
//-------------
bool
DetType::operator==( const DetType& otherType ) const
{
  bool answer = false;

  if ( _itype == otherType._itype &&
       _tname == otherType._tname ) answer = true;

  return answer;
}
  
//-------------
// Selectors --
//-------------
void
DetType::print( ostream& os ) const 
{
  os << "Detector type # " << _itype << " " << _tname;
}
    
const std::vector< TypeCoord* >*
DetType::pointStore() const
{
  return _pointStore;
}

const std::vector< TypeCoord* >*
DetType::outline() const
{
  return _outline;
}

//-------------
// Modifiers --
//-------------

//		-----------------------------------------------
// 		-- Static Data & Function Member Definitions --
//		-----------------------------------------------

//		-------------------------------------------
// 		-- Protected Function Member Definitions --
//		-------------------------------------------
std::vector< TypeCoord* >*
DetType::myOutline()
{
  return _outline;
}

std::vector< TypeCoord* >*
DetType::myPointStore()
{
  return _pointStore;
}


//		-----------------------------------------
// 		-- Private Function Member Definitions --
//		-----------------------------------------

//		-----------------------------------
// 		-- Internal Function Definitions --
//		-----------------------------------

//

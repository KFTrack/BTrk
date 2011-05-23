//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DetVolSideIntersection.cc,v 1.3 2000/11/07 18:19:09 elmer Exp $
//
// Description:
//	Class DetVolSideIntersection
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	Stephen J. Gowdy	        Originator
//	Gautier Hamel de Monchenault	Originator
//
// Copyright Information:
//	Copyright (C) 1997	University of Edinburgh
//	Copyright (C) 1997	CEA - Centre d'Etude de Saclay
//
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"
#include "ErrLogger/ErrLog.hh"

//-----------------------
// This Class's Header --
//-----------------------
#include "DetectorModel/DetVolSideIntersection.hh"

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

//		----------------------------------------
// 		-- Public Function Member Definitions --
//		----------------------------------------

//----------------
// Constructors --
//----------------
DetVolSideIntersection::DetVolSideIntersection()
  : _pathlen(0.0),
    _side(0)
{
}

DetVolSideIntersection::DetVolSideIntersection
                                  ( const DetVolSideIntersection& other ) 
  : _pathlen(other._pathlen),
    _side(other._side),
    _coord(other._coord)
{
}

DetVolSideIntersection::DetVolSideIntersection( double len, int i,
						const TwoDCoord& c )
  : _pathlen(len),
    _side(i),
    _coord(c)
{
}

//--------------
// Destructor --
//--------------
DetVolSideIntersection::~DetVolSideIntersection()
{
}

//-------------
// Methods   --
//-------------
    
//-------------
// Operators --
//-------------
bool
DetVolSideIntersection::operator==( const DetVolSideIntersection& test ) const
{
  return _side == test._side; 
}
    
bool 
DetVolSideIntersection::operator<( const DetVolSideIntersection& test )  const
{
  return _pathlen < test._pathlen; 
}

bool
DetVolSideIntersection::operator>( const DetVolSideIntersection& test )  const
{
  return _pathlen > test._pathlen; 
}  

DetVolSideIntersection&
DetVolSideIntersection::operator=( const DetVolSideIntersection& copy )
{
  _side    = copy._side;
  _pathlen = copy._pathlen; 
  _coord   = copy._coord;

  return *this;
}

//-------------
// Selectors --
//-------------
    
//-------------
// Modifiers --
//-------------

//		-----------------------------------------------
// 		-- Static Data & Function Member Definitions --
//		-----------------------------------------------

//		-------------------------------------------
// 		-- Protected Function Member Definitions --
//		-------------------------------------------

//		-----------------------------------------
// 		-- Private Function Member Definitions --
//		-----------------------------------------

//		-----------------------------------
// 		-- Internal Function Definitions --
//		-----------------------------------

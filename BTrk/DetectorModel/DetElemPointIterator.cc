//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DetElemPointIterator.cc,v 1.7 2002/12/30 15:44:28 dbrown Exp $
//
// Description:
//	Class DetElemPointIterator
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	Gautier Hamel de Monchenault
//      Stephen J. Gowdy
//
// Copyright Information:
//	Copyright (C) 1997	Lawrence Berkeley Lab
//
//------------------------------------------------------------------------
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/BaBar/ErrLog.hh"

//-----------------------
// This Class's Header --
//-----------------------
#include "BTrk/DetectorModel/DetElemPointIterator.hh"

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "BTrk/BbrGeom/HepPoint.h"
#include "BTrk/DetectorModel/DetElem.hh"


//-----------------------------------------------------------------------
// Local Macros, Typedefs, Structures, Unions and Forward Declarations --
//-----------------------------------------------------------------------

//		----------------------------------------
// 		-- Public Function Member Definitions --
//		----------------------------------------

//----------------
// Constructors --
//----------------
DetElemPointIterator::DetElemPointIterator( const DetElem& elem )
  :  _elem( &elem ), _vect( _elem->detectorType()->outline()), _current( 0 )  
{
}

//--------------
// Destructor --
//--------------
DetElemPointIterator::~DetElemPointIterator() 
{
}

//-------------
// Methods   --
//-------------
Action
DetElemPointIterator::next( HepPoint& point ) 
{
  if( _current == _vect->size() )  return CloseShape;
  point = HepPoint(0.,0.,0.);
  TypeCoord* ptr = (*_vect)[ _current++ ];
  if( ptr == 0 ) {
    if( next( point ) == Continue ) 
      return CloseLine;
    else
      return CloseShape;
  }
  point = _elem->coordToPoint( ptr );
  return Continue;
}
    
//-------------
// Operators --
//-------------
    
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





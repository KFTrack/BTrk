//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DetVolumeType.cc,v 1.10 2002/12/30 15:44:30 dbrown Exp $
//
// Description:
//	Class DetVolumeType. See header for information.
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
#include "DetectorModel/DetVolumeType.hh"

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
#include "DetectorModel/DetSurface.hh"
#include <vector>

//-----------------------------------------------------------------------
// Local Macros, Typedefs, Structures, Unions and Forward Declarations --
//-----------------------------------------------------------------------

//		----------------------------------------
// 		-- Public Function Member Definitions --
//		----------------------------------------

//----------------
// Constructors --
//----------------
DetVolumeType::DetVolumeType( const char* name, int id )
  : DetType( name, id ),
    _sides( new std::vector< DetSurface* > ),
    _sideCorners( new std::vector< std::vector< SurfacePoint* > >( 0 ) )
{
  
}

//--------------
// Destructor --
//--------------
DetVolumeType::~DetVolumeType()
{
  // Delete the side corners and their lists
  size_t i=0, nLists = _sideCorners->size();
  for ( i=0; i<nLists; i++ )
    {
      size_t j=0, nPoints = (*_sideCorners)[i].size();
      for ( j=0; j<nPoints; j++ ) delete (*_sideCorners)[i][j];
    }
  delete _sideCorners;
  // Delete the surfaces
  size_t nSides = _sides->size();
  for ( i=0; i<nSides; i++ ) delete (*_sides)[i];
  delete _sides;
}

//-------------
// Methods   --
//-------------
    
//-------------
// Operators --
//-------------
    
//-------------
// Selectors --
//-------------
bool
DetVolumeType::physicalMaterial(const TypeCoord*) const
{
  return true;
}

bool
DetVolumeType::insideLimitsOf( int side,
			       const SurfacePoint& thisPoint ) const
{
  size_t i=0, numberOfCorners = (*_sideCorners)[side].size();

  // Loop over all the edges of this side and see if the point is inside
  //   (inside means the same side as the origin)
  for ( i=0; i < ( numberOfCorners - 1 ); i++ )
    if ( ! insideLine( thisPoint, *(*_sideCorners)[side][i],
		       *(*_sideCorners)[side][i+1] ) ) return false;
  // To check the last edge has a special case
  if ( ! insideLine( thisPoint, *(*_sideCorners)[side][numberOfCorners-1],
		     *(*_sideCorners)[side][0] ) ) return false;

  return true;
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

//		-----------------------------------------
// 		-- Private Function Member Definitions --
//		-----------------------------------------

//		-----------------------------------
// 		-- Internal Function Definitions --
//		-----------------------------------

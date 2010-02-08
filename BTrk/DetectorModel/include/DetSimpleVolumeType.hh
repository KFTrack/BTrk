//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DetSimpleVolumeType.hh,v 1.6 2002/12/30 15:44:29 dbrown Exp $
//
// Description:
//	Class DetSimpleVolumeType. This class describe the type of a volume
//                                 made from planes.
//
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

#ifndef DETSIMPLEVOLUMETYPE_HH
#define DETSIMPLEVOLUMETYPE_HH

//-------------
// C Headers --
//-------------
extern "C" {
}

//---------------
// C++ Headers --
//---------------

//----------------------
// Base Class Headers --
//----------------------
#include "DetectorModel/DetVolumeType.hh"

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include <vector>

//------------------------------------
// Collaborating Class Declarations --
//------------------------------------
class DetMaterial;
class HepPoint;

//		---------------------
// 		-- Class Interface --
//		---------------------

class DetSimpleVolumeType : public DetVolumeType
{

//--------------------
// Declarations     --
//--------------------

  // Typedefs, consts, and enums
  enum { near, far };

//--------------------
// Instance Members --
//--------------------

public:

  // Constructors
  DetSimpleVolumeType( const char* typeName, int typeId,
		       const DetMaterial* theMaterial,
		       std::vector< HepPoint >& theCorners );

  // Copy Constructor

  // Destructor
  virtual ~DetSimpleVolumeType();

  // Operators

  // Selectors (const)
  virtual const DetMaterial& material( const TypeCoord* here ) const
    { return *_material; }

  // Modifiers

protected:

  // Helper functions
  virtual bool insideLine( const SurfacePoint& thisPoint,
			   const SurfacePoint& p1,
			   const SurfacePoint& p2 ) const;
  

private:

  // Friends

  // Data members

//------------------
// Static Members --
//------------------

public:

  // Selectors (const)

  // Modifiers

private:

  // Data members
  const DetMaterial* _material;

};

//#include "DetectorModel/DetSimpleVolumeType.icc"

#endif // DETSIMPLEVOLUMETYPE_HH

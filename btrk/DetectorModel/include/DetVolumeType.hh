//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DetVolumeType.hh,v 1.10 2002/12/30 15:44:30 dbrown Exp $
//
// Description:
//	Class DetVolumeType. This class describe the type of a volume
//                           made from surfaces.
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

#ifndef DETVOLUMETYPE_HH
#define DETVOLUMETYPE_HH

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
#include "DetectorModel/DetType.hh"

//-------------------------------
// Collaborating Class Headers --
//-------------------------------

#include <vector>

//------------------------------------
// Collaborating Class Declarations --
//------------------------------------
class DetSurface;
class SurfacePoint;

//		---------------------
// 		-- Class Interface --
//		---------------------

class DetVolumeType : public DetType
{
//--------------------
// Declarations     --
//--------------------

  // Typedefs, consts, and enums

//--------------------
// Instance Members --
//--------------------

public:

  // Constructors
  DetVolumeType( const char* name, int id );

  // Copy Constructor

  // Destructor
  virtual ~DetVolumeType();

  // Operators
  
  // Selectors (const)
  virtual bool physicalMaterial(const TypeCoord*) const;
  virtual bool insideLimitsOf( int aSide,const SurfacePoint& thisPoint ) 
    const;
  virtual const std::vector< DetSurface* >* sides() const
    { return _sides; }
  virtual const std::vector< std::vector< SurfacePoint* > >* sideCorners() 
    const
    { return _sideCorners; }

  // Modifiers

protected:

  // Helper functions
  virtual bool insideLine( const SurfacePoint& thisPoint,
			   const SurfacePoint& p1,
			   const SurfacePoint& p2 ) const = 0;
  virtual std::vector< DetSurface* >* mySides() 
    { return _sides; }
  virtual std::vector< std::vector< SurfacePoint* > >* mySideCorners()
    { return _sideCorners; }
private:

  // Friends
  friend class DetVolumeElem;

  // Data members
  std::vector< DetSurface* >* _sides;
  std::vector< std::vector< SurfacePoint* > >* _sideCorners;

//------------------
// Static Members --
//------------------

public:

  // Selectors (const)

  // Modifiers

private:

// prohibit
  DetVolumeType& operator = (const DetVolumeType&);
  DetVolumeType(const DetVolumeType&);

};


#endif // DETVOLUMETYPE_HH

//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DetVolSideIntersection.hh,v 1.1 1997/05/09 18:08:21 gowdy Exp $
//
// Description:
//	Class DetVolSideIntersection. This class describe the intersection
//                                    of a trajectory on a side of a 
//                                    DetVolume.
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

#ifndef DETVOLSIDEINTERSECTION_HH
#define DETVOLSIDEINTERSECTION_HH

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

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "BTrk/DetectorModel/DetType.hh"

//------------------------------------
// Collaborating Class Declarations --
//------------------------------------

//		---------------------
// 		-- Class Interface --
//		---------------------

class DetVolSideIntersection
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
  DetVolSideIntersection();
  DetVolSideIntersection( const DetVolSideIntersection& otherIntersection );
  DetVolSideIntersection( double pathLength, int nSide,
			  const TwoDCoord& theCoords );

  // Copy Constructor

  // Destructor
  virtual ~DetVolSideIntersection();

  // Operators
  bool operator==( const DetVolSideIntersection& otherInter ) const;
  bool operator<( const DetVolSideIntersection& otherInter ) const;  
  bool operator>( const DetVolSideIntersection& otherInter ) const;
  DetVolSideIntersection& operator=( const DetVolSideIntersection& toAssign );

  // Selectors (const)
  virtual double pathLength() const { return _pathlen; }
  virtual int side() const { return _side; }

  // Modifiers

protected:

  // Helper functions

private:

  // Friends

  // Data members
  double    _pathlen; // path length at intersection
  int       _side;    // side number
  TwoDCoord _coord;   // coordinates at intersection

//------------------
// Static Members --
//------------------

public:

  // Selectors (const)

  // Modifiers

private:

  // Data members

};

//#include "BTrk/DetectorModel/DetVolSideIntersection.icc"

#endif // DETVOLSIDEINTERSECTION_HH

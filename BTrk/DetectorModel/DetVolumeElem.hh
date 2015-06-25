//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DetVolumeElem.hh,v 1.13 2004/08/06 05:58:33 bartoldu Exp $
//
// Description:
//	DetVolumeElem Class - This class represents a volume made from
//                        surfaces.
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	Gautier Hamel de Monchenault - CEN Saclay & Lawrence Berkeley Lab
//      Stephen J. Gowdy              University of Edinburgh
//
// History (add to end):
//      Gautier   Sep 16, 1996  - creation
//      Stephen   Feb 20, 1997  - adopt for EmcXtal (was DrcBar)
//      Both      May 08, 1997  - adopt for DetVolumeElem
//
// Copyright Information:
//	Copyright (C) 1997		University of Edinburghatory
//	Copyright (C) 1997		Lawrence Berkeley Laboratory
//	CopPMT   yright (C) 1997	       CEA - Centre d'Etude de Saclay
//
//------------------------------------------------------------------------
#ifndef DETVOLUMEELEM_HH
#define DETVOLUMEELEM_HH

//-------------
// C Headers --
//-------------
extern "C" {
}

//---------------
// C++ Headers --
//---------------
#include <vector>
#include <algorithm>

//----------------------
// Base Class Headers --
//----------------------
#include "BTrk/DetectorModel/DetElem.hh"

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "BTrk/DetectorModel/DetVolSideIntersection.hh"
//------------------------------------
// Collaborating Class Declarations --
//------------------------------------
class DetAlignElem;
class DetIntersection;
class DetSurface;
class DetVolumeType;
class HepTransformation;
class Trajectory;
#include <iosfwd>

//		---------------------
// 		-- Class Interface --
//		---------------------

class DetVolumeElem : public DetElem
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
  DetVolumeElem( DetVolumeType* itsType, const char* name, int id,
		 const HepTransformation& theAlignment );

  // Destructor
  virtual ~DetVolumeElem();

  // Operators

  // Selectors (const)
  virtual int intersect( const Trajectory*, DetIntersection& ) const;
  virtual void sideIntersect( const Trajectory*, 
			      std::vector< DetVolSideIntersection >&,
			      double flightDistance,
			      double* arrayOfRange ) const;
  virtual void physicalOutline(std::vector<HepPoint>&) const;
  virtual void print(std::ostream& o) const;
  virtual void printAll(std::ostream& o) const;
  virtual const std::vector< DetSurface* >* sides() const
      { return _sides; }

  // Modifiers
  virtual void createCache();
  virtual void updateCache();

protected:

  // Helper functions
  virtual HepPoint coordToPoint( const TypeCoord* aCoord ) const;

private:
  std::vector< DetSurface* >* _sides;
// prohibit
  DetVolumeElem& operator = (const DetVolumeElem&);
  DetVolumeElem(const DetVolumeElem&);
};

std::ostream& operator<<( std::ostream& o, const DetVolumeElem& );

#endif // DETVOLUMEELEM_HH

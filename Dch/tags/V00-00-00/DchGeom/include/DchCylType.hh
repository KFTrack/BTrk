//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DchCylType.hh 123 2010-04-29 14:41:45Z stroili $
//
// Description:
//	Class DchCylType
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	R. Stroili		originator
//	
//
// Copyright Information:
//	Copyright (C) 1997	INFN-Pd
//
//------------------------------------------------------------------------

#ifndef DCHCYLTYPE_HH
#define DCHCYLTYPE_HH

//-------------
// C Headers --
//-------------

//---------------
// C++ Headers --
//---------------

//----------------------
// Base Class Headers --
//----------------------
#include "DetectorModel/DetSurfaceType.hh"

//-------------------------------
// Collaborating Class Headers --
//-------------------------------

//------------------------------------
// Collaborating Class Declarations --
//------------------------------------
class DchSimpleCyl;

//		---------------------
// 		-- Class Interface --
//		---------------------

class DchCylType : public DetSurfaceType {

//--------------------
// Declarations     --
//--------------------

    // Typedefs, consts, and enums

//--------------------
// Instance Members --
//--------------------

public:

  // Constructors
  DchCylType(const char* name, double radii[2], double length,
	     const DetMaterial*, int);

  DchCylType(const char* name, const DchSimpleCyl* cyl, const DetMaterial*, int);

  // Destructor
  virtual ~DchCylType( );
  
  // Operators
  //
  //  real versions of the virtual functions
  //
  bool physicalMaterial(const TypeCoord*) const ;
  //   bool activeMaterial(const TypeCoord* coord) const {return false;}
  //   const TypeCoord* physicalOutline(int&) const;
  //   const TypeCoord* activeOutline(int& ncorn) const {ncorn = 0; return 0;}
  const DetMaterial& material(const TypeCoord*) const;
  double thickness(const TwoDCoord*) const;
  double inRad(void) { return _radius; }
  double outRad(void) { return _radius + _thick;}
  double length(void) { return _length; }
  //    DchCylType&       operator= ( const DchCylType& );
  //    virtual Boolean operator==( const DchCylType& ) const;
  //            Boolean operator!=( const DchCylType& ) const;
  
  // Selectors (const)
  
  // Modifiers
  friend class DchGasElem;

private:

  // Copy Constructor
  DchCylType( const DchCylType& );
  DchCylType& operator= ( const DchCylType& );

  // Friends
  
  // Data members
  double _thick;     // tube wall thickness
  double _length;    // tube wall length
  double _radius;    // tube wall inner radius
  //
  //  Approximate outline, in 2-D type coordinates
  //
  TwoDCoord* _perimeter;
  //
  //  Material description.  This is a hack, ultimately this should point
  //  back to a database object
  //
  const DetMaterial* _thematerial; // material of this type

};

#endif

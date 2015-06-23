//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DchHyperType.hh 123 2010-04-29 14:41:45Z stroili $
//
// Description:
//	Class DchHyperType
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	<Author1>		<originator/contributor etc.>
//	<Author2>		<originator/contributor etc.>
//
// Copyright Information:
//	Copyright (C) 1994	<Institution>
//
//------------------------------------------------------------------------

#ifndef DCHHYPERTYPE_HH
#define DCHHYPERTYPE_HH

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
#include "DetectorModel/DetSurfaceType.hh"

//-------------------------------
// Collaborating Class Headers --
//-------------------------------

//------------------------------------
// Collaborating Class Declarations --
//------------------------------------

//		---------------------
// 		-- Class Interface --
//		---------------------

class DchHyperType : public DetSurfaceType {

  //--------------------
  // Declarations     --
  //--------------------


  //--------------------
  // Instance Members --
  //--------------------

public:

  // Constructors
  DchHyperType(const char*, const DetMaterial*, int, double, double);

  // Copy Constructor

  // Destructor
  virtual
  ~DchHyperType();

  // Operators

  //    virtual Boolean operator==( const DchHyperType& ) const;
  //            Boolean operator!=( const DchHyperType& ) const;

  // Selectors (const)
  bool
  physicalMaterial(const TypeCoord*) const;
  //   virtual bool activeMaterial(const TypeCoord*) const { return false; }
  //   const TypeCoord* activeOutline(int& ncorn) const {ncorn = 0; return 0;}
  const DetMaterial&
  material(const TypeCoord*) const;
  double
  thickness(const TwoDCoord* point) const;
  double
      effectiveThickness(const TwoDCoord* point, double dircos1, double dircos2) const;
  double
  getTwist()
  {
    return _twist;
  }

private:

  // Data members
  double _twist;
  double _length;
  TwoDCoord* _perimeter;
  const DetMaterial* _thematerial; // material of this type

  DchHyperType(const DchHyperType&);
  DchHyperType&
  operator=(const DchHyperType&);

};

#endif

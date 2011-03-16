//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DchCellPlaneType.hh 123 2010-04-29 14:41:45Z stroili $
//
// Description:
//	Class DchCellPlaneType.
//      Do not use this for DchCellPlaneTyped class (foo<T>).  
//      use DchCellPlaneTypeDchCellPlaneType.hh
//      instead.
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

#ifndef DCHCELLPLANETYPE_HH
#define DCHCELLPLANETYPE_HH

//-------------
// C Headers --
//-------------
extern "C" {
}

//---------------
// C++ Headers --
//---------------
#include <vector>

//----------------------
// Base Class Headers --
//----------------------
#include "DetectorModel/DetSurfaceType.hh"

//		---------------------
// 		-- Class Interface --
//		---------------------

class DchCellPlaneType : public DetSurfaceType {

public:

  // Constructors
  DchCellPlaneType(const char *, int id, std::vector<TwoDCoord>& corners);

    // Copy Constructor
//    DchCellPlaneType( const DchCellPlaneType& );

    // Destructor
  virtual ~DchCellPlaneType( );

//  real versions of the virtual functions
//
  bool physicalMaterial(const TypeCoord*) const;
  double thickness(const TwoDCoord*) const;
  const DetMaterial& material(const TypeCoord*) const;
    // Operators
    
//    DchCellPlaneType&       operator= ( const DchCellPlaneType& );
//    virtual int operator==( const DchCellPlaneType& ) const;
//            int operator!=( const DchCellPlaneType& ) const;

private:

  // Data members
  TwoDCoord _incorner;
  TwoDCoord _outcorner; // describe smallest enclosing rectangle
  double _thick; // describe the thickness
  const DetMaterial* _material;
};

#endif // DCHCELLPLANETYPE_HH

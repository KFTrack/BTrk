//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DchVolType.hh 123 2010-04-29 14:41:45Z stroili $
//
// Description:
//	Class DchVolType
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	R. Stroili		originator
//	
//
// Copyright Information:
//	Copyright (C) 1997	INFN - Pd
//
//------------------------------------------------------------------------

#ifndef DCHVOLTYPE_HH
#define DCHVOLTYPE_HH

//-------------
// C Headers --
//-------------

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

//------------------------------------
// Collaborating Class Declarations --
//------------------------------------
class HepPoint;

//		---------------------
// 		-- Class Interface --
//		---------------------

class DchVolType : public DetVolumeType {

//--------------------
// Declarations     --
//--------------------

  // Typedefs, consts, and enums

//--------------------
// Instance Members --
//--------------------

public:

  // Constructors
  DchVolType(const char* name, double rmin, double rmax, double zmin, double zmax,
	     const DetMaterial*, int);


  // Destructor
  virtual ~DchVolType( );

  // Operators
    
//    virtual Boolean operator==( const DchVolType& ) const;
//            Boolean operator!=( const DchVolType& ) const;

  // Selectors (const)
  const DetMaterial& material( const TypeCoord* here ) const;
  bool insideLimitsOf( int aSide,const SurfacePoint& thisPoint ) const;
  
  // Modifiers
  void setDebug(bool deb) { _debug = deb; }

  // Accessors
  double rmin(void)  const { return _rmin;  }
  double rmax(void)  const { return _rmax;  }
  double zmin(void)  const { return _zmin;  }
  double zmax(void)  const { return _zmax;  }
  bool   debug(void) const { return _debug; }

protected:

  // Helper functions
  bool insideLine( const SurfacePoint& thisPoint,
		   const SurfacePoint& p1,
		   const SurfacePoint& p2 ) const;

private:

  // Friends
  friend class DchVolElem;

  // Data members
  double _rmin;
  double _rmax;
  double _zmin;
  double _zmax;
  const DetMaterial* _theMaterial; // material of this type
  bool _debug;

  // Copy Constructor
  DchVolType( const DchVolType& );
  DchVolType&       operator= ( const DchVolType& );
};

#endif

//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DchPhiSegmVolType.hh 123 2010-04-29 14:41:45Z stroili $
//
// Description:
//	Class DchPhiSegmVolType
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

#ifndef DCHPHISEGMVOLTYPE_HH
#define DCHPHISEGMVOLTYPE_HH

//-------------
// C Headers --
//-------------

//---------------
// C++ Headers --
//---------------

//----------------------
// Base Class Headers --
//----------------------
#include "DchGeom/DchVolType.hh"

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

class DchPhiSegmVolType : public DchVolType {

//--------------------
// Declarations     --
//--------------------

  // Typedefs, consts, and enums

//--------------------
// Instance Members --
//--------------------

public:

  // Constructors
  DchPhiSegmVolType(const char* name, double rmin, double rmax, double zmin, double zmax,
                  double phi0, double solidDeltaPhi, double hollowDeltaPhi,
	     const DetMaterial*, int);


  // Destructor
  virtual ~DchPhiSegmVolType( );

  // Operators
    
//    virtual Boolean operator==( const DchPhiSegmVolType& ) const;
//            Boolean operator!=( const DchPhiSegmVolType& ) const;

  // Selectors (const)
//  bool insideLimitsOf( int aSide,const SurfacePoint& thisPoint ) const;
  
  double phi0(void)  const { return _phi0;  }
  double solidDeltaPhi(void)  const { return _solidDeltaPhi;  }
  double hollowDeltaPhi(void)  const { return _hollowDeltaPhi;  }
  double totDeltaPhi(void) const { return _deltaPhi;  }
  double maxNRotation(void) const { return _nRot;  }

protected:

  // Helper functions
//  bool insideLine( const SurfacePoint& thisPoint,
//		   const SurfacePoint& p1,
//		   const SurfacePoint& p2 ) const;

private:

  // Friends
  friend class DchVolElem;

  // Data members
  double _phi0;
  double _solidDeltaPhi;
  double _hollowDeltaPhi;
  double _deltaPhi;
  double _nRot;

  // Copy Constructor
  DchPhiSegmVolType( const DchPhiSegmVolType& );
  DchPhiSegmVolType&       operator= ( const DchPhiSegmVolType& );
};

#endif

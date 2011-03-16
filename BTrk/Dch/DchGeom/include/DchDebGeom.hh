//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DchDebGeom.hh 91 2010-01-14 12:37:23Z stroili $
//
// Description:
//	Class DchDebGeom
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

#ifndef DCHDEBGEOM_HH
#define DCHDEBGEOM_HH

//-------------
// C Headers --
//-------------

//---------------
// C++ Headers --
//---------------

//----------------------
// Base Class Headers --
//----------------------
#include "Framework/AppModule.hh"

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "Framework/AbsParmBool.hh"

//------------------------------------
// Collaborating Class Declarations --
//------------------------------------

//		---------------------
// 		-- Class Interface --
//		---------------------

class DchDebGeom : public AppModule {

//--------------------
// Declarations     --
//--------------------

    // Typedefs, consts, and enums

//--------------------
// Instance Members --
//--------------------

public:

    // Constructors
  DchDebGeom( const char* const theName, const char* const theDescription );

    // Copy Constructor
//    DchDebGeom( const DchDebGeom& );

  // Destructor
  virtual ~DchDebGeom( );

  // Operators
  virtual AppResult           beginJob( AbsEvent* anEvent );
  // dummies 
  virtual AppResult           endJob( AbsEvent* anEvent );
    
//    DchDebGeom&       operator= ( const DchDebGeom& );
//    virtual Boolean operator==( const DchDebGeom& ) const;
//            Boolean operator!=( const DchDebGeom& ) const;

    // Selectors (const)

    // Modifiers

protected:

    // Helper functions

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

};


#endif

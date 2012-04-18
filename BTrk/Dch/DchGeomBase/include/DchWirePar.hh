//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DchWirePar.hh 92 2010-01-14 12:38:30Z stroili $
//
// Description:
//	Class DchWirePar.
//      class used to store basic wire information in the database
//      
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	R. Stroili		originator
//	
//
// Copyright Information:
//	Copyright (C) 1998	INFN & Padova University
//
//------------------------------------------------------------------------

#ifndef DCHWIREPAR_HH
#define DCHWIREPAR_HH

//-------------
// C Headers --
//-------------

//---------------
// C++ Headers --
//---------------
#include <iostream>

//----------------------
// Base Class Headers --
//----------------------

//-------------------------------
// Collaborating Class Headers --
//-------------------------------

//------------------------------------
// Collaborating Class Declarations --
//------------------------------------

//		---------------------
// 		-- Class Interface --
//		---------------------

class DchWirePar {

//--------------------
// Declarations     --
//--------------------

    // Typedefs, consts, and enums

//--------------------
// Instance Members --
//--------------------

public:

  // Constructors
  DchWirePar();
  DchWirePar( int id );
  DchWirePar( int id, double radR, double phiR, double zR, double radF, 
	   double phiF, double zF, double sag=0., double diameter=0. );
  DchWirePar( int id, double radR, double phiR, double zR, double length,
	   double twist, double sag=0., double diameter=0. );
  DchWirePar( const DchWirePar& );       // Copy Constructor

  // Destructor
  virtual ~DchWirePar( );

  // Operators
  DchWirePar& operator= ( const DchWirePar& );  // Assignment op
    
  int operator==( const DchWirePar& ) const;
//            int operator!=( const DchWirePar& ) const;

  void printAll(std::ostream& o=std::cout) const ;

  // Selectors (const)
  int getID(void) const { return _ID; }
  double getRadRear(void) const { return _radRear; }
  double getPhiRear(void) const { return _phiRear; }
  double getZRear(void) const { return _zRear; }
  double getRadForw(void) const { return _radForw; }
  double getPhiForw(void) const { return _phiForw; }
  double getZForw(void) const { return _zForw; }
  double getSag(void) const { return _sag; }
  double getDiameter(void) const { return _diameter; }

  double twist(void) const { return _phiForw - _phiRear; }
  double zLength(void) const { return _zForw - _zRear; }
  double stereo(void) const ;


private:
  friend class DchGDchCmpr;
//  friend class DchWireAlignCmpr;
  friend class DchPlateDeflCmpr;
  // Data members
  int    _ID;                    // wire identifier
  double _radRear;               // radius on rear end plate
  double _phiRear;               // phi angle on rear end plate
  double _zRear;                 // z position on rear end plate
  double _radForw;               // radius on forward end plate
  double _phiForw;               // phi angle on forward end plate
  double _zForw;                 // z position on forward end plate
  double _sag;                   // sag of wire
  double _diameter;              // wire diameter
};

std::ostream&  operator<<(std::ostream& o, DchWirePar& cell);

#endif // DCHWIREPAR_HH

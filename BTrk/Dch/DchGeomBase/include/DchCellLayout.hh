//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DchCellLayout.hh 92 2010-01-14 12:38:30Z stroili $
//
// Description:
//	Class DchCellLayout.
//      Class used to store Dch cell layout in the database
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

#ifndef DCHCELLLAYOUT_HH
#define DCHCELLLAYOUT_HH

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "DchGeomBase/DchWirePar.hh"

//		---------------------
// 		-- Class Interface --
//		---------------------

#define MAXWIRES 7

class DchCellLayout {

//--------------------
// Declarations     --
//--------------------

    // Typedefs, consts, and enums

//--------------------
// Instance Members --
//--------------------

public:

    // Constructors
  DchCellLayout(void);
  DchCellLayout(int id, int nwires, const DchWirePar* wires);

  // Destructor
  virtual ~DchCellLayout( );

  // Operators
  DchCellLayout& operator= ( const DchCellLayout& );  // Assignment op
    
  // Selectors (const)
  int ID(void) const { return _ID; }
  int nofWires(void) const { return _nofWires; }
  const DchWirePar& fWire(int i) const { return _fieldWires[i]; }

private:
  friend class DchGDchCmpr;

  // Data members
  int _ID;                           // cell ID (?)
  int _nofWires;                     // number of field wires used to describe 
                                     // the cell
  DchWirePar _fieldWires[MAXWIRES];  // field wires describing cell layout

};

#endif // DCHCELLLAYOUT_HH

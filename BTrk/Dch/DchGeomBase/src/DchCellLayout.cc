//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DchCellLayout.cc 92 2010-01-14 12:38:30Z stroili $
//
// Description:
//	Class DchCellLayout
//      Dch cell layout using field wires
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
#include "BaBar/BaBar.hh"

//-----------------------
// This Class's Header --
//-----------------------
#include "DchGeomBase/DchCellLayout.hh"

//----------------
// Constructors --
//----------------
DchCellLayout::DchCellLayout()
  : _ID(0), _nofWires(0)
{
}

DchCellLayout::DchCellLayout( int id, int nofWires, const DchWirePar* wires )
  : _ID(id), _nofWires(nofWires)
{
  for (int i=0; i<nofWires; i++) {
    _fieldWires[i] = wires[i];
  }
}

//--------------
// Destructor --
//--------------
DchCellLayout::~DchCellLayout()
{}

//-------------
// Operators --
//-------------
DchCellLayout&       
DchCellLayout::operator= ( const DchCellLayout& cell)
{
  if (&cell == this) return *this;
  
  _ID =  cell.ID();
  _nofWires = cell.nofWires();
  for (int i=0; i<_nofWires; i++) {
    _fieldWires[i] = cell.fWire(i);
  }
  return *this;
}
    

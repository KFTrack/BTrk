//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DchSimpleCyl.cc 92 2010-01-14 12:38:30Z stroili $
//
// Description:
//	Class DchSimpleCyl
//      Do not use this for DchSimpleCyld class (foo<T>).  use DchSimpleCylDchSimpleCyl.hh
//      instead.
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
#include "DchGeomBase/DchSimpleCyl.hh"

//-------------
// C Headers --
//-------------
extern "C" {
}

//----------------
// Constructors --
//----------------
DchSimpleCyl::DchSimpleCyl()
  : _innerRadius(0.), _outerRadius(0.), _length(0.), _isPhiSegmented(false)
{
}

DchSimpleCyl::DchSimpleCyl( double inR, double outR, double length, bool isPhiSegmented )
  : _innerRadius(inR), _outerRadius(outR), _length(length), _isPhiSegmented(isPhiSegmented)
{
}

DchSimpleCyl::DchSimpleCyl( const DchSimpleCyl& c ) 
{
  _innerRadius = c.getInnerRadius();
  _outerRadius = c.getOuterRadius();
  _length = c.getLength();
  _isPhiSegmented = c.isPhiSegmented();
}

//--------------
// Destructor --
//--------------
DchSimpleCyl::~DchSimpleCyl() 
{}

//-------------
// Operators --
//-------------
DchSimpleCyl& 
DchSimpleCyl::operator= ( const DchSimpleCyl& c ) 
{
  _innerRadius = c.getInnerRadius();
  _outerRadius = c.getOuterRadius();
  _length = c.getLength();
  _isPhiSegmented = c.isPhiSegmented();

  return *this;
}
    

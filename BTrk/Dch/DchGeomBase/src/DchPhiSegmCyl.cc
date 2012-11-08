//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id:  $
//
// Description:
//	Class DchPhiSegmCyl
//      Do not use this for DchPhiSegmCyl class (foo<T>).  use DchPhiSegmCyl.hh
//      instead.
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	
//
//
// Copyright Information:
//
//
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"

//-----------------------
// This Class's Header --
//-----------------------
#include "DchGeomBase/DchPhiSegmCyl.hh"

//-------------
// C Headers --
//-------------
extern "C" {
}

//----------------
// Constructors --
//----------------
DchPhiSegmCyl::DchPhiSegmCyl()
  : DchSimpleCyl(), _phi0(0.), _solidDeltaPhi(0.), _hollowDeltaPhi(0.)
{
}

DchPhiSegmCyl::DchPhiSegmCyl( double inR, double outR, double length)
  : DchSimpleCyl(inR, outR, length, false),
    _phi0(0.), _solidDeltaPhi(0.), _hollowDeltaPhi(0.)
{
}

DchPhiSegmCyl::DchPhiSegmCyl( double inR, double outR, double length,
                double phi0, double solidDeltaPhi, double hollowDeltaPhi )
  : DchSimpleCyl(inR, outR, length, true),
    _phi0(phi0), _solidDeltaPhi(solidDeltaPhi), _hollowDeltaPhi(hollowDeltaPhi)
{
}

DchPhiSegmCyl::DchPhiSegmCyl( const DchPhiSegmCyl& c )
{
  _innerRadius = c.getInnerRadius();
  _outerRadius = c.getOuterRadius();
  _length = c.getLength();
  _isPhiSegmented = c.isPhiSegmented();
  _phi0 = c.getPhi0();
  _solidDeltaPhi = c.getSolidDeltaPhi();
  _hollowDeltaPhi = c.getHollowDeltaPhi();
}

//--------------
// Destructor --
//--------------
DchPhiSegmCyl::~DchPhiSegmCyl()
{}

//-------------
// Operators --
//-------------
DchPhiSegmCyl&
DchPhiSegmCyl::operator= ( const DchPhiSegmCyl& c )
{
  _innerRadius = c.getInnerRadius();
  _outerRadius = c.getOuterRadius();
  _length = c.getLength();
  _isPhiSegmented = c.isPhiSegmented();
  _phi0 = c.getPhi0();
  _solidDeltaPhi = c.getSolidDeltaPhi();
  _hollowDeltaPhi = c.getHollowDeltaPhi();

  return *this;
}
    

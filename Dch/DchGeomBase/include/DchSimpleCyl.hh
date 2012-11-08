//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DchSimpleCyl.hh 92 2010-01-14 12:38:30Z stroili $
//
// Description:
//	Class DchSimpleCyl.
//      Do not use this for DchSimpleCyld class (foo<T>).  use DchSimpleCylDchSimpleCyl.hh
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
//	Copyright (C) 1997	<Institution>
//
//------------------------------------------------------------------------

#ifndef DCHSIMPLECYL_HH
#define DCHSIMPLECYL_HH

//-------------
// C Headers --
//-------------
extern "C" {
}

//		---------------------
// 		-- Class Interface --
//		---------------------

class DchSimpleCyl {

public:

  // Constructors
  DchSimpleCyl();
  DchSimpleCyl( double inR, double outR, double length, bool isPhiSegmented=false );
  DchSimpleCyl( const DchSimpleCyl& );       // Copy Constructor

  // Destructor
  virtual ~DchSimpleCyl( );

  // Operators
  DchSimpleCyl& operator= ( const DchSimpleCyl& );  // Assignment op
  
  // Selectors (const)
  double getInnerRadius(void) const { return _innerRadius; }
  double getOuterRadius(void) const { return _outerRadius; }
  double getLength(void) const { return _length; }
  bool   isPhiSegmented(void) const { return _isPhiSegmented; }

protected:

  friend class DchGDchCmpr;
  //Note: if your class needs a copy constructor or an assignment operator, 
  //  make one of the following public and implement it.

  // Data members
  double _innerRadius;         // inner radius
  double _outerRadius;         // outer radius
  double _length;              // cylinder length
  bool   _isPhiSegmented;

};

#endif // DCHSIMPLECYL_HH

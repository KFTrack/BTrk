//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id:  $
//
// Description:
//	Class DchPhiSegmCyl.
//      Do not use this for DchPhiSegmCyld class (foo<T>).  use DchPhiSegmCylDchPhiSegmCyl.hh
//      instead.
//
// Environment:
//
//
// Author List:
//	<Author1>		<originator/contributor etc.>
//	<Author2>		<originator/contributor etc.>
//
// Copyright Information:
//	Copyright (C) 1997	<Institution>
//
//------------------------------------------------------------------------

#ifndef DCHPHISEGMCYL_HH
#define DCHPHISEGMCYL_HH

#include "DchGeomBase/DchSimpleCyl.hh"

//-------------
// C Headers --
//-------------
extern "C" {
}

//		---------------------
// 		-- Class Interface --
//		---------------------

class DchPhiSegmCyl : public DchSimpleCyl {

public:

  // Constructors
  DchPhiSegmCyl();
  DchPhiSegmCyl( double inR, double outR, double length );
  DchPhiSegmCyl( double inR, double outR, double length,
                 double phi0, double solidDeltaPhi, double hollowDeltaPhi );
  DchPhiSegmCyl( const DchPhiSegmCyl& );       // Copy Constructor

  // Destructor
  virtual ~DchPhiSegmCyl( );

  // Operators
  DchPhiSegmCyl& operator= ( const DchPhiSegmCyl& );  // Assignment op
  
  // Selectors (const)
  double getPhi0(void) const { return _phi0; }
  double getSolidDeltaPhi(void) const { return _solidDeltaPhi; }
  double getHollowDeltaPhi(void) const { return _hollowDeltaPhi; }

private:

  friend class DchGDchCmpr;
  //Note: if your class needs a copy constructor or an assignment operator, 
  //  make one of the following public and implement it.

  // Data members
  double _phi0;           // phi0
  double _solidDeltaPhi;  // delta phi of the solid part
  double _hollowDeltaPhi; // delta phi of the hollow parts

};

#endif // DCHSIMPLECYL_HH

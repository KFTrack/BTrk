//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DchMCEff.hh 88 2010-01-14 12:32:57Z stroili $
//
// Description:
//	Class DchMCEff
//      Base class for DchClusEff, DchConstEff, DchDBEff
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	S. Sen		originator
//	M. van Hoek     cleanup
//	
//
// Copyright Information:
//	Copyright (C) 1998	University of Colorado
//
//------------------------------------------------------------------------

#ifndef DCHMCEFF_HH
#define DCHMCEFF_HH

class DchGHit;

class DchMCEff {
public:

 DchMCEff();
 virtual  ~DchMCEff();

 // raison d'etre...
 virtual bool keepHit(unsigned layer, unsigned wire, double pathLength, double Q) const = 0;
 virtual const char* name() const = 0;

protected:

private:

  //Note: if your class needs a copy constructor or an assignment operator, 
  //  make one of the following public and implement it.
  DchMCEff( const DchMCEff& );       // Copy Constructor
  DchMCEff&       operator= ( const DchMCEff& );  // Assignment op
};
#endif // DCHEFF_HH

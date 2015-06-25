//==========================================================================
// File and Version Information:
// 	$Id: BbrDoubleErr.cc 491 2010-01-13 16:59:16Z stroili $
//
//--------------------------------------------------------------------------
// Description:
//	See BbrDoubleErr.hh
//
//--------------------------------------------------------------------------
// Sample User Code:
//
//--------------------------------------------------------------------------
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
//--------------------------------------------------------------------------
// Author List:
//	Abi Soffer              (Original author)
//
//--------------------------------------------------------------------------
// Copyright Information:
//	Copyright (C) 1998	Colorado State University
//
//==========================================================================

//----------------
// BaBar Header --
//----------------
#include "BTrk/BaBar/BaBar.hh"

//-----------------------
// This Class's Header --
//-----------------------

#include "BTrk/BbrGeom/BbrDoubleErr.hh"

//-------------
// C Headers --
//-------------
extern "C" {
}

//---------------
// C++ Headers --
//---------------

#include <math.h>
#include <iostream>
using std::ostream;

//-------------------------------
// Collaborating Class Headers --
//-------------------------------

//-----------------------------------------------------------------------
// Local Macros, Typedefs, Structures, Unions and Forward Declarations --
//-----------------------------------------------------------------------

static const char rscid[] = "$Id: BbrDoubleErr.cc 491 2010-01-13 16:59:16Z stroili $";

//-----------------------------------------------
//-- Static Data & Function Member Definitions --
//-----------------------------------------------

//----------------------------------------
//-- Public Function Member Definitions --
//----------------------------------------

//----------------
// Constructors --
//----------------
// Copy Constructor:

BbrDoubleErr::BbrDoubleErr(const BbrDoubleErr & original) {
  *this = original;
}

//-----------------------
// Assignment Operator --
//-----------------------

BbrDoubleErr & BbrDoubleErr::operator = (const BbrDoubleErr & original) {
  _value = original._value;
  _covariance = original._covariance;
  return *this;
}
    
//-------------
// Accessors --
//-------------
    
double 
BbrDoubleErr::determineChisq(double ref) const {
  if (0 >= _covariance){
    return BbrError::chisqUndef;
  }
  else{
    double diff = (ref - _value);
    return (diff * diff) / _covariance;
  }
}

//--------------------------------------------------------------------
// unary -:
BbrDoubleErr BbrDoubleErr::operator-() {
  return BbrDoubleErr(-_value, _covariance);
}

//--------------------------------------------------------------------
BbrDoubleErr BbrDoubleErr::operator - (const BbrDoubleErr & bde) {
  return BbrDoubleErr(_value - bde._value, _covariance - bde._covariance);
}

//--------------------------------------------------------------------
BbrDoubleErr BbrDoubleErr::operator + (const BbrDoubleErr & bde) {
  return BbrDoubleErr(_value + bde._value, _covariance + bde._covariance);
}

//--------------------------------------------------------------------
BbrDoubleErr BbrDoubleErr::operator * (const BbrDoubleErr & bde) {
  return BbrDoubleErr(_value * bde._value, 
		      _covariance * bde._value * bde._value +
		      bde._covariance * _value * _value);
}

//--------------------------------------------------------------------
BbrDoubleErr BbrDoubleErr::operator / (const BbrDoubleErr & bde) {
  register double bde2 = bde._value * bde._value;

  return BbrDoubleErr(_value / bde._value,
		      _covariance / bde2 +
		      bde._covariance * _value * _value / (bde2 * bde2));
}		      

//-------------
// Modifiers --
//-------------
 
BbrDoubleErr & BbrDoubleErr::operator += (const BbrDoubleErr & bde) {
  *this = *this + bde;
  return *this;
}

//-----------------------------------------------------------------------
BbrDoubleErr & BbrDoubleErr::operator -= (const BbrDoubleErr & bde) {
  *this = *this - bde;
  return *this;
}

//-----------------------------------------------------------------------
BbrDoubleErr & BbrDoubleErr::operator *= (const BbrDoubleErr & bde) {
  *this = *this * bde;
  return *this;
}

//-----------------------------------------------------------------------
BbrDoubleErr & BbrDoubleErr::operator /= (const BbrDoubleErr & bde) {
  *this = *this / bde;
  return *this;
}


//-----------
// Globals --
//-----------

ostream & operator<<(ostream & stream, const BbrDoubleErr & bde) {
  stream << "value: " << bde.value() 
	 << " covariance: " << bde.covariance();

  return stream;
}
    

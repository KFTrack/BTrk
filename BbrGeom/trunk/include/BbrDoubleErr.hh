//==========================================================================
// File and Version Information:
// 	$Id: BbrDoubleErr.hh 491 2010-01-13 16:59:16Z stroili $
//
//--------------------------------------------------------------------------
// Description:
// Class BbrDoubleErr holds a double and its error squared,
// the equivalent of a BbrVectorErr with one dimension.
//
// determineChisq(double ref) returns (ref-value())^2 / covariance().
// If the covariance() <= 0 then returns BbrError::chisqUndef.
//
// The mathematical operators can be used to correctly take into account 
// error propagation through simple mathematical operations. Beware of 
// operations whose covariance may not be the same as the one which 
// results from a sequence of simple operations. For example, a*a will give
// the wrong covariance, since the algorithm assumes that the two arguments
// of the operator* are independent.
// 
//--------------------------------------------------------------------------
// Collaborating classes:
//
//--------------------------------------------------------------------------
// Sample User Code:
//
//--------------------------------------------------------------------------
// Compiler Default Member Functions:
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

#ifndef BBRDOUBLEERR_HH
#define BBRDOUBLEERR_HH

//-------------
// C Headers --
//-------------
extern "C" {
}

//---------------
// C++ Headers --
//---------------

//-----------------
// BaBar Headers --
//-----------------
#include "BaBar/BaBar.hh"

//----------------------
// Base Class Headers --
//----------------------

#include "BbrGeom/BbrError.hh"       // for chisqUndef only

//-------------------------------
// Collaborating Class Headers --
//-------------------------------

//--------------------------------------------
// Collaborating Class Forward Declarations --
//--------------------------------------------

#include <iosfwd>

//		---------------------
// 		-- Class Interface --
//		---------------------

class BbrDoubleErr {
public:

  // Constructors: The _covariance and the _value are 0.0 unless specified:
  BbrDoubleErr() :                       _value(0.0), _covariance(0.0) {}
  BbrDoubleErr(double val) :             _value(val), _covariance(0.0) {}
  BbrDoubleErr(double val, double cov) : _value(val), _covariance(cov) {}

  // Copy Constructor
  BbrDoubleErr(const BbrDoubleErr &);

  // Destructor
  virtual ~BbrDoubleErr() {}

  // Assignment operator:
  BbrDoubleErr & operator = (const BbrDoubleErr &);

  // Accessors (const)
  double value()                    const {return _value;}
  double covariance()               const {return _covariance;}
  double determineChisq(double ref) const;

  BbrDoubleErr operator-();   // value() -> -value(), covariance() unaffected
  BbrDoubleErr operator - (const BbrDoubleErr &);
  BbrDoubleErr operator + (const BbrDoubleErr &);

  // NOTE: (a * b).covariance() is 
  // b^2 * a.covariance() + a^2 * b.covariance()
  BbrDoubleErr operator * (const BbrDoubleErr &);

  // NOTE: (a / b).covariance() is 
  // a.covariance() / b^2 + b.covariance() * a^2 / b^4
  BbrDoubleErr operator / (const BbrDoubleErr &);

  // modifiers:
  void set_value(double val)      {_value = val;}
  void set_covariance(double cov) {_covariance = cov;}

  BbrDoubleErr & operator += (const BbrDoubleErr &);
  BbrDoubleErr & operator -= (const BbrDoubleErr &);
  BbrDoubleErr & operator *= (const BbrDoubleErr &);
  BbrDoubleErr & operator /= (const BbrDoubleErr &);

  // needed for RWTValOrderedVector
  bool           operator == (const BbrDoubleErr & other) const 
      { return (_value == other._value && _covariance == other._covariance); }
  bool           operator <   (const BbrDoubleErr & other) const 
      { return (_value < other._value); }


private:
  // Data members
  double _value;
  double _covariance;
};

// globals:

std::ostream & operator<<(std::ostream & stream, const BbrDoubleErr & bde);

#endif

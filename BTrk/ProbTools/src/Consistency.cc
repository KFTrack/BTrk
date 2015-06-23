//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: Consistency.cc 458 2010-01-15 11:37:35Z stroili $
//
// Description:
//	Class Consistency; encapsulates statistics
//
//      See header for more info.  This implementation
//      is not complete
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	Bob Jacobsen		Original Author
//
//      Dec 2007, Alexandre Telnov (Princeton)
//            Added member _logLikelihood and accessor logLikelihood() 
//
// Copyright Information:
//	Copyright (C) 1995	Lawrence Berkeley Laboratory
//
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"

//------------------
// This C Headers --
//------------------
extern "C" {
#include <math.h>
#include <assert.h>
#include <float.h> 
}

//--------------------
// This C++ Headers --
//--------------------
#include <iostream>
#include <iomanip>

//-----------------------
// This Class's Header --
//-----------------------
#include "ProbTools/Consistency.hh"
using std::endl;
using std::ios;
using std::ostream;
using std::setiosflags;
using std::setprecision;
using std::setw;


//              ----------------------------------------
//              -- Public Function Member Definitions --
//              ----------------------------------------


Consistency::Consistency() 
  : _stat(noMeasure), _value(0), _likelihood(1.), _logLikelihood(0.), _sign(unknown)
{
}

Consistency::Consistency(double v, double l)
  : _stat(OK), _value(v), _likelihood(l), _logLikelihood(likeLog(l)), _sign(unknown)
{
}

Consistency::Consistency(const Consistency& rhs) 
  : _stat(rhs._stat), _value(rhs._value)
  , _likelihood(rhs._likelihood), _logLikelihood(rhs._logLikelihood), _sign(rhs._sign)
{
}

// operators
Consistency& 
Consistency::operator=(const Consistency& rhs)
{
  if(this != &rhs){
    _stat=rhs._stat;
    _value=rhs._value;
    _likelihood=rhs._likelihood;
    _logLikelihood=rhs._logLikelihood;
    _sign=rhs._sign;
  }
  return *this;
}

bool
Consistency::operator==(const Consistency& rhs) const
{
  bool answer = false;
  if (_stat       == rhs._stat &&
      _value      == rhs._value &&
//      _likelihood == rhs._likelihood &&
      ( ( _likelihood > 5*DBL_MIN && _likelihood == rhs._likelihood ) // DBL_MIN = 2.2250738585072014e-308 
	|| ( _likelihood <= 5*DBL_MIN && _logLikelihood == rhs._logLikelihood ) ) &&
      _sign       == rhs._sign) {
    answer = true;
  }

  return answer;
}

bool
Consistency::operator<(const Consistency& rhs) const
{
  //
  // this is not unique...
  //

  bool answer = false;
  if (_stat == OK && rhs._stat == OK ) {
    //
    // Dare to use consistency for comparison. 
    // Likelihood may give a different answer
    //
    if (_value < rhs._value) {
      answer = true;
    }
  } else {
    if (rhs._stat == OK) {
      answer = true;
    }
  }

  return answer;
}

bool
Consistency::operator>(const Consistency& rhs) const
{
  // ghm
  if( *this==rhs ) return false;
  return( ! (*this<rhs) );
}


const ConsistencySet* 
Consistency::genealogy() const
{
  return 0;
}


void
Consistency::print(ostream& os) const
{
  os << "Likelihood = ";
  os << setiosflags(ios::fixed) << setw(7) << setprecision(4);
  os << _likelihood;
  os << "log(Likelihood) = ";
  os << setiosflags(ios::fixed) << setw(7) << setprecision(4);
  os << _logLikelihood;
  os << "\t SignificanceLevel = ";
  os << setiosflags(ios::fixed) << setw(7) << setprecision(4);
  os << _value;
  os << "\t Status ";
  switch( status() ) {
  case Consistency::OK  :
    os << "OK"; 
    break;
  case  Consistency::noMeasure  :
    os << "noMeasure"; 
    break;
  case  Consistency::underFlow  :
    os << "underFlow"; 
    break;
  case  Consistency::unPhysical  :
    os << "unPhysical"; 
    break;
  default:
    os << "unknown";
  }
  os << endl;
}

//------------------
// Static methods --
//------------------

const Consistency&
Consistency::badMeasurement() 
{
  static Consistency bad(1.,1.);
  bad.setStatus(noMeasure);

  return bad;
}

// log(like) if argument is positive, -999 if 0, -9999 if negative 
double Consistency::likeLog(double like) {
  if (like>0.) return log(like);
  if (like==0.) return -999.;
  return -9999.;
}

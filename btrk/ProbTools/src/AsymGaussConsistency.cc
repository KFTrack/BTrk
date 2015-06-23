//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: AsymGaussConsistency.cc 458 2010-01-15 11:37:35Z stroili $
//
// Description:
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Yury Kolomensky
//      Alexandre Telnov, December 2007: 
//         fix bugs in setting the underflow status; add logLikelihood 
//
// Copyright Information:
//	Copyright (C) 1998       Caltech
//      2007  Princeton University
//
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"

//-----------------------
// This Class's Header --
//-----------------------
#include "BaBar/Constants.hh"
#include "ProbTools/AsymGaussConsistency.hh"

#include <assert.h>
#include <math.h> 
#include <float.h>

extern "C" double erfc(double);

//		----------------------------------------
// 		-- Public Function Member Definitions --
//		----------------------------------------


AsymGaussConsistency::AsymGaussConsistency( double delta, 
					    double sigmaMinus,
					    double sigmaPlus) 
  : _delta(delta)
  , _sigmaMinus(sigmaMinus) 
  , _sigmaPlus(sigmaPlus) 
{
  _value = 0.; _likelihood = 0.; _logLikelihood = -999.;
  if (delta < 0.) {
    _sign = Consistency::left;
  }
  if (delta > 0.) {
    _sign = Consistency::right;
  }
  if ( calc() ) {
    setStatus(Consistency::OK);
  } else {
    setStatus(Consistency::underFlow);
  }
}

bool AsymGaussConsistency::calc() {
  double sigma = _sigmaMinus;
  if ( _delta > 0 ) {
    sigma = _sigmaPlus;
  }

  if ( sigma == 0 ) return false;
  const double arg = _delta/sigma;

  _value = erfc(fabs(arg*M_SQRT1_2));

  const double arg2 = -0.5*arg*arg;
  const double norm = 0.5*sqrt(Constants::twoPi)*(_sigmaPlus+_sigmaMinus);

  _likelihood = exp(arg2)/norm;
  _logLikelihood = arg2 - log(norm);

// this sets the underflow flag in computing the likelihood
// _logLikelihood is correct in any case
  if ( arg2 < 0.71*DBL_MIN_EXP || _logLikelihood < 0.71*DBL_MIN_EXP ) return false; 
  
  return true;

}

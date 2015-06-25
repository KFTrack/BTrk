//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: GaussConsistency.cc 458 2010-01-15 11:37:35Z stroili $
//
// Description:
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	Bob Jacobsen, Ed Iskander
//      Paul Harrison:            Added sign of consistency
//      Alexandre Telnov, December 2007: 
//         add logLikelihood and method "double logErfC(double x)",
//         fix bugs in setting the underflow status  
//
// Copyright Information:
//	Copyright (C) 1996	Lawrence Berkeley Laboratory
//      2007  Princeton University
//
//------------------------------------------------------------------------

//-----------------------
// This Class's Header --
//-----------------------
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/BaBar/Constants.hh"
#include "BTrk/ProbTools/GaussConsistency.hh"

#include <assert.h>
#include <math.h> 
#include <float.h> 

extern "C" double erfc(double);

//		----------------------------------------
// 		-- Public Function Member Definitions --
//		----------------------------------------


GaussConsistency::GaussConsistency( double delta, double sigma) :
  _delta(delta), _sigma(sigma) 
{ 
  _value = 0.; _likelihood = 0.; _logLikelihood = -999.;
  if (delta < 0.) {
    _sign = Consistency::left;
  }
  if (delta > 0.) {
    _sign = Consistency::right;
  }
  if( calc() )
    setStatus(Consistency::OK);
  else
    setStatus(Consistency::underFlow); 
}

// Added by Alexandre Telnov, borrowed from RooFitModels/RooGExpModel.cc
double GaussConsistency::logErfC(double x)
{
// Approximated log(erfc(x))
  double t,z,ans;
  z=fabs(x);
  t=1.0/(1.0+0.5*z);
  
  if(x >= 0.0)
    ans=log(t)+(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+t*(-0.18628806+
        t*(0.27886807+t*(-1.13520398+t*(1.48851587+t*(-0.82215223+t*0.17087277)))))))));
  else
    ans=log(2.0-t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+t*(-0.18628806+
        t*(0.27886807+t*(-1.13520398+t*(1.48851587+t*(-0.82215223+t*0.17087277))))))))));

  return ans;
}

bool 
GaussConsistency::calc() 
{
  // Gautier 01/28/99 : add protection against unphysical values of arg
  // 
  // Alexandre Telnov, Dec 2007: comparing arg2 vs DBL_MIN_EXP == -1021, 
  // the way it was done until now, is incorrect. In the IEEE standard, 
  // exp(x) is a denormalized number if -745 <= x <= log(2)*DBL_MIN_EXP == -708,
  // zero if x < -745. Let's allow x down to -725, which is roughly half-way
  // through denormalized floating-point numbers. Please note that now we
  // calculate the likelihood even when it is a "very borderline" denormalized
  // number - but such marginal cases are covered by the underFlow flag.  
  if( _sigma==0. ) return false;
  const double arg = _delta/_sigma;

  _value = erfc(fabs(arg*M_SQRT1_2));

  const double arg2 = -0.5*arg*arg;
  const double norm = _sigma*sqrt(Constants::twoPi);

  _likelihood = exp(arg2)/norm;
  _logLikelihood = arg2 - log(norm);

// this sets the underflow flag in computing the likelihood
// _logLikelihood is correct in any reasonable case
  if ( arg2 < 0.71*DBL_MIN_EXP || _logLikelihood < 0.71*DBL_MIN_EXP ) return false; 

  return true;

}







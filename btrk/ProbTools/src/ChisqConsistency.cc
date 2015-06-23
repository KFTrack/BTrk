//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: ChisqConsistency.cc 560 2010-03-29 23:44:43Z brownd $
//
// Description:
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	Bob Jacobsen, Ed Iskander
//      Alexandre Telnov, July 2007: 
//         fix bugs in setting the underflow status; add logLikelihood 
//
// Copyright Information:
//	Copyright (C) 1996	Lawrence Berkeley Laboratory
//      2007  Princeton University
//
//------------------------------------------------------------------------

//-----------------------
// This Class's Header --
//-----------------------
#include "BaBar/BaBar.hh"
#include "ProbTools/ChisqConsistency.hh"
#include "TMath.h"

#include <iostream>
#include <math.h>
#include <float.h>

#include "ProbTools/NumRecipes.hh"

#include "ErrLogger/ErrLog.hh"

// prototype the cernlib function
//extern "C" {
//  float chisin_(const float&, const int&);
//}

ChisqConsistency::ChisqConsistency() :
  _chisq(-1.0), _nDof(0)
{}

ChisqConsistency::ChisqConsistency(double chisq, double nDof) :
   _chisq(chisq), _nDof(nDof)
{
  double z2 = 0.5*_chisq;
  double n2 = 0.5*_nDof;
  
  if (n2<=0 || z2<0) { 
//    ErrMsg(warning) << " Got unphysical values: chisq = " << chisq 
//		    << " #dof = " << nDof << endmsg;
    _value=0;
    _likelihood=0; _logLikelihood = -999.;
    setStatus(Consistency::unPhysical);
    return;
  }
  setStatus(Consistency::OK);
  
// given that n2>0 && z2>=0, gammq will NOT abort
  _value = NumRecipes::gammq(n2,z2);
  
  if (_chisq==0) { 
    _likelihood=1; _logLikelihood=0.;
  } else {
    const double loglike=(n2-1)*log(z2)-z2-NumRecipes::gammln(n2);
    _logLikelihood = loglike+log(0.5); // correct even if underFlow
    if ( loglike < 0.71*(DBL_MIN_EXP+1) ) {
    //_likelihood = 0; 
      _likelihood = 0.5*exp(loglike); // calculate but set underFlow flag
      setStatus(Consistency::underFlow);
    } else {
      _likelihood = 0.5*exp(loglike);
    }
  }
}



ChisqConsistency::ChisqConsistency(unsigned nDof, double prob) :
  _nDof(nDof)
{
  if(prob >= 0.0|| prob <= 1.0 || nDof < 0)
    _value = prob;
  else {
    ErrMsg(warning) << " Got unphysical values: prob = " << prob 
		    << " #dof = " << nDof << endmsg;
    _value=0;
    _likelihood=0; _logLikelihood = -999.;
    setStatus(Consistency::unPhysical);
    return;
  }
  setStatus(Consistency::OK);
  if(prob != 1.0){
// use the cernlib function to get chisq.  Note the funny convention on prob!  
    if(_value < FLT_EPSILON)
      _chisq = log(double(FLT_MAX));
    else {
      float value = 1.0f-float(_value);
      int ndof = nDof;
//      _chisq = chisin_(value,ndof);
      _chisq = TMath::ChisquareQuantile(value, ndof);
    }

// use the same algorithm as above to get loglikelihood
    double z2 = 0.5*_chisq;
    double n2 = 0.5*_nDof;
    if (_chisq==0) { 
      _likelihood=1; _logLikelihood=0.;
    } else {
      double loglike=(n2-1)*log(z2)-z2-NumRecipes::gammln(n2);
      _logLikelihood = loglike+log(0.5); // correct even if underFlow
      if ( loglike < 0.71*(DBL_MIN_EXP+1) ) {
      // _likelihood = 0; 
	_likelihood = 0.5*exp(loglike); // calculate anyway but set underFlow flag
        setStatus(Consistency::underFlow);
      } else {
	_likelihood = 0.5*exp(loglike);
      }
    }
  }
}


ChisqConsistency::ChisqConsistency(const ChisqConsistency& other) :
  Consistency(other), _chisq(other._chisq), _nDof(other._nDof)
{}

ChisqConsistency&
ChisqConsistency::operator =(const ChisqConsistency& other) {
  if(this != &other){
    Consistency::operator =(other);
    _chisq = other._chisq;
    _nDof = other._nDof;
  }
  return *this;
}

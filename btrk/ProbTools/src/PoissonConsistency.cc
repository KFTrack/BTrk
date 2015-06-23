//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: PoissonConsistency.cc 458 2010-01-15 11:37:35Z stroili $
//
// Description:
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	Bob Jacobsen, Ed Iskander
//      Paul Harrison:              Added sign of consistency
//      Yury Kolomensky             Corrected calculation of significance level
//      Alexandre Telnov, July 2007: 
//         fix bugs in setting the underflow status; add logLikelihood 
//
// Copyright Information:
//	Copyright (C) 1996	Lawrence Berkeley Laboratory
//      2007  Princeton University
//
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"

//-----------------------
// This Class's Header --
//-----------------------
#include "ProbTools/PoissonConsistency.hh"

//-------------
// C Headers --
//-------------
extern "C" {
#include <assert.h>
#include <math.h> 
#include <float.h> 
#include <stdlib.h>
}

//---------------
// C++ Headers --
//---------------
#include <iostream>

//-------------------------------
// Collaborating Class Headers --
//-------------------------------
#include "ProbTools/NumRecipes.hh"
#include "ProbTools/GaussConsistency.hh"

//-----------------------------------------------------------------------
// Local Macros, Typedefs, Structures, Unions and Forward Declarations --
//-----------------------------------------------------------------------
//static double EPS_CONVERGE = 0.01;

//		----------------------------------------
// 		-- Public Function Member Definitions --
//		----------------------------------------


PoissonConsistency::PoissonConsistency(int nObs, double mu) 
  : _nObs(nObs), _mu(mu)  
{ 
  if ( _mu > 0 && _nObs >= 0 ) { // idiot check
    if (nObs < mu) {
      _sign = Consistency::left;
    }
    if (nObs > mu) {
      _sign = Consistency::right;
    }
    calc();
  } else {
    _value = _likelihood = 0; _logLikelihood = -999.;
    setStatus(Consistency::underFlow);
  }
}

//		-------------------------------------------
// 		-- Protected Function Member Definitions --
//		-------------------------------------------

void 
PoissonConsistency::calc() 
{
  double logMu = log( _mu ) ;

  // get the likelihood - P(_nObs;mu)
  _logLikelihood = getLogLikelihood( _nObs, _mu, logMu ) ;
  _likelihood = exp(_logLikelihood);
  if (_logLikelihood < 0.71*DBL_MIN_EXP) setStatus(Consistency::underFlow); 

  if ( _likelihood == 0 ) {
    // looks like we are too far from peak
    _value = 0 ;
    return ; 
  }

  // the number can be adjusted
  if ( fabs( _nObs - _mu ) < 50 ) {

    calcDirectSum() ;

  } else if ( _mu > 200 ) {

    // gaussian is good for very big numbers only
    calcGauss() ;

  } else {

    // Gauss approximation is NOT good enough
    calcTails( logMu );

  }
   
  setStatus(Consistency::OK);
}


double
PoissonConsistency::getLikelihood( int n, double mu, double logMu ) const
{
  return exp(-mu+n*logMu-NumRecipes::gammln(n+1));
}

double
PoissonConsistency::getLogLikelihood( int n, double mu, double logMu ) const
{
  return -mu+n*logMu-NumRecipes::gammln(n+1);
}

void 
PoissonConsistency::calcGauss()
{  

  GaussConsistency cons(_nObs-_mu,sqrt(_mu));
  _value = cons.significanceLevel();

}
  
  
void
PoissonConsistency::calcDirectSum()
{

  // simply count all P(i) which > _likelihood
  // the algorithm is rather simple, but may be slow for big numbers,
  // so I limit it only for small differences between _nObs and _mu
  
  double P_current = _likelihood ;
  double significance = 1. ;
  int n = _nObs ;
  if ( _nObs <= _mu ) {
    double P_next = P_current * _mu / (n+1) ;
    while ( P_next > _likelihood ) {
      significance -= P_next ;
      P_current = P_next ;
      n ++ ;
      P_next = P_current * _mu / (n+1) ;
    }
  } else {                              // no, need go down
    while ( n > 0 ) {                  // at max to zero
      double P_next = P_current / _mu * n ;     // get P(I-1)
      if ( P_next > _likelihood ) {
	significance -= P_next ;
	P_current = P_next ;
	n -- ;
      } else {
	break ;                        // finish at first which is <= P(N)
      }
    }
  }
  
  _value = significance ;
}


void 
PoissonConsistency::calcTails( double logMu )
{
  
  // find the bin on the other side of peak - we need the bin which has P(i) > P(_nObs)
  int otherN = int(2*_mu - _nObs) ;
  if ( otherN < 0 ) otherN = 0 ; 
  
  // Now wander around to find exact position
  double otherL = getLikelihood( otherN, _mu, logMu ) ;
  
  int direction = 0 ;
  if ( otherN > _mu ) {
    if ( otherL < _likelihood ) {
      direction = -1 ;
    } else {
      direction = +1 ;
    } 
  } else { 
    if ( otherL < _likelihood ) {
      direction = +1 ;
    } else {
      direction = -1 ;
    } 
  }
  
  
  // we need to find such N that P(N) > _ likelihood and 
  // ( P(N+1) <= _likelihood or P(N-1) <= _likelihood )
  while ( true ) {
    int nextN = otherN + direction ;
    if ( nextN < 0 ) {
      // this can happen only if we are at 0 already but need to go 
      // even futher - we can't - stop here
      break ;
    }
    double nextL = direction > 0 ? otherL * _mu / (otherN+1) : otherL / _mu * otherN ;
    
    if ( otherL <= _likelihood && nextL > _likelihood ) {
      // next bin is OK
      otherN = nextN ;
      break ;
    } else if ( otherL > _likelihood && nextL <= _likelihood ) {
      // this bin is OK
      break ;
    } else {
      // nope, continue then
      otherN = nextN ;
      otherL = nextL ;
    }
  } 
  
  // now we should be at the good bin
  if ( otherN > _nObs ) {
    double thisP = NumRecipes::gammq ( _nObs+1, _mu ) ;
    double otherP = NumRecipes::gammp ( otherN+1, _mu ) ;
    _value = otherP + thisP ;
  } else {
    double thisP = NumRecipes::gammp ( _nObs, _mu ) ;
    double otherP = 0 ;
    if ( otherN > 0 ) otherP = NumRecipes::gammq ( otherN, _mu ) ;
    _value = otherP + thisP ;
  } 
  
}

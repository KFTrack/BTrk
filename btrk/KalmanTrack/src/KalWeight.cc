// ------------------------------------------------------------------------------
// File and Version Information:
//      $Id: KalWeight.cc,v 1.11 2008/03/17 13:14:30 brownd Exp $
//
//  Description: KalWeight
//
// Copyright Information:
//	Copyright (C) 1997	Lawrence Berkeley Laboratory
//
//  Authors: Dave Brown, 8/22/97
//------------------------------------------------------------------------------
#include "BaBar/BaBar.hh"

#include "KalmanTrack/KalWeight.hh"
#include "KalmanTrack/KalParams.hh"
#include <iostream>
using std::endl;
using std::ostream;


KalWeight::KalWeight() : _status(1)
{;}

KalWeight::KalWeight(int ndim) : 
  _wvector(ndim,0),_wcov(ndim,0),_status(1)
{;}

KalWeight::KalWeight(const HepVector& pvec, const HepSymMatrix& pcov) : 
  _wvector(pvec),_wcov(pcov),_status(0)
{;}

KalWeight::KalWeight(const KalWeight& other) :
  _wvector(other._wvector),_wcov(other._wcov),_status(other._status)
{
#ifdef KALDEBUG	
  double wdet = _wcov.determinant();
  if(isnan(wdet)){
	std::cout << "bad copy input, status = " << _status << endl;
	this->print(std::cout);
  }
#endif
}

KalWeight::KalWeight(const KalParams& params) : _status(params.status()){
// failsafe construction, to insure at least dimensional correctness
  if(matrixOK())
    _wcov = params.covarianceMatrix().inverse(_status);
  else
    _wcov = HepSymMatrix(params.parameterVector().num_row(),1);
  if(matrixOK())
    _wvector = _wcov*params.parameterVector();
  else
    _wvector = HepVector(params.parameterVector().num_row(),0);

#ifdef KALDEBUG
  double wdet = _wcov.determinant();
  double pdet = params.covarianceMatrix().determinant();
  if(isnan(wdet) || isnan(pdet)){
	std::cout << "bad inversion, status = " << _status << endl;
	params.print(std::cout);
	this->print(std::cout);
  }	
#endif
}

KalWeight& 
KalWeight::operator = (const KalParams& params) {
// failsafe construction, to insure at least dimensional correctness
  _status = params.status();
  if(matrixOK())
    _wcov = params.covarianceMatrix().inverse(_status);
  else
    _wcov = HepSymMatrix(params.parameterVector().num_row(),1);
  if(matrixOK())
    _wvector = _wcov*params.parameterVector();
  else
    _wvector = HepVector(params.parameterVector().num_row(),0);

#ifdef KALDEBUG
  double wdet = _wcov.determinant();
  double pdet = params.covarianceMatrix().determinant();
  if(isnan(wdet) || isnan(pdet)){
	std::cout << "bad inversion, status = " << _status << endl;
	params.print(std::cout);
	this->print(std::cout);
  }
#endif
  return *this;
}
   

KalWeight& KalWeight::operator = (const KalWeight& other){
  if (&other != this) {
    _wvector = other._wvector;
    _wcov = other._wcov;
    _status = other._status;
  }

#ifdef KALDEBUG
  double wdet = _wcov.determinant();
  if(isnan(wdet)){
  	std::cout << "bad input, status = " << _status << endl;
 	this->print(std::cout);
  }
#endif

  return *this;
}

KalWeight&
KalWeight::operator += (const KalWeight& other) {
  if(matrixOK() && other.matrixOK()){
    _wvector += other._wvector;
    _wcov += other._wcov;
  } else if(other.matrixOK()){
    *this = other;
  }
#ifdef KALDEBUG
  double wdet = _wcov.determinant();
  if(isnan(wdet)){
	std::cout << "bad addition, status = " << _status << endl;
	this->print(std::cout);
  }
#endif

  return *this;
}

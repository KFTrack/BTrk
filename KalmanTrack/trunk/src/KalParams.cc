// ------------------------------------------------------------------------------
// File and Version Information:
//      $Id: KalParams.cc,v 1.22 2008/03/17 13:14:30 brownd Exp $
//
//  Description: KalParams
//
// Copyright Information:
//	Copyright (C) 1997	Lawrence Berkeley Laboratory
//
//  Authors: Dave Brown, 8/22/97
//------------------------------------------------------------------------------
#include "BaBar/BaBar.hh"

#include "KalmanTrack/KalParams.hh"
#include "KalmanTrack/KalWeight.hh"
#include "TrkBase/TrkParams.hh"
#include <iostream>
using std::endl;
using std::ostream;

KalParams::KalParams() : _status(1)
{;}

KalParams::KalParams(int ndim) : _pvector(ndim,0),_pcov(ndim,0),
  _status(1)
{;}

KalParams::KalParams(const HepVector& pvec, const HepSymMatrix& pcov) : 
  _pvector(pvec),_pcov(pcov),_status(0)
{;}

KalParams::KalParams(const KalParams& other) :
  _pvector(other._pvector),_pcov(other._pcov),_status(other._status)
{;}

KalParams::KalParams(const KalWeight& weight) : _status(weight.status()){
// failsafe construction, to insure at least dimensional correctness
  if(matrixOK())
    _pcov = weight.weightMatrix().inverse(_status);
  else
    _pcov = HepSymMatrix(weight.weightVector().num_row(),1);
  if(matrixOK())
    _pvector = _pcov*weight.weightVector();
  else
    _pvector = HepVector(weight.weightVector().num_row(),0);

#ifdef KALDEBUG
  double wdet = weight.weightMatrix().determinant();
  double pdet = _pcov.determinant();
  if(isnan(wdet) || isnan(pdet)){
	std::cout << "bad inversion, status = " << _status << endl;
	weight.print(std::cout);
	this->print(std::cout);
  }
#endif

}

KalParams::KalParams(const TrkParams& tpar) :
  _pvector(tpar.parameter()),_pcov(tpar.covariance()),_status(0)
{
}

void
KalParams::diagonalize() {
  if(matrixOK()){
    int nrow = _pcov.num_row();
    for(int irow=0;irow<nrow;irow++)
      for(int icol=0;icol<irow;icol++)
	_pcov.fast(irow+1,icol+1) = 0.0;
  }
}

KalParams&    
KalParams::operator = (const KalWeight& weight) {
// failsafe construction, to insure at least dimensional correctness
  _status = weight.status();
  if(matrixOK())
    _pcov = weight.weightMatrix().inverse(_status);
  else
    _pcov = HepSymMatrix(weight.weightVector().num_row(),1);
  if(matrixOK())
    _pvector = _pcov*weight.weightVector();
  else
    _pvector = HepVector(weight.weightVector().num_row(),0);

#ifdef KALDEBUG
  double wdet = weight.weightMatrix().determinant();
  double pdet = _pcov.determinant();
  if(isnan(wdet) || isnan(pdet)){
	std::cout << "bad inversion, status = " << _status << endl;
	weight.print(std::cout);
	this->print(std::cout);
  }
#endif

  return *this;
}

double
KalParams::chisq(const KalParams& other,bool* tparams) const {
  double chisq(-1.0);
  if(matrixOK() && other.matrixOK()){
// compute parameter difference and covariance sum
    HepVector pdiff = _pvector - other._pvector;
    HepSymMatrix csum(_pcov + other._pcov);
    if(tparams == 0 ){
      int ierr(0);
      csum.invert(ierr);
      if(ierr == 0)
	chisq = csum.similarity(pdiff);
    } else {
// mask off parameters as specified
// We must reduce the dimensionality to be able to invert the combined
// error matrix.  First, compute a reduction matrix
      unsigned npar = pdiff.num_row();
      unsigned ndof(0);
      for(unsigned ipar=0;ipar<npar;++ipar)
	if(tparams[ipar])++ndof;
      if(ndof > 0){
	HepMatrix reduce(ndof,npar,0);
	unsigned idof(0);
	for(unsigned ipar=0;ipar<npar;ipar++){
	  if(tparams[ipar]){
	    reduce[idof][ipar] = 1.0;
	    idof++;
	  }
	}
// project out the unconstrained parameters
	HepSymMatrix rcsum = csum.similarity(reduce);
	HepVector rpdiff = reduce*pdiff;
// invert the reduced covariance and caclculate the chisq
	int ierr;
	rcsum.invert(ierr);
	if(ierr == 0)
	  chisq = rcsum.similarity(rpdiff);
      } else
// if there are no real DOFs, the chisq = 0;
	chisq = 0.0;
    }
  }
  return chisq; 
}


void
KalParams::print(ostream& os) const {
  static const char* pnames[5] = {"d0","phi0","omega","z0","tandip"};
  unsigned nparams = _pvector.num_row();
  for(unsigned iparam=0;iparam<nparams;iparam++){
    os << "Parameter " << pnames[iparam] << " = " 
       << _pvector[iparam] << " +- "
       << sqrt(_pcov[iparam][iparam]) << endl;
  }
}

// inversion function
KalParams&
KalParams::invert(const TrkSimpTraj* traj) {
  if(matrixOK()){
  // Get the TrkParams equiv of the KalParams
    TrkParams params = trackParameters();
  // Let the trajectory handle the parameter inversion
    std::vector<bool> flags(nPar(),false);
    traj->invertParams(&params, flags);
  // Put results into _pvector
    _pvector = params.parameter();
  // Invert covariance matrix, too.
    for (int ipar = 0; ipar < nPar(); ipar++) {
      for (int jpar = 0; jpar < ipar; jpar++) {
        if ( (flags[ipar] && !flags[jpar]) ||
             (!flags[ipar] && flags[jpar]) )
          _pcov.fast(ipar+1, jpar+1) *= -1.0;
      }
    }
  }
  return *this;
}

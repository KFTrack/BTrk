//--------------------------------------------------------------------------
//
// Description:
//  Code for the NeutParams neutral track parameterization class
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author(s): Justin Albert, Valery Miftahov (based on HelixParams.cc by
//                                            Dave Brown)
//
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"
#include "TrkBase/NeutParams.hh"
#include "TrkBase/TrkExchangePar.hh"
#include "ErrLogger/ErrLog.hh"

// construct from vector and covariance matrix
//------------------------------------------------------------------------
NeutParams::NeutParams(const HepVector& pvec,const HepSymMatrix& pcov) :
  TrkParams(pvec,pcov){
//------------------------------------------------------------------------

//  Make sure the dimensions of the input matrix and vector are correct

    if( pvec.num_row() != _nneutprm ||
	pcov.num_row() != _nneutprm ){
      ErrMsg(error) << 
	"NeutParams: incorrect constructor vector/matrix dimension" << endmsg;
      parameter() = HepVector(_nneutprm,0);
      covariance() = HepSymMatrix(_nneutprm,0);
    }
}

//  Construct from the fit parameters directly
//------------------------------------------------------------------------
NeutParams::NeutParams(double n_d0,double n_phi0,double n_p,double n_z0,
			 double n_tanDip,double n_s0) : 
//------------------------------------------------------------------------
  TrkParams(_nneutprm) {
    d0() = n_d0;
    phi0() = n_phi0;
    p() = n_p;
    z0() = n_z0;
    s0() = n_s0;
    tanDip() = n_tanDip;
}


//  Copy constructor
//------------------------------------------------------------------------
NeutParams::NeutParams(const NeutParams& old) :
//------------------------------------------------------------------------
  TrkParams(old){
}
//------------------------------------------------------------------------
double 
NeutParams::sinPhi0() const {
//------------------------------------------------------------------------
  return sin(parameter()[_phi0]); 
}

//------------------------------------------------------------------------
double 
NeutParams::cosPhi0() const {
//------------------------------------------------------------------------
  return cos(parameter()[_phi0]); 
}

//------------------------------------------------------------------------
double 
NeutParams::arcRatio() const {
//------------------------------------------------------------------------
  return sqrt(1. + parameter()[_tanDip] * parameter()[_tanDip]); 
}

//------------------------------------------------------------------------
NeutParams::~NeutParams()
{}
//------------------------------------------------------------------------


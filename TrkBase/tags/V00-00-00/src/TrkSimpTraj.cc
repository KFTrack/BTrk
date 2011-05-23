//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkSimpTraj.cc,v 1.31 2004/11/29 22:13:33 brownd Exp $
//
// Description:
//     
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author(s): Steve Schaffner
//
//------------------------------------------------------------------------

#include "BaBar/BaBar.hh"
#include "BaBar/Constants.hh"
#include "BbrGeom/BbrAngle.hh"
#include "TrkBase/TrkSimpTraj.hh"
#include "TrkBase/TrkParams.hh"
#include "CLHEP/Matrix/SymMatrix.h"
#include "TrkBase/TrkPocaXY.hh"
#include "ErrLogger/ErrLog.hh"
using std::ostream;

// statics
HepPoint TrkSimpTraj::_theOrigin(0.0,0.0,0.0);
//Constructors
//----------------------------------------------------------------------------
TrkSimpTraj::TrkSimpTraj(const HepVector& params, const HepSymMatrix& cov, 
                         const double lowlim,const double hilim,
			 const HepPoint& refpoint) :
//----------------------------------------------------------------------------
  TrkDifTraj(lowlim, hilim),_dtparams(params, cov),_refpoint(refpoint)
{;}

TrkSimpTraj::TrkSimpTraj(const TrkParams& params,
                         const double lowlim,const double hilim,
			 const HepPoint& refpoint) :
//----------------------------------------------------------------------------
  TrkDifTraj(lowlim, hilim),_dtparams(params),_refpoint(refpoint)
{;}

//----------------------------------------------------------------------------
TrkSimpTraj::TrkSimpTraj(const TrkSimpTraj& other) :
//----------------------------------------------------------------------------
  TrkDifTraj(other.lowRange(),other.hiRange()),
  _dtparams(other._dtparams),
  _refpoint(other._refpoint)
{}

//----------------------------------------------------------------------------
TrkSimpTraj::~TrkSimpTraj()
//----------------------------------------------------------------------------
{ }

//----------------------------------------------------------------------------
const TrkSimpTraj*
TrkSimpTraj::localTrajectory(double fltLen, double& localFlt) const {
//----------------------------------------------------------------------------
  localFlt = fltLen;
  return this;
}

//
//  The following is the only useful function which can be implemented at this level
//
//----------------------------------------------------------------------------
void
TrkSimpTraj::changePoint(const HepPoint& newpoint,double& fltlen) {
//----------------------------------------------------------------------------
  if(newpoint != _refpoint){
//  find POCA to the new point
    TrkPocaXY endpoca(*this,fltlen,newpoint);
    if(endpoca.status().failure()){
      ErrMsg(error) << "poca failure changing reference point" << endmsg;
      return;
    } else {
// find POCA to the old point: this is used to reset the flightlength
      TrkPocaXY oldpoca(*this,0.0,_refpoint);
      double oldfltlen = oldpoca.flt1();
// update flight length
      fltlen = endpoca.flt1();
//  Get the translation function
      TranslateParams pfunc = paramFunction();
//  Use it on the SimpTraj parameters
      pfunc(_refpoint,newpoint,
	    parameters()->parameter(),parameters()->covariance(),
	    _dtparams.parameter(),_dtparams.covariance(),
	    fltlen);
      _refpoint = newpoint;
// update the flight range to correspond to the same range in space as before
      double newrange[2];
      newrange[0] = lowRange() - fltlen + oldfltlen;
      newrange[1] = hiRange() - fltlen + oldfltlen;
      setFlightRange(newrange);
    }
  }
  return;
}
//----------------------------------------------------------------------------
void 
TrkSimpTraj::printAll(ostream& os) const {
//----------------------------------------------------------------------------
    os << "Simple ";
    Trajectory::printAll(os);
    os << "SimpTraj parameter vector = "
       << _dtparams.parameter();
    os << "  and covariance matrix =  " 
       << _dtparams.covariance();
  }
//----------------------------------------------------------------------------
void 
TrkSimpTraj::print(ostream& os) const {
//----------------------------------------------------------------------------
  os << "Simple ";
  Trajectory::print(os);
}
// Inversion.  This changes the flightrange from A->B to (-B)->(-A), as well
// adjusting the parameters so that the same arc in space is described (just
// going the opposite direction).  This is now implemented fully generically.
TrkSimpTraj&
TrkSimpTraj::invert()
{
  // Invert parameters
  std::vector<bool> flags(parameters()->nPar(),false);
  invertParams(parameters(), flags);
  // loop over parameters and invert covariance matrix
  for(int iparam=0;iparam<parameters()->nPar();iparam++){
    bool iinvert = flags[iparam];
// do covariance cross-terms too
    for(int jparam=iparam+1;jparam<parameters()->nPar();jparam++){
      bool jinvert = flags[jparam];
      if( (iinvert && !jinvert) || (!iinvert && jinvert) ) {
// cross-terms change sign
        parameters()->covariance()[iparam][jparam] *= -1.0;
      }
    }
  }
// invert the flightlength
  double range[2];
  range[0] = -hiRange();
  range[1] = -lowRange();
  setFlightRange(range);
// done
  return *this;
}


bool 
TrkSimpTraj::operator==(const TrkSimpTraj& x) const
{
     if (lowRange()!=x.lowRange() || hiRange()!=x.hiRange()) return false;
     const HepVector &m=_dtparams.parameter(); 
     unsigned int mp=m.num_row();
     const HepVector &n=x._dtparams.parameter(); 
     unsigned int np=n.num_row();
     if (np!=mp) return false;
     for(unsigned i=0;i<np;++i){
	if(m[i] != n[i]) return false;
     }
     return _refpoint==x._refpoint;
}

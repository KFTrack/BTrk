//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkExchangePar.hh,v 1.17 2004/08/06 06:31:40 bartoldu Exp $
//
// Description: Class to pass around a minimal set of track parameters from 
//   one class to another.  It has no functionality.
//     
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Authors: Steve Schaffner
//
// Revision History:
//	20020417  M. Kelsey -- Add print(), printAll(), and operator<<
//------------------------------------------------------------------------
#ifndef TRKEXCHANGEPAR_HH
#define TRKEXCHANGEPAR_HH
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/SymMatrix.h"

#include <iosfwd>

// Class interface //
class TrkExchangePar {
public:
  enum {nParam = 5};
  enum {ex_d0, ex_phi0, ex_omega, ex_z0, ex_tanDip};

  TrkExchangePar(const HepVector&, const HepMatrix&);
  TrkExchangePar(const HepVector&, const HepSymMatrix&);
  TrkExchangePar(const HepVector&);
  TrkExchangePar(double d0In, double phi0In, double omegaIn, 
		 double z0In, double tanDipIn);
  ~TrkExchangePar();
  
  double d0() const                              {return paramVec[ex_d0];}
  double phi0() const                            {return paramVec[ex_phi0];}
  double omega() const                           {return paramVec[ex_omega];}
  double z0() const                              {return paramVec[ex_z0];}
  double tanDip() const                          {return paramVec[ex_tanDip];}

  const HepVector& params() const                {return paramVec;}
  HepVector& params()                            {return paramVec;}
  const HepSymMatrix& covariance() const            {return paramErr;}
  HepSymMatrix& covariance()                        {return paramErr;}

  void setD0(double in)                          {paramVec[ex_d0] = in;} 
  void setPhi0(double in)                        {paramVec[ex_phi0] = in;} 
  void setOmega(double in)                       {paramVec[ex_omega] = in;} 
  void setZ0(double in)                          {paramVec[ex_z0] = in;} 
  void setTanDip(double in)                      {paramVec[ex_tanDip] = in;} 
  void setError(const HepSymMatrix& in)          {paramErr = in;}

  void print(std::ostream& o) const;		// Print parameters on one line
  void printAll(std::ostream& o) const;	// Print parameters and error matrix

private:	
  HepVector paramVec;
  HepSymMatrix paramErr;

};

// Output operator, useful for debugging
std::ostream& operator<<(std::ostream& o, const TrkExchangePar& helix);

#endif

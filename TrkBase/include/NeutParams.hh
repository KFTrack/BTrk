//  Justin Albert and Valery Miftahov 10/27/97

#ifndef NEUTPARAMS_HH
#define NEUTPARAMS_HH

#include <stdio.h>
#include <iostream>
#include <math.h>

#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/SymMatrix.h"

#include "TrkBase/TrkParams.hh"

class NeutParams : public TrkParams {

friend class NeutTraj;

public:

//  Constructors and destructor
//
  NeutParams(const HepVector&,const HepSymMatrix&);
  NeutParams(double,double,double,double,double,double);
  NeutParams(const NeutParams& old);
  ~NeutParams();

//  Define the parameter meaning and order by an enum
  enum {_d0, _phi0, _p, _z0, _tanDip, _s0,_nneutprm = 6};

//  access
//
  double& d0()                               { return parameter()[_d0]; }
  double& phi0()                             { return parameter()[_phi0]; }
  double& p()                                { return parameter()[_p]; }
  double& z0()                               { return parameter()[_z0]; }
  double& tanDip()                           { return parameter()[_tanDip]; }
  double& s0()                               { return parameter()[_s0]; }

  double d0() const                          { return parameter()[_d0]; }
  double phi0() const                        { return parameter()[_phi0]; }
  double p() const                           { return parameter()[_p]; }
  double z0() const                          { return parameter()[_z0]; }
  double tanDip() const                      { return parameter()[_tanDip]; }
  double s0()  const                         { return parameter()[_s0]; }

  double sinPhi0() const;                    
  double cosPhi0() const;
  double arcRatio() const;  // = fltLen / 2-d arclen
private:
};
#endif


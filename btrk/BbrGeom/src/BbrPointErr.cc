//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: BbrPointErr.cc 491 2010-01-13 16:59:16Z stroili $
//
// Description:
//      Class BbrPointErr
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Forest Rouse            February 1996
//      Victoria Novotny        August   1996
//
// Copyright Information:
//      Copyright (C) 1996      U.C. Davis
//
//------------------------------------------------------------------------
// File BbrPointErr.cc
// Source file for class BbrPointErr
//
// For advice, input, or any questions then please contact either
// Bob Jacobsen <Bob_Jacobsen@lbl.gov> or
// Forest Rouse <rouse@ucdhep.ucdavis.edu>
//
// =====================================================================
// Name           Change description
// Date
// Version
// =====================================================================
#include "BaBar/BaBar.hh"

#include <stdio.h>
#include <float.h>
#include "difAlgebra/DifNumber.hh"
#include "difAlgebra/DifArray.hh"
#include "BbrGeom/BbrPointErr.hh"
using std::istream;
using std::ostream;

BbrPointErr operator + (const BbrPointErr& v, const BbrVectorErr& w){
    BbrPointErr pe(HepPoint(v.x()+w.x(),v.y()+w.y(),v.z()+w.z()),
  	   (v.covMatrix()+w.covMatrix()));
    return pe;
  }

BbrPointErr operator + (const BbrVectorErr& w, const BbrPointErr& v){
    BbrPointErr pe(HepPoint(v.x()+w.x(),v.y()+w.y(),v.z()+w.z()), 
                   (v.covMatrix()+w.covMatrix()));
    return pe;
}

BbrPointErr operator - (const BbrPointErr& v, const BbrVectorErr& w){
    BbrPointErr pe(HepPoint(v.x()-w.x(),v.y()-w.y(),v.z()-w.z()), 
                   (v.covMatrix()+w.covMatrix()));
    return pe;
}

BbrVectorErr operator - (const BbrPointErr& v, const BbrPointErr& w){
    BbrVectorErr ve(Hep3Vector(v.x()-w.x(),v.y()-w.y(),v.z()-w.z()),
		    (v.covMatrix()+w.covMatrix()));
    return ve;
}

ostream & operator<<(ostream & stream, const BbrPointErr & verr) {
  stream << (const HepPoint&)verr
	 << ", " << verr.covMatrix();
  
  return stream;
}

istream & operator>>(istream & stream, BbrPointErr & verr) {
  BbrError mat(verr.SIZE);
  stream >> (HepPoint&)verr >> mat;
  verr.setCovMatrix(mat);
  
  return stream;
}

BbrError BbrPointErr::covRTPMatrix() const{
  // protect against 0's
  double xv = x()==0 ?  FLT_MIN : x();
  double yv = y()==0 ?  FLT_MIN : y();
  double zv = z()==0 ?  FLT_MIN : z();
  DifNumber xDF(xv,X+1,3), yDF(yv,Y+1,3), zDF(zv,Z+1,3);
  DifArray pars(3,NUM_PCOORDINATES);
  pars[Rho]   =  sqrt(xDF*xDF + yDF*yDF + zDF*zDF);
  pars[Phi]   = atan2(yDF,xDF);
  pars[Theta] = acos(zDF/pars[Rho]);
  return covMatrix().similarity(pars.jacobian());
}

BbrError BbrPointErr::covRZPMatrix() const{
  // protect against 0's
  double xv = x()==0 ?  FLT_MIN : x();
  double yv = y()==0 ?  FLT_MIN : y();
  double zv = z()==0 ?  FLT_MIN : z();
  DifNumber xDF(xv,X+1,3), yDF(yv,Y+1,3), zDF(zv,Z+1,3);
  DifArray pars(3,NUM_CCOORDINATES);
  pars[C_Rho]   =  sqrt(xDF*xDF + yDF*yDF );
  pars[C_Phi]   = atan2(yDF,xDF);
  pars[C_Zeta] =  zDF;
  return covMatrix().similarity(pars.jacobian());
}


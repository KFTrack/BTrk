//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: BbrVectorErr.cc 491 2010-01-13 16:59:16Z stroili $
//
// Description:
//      Class BbrVectorErr
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
// 2002/01/19, Sasha Telnov: added operator * (scaling by a real number) 
//
//------------------------------------------------------------------------
// File BbrVectorErr.cc
// Source file for class BbrVectorErr
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

#include "BTrk/BaBar/BaBar.hh"
#include <float.h>
#include <iostream>
#include "BTrk/BbrGeom/BbrVectorErr.hh"
#include "BTrk/difAlgebra/DifNumber.hh"
#include "BTrk/difAlgebra/DifArray.hh"
using std::istream;
using std::ostream;
using namespace CLHEP;

//void BbrVectorErr::printOn(ostream& out) const
//{
//    out << "Bbr Vector and Covariance Matrix:" << endl;
//    BbrVector::printOn(out);
//    out << endl;
//    out << _covMatrix;
//}
 
BbrError BbrVectorErr::covRTPMatrix() const{
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

BbrError BbrVectorErr::covRZPMatrix() const{
  // protect against 0's
  double xv = x()==0 ?  FLT_MIN : x();
  double yv = y()==0 ?  FLT_MIN : y();
  double zv = z()==0 ?  FLT_MIN : z();
  DifNumber xDF(xv,X+1,3), yDF(yv,Y+1,3), zDF(zv,Z+1,3);
  DifArray pars(3,NUM_CCOORDINATES);
  pars[C_Rho]   =  sqrt(xDF*xDF + yDF*yDF);
  pars[C_Phi]   = atan2(yDF,xDF);
  pars[C_Zeta] = zDF;
  return covMatrix().similarity(pars.jacobian());
}

double BbrVectorErr::determineChisq(const Hep3Vector& refVector) const
{
   HepVector temp(NUM_COORDINATES, 0);
   temp[0] = refVector.x()-this->x();
   temp[1] = refVector.y()-this->y();
   temp[2] = refVector.z()-this->z();
   return _covMatrix.determineChisq(temp);
}


BbrVectorErr operator + (const BbrVectorErr& v, const BbrVectorErr& w){
    BbrVectorErr ve(Hep3Vector(v.x()+w.x(),v.y()+w.y(),v.z()+w.z()),
		    (v.covMatrix()+w.covMatrix()));
    return ve;
}

BbrVectorErr operator - (const BbrVectorErr& v, const BbrVectorErr& w){
    BbrVectorErr ve(Hep3Vector(v.x()-w.x(),v.y()-w.y(),v.z()-w.z()), 
                    (v.covMatrix()+w.covMatrix()));
    return ve;
}

// Added by Sasha Telnov
BbrVectorErr operator * (const BbrVectorErr & p, HepDouble a) {
  BbrVectorErr ve(Hep3Vector(a*p.x(), a*p.y(), a*p.z()),
                  (p.covMatrix()*a*a));
  return ve;
}

BbrVectorErr operator * (HepDouble a, const BbrVectorErr & p) {
  BbrVectorErr ve(Hep3Vector(a*p.x(), a*p.y(), a*p.z()),
                  (p.covMatrix()*a*a));
  return ve;
}


ostream & operator<<(ostream & stream, const BbrVectorErr & verr) {
  stream << (const Hep3Vector&)verr
	 << ", " << verr.covMatrix();
  
  return stream;
}

istream & operator>>(istream & stream, BbrVectorErr & verr) {
  BbrError mat(verr.SIZE);
  stream >> (Hep3Vector&)verr >> mat;
  verr.setCovMatrix(mat);
  
  return stream;
}


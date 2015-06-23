//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DifRotation.cc 501 2010-01-14 12:46:50Z stroili $
//
// Description:
//	Class Implementation for |DifRotation| 
//      Rotates things
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	A. Snyder
//
// Copyright Information:
//	Copyright (C) 1996	SLAC
//
//------------------------------------------------------------------------

#include "BaBar/BaBar.hh"

#include "difAlgebra/DifRotation.hh"

DifRotation::DifRotation()
  :_xnew(1.0,0.0,0.0),_ynew(0.0,1.0,0.0),_znew(0.0,0.0,1.0)
{
}

DifRotation::DifRotation
(const DifNumber& alpha,const DifNumber& beta,const DifNumber& gamma)
{
  
  //active transformation  - move the vector (reverse sines)
  DifNumber calpha=cos(alpha); DifNumber salpha=-sin(alpha);
  DifNumber cbeta=cos(beta); DifNumber sbeta=-sin(beta);
  DifNumber cgamma=cos(gamma); DifNumber sgamma=-sin(gamma);

 _xnew.x=cbeta*calpha*cgamma-salpha*sgamma;
 _xnew.y=cbeta*salpha*cgamma+calpha*sgamma;
 _xnew.z=-sbeta*cgamma;

 _ynew.x=-cbeta*calpha*sgamma-salpha*cgamma;
 _ynew.y=-cbeta*salpha*sgamma+calpha*cgamma;
 _ynew.z=sbeta*sgamma;

 _znew.x=sbeta*calpha;
 _znew.y=sbeta*salpha;
 _znew.z=cbeta;
}
DifRotation::DifRotation
(double alpha,double beta,double gamma)
{
  //active transformation  - move the vector (reverse sines)
  DifNumber calpha(cos(alpha)); DifNumber salpha(-sin(alpha));
  DifNumber cbeta(cos(beta)); DifNumber sbeta(-sin(beta));
  DifNumber cgamma(cos(gamma)); DifNumber sgamma(-sin(gamma));

 _xnew.x=cbeta*calpha*cgamma-salpha*sgamma;
 _xnew.y=cbeta*salpha*cgamma+calpha*sgamma;
 _xnew.z=-sbeta*cgamma;

 _ynew.x=-cbeta*calpha*sgamma-salpha*cgamma;
 _ynew.y=-cbeta*salpha*sgamma+calpha*cgamma;
 _ynew.z=sbeta*sgamma;

 _znew.x=sbeta*calpha;
 _znew.y=sbeta*salpha;
 _znew.z=cbeta;
}



DifRotation::DifRotation
(const DifVector& xp,const DifVector& yp,const DifVector &zp)
  :_xnew(xp),_ynew(yp),_znew(zp)
{}				

DifRotation::DifRotation
(const DifVector& xp,const DifVector& yp) 
  :_xnew(xp),_ynew(yp),_znew()
{
  _znew=cross(_xnew,_ynew);
}

int DifRotation::fail()const {
  return 0;
}

void DifRotation::rotate(DifVector &v)const {
  DifNumber xcomp=xnew()*v;
  DifNumber ycomp=ynew()*v;
  DifNumber zcomp=znew()*v;
  v.x=xcomp; v.y=ycomp; v.z=zcomp;
}

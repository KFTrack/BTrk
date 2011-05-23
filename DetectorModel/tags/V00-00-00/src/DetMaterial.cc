// ------------------------------------------------------------------------------
// File and Version Information:
//      $Id: DetMaterial.cc,v 1.35 2008/09/08 22:04:29 brownd Exp $
//
//  Description:
//  This class defines the standard values and formulas for particles interacting
//  with matter.  It is used in the DetectorModel description heirarchy.
//
//
// Copyright Information:
//	Copyright (C) 1996	Lawrence Berkeley Laboratory
//
//  Authors: Dave Brown, 11/21/96
//           Matthias Steinke 04/09/99
//------------------------------------------------------------------------------
#include "BaBar/BaBar.hh"
#include "DetectorModel/DetMaterial.hh"
#include <iostream>
#include <cfloat>
#include <string>
#include <vector>

#include "CLHEP/Units/PhysicalConstants.h"
using namespace CLHEP;

using std::endl;
using std::ostream;

//
// setup the static members
//
//  The multiple scattering constant has an effective thickness
//  for the logrithmic term in the scattering formula corresponding to
//  the average material interior to the tracking volume at dip=~30 degrees
//

const double DetMaterial::_msmom = 12.2*MeV;
const double DetMaterial::_dgev = 0.153536*MeV*cm*cm;
const double bg2lim = 0.0169 , taulim = 8.4146e-3 ;
const double twoln10 = 2.0*log(10.);

const double betapower = 1.667; // most recent PDG gives beta^-5/3 as dE/dx
const int maxnstep = 10; // maximum number of steps through a single material
//
//  Constructor
//

DetMaterial::DetMaterial(const char* detMatName, const DetMtrProp* detMtrProp):
  _name(detMatName),
  _za(detMtrProp->getZ()/detMtrProp->getA()),
  _zeff(detMtrProp->getZ()),
  _aeff(detMtrProp->getA()),
  _radthick(detMtrProp->getRadLength()/cm/cm),
  _intLength(detMtrProp->getIntLength()/detMtrProp->getDensity()),
  _meanion(2.*log(detMtrProp->getMeanExciEnergy()*1.0e6)),
  _eexc(detMtrProp->getMeanExciEnergy()),
  _x0(detMtrProp->getX0density()),
  _x1(detMtrProp->getX1density()),
  _delta0(detMtrProp->getDEdxFactor()),
  _afactor(detMtrProp->getAdensity()),
  _mpower(detMtrProp->getMdensity()),
  _bigc(detMtrProp->getCdensity()),
  _density(detMtrProp->getDensity()/cm/cm/cm),
  _noem(detMtrProp->getNumberOfElements()),
  _taul(detMtrProp->getTaul())
{
  _shellCorrectionVector = 
    new std::vector< double >(detMtrProp->getShellCorrectionVector());
  _vecNbOfAtomsPerVolume = 
    new std::vector< double >(detMtrProp->getVecNbOfAtomsPerVolume());
  _vecTau0 = 
    new std::vector< double >(detMtrProp->getVecTau0());
  _vecAlow = new std::vector< double >(detMtrProp->getVecAlow());
  _vecBlow = new std::vector< double >(detMtrProp->getVecBlow());
  _vecClow = new std::vector< double >(detMtrProp->getVecClow());
  _vecZ = new std::vector< double >(detMtrProp->getVecZ());
}

DetMaterial::~DetMaterial()
{
  delete _shellCorrectionVector;
  delete _vecNbOfAtomsPerVolume;
  delete _vecTau0;
  delete _vecAlow;
  delete _vecBlow;
  delete _vecClow;
  delete _vecZ;
}

//
//  Multiple scattering static function (from PDG 96)
//
double
DetMaterial::scatterAngleRMS(double mom,
					 double pathlen,double mass) const {
  if(mom>0.0){
    double radfrac = _density*fabs(pathlen)/_radthick;
//  The logrithmic term is non-local, and so shouldn't be used.  An effective
//  term is included in _msmom, see above
    return _msmom*sqrt(radfrac)/(mom*particleBeta(mom,mass));
  } else
    return 1.0; // 'infinite' scattering
}


/********************** New Routines **************************/ 


double
DetMaterial::dEdx(double mom,dedxtype type,double mass) const {
  if(mom>0.0){
    double Eexc2 = _eexc*_eexc ;
  
  // New energy loss implementation
  
    double Tmax,gamma2,beta2,bg2,rcut,delta,x,sh,dedx ;
    double beta  = particleBeta(mom,mass) ;
    double gamma = particleGamma(mom,mass) ;
    double tau = gamma-1. ;
  
  // high energy part , Bethe-Bloch formula 
  
    beta2 = beta*beta ;
    gamma2 = gamma*gamma ;
    bg2 = beta2*gamma2 ;
  
    double RateMass = Pdt::mass(PdtPid::electron) / mass;
  
    Tmax = 2.*electron_mass_c2*bg2
      /(1.+2.*gamma*RateMass+RateMass*RateMass) ;

    const double cutOffEnergy = 50.;
    rcut =  ( cutOffEnergy< Tmax) ? cutOffEnergy/Tmax : 1;

    dedx = log(2.*electron_mass_c2*bg2*Tmax/Eexc2);
    if(type == loss)
      dedx -= 2.*beta2;
    else
      dedx += log(rcut)-(1.+rcut)*beta2;
  
// density correction 
    x = log(bg2)/twoln10 ;
    if ( x < _x0 ) {
      if(_delta0 > 0) {
	delta = _delta0*pow(10.0,2*(x-_x0));
      }
      else {
	delta = 0.;
      }
    } else {
      delta = twoln10*x - _bigc;
      if ( x < _x1 )
	delta += _afactor * pow((_x1 - x), _mpower);
    } 
  
  // shell correction          
    if ( bg2 > bg2lim ) {
      sh = 0. ;      
      x = 1. ;
      for (int k=0; k<=2; k++) {
	x *= bg2 ;
	sh += (*_shellCorrectionVector)[k]/x;
      }
    }
    else {
      sh = 0. ;      
      x = 1. ;
      for (int k=0; k<2; k++) {
	x *= bg2lim ;
	sh += (*_shellCorrectionVector)[k]/x;
      }
      sh *= log(tau/_taul)/log(taulim/_taul);
    }
    dedx -= delta + sh ;
    dedx *= -_dgev*_density*_za / beta2 ;
    return dedx;
  } else
    return 0.0;
}



double 
DetMaterial::energyLoss(double mom, double pathlen,double mass) const {
// make sure we take positive lengths!
  pathlen = fabs(pathlen);
  double dedx = dEdx(mom,loss,mass);
// see how far I can step, within tolerance, given this energy loss
  double maxstep = maxStepdEdx(mom,mass,dedx);
// if this is larger than my path, I'm done
  if(maxstep>pathlen){
    return dedx*pathlen;
  } else {
// subdivide the material
    unsigned nstep = std::min(int(pathlen/maxstep) + 1,maxnstep);
    double step = pathlen/nstep;
    double energy = particleEnergy(mom,mass);
    double deltae = step*dedx;
    double newenergy(energy+deltae);
    double eloss(deltae);
    for(unsigned istep=0;istep<nstep-1;istep++){
      if(newenergy>mass){
// compute the new dedx given the new momentum
        double newmom = particleMomentum(newenergy,mass);
        deltae = step*dEdx(newmom,loss,mass);
// compute the loss in this step
        eloss += deltae;
        newenergy += deltae;
      } else {
// lost all kinetic energy; stop
        eloss = mass-energy;
        break;
      }
    }
    return eloss;
  }
}  


double 
DetMaterial::energyGain(double mom, double pathlen, double mass) const {
  // make sure we take positive lengths!
  pathlen = fabs(pathlen);
  double dedx = dEdx(mom,loss,mass);
// see how far I can step, within tolerance, given this energy loss
  double maxstep = maxStepdEdx(mom,mass,dedx);
// if this is larger than my path, I'm done
  if(maxstep>pathlen){
    return -dedx*pathlen;
  } else {
// subdivide the material
    unsigned nstep = std::min(int(pathlen/maxstep) + 1,maxnstep);
    double step = pathlen/nstep;
    double energy = particleEnergy(mom,mass);
    double deltae = -step*dedx;
// move to the middle of the slice of material
    double newenergy(energy+deltae);
    double egain(deltae);
    for(unsigned istep=0;istep<nstep-1;istep++){
// compute the new dedx given the new momentum
      double newmom = particleMomentum(newenergy,mass);
      double deltae = -step*dEdx(newmom,loss,mass);
      egain += deltae;
      newenergy += deltae;
    }
    return egain;
  }
}  
//
// calculate the energy deposited in an absorber. That's similiar to 
// energyLoss, but the delta electron correction in the Bethe Bloch is
// switched on, and there is no Bremsstrahlung 
//
double 
DetMaterial::energyDeposit(double mom, double pathlen, double mass) const {
  double dedx = dEdx(mom,deposit,mass);
  return dedx*fabs(pathlen);
}  

/********************** end of New Routines **************************/ 

//
//  RMS of energy loss.  This is a gaussian approximation, stolen from
//  Geant3 (see Phys332)
//
double
DetMaterial::energyLossRMS(double mom,double pathlen,double mass) const {
  const double minkappa(1.0e-4);
  double beta = particleBeta(mom,mass);
  double emax = _emax(mom,mass);
  double xi = _xi(beta,fabs(pathlen));
  double kappa = xi/emax;
  double gam = sqrt(1.0-0.5*sqr(beta));;
// formula comes from GFLUCT.F in gphys dnb Jun 4 2004
//
// this formula seriously overestimates the rms when kappa<0.001
// This only really affects electrons
// as for heavier particles resolution effects already dominate when we get to
// this range.  I'll truncate
  if(kappa < minkappa)kappa = minkappa;
  double elossrms = xi*sqrt(gam/kappa);
//  cout << "beta = " << beta
//       << " emax = " << emax
//       << " xi = " << xi
//       << " kappa = " << kappa
//       << " elossrms = " << elossrms << endl;
  return elossrms;
}

//
//  Functions needed for energy loss calculation, see reference above
//
double
DetMaterial::_emax(double mom,double mass){
	static double emass =Pdt::mass(PdtPid::electron);
  double beta = particleBeta(mom,mass);
  double gamma = particleGamma(mom,mass);
	double mratio = emass/mass;
  double emax = 2*emass*sqr(beta)*sqr(gamma)/
    (1+2*gamma*mratio + sqr(mratio));
  if(mass <= emass)
    emax *= 0.5;
  return emax;
}

double
DetMaterial::_xi(double beta,double pathlen) const{
  return _dgev*_za*_density*fabs(pathlen)/sqr(beta);
}

void
DetMaterial::print(ostream& os) const {
  os << "Material " << _name << endl;
}

void
DetMaterial::printAll(ostream& os) const {
  os << "Material " << _name << " has properties : " << endl
  << "  Effective Z = " << _zeff << endl
  << "  Effective A = " << _aeff << endl
  << "  Density (g/cm^3) = " << _density*cm*cm*cm << endl
  << "  Radiation Length (g/cm^2) = " << _radthick*cm*cm<< endl
  << "  Interaction Length (g/cm^2) = " << _intLength << endl
//   << "  Mean Ionization energy (MeV) = " << _meanion << endl
  << "  Mean Ionization energy (MeV) = " << _eexc << endl;
}

double
DetMaterial::maxStepdEdx(double mom,double mass,double dEdx,double tol) {
// compute betagamma at entrance
  double betagamma = particleBetaGamma(mom,mass);
  double energy = particleEnergy(mom,mass);
// basic calculation, based on constant dE/dx
  if(dEdx<0.0){
    double maxstep = -tol*energy/dEdx;
// Modify for steep rise at low momentum
    if(betagamma<2.0)
      maxstep *= sqr(betagamma)/betapower;
    return maxstep;
  }
  else
    return 1.0e6; // large step 
}

// ------------------------------------------------------------------------------
// File and Version Information:
//      $Id: DetMaterial.hh,v 1.35 2008/07/21 23:33:11 brownd Exp $
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
//------------------------------------------------------------------------------
#ifndef DETECTORMATERIAL_HH 
#define DETECTORMATERIAL_HH


#include <math.h>
#include <algorithm>
//
//  Babar includes
//
#include "BaBar/Constants.hh"
#include "BaBar/PdtPid.hh"
#include "BaBar/Pdt.hh"
#include "DetectorModel/DetMtrProp.hh"
#include <iostream>
#include <string>
#include <vector>

class DetMaterial{
public:
  enum dedxtype {loss=0,deposit};
//
//  Constructor
  // new style
  DetMaterial(const char* detName, const DetMtrProp* detMtrProp);

  ~DetMaterial();
//
//  Access
//
  const std::string* name() const {return &_name;}

//
//  comparison function
//
  bool operator == (const DetMaterial& other) const {
    return _name == other._name; }
			
// DNB 3/13/00.  added energyGain to pair with energyLoss, to account
// for cases where the momentum is known _after_ traversal through the
// material as opposed to before.  Both functions are now implemented
// to subdivide the material into steps if it's too thick (dE/dx calculated
// using the exit momentum differs more than 1% of dE/dx for entrance
// momentum).
  double energyLoss(double mom,double pathlen,
		PdtPid::PidType ipart=PdtPid::pion) const {
			return energyLoss(mom,pathlen,Pdt::mass(ipart)); }
	double energyLoss(double mom,double pathlen,double mass) const;

  double energyGain(double mom,double pathlen,
		PdtPid::PidType ipart=PdtPid::pion) const {
			return energyGain(mom,pathlen,Pdt::mass(ipart));	}
	double energyGain(double mom,double pathlen, double mass) const;

  double energyDeposit(double mom, double pathlen,
		PdtPid::PidType ipart = PdtPid::pion) const {
			return energyDeposit(mom,pathlen,Pdt::mass(ipart));	}
	double energyDeposit(double mom, double pathlen,double mass) const;

// raw dE/dx function, used by the above
  double dEdx(double mom,dedxtype type=loss,
		PdtPid::PidType ipart=PdtPid::pion) const {
			return dEdx(mom,type,Pdt::mass(ipart));	}
		
	double dEdx(double mom,dedxtype type,double mass) const;
// 'error' on energy loss (or gain)
  double energyLossRMS(double mom,double pathlen,
		PdtPid::PidType ipart = PdtPid::pion) const {
			return energyLossRMS(mom,pathlen,Pdt::mass(ipart));	}

	double energyLossRMS(double mom,double pathlen,double mass) const;
// functions to describe scattering
// highland formula (without log term)
  double highlandSigma(double mom,double pathlen,
  PdtPid::PidType ipart=PdtPid::pion) const {
    return highlandSigma(mom,pathlen,Pdt::mass(ipart));
  }
  double highlandSigma(double mom,double pathlen, double mass) const;
// average number of singla scatters
  double nSingleScatter(double mom,double pathlen,
  PdtPid::PidType ipart=PdtPid::pion) const {
    return nSingleScatter(mom,pathlen,Pdt::mass(ipart));
  }
  double nSingleScatter(double mom,double pathlen, double mass) const;
// terms used in first-principles single scattering model
  double aParam(double mom) const { return 2.66e-6*pow(_zeff,0.33333333333333)/mom; }
  double bParam(double mom) const { return    0.14/(mom*pow(_aeff,0.33333333333333)); }
//
// Single Gaussian approximation, used in Kalman filtering
  double scatterAngleRMS(double mom,double pathlen,
  PdtPid::PidType ipart = PdtPid::pion) const {
    return scatterAngleRMS(mom,pathlen,Pdt::mass(ipart));
  }
  double scatterAngleRMS(double mom,double pathlen,double mass) const;
//
//  Generic kinematic functions
//
  static double particleEnergy(double mom,PdtPid::PidType ipart)  {
    return particleEnergy(mom,Pdt::mass(ipart)); }
	static double particleEnergy(double mom,double mass) {
		return sqrt(pow(mom,2)+pow(mass,2)); }

  static double particleMomentum(double energy,PdtPid::PidType ipart)  {
    return particleMomentum(energy,Pdt::mass(ipart)); }
	static double particleMomentum(double energy,double mass) {
		return sqrt(particleMom2(energy,mass)); }

  static double particleMom2(double energy,PdtPid::PidType ipart)  {
    return particleMom2(energy,Pdt::mass(ipart)); }
	static double particleMom2(double energy,double mass)  {
		return std::max(pow(energy,2)-pow(mass,2),0.0); }

  static double particleBeta(double mom,PdtPid::PidType ipart) {
    return particleBeta(mom,Pdt::mass(ipart)); }
	static double particleBeta(double mom,double mass) {
	  return mom/particleEnergy(mom,mass); }

  static double particleGamma(double mom,PdtPid::PidType ipart)  {
    return particleGamma(mom,Pdt::mass(ipart)); }
	static double particleGamma(double mom,double mass)  {
	  return particleEnergy(mom,mass)/mass; }
	
	static double particleKinEnergy(double E, PdtPid::PidType ipart)  {
	  return particleKinEnergy(E,Pdt::mass(ipart)); }
	static double particleKinEnergy(double E, double mass)  {
	  return E - mass; }

  static double particleBetaGamma(double mom,PdtPid::PidType ipart) {
    return particleBetaGamma(mom,Pdt::mass(ipart)); }
	static double particleBetaGamma(double mom,double mass) {
	  return mom/mass; }

//
//  functions used to compute energy loss
//
    static double eloss_emax(double mom,double mass) ;
    double eloss_xi(double beta,double pathlen) const;
  double kappa(double mom,double pathlen,double mass) const { 
    return eloss_xi(particleBeta(mom,mass),pathlen)/eloss_emax(mom,mass);}
//
// return the maximum step one can make through this material
// for a given momentum and particle type without dE/dx changing
// by more than the given tolerance (fraction).  This is an _approximate_ 
// function, based on a crude model of dE/dx.
  static double maxStepdEdx(double mom,double mass,
			    double dEdx,double tol=0.05);
private:
//
//  Constants used in material calculations
//
  static double _msmom; // constant in Highland scattering formula
  static double _minkappa; // ionization randomization parameter
  static double _dgev; // energy characterizing energy loss
  static const double _alpha; // fine structure constant
  static double _scatterfrac; // fraction of scattering distribution to include in RMS

//
//  Specific data for this material
//
  std::string _name;
  double _za; // ratio atomic number to atomic weight
  double _zeff; // effective Z of our material
  double _aeff; // effective Z of our material
  double _radthick; // radiation thickness in g/cm**2
  double _intLength; // ineraction length from MatMtrObj in g/cm**2
  double _meanion; // mean ionization energy loss
  double _eexc; // mean ionization energy loss for new e_loss routine
  double _x0; /*  The following specify parameters for energy loss. see
		  Sternheimer etal,'Atomic Data and
		  Nuclear Data Tables', 1984 (40) 267 */
  double _x1;
  double _delta0; 
  double _afactor;
  double _mpower;
  double _bigc;
  double _density;

  int _noem;
  std::vector< double >* _shellCorrectionVector;
  std::vector< double >* _vecNbOfAtomsPerVolume;
  std::vector< double >* _vecTau0;
  std::vector< double >* _vecAlow;
  std::vector< double >* _vecBlow;
  std::vector< double >* _vecClow;
  std::vector< double >* _vecZ;
  double _taul;
  
// cached values to speed calculations
  double _invx0;
  double _nbar;
  double _chic2;
  double _chia2_1;
  double _chia2_2;

public:
// baseic accessors
  double ZA()const {return _za;}
	double zeff() const { return _zeff;}
	double aeff() const { return _aeff;}
  double radiationLength()const {return _radthick;}
  double intLength()const {return _intLength;}
  double meanIon()const {return _meanion;}
	double eexc() const { return _eexc; }
  double X0()const {return _x0;}
  double X1()const {return _x1;}
  double delta0()const {return _delta0;}
  double aFactor()const {return _afactor;}
  double mPower()const {return _mpower;}
  double bigC()const {return _bigc;}
  double density()const {return _density;}
  double inverseX0() const { return _invx0; }
// returns fraction of radiation lengths traversed for a given
// physical distance through this material
  double radiationFraction(double pathlen) const {
    return _density*pathlen/_radthick; }
  void print(std::ostream& os) const;
  void printAll(std::ostream& os ) const;

// parameters used in ionization energy loss
  static double energyLossScale() { return _dgev; }
  static void setEnergyLossScale(double dgev) { _dgev = dgev; }
// parameters used in ionization energy loss randomization
  static double minKappa() { return _minkappa; }
  static void setMinimumKappa(double minkappa) { _minkappa = minkappa; }
// scattering parameter
  static double scatterFraction() { return _scatterfrac;}
  static void setScatterFraction(double scatterfrac) {_scatterfrac = scatterfrac;}
};
#endif


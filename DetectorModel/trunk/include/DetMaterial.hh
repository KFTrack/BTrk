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
#include "PDT/Pdt.hh"
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
//
//  Functions; Multiple scattering RMS,energy loss, and loss RMS
//  (pion is the default particle).
//
  double scatterAngleRMS(double mom,double pathlen,
			 PdtPid::PidType ipart = PdtPid::pion) const {
				return scatterAngleRMS(mom,pathlen,Pdt::mass(ipart));	}
	double scatterAngleRMS(double mom,double pathlen,double mass) const;
			
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

private:
//
//  Constants used in material calculations
//
  static const double _msmom; // momentum characterizing multiple scattering
  static const double _dgev; // energy characterizing energy loss
//
//  primitive functions; many of these are static since they
//  don't depend on the data members
//
  static double _gamma(double beta)  {
    return 1.0/sqrt(1-pow(beta,2)); }
  static double _emax(double mom,double mass) ;
  double _xi(double beta,double pathlen) const;
//
// return the maximum step one can make through this material
// for a given momentum and particle type without dE/dx changing
// by more than the given tolerance (fraction).  This is an _approximate_ 
// function, based on a crude model of dE/dx.
  static double maxStepdEdx(double mom,double mass,
			    double dEdx,double tol=0.05);
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
// returns fraction of radiation lengths traversed for a given
// physical distance through this material
  double radiationFraction(double pathlen) const {
    return _density*pathlen/_radthick; }
  void print(std::ostream& os) const;
  void printAll(std::ostream& os ) const;

  static double scatterMomentum() { return _msmom; }
//  static void setScatterMomentum(double mom) { _msmom = mom; }
    
};
#endif


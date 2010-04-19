// $Id: BFieldIntegrator.hh 497 2010-01-14 09:06:53Z stroili $
// Description:	Class Header for |BFieldIntegrator|
//              Do integral B*dl along track
// Author List:A. Snyder, Copyright (C) 1998	SLAC
#ifndef BFieldIntegrator_HH
#define BFieldIntegrator_HH

#include "BaBar/BaBar.hh"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Matrix/Vector.h"
#include <vector>

class Trajectory;
class BField;



// class interface //
class BFieldIntegrator {

public:

  //construct
  BFieldIntegrator
  (const BField &bField);	// full field

public:

  //destroy
  virtual ~BFieldIntegrator();

public:

  //integrate
  Hep3Vector deltaMomentum
  (const Trajectory *traj,	// trajectory to integrate
   double range[2])const;	// integration range

  Hep3Vector deltaMomentum	
  (const Trajectory *traj,	// trajectory to integrate
   double slo,			// lower integration limit
   double sup)const;		// upper integration limit

  void 
  divideRange	
  (const Trajectory *traj,	// trajectoyr to divide
   double range[2],		// range of interest
   std::vector<double>&)const;	// STL list of site positions
 
  void 
  divideRange
  (const Trajectory *traj,	// trajector to divide
   double slo,			// lower end of range
   double sup,			// upper end of range
   std::vector<double>& posList)const;// list of site positions

public:
  
  //set
  void setTolerance(double val) {_tolerance=val;}
  void setPathMin(double val) {_pathMin=val;}
  void setStepFrac(double val) {_stepFrac=val;}
  void setStepCeiling(double val) {_stepCeiling=val;}
  void setDivTolerance(double val) {_divTolerance=val;}
  void setDivPathMin(double val) {_divPathMin=val;}
  void setDivStepFrac(double val) {_divStepFrac=val;}
  void setDivStepCeiling(double val) {_divStepCeiling=val;}

  //access
  double tolerance()const {return _tolerance;}
  double pathMin()const {return _pathMin;}
  double stepFrac() const {return _stepFrac;}
  double stepCeiling() const {return _stepCeiling;}
  double divTolerance()const {return _divTolerance;}
  double divPathMin()const {return _divPathMin;}
  double divStepFrac() const {return _divStepFrac;}
  double divStepCeiling() const {return _divStepCeiling;}

protected:


private:
  //data

  const BField &_field;		// full field
  Hep3Vector _bNominal;		// nominal field
  double _tolerance;		// tolerence of field variation
  double _pathMin;		// minimum path size to split
  double _stepFrac;	        // fraction of radius for step
  double _stepCeiling;	        // maximum step size
  double _divTolerance;		// tolerence of field variation
  double _divPathMin;		// minimum path size to split
  double _divStepFrac;	        // fraction of radius for step
  double _divStepCeiling;	// maximum step size

public:
  //functions
   
  const BField &field()const {return _field;}
  double bNominal()const {return _bNominal.z();}
  const Hep3Vector& nominalField() const { return _bNominal; }
  
private:
  // Preempt copy constructor and operator=
  BFieldIntegrator&   operator= (const BFieldIntegrator&);
  BFieldIntegrator(const BFieldIntegrator &);
};
#endif




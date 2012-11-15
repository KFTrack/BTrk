// $Id: BFieldIntegrator.hh 830 2011-04-09 07:14:58Z brownd $
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

struct BFieldIntConfig {
  double _maxRange; // maximum bfield integration range
  double _intTolerance;		// tolerence of field variation
  double _intPathMin;		// minimum path size to split
  double _divTolerance;		// tolerence of field variation
  double _divPathMin;		// minimum path size to split
  double _divStepCeiling;	// maximum step size
};

struct BFieldIntRange{
  double _slo;
  double _smid;
  double _shi;
  double range() const { return _shi-_slo; }
  void invert() { double temp = _slo; _slo = _shi; _shi = temp; }
  BFieldIntRange(double slo, double shi) : _slo(slo), _smid(0.5*(slo+shi)),_shi(shi) {}
};

// class interface //
class BFieldIntegrator {

public:

  //construct
  BFieldIntegrator(const BField &bField,BFieldIntConfig const& config);
  BFieldIntegrator(const BField &bField,BFieldIntConfig const& config,const Hep3Vector& bnom);

  //destroy
  virtual ~BFieldIntegrator();

  //integrate
  Hep3Vector deltaMomentum
  (const Trajectory *traj,BFieldIntRange const& range) const;
// subdivide a range
  void divideRange(const Trajectory *traj,	// trajectoyr to divide
    BFieldIntRange const& range, // initial range
   std::vector<BFieldIntRange>& rdiv)const; // subdivided ranges
 // accessors 
  const BFieldIntConfig& config() const { return _config; }
  const BField &field()const {return _field;}
  double bNominal()const {return _bNominal.z();}
  const Hep3Vector& nominalField() const { return _bNominal; }

private:
  //data

  const BField& _field;		// full field
  BFieldIntConfig _config;
  Hep3Vector _bNominal;		// nominal field
  // Preempt copy constructor and operator=
  BFieldIntegrator&   operator= (const BFieldIntegrator&);
  BFieldIntegrator(const BFieldIntegrator &);
};
#endif




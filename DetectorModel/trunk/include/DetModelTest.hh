//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DetModelTest.hh,v 1.9 2003/05/02 16:42:13 brownd Exp $
//
// Description:
//      Test module for generic DetectorModels.  This intersects a
//      set of straight-line trajectories to the (input) detector model,
//      and stores the resulting intersections in an ntuple.
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Dave Brown   12/18/97
//
// Copyright Information:
//	Copyright (C) 1997    LBL
//
//------------------------------------------------------------------------

#ifndef DETMODELTEST_HH
#define DETMODELTEST_HH

class HepTupleManager;
class HepTuple;
class DetSet;
class TrkSimpTraj;
class Trajectory;
class HepHistogram;
#include "PDT/PdtPid.hh"

class DetModelTest {
public:
  enum trajtype { line=0,helix=1};
// The test occurs on construction of this object
  DetModelTest(const DetSet& testset,HepTupleManager* tmanager,
	       double range, // defines the trajectory length
	       int nct,double ctrange[2],
	       int nphi,double phirange[2],
	       PdtPid::PidType hypo=PdtPid::pion,
	       bool validateNextIntersection=false);
// test using a helix trajectory
  DetModelTest(const DetSet& testset,HepTupleManager* tmanager,
	       double range, // defines the trajectory length
	       int nct,double ctrange[2],
	       int nphi,double phirange[2],
	       double curvature,
	       PdtPid::PidType hypo=PdtPid::pion,
	       bool validNextIntersection=false);
// generic constructor: should supersceede the above
  DetModelTest(const DetSet& testset,HepTupleManager* tmanager,
	       double range, // defines the trajectory length
	       int nct,double ctrange[2],
	       int nphi,double phirange[2],
	       double curvmom, // curvature or momentum, depending on the following
	       trajtype type,
	       PdtPid::PidType hypo=PdtPid::pion,
	       bool validateNextIntersection=false);
  DetModelTest(const DetSet& testset,HepTupleManager* tmanager,
	       double range[2], // defines the trajectory length
	       int nct,double ctrange[2],
	       int nphi,double phirange[2],
	       double curvmom, // curvature or momentum, depending on the following
	       trajtype type,
	       PdtPid::PidType hypo=PdtPid::pion,
	       bool validateNextIntersection=false);
  virtual ~DetModelTest(); // might someone want to subclass?
// return the output heptuple
  HepTuple* intersectionTuple() { return _intersectntup; }
private:
  HepTuple* _intersectntup; // intersection test results ntuple
  HepHistogram* _scat; // profiles of material affects vs flightlength
  HepHistogram* _radlen; 
  HepHistogram* _radlensum;
  HepHistogram* _de;
  HepHistogram* _desum;
  HepHistogram* _count;
  TrkSimpTraj* _traj; // trajectory that gets moved around
  double _range[2];
  bool _validateNext; // whether to tests next intersection or not
  double _mom; // initial momentum
  trajtype _type;
  PdtPid::PidType _hypo; // particle type
  void intersect(const DetSet& testset,
		 int nct,double ctrange[2],
		 int nphi,double phirange[2]);
  void rotate(double,double);
  static unsigned _nprofbins;
  static void integrate(const HepHistogram* in,
			unsigned norm,
			HepHistogram* out);
};


#endif




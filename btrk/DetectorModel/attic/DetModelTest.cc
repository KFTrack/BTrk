//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: DetModelTest.cc,v 1.28 2004/08/06 05:58:31 bartoldu Exp $
//
// Description:
//      class DetModelTest
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
#include "BaBar/BaBar.hh"
#include "ErrLogger/ErrLog.hh"
#include "DetectorModel/DetModelTest.hh"
#include "DetectorModel/DetSet.hh"
#include "DetectorModel/DetMaterial.hh"
#include "HepTuple/TupleManager.h"
#include "HepTuple/Tuple.h"
#include "TrkBase/TrkExchangePar.hh"
#include "TrkFitter/TrkDifLineTraj.hh"
#include "TrkBase/TrkMomCalculator.hh"
#include "TrkBase/HelixTraj.hh"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Geometry/Transformation.h"
#include "CLHEP/Utilities/CLHEP.h"
#include "GenEnv/GenEnv.hh"
#include "TrkEnv/TrkEnv.hh"
#include "AbsEnv/AbsEnv.hh"

#include <vector>
using std::endl;
using std::ostream;

#define MATCH 1.0e-5
unsigned DetModelTest::_nprofbins(1000); // number of scatter-plot bins
//
//  All action takes place in the constructor
//
DetModelTest::DetModelTest(const DetSet& testset,HepTupleManager* tmanager,
			   double range,
			   int nct,double ctrange[2],
			   int nphi,double phirange[2],
			   double curvature,
			   PdtPid::PidType hypo,
			   bool validateNext) :
  _intersectntup(0), _traj(0), _type(helix), _hypo(hypo),_validateNext(validateNext) {
  _range[0] = 0.0;
  _range[1] = range;
//  Book the tuple
  _intersectntup = tmanager->ntuple("DetModelTest");
  _radlen = tmanager->histogram("RadLen vs fltlen",_nprofbins,_range[0],_range[1]);
  _scat = tmanager->histogram("Scattering angle vs fltlen",_nprofbins,_range[0],_range[1]);
  _de = tmanager->histogram("Energy loss vs fltlen",_nprofbins,_range[0],_range[1]);
  _radlensum = tmanager->histogram("RadLen integral vs fltlen",_nprofbins,_range[0],_range[1]);
  _desum = tmanager->histogram("Energy loss integral vs fltlen",_nprofbins,_range[0],_range[1]);
  _count = tmanager->histogram("fltlen",_nprofbins,_range[0],_range[1]);
// helix traj
  _traj = new HelixTraj(TrkExchangePar(0.0,0.0,curvature,0.0,0.0));
// momentum
  const BField* bfield = gblEnv->getTrk()->magneticField();
  Hep3Vector momentum = TrkMomCalculator::vecMom(*_traj,*bfield,0.0);
  _mom = momentum.mag();
  intersect(testset,nct,ctrange,nphi,phirange);
}

DetModelTest::DetModelTest(const DetSet& testset,HepTupleManager* tmanager,
			   double range,
			   int nct,double ctrange[2],
			   int nphi,double phirange[2],
			   PdtPid::PidType hypo,
			   bool validateNext) :
  _intersectntup(0), _traj(0), _type(line),_hypo(hypo),_validateNext(validateNext) {
  _range[0] = 0.0;
  _range[1] = range;
//  Book the tuple
  _intersectntup = tmanager->ntuple("DetModelTest");
  _radlen = tmanager->histogram("RadLen vs fltlen",_nprofbins,_range[0],_range[1]);
  _scat = tmanager->histogram("Scattering angle vs fltlen",_nprofbins,_range[0],_range[1]);
  _de = tmanager->histogram("Energy loss vs fltlen",_nprofbins,_range[0],_range[1]);
  _radlensum = tmanager->histogram("RadLen integral vs fltlen",_nprofbins,_range[0],_range[1]);
  _desum = tmanager->histogram("Energy loss integral vs fltlen",_nprofbins,_range[0],_range[1]);
  _count = tmanager->histogram("fltlen",_nprofbins,_range[0],_range[1]);

//  vector for intersections
  std::vector<DetIntersection> tlist;
//  Build a LineTraj
  _traj =  new TrkDifLineTraj(TrkExchangePar(0.0,0.0,1.0,0.0,0.0));
// set momentum by hand
  _mom = 1.0;
  intersect(testset,nct,ctrange,nphi,phirange);
}

DetModelTest::DetModelTest(const DetSet& testset,HepTupleManager* tmanager,
			   double range,
			   int nct,double ctrange[2],
			   int nphi,double phirange[2],
			   double curvmom,
			   trajtype type,
			   PdtPid::PidType hypo,
			   bool validateNext) :
  _intersectntup(0), _traj(0), _type(type),_hypo(hypo),_validateNext(validateNext) {
  _range[0] = 0.0;
  _range[1] = range;
//  Book the tuple
  _intersectntup = tmanager->ntuple("DetModelTest");
  _radlen = tmanager->histogram("RadLen vs fltlen",_nprofbins,_range[0],_range[1]);
  _scat = tmanager->histogram("Scattering angle vs fltlen",_nprofbins,_range[0],_range[1]);
  _de = tmanager->histogram("Energy loss vs fltlen",_nprofbins,_range[0],_range[1]);
  _radlensum = tmanager->histogram("RadLen integral vs fltlen",_nprofbins,_range[0],_range[1]);
  _desum = tmanager->histogram("Energy loss integral vs fltlen",_nprofbins,_range[0],_range[1]);
  _count = tmanager->histogram("fltlen",_nprofbins,_range[0],_range[1]);

  if(_type == helix){
    _traj = new HelixTraj(TrkExchangePar(0.0,0.0,curvmom,0.0,0.0));
// momentum
    const BField* bfield =  gblEnv->getTrk()->magneticField();
    Hep3Vector momentum = TrkMomCalculator::vecMom(*_traj,*bfield,0.0);
    _mom = momentum.mag();
  } else {
    _traj =  new TrkDifLineTraj(TrkExchangePar(0.0,0.0,1.0,0.0,0.0));
    _mom = curvmom;
  }
  intersect(testset,nct,ctrange,nphi,phirange);
}


DetModelTest::DetModelTest(const DetSet& testset,HepTupleManager* tmanager,
			   double range[2],
			   int nct,double ctrange[2],
			   int nphi,double phirange[2],
			   double curvmom,
			   trajtype type,
			   PdtPid::PidType hypo,
			   bool validateNext) :
  _intersectntup(0), _traj(0), _type(type), _hypo(hypo),_validateNext(validateNext) {
  _range[0] = range[0];
  _range[1] = range[1];
//  Book the tuple
  _intersectntup = tmanager->ntuple("DetModelTest");
  _radlen = tmanager->histogram("RadLen vs fltlen",_nprofbins,_range[0],_range[1]);
  _scat = tmanager->histogram("Scattering angle vs fltlen",_nprofbins,_range[0],_range[1]);
  _de = tmanager->histogram("Energy loss vs fltlen",_nprofbins,_range[0],_range[1]);
  _radlensum = tmanager->histogram("RadLen integral vs fltlen",_nprofbins,_range[0],_range[1]);
  _desum = tmanager->histogram("Energy loss integral vs fltlen",_nprofbins,_range[0],_range[1]);
  _count = tmanager->histogram("fltlen",_nprofbins,_range[0],_range[1]);

  if(_type == helix){
    _traj = new HelixTraj(TrkExchangePar(0.0,0.0,curvmom,0.0,0.0));
// momentum
    const BField* bfield =  gblEnv->getTrk()->magneticField();
    Hep3Vector momentum = TrkMomCalculator::vecMom(*_traj,*bfield,0.0);
    _mom = momentum.mag();
  } else {
    _traj =  new TrkDifLineTraj(TrkExchangePar(0.0,0.0,1.0,0.0,0.0));
    _mom = curvmom;
  }
  intersect(testset,nct,ctrange,nphi,phirange);
}


void
DetModelTest::rotate(double phi, double tanlambda) {
  static HepPoint origin(0,0,0); // all the line trajs start at the origin
  if(_type == line) {
    TrkDifLineTraj* ltraj = (TrkDifLineTraj*)_traj;
    *ltraj = TrkDifLineTraj(TrkExchangePar(0.0,phi,1.0,0.0,tanlambda));
  } else if(_type == helix) {
    HelixTraj* htraj = (HelixTraj*)_traj;
    *htraj = HelixTraj(TrkExchangePar(0.0,phi,htraj->omega(),0.0,tanlambda));
  }
  _traj->setFlightRange(_range);
}

void
DetModelTest::intersect(const DetSet& testset,
			int nct,double ctrange[2],
			int nphi,double phirange[2]) {
//   list for intersections
  std::vector<DetIntersection> tlist;
  int iphi,ict;
  double phi,costheta,sintheta,tanl;
  double ctstep,pstep;
  if(nct>1)
    ctstep = (ctrange[1]-ctrange[0])/(nct-1);
  else
    ctstep = (ctrange[1]-ctrange[0]);
  if(nphi>1)
    pstep = (phirange[1]-phirange[0])/(nphi-1);
  else
    pstep = (phirange[1]-phirange[0]);
  for(ict=0;ict<nct;ict++){
    costheta = ctrange[0] + ict*ctstep;
    sintheta = sqrt(1-costheta*costheta);
    tanl = (sintheta!= 0.0) ? costheta/sintheta : 1000*costheta;
    for(iphi=0;iphi<nphi;iphi++){
      phi = phirange[0] + iphi*pstep;
// rotate the trajectory
      rotate(phi,tanl);
//  Intersect the trajectory with the DetectorSet
      testset.intersection(tlist,_traj, 0, true);
//  Loop over the intersections and store them in the ntuple
      int ielem;
      DetIntersection nextinter;
      unsigned nelems=tlist.size();
      const DetElem* prevelem(0);
      double prevlen (0.0);
      double radlensum,desum,dpsum;
      radlensum = desum = dpsum = 0.0;
      for(ielem=0;ielem<nelems;ielem++){
	DetIntersection& dinter(tlist[ielem]);
	const DetElem* delem = dinter.delem;
// check for duplicates
	if(delem == prevelem &&
	   dinter.pathlen == prevlen) {
	  ostream& os = ErrMsg(error);
	  os << "DetModelTest Error: intersection appears multiple times for element ";
	  delem->printAll(os);
	  os << "At flightlength " << dinter.pathlen << endmsg;
	}
	prevelem = delem;
	prevlen = dinter.pathlen;
// test nextIntersection if desired
	if(_validateNext){
	  bool intersected;
	  if(ielem == 0)
	    intersected = testset.firstIntersection(_traj,nextinter);
	  else
	    intersected = testset.nextIntersection(_traj,nextinter,&nextinter);
	  if(! intersected ||
	     fabs(nextinter.pathlen-dinter.pathlen)> MATCH ){
	    ostream& os = ErrMsg(error);
	    os << "DetModelTest: nextIntersection mismatch" << endl;
	    os << "Intersection element = ";
	    delem->print(os);
	    if(intersected){
	      os << "NextIntersection element = ";
	      nextinter.delem->print(os);
	    } else
	      os << "No nextIntersection element" << endmsg;
	  } else if(delem != nextinter.delem){
	    ostream& os = ErrMsg(error);
	    os << "Warning: NextIntersection elements out of order" << endl
			  << "Intersection path = " << dinter.pathlen 
	       << ", element = " << endl;
	    delem->print(os);
	    os << "NextIntersection path = " << nextinter.pathlen 
	       << ", element = " << endl;
	    nextinter.delem->print(os);
	    os << endmsg;
	  }
	}
// fill the ntuple
	_intersectntup->column("costheta",float(costheta));
	_intersectntup->column("phi",float(phi));
	_intersectntup->column("pathlen",float(dinter.pathlen));
	_intersectntup->column("range",float(dinter.pathrange[1]-
			       dinter.pathrange[0]));

	HepPoint enterpoint = dinter.trajet->position(dinter.pathrange[0]);
	HepPoint exitpoint = dinter.trajet->position(dinter.pathrange[1]);

	_intersectntup->column("enterx",float(enterpoint.x()));
	_intersectntup->column("entery",float(enterpoint.y()));
	_intersectntup->column("enterz",float(enterpoint.z()));

	_intersectntup->column("exitx",float(exitpoint.x()));
	_intersectntup->column("exity",float(exitpoint.y()));
	_intersectntup->column("exitz",float(exitpoint.z()));

	HepPoint spoint = dinter.trajet->position(dinter.pathlen);
	HepPoint lpoint = delem->transform().transTo(spoint);
	
	_intersectntup->column("localx",float(lpoint.x()));
	_intersectntup->column("localy",float(lpoint.y()));
	_intersectntup->column("localz",float(lpoint.z()));

	_intersectntup->column("elemid",delem->elementNumber());
	_intersectntup->column("typeid",delem->detectorType()->
						typeNumber());
//
// fill profile histograms
	const DetMaterial& mat = delem->material(dinter);
	double dist = dinter.pathrange[1]-dinter.pathrange[0];
	double radlen = mat.radiationFraction(dist);
	double scat = mat.scatterAngleRMS(_mom,dist,_hypo);
	double de = fabs(mat.energyLoss(_mom,dist,_hypo));
	double dp = DetMaterial::particleMomentum
	  (DetMaterial::particleEnergy(_mom,_hypo)+de,_hypo) - _mom;
	_radlen->accumulate(dinter.pathlen,radlen);
	_scat->accumulate(dinter.pathlen,scat);
	_de->accumulate(dinter.pathlen,de);
	_count->accumulate(dinter.pathlen);
	radlensum += radlen;
	desum += de;
	dpsum += dp;

	_intersectntup->column("radlen",float(radlen));
	_intersectntup->column("scat",float(scat));
	_intersectntup->column("de",float(de));
	_intersectntup->column("dp",float(dp));
	_intersectntup->column("radlensum",float(radlensum));
	_intersectntup->column("desum",float(desum));
	_intersectntup->column("dpsum",float(dpsum));
	_intersectntup->dumpData();
	
      }
    }
  }
// integrate the profile histograms
  integrate(_radlen,nct*nphi,_radlensum);
  integrate(_de,nct*nphi,_desum);
}

DetModelTest::~DetModelTest()
{
  delete _traj;
}

//utility function
void
DetModelTest::integrate(const HepHistogram* input,
			unsigned norm,
			HepHistogram* output) {
  assert(input->getNbins(HepHistogram::dimX) == input->getNbins(HepHistogram::dimX) &&
	 input->getNbins(HepHistogram::dimY) == input->getNbins(HepHistogram::dimY) &&
	 input->getLow(HepHistogram::dimX) == input->getLow(HepHistogram::dimX) &&
	 input->getHigh(HepHistogram::dimX) == input->getHigh(HepHistogram::dimX) );
  output->reset();
  unsigned nbins = input->getNbins(HepHistogram::dimX);
  double sum(0.0);
  double binsize = (input->getHigh(HepHistogram::dimX) -
		    input->getLow(HepHistogram::dimX))/float(nbins);
  double center = input->getLow(HepHistogram::dimX)+binsize/2.0;
  for (unsigned ibin=0;ibin<nbins;ibin++){
    sum += input->getContents(ibin)/float(norm);
    output->accumulate(center,sum);
    center += binsize;
  }
}






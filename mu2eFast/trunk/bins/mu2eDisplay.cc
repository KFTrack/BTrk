/* testOneTrack
*
	* Generate one track, run it through KalmanTrack, and analyze the results.
*/


/*
	* Notes to myself:
*
	* Potential methods of comparison for a single simulated/reconstructed trajectory pair:
* - compare helix parameters at beginning and end of track
	* - compute absolute difference in L^1 or L^2 norm (i.e. integrate |X_g - X_r|^2)
	* - look at chi-square statistics for goodness-of-fit
	* - visually display the trajectories (along with the generated trajectory)
*/

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include <string.h>

#include <Gtypes.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>

#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Vector.h"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Random/Random.h"
#include "CLHEP/Random/RanecuEngine.h"
//#include "CLHEP/config/TemplateFunctions.h"
#include "KalmanTrack/KalRep.hh"
#include "KalmanTrack/KalInterface.hh"
#include "PDT/Pdt.hh"
#include "PDT/PdtEntry.hh"
#include "PDT/PdtPid.hh"
#include "TrajGeom/TrkLineTraj.hh"
#include "TrkBase/TrkDifPieceTraj.hh"
#include "TrkBase/HelixTraj.hh"
#include "TrkBase/TrkDifPieceTraj.hh"
#include "TrkBase/TrkMomCalculator.hh"
#include "TrkBase/TrkPoca.hh"
#include "TrkBase/TrkSimpTraj.hh"
#include "TrkBase/TrkHelixUtils.hh"
#include "TrkBase/TrkRecoTrk.hh"
#include "TrkBase/TrkExchangePar.hh"
#include "BbrGeom/TrkPieceTraj.hh"

#include "PacEnv/PacConfig.hh"
#include "PacDetector/PacCylDetector.hh"
#include "PacEnv/PacBuildEnv.hh"
#include "PacSim/PacSimulate.hh"
#include "PacTrk/PacReconstructTrk.hh"
#include "PacSim/PacSimTrack.hh"
#include "PacSim/PacSimHit.hh"
#include "PacGeom/PacHelix.hh"
#include "PacGeom/PacPieceTraj.hh"
#include "PacEnv/pstring.hh"
#include "PacDisplay/PacEvtDisplay.hh"

#include "ProxyDict/Ifd.hh"
#include "ProxyDict/IfdDataProxyUnowned.hh"
#include "AbsEnv/AbsEnv.hh"
#include "DetectorModel/DetSet.hh"
#include "BField/BField.hh"

#include "Framework/AppFileName.hh"
#include <TParticle.h>
#include <TParticlePDG.h>

#include "TF1.h"


	using namespace std;


/*******************************************************************************
* Utility functions
	*******************************************************************************/


double IntegratedDifference(const Trajectory* A, const Trajectory* B, double smin = -1, double smax = -1, int M = 100) {
		/* Determine the common range of validity of the two trajectories */
	if(smin == -1)
		smin = fmax(A->lowRange(), B->lowRange());
	if(smax == -1)
		smax = fmin(A->hiRange(), B->hiRange());
	if(smin > smax) {
		cout << "IntegratedDifference: range mismatch, (smin, smax) = (" << smin << ", " << smax << ")" << endl;
		return 0;
	}

		/* Perform the integration numerically using the trapezoid rule */
	double ds = (smax - smin)/M;
	double integral = 0;
	integral += ds * 0.5 * (A->position(smin) - B->position(smin)).mag();
	for(int m = 1; m <= M-1; m++) {
		double s = smin + m * ds;
		integral += ds * (A->position(s) - B->position(s)).mag();
	}
	integral += ds * 0.5 * (A->position(smax) - B->position(smax)).mag();

	return integral;
}

void SaveDifferencePlot(const char* name, const Trajectory* A, const Trajectory* B, double smin = -1, double smax = -1, int M = 100) {
		/* Determine the common range of validity of the two trajectories */
	if(smin == -1)
		smin = fmax(A->lowRange(), B->lowRange());
	if(smax == -1)
		smax = fmin(A->hiRange(), B->hiRange());
	if(smin > smax) {
		cout << "SaveDifferencePlot: range mismatch, (smin, smax) = (" << smin << ", " << smax << ")" << endl;
		return;
	}

	double s[M+1];
	for(int m = 1; m <= M; m++)
		s[m] = smin + m * (smax - smin)/M;

	double diff[M+1];
	for(int m = 1; m <= M; m++)
		diff[m] = (A->position(s[m]) - B->position(s[m])).mag();

	TGraph g(M+1, s, diff);
	g.Write(name);
}

/*******************************************************************************
* Main
	*******************************************************************************/

int main(int argc, char* argv[]) {

    // get configuration
    //
  if(argc <= 1){
    cout << "reading standard config file" << endl;
    gconfig.parsefile(AppFileName("mu2eFast/mu2e.xml").pathname().c_str());
    gconfig.parsefile(AppFileName("mu2eFast/mu2e_display.xml").pathname().c_str());
  } else {
    cout << "reading config from arguments ";
    for(unsigned iarg=1;iarg<argc;iarg++){
      cout << argv[iarg] << " ";
    }
    cout << endl;  
    gconfig.parseargs(argc, argv);
	}
//  gconfig.dump();
  
// build environment.  This creates the BField, PID tables, materials, etc.
	PacBuildEnv penv;
	penv.buildCore();
	penv.buildTrk();
	
	// get back the bfield
	const BField* bfield = Ifd<BField>::get(gblPEnv,"Default");

// build detector
	PacCylDetector* detector = new PacCylDetector();

	/* put the DetectorModel into the event */
  IfdDataProxyUnowned<DetSet>* dsproxy= new IfdDataProxyUnowned<DetSet>(const_cast<DetSet*>(detector->detectorModel()));
	if(! Ifd< DetSet >::put( gblPEnv,dsproxy,"Tracking Set"))
		cout << "Can't put Detector Set" << endl;
		
	if(! Ifd< PacDetector >::put( gblPEnv,detector,"Tracking Det"))
  	cout << "Can't put Detector" << endl;

// if requested, setup display
  PacEvtDisplay display;
  bool disptrack = gconfig.getbool("displaytrack");
  bool printtrack = gconfig.getbool("printtrack");
	if(disptrack){
    display.init(gconfig.getcstr("displayfile"),gconfig.getint("displayresolution"));
    display.drawDetector();
    display.drawVoxels();
    display.reset();
	}

		/* Get helix parameters */
	HepPoint initpos = gconfig.getHepPoint("initpos");
	Hep3Vector initmom = gconfig.getHep3Vector("initmom");
	PdtPdg::PdgType pdgid = (PdtPdg::PdgType)gconfig.getint("PdtPdg",13);
	PdtEntry* pdt = Pdt::lookup(pdgid);

// create particle

  TParticle* part = new TParticle();
  part->SetPdgCode(pdgid);
  double mass = pdt->mass();
  double energy = sqrt(mass*mass+initmom.mag2());
  part->SetMomentum(initmom.x(),initmom.y(),initmom.z(),energy);
  part->SetProductionVertex(initpos.x(),initpos.y(),initpos.z(),0.0);
  part->SetWeight(1.0); // all particles have same weight
  part->SetStatusCode(1);

  cout << "Simulating particle type ";
  pdt->printOn(cout);
 	
	double q = pdt->charge();
		/* Generate the track */
		
	PacSimulate sim(bfield,detector);
	HepRandomEngine* engine = new RanecuEngine();
	HepRandom::setTheEngine(engine);
	sim.setRandomEngine(engine);
	detector->setRandomEngine(engine);
	PacSimTrack* simtrk = sim.simulateParticle(part);
	const PacPieceTraj* simtraj = simtrk->getTraj();
	HepVector gparams(5);
	double fltlen;
	TrkHelixUtils::helixFromMom(gparams,fltlen,initpos,initmom,q,*bfield);
	PacHelix gentraj(gparams,fltlen,simtraj->hiRange());
  if(printtrack){
    cout << "Initial parameters " << gparams << endl;
    const std::vector<PacSimHit>& hits = simtrk->getHitList();
    cout << "[SimHit]" << endl;
    for(int ihit = 0; ihit < hits.size(); ihit++){
      cout << " SimHit " << ihit << ": " << hits[ihit].position()
        << " effect " << PacSimHit::detEffectName(hits[ihit].detEffect()) 
        << " global len " << hits[ihit].globalFlight()
        << " dmom " << hits[ihit].momChange();
      if(hits[ihit].detIntersection().delem != 0){
        cout << " element " << hits[ihit].detIntersection().delem->elementName() 
          << " " <<  hits[ihit].detIntersection().delem->elementNumber();
      }
      cout << endl;
    }
  }
  if(disptrack){
     display.drawParticle(part,simtrk->lastHit()->globalFlight(),bfield);
     display.drawSimTrack(simtrk);
     display.drawSimHits(simtrk,0);
   }
  
		/* Reconstruct the track with KalmanTrack (using the list of hits) */
	PacReconstructTrk trackreco(bfield,penv.getKalContext());
	trackreco.setRandomEngine(engine);
	std::vector<const PacSimTrack*> strks;
  strks.push_back(simtrk);
  trackreco.makeTracks(strks);
	const TrkRecoTrk* trk = trackreco.findTrack(simtrk);
	if(trk != 0 && trk->status() != 0 && trk->status()->fitCurrent() ){

		KalInterface kinter;
		trk->attach(kinter,trk->defaultType());
		const KalRep* kalrep = kinter.kalmanRep();
		const TrkDifPieceTraj& recotraj = kalrep->pieceTraj();
//		const TrkHotList* hotlist = kalrep->hotList();
		TrkExchangePar helix = kalrep->helix(0);
		HepVector recoparams = helix.params();

		
		TrkLineTraj zaxis(HepPoint(0, 0, -10), Hep3Vector(0, 0, 1), 20);
		TrkPoca generatedpoca(gentraj, 0, zaxis, 10, 1e-12);
		TrkPoca genpoca(*simtraj, 0, zaxis, 10, 1e-12);
		TrkPoca recopoca(recotraj, 0, zaxis, 10, 1e-12);

    if(printtrack){
      cout << "  initial helix parameters = " << gparams << endl;
      cout << "[Generated trajectory (green)]" << endl;
      cout << "  range = " << gentraj.lowRange() << ", " << gentraj.hiRange() << endl;
      cout << "  z-axis poca: s = " << generatedpoca.flt1() << ", distance = " << generatedpoca.doca() << endl;
      gentraj.printAll(cout);

      cout << "[Simulated trajectory (blue)]" << endl;
      cout << "  range = " << simtraj->lowRange() << ", " << simtraj->hiRange() << endl;
      cout << "  z-axis poca: s = " << genpoca.flt1() << ", distance = " << genpoca.doca() << endl;
      simtraj->printAll(cout);

      cout << "[Reconstructed trajectory (red)]" << endl;
      cout << "  range = " << recotraj.lowRange() << ", " << recotraj.hiRange() << endl;
      cout << "  z-axis poca: s = " << recopoca.flt1() << ", distance = " << recopoca.doca() << endl;
      cout << "  reco helix parameters = " << recoparams << endl;
      recotraj.printAll(cout);

//      kalrep->printAll(cout);

    /* Check that the initial momenta match */
      double localflight;
      const TrkSimpTraj* inithelix = recotraj.localTrajectory(0, localflight);
      const TrkGeomTraj* first_segment = simtraj->localTrajectory(0);
      Hep3Vector recoinitmom = TrkMomCalculator::vecMom(*inithelix, *bfield, localflight);
      cout << "[Initial momentum]" << endl;
      cout << "  generated         = " << initmom << endl;
//		cout << "  simulated     = " << first_segment.momentum(0, bfield->bFieldZ()) << endl;
      cout << "  reconstructed = " << recoinitmom << endl;

      cout << "[Initial position]" << endl;
      cout << "  generated         = " << initpos << endl;
      cout << "  simulated     = " << first_segment->position(0) << endl;
      cout << "  reconstructed = " << recotraj.position(0) << endl;

      cout << "[Fit chi-square]" << endl;
      cout << "  chisq/ndof = " << kalrep->chisq() << "/" << kalrep->nDof() << endl;

      cout << "[Integrated differences]" << endl;
      cout << "  simulated - generated = " << IntegratedDifference(simtraj, &gentraj, 25, 75) << endl;
      cout << "  simulated - reconstructed = " << IntegratedDifference(simtraj, &recotraj) << endl;
      cout << "  generated - reconstructed = " << IntegratedDifference(&gentraj, &recotraj) << endl;
    }

    if(disptrack)
      display.drawRecTrack(trk);
	}
  if(disptrack){
    display.fillTrees();
    display.finalize();
  }  
	return 0;
}

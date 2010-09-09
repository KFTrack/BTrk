/* testTrackReconstruction
*
* Generate several tracks, fit them with KalmanTrack, and analyze the results.
*/

#include <cmath>
#include <iostream>
#include <string>
#include <cstdio>
#include <vector>
#include <time.h>

#include <TTree.h>
#include <TBranch.h>
#include <Gtypes.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>
#include <TRandom3.h>

#include "BField/BFieldFixed.hh"
#include "BField/BFieldIntegrator.hh"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/config/TemplateFunctions.h"
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
#include "DetectorModel/DetSet.hh"
#include "DetectorModel/DetMaterial.hh"
#include "DetectorModel/DetSurfaceElem.hh"

#include "PacEnv/PacConfig.hh"
#include "PacDetector/PacCylDetector.hh"
#include "PacEnv/PacBuildEnv.hh"
#include "PacSim/PacSimulate.hh"
#include "PacTrk/PacReconstructTrk.hh"
#include "PacTrk/PacTrkSimHotMap.hh"
#include "PacSim/PacSimTrack.hh"
#include "PacSim/PacSimHit.hh"
#include "PacSim/PacShowerInfo.hh"
#include "PacGeom/PacHelix.hh"
#include "PacGeom/PacPieceTraj.hh"
#include "PacGeom/PacMeasurement.hh"
#include "PacDisplay/PacEvtDisplay.hh"
#include "mu2eFast/PacSimHitInfo.rdl"
#include "mu2eFast/TrajDiff.rdl"
#include "mu2eFast/PacSimTrkSummary.rdl"
#include "mu2eFast/mu2eDSField.hh"


#include "ProxyDict/Ifd.hh"
#include "ProxyDict/IfdDataProxyUnowned.hh"

#include "AbsEnv/AbsEnv.hh"
#include "BField/BField.hh"
#include "CLHEP/Random/Random.h"
#include "CLHEP/Random/RanecuEngine.h"
#include "CLHEP/Vector/ThreeVector.h"

#include "G3Data/GVertex.hh"
#include "G3Data/GTrack.hh"
#include "GUtils/GTrkUtils.hh"

#include "Framework/AppFileName.hh"

using namespace std;

#define RNGSEED 9082459

void fillSimHitInfo(const PacSimTrack* strk, std::vector<PacSimHitInfo>& sinfo);
void fillTrajDiff(const PacSimTrack* strk, const TrkDifPieceTraj& ptraj, std::vector<TrajDiff>& tdiff);
void fillSimTrkSummary(const PacSimTrack* strk, PacSimTrkSummary& ssum);

PacTrkSimHotMap simHotMap; // used to access nasty statics  

// field integral test stuff
BField* mecofield(0);
BFieldIntegrator* fieldint(0);


int main(int argc, char* argv[]) {
  gconfig.verbose(true);
  if(argc <= 1){
    gconfig.parsefile(AppFileName("mu2eFast/mu2e.xml").pathname().c_str());
    gconfig.parsefile(AppFileName("mu2eFast/mu2e_test.xml").pathname().c_str());
  }
  gconfig.parseargs(argc, argv);

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
  cout << "Built PacCylDetector " << endl;
    
  // if requested, setup display
  PacEvtDisplay display;
  bool disptrack = gconfig.getbool("displaytrack");
	if(disptrack){
    display.init(gconfig.getcstr("displayfile"),gconfig.getint("displayresolution"));
    display.drawDetector();
    display.reset();
	}
  
  // config parameters
  bool hittuple = gconfig.getbool("hittuple",false);
  bool trajdiff = gconfig.getbool("trajdiff",false);
//  bool calodiff = gconfig.getbool("calodiff",false);

    // Read helix generation parameters 
  double p_min = double(gconfig["p_min"]);
  double p_max = double(gconfig["p_max"]);
  double cost_min = double(gconfig["cost_min"]);
  double cost_max = double(gconfig["cost_max"]);
  double r0_mean = double(gconfig["r0_mean"]);
  double r0_sigma = double(gconfig["r0_sigma"]);
  double z0_mean = double(gconfig["z0_mean"]);
  double z0_sigma = double(gconfig["z0_sigma"]);
  
    // Prepare to construct tracks 
  int numtracks = int(gconfig["numtracks"]);
  const int rndseed = gconfig.getint("rndseed", 0);
  TRandom3 rng(RNGSEED + rndseed);

  //Create File		test.root
  const char* outputfile = gconfig.get("outputfile", "mu2e_test.root");
  TFile file(outputfile,"RECREATE");
  //Create Tree to store track info
  TTree* trackT = new TTree("tracks","Tracks");
  //Variables to store track information
  Int_t    itrack;
  Int_t    trknum;
  Float_t sim_d0;
  Float_t sim_phi0;
  Float_t sim_omega;
  Float_t sim_z0;
  Float_t sim_tandip;	
  Float_t sim_lowrange;
  Float_t sim_hirange;
  Float_t sim_poca;
  Float_t sim_doca;
  Float_t sim_mom_pt;
  Float_t sim_mom_z;
  Float_t sim_mom_mag;
  Float_t sim_mom_cost;
  Float_t sim_mom_phi;
  Float_t sim_inipos_x;
  Float_t sim_inipos_y;
  Float_t sim_inipos_z;
  
  Float_t rec_chisqr;
  Float_t rec_fitprob;
  Float_t rec_lowrange;
  Float_t rec_hirange;
  Float_t rec_poca;
  Float_t rec_doca;
  Float_t rec_mom_z;
  Float_t rec_mom_mag;
  Float_t rec_mom_pt;
  Float_t rec_mom_err;
  Float_t rec_mom_pterr;
  Int_t rec_ndof;
  Int_t rec_nsvt;
  Int_t rec_ndch;
  Int_t rec_nactive;
  Int_t rec_nhit;
  Float_t rec_d0;
  Float_t rec_phi0;
  Float_t rec_omega;
  Float_t rec_z0;
  Float_t rec_tandip;
  Float_t reccov_d0;
  Float_t reccov_phi0;
  Float_t reccov_omega;
  Float_t reccov_z0;
  Float_t reccov_tandip;
  Float_t pull_d0;
  Float_t pull_phi0;
  Float_t pull_omega;
  Float_t pull_omegaabs;
  Float_t pull_z0;
  Float_t pull_tandip;
  
  std::vector<PacSimHitInfo> sinfo;
  std::vector<TrajDiff> tdiff;
//  std::vector<TrajDiff> cdiff;
  PacSimTrkSummary ssum;
  
  //Create TBranch to store track info
  trackT->Branch("itrack",&itrack,"itrack/I");
  trackT->Branch("trknum",&trknum,"trknum/I");
  trackT->Branch("sim_d0",&sim_d0,"sim_d0/F");
  trackT->Branch("sim_phi0",&sim_phi0,"sim_phi0/F");
  trackT->Branch("sim_omega",&sim_omega,"sim_omega/F");
  trackT->Branch("sim_z0",&sim_z0,"sim_z0/F");
  trackT->Branch("sim_tandip",&sim_tandip,"sim_tandip/F");
  trackT->Branch("sim_lowrange",&sim_lowrange,"sim_lowrange/F");
  trackT->Branch("sim_hirange",&sim_hirange,"sim_hirange/F");
  trackT->Branch("sim_poca",&sim_poca,"sim_poca/F");
  trackT->Branch("sim_doca",&sim_doca,"sim_doca/F");
  trackT->Branch("sim_lowrange",&sim_lowrange,"sim_lowrange/F");
  trackT->Branch("sim_hirange",&sim_hirange,"sim_hirange/F");
  trackT->Branch("sim_poca",&sim_poca,"sim_poca/F");
  trackT->Branch("sim_doca",&sim_doca,"sim_doca/F");
  trackT->Branch("sim_mom_z",&sim_mom_z,"sim_mom_z/F");
  trackT->Branch("sim_mom_mag",&sim_mom_mag,"sim_mom_mag/F");
  trackT->Branch("sim_mom_cost",&sim_mom_cost,"sim_mom_cost/F");
  trackT->Branch("sim_mom_phi",&sim_mom_phi,"sim_mom_phi/F");
  trackT->Branch("sim_mom_pt",&sim_mom_pt,"sim_mom_pt/F");
  trackT->Branch("sim_inipos_x",&sim_inipos_x,"sim_inipos_x/F");
  trackT->Branch("sim_inipos_y",&sim_inipos_y,"sim_inipos_y/F");
  trackT->Branch("sim_inipos_z",&sim_inipos_z,"sim_inipos_z/F");

  trackT->Branch("rec_d0",&rec_d0,"rec_d0/F");
  trackT->Branch("rec_phi0",&rec_phi0,"rec_phi0/F");
  trackT->Branch("rec_omega",&rec_omega,"rec_omega/F");
  trackT->Branch("rec_z0",&rec_z0,"rec_z0/F");
  trackT->Branch("rec_tandip",&rec_tandip,"rec_tandip/F");
  trackT->Branch("reccov_d0",&reccov_d0,"reccov_d0/F");
  trackT->Branch("reccov_phi0",&reccov_phi0,"reccov_phi0/F");
  trackT->Branch("reccov_omega",&reccov_omega,"reccov_omega/F");
  trackT->Branch("reccov_z0",&reccov_z0,"reccov_z0/F");
  trackT->Branch("reccov_tandip",&reccov_tandip,"reccov_tandip/F");  
  trackT->Branch("rec_lowrange",&rec_lowrange,"rec_lowrange/F");
  trackT->Branch("rec_hirange",&rec_hirange,"rec_hirange/F");
  trackT->Branch("rec_poca",&rec_poca,"rec_poca/F");
  trackT->Branch("rec_doca",&rec_doca,"rec_doca/F");
  trackT->Branch("rec_fitprob",&rec_fitprob,"rec_fitprob/F");
  trackT->Branch("rec_chisqr",&rec_chisqr,"rec_chisqr/F");
  trackT->Branch("rec_ndof",&rec_ndof,"rec_ndof/I");
  trackT->Branch("rec_nsvt",&rec_nsvt,"rec_nsvt/I");
  trackT->Branch("rec_ndch",&rec_ndch,"rec_ndch/I");
  trackT->Branch("rec_nactive",&rec_nactive,"rec_nactive/I");
  trackT->Branch("rec_nhit",&rec_nhit,"rec_nhit/I");
  trackT->Branch("rec_mom_pt",&rec_mom_pt,"rec_mom_pt/F");
  trackT->Branch("rec_mom_err",&rec_mom_err,"rec_mom_err/F");
  trackT->Branch("rec_mom_pterr",&rec_mom_pterr,"rec_mom_pterr/F");
  trackT->Branch("rec_mom_z",&rec_mom_z,"rec_mom_z/F");
  trackT->Branch("rec_mom_mag",&rec_mom_mag,"rec_mom_mag/F");
  
  trackT->Branch("pull_d0",&pull_d0,"pull_d0/F");
  trackT->Branch("pull_phi0",&pull_phi0,"pull_phi0/F");
  trackT->Branch("pull_omega",&pull_omega,"pull_omega/F");
  trackT->Branch("pull_omegaabs",&pull_omegaabs,"pull_omegaabs/F");
  trackT->Branch("pull_z0",&pull_z0,"pull_z0/F");
  trackT->Branch("pull_tandip",&pull_tandip,"pull_tandip/F");
  
// branch for individual simhit info (TClonesArray)
  if(hittuple)
    trackT->Branch("simhit",&sinfo);
// test of trajectory differences
  if(trajdiff){
    trackT->Branch("trajdiff",&tdiff);
    double dfactor = gconfig.getfloat("distortionfactor",1.0);
    std::string fmap = gconfig.get("fieldmap");
    mecofield = new mu2eDSField(fmap,dfactor);
    assert(mecofield != 0);
    fieldint = new BFieldIntegrator(*mecofield);
  }
//  if(calodiff)
//    trackT->Branch("calodiff",&cdiff);
// branch for simtrack summary
  trackT->Branch("simtrk",&ssum.nsimhit,PacSimTrkSummary::rootnames());    

  int tracknum = 1;
  const int printfreq = gconfig.getint("printfreq", 100);

  PacSimulate sim(bfield,detector);
  HepRandomEngine* engine = new RanecuEngine();
  HepRandom::setTheEngine(engine);
  sim.setRandomEngine(engine);
  detector->setRandomEngine(engine);

  PacReconstructTrk trackreco(bfield,penv.getKalContext());

  for(int itrk = 0; itrk < numtracks; itrk++) {
    if(0 == (itrk+1)%printfreq) {
      printf("Count: %i \n",itrk+1);
    }
    itrack = itrk;
// must clear the nasty statics
    simHotMap.Clear();
// clear vectors
    sinfo.clear();
    tdiff.clear();
//    cdiff.clear();
// Generate track parameters; first, origin vertex
    double posphi =  rng.Uniform(0, 2*M_PI);
    double dx = rng.Gaus(0, r0_sigma);
    double dy = rng.Gaus(0, r0_sigma);
    double x = r0_mean*cos(posphi)+dx;
    double y = r0_mean*sin(posphi)+dy;
    double z = rng.Gaus(z0_mean, z0_sigma);
    HepPoint pos(x, y, z);
// now momentum
    double mom	= fabs(rng.Uniform(p_min, p_max));                // transverse momentum
    double phi	= rng.Uniform(0, 2*M_PI);
    double cost = rng.Uniform(cost_min,cost_max);
    double pz	= mom*cost;                // longitudinal momentum
    double pt = mom*sqrt(1.0-cost*cost);
    Hep3Vector momvec(pt*cos(phi), pt*sin(phi), pz);		
    double flightlen(0.0);

    //Create Initial Track
    GVertex gvtx;
    gvtx.setPosition(pos);
    GTrack gtrk;
    PdtPdg::PdgType pdgid = (PdtPdg::PdgType)gconfig.getint("PdtPdg",13);
    PdtEntry* pdt = Pdt::lookup(pdgid);
    double q = pdt->charge();
    gtrk.setVertex(&gvtx);
    gtrk.setPDT(pdt);
    HepLorentzVector p4; p4.setVectM(momvec,pdt->mass());
    gtrk.setP4(p4);

    //Simulate Track through detectors
    PacSimTrack* simtrk = sim.simulateGTrack(&gtrk);

    // global information about simtrk
    fillSimTrkSummary(simtrk,ssum);
    // Timing information
    
    if(disptrack){
       display.reset();
       display.drawGTrack(&gtrk,simtrk->lastHit()->globalFlight(),bfield);
       display.drawSimTrack(simtrk);
       display.drawSimHits(simtrk,0);
     }
        
    const PacPieceTraj* simtraj = simtrk->getTraj();

    //Fill initial parameters
    HepVector simparams(5);
    TrkHelixUtils::helixFromMom(simparams,flightlen,pos,momvec,q,*bfield);

    //Generated Track
    PacHelix gentraj(simparams,flightlen,simtraj->hiRange());

    TrkLineTraj zaxis(HepPoint(0, 0, -10), Hep3Vector(0, 0, 1), 20);
    TrkPoca genpoca(gentraj, 0, zaxis, 10, 1e-12);
    TrkPoca simpoca(*simtraj, 0, zaxis, 10, 1e-12);

      //Store Momentum and Position
    sim_mom_z	= momvec.z();
    sim_mom_mag	= mom;
    sim_inipos_x	= pos.x();
    sim_inipos_y	= pos.y();
    sim_inipos_z	= pos.z();
    
    sim_mom_cost = cost;
    sim_mom_phi = phi;
    sim_mom_pt = pt;

      //Store Parameters
    sim_d0		= simparams(1);
    sim_phi0	= simparams(2);
    sim_omega	= simparams(3);
    sim_z0		= simparams(4);
    sim_tandip	= simparams(5);

      //Store Range of trajectory
    sim_lowrange	= simtraj->lowRange();
    sim_hirange		= simtraj->hiRange();
    sim_poca		= simpoca.flt1();
    sim_doca		= simpoca.doca();

      // Reconstruct the track with KalmanTrack (using the list of hits) 
    TrkRecoTrk* trk = trackreco.makeTrack(simtrk);
    if(trk != 0){
        //Get Reconstructed Track data
      KalInterface kinter;
      trk->attach(kinter,penv.getKalContext()->defaultType());
      const KalRep* kalrep = kinter.kalmanRep();
      const TrkDifPieceTraj& recotraj = kalrep->pieceTraj();
        // find POCA to true production Point
      TrkPoca recpoca(recotraj,0.0,zaxis,10, 1e-12);
      double fltlen(0.0);
      if(recpoca.status().success())
        fltlen = recpoca.flt1();
      TrkExchangePar helix = kalrep->helix(fltlen);
      HepVector recoparams = helix.params();
      HepSymMatrix recocovar = helix.covariance();

        // Get initial momenta of reconstructed track
//      double localflight;
//      const TrkSimpTraj* inithelix = recotraj.localTrajectory(fltlen, localflight);
      TrkPoca recopoca(recotraj, 0, zaxis, 10, 1e-12);

      rec_lowrange	= recotraj.lowRange();
      rec_hirange		= recotraj.hiRange();
      rec_poca		= recopoca.flt1();
      rec_doca		= recopoca.doca();

      rec_d0		= recoparams(1);
      rec_phi0	= recoparams(2);
      rec_omega	= recoparams(3);
      rec_z0		= recoparams(4);
      rec_tandip	= recoparams(5);
      reccov_d0		= sqrt(recocovar.fast(1,1));
      reccov_phi0		= sqrt(recocovar.fast(2,2));
      reccov_omega	= sqrt(recocovar.fast(3,3));
      reccov_z0		= sqrt(recocovar.fast(4,4));
      reccov_tandip	= sqrt(recocovar.fast(5,5));	

      Hep3Vector recoinitmom = kalrep->momentum(0.0);
      BbrVectorErr momerr = kalrep->momentumErr(0.0);
      rec_mom_pt	= recoinitmom.perp();
      rec_mom_z	= recoinitmom.z();
      rec_mom_mag	= recoinitmom.mag();
      Hep3Vector momdir = recoinitmom.unit();
      HepVector momvec(3);
      for(int icor=0;icor<3;icor++)
        momvec[icor] = momdir[icor];
      rec_mom_err = sqrt(momerr.covMatrix().similarity(momvec));
      Hep3Vector ptdir = Hep3Vector(momdir.x(),momdir.y(),0.0).unit();
      for(int icor=0;icor<3;icor++)
        momvec[icor] = ptdir[icor];
      rec_mom_pterr = sqrt(momerr.covMatrix().similarity(momvec));

      //Pull Calculation
      HepVector pull(5);
      trknum = tracknum;
      for(int i = 1; i <= 5; i++) {
        pull(i) = (recoparams(i) - simparams(i)) / sqrt(recocovar(i,i));
      }
      //Store Pull information
      pull_d0		= pull(1);
      pull_phi0	= pull(2);
      pull_omega	= pull(3);
      pull_omegaabs	= (abs(recoparams(3)) - abs(simparams(3))) / sqrt(recocovar(3,3));
      pull_z0		= pull(4);
      pull_tandip	= pull(5);

      //Store fit information
      rec_chisqr	= kalrep->chisq();
      rec_ndof	= kalrep->nDof();
      rec_fitprob = kalrep->chisqConsistency().significanceLevel();
      
      rec_nsvt	= kalrep->hotList()->nSvt();
      rec_ndch	= kalrep->hotList()->nDch();
      rec_nactive = kalrep->hotList()->nActive();
      rec_nhit = kalrep->hotList()->nHit();
      
      // test of position difference between
      if(trajdiff)fillTrajDiff(simtrk,recotraj,tdiff);
      
      if(disptrack)
        display.drawRecTrack(trk);
      
    } else {
// no track: fill with dummy parameters
      rec_lowrange	= -100;
      rec_hirange		= -100;
      rec_poca		= -100;
      rec_doca		= -100;
      rec_d0		= -100;
      rec_phi0	= -100;
      rec_omega	= -100;
      rec_z0		= -100;
      rec_tandip	= -100;
      reccov_d0		= -100;
      reccov_phi0		= -100;
      reccov_omega	= -100;
      reccov_z0		= -100;
      reccov_tandip	= -100;
      rec_mom_pt	= -100;
      rec_mom_z	= -100;
      rec_mom_mag	= -100;
      rec_mom_err	= -100;
      rec_mom_pterr	= -100;
      pull_d0		= -100;
      pull_phi0	= -100;
      pull_omega	= -100;
      pull_omegaabs	= -100;
      pull_z0		= -100;
      pull_tandip	= -100;
      rec_chisqr	= -100;
      rec_fitprob	= -100;
      rec_ndof	= -100;
      rec_nsvt	= -100;
      rec_ndch	= -100;
      rec_nactive	= -100;
      rec_nhit	= -100;
      trknum = -100;
    }
    if(hittuple)fillSimHitInfo(simtrk, sinfo);
    
// cleanup
    delete trk;
    delete simtrk;
    trackT->Fill();
    
    if(disptrack)
      display.fillTrees();
    tracknum++;
  }//int m loop

  //Write track info to file
  trackT->Write();
//  delete trackT;
  file.Close();
  if(disptrack)display.finalize();
  cout << endl;
  return 0;
}

void
fillSimHitInfo(const PacSimTrack* strk, std::vector<PacSimHitInfo>& svec) {
  const std::vector<PacSimHit>& shs = strk->getHitList();
// loop over the simhits
  double radlenint(0);
  double intlenint(0);
  PacSimHitInfo sinfo;
  for(int ish=0;ish<shs.size();ish++){
    sinfo.shi = ish;
    const PacSimHit& sh = shs[ish];
    sinfo.shgloblen = sh.globalFlight();
    sinfo.sheffect = sh.detEffect();
    HepPoint pos = sh.position();
    Hep3Vector sdir = sh.momentumIn().unit();
    sinfo.shx = pos.x();
    sinfo.shy = pos.y();
    sinfo.shz = pos.z();
    sinfo.shmomin = sh.momentumIn().mag();
    sinfo.shmomout = sh.momentumOut().mag();
    if(sh.showerInfo() != 0){
      sinfo.sheinfrac = sh.showerInfo()->fractionIn();
    } else {
      sinfo.sheinfrac = -1.;
    }
    sinfo.shtime =  sh.time();
    const DetIntersection& dinter = sh.detIntersection();
    double pathlen = dinter.pathLength();
    sinfo.shpathlen = pathlen;
    const DetElem* delem = dinter.delem;
    const DetSurfaceElem* selem = dynamic_cast<const DetSurfaceElem*>(delem);
    Hep3Vector snorm(0,0,0);
    if(selem != 0){
      selem->surface()->normTo(pos,snorm);
    }
    sinfo.sdot = snorm.dot(sdir);
    const DetMaterial* mat(0);
    if(delem != 0){
      sinfo.shelemnum = delem->elementNumber();
      sinfo.shtypenum = delem->detectorType()->typeNumber();
      mat = &(delem->material(dinter));
      // NA
      const PacDetElem* pelem = dynamic_cast<const PacDetElem *>(delem);
      if( pelem != 0 && pelem->measurement()!= 0 ) {
        sinfo.shmeastype =  (int)pelem->measurement()->measurementType();
      } else {
        sinfo.shmeastype = -1;
      }
    } else {
      sinfo.shelemnum = -1;
      sinfo.shtypenum = -1;
      sinfo.shmeastype = -1;
    }
    double radlen(0.0);
    double intlen(0.0);
    if(mat != 0){
      radlen = mat->radiationFraction(pathlen);
      intlen = pathlen/mat->intLength();
    }
    sinfo.shradlen = radlen;
    sinfo.shintlen = intlen;
    radlenint += radlen;
    intlenint += intlen;
    sinfo.shradlenint = radlenint;
    sinfo.shintlenint = intlenint;
// look for HOT info
    std::vector<const TrkHitOnTrk*>hots = simHotMap.getHots(&sh);
    sinfo.shnhot = hots.size();
// save one entry/hot, or just 1 entry if there are no hots
    if(hots.size() > 0){
      for(unsigned ihot=0;ihot<hots.size();ihot++){
        const TrkHitOnTrk* hot= hots[ihot];
        sinfo.hview = hot->whatView();
        sinfo.hlay = hot->layerNumber();
        HepPoint trkpoint = hot->trkTraj()->position(hot->fltLen());
        Hep3Vector trkdir = hot->trkTraj()->direction(hot->fltLen());
        HepPoint hitpoint = hot->hitTraj()->position(hot->hitLen());
        Hep3Vector hitdir = hot->hitTraj()->direction(hot->hitLen());
        Hep3Vector pocadir = hitdir.cross(sdir);
        sinfo.resid = trkpoint.distanceTo(hitpoint);
        sinfo.hresid = (hitpoint - pos).dot(pocadir);
        // find planar intersection of track
        double sval = Hep3Vector(pos-trkpoint).dot(snorm)/trkdir.dot(snorm);
        Hep3Vector mdir = hitdir.cross(snorm).unit();
        sinfo.tresid = Hep3Vector(pos - (trkpoint+sval*trkdir)).dot(mdir);
        sinfo.mdot = pocadir.dot(snorm);
        sinfo.herr = hot->hitRms();
        svec.push_back(sinfo);
      }
    } else {
      sinfo.hview = -1;
      sinfo.hlay = -1;
      sinfo.resid = -1.;
      sinfo.tresid = -1.;
      sinfo.hresid = -1.;
      sinfo.mdot = -1.;
      sinfo.herr = -1.;
      svec.push_back(sinfo);
    }
  }
}

void
fillSimTrkSummary(const PacSimTrack* strk, PacSimTrkSummary& ssum) {
// initialize
  ssum = PacSimTrkSummary();
  const std::vector<PacSimHit>& shs = strk->getHitList();
  ssum.dmom = shs.back().momentumOut().mag() - shs.front().momentumIn().mag();
  Hep3Vector ddir = shs.back().momentumOut().unit() - shs.front().momentumIn().unit();
  ssum.ddir = ddir.mag();
  ssum.ddirphi= ddir.phi();
  ssum.pathlen = shs.back().globalFlight() - shs.front().globalFlight();
  ssum.nsimhit = shs.size();
  
  for(int ish=0;ish<shs.size();ish++){
    const PacSimHit& sh = shs[ish];
    if(sh.detEffect() == PacSimHit::normal)
      ssum.nscatter++;
    else if(sh.detEffect() == PacSimHit::brems)
      ssum.nbrems++;
    else if(sh.detEffect() == PacSimHit::shower)
      ssum.nshower++;
    else
      ssum.nother++;
    const DetIntersection& dinter = sh.detIntersection();
    const DetElem* delem = dinter.delem;
    if(delem != 0){
      if(delem->elementName() == "straw")
        ssum.nstraw++;
      else if(delem->elementName() == "gas")
        ssum.ngas++;
      else if(delem->elementName() == "wire")
        ssum.nwire++;
      const PacDetElem* pelem = dynamic_cast<const PacDetElem *>(delem);
      if( pelem != 0 && pelem->measurement()!= 0 ) {
        if(delem->elementName() == "gas")
          ssum.nwiremeas++;
        else
          ssum.npadmeas++;
      }
      const DetMaterial* mat = &(delem->material(dinter));
      if(mat != 0){
        double pathlen = dinter.pathLength();
        ssum.radlenint += mat->radiationFraction(pathlen);
      }
    }
  }
}


void
fillTrajDiff(const PacSimTrack* strk, const TrkDifPieceTraj& ptraj, std::vector<TrajDiff>& tdiff) {
  static const unsigned npts(20);
  const std::vector<PacSimHit>& shs = strk->getHitList();
  const PacPieceTraj* straj = strk->getTraj();
// loop over pairs of measurement simhits
  TrajDiff td;
  for(int ish=0;ish<shs.size();ish++){
    const PacSimHit& sh1 = shs[ish];
    const PacDetElem* pelem1 = dynamic_cast<const PacDetElem *>(sh1.detIntersection().delem);
    if( pelem1 != 0 && pelem1->measurement()!= 0 ) {
// find the next measurement
      for(int jsh=ish+1;jsh<shs.size();jsh++){
        const PacSimHit& sh2 = shs[jsh];
        const PacDetElem* pelem2 = dynamic_cast<const PacDetElem *>(sh2.detIntersection().delem);
        if( pelem2 != 0 && pelem2->measurement()!= 0 ) {
          td.shi = ish;
          td.shj = jsh;
          td.startglen = sh1.globalFlight();
          td.endglen = sh2.globalFlight();
          td.startrho = sh1.position().perp();
          td.startz = sh1.position().z();
          td.endrho = sh2.position().perp();
          td.endz = sh2.position().z();
          td.ddiff = 0.0;
          double step = (td.endglen-td.startglen)/(npts-1);
          Hep3Vector startdir = ptraj.direction(sh1.globalFlight());
          Hep3Vector enddir = ptraj.direction(sh2.globalFlight());
          Hep3Vector average = (startdir + enddir).unit();
          for(unsigned ipt=0;ipt<npts;ipt++){
            double glen = td.startglen+ipt*step;
            HepPoint spt = straj->position(glen);
            TrkPoca tpoca(ptraj, glen, spt, 1e-12);
            HepPoint rpt = ptraj.position(tpoca.flt1());
            Hep3Vector diff = rpt - spt;
            td.ddiff += diff.mag();
          }
          td.ddiff /= npts;
      // integrate the meco field over this range to compute dp    
          Hep3Vector tdp = fieldint->deltaMomentum(straj,td.startglen,td.endglen);
          td.truedp = tdp.mag();
          
          HepPoint tstart = straj->position(td.startglen);
          HepPoint tend = straj->position(td.endglen);
          double pstart=td.startglen;
          TrkPoca tspoca(ptraj, pstart, tstart, 1e-12);
          double pend = td.endglen;
          TrkPoca tepoca(ptraj, pend, tend, 1e-12);
          Hep3Vector rdp = fieldint->deltaMomentum(&ptraj,tspoca.flt1(),tepoca.flt1());
          td.recodp = rdp.mag();
          td.deltadp = (tdp-rdp).mag();
          td.deltadpmag = (tdp-rdp).dot(average);
          tdiff.push_back(td);
          break;
        }
      }
    }
  }
}


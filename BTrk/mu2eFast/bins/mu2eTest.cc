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
#include <TParticle.h>
#include <TParticlePDG.h>
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
#include "PacTrk/PacHitOnTrk.hh"
#include "PacSim/PacSimTrack.hh"
#include "PacSim/PacSimHit.hh"
#include "PacSim/PacWriteParticles.hh"
#include "PacSim/PacShowerInfo.hh"
#include "PacGeom/PacHelix.hh"
#include "PacGeom/PacPieceTraj.hh"
#include "PacGeom/PacMeasurement.hh"
#include "PacDisplay/PacEvtDisplay.hh"
#include "mu2eFast/PacSimHitInfo.rdl"
#include "mu2eFast/TrajDiff.rdl"
#include "mu2eFast/BDiff.rdl"
#include "mu2eFast/PacSimTrkSummary.rdl"
#include "mu2eFast/mu2eDSField.hh"

#include "ProxyDict/Ifd.hh"
#include "ProxyDict/IfdDataProxyUnowned.hh"

#include "AbsEnv/AbsEnv.hh"
#include "BField/BField.hh"
#include "CLHEP/Random/Random.h"
#include "CLHEP/Random/RanecuEngine.h"
#include "CLHEP/Vector/ThreeVector.h"

#include "mu2eFast/Mu2eSimpleInput.hh"
#include "mu2eFast/Mu2eRootInput.hh"
#include "mu2eFast/Mu2eBkgInput.hh"
#include "mu2eFast/Mu2eBeamInput.hh"
#include "mu2eFast/Mu2eTargetInput.hh"

#include "Framework/AppFileName.hh"

using namespace std;

const size_t ncount(5);
struct HitCount{
  Int_t nhit[ncount];
  Int_t nhit_ge[ncount];
  Int_t nstation;
  Int_t ndlayer;
  Int_t nabs;
  Int_t ntar;
  HitCount() { reset();}
  void reset(){
    nstation = ndlayer = 0;
    for(unsigned icount=0;icount<ncount;icount++){
      nhit[icount] = 0;
      nhit_ge[icount] = 0;
    }
    nabs = 0;
    ntar = 0;
  }
  void increment(unsigned nhits){
    nstation++;
    if(nhits<ncount)nhit[nhits]++;
    for(unsigned icount=0;icount <= std::min(nhits,(unsigned)ncount-1);icount++)
      nhit_ge[icount]++;
  }
};

// track selector
class MyTrkSel : public SimTrkSel {
public:
  MyTrkSel(bool invert) : _invert(invert) {}
  void add(const PacSimTrack* strk) { _mytrks.push_back(strk); }
  virtual bool select(const PacSimTrack* strk) const {
    std::vector<const PacSimTrack*>::const_iterator ifnd = std::find(_mytrks.begin(),_mytrks.end(),strk);
    bool found = ifnd != _mytrks.end();
    return _invert ^ found;
  }
  void clear() { _mytrks.clear();}
private:
  std::vector<const PacSimTrack*> _mytrks;
  bool _invert;
};

void fillSimHitInfo(const PacSimTrack* strk, std::vector<PacSimHitInfo>& sinfo);
void fillTrajDiff(const PacSimTrack* strk, const TrkDifPieceTraj& ptraj, std::vector<TrajDiff>& tdiff,
  std::vector<BDiff>& bdiff,PacSimTrkSummary& ssum);
void fillSimTrkSummary(const PacSimTrack* strk, PacSimTrkSummary& ssum);
void countHits(const PacSimTrack* strk,  HitCount& icount);
const PacSimHit* findFirstHit(const PacSimTrack* strk);
void createSim(const PacSimulate& sim,Mu2eEvent& event);


// field integral test stuff
BField* mecofield(0);
BFieldIntegrator* fieldint(0);

PacReconstructTrk* trackreco(0);

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

// get back the field
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
//  bool verbose = gconfig.getbool("verbose",false);
  bool hittuple = gconfig.getbool("hittuple",false);
  bool trajdiff = gconfig.getbool("trajdiff",false);
  bool writeallsim = gconfig.getbool("writeallsim",false);
//  bool calodiff = gconfig.getbool("calodiff",false);

// if requested, setup particle output
  PacWriteParticles writer;

  //Create tuple output
  const char* outputfile = gconfig.get("outputfile", "mu2e_test.root");
  TFile* file = new TFile(outputfile,"RECREATE");
  assert(file != 0);
  //Create Tree to store track info
  TTree* trackT = new TTree("tracks","Tracks");
  //Variables to store track information
  Int_t    itrack;
  Int_t    trknum;
  
  Int_t bkg_ntrks, bkg_nhits;
  Int_t sim_nzero, sim_nsingle, sim_ndouble, sim_ntriple, sim_nquad;
  Int_t sim_nzero_ge, sim_nsingle_ge, sim_ndouble_ge, sim_ntriple_ge, sim_nquad_ge;
  Int_t sim_nstation, sim_ndlayer, sim_nabsorber, sim_ntarget;
  Int_t sim_pdgid;
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
  
  Float_t simt_d0;
  Float_t simt_phi0;
  Float_t simt_omega;
  Float_t simt_z0;
  Float_t simt_tandip;
  Float_t simt_pos_x;
  Float_t simt_pos_y;
  Float_t simt_pos_z;
  Float_t simt_mom_mag;
  Float_t simt_mom_cost;
  Float_t simt_mom_phi;
  Float_t simt_mom_pt;
  
  Int_t rec_pdgid;
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
  Int_t rec_nmerged;
  Int_t rec_nmergeda;
  Int_t rec_nshadowed;
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
  std::vector<BDiff> bdiff;
  PacSimTrkSummary ssum;
  
  //Create TBranch to store track info
  trackT->Branch("itrack",&itrack,"itrack/I");
  trackT->Branch("trknum",&trknum,"trknum/I");
  trackT->Branch("bkg_ntrks",&bkg_ntrks,"bkg_ntrks/I");  
  trackT->Branch("bkg_nhits",&bkg_nhits,"bkg_nhits/I");  
  trackT->Branch("sim_pdgid",&sim_pdgid,"sim_pdgid/I");
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
  
  trackT->Branch("sim_nzero",&sim_nzero,"sim_nzero/I");
  trackT->Branch("sim_nsingle",&sim_nsingle,"sim_nsingle/I");
  trackT->Branch("sim_ndouble",&sim_ndouble,"sim_ndouble/I");
  trackT->Branch("sim_ntriple",&sim_ntriple,"sim_ntriple/I");
  trackT->Branch("sim_nquad",&sim_nquad,"sim_nquad/I");
  trackT->Branch("sim_nzero_ge",&sim_nzero_ge,"sim_nzero_ge/I");
  trackT->Branch("sim_nsingle_ge",&sim_nsingle_ge,"sim_nsingle_ge/I");
  trackT->Branch("sim_ndouble_ge",&sim_ndouble_ge,"sim_ndouble_ge/I");
  trackT->Branch("sim_ntriple_ge",&sim_ntriple_ge,"sim_ntriple_ge/I");
  trackT->Branch("sim_nquad_ge",&sim_nquad_ge,"sim_nquad_ge/I");
  trackT->Branch("sim_nstation",&sim_nstation,"sim_nstation/I");
  trackT->Branch("sim_ndlayer",&sim_ndlayer,"sim_ndlayer/I");
  trackT->Branch("sim_nabsorber",&sim_nabsorber,"sim_nabsorber/I");
  trackT->Branch("sim_ntarget",&sim_ntarget,"sim_ntarget/I");
    
  trackT->Branch("simt_d0",&simt_d0,"simt_d0/F");
  trackT->Branch("simt_phi0",&simt_phi0,"simt_phi0/F");
  trackT->Branch("simt_omega",&simt_omega,"simt_omega/F");
  trackT->Branch("simt_z0",&simt_z0,"simt_z0/F");
  trackT->Branch("simt_tandip",&simt_tandip,"simt_tandip/F");
  trackT->Branch("simt_mom_mag",&simt_mom_mag,"simt_mom_mag/F");
  trackT->Branch("simt_mom_cost",&simt_mom_cost,"simt_mom_cost/F");
  trackT->Branch("simt_mom_phi",&simt_mom_phi,"simt_mom_phi/F");
  trackT->Branch("simt_mom_pt",&simt_mom_pt,"simt_mom_pt/F");
  trackT->Branch("simt_pos_x",&simt_pos_x,"simt_pos_x/F");
  trackT->Branch("simt_pos_y",&simt_pos_y,"simt_pos_y/F");
  trackT->Branch("simt_pos_z",&simt_pos_z,"simt_pos_z/F");

  trackT->Branch("rec_pdgid",&rec_pdgid,"rec_pdgid/I");
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
  trackT->Branch("rec_nmerged",&rec_nmerged,"rec_nmerged/I");
  trackT->Branch("rec_nmergeda",&rec_nmergeda,"rec_nmergeda/I");
  trackT->Branch("rec_nshadowed",&rec_nshadowed,"rec_nshadowed/I");
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
    trackT->Branch("bdiff",&bdiff);
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

// simulation and reconstruction factories

  PacSimulate sim(bfield,detector);
  HepRandomEngine* engine = new RanecuEngine();
  HepRandom::setTheEngine(engine);
  sim.setRandomEngine(engine);
  detector->setRandomEngine(engine);
  trackreco = new PacReconstructTrk(bfield,penv.getKalContext());
  trackreco->setRandomEngine(engine);  

  Mu2eInput* input(0);
  // input specification; check for input file first
  if(gconfig.has("RootFile.inputfile")){
    PacConfig rootconfig = gconfig.getconfig("RootFile.");
    input = new Mu2eRootInput(rootconfig);
  } else if (gconfig.has("TargetInput.nevents")){
    PacConfig targetconfig = gconfig.getconfig("TargetInput.");
    input = new Mu2eTargetInput(targetconfig,detector);
  } else if (gconfig.has("BeamInput.nevents")){
    PacConfig beamconfig = gconfig.getconfig("BeamInput.");
    input = new Mu2eBeamInput(beamconfig,&sim);
  } else if(gconfig.has("SimpleInput.nevents")){
    PacConfig simpleconfig = gconfig.getconfig("SimpleInput.");
    input = new Mu2eSimpleInput(simpleconfig);
  } else {
    std::cerr << "No input specified: aborting" << std::endl;
    return 1;
  }

  Mu2eInput* bkginput(0);
  MyTrkSel* strksel = 0;
  if(gconfig.has("BkgInput.inputfile")){
// we have backgrounds to merge with signal.  Create and configure the input object
    PacConfig bkgconfig = gconfig.getconfig("BkgInput.");
    bkginput = new Mu2eBkgInput(bkgconfig);
    bool fitbkg = bkgconfig.getbool("fitbkg",false);
// reverse logic
    strksel = new MyTrkSel(!fitbkg);
  }
  
  const int printfreq = gconfig.getint("printfreq", 100);
  unsigned nevt(0);
  Mu2eEvent event;
  Mu2eEvent bkgevt;
  bool goodevent;
  while(goodevent = input->nextEvent(event)){
    if(0 == (nevt+1)%printfreq) {
      printf("Count: %i \n",nevt+1);
    }
    nevt++;
// if bkg input exists, merge backgrounds with this event
    if(bkginput != 0 && bkginput->nextEvent(bkgevt)){
      createSim(sim,bkgevt);
// create background selector
      strksel->clear();
// count background
      bkg_ntrks = bkgevt._strks.size();
      bkg_nhits = 0;
      for(unsigned istrk=0;istrk<bkgevt._strks.size();istrk++){
        strksel->add(bkgevt._strks[istrk]);
        const std::vector<PacSimHit>& shs = bkgevt._strks[istrk]->getHitList();
        for(int ish=0;ish<shs.size();ish++){
          const PacSimHit& sh = shs[ish];
          const DetElem* delem = sh.detIntersection().delem;
          const PacDetElem* pelem = dynamic_cast<const PacDetElem *>(delem);
          if( pelem != 0 && pelem->measurementDevices().size()!= 0 )        
            bkg_nhits++;
        }
      }
    } else {
      bkg_ntrks = 0;
      bkg_nhits = 0; 
    }
// simulate signal particle
    createSim(sim,event);
// combine bkg tracks
    event._strks.insert(event._strks.end(),bkgevt._strks.begin(),bkgevt._strks.end());
// create reco tracks
    trackreco->makeTracks(event._strks,strksel);
    for(unsigned istrk=0;istrk<event._strks.size();istrk++){
//      if(verbose)strks[istrk]->print();
// clear vectors
      sinfo.clear();
      tdiff.clear();
//    cdiff.clear();
      const PacSimTrack* simtrk = event._strks[istrk];
      // Find the reconstructed track
      const TrkRecoTrk* trk = trackreco->findTrack(simtrk);
      if(writeallsim || ( trk!= 0 &&  trk->status() != 0 && trk->status()->fitCurrent() ) ){
        const TParticle* tpart = simtrk->getTParticle();
    // global information about simtrk
        fillSimTrkSummary(simtrk,ssum);
    // Timing information

        if(disptrack){
          display.drawParticle(simtrk->getTParticle(),simtrk->lastHit()->globalFlight(),bfield);
          display.drawSimTrack(simtrk);
          display.drawSimHits(simtrk,0);
        }

        const PacPieceTraj* simtraj = simtrk->getTraj();

    //Fill initial parameters
        HepVector simparams(5);
        double flightlen(0.0);
        HepPoint gpos(tpart->Vx(),tpart->Vy(),tpart->Vz());
        Hep3Vector gmom(tpart->Px(),tpart->Py(),tpart->Pz());
        double gcharge = (const_cast<TParticle*>(tpart))->GetPDG()->Charge()/3.0;
        TrkHelixUtils::helixFromMom(simparams,flightlen,gpos,gmom,gcharge,*bfield);

    //Generated Track
        PacHelix gentraj(simparams,flightlen,max(flightlen+10,simtraj->hiRange()));

        TrkLineTraj zaxis(HepPoint(0, 0, -10), Hep3Vector(0, 0, 1), 20);
        TrkPoca genpoca(gentraj, 0, zaxis, 10, 1e-12);
        TrkPoca simpoca(*simtraj, 0, zaxis, 10, 1e-12);

        sim_pdgid = simtrk->pdgId();

      //Store initial Momentum and Position
        const PacSimHit& fhit = simtrk->getHitList()[0];
        Hep3Vector momvec = fhit.momentumIn();
        sim_mom_z	= momvec.z();
        sim_mom_mag	= momvec.mag();
        sim_inipos_x	= fhit.position().x();
        sim_inipos_y	= fhit.position().y();
        sim_inipos_z	= fhit.position().z();

        sim_mom_cost = momvec.cosTheta();
        sim_mom_phi = momvec.phi();
        sim_mom_pt = momvec.perp();

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

    // count the number of measurements per station.
        HitCount hcount;
        countHits(simtrk,hcount);  
        sim_nzero = hcount.nhit[0];
        sim_nzero_ge = hcount.nhit_ge[0];
        sim_nsingle = hcount.nhit[1];
        sim_nsingle_ge = hcount.nhit_ge[1];
        sim_ndouble = hcount.nhit[2];
        sim_ndouble_ge = hcount.nhit_ge[2];
        sim_ntriple = hcount.nhit[3];
        sim_ntriple_ge = hcount.nhit_ge[3];
        sim_nquad = hcount.nhit[4];
        sim_nquad_ge = hcount.nhit_ge[4];
        sim_nstation = hcount.nstation;
        sim_ndlayer = hcount.ndlayer;
        sim_ntarget = hcount.ntar;
        sim_nabsorber = hcount.nabs;
        // look for parameters at the tracker
        
        simt_d0		= -100.0;
        simt_phi0	= -100.0;
        simt_omega	= -100.0;
        simt_z0		= -100.0;
        simt_tandip	= -100.0;
        simt_pos_x	= -100.0;
        simt_pos_y	= -100.0;
        simt_pos_z	= -100.0;
        simt_mom_mag = -100.0;
        simt_mom_cost = -100.0;
        simt_mom_phi = -100.0;
        simt_mom_pt = -100.0;        
        
        const std::vector<PacSimHit>& shs = simtrk->getHitList();
        for(int ish=0;ish<shs.size();ish++){
          const PacSimHit& sh = shs[ish];
          const DetElem* delem = sh.detIntersection().delem;
          if(delem != 0 && delem->elementNumber() == 2){
            simt_pos_x	= sh.position().x();
            simt_pos_y	= sh.position().y();
            simt_pos_z	= sh.position().z();
            simt_mom_mag = sh.momentumIn().mag();
            simt_mom_cost = sh.momentumIn().cosTheta();
            simt_mom_phi = sh.momentumIn().phi();
            simt_mom_pt = sh.momentumIn().perp();

            HepVector simtparams(5);
            double flightlen(0.0);
            TrkHelixUtils::helixFromMom(simtparams,flightlen,
              sh.position(),sh.momentumIn(),gcharge,*bfield);
            //Store Parameters
            simt_d0		= simtparams(1);
            simt_phi0	= simtparams(2);
            simt_omega	= simtparams(3);
            simt_z0		= simtparams(4);
            simt_tandip	= simtparams(5);
            break;
          }
        }
// look for reconstructed tracks
        if(trk != 0 && trk->status() != 0 && trk->status()->fitCurrent() ){
        //Get Reconstructed Track data
          KalInterface kinter;
          trk->attach(kinter,trk->defaultType());
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

          const PdtEntry* pdt = Pdt::lookup(kalrep->particleType(),kalrep->charge());
          if(pdt != 0)
            rec_pdgid = pdt->pdgId();
          else
            rec_pdgid = -1000;
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
          trknum = trk->id();
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
// count merged hits
          rec_nmerged = 0;
          rec_nmergeda = 0;
          rec_nshadowed = 0;
          for (TrkHotList::hot_iterator i = kalrep->hotList()->begin();i!=kalrep->hotList()->end();++i) {
            if (i->usability() >= 11)rec_nmerged++;
            if (i->usability() >= 11 && i->isActive())rec_nmergeda++;
            if (i->usability() <= -10 && !i->isActive())rec_nshadowed++;
          }

      // test of position difference between
          if(trajdiff)fillTrajDiff(simtrk,recotraj,tdiff,bdiff,ssum);

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
          rec_nmerged	= -100;
          rec_nhit	= -100;
          trknum = -100;
        }
    // if we're writing particles, test this one
        if(writer.isActive()){
          const PacSimHit* writehit = findFirstHit(simtrk);
          if(writehit != 0){
            writer.writeParticle(*writehit);
          }
        }
        if(hittuple)fillSimHitInfo(simtrk, sinfo);
        trackT->Fill();
      }
    }
// write display
    if(disptrack){
      display.fillTrees();
      display.reset();
    }
// cleanup this event
    for(unsigned istrk=0;istrk<event._strks.size();istrk++){
      PacSimTrack* simtrk = const_cast<PacSimTrack*>(event._strks[istrk]);
      TrkRecoTrk* trk = const_cast<TrkRecoTrk*>(trackreco->findTrack(simtrk));
      delete simtrk;
      delete trk;
    }
    event._strks.clear();
    bkgevt._strks.clear();
// fill the tree every event
    if(writer.isActive()){
      writer.fillTree();
    }
  }
// close input
  delete input;
//Write track info to file
  TFile* file2  = trackT->GetCurrentFile();
//  trackT->Write();
//  delete trackT;
  file2->Write();
  file2->Close();
  if(disptrack)display.finalize();
  cout << endl;
  return 0;
}

const PacSimHit*
findFirstHit(const PacSimTrack* strk){
  const PacSimHit* retval(0);
  const std::vector<PacSimHit>& shs = strk->getHitList();
  for(int ish=0;ish<shs.size();ish++){
    const PacSimHit& sh = shs[ish];
    const DetElem* delem = sh.detIntersection().delem;
    const PacDetElem* pelem = dynamic_cast<const PacDetElem *>(delem);
    if( pelem != 0 && pelem->measurementDevices().size()!= 0 ) {
      retval = &sh;
      break;
    }
  }
  return retval;
}

void
countHits(const PacSimTrack* strk, HitCount& count) {
// storage for station hit counts;
  std::map<unsigned,unsigned> stations;
  std::vector<unsigned> elements;
// loop over the simhits
  const std::vector<PacSimHit>& shs = strk->getHitList();
  for(int ish=0;ish<shs.size();ish++){
    const PacSimHit& sh = shs[ish];
    const DetElem* delem = sh.detIntersection().delem;
    const PacDetElem* pelem = dynamic_cast<const PacDetElem *>(delem);
    if( pelem != 0 && pelem->measurementDevices().size()!= 0 ) {
// don't count hits with identical layer numbers
      if(std::find(elements.begin(),elements.end(),delem->elementNumber()) == elements.end()){
// extract the 'station number' from the element ID
        div_t idiv = div(delem->elementNumber(),100);
        unsigned station = idiv.rem;
        stations[station]++;
        elements.push_back(delem->elementNumber());
      } else {
// double layer: increment taht
        count.ndlayer++;
      }
    } else if(delem != 0){
      if(delem->elementNumber() == 0)
        count.nabs++;
      else if(delem->elementNumber() ==1 )
        count.ntar++;
    }
  }
// fill the struct
  for(std::map<unsigned,unsigned>::const_iterator istat=stations.begin();istat!= stations.end();istat++){
    count.increment(istat->second);
  }
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
    sinfo.shdel = sqrt(std::max(0.0,1.0 -sh.momentumIn().unit().dot(sh.momentumOut().unit())));
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
      const PacDetElem* pelem = sh.detElem();
      if( pelem != 0 && pelem->measurementDevices().size()!= 0 ) {
        sinfo.shmeastype =  (int)pelem->measurementDevices()[0]->measurementType();
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
    std::vector<const PacHitOnTrk*>hots = trackreco->simHotMap().getHots(&sh);
    sinfo.shnhot = hots.size();
// save one entry/hot, or just 1 entry if there are no hots
    if(hots.size() > 0 && hots[0]->getParentRep() != 0){
// find track trajectory
      const TrkDifTraj& ttraj = hots[0]->getParentRep()->traj();
// multiple hots/simhit not currently handled: fixme!!!
      for(unsigned ihot=0;ihot<hots.size();ihot++){
        const PacHitOnTrk* hot= dynamic_cast<const PacHitOnTrk*>(hots[ihot]);
        if(hot->trkTraj() != 0){
          sinfo.ihot = ihot;
          sinfo.hview = hot->whatView();
          sinfo.hlay = hot->layerNumber();
          sinfo.active = hot->isActive();
          
          HepPoint trkpoint = ttraj.position(hot->fltLen());
          Hep3Vector trkdir = ttraj.direction(hot->fltLen());
          HepPoint hitpoint = hot->hitTraj()->position(hot->hitLen());
          Hep3Vector hitdir = hot->hitTraj()->direction(hot->hitLen());
          Hep3Vector pocadir = hitdir.cross(sdir);
          sinfo.resid = trkpoint.distanceTo(hitpoint);
//        sinfo.hresid = (hitpoint - pos).dot(pocadir);
          sinfo.hresid = hitpoint.distanceTo(pos);
        // find planar intersection of track
          double sval = Hep3Vector(pos-trkpoint).dot(snorm)/trkdir.dot(snorm);
          Hep3Vector mdir = hitdir.cross(snorm).unit();
          sinfo.tresid = Hep3Vector(pos - (trkpoint+sval*trkdir)).dot(mdir);
          sinfo.mdot = pocadir.dot(snorm);
          sinfo.herr = hot->hitRms();
          sinfo.serr = hot->hitInfo()._sres;
          double xresid,xerr;
          hot->resid(xresid,xerr,true);
          sinfo.xresid = xresid;
          sinfo.xerr = xerr;
          svec.push_back(sinfo);
        }
      }
    } else {
      sinfo.hview = -1;
      sinfo.hlay = -1;
      sinfo.active = false;
      sinfo.resid = -1.;
      sinfo.tresid = -1.;
      sinfo.hresid = -1.;
      sinfo.mdot = -1.;
      sinfo.herr = -1.;
      sinfo.serr = -1.;
      sinfo.xresid = -1;
      sinfo.xerr = -1;
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
  ssum.ifirsthit = -1;
  ssum.ilasthit = -1;
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
    double pathlen = dinter.pathLength();
    const DetElem* delem = dinter.delem;
    if(delem != 0){
      double density(-1.);
      const DetMaterial* mat = &(delem->material(dinter));
      if(mat != 0){
        density = mat->density();
        ssum.radlenint += mat->radiationFraction(pathlen);
      }  
      const PacDetElem* pelem = dynamic_cast<const PacDetElem *>(delem);
      assert(pelem != 0);
      if(density < 0.01)
        ssum.ngas++;
      else if(density > 2.8)
        ssum.nwire++;
      else
        ssum.nstraw++;
      if( pelem->measurementDevices().size()!= 0 ) {
        if(delem->elementName().find("cath")!=string::npos)
          ssum.npadmeas++;
        else
          ssum.nwiremeas++;
      }
    }
// record the last hit
    const PacDetElem* pelem = sh.detElem();
    if( pelem != 0 && pelem->measurementDevices().size()!= 0){
      for(std::vector<const PacMeasurement*>::const_iterator imdev = pelem->measurementDevices().begin();
      imdev != pelem->measurementDevices().end();imdev++){
        const PacMeasurement* meas = *imdev;
        if( meas->measurementType() == PacMeasurement::TrkHit && sh.detIntersection().delem->elementNumber()< 5000){
          if(ssum.ifirsthit<0)ssum.ifirsthit=ish;
          ssum.ilasthit = ish;
        }
      }
    }
  }
}


void
fillTrajDiff(const PacSimTrack* strk, const TrkDifPieceTraj& ptraj,
  std::vector<TrajDiff>& tdiff,std::vector<BDiff>& bdiff,PacSimTrkSummary& ssum) {
  static const unsigned npts(20);
  const std::vector<PacSimHit>& shs = strk->getHitList();
  const PacPieceTraj* straj = strk->getTraj();
// loop over pairs of measurement simhits
  TrajDiff td;
  BDiff bd;
  Hep3Vector binttru;
  Hep3Vector bintrec;
  for(int ish=0;ish<shs.size();ish++){
    const PacSimHit& sh1 = shs[ish];
    const PacDetElem* pelem1 = dynamic_cast<const PacDetElem *>(sh1.detIntersection().delem);
    if( pelem1 != 0 && pelem1->measurementDevices().size()!= 0 ) {
// find the next measurement
      for(int jsh=ish+1;jsh<shs.size();jsh++){
        const PacSimHit& sh2 = shs[jsh];
        const PacDetElem* pelem2 = dynamic_cast<const PacDetElem *>(sh2.detIntersection().delem);
        if( pelem2 != 0 && pelem2->measurementDevices().size()!= 0 ) {
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
          for(unsigned ipt=0;ipt<npts;ipt++){
            double glen = td.startglen+ipt*step;
            HepPoint spt = straj->position(glen);
            TrkPoca tpoca(ptraj, glen, spt, 1e-12);
            HepPoint rpt = ptraj.position(tpoca.flt1());
            Hep3Vector diff = rpt - spt;
            Hep3Vector dir = straj->direction(glen);
            Hep3Vector dtrans = diff - dir*diff.dot(dir);
            td.ddiff += dtrans.mag();
            bd.tx = spt.x();
            bd.ty = spt.y();
            bd.tz = spt.z();
            bd.rx = rpt.x();
            bd.ry = rpt.y();
            bd.rz = rpt.z();
            Hep3Vector tb = mecofield->bFieldVect(spt);
            bd.tbx = tb.x();
            bd.tby = tb.y();
            bd.tbz = tb.z();
            Hep3Vector rb = mecofield->bFieldVect(rpt);
            bd.tbx = rb.x();
            bd.tby = rb.y();
            bd.tbz = rb.z();
            bdiff.push_back(bd);
          }
          td.ddiff /= npts;
      // integrate the meco field over this range to compute dp    
          Hep3Vector tdp = fieldint->deltaMomentum(straj,td.startglen,td.endglen);
          td.truedp = tdp.mag();
          td.truedpp = tdp.perp();
          
          HepPoint tstart = straj->position(td.startglen);
          HepPoint tend = straj->position(td.endglen);
          double pstart=td.startglen;
          TrkPoca tspoca(ptraj, pstart, tstart, 1e-12);
          double pend = td.endglen;
          TrkPoca tepoca(ptraj, pend, tend, 1e-12);
          Hep3Vector rdp = fieldint->deltaMomentum(&ptraj,tspoca.flt1(),tepoca.flt1());
          td.recodp = rdp.mag();
          td.recodpp = rdp.perp();
          td.deltadp = (tdp-rdp).mag();
          tdiff.push_back(td);
          
          binttru += tdp;
          bintrec += rdp;
          break;
        }
      }
    }
  }
  ssum.binttru = binttru.mag();
  ssum.bintrec = bintrec.mag();
}

void createSim(const PacSimulate& sim, Mu2eEvent& event) {
  for(std::vector<TParticle*>::const_iterator ipar = event._particles.begin();ipar != event._particles.end();ipar++){
    TParticle* part = *ipar;
    PacSimTrack* simtrk = sim.simulateParticle(part);
    if(simtrk != 0)event._strks.push_back(simtrk);
  }
}


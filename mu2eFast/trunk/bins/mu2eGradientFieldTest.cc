//
//specific to this code, taken from framework

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Geometry/HepPoint.h"
#include "mu2eFast/mu2eGradientField.hh"
#include "Framework/AppFileName.hh"
#include "PacEnv/PacConfig.hh"
#include "PacEnv/PacBuildEnv.hh"
#include "mu2eFast/Mu2eSimpleInput.hh"
#include "PacDisplay/PacEvtDisplay.hh"
#include "PacGeom/PacHelix.hh"
#include "PacGeom/PacPieceTraj.hh"
#include "PacSim/PacSimTrack.hh"
#include "PacSim/PacSimHit.hh"
#include "TrkBase/TrkHelixUtils.hh"
#include "BField/BFieldIntegrator.hh"
#include "Pdt/Pdt.hh"
#include "Pdt/PdtEntry.hh"
#include <iostream>
#include "TParticlePDG.h"

double
  howfar(const BField& field,double mom, double charge, PacHelix* helix,  HepPoint& endpos, Hep3Vector& endmom,double fltlen, double step, double tol);

void
  delta(const BField& field,double mom, double charge, const PacHelix* helix, double fltlen, double step,Hep3Vector& delx);


int main(int argc, char* argv[]) {
  gconfig.verbose(true);
  if(argc <= 1){
    gconfig.parsefile(AppFileName("mu2eFast/testGradient.xml").pathname().c_str());
  }
  gconfig.parseargs(argc, argv);
  PacBuildEnv penv;
  penv.buildCore();

// gradient is 0.3T/m, cenetered on the target, from 2 to 1 tesla  
  gconfig.verbose(true);
  gconfig.parseargs(argc, argv);
  double b0 = gconfig.getfloat("b0",2.0);
  double z0 = gconfig.getfloat("z0",-596);
  double b1 = gconfig.getfloat("b1",1.0);
  double z1 = gconfig.getfloat("z1",-263);
  double step = gconfig.getfloat("step",1.0);
  double tol = gconfig.getfloat("tol",1e-3);

  mu2eGradientField gfield(b0,z0,b1,z1,0.9);

  double bnom = gfield.bFieldNominal();
  std::cout <<" nominal field = " << bnom << std::endl;

  std::cout << "z scan " << std::endl;
  double x = 0.;
  double y = 0.;
  for (int i = -700; i <=-200; i+=10){
    double z = i;
    for(int ir=0;ir<10;ir++){
      x = ir*10.0;
      HepPoint testpoint = HepPoint(x,y,z);
      Hep3Vector bfield = gfield.bFieldVect(testpoint);
      std::cout << "testpoint is " << testpoint << " and BField is " << bfield << std::endl;
    }
  }

  PacEvtDisplay display;
  display.init(gconfig.getcstr("displayfile"),gconfig.getint("displayresolution"));
  display.reset();

// now create a particle and track it
  PacConfig simpleconfig = gconfig.getconfig("SimpleInput.");
  Mu2eSimpleInput* input = new Mu2eSimpleInput(simpleconfig);
  Mu2eEvent event;
  unsigned nevt(0);
  bool goodevent;
  while(goodevent = input->nextEvent(event)){
    printf("Count: %i \n",nevt+1);
    nevt++;
    for(std::vector<TParticle*>::const_iterator ipar = event._particles.begin();ipar != event._particles.end();ipar++){
      TParticle* part = *ipar;
// convert to a helix
      HepVector hpars(5,0);
      HepPoint pos(part->Vx(),part->Vy(),part->Vz());
      Hep3Vector mom(part->Px(),part->Py(),part->Pz());
      double fmax(1000);
      GTrack* gtrk = new GTrack;
      gtrk->setP4(HepLorentzVector( part->Px(),
        part->Py(),
        part->Pz(),
        part->Energy()));
      PdtEntry* pdt = Pdt::lookup((PdtPdg::PdgType)part->GetPdgCode());
      gtrk->setPDT( pdt );      
      PacSimTrack* strk = new PacSimTrack(gtrk);
      double localflt(0);
      double charge = part->GetPDG()->Charge()/3.0;
      TrkHelixUtils::helixFromMom(hpars,localflt,pos,mom,charge,gfield.bFieldVect(pos).z());
      PacHelix* hpart = new PacHelix(hpars,localflt,localflt+fmax);
      DetIntersection dinter;
      dinter.trajet = hpart;
      dinter.pathlen = localflt;
      PacSimHit origin(strk,dinter,pos,mom,mom,PacSimHit::creation);
      strk->addHit(origin);
      strk->addTraj(hpart); 
      while(strk->lastHit()->globalFlight() < fmax){
        Hep3Vector oldmom = mom;
        double dflt = howfar(gfield,mom.mag(),charge,hpart,pos,mom,localflt,step,tol);
        dinter.trajet = hpart;
        dinter.pathlen = localflt+dflt;
        PacSimHit bend(strk,dinter,pos,oldmom,mom,PacSimHit::bend);
        strk->addHit(bend);
        TrkHelixUtils::helixFromMom(hpars,localflt,pos,mom,charge,gfield.bFieldVect(pos).z());
        hpart = new PacHelix(hpars,localflt,localflt+fmax);
        strk->addTraj(hpart);
        localflt += dflt;
      }
      strk->finalize();
      display.drawSimTrack(strk);      
    }
    display.fillTrees();
    display.reset();
  }
  display.finalize();
  return 0;
}

double
howfar(const BField& field,double mom, double charge, PacHelix* helix, HepPoint& endpos, Hep3Vector& endmom, double fltlen, double step, double tol) {
  static const double smax(1000.);
  static const double smalln(10);
  double range[2] = {fltlen,fltlen+step};
  Hep3Vector delx(0,0,0);
  delta(field,mom,charge,helix,fltlen,step,delx);
  double dx = delx.mag();
  if(dx == 0.0) {
    range[1] = fltlen+smax;
  } else {
    double nstep = sqrt(tol/dx);
    if(nstep*step > smax){
      range[1] = fltlen+smax;
// for fractional steps or small steps, return the prediction
    } else if(nstep < smalln){
      range[1] = fltlen + nstep*step;
    } else {
// otherwise, step until we reach tolerance, or the step limit
      unsigned istep(1);
      while(dx < tol && istep < nstep && range[1] < smax){
        double newflt = fltlen+istep*step;
        delta(field,mom,charge,helix,newflt,step,delx);
         dx = delx.mag();
         if(dx < tol)range[1] = newflt;
        istep++;
      }
    }
  }
// reset the range for this traj
  helix->setFlightRange(range);
// compute the momentum change over this flight range, and correct
  BFieldIntegrator ib(field);  
  Hep3Vector dmom = charge*ib.deltaMomentum(helix,range);
  Hep3Vector enddir;
  helix->getInfo(range[1],endpos,enddir);
  endmom = mom*(enddir+dmom/mom).unit();
  return range[1] - range[0];
}


void
delta(const BField& field,double mom, double charge, const PacHelix* helix, double fltlen, double step,Hep3Vector& delx) {
// get the position at this flightlength
  HepPoint pos0;
  Hep3Vector dir0;
  helix->getInfo(fltlen,pos0,dir0);
// get the field at this point
  Hep3Vector bvect0 = field.bFieldVect(pos0);
// Strength at this point is taken as the nominal
  Hep3Vector bnom(0.0,0.0,bvect0.z());
  Hep3Vector delb0 = bvect0-bnom;
// move ahead the estimated step
  HepPoint pos1;
  Hep3Vector dir1;
  helix->getInfo(fltlen+step,pos1,dir1);
  Hep3Vector bvect1 = field.bFieldVect(pos1);
  Hep3Vector delb1 = bvect1-bnom;
// compute deviation
  Hep3Vector d0 = dir0.cross(delb0);
  Hep3Vector d1 = dir1.cross(delb1);
  delx += charge*(0.5*step*step/mom)*(d0 + d1);
}
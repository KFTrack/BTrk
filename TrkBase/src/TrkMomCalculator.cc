//--------------------------------------------------------------------------
// File and Version Information:
// $Id: TrkMomCalculator.cc,v 1.7 2005/05/14 01:07:26 brownd Exp $
// Description:  See the .hh file
//
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author(s): Justin Albert, Steve Schaffner
// Add functions needed for Vertexing: Mario Bondioli, Riccardo Faccini, Eugenio Paoloni.
// Add fast Helix->PX methods: Aaron Roodman
//
//------------------------------------------------------------------------

#include "BaBar/BaBar.hh"
#include "TrkBase/TrkSimpTraj.hh"
#include "TrkBase/TrkMomCalculator.hh"
#include "TrkBase/TrkMomVisitor.hh"
#include "TrkBase/HelixParams.hh"
#include "TrkBase/NeutParams.hh"
#include "BField/BField.hh"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Geometry/HepPoint.h"
#include "difAlgebra/DifPoint.hh"
#include "difAlgebra/DifNumber.hh"
#include "difAlgebra/DifVector.hh"
#include "difAlgebra/DifIndepPar.hh"
#include "BbrGeom/BbrVectorErr.hh"
#include "TrkBase/HelixTraj.hh"
#include "ErrLogger/ErrLog.hh"
using std::endl;

//------------------------------------------------------------------------
TrkMomCalculator::~TrkMomCalculator() {
//------------------------------------------------------------------------
}

//------------------------------------------------------------------------
TrkMomCalculator::TrkMomCalculator() {
//------------------------------------------------------------------------
}

//------------------------------------------------------------------------
double
TrkMomCalculator::ptMom(const TrkSimpTraj& theTraj, const BField&
                        theField, double fltlen) {
//------------------------------------------------------------------------

  TrkMomVisitor theVisitor(theTraj);

  if (theVisitor.helix() != 0 || theVisitor.circle() != 0) {

// treat as curve, calculate Pt accordingly (see calcCurvMom function
// below)
    return calcCurvPtMom(theTraj.direction(fltlen),
                         theTraj.curvature(fltlen), theField);

  } else if (theVisitor.neut() != 0) {

// treat as neutral particle (with the ptot as a parameter)
     double sindip = theTraj.direction(fltlen).z();
     double arg = 1.0-sindip*sindip;
     if (arg < 0.0) arg = 0.0;
     double cosdip = sqrt(arg);
     double ptot = theTraj.parameters()->parameter()[NeutParams::_p];

     return cosdip * ptot; 

  } else {

// particle must be a plain line--no way to calculate momentum
    return 999999.99;

  }
  
}

//------------------------------------------------------------------------
Hep3Vector
TrkMomCalculator::vecMom(const TrkSimpTraj& theTraj, const BField&
                         theField, double fltlen) {
//------------------------------------------------------------------------

  TrkMomVisitor theVisitor(theTraj);

  if (theVisitor.helix() != 0 || theVisitor.circle() != 0) {

// treat as curve, calculate VecMom accordingly (see calcCurvMom function
// below)
    return calcCurvVecMom(theTraj.direction(fltlen),
                          theTraj.curvature(fltlen), theField);

  } else if (theVisitor.neut() != 0) {

// treat as neutral particle (with the pt as a parameter)
     Hep3Vector theMom = theTraj.direction(fltlen);
     theMom.setMag(theTraj.parameters()->parameter()[NeutParams::_p]);
     return theMom;

  } else {

// particle must be a plain line--no way to calculate momentum
    return Hep3Vector(999, 999, 999);

  }
}

//------------------------------------------------------------------------
BbrVectorErr
TrkMomCalculator::errMom(const TrkSimpTraj& theTraj, const BField&
                         theField, double fltlen) {
//------------------------------------------------------------------------

  TrkMomVisitor theVisitor(theTraj);

  if (theVisitor.helix() != 0 || theVisitor.circle() != 0) {

// treat as curve, calculate ErrMom accordingly (see calcCurvMom function
// below)
    return calcCurvErrMom(theTraj, theField, fltlen);

  } else if (theVisitor.neut() != 0) {

// treat as neutral particle, same as curve in this case
     return calcNeutErrMom(theTraj, theField, fltlen);

  } else {

// particle must be a plain line--no way to calculate momentum or err
// The matrix is initialized to zero (see BbrError constructor)
    BbrError theErr(3); 
    return BbrVectorErr(Hep3Vector(999, 999, 999), theErr);
  }
}

//------------------------------------------------------------------------
int
TrkMomCalculator::charge(const TrkSimpTraj& theTraj, const BField&
                         theField, double fltlen) {
//------------------------------------------------------------------------

  TrkMomVisitor theVisitor(theTraj);

  if (theVisitor.helix() != 0 || theVisitor.circle() != 0) {

// treat as curve, calculate Pt accordingly (see calcCurvMom function
// below)

    bool plus = false;
    HelixTraj::ParIndex parInd = HelixTraj::omegaIndex;
    int index = parInd;
    plus = (theField.bFieldNominal() < 0.0 &&
            theTraj.parameters()->parameter()[index] > 0.0)
            || 
           (theField.bFieldNominal() > 0.0 &&
            theTraj.parameters()->parameter()[index] < 
            0.0);
    return ( plus ? 1 : -1 );

//    return calcCurvCharge(
//                          theTraj.direction(fltlen),
//                          theTraj.curvature(fltlen), theField);

  } else if (theVisitor.neut() != 0) {

// treat as neutral particle, so charge is zero
     return 0;

  } else {

// particle must be a plain line--take charge as zero
    return 0;
  }
}

//------------------------------------------------------------------------
HepMatrix
TrkMomCalculator::posmomCov(const TrkSimpTraj& theTraj,const BField& theField,
			    double fltlen) {
//------------------------------------------------------------------------

  TrkMomVisitor theVisitor(theTraj);

  if (theVisitor.helix() != 0 || theVisitor.circle() != 0) {

// treat as curve, calculate ErrMom accordingly (see calcCurvMom function
// below)
    return calcCurvPosmomCov(theTraj, theField, fltlen);

  } else if (theVisitor.neut() != 0) {

// treat as neutral particle, same as curve in this case
    return calcNeutPosmomCov(theTraj, theField, fltlen);

  } else {

// particle must be a plain line--no way to calculate momentum or err
// The matrix is initialized to zero (see BbrError constructor)
    return HepMatrix(3,3,0);
  }
}

//------------------------------------------------------------------------
void
TrkMomCalculator::getAllCovs(const TrkSimpTraj& theTraj,
			     const BField& theField,
			     double fltlen,
			     HepSymMatrix& xxCov,
			     HepSymMatrix& ppCov,
			     HepMatrix&    xpCov) {
//------------------------------------------------------------------------

  TrkMomVisitor theVisitor(theTraj);

  if (theVisitor.helix() != 0) {

    // fast inline calculation...
    calcCurvAllCovs(theTraj,theField,fltlen,xxCov,ppCov,xpCov);

  } else if (theVisitor.circle() != 0){

    // treat as curve, calculate ErrMom accordingly (see calcCurvMom function
    // below)
    calcCurvAllCovsOLD(theTraj,theField,fltlen,xxCov,ppCov,xpCov);

  } else if (theVisitor.neut() != 0) {
    
    // treat as neutral particle, same as curve in this case
    calcNeutAllCovs(theTraj,theField,fltlen,xxCov,ppCov,xpCov);

  } else {

    // particle must be a plain line--no way to calculate momentum or err
    // The matrix is initialized to zero (see BbrError constructor)
    int i,j;
    for(i=0;i<3;i++)
      for(j=i;j<3;j++)
	{
	  xxCov[i][j]=0;
	  ppCov[i][j]=0;
	  xpCov[i][j]=0;
	  xpCov[j][i]=0;
	}
  }

}

//------------------------------------------------------------------------
void
TrkMomCalculator::getAllWeights(const TrkSimpTraj& theTraj,
				const BField& theField,
				double fltlen,
				HepVector& pos,
				HepVector& mom,
				HepSymMatrix& xxWeight,
				HepSymMatrix& ppWeight,
				HepMatrix&    xpWeight) {
//------------------------------------------------------------------------

  TrkMomVisitor theVisitor(theTraj);

  if (theVisitor.helix() != 0){

    // treat as curve, calculate ErrMom accordingly (see calcCurvMom function
    // below)
    calcCurvAllWeights(theTraj,theField,fltlen,
		       pos,mom,xxWeight,ppWeight,xpWeight);

  } else if (theVisitor.circle() != 0) {

    calcCurvAllWeightsOLD(theTraj,theField,fltlen,
		       pos,mom,xxWeight,ppWeight,xpWeight);

  } else if (theVisitor.neut() != 0) {
    
    // treat as neutral particle, same as curve in this case
    calcNeutAllWeights(theTraj,theField,fltlen,
		       pos,mom,xxWeight,ppWeight,xpWeight);

  } else {

    // particle must be a plain line--no way to calculate momentum or err
    // temporary: initialize everything to 0 
    int i,j;
    for(i=0;i<3;i++)
      for(j=i;j<3;j++) {
	xxWeight[i][j]=0;
	ppWeight[i][j]=0;
	xpWeight[i][j]=0;
	xpWeight[j][i]=0;
      }

    for(i=0;i<3;i++) {
      pos[i]=0;
      mom[i]=0;
    }
  }
}

//------------------------------------------------------------------------
double
TrkMomCalculator::calcCurvPtMom(const Hep3Vector& direction,
                                double curvature,
                                const BField& theField) {
//------------------------------------------------------------------------
  Hep3Vector theMomVec = calcCurvVecMom(direction, curvature,
                                        theField);
  return theMomVec.perp();

}

//------------------------------------------------------------------------
Hep3Vector
TrkMomCalculator::calcCurvVecMom(const Hep3Vector& direction,
                                 double curvature,
                                 const BField& theField) {
//------------------------------------------------------------------------
  double sindip = direction.z();
  double arg = 1.0-sindip*sindip;
  if (arg < 0.0) arg = 0.0;
  double cosdip = sqrt(arg);

  double momMag =
           fabs( BField::mmTeslaToMeVc*theField.bFieldNominal()*
		cosdip/curvature );

  Hep3Vector momVec = direction;
  momVec.setMag(momMag);

  return momVec;

}

DifNumber
TrkMomCalculator::momMag(const HelixTraj& theTraj,
			 const BField& theField) {
  DifNumber d0Df(theTraj.d0(), HelixTraj::d0Index+1, HelixTraj::NHLXPRM);
  DifNumber phi0Df(theTraj.phi0(), HelixTraj::phi0Index+1, HelixTraj::NHLXPRM);
  DifNumber omegaDf(theTraj.omega(), HelixTraj::omegaIndex+1, HelixTraj::NHLXPRM);
  DifNumber z0Df(theTraj.z0(), HelixTraj::z0Index+1, HelixTraj::NHLXPRM);
  DifNumber tanDipDf(theTraj.tanDip(), HelixTraj::tanDipIndex+1, HelixTraj::NHLXPRM);
  phi0Df.setIndepPar( theTraj.parameters() );
  d0Df.setIndepPar( theTraj.parameters() );
  z0Df.setIndepPar( theTraj.parameters() );
  tanDipDf.setIndepPar( theTraj.parameters() );
  omegaDf.setIndepPar( theTraj.parameters() );

  DifNumber dipDf = atan(tanDipDf);
  DifNumber cosdip = cos(dipDf);
  DifNumber momMag =
    BField::mmTeslaToMeVc*theField.bFieldNominal() / (cosdip * omegaDf);
  momMag.absolute();
  return momMag;
}

//------------------------------------------------------------------------
BbrVectorErr
TrkMomCalculator::calcCurvErrMom(const TrkSimpTraj& theTraj,
                                 const BField& theField,
                                 double fltlen) {
//------------------------------------------------------------------------

  DifPoint  PosDif;
  DifVector DirDif;
  DifVector delDirDif;
  DifVector MomDif(0., 0., 0.);

  theTraj.getDFInfo(fltlen, PosDif, DirDif, delDirDif);
  if (delDirDif.length() == 0.) {   
  }
  else {
    DifNumber sindip=DirDif.z;
    DifNumber arg = 1.0-sindip*sindip;
    if (arg.number() < 0.0) {arg.setNumber(0.0);}
    DifNumber cosdip = sqrt(arg);
    
    DifNumber momMag =
      BField::mmTeslaToMeVc*theField.bFieldNominal()*cosdip /
      delDirDif.length();
    momMag.absolute();

    MomDif = DirDif * momMag;

    if (ErrLogging(debugging) && 0) {
      HepMatrix e1 = DirDif.errorMatrix(MomDif.x.indepPar()->covariance());
      double e2 = momMag.error(MomDif.x.indepPar()->covariance());
      HepMatrix e3 = MomDif.errorMatrix(MomDif.x.indepPar()->covariance());
      const HepMatrix& c1 = MomDif.x.indepPar()->covariance();
      
      ErrMsg(debugging) << endl << endl << "Start" << endl
      << "param cov: " << endl
      << c1(1,1) << endl
      << c1(2,1) << "  " << c1(2,2) << endl
      << c1(3,1) << "  " << c1(3,2) << "  " << c1(3,3) << endl
      << c1(4,1) << "  " << c1(4,2) << "  " << c1(4,3) 
      << "  " << c1(4,4) << endl
      << c1(5,1) << "  " << c1(5,2) << "  " << c1(5,3) 
      << "  " << c1(5,4) << "  " << c1(5,5) << endl
      
      << "Dir " << e1.num_row() << " " << e1.num_col() << endl
      << "output:" << endl
      << e1(1,1) << endl
      << e1(2,1) << "  " << e1(2,2) << endl
      << e1(3,1) << "  " << e1(3,2) << "  " << e1(3,3) << endl
      << "deriv: " << endl
      << DirDif.x.derivatives() << endl
      << endl
      
      << "Mag " << endl
      << "output:" << endl
      << e2 << endl
      << "deriv: " << endl
      << momMag.derivatives() << endl
      << endl
      
      << "Momdif " << e3.num_row() << " " << e3.num_col() << endl
      << "output:" << endl
      << e3(1,1) << endl
      << e3(2,1) << "  " << e3(2,2) << endl
      << e3(3,1) << "  " << e3(3,2) << "  " << e3(3,3) << endl
      << "deriv: " << endl
      << MomDif.x.derivatives() << endl
      << endl
      
      << "End" << endl << endmsg; 
    }
  }
  BbrError  symErr(MomDif.errorMatrix(
                              MomDif.x.indepPar()->covariance()));
  return BbrVectorErr(Hep3Vector(MomDif.x.number(), MomDif.y.number(),
                    MomDif.z.number()), symErr);
}

//------------------------------------------------------------------------
BbrVectorErr
TrkMomCalculator::calcNeutErrMom(const TrkSimpTraj& theTraj,
                                 const BField& theField,
                                 double fltlen) {
//------------------------------------------------------------------------

  DifPoint  posDif;
  DifVector dirDif;
  DifVector delDirDif;
 
  theTraj.getDFInfo(fltlen, posDif, dirDif, delDirDif);

// set the momentum's direction, and then its magnitude
  DifVector momDif = dirDif;

// The + 1 is necessary because DifIndepPar starts counting at 1, whereas
// the enum for NeutParams starts at 0.  Stupid, ain't it?
  momDif *= theTraj.parameters()->difPar(NeutParams::_p + 1);

// now create an error matrix
  BbrError  symErr(momDif.errorMatrix(
                              momDif.x.indepPar()->covariance()));

// get the result in the correct form
  return BbrVectorErr(Hep3Vector(momDif.x.number(), momDif.y.number(),
                    momDif.z.number()), symErr);
}

// The following functions may be used at a later date (if and when we
// ever want to drop the assumption that all recon'ed charges are 0, +1,
// or -1):

//------------------------------------------------------------------------
int
TrkMomCalculator::calcCurvCharge(const Hep3Vector& direction,
                                 double curvature,
                                 const BField& theField) {
//------------------------------------------------------------------------

   Hep3Vector momVec = 
                calcCurvVecMom(direction, curvature, theField);

   if (theField.bFieldNominal() > 0.) {
     return -nearestInt(momVec.mag() * curvature / theField.bFieldNominal());
   } else {
     return nearestInt(momVec.mag() * curvature / theField.bFieldNominal());
   }
}                        

//------------------------------------------------------------------------
int
TrkMomCalculator::nearestInt(double floatPt) {
//------------------------------------------------------------------------

  if (floatPt > 0.) {
    return (int)(floatPt+0.5);
  } else {
    return (int)(floatPt-0.5);
  }
}

/*
 * These function were added on 7/18/98 to accomplish 
 * vertexing interface (M.Bondioli)
 */

//------------------------------------------------------------------------
HepMatrix
TrkMomCalculator::calcCurvPosmomCov(const TrkSimpTraj& theTraj,
				    const BField& theField,
				    double fltlen) {
//------------------------------------------------------------------------

  DifPoint  PosDif;
  DifVector DirDif;
  DifVector delDirDif;
  DifVector MomDif(0., 0., 0.);

  theTraj.getDFInfo(fltlen, PosDif, DirDif, delDirDif);
  if (delDirDif.length() == 0.) {   
  }
  else {
    DifNumber sindip=DirDif.z;
    DifNumber arg = 1.0-sindip*sindip;
    if (arg.number() < 0.0) {arg.setNumber(0.0);}
    DifNumber cosdip = sqrt(arg);


    DifNumber momMag =
      BField::mmTeslaToMeVc*theField.bFieldNominal()*cosdip /
      delDirDif.length();
    momMag.absolute();

    MomDif = DirDif * momMag;

  }

  // computes the correlation among position and momentum
  HepMatrix xpCov(3,3);
  const HepSymMatrix& nnCov=MomDif.x.indepPar()->covariance();
  xpCov(1,1)=correlation(PosDif.x,MomDif.x,nnCov);
  xpCov(1,2)=correlation(PosDif.x,MomDif.y,nnCov);
  xpCov(1,3)=correlation(PosDif.x,MomDif.z,nnCov);
  xpCov(2,1)=correlation(PosDif.y,MomDif.x,nnCov);
  xpCov(2,2)=correlation(PosDif.y,MomDif.y,nnCov);
  xpCov(2,3)=correlation(PosDif.y,MomDif.z,nnCov);
  xpCov(3,1)=correlation(PosDif.z,MomDif.x,nnCov);
  xpCov(3,2)=correlation(PosDif.z,MomDif.y,nnCov);
  xpCov(3,3)=correlation(PosDif.z,MomDif.z,nnCov);

  return xpCov;
}

//------------------------------------------------------------------------
HepMatrix
TrkMomCalculator::calcNeutPosmomCov(const TrkSimpTraj& theTraj,
				    const BField& theField,
				    double fltlen) {
//------------------------------------------------------------------------

  DifPoint  PosDif;
  DifVector DirDif;
  DifVector delDirDif;
 
  theTraj.getDFInfo(fltlen, PosDif, DirDif, delDirDif);

  // set the momentum's direction, and then its magnitude
  DifVector MomDif = DirDif;

  // The + 1 is necessary because DifIndepPar starts counting at 1, whereas
  // the enum for NeutParams starts at 0.  Stupid, ain't it?
  MomDif *= theTraj.parameters()->difPar(NeutParams::_p + 1);

  // computes the correlation among position and momentum
  HepMatrix xpCov(3,3);
  const HepSymMatrix& nnCov=MomDif.x.indepPar()->covariance();
  xpCov(1,1)=correlation(PosDif.x,MomDif.x,nnCov);
  xpCov(1,2)=correlation(PosDif.x,MomDif.y,nnCov);
  xpCov(1,3)=correlation(PosDif.x,MomDif.z,nnCov);
  xpCov(2,1)=correlation(PosDif.y,MomDif.x,nnCov);
  xpCov(2,2)=correlation(PosDif.y,MomDif.y,nnCov);
  xpCov(2,3)=correlation(PosDif.y,MomDif.z,nnCov);
  xpCov(3,1)=correlation(PosDif.z,MomDif.x,nnCov);
  xpCov(3,2)=correlation(PosDif.z,MomDif.y,nnCov);
  xpCov(3,3)=correlation(PosDif.z,MomDif.z,nnCov);

  return HepMatrix(3,3,0);
}

bool TrkMomCalculator::weightToCov(const HepSymMatrix& inXX,const HepSymMatrix& inPP,const HepMatrix& inXP,
				      HepSymMatrix& outXX,HepSymMatrix& outPP,HepMatrix& outXP){
  assert(inXX.num_row()==outXX.num_row());
  assert(inPP.num_row()==outPP.num_row());
  assert(inXP.num_row()==outXP.num_row());
  assert(inXP.num_col()==outXP.num_col());
  assert(inXX.num_row()==inXP.num_row());
  assert(inPP.num_row()==inXP.num_col());
  int status;
  HepSymMatrix aInv = inXX.inverse(status);
  if(status)return false;
  HepSymMatrix beta = inPP-aInv.similarityT(inXP);
  outPP = beta.inverse(status);
  if(status)return false;
  outXP = -aInv*inXP*outPP;
  HepMatrix alpha(aInv-aInv*inXP*outXP.T());
  outXX.assign(alpha);
  return true;
}

//------------------------------------------------------------------------
void
TrkMomCalculator::calcCurvAllCovs(const TrkSimpTraj& theTraj,
				  const BField& theField,
				  double fltlen,
				  HepSymMatrix& xxCov,
				  HepSymMatrix& ppCov,
				  HepMatrix&    xpCov) {
//------------------------------------------------------------------------

  const HepVector& v = theTraj.parameters()->parameter();
  const HepSymMatrix& m = theTraj.parameters()->covariance();

  // fast inlined conversion from Helix to PX representation

  // track parameters
  double s = v[HelixParams::tanDipIndex];
  double omega = v[HelixParams::omegaIndex];
  double d0 = v[HelixParams::d0Index];
  //double z0 = v[HelixParams::z0Index];
  double phi0 = v[HelixParams::phi0Index];
  double cosDip = 1/sqrt(1.0+s*s);
  double l = fltlen*cosDip ;

  // calculate some sin,cos etc..
  double dlds = -fltlen*s*(cosDip*cosDip*cosDip) ;
  double phi = phi0 + omega*l;
  double sinphi0 = sin(phi0);
  double cosphi0 = cos(phi0);
  double sinphi = sin(phi);
  double cosphi = cos(phi);
  double pt = fabs( BField::mmTeslaToMeVc * theField.bFieldNominal() / omega );
  double r = 1.0/omega;
  
  // Calculate derivatives for Jacobian matrix
  double d_x_d0 = -sinphi0;
  double d_x_phi0 = r*cosphi - (r+d0)*cosphi0;
  double d_x_omega = -r*r*sinphi + r*r*sinphi0 + l*r*cosphi;
  double d_x_tanDip = cosphi*dlds;
  double d_y_d0 = cosphi0;
  double d_y_phi0 = r*sinphi - (r+d0)*sinphi0;
  double d_y_omega = r*r*cosphi - r*r*cosphi0 + l*r*sinphi;
  double d_y_tanDip = sinphi*dlds;
  double d_z_z0 = 1.0;
  double d_z_tanDip = l + dlds*s;
  double d_px_phi0 = -pt*sinphi;
  double d_px_omega = -pt*cosphi/omega - pt*l*sinphi;
  double d_px_tanDip = -pt*omega*sinphi*dlds;
  double d_py_phi0 = pt*cosphi;
  double d_py_omega = -pt*sinphi/omega + pt*l*cosphi;
  double d_py_tanDip = pt*omega*cosphi*dlds;
  double d_pz_omega = -pt*s/omega;
  double d_pz_tanDip = pt;

  // Fill temporary variables for m
  double m_d0_d0 =  m[HelixParams::d0Index][HelixParams::d0Index];
  double m_phi0_d0 =  m[HelixParams::phi0Index][HelixParams::d0Index];
  double m_phi0_phi0 =  m[HelixParams::phi0Index][HelixParams::phi0Index];
  double m_omega_d0 =  m[HelixParams::omegaIndex][HelixParams::d0Index];
  double m_omega_phi0 =  m[HelixParams::omegaIndex][HelixParams::phi0Index];
  double m_omega_omega =  m[HelixParams::omegaIndex][HelixParams::omegaIndex];
  double m_z0_d0 =  m[HelixParams::z0Index][HelixParams::d0Index];
  double m_z0_phi0 =  m[HelixParams::z0Index][HelixParams::phi0Index];
  double m_z0_omega =  m[HelixParams::z0Index][HelixParams::omegaIndex];
  double m_z0_z0 =  m[HelixParams::z0Index][HelixParams::z0Index];
  double m_tanDip_d0 =  m[HelixParams::tanDipIndex][HelixParams::d0Index];
  double m_tanDip_phi0 =  m[HelixParams::tanDipIndex][HelixParams::phi0Index];
  double m_tanDip_omega =  m[HelixParams::tanDipIndex][HelixParams::omegaIndex];
  double m_tanDip_z0 =  m[HelixParams::tanDipIndex][HelixParams::z0Index];
  double m_tanDip_tanDip =  m[HelixParams::tanDipIndex][HelixParams::tanDipIndex];

  // Fill xxCov,xpCov,ppCov - nb. this code generated by script writecov.py
  xxCov(1,1) = 
    d_x_d0* (  d_x_d0*m_d0_d0 + d_x_phi0*m_phi0_d0 + d_x_omega*m_omega_d0 + d_x_tanDip*m_tanDip_d0  )   + 
    d_x_phi0* (  d_x_d0*m_phi0_d0 + d_x_phi0*m_phi0_phi0 + d_x_omega*m_omega_phi0 + d_x_tanDip*m_tanDip_phi0  )   + 
    d_x_omega* (  d_x_d0*m_omega_d0 + d_x_phi0*m_omega_phi0 + d_x_omega*m_omega_omega + d_x_tanDip*m_tanDip_omega  )   + 
    d_x_tanDip* (  d_x_d0*m_tanDip_d0 + d_x_phi0*m_tanDip_phi0 + d_x_omega*m_tanDip_omega + d_x_tanDip*m_tanDip_tanDip  )   ; 
  xxCov(2,1) = 
    d_y_d0* (  d_x_d0*m_d0_d0 + d_x_phi0*m_phi0_d0 + d_x_omega*m_omega_d0 + d_x_tanDip*m_tanDip_d0  )   + 
    d_y_phi0* (  d_x_d0*m_phi0_d0 + d_x_phi0*m_phi0_phi0 + d_x_omega*m_omega_phi0 + d_x_tanDip*m_tanDip_phi0  )   + 
    d_y_omega* (  d_x_d0*m_omega_d0 + d_x_phi0*m_omega_phi0 + d_x_omega*m_omega_omega + d_x_tanDip*m_tanDip_omega  )   + 
    d_y_tanDip* (  d_x_d0*m_tanDip_d0 + d_x_phi0*m_tanDip_phi0 + d_x_omega*m_tanDip_omega + d_x_tanDip*m_tanDip_tanDip  )   ; 
  xxCov(2,2) = 
    d_y_d0* (  d_y_d0*m_d0_d0 + d_y_phi0*m_phi0_d0 + d_y_omega*m_omega_d0 + d_y_tanDip*m_tanDip_d0  )   + 
    d_y_phi0* (  d_y_d0*m_phi0_d0 + d_y_phi0*m_phi0_phi0 + d_y_omega*m_omega_phi0 + d_y_tanDip*m_tanDip_phi0  )   + 
    d_y_omega* (  d_y_d0*m_omega_d0 + d_y_phi0*m_omega_phi0 + d_y_omega*m_omega_omega + d_y_tanDip*m_tanDip_omega  )   + 
    d_y_tanDip* (  d_y_d0*m_tanDip_d0 + d_y_phi0*m_tanDip_phi0 + d_y_omega*m_tanDip_omega + d_y_tanDip*m_tanDip_tanDip  )   ; 
  xxCov(3,1) = 
    d_z_z0* (  d_x_d0*m_z0_d0 + d_x_phi0*m_z0_phi0 + d_x_omega*m_z0_omega + d_x_tanDip*m_tanDip_z0  )   + 
    d_z_tanDip* (  d_x_d0*m_tanDip_d0 + d_x_phi0*m_tanDip_phi0 + d_x_omega*m_tanDip_omega + d_x_tanDip*m_tanDip_tanDip  )   ; 
  xxCov(3,2) = 
    d_z_z0* (  d_y_d0*m_z0_d0 + d_y_phi0*m_z0_phi0 + d_y_omega*m_z0_omega + d_y_tanDip*m_tanDip_z0  )   + 
    d_z_tanDip* (  d_y_d0*m_tanDip_d0 + d_y_phi0*m_tanDip_phi0 + d_y_omega*m_tanDip_omega + d_y_tanDip*m_tanDip_tanDip  )   ; 
  xxCov(3,3) = 
    d_z_z0* (  d_z_z0*m_z0_z0 + d_z_tanDip*m_tanDip_z0  )   + 
    d_z_tanDip* (  d_z_z0*m_tanDip_z0 + d_z_tanDip*m_tanDip_tanDip  )   ; 

  xpCov(1,1) = 
    d_px_phi0* (  d_x_d0*m_phi0_d0 + d_x_phi0*m_phi0_phi0 + d_x_omega*m_omega_phi0 + d_x_tanDip*m_tanDip_phi0  )   + 
    d_px_omega* (  d_x_d0*m_omega_d0 + d_x_phi0*m_omega_phi0 + d_x_omega*m_omega_omega + d_x_tanDip*m_tanDip_omega  )   + 
    d_px_tanDip* (  d_x_d0*m_tanDip_d0 + d_x_phi0*m_tanDip_phi0 + d_x_omega*m_tanDip_omega + d_x_tanDip*m_tanDip_tanDip  )   ; 
  xpCov(2,1) = 
    d_px_phi0* (  d_y_d0*m_phi0_d0 + d_y_phi0*m_phi0_phi0 + d_y_omega*m_omega_phi0 + d_y_tanDip*m_tanDip_phi0  )   + 
    d_px_omega* (  d_y_d0*m_omega_d0 + d_y_phi0*m_omega_phi0 + d_y_omega*m_omega_omega + d_y_tanDip*m_tanDip_omega  )   + 
    d_px_tanDip* (  d_y_d0*m_tanDip_d0 + d_y_phi0*m_tanDip_phi0 + d_y_omega*m_tanDip_omega + d_y_tanDip*m_tanDip_tanDip  )   ; 
  xpCov(3,1) = 
    d_px_phi0* (  d_z_z0*m_z0_phi0 + d_z_tanDip*m_tanDip_phi0  )   + 
    d_px_omega* (  d_z_z0*m_z0_omega + d_z_tanDip*m_tanDip_omega  )   + 
    d_px_tanDip* (  d_z_z0*m_tanDip_z0 + d_z_tanDip*m_tanDip_tanDip  )   ; 
  xpCov(1,2) = 
    d_py_phi0* (  d_x_d0*m_phi0_d0 + d_x_phi0*m_phi0_phi0 + d_x_omega*m_omega_phi0 + d_x_tanDip*m_tanDip_phi0  )   + 
    d_py_omega* (  d_x_d0*m_omega_d0 + d_x_phi0*m_omega_phi0 + d_x_omega*m_omega_omega + d_x_tanDip*m_tanDip_omega  )   + 
    d_py_tanDip* (  d_x_d0*m_tanDip_d0 + d_x_phi0*m_tanDip_phi0 + d_x_omega*m_tanDip_omega + d_x_tanDip*m_tanDip_tanDip  )   ; 
  xpCov(2,2) = 
    d_py_phi0* (  d_y_d0*m_phi0_d0 + d_y_phi0*m_phi0_phi0 + d_y_omega*m_omega_phi0 + d_y_tanDip*m_tanDip_phi0  )   + 
    d_py_omega* (  d_y_d0*m_omega_d0 + d_y_phi0*m_omega_phi0 + d_y_omega*m_omega_omega + d_y_tanDip*m_tanDip_omega  )   + 
    d_py_tanDip* (  d_y_d0*m_tanDip_d0 + d_y_phi0*m_tanDip_phi0 + d_y_omega*m_tanDip_omega + d_y_tanDip*m_tanDip_tanDip  )   ; 
  xpCov(3,2) = 
    d_py_phi0* (  d_z_z0*m_z0_phi0 + d_z_tanDip*m_tanDip_phi0  )   + 
    d_py_omega* (  d_z_z0*m_z0_omega + d_z_tanDip*m_tanDip_omega  )   + 
    d_py_tanDip* (  d_z_z0*m_tanDip_z0 + d_z_tanDip*m_tanDip_tanDip  )   ; 
  xpCov(1,3) = 
    d_pz_omega* (  d_x_d0*m_omega_d0 + d_x_phi0*m_omega_phi0 + d_x_omega*m_omega_omega + d_x_tanDip*m_tanDip_omega  )   + 
    d_pz_tanDip* (  d_x_d0*m_tanDip_d0 + d_x_phi0*m_tanDip_phi0 + d_x_omega*m_tanDip_omega + d_x_tanDip*m_tanDip_tanDip  )   ; 
  xpCov(2,3) = 
    d_pz_omega* (  d_y_d0*m_omega_d0 + d_y_phi0*m_omega_phi0 + d_y_omega*m_omega_omega + d_y_tanDip*m_tanDip_omega  )   + 
    d_pz_tanDip* (  d_y_d0*m_tanDip_d0 + d_y_phi0*m_tanDip_phi0 + d_y_omega*m_tanDip_omega + d_y_tanDip*m_tanDip_tanDip  )   ; 
  xpCov(3,3) = 
    d_pz_omega* (  d_z_z0*m_z0_omega + d_z_tanDip*m_tanDip_omega  )   + 
    d_pz_tanDip* (  d_z_z0*m_tanDip_z0 + d_z_tanDip*m_tanDip_tanDip  )   ; 

  ppCov(1,1) = 
    d_px_phi0* (  d_px_phi0*m_phi0_phi0 + d_px_omega*m_omega_phi0 + d_px_tanDip*m_tanDip_phi0  )   + 
    d_px_omega* (  d_px_phi0*m_omega_phi0 + d_px_omega*m_omega_omega + d_px_tanDip*m_tanDip_omega  )   + 
    d_px_tanDip* (  d_px_phi0*m_tanDip_phi0 + d_px_omega*m_tanDip_omega + d_px_tanDip*m_tanDip_tanDip  )   ; 
  ppCov(2,1) = 
    d_py_phi0* (  d_px_phi0*m_phi0_phi0 + d_px_omega*m_omega_phi0 + d_px_tanDip*m_tanDip_phi0  )   + 
    d_py_omega* (  d_px_phi0*m_omega_phi0 + d_px_omega*m_omega_omega + d_px_tanDip*m_tanDip_omega  )   + 
    d_py_tanDip* (  d_px_phi0*m_tanDip_phi0 + d_px_omega*m_tanDip_omega + d_px_tanDip*m_tanDip_tanDip  )   ; 
  ppCov(2,2) = 
    d_py_phi0* (  d_py_phi0*m_phi0_phi0 + d_py_omega*m_omega_phi0 + d_py_tanDip*m_tanDip_phi0  )   + 
    d_py_omega* (  d_py_phi0*m_omega_phi0 + d_py_omega*m_omega_omega + d_py_tanDip*m_tanDip_omega  )   + 
    d_py_tanDip* (  d_py_phi0*m_tanDip_phi0 + d_py_omega*m_tanDip_omega + d_py_tanDip*m_tanDip_tanDip  )   ; 
  ppCov(3,1) = 
    d_pz_omega* (  d_px_phi0*m_omega_phi0 + d_px_omega*m_omega_omega + d_px_tanDip*m_tanDip_omega  )   + 
    d_pz_tanDip* (  d_px_phi0*m_tanDip_phi0 + d_px_omega*m_tanDip_omega + d_px_tanDip*m_tanDip_tanDip  )   ; 
  ppCov(3,2) = 
    d_pz_omega* (  d_py_phi0*m_omega_phi0 + d_py_omega*m_omega_omega + d_py_tanDip*m_tanDip_omega  )   + 
    d_pz_tanDip* (  d_py_phi0*m_tanDip_phi0 + d_py_omega*m_tanDip_omega + d_py_tanDip*m_tanDip_tanDip  )   ; 
  ppCov(3,3) = 
    d_pz_omega* (  d_pz_omega*m_omega_omega + d_pz_tanDip*m_tanDip_omega  )   + 
    d_pz_tanDip* (  d_pz_omega*m_tanDip_omega + d_pz_tanDip*m_tanDip_tanDip  )   ; 
  

}

//------------------------------------------------------------------------
void
TrkMomCalculator::calcCurvAllCovsOLD(const TrkSimpTraj& theTraj,
				  const BField& theField,
				  double fltlen,
				  HepSymMatrix& xxCov,
				  HepSymMatrix& ppCov,
				  HepMatrix&    xpCov) {
//------------------------------------------------------------------------
  
  DifPoint  PosDif;
  DifVector DirDif;
  DifVector delDirDif;
  DifVector MomDif(0., 0., 0.);

  theTraj.getDFInfo(fltlen, PosDif, DirDif, delDirDif);
  if (delDirDif.length() != 0) {   
    DifNumber sindip=DirDif.z;
    DifNumber arg = 1.0-sindip*sindip;
    if (arg.number() < 0.0) {arg.setNumber(0.0);}
    DifNumber cosdip = sqrt(arg);

    DifNumber momMag =
      BField::mmTeslaToMeVc*theField.bFieldNominal()*cosdip /
      delDirDif.length();
    momMag.absolute();

    MomDif = DirDif * momMag;

  }

  // computes position covariances
  
  xxCov.assign(PosDif.errorMatrix(PosDif.x.indepPar()->covariance()));

  // computes momentum covariances
  ppCov.assign(MomDif.errorMatrix(MomDif.x.indepPar()->covariance()));

  // computes correlations
  const HepSymMatrix& nnCov=MomDif.x.indepPar()->covariance();
  xpCov(1,1)=correlation(PosDif.x,MomDif.x,nnCov);
  xpCov(1,2)=correlation(PosDif.x,MomDif.y,nnCov);
  xpCov(1,3)=correlation(PosDif.x,MomDif.z,nnCov);
  xpCov(2,1)=correlation(PosDif.y,MomDif.x,nnCov);
  xpCov(2,2)=correlation(PosDif.y,MomDif.y,nnCov);
  xpCov(2,3)=correlation(PosDif.y,MomDif.z,nnCov);
  xpCov(3,1)=correlation(PosDif.z,MomDif.x,nnCov);
  xpCov(3,2)=correlation(PosDif.z,MomDif.y,nnCov);
  xpCov(3,3)=correlation(PosDif.z,MomDif.z,nnCov);
  
}

//------------------------------------------------------------------------
void
TrkMomCalculator::calcNeutAllCovs(const TrkSimpTraj& theTraj,
				  const BField& theField,
				  double fltlen,
				  HepSymMatrix& xxCov,
				  HepSymMatrix& ppCov,
				  HepMatrix&    xpCov) {
//------------------------------------------------------------------------
  HepVector p0(3),x0(3);
  calcNeutAllCovs(theTraj,theField,fltlen,x0,p0,xxCov,ppCov,xpCov);
  
}

//------------------------------------------------------------------------
void
TrkMomCalculator::calcNeutAllCovs(const TrkSimpTraj& theTraj,
				  const BField& theField,
				  double fltlen,
				  HepVector& x0,HepVector& p0,
				  HepSymMatrix& xxCov,
				  HepSymMatrix& ppCov,
				  HepMatrix&    xpCov) {
//------------------------------------------------------------------------

  DifPoint  PosDif;
  DifVector DirDif;
  DifVector delDirDif;
 
  theTraj.getDFInfo(fltlen, PosDif, DirDif, delDirDif);

  // set the momentum's direction, and then its magnitude
  DifVector MomDif = DirDif;

  // The + 1 is necessary because DifIndepPar starts counting at 1, whereas
  // the enum for NeutParams starts at 0.  Stupid, ain't it?
  MomDif *= theTraj.parameters()->difPar(NeutParams::_p + 1);

  // fill momenta and positions
  p0[0] = MomDif.x.number(); 
  p0[1] = MomDif.y.number(); 
  p0[2] = MomDif.z.number(); 
  x0[0] = PosDif.x.number(); 
  x0[1] = PosDif.y.number(); 
  x0[2] = PosDif.z.number(); 

  // computes momentum covariances
  xxCov.assign(PosDif.errorMatrix(PosDif.x.indepPar()->covariance()));

  // computes momentum covariances
  ppCov.assign(MomDif.errorMatrix(MomDif.x.indepPar()->covariance()));

  // computes correlations
  const HepSymMatrix& nnCov=MomDif.x.indepPar()->covariance();
  xpCov(1,1)=correlation(PosDif.x,MomDif.x,nnCov);
  xpCov(1,2)=correlation(PosDif.x,MomDif.y,nnCov);
  xpCov(1,3)=correlation(PosDif.x,MomDif.z,nnCov);
  xpCov(2,1)=correlation(PosDif.y,MomDif.x,nnCov);
  xpCov(2,2)=correlation(PosDif.y,MomDif.y,nnCov);
  xpCov(2,3)=correlation(PosDif.y,MomDif.z,nnCov);
  xpCov(3,1)=correlation(PosDif.z,MomDif.x,nnCov);
  xpCov(3,2)=correlation(PosDif.z,MomDif.y,nnCov);
  xpCov(3,3)=correlation(PosDif.z,MomDif.z,nnCov);

}

//------------------------------------------------------------------------
void
TrkMomCalculator::calcCurvAllWeights(const TrkSimpTraj& theTraj,
				     const BField& theField,
				     double fltlen,
				     HepVector& pos,
				     HepVector& mom,
				     HepSymMatrix& xxWeight,
				     HepSymMatrix& ppWeight,
				     HepMatrix&    xpWeight) {
//------------------------------------------------------------------------
  const HepVector&    v = theTraj.parameters()->parameter();
  const HepSymMatrix& w = theTraj.parameters()->weightMatrix();

  // fast inlined conversion from Helix to PX representation

  // track parameters
  double s = v[HelixParams::tanDipIndex];
  double omega = v[HelixParams::omegaIndex];
  double d0 = v[HelixParams::d0Index];
  double z0 = v[HelixParams::z0Index];
  double phi0 = v[HelixParams::phi0Index];
  double l = fltlen / sqrt(1.0+s*s) ;

  double phi = phi0 + omega*l;
  double sinphi0 = sin(phi0);
  double cosphi0 = cos(phi0);
  double sinphi = sin(phi);
  double cosphi = cos(phi);
  double C = BField::mmTeslaToMeVc * theField.bFieldNominal();
  double q(1.0); 
  -omega>0 ? q=1.0: q=-1.0;
  double qC = q*C;
  double pt = fabs( -qC / omega );
  double r = 1.0/omega;

  // calculate position and momentum
  pos(1) = r*sinphi - (r + d0)*sinphi0;
  pos(2) = -r*cosphi + (r + d0)*cosphi0;
  pos(3) = z0 + l*s;
  mom(1) = pt*cosphi;
  mom(2) = pt*sinphi;
  mom(3) = pt*s;
 
  // inverse of the jacobian matrix - see helix.mws Maple worksheet

  // protect against divide by 0.
  if ( (1+d0*omega)==0.0 ){
    calcCurvAllWeightsOLD(theTraj,theField,fltlen,
			  pos,mom,xxWeight,ppWeight,xpWeight);
    return;
  }

  double dinv_x_d0 = -sinphi0;
  double dinv_x_phi0 = -omega * cosphi0 / (1 + d0 * omega);
  double dinv_x_z0 = -s * cosphi0 / (1 + d0 * omega);
  
  double dinv_y_d0 = cosphi0;
  double dinv_y_phi0 = -omega * sinphi0 / (1 + d0 * omega);
  double dinv_y_z0 = -s * sinphi0 / (1 + d0 * omega);
  
  double dinv_z_z0 = 1;
  
  double dinv_px_d0 = (cosphi - cosphi0) / qC;
  double dinv_px_phi0 = omega * sinphi0 / (1 + d0 * omega) / qC;
  double dinv_px_omega = omega * omega * cosphi / qC;
  double dinv_px_z0 = -s * ( sinphi*(1 + d0 * omega) - sinphi0 ) / (qC * (1 + d0 * omega));
  double dinv_px_tanDip = omega * cosphi * s / qC;
  
  double dinv_py_d0 = (sinphi - sinphi0) / qC;
  double dinv_py_phi0 = -omega * cosphi0 / (1 + d0 * omega) / qC;
  double dinv_py_omega = omega * omega * sinphi / qC;
  double dinv_py_z0 = s * ( cosphi*(1 + d0 * omega) - cosphi0) / (qC * (1 + d0 * omega));
  double dinv_py_tanDip = omega * sinphi * s / qC;
  
  double dinv_pz_z0 = l * omega / qC;
  double dinv_pz_tanDip = -omega / qC;

  // local variables for the weight matrix
  double w_d0_d0 =  w[HelixParams::d0Index][HelixParams::d0Index];
  double w_phi0_d0 =  w[HelixParams::phi0Index][HelixParams::d0Index];
  double w_phi0_phi0 =  w[HelixParams::phi0Index][HelixParams::phi0Index];
  double w_omega_d0 =  w[HelixParams::omegaIndex][HelixParams::d0Index];
  double w_omega_phi0 =  w[HelixParams::omegaIndex][HelixParams::phi0Index];
  double w_omega_omega =  w[HelixParams::omegaIndex][HelixParams::omegaIndex];
  double w_z0_d0 =  w[HelixParams::z0Index][HelixParams::d0Index];
  double w_z0_phi0 =  w[HelixParams::z0Index][HelixParams::phi0Index];
  double w_z0_omega =  w[HelixParams::z0Index][HelixParams::omegaIndex];
  double w_z0_z0 =  w[HelixParams::z0Index][HelixParams::z0Index];
  double w_tanDip_d0 =  w[HelixParams::tanDipIndex][HelixParams::d0Index];
  double w_tanDip_phi0 =  w[HelixParams::tanDipIndex][HelixParams::phi0Index];
  double w_tanDip_omega =  w[HelixParams::tanDipIndex][HelixParams::omegaIndex];
  double w_tanDip_z0 =  w[HelixParams::tanDipIndex][HelixParams::z0Index];
  double w_tanDip_tanDip =  w[HelixParams::tanDipIndex][HelixParams::tanDipIndex];

  // calculate the Weight matrix from dinv (the indices are xpWeight(ip,ix) ) (used writewgt.py script)
  xxWeight(1,1) = 
    dinv_x_d0* ( dinv_x_d0*w_d0_d0 +dinv_x_phi0*w_phi0_d0 +dinv_x_z0*w_z0_d0  )   + 
    dinv_x_phi0* ( dinv_x_d0*w_phi0_d0 +dinv_x_phi0*w_phi0_phi0 +dinv_x_z0*w_z0_phi0  )   + 
    dinv_x_z0* ( dinv_x_d0*w_z0_d0 +dinv_x_phi0*w_z0_phi0 +dinv_x_z0*w_z0_z0  )   ; 
  xxWeight(2,1) = 
    dinv_y_d0* ( dinv_x_d0*w_d0_d0 +dinv_x_phi0*w_phi0_d0 +dinv_x_z0*w_z0_d0  )   + 
    dinv_y_phi0* ( dinv_x_d0*w_phi0_d0 +dinv_x_phi0*w_phi0_phi0 +dinv_x_z0*w_z0_phi0  )   + 
    dinv_y_z0* ( dinv_x_d0*w_z0_d0 +dinv_x_phi0*w_z0_phi0 +dinv_x_z0*w_z0_z0  )   ; 
  xxWeight(2,2) = 
    dinv_y_d0* ( dinv_y_d0*w_d0_d0 +dinv_y_phi0*w_phi0_d0 +dinv_y_z0*w_z0_d0  )   + 
    dinv_y_phi0* ( dinv_y_d0*w_phi0_d0 +dinv_y_phi0*w_phi0_phi0 +dinv_y_z0*w_z0_phi0  )   + 
    dinv_y_z0* ( dinv_y_d0*w_z0_d0 +dinv_y_phi0*w_z0_phi0 +dinv_y_z0*w_z0_z0  )   ; 
  xxWeight(3,1) = 
    dinv_z_z0* ( dinv_x_d0*w_z0_d0 +dinv_x_phi0*w_z0_phi0 +dinv_x_z0*w_z0_z0  )   ; 
  xxWeight(3,2) = 
    dinv_z_z0* ( dinv_y_d0*w_z0_d0 +dinv_y_phi0*w_z0_phi0 +dinv_y_z0*w_z0_z0  )   ; 
  xxWeight(3,3) = 
    dinv_z_z0* ( dinv_z_z0*w_z0_z0  )   ; 

  xpWeight(1,1) = 
    dinv_px_d0* ( dinv_x_d0*w_d0_d0 +dinv_x_phi0*w_phi0_d0 +dinv_x_z0*w_z0_d0  )   + 
    dinv_px_phi0* ( dinv_x_d0*w_phi0_d0 +dinv_x_phi0*w_phi0_phi0 +dinv_x_z0*w_z0_phi0  )   + 
    dinv_px_omega* ( dinv_x_d0*w_omega_d0 +dinv_x_phi0*w_omega_phi0 +dinv_x_z0*w_z0_omega  )   + 
    dinv_px_z0* ( dinv_x_d0*w_z0_d0 +dinv_x_phi0*w_z0_phi0 +dinv_x_z0*w_z0_z0  )   + 
    dinv_px_tanDip* ( dinv_x_d0*w_tanDip_d0 +dinv_x_phi0*w_tanDip_phi0 +dinv_x_z0*w_tanDip_z0  )   ; 
  xpWeight(2,1) = 
    dinv_px_d0* ( dinv_y_d0*w_d0_d0 +dinv_y_phi0*w_phi0_d0 +dinv_y_z0*w_z0_d0  )   + 
    dinv_px_phi0* ( dinv_y_d0*w_phi0_d0 +dinv_y_phi0*w_phi0_phi0 +dinv_y_z0*w_z0_phi0  )   + 
    dinv_px_omega* ( dinv_y_d0*w_omega_d0 +dinv_y_phi0*w_omega_phi0 +dinv_y_z0*w_z0_omega  )   + 
    dinv_px_z0* ( dinv_y_d0*w_z0_d0 +dinv_y_phi0*w_z0_phi0 +dinv_y_z0*w_z0_z0  )   + 
    dinv_px_tanDip* ( dinv_y_d0*w_tanDip_d0 +dinv_y_phi0*w_tanDip_phi0 +dinv_y_z0*w_tanDip_z0  )   ; 
  xpWeight(3,1) = 
    dinv_px_d0* ( dinv_z_z0*w_z0_d0  )   + 
    dinv_px_phi0* ( dinv_z_z0*w_z0_phi0  )   + 
    dinv_px_omega* ( dinv_z_z0*w_z0_omega  )   + 
    dinv_px_z0* ( dinv_z_z0*w_z0_z0  )   + 
    dinv_px_tanDip* ( dinv_z_z0*w_tanDip_z0  )   ; 
  xpWeight(1,2) = 
    dinv_py_d0* ( dinv_x_d0*w_d0_d0 +dinv_x_phi0*w_phi0_d0 +dinv_x_z0*w_z0_d0  )   + 
    dinv_py_phi0* ( dinv_x_d0*w_phi0_d0 +dinv_x_phi0*w_phi0_phi0 +dinv_x_z0*w_z0_phi0  )   + 
    dinv_py_omega* ( dinv_x_d0*w_omega_d0 +dinv_x_phi0*w_omega_phi0 +dinv_x_z0*w_z0_omega  )   + 
    dinv_py_z0* ( dinv_x_d0*w_z0_d0 +dinv_x_phi0*w_z0_phi0 +dinv_x_z0*w_z0_z0  )   + 
    dinv_py_tanDip* ( dinv_x_d0*w_tanDip_d0 +dinv_x_phi0*w_tanDip_phi0 +dinv_x_z0*w_tanDip_z0  )   ; 
  xpWeight(2,2) = 
    dinv_py_d0* ( dinv_y_d0*w_d0_d0 +dinv_y_phi0*w_phi0_d0 +dinv_y_z0*w_z0_d0  )   + 
    dinv_py_phi0* ( dinv_y_d0*w_phi0_d0 +dinv_y_phi0*w_phi0_phi0 +dinv_y_z0*w_z0_phi0  )   + 
    dinv_py_omega* ( dinv_y_d0*w_omega_d0 +dinv_y_phi0*w_omega_phi0 +dinv_y_z0*w_z0_omega  )   + 
    dinv_py_z0* ( dinv_y_d0*w_z0_d0 +dinv_y_phi0*w_z0_phi0 +dinv_y_z0*w_z0_z0  )   + 
    dinv_py_tanDip* ( dinv_y_d0*w_tanDip_d0 +dinv_y_phi0*w_tanDip_phi0 +dinv_y_z0*w_tanDip_z0  )   ; 
  xpWeight(3,2) = 
    dinv_py_d0* ( dinv_z_z0*w_z0_d0  )   + 
    dinv_py_phi0* ( dinv_z_z0*w_z0_phi0  )   + 
    dinv_py_omega* ( dinv_z_z0*w_z0_omega  )   + 
    dinv_py_z0* ( dinv_z_z0*w_z0_z0  )   + 
    dinv_py_tanDip* ( dinv_z_z0*w_tanDip_z0  )   ; 
  xpWeight(1,3) = 
    dinv_pz_z0* ( dinv_x_d0*w_z0_d0 +dinv_x_phi0*w_z0_phi0 +dinv_x_z0*w_z0_z0  )   + 
    dinv_pz_tanDip* ( dinv_x_d0*w_tanDip_d0 +dinv_x_phi0*w_tanDip_phi0 +dinv_x_z0*w_tanDip_z0  )   ; 
  xpWeight(2,3) = 
    dinv_pz_z0* ( dinv_y_d0*w_z0_d0 +dinv_y_phi0*w_z0_phi0 +dinv_y_z0*w_z0_z0  )   + 
    dinv_pz_tanDip* ( dinv_y_d0*w_tanDip_d0 +dinv_y_phi0*w_tanDip_phi0 +dinv_y_z0*w_tanDip_z0  )   ; 
  xpWeight(3,3) = 
    dinv_pz_z0* ( dinv_z_z0*w_z0_z0  )   + 
    dinv_pz_tanDip* ( dinv_z_z0*w_tanDip_z0  )   ; 


  ppWeight(1,1) = 
    dinv_px_d0* ( dinv_px_d0*w_d0_d0 +dinv_px_phi0*w_phi0_d0 +dinv_px_omega*w_omega_d0 +dinv_px_z0*w_z0_d0 +dinv_px_tanDip*w_tanDip_d0  )   + 
    dinv_px_phi0* ( dinv_px_d0*w_phi0_d0 +dinv_px_phi0*w_phi0_phi0 +dinv_px_omega*w_omega_phi0 +dinv_px_z0*w_z0_phi0 +dinv_px_tanDip*w_tanDip_phi0  )   + 
    dinv_px_omega* ( dinv_px_d0*w_omega_d0 +dinv_px_phi0*w_omega_phi0 +dinv_px_omega*w_omega_omega +dinv_px_z0*w_z0_omega +dinv_px_tanDip*w_tanDip_omega  )   + 
    dinv_px_z0* ( dinv_px_d0*w_z0_d0 +dinv_px_phi0*w_z0_phi0 +dinv_px_omega*w_z0_omega +dinv_px_z0*w_z0_z0 +dinv_px_tanDip*w_tanDip_z0  )   + 
    dinv_px_tanDip* ( dinv_px_d0*w_tanDip_d0 +dinv_px_phi0*w_tanDip_phi0 +dinv_px_omega*w_tanDip_omega +dinv_px_z0*w_tanDip_z0 +dinv_px_tanDip*w_tanDip_tanDip  )   ; 
  ppWeight(2,1) = 
    dinv_py_d0* ( dinv_px_d0*w_d0_d0 +dinv_px_phi0*w_phi0_d0 +dinv_px_omega*w_omega_d0 +dinv_px_z0*w_z0_d0 +dinv_px_tanDip*w_tanDip_d0  )   + 
    dinv_py_phi0* ( dinv_px_d0*w_phi0_d0 +dinv_px_phi0*w_phi0_phi0 +dinv_px_omega*w_omega_phi0 +dinv_px_z0*w_z0_phi0 +dinv_px_tanDip*w_tanDip_phi0  )   + 
    dinv_py_omega* ( dinv_px_d0*w_omega_d0 +dinv_px_phi0*w_omega_phi0 +dinv_px_omega*w_omega_omega +dinv_px_z0*w_z0_omega +dinv_px_tanDip*w_tanDip_omega  )   + 
    dinv_py_z0* ( dinv_px_d0*w_z0_d0 +dinv_px_phi0*w_z0_phi0 +dinv_px_omega*w_z0_omega +dinv_px_z0*w_z0_z0 +dinv_px_tanDip*w_tanDip_z0  )   + 
    dinv_py_tanDip* ( dinv_px_d0*w_tanDip_d0 +dinv_px_phi0*w_tanDip_phi0 +dinv_px_omega*w_tanDip_omega +dinv_px_z0*w_tanDip_z0 +dinv_px_tanDip*w_tanDip_tanDip  )   ; 
  ppWeight(2,2) = 
    dinv_py_d0* ( dinv_py_d0*w_d0_d0 +dinv_py_phi0*w_phi0_d0 +dinv_py_omega*w_omega_d0 +dinv_py_z0*w_z0_d0 +dinv_py_tanDip*w_tanDip_d0  )   + 
    dinv_py_phi0* ( dinv_py_d0*w_phi0_d0 +dinv_py_phi0*w_phi0_phi0 +dinv_py_omega*w_omega_phi0 +dinv_py_z0*w_z0_phi0 +dinv_py_tanDip*w_tanDip_phi0  )   + 
    dinv_py_omega* ( dinv_py_d0*w_omega_d0 +dinv_py_phi0*w_omega_phi0 +dinv_py_omega*w_omega_omega +dinv_py_z0*w_z0_omega +dinv_py_tanDip*w_tanDip_omega  )   + 
    dinv_py_z0* ( dinv_py_d0*w_z0_d0 +dinv_py_phi0*w_z0_phi0 +dinv_py_omega*w_z0_omega +dinv_py_z0*w_z0_z0 +dinv_py_tanDip*w_tanDip_z0  )   + 
    dinv_py_tanDip* ( dinv_py_d0*w_tanDip_d0 +dinv_py_phi0*w_tanDip_phi0 +dinv_py_omega*w_tanDip_omega +dinv_py_z0*w_tanDip_z0 +dinv_py_tanDip*w_tanDip_tanDip  )   ; 
  ppWeight(3,1) = 
    dinv_pz_z0* ( dinv_px_d0*w_z0_d0 +dinv_px_phi0*w_z0_phi0 +dinv_px_omega*w_z0_omega +dinv_px_z0*w_z0_z0 +dinv_px_tanDip*w_tanDip_z0  )   + 
    dinv_pz_tanDip* ( dinv_px_d0*w_tanDip_d0 +dinv_px_phi0*w_tanDip_phi0 +dinv_px_omega*w_tanDip_omega +dinv_px_z0*w_tanDip_z0 +dinv_px_tanDip*w_tanDip_tanDip  )   ; 
  ppWeight(3,2) = 
    dinv_pz_z0* ( dinv_py_d0*w_z0_d0 +dinv_py_phi0*w_z0_phi0 +dinv_py_omega*w_z0_omega +dinv_py_z0*w_z0_z0 +dinv_py_tanDip*w_tanDip_z0  )   + 
    dinv_pz_tanDip* ( dinv_py_d0*w_tanDip_d0 +dinv_py_phi0*w_tanDip_phi0 +dinv_py_omega*w_tanDip_omega +dinv_py_z0*w_tanDip_z0 +dinv_py_tanDip*w_tanDip_tanDip  )   ; 
  ppWeight(3,3) = 
    dinv_pz_z0* ( dinv_pz_z0*w_z0_z0 +dinv_pz_tanDip*w_tanDip_z0  )   + 
    dinv_pz_tanDip* ( dinv_pz_z0*w_tanDip_z0 +dinv_pz_tanDip*w_tanDip_tanDip  )   ; 

}

//------------------------------------------------------------------------
void
TrkMomCalculator::calcCurvAllWeightsOLD(const TrkSimpTraj& theTraj,
				     const BField& theField,
				     double fltlen,
				     HepVector& pos,
				     HepVector& mom,
				     HepSymMatrix& xxWeight,
				     HepSymMatrix& ppWeight,
				     HepMatrix&    xpWeight) {
//------------------------------------------------------------------------

  DifPoint  PosDif;
  DifVector DirDif;
  DifVector delDirDif;
  DifNumber momMag;
  DifVector MomDif;

  theTraj.getDFInfo(fltlen, PosDif, DirDif, delDirDif);
  if (delDirDif.length() != 0) {   
    DifNumber sindip=DirDif.z;
    DifNumber arg = 1.0-sindip*sindip;
    if (arg.number() < 0.0) {arg.setNumber(0.0);}
    DifNumber cosdip = sqrt(arg);
    
    momMag =
      BField::mmTeslaToMeVc*theField.bFieldNominal()*cosdip /
      delDirDif.length();
    momMag.absolute();

    MomDif = DirDif * momMag;

  }

  /*
   * start computing the inverse of the Jacobian needed:
   *
   *      D(x,p)
   *      ------
   *       D(n)
   */
  HepMatrix Jx_n(PosDif.jacobian());
  HepMatrix Jp_n(MomDif.jacobian());

  int          i,j;
  HepMatrix    Jxp_ns(6,6);

  for(i=0;i<3;i++)
    for(j=0;j<5;j++)
      {
	Jxp_ns[i  ][j]=Jx_n[i][j];
	Jxp_ns[i+3][j]=Jp_n[i][j];
      }

  /*
   * now we need: 
   *
   *     D(x,p)
   *     ------
   *      D(s)
   *
   */

  Jxp_ns[0][5]=DirDif.x.number();
  Jxp_ns[1][5]=DirDif.y.number();
  Jxp_ns[2][5]=DirDif.z.number();

  Jxp_ns[3][5]=momMag.number()*delDirDif.x.number();
  Jxp_ns[4][5]=momMag.number()*delDirDif.y.number();
  Jxp_ns[5][5]=momMag.number()*delDirDif.z.number();

  /*
   * with an inversion we obtain
   *     D(n)
   *   --------
   *    D(x,p)
   *
   */
  int invStatus;
  
  Jxp_ns.invert(invStatus);
  
  HepMatrix Jn_x(5,3);
  HepMatrix Jn_p(5,3);
  
  for(i=0;i<5;i++)
    for(j=0;j<3;j++)
      {
	Jn_x[i][j]=Jxp_ns[i][j  ];
	Jn_p[i][j]=Jxp_ns[i][j+3];
      }
  // this is the weight matrix for the helix parameters
  HepSymMatrix Wnn(PosDif.x.indepPar()->covariance());

  Wnn.invert(invStatus);
  /*
   * now we have the weight matrices
   *
   */
  xxWeight   = Wnn.similarityT(Jn_x);
  ppWeight   = Wnn.similarityT(Jn_p);
  xpWeight   = Jn_x.T()*Wnn*Jn_p;

  pos[0]=PosDif.x.number();
  pos[1]=PosDif.y.number();
  pos[2]=PosDif.z.number();

  mom[0]=MomDif.x.number();
  mom[1]=MomDif.y.number();
  mom[2]=MomDif.z.number();

}

//------------------------------------------------------------------------
void
TrkMomCalculator::calcNeutAllWeights(const TrkSimpTraj& theTraj,
				     const BField& theField,
				     double fltlen,
				     HepVector& pos,
				     HepVector& mom,
				     HepSymMatrix& xxWeight,
				     HepSymMatrix& ppWeight,
				     HepMatrix&    xpWeight) {
//------------------------------------------------------------------------
  DifPoint  PosDif;
  DifVector DirDif;
  DifVector delDirDif;
  DifNumber momMag;

  theTraj.getDFInfo(fltlen, PosDif, DirDif, delDirDif);

  // set the momentum's direction, and then its magnitude
  DifVector MomDif = DirDif;

  MomDif *= theTraj.parameters()->difPar(NeutParams::_p + 1);

  HepMatrix Jx_n(PosDif.jacobian());
  HepMatrix Jp_n(MomDif.jacobian());

  int          i,j;
  HepMatrix    Jxp_ns(6,6);

  for(i=0;i<3;i++)
    for(j=0;j<6;j++)
      {
	Jxp_ns[i  ][j]=Jx_n[i][j];
	Jxp_ns[i+3][j]=Jp_n[i][j];
      }
  int invStatus;
  
  Jxp_ns.invert(invStatus);
  
  HepMatrix Jn_x(5,3);
  HepMatrix Jn_p(5,3);
  
  for(i=0;i<5;i++)
    for(j=0;j<3;j++)
      {
	Jn_x[i][j]=Jxp_ns[i][j  ];
	Jn_p[i][j]=Jxp_ns[i][j+3];
      }
  // this is the weight matrix for the helix parameters
  HepSymMatrix Wnn(PosDif.x.indepPar()->covariance().sub(1,5));

  Wnn.invert(invStatus);
  xxWeight   = Wnn.similarityT(Jn_x);
  ppWeight   = Wnn.similarityT(Jn_p);
  xpWeight   = Jn_x.T()*Wnn*Jn_p;

  pos[0]=PosDif.x.number();
  pos[1]=PosDif.y.number();
  pos[2]=PosDif.z.number();

  mom[0]=MomDif.x.number();
  mom[1]=MomDif.y.number();
  mom[2]=MomDif.z.number();
}

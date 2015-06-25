// $Id: HelixTraj.hh,v 1.7 2007/09/24 21:56:26 gapon Exp $
//
// Trajectory subclass that implements a helix simple along the z axis
// Author : Gautier Hamel de Monchenault, from TrkParms and other similar classes
//
#ifndef HELIXTRAJ_HH
#define HELIXTRAJ_HH
#include "BTrk/TrkBase/TrkSimpTraj.hh"

class TrkVisitor;
#include "CLHEP/Matrix/Vector.h"
class HelixParams;

class HelixTraj : public TrkSimpTraj {
public:
//  Define the parameters 
  enum ParIndex {d0Index=0, phi0Index, omegaIndex, z0Index, tanDipIndex, NHLXPRM};
  HelixTraj(const CLHEP::HepVector&,const CLHEP::HepSymMatrix&, double lowlim=-99999.,
	    double hilim=99999., const HepPoint& refpoint = _theOrigin);
  HelixTraj(const HelixParams&, double lowlim=-99999.,
	    double hilim=99999., const HepPoint& refpoint = _theOrigin);
  HelixTraj(const TrkParams&, double lowlim=-99999.,
	    double hilim=99999., const HepPoint& refpoint = _theOrigin);
  HelixTraj( const HelixTraj&  );  
  HelixTraj* clone() const;

  virtual ~HelixTraj();

// operators
  HelixTraj& operator=(const HelixTraj&);
  virtual HepPoint   position(double fltLen)  const;
  virtual CLHEP::Hep3Vector direction(double fltLen) const;
  virtual CLHEP::Hep3Vector delDirect(double)        const;
  virtual void       getInfo(double fltLen, HepPoint& pos, 
			     CLHEP::Hep3Vector& dir) const;
  virtual void       getInfo(double fltLen, HepPoint& , CLHEP::Hep3Vector& dir, 
			     CLHEP::Hep3Vector& delDir) const;
  virtual void       getDFInfo(double fltLen, DifPoint& , DifVector& dir, 
			       DifVector& delDir) const;
  virtual void       getDFInfo2(double fltLen, DifPoint& pos, DifVector& 
			       dir) const;

  // How far can you go using given approximation before error > tolerance, 
  //   in direction pathDir?
  virtual double distTo1stError(double s, double tol, int pathDir) const;
  virtual double distTo2ndError(double s, double tol, int pathDir) const;


//  Real versions of the base class derrivative functions
  double curvature( double fltLen) const;
  CLHEP::HepMatrix derivDeflect(double fltlen,deflectDirection) const;
  CLHEP::HepMatrix derivDisplace(double fltlen,deflectDirection idir) const;
  CLHEP::HepMatrix derivPFract(double fltlen) const;
//  PointTraj functions
  TranslateParams paramFunction() const { return HelixTraj::paramFunc; }
  // Invert the parameters.  Returns true in flags if the inversion
  // requires a change of sign in the covariance matrix.
  void invertParams(TrkParams* params, std::vector<bool>& flags) const;
  // Helix-specific accessors
  double d0() const {return parameters()->parameter()[d0Index];}
  double phi0() const;
  double omega() const {return parameters()->parameter()[omegaIndex]; }
  double z0() const {return parameters()->parameter()[z0Index]; }
  double tanDip() const {  return parameters()->parameter()[tanDipIndex]; }

  //--------------------------------------------------
  // Visitor access
  //--------------------------------------------------

  virtual void visitAccept(TrkVisitor* vis) const;

  virtual void               print(std::ostream& os) const;
  virtual                    void printAll(std::ostream& os) const;

  //  double x( const double& ) const;
  //  double y( const double& ) const;
  double z( const double& ) const;
  double zFlight(double zpos) const;
  double dip() const {return atan(tanDip());}
  double cosDip() const {return 1./sqrt(1.+(tanDip()*tanDip())); }
  double sinDip() const {return tanDip()*cosDip(); }
  double translen(const double& f) const {return cosDip()*f;}
  double arc( const double& f) const {return translen(f)*omega();}
  double angle(const double& f) const;
private:
// the real point translation function
  static void paramFunc(const HepPoint& oldpoint,const HepPoint& newpoint,
			const CLHEP::HepVector& oldpar,const CLHEP::HepSymMatrix& oldcov,
			CLHEP::HepVector& newpar,CLHEP::HepSymMatrix& newcov,
			double fltlen);
};
#endif


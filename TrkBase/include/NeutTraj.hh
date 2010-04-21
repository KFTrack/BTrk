//
// Trajectory subclass that implements a neutral simple
// Author : Justin Albert and Valery Miftahov, from HelixTraj
//
#ifndef NEUTTRAJ_HH
#define NEUTTRAJ_HH
#include "TrkBase/TrkSimpTraj.hh"
#include "TrkBase/NeutParams.hh"
#include <algorithm>

class TrkVisitor;
class NeutParams;

class NeutTraj : public TrkSimpTraj
{
public:
  NeutTraj( const NeutParams&, double lowlim=-99999.,double hilim=99999.,
            const HepPoint& refpoint = _theOrigin);
  NeutTraj( const NeutTraj&  );   // copy ctor
  NeutTraj* clone() const;

// destructor
  virtual ~NeutTraj();

// operators
  NeutTraj& operator=(const NeutTraj&);
// needed implementations for intersection with a Surface
  virtual HepPoint   position( double fltLen) const;
  virtual Hep3Vector direction( double fltLen) const;
  virtual Hep3Vector delDirect( double ) const;
  void  getInfo( double fltLen, HepPoint& pos, Hep3Vector& dir ) const;
  void  getInfo( double fltLen, HepPoint& , Hep3Vector& dir,
			Hep3Vector& delDir ) const;
  void  getDFInfo( double fltLen, DifPoint& , DifVector& dir,
			DifVector& delDir ) const;

  // How far can you go using given approximation before error > tolerance,
  //   in direction pathDir?
  virtual double distTo1stError(double s, double tol, int pathDir) const;
  virtual double distTo2ndError(double s, double tol, int pathDir) const;

  // Expose the parameters for user manipulation
  NeutParams& params()                                    {return *_np();}
  const NeutParams& params() const                        {return *_np();}
//
//  Real versions of the base class derivative functions
//
  double curvature( double fltLen) const                  {return 0.;}
  HepMatrix derivDeflect(double fltlen,deflectDirection) const;
  HepMatrix derivDisplace(double fltlen,deflectDirection idir) const;
  HepMatrix derivPFract(double fltlen) const;
//  PointTraj functions
  TranslateParams paramFunction() const { return NeutTraj::paramFunc; }
  // Invert the parameter.  Returns true in flags if the inversion
  //requires a change of sign in the covariance matrix.
  void invertParams(TrkParams* newparams, std::vector<bool>& flags) const;

  TrkSimpTraj& invert();

  //--------------------------------------------------
  // Visitor access
  //--------------------------------------------------

  virtual void visitAccept(TrkVisitor* vis) const;


private:


//
//  Private functions (data members are part of the base class)
//
  double x( const double& ) const;
  double y( const double& ) const;
  double z( const double& ) const;
  inline NeutParams* _np() const {return (NeutParams*) &_dtparams; }
  inline double d0() const {return _np()->d0(); }
         double phi0() const;
  inline double p() const {return _np()->p(); }
  inline double z0() const {return _np()->z0(); }
  inline double s0() const {return _np()->s0(); }
  inline double tanDip() const {  return _np()->tanDip(); }
  inline double dip() const {return atan(tanDip());}
  inline double cosDip() const {return 1./sqrt(1.+std::sqr(tanDip())); }
  inline double sinDip() const {return tanDip()*cosDip(); }
  inline double translen(const double& f) const {return cosDip()*f;}
// the real point translation function
  static void paramFunc(const HepPoint& oldpoint,const HepPoint& newpoint,
                        const HepVector& oldpar,const HepSymMatrix& oldcov,
                        HepVector& newpar,HepSymMatrix& newcov,
			double fltlen);
};
#endif

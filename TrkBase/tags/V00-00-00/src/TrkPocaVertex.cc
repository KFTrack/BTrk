//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkPocaVertex.cc,v 1.9 2007/05/18 03:39:52 kelsey Exp $
//
// Description:
//
// Utility to vertex two tracks, or a track and a point, using TrkPoca.
//
// I'll put the algebra for deriving chisquares and errors here. A
// track is a line, yet we'll represent it as a point and a
// corresponding error. Given the poca of the two tracks, we would
// like to calculate the 'minimum chisquare' position based on a
// weighted average of the two track positions. However, since the
// track is a line, the covariance matrix of the point cannot be
// inverted. What we do is to make a projection in the plane
// perpendicular to the track. The chisquare contribution can be
// written as:
//
//    chisquare = ( P (xvtx-xtrk) )^T (P^T Vtrk P)^-1 P (xvtx-xtrk)
//
// where P is the appropriate projection matrix. Note that this
// is equivalent to
//
//    chisquare = (xvtx-xtrk)^T   P^T(P^T Vtrk P)^-1 P    (xvtx-xtrk)
// 
// where the middle term is now the weight matrix 'W'. We can obtain
// the best position from the weighted average of the two tracks
//
//    xvtx = (W1 + W2)^-1 * ( W1 xtrk1 + W2 xtrk2 )             [1]
//    Cvtx = (W1 + W2)^-1                                       [2]
//  
// As you can see below we will use this calculation only for
// extracting the cov matrix and not the actual best position. The
// reason is that this calculation is unstable if the two tracks are
// nearly parallel, as is the case for conversions. Therefore, we will
// actually rely on TrkPoca itself, to find the place where the tracks
// are closest and obtain the best vertex position from a weighted
// average of the two points _along_ the residual only: The expression
// is in fact the same as [1] above, except that the weights are
// just numbers, not projection matrices.
// 
// Because the chisquare of the trk-trk vertex has actually only 1
// degree of freedom, we can also write it as
//
//    chisquare =      doca^2/(var1 + var2)
//
// where doca is just doca = |residual| = |xtrk-xtrk2|. The
// 1-dimensional weights are given by the projection of the covariance
// matrix on the residual vector. This procedure gives us both the
// chisquare and the best vertex, but not the covariance matrix. For
// the latter, we need the full calculation.
//
// For the trk-to-point vertex we cannot apply this trick since the
// chisquare has two degrees of freedom. Therefore, we use the full 3D
// calculation. In principle, we can project out the two dimensions
// only, to speed things a bit up. That is for later.
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author(s): Wouter Hulsbergen
//
//------------------------------------------------------------------------

#include "BaBar/BaBar.hh"
#include "TrkBase/TrkPocaVertex.hh"
#include "TrkBase/TrkFit.hh"
#include "TrkBase/TrkRecoTrk.hh"
#include "TrkBase/TrkDifTraj.hh"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "PDT/Pdt.hh"
#include <math.h>

template<class T>
HepVector 
TrkPocaVertex::convert(const T& vec) 
{
  HepVector rc(3) ;
  rc(1) = vec.x() ;
  rc(2) = vec.y() ;
  rc(3) = vec.z() ;
  return rc ;
}

HepVector
TrkPocaVertex::crossproduct(const HepVector& v1,const HepVector& v2) 
{
  HepVector rc(3) ;
  rc(1) = v1(2)*v2(3) - v1(3)*v2(2) ;
  rc(2) = v1(3)*v2(1) - v1(1)*v2(3) ;
  rc(3) = v1(1)*v2(2) - v1(2)*v2(1) ;
  return rc ;
} 

TrkPocaVertex::TrkPocaVertex(const TrkRecoTrk* trk1, double flt1, PdtPid::PidType pid1,
			     const TrkRecoTrk* trk2, double flt2, PdtPid::PidType pid2,
			     const PdtEntry* vtype,
			     double precision) : TrkVertex(vtype)
{
  _usedtracks.push_back(trk1);
  _usedtracks.push_back(trk2);
  const TrkFit* fit1 = trk1->fitResult(pid1);
  const TrkFit* fit2 = trk2->fitResult(pid2);
  if(fit1 != 0 && fit2 != 0)
    _poca = TrkPoca(fit1->traj(), flt1, fit2->traj(), flt2, precision );
  _info[trk1] = TrkVtxInfo(_poca.doca(),_poca.flt1(),pid1);
  _info[trk2] = TrkVtxInfo(_poca.doca(),_poca.flt2(),pid2);
  _status = _poca.status();
  if(_status.success()) {
    const TrkAbsFit* trkfit[2]  ;
    double  flt[2], mass[2] ;
    trkfit[0] = fit1 ; flt[0] = _poca.flt1() ; mass[0] = Pdt::lookup(pid1)->mass();
    trkfit[1] = fit2 ; flt[1] = _poca.flt2() ; mass[1] = Pdt::lookup(pid2)->mass();
 
    HepSymMatrix cov[2] ;
    HepVector    pos[2] ;
    HepVector    dir[2] ;
    double       vardoca[2] ;

    for(int i=0; i<2; ++i) {
      BbrPointErr position = trkfit[i]->positionErr(flt[i]) ;
      //std::cout << "pos at poca: " << position.x() << " " << position.y() << " " << position.z() << std::endl; 
      cov[i] = position.covMatrix() ;
      pos[i] = convert<HepPoint>(position) ;
      BbrVectorErr momentum = trkfit[i]->momentumErr(flt[i] ) ;
      dir[i] = convert<Hep3Vector>(momentum.unit()) ;
      _p4 += BbrLorentzVectorErr( momentum, mass[i] ) ; 
    }

    // calculate the chisquare (a bit of hocus pocus)
    HepVector residual = pos[1] - pos[0] ;
    double doca = residual.norm() ;
    // in principle, the 'jacobian' is along the residual. however, this
    // might be a bit unstable for small residuals. so, we could take
    // the cross product of the momentum vectors ... which is unsuitable
    // for parallel tracks, like conversions. Let's take the residual after all.
    HepVector jacobian   = residual ;
    if( fabs(doca) < 1e-12) jacobian = crossproduct( dir[0], dir[1] ) ;    
    jacobian /= jacobian.norm() ;
    
    for(int i=0; i<2; ++i) vardoca[i] = cov[i].similarity(jacobian) ;
    setConsistency(doca*doca/(vardoca[0]+vardoca[1]),1);

    // the new position is simply the weighted average
    HepVector vertexpos = (pos[0]/vardoca[0]+pos[1]/vardoca[1])/(1/vardoca[0] + 1/vardoca[1]) ;
    
    // the covariance matrix is a bit more expensive 
    HepSymMatrix weight[2] ;
    int ierr(0) ;
    for(int i=0; i<2 && !ierr; ++i) {
      // create the projection matrix. it has two rows: the first is
      // the unit vector along the residual (called jacobian). the
      // second is the cross product with the direction.
      HepVector secondrow = crossproduct(jacobian,dir[i]) ;
      secondrow /= secondrow.norm() ;
      HepMatrix P(2,3) ;
      for(int col=1; col<=3; ++col) {
	P(1,col) = jacobian(col) ;
	P(2,col) = secondrow(col) ;
      }
      HepSymMatrix thiscov = cov[i].similarity(P) ;
      thiscov.invert(ierr) ;
      weight[i] = thiscov.similarityT(P) ;
      if(ierr) _status.setFailure(31,"inversion error") ;
    }
    if(_status.success()){
      HepSymMatrix vertexcov = (weight[0] + weight[1]).inverse(ierr) ;
      if(ierr==0) {
    // This 
    //
    //    HepVector vertexpos = vertexcov * ( weight[0]*pos[0] + weight[1]*pos[1]) ;
    // 
    // is the position a vertex fit would give us, if it would treat
    // the tracks as straight lines. Because this is unsuitable for
    // conversions, we will not use it.
    
        _position.setX(vertexpos(1)) ;
        _position.setY(vertexpos(2)) ;
        _position.setZ(vertexpos(3)) ;
        _position.setCovMatrix( vertexcov ) ;
      } else
        _status.setFailure(31,"inversion error") ;
    }
  }
}

TrkPocaVertex::TrkPocaVertex(const TrkRecoTrk* trk, double flt, PdtPid::PidType pid,
			     const BbrPointErr& point, 
			     const PdtEntry* vtype, double precision) : TrkVertex(vtype)
{
  _usedtracks.push_back(trk);  
  const TrkFit* fit = trk->fitResult(pid);
  if(fit != 0)
    _poca = TrkPoca(fit->traj(), flt, point, precision );
  _info[trk] = TrkVtxInfo(_poca.doca(),_poca.flt1(),pid);
  _status = _poca.status();
  if(_status.success()) {

    HepSymMatrix cov[2] ;
    HepVector    pos[2] ;
    double mass = Pdt::lookup(pid)->mass();

    BbrPointErr position = fit->positionErr(_poca.flt1()) ;
    cov[0] = position.covMatrix() ;
    pos[0] = convert<HepPoint>(position) ;
    BbrVectorErr momentum = fit->momentumErr(_poca.flt1() ) ;
    HepVector dir = convert<Hep3Vector>(momentum.unit()) ;
    _p4 += BbrLorentzVectorErr( momentum, mass ) ; 
    
    pos[1] = convert<HepPoint>(point) ;
    cov[1] = point.covMatrix() ;

    // calculate one vector perpendicular to the track. (call it
    // jacobian, in analogy with the trk-trk case above). 
    HepVector residual = pos[1] - pos[0] ;
    HepVector jacobian = residual ;
    const double mindoca = 1e-12 ;
    double norm ;
    if( fabs(norm = jacobian.norm())<mindoca ) {
      // construct the jacobian from any vector perpendicular to the track
      HepVector zaxis(3,0) ; zaxis(3) = 1 ;
      jacobian = crossproduct(zaxis,dir) ;
    }
    jacobian /= jacobian.norm() ;

    // the covariance matrix is a bit more expensive.
    HepSymMatrix weight[2] ;
    int ierr(0);
    // first the track
    HepVector secondrow = crossproduct(jacobian,dir) ;
    secondrow /= secondrow.norm() ;
    HepMatrix P(2,3) ;
    for(int col=1; col<=3; ++col) {
      P(1,col) = jacobian(col) ;
      P(2,col) = secondrow(col) ;
    }
    HepSymMatrix thiscov = cov[0].similarity(P) ;
    thiscov.invert(ierr) ;
    if(ierr==0){
      weight[0] = thiscov.similarityT(P) ;
    // now the point
      weight[1] = cov[1].inverse(ierr) ;
      if(ierr==0){
        // take the weighted average
        HepSymMatrix vertexcov = (weight[0] + weight[1]).inverse(ierr) ;
        if(ierr==0){
          HepVector vertexpos = vertexcov * ( weight[0]*pos[0] + weight[1]*pos[1]) ;
      // we must be able to calculate this more efficiently ... :-(
          double chisq = weight[0].similarity(vertexpos-pos[0]) + weight[1].similarity(vertexpos-pos[1]);
          if(!isnan(chisq)){
            setConsistency(chisq,2);
            _position.setX(vertexpos(1)) ;
            _position.setY(vertexpos(2)) ;
            _position.setZ(vertexpos(3)) ;
            _position.setCovMatrix( vertexcov );
          } else
            _status.setFailure(32,"nan");
        } else
          _status.setFailure(31,"inversion error") ;
      } else
        _status.setFailure(31,"inversion error") ;
    } else 
      _status.setFailure(31,"inversion error") ;
  }
}

TrkPocaVertex::TrkPocaVertex(const TrkPocaVertex& other) :
  TrkVertex(other),_poca(other._poca)
{}

TrkPocaVertex::TrkPocaVertex()
{}

TrkPocaVertex&
TrkPocaVertex::operator = (const TrkPocaVertex& other) {
  if(this != &other){
    TrkVertex::operator = (other);
    _poca = other._poca;
  }
  return *this;
}

TrkPocaVertex*
TrkPocaVertex::clone() const {
  return new TrkPocaVertex(*this);
}

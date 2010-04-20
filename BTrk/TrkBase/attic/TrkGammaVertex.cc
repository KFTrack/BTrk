//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkGammaVertex.cc,v 1.26 2008/05/05 16:48:02 fransham Exp $
//
// Description:
//      A vertexing class specialized in finding photon conversions.
//      The algorithm combines POCA and parallelism in a least-squares fit.
//      Intended for use in the photon conversion list at the Beta-level
//      (SimpleComposition) and at the track-level as a V0 finder (TrkFixup).
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Tomohiko Tanabe, David Brown, LBNL (12/14/2006)
//
// Copyright Information:
//      Copyright (C) 2006              Lawrence Berkeley Laboratory
//
//------------------------------------------------------------------------
#include "BaBar/BaBar.hh"
#include "BaBar/Constants.hh"
#include "TrkBase/TrkGammaVertex.hh"
#include "CLHEP/Geometry/HepPoint.h"
#include "ErrLogger/ErrLog.hh"
#include "TrkBase/TrkDifTraj.hh"
#include "TrkBase/TrkSimpTraj.hh"
#include "TrkBase/TrkFit.hh"
#include "TrkBase/TrkPocaXY.hh"
#include "TrkBase/TrkRecoTrk.hh"
#include "TrkBase/TrkExchangePar.hh"
#include "TrkBase/HelixTraj.hh"
#include "PDT/Pdt.hh"
#include "BField/BField.hh"
#include "HepTuple/Histogram.h"
#include "HepTuple/Tuple.h"
#include "HepTuple/TupleManager.h"

// statics
const double TrkGammaVertex::_maxdoca(2.0);
const double TrkGammaVertex::_maxddip(0.25);
const double TrkGammaVertex::_awconv(1.0/sqr(0.0062));
const double TrkGammaVertex::_awdalitz(1.0/sqr(0.057));


// cache PdtEntry for fast lookups
const PdtEntry *
TrkGammaVertex::getGammaPdt() 
{
  static PdtEntry *gpdt(0);
  if (!gpdt) gpdt = Pdt::lookup(PdtPdg::gamma);
  return gpdt;
}

const PdtEntry *
TrkGammaVertex::getElectronPdt() 
{
  static PdtEntry *epdt(0);
  if (!epdt) epdt = Pdt::lookup(PdtPdg::e_minus);
  return epdt;
}

// Vertex two reco tracks. 
TrkGammaVertex::TrkGammaVertex(const TrkRecoTrk* trk1, const TrkRecoTrk* trk2,
    enum AlgoMode algomode):
  TrkVertex(getGammaPdt()),
  _flt1(0), _flt2(0), _doca(0), _algomode(algomode)
{
  _htrajs[0] = 0;
  _htrajs[1] = 0;

  if (trk1 != 0 &&  trk2 != 0) 
    vertexIt(trk1,trk2);
  else
    _status.setFailure(91,"null track argument");
}

TrkGammaVertex::TrkGammaVertex(const TrkVertex* vtx, enum AlgoMode algomode) :
  TrkVertex(getGammaPdt()),
  _flt1(0), _flt2(0), _doca(0), _algomode(algomode)
{
  _htrajs[0] = 0;
  _htrajs[1] = 0;

  const trkcontainer& trks = vtx->trackList();
  if(trks.size() == 2)
    vertexIt(trks[0],trks[1]);
  else
    _status.setFailure(92,"wrong # tracks");
}

void
TrkGammaVertex::vertexIt(
    const TrkRecoTrk* const & trk1,
    const TrkRecoTrk* const & trk2){
  _usedtracks.push_back(trk1);
  _usedtracks.push_back(trk2);

  const TrkFit* fit1 = trk1->fitResult(PdtPid::electron);
  const TrkFit* fit2 = trk2->fitResult(PdtPid::electron);

// apply pre-selection

  if (fit1 && fit2 && fitSelect(fit1,fit2)) {

// find transverse parallel point.  This overwrites the flightlengths, and
// finds the appropriate local trajectory.
    const TrkDifTraj* trajs[2];
    trajs[0] = &fit1->traj();
    trajs[1] = &fit2->traj();
    double flt[2] = {0.0,0.0};
    HepPoint pp;
    _status = transParallel(trajs,flt,pp);
    _pp = pp;
    if(_status.success()){
      double lflt[2];
      static HepPoint origin(0,0,0);

      for (int itraj=0;itraj<2; itraj++) {
        const HelixTraj* htraj = static_cast<const HelixTraj*>(trajs[itraj]->localTrajectory(flt[itraj],lflt[itraj]));
        if(0 == htraj){
          _status = TrkErrCode(TrkErrCode::fail,101,"Not helix");
          return;
        }
// clone the helix
        _htrajs[itraj] = htraj->clone();
        assert(_htrajs[itraj] != 0);
// move the parameters and covariance to the parallel point.  This makes
// the math much easier
        TrkParams* params = _htrajs[itraj]->parameters();
        TranslateParams pfunc = _htrajs[itraj]->paramFunction();
        pfunc(origin,pp,
              params->parameter(),params->covariance(),
              params->parameter(),params->covariance(),lflt[itraj]);
        _htrajs[itraj]->setPoint(pp);
      }
// test parameters
      double deltaphi0 = _htrajs[1]->parameters()->parameter()[HelixTraj::phi0Index]-
	_htrajs[0]->parameters()->parameter()[HelixTraj::phi0Index];
      if(fabs(deltaphi0)>Constants::pi){
	if(deltaphi0>0)
	  deltaphi0 -= Constants::twoPi;
	else
	  deltaphi0 += Constants::twoPi;
      }
      if(fabs(deltaphi0)>0.0001){
	ErrMsg(error) << "parallel point not parallel! " << endmsg;
	_status = TrkErrCode(TrkErrCode::fail,105,"Not parallel");
      }
// minimize; this refines the parallel position taking into account all the track information.
// It also computes the covariance matrices.  There are 10 parameters: 3 track parameters
// for each track (phi0 drops out) plus the transverse flight length (tf).
// We  work with the average dip angle (dip) and the difference between dip angles (eta),
// related as td1 = td + eta, td2 = td - eta.
      HepVector gpar(9,0);
      HepSymMatrix gcov(9,0);
      double chisq(0.0);
      _status = minimize(_htrajs,gpar,gcov,chisq);
      _gpar = gpar;
      _gcov = gcov;

      if (_status.success()) {
// convert these back to physical space.
        findVertex(trk1->bField(),gpar,gcov,pp);
// fail the rare cases when the position covariance matrix has negative determinant
        HepSymMatrix xcov = _position.covMatrix();
        if (xcov.determinant()<0.) {
          _status.setFailure(205,"inversion error");
        } else {
// set info; these are not exactly correct, as the physics model of the
// conversion re-defines the track parameters.

// test position at updated point; should be the same for both tracks
          double deltal = gpar[d0]*gpar[dphi0];
          double l1 = gpar[tf] + 0.5*deltal;
          double l2 = gpar[tf] - 0.5*deltal;
          double f1 = l1/cos(gpar[dip]+0.5*gpar[ddip]);
          double f2 = l2/cos(gpar[dip]-0.5*gpar[ddip]);
          HepPoint p1 = _htrajs[0]->position(f1);
          HepPoint p2 = _htrajs[1]->position(f2);
// tests
          Hep3Vector diff(p1.x()-p2.x(),p1.y()-p2.y(),p1.z()-p2.z());
// override the position
          Hep3Vector average(p1.x()+p2.x(),p1.y()+p2.y(),p1.z()+p2.z());
          average *= 0.5;
          Hep3Vector vpos(_position.x(),_position.y(),_position.z());
          Hep3Vector ddiff = vpos - average;
          BbrPointErr newpos(HepPoint(average.x(),average.y(),average.z()),
              _position.covMatrix());        
          _position = newpos; 
// update flight parameters (WRT real origin)
          _flt1 = flt[0]+f1;
          _flt2 = flt[1]+f2;
          _doca = p1.distanceTo(p2);
          _info[trk1] = TrkVtxInfo(_doca,_flt1,PdtPid::electron);
          _info[trk2] = TrkVtxInfo(_doca,_flt2,PdtPid::electron);

          setConsistency(std::max(chisq,0.),3);
// still to do: deal correctly with updated daughters (= new parameters)  FIXME!!!
        }
      }
    }
  } else {
    _status.setFailure(107,"tracks fail preselection");
  }
}

TrkGammaVertex::TrkGammaVertex()
{}

TrkGammaVertex::TrkGammaVertex(const TrkGammaVertex& other) :
  TrkVertex(other),
  _flt1(other.flt1()), _flt2(other.flt2()),
  _doca(other._doca),_pp(other._pp),
  _algomode(other._algomode),
  _vxp(other._vxp)
{
  for(int itraj=0;itraj<2;itraj++){
    _htrajs[itraj] = 0;
    if(other._htrajs[itraj] != 0)
      _htrajs[itraj] = other._htrajs[itraj]->clone();
  }
}

TrkGammaVertex::~TrkGammaVertex()
{
  // clean up.
  delete _htrajs[0];
  delete _htrajs[1];
}

TrkGammaVertex&
TrkGammaVertex::operator = (const TrkGammaVertex& other)
{
  if (this != &other) {
    TrkVertex::operator=(other);
    _flt1 = other.flt1();
    _flt2 = other.flt2();
    for(int itraj=0;itraj<2;itraj++){
      delete _htrajs[itraj];
      _htrajs[itraj] = 0;
      if(other._htrajs[itraj] != 0)
        _htrajs[itraj] = other._htrajs[itraj]->clone();
    }
    _doca = other._doca;
    _pp = other._pp;
    _algomode = other._algomode;
    _vxp = other._vxp;
  }
  return *this;
}

TrkGammaVertex*
TrkGammaVertex::clone() const
{
  return new TrkGammaVertex(*this);
}

TrkErrCode
TrkGammaVertex::transParallel(const TrkDifTraj** intrajs, double* flt, HepPoint& pp) const {
// start by getting the local helix interpretation of the fit
  static int ntrajs(2);
  Hep3Vector center[2];
  double rad[2];
  const HelixTraj* trajs[2] = {0,0};
// iterate to find the right local trajectories
  int iter(0);
  static int maxniter(5);
  while(iter<maxniter){
    bool sametrajs(true);
    for(int itraj=0;itraj<ntrajs;itraj++){
// update flightlength and local trajectory
      double lflt;
      const TrkSimpTraj* ltraj = intrajs[itraj]->localTrajectory(flt[itraj],lflt);
      const HelixTraj* htraj = static_cast<const HelixTraj*>(ltraj);
      if(0 == htraj)
        return TrkErrCode(TrkErrCode::fail,101,"Not helix");
      sametrajs &= (htraj == trajs[itraj]);
      trajs[itraj] = htraj;
    }
// escape the iteration loop if we've found the same local trajectories as last time
    if(sametrajs)break;
    iter++;
    for(int itraj=0;itraj<ntrajs;itraj++){
// find the center points in 2-d
      const HelixTraj* htraj = trajs[itraj];
      rad[itraj] = 1.0/htraj->omega();
      double traj_phi0 = htraj->phi0();
      double cphi0 = cos(traj_phi0);
      double sphi0 = sin(traj_phi0);
      center[itraj] = Hep3Vector(-(htraj->d0()+rad[itraj])*sphi0,
				 (htraj->d0()+rad[itraj])*cphi0,0.0);
    }
    Hep3Vector delc;
    if(trajs[0]->omega() > 0)
      delc = center[1]-center[0];
    else
      delc = center[0]-center[1];
// test
    double doca = fabs(delc.mag() - fabs(rad[0]) - fabs(rad[1]));
// remove garbage
    if(doca > _maxdoca)
      return TrkErrCode(TrkErrCode::fail,103,"Bad Doca");

    double phi = atan2(delc.x(),-delc.y());
    double tflt[2];
    Hep3Vector pos[2];
// resolve ambiguity
    for(int itraj=0;itraj<ntrajs;itraj++){
      const HelixTraj* htraj = trajs[itraj];
      double tphi = phi;
      if(fabs(tphi-htraj->phi0()) > Constants::pi){
        if(tphi-htraj->phi0() > 0 ) 
          tphi -= Constants::twoPi;
        else
          tphi += Constants::twoPi;
      }
// compute TRANSVERSE flightlength and space flightlength
      tflt[itraj] = (tphi - htraj->phi0() )*rad[itraj];
      flt[itraj] = tflt[itraj]*sqrt(1.0+sqr(htraj->tanDip()));
      pos[itraj] = Hep3Vector(rad[itraj]*sin(tphi) - (rad[itraj]+htraj->d0())*sin(htraj->phi0()),
			      -rad[itraj]*cos(tphi) + (rad[itraj]+htraj->d0())*cos(htraj->phi0()),
			      0.0);
    }
    double doca2 = (pos[0]- pos[1]).mag();
    if(fabs(doca-doca2) > 0.001)
      return TrkErrCode(TrkErrCode::fail,102,"Ambiguity failure");
// fill in the Z coordinate of the points
    for(int itraj=0;itraj<ntrajs;itraj++){
      const HelixTraj* htraj = trajs[itraj];
      pos[itraj].setZ(htraj->z0() + tflt[itraj]*htraj->tanDip());
    }
// reset to the new point.  Choose the midpoint between the 2 trajectories at parallel
    pp = HepPoint(0.0,0.0,0.0);
    pp += 0.5*(pos[0] + pos[1]);
  }
  return TrkErrCode(TrkErrCode::succeed,iter);
}

bool
TrkGammaVertex::fitSelect(const TrkFit* fit1,const TrkFit* fit2) const {
  bool retval = fit1 != 0 && fit2 != 0;
  if(retval){
// check dip difference (same as theta difference)
    double dtheta = fit1->direction(0.0).theta() - fit2->direction(0.0).theta();
    retval = fabs(dtheta) <  _maxddip;
  }
  return retval;
}

TrkErrCode
TrkGammaVertex::minimize(HelixTraj** htrajs,HepVector& gpar,HepSymMatrix& gcov,double& chisq) const {
// compute initial chisquared (this comes only from parameters)
  HepVector pmeas[2];
  HepSymMatrix twt[2];

// compute the parameter weight matrices for the 2 tracks.
  for(int itraj=0;itraj<2;itraj++){
    pmeas[itraj] = htrajs[itraj]->parameters()->parameter();
    twt[itraj] = htrajs[itraj]->parameters()->covariance();
    int ifail;
    twt[itraj].invert(ifail);
    if(ifail)return TrkErrCode(TrkErrCode::fail,205,"inversion error");
  }
 
// initialize parallel point parameters.  Most are 0
  HepVector gpar0(9,0);
// omega comes from the tracks directly
  gpar0[om1] = pmeas[0][HelixTraj::omegaIndex];
  gpar0[om2] = pmeas[1][HelixTraj::omegaIndex];
// average phi for phi0.  Be sure to treat WRAPPING!!!
  double phi01 = pmeas[0][HelixTraj::phi0Index];
  double phi02 = pmeas[1][HelixTraj::phi0Index];

  if(fabs(phi01-phi02) > Constants::pi){
    if(phi01 > phi02) 
      phi01 -= Constants::twoPi;
    else
      phi01 += Constants::twoPi;
  }
  gpar0[phi0] = 0.5*(phi01+phi02);
// average dip: no wrapping here
  gpar0[dip] = 0.5*(atan(pmeas[0][HelixTraj::tanDipIndex])+atan(pmeas[1][HelixTraj::tanDipIndex]));
  HepVector tpar0[2];
  tpar0[0] = pmeas[0];
  tpar0[1] = pmeas[1];

  fillParameters(gpar0,tpar0);
  double chisq0(0.0);
  HepVector deltaw[2];
  for(int itraj=0;itraj<2;itraj++){
    HepVector tdeltap = tpar0[itraj]-pmeas[itraj];
    chisq0 += twt[itraj].similarity(tdeltap);
    deltaw[itraj] = twt[itraj]*tdeltap;
  }

// derivatives
  HepMatrix dPdQ[2];
  //std::vector<HepSymMatrix> d2PdQ2[2];
  HepSymMatrix d2PdQ2[2][HelixTraj::NHLXPRM];
  HepMatrix null1(HelixTraj::NHLXPRM,ngind,0);
  HepSymMatrix null2(ngind,0);
  for(int itraj=0;itraj<2;itraj++){
    dPdQ[itraj] = null1; 
    for(int ipar=0;ipar<HelixTraj::NHLXPRM;ipar++) {
      d2PdQ2[itraj][ipar] = null2;
    }
  }
// pre-compute some common values
  double cosdip = cos(gpar0[dip]);
  double invcosdip = 1.0/cosdip;
  double invcosdip2 = sqr(invcosdip);
  double tandip = tan(gpar0[dip]);
  double tdipinvcdip2 = tandip*invcosdip2;
  double domega = gpar0[om1]-gpar0[om2];

// d0
  dPdQ[0][HelixTraj::d0Index][d0] = 1.0;
  dPdQ[1][HelixTraj::d0Index][d0] = 1.0;
  
  d2PdQ2[0][HelixTraj::d0Index][dphi0][tf] = -0.5;
  d2PdQ2[1][HelixTraj::d0Index][dphi0][tf] = 0.5;
  
  d2PdQ2[0][HelixTraj::d0Index][tf][tf] = -0.5*domega;
  d2PdQ2[1][HelixTraj::d0Index][tf][tf] = 0.5*domega;
//phi0
  dPdQ[0][HelixTraj::phi0Index][phi0] = 1.0;
  dPdQ[1][HelixTraj::phi0Index][phi0] = 1.0;
  
  dPdQ[0][HelixTraj::phi0Index][dphi0] = 0.5;
  dPdQ[1][HelixTraj::phi0Index][dphi0] = -0.5;
//omega
  dPdQ[0][HelixTraj::omegaIndex][om1] = 1.0;
  dPdQ[1][HelixTraj::omegaIndex][om2] = 1.0;
//z0
  dPdQ[0][HelixTraj::z0Index][z0] = 1.0;
  dPdQ[1][HelixTraj::z0Index][z0] = 1.0;
  
  d2PdQ2[0][HelixTraj::z0Index][ddip][tf] = -0.5*invcosdip2;
  d2PdQ2[1][HelixTraj::z0Index][ddip][tf] = 0.5*invcosdip2;

  d2PdQ2[0][HelixTraj::z0Index][d0][dphi0] = -0.5*tandip;
  d2PdQ2[1][HelixTraj::z0Index][d0][dphi0] = 0.5*tandip;

//tandip
    
  dPdQ[0][HelixTraj::tanDipIndex][dip] = invcosdip2;
  dPdQ[1][HelixTraj::tanDipIndex][dip] = invcosdip2;
  
  dPdQ[0][HelixTraj::tanDipIndex][ddip] = 0.5*invcosdip2;
  dPdQ[1][HelixTraj::tanDipIndex][ddip] = -0.5*invcosdip2;

  d2PdQ2[0][HelixTraj::tanDipIndex][dip][dip] = 2*tdipinvcdip2;
  d2PdQ2[1][HelixTraj::tanDipIndex][dip][dip] = 2*tdipinvcdip2;
  
  d2PdQ2[0][HelixTraj::tanDipIndex][ddip][ddip] = 0.5*tdipinvcdip2;
  d2PdQ2[1][HelixTraj::tanDipIndex][ddip][ddip] = 0.5*tdipinvcdip2;
  
  d2PdQ2[0][HelixTraj::tanDipIndex][dip][ddip] = tdipinvcdip2;
  d2PdQ2[1][HelixTraj::tanDipIndex][dip][ddip] = -tdipinvcdip2;

// angular terms.  Only 1st order terms as no measurement
// first phi

  HepVector dphidQ(ngind,0);

  dphidQ[dphi0] = cosdip;
  dphidQ[tf] =  cosdip*domega;

// theta

  HepVector dthetadQ(ngind,0);
  dthetadQ[ddip]=1.0;

// now minimize chisq; linear terms
  HepVector beta(ngind,0);
  HepSymMatrix gamma(ngind,0);
  for(int itraj=0;itraj<2;itraj++){
    HepMatrix dpdq = dPdQ[itraj].T();
    beta += dpdq*deltaw[itraj];
    gamma += twt[itraj].similarity(dpdq);
// now 2nd order terms.  The loop is because there's no 3rd rank tensor in CLHEP.
    for(int ipar=0;ipar<HelixTraj::NHLXPRM;ipar++) {
      gamma += d2PdQ2[itraj][ipar]*deltaw[itraj][ipar];
    }
  }
// Add angular difference terms; there is no contribution to beta since the
// nominal angular difference = 0
  double aweight(0);
  if (_algomode == Conversion)
    aweight = _awconv;
  else
    aweight = _awdalitz;

  gamma += aweight*vT_times_v(dphidQ);
  gamma += aweight*vT_times_v(dthetadQ);

  /*
  // rescale gamma first
  const double scale[ngind] = {
    3e4, 3e4, 5e2, 1e3, 1e2, 1e2, 2e3, 2e3, 1e2
  };
  for (int i=0; i<ngind; ++i) {
    for (int j=i; j<ngind; ++j) {
      gamma[i][j] *= scale[i]*scale[j];
    }
  }
  */
  /* original values for chisq scan
  static double delta[ngind] = {
    0.0004, //om1
    0.0004, //om2
    0.002, //phi0
    0.002, //dip
    0.005, //d0
    0.01, //z0
    0.002, //dphi0
    0.002, //ddip
    0.1    //tf
  };
  */
  // "hand scan" of second derivatives
  static double delta[ngind] = {
    0.00004, //om1
    0.00004, //om2
    0.0002, //phi0
    0.0002, //dip
    0.0005, //d0
    0.001, //z0
    0.0002, //dphi0
    0.0002, //ddip
    0.01    //tf
  };

  HepVector secondDerivs(ngind);
  const int nSample = 5;
  double sampleX[nSample];
  double sampleY[nSample];

  for (int iind=0; iind<ngind; ++iind) {
    for (int j=0; j<nSample; ++j) {
      double n = j-(nSample-1.)/2.;
      double chisq(0.);
      HepVector gpar = gpar0;
      gpar[iind] += n*delta[iind];
      HepVector tpar[2];
      tpar[0] = pmeas[0];
      tpar[1] = pmeas[1];
      fillParameters(gpar,tpar);
      for(int itraj=0;itraj<2;itraj++){
        HepVector tdeltap = tpar[itraj]-pmeas[itraj];
        chisq += twt[itraj].similarity(tdeltap);
      }
      double dphi = gpar[dphi0] + gpar[tf]*(gpar[om1]-gpar[om2])
        + 0.5*gpar[d0]*gpar[dphi0]*(gpar[om1]+gpar[om2]);

      chisq += sqr(dphi*cos(gpar[dip]))*aweight;
      chisq += sqr(gpar[ddip])*aweight;

      sampleX[j] = gpar[iind]-gpar0[iind];
      sampleY[j] = chisq;
    }
    double sumXsq(0.);
    double sumX4(0.);
    double sumY(0.);
    double sumXsqY(0.);
    for (int j=0; j<nSample; ++j) {
      double Xsq = sampleX[j]*sampleX[j];
      sumXsq += Xsq;
      sumXsqY += Xsq*sampleY[j];
      sumX4 += Xsq*Xsq;
      sumY += sampleY[j];
    }
    sumXsq /= (double)nSample;
    sumX4 /= (double)nSample;
    sumY /= (double)nSample;
    sumXsqY /= (double)nSample;
    HepMatrix mScan(2,2);
    mScan(1,1) = 1.;
    mScan(1,2) = sumXsq;
    mScan(2,1) = sumXsq;
    mScan(2,2) = sumX4;
    HepVector mVec(2);
    mVec(1) = sumY;
    mVec(2) = sumXsqY;
    int ifail(0);
    mScan.invert(ifail);
    if(ifail)return TrkErrCode(TrkErrCode::fail,205,"inversion error");
    HepVector resVec = mScan*mVec;
    secondDerivs[iind] = resVec(2);
  }

  bool saddlePoint(false); // flag to check if chisq is at a saddle point in tf
  if (gamma[tf][tf]<=0. || gamma[tf][tf] > 100.) {
    saddlePoint = true;
  } else {
    // check for weird correlations between tf and other variables
    for (int iind=0; iind<tf; ++iind) {
      HepSymMatrix msub(2);
      msub(1,1) = gamma[iind][iind];
      msub(1,2) = gamma[iind][tf];
      msub(2,2) = gamma[tf][tf];
      saddlePoint = msub.determinant()<=0;
      if (saddlePoint) break;
    }
  }

  if (!saddlePoint) {
    // check for weird correlations between tf, dphi0, and ddip
    HepSymMatrix msub(3);
    msub(1,1) = gamma[dphi0][dphi0];
    msub(1,2) = gamma[dphi0][ddip];
    msub(1,3) = gamma[dphi0][tf];
    msub(2,2) = gamma[ddip][ddip];
    msub(2,3) = gamma[ddip][tf];
    msub(3,3) = gamma[tf][tf];
    saddlePoint = msub.determinant()<=0;
  }

  if (saddlePoint) {
    //ErrMsg(warning) << "saddle point in tf" << endmsg;
    for (int i=0; i<ngind; ++i) {
      gamma[tf][i] = 0.;
    }
    // set tf error squared to (0.5cm)^2
    gamma[tf][tf] = 4.; // = 1/(0.5)^2
    // do not update the longitudinal position
    beta[tf]=0.;
  }

// now solve
  int ifail(0);
  gcov = gamma.inverse(ifail);
  if(ifail)return TrkErrCode(TrkErrCode::fail,205,"inversion error");

// check for saddle points in tf
  /*
  bool saddlePoint(false); // flag to check if chisq is at a saddle point in tf
  if (gcov[tf][tf] < 0.001) {
    //ErrMsg(warning) << "saddle point in tf" << endmsg;
    saddlePoint = true;
    for (int i=0; i<ngind; ++i) {
      gcov[tf][i] = 0;
    }
    // set tf error squared to (0.5cm)^(-2)
    gcov[tf][tf] = 0.25;
  }
  */

  HepVector delgpar = gcov*beta;
  gpar = gpar0-delgpar;

// compute chisq
  chisq = 0.0;
  HepVector tpar[2];
  tpar[0] = pmeas[0];
  tpar[1] = pmeas[1];

  fillParameters(gpar,tpar);
  for(int itraj=0;itraj<2;itraj++){
    HepVector tdeltap = tpar[itraj]-pmeas[itraj];
    chisq += twt[itraj].similarity(tdeltap);
  }
// angular terms too.  Be careful to compute phi difference including ALL TERMS!
  double dphi = gpar[dphi0] + gpar[tf]*(gpar[om1]-gpar[om2]) 
    + 0.5*gpar[d0]*gpar[dphi0]*(gpar[om1]+gpar[om2]);

  chisq += sqr(dphi*cos(gpar[dip]))*aweight;
  chisq += sqr(gpar[ddip])*aweight;

// update trajectory parameters

  htrajs[0]->parameters()->parameter() = tpar[0];
  htrajs[1]->parameters()->parameter() = tpar[1];

// check: linear update to chisq0
//  double chisq1 = chisq0 -gcov.similarity(beta);

  return TrkErrCode(TrkErrCode::succeed,1);
}

void
TrkGammaVertex::fillParameters(const HepVector& gpar, HepVector* trkpars) const {

  double tandip = tan(gpar[dip]);
  double cosdip = cos(gpar[dip]);
  double cosdip2 = cosdip*cosdip;
  double delom = gpar[om1] - gpar[om2];
  double deltad0 = -(gpar[dphi0]+0.5*delom*gpar[tf])*gpar[tf];
  double deltaz0 = -gpar[ddip]*gpar[tf]/cosdip2 - gpar[dphi0]*gpar[d0]*tandip;

  trkpars[0][HelixTraj::d0Index] = gpar[d0] + 0.5*deltad0;
  trkpars[1][HelixTraj::d0Index] = gpar[d0] - 0.5*deltad0;
  
// check for wrapping on phi.  This relies on the vectors having been initialized with the 
// right-wrap parameters
  
  double phi1 = gpar[phi0] + 0.5*gpar[dphi0];
  double phi2 = gpar[phi0] - 0.5*gpar[dphi0];
  if(abs(phi1 - trkpars[0][HelixTraj::phi0Index]) < Constants::pi)
      trkpars[0][HelixTraj::phi0Index] = phi1; 
  else {
    if(phi1 > trkpars[0][HelixTraj::phi0Index]) 
      trkpars[0][HelixTraj::phi0Index] = phi1 - Constants::twoPi;
    else
      trkpars[0][HelixTraj::phi0Index] = phi1 + Constants::twoPi;
  }
  
  if(abs(phi2- trkpars[1][HelixTraj::phi0Index]) < Constants::pi)
    trkpars[1][HelixTraj::phi0Index] = phi2;
  else {
    if(phi2 > trkpars[1][HelixTraj::phi0Index])
      trkpars[1][HelixTraj::phi0Index] = phi2 - Constants::twoPi;
    else  
      trkpars[1][HelixTraj::phi0Index] = phi2 + Constants::twoPi;
  } 

  trkpars[0][HelixTraj::omegaIndex] = gpar[om1];
  trkpars[1][HelixTraj::omegaIndex] = gpar[om2];
  
  trkpars[0][HelixTraj::z0Index] = gpar[z0] + 0.5*deltaz0;
  trkpars[1][HelixTraj::z0Index] = gpar[z0] - 0.5*deltaz0;

  trkpars[0][HelixTraj::tanDipIndex] = tan(gpar[dip] + 0.5*gpar[ddip]);
  trkpars[1][HelixTraj::tanDipIndex] = tan(gpar[dip] - 0.5*gpar[ddip]);

}


void
TrkGammaVertex::findVertex(const BField& bfield, const HepVector& gpar,
    const HepSymMatrix& gcov, const HepPoint& pp)
{
// compute the change in position
  double cp0 = cos(gpar[phi0]);
  double sp0 = sin(gpar[phi0]);
  double td = tan(gpar[dip]);
  double invcdip2 = 1.0+sqr(td);
  Hep3Vector dr(cp0*gpar[tf]-sp0*gpar[d0],
		sp0*gpar[tf]+cp0*gpar[d0],
		gpar[z0]+td*gpar[tf]);
// don't update point for Dalitz
  HepPoint newr = pp;
  if(_algomode == Conversion)
    newr += dr;
// derivatives of this
  HepMatrix drdq(3,ngind,0);

  drdq[xind][phi0] = -sp0*gpar[tf] - cp0*gpar[d0];
  drdq[yind][phi0] = cp0*gpar[tf] - sp0*gpar[d0];

  drdq[zind][dip] = invcdip2*gpar[tf];

  drdq[xind][tf] = cp0;
  drdq[yind][tf] = sp0;
  drdq[zind][tf] = td;

  drdq[xind][d0] = -sp0;
  drdq[yind][d0] = cp0;

  drdq[zind][z0] = 1;

// propagate the errors
  HepSymMatrix sigr = gcov.similarity(drdq);
// put it together
  _position = BbrPointErr(newr,sigr);

// now momentum.

  double R[2] = {1.0/fabs(gpar[om1]),1.0/fabs(gpar[om2])};
  double dsign[2] = {1.0,-1.0};
  double bf = BField::cmTeslaToGeVc*bfield.bFieldNominal();
  double emass2 = sqr(getElectronPdt()->mass());
  HepLorentzVector g4;
  HepSymMatrix sigp(3,0);
  int om[2] = {om1,om2};
  HepMatrix dpdq(3,ngind,0);
  for(unsigned itraj=0;itraj<2;itraj++){
// pt is given just by the curvature
    double pt = bf*R[itraj];
    double tflt = gpar[tf] + 0.5*gpar[d0]*gpar[dphi0]*dsign[itraj];
    double phi = gpar[phi0] + 0.5*gpar[dphi0]*dsign[itraj] + tflt*gpar[om[itraj]];
    double dipangle = gpar[dip] + 0.5*gpar[ddip]*dsign[itraj];
    double cosphi = cos(phi);
    double sinphi = sin(phi);
    double tandip= tan(dipangle);
    double invcosdip2 = 1.0+sqr(tandip);
    Hep3Vector ep = pt*Hep3Vector(cosphi,sinphi,tandip);
    double ee = sqrt(ep.mag2() + emass2);
    g4 += HepLorentzVector(ep,ee);
// now derivatives

    dpdq[xind][om[itraj]] = -ep[xind]/gpar[om[itraj]] - pt*sinphi*tflt;
    dpdq[yind][om[itraj]] = -ep[yind]/gpar[om[itraj]] + pt*cosphi*tflt;
    dpdq[zind][om[itraj]] = -ep[zind]/gpar[om[itraj]];

    dpdq[xind][phi0] += -pt*sinphi;
    dpdq[yind][phi0] += pt*cosphi;

    dpdq[xind][dphi0] += -0.5*pt*sinphi*dsign[itraj]*(1.0+gpar[d0]*gpar[om[itraj]]);
    dpdq[yind][dphi0] += 0.5*pt*cosphi*dsign[itraj]*(1.0+gpar[d0]*gpar[om[itraj]]);

    dpdq[xind][tf] += -pt*sinphi*gpar[om[itraj]];
    dpdq[yind][tf] += pt*cosphi*gpar[om[itraj]];

    dpdq[xind][d0] += -0.5*dsign[itraj]*pt*sinphi*gpar[om[itraj]]*gpar[dphi0];
    dpdq[yind][d0] += 0.5*dsign[itraj]*pt*cosphi*gpar[om[itraj]]*gpar[dphi0];

    dpdq[zind][dip] += pt*invcosdip2;
    dpdq[zind][ddip] += 0.5*dsign[itraj]*pt*invcdip2;

  }
// propagate the errors
  
  sigp = gcov.similarity(dpdq);

  _p4 = BbrLorentzVectorErr(g4,sigp);

  _vxp = drdq*gcov*dpdq.T();  

}

// do a scan of the chisq with respect to the supplied true position
const HepHistogram*
TrkGammaVertex::fillChisq(
    HepTupleManager* manager, const char* name,
    HepTuple* tuple, const HepPoint& tpos,
    int gind1, int gind2
    ) const {
  assert(_status.success());
  assert(manager);
  assert(tuple);
  assert(_htrajs[0]);
  assert(_htrajs[1]);

  Hep3Vector delpos = tpos-_pp;

// compute initial chisquared (this comes only from parameters)
  HepVector pmeas[2];
  HepSymMatrix twt[2];

// compute the parameter weight matrices for the 2 tracks.
  for(int itraj=0;itraj<2;itraj++){
    pmeas[itraj] = _htrajs[itraj]->parameters()->parameter();
    twt[itraj] = _htrajs[itraj]->parameters()->covariance();
    int ifail;
    twt[itraj].invert(ifail);
    if(ifail)return 0;
  }

// initialize parallel point parameters.  Most are 0
  HepVector gpar0(9,0);
// omega comes from the tracks directly
  gpar0[om1] = pmeas[0][HelixTraj::omegaIndex];
  gpar0[om2] = pmeas[1][HelixTraj::omegaIndex];
// average phi for phi0.  Be sure to treat WRAPPING!!!
  double phi01 = pmeas[0][HelixTraj::phi0Index];
  double phi02 = pmeas[1][HelixTraj::phi0Index];

  if(fabs(phi01-phi02) > Constants::pi){
    if(phi01 > phi02) 
      phi01 -= Constants::twoPi;
    else
      phi01 += Constants::twoPi;
  }
  gpar0[phi0] = 0.5*(phi01+phi02);
// average dip: no wrapping here
  gpar0[dip] = 0.5*(atan(pmeas[0][HelixTraj::tanDipIndex])+atan(pmeas[1][HelixTraj::tanDipIndex]));
  if (tuple) {
    double cdip = cos(gpar0[dip]);
    double sdip = sin(gpar0[dip]);
    double cphi = cos(gpar0[phi0]);
    double sphi = sin(gpar0[phi0]);
    Hep3Vector drdL(cphi*cdip,sphi*cdip,sdip);
    Hep3Vector drdL2(cphi,sphi,0);
    double trueL = drdL.dot(delpos);
    double trueL2 = drdL2.dot(delpos);

    tuple->column("L",(float)trueL);
    tuple->column("L2",(float)trueL2);
    tuple->column("om1",(float)gpar0[om1]);
    tuple->column("om2",(float)gpar0[om2]);
    tuple->column("phi0",(float)gpar0[phi0]);
    tuple->column("dip",(float)gpar0[dip]);
    tuple->column("d0",(float)gpar0[d0]);
    tuple->column("z0",(float)gpar0[z0]);
    tuple->column("dphi0",(float)gpar0[dphi0]);
    tuple->column("ddip",(float)gpar0[ddip]);
    tuple->column("tf",(float)gpar0[tf]);
    tuple->dumpData();
  }

  static double delta[ngind] = {
    0.0004, //om1
    0.0004, //om2
    0.002, //phi0
    0.002, //dip
    0.005, //d0
    0.01, //z0
    0.002, //dphi0
    0.002, //ddip
    0.1    //tf
  };

  const static int nWidth = 50;
  double min1 = gpar0[gind1]-nWidth*delta[gind1];
  double max1 = gpar0[gind1]+nWidth*delta[gind1];
  double min2 = gpar0[gind2]-nWidth*delta[gind2];
  double max2 = gpar0[gind2]+nWidth*delta[gind2];

  HepHistogram* hist = manager->histogram(name,101,min1,max1,101,min2,max2);

  double aweight(0);
  if (_algomode == Conversion)
    aweight = _awconv;
  else
    aweight = _awdalitz;

  for (int i=-nWidth; i<=nWidth; ++i) {
    for (int j=-nWidth; j<=nWidth; ++j) {
      double chisq = 0.0;
      HepVector gpar = gpar0;
      gpar[gind1] += i*delta[gind1];
      gpar[gind2] += j*delta[gind2];

      HepVector tpar[2];
      tpar[0] = pmeas[0];
      tpar[1] = pmeas[1];

      fillParameters(gpar,tpar);
      for(int itraj=0;itraj<2;itraj++){
        HepVector tdeltap = tpar[itraj]-pmeas[itraj];
        chisq += twt[itraj].similarity(tdeltap);
      }
// angular terms too.  Be careful to compute phi difference including ALL TERMS!
      double dphi = gpar[dphi0] + gpar[tf]*(gpar[om1]-gpar[om2])
        + 0.5*gpar[d0]*gpar[dphi0]*(gpar[om1]+gpar[om2]);

      chisq += sqr(dphi*cos(gpar[dip]))*aweight;
      chisq += sqr(gpar[ddip])*aweight;

      hist->accumulate2(gpar[gind1],gpar[gind2],chisq);
    }
  }
  return hist;
}

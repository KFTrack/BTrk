// ------------------------------------------------------------------------------
// File and Version Information:
//      $Id: KalParams.hh,v 1.25 2006/04/24 18:53:06 brownd Exp $
//
//  Description: KalParams
//  define the parameter aspect of a KalmanSite: this is just a parameter vector
//  and a covariance matrix
//
// Copyright Information:
//	Copyright (C) 1997	Lawrence Berkeley Laboratory
//
//  Authors: Dave Brown, 8/22/97
//------------------------------------------------------------------------------
#ifndef KALPARAMS_HH
#define KALPARAMS_HH
#include <iostream>
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Matrix/Vector.h"
#include "TrkBase/TrkSimpTraj.hh"
class KalWeight;
class TrkParams;
class KalParams{
public:
  KalParams();
  KalParams(int); // specify the dimension
  KalParams(const HepVector&, const HepSymMatrix&);
  KalParams(const KalWeight&);
  KalParams(const KalParams&);
  KalParams(const TrkParams&); // construct from track parameters
  ~KalParams(){;}
  int nPar() const { return _pvector.num_row(); }
  const HepVector& parameterVector() const { return _pvector; }
  const HepSymMatrix& covarianceMatrix() const { return _pcov; }

  HepVector& parameterVector() { return _pvector; }
  HepSymMatrix& covarianceMatrix() { return _pcov; }

  KalParams& operator = (const KalParams& other) {
    if (&other != this){
      _pvector = other._pvector;
      _pcov = other._pcov;
      _status = other._status;
    }
    return *this; }
  KalParams& operator = (const KalWeight& other);
// addition operation, used to add the effect of 2 parameters
  KalParams& operator += (const KalParams& other) {
    if(matrixOK() && other.matrixOK()){
      _pvector += other._pvector;
      _pcov += other._pcov;
    } else
      _status = -1;
    return *this; }
// subtraction operation, used to _remove_ the effect of one parameter from the other
  KalParams& operator -= (const KalParams& other) {
    if(matrixOK() && other.matrixOK()){
      _pvector -= other._pvector;
      _pcov -= other._pcov;
    } else
      _status = -1;
    return *this; }
// Same for HepVector (changes paremeters, not covariance)
  KalParams& operator += (const HepVector& transvec) {
    if(matrixOK())
      _pvector += transvec;
    return *this; }
  KalParams& operator -= (const HepVector& transvec) {
    if(matrixOK())
      _pvector -= transvec;
    return *this; }
// define negation as applying only to the vector, not covariance
  KalParams& operator -() {
    _pvector = -_pvector;
    return *this;
  }
// multiplication by a float blows up the covariance matrix
  KalParams& operator *= (float factor) {
    if(matrixOK())
      _pcov *= factor;
    return *this; }
// explicitly name some functions used in processing
  void addEffect(const KalParams& other) {
    if(matrixOK() && other.matrixOK()){
      _pvector += other._pvector;
      _pcov += other._pcov;
    } else
      _status = -1;
  }
//
  void subEffect(const KalParams& other) {
    if(matrixOK() && other.matrixOK()){
      _pvector -= other._pvector;
      _pcov += other._pcov;
    } else
      _status = -1;
  }
// test matrix status after inversion
  bool matrixOK() const { return _status == 0; }
  int status() const { return _status; }
// invert and overwrite
  KalParams& invert(const TrkSimpTraj* traj);

// translation function: note that only parameters (not weights) can be
// spatially translated
  KalParams& translatePoint(TranslateParams tfunc,
			    const HepPoint& oldpoint,
			    const HepPoint& newpoint,
			    double fltlen){
    tfunc(oldpoint,newpoint,_pvector,_pcov,_pvector,_pcov,fltlen);
    return *this; }
// compare 2 sites in terms of a chisq, optionally only for a subset of parameters
  double chisq(const KalParams& other,bool* tparams=0) const;
// diagonalize the covariance matrix
  void diagonalize();
// Return equivalent Track parameters (yes, by value)
  TrkParams trackParameters() const {
    return TrkParams(_pvector,_pcov); }
// printout
  void printAll(std::ostream& os=std::cout) const {
    print(os);
  }
  void print(std::ostream& os=std::cout) const;
// invalidate to force re-computation
  void invalidate() { _status = 1;}
private:
  HepVector _pvector;
  HepSymMatrix _pcov;
  int _status;
};
#endif

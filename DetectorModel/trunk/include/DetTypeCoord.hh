// ------------------------------------------------------------------------------
// File and Version Information:
//      $Id: DetTypeCoord.hh,v 1.4 2002/12/30 15:44:30 dbrown Exp $
//  Description:
//  Trivial classes for define 2 and 3 dimensional (abstract) coordinates
//
// Copyright Information:
//	Copyright (C) 1996	Lawrence Berkeley Laboratory
//
//  Authors: Dave Brown, 11/5/96
//------------------------------------------------------------------------------
//
#ifndef TYPECOORD_HH
#define TYPECOORD_HH
#include <math.h>
#include "BaBar/BaBar.hh"
#include "ErrLogger/ErrLog.hh"
struct TypeCoord{
  virtual double & operator [](int) = 0; // return an element
  virtual const double & operator [](int) const = 0; // return an element
  virtual const double* array() const = 0; // return the array
  virtual double* array() = 0; // return the array
  virtual TypeCoord* copyOf() = 0; // copy function
  double magnitude() const {
    return sqrt(magnitudesqr());
  }  // magnitude function
  virtual double magnitudesqr() const = 0; // magnitude squared
//  more functions
  virtual bool operator == (const TypeCoord& other) const = 0;
  virtual int nDimensions() const = 0;
  virtual ~TypeCoord() {;}
};
//
//  Specific coordinates for 2-D and 3-D types
//
struct ThreeDCoord: public TypeCoord{
  ThreeDCoord(const double& u,const double& v,const double& w){
    _uvw[0] = u; _uvw[1] = v; _uvw[2] = w;}
  ThreeDCoord(double *uvw){
    _uvw[0] = uvw[0]; _uvw[1] = uvw[1]; _uvw[2] = uvw[2];}
  ThreeDCoord(){;}
  ThreeDCoord(const ThreeDCoord& other){
	_uvw[0] = other._uvw[0];
	_uvw[1] = other._uvw[1];
	_uvw[2] = other._uvw[2];       }

  double _uvw[3];
  double & operator [](int icoord){
    return _uvw[icoord]; }
  const double & operator [](int icoord) const {
    return _uvw[icoord]; }
  double* array(){return _uvw;}
  const double* array() const{return _uvw;}
  TypeCoord* copyOf(){
    return (TypeCoord*) new ThreeDCoord(*this);}
  double magnitudesqr() const {
    return _uvw[0]*_uvw[0]+_uvw[1]*_uvw[1]+_uvw[2]*_uvw[2];}
  int nDimensions() const { return 3; }
  bool operator == (const TypeCoord& other) const {
    return other.nDimensions() == 3 && _uvw[0] == other[0]
      && _uvw[1] == other[1] &&  _uvw[2] == other[2]; }
  ThreeDCoord& operator = (const ThreeDCoord& other) {
    _uvw[0] = other._uvw[0];
    _uvw[1] = other._uvw[1];
    _uvw[2] = other._uvw[2];
    return *this;
  }
};
///
struct TwoDCoord: public TypeCoord{
  TwoDCoord(const double& u,const double& v){
    _uv[0] = u; _uv[1] = v;}
  TwoDCoord(double *uv){
    _uv[0] = uv[0]; _uv[1] = uv[1];}
  TwoDCoord(){;}
  TwoDCoord(const TwoDCoord& other){
	_uv[0] = other._uv[0];
	_uv[1] = other._uv[1]; }
  double _uv[2];
  double & operator [](int icoord){
    return _uv[icoord]; }
  const double & operator [](int icoord) const {
    return _uv[icoord]; }
  double* array(){return _uv;}
  const double* array() const {return _uv;}
  TypeCoord* copyOf(){
    return (TypeCoord*) new TwoDCoord(*this);}
  double magnitudesqr() const {
    return _uv[0]*_uv[0]+_uv[1]*_uv[1];}
  int nDimensions() const { return 2; }
  bool operator == (const TypeCoord& other) const {
    return other.nDimensions() == 2 && _uv[0] == other[0] &&
      _uv[1] == other[1]; }
  TwoDCoord& operator = (const TwoDCoord& other) {
    _uv[0] = other._uv[0];
    _uv[1] = other._uv[1];
    return *this;
  }
};
#endif

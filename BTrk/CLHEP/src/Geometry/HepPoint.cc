// -*- C++ -*-
// $Id: HepPoint.cc 379 2009-07-06 11:09:56Z stroili $
// ---------------------------------------------------------------------------
//
// This file is part of the CLHEP - a Class Library for High Energy Physics.
//
// This is the implementation of the HepPoint class.
//
// .SS History
// Author:Victoria Novotny(LBL)

#ifdef __GNUC__
#pragma implementation
#endif

#include <cstdio>

#include "CLHEP/config/sstream.h"
#include <ctype.h>

#include "CLHEP/config/TemplateFunctions.h"

#include "CLHEP/Geometry/HepPoint.h"
#include "CLHEP/Vector/ThreeVector.h"
// #include "CLHEP/String/Strings.h"
#include <string>
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Geometry/Translation.h"
#include "CLHEP/Geometry/Transformation.h"

HepPoint::HepPoint(double x, double y, double z)
    : dx(x), dy(y),dz(z){}

HepPoint::HepPoint(const HepPoint & p) : dx(p.x()), dy(p.y()), dz(p.z()) {}


HepPoint & HepPoint::operator = (const HepPoint & p) {
    dx = p.x();
    dy = p.y();
    dz = p.z();
    return *this;
}

HepPoint HepPoint::operator - () const {
    HepPoint q(-dx,-dy,-dz);
    return q;
}

HepPoint & HepPoint::operator += (const Hep3Vector & p){
    dx += p.x();
    dy += p.y();
    dz += p.z();
    return *this;
}

HepPoint & HepPoint::operator -= (const Hep3Vector & p){
    dx -= p.x();
    dy -= p.y();
    dz -= p.z();
    return *this;
}

double HepPoint::mag2() const {
    return dx*dx + dy*dy + dz*dz;
}

double HepPoint::mag() const {
    return sqrt(mag2());
}

double HepPoint::phi() const{
    return dx == 0.0 && dy == 0.0 ? 0.0 : atan2 (dy, dx);
}

double HepPoint::theta() const{
    return dx == 0.0 && dy == 0.0 && dz == 0.0 ? 0.0 :
           atan2(sqrt(dx*dx+dy*dy),dz);
}

double HepPoint::cosTheta() const {
    float ptot = mag();
    return ptot == 0.0 ? 1.0 : dz/ptot;
}

double HepPoint::distanceTo(const HepPoint& other) const {
  return sqrt( sqr(dx - other.dx) +
	       sqr(dy - other.dy) +
	       sqr(dz - other.dz));
}

HepPoint & HepPoint::operator *= (const HepRotation & m){
    *this = m * (*this);
    return *this;
}

HepPoint & HepPoint::transform(const HepRotation & m){
    *this = m * (*this);
    return *this;
}

HepPoint & HepPoint::transform(const HepTranslation & m) {
    *this += m.trans_vec();
    return *this;
}

HepPoint & HepPoint::transform(const HepTransformation & m) {
    transform(m.rot_mat()).transform(m.tra_vec());
    return *this;
}

void HepPoint::rotate(float psi, const Hep3Vector & axis){
    HepRotation trans;
    trans.rotate(psi, axis);
    operator*=(trans);
    return;
}

void HepPoint::rotateX(float psi){
    double sp = sin(psi);
    double cp = cos(psi);
    float py = dy;
    float pz = dz;
    dy = cp * py - sp * pz;
    dz = sp * py + cp * pz;
    return;
}

void HepPoint::rotateY(float theta){
    double st = sin(theta);
    double ct = cos(theta);
    float px = dx;
    float pz = dz;
    dx = ct * px + st * pz;
    dz = -st * px + ct * pz;
    return;
}

void HepPoint::rotateZ(float phi){
    double sp = sin(phi);
    double cp = cos(phi);
    float px = dx;
    float py = dy;
    dx = cp * px - sp * py;
    dy = sp * px + cp * py;
    return;
}

void HepPoint::setMag(double ma){
    double th = theta();
    double ph = phi();

    setX(ma*sin(th)*cos(ph));
    setY(ma*sin(th)*sin(ph));
    setZ(ma*cos(th));
}

void HepPoint::setTheta(double th){
    double ma = mag();
    double ph = phi();

    setX(ma*sin(th)*cos(ph));
    setY(ma*sin(th)*sin(ph));
    setZ(ma*cos(th));
}


void HepPoint::setPhi(double ph){
    double ma = mag();
    double th = theta();

    setX(ma*sin(th)*cos(ph));
    setY(ma*sin(th)*sin(ph));
    setZ(ma*cos(th));
}

int operator == (const HepPoint & p, const HepPoint & q){
  if ((q.x() == p.x()) && (q.y() == p.y()) && (q.z() == p.z())) {return 1;}
  return 0;
}

int operator != (const HepPoint & p, const HepPoint & q){
  if ((q.x() != p.x()) || (q.y() != p.y()) || (q.z() != p.z())) {return 1;}
  return 0;
}

std::ostream & operator << (std::ostream & s, const HepPoint & q) {
  return s << "(" << q.x() << "," << q.y() << "," << q.z() << ")";
}

std::istream & operator >> (std::istream & s, HepPoint & q) {
  // Parse:
  std::string cleanString;
  char c;
  int numRead = 0;
  bool previosIsSpace = false;

  while(EOF != (c = s.get()) && c != ')' && s.good()) {
    // ignore initial bracket or spaces
    if (!cleanString.empty() || ('(' != c && 0 == isspace(c))){
      if (',' == c){          // replace ',' with ' '
	c = ' ';
      }
      cleanString += c;     // take c
      if (isspace(c)){      // count # of white spaces (excluding initials)
	if (false == previosIsSpace) {
	  ++numRead;
	}
	if (numRead >= q.SIZE){  // read enough
	  break;
	}
      }
      else {
	previosIsSpace = false;
      }
    }
  }

  std::istringstream newStream(cleanString.c_str());
  double xx(q.x()), yy(q.y()), zz(q.z());
  newStream >> xx >> yy >> zz;
  q.setX(xx);
  q.setY(yy);
  q.setZ(zz);

  return s;
}
HepPoint operator + (const HepPoint & q, const Hep3Vector & p) {
  HepPoint r(q.x()+p.x(),q.y()+p.y(),q.z()+p.z());
  return r;
}

HepPoint operator + (const Hep3Vector & q, const HepPoint & p) {
  HepPoint r(q.x()+p.x(),q.y()+p.y(),q.z()+p.z());
  return r;
}

Hep3Vector operator - (const HepPoint & q, const HepPoint & p){
  Hep3Vector r(q.x()-p.x(),q.y()-p.y(),q.z()-p.z());
  return r;
}

HepPoint operator - (const HepPoint & q, const Hep3Vector & p) {
  HepPoint r(q.x()-p.x(),q.y()-p.y(),q.z()-p.z());
  return r;
}

HepPoint operator * (const HepRotation & r, const HepPoint & p) {
  return HepPoint(r.xx()*p.x() + r.xy()*p.y() + r.xz()*p.z(),
                  r.yx()*p.x() + r.yy()*p.y() + r.yz()*p.z(),
                  r.zx()*p.x() + r.zy()*p.y() + r.zz()*p.z());
}


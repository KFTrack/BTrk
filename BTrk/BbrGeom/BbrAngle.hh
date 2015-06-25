#ifndef BBRANGLE_HH
#define BBRANGLE_HH
////////////////////////////////
//
//  Angle class 
//
//  Luca Lista
//  06 Jan 1996
//
////////////////////////////////
//
// BaBar angles are in radians, and degress should inly be used
// when absolutely necessary. Automatic conversions to and from
// the radians form are provided, but you have to manually
// go to and from degrees
//
// By convention, angles are represented as (-pi, pi]
//
// 14 Oct 1996  Ed Frank    Replace all instances of *this with _phi.
//                          This avoids formation of a termporary.  I don't
//                          think it affects code longevity.
//  2 Nov 1996  M. Haire    changed input classes to references (HP-gcc2.7 bug)   
//
//  5 May 1997  L.Lista     * BaBar method name convention applied 
//                          * fixed rad to degrees conversion
//                          * added the metod arc(double)
// 22 Sep 1997  L.Lista     Removed pi, twoPi, toDegree
//                          now taken from Constants
//
#include "BTrk/BaBar/Constants.hh"
#include <string>

class BbrAngle
{
public:
  inline BbrAngle();
  inline BbrAngle(const double);
  inline ~BbrAngle();

  inline operator double() const  { return _phi;};   // automatic conversion to double

  inline double rad() const;  
  inline double deg() const;  
  // convention : returns value in [-180, 180]

  inline double posRad() const;
  inline double posDeg() const;
  // returns 0 to 360.  This is not the BaBar convention, and should not 
  // be used for handing off values to anyone else's code

  inline void setRad(const double);
  inline void setDeg(const double);

  inline double arc(double radius) const;

  inline void setSector(int n, int max);
  inline void setSector(int n, int max, BbrAngle phi_0);

  inline BbrAngle operator + (const BbrAngle&) const;
  inline BbrAngle operator + (const double) const;  //assumes double in radians
  inline BbrAngle operator - (const BbrAngle&) const;
  inline BbrAngle operator - (const double) const;  //assumes double in radians
  inline BbrAngle operator * (const double) const;
  inline BbrAngle operator / (const double) const;
  inline friend BbrAngle operator * (const double, const BbrAngle&);

  inline void operator = (const BbrAngle);
  inline void operator += (BbrAngle);
  inline void operator -= (BbrAngle);
  inline void operator += (double);  //assumes double in radians
  inline void operator -= (double);  //assumes double in radians
  inline void operator *= (double);
  inline void operator /= (double);

  // note : > and < should have no well defined meaning ?

  inline int sector(int max);
  inline int sector(int max, BbrAngle phi_0); 
   // convention : returns values [1..max]
  std::string degString() const;
  inline friend double sin(const BbrAngle);
  inline friend double cos(const BbrAngle);
  inline friend double tan(const BbrAngle);

  // class static constants defined in .cc file
  // these are generally available as Constants::pi, Constants::twoPi, etc,
  // and once the BTrk/BaBar/Constants class is in a release they should be
  /// used instead.

  static const double pi;
  static const double twoPi;

  // old names, forwarded for migration BobJ May 21 97
  inline double Rad() const { return rad(); }
  inline double Deg() const { return deg(); }
  inline std::string DegString() const { return degString(); }
  inline int Sector(int max) { return sector(max); }
  inline int Sector(int max, BbrAngle phi_0) { return sector(max, phi_0); }

protected:
  double _phi;

  inline static double normalize(double);

  static const double toDegrees;
  static  const std::string degChar, deg1Char, deg2Char;

};

//
// Methods for BbrAngle
//

inline double BbrAngle::normalize(double angle) {
  if (angle < - Constants::pi) {
    angle += Constants::twoPi;
    if (angle < - Constants::pi) angle = fmod(angle+ Constants::pi, Constants::twoPi) + Constants::pi; 
  }
  else if (angle > Constants::pi) {
    angle -= Constants::twoPi;
    if (angle > Constants::pi) angle = fmod(angle+Constants::pi, Constants::twoPi) - Constants::pi;
  }
  return angle;
}

inline BbrAngle::BbrAngle() : _phi(0)
{ }

inline BbrAngle::BbrAngle(const double phi) : _phi(normalize(phi))
{}

inline BbrAngle::~BbrAngle() {}

inline double BbrAngle::rad() const
{ return _phi; }

inline double BbrAngle::deg() const
{ return _phi *  Constants::radToDegrees; }

inline double BbrAngle::posRad() const
{ 
  if (_phi >= 0.0) return _phi; 
  else return _phi + Constants::twoPi;
}

inline double BbrAngle::posDeg() const
{ return posRad() *  Constants::radToDegrees; }

inline void BbrAngle::setRad(const double phi)
{ _phi = normalize(phi); }

inline void BbrAngle::setDeg(const double phi)
{ setRad(phi / Constants::radToDegrees); }

inline double BbrAngle::arc(double radius) const
{ return radius * rad(); }

inline void BbrAngle::setSector(int n, int max)
{ setRad((n + 0.5) * Constants::twoPi / max); }

inline void BbrAngle::setSector(int n, int max, BbrAngle phi_0)
{ setRad((n + 0.5) * Constants::twoPi / max + phi_0._phi); }

inline BbrAngle BbrAngle::operator + (const BbrAngle &a) const
{ return BbrAngle(_phi + a._phi); }

inline BbrAngle BbrAngle::operator + (const double a) const
{ return BbrAngle(_phi + a); }

inline BbrAngle BbrAngle::operator - (const BbrAngle &a) const
{ return BbrAngle(_phi - a._phi); }

inline BbrAngle BbrAngle::operator - (const double a) const
{ return BbrAngle(_phi - a); }

inline BbrAngle BbrAngle::operator * (const double x) const
{ return BbrAngle(_phi * x); }

inline BbrAngle BbrAngle::operator / (const double x) const
{ return BbrAngle(_phi / x); }

inline BbrAngle operator * (const double x, const BbrAngle &a)
{ return BbrAngle(a * x); }

inline void BbrAngle::operator = (const BbrAngle a)
{ _phi = normalize(a._phi); 
}

inline void BbrAngle::operator += (BbrAngle a)
{ _phi = normalize(_phi + a._phi );
}

inline void BbrAngle::operator += (double a)
{ _phi = normalize(_phi + a); 
}

inline void BbrAngle::operator -= (BbrAngle a)
{ _phi = normalize(_phi - a._phi ); 
}

inline void BbrAngle::operator -= (double a)
{ _phi = normalize(_phi - a); 
}

inline void BbrAngle::operator *= (double x)
{ _phi = normalize(_phi*x); 
}

inline void BbrAngle::operator /= (double x)
{ _phi = normalize(_phi/x); 
}

inline int BbrAngle::sector(int max)
{ 
  double phi = _phi;
  if (phi < 0) phi += Constants::twoPi;
  int tmp=( int (phi / (Constants::twoPi / max) ) + 1 );
  //tmp can be bigger than max if _phi is -0.0 ryd, 020306
  if (tmp>max) return max;
  return tmp;
}

inline int BbrAngle::sector(int max, BbrAngle phi_0)
{ 
  BbrAngle t( _phi - phi_0._phi);
  return t.sector(max);
}

inline double sin(const BbrAngle a)
{ return sin(a._phi); }

inline double cos(const BbrAngle a)
{ return cos(a._phi); }

inline double tan(const BbrAngle a)
{ return tan(a._phi); }

#endif



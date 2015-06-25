//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: Consistency.hh 458 2010-01-15 11:37:35Z stroili $
//
// Description:
//	The Consistency class is a top level parent for various
//      statistical classes
//
//      This class also defines an enum ConsistentStatus for representing
//      various error conditions
//
//      Return an object of this class from a chisquared calculation, for
//      example, so that probability and chisquared can be treated
//      directly.  Future subclasses might implement other statistics
//      (Poisson, etc), and this class could use operations to combine, 
//      print, etc. (Or other statistics could be added to this class, 
//      allowing it to convert between different measures; the uses
//      will dictate the structure)
//
//      Note: implementation of this class is minimal at present; please
//      read the comments!
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Bob Jacobsen 
//
// History :
//      April 98 - Gautier Hamel de Monchenault
//                 o  introduce genealogy methods
//                 o  add significanceLevel() accessor which returns _value
//                 o  comment out the automatic translation to double
//                 o  add 'const' to operators == and <
//
//      Dec 2007, Alexandre Telnov (Princeton)
//            Added member _logLikelihood and accessor logLikelihood() 
//
// Copyright Information:
//	Copyright (C) 1995, 1996
//
//------------------------------------------------------------------------
#ifndef CONSISTENCY_HH
#define CONSISTENCY_HH

#include "BTrk/BaBar/BaBar.hh"

#include <iosfwd>
class ConsistencySet;

//-----------------
// BaBar Headers --
//-----------------

class Consistency {

  public:

  // default constructor; sets internal state to noMeasure and value to zero
  Consistency();

  // default value for likelihood is for backward compatibility
  Consistency( double consistency, double likelihood=0.);

  // copy constructor
  Consistency(const Consistency& rhs);

  virtual ~Consistency() {}

  // assignment
  Consistency& operator= (const Consistency& rhs);

  // equality -- this will never work right since one has to compare 
  // two doubles. But RW sometimes requires this. 
  bool operator==(const Consistency& rhs) const;
  // less operator is better in that regard, but it's not unique
  bool operator < (const Consistency& rhs) const; 
  // greater than operator 
  bool operator > (const Consistency& rhs) const; 

  // members for returning 
  // the statistical significance level
  // and the likelihood 
  // of the represented observation; 
  // subclasses are responsible
  // for keeping the data member up to date
  double significanceLevel() const { return _value; }
  double likelihood() const { return  _likelihood;}
  double logLikelihood() const { return  _logLikelihood;}

  // Old function kept in the interface for back compatibility.
  // Avoid to use it - use significanceLevel() instead.
  // Will disappear eventually
  double consistency() const { return _value; }

  // represent the various ways that a statistical measure could be bad
  enum ConsistentStatus { OK=0, noMeasure, underFlow, unPhysical };

  // OK          just that; usable statistical data
  // noMeasure   no statistical statement possible; DOF == 0, for example.
  //             may or may not be consistent
  // underFlow   probability too small to calculate with; may or may not want
  //             to treat as zero, but value may have been treated by machine
  //             underflow handling
  // unPhysical  because of a "can't happen" condition, probability is
  //             identically equal to zero;

  // return the status
  ConsistentStatus status() const {return _stat;}
  void setStatus(ConsistentStatus s ) { _stat = s; }

  // represent whether the measured value was "left" or "right" of
  // the expected.  This is generally only useful if there is more
  // than one provider of consistencies in a given situation (e.g.
  // multiple systems doing PID), and they have a convention for 
  // what these mean
  enum ConsistentSign { left=-1, unknown=0, right=1 };

  // return the sign
  ConsistentSign sign() const {return _sign; }
  void setSign(ConsistentSign s ) { _sign = s; }

  // use as a double returns the consistency value; this member is
  // intended as just a convenience, and may be removed if it makes
  // trouble.

  // genealogy 
  virtual const ConsistencySet* genealogy() const;

  // print
  virtual void print(std::ostream& ) const;

// log(like) if argument is positive, -999 if 0, -9999 if negative 
// (Yes, some BaBar code sets likelihood to a negative value!)
  static double likeLog(double like);

protected:

  ConsistentStatus _stat;
  double _value;        // value of the consistency
  double _likelihood;   // value of the likelihood 
  double _logLikelihood;   // value of the likelihood's natural logarithm
  ConsistentSign _sign;
  

  //------------------
  // Static methods --
  //------------------

public:
  static const Consistency& badMeasurement();

};

#endif





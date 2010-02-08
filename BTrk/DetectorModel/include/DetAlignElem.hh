// ------------------------------------------------------------------------------
// File and Version Information:
//      $Id: DetAlignElem.hh,v 1.20 2007/09/20 23:01:42 gapon Exp $
//
//  Description:  Define the alignment parameters for a detector element.  These
//  consist of a translation vector plus a set of alignment angles (sequential
//  rotations about x,y,z axes).  Note that the coordinate system is NOT FIXED
//  to the BaBar coordinate system, but could be local to a part of BaBar.
//  Access in the form of a HepTransform is provided for the DetElem
//  functions.
//
// Copyright Information:
//	Copyright (C) 1996	Lawrence Berkeley Laboratory
//
//  Authors: Dave Brown, 12/5/96
//------------------------------------------------------------------------------
#ifndef DetAlignElem_HH
#define DetAlignElem_HH
#include "BaBar/BaBarODMGTypes.h" // Objy-safe types
#include <string>
#include <iostream>
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Vector/ThreeVector.h"
class HepTransformation;
#include "CLHEP/Vector/Rotation.h"
class AlignAngles;
#define NAMELENGTH 20
//
//  Define the class
//
class DetAlignElem{
public:
  enum {alndx,alndy,alndz,alnrx,alnry,alnrz}; // define the components of the vector
//  construct from a displacement and angle set (error set to 0)
  DetAlignElem(const Hep3Vector&,const AlignAngles&,const char*,int);
//  construct from a HepTransformation (error set to 0)
  DetAlignElem(const HepTransformation &, const char*, int);
//  construct from a parameter vector and matrix
  DetAlignElem(const HepVector&,const HepSymMatrix&,const char*,int);
  DetAlignElem(); // dummy constructor
  DetAlignElem(const char*,int); // zero parameter constructor for a particular element
  DetAlignElem(const DetAlignElem&); // copy constructor
  virtual ~DetAlignElem(){;}
//
//  operators
//
  DetAlignElem& operator = (const DetAlignElem&);
  DetAlignElem& operator *= (const DetAlignElem&);
  bool operator == (const DetAlignElem& other) const;

//
//  Access
//
  void print(std::ostream& out) const;
  Hep3Vector displacement() const;
  AlignAngles angles() const;
  HepTransformation transform() const;
  HepVector parameterVector() const ; // return parameter vector as a HepVector
  HepSymMatrix parameterCovariance() const ; // same for covariance
  HepTransformation inverseTransform() const;
  int elementNumber() const {return _ielem; } // element number
  std::string elementName() const;
  std::string elName() const;
// modifier
  void setCovariance(const HepSymMatrix& covar);
private:
//
//  The data members are declared as fixxed-length double arrays, so that they
//  can be introduced into objectivity without having to write a ddl.  The
//  intention is that specific subdetector alignment will define database
//  classes which use these to describe the alignment value (see SVTGlobalAlign
//  and SVTWaferAlign in the SvtTrack package).
//
  d_Long _ielem; // identify the element
  d_Char _ename[NAMELENGTH]; // name, must be fixed length for Objy storage
  d_Double _parvec[6]; // vector of parameters
  d_Double _parcov[21]; // upper diagonal covariance matrix

  float angle( float sina, float cosa );
  
  AlignAngles alignAngles( const HepRotation & r );

  friend class  DchPlateAligCmpr;
  friend class  DchGlobalAlignCmpr;
  friend class  EmcGlobalAlignCmpr;
  friend class  EmcGlobalAlignCmpr_001;
  friend class SvtWaferAlignCmpr;
  friend class SvtGlobalAlignCmpr;
  friend class  DrcAlignElementCmpr; 
  friend class DetAlignElemCdbR;
  friend class TestDetectorModelCdb;
  friend bool testCdb(const DetAlignElem*, const DetAlignElem*);

};
#endif

// ------------------------------------------------------------------------------
// File and Version Information:
//      $Id: DetElemSurfRange.hh,v 1.3 1999/11/30 21:42:35 brownd Exp $
//
//  Description:
//     Trivial class to keep track of the out-of-surface range
//     subtended by an element in a surface set.
//
// Copyright Information:
//	Copyright (C) 1999	Lawrence Berkeley Laboratory
//
//  Authors: Dave Brown, 4/19/99
// ------------------------------------------------------------------------------
#ifndef DETELEMSURFRANGE_HH
#define DETELEMSURFRANGE_HH

class DetElem;

class DetElemSurfRange {
public:
  DetElemSurfRange() : _elem(0), _mindist(0.0), _maxdist(-1.0)
    {;}
  DetElemSurfRange(const DetElem* elem,double mindist,double maxdist) :
    _elem(elem), _mindist(mindist), _maxdist(maxdist)
    {;}
  DetElemSurfRange(const DetElemSurfRange& other) :
    _elem(other._elem), _mindist(other._mindist), _maxdist(other._maxdist)
    {;}
  ~DetElemSurfRange(){;}
  DetElemSurfRange& operator = (const DetElemSurfRange& other) {
    if (&other != this){
      _elem = other._elem;
      _mindist = other._mindist;
      _maxdist = other._maxdist;
    }
    return *this;
  }
// interface
  const DetElem* element() const { return _elem; }
  const double& minDist() const { return _mindist; }
  const double& maxDist() const { return _maxdist; }
// needed for RW sorting
  bool operator == (const DetElemSurfRange& other) const {
    return other._elem == _elem; }
  bool operator < ( const DetElemSurfRange& other) const {
    return _elem < other._elem; }
private:
  const DetElem* _elem;
  double _mindist;
  double _maxdist;
};

#endif

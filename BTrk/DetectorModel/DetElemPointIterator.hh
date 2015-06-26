// ------------------------------------------------------------------------------
// File and Version Information:
//      $Id: DetElemPointIterator.hh,v 1.2 2002/12/30 15:44:28 dbrown Exp $
//
//  Description:
//
// Copyright Information:
//	Copyright (C) 1997	Lawrence Berkeley Laboratory
//
//  Authors: Gautier Hamel de Monchenault and Stephen Gowdy
// ------------------------------------------------------------------------------
#ifndef DETELEMPOINTITERATOR_HH
#define DETELEMPOINTITERATOR_HH

class DetElem;
class HepPoint;
class TypeCoord;
#include <vector>

enum Action { Continue=0, CloseLine, CloseShape };

class DetElemPointIterator {

public :

// construct the iterator 
  DetElemPointIterator( const DetElem& );

// destructor
  virtual ~DetElemPointIterator();

// iterator function - 
//    returns the next action to perform and passes the next point
  Action next( HepPoint& );

private:

  const DetElem* _elem;
  const std::vector<TypeCoord*>* _vect;
  unsigned _current;

};

#endif

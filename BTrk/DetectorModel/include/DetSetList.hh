//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: DetSetList.hh,v 1.1 2004/12/14 07:10:18 bartoldu Exp $
//
// Description: 
//      Container for detector sets.
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Rainer Bartoldus
//
// Copyright Information:
//      Copyright (C) 2004      Stanford Linear Accelerator Center
//
//--------------------------------------------------------------------------
#ifndef DETSETLIST_HH
#define DETSETLIST_HH

//---------------
// C++ Headers --
//---------------
//#include <list>
#include <vector>

//------------------------------------
// Collaborating Class Declarations --
//------------------------------------
class DetSet;

//typedef std::list<DetSet*> DetSetList;
typedef std::vector<DetSet*> DetSetList;

#endif // DETSETLIST_HH

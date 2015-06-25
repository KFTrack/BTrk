// -----------------------------------------------------------------------
// File and Version Information:
// $Id: TrkDirection.cc,v 1.2 2004/08/06 06:31:39 bartoldu Exp $
//
// This defines an output operator to write the strings "trkIn" and "trkOut"
// for the trkDirection enumeration values.
//
// Michael Kelsey <kelsey@slac.stanford.edu>
//
// Copyright Information:
//	Copyright (C) 1999  Princeton University
//
//------------------------------------------------------------------------

#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/TrkBase/TrkDirection.hh"
#include <iostream>
using std::ostream;

ostream& operator<<(ostream& os, const trkDirection& dir) {
  if (dir==trkIn) return os << "trkIn";
  if (dir==trkOut) return os << "trkOut";
  return os << (int)dir;
}

//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkKalComposite.hh,v 1.3 2002/10/08 18:42:33 brownd Exp $
//
// Description:
//      class TrkKalComposite.  Transient interface base class for persistent
//      form of Kalman track list.
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Copyright Information:
//	Copyright (C) 2000	Lawrence Berkeley Laboratory
//
// Author List:
//      Dave Brown 10/31/00
//------------------------------------------------------------------------
#ifndef TRKKALCOMPOSITE_HH
#define TRKKALCOMPOSITE_HH

class TrkKalComposite {
public:
// only default constructor, as there can be NO DATA MEMBERS in this class
  TrkKalComposite();
// pure virtual destructor, as this is a pure interface class
  virtual ~TrkKalComposite() = 0;
// a TrkKalComposite knows how many svt, dch hots it has
//  These are the total for all tracks combined
  virtual unsigned nSvt() const = 0;
  virtual unsigned nDch() const = 0;
// a TrkKalComposite knows how many fit results it has
  virtual unsigned nFit() const = 0;
private:
};

#endif

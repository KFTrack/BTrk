//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkOutAvgTimeCalculator.hh,v 1.1 2001/07/23 22:25:25 raven Exp $
//
// Description:
//   class TrkOutAvgTimeCalculater.  A simple implementation of TrkTimeCalculator
//   that takes the weighted average of all selected hot's times assuming a
//   Gaussian error, but rejecting outliers based on their pull.
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Copyright Information:
//	Copyright (C) 2001	UC, San Diego
//
// Author List:
//      Gerhard Raven 7/19/01
//------------------------------------------------------------------------

#ifndef TRKOUTAVGTIMECALCULATOR_HH
#define TRKOUTAVGTIMECALCULATOR_HH

#include "TrkBase/TrkTimeCalculator.hh"
#include <assert.h>
#include <math.h>

class TrkOutAvgTimeCalculator : public TrkTimeCalculator {
public:
  TrkOutAvgTimeCalculator(const TrkHotSelector& selector,double maxpull);
  virtual ~TrkOutAvgTimeCalculator();
// The Function
  virtual bool trackTime(const TrkRecoTrk& trk,
                         double& time, double& timeerr,
                         int& nHotsUsed) const;
private:
// disallow
  TrkOutAvgTimeCalculator(const TrkOutAvgTimeCalculator&);
  TrkOutAvgTimeCalculator& operator = (const TrkOutAvgTimeCalculator&);

  double _maxpull;

  class ws {
      public:
          ws(double x=0,double w=0):_w(w),_wx(w*x),_n(w>0?1:0) { assert(!(w<0));}
          ws& operator+=(const ws& x) { _wx += x._wx; _w+=x._w; _n+=x._n; return *this;}
          ws& operator-=(const ws& x) { _wx -= x._wx; _w-=x._w; _n-=x._n; return *this;}
          bool operator==(const ws& x) const { return _w==x._w && _wx==x._wx && _n==x._n;}
          bool isPhysical() const { return _w>0;}
          double mean() const {return _wx/_w;}
          double sigma() const { return double(1)/sqrt(_w);}
          double pull(const ws& x) const { return (mean()-x.mean())/sqrt(sigma2()+x.sigma2()); }
          unsigned n() const { return _n;}
      private:
          double sigma2() const { return double(1)/_w;}
          double _w,_wx;
          unsigned _n;
  };
};
#endif

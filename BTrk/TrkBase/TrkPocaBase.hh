//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkPocaBase.hh,v 1.11 2006/03/25 15:15:56 brownd Exp $
//
// Description:
//   Base class for various Poca classes; holds infrastructure, and one 
//   common algorithm.  Ctor and dtor protected, to prevent instantiation.
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author(s): Steve Schaffner
//------------------------------------------------------------------------
#ifndef TRKPOCABASE_HH
#define TRKPOCABASE_HH

#include "BTrk/TrkBase/TrkErrCode.hh"
class Trajectory;
class HepPoint;

// struct for passing around associated traj information
struct TrkPocaTraj{
  TrkPocaTraj(const Trajectory& traj, double flt, bool rflt) :
    _traj(traj),_flt(flt),_rflt(rflt) {}
  const Trajectory& _traj;
  double _flt;
  bool _rflt;
};

// Class interface
class TrkPocaBase {

protected:
  TrkPocaBase(double flt1, double flt2, double precision);
  TrkPocaBase(double flt1, double precision);
  TrkPocaBase();
  TrkPocaBase(const TrkPocaBase& other);
  TrkPocaBase& operator = (const TrkPocaBase& other);

public:
  inline const TrkErrCode& status() const;     // did the calculation succeed?
  inline double flt1() const;                  // path length on traj 1 @ poca
  inline double flt2() const;
  inline double precision();                   // In case anyone wants to know
  virtual ~TrkPocaBase();
  virtual double doca() const = 0;
  virtual TrkPocaBase* clone() const = 0;  // returns ownership

protected:
  double _precision;
  double _flt1;
  double _flt2;
  TrkErrCode _status;

  void minimize(TrkPocaTraj& traj1, TrkPocaTraj& traj2);
  void minimize(TrkPocaTraj& traj1, const HepPoint& pt);

  void stepTowardPoca(TrkPocaTraj& traj1,TrkPocaTraj& traj2);
  void stepToPointPoca(TrkPocaTraj& traj, const HepPoint& pt);

  static void setFlt(double flt,TrkPocaTraj& traj);

  static double _maxDist;                   // step > maxDist => parallel
  static int _maxTry;
  static double _extrapToler;            // error allowed in picking step size
};

// Inlined functions
double
TrkPocaBase::precision()                 {return _precision;}

const TrkErrCode&
TrkPocaBase::status() const              {return _status;}

double
TrkPocaBase::flt1() const                {return _flt1;}

double
TrkPocaBase::flt2() const                {return _flt2;}

#endif

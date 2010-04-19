// ------------------------------------------------------------------------------
// File and Version Information:
//      $Id: KalPairSite.hh,v 1.21 2004/08/06 06:12:48 bartoldu Exp $
//
//    class KalPairSite.  Special site to define the constraint connecting two
//    tracks at a point.  This is implemented as a vertex constraint (they come
//    from the same point in space) and a momentum constraint (from pre-knowledge
//    of the total momentum before decay).  The error matrix can express
//    uncertainty about the momentum constraint.
//
// Copyright Information:
//	Copyright (C) 1998	Lawrence Berkeley Laboratory
//
//  Authors: Dave Brown, 3/21/98
//------------------------------------------------------------------------------
//
#ifndef KALPAIRSITE_HH
#define KALPAIRSITE_HH

#include "KalmanTrack/KalConstraint.hh"
#include "KalmanTrack/KalPairConduit.hh"

class KalPairConduit;

class KalPairSite : public KalConstraint {
public:
// contruct from momentum constraint, vertex constraint is implicit
  KalPairSite(const TrkDifPieceTraj* ptraj, const TrkParams& constraint,
	      bool* constrainParams, double fltlen, KalPairConduit* conduit);
  virtual ~KalPairSite();
// KalSite functions
  KalPairSite* clone(const KalRep*) const;
  bool update(const TrkDifPieceTraj*,double);
  bool process(const KalSite* prevsite, trkDirection idir); 
  void printAll(std::ostream& os=std::cout) const;

private:
  // Preempt
  KalPairSite(const KalPairSite&);

  KalPairConduit* _conduit;

  void changeConstraint(const KalParams& constraint);

};

#endif

//--------------------------------------------------------------------------
// File and Version Information:
//
// Description:  TrkVisitor is an implementation of the Visitor pattern
//      for use in calculating various things for different types of
//      trajectories. (See the book _Design Patterns_ or one of the
//      authors for a definition of the Visitor pattern.)  It is the 
//      abstract base class for visitors such as MomVisitor (the
//      momentum visitor) and others.
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author(s): Justin Albert, Steve Schaffner
//
//------------------------------------------------------------------------

#ifndef TRKVISITOR_HH
#define TRKVISITOR_HH

class HelixTraj;
class TrkCircleTraj;
class NeutTraj;
class TrkDifLineTraj;

// Class interface //
class TrkVisitor {

public:

  TrkVisitor();
  virtual ~TrkVisitor();

  //********************************
  //The visitor functions:
  //********************************

  virtual void trkVisitHelixTraj(const HelixTraj*) = 0;
  virtual void trkVisitCircleTraj(const TrkCircleTraj*) = 0;
  virtual void trkVisitNeutTraj(const NeutTraj*) = 0;
  virtual void trkVisitLineTraj(const TrkDifLineTraj*) = 0;


};

#endif

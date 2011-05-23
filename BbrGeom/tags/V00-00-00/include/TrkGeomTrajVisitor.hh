//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkGeomTrajVisitor.hh 491 2010-01-13 16:59:16Z stroili $
//
// Description:
//     Part of Visitor pattern for TrkGeomTraj class hierarchy.  This is 
//     the base class for any Visitor that will visit a TrkGeomTraj.  As 
//     new GeomTraj's are introduced, they must have an appropriate 
//     visitXXXX function added here and in all derived Visitors.
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author(s): Steve Schaffner
//
//------------------------------------------------------------------------

#ifndef TRKGEOMTRAJVISITOR_HH
#define TRKGEOMTRAJVISITOR_HH

class TrkLineTraj;
class TrkPieceLineTraj;
class TrkParabolaTraj;
class BbrHelixTraj;

// Class interface //
class TrkGeomTrajVisitor {

public:
  TrkGeomTrajVisitor();
  virtual ~TrkGeomTrajVisitor();

  virtual void visitLine(const TrkLineTraj*) = 0;
  //  virtual void visitParabola(const TrkParabolaTraj*) = 0;
  virtual void visitPieceLine(const TrkPieceLineTraj*) = 0;
  virtual void visitHelix(const BbrHelixTraj*) = 0;
  
private:	

  // Preempt
  TrkGeomTrajVisitor&   operator= (const TrkGeomTrajVisitor&);
  TrkGeomTrajVisitor(const TrkGeomTrajVisitor &);
};

#endif








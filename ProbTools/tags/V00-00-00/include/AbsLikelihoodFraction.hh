//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: AbsLikelihoodFraction.hh 458 2010-01-15 11:37:35Z stroili $
//
// Description:
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	Gautier Hamel de Monchenault - CEN Saclay & Lawrence Berkeley Lab
//
// History (add to end):
//      Gautier   Apr, 1998  - creation
//
// Copyright Information:
//
//------------------------------------------------------------------------

#ifndef AbsLikelihoodFraction_hh
#define AbsLikelihoodFraction_hh

#include <sys/types.h>

class AbsLikelihoodFraction
{

public:

  // no constructor, it's a pure interface
  virtual ~AbsLikelihoodFraction();
  
  // Enum for likelihood fractions.  OK is what it says.
  // Underflow says that the hypothesis is extremely unlikely
  // Overflow says that the sum of all other hypotheses is very small.
  // IllDefined says that the hypothesis is extremely unlikely 
  // *and* the sum of all other hypotheses is extremely unlikely.
  // Extremely unlikely is 10^-100 or less.
  
  enum LikelihoodStatus{OK=0, underflow, overflow, illDefined};

  // pure virtual function declarations
  virtual size_t        nHypos()              const =0;
  virtual const char* hypoName( size_t hypo ) const =0;
  virtual double  aPrioriProba( size_t hypo ) const =0; // P(H=h)
  virtual bool         isValid( size_t hypo ) const =0;
  virtual double    likelihood( size_t hypo ) const =0; // L(h|x)=f(x|h)

  // public functions implemented in the base class
  double conditionalPdf( size_t hypo ) const; // f(x|h)=L(h|x)
  double       jointPdf( size_t hypo ) const; // f(x,h)=P(H=h)*f(x|h)
  double    marginalPdf() const;              // f(x)=Sum_h{ f(x,h) }

  double likelihoodFraction( size_t hypo) const; 
  // f(x) = L(j|x)/Sum_h{L(x|h)} for h not equal to i
  double logLikelihoodFraction( size_t hypo) const; 
  // f(x) = ln[L(j|x)] - ln[Sum_h{L(x|h)}] for h not equal to i

  double inclusiveLikelihoodFraction( size_t hypo) const; 
  // f(x) = L(j|x)/Sum_h{L(x|h)} for all h
  // This is just op() below.
  double logInclusiveLikelihoodFraction( size_t hypo) const; 
  // f(x) = ln[L(j|x)] - ln[Sum_h{L(x|h)}] for all h

  // Each of the above sets the _lastCalcStatus.  This should 
  // be checked for the validity of the result.  All of the above 
  // include safeguards for FPEs.
  LikelihoodStatus statusOfLastCalc() const {return _lastCalcStatus;}

  // operators
  virtual double operator()( size_t hypo ) const;
  virtual double operator[]( const char* hypoName ) const;

protected:
  
  // helper functions
  virtual bool getIndex_( size_t& index, const char* hypoName ) const;

  
private:
  
  AbsLikelihoodFraction *myself() const 
    {return (AbsLikelihoodFraction*)this;}

  void setStatus(LikelihoodStatus);

  LikelihoodStatus combinedStatus(LikelihoodStatus status1, 
				  LikelihoodStatus status2) const;
  // status1 = numerator status 
  // status2 = denominator status
  // OK / OK = OK
  // OK / underflow = overflow
  // underflow / OK = underflow
  // OK / overflow = underflow
  // overflow / OK = overflow
  // status1 = status2 and status1 is not OK ==> illDefined
  // status1 == illDefined || status2 == illDefined ==> illDefined
  // (shouldn't happen.  illDefined should only happen as a result 
  // of this combination!)
  LikelihoodStatus _lastCalcStatus;

};

#endif


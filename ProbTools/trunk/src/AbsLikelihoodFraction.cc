//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: AbsLikelihoodFraction.cc 458 2010-01-15 11:37:35Z stroili $
//
// Description:
//	Class AbsLikelihoodFraction
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Gautier Hamel de Monchenault, LBNL & Saclay
//
// History (add to end):
//
// Copyright Information:
//
//------------------------------------------------------------------------

//----------------
// BaBar Header --
//----------------
#include "BaBar/BaBar.hh"
#include <iostream>
//-----------------------
// This Class's Header --
//-----------------------
#include "ProbTools/AbsLikelihoodFraction.hh"

//---------------
// C++ Headers --
//---------------
#include <assert.h>
#include <string>
#include <math.h>
//-------------------------------
// Collaborating Class Headers --
//-------------------------------

//-----------------------------------------------------------------------
// Local Macros, Typedefs, Structures, Unions and Forward Declarations --
//-----------------------------------------------------------------------

static const char rscid[] = "$Id: AbsLikelihoodFraction.cc 458 2010-01-15 11:37:35Z stroili $";

//		-----------------------------------------------
// 		-- Static Data & Function Member Definitions --
//		-----------------------------------------------

//		----------------------------------------
// 		-- Public Function Member Definitions --
//		----------------------------------------


AbsLikelihoodFraction::~AbsLikelihoodFraction()
{}

double 
AbsLikelihoodFraction::operator()( size_t hypo ) const
{
  double sum = marginalPdf();  
  LikelihoodStatus status = statusOfLastCalc();
  if( sum>0. ){
    double result = jointPdf( hypo )/sum;
    LikelihoodStatus status2 = statusOfLastCalc();
    myself()->setStatus(combinedStatus(status2, status));
    return result;
  }
  return 0.;
}

double 
AbsLikelihoodFraction::operator[]( const char* hypoName ) const
{
  size_t index;
  assert( getIndex_(index,hypoName) );
  return (*this)(index);
}

double 
AbsLikelihoodFraction::conditionalPdf( size_t hypo ) const
{
  return likelihood( hypo );
}


double 
AbsLikelihoodFraction::jointPdf( size_t hypo ) const
{
  double result = aPrioriProba( hypo )*conditionalPdf( hypo );
  if (result > 1.0E-100 & result < 1.0E100) {
    myself()->setStatus(AbsLikelihoodFraction::OK);
    return result;
  }
  if (result < 1.0E-100) {
    myself()->setStatus(AbsLikelihoodFraction::underflow);
    result = 1.0e-100;
  } else {
    myself()->setStatus(AbsLikelihoodFraction::overflow);
    result = 1.0e100;
  }
  return result;
}

double 
AbsLikelihoodFraction::marginalPdf( ) const
{
  double sum(0);
  for( size_t hypo=0; hypo<nHypos(); hypo++ )
    sum += jointPdf( hypo );
  if (sum > 1.0E-100 & sum < 1.0E100) {
    myself()->setStatus(AbsLikelihoodFraction::OK);
    return sum;
  }
  if (sum < 1.0E-100) {
    myself()->setStatus(AbsLikelihoodFraction::underflow);
    sum = 1.0e-100;
  } else {
    myself()->setStatus(AbsLikelihoodFraction::overflow);
    sum = 1.0e100;
  }
  return sum;
}

double 
AbsLikelihoodFraction::likelihoodFraction(size_t hypo) const{
  

  double sum = marginalPdf();  
  LikelihoodStatus status = statusOfLastCalc();
  double thisHypo = jointPdf(hypo);
  LikelihoodStatus status2 = statusOfLastCalc();
  sum -= thisHypo;  
  if( sum>0. ){
    double result = thisHypo/sum;
    myself()->setStatus(combinedStatus(status2, status));
    return result;
  }
  return 0.;
  
}

double 
AbsLikelihoodFraction::logLikelihoodFraction(size_t hypo) const{
  double result = likelihoodFraction(hypo);

  // This _really_ shouldn't happen.  
  if (result==0.0) {
    myself()->setStatus(AbsLikelihoodFraction::illDefined);
    return log(1.0e-100);
  }

  if (statusOfLastCalc() == AbsLikelihoodFraction::OK) return log(result);
  if (statusOfLastCalc() == AbsLikelihoodFraction::overflow) return log(1.0e100);
  return log(1.0e-100);
}

double 
AbsLikelihoodFraction::inclusiveLikelihoodFraction(size_t hypo) const{
  return (*this)(hypo);
}

double 
AbsLikelihoodFraction::logInclusiveLikelihoodFraction(size_t hypo) const{
  double result = inclusiveLikelihoodFraction(hypo);

  // This _really_ shouldn't happen.  
  if (result==0.0) {
    myself()->setStatus(AbsLikelihoodFraction::illDefined);
    return log(1.0e-100);
  }


  if (statusOfLastCalc() == AbsLikelihoodFraction::OK) return log(result);
  if (statusOfLastCalc() == AbsLikelihoodFraction::overflow) return log(1.0e100);
  return log(1.0e-100);
}

//		-----------------------------------
// 		-- Internal Function Definitions --
//		-----------------------------------


bool
AbsLikelihoodFraction::getIndex_( size_t& index, const char* aName ) const
{
  std::string theString(aName);
  for( size_t hypo=0; hypo<nHypos(); hypo++ )
    {
      std::string theHypoName (hypoName( hypo ) );
      if( theString == theHypoName ) 
	{
	  index=hypo;
	  return true;
	}
    }
  return false;
}


AbsLikelihoodFraction::LikelihoodStatus 
AbsLikelihoodFraction::combinedStatus(AbsLikelihoodFraction::LikelihoodStatus status1, 
				      AbsLikelihoodFraction::LikelihoodStatus status2) const{
  
  // Start with the most likely combo to save time
  
  if (status1 == status2) {

    if (status1 == AbsLikelihoodFraction::OK) 
      return AbsLikelihoodFraction::OK;
    else
      return AbsLikelihoodFraction::illDefined;

  } else if (status1 == AbsLikelihoodFraction::OK){

    if (status2 == AbsLikelihoodFraction::underflow) 
      return AbsLikelihoodFraction::overflow;

    if (status2 == AbsLikelihoodFraction::overflow) 
      return AbsLikelihoodFraction::underflow;

  } else if (status2 == AbsLikelihoodFraction::OK
	     && status1 != AbsLikelihoodFraction::illDefined){

    return status1;

  }
  return AbsLikelihoodFraction::illDefined;
}

void
AbsLikelihoodFraction::setStatus(AbsLikelihoodFraction::LikelihoodStatus status){
  _lastCalcStatus = status;

}

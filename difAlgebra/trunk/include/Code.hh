//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: Code.hh 501 2010-01-14 12:46:50Z stroili $
//
// Description:
//	Class Header for |Code|
//      Tell about success and failure
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	A. Snyder
//
// Copyright Information:
//	Copyright (C) 1996	SLAC
//
//------------------------------------------------------------------------

#ifndef Code_HH
#define Code_HH

#include <assert.h>
#include <stdlib.h>

#include <iosfwd>

class Code {

public:

  //constructors
  
  //default to success, default success code is 1
  Code(int s=1,int f=0):_fail(0),_success(0) 
  {
    if(f==0) {setSuccess(s);}
    else if(s==0) {setFail(f);}
  }

  //copy
  Code(const Code &c) 
    :_fail(c.fail()),_success(c.success())
  {}

  //access

  inline int fail()const {return _fail;}
  inline int success()const {return _success;}

  //set
  inline void setFail(int i) 
  {assert(i); _fail=i; _success=0;}
  inline void setSuccess(int i) 
  {assert(i); _success=i; _fail=0;}


private:

  //data

  int _fail;			// failure code
  int _success;			// success code

};

#endif

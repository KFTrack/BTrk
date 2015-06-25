//--------------------------------------------------------------------------
// File and Version Information:
// 	$Id: NumRecipes.cc 458 2010-01-15 11:37:35Z stroili $
//
// Description:
//
// Environment:
//	Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//	Bob Jacobsen, Ed Iskander
//
// Copyright Information:
//	Copyright (C) 1996	Lawrence Berkeley Laboratory
//
//------------------------------------------------------------------------
#include "BTrk/BaBar/BaBar.hh"

//-----------------------
// This Class's Header --
//-----------------------
#include "BTrk/ProbTools/NumRecipes.hh"

//-------------
// C Headers --
//-------------
extern "C" {
#include <assert.h>
#include <math.h> 
#include <stdlib.h>
}

//---------------
// C++ Headers --
//---------------
#include <iostream>
#include "BTrk/BaBar/ErrLog.hh"
using std::endl;

double NumRecipes::gammln(double xx)
{
   double x,y,tmp,ser;
   static double cof[6]={76.18009172947146,-86.50532032941677,
       24.01409824083091,-1.231739572450155,
       0.1208650973866179e-2,-0.5395239384953e-5};
   int j;

   y=x=xx;
   tmp=x+5.5;
   tmp -= (x+0.5)*log(tmp);
   ser=1.000000000190015;
   for (j=0;j<=5;j++) ser += cof[j]/++y;
   return -tmp+log(2.5066282746310005*ser/x);
}

double 
NumRecipes::gammp(double a, double x)
{
   double gamser,gammcf,gln;

   if (x < 0.0 || a <= 0.0) ErrMsg(error) <<" Invalid arguments in routine gammp x=" << x << " a=" << a << endl;
   if (x < (a+1.0)) {
      gser(&gamser,a,x,&gln);
      return gamser;
   } else {
      gcf(&gammcf,a,x,&gln);
      return 1.0-gammcf;
   }
}

double 
NumRecipes::gammq(double a, double x)
{
   double gamser,gammcf,gln;

   if (x < 0.0 || a <= 0.0) recipesErr(" Invalid arguments in routine GAMMQ");
   if (x < (a+1.0)) {
      gser(&gamser,a,x,&gln);
      return 1.0-gamser;
   } else {
      gcf(&gammcf,a,x,&gln);
      return gammcf;
   }
}

void NumRecipes::gcf(double* gammcf, double a, double x, double* gln)
   {
     int n;
     double gold=0.0,g,fac=1.0,b1=1.0;
     double b0=0.0,anf,ana,an,a1,a0=1.0;

     *gln=gammln(a);
     a1=x;
     for (n=1;n<=NUMREC_ITMAX;n++) {
       an=(double) n;
       ana=an-a;
       a0=(a1+a0*ana)*fac;
       b0=(b1+b0*ana)*fac;
       anf=an*fac;
       a1=x*a0+anf*a1;
       b1=x*b0+anf*b1;
       if (a1) {
         fac=1.0/a1;
         g=b1*fac;
         if (fabs((g-gold)/g) < NUMREC_EPS) {
            *gammcf=exp(-x+a*log(x)-(*gln))*g;
            return;
         }
         gold=g;
       }
     }
     recipesErr(" a too large, NUMREC_ITMAX too small in routine GCF");
   }

void NumRecipes::gser(double* gamser, double a, double x, double* gln)
{
   int n;
   double sum,del,ap;

   *gln=gammln(a);
   if (x <= 0.0) {
      if (x < 0.0) recipesErr(" x less than 0 in routine GSER");
      *gamser=0.0;
      return;
   } else {
      ap=a;
      del=sum=1.0/a;
      for (n=1;n<=NUMREC_ITMAX;n++) {
         ap += 1.0;
         del *= x/ap;
         sum += del;
         if (fabs(del) < fabs(sum)*NUMREC_EPS) {
            *gamser=sum*exp(-x+a*log(x)-(*gln));
            return;
         }
      }
      recipesErr(" a too large, NUMREC_ITMAX too small in routine GSER");
      return;
   }
}

void NumRecipes::recipesErr(const char* c)
{
   ErrMsg(fatal) << " Numerical Recipes run-time error...\n" << c 
                 << "\n ...now exiting to system..." << endmsg;
}


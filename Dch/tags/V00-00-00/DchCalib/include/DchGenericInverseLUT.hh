//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: DchGenericInverseLUT.hh 88 2010-01-14 12:32:57Z stroili $
//
// Description:
//      Class DchGenericInverseLut:
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Gerhard Raven           8/05/98
//
// Copyright Information:
//      Copyright (C) 1998      University of California, San Diego
//
//------------------------------------------------------------------------
#ifndef DCHGENERICINVERSELUT_HH
#define DCHGENERICINVERSELUT_HH

#include <list>
#include <functional>
#include <iostream>
class DchCalibFun;
#include "CLHEP/Matrix/Vector.h"

class DchGenericInverseLUT
{
public:
        DchGenericInverseLUT(const DchCalibFun& fun,HepVector &x, 
                             double mintime, double maxtime, double tolerance);
        double value(double x) const { return interpol(x,findInter(x)); }
        double fieldValue(unsigned i) const;
        unsigned entries() const { return _tbl.size();}
private:
        class less_d;
        class greater_d;
        friend class less_d;
        friend class greater_d;
        class doublet {
        public:
                doublet(): t(0),d(0) {;}
                doublet(double xt,double xd): t(xt),d(xd) {;}
                double t,d;
                bool operator<(const doublet& x) const { return t<x.t;}
                bool operator>(const doublet& x) const { return !operator<(x);}
                bool operator==(const doublet& x) const { return t==x.t && d==x.d; }
                bool operator!=(const doublet& x) const { return !operator==(x);}
                std::ostream& print(std::ostream& o) const { return o << "(" << t << "," << d << ")"; }
        };
        std::list<doublet> _tbl;
        std::list<doublet>::iterator findInter(double) const;
        double interpol(double d,std::list<doublet>::iterator l) const
        {
                const doublet &a=*l; const doublet &b=*(++l);
                double delta=b.d-a.d;
                return delta==0?b.t:a.t+(d-a.d)*(b.t-a.t)/delta;
        }
        class less_d : public std::unary_function<doublet,bool>
        {
             public:
                less_d(double m) : _m(m) {};
                bool operator()(const doublet& x) const { return x.d < _m; }
             private:
                double _m;
        };
        class greater_d : public std::unary_function<doublet,bool>
        {
             public:
                greater_d(double m) : _m(m) {};
                bool operator()(const doublet& x) const { return _m < x.d ; }
             private:
                double _m;
        };
};
#endif

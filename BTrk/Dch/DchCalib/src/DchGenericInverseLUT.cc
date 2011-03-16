//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: DchGenericInverseLUT.cc 88 2010-01-14 12:32:57Z stroili $
//
// Description:
//      Class DchGenericInverseLUT
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
#include "BaBar/BaBar.hh"

#include "DchCalib/DchGenericInverseLUT.hh"
#include "DchCalib/DchCalibFun.hh"
#include "CLHEP/Matrix/Vector.h"
#include "ErrLogger/ErrLog.hh"
#include <assert.h>
#include <iostream>
#include <algorithm>

// std::ostream& operator<<(std::ostream &o, const DchGenericInverseLUT::doublet& x) 
// { return x.print(o); };


DchGenericInverseLUT::DchGenericInverseLUT(const DchCalibFun& fun,HepVector &x, 
                     double mintime, double maxtime, double tolerance)
{
        x[0]=mintime; double y0=fun.value(x); _tbl.push_back(doublet(x[0],y0));
        x[0]=maxtime; double y1=fun.value(x); _tbl.push_back(doublet(x[0],y1>=y0?y1:y0));
        double timestep=tolerance/2;
        bool hasLocalMaximum=false;
RESTART:
        for (x[0]=mintime;x[0]<maxtime;x[0]+=timestep) {
                double d=fun.value(x);
                std::list<doublet>::iterator left=findInter(d);
                double t=interpol(d,left);
                if ( fabs(x[0]-t)<tolerance) continue;

                mintime=left->t; // don't need to recheck left of current interval
                // bisect the interval, add node and restart...
                std::list<doublet>::iterator right=left; ++right;
                x[0]=(left->t+right->t)/2;
                d=fun.value(x);
                if (left->d < d && d < right->d)  {
                        _tbl.insert(right,doublet(x[0],d));
                        // cout << "insert loop  restart: " << x[0] << " , " << d << endl;
                } else {
                        hasLocalMaximum=true;
                        // oops: there is a maximum between _tbl[left].t and _tbl[left+1].t
                        // Find the maximum using a binary search
                        // FIXME: could be speed up by doing a golden rule search...
                        double t1=left->t; double t3=right->t; double t2=(t1+t3)/2;
                        double d1=left->d; double d3=right->d; x[0]=t2 ; double d2=fun.value(x);
                        // cout << " bracket max between " << t1 << " and " << t3 << endl; 
                        while (t3>t1+tolerance/10) {
                                // cout << "looking for max " << t1 << " " << t3 << endl;
                                // cout << "                " << d1 << " " << d3 << endl;
                                double t12=(t1+t2)/2; x[0]=t12; double d12=fun.value(x);
                                double t23=(t2+t3)/2; x[0]=t23; double d23=fun.value(x);
                                if (d12>d2 && d12>d23) { // recenter on d12
                                        t3=t2; d3=d2;
                                        t2=t12; d2=d12;
                                } else if (d2>d12 && d2>d23) { // recenter of d2
                                        t1=t12; d1=d12;
                                        t3=t23; d3=d23;
                                } else if ( d23>d&&d23>d12 ) { // recenter on d23
                                        t1=t2; d1=d2;
                                        t2=t23; d2=d23;
                                } else {
                                        ErrMsg(fatal) << " More than one maximum..." 
                                                      << endmsg;
                                }
                        }
                        double tmax=t3; double dmax=d3;
                        std::list<doublet>::iterator it =_tbl.insert(right,doublet(tmax,dmax));
                        // cout << "insert max: " << tmax << " , " << dmax << endl;

                        // eliminate all points to the right of the maximum which are smaller
                        _tbl.erase(remove_if(it,_tbl.end(),DchGenericInverseLUT::less_d(dmax)),_tbl.end());

                        // cout << "list after pruning" << endl;
                        // copy(_tbl.begin(),_tbl.end(),ostream_iterator<doublet>(cout));
                        // cout << endl;

                        // check if to the right of the maximum the function ever
                        // becomes larger. If so, continue there... otherwise add
                        // final point at same height as maximum.
                        for (x[0]=tmax;x[0]<maxtime;x[0]+=timestep) {
                                if (fun.value(x)>dmax) break;
                        }
                        _tbl.push_back(doublet(x[0],dmax));
                        mintime=x[0]+timestep;
                        // cout << "list at restart at " << mintime << endl;
                        // copy(_tbl.begin(),_tbl.end(),ostream_iterator<doublet>(cout));
                        // cout << endl;

GAPFILLER_RESTART:
                        // fixup from mintime until the time of the local max
                        for (x[0]=left->t;x[0]<tmax;x[0]+=timestep) {
                                double d=fun.value(x);
                                std::list<doublet>::iterator l=findInter(d);
                                double t=interpol(d,l);
                                if ( fabs(x[0]-t)<tolerance) continue;
                                std::list<doublet>::iterator r=l; ++r;
                                x[0]=(l->t+r->t)/2;
                                if (x[0]>tmax) continue;
                                d=fun.value(x);
                                assert(l->d < d && d <= r->d) ;
                                _tbl.insert(r,doublet(x[0],d));
                                // cout << "insert loop  minor restart: " << x[0] << " , " << d << endl;
                                goto GAPFILLER_RESTART;
                        }
                        // and now continue as usual at mintime...
                }
                goto RESTART;

        }
        // cout << "list finished" << endl;
        // for (unsigned i=0;i<_tbl.entries();++i)  cout << " ( " << _tbl[i].t << " , " << _tbl[i].d << " ) " ;
}


std::list<DchGenericInverseLUT::doublet>::iterator
DchGenericInverseLUT::findInter(double d) const
{
        std::list<doublet> &tbl =const_cast<DchGenericInverseLUT*>(this)->_tbl;
        if (_tbl.size()<2) return tbl.begin();

#if 0
        // can no longer do binary search, as the 'd' range 
        // is not neccessarily sorted anymore ;-(
        while (left+1<right) {
             unsigned pivot=left+(right-left)/2;
             ((d>_tbl[pivot].d)?left:right) = pivot;
        }
#endif
        std::list<doublet>::iterator right=find_if(tbl.begin(),tbl.end(),DchGenericInverseLUT::greater_d(d));
        if (right==tbl.end()) --right;
        std::list<doublet>::iterator left=right; --left;
        // cout << "[" << *left << " | " << d << " | " << *right << " ]" << endl;

        assert(left->d<=d||left==tbl.begin());
        assert(d<right->d||right==--tbl.end());
        assert(left->d<=right->d);
        return left;
}

double
DchGenericInverseLUT::fieldValue(unsigned i) const
{
        assert(i<2*_tbl.size());
        std::list<DchGenericInverseLUT::doublet>::const_iterator j(_tbl.begin());
        for (;i>1;i-=2) ++j;
        return i==0?j->t:j->d;
}

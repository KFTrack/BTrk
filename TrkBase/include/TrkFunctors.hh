//--------------------------------------------------------------------------
// File and Version Information:
//      $Id: TrkFunctors.hh,v 1.3 2004/09/30 16:59:41 hulsberg Exp $
//
// Description:
//      general, simple tracking predicates
//
//      All classes (or structs) here should inherit from either
//      unary_function< X, Z> or binary_function<X,Y,Z> and
//      provide Z operator()(const X&,const Y&) const.
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//      Gerhard Raven
//------------------------------------------------------------------------
#ifndef TRKFUNCTORS_HH
#define TRKFUNCTORS_HH
#include <functional>

#include "TrkBase/TrkHitOnTrk.hh"
#include "TrkBase/TrkErrCode.hh"
class TrkHitOnTrkUpdater;
class TrkRep;

namespace TrkBase { namespace Functors {



        class cloneHot : public std::unary_function<TrkHitOnTrk,TrkHitOnTrk*> {
        public:
                cloneHot(TrkRep *parentRep, const TrkDifTraj* trkTraj=0) : _r(parentRep), _t(trkTraj) {}

                TrkHitOnTrk *operator()(const TrkHitOnTrk& h) const 
                { return h.clone(_r,_t); }
        private:
                TrkRep *_r;
                const TrkDifTraj* _t;
        };

        template <class T> struct takeAddress : std::unary_function<T,T*> {
                T* operator()(T& t) { return &t; }
        };


        // Functors which are only accessible by inheriting from TrkHitOnTrkUpdater
        // FIXME: maybe these should live in the TrkHitOnTrkUpdater namespace instead???

        class updateMeasurement : public std::unary_function<TrkHitOnTrk,TrkErrCode> {
        public:
                TrkErrCode operator()(TrkHitOnTrk& h) const
                { return h.updateMeasurement(_t); }
        private:
                // Only TrkHitOnTrkUpdater can create one of these...
                friend class ::TrkHitOnTrkUpdater;
                updateMeasurement(const TrkDifTraj* traj=0) : _t(traj) {}
                const TrkDifTraj *_t;
        };

        class setActive : public std::unary_function<TrkHitOnTrk,void> {
        public:
                TrkHitOnTrk* operator()(TrkHitOnTrk& h) const
                { return h.setActive(_a); }
                TrkHitOnTrk* operator()(TrkHitOnTrk*& h) const
                { return h->setActive(_a); }
        private:
                // Only TrkHitOnTrkUpdater can create one of these...
                friend class ::TrkHitOnTrkUpdater;
                setActive(bool active) : _a(active){}
                bool _a;
        };

        class setParent : public std::unary_function<TrkHitOnTrk,void> {
        public:
                TrkHitOnTrk* operator()(TrkHitOnTrk& h) const
                { return h.setParent(_p); }
        private:
                // Only TrkHitOnTrkUpdater can create one of these...
                friend class ::TrkHitOnTrkUpdater;
                setParent(TrkRep* parent) : _p(parent){}
                TrkRep* _p;
        };
}};
#endif

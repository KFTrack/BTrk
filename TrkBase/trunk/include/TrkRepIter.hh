#ifndef TRKREPITER_HH
#define TRKREPITER_HH
#include "boost/shared_ptr.hpp"
#include <iterator>

//FIXME: only works on Linux:    class TrkHitOnTrkIter : public std::iterator_traits<TrkHitOnTrk *>
//FIXME: only works on Linux:    class TrkHitOnTrkIter : public std::random_access_iterator<const TrkHitOnTrk, ptrdiff_t>
//FIXME: only works on SunOS58:  class TrkHitOnTrkIter : public std::iterator<std::random_access_iterator_tag,TrkHitOnTrk>
//
//
//FIXME: new idea: give TrkRepIter a constructor which tells it 
//                 to iterate either of all reps, or unique ones only
//                 (i.e. give the iterator a memory of those reps it
//                 has already seeni, or search the vector 'in front'
//                 of the current position?)
//                 Furthermore, give a reference to the TrkRecoTrk vector
//                 of sharep pointers of reps...

class TrkRep;

class TrkRepIter
{
public:
        typedef std::random_access_iterator_tag  iterator_category;
        // typedef const TrkRep                     value_type;
        typedef TrkRep                           value_type;
        typedef ptrdiff_t                        difference_type;
        typedef value_type*                      pointer;
        typedef value_type&                      reference;

        TrkRepIter(const TrkRepIter& rhs) : _l(rhs._l),_i(rhs._i) {};

        pointer   get()        const { return (*_l)[_i]; }
        pointer   operator->() const { return get(); }
        reference operator*()  const { return *get(); }

        TrkRepIter& operator-=(int i) { _i-=i; return *this; }
        TrkRepIter& operator+=(int i) { _i+=i; return *this; }
        TrkRepIter& operator--()      { --_i; return *this; } // prefix --
        TrkRepIter& operator++()      { ++_i; return *this; } // prefix ++

        TrkRepIter  operator--(int)   { TrkRepIter x(*this); --_i; return x; } // postfix --
        TrkRepIter  operator++(int)   { TrkRepIter x(*this); ++_i; return x; } // postfix ++


        ptrdiff_t operator-(const TrkRepIter& i) const { return _i - i._i; }
        bool operator==(const TrkRepIter& i) const { return _l.get() == i._l.get()
                                                                && _i == i._i; }
        bool operator!=(const TrkRepIter& i) const { return !operator==(i); }

private:
        friend class TrkRecoTrk; // allow TrkRecoTrk to create TrkRepIters -- and nobody else!
        typedef std::vector<pointer> TRL;
        TrkRepIter(boost::shared_ptr<TRL> l,int i) : _l(l),_i(i) {};

        boost::shared_ptr< std::vector<pointer> > _l;
        unsigned _i;
};
#endif

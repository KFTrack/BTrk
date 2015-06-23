// $Id: TrkHitOnTrkIter.hh,v 1.3 2005/11/28 05:44:02 kelsey Exp $
//
// "Adaptor" class which takes a class T, and presents a _guaranteed_ interface to T
// This is done with templates instead of inheritance to avoid (virtual) function 
// call overhead, as these functions are going to be called from lots of loops, 
// and this way the interface gets inlined and (if everything goes well) no 
// overhead is incurred for using this iterator class compared to using 
// a 'raw' T instead... (but then you couldn't deal with different T's...)
//
// There is a similar example in Stroustroup, which turns an iterator into a 'checked'
// iterator
//
// For the above to work, an implementation specific to 'T' has to be provided.
// In order to allow 'T' to specify which implementation should be used, it should
// provide a typename T::iterator_implementation that we can use here...
//
// To allow iterators over const TrkHitOnTrk and TrkHitOnTrk, 'T' should provide
// a T::value_type which is either a const TrkHitOnTrk, or a TrkHitOnTrk
//
// BTW, the implementations tend to iterate over pointers to the valuetype,
// and this class derefences the implementation to iterator over the valuetype instead...
//
#ifndef TRKHITONTRKITER_HH
#define TRKHITONTRKITER_HH
#include <iterator>
#include <cstddef>

//FIXME: only works on Linux:    class TrkHitOnTrkIter : public std::iterator_traits<TrkHitOnTrk *>
//FIXME: only works on Linux:    class TrkHitOnTrkIter : public std::random_access_iterator<const TrkHitOnTrk, ptrdiff_t>
//FIXME: only works on SunOS58:  class TrkHitOnTrkIter : public std::iterator<std::random_access_iterator_tag,TrkHitOnTrk>

class TrkHitOnTrk;

template <class T>
class TrkHitOnTrkIter
{
        // by using T::iterator_value_type as our value_type, we can re-use the same
        // code for a non-const iterator as a const iterator -- all that is needed
        // is a slightly different 'traits' class T to be passed here; The underlying
        // real iterator (the iterator_implementation) can be the same regardless of
        // const or non-const...

public:
        typedef std::random_access_iterator_tag       iterator_category;
        typedef typename T::iterator_value_type       value_type;
        typedef std::ptrdiff_t                             difference_type;
        typedef value_type*                           pointer;
        typedef value_type&                           reference;

        typedef typename T::iterator_implementation   iterator_implementation;

        TrkHitOnTrkIter<T>() :_i() {}; // create an invalid iter...
        TrkHitOnTrkIter<T>(const TrkHitOnTrkIter<T>& i) : _i(i._i) { }
        TrkHitOnTrkIter<T>(const iterator_implementation& i) : _i(i) { }

        pointer get() const { return *_i; }  // this function (together with * and ->)  is one of the main 
                                             // reasons for this class:
                                             //   most (all?) underlying containers contain pointers to
                                             //   TrkHitOnTrk, and we need to double-dereference to
                                             //   create the illusion of something that iterates over
                                             //   (const) TrkHitOnTrk instead of (const) TrkHitOnTrk*

        pointer   operator->() const { return this->get(); }
        reference operator*() const  { return *this->get(); }


        // next: forward all usual random access iterator operations 
        //       to the underlying actual implementation...

        TrkHitOnTrkIter<T>& operator-=(int i) { _i-=i; return *this; }
        TrkHitOnTrkIter<T>& operator+=(int i) { _i+=i; return *this; }

        TrkHitOnTrkIter<T>  operator-(int i)   { TrkHitOnTrkIter<T> x(_i); x-=i; return x; }
        TrkHitOnTrkIter<T>  operator+(int i)   { TrkHitOnTrkIter<T> x(_i); x+=i; return x; }

        TrkHitOnTrkIter<T>& operator--()      { --_i; return *this; } // prefix --
        TrkHitOnTrkIter<T>& operator++()      { ++_i; return *this; } // prefix ++

        TrkHitOnTrkIter<T>  operator--(int)   { TrkHitOnTrkIter<T> x(_i); --_i; return x; } // postfix --
        TrkHitOnTrkIter<T>  operator++(int)   { TrkHitOnTrkIter<T> x(_i); ++_i; return x; } // postfix ++

        ptrdiff_t operator-(const TrkHitOnTrkIter<T>& i) const { return _i - i._i; }
        bool operator==(const TrkHitOnTrkIter<T>& i) const { return _i==i._i; }
        bool operator!=(const TrkHitOnTrkIter<T>& i) const { return !operator==(i); }
        bool operator< (const TrkHitOnTrkIter<T>& i) const { return _i<i._i;}
        bool operator>=(const TrkHitOnTrkIter<T>& i) const { return !operator<(i);}
        bool operator> (const TrkHitOnTrkIter<T>& i) const { return _i>i._i;}
        bool operator<=(const TrkHitOnTrkIter<T>& i) const { return !operator>(i);}

private:
        iterator_implementation _i;
};
#endif

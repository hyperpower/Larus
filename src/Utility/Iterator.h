//
//  Iterator.h
//  LarusX
//
//  Created by zhou on 16/12/14.
//  Copyright (c) 2014 zhou. All rights reserved.
//

#ifndef _Iterator_h
#define _Iterator_h

#include "../TypeDef.h"

namespace Larus {
    //Define the iterator tags=====================
    struct input_iterator_tag {};
    struct output_iterator_tag {};
    struct forward_iterator_tag : public input_iterator_tag {};
    struct bidirectional_iterator_tag : public forward_iterator_tag {};
    struct random_access_iterator_tag : public bidirectional_iterator_tag {};
    
    template <class _Tp, class _Distance>
    struct input_iterator {
        typedef input_iterator_tag iterator_category;
        typedef _Tp                value_type;
        typedef _Distance          difference_type;
        typedef _Tp*               pointer;
        typedef _Tp&               reference;
    };
    
    struct output_iterator {
        typedef output_iterator_tag iterator_category;
        typedef void                value_type;
        typedef void                difference_type;
        typedef void                pointer;
        typedef void                reference;
    };
    
    template <class _Tp, class _Distance>
    struct forward_iterator {
        typedef forward_iterator_tag iterator_category;
        typedef _Tp                  value_type;
        typedef _Distance            difference_type;
        typedef _Tp*                 pointer;
        typedef _Tp&                 reference;
    };
    
    
    template <class _Tp, class _Distance>
    struct bidirectional_iterator {
        typedef bidirectional_iterator_tag iterator_category;
        typedef _Tp                        value_type;
        typedef _Distance                  difference_type;
        typedef _Tp*                       pointer;
        typedef _Tp&                       reference;
    };
    
    template <class _Tp, class _Distance>
    struct random_access_iterator {
        typedef random_access_iterator_tag iterator_category;
        typedef _Tp                        value_type;
        typedef _Distance                  difference_type;
        typedef _Tp*                       pointer;
        typedef _Tp&                       reference;
    };
    
    //template of iterator=============================
    template <class _Category,
              class _Tp,
              class _Distance = size_t,
              class _Pointer = _Tp*,
              class _Reference = _Tp&>
    struct iterator {
        typedef _Category  iterator_category;
        typedef _Tp        value_type;
        typedef _Distance  difference_type;
        typedef _Pointer   pointer;
        typedef _Reference reference;
    };
    //iterator_traits=================================
    template <class _Iterator>
    struct iterator_traits {
        typedef typename _Iterator::iterator_category iterator_category;
        typedef typename _Iterator::value_type        value_type;
        typedef typename _Iterator::difference_type   difference_type;
        typedef typename _Iterator::pointer           pointer;
        typedef typename _Iterator::reference         reference;
    };
    // for native pointer
    template <class _Tp>
    struct iterator_traits<_Tp*> {
        typedef random_access_iterator_tag  iterator_category;
        typedef _Tp                         value_type;
        typedef size_t                   difference_type;
        typedef _Tp*                        pointer;
        typedef _Tp&                        reference;
    };
    // for const native pointer
    template <class _Tp>
    struct iterator_traits<const _Tp*> {
        typedef random_access_iterator_tag iterator_category;
        typedef _Tp                         value_type;
        typedef size_t                   difference_type;
        typedef const _Tp*                  pointer;
        typedef const _Tp&                  reference;
    };
    //=========================
    template <class _Iter>
    inline typename iterator_traits<_Iter>::iterator_category
    iterator_category(const _Iter&) {
        typedef typename iterator_traits<_Iter>::iterator_category _Category;
        return _Category();
    }
    
    template <class _Iter>
    inline typename iterator_traits<_Iter>::difference_type*
    distance_type(const _Iter&) {
        return static_cast<typename iterator_traits<_Iter>::difference_type*>(0);
    }
    
    template <class _Iter>
    inline typename iterator_traits<_Iter>::value_type*
    value_type(const _Iter& __i) {
        return static_cast<typename iterator_traits<_Iter>::value_type*>(0);
    }
    //distance ==========================
    template <class _InputIterator>
    inline typename iterator_traits<_InputIterator>::difference_type
    __distance(_InputIterator __first,
               _InputIterator __last,
               input_iterator_tag)
    {
        typename iterator_traits<_InputIterator>::difference_type __n = 0;
        while (__first != __last) {
            ++__first; ++__n;
        }
        return __n;
    }
    
    template <class _RandomAccessIterator>
    inline typename iterator_traits<_RandomAccessIterator>::difference_type
    __distance(_RandomAccessIterator __first,
               _RandomAccessIterator __last,
               random_access_iterator_tag) {
        return __last - __first;
    }
    
    template <class _InputIterator>
    inline typename iterator_traits<_InputIterator>::difference_type
    distance(_InputIterator __first,
             _InputIterator __last) {
        typedef typename iterator_traits<_InputIterator>::iterator_category
        _Category;
        
        return __distance(__first, __last, _Category());
    }
    //advance===============================
    template <class _InputIter, class _Distance>
    inline void __advance(_InputIter& __i, _Distance __n, input_iterator_tag) {
        while (__n--) ++__i;
    }
    template <class _BidirectionalIterator, class _Distance>
    inline void __advance(_BidirectionalIterator& __i, _Distance __n,
                          bidirectional_iterator_tag) {
        if (__n >= 0)
            while (__n--) ++__i;
        else
            while (__n++) --__i;
    }
    template <class _RandomAccessIterator, class _Distance>
    inline void __advance(_RandomAccessIterator& __i,
                          _Distance __n,
                          random_access_iterator_tag) {
        __i += __n;
    }
    template <class _InputIterator, class _Distance>
    inline void advance(_InputIterator& __i, _Distance __n) {
        __advance(__i, __n, iterator_category(__i));
    }
}

#endif

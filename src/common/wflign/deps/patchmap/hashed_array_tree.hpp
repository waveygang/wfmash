#ifndef HASHED_ARRAY_TREE_H
#define HASHED_ARRAY_TREE_H

#include <algorithm>
#include <bitset>
#include <cassert>
#include <cstdint>
#include <cstring>
#include <functional>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <boost/container/allocator.hpp>
#include <boost/interprocess/allocators/allocator.hpp>
#include <typeinfo>
#include <exception>
#include <memory>
#include <utility>
#include "wmath_forward.hpp"
#include "wmath_bits.hpp"
#include "wmath_hash.hpp"
#include "wmath_math.hpp"

namespace wmath{
  template<
    typename T,
    class alloc = std::allocator<T>
  >
  class hashed_array_tree{
    public:
      typedef typename alloc::size_type size_type;
      typedef alloc allocator_type;
    private:
      size_type n = 0;
      T** data = nullptr;
      allocator_type allocator;
      size_type log2rootsize_n;
      size_type static constexpr log2rootsize(const size_type& n) {
        return (digits<size_type>()-clz(n)+1)/2;
      }
      size_type static constexpr rootsize(const size_type& n) {
        return n?(size_type(1)<<((digits<size_type>()-clz(n)+1)/2)):0;
      }
      size_type static constexpr rootsize_from_log2(
          const size_type& log2rootsize_n){
        return size_type(1)<<log2rootsize_n;
      }
      tuple<size_type,size_type> constexpr resolve(
          const size_type& i,
          const size_type& n,
          const size_type& log2rootsize_n
          ) const {
        assert(n);
        assert(i<n);
        const size_type m = (size_type(1)<<log2rootsize_n)-1;
        return {i>>log2rootsize_n,i&m};
      }
      tuple<size_type,size_type> constexpr resolve(
          const size_type& i,
          const size_type& n) const {
        const size_type log2roosize_n = log2rootsize(n);
        return resolve(i,n,log2rootsize_n);
      }
    public:
      ~hashed_array_tree() {
        if (data!=nullptr) {
          for (size_type i=0;i!=rootsize(n);++i) {
            allocator.deallocate(data[i],rootsize(n));
          }
          delete[] data;
        }
      }
      inline const T& operator[](const size_type& i) const {
        const auto r = resolve(i,n,log2rootsize_n);
        //cerr << i << " " << n << " " << get<0>(r) << " " << get<1>(r) << endl;
        assert(data!=nullptr);
        assert(get<0>(r)<rootsize(n));
        assert(data[get<0>(r)]!=nullptr);
        assert(get<1>(r)<rootsize(n));
        return data[get<0>(r)][get<1>(r)];
      }
      inline T& operator[](const size_type& i) {
        const auto r = resolve(i,n,log2rootsize_n);
        //cerr << i << " " << n << " " << get<0>(r) << " " << get<1>(r) << endl;
        assert(data!=nullptr);
        assert(get<0>(r)<rootsize(n));
        assert(data[get<0>(r)]!=nullptr);
        assert(get<1>(r)<rootsize(n));
        return data[get<0>(r)][get<1>(r)];
      }
      template<bool is_const>
      class const_noconst_iterator {
        friend class hashed_array_tree;
        private:
          typename conditional<is_const,const T**,T**>::type data;
          size_type position;
          const size_type log2rootsize;
        public:
          typedef typename alloc::difference_type difference_type;
          typedef typename alloc::value_type value_type;
          typedef typename
            conditional<is_const,
                        const typename alloc::reference,
                              typename alloc::reference
                       >::type
            reference;
          typedef typename
            conditional<is_const,
                        const typename alloc::pointer,
                              typename alloc::pointer
                       >::type
            pointer;
          typedef std::random_access_iterator_tag iterator_category;
          const_noconst_iterator(
            typename conditional<is_const,
                                 const hashed_array_tree*,
                                       hashed_array_tree*
                                >::type tree
                                )
            :data(tree->data),
             position(0),
             log2rootsize(hashed_array_tree::log2rootsize(tree->n))
          {}
          const_noconst_iterator(
              typename conditional<is_const,
                                 const hashed_array_tree*,
                                       hashed_array_tree*
                                >::type tree,
              const size_type& i
                                )
            :data(tree->data),
             position(i),
             log2rootsize(hashed_array_tree::log2rootsize(tree->n))
          {}
          template<bool is_const_other>
          bool operator==(
              const const_noconst_iterator<is_const_other>& o) const {
            if (position!=o.position) return false;
            if (data!=o.data) return false;
            return true;
          }
          template<bool is_const_other>
          bool operator!=(
              const const_noconst_iterator<is_const_other>& o) const {
            return !((*this)==o);
          }
          template<bool is_const_other>
          bool operator< (
              const const_noconst_iterator<is_const_other>& o) const {
            if (position < o.position) return true;
            else return false;
          }
          template<bool is_const_other>
          bool operator> (
              const const_noconst_iterator<is_const_other>& o) const {
            if (position > o.position) return true;
            else return false;
          }
          template<bool is_const_other>
          bool operator<=(
              const const_noconst_iterator<is_const_other>& o) const {
            if (position <=o.position) return true;
            else return false;
          }
          template<bool is_const_other>
          bool operator>=(
              const const_noconst_iterator<is_const_other>& o) const {
            if (position >=o.position) return true;
            else return false;
          }
          const_noconst_iterator<is_const>& operator++(){   // prefix
            ++position;
            return *this;
          }
          const_noconst_iterator<is_const> operator++(int){ // postfix
            iterator pre(*this);
            operator++();
            return pre;
          }
          const_noconst_iterator<is_const>& operator--(){   // prefix
            --position;
            return *this;
          }
          const_noconst_iterator<is_const> operator--(int){ // postfix
            iterator pre(*this);
            operator--();
            return pre;
          }
          const_noconst_iterator<is_const>& operator+=(const size_type& n){
            position+=n;
            return *this;
          }
          const_noconst_iterator<is_const> operator+(const size_type& n) const {
            return (const_noconst_iterator<is_const>(*this)+=n);
          }
          friend const_noconst_iterator<is_const> operator+(
              const size_type& n,
              const const_noconst_iterator<is_const>& it){
            return (const_noconst_iterator<is_const>(it)+=n);
          }
          const_noconst_iterator<is_const>& operator-=(const size_type& n){
            position-=n;
            return *this;
          }
          const_noconst_iterator<is_const> operator-(const size_type& n) const{
            return (const_noconst_iterator<is_const>(*this)-=n);
          }
          reference operator*() {
            const size_type m = (size_type(1)<<log2rootsize)-1;
            return data[position>>log2rootsize][position&m];
          }
          pointer operator->() {
            const size_type m = (size_type(1)<<log2rootsize)-1;
            return data[position>>log2rootsize]+(position&m);
          }
          reference operator*() const {
            const size_type m = (size_type(1)<<log2rootsize)-1;
            return data[position>>log2rootsize][position&m];
          }
          pointer operator->() const {
            const size_type m = (size_type(1)<<log2rootsize)-1;
            return data[position>>log2rootsize]+(position&m);
          }
      };
      typedef const_noconst_iterator<false> iterator;
      typedef const_noconst_iterator<true>  const_iterator;
      iterator begin(){
        return iterator(this);
      }
      const_iterator begin() const {
        return const_iterator(this);
      }
      const_iterator cbegin() const {
        return const_iterator(this);
      }
      iterator end() {
        return iterator(this,n);
      }
      const_iterator end() const {
        return const_iterator(this,n);
      }
      const_iterator cend() const {
        return const_iterator(this,n);
      }
      size_type constexpr max_size() const {
        return numeric_limits<size_type>::max();
      }
      const size_type size() const {
        return n;
      }
      size_type constexpr static next_size(const size_type& n) {
        return (1+(n-1)/rootsize(n))*rootsize(n);
      }
      void inline resize(const size_type& m) {
        //cerr << "resize hashed_array_tree from " << n << " to " << m << endl; 
        if (m==n) return;
        const size_type rootsize_n = rootsize_from_log2(log2rootsize_n);
        const size_type log2rootsize_m = log2rootsize(m);
        const size_type rootsize_m = rootsize(m);
        if (m==0) {
          for (size_type i=0;i!=rootsize_n;++i){
            allocator.deallocate(data[i],rootsize_n);
          }
          delete[] data;
          data = nullptr;
          n = m;
          log2rootsize_n = log2rootsize(n);
          return;
        }
        //cout << rootsize_n << " " << rootsize_m << endl;
        if (rootsize_m!=rootsize_n) {
          //cerr << "new index needed" << endl;
          T** new_data = new T*[rootsize_m];
          //cout << "first new[]" << endl;
          //cout << "new rootsize = " << rootsize_m << endl;
          for (size_type i=0;i!=rootsize_m;++i) new_data[i] = nullptr;
          size_type old_leaf = 0;
          size_type old_leafposition = 0;
          size_type new_leaf = 0;
          size_type new_leafposition = 0;
          for (size_type i=0;;++i) {
            if (new_leafposition == 0) {
              new_data[new_leaf] = allocator.allocate(rootsize_m);
              //cout << "new leaf" << endl;
            }
            if (i==min(m,n)) break;
            //cerr << "copying values" << endl;
            new_data[new_leaf][new_leafposition] =
              std::move(data[old_leaf][old_leafposition]);
            if (i+1==min(m,n)) break;
            if (++new_leafposition==rootsize_m) {
              ++new_leaf;
              new_leafposition = 0;
            }
            if (++old_leafposition==rootsize_n) {
              allocator.deallocate(data[old_leaf],rootsize_n);
              ++old_leaf;
              old_leafposition = 0;
            }
          }
          //cout << new_leaf << " " << (((m-1)>>log2rootsize_m)+1) << endl;
          for (++new_leaf;
              new_leaf!=(((m-1)>>log2rootsize_m)+1);
               ++new_leaf) {
            new_data[new_leaf] = allocator.allocate(rootsize_m);
            //cout << "new leaf" << endl;
          }
          //cout << "this were supposed to be all the leaves" << endl;
          if (data!=nullptr) {
            for (;old_leaf!=(((n-1)>>log2rootsize_n)+1);++old_leaf) {
              allocator.deallocate(data[old_leaf],rootsize_n);
            } 
            delete[] data;
          }
          data = new_data;
          n = m;
          log2rootsize_n = log2rootsize(n);
          //cout << "I think I am done oO " << endl;
        } else {
          //cout << "only new leaves needed" << endl;
          if (m>n) {
            //cout << ((n-1)>>log2rootsize_n)+1 << " "
            //     << ((m-1)>>log2rootsize_n)+1 << endl;
            for (size_type i = n?((n-1)>>log2rootsize_n)+1:0;
                           i!=(m?((m-1)>>log2rootsize_n)+1:0);
                         ++i) {
              //cout << i << " " << rootsize_n << endl;
              assert(i<rootsize_n);
              data[i] = allocator.allocate(rootsize_n);
            }
          } else {
            for (size_type i = m?((m-1)>>log2rootsize_n)+1:0;
                           i!=(n?((n-1)>>log2rootsize_n)+1:0);
                         ++i) {
              assert(i<rootsize_n);
              allocator.deallocate(data[i],rootsize_n);
            }
          }
          n = m;
          log2rootsize_n = log2rootsize(n);
        }
      }
  };
}
#endif // HASHED_ARRAY_TREE_H

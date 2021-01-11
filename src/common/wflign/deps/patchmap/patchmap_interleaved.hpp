#ifndef ORDERED_PATCH_MAP_H
#define ORDERED_PATCH_MAP_H

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

namespace whash{
  bool constexpr VERBOSE_PATCHMAP = false;
  using std::allocator_traits;
  using std::array;
  using std::cerr;
  using std::conditional;
  using std::cout;
  using std::enable_if;
  using std::endl;
  using std::false_type;
  using std::get;
  using std::index_sequence;
  using std::index_sequence_for;
  using std::initializer_list;
  using std::integral_constant;
  using std::is_const;
  using std::is_fundamental;
  using std::is_same;
  using std::is_trivially_copyable;
  using std::numeric_limits;
  using std::pair;
  using std::setw;
  using std::true_type;
  using std::tuple;
  using std::swap;

  template<class T>
  double frac(const T& n){
    return n*pow(0.5,numeric_limits<T>::digits);
  }

  template<typename T>
  struct dummy_comp{ // dummy comparator for when we don't need a comparator
    constexpr bool operator()(const T&,const T&) const {return false;}
  };

  template <typename T>
  constexpr size_t digits(const T& n=0){
    return numeric_limits<T>::digits;
  }

  template<typename T>
  typename std::enable_if<std::is_unsigned<T>::value,tuple<T,T>>::type
  constexpr long_mul(const T& a, const T& b);

  // calculate a * b = r0r1
  template<typename T>
  typename std::enable_if<std::is_unsigned<T>::value,tuple<T,T>>::type
  constexpr long_mul(const T& a, const T& b){
    const T N  = digits<T>()/2;
    const T t0 = (a>>N)*(b>>N);
    const T t1 = ((a<<N)>>N)*(b>>N);
    const T t2 = (a>>N)*((b<<N)>>N);
    const T t3 = ((a<<N)>>N)*((b<<N)>>N);
    const T t4 = t3+(t1<<N);
    const T r1 = t4+(t2<<N);
    const T r0 = (r1<t4)+(t4<t3)+(t1>>N)+(t2>>N)+t0;
    return {r0,r1};
  }
  
#ifdef __SIZEOF_INT128__
  template<>
  tuple<uint64_t,uint64_t>
  constexpr long_mul(const uint64_t& a, const uint64_t& b){
    unsigned __int128 r = ((unsigned __int128)(a))*((unsigned __int128)(b));
    return {r>>64,r};
  }
#endif

  template<>
  tuple<uint8_t,uint8_t> constexpr long_mul(const uint8_t& a,const uint8_t& b){
    const int_fast16_t r = int_fast16_t(a)*int_fast16_t(b);
    return {uint8_t(r>>8),uint8_t(r)};
  }
 
  template<>
  tuple<uint16_t,uint16_t> constexpr long_mul(
      const uint16_t& a,
      const uint16_t& b){
    const int_fast32_t r = int_fast32_t(a)*int_fast32_t(b);
    return {uint16_t(r>>16),uint16_t(r)};
  }
  
  template<>
  tuple<uint32_t,uint32_t> constexpr long_mul(
      const uint32_t& a,
      const uint32_t& b){
    const int_fast64_t r = int_fast64_t(a)*int_fast64_t(b);
    return {uint32_t(r>>32),uint32_t(r)};
  }
  
  template <typename T>
  constexpr size_t popcount(const T n){
    size_t c=0;
    while(n) (n&=(n-1),++c);
    return c;
  }

  constexpr size_t popcount(const uint32_t n){
    return __builtin_popcountl(n);
  }

  constexpr size_t popcount(const uint64_t n){
    return __builtin_popcountll(n);
  }

  template <typename T>
  typename std::enable_if<std::is_unsigned<T>::value,T>::type
  constexpr distribute(const T& a); // mix the hash value good, clmul_circ
                                    // with odious integer is suitable

  uint8_t  constexpr inline distribute(const uint8_t& a){
    return (a+111)*97;
  }

  uint16_t constexpr inline distribute(const uint16_t& a){
    return (a+36690)*43581;
  }

  uint32_t const distribute(uint32_t a){
    const uint32_t  b = 0x55555555ul;
    const uint32_t c0 = 3107070805ul;
    const uint32_t c1 = 3061963241ul;
    a = (a^(a>>16))*b;
    a = a^(a>>16);
    return a;
  }

  uint64_t constexpr distribute(uint64_t a){
    const uint64_t  b =   0x5555555555555555ull;
    const uint64_t c0 = 16123805160827025777ull;
    const uint64_t c1 = 13834579444137454003ull;
    const uint64_t c2 = 14210505232527258663ull;
    a^=(a>>32); a*=b; // a^=(a>>32);
    return a;
  }

  template<typename,typename=void>
  struct is_injective : false_type {};

  template<typename T>
  struct is_injective<T,typename enable_if<T::is_injective::value>::type>
  : true_type {};

  template<typename,typename=void>
  struct has_std_hash : false_type {};

  template<typename T>
  struct has_std_hash<T,decltype(std::hash<T>()(std::declval<T>()),void())>
  : true_type {};

  template<typename T>
  typename
  enable_if<has_std_hash<T>::value&&(!is_fundamental<T>::value),size_t>::type
  constexpr hash(const T& v){
    return std::hash<T>()(v);
  }

  size_t constexpr hash() {
    return 0;
  }

  size_t constexpr hash(const size_t& seed,const size_t& n) {
    return seed^distribute(n);
  }

  template<class K>
  typename enable_if<(sizeof(K)>sizeof(size_t)
                   &&is_fundamental<K>::value),size_t>::type
  constexpr hash(const K& k){
    size_t h(k);
    const size_t n = sizeof(K)/sizeof(size_t);
    for (size_t i=sizeof(size_t);i<sizeof(K);i+=sizeof(size_t))
      h = hash(h,size_t(k>>(i*CHAR_BIT)));
    return h;
  }

  uint8_t constexpr hash(const uint8_t& v){
    return v;
  }

  uint8_t constexpr hash(const int8_t& v){
    return v;
  }

  uint16_t constexpr hash(const uint16_t& v){
    return v;
  }

  uint16_t constexpr hash(const int16_t& v){
    return v;
  }

  uint32_t constexpr hash(const uint32_t& v){
    return v;
  }

  uint32_t constexpr hash(const int32_t& v){
    return v;
  }

  uint64_t constexpr hash(const uint64_t& v){
    return v;
  }

  uint64_t constexpr hash(const int64_t& v){
    return v;
  }
 
  template <typename T,typename... Rest>
  size_t constexpr hash(const T& v,Rest... rest);

  uint16_t constexpr hash(const uint8_t& v0,const uint8_t& v1){
    return (uint16_t(v0)<<8)^uint16_t(v1);
  }

  uint32_t constexpr hash(const uint16_t& v0,const uint16_t& v1){
    return (uint32_t(v0)<<16)^(uint32_t(v1));
  }
  
  uint64_t constexpr hash(const uint32_t& v0,const uint32_t& v1){
    return (uint64_t(v0)<<32)^(uint64_t(v1));
  }
  
  template<typename T,size_t... I>
  size_t constexpr hash_tuple_impl(const T& t, index_sequence<I...>){
    return hash(std::get<I>(t)...);
  }

  template<typename... Ts>
  size_t constexpr hash(const tuple<Ts...>& t){
    return hash_tuple_impl(t,index_sequence_for<Ts...>{});
  }

  template<typename T,size_t n>
  size_t constexpr hash(const array<T,n> a){
    size_t h(0);
    for (size_t i=0;i!=n;++i) {
      if constexpr(sizeof(T)<=sizeof(size_t))
        if (i%(sizeof(size_t)/sizeof(T))==0) h = distribute(h); 
      h = hash(h,a[i]);
      if constexpr(sizeof(T)>sizeof(size_t)) h = distribute(h);
    }
    return h;
  }
  
  template <typename T, typename... Rest>
  size_t constexpr hash(const T& v, Rest... rest) {
    return hash(hash(v),hash(rest...));
  }

  template<class K,class enable = void>
  struct hash_functor{
    typedef typename false_type::type is_injective;
    size_t operator()(const K& k) const {
      return hash(k);
    }
  };
  
  template<class K>
  struct hash_functor<
    K,
    typename enable_if<is_fundamental<K>::value,void>::type
  >{
    // if size_t has at least as many digits as the hashed type the hash can
    // be injective and I will make it so
    typedef typename integral_constant<bool,sizeof(K)<=sizeof(size_t)>::type
      is_injective;
    auto constexpr operator()(const K& k) const {
      return hash(k);
    }
  };
  
  template<typename T,size_t n>
  struct hash_functor<
    array<T,n>,void>
  {
    // if size_t has at least as many digits as the hashed type the hash can
    // be injective and I will make it so
    typedef typename integral_constant<bool,n*sizeof(T)<=sizeof(size_t)>::type
      is_injective;
    auto constexpr operator()(const array<T,n>& k) const {
      return hash(k);
    }
  };
  
  template<typename T0,typename T1,typename T2>
  constexpr T0 clip(const T0& n,const T1& l,const T2& h){
    return n<l?l:n>h?h:n;
  }
  
  template <typename T>
  typename std::enable_if<std::is_unsigned<T>::value,T>::type
  constexpr clz(const T x,const T lower=0,const T upper=digits<T>()){
    return (upper-lower==T(1))?digits<T>()-upper:
      (x&(T(0)-T(1)<<((upper+lower)/2))?
           clz(x,(upper+lower)/2,upper):
           clz(x,lower,(upper+lower)/2));
  }
 
  template <typename T>
  typename std::enable_if<std::is_unsigned<T>::value,T>::type
  constexpr ctz(const T x,const T lower=0,const T upper=digits<T>()){
    return
      (upper-lower==T(1))?lower:(x&(T(0)-T(1)<<((upper+lower)/2))?
          ctz(x,(upper+lower)/2,upper):
          ctz(x,lower,(upper+lower)/2));
    // TODO
  }


  template <typename T>
  typename std::enable_if<std::is_unsigned<T>::value,T>::type
  constexpr log2(const T x,const T lower=0,const T upper=digits<T>()){
    return (upper-lower==T(1))?lower:(x&(T(0)-T(1)<<((upper+lower)/2))?
           log2(x,(upper+lower)/2,upper):
           log2(x,lower,(upper+lower)/2));
  }

#if __GNUC__ > 3 || __clang__
  uint32_t constexpr clz(const uint32_t x){
    return x==0?32:__builtin_clz(x);
  }
  
  uint32_t constexpr ctz(const uint32_t x){
    return x==0?32:__builtin_ctz(x);
  }
  
  uint64_t constexpr clz(const uint64_t x){
    return x==0?64:__builtin_clzll(x);
  }
  
  uint64_t constexpr ctz(const uint64_t x){
    return x==0?64:__builtin_ctzll(x);
  }

  uint32_t constexpr log2(const uint32_t x){
    return x==0?0:31-__builtin_clz(x);
  }
  
  uint64_t constexpr log2(const uint64_t x){
    return x==0?0:63-__builtin_clzll(x);
  }
#endif
  
  template <typename T,typename S>
  typename std::enable_if<std::is_unsigned<T>::value,T>::type
  constexpr shl(const T n, const S i){
    if ((i<digits<T>())&&(i>=0)) return n<<i;
    return 0;
  }
  
  template <typename T,typename S>
  typename std::enable_if<std::is_unsigned<T>::value,T>::type
  constexpr shr(const T n, const S i){
    if ((i<digits<T>())&&(i>=0)) return n>>i;
    return 0;
  }

  template<class key_type    = int,  // int is the default, why not
           class mapped_type = int,  // int is the default, why not
           class hash        = hash_functor<key_type>,
           class equal       = std::equal_to<key_type>,
           class comp        = typename conditional<
             is_injective<hash>::value,
             dummy_comp<key_type>,
             std::less<key_type>>::type,
           bool dynamic      = true
          >
  class patchmap{
    public:
      typedef typename conditional<is_same<mapped_type,void>::value,
                                  std::true_type,
                                  mapped_type>::type
                                  _mapped_type;
      typedef pair<key_type,_mapped_type>  value_type;
      typedef value_type* value_pointer;
      typedef value_type& reference;
      typedef const value_type& const_reference;
      typedef size_t size_type;
      typedef uint64_t flag_type;
      typedef typename std::result_of<hash(key_type)>::type hash_type;
      static const size_type stride = 64;
      struct bucket {
        value_type values[stride];
        flag_type flag;
      };
    private:
      size_type num_data;
      size_type datasize;
      bucket* data;
      comp  comparator;
      equal equator;
      hash  hasher;
      using uphold_iterator_validity = false_type;
      // size_type const inline masksize() const {
      //  return (datasize+digits<size_type>()-1)/digits<size_type>();
      //}
      template<typename T>
      const key_type& key_of(T&& value) const {
        if constexpr (is_same<void,mapped_type>::value) return value;
        else return value.first;
      }
      size_type inline map(
          const hash_type& h,
          const hash_type& n
          ) const {
        return get<0>(long_mul(h,n));
      }
      size_type inline map(const hash_type& h) const {
        return map(h,datasize);
      }
      size_type inline map_diff(
          const hash_type& h0,
          const hash_type& h1,
          const hash_type& n
          ) const {
        const auto lm = long_mul(hash_type(h0-h1),n);
        return get<0>(lm);
      }
      size_type inline map_diff(
          const hash_type& h0,
          const hash_type& h1
          ) const {
        return map_diff(h0,h1,datasize);
      }
      size_type inline map_diff_round(
          const hash_type& h0,
          const hash_type& h1,
          const hash_type& n
          ) const {
        const auto lm = long_mul(hash_type(h0-h1),n);
        return get<0>(lm)+(get<1>(lm)>((~hash_type(0))>>1));
      }
      size_type inline map_diff_round(
          const hash_type& h0,
          const hash_type& h1
          ) const {
        return map_diff_round(h0,h1,datasize);
      }
      hash_type inline order(const key_type& k) const {
        return distribute(hasher(k));
      }
      inline const value_type& get_const_value_at(const size_type& i) const {
        return data[i/stride].values[i%stride];
      }
      inline value_type& get_value_at(const size_type& i) {
        return data[i/stride].values[i%stride];
      }
      inline const key_type& get_const_key_at(const size_type& i) const {
        return data[i/stride].values[i%stride].first;
      }
      inline key_type& get_key_at(const size_type& i) {
        return data[i/stride].values[i%stride].first;
      }
      inline const mapped_type& get_const_mapped_at(const size_type& i) const {
        return data[i/stride].values[i%stride].second;
      }
      inline mapped_type& get_mapped_at(const size_type& i) {
        return data[i/stride].values[i%stride].second;
      }
      bool inline is_less(
          const key_type& a,
          const key_type& b,
          const hash_type& oa,
          const hash_type& ob
          ) const {
        if constexpr (is_injective<hash>::value){
          assert(equator(a,b)==(oa==ob));
          if (oa<ob) return true;
          else       return false;
        } else {
          if (oa<ob) return true;
          if (oa>ob) return false;
          return comparator(a,b);
        }
      }
      bool inline is_less(
          const key_type& a,
          const key_type& b,
          const hash_type& oa
          ) const {
        return is_less(a,b,oa,order(b));
      }
      bool inline is_less(const key_type& a,const key_type& b) const {
        return is_less(a,b,order(a),order(b));
      }
      bool inline is_more(
          const key_type& a,
          const key_type& b,
          const hash_type& oa,
          const hash_type& ob
          ) const {
        if constexpr (is_injective<hash>::value){
          assert(equator(a,b)==(oa==ob));
          if (oa>ob) return true;
          else       return false;
        } else {
          if (oa>ob) return true;
          if (oa<ob) return false;
          return !((comparator(a,b))||(equator(a,b)));
        }
      }
      bool inline is_more(
          const key_type& a,
          const key_type& b,
          const hash_type& oa
          ) const {
        return is_more(a,b,oa,order(b));
      }
      bool inline is_more(
          const key_type& a,
          const key_type& b
          ) const {
        return is_more(a,b,order(a),order(b));
      }
      bool inline is_set(const size_type& n) const {        
        const size_type i = n/stride;
        const size_type j = n%stride;
        return data[i].flag&(flag_type(1)<<j);
      }
      void inline set(const size_type& n) {
        const size_type i = n/stride;
        const size_type j = n%stride;
        data[i].flag|=flag_type(1)<<j;
      }
      void inline unset(const size_type& n) {
        const size_type i = n/stride;
        const size_type j = n%stride;
        data[i].flag&=(~flag_type(0))^(flag_type(1)<<j);
      }
      void inline swap_set(const size_type& i,const size_type& j){
        if (is_set(i)==is_set(j)) return;
        if (is_set(i)){
          set(j);
          unset(i);
        }else{
          set(i);
          unset(j);
        }
      }
      bool inline index_key_is_less(const size_type& i,const key_type& k) const{
        if (is_set(i)) return is_less(get_const_key_at(i),k);
        return i<map(order(k));
      }
      bool inline key_index_is_less(const key_type& k,const size_type& i) const{
        if (is_set(i)) return is_less(k,get_const_key_at(i));
        return map(order(k))<i;
      }
      bool inline index_index_is_less(const size_type& i,const size_type& j)
        const {
        assert(i<datasize);
        assert(j<datasize);
        if (is_set(i)&&is_set(j))
          return is_less(get_const_key_at(i),get_const_key_at(j));
        if (is_set(i)) return map(order(get_const_key_at(i)))<j;
        if (is_set(j)) return i<map(order(get_const_key_at(j)));
        return i<j;
      }
      bool inline index_index_is_more(const size_type& i,const size_type& j)
        const {
        return index_index_is_less(j,i);
      }
      size_type inline find_first() const {
        size_type i=0;
        while (i!=datasize) if (is_set(i)) return i;
        return ~size_type(0);
      }
      // search for free bucket in decreasing order
      size_type inline search_free_dec(size_type i) const {
        while(true) {
          if (!is_set(i)) return i;
          if (i--==0) return ~size_type(0);
        }
        return ~size_type(0);
      }
      // search for free bucket in increasing order
      size_type inline search_free_inc(size_type i) const {
        while(true) {
          if (!is_set(i)) return i;
          if (i++==0) return ~size_type(0);
        }
        return ~size_type(0);
      }
      // search for free bucket bidirectional
      size_type inline search_free_bidir(const size_type& i) const {
        if (!is_set(i)) return i;
        for (size_t j=1;;++j) {
          if (i+j==datasize) {
            for (;i-j!=~size_type(0);++j) if (!is_set(i-j)) return i-j;
            return ~size_type(0);
          }
          if (i-j==~size_type(0)) {
            for (;i+j!=datasize     ;++j) if (!is_set(i+j)) return i+j;
            return ~size_type(0);
          }
          if (!is_set(i+j)) return i+j;
          if (!is_set(i-j)) return i-j;
        }
      }
      size_type const inline reserve_node(
          const  key_type&   k,
          const hash_type&  ok,
          const size_type& mok
          ){
        if (!is_set(mok)) {
          set(mok);
          ++num_data;
          return mok;
        }
        const size_type j = search_free_bidir(mok);
        assert(j<datasize);
        assert(!is_set(j));
        set(j);
        ++num_data;
        size_type i = j;
        while(true){
          if (i==0) break;
          if (!is_set(i-1)) break;
          if (is_less(get_key_at(i-1),k,order(get_key_at(i-1)),ok)) break; 
          swap(get_value_at(i),get_value_at(i-1));
          --i;
        }
        if (i!=j) return i;
        while(true){
          if (i+1>=datasize) break;
          if (!is_set(i+1)) break;
          if (is_less(k,get_key_at(i+1),ok,order(get_key_at(i+1)))) break; 
          swap(get_value_at(i),get_value_at(i+1));
          ++i;
        }
        return i;
      }
      size_type inline reserve_node(
          const key_type&   k,
          const hash_type& ok){
        const hash_type mok = map(ok);
        return reserve_node(k,ok,mok);
      }
      size_type inline reserve_node(const key_type& k){
        const hash_type  ok = order( k);
        const size_type mok =   map(ok);
        assert(mok<datasize);
        return reserve_node(k,ok,mok);
      }
      
      size_type inline interpol(
          const hash_type& ok,
          const hash_type& olo,
          const hash_type& ohi,
          const size_type& lo,
          const size_type& hi
          ) const {
        auto lm             = long_mul(size_type(ok)-size_type(olo),hi-lo);
        const size_type n   = clz(get<0>(lm));
        const size_type m   = digits<size_type>()-n;
        const hash_type den = (size_type(ohi)-size_type(olo))>>m;
        const hash_type nom = (get<0>(lm)<<n)+(get<1>(lm)>>m);
        return lo+nom/den;
      }

      size_type inline find_node_interpol(
        const  key_type&   k,
        const hash_type&  ok,
        const size_type& mok,
              size_type   lo,
              hash_type  olo,
              bool is_set_lo,
              size_type   hi,
              size_type  ohi,
              bool is_set_hi
          ) const {
        assert(lo<=hi||datasize==0);
        size_type mi;
        //size_t i = 0;
        while(true) {
          //if (i++>16) cout << lo << " " << hi << endl;
          //cerr << lo << " " << hi << " " << datasize << endl;
          //cerr << frac(olo) << " " << frac(ohi) << endl;
          //cerr << frac(order(get_const_key_at(lo))) << " "
          //     << frac(order(get_const_key_at(hi))) << endl;
          if (hi-lo<8) {
            if (hi-lo<4) {
              if (hi-lo<2) {
                if (hi-lo<1) {
                  if (is_set(lo))
                    if (equator(k,get_const_key_at(lo))) return lo;
                  return ~size_type(0);
                } else {
                  if (is_set_lo&&is_set_hi)
                    return ~size_type(0);
                  if (is_set(lo)) if (equator(k,get_const_key_at(lo)))
                    return lo;
                  if (is_set(hi)) if (equator(k,get_const_key_at(hi)))
                    return hi;
                  return ~size_type(0);
                }
              } else {
                mi = lo + ((hi-lo)>>1);
              }
            } else {
              if (is_set_hi && is_set_lo) {
                mi = lo + ((hi-lo)>>1);
              } else if (is_set_lo) {
                mi = lo + ((hi-lo+2)>>2);
              } else if (is_set_hi) {
                mi = hi - ((hi-lo+2)>>2);
              } else {
                return ~size_type(0);
              } 
            }
          } else {
            if (is_set_hi && is_set_lo) {
              mi = interpol(ok,olo,ohi,lo,hi);
            } else if (is_set_lo) {
              const size_type st = map_diff(ok,olo);
              mi = lo+st<hi?lo+st:hi;
            } else if (is_set_hi) {
              const size_type st = map_diff(ohi,ok);
              mi = lo+st<hi?hi-st:lo;
            } else {
              return ~size_type(0);
            }
            mi = clip(mi,lo+1,hi-1);
          }
          if (!is_set(mi)) {
            if (mi<mok) {
              lo = mi;
              is_set_lo=false;
              continue;
            }
            if (mi>mok) {
              hi = mi;
              is_set_hi=false;
              continue;
            }
            return ~size_type(0);
          }
          if (equator(k,get_const_key_at(mi))) return mi;
          const hash_type omi = order(get_const_key_at(mi));
          if (ok<omi) {
            hi = mi;
            ohi = omi;
            is_set_hi = true;
            continue;
          }
          if (ok>omi) {
            lo = mi;
            olo = omi;
            is_set_lo = true;
            continue;
          }
          if constexpr (is_injective<hash>::value) {
            return ~size_type(0);
          } else {
            if (comparator(k,get_const_key_at(mi))) {
              hi = mi;
              ohi = omi;
              is_set_hi = true;
              continue;
            }
            if (comparator(get_const_key_at(mi),k)) {
              lo = mi;
              olo = omi;
              is_set_lo = true;
              continue;
            }
          }
          return ~size_type(0);
        }
        return ~size_t(0);
      }

      size_type inline find_node(
          const key_type &   k,
          const hash_type&  ok,
          const size_type& mok)
        const {
        assert((mok<datasize)||(datasize==0));
        if (datasize==0) return ~size_type(0);
        //cout << "find_node " << frac(ok) << " " << datasize << endl;
        if (!is_set(mok)) return ~size_type(0);
        if (equator(get_const_key_at(mok),k)) return mok;
        const hash_type omi = order(get_const_key_at(mok));
        if (omi<ok) {
          return find_node_interpol(k,ok,mok,
              mok       ,omi          ,true ,
              datasize-1,~size_type(0),false);
        } else {
          return find_node_interpol(k,ok,mok,
              0         ,0            ,false,
              mok       ,          omi,true );
        }
      }
      
      size_type const inline find_node(
          const  key_type&  k,
          const size_type& ok
          ) const { return find_node(k,ok,map(ok)); }
      
      size_type const inline find_node(const key_type& k)
      const { return find_node(k,order(k)); }

      size_type const inline find_node_bruteforce(const key_type& k) const {
        for (size_type i=0;i!=datasize;++i)
          if (is_set(i)) if (equator(get_const_key_at(i),k)) return i;
        return ~size_type(0);
      }

      template<typename map_type>
      typename conditional<is_const<map_type>::value,
        const mapped_type&,
              mapped_type&>::type
      static inline const_noconst_at(map_type& hashmap,const key_type& k) {
        size_type i = hashmap.find_node(k);
        if (i<hashmap.datasize){
          assert(hashmap.is_set(i));
          return hashmap.data[i/8].values[i%8].second;
        } else throw std::out_of_range(
            std::string(typeid(hashmap).name())
            +".const_noconst_at("+typeid(k).name()+" k)"
            +"key not found, array index "
            +std::to_string(i)+" out of bounds"
           );
      }
      void const resize_out_of_place(const size_type& n){
        //cout << "resize out of place" << endl;
        size_type old_datasize = n;
        bucket* old_data = new bucket[(old_datasize+stride-1)/stride]; 
        for (size_type i=0;i!=(old_datasize+stride-1)/stride;++i)
          old_data[i].flag=0;
        num_data = 0;
        swap(old_data,data);
        swap(old_datasize,datasize);
        for (size_type n=0;n<old_datasize;++n) {
          const size_type i = n/stride;
          const size_type j = n%stride;
          if (old_data[i].flag&(flag_type(1)<<j)) {
            const key_type key = old_data[i].values[j].first;
            const size_type l = reserve_node(key);
            get_value_at(l) = old_data[i].values[j];
            //cout << "copied " << frac(order(get_const_key_at(l))) << endl;
            set(l);
          }
        }
        delete[] old_data;
      }
    public:
      // constructor
      patchmap(const size_type& datasize = 0)
        :datasize(datasize)
      {
        num_data = 0;
        if (datasize) data = new bucket[(datasize+stride-1)/stride];
        else          data = nullptr;
        for (size_type i=0;i!=(datasize+stride-1)/stride;++i) data[i].flag=0;
      }
      ~patchmap(){                                  // destructor
        delete[] data;
      }
      patchmap(patchmap&& other) noexcept  // move constructor
      {
        data = nullptr;
        datasize = 0;
        swap(data,other.data);
        swap(datasize,other.datasize);
      }
      template<
        class key_type_other,
        class mapped_type_other,
        class hash_other,
        class equal_other,
        class comp_other,
        class alloc_other
              >
      inline patchmap& operator=                   // copy assignment
        (const patchmap<
           key_type_other,
           mapped_type_other,
           hash_other,
           equal_other,
           comp_other
           //alloc_other
         >& other)
      {
        typedef patchmap<
           key_type_other,
           mapped_type_other,
           hash_other,
           equal_other,
           comp_other
           //alloc_other
         > other_type;
        num_data = other.num_data;
        datasize = other.datasize;
        if (datasize) data = new bucket[(datasize+stride-1)/stride];
        else data = nullptr;
        for (size_type i=0;i!=(datasize+stride-1)/stride;++i) data[i].flag=0;
        if constexpr (
            is_same<hash , hash_other>::value
          &&is_same<equal,equal_other>::value
          &&is_same<comp , comp_other>::value
          ){
          if constexpr (
              is_trivially_copyable<value_type>::value
            &&is_same<value_type,typename other_type::value_type>::value)
            memcpy(reinterpret_cast<void*>(data),
                   reinterpret_cast<void*>(other.data),
                   datasize*sizeof(value_type));
          else for (size_type i=0;i!=datasize;++i)
            get_value_at(i)=other.get_value_at(i);;
        } else {
          for (auto it=other.begin();it!=other.end();++it) insert(*it);
        }
      }
      patchmap(const patchmap& other){
        num_data = other.num_data;
        datasize = other.datasize;
        if (datasize) data = new bucket[(datasize+stride-1)/stride];
        else data = nullptr;
        if constexpr (is_trivially_copyable<value_type>::value) {
            memcpy(reinterpret_cast<void*>(data),
                   reinterpret_cast<void*>(other.data),
                   (datasize+7)/8*sizeof(bucket));
        } else {
          for (size_type i=0;i!=(datasize+stride-1)/stride;++i) data[i].flag=0;
          for (size_type i=0;i!=datasize;++i)
            if (is_set(i)) get_value_at(i)=other.get_value_at(i);
        }
      }
      inline patchmap& operator=                   // copy assignment
        (const patchmap& other)
      {
        return *this = patchmap(other);
      }
      inline patchmap& operator=                   // move assignment
        (patchmap&& other)
        noexcept{
        swap(data,other.data);
        swap(datasize,other.datasize);
        return *this;
      }
      size_type erase(
          const  key_type&   k,
          const hash_type&  ok,
          const size_type& mok){
        size_type i = find_node(k,ok,mok);
        if (i>=datasize) return 0;
        //cout << "erasing " << frac(ok) << endl;
        //cout << "found at " << i << endl;
        const size_type j = i;
        while(true){
          if (i+1==datasize) break;
          if (!is_set(i+1)) break;
          if (map(order(get_const_key_at(i+1)))>i) break;
          swap(get_value_at(i),get_value_at(i+1));
          ++i;
        }
        if (i==j){
          while(true){
            if (i==0) break;
            if (!is_set(i-1)) break;
            if (map(order(get_const_key_at(i-1)))<i) break;
            swap(get_value_at(i),get_value_at(i-1));
            --i;
          }
        }
        unset(i);
        get_value_at(i)=value_type();
        --num_data;
        assert(num_data<datasize);
        return 1;
      }
      size_type erase(const  key_type& k,const size_type& ok){
        const hash_type hint = map(ok);
        return erase(k,ok,hint);
      }
      size_type erase(const key_type & k){
        const size_type ok = order(k);
        return erase(k,ok);
      }
      void inline clear(){
        for (size_type i=0;i!=datasize/stride;++i) {
          for (size_type j=0;j!=stride;++j) data[i].values[j] = value_type();
          data[i].flag=0;
        }
        num_data=0;
      }
      void const resize(const size_type& n){
        resize_out_of_place(n); return;
      }
      size_type inline size() const { return num_data; }
      size_type const test_size() const {
        size_type test = 0;
        for (size_type i=0;i!=datasize;++i) test += is_set(i);
        return test;
      }
      void inline ensure_size() {
        if constexpr (!dynamic) return;
#if defined PATCHMAP_EXPANSIVE
        //if (num_data*9<datasize*7) return; 
        if (num_data*5<datasize*4) return;
        if (datasize) resize(stride*(((12*datasize+6)/7+stride-1)/stride));
        //if (datasize) resize((7*datasize+2)/4);
        else resize(256);
        return;
#endif
        if (num_data*8 < datasize*7 ) return;
        size_type nextsize;
        if (datasize < 257) {
          if (datasize == 0) nextsize = 8;
          else nextsize = 2*datasize;
        } else {
          //nextsize = 50*datasize/31;
          //nextsize = 48*datasize/31;
          nextsize = 47*datasize/31;
          //nextsize = 45*datasize/31;
          nextsize = (nextsize+stride-1)/stride;
          nextsize*= stride;
        }
        resize(nextsize);
      }
      mapped_type& operator[](const key_type& k){
        const size_type i = find_node(k);
        if (i<datasize) return get_mapped_at(i);
        ensure_size();
        const size_type j = reserve_node(k);
        get_value_at(j) = {k,_mapped_type()};
        return get_mapped_at(j);
      }
      const mapped_type& operator[](const key_type& k) const {
        const size_type i = find_node(k);
        if (i<datasize) return get_const_mapped_at(i);
        else throw std::out_of_range(
            std::string(typeid(*this).name())
            +".const_noconst_at("+typeid(k).name()+" k)"
            +"key not found, array index "
            +std::to_string(i)+" out of bounds"
           );
      }
      mapped_type& at(const key_type& k){
        return const_noconst_at(*this,k);
      }
      const mapped_type& at(const key_type& k) const {
        return const_noconst_at(*this,k);
      }
      size_type const inline count(const key_type& k) const {
        return (find_node(k)<datasize);
      }
      template<class key_type_other,
               class mapped_type_other,
               class hash_other,
               class equal_other,
               class comp_other
               //class alloc_other
              >
      bool operator==(
          const patchmap<
            key_type_other,
            mapped_type_other,
            hash_other,
            equal_other,
            comp_other
            //alloc_other
            >& other)
      const {
        if (datasize!=other.datasize) return false;
        if constexpr (
            is_same<hash , hash_other>::value
          &&is_same<equal,equal_other>::value
          &&is_same<comp , comp_other>::value
          ){
          auto it0 = begin();
          auto it1 = other.begin();
          while (true){
            if (it0==end()) return true;
            if ((*it0)!=(*it1)) return false;
            ++it0;++it1;
          }
        } else {
          for (auto it=other.begin();it!=other.end();++it){
            if (nount(it->first)) if (at(it->first)==it->second) continue;
            return false;
          }
          return true;
        }
      }
      template<class key_type_other,
               class mapped_type_other,
               class hash_other,
               class equal_other,
               class comp_other
               //class alloc_other
              >
      bool operator!=(
          const patchmap<
            key_type_other,
            mapped_type_other,
            hash_other,
            equal_other,
            comp_other
            //alloc_other
            >& o)
      const{ return !((*this)==o); }
      equal key_eq() const{ // get key equivalence predicate
        return equal{};
      }
      comp key_comp() const{ // get key order predicate
        return comp{};
      }
      /*
      alloc get_allocator() const{
        return allocator;
      }
      */
      hash hash_function() const{ // get hash function
        return hash{};
      }  
      template<bool is_const>
      class const_noconst_iterator {
        friend class patchmap;
        public:
          size_type hint;
          key_type key;
          typename conditional<is_const,const patchmap&,
                                              patchmap&>::type map;
        private:
          void inline update_hint() {
            if constexpr (!uphold_iterator_validity::value) return;
            if (hint<map.datasize)
              if (equator(map.get_const_key_at(hint),key)) return;
            hint = map.find_node(key,hint);
            if (hint>=map.datasize) hint = ~size_type(0);
          }
          void inline unsafe_increment() { // assuming hint is valid
            while (++hint<map.datasize) if (map.is_set(hint)) return;
          }
          void inline unsafe_decrement(){ // assuming hint is valid
            while (--hint<map.datasize) if (map.is_set(hint)) return;
          }
        public:
          //typedef typename alloc::difference_type difference_type;
          typedef ptrdiff_t difference_type;
          typedef pair<key_type,_mapped_type> value_type;
          typedef typename
            conditional<is_const,const value_type&,value_type&>::type
            reference;
          typedef typename
            conditional<is_const,const value_type*,value_type*>::type
            pointer;
          typedef std::bidirectional_iterator_tag iterator_category;
          const_noconst_iterator(){
            //cout << "constructor 0" << endl;
          }
          const_noconst_iterator(
            const size_t& hint,
            typename conditional<is_const,
                                 const patchmap*,
                                       patchmap*
                                >::type map)
            :hint(hint),key(key_type{}),map(map){
            //cout << "constructor 1 " << hint << endl;
          }
          const_noconst_iterator(
            const size_t& hint,
            const key_type& key,
            typename conditional<is_const,
                                 const patchmap*,
                                       patchmap*
                                >::type map)
            :hint(hint),key(key),map(map) {
              //cout << "constructor 2 " << hint << endl;
          }
          ~const_noconst_iterator(){
          //cout << "destructor of const_noconst_iterator " << is_const << endl;
          //cout << hint << endl;
          }
          // copy constructor
          template<bool is_const_other>
          const_noconst_iterator(const const_noconst_iterator<is_const_other>& o)
          :hint(o.hint),key(o.key),map(o.map){
            //cout << "copy constructor" << endl;
          }
          // move constructor
          template<bool is_const_other>
          const_noconst_iterator(
              const_noconst_iterator<is_const_other>&& o) noexcept{
            //cout << "move constructor" << endl;
            swap(hint,o.hint);
            swap(key,o.key);
            swap(map,o.map);
          }
          // copy assignment
          template<bool is_const_other>
          const_noconst_iterator<is_const>& operator=(
              const const_noconst_iterator<is_const_other>& other){
            //cout << "copy assignment" << endl;
            return  (*this=const_noconst_iterator<is_const>(other));
          }
          template<bool is_const_other>
          bool operator==(
              const const_noconst_iterator<is_const_other>& o) const {
            if ((hint>=map.datasize)&&(o.hint>=o.map.datasize)) return true;
            if ((hint>=map.datasize)||(o.hint>=o.map.datasize)) return false;
            // comparisons are only valid for iterators on same container
            // if ((&map)!=(&o.map)) return false;
            if (key!=o.key) return false;
            return true;
          }
          template<bool is_const_other>
          bool operator!=(
              const const_noconst_iterator<is_const_other>& o) const{
            return !((*this)==o);
          }
          template<bool is_const_other>
          bool operator< (
              const const_noconst_iterator<is_const_other>& o) const{
            if ((o.hint<o.map.datasize)){
              if (hint<map.datasize){
                return comp(key,o.key);
              }else{
                return false;
              }
            } else {
              return false;
            }
          }
          template<bool is_const_other>
          bool operator> (
              const const_noconst_iterator<is_const_other>& o) const{
            if ((o.hint<o.map.datasize)){
              if (hint<map.datasize){
                return (!comp(key,o.key))&&(!equal(key,o.key));
              }else{
                return true;
              }
            } else {
              return false;
            }
          }
          template<bool is_const_other>
          bool operator<=(
              const const_noconst_iterator<is_const_other>& o) const{
            if ((o.hint<o.map.datasize)){
              if (hint<map.datasize){
                return comp(key,o.key)||equal(key,o.key);
              }else{
                return false;
              }
            } else {
              return true;
            }
          }
          template<bool is_const_other>
          bool operator>=(
              const const_noconst_iterator<is_const_other>& o) const{
            if ((o.hint<o.map.datasize)){
              if (hint<map.datasize){
               return !comp(key,o.key);
              }else{
                return true;
              }
            } else {
              return true;
            }
          }
          const_noconst_iterator<is_const>& operator++(){   // prefix
            update_hint();
            unsafe_increment();
            return *this;
          }
          const_noconst_iterator<is_const> operator++(int){ // postfix
            update_hint();
            iterator pre(*this);
            unsafe_increment();
            return pre;
          }
          const_noconst_iterator<is_const>& operator--(){   // prefix
            update_hint();
            unsafe_decrement();
            return *this;
          }
          const_noconst_iterator<is_const> operator--(int){ // postfix
            update_hint();
            iterator pre(*this);
            unsafe_decrement();
            return pre;
          }
          reference operator*() {
            update_hint();
            return map.data[hint/8].values[hint%8];
          }
          pointer operator->() {
            update_hint();
            return &(map.data[hint/8].values[hint%8]);
          }
          reference operator*() const {
            if (hint<map->datasize) if (equator(map.get_const_key_at(hint),key))
              return map.data[hint/stride].values[hint%stride];
            const size_type i = map.find_node(key);
            return map.data[i/stride].value[i%stride];
          }
          pointer operator->() const {
            return &(operator*());
          }
    };
    typedef const_noconst_iterator<false> iterator;
    typedef const_noconst_iterator<true>  const_iterator;    
    iterator begin() {
      const size_type i = find_first();
      return iterator(i,get_key_at(i),this);
    }
    const_iterator begin() const {
      const size_type i = find_first();
      return const_iterator(i,get_const_key_at(i),this);
    }
    const_iterator cbegin() const {
      const size_type i = find_first();
      return const_iterator(i,get_const_key_at(i),this);
    }
    iterator end() {
      const size_type i = find_first();
      return iterator(~size_type(0),this);
    }
    const_iterator end() const {
      return const_iterator(~size_type(0),this);
    }
    const_iterator cend() const {
      return const_iterator(~size_type(0),this);
    }
    // void swap(unpatchmap&); // TODO
    size_type max_size()         const {
      return std::numeric_limits<size_type>::max();
    }
    bool empty()                 const {return (num_data==0);}
    size_type bucket_count()     const {return datasize;}
    size_type max_bucket_count() const {
      return std::numeric_limits<size_type>::max();
    }
    void rehash(const size_type& n) { if (n>=size()) resize(n); }
    void reserve(const size_type& n){ if (3*n>=2*(size()+1)) resize(n*3/2); }
    pair<iterator,bool> insert ( const value_type& val ){
      const size_type i = find_node(key_of(val));
      if (i<datasize) return {iterator(i,key_of(val),this),false};
      ensure_size();
      const size_type j = reserve_node(key_of(val));
      get_value_at(j) = val;
      return {{j,key_of(val),this},true};
    }
    template <class P>
    pair<iterator,bool> insert ( P&& val ){
      const size_type i = find_node(key_of(val));
      if (i<datasize) return {iterator(i,key_of(val),this),false};
      ensure_size();
      const size_type j = reserve_node(key_of(val));
      if constexpr (is_same<void,mapped_type>::value)
        get_value_at(j) = {val,true_type{}};
      else
        get_value_at(j) = val;
      return {{j,key_of(val),this},true};
    }
    iterator insert ( const_iterator hint, const value_type& val ){
      return insert(val); // hint is useless
    }
    template <class P>
    iterator insert ( const_iterator hint, P&& val ) {
      return insert(val); // hint is useless
    }
    template <class InputIterator>
    void insert ( InputIterator first, InputIterator last ){
      for (auto it(first);it!=last;++it){
        insert(*it);
      }
    }
    void insert ( initializer_list<value_type> il ){
      insert(il.begin(),il.end());
    }
    template <class... Args>
    pair<iterator, bool> emplace ( Args&&... args ){
      insert(value_type(args...));
    }
    template <class... Args>
    iterator emplace_hint(const_iterator position,Args&&... args){
      insert(position,value_type(args...));
    }
    pair<iterator,iterator> equal_range(const key_type& k){
      const size_type i = find_node(k);
      if (i>=datasize) return {end(),end()};
      iterator lo(i,get_const_key_at(i),this);
      iterator hi(lo);
      ++hi;
      return {lo,hi};
    }
    pair<const_iterator,const_iterator>
    equal_range ( const key_type& k ) const {
      const size_type i = find_node(k);
      if (i>=datasize) return {cend(),cend()};
      iterator lo(i,get_const_key_at(i),this);
      iterator hi(lo);
      ++hi;
      return {lo,hi};
    }
    float load_factor() const noexcept{
      return float(num_data)/float(datasize);
    }
    float const max_load_factor() const noexcept {
      return 1.0;
    }
    template<bool is_const>
    iterator erase(const_noconst_iterator<is_const> position){
      iterator it(position);
      ++it;
      erase(position.key);//,position.hint);
      return it;
    }
    template<bool is_const>
    iterator erase(
        const_noconst_iterator<is_const> first,
        const_noconst_iterator<is_const> last){
      for (auto it=first;it!=last;it=erase(it));
    }
    [[deprecated(
        "disabled for performance reasons"
    )]] void max_load_factor(float z) {
      // m = number of elements ; n = number of buckets
      // n*load_factor >= m
      // n*mul_n >= m*mul_m
      size_type mul_n =  ceil(z*16);
      size_type mul_m = floor(z*16);
    }
  }; 

  /* TODO
  template<class K,
           class T,
           class hash=hash_functor<K>,
           class equal = std::equal_to<K>,
           class comp = std::less<K>,
           class A = std::allocator<std::pair<K,T>>
          >
  void swap(patchmap<K,K,hash,equal,comp,A>&,
            patchmap<K,K,hash,equal,comp,A>&);
  */
  
  template<class key_type    = int,  // int is the default, why not
           class mapped_type = int,  // int is the default, why not
           class hash        = hash_functor<key_type>,
           class equal       = std::equal_to<key_type>,
           class comp        = typename conditional<is_injective<hash>::value,
                                                    dummy_comp<key_type>,
                                                    std::less<key_type>>::type
           /*class alloc       = typename boost::container::allocator<
             typename conditional<
               std::is_same<mapped_type,void>::value,
               std::pair<key_type,true_type>,
               std::pair<key_type,mapped_type>
             >::type,2>*/
          >
  using static_patchmap =
    patchmap<key_type,mapped_type,hash,equal,comp,false>;
  
  template<class key_type,           // unordered_map has no default key_type
           class mapped_type,        // unordered_map has no default mapped_type
           class hash        = hash_functor<key_type>,
           class equal       = std::equal_to<key_type>,
           //class alloc       = typename // mapped_type must not be void
           //  boost::container::allocator<std::pair<key_type,mapped_type>,2>,
           class comp        = typename conditional<is_injective<hash>::value,
             dummy_comp<key_type>,typename std::less<key_type>::type>::type
          >
  using unordered_map =
    patchmap<key_type,mapped_type,hash,equal,comp>;
  
  template<class key_type,           // unordered_set has no default key_type
           class hash        = hash_functor<key_type>,
           class equal       = std::equal_to<key_type>,
           //class alloc       = typename
           //  boost::container::allocator<pair<key_type,true_type>,2>,
           class comp        = typename conditional<is_injective<hash>::value,
             dummy_comp<key_type>,typename std::less<key_type>::type>::type
          >
  using unordered_set =
    patchmap<key_type,void,hash,equal,comp>;
  
}
#endif // ORDERED_PATCH_MAP_H

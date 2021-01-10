#include <iomanip>
#include <iostream>
#include <limits>
#include <algorithm>
#include <iostream>

using std::swap;
using std::min;
using std::cerr;
using std::cout;
using std::endl;
using std::numeric_limits;
using std::get;
using std::setw;
using std::fixed;
using std::setprecision;

constexpr size_t VERBOSITY = 0;

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

class fibonacci_groth_policy{
  private:
    size_t n0 = 1;
  public:
    const size_t operator()(const size_t& n){
      size_t n1 = 0;
      if (n>n0){
        n1=n0+n;
        n0=n;
      } else{
        n0 =n+1;
        n1 =3*(n+2)/2;
      }
      return n1;
    }
  public:
    const bool operator()(
        const size_t& n_elem,
        const size_t& n_data){
      return n_elem*8 >= n_data*7;
    }
};

double frac(const size_t& n){
  return n*pow(0.5,numeric_limits<size_t>::digits);
}

class ordered_patch_map_v0{
  private:
    fibonacci_groth_policy grow;
    bool *       mask;
    size_t (*data)[2];
    size_t n_elem = 0;
    size_t n_data = 0;
  public:
    void print(){
      cerr << n_elem << " " << n_data << endl;
      for (size_t i=0;i!=n_data;++i){
        cout << fixed << setprecision(16);
        if (mask[i]) cout << setw(6) << i;
        else         cout << "      ";
        cout << setw(20) << frac(data[i][0])
             << setw(20) << frac(data[i][1]);
        if (mask[i]) cout << setw(8) << map(data[i][0])
                          << setw(8) << int(map(data[i][0]))-int(i);
        else         cout << setw(8) << i
                          << setw(8) << 0;
        cout << endl;
      }
      cout << endl;
    }
  private:
    const size_t map(const size_t& n) const {
      return get<0>(wmath::long_mul(n,n_data));
    }
    bool find_element(const size_t& key,size_t& n) const {
      if (VERBOSITY) cout << "find element " << frac(key) << endl;
      if (n_data==0) return false;
      const size_t patchsize = 1+(n_data+1)/(n_data+1-n_elem);
      size_t lower = 0;
                 n = map(key);
      size_t upper = n_data-1;
      while (true){
        if (VERBOSITY){  cout << setw( 8) << lower
             << setw( 8) << n
             << setw( 8) << upper
             << setw( 4) << mask[n]
             << setw(20) << frac(data[n][0])
             << setw(20) << frac(key)
             << setw(20) << map(key)
             << endl;
        }
        if (mask[n]) {
          if (data[n][0]==key) return true;
        }
        if (!(lower<upper)) break;
        bool is_less;
        if (mask[n]) is_less = data[n][0]>key;
        else         is_less = n>map(key); 
        if (is_less){
          if (VERBOSITY) cout << "is_less" << endl;
          if (n==0) break;
          upper = n-1;
          n    -= min((upper+2-lower)/2,patchsize);
        } else {
          if (VERBOSITY) cout << "is_more" << endl;
          if (n==n_data-1) break;
          lower = n+1;
          n    += min((upper+2-lower)/2,patchsize);
        }
      }
      return false;
    }
    size_t insert_element(const size_t& key){
      ++n_elem;
      size_t estim = map(key);
      for (size_t i=0;;++i){
        if (estim+i<n_data) if (!mask[estim+i]){
          estim+=i;
          break;
        }
        if (estim>=i)       if (!mask[estim-i]){
          estim-=i;
          break;
        }
        if ((estim<i)&&(estim+i>=n_data)) break;
      }
      mask[estim]=true;
      data[estim][0] = key;
      while (estim+1<n_data) {
        if (!mask[estim+1]) break;
        if (data[estim][0]<data[estim+1][0]) break;
        swap(data[estim][0],data[estim+1][0]);
        swap(data[estim][1],data[estim+1][1]);
        ++estim;
      }
      while (estim>0) {
        if (!mask[estim-1]) break;
        if (data[estim][0]>data[estim-1][0]) break;
        swap(data[estim][0],data[estim-1][0]);
        swap(data[estim][1],data[estim-1][1]);
        --estim;
      }
      return estim;
    }
  public:
    ordered_patch_map_v0(){
      mask = (bool*)malloc(0);
      data = (size_t(*)[2])malloc(0);
    }
    ~ordered_patch_map_v0(){
      free(mask);
      free(data);
    }
    size_t count(const size_t& key) const {
      size_t n;
      return find_element(key,n);
      if (find_element(key,n)) return mask[n];
      return false;
    }
    size_t& operator[](const size_t& n) {
      size_t f;
      if(find_element(n,f)){
        if (VERBOSITY) cout << "found at " << f << endl; 
        return data[f][1];
      } else {
        if (VERBOSITY) cout << "did not find " << frac(n) << endl; 
      }
      if (grow(n_elem,n_data)) resize(grow(n_data));
      return data[insert_element(n)][1];
    }
    size_t size() const {
      return n_elem;
    }
    void erase(const size_t& n){
      if (VERBOSITY) cout << "erase " << frac(n) << endl;
      size_t f;
      if (!find_element(n,f)) return;
      if (VERBOSITY) cout << "found at " << f << endl;
      //if (f<map(n)){
        while (f>0){
          if (!mask[f-1]) break;
          if (map(data[f-1][0])>f-1){
            if (VERBOSITY) cerr << "swap down" << endl;
            swap(data[f][0],data[f-1][0]);
            swap(data[f][1],data[f-1][1]);
          } else {
            break;
          }
          --f;
        }
      //}
      //if (f>map(n)){
        while (f+1<n_data){
          if (!mask[f+1]) break;
          if (map(data[f+1][0])<f+1){
            if (VERBOSITY) cerr << "swap up" << endl;
            swap(data[f][0],data[f+1][0]);
            swap(data[f][1],data[f+1][1]);
          } else {
            break;
          }
          ++f;
        }
      //}
      if (data[f][0]!=n){
        cout << "erase failed for " << frac(n) << endl;
      }
      mask[f]=false;
      --n_elem;
    }
    void resize(const size_t& n){
      if (VERBOSITY) cerr << "resize" << endl;
      mask = (bool*)realloc(mask,n);
      for (size_t i=n_data;i<n;++i) mask[i] = false;
      data = (size_t(*)[2])realloc(data,n*2*sizeof(size_t));
      const size_t old_size = n_data;
      n_data = n;
      for (size_t i=old_size-1;i!=~size_t(0);--i){
        if (!mask[i]) continue;
        size_t tmp[2];
        tmp[0] = data[i][0];
        tmp[1] = data[i][1];
        --n_elem;
        mask[i]=false;
        //erase(tmp[0]);
        if (VERBOSITY) print();
        data[insert_element(tmp[0])][1]=tmp[1];
        if (VERBOSITY) print();
        if (VERBOSITY) cout << n_elem << endl;
        //(*this)[tmp[0]] = tmp[1];
      }
      if (VERBOSITY) cerr << "resize lives" << endl;
    }
    bool is_in_order(){
      size_t last_key = 0;
      for (size_t i=0;i!=n_data;++i){
        if (mask[i]) if (data[i][0]<last_key) {
          return false;
          last_key = data[i][0];
        }
      }
      for (size_t i=1;i<n_data;++i)
        if (mask[i])
          if (map(data[i][0])<i)
            if (!mask[i-1]) return false;
      for (size_t i=0;i+1<n_data;++i)
        if (mask[i])
          if (map(data[i][0])>i)
            if (!mask[i+1]) return false;
      return true;
    }
};

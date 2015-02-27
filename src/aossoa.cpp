
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include "simd.h"

#ifndef FP32
using real_t = double;
#else
using real_t = float;
#endif

/* Base class describing the AoS */
template<typename T>
struct DataBaseT
{
  public:
    T a, b, c, d, e;
    DataBaseT() {}
    DataBaseT(
        const T &_a,
        const T &_b,
        const T &_c,
        const T &_d,
        const T &_e) : a(_a), b(_b), c(_c), d(_d), e(_e) {}
    DataBaseT(const T &x) : a(x), b(x), c(x), d(x), e(x) {}

  protected:
    template<typename TT> void copy(const TT val)
    {
      a = val.a;
      b = val.b;
      c = val.c;
      d = val.d;
      e = val.e;
    }
};

template<typename DataT>
struct DataRefT;
template<typename DataT>
struct DataIRefT;

template<typename T, typename Tref, typename Tscal = T, typename Tiref = Tref>
struct DataT : public DataBaseT<T>
{
  using type = T;
  using ref_t  = Tref;
  using scal_t = Tscal;
  using iref_t = Tiref;

  DataT() : DataBaseT<T>() {}
  DataT(
      const T &a,
      const T &b,
      const T &c,
      const T &d,
      const T &e) : DataBaseT<T>(a,b,c,d,e) {}
  DataT(const T   &a) : DataBaseT<T>(a) {}
  friend DataT operator+(const DataT &x, const DataT &y) 
  {
    return DataT(
        x.a + y.a,
        x.b + y.b,
        x.c + y.c,
        x.d + y.d,
        x.e + y.e);

  }
  friend DataT operator-(const DataT &x, const DataT &y) 
  {
    return DataT(
        x.a - y.a,
        x.b - y.b,
        x.c - y.c,
        x.d - y.d,
        x.e - y.e);

  }
  friend DataT operator*(const DataT &x, const DataT &y) 
  {
    return DataT(
        x.a * y.a,
        x.b * y.b,
        x.c * y.c,
        x.d * y.d,
        x.e * y.e);

  }
  friend DataT operator/(const DataT &x, const DataT &y) 
  {
    return DataT(
        x.a / y.a,
        x.b / y.b,
        x.c / y.c,
        x.d / y.d,
        x.e / y.e);

  }
  void operator=(const DataIRefT<DataT> p) { this->copy(p); }
  void operator=(const DataRefT<DataT>  p) { this->copy(p); }
  void operator=(const DataT           &p) { this->copy(p); }
};

template<typename DataT, typename Tref>
struct DataRefBaseT : public DataBaseT<Tref>
{
  DataRefBaseT(
      Tref a,
      Tref b,
      Tref c,
      Tref d,
      Tref e) : DataBaseT<Tref>(a,b,c,d,e) {}
  operator DataT() const  
  {
    return DataT(
        this->a,
        this->b,
        this->c,
        this->d,
        this->e);
  }
};

template<typename DataT>
struct DataRefT : public DataRefBaseT<DataT, typename DataT::ref_t>
{
  using Tref  = typename DataT::ref_t;
  using Tscal = typename DataT::scal_t;
  DataRefT(
      Tscal &a,
      Tscal &b,
      Tscal &c,
      Tscal &d,
      Tscal &e) : DataRefBaseT<DataT,Tref>(a,b,c,d,e) {}
  void operator=(const DataIRefT<DataT> p) { this->copy(p); }
  void operator=(const DataRefT<DataT>  p) { this->copy(p); }
  void operator=(const DataT           &p) { this->copy(p); }
};
template<typename DataT>
struct DataIRefT : public DataRefBaseT<DataT, typename DataT::iref_t>
{
  using Tscal = typename DataT::scal_t;
  using Tiref = typename DataT::iref_t;
  DataIRefT(
      Tscal &a,
      Tscal &b,
      Tscal &c,
      Tscal &d,
      Tscal &e,
      const Simd<int> &idx) : 
    DataRefBaseT<DataT, Tiref>(
        Tiref(a,idx),
        Tiref(b,idx),
        Tiref(c,idx),
        Tiref(d,idx),
        Tiref(e,idx)) {}
  void operator=(const DataIRefT<DataT> p) { this->copy(p); }
  void operator=(const DataRefT<DataT>  p) { this->copy(p); }
  void operator=(const DataT           &p) { this->copy(p); }
};
 
using Data         = DataT<real_t, real_t&>;
using DataRef      = DataRefT<Data>;
using SimdData     = DataT    <Simd<real_t>, SimdRefT<real_t>, real_t, SimdIRefT<real_t>>;
using SimdDataRef  = DataRefT <SimdData>;
using SimdDataIRef = DataIRefT<SimdData>;

template<template<typename> class T1, template<typename> class T2>
static SimdData operator+(const T1<SimdData> x, const T2<SimdData> y) 
{
  return static_cast<SimdData>(x) + static_cast<SimdData>(y);
}
template<template<typename> class T1, template<typename> class T2>
static SimdData operator-(const T1<SimdData> x, const T2<SimdData> y) 
{
  return static_cast<SimdData>(x) - static_cast<SimdData>(y);
}
template<template<typename> class T1, template<typename> class T2>
static SimdData operator*(const T1<SimdData> x, const T2<SimdData> y) 
{
  return static_cast<SimdData>(x) * static_cast<SimdData>(y);
}
template<template<typename> class T1, template<typename> class T2>
static SimdData operator/(const T1<SimdData> x, const T2<SimdData> y) 
{
  return static_cast<SimdData>(x) / static_cast<SimdData>(y);
}

/////////////

struct DataSoA
{
  std::vector<real_t> a, b, c, d, e;

  DataSoA(const int n) 
  {
    a.resize(n);
    b.resize(n);
    c.resize(n);
    d.resize(n);
    e.resize(n);
  }
  DataRef operator[](const int i) { return DataRef(a[i],b[i],c[i],d[i],e[i]); }
  void push_back(const Data &p)
  {
    a.push_back(p.a);
    b.push_back(p.b);
    c.push_back(p.c);
    d.push_back(p.d);
    e.push_back(p.e);
  }
  size_t size() const { return a.size(); }
};

struct DataSoA_simd
{
  using vector = std::vector<real_t>;
  vector &a, &b, &c, &d, &e;
  DataSoA_simd(DataSoA &data) :
    a(data.a),
    b(data.b),
    c(data.c),
    d(data.d),
    e(data.e) {}
  SimdDataRef operator[](const int i){ return SimdDataRef(a[i],b[i],c[i],d[i],e[i]); }
  SimdDataIRef operator[](const Simd<int> &i) 
  {
    return SimdDataIRef(a[0],b[0],c[0],d[0],e[0],i);
  }
};

int main(int argc, char * argv[])
{
  using std::cout;
  using std::endl;
  using std::sqrt;

  const int n = argc > 1 ? atoi(argv[1]) : 34;

  cout << n << endl;
 
  cout << "\nfilling array \n"; 

  DataSoA array(n);
#pragma ivdep
  for (int i = 0; i < n; i++)
  {
#if 0
    auto&& p = array[i];
    p.a = i;
    p.b = i*i;
    p.c = sqrt(i);
    p.d = -1;
    p.e = 2*i+1;
#else
    Data p;
    p.a = i;
    p.b = i*i;
    p.c = sqrt(i);
    p.d = 123;
    p.e = 2*i+1;
    array[i] = p;
#endif
  }

  cout << "\nprinting array \n"; 
#pragma ivdep
  for (int i= 0; i < n; i++)
  {
    const auto&& pp = array[i];
    cout 
      << pp.a  << " "
      << pp.b  << " "
      << pp.c  << " "
      << pp.d  << " "
      << pp.e  << " "
      << endl;
  }

  ///////
  std::vector<Data> darray(n);
  std::vector<Data> dres(n);
  asm("#test1a");
#pragma ivdep
  for (int i = 1; i < n-1; i++)
  {
    Data d0 = darray[i];
    Data d1 = darray[i-1];
    Data d2 = darray[i+1];
#if 0
    d0.a = d0.a*d1.b + d2.c;
    d0.b = d0.b*d0.d + d1.c;
    res[i] = (d0/d1 + d0.b) * d2 + d0.a;
#else
    dres[i] = d0/d1 + d2;
#endif
  }

  asm("#test2a");

  DataSoA res(n);
  asm("#test1");
#pragma ivdep
  for (int i = 1; i < n-1; i++)
  {
    Data d0 = array[i];
    Data d1 = array[i-1];
    Data d2 = array[i+1];
#if 0
    d0.a = d0.a*d1.b + d2.c;
    d0.b = d0.b*d0.d + d1.c;
    res[i] = (d0/d1 + d0.b) * d2 + d0.a;
#else
    res[i] = d0/d1 + d2;
#endif
  }

  asm("#test2");
  cout << "\nprinting res \n"; 
  
  for (int i= 0; i < n; i++)
  {
    const auto&& pp = res[i];
    cout 
      << pp.a  << " "
      << pp.b  << " "
      << pp.c  << " "
      << pp.d  << " "
      << pp.e  << " "
      << endl;
  }

  asm("#test1");

  for (int i= 0; i < n; i++)
    res[i] = 0;

  DataSoA_simd res_vec(res);
  DataSoA_simd array_vec(array);
  asm("#testvec1");
  for (int i = 1; i < n-1; i += Simd<real_t>::VLEN)
  {
    SimdData d0 = array_vec[i];
    SimdData d1 = array_vec[i-1];
    SimdData d2 = array_vec[i+1];
#if 0
    array_vec[i].a = d0.a*d1.b + d2.c;
    d0.b = d0.b*d0.d + d1.c;
    res_vec[i] = (d0/d1 + d0.b) * d2 + d0.a;
#else
    res_vec[i] = d0*d1 + d2;
#endif
  }
  asm("#testvec2");
  cout << "\nprinting resvec \n"; 
  
  for (int i= 0; i < n; i++)
  {
    const auto&& pp = res[i];
    cout 
      << pp.a  << " "
      << pp.b  << " "
      << pp.c  << " "
      << pp.d  << " "
      << pp.e  << " "
      << endl;
  }
  
  //////
  
  cout << "VLEN= " << Simd<real_t>::VLEN << endl;
  cout << "sizeof(Data):     " << sizeof(Data)     << " = " << sizeof(Data)/sizeof(real_t) << "x real_t" << endl;
  cout << "sizeof(SimdData): " << sizeof(SimdData) << " = " << sizeof(SimdData)/sizeof(real_t) << "x real_t" << endl;

  std::vector<int> idx_src(n), idx_dst(n);
  std::iota(idx_src.begin(), idx_src.end(), 0);
  std::random_shuffle(idx_src.begin(), idx_src.end());
  std::iota(idx_dst.begin(), idx_dst.end(), 0);
  std::random_shuffle(idx_dst.begin(), idx_dst.end());

#if 1
  asm("#test3");
  auto& res1 = res;
#pragma ivdep
  for (int i = 0; i < n; i++)
  {
    const auto src = (idx_src[i]);
    const auto dst = (idx_dst[i]);
    auto p  = array[src] + dres[src];
    p.a = p.b;
    auto res = p.a*p.a + p.b*p.b + p.c*p.c;
    res1[dst].a = res;
    res1[dst].b = res*res;
    res1[dst].d = res*res-res;
    res1[dst].e = res/(res*res-res);
    Data d = array[i];
    d.a = res;
    d.b = res*res;
    d.c = res/res + res;
    d.d = res*res*res-res;
    d.e = res/(res*res - res);
    res1[dst] = array[src]; // + res1[dst] + d;
  };
  asm("#test4");
#endif

#if 1
  asm("#testvec3");
  for (int i = 0; i < n; i += Simd<real_t>::VLEN)
  {
    const auto src = SimdRefT<int>(idx_src[i]);
    const auto dst = SimdRefT<int>(idx_dst[i]);
    auto p  = array_vec[src];
    p.a = p.b;
    auto res = p.a*p.a + p.b*p.b + p.c*p.c;
    res_vec[dst].a = res;
    res_vec[dst].b = res*res;
    res_vec[dst].d = res*res-res;
    res_vec[dst].e = res/(res*res-res);
    SimdData d = array_vec[i];
    d.a = res;
    d.b = res*res;
    d.c = res/res + res;
    d.d = res*res*res-res;
    d.e = res/(res*res - res);
    res_vec[dst] = res_vec[src] + res_vec[dst];
  }
  asm("#testvec4");
#endif


  return 0;
}

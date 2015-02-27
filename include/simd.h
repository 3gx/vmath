#pragma once

#include <immintrin.h>

template<typename T>
struct Simd;


template<>
struct Simd<double>
{
  private:
    __m256d v1, v2;
    Simd(const __m256d &val1, const __m256d &val2) : v1(val1), v2(val2) {}
  public:
    enum {VLEN = 8};

    Simd() {}

    static Simd vloadu(const double *ptr)
    {
      return Simd(_mm256_loadu_pd(ptr), _mm256_loadu_pd(ptr+4));
    }
    static void vstoreu(double *ptr, const Simd &x)
    {
      _mm256_storeu_pd(ptr  , x.v1);
      _mm256_storeu_pd(ptr+4, x.v2);
    }
    inline static Simd vgather(const double *base, const Simd<int> &idx);
    inline static void vscatter(double *base, const Simd<int> &idx, const Simd &x);

    friend Simd operator+(const Simd &x, const Simd &y)
    {
      return Simd(x.v1 + y.v1, x.v2 + y.v2);
    }
    friend Simd operator-(const Simd &x, const Simd &y)
    {
      return Simd(x.v1 - y.v1, x.v2 - y.v2);
    }
    friend Simd operator*(const Simd &x, const Simd &y)
    {
      return Simd(x.v1 * y.v1, x.v2 * y.v2);
    }
    friend Simd operator/(const Simd &x, const Simd &y)
    {
      return Simd(x.v1 / y.v1, x.v2 / y.v2);
    }
};

template<>
struct Simd<float>
{
  private:
    __m256 v;
    Simd(const __m256 &value) : v(value) {}

  public:
    enum {VLEN = 8};
    Simd() {}

    inline static Simd vgather(const float *base, const Simd<int> &idx);
    inline static void vscatter(float *base, const Simd<int> &idx, const Simd &x);

    static Simd vloadu(const float *ptr)
    {
      return Simd(_mm256_loadu_ps(ptr));
    }
    static void vstoreu(float *ptr, const Simd &x)
    {
      _mm256_storeu_ps(ptr, x.v);
    }

    friend Simd operator+(const Simd &x, const Simd &y)
    {
      return Simd(x.v + y.v);
    }
    friend Simd operator-(const Simd &x, const Simd &y)
    {
      return Simd(x.v - y.v);
    }
    friend Simd operator*(const Simd &x, const Simd &y)
    {
      return Simd(x.v * y.v);
    }
    friend Simd operator/(const Simd &x, const Simd &y)
    {
      return Simd(x.v / y.v);
    }
};

template<>
struct Simd<int>
{
  private:
    union
    {
      __m256i v;
      struct {__m128i v1, v2;};
    };
    Simd(const __m256i &value) : v(value) {}

  public:
    enum {VLEN = 8};
    __m256i getFull() const { return v; }
    __m128i getFirst() const {return v1;}
    __m128i getSecond() const {return v2;}
    
    Simd() {}
    static Simd vloadu(const int *ptr)
    {
      return Simd(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(ptr)));
    }
    static void vstoreu(int *ptr, const Simd &x)
    {
      _mm256_storeu_si256(reinterpret_cast<__m256i*>(ptr), x.v);
    }
    

  private:
    template<typename Op>
      static Simd compute(const Simd &x, const Simd &y, Op&& op)
      {
        int sx[VLEN], sy[VLEN], sz[VLEN];
        vstoreu(sx,x);
        vstoreu(sy,y);
        for (int i = 0; i < VLEN; i++)
          sz[i] = op(sx[i], sy[i]);
        return vloadu(sz);
      }

  public:


    friend Simd operator+(const Simd &x, const Simd &y)
    {
#ifdef __clang__
      return compute(x, y, [](int x, int y) { return x + y; });
#else
      return Simd(x.v + y.v);
#endif
    }
    friend Simd operator-(const Simd &x, const Simd &y)
    {
#ifdef __clang__
      return compute(x, y, [](int x, int y) { return x - y; });
#else
      return Simd(x.v - y.v);
#endif
    }
    friend Simd operator*(const Simd &x, const Simd &y)
    {
#ifdef __clang__
      return compute(x, y, [](int x, int y) { return x * y; });
#else
      return Simd(x.v * y.v);
#endif
    }
    friend Simd operator/(const Simd &x, const Simd &y)
    {
#ifdef __clang__
      return compute(x, y, [](int x, int y) { return x / y; });
#else
      return Simd(x.v / y.v);
#endif
    }
};
   
inline Simd<double> Simd<double>::vgather(const double *base, const Simd<int> &idx)
{
  return Simd<double>(
      _mm256_i32gather_pd(base, idx.getFirst(), /*scale*/ 8),
      _mm256_i32gather_pd(base, idx.getSecond(), /*scale*/ 8)
      );
}
inline void Simd<double>::vscatter(double *base, const Simd<int> &idx, const Simd<double> &x)
{
  auto index = reinterpret_cast<const    int*>(&idx);
  auto value = reinterpret_cast<const double*>(&  x);
  for (int i= 0; i < VLEN; i++)
    base[index[i]] = value[i];
}
inline Simd<float> Simd<float>::vgather(const float *base, const Simd<int> &idx)
{
  return Simd<float>(
      _mm256_i32gather_ps(base, idx.getFull(), /*scale*/ 4)
      );
}
inline void Simd<float>::vscatter(float *base, const Simd<int> &idx, const Simd<float> &x)
{
  auto index = reinterpret_cast<const   int*>(&idx);
  auto value = reinterpret_cast<const float*>(&  x);
  for (int i= 0; i < VLEN; i++)
    base[index[i]] = value[i];
}

/////////////

template<typename T>
struct SimdRefT
{
  private:
    T &ref;
  public:
    SimdRefT(T& v) : ref(v) {}
    operator Simd<T>() const 
    {
      return Simd<T>::vloadu(&ref);
    }
    SimdRefT operator=(const Simd<T> &value)
    {
      Simd<T>::vstoreu(&ref,value);
      return *this;
    }
    SimdRefT operator=(const SimdRefT &value)
    {
      Simd<T>::vstoreu(&ref,static_cast<Simd<T>>(value));
      return *this;
    }
};

template<typename T>
struct SimdIRefT
{
  private:
    T &base;
    const Simd<int> &idx;
  public:
    SimdIRefT(T &_base, const Simd<int> &_idx) : base(_base), idx(_idx) {}
    operator Simd<T>() const
    {
      return Simd<T>::vgather(&base,idx);
    }
    SimdIRefT operator=(const Simd<T> &value)
    {
      Simd<T>::vscatter(&base, idx, value);
      return *this;
    }
    SimdIRefT operator=(const SimdIRefT &value)
    {
      Simd<T>::vscatter(&base, idx, static_cast<Simd<T>>(value));
      return *this;
    }
};

template<template<typename> class Tref1, template<typename> class Tref2, typename T>
static inline Simd<T> operator*(const Tref1<T> x, const Tref2<T> y)
{
  return static_cast<Simd<T>>(x) * static_cast<Simd<T>>(y);
}

template<template<typename> class Tref1, template<typename> class Tref2, typename T>
static inline Simd<T> operator+(const Tref1<T> x, const Tref2<T> y)
{
  return static_cast<Simd<T>>(x) + static_cast<Simd<T>>(y);
}

template<template<typename> class Tref1, template<typename> class Tref2, typename T>
static inline Simd<T> operator-(const Tref1<T> x, const Tref2<T> y)
{
  return static_cast<Simd<T>>(x) - static_cast<Simd<T>>(y);
}

template<template<typename> class Tref1, template<typename> class Tref2, typename T>
static inline Simd<T> operator/(const Tref1<T> x, const Tref2<T> y)
{
  return static_cast<Simd<T>>(x) / static_cast<Simd<T>>(y);
}


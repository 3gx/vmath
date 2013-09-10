#pragma once

#ifndef __AVX__
#error "Please enable AVX support"
#endif

#define __out

#include <cassert>

namespace SIMD
{
  enum {WIDTH = 8};
  typedef float _v8sf  __attribute__((vector_size(32)));

  struct fvec 
  {
    _v8sf val;
    fvec() {}
    fvec(_v8sf _val) : val(_val) {};
    fvec(const float  f) : val((_v8sf){f,f,f,f,f,f,f,f}) {}
    fvec(const float  a, const float b, const float c, const float d, const float e, const float f, const float g, const float h) : 
      val((_v8sf){a,b,c,d,e,f,g,h}) {}
    fvec(const float *p, const bool aligned = true) : 
      val(aligned ? *(_v8sf*)p : __builtin_ia32_loadups256(p)) {}

    operator _v8sf() const {return val;}

    const fvec& operator =(const fvec a) {val  = a.val; return *this;}
    const fvec& operator+=(const fvec a) {val += a.val; return *this;}
    const fvec& operator-=(const fvec a) {val -= a.val; return *this;}
    const fvec& operator*=(const fvec a) {val *= a.val; return *this;}
    const fvec& operator/=(const fvec a) {val /= a.val; return *this;}
    const fvec operator -()  const {return fvec(-val);}
    const fvec operator +(const fvec a) const {return fvec(val + a.val);}
    const fvec operator -(const fvec a) const {return fvec(val - a.val);}
    const fvec operator *(const fvec a) const {return fvec(val * a.val);}
    const fvec operator /(const fvec a) const {return fvec(val / a.val);}

    operator bool() const {return true;}
    const fvec operator!=(const fvec a) const 
    {
      return fvec(__builtin_ia32_cmpps256(val, a.val, 28));
    }
#if 0
    const fvec operator==(const fvec a) const {
      return fvec(__builtin_ia32_cmpeqps(val, a.val));
    }

    const fvec operator<(const fvec &rhs) const{
      return __builtin_ia32_cmpltps(val, rhs.val);
    }
    const fvec operator<=(const fvec &rhs) const{
      return __builtin_ia32_cmpleps(val, rhs.val);
    }
    const fvec operator>(const fvec &rhs) const{
      return __builtin_ia32_cmpgtps(val, rhs.val);
    }
    const fvec operator>=(const fvec &rhs) const{
      return __builtin_ia32_cmpgeps(val, rhs.val);
    }
    const fvec operator|(const fvec &rhs) const{
      return __builtin_ia32_orps(val, rhs.val);
    }
    const fvec operator&(const fvec &rhs) const{
      return __builtin_ia32_andps(val, rhs.val);
    }
#endif


    const float operator[] (const int i) const {
      union {_v8sf v; float s[8];} test;
      test.v = val;
      return test.s[i];
    }


  };

  inline static void transpose(
      const fvec v[8], fvec t[8])
  {
#define _MM_SHUFFLE(z,y,x,w) ((z<<6) | (y<<4) | (x<<2) | w)
    const _v8sf __t0 = __builtin_ia32_unpcklps256(v[0], v[1]);
    const _v8sf __t1 = __builtin_ia32_unpckhps256(v[0], v[1]);
    const _v8sf __t2 = __builtin_ia32_unpcklps256(v[2], v[3]);
    const _v8sf __t3 = __builtin_ia32_unpckhps256(v[2], v[3]);
    const _v8sf __t4 = __builtin_ia32_unpcklps256(v[4], v[5]);
    const _v8sf __t5 = __builtin_ia32_unpckhps256(v[4], v[5]);
    const _v8sf __t6 = __builtin_ia32_unpcklps256(v[6], v[7]);
    const _v8sf __t7 = __builtin_ia32_unpckhps256(v[6], v[7]);
    const _v8sf __tt0 = __builtin_ia32_shufps256(__t0,__t2, _MM_SHUFFLE(1,0,1,0));
    const _v8sf __tt1 = __builtin_ia32_shufps256(__t0,__t2, _MM_SHUFFLE(3,2,3,2));
    const _v8sf __tt2 = __builtin_ia32_shufps256(__t1,__t3, _MM_SHUFFLE(1,0,1,0));
    const _v8sf __tt3 = __builtin_ia32_shufps256(__t1,__t3, _MM_SHUFFLE(3,2,3,2));
    const _v8sf __tt4 = __builtin_ia32_shufps256(__t4,__t6, _MM_SHUFFLE(1,0,1,0));
    const _v8sf __tt5 = __builtin_ia32_shufps256(__t4,__t6, _MM_SHUFFLE(3,2,3,2));
    const _v8sf __tt6 = __builtin_ia32_shufps256(__t5,__t7, _MM_SHUFFLE(1,0,1,0));
    const _v8sf __tt7 = __builtin_ia32_shufps256(__t5,__t7, _MM_SHUFFLE(3,2,3,2));
    t[0] = __builtin_ia32_vperm2f128_ps256(__tt0, __tt4, 0x20);
    t[1] = __builtin_ia32_vperm2f128_ps256(__tt1, __tt5, 0x20);
    t[2] = __builtin_ia32_vperm2f128_ps256(__tt2, __tt6, 0x20);
    t[3] = __builtin_ia32_vperm2f128_ps256(__tt3, __tt7, 0x20);
    t[4] = __builtin_ia32_vperm2f128_ps256(__tt0, __tt4, 0x31);
    t[5] = __builtin_ia32_vperm2f128_ps256(__tt1, __tt5, 0x31);
    t[6] = __builtin_ia32_vperm2f128_ps256(__tt2, __tt6, 0x31);
    t[7] = __builtin_ia32_vperm2f128_ps256(__tt3, __tt7, 0x31);
#undef _MM_SHUFFLE
  }

  template<const int N>
    inline static void transpose(fvec *d)
    {
      if (N == 1)
      {
        transpose(d, d);
      }
      else if (N == 2)
      {
        fvec m11[8] = {d[ 0], d[ 2], d[ 4], d[ 6], d[ 8], d[10], d[12], d[14]};
        fvec m12[8] = {d[ 1], d[ 3], d[ 5], d[ 7], d[ 9], d[11], d[13], d[15]};
        fvec m21[8] = {d[16], d[18], d[20], d[22], d[24], d[26], d[28], d[30]};
        fvec m22[8] = {d[17], d[19], d[21], d[23], d[25], d[27], d[29], d[31]};
        transpose(m11, m11);
        transpose(m12, m12);
        transpose(m21, m21);
        transpose(m22, m22);
        for (int ch = 0; ch < 16; ch += 2)
        {
          d[   ch  ] = m11[ch>>1];
          d[   ch+1] = m21[ch>>1];
          d[16+ch  ] = m12[ch>>1];
          d[16+ch+1] = m22[ch>>1];
        }
      }
      else
        assert(N == 1);
    }

  template<class T>
    inline static void sstore(T &ptr, const fvec val)
    {
      __builtin_ia32_movntps256((float*)&ptr, val);
    }

  template<const bool ALIGN, typename T>
    inline static fvec store(T &ptr, const fvec val)
    {
      if (ALIGN)
      {
        *(_v8sf*)&ptr = val;
      }
      else
      {
        __builtin_ia32_storeups256((float*)&ptr, val);
      }
      return val;
    }
  
  template<const bool ALIGN, typename T>
    inline static fvec load(const T &ptr)
    {
      return fvec((float*)&ptr, ALIGN);
    }

  template<const int N>
    inline static fvec broadcast(const fvec x)
    {
      const int NN = N & 3;
      const int mask = 
        NN == 0 ? 0 :
        NN == 1 ? 0x55 :
        NN == 2 ? 0xAA  : 0xFF;

      const _v8sf y = __builtin_ia32_shufps256(x, x, mask);

      return N < 4 ? 
        __builtin_ia32_vperm2f128_ps256(y, y, 0x20) : 
        __builtin_ia32_vperm2f128_ps256(y, y, 0x31);
    }

  typedef fvec real;

}
#define LDA(x)   (SIMD::load<true >(x))
#define LDU(x)   (SIMD::load<false>(x))
#define STU(x,y) (SIMD::store<false>(x, y))
#define STA(x)   (*(SIMD::real*)&x)


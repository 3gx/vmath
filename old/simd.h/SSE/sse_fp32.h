#pragma once

struct fvec 
{
  enum {WIDTH = 4};
  _v4sf val;
  fvec() {}
  fvec(_v4sf _val) : val(_val) {};
  fvec(const float  f) : val((_v4sf){f,f,f,f}) {}
  fvec(const float  a, const float b, const float c, const float d) : val((_v4sf){a,b,c,d}) {}
  fvec(const float *p, const bool aligned = true) : 
    val(aligned ? *(_v4sf*)p : __builtin_ia32_loadups(p)) {}

  operator _v4sf() const {return val;}

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

  const fvec operator!=(const fvec a) const 
  {
    return fvec(__builtin_ia32_cmpneqps(val, a.val));
  }
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

  operator bool() const {return true;}

  const float operator[] (const int i) const {
    switch(i) {
      case 0:	return __builtin_ia32_vec_ext_v4sf(val, 0);
      case 1:	return __builtin_ia32_vec_ext_v4sf(val, 1);
      default: assert(0);
    }
  }


};

inline static void transpose(
    const fvec v0, const fvec v1, const fvec v2, const fvec v3,
    fvec &t0, fvec &t1, fvec &t2, fvec &t3)
{
  t0 = __builtin_ia32_unpcklps(v0, v1);
  t1 = __builtin_ia32_unpckhps(v0, v1);

  const _v4sf a0 = __builtin_ia32_unpcklps(v0, v2);
  const _v4sf a1 = __builtin_ia32_unpckhps(v0, v2);
  const _v4sf a2 = __builtin_ia32_unpcklps(v1, v3);
  const _v4sf a3 = __builtin_ia32_unpckhps(v1, v3);

  t0 = __builtin_ia32_unpcklps(a0, a2);
  t1 = __builtin_ia32_unpckhps(a0, a2);
  t2 = __builtin_ia32_unpcklps(a1, a3);
  t3 = __builtin_ia32_unpckhps(a1, a3);
}

  template<const int N>
inline static void transpose(fvec *d)
{
  if (N == 1)
    transpose(d[0], d[1], d[2], d[3], d[0], d[1], d[2], d[3]);
  else if (N==2)
  {
    const fvec m11[4] = {d[0], d[ 2], d[ 4], d[ 6]};
    const fvec m12[4] = {d[1], d[ 3], d[ 5], d[ 7]};
    const fvec m21[4] = {d[8], d[10], d[12], d[14]};
    const fvec m22[4] = {d[9], d[11], d[13], d[15]};
    transpose(m11[0], m11[1], m11[2], m11[3], d[0], d[ 2], d[ 4], d[ 6]);
    transpose(m21[0], m21[1], m21[2], m21[3], d[1], d[ 3], d[ 5], d[ 7]);
    transpose(m12[0], m12[1], m12[2], m12[3], d[8], d[10], d[12], d[14]);
    transpose(m22[0], m22[1], m22[2], m22[3], d[9], d[11], d[13], d[15]);
  }
  else
    assert(N==1);
}

inline static void sstore(float &ptr, const fvec val)
{
  __builtin_ia32_movntps(&ptr, val);
}

  template<const bool ALIGN>
inline static fvec store(float &ptr, const fvec val)
{
  if (ALIGN)
  {
    *(_v4sf*)&ptr = val;
  }
  else
  {
    __builtin_ia32_storeups((float*)&ptr, val);
  }
  return val;
}

template<const bool ALIGN>
inline static fvec load(const float &ptr)
{
  return fvec(&ptr, ALIGN);
}

  template<const int N>
inline static fvec broadcast(const fvec x)
{
  return
    N == 0 ? __builtin_ia32_shufps(x,x, 0x00) : 
    N == 1 ? __builtin_ia32_shufps(x,x, 0x55) :
    N == 2 ? __builtin_ia32_shufps(x,x, 0xAA) :
    __builtin_ia32_shufps(x,x, 0xFF);
}


inline static fvec sqrt(const fvec x)
{
  return __builtin_ia32_sqrtps(x);
}

inline static fvec rsqrt(const fvec x)
{
#if 1
  const fvec y = __builtin_ia32_rsqrtps(x);
  return (fvec(-0.5f) * y) * (x*y*y + fvec(-3.0f));
#else
  return __builtin_ia32_rsqrtps(x);
#endif
}



#pragma once

struct dvec 
{
  enum {WIDTH = 2};
  _v2df val;
  dvec() {}
  dvec(_v2df _val) : val(_val) {};
  dvec(const double  f) : val((_v2df){f,f}) {}
  dvec(const double  a, const double b) : val((_v2df){a,b}) {}
  dvec(const double *p, const bool aligned = true) : 
    val(aligned ? *(_v2df*)p : __builtin_ia32_loadupd(p)) {}

  operator _v2df() const {return val;}

  const dvec& operator =(const dvec a) {val  = a.val; return *this;}
  const dvec& operator+=(const dvec a) {val += a.val; return *this;}
  const dvec& operator-=(const dvec a) {val -= a.val; return *this;}
  const dvec& operator*=(const dvec a) {val *= a.val; return *this;}
  const dvec& operator/=(const dvec a) {val /= a.val; return *this;}
  const dvec operator -()  const {return dvec(-val);}
  const dvec operator +(const dvec a) const {return dvec(val + a.val);}
  const dvec operator -(const dvec a) const {return dvec(val - a.val);}
  const dvec operator *(const dvec a) const {return dvec(val * a.val);}
  const dvec operator /(const dvec a) const {return dvec(val / a.val);}

  const dvec operator!=(const dvec a) const 
  {
    return dvec(__builtin_ia32_cmpneqpd(val, a.val));
  }
  const dvec operator==(const dvec a) const {
    return dvec(__builtin_ia32_cmpeqpd(val, a.val));
  }

  const dvec operator<(const dvec &rhs) const{
    return __builtin_ia32_cmpltpd(val, rhs.val);
  }
  const dvec operator<=(const dvec &rhs) const{
    return __builtin_ia32_cmplepd(val, rhs.val);
  }
  const dvec operator>(const dvec &rhs) const{
    return __builtin_ia32_cmpgtpd(val, rhs.val);
  }
  const dvec operator>=(const dvec &rhs) const{
    return __builtin_ia32_cmpgepd(val, rhs.val);
  }
  const dvec operator|(const dvec &rhs) const{
    return __builtin_ia32_orpd(val, rhs.val);
  }
  const dvec operator&(const dvec &rhs) const{
    return __builtin_ia32_andpd(val, rhs.val);
  }

  operator bool() const {return true;}

  const double operator[] (const int i) const {
    switch(i) {
      case 0:	return __builtin_ia32_vec_ext_v2df(val, 0);
      case 1:	return __builtin_ia32_vec_ext_v2df(val, 1);
      default: assert(0);
    }
  }


};

inline static void transpose(const dvec v0, const dvec v1, dvec &t0, dvec &t1)
{
  t0 = __builtin_ia32_unpcklpd(v0, v1);
  t1 = __builtin_ia32_unpckhpd(v0, v1);
}

  template<const int N>
inline static void transpose(dvec *d)
{
  dvec v1, v2;
  switch(N)
  {
    case 1:
      transpose(d[0], d[1], d[0], d[1]);
      break;
    case 2:
      transpose(d[0], d[2], d[0], d[2]);
      transpose(d[5], d[7], d[5], d[7]);
      transpose(d[1], d[3], v1, v2);
      transpose(d[4], d[6], d[1], d[3]);
      d[4] = v1;
      d[6] = v2;
      break;
    default:
      assert(N == 1);
  }
}


inline static void sstore(double &ptr, const dvec val)
{
  __builtin_ia32_movntpd(&ptr, val);
}

  template<const bool ALIGN>
inline static dvec store(double &ptr, const dvec val)
{
  if (ALIGN)
  {
    *(_v2df*)&ptr = val;
  }
  else
  {
    __builtin_ia32_storeupd(&ptr, val);
  }
  return val;
}

  template<const bool ALIGN>
inline static dvec load(const double &ptr)
{
  return dvec(&ptr, ALIGN);
}

  template<const int N>
inline static dvec broadcast(const dvec x)
{
  const int mask =
    N == 0 ? (0 << 0) + (0 << 0) : (1 << 0) + (1 << 1);
  return __builtin_ia32_shufpd(x, x, mask);
}

inline static dvec sqrt(const dvec x)
{
  return __builtin_ia32_sqrtpd(x);
}

inline static dvec rsqrt(const dvec x)
{
#if 1
  const _v2df y = __builtin_ia32_cvtps2pd(__builtin_ia32_rsqrtps(__builtin_ia32_cvtpd2ps(x)));
  return (dvec(-0.5) * y) * (x*y*y + dvec(-3.0));
#else
  const _v4sf y = __builtin_ia32_cvtpd2ps(x);
  const _v4sf z = __builtin_ia32_rsqrtps(y);
  return __builtin_ia32_cvtps2pd(z);
#endif
}




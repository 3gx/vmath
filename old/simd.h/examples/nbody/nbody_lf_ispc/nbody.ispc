#include "defs.h"

typedef real<3> vec3;
typedef real<4> vec4;

struct Predictor
{
  real mass;
  vec3 pos;
};

struct Force
{
  vec3 acc;
  real pot;
};

struct real4
{
  real x, y, z, w;
};

static inline 
real dot(const vec3 &a, const vec3 &b) {return a.x*b.x + a.y*b.y + a.z*b.z;}

static inline
real RSQRT(const real x) {
#if 0
  return 1.0/sqrt(x);
#else
  return rsqrt(x);
#endif
}

static inline 
void body_body_force(
    Force &f,
    const Predictor &pi, 
    const Predictor &pj, 
    const real eps2)
{
  const vec3 dr  = pj.pos - pi.pos;
  const real ds2 = dot(dr,dr) + eps2;

  const real  inv_ds  = RSQRT(ds2);
  const real  inv_ds2 = inv_ds*inv_ds;
  const real minv_ds  = inv_ds  * pj.mass;
  const real minv_ds3 = inv_ds2 * minv_ds;

  f.acc += minv_ds3 * dr;
  f.pot -= minv_ds;
}

export void compute_forces(
    uniform const int     n,
    uniform const real mass[],
    uniform const real posx[],
    uniform const real posy[],
    uniform const real posz[],
    uniform       real accx[],
    uniform       real accy[],
    uniform       real accz[],
    uniform       real gpot[],
    uniform const real eps2)
{
  foreach (i = 0 ... n)
  {
    Force fi;
    Predictor pi;
    pi.mass  = mass[i];
    pi.pos.x = posx[i];
    pi.pos.y = posy[i];
    pi.pos.z = posz[i];
    fi.acc.x = 0.0;
    fi.acc.y = 0.0;
    fi.acc.z = 0.0;
    fi.pot   = 0.0;

    for (uniform int j = 0; j < n; j++)
    {
      Predictor pj;
      pj.mass  = mass[j];
      pj.pos.x = posx[j];
      pj.pos.y = posy[j];
      pj.pos.z = posz[j];
      body_body_force(fi, pi, pj, eps2);
    }
    accx[i] = fi.acc.x;
    accy[i] = fi.acc.y;
    accz[i] = fi.acc.z;
    gpot[i] = fi.pot;
  }
}


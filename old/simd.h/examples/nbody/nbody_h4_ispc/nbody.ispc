#include "defs.h"

typedef real<3> vec3;
typedef real<4> vec4;


struct Predictor
{
  real mass;
  vec3 pos;
  vec3 vel;
  real iPad;
};

struct Force
{
  vec3 acc;
  vec3 jrk;
  real pot;
  real iPad;
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
  
  const vec3 dv = pj.vel - pi.vel;
  const real rv = dot(dr,dv);

  const real Jij = (-3.0) * rv * inv_ds2 * minv_ds3;

  f.jrk += minv_ds3*dv + Jij*dr;

}

export void compute_forces(
    uniform const int     n,
    uniform const real mass[],
    uniform const real posx[],
    uniform const real posy[],
    uniform const real posz[],
    uniform const real velx[],
    uniform const real vely[],
    uniform const real velz[],
    uniform       real accx[],
    uniform       real accy[],
    uniform       real accz[],
    uniform       real jrkx[],
    uniform       real jrky[],
    uniform       real jrkz[],
    uniform       real gpot[],
    uniform const real eps2)
{
  foreach (i = 0 ... n)
  {
    const real imass   = mass[i];
    const vec3 ipos    = {posx[i], posy[i], posz[i]};
    const vec3 ivel    = {velx[i], vely[i], velz[i]};
    const Predictor pi = {imass, ipos, ivel};
    Force fi;
    const vec3 null = {0.0, 0.0, 0.0};
    fi.acc = fi.jrk = null;
    fi.pot = 0.0;
    for (uniform int j = 0; j < n; j++)
    {
      const real jmass   = mass[j];
      const vec3 jpos    = {posx[j], posy[j], posz[j]};
      const vec3 jvel    = {velx[j], vely[j], velz[j]};
      const Predictor pj = {jmass, jpos, jvel};
      body_body_force(fi, pi, pj, eps2);
    }
    accx[i] = fi.acc.x;
    accy[i] = fi.acc.y;
    accz[i] = fi.acc.z;
    gpot[i] = fi.pot;

    jrkx[i] = fi.jrk.x;
    jrky[i] = fi.jrk.y;
    jrkz[i] = fi.jrk.z;
  }
}


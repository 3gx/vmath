#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <vector>
#include <cassert>

#if 0
#include "SSE/sse.h"
#else
#include "AVX/avx.h"
#endif

#if 1
typedef float       real;
typedef SIMD::fvec vreal;
#else
typedef double real;
typedef SIMD::dvec vreal;
#endif


struct Force
{
  vreal acc[3];
  vreal jrk[3];
  vreal pot;
};

struct Predictor
{
  vreal pos[3];
  vreal vel[3];
  vreal mass;
};

static inline 
vreal dot(const vreal a[3], const vreal b[3]) 
{
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static inline
vreal RSQRT(const vreal x) 
{
#if 0
  return vreal(1.0)/SIMD::sqrt(x);
#else
  return SIMD::rsqrt(x);
#endif
}

static inline 
void body_body_force(
    Force &f,
    const Predictor &pi, 
    const Predictor &pj, 
    const vreal eps2)
{
  const vreal dr[3] = {
    pj.pos[0] - pi.pos[0],
    pj.pos[1] - pi.pos[1],
    pj.pos[2] - pi.pos[2] };
  const vreal ds2 = dot(dr,dr) + eps2;

  const vreal  inv_ds  = RSQRT(ds2);
  const vreal  inv_ds2 = inv_ds*inv_ds;
  const vreal minv_ds  = inv_ds  * pj.mass;
  const vreal minv_ds3 = inv_ds2 * minv_ds;

  f.acc[0] += minv_ds3 * dr[0];
  f.acc[1] += minv_ds3 * dr[1];
  f.acc[2] += minv_ds3 * dr[2];
  f.pot -= minv_ds;
  
  const vreal dv[3] = {
    pj.vel[0] - pi.vel[0],
    pj.vel[1] - pi.vel[1],
    pj.vel[2] - pi.vel[2] };
  const vreal rv = dot(dr,dv);
  
  const vreal Jij = vreal(-3.0) * rv * inv_ds2 * minv_ds3;
  
  f.jrk[0] += minv_ds3*dv[0] + Jij*dr[0];
  f.jrk[1] += minv_ds3*dv[1] + Jij*dr[1];
  f.jrk[2] += minv_ds3*dv[2] + Jij*dr[2];
}



void compute_forces(
    const int     n,
    const real mass[],
    const real posx[],
    const real posy[],
    const real posz[],
    const real velx[],
    const real vely[],
    const real velz[],
    real accx[],
    real accy[],
    real accz[],
    real jrkx[],
    real jrky[],
    real jrkz[],
    real gpot[],
    const real seps2)
{
  const vreal eps2(seps2);
  for (int i = 0; i < n; i+= vreal::WIDTH)
  {
    Force fi;
    Predictor pi;
    pi.mass   = LDA(mass[i]);
    pi.pos[0] = LDA(posx[i]);
    pi.pos[1] = LDA(posy[i]);
    pi.pos[2] = LDA(posz[i]);
    pi.vel[0] = LDA(velx[i]);
    pi.vel[1] = LDA(vely[i]);
    pi.vel[2] = LDA(velz[i]);
    fi.acc[0] = 0.0;
    fi.acc[1] = 0.0;
    fi.acc[2] = 0.0;
    fi.jrk[0] = 0.0;
    fi.jrk[1] = 0.0;
    fi.jrk[2] = 0.0;
    fi.pot   = 0.0;

    for (int j = 0; j < n; j++)
    {
      Predictor pj;
      pj.mass   = mass[j];
      pj.pos[0] = posx[j];
      pj.pos[1] = posy[j];
      pj.pos[2] = posz[j];
      pj.vel[0] = velx[j];
      pj.vel[1] = vely[j];
      pj.vel[2] = velz[j];
      body_body_force(fi, pi, pj, eps2);
    }
    STA(accx[i]) = fi.acc[0];
    STA(accy[i]) = fi.acc[1];
    STA(accz[i]) = fi.acc[2];
    STA(jrkx[i]) = fi.jrk[0];
    STA(jrky[i]) = fi.jrk[1];
    STA(jrkz[i]) = fi.jrk[2];
    STA(gpot[i]) = fi.pot;
  }
}

#define FTINY 1e-10
#define FHUGE 1e+10

#include <sys/time.h>
inline double get_time() {
  struct timeval Tvalue;
  struct timezone dummy;

  gettimeofday(&Tvalue,&dummy);
  return ((double) Tvalue.tv_sec +1.e-6*((double) Tvalue.tv_usec));
}


struct real4 {
  real x, y, z, w;
  real4() {};
  real4(const real r) : x(r), y(r), z(r), w(r) {}
  real4(const real &_x, 
      const real &_y,
      const real &_z,
      const real &_w = 0) : x(_x), y(_y), z(_z), w(_w) {};
  const real abs2() const {return x*x + y*y + z*z;}
};

inline const real sqr(const real &x) {return x*x;}


#define BLOCKDIM 64

void forces(
    const int n,
    const real4 pos[],
    const real4 vel[],
    real4 acc[],
    real4 jrk[],
    const real eps2) 
{
  const int NALLOC = 8192+8;
  assert(n <= NALLOC);
  static real mass[NALLOC] __attribute__((aligned(64)));
  static real posx[NALLOC] __attribute__((aligned(64)));
  static real posy[NALLOC] __attribute__((aligned(64)));
  static real posz[NALLOC] __attribute__((aligned(64)));
  static real velx[NALLOC] __attribute__((aligned(64)));
  static real vely[NALLOC] __attribute__((aligned(64)));
  static real velz[NALLOC] __attribute__((aligned(64)));
  static real accx[NALLOC] __attribute__((aligned(64)));
  static real accy[NALLOC] __attribute__((aligned(64)));
  static real accz[NALLOC] __attribute__((aligned(64)));
  static real jrkx[NALLOC] __attribute__((aligned(64)));
  static real jrky[NALLOC] __attribute__((aligned(64)));
  static real jrkz[NALLOC] __attribute__((aligned(64)));
  static real gpot[NALLOC] __attribute__((aligned(64)));
  for (int i =0; i < n; i++)
  {
    mass[i] = pos[i].w;
    posx[i] = pos[i].x;
    posy[i] = pos[i].y;
    posz[i] = pos[i].z;
    velx[i] = vel[i].x;
    vely[i] = vel[i].y;
    velz[i] = vel[i].z;
  }
  compute_forces(
      n,
      &mass[0],
      &posx[0],
      &posy[0],
      &posz[0],
      &velx[0],
      &vely[0],
      &velz[0],
      &accx[0],
      &accy[0],
      &accz[0],
      &jrkx[0],
      &jrky[0],
      &jrkz[0],
      &gpot[0],
      eps2);
  for (int i = 0; i < n; i++)
  {
    acc[i].w = gpot[i];
    acc[i].x = accx[i];
    acc[i].y = accy[i];
    acc[i].z = accz[i];
    jrk[i].x = jrkx[i];
    jrk[i].y = jrky[i];
    jrk[i].z = jrkz[i];
  }
}

const real timestep(
    const int n,
    const real4 vel[],
    const real4 acc[],
    const real eta,
    const real eps2) 
{
  real dt_min = FHUGE;
  for (int i = 0; i < n; i++) {
    const real vabs = std::sqrt(vel[i].abs2());
    const real aabs = std::sqrt(acc[i].abs2());

    const real s = std::sqrt(eps2);
    const real idt1 = s/(vabs + FTINY);
    const real idt2 = std::sqrt(s/(aabs + FTINY));
    dt_min = std::min(dt_min, std::min(idt1, idt2));
  }
  dt_min *= eta;

  return dt_min;
}

void energy(
    const int n,
    const real4 pos[],
    const real4 vel[],
    const real4 acc[],
    real &Ekin, real &Epot) {
  Ekin = Epot = 0;
  for (int i = 0; i < n; i++) {
    Ekin += pos[i].w * vel[i].abs2() * 0.5;
    Epot += 0.5*pos[i].w * acc[i].w;
  }
}

const real iterate(
    const int n,
    real4 pos[],
    real4 vel[],
    real4 acc[],
    real4 jrk[],
    const real eta,
    const real eps2,
    const real dt) 
{
  std::vector<real4> acc0(n), jrk0(n);
  const real dt2 = dt*(1.0/2.0);
  const real dt3 = dt*(1.0/3.0);

  for (int i = 0; i < n; i++)
  {
    acc0[i] = acc[i];
    jrk0[i] = jrk[i];

    pos[i].x += dt*(vel[i].x + dt2*(acc[i].x + dt3*jrk[i].x));
    pos[i].y += dt*(vel[i].y + dt2*(acc[i].y + dt3*jrk[i].y));
    pos[i].z += dt*(vel[i].z + dt2*(acc[i].z + dt3*jrk[i].z));

    vel[i].x += dt*(acc[i].x + dt2*jrk[i].x);
    vel[i].y += dt*(acc[i].y + dt2*jrk[i].y);
    vel[i].z += dt*(acc[i].z + dt2*jrk[i].z);
  }

  forces(n, pos, vel, acc, jrk, eps2);

  if (dt > 0.0)
  {
    const real h    = dt*0.5;
    const real hinv = 1.0/h;
    const real f1   = 0.5*hinv*hinv;
    const real f2   = 3.0*hinv*f1;

    const real dt2  = dt *dt * (1.0/2.0);
    const real dt3  = dt2*dt * (1.0/3.0);
    const real dt4  = dt3*dt * (1.0/4.0);
    const real dt5  = dt4*dt * (1.0/5.0);

    for (int i = 0; i < n; i++)
    {

      /* compute snp & crk */

      const real4 Am(   acc[i].x - acc0[i].x,     acc[i].y - acc0[i].y,     acc[i].z - acc0[i].z);
      const real4 Jm(h*(jrk[i].x - jrk0[i].x), h*(jrk[i].y - jrk0[i].y), h*(jrk[i].z - jrk0[i].z));
      const real4 Jp(h*(jrk[i].x + jrk0[i].x), h*(jrk[i].y + jrk0[i].y), h*(jrk[i].z + jrk0[i].z));
      real4 snp(f1* Jm.x,         f1* Jm.y,         f1* Jm.z        );
      real4 crk(f2*(Jp.x - Am.x), f2*(Jp.y - Am.y), f2*(Jp.z - Am.z));

      snp.x -= h*crk.x;
      snp.y -= h*crk.y;
      snp.z -= h*crk.z;

      /* correct */

      pos[i].x += dt4*snp.x + dt5*crk.x;
      pos[i].y += dt4*snp.y + dt5*crk.y;
      pos[i].z += dt4*snp.z + dt5*crk.z;

      vel[i].x += dt3*snp.x + dt4*crk.x;
      vel[i].y += dt3*snp.y + dt4*crk.y;
      vel[i].z += dt3*snp.z + dt4*crk.z;
    }
  }

  return timestep(n, vel, acc, eta, eps2);
}


void integrate(
    const int n,
    real4 pos[],
    real4 vel[],
    const real eta,
    const real eps2,
    const real t_end) {

  static real4 acc[131072] __attribute__((aligned(32)));
  static real4 jrk[131072] __attribute__((aligned(32)));


  const double tin = get_time();
  forces(n, &pos[0], &vel[0], &acc[0], &jrk[0], eps2);
  const double fn = n;
  fprintf(stderr, " mean flop rate in %g sec [%g GFLOP/s]\n", get_time() - tin,
      fn*fn*60/(get_time() - tin)/1e9);

  real Epot0, Ekin0;
  energy(n, &pos[0], &vel[0], &acc[0], Ekin0, Epot0);
  const real Etot0 = Epot0 + Ekin0;
  fprintf(stderr, " E: %g %g %g \n", Epot0, Ekin0, Etot0);

  /////////

  real t_global = 0;
  double t0 = 0;
  int iter = 0;
  int ntime = 10;
  real dt = 0;
  real Epot, Ekin, Etot = Etot0;
  while (t_global < t_end) {
    if (iter % ntime == 0) 
      t0 = get_time();

    dt = iterate(n, &pos[0], &vel[0], &acc[0], &jrk[0], eta, eps2, dt);
    iter++;
    t_global += dt;

    const real Etot_pre = Etot;
    energy(n, &pos[0], &vel[0], &acc[0], Ekin, Epot);
    Etot = Ekin + Epot;

    if (iter % 1 == 0) {
      const real Etot = Ekin + Epot;
      fprintf(stderr, "iter= %d: t= %g  dt= %g Ekin= %g  Epot= %g  Etot= %g , dE = %g d(dE)= %g \n",
          iter, t_global, dt, Ekin, Epot, Etot, (Etot - Etot0)/std::abs(Etot0),
          (Etot - Etot_pre)/std::abs(Etot_pre)   );
    }

    if (iter % ntime == 0) {
      fprintf(stderr, " mean flop rate in %g sec [%g GFLOP/s]\n", get_time() - t0,
          fn*fn*60/(get_time() - t0)/1e9*ntime);
    }

  }
};

int main(int argc, char *argv[]) {
//  std::vector<real4> pos, vel;
  static real4 pos[131072] __attribute__((aligned(32)));
  static real4 vel[131072] __attribute__((aligned(32)));

  int nbodies = 131072;
#if 1
  nbodies = 32768;
#endif
  nbodies = 1024;
  nbodies = 2048;
  nbodies = 4096;
  nbodies = 8192;
  fprintf(stderr, "nbodies= %d\n", nbodies);
  const real R0 = 1;
  const real mp = 1.0/nbodies;
  for (int i = 0; i < nbodies; i++) {
    real xp, yp, zp, s2 = 2*R0;
    real vx, vy, vz;
    while (s2 > R0*R0) {
      xp = (1.0 - 2.0*drand48())*R0;
      yp = (1.0 - 2.0*drand48())*R0;
      zp = (1.0 - 2.0*drand48())*R0;
      s2 = xp*xp + yp*yp + zp*zp;
      vx = drand48() * 0.1;
      vy = drand48() * 0.1;
      vz = drand48() * 0.1;
    } 
    pos[i] = (real4(xp, yp, zp, mp));
    vel[i] = (real4(vx, vy, vz, -1));
  }
  const real eps2 = 4.0f/nbodies;
  real eta  = 0.01f;

  if (argc > 1) eta *= atof(argv[1]);
  fprintf(stderr, " eta= %g \n", eta);
  fprintf(stderr, " stargting ... \n");
  const real tend = 1.0;
  integrate(nbodies, pos, vel, eta, eps2, tend);

}

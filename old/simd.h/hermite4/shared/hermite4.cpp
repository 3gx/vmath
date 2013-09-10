/* Hermite4 N-body integrator */
/* Makino and Aarseth, 1992 */
/* http://adsabs.harvard.edu/abs/1992PASJ...44..141M and references there in*/
/* ISPC/ICC version with XeonPhi support */
/* v1.0, works with single/double precision, compiles with icpc and ispc */
/* (c) Evghenii Gaburov, 2013 */

#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <vector>
#include <cassert>
#include <omp.h>

#include <sys/time.h>
inline double rtc() {
  struct timeval Tvalue;
  struct timezone dummy;

  gettimeofday(&Tvalue,&dummy);
  return ((double) Tvalue.tv_sec +1.e-6*((double) Tvalue.tv_usec));
}

#include "templates.h"
#include "hermite4_ispc.h"

template<typename real>
static inline void body_body_force(
    real &ax, real &ay, real &az, real &gp,
    real &jx, real &jy, real &jz,
    const real xi, const real yi, const real zi,
    const real vxi, const real vyi, const real vzi,
    const real xj, const real yj, const real zj,
    const real vxj, const real vyj, const real vzj, const real mj)
{
  const real dx = xj - xi;
  const real dy = yj - yi;
  const real dz = zj - zi;

  const real ds2 = dx*dx + dy*dy + dz*dz;

  const real  inv_ds  = ds2 > 0.0 ? (real)(1.0f/sqrtf((float)ds2)) : 0.0;
  const real  inv_ds2 = inv_ds*inv_ds;
  const real minv_ds  = inv_ds  * mj;
  const real minv_ds3 = inv_ds2 * minv_ds;

  ax += minv_ds3 * dx;
  ay += minv_ds3 * dy;
  az += minv_ds3 * dz;
  gp -= minv_ds;       // 21

  const real dvx = vxj - vxi;
  const real dvy = vyj - vyi;
  const real dvz = vzj - vzi;
  const real rv  = dx*dvx + dy*dvy + dz*dvz;

  const real Jij = (-3.0) * (rv * inv_ds2 * minv_ds3);

  jx += minv_ds3*dvx + Jij*dx;
  jy += minv_ds3*dvy + Jij*dy;
  jz += minv_ds3*dvz + Jij*dz;  // 23
} // 44 flops/interaction

template<typename real>
static void compute_forces(
    const int     n,
    const int     n1, const int n2,
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
    real gpot[])
{
  for (int i = n1; i < n2; i++)
  {
    real ax = 0.0;
    real ay = 0.0;
    real az = 0.0;
    real gp = 0.0;
    real jx = 0.0;
    real jy = 0.0;
    real jz = 0.0;

    const real  xi = posx[i];
    const real  yi = posy[i];
    const real  zi = posz[i];
    const real vxi = velx[i];
    const real vyi = vely[i];
    const real vzi = velz[i];

#pragma ivdep 
    for (int j = 0; j < n; j++)
    {
      const real  xj = posx[j];
      const real  yj = posy[j];
      const real  zj = posz[j];
      const real vxj = velx[j];
      const real vyj = vely[j];
      const real vzj = velz[j];
      const real  mj = mass[j];
      body_body_force<real>(ax,ay,az,gp,jx,jy,jz,
          xi,yi,zi,vxi,vyi,vzi,
          xj,yj,zj,vxj,vyj,vzj,mj);
    }
    accx[i] = ax;
    accy[i] = ay;
    accz[i] = az;
    jrkx[i] = jx;
    jrky[i] = jy;
    jrkz[i] = jz;
    gpot[i] = gp;
  }
}

template<typename real>
struct Hermite4
{
  enum {PP_FLOP=44};
  const int n;
  const real eta;
  const bool use_ispc;
  real *g_mass, *g_gpot;
  real *g_posx, *g_posy, *g_posz;
  real *g_velx, *g_vely, *g_velz;
  real *g_accx, *g_accy, *g_accz;
  real *g_jrkx, *g_jrky, *g_jrkz;

  std::vector<real> accx0, accy0, accz0;
  std::vector<real> jrkx0, jrky0, jrkz0;

  Hermite4(const int _n = 8192, const real _eta = 0.1, const bool _ispc = true) : n(_n), eta(_eta), use_ispc(_ispc)
  {
    const int alignment = 64;
    const int size = n*sizeof(real);
    assert(posix_memalign((void**)&g_mass, alignment, size) == 0);
    assert(posix_memalign((void**)&g_gpot, alignment, size) == 0);
    assert(posix_memalign((void**)&g_posx, alignment, size) == 0);
    assert(posix_memalign((void**)&g_posy, alignment, size) == 0);
    assert(posix_memalign((void**)&g_posz, alignment, size) == 0);
    assert(posix_memalign((void**)&g_velx, alignment, size) == 0);
    assert(posix_memalign((void**)&g_vely, alignment, size) == 0);
    assert(posix_memalign((void**)&g_velz, alignment, size) == 0);
    assert(posix_memalign((void**)&g_accx, alignment, size) == 0);
    assert(posix_memalign((void**)&g_accy, alignment, size) == 0);
    assert(posix_memalign((void**)&g_accz, alignment, size) == 0);
    assert(posix_memalign((void**)&g_jrkx, alignment, size) == 0);
    assert(posix_memalign((void**)&g_jrky, alignment, size) == 0);
    assert(posix_memalign((void**)&g_jrkz, alignment, size) == 0);

    accx0.resize(n);
    accy0.resize(n);
    accz0.resize(n);
    jrkx0.resize(n);
    jrky0.resize(n);
    jrkz0.resize(n);

    printf("---Intializing nbody--- \n");
    int nthreads = 0;
#pragma omp parallel
#pragma omp critical
    nthreads++;
    printf("nthreads= %d\n", nthreads);

    const real R0 = 1;
    const real mp = 1.0/n;
    for (int i = 0; i < n; i++) {
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
      g_posx[i] = xp;
      g_posy[i] = yp;
      g_posz[i] = zp;
      g_velx[i] = vx;
      g_vely[i] = vy;
      g_velz[i] = vz;
      g_mass[i] = mp;
    }

  }

  ~Hermite4()
  {
    free(g_mass);
    free(g_gpot);
    free(g_posx);
    free(g_posy);
    free(g_posz);
    free(g_velx);
    free(g_vely);
    free(g_velz);
    free(g_accx);
    free(g_accy);
    free(g_accz);
    free(g_jrkx);
    free(g_jrky);
    free(g_jrkz);
  }

  void forces();


  real step(const real dt) 
  {
    const real dt2 = dt*real(1.0/2.0);
    const real dt3 = dt*real(1.0/3.0);

    real dt_min = HUGE;

#pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
      accx0[i] = g_accx[i];
      accy0[i] = g_accy[i];
      accz0[i] = g_accz[i];
      jrkx0[i] = g_jrkx[i];
      jrky0[i] = g_jrky[i];
      jrkz0[i] = g_jrkz[i];

      g_posx[i] += dt*(g_velx[i] + dt2*(g_accx[i] + dt3*g_jrkx[i]));
      g_posy[i] += dt*(g_vely[i] + dt2*(g_accy[i] + dt3*g_jrky[i]));
      g_posz[i] += dt*(g_velz[i] + dt2*(g_accz[i] + dt3*g_jrkz[i]));

      g_velx[i] += dt*(g_accx[i] + dt2*g_jrkx[i]);
      g_vely[i] += dt*(g_accy[i] + dt2*g_jrky[i]);
      g_velz[i] += dt*(g_accz[i] + dt2*g_jrkz[i]);
    }

    forces();

    if (dt > 0.0)
    {
      const real h    = dt*real(0.5);
      const real hinv = real(1.0)/h;
      const real f1   = real(0.5)*hinv*hinv;
      const real f2   = real(3.0)*hinv*f1;

      const real dt2  = dt *dt * real(1.0/2.0);
      const real dt3  = dt2*dt * real(1.0/3.0);
      const real dt4  = dt3*dt * real(1.0/4.0);
      const real dt5  = dt4*dt * real(1.0/5.0);

#pragma omp parallel 
      {
        real dt_min_loc = HUGE;
#pragma omp for
        for (int i = 0; i < n; i++)
        {
          /* compute snp & crk */

          const real Amx = g_accx[i] - accx0[i]; 
          const real Amy = g_accy[i] - accy0[i]; 
          const real Amz = g_accz[i] - accz0[i]; 

          const real Jmx = h*(g_jrkx[i] - jrkx0[i]);
          const real Jmy = h*(g_jrky[i] - jrky0[i]);
          const real Jmz = h*(g_jrkz[i] - jrkz0[i]);

          const real Jpx = h*(g_jrkx[i] + jrkx0[i]);
          const real Jpy = h*(g_jrky[i] + jrky0[i]);
          const real Jpz = h*(g_jrkz[i] + jrkz0[i]);


          real snpx = f1*Jmx;
          real snpy = f1*Jmy;
          real snpz = f1*Jmz;

          real crkx = f2*(Jpx - Amx);
          real crky = f2*(Jpy - Amy);
          real crkz = f2*(Jpz - Amz);

          snpx -= h*crkx;
          snpy -= h*crky;
          snpz -= h*crkz;

          /* correct */

          g_posx[i] += dt4*snpx + dt5*crkx;
          g_posy[i] += dt4*snpy + dt5*crky;
          g_posz[i] += dt4*snpz + dt5*crkz;

          g_velx[i] += dt3*snpx + dt4*crkx;
          g_vely[i] += dt3*snpy + dt4*crky;
          g_velz[i] += dt3*snpz + dt4*crkz;

          /* compute new timestep */

          const real s0 = g_accx[i]*g_accx[i] + g_accy[i]*g_accy[i] + g_accz[i]*g_accz[i];
          const real s1 = g_jrkx[i]*g_jrkx[i] + g_jrky[i]*g_jrky[i] + g_jrkz[i]*g_jrkz[i];
          const real s2 = snpx*snpx + snpy*snpy + snpz*snpz;
          const real s3 = crkx*crkx + crky*crky + crkz*crkz;

          const double u = std::sqrt(s0*s2) + s1;
          const double l = std::sqrt(s1*s3) + s2;
          assert(l > 0.0f);
          const real dt_loc = eta *std::sqrt(u/l);
          dt_min_loc = std::min(dt_min_loc, dt_loc);
        }
#pragma omp critical
        dt_min = std::min(dt_min, dt_min_loc);
      }
    }

    if (dt_min == HUGE) return dt;
    else return dt_min;
  }

  void energy(
      real &Ekin, real &Epot) 
  {
    real ekin = 0, epot = 0;

#pragma omp parallel for reduction(+:ekin,epot)
    for (int i = 0; i < n; i++) 
    {
      ekin += g_mass[i] * (g_velx[i]*g_velx[i] + g_vely[i]*g_vely[i] + g_velz[i]*g_velz[i]) * real(0.5f);
      epot += real(0.5f)*g_mass[i] * g_gpot[i];
    }
    Ekin = ekin;
    Epot = epot;
  }

  void integrate(const int niter, const real t_end = HUGE)
  {
    const double tin = rtc();
    forces();
    const double fn = n;
    printf(" mean flop rate in %g sec [%g GFLOP/s]\n", rtc() - tin,
        fn*fn*PP_FLOP/(rtc() - tin)/1e9);

    real Epot0, Ekin0;
    energy(Ekin0, Epot0);
    const real Etot0 = Epot0 + Ekin0;
    printf(" E: %g %g %g \n", Epot0, Ekin0, Etot0);

    /////////

    real t_global = 0;
    double t0 = 0;
    int iter = 0;
    int ntime = 10;
    real dt = 1.0/131072;
    real Epot, Ekin, Etot = Etot0;
    while (t_global < t_end) {
      if (iter % ntime == 0) 
        t0 = rtc();

      if (iter >= niter) return;

      dt = step(dt);
      iter++;
      t_global += dt;

      const real Etot_pre = Etot;
      energy(Ekin, Epot);
      Etot = Ekin + Epot;

      if (iter % 1 == 0) {
        const real Etot = Ekin + Epot;
        printf("iter= %d: t= %g  dt= %g Ekin= %g  Epot= %g  Etot= %g , dE = %g d(dE)= %g \n",
            iter, t_global, dt, Ekin, Epot, Etot, (Etot - Etot0)/std::abs(Etot0),
            (Etot - Etot_pre)/std::abs(Etot_pre)   );
      }

      if (iter % ntime == 0) {
        printf(" mean flop rate in %g sec [%g GFLOP/s]\n", rtc() - t0,
            fn*fn*PP_FLOP/(rtc() - t0)/1e9*ntime);
      }

      fflush(stdout);

    }

  }


};

  template<>
void Hermite4<float>::forces()
{
#pragma omp parallel
  {
    const int blockDim = 64;
    const int gridDim  = omp_get_num_threads();
    const int blockIdx = omp_get_thread_num();

    for (int i = blockIdx*blockDim; i < n; i += blockDim*gridDim)
    {
      if (use_ispc)
        TEMPLATE(ispc::compute_forces,float)(
            n,
            i, std::min(i + blockDim,n),
            &g_mass[0],
            &g_posx[0],
            &g_posy[0],
            &g_posz[0],
            &g_velx[0],
            &g_vely[0],
            &g_velz[0],
            &g_accx[0],
            &g_accy[0],
            &g_accz[0],
            &g_jrkx[0],
            &g_jrky[0],
            &g_jrkz[0],
            &g_gpot[0]);
      else
        compute_forces<float>(
            n,
            i, std::min(i + blockDim,n),
            &g_mass[0],
            &g_posx[0],
            &g_posy[0],
            &g_posz[0],
            &g_velx[0],
            &g_vely[0],
            &g_velz[0],
            &g_accx[0],
            &g_accy[0],
            &g_accz[0],
            &g_jrkx[0],
            &g_jrky[0],
            &g_jrkz[0],
            &g_gpot[0]);
    }
  }
}

  template<>
void Hermite4<double>::forces()
{
#pragma omp parallel
  {
    const int blockDim = 32;
    const int gridDim  = omp_get_num_threads();
    const int blockIdx = omp_get_thread_num();

    for (int i = blockIdx*blockDim; i < n; i += blockDim*gridDim)
    {
      if (use_ispc)
        TEMPLATE(ispc::compute_forces,double)(
            n,
            i, std::min(i + blockDim,n),
            &g_mass[0],
            &g_posx[0],
            &g_posy[0],
            &g_posz[0],
            &g_velx[0],
            &g_vely[0],
            &g_velz[0],
            &g_accx[0],
            &g_accy[0],
            &g_accz[0],
            &g_jrkx[0],
            &g_jrky[0],
            &g_jrkz[0],
            &g_gpot[0]);
      else
        compute_forces<double>(
            n,
            i, std::min(i + blockDim,n),
            &g_mass[0],
            &g_posx[0],
            &g_posy[0],
            &g_posz[0],
            &g_velx[0],
            &g_vely[0],
            &g_velz[0],
            &g_accx[0],
            &g_accy[0],
            &g_accz[0],
            &g_jrkx[0],
            &g_jrky[0],
            &g_jrkz[0],
            &g_gpot[0]);
    }
  }

}


  template<typename real>
void run(const int nbodies, const real eta, const int nstep, const bool use_ispc = true)
{
  Hermite4<real> h4(nbodies, eta, use_ispc);
  h4.integrate(nstep);
}

int main(int argc, char *argv[]) 
{
  if (argc <= 2)
  {
    printf("------------------- \n");
    printf("  Usage: %s fp64 icc [nbodies=8192] [nsteps=40] [eta=0.1] \n", argv[0]);
    printf("    fp64= 1 to use double precision, otherwise use single precision \n"); 
    printf("    icc=  1 to use icc_code, otherwise use ispc code \n");
    printf("------------------- \n");
    return -1;
  }

  const bool fp64 = atoi(argv[1]) == 1;
  const bool icc  = atoi(argv[2]) == 1;

  int nbodies = 8192;
  if (argc > 3) nbodies = atoi(argv[3]);

  int nstep = 40;
  if (argc > 4) nstep = atoi(argv[4]);

  float eta = 0.1;
  if (argc > 5) eta = atof(argv[5]);


  printf("precision: %s \n", fp64 ? "double" : "float"); 
  printf("code: %s \n", icc ? "icc" : "ispc");
  printf("nbodies= %d\n", nbodies);
  printf("nstep= %d\n", nstep);
  printf(" eta= %g \n", eta);

  if (fp64) run<double>(nbodies, eta, nstep, !icc);
  else      run<float >(nbodies, eta, nstep, !icc);

  return 0;
}


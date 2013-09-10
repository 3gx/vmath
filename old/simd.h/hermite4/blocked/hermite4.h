#include <cassert>
#include "pownth.h"
#include <omp.h>
#include "hermite4_ispc.h"

#include <sys/time.h>
static double wtime(){
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return double(tv.tv_sec) + 1.e-6*double(tv.tv_usec);
}

struct Force{
  dvec3 acc;
  dvec3 jrk;
};

struct Particle
{
  enum {
    order = 4,
    flops = 44
  };
  long  id;
  double mass;
  double tlast;
  double dt;
  dvec3 pos;
  dvec3 vel;
  dvec3 acc;
  dvec3 jrk;
  dvec3 snp; // now save all the derivatives
  dvec3 crk;
  double pad[2]; // 24 DP words

  Particle(
      const long    _id, 
      const double  _mass,
      const dvec3 & _pos,
      const dvec3 & _vel)
    : id(_id), mass(_mass), tlast(0.0), dt(0.0), 
    pos(_pos), vel(_vel), acc(0.0), jrk(0.0)
  {
    assert(sizeof(Particle) == 192);
    pad[0] = pad[1] = 0.0;
  }

  void dump(FILE *fp){
    fprintf(fp, "%8ld %A", id, mass);
    fprintf(fp, "   %A %A %A", pos.x, pos.y, pos.z);
    fprintf(fp, "   %A %A %A", vel.x, vel.y, vel.z);
    fprintf(fp, "   %A %A %A", acc.x, acc.y, acc.z);
    fprintf(fp, "   %A %A %A", jrk.x, jrk.y, jrk.z);
    fprintf(fp, "   %A %A %A", snp.x, snp.y, snp.z);
    fprintf(fp, "   %A %A %A", crk.x, crk.y, crk.z);
    fprintf(fp, "\n");
  }
  void restore(FILE *fp){
    int nread = 0;
    nread += fscanf(fp, "%ld %lA", &id, &mass);
    nread += fscanf(fp, "%lA %lA %lA", &pos.x, &pos.y, &pos.z);
    nread += fscanf(fp, "%lA %lA %lA", &vel.x, &vel.y, &vel.z);
    nread += fscanf(fp, "%lA %lA %lA", &acc.x, &acc.y, &acc.z);
    nread += fscanf(fp, "%lA %lA %lA", &jrk.x, &jrk.y, &jrk.z);
    nread += fscanf(fp, "%lA %lA %lA", &snp.x, &snp.y, &snp.z);
    nread += fscanf(fp, "%lA %lA %lA", &crk.x, &crk.y, &crk.z);
    assert(20 == nread);
  }

  void assign_force(const Force &f){
    acc = f.acc;
    jrk = f.jrk;
  }

  void init_dt(const double eta_s, const double dtmax){
    const double dtnat = eta_s * sqrt(acc.norm2() / jrk.norm2());
    dt = 0.25 * dtmax;
    while(dt > dtnat) dt *= 0.5;
  }

  static double aarseth_step_quant(
      const dvec3 &a0, 
      const dvec3 &a1, 
      const dvec3 &a2, 
      const dvec3 &a3, 
      const double etapow)
  {
    double s0 = a0.norm2();
    double s1 = a1.norm2();
    double s2 = a2.norm2();
    double s3 = a3.norm2();

    double u = sqrt(s0*s2) + s1;
    double l = sqrt(s1*s3) + s2;
    const double dtpow = etapow * (u/l);
    return pow_one_nth_quant<2>(dtpow);
  }
  static double aarseth_step(
      const dvec3 &a0, 
      const dvec3 &a1, 
      const dvec3 &a2, 
      const dvec3 &a3, 
      const double eta)
  {
    double s0 = a0.norm2();
    double s1 = a1.norm2();
    double s2 = a2.norm2();
    double s3 = a3.norm2();

    double u = sqrt(s0*s2) + s1;
    double l = sqrt(s1*s3) + s2;
    return eta * sqrt(u/l);
  }

  void recalc_dt(const double eta, const double dtmax){
    const double dta = aarseth_step(acc, jrk, snp, crk, eta);
    dt = dtmax;
    while(dt > dta) dt *= 0.5;
  }

  void correct(const Force &f, const double eta, const double etapow, const double dtlim){
    const double h = 0.5 * dt;
    const double hinv = 2.0/dt;

    const dvec3 Ap = (f.acc + acc);
    const dvec3 Am = (f.acc - acc);
    const dvec3 Jp = (f.jrk + jrk)*h;
    const dvec3 Jm = (f.jrk - jrk)*h;

    // do correct
    dvec3 vel1 = vel + h*(Ap - (1./3.)*Jm);
    pos += h*((vel + vel1) + h*(-(1./3.)*Am));
    vel = vel1;
    acc = f.acc;
    jrk = f.jrk;
    tlast += dt;

    // update timestep
    snp = (0.5*hinv*hinv) * Jm;
    crk = (1.5*hinv*hinv*hinv) * (Jp - Am);
    snp += h * crk;
#if 0
    const double dta = aarseth_step(acc, jrk, snp, crk, eta);
    dt = dtlim;
    while(dt > dta) dt *= 0.5;
#else
    const double dtq = aarseth_step_quant(acc, jrk, snp, crk, etapow);
    dt = dtq<dtlim ? dtq : dtlim;
#endif
  }

} __attribute__((aligned(64)));


struct Gravity{
  typedef Particle GParticle;
  struct GPredictor{
    dvec3  pos;
    double mass;
    dvec3  vel;
    long   id;

    GPredictor(const GParticle &p, const double tsys){
      const double dt = tsys - p.tlast;
      const double dt2 = (1./2.) * dt;
      const double dt3 = (1./3.) * dt;

      pos  = p.pos + dt * (p.vel + dt2 * (p.acc + dt3 * p.jrk));
      vel  = p.vel + dt * (p.acc + dt2 * (p.jrk));
      mass = p.mass;
      id   = p.id;
    }
  };

  struct GPredictor_SoA
  {
    double *x;
    double *y;
    double *z;
    double *vx;
    double *vy;
    double *vz;
    double *mass;
    GPredictor_SoA() : x(NULL), y(NULL), z(NULL), vx(NULL), vy(NULL), vz(NULL), mass(NULL) {}
    ~GPredictor_SoA() 
    {
      if (x == NULL) return;
      free(x);
      free(y);
      free(z);
      free(vx);
      free(vy);
      free(vz);
      free(mass);
    }
    void allocate(const int nbody) 
    {
      x    = ::allocate<double,64>(nbody);
      y    = ::allocate<double,64>(nbody);
      z    = ::allocate<double,64>(nbody);
      vx   = ::allocate<double,64>(nbody);
      vy   = ::allocate<double,64>(nbody);
      vz   = ::allocate<double,64>(nbody);
      mass = ::allocate<double,64>(nbody);
    }
  };

  struct GForce_SoA
  {
    double *ax;
    double *ay;
    double *az;
    double *jx;
    double *jy;
    double *jz;
    GForce_SoA() : ax(NULL), ay(NULL), az(NULL), jx(NULL), jy(NULL), jz(NULL) {}

    void allocate(const int nbody)
    {
      ax = ::allocate<double,64>(nbody);
      ay = ::allocate<double,64>(nbody);
      az = ::allocate<double,64>(nbody);
      jx = ::allocate<double,64>(nbody);
      jy = ::allocate<double,64>(nbody);
      jz = ::allocate<double,64>(nbody);
    }
    ~GForce_SoA() 
    {
      if (ax == NULL) return;
      free(ax);
      free(ay);
      free(az);
      free(jx);
      free(jy);
      free(jz);
    }
  };

  const int  nbody;
  GParticle  *ptcl;
  GPredictor *pred;
  GPredictor_SoA pred_SoA;
  GForce_SoA force_SoA;


  Gravity(const int _nbody) : nbody(_nbody) {
    ptcl = allocate<GParticle,  64> (nbody);
    pred = allocate<GPredictor, 64> (nbody);
    pred_SoA.allocate(nbody);
    force_SoA.allocate(nbody);
  }
  ~Gravity(){
    free(ptcl);
    free(pred);
  }

  void set_jp(const int addr, const Particle &p){
    ptcl[addr] = p;
  }

  void predict_all(const double tsys){
#pragma omp parallel for
    for(int i=0; i<nbody; i++)
    {
      pred[i] = GPredictor(ptcl[i], tsys);
      pred_SoA.   x[i] = pred[i].pos.x;
      pred_SoA.   y[i] = pred[i].pos.y;
      pred_SoA.   z[i] = pred[i].pos.z;
      pred_SoA.  vx[i] = pred[i].vel.x;
      pred_SoA.  vy[i] = pred[i].vel.y;
      pred_SoA.  vz[i] = pred[i].vel.z;
      pred_SoA.mass[i] = pred[i].mass;
    }
  }
 

  void calcforce(
      const int ibeg, const int iend,
      const int ioffset,
      const int jbeg, const int jend,
      const int dnj,
      const double eps2,
      const GPredictor_SoA &pred,
      double *gax,
      double *gay,
      double *gaz,
      double *gjx,
      double *gjy,
      double *gjz)
  {
    double *xptr = pred.x;
    double *yptr = pred.y;
    double *zptr = pred.z;
    double *vxptr = pred.vx;
    double *vyptr = pred.vy;
    double *vzptr = pred.vz;
    double *massptr = pred.mass;
    for (int i = ibeg; i < iend; i++)
    {
      double ax, ay, az;
      double jx, jy, jz;
      ax = ay = az = 0.0;
      jx = jy = jz = 0.0;

      for (int j = jbeg; j < jend; j += dnj)
      {
        const double  xi = xptr [i];
        const double  yi = yptr [i];
        const double  zi = zptr [i];
        const double vxi = vxptr[i];
        const double vyi = vyptr[i];
        const double vzi = vzptr[i];

        const double  xj = xptr[j];
        const double  yj = yptr[j];
        const double  zj = zptr[j];
        const double vxj = vxptr[j];
        const double vyj = vyptr[j];
        const double vzj = vzptr[j];
        const double  mj = massptr[j];

        const double dx = xj - xi;
        const double dy = yj - yi;
        const double dz = zj - zi;

        const double ds2 = dx*dx + dy*dy + dz*dz + eps2;

        const double inv_ds = 1.0f/sqrtf((float)ds2);
        const double  inv_ds2 = inv_ds  * inv_ds;
        const double minv_ds  = inv_ds  * mj;
        const double minv_ds3 = inv_ds2 * minv_ds;


        ax += minv_ds3 * dx;
        ay += minv_ds3 * dy;
        az += minv_ds3 * dz;

        const double dvx = vxj - vxi;
        const double dvy = vyj - vyi;
        const double dvz = vzj - vzi;
        const double rv  = dx*dvx + dy*dvy + dz*dvz;

        const double Jij = (-3.0) * (rv * inv_ds2 * minv_ds3);

        jx += minv_ds3*dvx + Jij*dx;
        jy += minv_ds3*dvy + Jij*dy;
        jz += minv_ds3*dvz + Jij*dz;
      }
      gax[i-ioffset] = ax;
      gay[i-ioffset] = ay;
      gaz[i-ioffset] = az;
      gjx[i-ioffset] = jx;
      gjy[i-ioffset] = jy;
      gjz[i-ioffset] = jz;
    }
  }

  template<int NTX, int NTYMAX, int IBLOCK>
    void calc_force_2d(
        const int    nact,
        const double eps2,
        const int nthreads,
        Force        force[])
    {
      const int ni = nact;
      const int nj = nbody;

      static double gax[NTYMAX][NTX][IBLOCK] __attribute__((aligned(64)));
      static double gay[NTYMAX][NTX][IBLOCK] __attribute__((aligned(64)));
      static double gaz[NTYMAX][NTX][IBLOCK] __attribute__((aligned(64)));
      static double gjx[NTYMAX][NTX][IBLOCK] __attribute__((aligned(64)));
      static double gjy[NTYMAX][NTX][IBLOCK] __attribute__((aligned(64)));
      static double gjz[NTYMAX][NTX][IBLOCK] __attribute__((aligned(64)));


      const int ntx = NTX;
      const int nty = nthreads/ntx;

      assert(nthreads%ntx == 0);
      assert(nty <= NTYMAX);

#pragma omp parallel num_threads(nthreads)
      {
        const int tid = omp_get_thread_num();

        const int tx = tid % ntx;
        const int ty = tid / ntx;

#if 1  /* this one appears to be faster */
        const int nj_per_ty = (nj+nty-1)/nty;
        const int jb =          ty * nj_per_ty;
        const int je = std::min(jb + nj_per_ty, nj);
        const int dnj = 1;
#else
        const int jb  = ty;
        const int je  = nj;
        const int dnj = nty;
#endif

        for (int iblock = 0; iblock < ni; iblock += ntx*IBLOCK)
        {
          const int ib = std::min(iblock + tx*IBLOCK,ni);
          const int ie = std::min(ib     +    IBLOCK,ni);
#if 1
          ispc::calcforce_on_first_nact(
              ib, ie, ib,
              jb, je, dnj,
              eps2,
              (ispc::Predictor*)&pred_SoA,
              &gax[ty][tx][0],
              &gay[ty][tx][0],
              &gaz[ty][tx][0],
              &gjx[ty][tx][0],
              &gjy[ty][tx][0],
              &gjz[ty][tx][0]);
#else
          calcforce(
              ib, ie, ib,
              jb, je, dnj,
              eps2,
              pred_SoA,
              &gax[ty][tx][0],
              &gay[ty][tx][0],
              &gaz[ty][tx][0],
              &gjx[ty][tx][0],
              &gjy[ty][tx][0],
              &gjz[ty][tx][0]);
#endif

#pragma omp barrier

          for (int i = ty; i < ie-ib; i += nty)
          {
            double ax = 0.0;
            double ay = 0.0;
            double az = 0.0;
            double jx = 0.0;
            double jy = 0.0;
            double jz = 0.0;
            for (int j = 0; j < nty; j++)
            {
              ax += gax[j][tx][i];
              ay += gay[j][tx][i];
              az += gaz[j][tx][i];
              jx += gjx[j][tx][i];
              jy += gjy[j][tx][i];
              jz += gjz[j][tx][i];
            }
            force[ib+i].acc = dvec3(ax,ay,az);
            force[ib+i].jrk = dvec3(jx,jy,jz);
          }

#pragma omp barrier
        }
      }
    }


  void calc_force_1d(
      const int    nact,
      const double eps2,
      Force        force[] )
  {
    const int ni = nact;
    const int nj = nbody;

#pragma omp parallel
    {
      const int blockDim = 16;
      const int gridDim  = omp_get_num_threads();
      const int blockIdx = omp_get_thread_num();

      for (int ib = blockIdx*blockDim; ib < ni; ib += blockDim*gridDim)
      {
        const int ie = std::min(ib+blockDim,ni);
#if 1
        ispc::calcforce_on_first_nact(
            ib, ie, 0,
            0,nj,1,
            eps2,
            (ispc::Predictor*)&pred_SoA,
            force_SoA.ax,
            force_SoA.ay,
            force_SoA.az,
            force_SoA.jx,
            force_SoA.jy,
            force_SoA.jz);
#else
        calcforce(
            ib, ie, 0,
            0,nj,1,
            eps2,
            pred_SoA,
            force_SoA.ax,
            force_SoA.ay,
            force_SoA.az,
            force_SoA.jx,
            force_SoA.jy,
            force_SoA.jz);
#endif


#pragma ivdep
        for (int i = ib; i < ie; i++)
        {
          force[i].acc = dvec3(force_SoA.ax[i], force_SoA.ay[i], force_SoA.az[i]);
          force[i].jrk = dvec3(force_SoA.jx[i], force_SoA.jy[i], force_SoA.jz[i]);
        }
      }
    }
  }


  void calc_force_on_first_nact(
      const int    nact,
      const double eps2,
      Force        force[] )
  {
    const int ni = nact;
    const int nj = nbody;

#if 1
    static int nthreads = 0;
    if (nthreads == 0)
    {
#pragma omp parallel
#pragma omp master
      nthreads = omp_get_num_threads();
    }

    const double t0 = wtime();

    const int NTX    =  4;   /* partion thread-pool into 2d grid of size NTX x NTY */
    const int IBLOCK = 16;   /* number of particles processed in an i-block, multiple of SIMD width */
    const int FUDGEF = 16;   /* fudge factor */
    if (nact < FUDGEF*NTX*IBLOCK)
    {
      const int NTYMAX = 256;
      calc_force_2d<NTX,NTYMAX,IBLOCK>(nact,eps2,nthreads,force);
    }
    else
    {
      calc_force_1d(nact,eps2,force);
    }

    const double t1 = wtime();
    if (ni == nj)
      printf("done in %g sec: ni= %d :: %g GFLOPs  \n", t1-t0, ni,
          44.0*ni*nj/(t1-t0)/1e9);
    //    assert(0);

#else  /* naive C++ code */

#pragma omp parallel for
    for(int i=0; i<ni; i++){
      const dvec3 posi = pred[i].pos;
      const dvec3 veli = pred[i].vel;
      dvec3 acc(0.0);
      dvec3 jrk(0.0);
      for(int j=0; j<nj; j++){
        const dvec3  dr    = pred[j].pos - posi;
        const dvec3  dv    = pred[j].vel - veli;
        const double r2    = eps2 + dr*dr;
        const double drdv  = dr*dv;
        const double invr  = 1.0f/sqrtf((float)r2);
        const double ri2   = invr*invr;
        const double mri3  = pred[j].mass * ri2 * invr;
        const double alpha = -3.0 * drdv * ri2;
        acc += mri3 * dr;
        jrk += mri3 * dv + alpha * (mri3 * dr);
      }
      force[i].acc = acc;
      force[i].jrk = jrk;
    }
#endif
  }

  void calc_potential(
      const double eps2,
      double       potbuf[] )
  {
    const int ni = nbody;
    const int nj = nbody;
#pragma omp parallel for
    for(int i=0; i<ni; i++){
      double pot = 0.0;
      const dvec3 posi = ptcl[i].pos;
      for(int j=0; j<nj; j++){
        // if(j == i) continue;
        const dvec3  posj = ptcl[j].pos;
        const dvec3  dr   = posj - posi;
        const double r2   = eps2 + dr*dr;
#if 0
        const double mj = (j != i) ? ptcl[j].mass : 0.0;
#else
        const double mj = (r2 > eps2) ? ptcl[j].mass : 0.0;
#endif
        pot -= mj *(1.0f / sqrtf((float)r2));
      }
      potbuf[i] = pot;
    }
  }
};

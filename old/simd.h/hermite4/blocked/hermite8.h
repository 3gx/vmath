#include <cassert>
#include <omp.h>
#include "pownth.h"
#include "hermite8_ispc.h"

#include <sys/time.h>
static double wtime(){
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return double(tv.tv_sec) + 1.e-6*double(tv.tv_usec);
}

struct Force{
	dvec3 acc;
	dvec3 jrk;
	dvec3 snp;
	dvec3 crk;
};

struct Particle{
	enum {
		order = 8,
		flops = 144,
	};
	long  id;
	double mass;
	double tlast;
	double dt;
	dvec3 pos;
	dvec3 vel;
	dvec3 acc;
	dvec3 jrk;
	dvec3 snp;
	dvec3 crk;
	dvec3 d4a;
	dvec3 d5a;
	dvec3 d6a; // save all the derivatives
	dvec3 d7a;
	double pad[2]; // 36 DP words

	Particle(
			const long    _id, 
			const double  _mass,
			const dvec3 & _pos,
			const dvec3 & _vel)
		: id(_id), mass(_mass), tlast(0.0), dt(0.0), 
		  pos(_pos), vel(_vel), 
		  acc(0.0), jrk(0.0), 
		  snp(0.0), crk(0.0),
		  d4a(0.0), d5a(0.0)
	{
		assert(sizeof(Particle) == 288);
	}

	void dump(FILE *fp){
		fprintf(fp, "%8ld %A", id, mass);
		fprintf(fp, "   %A %A %A", pos.x, pos.y, pos.z);
		fprintf(fp, "   %A %A %A", vel.x, vel.y, vel.z);
		fprintf(fp, "   %A %A %A", acc.x, acc.y, acc.z);
		fprintf(fp, "   %A %A %A", jrk.x, jrk.y, jrk.z);
		fprintf(fp, "   %A %A %A", snp.x, snp.y, snp.z);
		fprintf(fp, "   %A %A %A", crk.x, crk.y, crk.z);
		fprintf(fp, "   %A %A %A", d4a.x, d4a.y, d4a.z);
		fprintf(fp, "   %A %A %A", d5a.x, d5a.y, d5a.z);
		fprintf(fp, "   %A %A %A", d6a.x, d6a.y, d6a.z);
		fprintf(fp, "   %A %A %A", d7a.x, d7a.y, d7a.z);
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
		nread += fscanf(fp, "%lA %lA %lA", &d4a.x, &d4a.y, &d4a.z);
		nread += fscanf(fp, "%lA %lA %lA", &d5a.x, &d5a.y, &d5a.z);
		nread += fscanf(fp, "%lA %lA %lA", &d6a.x, &d6a.y, &d6a.z);
		nread += fscanf(fp, "%lA %lA %lA", &d7a.x, &d7a.y, &d7a.z);
		assert(32 == nread);
	}

	void assign_force(const Force &f){
		acc = f.acc;
		jrk = f.jrk;
		snp = f.snp;
		crk = f.crk;
	}

	void init_dt(const double eta_s, const double dtmax){
		double s0 = acc.norm2();
		double s1 = jrk.norm2();
		double s2 = snp.norm2();
		double s3 = crk.norm2();

		double u = sqrt(s0*s2) + s1;
		double l = sqrt(s1*s3) + s2;

		double dtnat =  eta_s * sqrt(u/l);
		dt = 0.25 * dtmax;
		while(dt > dtnat) dt *= 0.5;
	}

	static double aarseth_step_quant(
			const dvec3 &a0, 
			const dvec3 &a1, 
			const dvec3 &a2, 
			const dvec3 &a3, // not used 
			const dvec3 &a4, // not used
			const dvec3 &a5, 
			const dvec3 &a6, 
			const dvec3 &a7, 
			const double etapow)
	{
		double s0 = a0.norm2();
		double s1 = a1.norm2();
		double s2 = a2.norm2();

		double s5 = a5.norm2();
		double s6 = a6.norm2();
		double s7 = a7.norm2();

		double u = sqrt(s0*s2) + s1;
		double l = sqrt(s5*s7) + s6;
		const double dtpow = etapow * (u/l);
		return pow_one_nth_quant<10>(dtpow);
	}
	static double aarseth_step(
			const dvec3 &a0, 
			const dvec3 &a1, 
			const dvec3 &a2, 
			const dvec3 &a3, // not used 
			const dvec3 &a4, // not used
			const dvec3 &a5, 
			const dvec3 &a6, 
			const dvec3 &a7, 
			const double eta)
	{
		double s0 = a0.norm2();
		double s1 = a1.norm2();
		double s2 = a2.norm2();

		double s5 = a5.norm2();
		double s6 = a6.norm2();
		double s7 = a7.norm2();

		double u = sqrt(s0*s2) + s1;
		double l = sqrt(s5*s7) + s6;
		return eta * pow(u/l, 1./10.);
	}

	void recalc_dt(const double eta, const double dtmax){
		const double dta = aarseth_step(acc, jrk, snp, crk, d4a, d5a, d6a, d7a, eta);
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
		const dvec3 Sp = (f.snp + snp)*(h*h);
		const dvec3 Sm = (f.snp - snp)*(h*h);
		const dvec3 Cp = (f.crk + crk)*(h*h*h);
		const dvec3 Cm = (f.crk - crk)*(h*h*h);

		// do correct
		dvec3 vel1 = vel + h*(Ap +  (-3./7.)*Jm + (2./21)*Sp - (1./105.)*Cm);
		pos += h*((vel + vel1) + h*((-3./7.)*Am + (2./21)*Jp - (1./105.)*Sm));
		vel = vel1;
		acc = f.acc;
		jrk = f.jrk;
		snp = f.snp;
		crk = f.crk;
		tlast += dt;

		// taylor series
		double hinv2 = hinv*hinv; 
		double hinv3 = hinv2*hinv; 
		double hinv4 = hinv2*hinv2;
		double hinv5 = hinv2*hinv3;
		double hinv6 = hinv3*hinv3;
		double hinv7 = hinv3*hinv4;

		d4a = (hinv4 *   24./32.)*(         - 5.*Jm + 5.*Sp - Cm);
		d5a = (hinv5 *  120./32.)*( 21.*Am - 21.*Jp + 8.*Sm - Cp);
		d6a = (hinv6 *  720./32.)*(              Jm -    Sp + Cm/3.);
		d7a = (hinv7 * 5040./32.)*( -5.*Am +  5.*Jp - 2.*Sm + Cp/3.);

		double h2 = 0.5*h, h3 = (1./3.)*h;
		d4a += h*(d5a + h2*(d6a + h3*d7a));
		d5a += h*(d6a + h2*d7a);
		d6a += h*d7a;

		// update timestep
#if 0
		const double dta = aarseth_step(acc, jrk, snp, crk, d4a, d5a, d6a, d7a, eta);
		dt = dtlim;
		while(dt > dta) dt *= 0.5;
#else
		const double dtq = aarseth_step_quant(acc, jrk, snp, crk, d4a, d5a, d6a, d7a, etapow);
		dt = dtq<dtlim ? dtq : dtlim;
#endif
	}

} __attribute__ ((aligned(32)));

struct Gravity{
	typedef Particle GParticle;
	struct GPredictor{
		dvec3  pos;
		double mass;
		dvec3  vel;
		long   id;
		dvec3  acc;
		dvec3  jrk; // 14 DP

		GPredictor(const GParticle &p, const double tsys){
			const double dt = tsys - p.tlast;
			const double dt2 = (1./2.) * dt;
			const double dt3 = (1./3.) * dt;
			const double dt4 = (1./4.) * dt;
			const double dt5 = (1./5.) * dt;
			const double dt6 = (1./6.) * dt;
			const double dt7 = (1./7.) * dt;

			pos  = p.pos + dt * (p.vel + dt2 * (p.acc + dt3 * (p.jrk + dt4 * (p.snp + dt5 * (p.crk + dt6 * (p.d4a + dt7 * (p.d5a)))))));
			vel  = p.vel + dt * (p.acc + dt2 * (p.jrk + dt3 * (p.snp + dt4 * (p.crk + dt5 * (p.d4a + dt6 * (p.d5a))))));
			acc  = p.acc + dt * (p.jrk + dt2 * (p.snp + dt3 * (p.crk + dt4 * (p.d4a + dt5 * (p.d5a)))));
			jrk  = p.jrk + dt * (p.snp + dt2 * (p.crk + dt3 * (p.d4a + dt4 * (p.d5a))));
			mass = p.mass;
			id   = p.id;
		}
	};

  struct GPredictor_SoA
  {
    double* pos[4];  /* x, y, z, mass */
    double* vel[3];
    double* acc[3];
    double* jrk[3];
    GPredictor_SoA() 
    {
      pos[3] = NULL;
      for (int k = 0; k < 3; k++)
        pos[k] = vel[k] = acc[k] = jrk[k] = NULL;
    }
    ~GPredictor_SoA() 
    {
      if (pos[0] == NULL) return;
      for (int k = 0; k < 3; k++)
      {
        free(pos[k]);
        free(vel[k]);
        free(acc[k]);
        free(jrk[k]);
      }
      free(pos[3]);
    }
    void allocate(const int nbody) 
    {
      pos[3] = ::allocate<double,64>(nbody);
      for (int k = 0; k < 3; k++)
      {
        pos[k]  = ::allocate<double,64>(nbody);
        vel[k]  = ::allocate<double,64>(nbody);
        acc[k]  = ::allocate<double,64>(nbody);
        jrk[k]  = ::allocate<double,64>(nbody);
      }
    }
  };

  struct GForce_SoA
  {
    double* acc[3];
    double* jrk[3];
    double* snp[3];
    double* crk[3];
    GForce_SoA() 
    {
      for (int k = 0; k <3; k++)
        acc[k] = jrk[k] = snp[k] = crk[k] = NULL;
    }

    void allocate(const int nbody)
    {
      for (int k = 0;k < 3; k++)
      {
        acc[k] = ::allocate<double,64>(nbody);
        jrk[k] = ::allocate<double,64>(nbody);
        snp[k] = ::allocate<double,64>(nbody);
        crk[k] = ::allocate<double,64>(nbody);
      }
    }
    ~GForce_SoA() 
    {
      if (acc[0] == NULL) return;
      for (int k = 0;k < 3; k++)
      {
        free(acc[k]);
        free(jrk[k]);
        free(snp[k]);
      }
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
		for(int i=0; i<nbody; i++){
			pred[i] = GPredictor(ptcl[i], tsys);
      for (int k = 0; k < 3; k++)
      {
        pred_SoA.pos[k][i] = pred[i].pos[k];
        pred_SoA.vel[k][i] = pred[i].vel[k];
        pred_SoA.acc[k][i] = pred[i].acc[k];
        pred_SoA.jrk[k][i] = pred[i].jrk[k];
      }
      pred_SoA.pos[3][i] = pred[i].mass;
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

      static double gacc[3][NTYMAX][NTX][IBLOCK] __attribute__((aligned(64)));
      static double gjrk[3][NTYMAX][NTX][IBLOCK] __attribute__((aligned(64)));
      static double gsnp[3][NTYMAX][NTX][IBLOCK] __attribute__((aligned(64)));
      static double gcrk[3][NTYMAX][NTX][IBLOCK] __attribute__((aligned(64)));

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
          ispc::calcforce_on_first_nact(
              ib, ie, ib,
              jb, je, dnj,
              eps2,
              (ispc::Predictor*)&pred_SoA,
              &gacc[0][ty][tx][0],
              &gacc[1][ty][tx][0],
              &gacc[2][ty][tx][0],
              &gjrk[0][ty][tx][0],
              &gjrk[1][ty][tx][0],
              &gjrk[2][ty][tx][0],
              &gsnp[0][ty][tx][0],
              &gsnp[1][ty][tx][0],
              &gsnp[2][ty][tx][0],
              &gcrk[0][ty][tx][0],
              &gcrk[1][ty][tx][0],
              &gcrk[2][ty][tx][0]);

#pragma omp barrier

          for (int i = ty; i < ie-ib; i += nty)
          {
            double ax = 0.0;
            double ay = 0.0;
            double az = 0.0;
            double jx = 0.0;
            double jy = 0.0;
            double jz = 0.0;
            double sx = 0.0;
            double sy = 0.0;
            double sz = 0.0;
            double cx = 0.0;
            double cy = 0.0;
            double cz = 0.0;
            for (int j = 0; j < nty; j++)
            {
              ax += gacc[0][j][tx][i];
              ay += gacc[1][j][tx][i];
              az += gacc[2][j][tx][i];
              jx += gjrk[0][j][tx][i];
              jy += gjrk[1][j][tx][i];
              jz += gjrk[2][j][tx][i];
              sx += gsnp[0][j][tx][i];
              sy += gsnp[1][j][tx][i];
              sz += gsnp[2][j][tx][i];
              cx += gcrk[0][j][tx][i];
              cy += gcrk[1][j][tx][i];
              cz += gcrk[2][j][tx][i];
            }
            force[ib+i].acc = dvec3(ax,ay,az);
            force[ib+i].jrk = dvec3(jx,jy,jz);
            force[ib+i].snp = dvec3(sx,sy,sz);
            force[ib+i].crk = dvec3(cx,cy,cz);
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
        ispc::calcforce_on_first_nact(
            ib, ie, 0,
            0,nj,1,
            eps2,
            (ispc::Predictor*)&pred_SoA,
            &force_SoA.acc[0][0],
            &force_SoA.acc[1][0],
            &force_SoA.acc[2][0],
            &force_SoA.jrk[0][0],
            &force_SoA.jrk[1][0],
            &force_SoA.jrk[2][0],
            &force_SoA.snp[0][0],
            &force_SoA.snp[1][0],
            &force_SoA.snp[2][0],
            &force_SoA.crk[0][0],
            &force_SoA.crk[1][0],
            &force_SoA.crk[2][0]);


#pragma ivdep
        for (int i = ib; i < ie; i++)
        {
          force[i].acc = dvec3(force_SoA.acc[0][i], force_SoA.acc[1][i], force_SoA.acc[2][i]);
          force[i].jrk = dvec3(force_SoA.jrk[0][i], force_SoA.jrk[1][i], force_SoA.jrk[2][i]);
          force[i].snp = dvec3(force_SoA.snp[0][i], force_SoA.snp[1][i], force_SoA.snp[2][i]);
          force[i].crk = dvec3(force_SoA.crk[0][i], force_SoA.crk[1][i], force_SoA.crk[2][i]);
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
         126.0*ni*nj/(t1-t0)/1e9);
    //    assert(0);

#else  /* naive C++ code */

#pragma omp parallel for
    for(int i=0; i<ni; i++){
      const dvec3 posi = pred[i].pos;
      const dvec3 veli = pred[i].vel;
      const dvec3 acci = pred[i].acc;
      const dvec3 jrki = pred[i].jrk;
      dvec3 acc(0.0);
      dvec3 jrk(0.0);
      dvec3 snp(0.0);
      dvec3 crk(0.0);
      for(int j=0; j<nj; j++){
        const dvec3  dr    = pred[j].pos - posi;
        const dvec3  dv    = pred[j].vel - veli;
        const dvec3  da    = pred[j].acc - acci;
        const dvec3  dj    = pred[j].jrk - jrki;

        const double r2    = eps2 + dr*dr;
#if 0
        if(r2 == 0.0) continue;
#endif
        const double drdv  = dr*dv;
        const double dvdv  = dv*dv;
        const double drda  = dr*da;
        const double dvda  = dv*da;
        const double drdj  = dr*dj;

        const double ri2   = 1.0 / r2;
        const double mri3  = pred[j].mass * ri2 * sqrt(ri2);
        const double alpha = drdv * ri2;
        const double beta  = (dvdv + drda)*ri2 + alpha*alpha;
        const double gamma = (3.0*dvda + drdj)*ri2 + alpha*(3.0*beta - 4.0*alpha*alpha);

#if 0
        const dvec3 aij = mri3 * dr;
        const dvec3 jij = mri3 * dv + (-3.0*alpha) * aij;
        const dvec3 sij = mri3 * da + (-6.0*alpha) * jij + (-3.0*beta) * aij;
        const dvec3 cij = mri3 * dj + (-9.0*alpha) * sij + (-9.0*beta) * jij + (-3.0*gamma) * aij;

        acc += aij;
        jrk += jij;
        snp += sij;
        crk += cij;
#else
        acc += mri3 * dr;
        dvec3 tmp1 = dv + (-3.0*alpha) * dr;
        jrk += mri3 * tmp1;
        dvec3 tmp2 = da + (-6.0*alpha) * tmp1 + (-3.0*beta) * dr;
        snp += mri3 * tmp2;
        dvec3 tmp3 = dj + (-9.0*alpha) * tmp2 + (-9.0*beta) * tmp1 + (-3.0*gamma) * dr;
        crk += mri3 * tmp3;
#endif
      }
      force[i].acc = acc;
      force[i].jrk = jrk;
      force[i].snp = snp;
      force[i].crk = crk;
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
        const double mj = (j != i) ? ptcl[j].mass : 0.0;
        pot -= mj / sqrt(r2);
      }
      potbuf[i] = pot;
    }
  }
};

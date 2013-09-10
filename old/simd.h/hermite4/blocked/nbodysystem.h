#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cassert>
#include <algorithm>
#include <sys/time.h>

// #include <mpi.h>

#include "vector3.h"

// #define NOPROFILE

#ifdef NOPROFILE
struct Profile{
	enum{
		FORCE = 0,
		POT,
		PREDICT,
		CORRECT,
		SORT,
		SET_JP,
		IO,
		// don't touch the below
		TOTAL,
		MISC,
		NUM_ELEM,
	};
	void flush(){}
	void beg(const int elem, const bool reuse = false){}
	void end(const int elem){}
	void show(
			FILE *fp = stderr,
			const char *fmt = " %s : %e\n")
	{}

	static double wtime(){
#if 1
		return 0.0;
#else
		struct timeval tv;
		gettimeofday(&tv, NULL);
		return double(tv.tv_sec) + 1.e-6*double(tv.tv_usec);
#endif
	}
};
#else
struct Profile{
	enum{
		FORCE = 0,
		POT,
		PREDICT,
		CORRECT,
		SORT,
		SET_JP,
		IO,
		// don't touch the below
		TOTAL,
		MISC,
		NUM_ELEM,
	};
	static const char *name(const int i){
		static const char *strs[NUM_ELEM] = {
			"force   ",
			"pot     ",
			"predict ",
			"correct ",
			"sort    ",
			"set_jp  ",
			"I/O     ",
			"total   ",
			"misc    ",
		};
		return strs[i];
	}
	double time[NUM_ELEM];
	double tprev;
	
	void flush(){
		for(int i=0; i<NUM_ELEM; i++) time[i] = 0.0;
	}
	void beg(const int elem, const bool reuse = false){
		if(reuse) time[elem] -= tprev;
		else      time[elem] -= wtime();
	}
	void end(const int elem){
		// time[elem] += wtime();
		tprev = wtime();
		time[elem] += tprev;
	}
	void show(
			FILE *fp = stderr,
			const char *fmt = " %s : %e\n")
	{
		time[MISC] = time[TOTAL];
		for(int i=0; i<NUM_ELEM-2; i++){
			time[MISC] -= time[i];
		}
		for(int i=0; i<NUM_ELEM; i++){
			fprintf(fp, fmt, name(i), time[i]);
		}
	}

	static double wtime(){
#ifdef __HPC_ACE__
#  if 1
		return 1.e-6 * __gettod();
#  else
		unsigned long x;
		asm volatile ("rd %%tick, %0" : "=r" (x));
		return (1.0/2.0e9) * (double)x;
#  endif
#else
		struct timeval tv;
		gettimeofday(&tv, NULL);
		return double(tv.tv_sec) + 1.e-6*double(tv.tv_usec);
#endif
	}
};
#endif

struct RAII_Timer{
	const char *name, *fmt;
	FILE *fp;
	double tbeg, tend;
	RAII_Timer(
			const char *_name = "timer",
			FILE       *_fp   = stdout,
			const char *_fmt  = "%s: %e sec\n")
		: name(_name), fmt(_fmt), fp(_fp)
	{
		tbeg = Profile::wtime();
	}
	~RAII_Timer(){
		tend = Profile::wtime();
		fprintf(fp, fmt, name, tend-tbeg);
	}
};

#if 0
template <typename T, size_t align>
T *allocate(size_t num_elem){
	void *ptr;
	int ret = posix_memalign(&ptr, align, num_elem * sizeof(T));
	assert(0 == ret);
	assert(ptr);
	return (T *)ptr;
}
#endif

struct CmpPtcl_dt{
	bool operator()(const Particle &p1, const Particle &p2) const {
		return (p1.dt < p2.dt);
	}
};
struct CmpPtcl_id{
	bool operator()(const Particle &p1, const Particle &p2) const {
		return (p1.id < p2.id);
	}
};

struct Split_dt{
	const double dt;
	Split_dt(const double _dt) : dt(_dt) {}
	bool operator()(const Particle &p){
		return (p.dt < dt);
	}
};

struct NbodySystem{
	// paramters
	long   nbody;
	double eps2;
	double eta, etapow, eta_s;
	double tsys, tend;
	double dtmax;
	// counters
	double init_energy, prev_energy;
	long   num_step, num_bstep;
	long   num_step_tot, num_bstep_tot;
	// buffers
	Particle *ptcl;
	Force    *force;
	double   *potbuf;
	// int      *ilist;
	// objects
	Gravity  *gravity;
	Profile   prof;


	NbodySystem(){
		nbody = 0;
		num_step      = 0;
		num_bstep     = 0;
		num_step_tot  = 0;
		num_bstep_tot = 0;
		dtmax = 0.0;
		eta = etapow = eta_s = 0.0;
		ptcl   = NULL;
		force  = NULL;
		potbuf = NULL;
		// ilist = NULL;
		// schedule = NULL;
		gravity  = NULL;
	}
	~NbodySystem(){
		release_resources();
	}

	void allocate_resources(){
		fprintf(stderr, "allocate resources for n = %ld\n", nbody);
		assert(nbody > 0);

		release_resources();

		ptcl   = allocate<Particle, 64> (nbody);
		force  = allocate<Force,    64> (nbody);
		potbuf = allocate<double,   64> (nbody);
		// ilist = allocate<int,      64> (nbody);

		// schedule = new Schedule(dtmax);
		gravity  = new Gravity(nbody);
	}
	void release_resources(){
		free(ptcl  ); ptcl   = NULL;
		free(force ); force  = NULL;
		free(potbuf); potbuf = NULL;
		// free(ilist); ilist = NULL;

		// delete schedule; schedule = NULL;
		delete gravity;  gravity  = NULL;
	}
#if 0
	void read_snapshot_berczic(const char *filename){
		FILE *fp = fopen(filename, "r");
		assert(fp);
		int snapid;
		int nread = fscanf(fp, "%d %ld %lf", &snapid, &nbody,  &tsys);
		assert(3 == nread);
		fprintf(stderr, "read snapshot, n = %ld\n", nbody);

		allocate_resources();

		for(int i=0; i<nbody; i++){
			int id;
			double mass;
			dvec3 pos, vel;
			nread = fscanf(fp, "%d %lf %lf %lf %lf %lf %lf %lf",
					&id, &mass, &pos.x, &pos.y, &pos.z, &vel.x, &vel.y, &vel.z);
			assert(8 == nread);
			// The constructer flushes higher derivatives
			ptcl[i] = Particle(id, mass, pos, vel);
#ifdef EIGHTHDD
			const dvec3 off(sqrt(3.0) * pow(2.0, 32.0));
			ptcl[i].posL = off;
			dd_twosum(ptcl[i].posH, ptcl[i].posL);
#endif
		}

		fclose(fp);
		fprintf(stderr, "read snapshot, done\n");
	}
#endif
	void write_snapshot_masaki(const char *filename, const int snapid){
		FILE *fp = fopen(filename, "w");
		assert(fp);

		fprintf(fp, "%d\n%ld\n%f\n", snapid, nbody, tsys);
		for(int i=0; i<nbody; i++){
			const Particle &p = ptcl[i];
#ifndef EIGTHDD
			const dvec3 pos = p.pos;
#else
			const dvec3 pos = p.posH + p.posL;
#endif
			fprintf(fp, "%6ld    %A    %A %A %A    %A %A %A    %A\n",
					p.id, p.mass, 
					pos.x, pos.y, pos.z,
					p.vel.x, p.vel.y, p.vel.z,
					potbuf[i]);
		}
		fclose(fp);
	}
	void read_snapshot_masaki(const char *filename, int &snapid){
		int nread;
		FILE *fp = fopen(filename, "r");
		assert(fp);

		nread = fscanf(fp, "%d %ld %lf", &snapid, &nbody, &tsys);
		assert(3 == nread);
		fprintf(stderr, "read snapshot (form M.I.), n = %ld, t = %f\n", nbody, tsys);

		allocate_resources();

		for(int i=0; i<nbody; i++){
			long id;
			double mass, pot;
			dvec3 pos, vel;
			nread = fscanf(fp, "%ld %lA %lA %lA %lA %lA %lA %lA %lA", 
					&id, &mass,
					&pos.x, &pos.y, &pos.z,
					&vel.x, &vel.y, &vel.z,
					&pot);
			assert(9 == nread);
			ptcl[i] = Particle(id, mass, pos, vel);
		}
		fclose(fp);
	}

	void dump(const char *filename, const int dumpid){
		FILE *fp = fopen(filename, "w");
		assert(fp);

#if 0
		sort_ptcl_by_id(nbody);
#endif

		fprintf(fp, "%d %d\n%ld\n%20.16f\n", Particle::order, dumpid, nbody, tsys);
		fprintf(fp, "%A %A\n", init_energy, prev_energy);
		fprintf(fp, "%ld %ld\n", num_step_tot, num_bstep_tot);
		for(int i=0; i<nbody; i++){
			ptcl[i].dump(fp);
		}
		fclose(fp);

#if 0
		sort_ptcl(nbody, dtmax);
		for(int i=0; i<nbody; i++){
			gravity->set_jp(i, ptcl[i]);
		}
#endif
	}
	void restore(const char *filename, int &dumpid){
		FILE *fp = fopen(filename, "r");
		assert(fp);
		int order;
		assert(4 == fscanf(fp, "%d %d %ld %lf", &order, &dumpid, &nbody, &tsys));
		assert(2 == fscanf(fp, "%lA %lA", &init_energy, &prev_energy));
		assert(2 == fscanf(fp, "%ld %ld", &num_step_tot, &num_bstep_tot));
		assert(order == Particle::order);
		fprintf(stderr, "read dump file : <%s>, order=%d, nbody=%ld, tsys=%f\n", 
				filename, order, nbody, tsys);

		allocate_resources();

		for(int i=0; i<nbody; i++){
			ptcl[i].restore(fp);
			ptcl[i].tlast = tsys;
		}
		fclose(fp);
	}
	void reinitialize_steps(const double _eta, const double _dtmax){
		eta   = _eta;
		dtmax = _dtmax;
		for(int i=0; i<nbody; i++){
			ptcl[i].recalc_dt(eta, dtmax);
		}

		sort_ptcl(nbody, dtmax);
		for(int i=0; i<nbody; i++){
			gravity->set_jp(i, ptcl[i]);
		}
	}

	void set_run_params(
			const double _tend,
			const double _dtmax,
			const double _eps2,
			const double _eta,
			const double _eta_s)
	{
		tend  = _tend;
		dtmax = _dtmax;
		eps2  = _eps2;
		eta   = _eta;
		eta_s = _eta_s;

		const int ipow = 2*(Particle::order - 3);
		etapow = pow(eta, double(ipow));

		// if(schedule) delete schedule;
		assert(dtmax > 0.0);
		// schedule = new Schedule(dtmax);
	}
	void set_thermal_eps(const double a){
		const double eps = a / nbody;
		fprintf(stderr, "eps is set to %e\n", eps);
		eps2 = eps * eps;
	}

	void sort_ptcl(const int nact, const double dtlim){
		// this version is faster than the version with std::sort
#if 0
		sort_ptcl(nact);
#else
		Particle *beg = ptcl;
		Particle *end = ptcl+nact;
		double dt = dtlim;
		while(end != beg){
			end = std::partition(beg, end, Split_dt(dt));
			dt *= 0.5;
		}
#endif
	}

	void sort_ptcl(const int nact){
		std::sort(ptcl+0, ptcl+nact, CmpPtcl_dt());
	}

	void sort_ptcl_by_id(const int nact){
		std::sort(ptcl+0, ptcl+nact, CmpPtcl_id());
	}

	void init_step(){
		const int n = nbody;
		for(int i=0; i<n; i++){
			ptcl[i].tlast = tsys;
			gravity->set_jp(i, ptcl[i]);
			// ilist[i] = i;
		}
		init_energy = prev_energy = calc_energy(true);
		if(Particle::order > 4){ 
			puts("INITIAL FORCE");
			gravity->predict_all(tsys);
			gravity->calc_force_on_first_nact(nbody, eps2, force);
			for(int i=0; i<n; i++){
				ptcl[i].assign_force(force[i]);
				gravity->set_jp(i, ptcl[i]);
			}
		}
		prof.beg(Profile::PREDICT);
		gravity->predict_all(tsys);
		prof.end(Profile::PREDICT);

		prof.beg(Profile::FORCE);
		gravity->calc_force_on_first_nact(nbody, eps2, force);
		prof.end(Profile::FORCE);
		for(int i=0; i<n; i++){
			ptcl[i].assign_force(force[i]);
			ptcl[i].init_dt(eta_s, dtmax);
		}
		// sort_ptcl(nbody);
		sort_ptcl(nbody, dtmax);
		for(int i=0; i<n; i++){
#if 0
			printf("%4ld, %A, %f\n", 
					ptcl[i].id, ptcl[i].dt, ptcl[i].acc.abs()/ptcl[i].jrk.abs());
#endif
			gravity->set_jp(i, ptcl[i]);
		}
		init_energy = calc_energy(true);
	}

	double calc_energy(bool print=false){
		prof.beg(Profile::POT);
		gravity->calc_potential(eps2, potbuf);
		prof.end(Profile::POT);
		double ke = 0.0;
		double pe = 0.0;
		for(int i=0; i<nbody; i++){
			ke += ptcl[i].mass * ptcl[i].vel.norm2();
			pe += ptcl[i].mass * potbuf[i];
		}
		ke *= 0.5;
		pe *= 0.5;
		if(print){
			fprintf(stderr, "ke = %24.16e, pe = %24.16e, e = %24.16e\n",
					ke, pe, ke+pe);
		}
		return ke+pe;
	}

	double calc_dtlim(const double tnext) const {
		double dtlim = dtmax;
		double s = tnext / dtmax;
		while(s != double(int(s))){
			s *= 2.0;
			dtlim *= 0.5;
			assert(dtlim >= 1.0/(1LL<<32));
		}
		return dtlim;
	}

	int count_nact(const double tnext) const {
		int nact;
		for(nact = 0; nact<nbody; nact++){
			const Particle &p = ptcl[nact];
			if(p.dt + p.tlast != tnext) break;
		}
		return nact;
	}

	__attribute__((noinline))
	void integrate_one_block(){
		const double tnext = ptcl[0].tlast + ptcl[0].dt;
		const double dtlim = calc_dtlim(tnext);
		const int    nact  = count_nact(tnext);
#if 0
		printf("t = %f, nact = %6d, dtlim = %A\n", tsys, nact, dtlim);
#endif
	  prof.beg(Profile::PREDICT);
		gravity->predict_all(tnext);
	  prof.end(Profile::PREDICT);

	  prof.beg(Profile::FORCE, true);
		gravity->calc_force_on_first_nact(nact, eps2, force);
	  prof.end(Profile::FORCE);

	  prof.beg(Profile::CORRECT);
#pragma omp parallel for
		for(int i=0; i<nact; i++){
			ptcl[i].correct(force[i], eta, etapow, dtlim);
		}
	  prof.end(Profile::CORRECT);

	  prof.beg(Profile::SORT, true);
		// sort_ptcl(nact);
		sort_ptcl(nact, dtlim);
	  prof.end(Profile::SORT);

	  prof.beg(Profile::SET_JP, true);
		for(int i=0; i<nact; i++){
			gravity->set_jp(i, ptcl[i]);
		}
	  prof.end(Profile::SET_JP);

		num_step += nact;
		num_bstep++;

		tsys = tnext;
	}

	void integrate_one_dtmax(){
		const double tt = tsys + dtmax;
		while(tsys < tt){
			integrate_one_block();
		}
		assert(tsys == tt);
	}

	void integrate(const double tcrit){
		prof.flush();
		prof.beg(Profile::TOTAL);
        while(tsys < tcrit){
			double t0 = Profile::wtime();
			integrate_one_dtmax();
			double t1 = Profile::wtime();
			print_stat(t1-t0);
		}
		prof.end(Profile::TOTAL);
		prof.show();
	}

	void print_stat(const double wtime, FILE *fp = stdout){
		const double energy = calc_energy(false);
		const double de_glo =((init_energy - energy) / init_energy);
		const double de_loc =((prev_energy - energy) / init_energy);
		num_step_tot  += num_step;
		num_bstep_tot += num_bstep;
		const double nact_avr = double(num_step) / double(num_bstep);
		const double Gflops = (((1.0e-9 * num_step) * nbody) * Particle::flops) / wtime;
		prof.beg(Profile::IO);
#if 0
		fprintf(fp, "t = %f\n", tsys);
		fprintf(fp, " steps: %ld %ld %ld %ld\n", num_bstep, num_step, num_bstep_tot, num_step_tot);
		fprintf(fp, " nact : %f\n", nact_avr);
		fprintf(fp, " de(local/global) : %+e %+e\n", de_loc, de_glo);
		fprintf(fp, " %f sec, %f Gflops\n", wtime, Gflops);
#else
		fprintf(stderr, "t = %f\n", tsys);
		fprintf(stderr, " steps: %ld %ld %ld %ld\n", num_bstep, num_step, num_bstep_tot, num_step_tot);
		fprintf(stderr, " nact : %f\n", nact_avr);
		fprintf(stderr, " de(local/global) : %+e %+e\n", de_loc, de_glo);
		fprintf(stderr, " %f sec, %f Gflops\n", wtime, Gflops);

        fprintf(fp, "%f %6ld %6ld %10ld %10ld  %8.2f  %+e  %+e  %f  %f \n", 
                tsys, 
                num_bstep, num_step, 
                num_bstep_tot, num_step_tot, 
                nact_avr, 
                de_loc, de_glo, 
                wtime, Gflops);
		prof.end(Profile::IO);
#endif 
		prev_energy = energy;
		num_step  = 0;
		num_bstep = 0;
	}
};


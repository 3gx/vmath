#include <cstdio>
#include <cassert>
#include <cassert>
#include "vector3.h"

#ifdef __SSE__
#warning SSE is available
typedef double v2df __attribute__((vector_size(16)));
#endif
#ifdef __AVX__
#warning AVX is available
typedef double v4df __attribute__((vector_size(32)));
#endif
#ifdef __HPC_ACE__
#warning HPC-ACE is available
#include <fjcex.h>
#endif

template <typename T, size_t align>
T *allocate(size_t num_elem){
	void *ptr;
	int ret = posix_memalign(&ptr, align, num_elem * sizeof(T));
	assert(0 == ret);
	assert(ptr);
	return (T *)ptr;
}

#if   defined FOURTH
#  include "hermite4.h"
#elif defined SIXTH
#  include "hermite6.h"
#elif defined EIGHTH
#  include "hermite8.h"
#elif defined EIGHTHDD
#  include "hermite8dd.h"
#else
#  error
#endif

#include "nbodysystem.h"

#if 0
int main(){
	NbodySystem sys;
	const double tend = 1.0;
    //const double dtmax = 1.0/16.0;
	sys.read_snapshot_berczic("data.inp");
#if   defined FOURTH
	const double eta   = 0.1;
	const double eta_s = 0.01;
#elif defined SIXTH
	const double eta   = 0.5;
	const double eta_s = 0.1;
#elif defined EIGHTH
	const double eta   = 0.75;
	const double eta_s = 0.1;
#elif defined EIGHTHDD
	const double eta   = 0.75;
	const double eta_s = 0.1;
#endif
	sys.set_run_params(tend, 1./16., 0.0, eta, eta_s);

#if defined EIGHTHDD
	const double kt = 4.0;
#else
	const double kt = 4.0;
#endif
	sys.set_thermal_eps(kt);

	sys.init_step();
	// sys.integrate(0.5);
	// sys.dump("dump.dat");
	sys.integrate(1.0);

	if(0){
		puts("");
		NbodySystem sys2;
		sys2.restore("dump.dat");
		sys2.set_run_params(tend, 1./16., 0.0, eta, eta_s);
		sys2.set_thermal_eps(kt);
		sys2.reinitialize_steps(eta, 1./16.);
		sys2.integrate(1.0);
	}

#if defined EIGHTHDD
	for(int i=0; i<10; i++){
		const Particle &p = sys.ptcl[i];
		printf("%ld: \n", p.id);
		printf(" x: %+A + %+A\n", p.posH.x, p.posL.x);
		printf(" y: %+A + %+A\n", p.posH.y, p.posL.y);
		printf(" z: %+A + %+A\n", p.posH.z, p.posL.z);
	}
#endif

	// sys.write_snapshot_masaki("masaki.out");
	// sys.read_snapshot_masaki("masaki.out");

	return 0;
}
#else
#include "parameter.h"
namespace std{
	extern istream cin;
	extern ostream cout;
}
int main(){
	Parameter param;
	param.read(std::cin);
	param.print(std::cout);
	param.assert_check(Particle::order);

	int snapid = param.snapid;
	int dumpid = param.dumpid;


	NbodySystem sys;
	const double eps2 = param.eps * param.eps;
	sys.set_run_params(param.tend, param.dtmax, eps2, param.eta, param.eta_s);
	if(param.snapin.length() > 0){
		sys.read_snapshot_masaki(param.snapin.c_str(), snapid);
		if(param.kt_for_eps >= 0.0){
			sys.set_thermal_eps(param.kt_for_eps);
		}
		sys.init_step();
	}
	else if(param.dumpin.length() > 0){
		sys.restore(param.dumpin.c_str(), dumpid);
		if(param.kt_for_eps >= 0.0){
			sys.set_thermal_eps(param.kt_for_eps);
		}
		sys.reinitialize_steps(param.eta, param.dtmax);
	}else{
		assert(!"xxx");
	}

	// sys.integrate(param.tend);

	double tnext_log  = sys.tsys + param.log_interval;
	double tnext_snap = sys.tsys + param.snap_interval;
	double tnext_dump = sys.tsys + param.dump_interval;

	double wt_log, wt_snap, wt_dump;
	wt_log = wt_snap = wt_dump = sys.prof.wtime();

	sys.prof.flush();
	sys.prof.beg(Profile::TOTAL);
	while(sys.tsys < param.tend){
		sys.integrate_one_dtmax();

		if(sys.tsys >= tnext_log){
			tnext_log += param.log_interval;
			double wt = sys.prof.wtime();
			sys.print_stat(wt - wt_log);
			wt_log = wt;
		}

		if(sys.tsys >= tnext_snap){
			tnext_snap += param.snap_interval;
			snapid++;
			// fprintf(stderr, "output snapshot [%d] (dummy)\n", snapid);
			if(param.snapout_base.length() > 0){
				static char fname[256];
				snprintf(fname, 256, param.snapout_base.c_str(), snapid);
				fprintf(stderr, "snap file : %s\n", fname);
				sys.write_snapshot_masaki(fname, snapid);
			}
		}

		if(sys.tsys >= tnext_dump){
			tnext_dump += param.dump_interval;
			dumpid++;

			sys.prof.end(Profile::TOTAL);
			sys.prof.show();

			// fprintf(stderr, "output dumpfile [%d] (dummy)\n", dumpid);
			if(param.dumpout_base.length() > 0){
				static char fname[256];
				snprintf(fname, 256, param.dumpout_base.c_str(), dumpid);
				fprintf(stderr, "dump file : %s\n", fname);
				sys.dump(fname, dumpid);
			}

			sys.prof.flush();
			sys.prof.beg(Profile::TOTAL);
		}
	}
}
#endif

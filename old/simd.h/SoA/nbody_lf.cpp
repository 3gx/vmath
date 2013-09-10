#include <cstdio>
#include <iostream>
#include "SSE/sse.h"
#include "particle.h"

typedef double real;
typedef Vector3<real> vec3s;

template<typename T>
inline T min(const T a, const T b) {return std::min(a,b); }// a < b ? a : b; }

struct NBodyLF
{
  ParticleData<real> ptcl;
};

template<int start, int end>
struct static_for
{
  template<typename Func>
  static void eval(Func func)
  {
    func(start);
    static_for<start+1, end>::template eval<Func>(func);
  }
};

template<int end>
struct static_for<end, end>
{
  template<typename Func>
  static void eval(Func func) {}
};

int main(int argc, char * argv[])
{
  int n = 1024;
  NBodyLF s;  
  s.ptcl.resize(n);
  for (int i = 0; i < n; i++)
  {
    vec3s pos(drand48(), drand48(), drand48());
    vec3s vel(drand48(), drand48(), drand48());
    real mass = 1.0/n;

    s.ptcl[i].pos  = pos;
    s.ptcl[i].vel  = vel;
    s.ptcl[i].mass = mass;
  };

  real mTotal = 0.0;
  for (int i = 0; i < n; i++)
    mTotal += s.ptcl[i].pos.x;
  fprintf(stderr, "mTotal= %g\n", mTotal);

  vec3s myvec(argc,2,3);
  fprintf(stderr, "myvec= %g %g %g \n", myvec[0], myvec[1], myvec[2]);
  asm("#begin-A");
  static_for<0,2>::eval([&] (const int i) 
      {
      myvec[i] += myvec[i]*myvec[min(i,2)];
#if 0
      switch(i)
      {
      case 0: asm("#case0"); myvec[0] += myvec[0]*myvec[0]; break;
      case 1: asm("#case1"); break;
      case 2: asm("#case2"); break;
      case 3: asm("#case3"); break;
      default: assert(0);
      };
#endif
//      fprintf(stderr, " i= %d : v= %g  n= %d\n", i, myvec[i], n);
      } );

  asm("#end-A");
  asm("#begin-B");
#define OP(i) myvec[i] += myvec[i]*myvec[min(i+1,2)]
  OP(0);
  OP(1);
  OP(2);
  asm("#end-B");
  
  fprintf(stderr, ">>mTotal= %g\n", mTotal);
  fprintf(stderr, "myvec= %g %g %g \n", myvec[0], myvec[1], myvec[2]);

  return 0;
};

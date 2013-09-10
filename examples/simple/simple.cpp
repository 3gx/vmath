#include <cstdio>
#include <cstdlib>
#include <vmath.h>

// Include the header file that the ispc compiler generates
//
//

struct simple
{
  simple(const double *vin, double *vout, const int count) 
  {
    mainloop m(vin,vout);
    foreach::loop(0,count,m);
  }

  struct mainloop
  {
    mainloop(const double *_vin, double *_vout) : vin(_vin), vout(_vout) {}
    const double *vin;
    double *vout;

    VECTORLOOPFUNCTION(vreal,vint,vmask,index)
    {
      printf("index= %d sizof(vreal)= %d\n", index, sizeof(vreal));
      vreal v = vload<vreal>(&vin[index]);

      // Do an arbitrary little computation, but at least make the
      // computation dependent on the value being processed
      cif (v < 9.0,
          IF  (_mask) { _mask(v, v*v); },
          ELSE(_mask) { _mask(v,__sqrt(v)); });

      // And write the result to the output array.
      vstore(&vout[index],v);
    }
  };
};



int main(int argc, char * argv[])
{
  double vin[23], vout[23];

  // Initialize input buffer
  for (int i = 0; i < 23; ++i)
    vin[i] = (double)i;

  // Call simple() function from simple.ispc file
  simple(vin, vout, 23);

  // Print results
  for (int i = 0; i < 23; ++i)
    printf("%2d: simple(%f) = %f\n", i, vin[i], vout[i]);

  return 0;
}

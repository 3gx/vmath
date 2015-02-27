/* Author: Intel */
/* vmath port: Evghenii Gaburov */

#include <cstdio>
#include <algorithm>
#include <rtc.h>
#include <vmath.h>

struct mandelbrot
{
  const double x0, y0;
  const double x1, y1;
  const int width, height;
  const int maxIterations;
  double *output;
  double dx,dy;
  int j;

  mandelbrot(
      const double _x0, const double _y0,
      const double _x1, const double _y1,
      const int _width, const int _height,
      const int _maxIterations,
      double *_output) :
    x0(_x0), y0(_y0),
    x1(_x1), y1(_y1),
    width(_width), height(_height),
    maxIterations(_maxIterations),
    output(_output)
  {
    const double dx = (x1 - x0) / width;
    const double dy = (y1 - y0) / height;

    for (int j = 0; j < height; j++)
    {
      vinteger vj (j);
      vinteger vi = vspan<vinteger>(0);
      const int nvec = vdouble::round(width);
      for (int i = 0; i < nvec; i += vdouble::SIZE, vi += vdouble::SIZE)
      {
        const vdouble x = x0 + vdouble(vi) * dx;
        const vdouble y = y0 + vdouble(vj) * dy;
        const int index = j * width + i;
        const vdouble pts = mandel<vdouble,vinteger,vboolean>(x,y,maxIterations);
        vstore(&output[index],pts);
      }
      for (int i = nvec; i < width; i++)
      {
        const double x = x0 + double(i) * dx;
        const double y = y0 + double(j) * dy;
        const int index = j * width + i;
        const double pts = mandel<double,int,bool>(x,y,maxIterations);
        vstore(&output[index],pts);
      }
    }
  }
  template<typename vreal, typename vint, typename vmask>
    vreal mandel(const vreal c_re, const vreal c_im, const int count) 
    {
      vreal z_re = c_re;
      vreal z_im = c_im;
      vreal counter(0.0);

      vmask break_flag = vfalse<vmask>();
      for (int i = 0; i < count; ++i) 
      {
        break_flag = break_flag | __square(z_re) + __square(z_im) > 4.0;
        if (all(break_flag)) break;

        cif (!break_flag, 
            IF(_mask)
            {
              const vreal new_re = z_re*z_re - z_im*z_im;
              const vreal new_im = 2.0 *z_re * z_im;
              _mask(z_re , c_re + new_re);
              _mask(z_im , c_im + new_im);
              _mask(counter, counter + 1.0);
            });
      }

      return counter;
    }

};

extern void mandelbrot_serial(double x0, double y0, double x1, double y1,
    int width, int height, int maxIterations,
    double output[]);

/* Write a PPM image file with the image of the Mandelbrot set */
static void
writePPM(int *buf, int width, int height, const char *fn) {
  FILE *fp = fopen(fn, "wb");
  fprintf(fp, "P6\n");
  fprintf(fp, "%d %d\n", width, height);
  fprintf(fp, "255\n");
  for (int i = 0; i < width*height; ++i) {
    // Map the iteration count to colors by just alternating between
    // two greys.
    char c = (buf[i] & 0x1) ? 240 : 20;
    for (int j = 0; j < 3; ++j)
      fputc(c, fp);
  }
  fclose(fp);
  printf("Wrote image file %s\n", fn);
}


int main() {
  unsigned int width = 768;
  unsigned int height = 512;
  float x0 = -2;
  float x1 = 1;
  float y0 = -1;
  float y1 = 1;

  int maxIterations = 256;
  double *buf = new double[width*height];

  //
  // Compute the image using the ispc implementation; report the minimum
  // time of three runs.
  //
  double minISPC = 1e30;
  for (int i = 0; i < 3; ++i) {
    const double t0 = rtc();
    mandelbrot(x0, y0, x1, y1, width, height, maxIterations, buf);
    const double dt = rtc() - t0;
    minISPC = std::min(minISPC, dt);
  }

  printf("[mandelbrot ispc]:\t\t[%.3f] million cycles\n", minISPC);
  int *ibuf = new int[width*height];
  for (int i = 0; i < width*height;i++) ibuf[i] = (int)buf[i];
  writePPM(ibuf, width, height, "mandelbrot-ispc.ppm");

  // Clear out the buffer
  for (unsigned int i = 0; i < width * height; ++i)
    buf[i] = 0;

  // 
  // And run the serial implementation 3 times, again reporting the
  // minimum time.
  //
  double minSerial = 1e30;
  for (int i = 0; i < 3; ++i) {
    const double t0 = rtc();
    mandelbrot_serial(x0, y0, x1, y1, width, height, maxIterations, buf);
    const double dt = rtc() - t0;
    minSerial = std::min(minSerial, dt);
  }

  printf("[mandelbrot serial]:\t\t[%.3f] million cycles\n", minSerial);
  for (int i = 0; i < width*height;i++) ibuf[i] = (int)buf[i];
  writePPM(ibuf, width, height, "mandelbrot-serial.ppm");

  printf("\t\t\t\t(%.2fx speedup from ISPC)\n", minSerial/minISPC);

  return 0;
}

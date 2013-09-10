#pragma once

#include <cstdio>
#include <cstdlib>
#include <cassert>

template<typename T, const int NX, const int NY>
struct FArray2D
{
  private:
    int _xmin;
    int _ymin;
    T __restrict__ (*_data)[NX];

  public:

    FArray2D(
        const T  *data = NULL,
        const int xmin = 0,
        const int ymin = 0) : _xmin(xmin), _ymin(ymin)
  {
      set_data(data);
  }

    void set_data(const T *data)
    {
      _data = (T (*)[NX])(data - _ymin*NX - _xmin);
    }

    const T& operator()(const int i, const int j) const { return _data[j][i]; }
    T& operator()(const int i, const int j) { return _data[j][i]; }
    
    operator T*() {return &_data[0][0];}

    int xmin() const { return  _xmin; }
    int ymin() const { return  _ymin; }
};

template<typename T, const int NX>
struct FArray1D
{
  private:
    int _xmin;
    T *_data;

  public:

    FArray1D(
        const T  *data = NULL,
        const int xmin = 0) : _xmin(xmin)
  {
    set_data(data);
  }

    void set_data(const T *data)
    {
      _data = (T*)(data - _xmin);
    }

    const T& operator()(const int i) const { return _data[i]; }
    T& operator()(const int i) { return _data[i]; }

    operator T*() {return _data;}

    int xmin() const { return _xmin; }
};


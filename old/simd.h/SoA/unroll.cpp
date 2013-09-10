#include <iostream>
#include <algorithm>
#include <cassert>


template<int start, int end, int step = 1>
struct static_for
{
  template<typename Func>
    static void eval(Func func)
    {
      func(start);
      static_for<start+step, end>::template eval<Func>(func);
    }
};

template<int end>
struct static_for<end, end>
{
  template<typename Func>
    static void eval(Func func) {}
};



int main(int argc, char * argv[] )
{
  int n = argc;

  constexpr int UNROLL = 8;
  for (int i = 0; i < n; i + UNROLL)
    static_unroll<UNROLL>([&] (const int i) 
        {
        fprintf(stderr ,"i= %d n= %d\n", i, n);
        }
        );

  return 0;
};

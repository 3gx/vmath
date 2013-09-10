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

  template<const int N>
void print()
{
  fprintf(stderr, " number= %d \n", N);
}



int main(int argc, char * argv[] )
{
  int n = argc;

  static_for<0,8,2>::eval([&] (int i) 
      {
#if 0
    print<i>();
#else
    fprintf(stderr, " i= %d  n= %d\n", i, n);
#endif
      } );

#if 0 
  print<10>();
  print1(10);
#else
  print<10>();
#endif
  
  return 0;
};

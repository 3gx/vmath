#include <iostream>
#include <algorithm>
#include <cassert>

template<int N>
struct static_unroll
{
  template<typename Lambda>
    static void eval(const int i, const Lambda& f)
    {
      static_unroll<N-1>::eval(i, f);
      f(i+N);
    }
};

template<>
struct static_unroll<0>
{
  template<typename Lambda>
    static void eval(const int i, const Lambda& f)
    {
      f(i);
    }
};

template<int UNROLL>
struct dynamic_for
{
  template<typename Lambda>
    static void eval(const int beg, const int end, const Lambda &f)
    {
      for (int i = beg; i < end; i += UNROLL)
        static_unroll<UNROLL>::eval(i, f);
    }
};


int main(int argc, char * argv[] )
{
  assert(argc > 1);
  int n = atoi(argv[1]);

  constexpr int UNROLL = 8;
  for (int i = 0; i < n; i += UNROLL)
    static_unroll<UNROLL>::eval(i, [&] (const int i) 
        {
        fprintf(stderr ," static_unroll: i= %d n= %d\n", i, n);
        }
        );

  dynamic_for<UNROLL>::eval(0, n, [&] (const int i)
      {
      fprintf(stderr ," dynmic_for: i= %d n= %d\n", i, n);
      });
  
};

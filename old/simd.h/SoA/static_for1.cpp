#include <array>
#include <iostream>

template<int N>
inline void print()
{
};

/******** static for **********/

template<int end, int step, int N, int... Ts>
struct static_for_unroll
{
  using type = typename static_for_unroll<end, step, N+step, N, Ts...>::type;
};

template<int end, int step, int... Ts>
struct static_for_unroll<end, step, end, Ts...>
{
  using type = static_for_unroll<end, step, end, Ts...>;
  static void eval()
  {
    constexpr int n = sizeof...(Ts);
    constexpr std::array<int, sizeof...(Ts)> data({{Ts...}});
    fprintf(stderr, " size= %d\n", n);
    for (auto i : data)
      fprintf(stderr, "i= %d \n", i);
  };
};

template<int beg, int end, int step = 1>
struct static_for
{
  static void eval()
  {
    using type = typename static_for_unroll<end, step, beg>::type;
    type::eval();
  };
};


/******* static unroll ********/


int main(int argc , char * argv[])
{
  static_for<0,10,1>::eval();
  static_for<0,10,2>::eval();
  static_for<0,12,3>::eval();
  return 0;
};



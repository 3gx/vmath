#include <array>
#include <iostream>

constexpr int f(int i) { return 2 * i; }

template<int I>
void print()
{
  fprintf(stderr, " >>> print<I>= %d \n", I);
}

void print1(const int I)
{
#if 0
  fprintf(stderr, " >>> print1(I)= %d \n", I);
#else
  print<I>();
#endif

}

template<int N, int... NS>
struct static_forX
{
  static void eval()
  {
    fprintf(stderr, "N= %d  sizeof(NS...)= %d\n", N, sizeof...(NS));
#if 0
    print<f(N)>();
#else
    print1(f(N));
#endif
    static_forX<NS...>::eval();
  }
};

template<>
struct static_forX<8>
{
  static void eval()
  {
    fprintf(stderr, " -- specialization -- \n");
  }
};


#if 0
  template<int N, int... NS>
inline void print()
{
  fprintf(stderr, " number= %d \n", N);
  print<NS...>();
}
#endif

#if 0
template<int... NS>
inline void print<0>() {}
#endif



template <int M, int N, int... Ts>
struct t { 
  using type = typename t<M, N - 1, Ts..., M+1 - N>::type; 
  t() {fprintf(stderr, " val= %d\n", M);}
};

template <int M, int... Ts>
struct t<M, 0u, Ts...>
{
  using type = t<M, 0u, Ts...>;
#if 0
  static std::array<int, sizeof...(Ts)> apply() { return {{f(Ts)...}}; }
#else
  static void apply() 
  {
    constexpr int n = sizeof...(Ts);
    static_forX<Ts...>::eval();
//    print<Ts...>();
//    constexpr std::array<int, sizeof...(Ts)> data({{Ts...}});
 //   print<data[0]>();
    #if 0
    for (int i = 0; i < n; i++)
      fprintf(stderr, "n= %d\n", data[i]);
#endif
  };
#endif
};


int main()
{
  using v = typename t<8, 8>::type;
  v::apply();
#if 0
  fprintf(stderr, "size= %d\n", x.size());
#if 0
  for (int i = 0; i < x.size(); i++)
  {
    fprintf(stderr, "x[%2d]= %d\n", i, x[i]);
  }
#endif
#endif
}

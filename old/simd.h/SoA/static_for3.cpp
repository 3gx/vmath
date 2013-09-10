#include <array>
#include <iostream>
#include <functional>
#include <vector>

void execute(const std::vector<std::function<void()>> &fs)
{
  for (auto& f: fs)
    f();
};

void plain_old_func()
{
  std::cerr << "I'm old plain function\n";
};

struct functor
{
  static void apply()
  {
    std::cerr << "I'm a functor \n";
  };
};

template<typename T>
void printX(T x)
{
  fprintf(stderr, "test= %d\n", x);
}

template<int N>
void print()
{
  fprintf(stderr, "test= %d\n", N);
};

template<int end, int step, int N, int... Ts>
struct static_for_unroll
{
  using type = typename static_for_unroll<end, step, N+step, N, Ts...>::type;
};

template<int end, int step, int... Ts>
struct static_for_unroll<end, step, end, Ts...>
{
  using type = static_for_unroll<end, step, end, Ts...>;
  static constexpr std::array<int, sizeof...(Ts)> eval()
  {
    return std::array<int, sizeof...(Ts)>({{Ts...}});
  };
};

template<int beg, int end, int step = 1>
struct static_for
{
  using type = typename static_for_unroll<end, step, beg>::type;
};

template<typename T, T t, int beg = 0>
struct print_loop
{
};

void do_test(const int n)
{
  constexpr auto tst = static_for<0,20>::type::eval();
  constexpr int nn = tst.size();
//  static_for<tst[0], tst[nn-1]> t;
  print<tst[(nn>>1)]>();
//  print<tst[n]>();

}


int main(int argc, char* argv[] )
{
  std::vector<std::function<void()>> x;
  x.push_back(plain_old_func);
  x.push_back(functor::apply);
  constexpr int i = 8;
  x.push_back([&] ()
      {
      std::cerr << "I'm a lambda, size= " <<  argc << std::endl;
      print<i>();
      });

  execute(x);

  printX<int>(9);


#if 0
  struct t
  {
    template<int N>
      static void doe()
      {
        fprintf(stderr,  "doe= %d\n", N);
      };
  };
#endif
  
  do_test(10);



//  contstexpr std::array<std::function<void()>> funcx({print<0>, print<1>, print<2>});
  return 0;
};

#include <iostream>
#include <algorithm>



  template <typename T, T t>
T  tfunc()
{
  return t + 10;
}

  template <typename T, T t>
constexpr T  func()
{
  return tfunc<T, t>();
}

#define FUNC(a)  func<decltype(a),a>()


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

/*********/

#if 0
template<int start, int end>
struct static_for1
{
  template<typename Func>
    static void eval(Func<int> func)
    {
      func<start>();
      static_for1<start+1, end>::template eval<Func>(func);
    }
};

template<int end>
struct static_for1<end, end>
{
  template<typename Func>
    static void eval(Func<int> func) {}
};
#endif

/*********/

  template<const int n>
void print()
{
  fprintf(stderr, " number= %d \n", n);
}

#if 0
template<const int N>
void loop()
{
  for (int i = 0; i < N; i++)
    print<i>();
}
#endif


#if 0
void print1(constexpr int N)
{
  print<N>();
}
#endif

template <typename F>
struct Foo
{
  static void eval()
  {

  };
};

int main(int argc, char * argv[] )
{
  int n = argc;

  static_for<0,4>::eval([&] (int i) 
      {
  //    print<i>();
      fprintf(stderr, " i= %d  n= %d\n", i, n);
      } );

#if 0 
  print<10>();
  print1(10);
#else
  print<10>();
#endif

  std::cout << FUNC(10) << std::endl;


  return 0;
};

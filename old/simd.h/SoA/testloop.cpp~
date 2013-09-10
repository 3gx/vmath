#include <cstdio>

// Loop<> implementation

// primary template
template < int start, int end, typename T >
struct Loop
{
    static void eval()     {
        const T v(start);                    // do work
        Loop<start+1,end, T>::eval() ;   // increment argument and 				         // recurse
    }
};
// specialization
template <int end, typename T >
struct Loop <end, end, T>
{
    static void eval()    
    {
    }
};
  



struct Test 
{

  template<int N>
    struct T1
    {
    };

  Test(const int N) {fprintf(stderr, "  cnt = %d\n", N);}
}; 

int main(int argc, char * argv[])
{

  Loop<0,10,Test>::eval();


  return 0;
}

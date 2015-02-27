#include <simd.h>
#include <cmath>
#include <array>
#include <gtest/gtest.h>

template<typename T>
class SimdMathTest : public testing::Test
{
  protected:
    using stype = T;
    using vtype = Simd<T>;
    using vref  = SimdRefT<T>;
    using iref  = SimdIRefT<T>;

    std::array<stype,vtype::VLEN> rhs, lhs, res;
    vref vrhs, vlhs, vres;

    SimdMathTest() : vrhs(rhs[0]), vlhs(lhs[0]), vres(res[0])  {}
    virtual ~SimdMathTest() {}

    virtual void SetUp()
    {
      for (int i = 0; i < vtype::VLEN; i++)
      {
        rhs[i] = static_cast<stype>(ceil(drand48()*(1<<16)));
        lhs[i] = static_cast<stype>(ceil(drand48()*(1<<16)));
      }
    }

    virtual void TearDown()
    {
    }

    void test_add()
    {
      vres = vrhs + vlhs;
      for (int i = 0; i < vtype::VLEN; i++)
        ASSERT_EQ(res[i], rhs[i] + lhs[i]);
    }

    void test_sub()
    {
      vres = vrhs - vlhs;
      for (int i = 0; i < vtype::VLEN; i++)
        ASSERT_EQ(res[i], rhs[i] - lhs[i]);
    }
    
    void test_mul()
    {
      vres = vrhs * vlhs;
      for (int i = 0; i < vtype::VLEN; i++)
        ASSERT_EQ(res[i], rhs[i] * lhs[i]);
    }

    void test_div()
    {
      vres = vrhs / vlhs;
      for (int i = 0; i < vtype::VLEN; i++)
        ASSERT_EQ(res[i], rhs[i] / lhs[i]);
    }

};

namespace 
{
  using SimdMathDouble= SimdMathTest<double>;
  TEST_F(SimdMathDouble, ADD)
  {
    test_add();
  }
  TEST_F(SimdMathDouble, SUB)
  {
    test_sub();
  }
  TEST_F(SimdMathDouble, MUL)
  {
    test_mul();
  }
  TEST_F(SimdMathDouble, DIV)
  {
    test_div();
  }
}

namespace 
{
  using SimdMathFloat= SimdMathTest<float>;
  TEST_F(SimdMathFloat, ADD)
  {
    test_add();
  }
  TEST_F(SimdMathFloat, SUB)
  {
    test_sub();
  }
  TEST_F(SimdMathFloat, MUL)
  {
    test_mul();
  }
  TEST_F(SimdMathFloat, DIV)
  {
    test_div();
  }
}

namespace 
{
  using SimdMathInt= SimdMathTest<int>;
  TEST_F(SimdMathInt, ADD)
  {
    test_add();
  }
  TEST_F(SimdMathInt, SUB)
  {
    test_sub();
  }
  TEST_F(SimdMathInt, MUL)
  {
    test_mul();
  }
  TEST_F(SimdMathInt, DIV)
  {
    test_div();
  }
}




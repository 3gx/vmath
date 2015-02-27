#include <simd.h>
#include <cmath>
#include <array>
#include <gtest/gtest.h>

template<typename T>
class SimdTest : public testing::Test
{
  protected:
    using stype = T;
    using vtype = Simd<T>;
    using vref  = SimdRefT<T>;
    using iref  = SimdIRefT<T>;

    std::array<stype,vtype::VLEN> rhs, lhs, res;
    vref vrhs, vlhs, vres;

    SimdTest() : vrhs(rhs[0]), vlhs(lhs[0]), vres(res[0])  {}
    virtual ~SimdTest() {}

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
  using SimdDouble= SimdTest<double>;
  TEST_F(SimdDouble, ADD)
  {
    test_add();
  }
  TEST_F(SimdDouble, SUB)
  {
    test_sub();
  }
  TEST_F(SimdDouble, MUL)
  {
    test_mul();
  }
  TEST_F(SimdDouble, DIV)
  {
    test_div();
  }
}

namespace 
{
  using SimdFloat= SimdTest<float>;
  TEST_F(SimdFloat, ADD)
  {
    test_add();
  }
  TEST_F(SimdFloat, SUB)
  {
    test_sub();
  }
  TEST_F(SimdFloat, MUL)
  {
    test_mul();
  }
  TEST_F(SimdFloat, DIV)
  {
    test_div();
  }
}

namespace 
{
  using SimdInt= SimdTest<int>;
  TEST_F(SimdInt, ADD)
  {
    test_add();
  }
  TEST_F(SimdInt, SUB)
  {
    test_sub();
  }
  TEST_F(SimdInt, MUL)
  {
    test_mul();
  }
  TEST_F(SimdInt, DIV)
  {
    test_div();
  }
}




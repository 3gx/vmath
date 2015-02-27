#include <simd.h>
#include <cmath>
#include <array>
#include <gtest/gtest.h>

template<typename T>
class SimdTest : public testing::Test
{
  protected:
    using simd = Simd<T>;
    using vref = SimdRefT<T>;

    std::array<T,simd::VLEN> rhs, lhs, res;
    vref vrhs, vlhs, vres;

    SimdTest() : vrhs(rhs[0]), vlhs(lhs[0]), vres(res[0])  {}
    virtual ~SimdTest() {}

    virtual void SetUp()
    {
      for (int i = 0; i < simd::VLEN; i++)
      {
        rhs[i] = static_cast<T>(ceil(drand48()*(1<<16)));
        lhs[i] = static_cast<T>(ceil(drand48()*(1<<16)));
      }
    }

    virtual void TearDown()
    {
    }

};

namespace 
{
  using SimdDouble= SimdTest<double>;
  TEST_F(SimdDouble, ADD)
  {
    vres = vrhs + vlhs;
    for (int i = 0; i < simd::VLEN; i++)
      ASSERT_EQ(res[i], rhs[i] + lhs[i]);
  }
  TEST_F(SimdDouble, SUB)
  {
    vres = vrhs - vlhs;
    for (int i = 0; i < simd::VLEN; i++)
      ASSERT_EQ(res[i], rhs[i] - lhs[i]);
  }
  TEST_F(SimdDouble, MUL)
  {
    vres = vrhs * vlhs;
    for (int i = 0; i < simd::VLEN; i++)
      ASSERT_EQ(res[i], rhs[i] * lhs[i]);
  }
  TEST_F(SimdDouble, DIV)
  {
    vres = vrhs / vlhs;
    for (int i = 0; i < simd::VLEN; i++)
      ASSERT_EQ(res[i], rhs[i] / lhs[i]);
  }
}

namespace 
{
  using SimdFloat= SimdTest<float>;
  TEST_F(SimdFloat, ADD)
  {
    vres = vrhs + vlhs;
    for (int i = 0; i < simd::VLEN; i++)
      ASSERT_EQ(res[i], rhs[i] + lhs[i]);
  }
  TEST_F(SimdFloat, SUB)
  {
    vres = vrhs - vlhs;
    for (int i = 0; i < simd::VLEN; i++)
      ASSERT_EQ(res[i], rhs[i] - lhs[i]);
  }
  TEST_F(SimdFloat, MUL)
  {
    vres = vrhs * vlhs;
    for (int i = 0; i < simd::VLEN; i++)
      ASSERT_EQ(res[i], rhs[i] * lhs[i]);
  }
  TEST_F(SimdFloat, DIV)
  {
    vres = vrhs / vlhs;
    for (int i = 0; i < simd::VLEN; i++)
      ASSERT_EQ(res[i], rhs[i] / lhs[i]);
  }
}

namespace 
{
  using SimdInt= SimdTest<int>;
  TEST_F(SimdInt, ADD)
  {
    vres = vrhs + vlhs;
    for (int i = 0; i < simd::VLEN; i++)
      ASSERT_EQ(res[i], rhs[i] + lhs[i]);
  }
  TEST_F(SimdInt, SUB)
  {
    vres = vrhs - vlhs;
    for (int i = 0; i < simd::VLEN; i++)
      EXPECT_EQ(res[i], rhs[i] - lhs[i]);
  }
  TEST_F(SimdInt, MUL)
  {
    vres = vrhs * vlhs;
    for (int i = 0; i < simd::VLEN; i++)
      ASSERT_EQ(res[i], rhs[i] * lhs[i]);
  }
  TEST_F(SimdInt, DIV)
  {
    vres = vrhs / vlhs;
    for (int i = 0; i < simd::VLEN; i++)
      ASSERT_EQ(res[i], rhs[i] / lhs[i]);
  }
}




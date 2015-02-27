#include <gtest/gtest.h>

#include <simd.h>
#include <numeric>
#include <cmath>
#include <array>
#include <algorithm>
#include <gtest/gtest.h>
#include <random>

template<typename T> 
class SimdMemTest : public testing::Test
{
  protected:
    using stype = T;
    using vtype = Simd<T>;
    using vref  = SimdRefT<T>;
    using iref  = SimdIRefT<T>;

    enum {NVEL=32, NSEL = NVEL*vtype::VLEN};
    std::array<stype, NSEL> data;

    SimdMemTest() {}
    virtual ~SimdMemTest() {}

    virtual void SetUp()
    {
      std::random_device rd;
      std::mt19937 gen(rd());
      std::uniform_int_distribution<> dis(1, 1<<16);
      for (int i = 0; i < NSEL; i++)
        data[i] = static_cast<stype>(dis(gen)); 
    }
    
    virtual void TearDown()
    {
    }

    void test_vref()
    {
      // test vload from simd reference to vtype
      vtype vLoad[NVEL];
      for (int i = 0; i < NVEL; i++)
        vLoad[i] = vref(data[i*vtype::VLEN]);

      auto sLoad = (stype*)vLoad;
      for (int i = 0; i < NSEL; i++)
        ASSERT_EQ(sLoad[i], data[i]);

      // test vstore from vtype to simd ref
      std::array<stype,NSEL> store;
      for (int i = 0; i < NVEL; i++)
      {
        vref vStore(store[i*vtype::VLEN]);
        vStore = vLoad[i];
      }

      for (int i = 0; i < NSEL; i++)
      {
        ASSERT_EQ(sLoad[i], store[i]);
        store[i] = 0;
      }

      // test vstore from simd ref to simd ref
      for (int i = 0; i < NVEL; i++)
      {
        vref vStore(store[i*vtype::VLEN]);
        vref vLoad (data [i*vtype::VLEN]);
        vStore = vLoad;
      }

      for (int i = 0; i < NSEL; i++)
        ASSERT_EQ(sLoad[i], store[i]);
    }

    void test_iref()
    {
      std::array<int,NSEL> src, dst;

      std::iota(src.begin(), src.end(), 0);
      std::iota(dst.begin(), dst.end(), 0);
      std::random_shuffle(src.begin(), src.end());
      std::random_shuffle(dst.begin(), dst.end());

      std::array<stype,NSEL> goldGatherStream, goldStreamScatter, goldGatherScatter;
      for (int i = 0; i < NSEL; i++)
      {
        goldGatherStream     [i]  = data[src[i]];
        goldStreamScatter[dst[i]] = data    [i] ;
        goldGatherScatter[dst[i]] = data[src[i]];
      }

      std::array<stype,NSEL> sGatherStream, sStreamScatter, sGatherScatter;
      for (int i = 0; i < NVEL; i++)
      {
        const auto vsrc = SimdRefT<int>(src[i*vtype::VLEN]);
        const auto vdst = SimdRefT<int>(dst[i*vtype::VLEN]);
        const auto gather = iref(data[0],vsrc);
        auto res1 = vref(sGatherStream[i*vtype::VLEN]);
        auto res2 = iref(sStreamScatter[0],vdst);
        auto res3 = iref(sGatherScatter[0],vdst);

        res1 = gather;
        res2 = vref(data[i*vtype::VLEN]);
        res3 = gather;
      }

      for (int i = 0; i < NSEL; i++)
      {
        ASSERT_EQ(sGatherStream [i], goldGatherStream [i]);
        ASSERT_EQ(sStreamScatter[i], goldStreamScatter[i]);
        ASSERT_EQ(sGatherScatter[i], goldGatherScatter[i]);
      }

    }
};

namespace
{
  using SimdMemDouble = SimdMemTest<double>;
  TEST_F(SimdMemDouble,VREF)
  {
    test_vref();
  }

  TEST_F(SimdMemDouble,IREF)
  {
    test_iref();
  }
}

namespace
{
  using SimdMemFloat = SimdMemTest<float>;
  TEST_F(SimdMemFloat,VREF)
  {
    test_vref();
  }

  TEST_F(SimdMemFloat,IREF)
  {
    test_iref();
  }
}

namespace
{
  using SimdMemInt = SimdMemTest<int>;
  TEST_F(SimdMemInt,VREF)
  {
    test_vref();
  }

  TEST_F(SimdMemInt,IREF)
  {
    test_iref();
  }
}


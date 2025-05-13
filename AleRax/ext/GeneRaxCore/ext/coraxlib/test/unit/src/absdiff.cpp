#include <gtest/gtest.h>
#include <corax/util/absdiff.h>

TEST(Util, absdiff)
{
  EXPECT_EQ((absdiff(0, 10)), 10);
  EXPECT_EQ((absdiff(10, 0)), 10);
  EXPECT_EQ((absdiff(std::numeric_limits<unsigned int>::max(), 1)),
            std::numeric_limits<unsigned int>::max() - 1);
  EXPECT_EQ((absdiff(1, std::numeric_limits<unsigned int>::max())),
            std::numeric_limits<unsigned int>::max() - 1);
  EXPECT_EQ((absdiff(42, 10)), 32);
  EXPECT_EQ((absdiff(10, 42)), 32);
  EXPECT_EQ((absdiff(1337, 1337)), 0);

  for (unsigned int low = 1834; low < 4312; low++)
  {
    for (unsigned int high = low; high < 4312; high++)
    {
      assert(high >= low);
      auto diff = high - low;
      EXPECT_EQ(absdiff(low, high), diff);
      EXPECT_EQ(absdiff(high, low), diff);
    }
  }
}

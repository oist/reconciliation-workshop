#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "corax/core/corax_assert.hpp"

using namespace ::testing;

#define EXPECT_CORAX_ASSERT_FAILS(CODE, FAILURE_MESSAGE)                       \
  EXPECT_EXIT({ CODE; }, testing::KilledBySignal(SIGABRT), FAILURE_MESSAGE);

#define ASSERT_CORAX_ASSERT_FAILS(CODE, FAILURE_MESSAGE)                       \
  ASSERT_EXIT({ CODE; }, testing::KilledBySignal(SIGABRT), FAILURE_MESSAGE);

#define EXPECT_CORAX_ASSERT_THROWS(CODE)                                       \
  EXPECT_THROW({ CODE; }, corax_assert_failed_exception);

#define ASSERT_CORAX_ASSERT_THROWS(CODE)                                       \
  EXPECT_THROW({ CODE; }, corax_assert_failed_exception);

// The tests in this file should be run once in Debug and once on Release mode.
// They should also be run once with CORAX_ASSERT_USE_EXCEPTIONS defined and
// once without.

TEST(CoraxAssertTest, corax_assert)
{
  corax_assert(true, "This should never fail.");
  corax_assert(1 < 2, "This should never fail.");
  corax_assert(true != false, "This should never fail.");
  corax_assert(true || (false && true), "This should never fail.");
  corax_assert(false || (true && true), "This should never fail.");

#ifndef NDEBUG
  EXPECT_CORAX_ASSERT_FAILS(
      { corax_assert(false, "This should always fail."); },
      "This should always fail.");

  EXPECT_CORAX_ASSERT_FAILS(
      { corax_assert(3 < 2, "This should always fail."); },
      "This should always fail.");

  EXPECT_CORAX_ASSERT_FAILS(
      { corax_assert(false != false, "This should always fail."); },
      "This should always fail.");

  EXPECT_CORAX_ASSERT_FAILS(
      { corax_assert(false || (false && false), "This should always fail."); },
      "This should always fail.");

  EXPECT_CORAX_ASSERT_FAILS(
      { corax_assert(true && (false || false), "This should always fail."); },
      "This should always fail.");

  EXPECT_CORAX_ASSERT_FAILS({ corax_assert(false); },
                            "Assertion 'false' failed");
#endif
}

TEST(CoraxAssertTest, corax_always_assert)
{
  corax_always_assert(true, "This should never fail.");
  corax_always_assert(1 < 2, "This should never fail.");
  corax_always_assert(true != false, "This should never fail.");
  corax_always_assert(true || (false && true), "This should never fail.");
  corax_always_assert(false || (true && true), "This should never fail.");

#ifndef CORAX_ASSERT_USE_EXCEPTIONS
  EXPECT_CORAX_ASSERT_FAILS(
      { corax_always_assert(false, "This should always fail."); },
      "This should always fail.");

  EXPECT_CORAX_ASSERT_FAILS(
      { corax_always_assert(3 < 2, "This should always fail."); },
      "This should always fail.");

  EXPECT_CORAX_ASSERT_FAILS(
      { corax_always_assert(false != false, "This should always fail."); },
      "This should always fail.");

  EXPECT_CORAX_ASSERT_FAILS(
      {
        corax_always_assert(false || (false && false),
                            "This should always fail.");
      },
      "This should always fail.");

  EXPECT_CORAX_ASSERT_FAILS(
      {
        corax_always_assert(true && (false || false),
                            "This should always fail.");
      },
      "This should always fail.");
#else
  EXPECT_CORAX_ASSERT_THROWS(
      { corax_always_assert(false, "This should always fail."); });

  EXPECT_CORAX_ASSERT_THROWS(
      { corax_always_assert(3 < 2, "This should always fail."); });

  EXPECT_CORAX_ASSERT_THROWS(
      { corax_always_assert(false != false, "This should always fail."); });

  EXPECT_CORAX_ASSERT_THROWS({
    corax_always_assert(false || (false && false), "This should always fail.");
  });

  EXPECT_CORAX_ASSERT_THROWS({
    corax_always_assert(true && (false || false), "This should always fail.");
  });
#endif
}

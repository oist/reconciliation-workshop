#include "environment.hpp"
#include <gtest/gtest.h>

CoraxlibEnvironment* env = nullptr;

int main(int argc, char **argv) {
  env = new CoraxlibEnvironment();
  ::testing::InitGoogleTest(&argc, argv);
  ::testing::AddGlobalTestEnvironment(env);
  return RUN_ALL_TESTS();
}

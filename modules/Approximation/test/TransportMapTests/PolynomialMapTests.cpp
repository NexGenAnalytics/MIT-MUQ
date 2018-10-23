#include <gtest/gtest.h>

#include "MUQ/Approximation/TransportMaps/PolynomialMap.h"

using namespace muq::Approximation;

TEST(PolynomialMap, BasicTest) {
  auto map = std::make_shared<PolynomialMap>();
}

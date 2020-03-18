#include "TestingUtility.h"

#include "Constants.h"

using namespace pfc;


template <typename Real>
class ConstantsTest : public ::testing::Test {
};

typedef ::testing::Types<float, double> types;
TYPED_TEST_CASE(ConstantsTest, types);

TYPED_TEST(ConstantsTest, Values)
{
    // Check that all members of Constants<T> can be instantiated for all T in types
    // and are the same as casted corresponding values of constants namespace
    ASSERT_EQ(static_cast<TypeParam>(constants::pi), Constants<TypeParam>::pi());
    ASSERT_EQ(static_cast<TypeParam>(constants::c), Constants<TypeParam>::c());
    ASSERT_EQ(static_cast<TypeParam>(constants::lightVelocity), Constants<TypeParam>::lightVelocity());
    ASSERT_EQ(static_cast<TypeParam>(constants::electronCharge), Constants<TypeParam>::electronCharge());
    ASSERT_EQ(static_cast<TypeParam>(constants::electronMass), Constants<TypeParam>::electronMass());
    ASSERT_EQ(static_cast<TypeParam>(constants::protonMass), Constants<TypeParam>::protonMass());
    ASSERT_EQ(static_cast<TypeParam>(constants::planck), Constants<TypeParam>::planck());
    ASSERT_EQ(static_cast<TypeParam>(constants::eV), Constants<TypeParam>::eV());
    ASSERT_EQ(static_cast<TypeParam>(constants::meV), Constants<TypeParam>::meV());
}

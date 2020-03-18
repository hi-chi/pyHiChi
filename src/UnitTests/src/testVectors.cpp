#include "TestingUtility.h"

#include "Vectors.h"

#include <iostream>
#include <sstream>

using namespace pfc;


TEST(VectorsTest, Int3DefaultConstructor)
{
    Int3 zero;
    ASSERT_EQ(zero.x, 0);
    ASSERT_EQ(zero.y, 0);
    ASSERT_EQ(zero.z, 0);
}

TEST(VectorsTest, Int3Constructor)
{
    Int3 v(3, 5, -2);
    ASSERT_EQ(v.x, 3);
    ASSERT_EQ(v.y, 5);
    ASSERT_EQ(v.z, -2);
}

TEST(VectorsTest, Int3IndexAccess)
{
    Int3 v(-6, 1, 5);
    ASSERT_EQ(v.x, v[0]);
    ASSERT_EQ(v.y, v[1]);
    ASSERT_EQ(v.z, v[2]);
    const Int3 u(7, -13, -4);
    ASSERT_EQ(u.x, u[0]);
    ASSERT_EQ(u.y, u[1]);
    ASSERT_EQ(u.z, u[2]);
}

TEST(VectorsTest, Int3Volume)
{
    Int3 v(11, 15, 2);
    ASSERT_EQ(v.volume(), 11 * 15 * 2);
}

TEST(VectorsTest, Int3Sum)
{
    Int3 v1(3, 5, -2), v2(6, -8, 2);
    Int3 res = v1 + v2;
    ASSERT_EQ_INT3(res, Int3(9, -3, 0));
}

TEST(VectorsTest, Int3SumAssignment)
{
    Int3 v1(3, 5, -2), v2(6, -8, 2);
    v1 += v2;
    ASSERT_EQ_INT3(v1, Int3(9, -3, 0));
}

TEST(VectorsTest, Int3Difference)
{
    Int3 v1(3, 5, -2), v2(6, -8, 2);
    Int3 res = v1 - v2;
    ASSERT_EQ_INT3(res, Int3(-3, 13, -4));
}

TEST(VectorsTest, Int3DifferenceAssignment)
{
    Int3 v1(3, 5, -2), v2(6, -8, 2);
    v1 -= v2;
    ASSERT_EQ_INT3(v1, Int3(-3, 13, -4));
}

TEST(VectorsTest, Int3Product)
{
    Int3 v1(3, 5, -2), v2(6, -8, 2);
    Int3 res = v1 * v2;
    ASSERT_EQ_INT3(res, Int3(18, -40, -4));
}

TEST(VectorsTest, Int3ProductAssignment)
{
    Int3 v1(3, 5, -2), v2(6, -8, 2);
    v1 *= v2;
    ASSERT_EQ_INT3(v1, Int3(18, -40, -4));
}

TEST(VectorsTest, Int3ProductScalar)
{
    Int3 v1(3, 5, -2);
    Int3 res = v1 * 5;
    ASSERT_EQ_INT3(res, Int3(15, 25, -10));
}

TEST(VectorsTest, Int3ProductScalarAssignment)
{
    Int3 v1(3, 5, -2);
    v1 *= 5;
    ASSERT_EQ_INT3(v1, Int3(15, 25, -10));
}

TEST(VectorsTest, Int3Quotient)
{
    Int3 v1(13, 25, 21), v2(6, 8, 2);
    Int3 res = v1 / v2;
    ASSERT_EQ_INT3(res, Int3(2, 3, 10));
}

TEST(VectorsTest, Int3QuotientAssignment)
{
    Int3 v1(13, 25, 21), v2(6, 8, 2);
    v1 /= v2;
    ASSERT_EQ_INT3(v1, Int3(2, 3, 10));
}

TEST(VectorsTest, Int3QuotientScalar)
{
    Int3 v1(13, 25, 21);
    Int3 res = v1 / 4;
    ASSERT_EQ_INT3(res, Int3(3, 6, 5));
}

TEST(VectorsTest, Int3QuotientScalarAssignment)
{
    Int3 v1(13, 25, 21);
    v1 /= 4;
    ASSERT_EQ_INT3(v1, Int3(3, 6, 5));
}

TEST(VectorsTest, Int3Equals)
{
    Int3 v1(13, 25, 21), v2(6, 8, 2);
    ASSERT_FALSE(v1 == v2);
    ASSERT_TRUE((v1 == v1) && (v2 == v2));
}

TEST(VectorsTest, Int3NotEquals)
{
    Int3 v1(13, 25, 21), v2(6, 8, 2);
    ASSERT_TRUE(v1 != v2);
    ASSERT_FALSE((v1 != v1) || (v2 != v2));
}

TEST(VectorsTest, Int3Lesser)
{
    Int3 v1(13, 25, 21), v2(6, 8, 2), v3(14, 6, 1);
    ASSERT_FALSE(v1 < v2);
    ASSERT_TRUE(v2 < v1);
    ASSERT_FALSE(v1 < v3);
    ASSERT_FALSE(v3 < v1);
    ASSERT_FALSE(v2 < v3);
    ASSERT_FALSE(v3 < v2);
}

TEST(VectorsTest, Int3LesserOrEqual)
{
    Int3 v1(13, 25, 21), v2(6, 8, 2), v3(14, 25, 1);
    ASSERT_FALSE(v1 <= v2);
    ASSERT_TRUE(v2 <= v1);
    ASSERT_FALSE(v1 <= v3);
    ASSERT_FALSE(v3 <= v1);
    ASSERT_FALSE(v2 <= v3);
    ASSERT_FALSE(v3 <= v2);
}

TEST(VectorsTest, Int3Greater)
{
    Int3 v1(13, 25, 21), v2(6, 8, 2), v3(14, 6, 1);
    ASSERT_TRUE(v1 > v2);
    ASSERT_FALSE(v2 > v1);
    ASSERT_FALSE(v1 > v3);
    ASSERT_FALSE(v3 > v1);
    ASSERT_FALSE(v2 > v3);
    ASSERT_FALSE(v3 > v2);
}

TEST(VectorsTest, Int3GreaterOrEqual)
{
    Int3 v1(13, 25, 21), v2(6, 8, 2), v3(14, 8, 1);
    ASSERT_TRUE(v1 >= v2);
    ASSERT_FALSE(v2 >= v1);
    ASSERT_FALSE(v1 >= v3);
    ASSERT_FALSE(v3 >= v1);
    ASSERT_FALSE(v2 >= v3);
    ASSERT_FALSE(v3 >= v2);
}

TEST(VectorsTest, Int3OperatorOutput)
{
    std::ostringstream os;
    os << Int3(2, -1, 0);
    ASSERT_EQ("(2, -1, 0)", os.str());
}


TEST(VectorsTest, FP3DefaultConstructor)
{
    FP3 zero;
    ASSERT_EQ(zero.x, (FP)0);
    ASSERT_EQ(zero.y, (FP)0);
    ASSERT_EQ(zero.z, (FP)0);
}

TEST(VectorsTest, FP3Constructor)
{
    FP3 v(3.1, 5.5, -2.6);
    ASSERT_EQ(v.x, (FP)3.1);
    ASSERT_EQ(v.y, (FP)5.5);
    ASSERT_EQ(v.z, (FP)-2.6);
}

TEST(VectorsTest, FP3ConstructorInt3)
{
    Int3 v1(3, 5, -1);
    FP3 v2(v1);
    ASSERT_EQ(v1.x, v2.x);
    ASSERT_EQ(v1.y, v2.y);
    ASSERT_EQ(v1.z, v2.z);
}

TEST(VectorsTest, FP3IndexAccess)
{
    FP3 v(3.1, 5.5, -2.6);
    ASSERT_EQ(v.x, v[0]);
    ASSERT_EQ(v.y, v[1]);
    ASSERT_EQ(v.z, v[2]);
    const FP3 u(3.1, 5.5, -2.6);
    ASSERT_EQ(u.x, u[0]);
    ASSERT_EQ(u.y, u[1]);
    ASSERT_EQ(u.z, u[2]);
}

TEST(VectorsTest, FP3Volume)
{
    FP3 v(11, 15, 2);
    ASSERT_EQ(v.volume(), 11 * 15 * 2);
}

TEST(VectorsTest, FP3Norm)
{
    FP3 v(2, 3, 6);
    ASSERT_EQ(v.norm(), (FP)7);
}

TEST(VectorsTest, FP3Norm2)
{
    FP3 v(2, 3, 6);
    ASSERT_EQ(v.norm2(), (FP)49);
}

TEST(VectorsTest, FP3Sum)
{
    FP3 v1(3, 5, -2), v2(6, -8, 2);
    FP3 res = v1 + v2;
    ASSERT_EQ_FP3(res, FP3(9, -3, 0));
}

TEST(VectorsTest, FP3SumAssignment)
{
    FP3 v1(3, 5, -2), v2(6, -8, 2);
    v1 += v2;
    ASSERT_EQ_FP3(v1, FP3(9, -3, 0));
}

TEST(VectorsTest, FP3Difference)
{
    FP3 v1(3, 5, -2), v2(6, -8, 2);
    FP3 res = v1 - v2;
    ASSERT_EQ_FP3(res, FP3(-3, 13, -4));
}

TEST(VectorsTest, FP3DifferenceAssignment)
{
    FP3 v1(3, 5, -2), v2(6, -8, 2);
    v1 -= v2;
    ASSERT_EQ_FP3(v1, FP3(-3, 13, -4));
}

TEST(VectorsTest, FP3Product)
{
    FP3 v1(3, 5, -2), v2(6, -8, 2);
    FP3 res = v1 * v2;
    ASSERT_EQ_FP3(res, FP3(18, -40, -4));
}

TEST(VectorsTest, FP3ProductAssignment)
{
    FP3 v1(3, 5, -2), v2(6, -8, 2);
    v1 *= v2;
    ASSERT_EQ_FP3(v1, FP3(18, -40, -4));
}

TEST(VectorsTest, FP3ProductInt3)
{
    FP3 v1(3, 5, -2);
    Int3 v2(6, -8, 2);
    FP3 res = v1 * FP3(v2);
    ASSERT_EQ_FP3(res, FP3(18, -40, -4));
}

TEST(VectorsTest, FP3ProductInt3Assignment)
{
    FP3 v1(3, 5, -2);
    Int3 v2(6, -8, 2);
    v1 *= FP3(v2);
    ASSERT_EQ_FP3(v1, FP3(18, -40, -4));
}

TEST(VectorsTest, FP3ProductScalar)
{
    FP3 v1(3, 5, -2);
    FP3 res1 = v1 * 5.0;
    FP3 res2 = 5.0 * v1;
    ASSERT_EQ_FP3(res1, FP3(15, 25, -10));
    ASSERT_EQ_FP3(res2, FP3(15, 25, -10));
}

TEST(VectorsTest, FP3ProductScalarAssignment)
{
    FP3 v1(3, 5, -2);
    v1 *= 5.0;
    ASSERT_EQ_FP3(v1, FP3(15, 25, -10));
}

TEST(VectorsTest, FP3Quotient)
{
    FP3 v1(12, 24, 21), v2(6, 8, 3);
    FP3 res = v1 / v2;
    ASSERT_EQ_FP3(res, FP3(2, 3, 7));
}

TEST(VectorsTest, FP3QuotientAssignment)
{
    FP3 v1(12, 24, 21), v2(6, 8, 3);
    v1 /= v2;
    ASSERT_EQ_FP3(v1, FP3(2, 3, 7));
}

TEST(VectorsTest, FP3QuotientInt3)
{
    FP3 v1(12, 24, 21);
    Int3 v2(6, 8, 3);
    FP3 res = v1 / FP3(v2);
    ASSERT_EQ_FP3(res, FP3(2, 3, 7));
}

TEST(VectorsTest, FP3QuotientInt3Assignment)
{
    FP3 v1(12, 24, 21);
    Int3 v2(6, 8, 3);
    v1 /= FP3(v2);
    ASSERT_EQ_FP3(v1, FP3(2, 3, 7));
}

TEST(VectorsTest, FP3QuotientScalar)
{
    FP3 v1(12, -24, 22);
    FP3 res = v1 / 2.0;
    ASSERT_EQ_FP3(res, FP3(6, -12, 11));
}

TEST(VectorsTest, FP3QuotientScalarAssignment)
{
    FP3 v1(12, -24, 22);
    v1 /= 2.0;
    ASSERT_EQ_FP3(v1, FP3(6, -12, 11));
}

TEST(VectorsTest, FP3Equals)
{
    FP3 v1(13, 25, 21), v2(6, 8, 2);
    ASSERT_FALSE(v1 == v2);
    ASSERT_TRUE((v1 == v1) && (v2 == v2));
}

TEST(VectorsTest, FP3NotEquals)
{
    FP3 v1(13, 25, 21), v2(6, 8, 2);
    ASSERT_TRUE(v1 != v2);
    ASSERT_FALSE((v1 != v1) || (v2 != v2));
}

TEST(VectorsTest, FP3Lesser)
{
    FP3 v1(13, 25, 21), v2(6, 8, 2), v3(14, 6, 1);
    ASSERT_FALSE(v1 < v2);
    ASSERT_TRUE(v2 < v1);
    ASSERT_FALSE(v1 < v3);
    ASSERT_FALSE(v3 < v1);
    ASSERT_FALSE(v2 < v3);
    ASSERT_FALSE(v3 < v2);
}

TEST(VectorsTest, FP3LesserOrEqual)
{
    FP3 v1(13, 25, 21), v2(6, 8, 2), v3(14, 25, 1);
    ASSERT_FALSE(v1 <= v2);
    ASSERT_TRUE(v2 <= v1);
    ASSERT_FALSE(v1 <= v3);
    ASSERT_FALSE(v3 <= v1);
    ASSERT_FALSE(v2 <= v3);
    ASSERT_FALSE(v3 <= v2);
}

TEST(VectorsTest, FP3Greater)
{
    FP3 v1(13, 25, 21), v2(6, 8, 2), v3(14, 6, 1);
    ASSERT_TRUE(v1 > v2);
    ASSERT_FALSE(v2 > v1);
    ASSERT_FALSE(v1 > v3);
    ASSERT_FALSE(v3 > v1);
    ASSERT_FALSE(v2 > v3);
    ASSERT_FALSE(v3 > v2);
}

TEST(VectorsTest, FP3GreaterOrEqual)
{
    FP3 v1(13, 25, 21), v2(6, 8, 2), v3(14, 8, 1);
    ASSERT_TRUE(v1 >= v2);
    ASSERT_FALSE(v2 >= v1);
    ASSERT_FALSE(v1 >= v3);
    ASSERT_FALSE(v3 >= v1);
    ASSERT_FALSE(v2 >= v3);
    ASSERT_FALSE(v3 >= v2);
}

TEST(VectorsTest, FP3Floor)
{
    FP3 v(5.2, 7.5, -1.4);
    ASSERT_EQ_INT3(floor(v), Int3(5, 7, -2));
}

TEST(VectorsTest, FP3Truncate)
{
    FP3 v(5.2, 7.5, -1.4);
    ASSERT_EQ_INT3(truncate(v), Int3(5, 7, -1));
}

TEST(VectorsTest, FP3Dot)
{
    FP3 v1(5, 7, -1), v2(1, -2, 3);
    ASSERT_EQ(dot(v1, v2), 5 * 1 + 7 * (-2) + (-1) * 3);
}

TEST(VectorsTest, FP3Cross)
{
    FP3 v1(5, 7, -1), v2(1, -2, 3);
    ASSERT_EQ_FP3(cross(v1, v2), FP3(19, -16, -17));
}

TEST(VectorsTest, FP3Dist)
{
    FP3 v1(5, 7, -1), v2(3, 4, -7);
    ASSERT_EQ(dist(v1, v2), 7);
}

TEST(VectorsTest, FP3OperatorOutput)
{
    std::ostringstream os;
    os << FP3(2, -0.25, 1.5);
    ASSERT_EQ("(2, -0.25, 1.5)", os.str());
}

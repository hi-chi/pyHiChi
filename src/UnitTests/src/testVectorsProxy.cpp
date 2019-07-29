#include "TestingUtility.h"

#include "Vectors.h"
#include "VectorsProxy.h"

#include <iostream>
#include <sstream>

using namespace pfc;


TEST(VectorsProxyTest, Int3FullConstructor)
{
	int x = 3, y = 5, z = -2;
	Int3Proxy vP(x, y, z);
	ASSERT_EQ(vP.x, 3);
	ASSERT_EQ(vP.y, 5);
	ASSERT_EQ(vP.z, -2);
}

TEST(VectorsProxyTest, Int3Constructor)
{
	Int3 v(3, 5, -2);
	Int3 vP(v);
	ASSERT_EQ(v.x, 3);
	ASSERT_EQ(v.y, 5);
	ASSERT_EQ(v.z, -2);
	ASSERT_EQ(vP.x, 3);
	ASSERT_EQ(vP.y, 5);
	ASSERT_EQ(vP.z, -2);
}

TEST(VectorsProxyTest, Int3IndexAccess)
{
	Int3 v(-6, 1, 5);
	Int3 vP(v);
	ASSERT_EQ(vP.x, vP[0]);
	ASSERT_EQ(vP.y, vP[1]);
	ASSERT_EQ(vP.z, vP[2]);
	ASSERT_EQ(v.x, vP[0]);
	ASSERT_EQ(v.y, vP[1]);
	ASSERT_EQ(v.z, vP[2]);
	ASSERT_EQ(vP.x, v[0]);
	ASSERT_EQ(vP.y, v[1]);
	ASSERT_EQ(vP.z, v[2]);
	Int3 u(7, -13, -4);
	const Int3Proxy uP(u);
	ASSERT_EQ(u.x, uP[0]);
	ASSERT_EQ(u.y, uP[1]);
	ASSERT_EQ(u.z, uP[2]);
	ASSERT_EQ(uP.x, u[0]);
	ASSERT_EQ(uP.y, u[1]);
	ASSERT_EQ(uP.z, u[2]);
	ASSERT_EQ(uP.x, uP[0]);
	ASSERT_EQ(uP.y, uP[1]);
	ASSERT_EQ(uP.z, uP[2]);
}

TEST(VectorsProxyTest, Int3Volume)
{
	Int3 v(11, 15, 2);
	Int3Proxy vP(v);
	ASSERT_EQ(vP.volume(), 11 * 15 * 2);
}

TEST(VectorsProxyTest, Int3Sum)
{
	Int3 v1(3, 5, -2), v2(6, -8, 2);
	Int3Proxy v1P(v1), v2P(v2);

	Int3 res1 = v1 + v2;
	Int3 res2 = Int3Proxy(v1) + v2P;
	Int3 res3 = v1P + Int3Proxy(v2);
	Int3 res4 = v1P + v2P;
	ASSERT_EQ_INT3(res1, Int3(9, -3, 0));
	ASSERT_EQ_INT3(res2, res1);
	ASSERT_EQ_INT3(res3, res1);
	ASSERT_EQ_INT3(res4, res1);
}

TEST(VectorsProxyTest, Int3SumAssignment)
{
	Int3 v11(3, 5, -2), v2(6, -8, 2);
	Int3 v12(v11), v13(v11), v14(v11);
	Int3Proxy v2P(v2), v12P(v12), v13P(v13), v14P(v14);

	v11 += v2;
	v12P += v2P;
	v13P += Int3Proxy(v2);

	ASSERT_EQ_INT3(v11, Int3(9, -3, 0));
	ASSERT_EQ_INT3(v12P, v11);
	ASSERT_EQ_INT3(v12P, v12);
	ASSERT_EQ_INT3(v13P, v11);
}

TEST(VectorsProxyTest, Int3Difference)
{
	Int3 v1(3, 5, -2), v2(6, -8, 2);
	Int3Proxy v1P(v1), v2P(v2);
	Int3 res1 = v1 - v2;
	Int3 res2 = Int3Proxy(v1) - v2P;
	Int3 res3 = v1P - Int3Proxy(v2);
	Int3 res4 = v1P - v2P;
	ASSERT_EQ_INT3(res1, Int3(-3, 13, -4));
	ASSERT_EQ_INT3(res2, res1);
	ASSERT_EQ_INT3(res3, res1);
	ASSERT_EQ_INT3(res4, res1);
}

TEST(VectorsProxyTest, Int3DifferenceAssignment)
{
	Int3 v11(3, 5, -2), v2(6, -8, 2);
	Int3 v12(v11), v13(v11), v14(v11);
	Int3Proxy v2P(v2), v12P(v12), v13P(v13), v14P(v14);

	v11 -= v2;
	v12P -= v2P;
	v13P -= Int3Proxy(v2);

	ASSERT_EQ_INT3(v11, Int3(-3, 13, -4));
	ASSERT_EQ_INT3(v12P, v11);
	ASSERT_EQ_INT3(v12P, v12);
	ASSERT_EQ_INT3(v13P, v11);
}

TEST(VectorsProxyTest, Int3Product)
{
	Int3 v1(3, 5, -2), v2(6, -8, 2);
	Int3Proxy v1P(v1), v2P(v2);

	Int3 res1 = v1 * v2;
	Int3 res2 = Int3Proxy(v1) * v2P;
	Int3 res3 = v1P * Int3Proxy(v2);
	Int3 res4 = v1P * v2P;

	ASSERT_EQ_INT3(res1, Int3(18, -40, -4));
	ASSERT_EQ_INT3(res2, res1);
	ASSERT_EQ_INT3(res3, res1);
	ASSERT_EQ_INT3(res4, res1);
}

TEST(VectorsProxyTest, Int3ProductAssignment)
{
	Int3 v11(3, 5, -2), v2(6, -8, 2);
	Int3 v12(v11), v13(v11), v14(v11);
	Int3Proxy v2P(v2), v12P(v12), v13P(v13), v14P(v14);

	v11 *= v2;
	v12P *= v2P;
	v13P *= Int3Proxy(v2);

	ASSERT_EQ_INT3(v11, Int3(18, -40, -4));
	ASSERT_EQ_INT3(v12P, v11);
	ASSERT_EQ_INT3(v12P, v12);
	ASSERT_EQ_INT3(v13P, v11);
}

TEST(VectorsProxyTest, Int3ProductScalar)
{
	Int3 v1(3, 5, -2);
	Int3Proxy v1P(v1);

	Int3 res1 = v1 * 5;
	Int3 res2 = v1P * 5;

	ASSERT_EQ_INT3(res1, Int3(15, 25, -10));
	ASSERT_EQ_INT3(res2, res1);
}

TEST(VectorsProxyTest, Int3ProductScalarAssignment)
{
	Int3 v1(3, 5, -2);
	Int3 v2(v1);
	Int3Proxy v2P(v2);

	v1 *= 5;
	v2P *= 5;

	ASSERT_EQ_INT3(v1, Int3(15, 25, -10));
	ASSERT_EQ_INT3(v2, v1);
	ASSERT_EQ_INT3(v2P, v1);
}

TEST(VectorsProxyTest, Int3Quotient)
{
	Int3 v1(13, 25, 21), v2(6, 8, 2);
	Int3Proxy v1P(v1), v2P(v2);

	Int3 res1 = v1 / v2;
	Int3 res2 = Int3Proxy(v1) / v2P;
	Int3 res3 = v1P / Int3Proxy(v2);
	Int3 res4 = v1P / v2P;

	ASSERT_EQ_INT3(res1, Int3(2, 3, 10));
	ASSERT_EQ_INT3(res2, res1);
	ASSERT_EQ_INT3(res3, res1);
	ASSERT_EQ_INT3(res4, res1);
}

TEST(VectorsProxyTest, Int3QuotientAssignment)
{
	Int3 v11(13, 25, 21), v2(6, 8, 2);
	Int3 v12(v11), v13(v11), v14(v11);
	Int3Proxy v2P(v2), v12P(v12), v13P(v13), v14P(v14);

	v11 /= v2;
	v12P /= v2P;
	v13P /= Int3Proxy(v2);

	ASSERT_EQ_INT3(v11, Int3(2, 3, 10));
	ASSERT_EQ_INT3(v12P, v11);
	ASSERT_EQ_INT3(v12P, v12);
	ASSERT_EQ_INT3(v13P, v11);
}

TEST(VectorsProxyTest, Int3QuotientScalar)
{
	Int3 v1(13, 25, 21);
	Int3Proxy v1P(v1);

	Int3 res1 = v1 / 4;
	Int3 res2 = v1P / 4;
	
	ASSERT_EQ_INT3(res1, Int3(3, 6, 5));
	ASSERT_EQ_INT3(res2, res1);
}

TEST(VectorsProxyTest, Int3QuotientScalarAssignment)
{
	Int3 v1(13, 25, 21);
	Int3 v2(v1);
	Int3Proxy v2P(v2);
	
	v1 /= 4;
	v2P /= 4;
	
	ASSERT_EQ_INT3(v1, Int3(3, 6, 5));
	ASSERT_EQ_INT3(v1, v2);
	ASSERT_EQ_INT3(v1, v2P);
}

TEST(VectorsProxyTest, Int3Equals)
{
	Int3 v1(13, 25, 21), v2(6, 8, 2);
	Int3Proxy v1P(v1), v2P(v2);
	
	ASSERT_FALSE((v1P == v2P));
	ASSERT_TRUE((v1P == v1P) && (v2P == v2P));
}

TEST(VectorsProxyTest, Int3NotEquals)
{
	Int3 v1(13, 25, 21), v2(6, 8, 2);
	Int3Proxy v1P(v1), v2P(v2);
	
	ASSERT_TRUE((v1P != v2P));
	ASSERT_FALSE((v1P != v1P) || (v2P != v2P));
}

TEST(VectorsProxyTest, Int3Lesser)
{
	Int3 v1(13, 25, 21), v2(6, 8, 2), v3(14, 6, 1);
	Int3Proxy v1P(v1), v2P(v2), v3P(v3);
	
	ASSERT_FALSE(v1P < v2P);
	ASSERT_TRUE(v2P < v1P);
	ASSERT_FALSE(v1P < v3P);
	ASSERT_FALSE(v3P < v1P);
	ASSERT_FALSE(v2P < v3P);
	ASSERT_FALSE(v3P < v2P);
}

TEST(VectorsProxyTest, Int3LesserOrEqual)
{
	Int3 v1(13, 25, 21), v2(6, 8, 2), v3(14, 25, 1);
	Int3Proxy v1P(v1), v2P(v2), v3P(v3);
	
	ASSERT_FALSE(v1P <= v2P);
	ASSERT_TRUE(v2P <= v1P);
	ASSERT_FALSE(v1P <= v3P);
	ASSERT_FALSE(v3P <= v1P);
	ASSERT_FALSE(v2P <= v3P);
	ASSERT_FALSE(v3P <= v2P);
}

TEST(VectorsProxyTest, Int3Greater)
{
	Int3 v1(13, 25, 21), v2(6, 8, 2), v3(14, 6, 1);
	Int3Proxy v1P(v1), v2P(v2), v3P(v3);
	
	ASSERT_TRUE(v1P > v2P);
	ASSERT_FALSE(v2P > v1P);
	ASSERT_FALSE(v1P > v3P);
	ASSERT_FALSE(v3P > v1P);
	ASSERT_FALSE(v2P > v3P);
	ASSERT_FALSE(v3P > v2P);
}

TEST(VectorsPoxyTest, Int3GreaterOrEqual)
{
	Int3 v1(13, 25, 21), v2(6, 8, 2), v3(14, 8, 1);
	Int3Proxy v1P(v1), v2P(v2), v3P(v3);
	
	ASSERT_TRUE(v1P >= v2P);
	ASSERT_FALSE(v2P >= v1P);
	ASSERT_FALSE(v1P >= v3P);
	ASSERT_FALSE(v3P >= v1P);
	ASSERT_FALSE(v2P >= v3P);
	ASSERT_FALSE(v3P >= v2P);
}

TEST(VectorsProxyTest, Int3OperatorOutput)
{
	std::ostringstream os;
	Int3 v(2, -1, 0);
	Int3Proxy vP(v);
	
	os << vP;
	
	ASSERT_EQ("(2, -1, 0)", os.str());
}


TEST(VectorsProxyTest, FP3FullConstructor)
{
	FP x = 3.1, y = 5.5, z = -2.6;
	FP3Proxy vP(x, y, z);
	
	ASSERT_EQ(vP.x, 3.1);
	ASSERT_EQ(vP.y, 5.5);
	ASSERT_EQ(vP.z, -2.6);
}

TEST(VectorsProxyTest, FP3Constructor)
{
	FP3 v(3.1, 5.5, -2.6);
	FP3 vP(v);

	ASSERT_EQ(vP.x, (FP)3.1);
	ASSERT_EQ(vP.y, (FP)5.5);
	ASSERT_EQ(vP.z, (FP)-2.6);
}

TEST(VectorsProxyTest, FP3IndexAccess)
{
	FP3 v(3.1, 5.5, -2.6);
	FP3 vP(v);

	ASSERT_EQ(vP.x, vP[0]);
	ASSERT_EQ(vP.y, vP[1]);
	ASSERT_EQ(vP.z, vP[2]);
	ASSERT_EQ(v.x, vP[0]);
	ASSERT_EQ(v.y, vP[1]);
	ASSERT_EQ(v.z, vP[2]);
	ASSERT_EQ(vP.x, v[0]);
	ASSERT_EQ(vP.y, v[1]);
	ASSERT_EQ(vP.z, v[2]);
	
	FP3 u(3.1, 5.5, -2.6);
	const FP3Proxy uP(u);
	
	ASSERT_EQ(u.x, uP[0]);
	ASSERT_EQ(u.y, uP[1]);
	ASSERT_EQ(u.z, uP[2]);
	ASSERT_EQ(uP.x, u[0]);
	ASSERT_EQ(uP.y, u[1]);
	ASSERT_EQ(uP.z, u[2]);
	ASSERT_EQ(uP.x, uP[0]);
	ASSERT_EQ(uP.y, uP[1]);
	ASSERT_EQ(uP.z, uP[2]);
}

TEST(VectorsProxyTest, FP3Volume)
{
	FP3 v(11, 15, 2);
	FP3Proxy vP(v);
	
	ASSERT_EQ(vP.volume(), 11 * 15 * 2);
}

TEST(VectorsProxyTest, FP3Norm)
{
	FP3 v(2, 3, 6);
	FP3Proxy vP(v);
	
	ASSERT_EQ(vP.norm(), (FP)7);
}

TEST(VectorsProxyTest, FP3Norm2)
{
	FP3 v(2, 3, 6);
	FP3Proxy vP(v);
	
	ASSERT_EQ(vP.norm2(), (FP)49);
}

TEST(VectorsProxyTest, FP3Sum)
{
	FP3 v1(3, 5, -2), v2(6, -8, 2);
	FP3Proxy v1P(v1), v2P(v2);

	FP3 res1 = v1 + v2;
	FP3 res2 = FP3Proxy(v1) + v2P;
	FP3 res3 = v1P + FP3Proxy(v2);
	FP3 res4 = v1P + v2P;
	
	ASSERT_EQ_FP3(res1, FP3(9, -3, 0));
	ASSERT_EQ_FP3(res2, res1);
	ASSERT_EQ_FP3(res3, res1);
	ASSERT_EQ_FP3(res4, res1);
}

TEST(VectorsProxyTest, FP3SumAssignment)
{
	Int3 v11(3, 5, -2), v2(6, -8, 2);
	FP3 v12(v11), v13(v11), v2f(v2);
	FP3Proxy v2P(v2f), v12P(v12), v13P(v13);

	v11 += v2;
	v12P += v2P;
	v13P += FP3Proxy(v2f);

	ASSERT_EQ_FP3(v11, FP3(9, -3, 0));
	ASSERT_EQ_FP3(v12P, v11);
	ASSERT_EQ_FP3(v12P, v12);
	ASSERT_EQ_FP3(v13P, v11);
}

TEST(VectorsProxyTest, FP3Difference)
{
	FP3 v1(3, 5, -2), v2(6, -8, 2);
	FP3Proxy v1P(v1), v2P(v2);
	
	FP3 res1 = v1 - v2;
	FP3 res2 = FP3Proxy(v1) - v2P;
	FP3 res3 = v1P - FP3Proxy(v2);
	FP3 res4 = v1P - v2P;
	
	ASSERT_EQ_FP3(res1, FP3(-3, 13, -4));
	ASSERT_EQ_FP3(res2, res1);
	ASSERT_EQ_FP3(res3, res1);
	ASSERT_EQ_FP3(res4, res1);
}

TEST(VectorsProxyTest, FP3DifferenceAssignment)
{
	FP3 v11(3, 5, -2), v2(6, -8, 2);
	FP3 v12(v11), v13(v11);
	FP3Proxy v2P(v2), v12P(v12), v13P(v13);

	v11 -= v2;
	v12P -= v2P;
	v13P -= FP3Proxy(v2);

	ASSERT_EQ_FP3(v11, FP3(-3, 13, -4));
	ASSERT_EQ_FP3(v12P, v11);
	ASSERT_EQ_FP3(v12P, v12);
	ASSERT_EQ_FP3(v13P, v11);
}

TEST(VectorsProxyTest, FP3Product)
{
	FP3 v1(3, 5, -2), v2(6, -8, 2);
	FP3Proxy v1P(v1), v2P(v2);
	
	FP3 res1 = v1 * v2;
	FP3 res2 = FP3Proxy(v1) * v2P;
	FP3 res3 = v1P * FP3Proxy(v2);
	FP3 res4 = v1P * v2P;
	
	ASSERT_EQ_FP3(res1, FP3(18, -40, -4));
	ASSERT_EQ_FP3(res2, res1);
	ASSERT_EQ_FP3(res3, res1);
	ASSERT_EQ_FP3(res4, res1);
}

TEST(VectorsProxyTest, FP3ProductAssignment)
{
	FP3 v11(3, 5, -2), v2(6, -8, 2);
	FP3 v12(v11), v13(v11);
	FP3Proxy v2P(v2), v12P(v12), v13P(v13);

	v11 *= v2;
	v12P *= v2P;
	v13P *= FP3Proxy(v2);

	ASSERT_EQ_FP3(v11, FP3(18, -40, -4));
	ASSERT_EQ_FP3(v12P, v11);
	ASSERT_EQ_FP3(v12P, v12);
	ASSERT_EQ_FP3(v13P, v11);
}

TEST(VectorsProxyTest, FP3ProductScalar)
{
	FP3 v1(3, 5, -2);
	FP3Proxy v1P(v1);
	
	FP3 res1 = v1P * 5.0;
	FP3 res2 = 5.0 * v1P;
	
	ASSERT_EQ_FP3(res1, FP3(15, 25, -10));
	ASSERT_EQ_FP3(res2, FP3(15, 25, -10));
}

TEST(VectorsProxyTest, FP3ProductScalarAssignment)
{
	FP3 v1(3, 5, -2);
	FP3Proxy v1P(v1);
	
	v1P *= 5.0;
	
	ASSERT_EQ_FP3(v1, FP3(15, 25, -10));
	ASSERT_EQ_FP3(v1P, FP3(15, 25, -10));
}

TEST(VectorsProxyTest, FP3Quotient)
{
	FP3 v1(12, 24, 21), v2(6, 8, 3);
	FP3Proxy v1P(v1), v2P(v2);
	
	FP3 res = v1P / v2P;
	
	ASSERT_EQ_FP3(res, FP3(2, 3, 7));
}

TEST(VectorsProxyTest, FP3QuotientAssignment)
{
	FP3 v1(12, 24, 21), v2(6, 8, 3);
	FP3Proxy v1P(v1), v2P(v2);
	
	v1P /= v2P;
	
	ASSERT_EQ_FP3(v1, FP3(2, 3, 7));
	ASSERT_EQ_FP3(v1P, v1);
}


TEST(VectorsProxyTest, FP3QuotientScalar)
{
	FP3 v1(12, -24, 22);
	FP3Proxy v1P(v1);
	
	FP3 res = v1P / 2.0;
	
	ASSERT_EQ_FP3(res, FP3(6, -12, 11));
}

TEST(VectorsProxyTest, FP3QuotientScalarAssignment)
{
	FP3 v1(12, -24, 22);
	FP3Proxy v1P(v1);
	
	v1P /= 2.0;
	
	ASSERT_EQ_FP3(v1, FP3(6, -12, 11));
	ASSERT_EQ_FP3(v1P, v1);
}

TEST(VectorsProxyTest, FP3Equals)
{
	FP3 v1(13, 25, 21), v2(6, 8, 2);
	FP3Proxy v1P(v1), v2P(v2);
	
	ASSERT_FALSE(v1P == v2P);
	ASSERT_TRUE((v1P == v1P) && (v2P == v2P));
}

TEST(VectorsProxyTest, FP3NotEquals)
{
	FP3 v1(13, 25, 21), v2(6, 8, 2);
	FP3Proxy v1P(v1), v2P(v2);
	
	ASSERT_TRUE(v1P != v2P);
	ASSERT_FALSE((v1P != v1P) || (v2P != v2P));
}

TEST(VectorsProxyTest, FP3Lesser)
{
	FP3 v1(13, 25, 21), v2(6, 8, 2), v3(14, 6, 1);
	FP3Proxy v1P(v1), v2P(v2), v3P(v3);
	
	ASSERT_FALSE(v1P < v2P);
	ASSERT_TRUE(v2P < v1P);
	ASSERT_FALSE(v1P < v3P);
	ASSERT_FALSE(v3P < v1P);
	ASSERT_FALSE(v2P < v3P);
	ASSERT_FALSE(v3P < v2P);
}

TEST(VectorsProxyTest, FP3LesserOrEqual)
{
	FP3 v1(13, 25, 21), v2(6, 8, 2), v3(14, 25, 1);
	FP3Proxy v1P(v1), v2P(v2), v3P(v3);
	
	ASSERT_FALSE(v1P <= v2P);
	ASSERT_TRUE(v2P <= v1P);
	ASSERT_FALSE(v1P <= v3P);
	ASSERT_FALSE(v3P <= v1P);
	ASSERT_FALSE(v2P <= v3P);
	ASSERT_FALSE(v3P <= v2P);
}

TEST(VectorsProxyTest, FP3Greater)
{
	FP3 v1(13, 25, 21), v2(6, 8, 2), v3(14, 6, 1);
	FP3Proxy v1P(v1), v2P(v2), v3P(v3);
	
	ASSERT_TRUE(v1P > v2P);
	ASSERT_FALSE(v2P > v1P);
	ASSERT_FALSE(v1P > v3P);
	ASSERT_FALSE(v3P > v1P);
	ASSERT_FALSE(v2P > v3P);
	ASSERT_FALSE(v3P > v2P);
}

TEST(VectorsProxyTest, FP3GreaterOrEqual)
{
	FP3 v1(13, 25, 21), v2(6, 8, 2), v3(14, 8, 1);
	FP3Proxy v1P(v1), v2P(v2), v3P(v3);
	
	ASSERT_TRUE(v1P >= v2P);
	ASSERT_FALSE(v2P >= v1P);
	ASSERT_FALSE(v1P >= v3P);
	ASSERT_FALSE(v3P >= v1P);
	ASSERT_FALSE(v2P >= v3P);
	ASSERT_FALSE(v3P >= v2P);
}

TEST(VectorsProxyTest, FP3Floor)
{
	FP3 v(5.2, 7.5, -1.4);
	FP3Proxy vP(v);
	
	FP3 res = floor(vP);
	
	ASSERT_EQ_INT3(res, Int3(5, 7, -2));
}

TEST(VectorsProxyTest, FP3Truncate)
{
	FP3 v(5.2, 7.5, -1.4);
	FP3Proxy vP(v);
	
	FP3 res = truncate(vP);
	
	ASSERT_EQ_INT3(res, Int3(5, 7, -1));
}

TEST(VectorsProxyTest, FP3Dot)
{
	FP3 v1(5, 7, -1), v2(1, -2, 3);
	FP3Proxy v1P(v1), v2P(v2);
	
	ASSERT_EQ(dot(v1P, v2P), 5 * 1 + 7 * (-2) + (-1) * 3);
}

TEST(VectorsProxyTest, FP3Cross)
{
	FP3 v1(5, 7, -1), v2(1, -2, 3);
	FP3Proxy v1P(v1), v2P(v2);
	
	ASSERT_EQ_FP3(cross(v1P, v2P), FP3(19, -16, -17));
}

TEST(VectorsProxyTest, FP3Dist)
{
	FP3 v1(5, 7, -1), v2(3, 4, -7);
	FP3Proxy v1P(v1), v2P(v2);
	
	ASSERT_EQ(dist(v1P, v2P), 7);
}

TEST(VectorsProxyTest, FP3OperatorOutput)
{
	std::ostringstream os;
	FP3 v(2, -0.25, 1.5);
	FP3Proxy v1P(v);
	
	os << v1P;
	
	ASSERT_EQ("(2, -0.25, 1.5)", os.str());
}

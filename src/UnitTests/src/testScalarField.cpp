#include "TestingUtility.h"

#include "ScalarField.h"

template <class Data>
class ScalarFieldTest : public BaseFixture {
public:
    typedef Data data;
    typedef ScalarField<Data> ScalarFieldType;
protected:
    virtual void SetUp() {
        BaseFixture::SetUp();
        maxAbsoluteError = (FP)1e-4;
        maxRelativeError = (FP)0.001;
    }

    ScalarFieldType createScalarField(const Int3& size) {
        ScalarFieldType f(size);
        for (int i = 0; i < size.x; i++)
            for (int j = 0; j < size.y; j++)
                for (int k = 0; k < size.z; k++)
                    f(i, j, k) = k + (j + i * size.y) * size.z;
        return f;
    }
};

typedef ::testing::Types<FP> types;
TYPED_TEST_CASE(ScalarFieldTest, types);

TYPED_TEST(ScalarFieldTest, Constructor) {
    typedef typename ScalarFieldTest<TypeParam>::ScalarFieldType ScalarField;
    ScalarField f(Int3(3, 4, 5));
}

TYPED_TEST(ScalarFieldTest, CopyConstructor) {
    typedef typename ScalarFieldTest<TypeParam>::ScalarFieldType ScalarField;
    ScalarField f(Int3(3, 4, 5));
    ScalarField g(f);
    g(0, 0, 0) = 1.0;
}

TYPED_TEST(ScalarFieldTest, Assignment) {
    typedef typename ScalarFieldTest<TypeParam>::ScalarFieldType ScalarField;
    ScalarField f(Int3(3, 4, 5)), g(Int3(1, 3, 2));
    g = f;
    g(0, 0, 0) = 1.0;

}

TYPED_TEST(ScalarFieldTest, IndexAccess) {
    typedef typename ScalarFieldTest<TypeParam>::ScalarFieldType ScalarField;
    Int3 size(5, 3, 8);
    ScalarField f1(this->createScalarField(size));
    const ScalarField f2(this->createScalarField(size));
    for (int i = 0; i < size.x; i++)
        for (int j = 0; j < size.y; j++)
            for (int k = 0; k < size.z; k++) {
                ASSERT_EQ(f1(i, j, k), k + (j + i * size.y) * size.z);
                ASSERT_EQ(f1(i, j, k), f2(i, j, k));
                ASSERT_EQ(f1(i, j, k), f1(Int3(i, j, k)));
                ASSERT_EQ(f1(i, j, k), f2(Int3(i, j, k)));
            }
}

/*TYPED_TEST(ScalarFieldTest, Zeroize) {
    typedef typename ScalarFieldTest<TypeParam>::ScalarFieldType ScalarField;
    Int3 size(5, 3, 8);
    ScalarField f(this->createScalarField(size));
    f(0, 0, 0) = 1;
    f.zeroize();
    for (int i = 0; i < size.x; i++)
        for (int j = 0; j < size.y; j++)
            for (int k = 0; k < size.z; k++)
                ASSERT_EQ(f(i, j, k), 0);
}*/

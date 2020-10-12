#include "TestingUtility.h"

int main(int argc, char **argv)
{
    omp_set_num_threads(1);

    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();

    return result;
}
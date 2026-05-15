#include "mpi.h"
#include "gtest/gtest.h"
#include <iostream>
#include <sstream>

class MpiTestListener : public testing::TestEventListener 
{
public:
    explicit MpiTestListener(testing::TestEventListener* default_listener, int rank) 
        : default_listener_(default_listener), rank_(rank) {}
    
    virtual void OnTestProgramStart(const testing::UnitTest& unit_test) override 
    {
        std::cout << "[RANK " << rank_ << "] ";
        default_listener_->OnTestProgramStart(unit_test);
    }
    
    virtual void OnTestIterationStart(const testing::UnitTest& unit_test, int iteration) override 
    {
        std::cout << "[RANK " << rank_ << "] ";
        default_listener_->OnTestIterationStart(unit_test, iteration);
    }
    
    virtual void OnEnvironmentsSetUpStart(const testing::UnitTest& unit_test) override 
    {
        std::cout << "[RANK " << rank_ << "] ";
        default_listener_->OnEnvironmentsSetUpStart(unit_test);
    }
    
    virtual void OnEnvironmentsSetUpEnd(const testing::UnitTest& unit_test) override 
    {
        std::cout << "[RANK " << rank_ << "] ";
        default_listener_->OnEnvironmentsSetUpEnd(unit_test);
    }
    
    virtual void OnTestCaseStart(const testing::TestCase& test_case) override 
    {
        std::cout << "[RANK " << rank_ << "] ";
        default_listener_->OnTestCaseStart(test_case);
    }
    
    virtual void OnTestStart(const testing::TestInfo& test_info) override 
    {
        std::cout << "[RANK " << rank_ << "] ";
        default_listener_->OnTestStart(test_info);
    }
    
    virtual void OnTestPartResult(const testing::TestPartResult& result) override 
    {
        std::cout << "[RANK " << rank_ << "] ";
        default_listener_->OnTestPartResult(result);
    }
    
    virtual void OnTestEnd(const testing::TestInfo& test_info) override 
    {
        std::cout << "[RANK " << rank_ << "] ";
        default_listener_->OnTestEnd(test_info);
    }
    
    virtual void OnTestCaseEnd(const testing::TestCase& test_case) override 
    {
        std::cout << "[RANK " << rank_ << "] ";
        default_listener_->OnTestCaseEnd(test_case);
    }
    
    virtual void OnEnvironmentsTearDownStart(const testing::UnitTest& unit_test) override 
    {
        std::cout << "[RANK " << rank_ << "] ";
        default_listener_->OnEnvironmentsTearDownStart(unit_test);
    }
    
    virtual void OnEnvironmentsTearDownEnd(const testing::UnitTest& unit_test) override 
    {
        std::cout << "[RANK " << rank_ << "] ";
        default_listener_->OnEnvironmentsTearDownEnd(unit_test);
    }
    
    virtual void OnTestIterationEnd(const testing::UnitTest& unit_test, int iteration) override 
    {
        std::cout << "[RANK " << rank_ << "] ";
        default_listener_->OnTestIterationEnd(unit_test, iteration);
    }
    
    virtual void OnTestProgramEnd(const testing::UnitTest& unit_test) override 
    {
        std::cout << "[RANK " << rank_ << "] ";
        default_listener_->OnTestProgramEnd(unit_test);
    }
    
private:
    testing::TestEventListener* default_listener_;
    int rank_;
};

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    testing::InitGoogleTest(&argc, argv);
    
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    testing::TestEventListeners& listeners = testing::UnitTest::GetInstance()->listeners();
    auto* default_printer = listeners.Release(listeners.default_result_printer());
    
    listeners.Append(new MpiTestListener(default_printer, rank));
    
    int result = RUN_ALL_TESTS();

    MPI_Barrier(MPI_COMM_WORLD);
    for (int i = 0; i < size; i++)
    {
        if (rank == i)
            std::cout << " Result: " << ((result == 0) ? "[PASSED]" : "[FAILED]") << ", Num Failed: " << result << std::endl;
        MPI_Barrier(MPI_COMM_WORLD);
    }
    
    MPI_Finalize();

    return result;
}
#ifndef TEST_LINALG_HPP
#define TEST_LINALG_HPP

#include "TestCase.hpp"
#include "TestCaseContainer.hpp"

#include "types.hpp"

#include <limits>

namespace vaspml
{

Real const defaultTolerance = 100.0 * std::numeric_limits<Real>::epsilon();

struct TestCaseLinalg_l2Norm : public TestCase
{
    Real     tolerance;
    Vec1Real norm;
    Vec2Real vector;

    TestCaseLinalg_l2Norm(String name) : TestCase(name), tolerance(defaultTolerance) {}
};

template<>
void TestCaseContainer<TestCaseLinalg_l2Norm>::setup()
{
    TestCaseLinalg_l2Norm* tc = nullptr;
    testCases.push_back(TestCaseLinalg_l2Norm("TestCaseLinalg_l2Norm"));
    tc = &(testCases.back());

    tc->vector.resize(2);
    tc->norm.resize(2);
    tc->vector[0].resize(10);
    for (Size i = 0; i < tc->vector[0].size(); i++) { tc->vector[0][i] = (Real)i; }
    tc->vector[1].resize(15);
    for (Size i = 0; i < tc->vector[1].size(); i++) { tc->vector[1][i] = (Real)i * i; }
    tc->norm[0] = 16.881943016134134;
    tc->norm[1] = 357.3331778606627;
}

struct TestCaseLinalg_scaleVector : public TestCase
{
    Real     tolerance;
    Real     scalar;
    Vec1Real vector;
    Vec1Real target;

    TestCaseLinalg_scaleVector(String name) : TestCase(name), tolerance(defaultTolerance) {}
};

template<>
void TestCaseContainer<TestCaseLinalg_scaleVector>::setup()
{

    TestCaseLinalg_scaleVector* tc = nullptr;
    testCases.push_back(TestCaseLinalg_scaleVector("TestCaseLinalg_scaleVector"));
    tc = &(testCases.back());

    tc->scalar = 5.5;
    tc->vector.resize(20);
    tc->target.resize(20);
    for (Size i = 0; i < 20; i++)
    {
        tc->vector[i] = (Real)i * i;
        tc->target[i] = tc->vector[i] * tc->scalar;
    }
}

struct TestCaseLinalg_matMul : public TestCase
{
    Real tolerance;

    Vec1Real MatrixA;
    Vec1Real MatrixB;
    Vec1Real MatrixC;

    Int n;
    Int m;
    Int k;

    TestCaseLinalg_matMul(String name) : TestCase(name), tolerance(defaultTolerance) {}
};

template<>
void TestCaseContainer<TestCaseLinalg_matMul>::setup()
{

    TestCaseLinalg_matMul* tc = nullptr;
    testCases.push_back(TestCaseLinalg_matMul("TestCaseLinalg_matMul"));
    tc = &(testCases.back());

    tc->m = 5;
    tc->k = 3;
    tc->n = 4;

    tc->MatrixA.resize(15);
    tc->MatrixB.resize(12);

    tc->MatrixA[0] = 1;
    tc->MatrixA[1] = 2;
    tc->MatrixA[2] = 3;
    tc->MatrixA[3] = 4;
    tc->MatrixA[4] = 5;
    tc->MatrixA[5] = 6;
    tc->MatrixA[6] = 7;
    tc->MatrixA[7] = 8;
    tc->MatrixA[8] = 9;
    tc->MatrixA[9] = 10;
    tc->MatrixA[10] = 11;
    tc->MatrixA[11] = 12;
    tc->MatrixA[12] = 13;
    tc->MatrixA[13] = 14;
    tc->MatrixA[14] = 15;

    tc->MatrixB[0] = 1;
    tc->MatrixB[1] = 4;
    tc->MatrixB[2] = 7;
    tc->MatrixB[3] = 10;
    tc->MatrixB[4] = 2;
    tc->MatrixB[5] = 5;
    tc->MatrixB[6] = 8;
    tc->MatrixB[7] = 11;
    tc->MatrixB[8] = 3;
    tc->MatrixB[9] = 6;
    tc->MatrixB[10] = 9;
    tc->MatrixB[11] = 12;

    tc->MatrixC.resize(20);

    tc->MatrixC[0] = 14;
    tc->MatrixC[1] = 32;
    tc->MatrixC[2] = 50;
    tc->MatrixC[3] = 68;
    tc->MatrixC[4] = 32;
    tc->MatrixC[5] = 77;
    tc->MatrixC[6] = 122;
    tc->MatrixC[7] = 167;
    tc->MatrixC[8] = 50;
    tc->MatrixC[9] = 122;
    tc->MatrixC[10] = 194;
    tc->MatrixC[11] = 266;
    tc->MatrixC[12] = 68;
    tc->MatrixC[13] = 167;
    tc->MatrixC[14] = 266;
    tc->MatrixC[15] = 365;
    tc->MatrixC[16] = 86;
    tc->MatrixC[17] = 212;
    tc->MatrixC[18] = 338;
    tc->MatrixC[19] = 464;
}

struct TestCaseLinalg_scaleVectorPlusVector : public TestCase
{

    Real     tolerance;
    Real     scalar;
    Vec1Real scaleVector;
    Vec1Real addVector;
    Vec1Real target;

    TestCaseLinalg_scaleVectorPlusVector(String name) : TestCase(name), tolerance(defaultTolerance)
    {}
};

template<>
void TestCaseContainer<TestCaseLinalg_scaleVectorPlusVector>::setup()
{

    TestCaseLinalg_scaleVectorPlusVector* tc = nullptr;
    testCases.push_back(
        TestCaseLinalg_scaleVectorPlusVector("TestCaseLinalg_scaleVectorPlusVector"));
    tc = &(testCases.back());

    tc->scalar = (Real)3.5;
    tc->scaleVector.resize(5);
    tc->scaleVector[0] = (Real)3.5;
    tc->scaleVector[1] = (Real)3.5;
    tc->scaleVector[2] = (Real)3.5;
    tc->scaleVector[3] = (Real)3.5;
    tc->scaleVector[4] = (Real)3.5;

    tc->addVector.resize(5);
    tc->addVector[0] = (Real)2.0;
    tc->addVector[1] = (Real)2.0;
    tc->addVector[2] = (Real)2.0;
    tc->addVector[3] = (Real)2.0;
    tc->addVector[4] = (Real)2.0;

    tc->target.resize(5);
    tc->target[0] = (Real)14.25;
    tc->target[1] = (Real)14.25;
    tc->target[2] = (Real)14.25;
    tc->target[3] = (Real)14.25;
    tc->target[4] = (Real)14.25;
}

struct TestCaseLinalg_dotProduct : public TestCase
{

    Real     tolerance;
    Vec1Real vectorA;
    Vec1Real vectorB;
    Real     target;

    TestCaseLinalg_dotProduct(String name) : TestCase(name), tolerance(defaultTolerance) {}
};

template<>
void TestCaseContainer<TestCaseLinalg_dotProduct>::setup()
{

    TestCaseLinalg_dotProduct* tc = nullptr;
    testCases.push_back(TestCaseLinalg_dotProduct("TestCaseLinalg_dotProduct"));
    tc = &(testCases.back());

    tc->vectorA.resize(6);
    tc->vectorA[0] = (Real)1;
    tc->vectorA[1] = (Real)2;
    tc->vectorA[2] = (Real)3;
    tc->vectorA[3] = (Real)4;
    tc->vectorA[4] = (Real)5;
    tc->vectorA[5] = (Real)6;

    tc->vectorB.resize(6);
    tc->vectorB[0] = (Real)2;
    tc->vectorB[1] = (Real)4;
    tc->vectorB[2] = (Real)6;
    tc->vectorB[3] = (Real)8;
    tc->vectorB[4] = (Real)10;
    tc->vectorB[5] = (Real)12;

    tc->target = (Real)182;
}



struct TestCaseLinalg_solveLeastSquaresSVD : public TestCase
{

    Real     tolerance;
    Vec1Real matrix;
    Vec1Real b;
    Vec1Real x;
    Vec1Real target;
    Int m;
    Int n;

    TestCaseLinalg_solveLeastSquaresSVD(String name) : TestCase(name), tolerance(defaultTolerance)
    {}
};

template<>
void TestCaseContainer<TestCaseLinalg_solveLeastSquaresSVD>::setup()
{

    TestCaseLinalg_solveLeastSquaresSVD* tc = nullptr;
    testCases.push_back(TestCaseLinalg_solveLeastSquaresSVD("TestCaseLinalg_solveLeastSquaresSVD"));
    tc = &(testCases.back());
    Vec2Real A = {
                 {1, 2},
                 {3, 4},
                 {5, 6}
                 };
    tc->m = A.size();
    tc->n = A[0].size();
    for ( Int i = 0; i < tc->m; i++ ) {
        for ( Int j = 0; j < tc->n; j++ ) {
            tc->matrix.push_back( A[i][j] );
        }
    }
    tc->b = {7.0, 8.0, 9.0};
    tc->target = {-6.0, 6.5};
}


struct TestCaseLinalg_determinants : public TestCase
{

    Real     tolerance;
    Vec1Real matrixSquare;
    Vec1Real matrixRect;
    Vec1Real target;
    Int n0;
    Int n1;

    TestCaseLinalg_determinants(String name) : TestCase(name), tolerance(defaultTolerance)
    {}
};

template<>
void TestCaseContainer<TestCaseLinalg_determinants>::setup()
{

    TestCaseLinalg_determinants* tc = nullptr;
    testCases.push_back(TestCaseLinalg_determinants("TestCaseLinalg_determinants"));
    tc = &(testCases.back());
    tc->matrixSquare = {1,2,3,4};
    tc->matrixRect   = {1,2,3,4,5,6};
    tc->n0 = 2;
    tc->n1 = 3;
    tc->target.push_back( -2.0 );
    tc->target.push_back( 54.0 );
}

} //namespace vaspml
#endif
